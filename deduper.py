#!/usr/bin/env python3
"""The purpose of this script is to remove PCR duplicates from single end sequencing sam file
based on strandeness, UNIs, leftmost position adjastment from CIGAR string"""
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="global variables to set")
    parser.add_argument("-u", "--umis", help="provide a file with umis", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="absolute path to sorted bam", required=True, type=str)
    parser.add_argument("-d", "--deduped", help="path to duplicates removed file", required=False, type=str)
    #parser.add_argument("-h", "--help")
    return parser.parse_args()


def create_umi_set(umis) -> set:

    """Takes file containing UMI list, adds UMIs into a set, returns a set of UMIs"""
    
    umi_set = set()
    with open(umis, 'r') as fh:
        for umi in fh:
            umi = umi.strip('\n')
            umi_set.add(umi)
    return umi_set


def check_strand(flag:int) -> bool:

    """Takes a bitwise flag, determines if strand is reverse or forward returns True if reverse,
    False if forward"""

    return((flag & 16) == 16)


def parse_cigar(cigar_str:str) -> list:

    """Takes cigar string and parses it into tuple of a character and a number that goes with this
      character. Returns a list of tuples: first element is number, second is character"""
    
    parsed = []
    parse = re.findall('\d+[M|N|D|S]', cigar_str)
    for piece in parse:
        match = re.match('(\d+)([M|N|D|S])', piece)
        result = (int(match.group(1)), match.group(2))
        parsed.append(result)
    return parsed

def adjust_position(parsed_cigar: list, direction: bool, lm_position) -> int:
    
    """Takes parsed cigar string, strand direction and leftmost position. Adjusts the leftmost position
    and returns it"""

    if direction == True:
        if parsed_cigar[0][1] == 'S':
            parsed_cigar.pop(0)
        for piece in parsed_cigar:
            lm_position += piece[0]
        #lm_position = lm_position - 1
    else:
        if parsed_cigar[0][1] == 'S':
            lm_position = lm_position - parsed_cigar[0][0]

    return lm_position


if __name__ == "__main__":
    args = get_args()
    umis: str = args.umis # hold a link to the umi file
    sam_f: str = args.outfile # holds path to a sorted sam file
    dedup_f: str = args.deduped # holds path for writing into duplicate removed file
    track_uniq_dict: dict = {} # holds unique chr, leftmost position, strand direction
    seen_chr: str  = '' # save the number of chromosome currently processed
    umi_set = create_umi_set(umis) # creates a set of UMIs 

    with open(sam_f, 'r') as sm, open(dedup_f, 'w') as dp:
        for line in sm:
            if line.startswith('@'):
                dp.write(line)
            else:
                spl_line = line.split('\t')
                umi = spl_line[0].split(':')[-1]
                
                # Chek if a UMI is in UMI set, if it is, proceed with read comparison
                if umi in umi_set: 
                    chrom = spl_line[2]
                    bwf = int(spl_line[1])
                    lm_pos = int(spl_line[3])
                    cigar = spl_line[5]
                    strand_dir = check_strand(bwf)
                    parsed_cigar = parse_cigar(cigar)

                    # adjast left most position of the read
                    adj_pos = adjust_position(parsed_cigar, strand_dir, lm_pos)
                    if chrom != seen_chr:
                        seen_chr = chrom
                        track_uniq_dict.clear()
                        
                        
                    
                    if (umi, chrom, strand_dir, adj_pos) in track_uniq_dict:
                        track_uniq_dict[(umi,chrom, strand_dir, adj_pos)]+=1

                    else:
                        track_uniq_dict[(umi, chrom, strand_dir, adj_pos)] = 1
                        dp.write(line)

                
