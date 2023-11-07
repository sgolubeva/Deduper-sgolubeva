#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp               #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB  
outfile=/projects/bgmp/sgolubev/bioinfo/Bi624/Deduper-sgolubeva/sorted_C1_SE_uniqAlign.sam 
deduped_file=/projects/bgmp/sgolubev/bioinfo/Bi624/Deduper-sgolubeva/C1_SE_dedup.sam
summary_file=/projects/bgmp/sgolubev/bioinfo/Bi624/Deduper-sgolubeva/C1_SE_summary.tsv
duplicate_file=/projects/bgmp/sgolubev/bioinfo/Bi624/Deduper-sgolubeva/C1_SE_duplicates.sam

/usr/bin/time -v python golubeva_deduper.py -u STL96.txt -o $outfile -d $deduped_file -s $summary_file -p $duplicate_file