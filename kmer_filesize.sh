#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2048M

output="logs/avg_kmer_filesize.csv"

> "$output"

for kdir in xeng_data/k*/; do
    kname=$(basename "$kdir")

    avg=$(
        find "$kdir" -type f \
            \( -name "*-host.fq.gz" -o -name "*-graft.fq.gz" \) \
            -exec stat -c%s {} + |
        awk '{sum+=$1; count++} END {if(count>0) print (sum/count)/1024/1024/1024; else print 0}'
    )

    echo "$kname,$avg" >> "$output"
done
