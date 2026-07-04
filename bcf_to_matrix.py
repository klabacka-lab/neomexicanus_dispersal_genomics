import csv
import sys

infile = snakemake.input.tsv
outfile = snakemake.output.matrix

with open(infile) as f, open(outfile, "w", newline="") as out:

    reader = csv.reader(f, delimiter="\t")
    writer = csv.writer(out)

    header = next(reader)

    new_header =["POS"] + header[2:]
    writer.writerow(new_header)

    for row in reader:
        gts = row[2:]
        snp = f"{row[0]}:{row[1]}"

        out_row =[snp]

        all_zero = True
        all_one = True

        for gt in gts:
            if gt == "./." or gt == ".":
                val = 0
            elif gt.startswith("0/0"):
                val = 0
            else:
                val = 1

            out_row.append(val)

            if val == 1:
                all_zero = False
            else:
                all_one = False

        if all_zero or all_one:
            continue
    
        writer.writerow(out_row)
