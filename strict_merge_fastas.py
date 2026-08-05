import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <host_fasta> <graft_fasta> <output_fasta>")
        sys.exit(1)

    host_fasta = sys.argv[1]
    graft_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    host = SeqIO.to_dict(SeqIO.parse(host_fasta, "fasta"))
    graft = SeqIO.to_dict(SeqIO.parse(graft_fasta, "fasta"))

    if set(host.keys()) != set(graft.keys()):
        missing_host = set(graft.keys()) - set(host.keys())
        missing_graft = set(host.keys()) - set(graft.keys())

        raise ValueError(
            f"Host and graft contain different sample IDs.\n"
            f"Only in graft: {sorted(missing_host)}\n"
            f"Only in host: {sorted(missing_graft)}"
        )

    combined_records = []

    for sample in sorted(host.keys()):
        combined_records.append(
            SeqRecord(
                host[sample].seq + graft[sample].seq,
                id=sample,
                description=""
            )
        )

    SeqIO.write(combined_records, output_fasta, "fasta")

    print("fastas merged")


if __name__ == "__main__":
    main()