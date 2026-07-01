from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

host_fasta = snakemake.input.host
graft_fasta = snakemake.input.graft

output_fasta = snakemake.output[0]

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
