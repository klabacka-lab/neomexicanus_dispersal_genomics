import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pca = pd.read_csv(snakemake.input.eigenvec, sep=r"\s+", header=None)

pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca.shape[1]-1)]

eig = np.loadtxt(snakemake.input.eigenval)

eig = eig[eig > 0]

pct_var = eig / eig.sum() * 100

plt.figure(figsize=(8,6))

plt.scatter(pca["PC1"], pca["PC2"])

for _, row in pca.iterrows():
    plt.annotate(str(row["IID"]), (row["PC1"], row["PC2"]))

plt.xlabel(f"PC1 ({pct_var[0]:.1f}%)")
plt.ylabel(f"PC2 ({pct_var[1]:.1f}%)")
plt.title(f"Merged PCA")

plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=300)
plt.close()
