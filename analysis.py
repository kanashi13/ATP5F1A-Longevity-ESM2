import esm
import torch
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import dendrogram, linkage
import os

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()

fasta_folder = "FASTA"
species_data = []

for filename in os.listdir(fasta_folder):
    if filename.endswith(".fasta"):
        filepath = os.path.join(fasta_folder, filename)
        record = SeqIO.read(filepath, "fasta")
        sequence = str(record.seq)
        species_name = filename.replace(".fasta", "")
        species_data.append({"species": species_name, "sequence": sequence})

embeddings = []
names = []

for item in species_data:
    seq = item["sequence"]
    name = item["species"]
    tokens = alphabet.encode(seq)
    tokens = torch.tensor([tokens])
    with torch.no_grad():
        results = model(tokens, repr_layers=[33])
    token_repr = results["representations"][33][0][1:-1]
    avg_embedding = token_repr.mean(0).numpy()
    embeddings.append(avg_embedding)
    names.append(name)
    print(f"{name} done (len={len(seq)})")

X = np.array(embeddings)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X)

lifespan_map = {
    "Mus musculus": ("Low", "red"),
    "Rattus norvegicus": ("Low", "red"),
    "Mesocricetus auratus": ("Low", "red"),
    "Canis lupus familiaris": ("Medium", "blue"),
    "Homo sapiens": ("High", "green"),
    "Pan troglodytes": ("Medium", "blue"),
    "Macaca mulatta": ("Medium", "blue"),
    "Equus caballus": ("Medium", "blue"),
    "Heterocephalus glaber": ("Medium", "blue"),
    "Loxodonta africana": ("High", "green"),
    "Myotis lucifugus": ("Medium", "blue")
}

colors = []
for name in names:
    if name in lifespan_map:
        _, c = lifespan_map[name]
        colors.append(c)
    else:
        colors.append("gray")

plt.figure(figsize=(12, 10))
for i, name in enumerate(names):
    plt.scatter(pca_result[i, 0], pca_result[i, 1], c=colors[i], s=120, edgecolors='k', alpha=0.8)
    plt.annotate(name, (pca_result[i, 0], pca_result[i, 1]), fontsize=9, ha='center', va='bottom')

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', edgecolor='k', label='Low (<10 yr)'),
    Patch(facecolor='blue', edgecolor='k', label='Medium (10-60 yr)'),
    Patch(facecolor='green', edgecolor='k', label='High (>60 yr)')
]
plt.legend(handles=legend_elements)
plt.title("ATP5F1A Protein Space (ESM2 + PCA)")
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("ATP5A_PCA.png", dpi=150)
plt.show()

cos_sim = cosine_similarity(X)
cos_dist = 1 - cos_sim

dist_df = pd.DataFrame(cos_dist, index=names, columns=names)

plt.figure(figsize=(10, 8))
sns.heatmap(dist_df, annot=True, fmt=".3f", cmap="YlOrRd", linewidths=0.5)
plt.title("Cosine Distance Matrix of ATP5F1A Embeddings")
plt.tight_layout()
plt.savefig("ATP5A_heatmap.png", dpi=150)
plt.show()

Z = linkage(cos_dist, method='average')

plt.figure(figsize=(10, 5))
dendrogram(Z, labels=names, leaf_rotation=45, leaf_font_size=10)
plt.title("Hierarchical Clustering of ATP5F1A Embeddings (Cosine Distance)")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig("ATP5A_dendrogram.png", dpi=150)
plt.show()