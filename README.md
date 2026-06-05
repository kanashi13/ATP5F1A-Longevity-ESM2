# ATP5F1A-Longevity-ESM2

**Exploring lifespan-associated evolutionary patterns in ATP synthase subunit ALPHA using ESM2 protein language model.**

---

### Part of the **ATP Synthase Longevity Series**

This project investigates whether evolutionary variations in **ATP5F1A** (the alpha subunit of mitochondrial ATP synthase) correlate with species lifespan across mammals.

## Project Overview

**ATP5F1A** is a core component of the mitochondrial ATP synthase (Complex V), responsible for ATP production. Due to its critical role, it is expected to be highly conserved.

### Species Analyzed (11 Mammals)

- **Homo sapiens** (Human) — 80 years
- **Pan troglodytes** (Chimpanzee) — 60 years
- **Loxodonta africana** (African savanna elephant) — 65 years
- **Heterocephalus glaber** (Naked mole-rat) — 37 years
- **Macaca mulatta** (Rhesus macaque) — 40 years
- **Equus caballus** (Horse) — 40 years
- **Myotis lucifugus** (Little brown bat) — 34 years *(used instead of Brandt's bat due to sequence availability)*
- **Canis lupus familiaris** (Dog) — 9 years
- **Mesocricetus auratus** (Golden hamster) — 4 years
- **Rattus norvegicus** (Norway rat) — 4 years
- **Mus musculus** (House mouse) — 4 years

## Methods

- Protein sequence retrieval from NCBI/UniProt
- Multiple Sequence Alignment using **Clustal Omega**
- Protein embeddings generated with **ESM2** (`esm2_t33_650M_UR50D`)
- Analysis: PCA, Cosine Distance, Hierarchical Clustering
- Variable sites identification

## Results

- **High Conservation**: Only \~35 variable positions out of 553 amino acids (\~6.3% variability).
- **No Strong Lifespan Signal**: Sequence variations are primarily **phylogenetic** (mainly distinguishing rodents from other mammals) rather than associated with longevity.
- Long-lived species (Human, Elephant, Naked mole-rat) do not share distinctive amino acid signatures in ATP5F1A.
- Clustering patterns (PCA, Heatmap, Dendrogram) reflect evolutionary relationships more than lifespan categories.

**Conclusion**:  
**ATP5F1A does not exhibit significant evolutionary signals related to species lifespan.** The protein is under strong purifying selection, consistent with its essential function in energy production.

---

## Repository Structure

- `FASTA/` — Individual and combined protein sequences
- `figures/` — PCA plot, Cosine Distance Heatmap, Dendrogram
- `results/` — Variable positions table and final conclusions
- `analysis.py` — Main analysis script
- `variable_positions.xlsx` — All variable sites

## How to Reproduce

```bash
pip install -r requirements.txt
python analysis.py
