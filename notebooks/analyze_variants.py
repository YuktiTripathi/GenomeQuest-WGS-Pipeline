#!/usr/bin/env python3
"""
GenomeQuest WGS Pipeline — Variant Analysis
Analyzes annotated variants from E. coli LTEE experiment (SRR2584866)
aligned to REL606 ancestor reference genome.
"""

import cyvcf2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
import csv
import os

os.makedirs('notebooks/figures', exist_ok=True)

# ── Parse annotated SNPs ──
snp_path = 'annotation/annotated_snps_rel606.vcf'
indel_path = 'annotation/annotated_indels_rel606.vcf'

effects, impacts, genes, positions = [], [], [], []

for variant in cyvcf2.VCF(snp_path):
    if variant.FILTER is not None:
        continue
    ann_raw = variant.INFO.get('ANN', None)
    if ann_raw:
        first_ann = str(ann_raw).split(',')[0]
        fields = first_ann.split('|')
        effects.append(fields[1].strip() if len(fields) > 1 else 'unknown')
        impacts.append(fields[2].strip() if len(fields) > 2 else 'unknown')
        genes.append(fields[3].strip() if len(fields) > 3 else 'unknown')
        positions.append(variant.POS)

# ── Figure 1: Effect types + Impact distribution ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('E. coli LTEE SRR2584866 — Variant Effects vs REL606 Ancestor',
             fontsize=13, fontweight='bold')

effect_counts = Counter(effects).most_common(8)
axes[0].barh([e[0].replace('_', ' ') for e in effect_counts],
             [e[1] for e in effect_counts], color='steelblue')
axes[0].set_xlabel('Count')
axes[0].set_title('Variant Effect Types (SNPs)')
axes[0].invert_yaxis()

impact_counts = Counter(impacts)
colors = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e',
          'LOW': '#2ca02c', 'MODIFIER': '#1f77b4'}
labels = list(impact_counts.keys())
axes[1].bar(labels, [impact_counts[l] for l in labels],
            color=[colors.get(l, 'gray') for l in labels])
axes[1].set_ylabel('Count')
axes[1].set_title('Variant Impact Distribution')
legend = [mpatches.Patch(color=v, label=k) for k, v in colors.items()]
axes[1].legend(handles=legend, fontsize=9)

plt.tight_layout()
plt.savefig('notebooks/figures/variant_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: notebooks/figures/variant_summary.png")

# ── Figure 2: Top mutated genes ──
gene_counts = Counter(genes).most_common(15)
fig, ax = plt.subplots(figsize=(10, 6))
ax.barh([g[0] for g in gene_counts], [g[1] for g in gene_counts], color='coral')
ax.set_xlabel('Number of SNPs')
ax.set_title('Top 15 Most Mutated Genes\nE. coli LTEE SRR2584866 vs REL606')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('notebooks/figures/top_genes.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: notebooks/figures/top_genes.png")

# ── Figure 3: Genomic position of variants ──
fig, ax = plt.subplots(figsize=(12, 3))
high_pos = [p for p, i in zip(positions, impacts) if i == 'HIGH']
mod_pos  = [p for p, i in zip(positions, impacts) if i == 'MODERATE']
low_pos  = [p for p, i in zip(positions, impacts) if i in ('LOW', 'MODIFIER')]

ax.scatter(low_pos,  [1]*len(low_pos),  c='#2ca02c', s=5, alpha=0.4, label='LOW/MODIFIER')
ax.scatter(mod_pos,  [1]*len(mod_pos),  c='#ff7f0e', s=8, alpha=0.7, label='MODERATE')
ax.scatter(high_pos, [1]*len(high_pos), c='#d62728', s=12, alpha=1.0, label='HIGH')
ax.set_xlabel('Genomic Position (bp)')
ax.set_title('Variant Distribution Across REL606 Genome (4.63 Mb)')
ax.set_xlim(0, 4629812)
ax.set_yticks([])
ax.legend(loc='upper right')
plt.tight_layout()
plt.savefig('notebooks/figures/genomic_distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: notebooks/figures/genomic_distribution.png")

# ── Export summary CSV ──
with open('annotation/variant_summary.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['position', 'gene', 'effect', 'impact'])
    for pos, gene, eff, imp in zip(positions, genes, effects, impacts):
        writer.writerow([pos, gene, eff, imp])
print("Saved: annotation/variant_summary.csv")

# ── Print summary ──
print(f"\n{'='*45}")
print(f"  VARIANT ANALYSIS SUMMARY")
print(f"{'='*45}")
print(f"  Total PASS SNPs:       {len(effects)}")
print(f"  HIGH impact:           {impacts.count('HIGH')}")
print(f"  MODERATE impact:       {impacts.count('MODERATE')}")
print(f"  LOW impact:            {impacts.count('LOW')}")
print(f"  MODIFIER:              {impacts.count('MODIFIER')}")
print(f"\n  Top 10 mutated genes:")
for gene, count in Counter(genes).most_common(10):
    print(f"    {gene:<30} {count} SNPs")
print(f"{'='*45}")
