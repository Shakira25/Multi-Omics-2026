<!-- Cloned 2026-11-12 from https://elsherbini.github.io/data-science-for-biology/materials/2-workshop2/5-enrichment/slides.html
Offline page: ../../../../website/materials/2-workshop2/5-enrichment/slides.html -->

# Sets and Enrichment

Workshop 2, Module 5

---

## Yesterday, preterm against term gave you this

1,698 genes with an adjusted p-value below 0.05.

This is the gene table tested against birth outcome, from the end of yesterday’s second exercise, not the anaerobe-against-*Lactobacillus* test with its 1,815 genes.

What is the finding?

---

## A gene list is hard to report

- You cannot read 1,698 names
- You cannot put them in an abstract
- Two groups with the same data will hand you two different lists

You need to say what the genes have **in common**.

---

## Cut-off tests and ranked tests

- **Cut, then count.** Take your significant genes and a predefined set, and ask whether the overlap is bigger than chance (a hypergeometric test). It is simple, and it ignores every gene below your cutoff.
- **Rank, then look.** Order every gene you tested from most up to most down, and ask whether a set’s members pile up at one end. This needs no cutoff. It is called gene set enrichment analysis (GSEA).

We use GSEA, because with the first approach the answer depends on a threshold you picked arbitrarily.

---

## A gene set is a named list

Someone curated it, from papers or experiments.

We use the **Hallmark** collection: 50 sets of about 200 genes, one for each broad biological program.

`HALLMARK_INFLAMMATORY_RESPONSE`
 `HALLMARK_TNFA_SIGNALING_VIA_NFKB`
 `HALLMARK_INTERFERON_GAMMA_RESPONSE`

The first 40 of the 195 `HALLMARK_TNFA_SIGNALING_VIA_NFKB` genes in our table:

PDE4B NFIL3 PTX3 CEBPB IL6 BTG2 FOS LITAF TRIB1 DUSP1 ATP2B1 NFAT5 BCL2A1 NAMPT EGR1 KDM6B PLAUR PTGS2 IL1B PLEK SLC16A6 SLC2A3 CCNL1 ETS2 CD44 SOCS3 SGK1 TLR2 PTPRE G0S2 TNFAIP3 SERPINE1 NFE2L2 BCL3 SNN DENND5A AREG IFNGR2 PFKFB3 NFKBIA

---

## Rank by the test statistic

```r
# res is your DESeq2 result table
stat <- res$stat
names(stat) <- rownames(res)
stat <- stat[!is.na(stat)]
```

- The p-value has no sign, so A and B look identical.
- The fold change alone puts noisy gene C at the top.
- The statistic is the fold change divided by its standard error, so it carries both size and confidence.
   Toy numbers

| gene | log2 FC | std. error | p-value | statistic |
|---|---|---|---|---|
| A | +3.0 | 0.5 | 2e-9 | +6.0 |
| B | -3.0 | 0.5 | 2e-9 | -6.0 |
| C | +6.0 | 4.0 | 0.13 | +1.5 |

---

## Walk down the list and keep score

Each tick is a set member. The score climbs at each hit and drifts down between them.

---

## The number you report

The highest point the score reaches is the enrichment score. Shuffle the gene labels many times to find out how surprising that peak is.

```r
library(fgsea)
library(msigdbr)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
sets <- split(hallmark$ensembl_gene, hallmark$gs_name)

fgsea(sets, stat, minSize = 15, maxSize = 500)
```

`NES` is the peak, normalised so sets of different sizes are comparable. Its sign tells you which end of the ranking the set sits at.

---

## Preterm versus term

---

## Reading the top sets

| set | size | NES | padj |
|---|---|---|---|
| Tnfa signaling via nfkb | 195 | 3.46 | 2.1e-39 |
| Inflammatory response | 184 | 3.32 | 1.2e-33 |
| Heme metabolism | 188 | -2.84 | 1.6e-20 |
| Interferon gamma response | 193 | 2.39 | 1.9e-12 |
| Complement | 185 | 2.35 | 6.4e-12 |
| Il6 jak stat3 signaling | 83 | 2.64 | 1.0e-10 |

Inflammatory signalling, tumour necrosis factor, interferon and complement are all shifted up in the people who delivered early.

You could put that sentence in an abstract, and it comes from the same 13,738 statistics that gave you the list of 1,698 genes.

---

## Quick check

---

## Now try the community contrast

Run the same analysis on yesterday’s other contrast, anaerobe versus Lactobacillus.

| set | size | NES | padj |
|---|---|---|---|
| Estrogen response early | 189 | -2.06 | 2.3e-08 |
| Tnfa signaling via nfkb | 195 | -1.87 | 4.7e-06 |
| Estrogen response late | 184 | -1.76 | 1.0e-04 |
| P53 pathway | 190 | -1.70 | 2.0e-04 |
| Hypoxia | 177 | -1.74 | 3.0e-04 |
| Il2 stat5 signaling | 190 | -1.64 | 6.5e-04 |

18 of 50 sets are significant, and 17 of those point the same way. Does that look like biology to you?

---

## What differs between these swabs

Aldunate et al., *Front. Physiol.* 2015, doi:10.3389/fphys.2015.00164

- The two communities sit on different epithelium, so their swabs pick up a different mix of cells. Bulk RNA-seq adds all of those cells together.
- Why does this make it hard to interpret the gene expression data?

---

## Checking that the statistics are centred

Is the distribution of test statistics centred?

| ranking | share of all genes positive |
|---|---|
| preterm vs term | 46.1% |
| anaerobe vs Lactobacillus | 55.0% |

Both are near half, so by this check both rankings look fine.

The check is wrong, because GSEA does not weight genes equally. It weights each gene by the size of its statistic.

---

## Look at the strong effects only

---

## Share positive at three levels

| ranking | all genes | significant genes | 500 strongest |
|---|---|---|---|
| preterm vs term | 46.1% | 49.4% | 49.6% |
| anaerobe vs Lactobacillus | 55.0% | 30.5% | 4.4% |

- Outcome: balanced at every level, so nothing pushes the sets one way.
- Community: balanced overall, but only 4.4% of the 500 strongest effects are positive.

Every large gene set will be dragged toward the Lactobacillus end.

---

## The sets follow the tail

| ranking | significant sets | pointing one way |
|---|---|---|
| preterm vs term | 26 | 23 |
| anaerobe vs Lactobacillus | 18 | 17 |

- Outcome: the strong effects are balanced, and the sets still agree with each other. Here the agreement comes from the biology.
- Community: the strong effects nearly all point one way, and so do the sets. The lopsided ranking alone would explain the direction.

---

## The individual genes still stand

The individual genes from the community contrast are still real. Mucins and filaggrin were up in Lactobacillus-dominant samples by a factor of eight, with adjusted p-values near ten to the minus eighteen.

The set names are harder to trust.

Go back to the question the experiment asked. If it was about immune activation, measured with RNA-seq on a mixed swab, tread carefully. The data is still telling you something about the biology. Deciding what it means is your job.

---

## Check the tail every time

Before you believe an enrichment, sort your genes by the size of their effect and look at the top few hundred.

If those nearly all point one way, every large set will lean that way. Ask what kind of genes they are before you believe the set names.

Our community ranking was centred overall and still failed this check.

---

## Limits of any gene set collection

- Sets overlap. TNFA signalling and inflammatory response share many genes, so their p-values are not independent and you cannot count them as separate findings.
- People curated them, mostly from experiments in cell lines and mostly not from cervical tissue.
- If nobody has curated a set for the biology you study, enrichment cannot find it.

---

## Quick check

---

## The same idea, on bacteria

A gene set is a label that groups features.

Species already carry labels that group them: **genus** and **family**.

---

## Rank species, test taxonomy

```r
# the same fgsea call; only the sets changed
genus_sets <- split(species_results$feature, species_results$Genus)
fgsea(genus_sets, species_stat, minSize = 2, maxSize = 50)
```

| genus | size | NES | padj |
|---|---|---|---|
| Gardnerella | 3 | 1.58 | 0.12 |
| Megasphaera | 2 | -1.26 | 0.62 |
| Mobiluncus | 2 | -1.30 | 0.62 |
| Sneathia | 2 | 1.18 | 0.80 |
| Actinomyces | 2 | -0.88 | 0.84 |
| Corynebacterium | 4 | -1.05 | 0.84 |
| Dialister | 2 | 0.82 | 0.84 |
| Porphyromonas | 2 | -1.02 | 0.84 |
| Prevotella | 12 | 0.79 | 0.84 |

*Gardnerella* leads, but only one of its three species ranks near the top, and at padj 0.12 it does not pass. With no species significant on its own, a set has little to gather.

---

## Set sizes on the two tables

|  | genes | species |
|---|---|---|
| features ranked | 13,738 | 53 |
| sets tested | 50 | 10 |
| typical set size | about 190 | 2 to 12 |
| sets significant | 26 | 0 |

Set-level analysis pays off when you have thousands of features and hundreds of curated groups. With 53 species in ten genera, you are better off reading the species list.

---

## Use the tree instead

microViz fits the same model to every taxon at every rank, from phyla near the centre to the 53 species at the edge.

- Colour: brown is higher in preterm, teal is higher in term
- Dot size: share of samples that carry the taxon

You can see a whole branch lean one way. *Lactobacillus* leans term and *Gardnerella* leans preterm, but only *Fenollaria* reaches p < 0.05, and that is before correcting for testing every node.

microViz: modelling taxa

---

## Choosing a method for your table

The code was identical on genes and on species, but it only told us something new about the genes.

Pick the method for how many features your table has and how they group, rather than for whether it came from RNA-seq or 16S.

---

## In the exercise

1. Ask Pi to rank your genes and run fgsea against Hallmark, then check that the identifiers matched
2. Let Pi interpret the top sets, then count how many genes they share
3. Check your ranking for imbalance before you believe the result
4. Compare against a cutoff-based test and a second set collection
5. Build taxonomic sets and run the same function on 53 species

---

## Your turn

Open `notebooks/05-enrichment.qmd`.

back to module
