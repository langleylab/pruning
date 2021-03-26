# Pruning pathway analysis results

[Tutorial](https://github.com/langleylab/pruning/blob/main/gsea_pruning.md) on how to prune FGSEA results according to the GO graph in R.

[Tutorial](https://github.com/langleylab/pruning/blob/main/gost_pruning.md) on how to prune g:Profiler GO ORA results according to the GO graph in R.

## TL;DR

FGSEA performed on GO gene sets overcorrects p values and reports non-specific terms because it treats every Gene Ontology category independently. In this tutorial we propose a simple approach to pruning GSEA results making use of the GO graph and show how terms can be "recovered" and give better biological insight.
Similarly, GO over-representation analyses performed by gprofiler2 shows repetitive terms, even thouth the FDR correction does respect the hierarchical structure. In a similar tutorial we propose the same approach to show specific terms.

## Small introduction

When doing differential gene expression analysis in RNA-seq experiments, it is often useful to perform **pathway enrichment analyses** on the results to have an overview of which biological processes and functions are being affected in the comparison between groups of samples. There are many methods to do so, and the two most popular approaches fall under two broad categories: 

- over-representation analyses (**ORA**)
- gene set enrichment analyses (**GSEA**)

There are several important differences between these two approaches. 

## ORA: testing for over-representation 

**ORA** tests whether a set of interesting genes (i.e. differentially expressed above an arbitrary threshold of the test statistic) is represented in pre-determined gene sets more than we would expect by chance, that is, the interesting genes are over-represented in the set. 
Expected values are calculated by taking into account the size of the interesting set, the size of the gene-sets, the size of their overlap and the size of the background set, that is the universe of genes under consideration whether they're differentially expressed or not, and ORA tests for the significance of this over-representation using a statistical test such as Pearson's χ^2 test, hypergeometric test or Fisher's exact test. 

The idea is simple: if we imagine a background set of 100 genes - the genes we are considering regardless of their status as interesting or not, i.e. our "universe" - and a set of 20 differentially expressed genes, then we have 20/100 = 20%. This means that if we pick genes at random from our background there is a 20% probability that they will be differentially expressed. This also applies to gene sets: if a gene set is made of 10 genes, and all its genes are included the bakcground, 20% of the genesets (meaning 2 genes) will be found by chance in this scenario. Let's consider the case in which, for this particular gene set, we have 4 differentially expressed genes overlapping with the gene set. This is 4/2 = 2 times as many genes as you would expect to find in the geneset by chance. This value of 2 is the over-representation value. How do we determine the statistical significance of this result?

It is useful to represent our 2x enrichment case in a 2x2 contingency table:

| Type | belongs to gene set | does not belong to gene set | total |
| :---: | :---: | :---: | :---: |
| **is DE** | 4 | 16 | 20 |
| **is not DE** | 6 | 74 | 80 |
| **total** | 10 | 90 | 100 |

The question we are asking is whether those 4 genes are found "by chance", i.e. whether the gene set is affected in the differential expression. We can slightly reframe the question as: what is the probability of finding an uneven distribution of DE genes within the gene set, under the null hypothesis that DE genes are found in the gene set by chance? We rewrite the table using symbols for simplicity:

| Type | belongs to gene set | does not belong to gene set | total |
| :---: | :---: | :---: | :---: |
| **is DE** | A | B | A + B |
| **is not DE** | C | D | C + D |
| **total** | A + C | B + D | A + B + C + D |

If we know the totals, A follows the hypergeometric distribution with A + B successes (DE genes) and C + D failures (non-DE genes). If we apply Fisher's exact test, which is one of the popular options for ORA, the probability of observing these values for A, B, C, D is:

<a href="https://www.codecogs.com/eqnedit.php?latex=p&space;=&space;\frac{(A&space;&plus;&space;B)!(A&space;&plus;&space;C)!(C&space;&plus;&space;D)!(B&space;&plus;&space;D)!}{A!&space;B!&space;C!&space;D!&space;(A&plus;B&plus;C&plus;D)!}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?p&space;=&space;\frac{(A&space;&plus;&space;B)!(A&space;&plus;&space;C)!(C&space;&plus;&space;D)!(B&space;&plus;&space;D)!}{A!&space;B!&space;C!&space;D!&space;(A&plus;B&plus;C&plus;D)!}" title="p = \frac{(A + B)!(A + C)!(C + D)!(B + D)!}{A! B! C! D! (A+B+C+D)!}" /></a>

Plugging in the numbers, we find that p = (20! * 10! * 80! * 90!)/(4! * 16! * 6! * 74! * 100!) = 0.0841073, meaning we have about an 8.4% probability of finding an overrepresentation as high as that by chance given these quantities. This is the p-value of the test for our observed over-representation. 

### ORA: choosing the background

The above equation highlights an important aspect of ORA, which is the choice of background. Since the background (A+B+C+D) quantity is included in the denominator, the bigger the background the lower the probability of observing that over-representation by chance, for any given set and all things being equal. This has important consequences for biological interpretation of results.

Suppose we treat fibroblasts with transforming growth factor beta (TGFB), a widely used model of fibrosis, and we use DMSO (in which TGFB is dissolve) as a control. We then perform RNA-seq and differential gene expression analysis. We expect the differentially expressed genes to show an over-representation in gene sets involved in collagen production and cell replication, whereas other gene sets (e.g. voltage dependent channels) should not be significantly over-represented, since TGFB does not regulate voltage dependent channels. What is a reasonable choice of background? We want to know what our differentially expressed genes are doing, so the discriminating feature - whether a gene goes on the first or the second row of the table above - is whether they are _differentially_ expressed. So we only consider _expressed_ genes as our background: non-expressed genes do not have a chance to be differentially expressed, therefore they would never appear in the contingency table. If we included all genes as our background, we would decrease the p value for all gene sets tested. This may lead to false positive results (or rather, in many cases, uninformative results). 
Consider the case in which you are testing several gene sets for the experiment below, for instance using the Gene Ontology sets. Using the whole set of genes as a background many significant results will pop up, mostly related to uninteresting aspects of biology, such as finding out that your fibroblasts are, indeeed, fibroblasts. Low p values would also appear for sets in which only one gene is differentially expressed, which is most likely due to chance.

Let's instead consider another kind of experiment: we gather blood samples from a cohort of patients with a new, unexplained congenital disease, and from their unaffected relatives. We perform whole exome sequencing on the DNA from these samples and we run a variant calling pipeline to identify mutations, and we find hundreds of genes affected. We want to know what pathways are being affected by these mutations, so we run an over-representation analysis on the mutated genes. What is in this case a reasonable choice of background? We no longer have the distinction between differentially expressed/unchanged, and we assume (somewhat problematically) that all genes can be mutated with the same probability. Moreover, we don't know a priori when an where these mutations affect development and/or function. In this case it makes sense to use the whole set of genes as a background.

### ORA: using the Gene Ontology

As mentioned above, we can run ORA using Gene Ontology sets. Each over-representation test is carried independently on each set, so the usual multiple-test correction strategies apply. However, GO sets are not independent: they exist in a hierarchical structure in which smaller, more specific sets - children - are contained within one (or more) broader sets - parents. Sets relate to each other in structures called Directed Acyclic Graphs (DAGs), in which sets are nodes and their hierarchical relation is determined by the edges going from more specific/smaller to more general/larger sets. GO has 3 DAGs: Biological Processes (BP), Molecular Function (MF), and Cellular Compartment (CC). So we cannot applying multiple test correction on the results from a GO ORA, as sets are not independent, and we are not really testing whether "limb development" is over-represented as opposed to the more specific "arm development". Moreover, the multiple test correction procedures will lead us to reject terms just based on the sheer number of sets which affects the amount of correction, resulting in many false negeatives.

Several ORA tools deal with this fact by simply removing terms that have enriched children, so that only the most specific, significant terms are left. Other strategies can also be used where the hierarchical filtering is not as strong, and weights are given to each gene set based on their position on the DAG. Whatever the strategy, the final result is that of a more interpretable list of pathways, in which the number of false negatives is greatly reduced.

In the tutorial we provide a solution to filter the results of a popular ORA tool, 'gprofiler2', which does not currently perform this kind of pruning.


## GSEA 

A different approach altogether, **Gene Set Enrichment Analysis** uses the whole set of genes and their associated test statistics (e.g. log fold change in a differential expression analysis) to identify interesting patterns in the data. The most frequently used GSEA strategy is a pre-ranked GSEA, where genes are ranked according to the test statistic. In this strategy, the algorithm asks whether there is a significant enrichment of members of a gene set in the extremes (top or bottom) of the rank. Going back to the previous case of TGFB-stimulated fibroblasts, we can expect an up-regulation of many genes related to collagen secretion. If we rank genes by their fold change, we expect to find these genes roughly in the top part of the rank. Hedgehog-related genes, instead, may be scattered across the list, sitting mostly in the central part (where fold changes approximate zero). 

Greatly simplifying, GSEA uses a permutation-based test in which gene labels are shuffled across ranks several times. This creates a null distribution against which the actual rank distribution can be tested for every gene set. The number of permutations in which genes in the gene set can be found in a similar position in the shuffled ranks determines the (empirical) p value. Say the first 10 DE genes all belong to collagen secretion. How many times do collagen secretion genes appear in the top 10 positions if we randomize the gene names 10000 times? Most likely, very rarely (although the number of permutations gives us a lower bound on p values). 

So GSEA does not require to select a background, and does not require to threshold genes based on a statistic. Moreover, considering that enrichment can be either on the top or bottom of the ranks, it also gives a nice directionality to the result. 

### GSEA: using the GO gene sets

A natural source of gene sets for GSEA is, once again, Gene Ontology. However, we have the same issues as before: nominal p values will be corrected for results on the whole DAG, and we will get many non-specific significant terms. However, unlike ORA tools, GSEA tools do not have approaches to deal with the DAG structure of Gene Ontology. In the tutorial we show a possible approach and discuss its limitations.
