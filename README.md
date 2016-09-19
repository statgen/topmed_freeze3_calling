TOPMed Freeze 3a Variant Calling Pipeline
=========================================

Overview of this repository
----------------------------

This repository is intended to provide a copy of software tools used for producing TOPMed Year 1 Freeze 3a variant calls and genotypes with a comprehensive documentation that allows investigators to understand the methods and reproduce the variant calls from the same set of aligned sequence reads.

This repository reflects specific versions of software tools that are under active development in the Center for Statistical Genetics (CSG). Most of the latest version of these software tools can be accessed at http://github.com/atks/vt or http://github.com/hyunminkang/apigenome , and this repository is focused on a freeze of software tools that can reproduce a variant calls compatible to TOPMed Freeze 3a.


Outline of the variant calling procedure
----------------------------------------

Our ``GotCloud vt`` pipeline detects and genotype variants from a list of aligned sequence reads. Specifically, the pipeline consist of the following six key steps. Most of these procedure will be integrated into ``GotCloud`` software package later this year. 

1. **Variant detection** : For each sequenced genome (in BAM/CRAMs), candidate variants are detected by ``vt discover2`` software tools, separated by each chromosome. The candidate variants are normalized by ``vt normalize`` algorithm. 
2. **Variant consolidation** : For each chromosome, the called variant sites are merged across the genomes, accounting for overlap of variants between genomes, using ``vt consolidate`` software tool.
3. **Genotype and feature collection** : For each 100kb chunk of genome, the genotyping module implemented in ``vt joint_genotype_sequential`` collects individual genotypes and variant features across the merged sites by iterating each sequence genome focusing on the selected region.  
4. **Identifying related and duplicate subjects** : Using ``king`` software and customized scripts included in this repository, we infer pedigree of nuclear families containing related and duplicated samples.
5. **Variant filtering** : We use the inferred pedigree of related and duplicated samples to calculate the Mendlian consistency statistics using ``vt milk-filter``, and train variant classifier using Support Vector Machine (SVM) implemented in the ``libsvm`` software package.
6. **Variant annotation and post-processing** : We use ``snpEff`` software tools to annotate each variants and produce a filtered version of variant calls using a customized script.


Variant Detection
-----------------
Variant detection from each sequence (ang aligned) genome is performed by ``vt discover2`` software tool. The script ``step-1-detect-variant.pl`` provide a mean to automate the variant detection across a large number of sequence genome.

The variant detection algorithm consider a variant as a potential candidate variant if there exists a mismatch between the aligned sequence reads and the reference genome. Because such a mismatch can easily occur by random errors, only potential candidate variants passing the following criteria are considered to be ***candidate variants*** in the next steps.

1. At least two identical evidence of variants must be observed from aligned sequence reads. 
  1. Each individual evidence will be normalized using the normalization algorithm implemented in ``vt normalize`` software tools.
  1. Only evidence on the reads with mapping quality 20 or greater will be considered.
  1. Duplicate reads, QC-passed reads, supplementary reads, secondary reads will be ignored. 
  1. Evidence of variant within overlapping fragments of read pairs will not be double counted. Either end of the overlapping read pair will be soft-clipped using ``bam clipOverlap`` software tool.  
1. Assuming per-sample heterozygosity of 0.1%, the posterior probability of having variant at the position should be greater than 50%. The method is equivalent to the `glfSingle` model described in http://www.ncbi.nlm.nih.gov/pubmed/25884587

 

