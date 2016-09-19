TOPMed Freeze 3a Variant Calling Pipeline
=========================================

Overview of this repository
----------------------------

This repository is intended to provide a copy of software tools used for producing TOPMed Year 1 Freeze 3a variant calls and genotypes with a comprehensive documentation that allows investigators to understand the methods and reproduce the variant calls from the same set of aligned sequence reads.

This repository reflects specific versions of software tools that are under active development in the Center for Statistical Genetics (CSG). Most of the latest version of these software tools can be accessed at http://github.com/atks/vt or http://github.com/hyunminkang/apigenome , and this repository is focused on a freeze of software tools that can reproduce a variant calls compatible to TOPMed Freeze 3a.


Outline of the variant calling procedure
----------------------------------------

Our ``GotCloud vt`` pipeline detects and genotype variants from a list of aligned sequence reads. Specifically, the pipeline consist of the following six key steps. Most of these procedure will be integrated into ``GotCloud`` software package later this year. 

1. Variant detection : For each sequenced genome (in BAM/CRAMs), candidate variants are detected by ``vt discover2`` software tools, separated by each chromosome. The candidate variants are normalized by ``vt normalize`` algorithm. 
2. Variant consolidation : For each chromosome, the called variant sites are merged across the genomes, accounting for overlap of variants between genomes, using ``vt consolidate`` software tool.
3. Genotype and feature collection : For each 100kb chunk of genome, the genotyping module implemented in ``vt joint_genotype_sequential`` collects individual genotypes and variant features across the merged sites by iterating each sequence genome focusing on the selected region.  
4. Identifying related and duplicate subjects : Using ``king`` software and customized scripts included in this repository, we infer pedigree of nuclear families containing related and duplicated samples.
5. Variant filtering : We use the inferred pedigree of related and duplicated samples to calculate the Mendlian consistency statistics using ``vt milk-filter``, and train variant classifier using Support Vector Machine (SVM) implemented in the ``libsvm`` software package.
6. Variant annotation and post processing : We use ``snpEff`` software tools to annotate each variants and produce a filtered version of variant calls using a customized script.
