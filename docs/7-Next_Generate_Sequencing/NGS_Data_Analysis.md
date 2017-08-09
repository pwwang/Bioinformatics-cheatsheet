# NGS Data Analysis
<!-- toc -->

## Read Simulation
- **[dwgsim](https://github.com/nh13/DWGSIM)**: Whole genome simulator.
- **[ReadSim](https://sourceforge.net/p/readsim/wiki/Home/)**: ReadSim is a fast and simple reads simulator to target long reads such as PacBio or Nanopore.
- **[simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/)**: simNGS is software for simulating observations from Illumina sequencing machines using the statistical models behind the AYB base-calling software.
- **[wgsim](https://github.com/lh3/wgsim)**: Wgsim is a small tool for simulating sequence reads from a reference genome.

## Read Trimming
- **[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)**: A flexible read trimming tool for Illumina NGS data.
- **[Sickle](https://github.com/najoshi/sickle)**: A windowed adaptive trimming tool for FASTQ files using quality.
- **[famas](https://github.com/andreas-wilm/famas)**: Yet another program for FastQ massaging with features: Quality- and length-based trimming, Random sampling, Splitting into multiple files, Order checking for paired-end files, Native gzip support.

## De-Duplication
- **[PiCard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)**: This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
- **[sambamba-markdup](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html)**: Find duplicate reads in BAM file.

## Alignment
- **[bwa](http://bio-bwa.sourceforge.net/)**: BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.
- **[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**: Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
- **[ABRA](https://github.com/mozack/abra)**: ABRA is a realigner for next generation sequencing data. It uses localized assembly and global realignment to align reads more accurately, thus improving downstream analysis (detection of indels and complex variants in particular). [@Ref](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4173014/)
- **[NextGenMap](http://cibiv.github.io/NextGenMap/)**: NextGenMap (NGM) is a flexible and fast read mapping program that is more than twice as fast as BWA, while achieving a mapping sensitivity similar to Stampy or Bowtie2. [@Ref](http://bioinformatics.oxfordjournals.org/content/29/21/2790.long)

## Quality Control
- **[ClinQC](https://sourceforge.net/projects/clinqc/)**: ClinQC is an integrated and user-friendly pipeline for quality control, filtering and trimming of Sanger and NGS sequencing data for hundred to thousands of samples/patients in a single run in clinical research.
- **[NGS QC Toolkit](http://www.nipgr.res.in/ngsqctoolkit.html)**: NGS QC Toolkit: A toolkit for the quality control (QC) of next generation sequencing (NGS) data.

## Peak Calling/Differential Peak Calling
__Peak calling__ is a computational method used to identify areas in a genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing or MeDIP-seq experiment. These areas are those where a protein interacts with DNA.  
__Differential peak calling__ is about identifying significant differences in two ChIP-seq signals.
- **[MACS](http://liulab.dfci.harvard.edu/MACS/)**: Model-based analysis of ChIP-seq (MACS) is a computational algorithm that identifies genome-wide locations of transcription/chromatin factor binding or histone modification from ChIP-seq data. 
- **[DBChIP](http://pages.cs.wisc.edu/~kliang/DBChIP/)**: detects differentially bound sharp binding sites across multiple conditions, with or without matching control samples.
- **[MAnorm](http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm.htm)**: a robust model for quantitative comparison of ChIP-Seq data sets.
- **[THOR](http://www.regulatory-genomics.org/thor-2/basic-intrstruction/)**: Differential peak calling of ChIP-seq signals with replicates. [@Ref](http://nar.oxfordjournals.org/content/early/2016/08/01/nar.gkw680.abstract)
- **[ODIN](http://www.regulatory-genomics.org/odin-2/basic-introduction/)**: ODIN is an HMM-based approach to detect and analyse differential peaks in pairs of ChIP-seq data. ODIN performs genomic signal processing, peak calling and p-value calculation in an integrated framework. [@Ref](http://bioinformatics.oxfordjournals.org/content/30/24/3467)
- **[MMDiff](https://www.bioconductor.org/packages/release/bioc/html/MMDiff.html)**: This package detects statistically significant difference between read enrichment profiles in different ChIP-Seq samples. [@Ref](bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-826)

## RNA-seq data analyses
- **[RSEM](http://deweylab.github.io/RSEM/)**: accurate transcript quantification from RNA-Seq data with or without a reference genome[@Ref](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)
- **[Limma](https://bioconductor.org/packages/release/bioc/html/limma.html)**: Linear Models for Microarray and RNA-Seq Data.
- **[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)**: Empirical Analysis of Digital Gene Expression Data in R.
- **[DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)**: Differential gene expression analysis based on the negative binomial distribution.
- **[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)**: Transcriptome assembly and differential expression analysis for RNA-Seq.
- **[MISO](https://miso.readthedocs.io/en/fastmiso/)**: MISO (Mixture-of-Isoforms) is a probabilistic framework that quantitates the expression level of alternatively spliced genes from RNA-Seq data, and identifies differentially regulated isoforms or exons across samples. [@Ref](http://www.nature.com/nmeth/journal/v7/n12/full/nmeth.1528.html)

## (Capture) Hi-C data analyses
- **[CHiCAGO](http://regulatorygenomicsgroup.org/chicago)**: CHiCAGO is a set of tools for calling significant interactions in Capture HiC data, such as Promoter Capture HiC. [@Ref](https://www.ncbi.nlm.nih.gov/pubmed/27306882)

## Chromatin status data analysis
- **[CENTIPEDE](http://centipede.uchicago.edu/)**: CENTIPEDE applies a hierarchical Bayesian mixture model to infer regions of the genome that are bound by particular transcription factors. [@Ref](http://genome.cshlp.org/content/21/3/447)

## Variant Calling
- **[FaSD](http://jjwanglab.org/fasd)**: a fast and accurate single-nucleotide polymorphism detection program that uses a binomial distribution-based algorithm and a mutation probability.
- **[SOAPsnp](http://soap.genomics.org.cn/soapsnp.html)**: SOAPsnp uses a method based on Bayes’ theorem (the reverse probability model) to call consensus genotype by carefully considering the data quality, alignment, and recurring experimental errors.
- **[SNVmix](http://compbio.bccrc.ca/software/snvmix/)**: SNVMix is designed to detect single nucleotide variants from next generation sequencing data. 
- **[CNVnator](https://github.com/abyzovlab/CNVnator)**: a tool for CNV discovery and genotyping from depth-of-coverage by mapped reads.
- **[bcftools](https://samtools.github.io/bcftools/bcftools.html)**: utilities for variant calling and manipulating VCFs and BCFs.
- **[GATK](https://software.broadinstitute.org/gatk/)**: Genome Analysis Toolkit offers a wide variety of tools with a primary focus on variant discovery and genotyping.
- **[Bambino](https://cgwb.nci.nih.gov/goldenPath/bamview/documentation/index.html)**: a variant detector and alignment viewer for next-generation sequencing data in the SAM/BAM format.
- **[CONSERTING](http://www.stjuderesearch.org/site/lab/zhang)**: integrating copy-number analysis with structural-variation detection. [@Ref](http://www.nature.com/nmeth/journal/v12/n6/full/nmeth.3394.html)
- **[CREST](http://www.stjuderesearch.org/site/lab/zhang)**: CREST (Clipping Reveals Structure) is a new algorithm for detecting genomic structural variations at base-pair resolution using next-generation sequencing data. [@Ref](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.1628.html)
- **[Control-FREEC](http://boevalab.com/FREEC/)**: Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH) using deep-sequencing data.
- **[HMMcopy](http://compbio.bccrc.ca/software/hmmcopy/)**: Copy number prediction with correction for GC and mappability bias for HTS data
- **[SegSeq](http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=182)**: an algorithm to identify chromosomal breakpoints using massively parallel sequence data.
- **[CNV-seq](http://tiger.dbs.nus.edu.sg/cnv-seq/)**: a new method to detect copy number variation using high-throughput sequencing.
- **[BICseq2](http://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.html)**: BICseq2 is an algorithm developed for the normalization of  high-throughput sequencing (HTS) data and detection of copy number variations (CNV) in the genome. BICseq2 can be used for detecting CNVs with or without a control genome.
- **[MuSE](http://bioinformatics.mdanderson.org/main/MuSE)**: a novel approach to mutation calling based on the F81 Markov substitution model for molecular evolution, which models the evolution of the reference allele to the allelic composition of the matched tumor and normal tissue at each genomic locus. [@Ref](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1029-6)
- **[VarScan](http://dkoboldt.github.io/varscan/)**: a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data.
- **[Pindel](https://github.com/genome/pindel)**: Pindel can detect breakpoints of large deletions, medium sized insertions, inversions, tandem duplications and other structural variants at single-based resolution from next-gen sequence data.
- **COPS**: A Sensitive and Accurate Tool for Detecting Somatic Copy Number Alterations Using Short-Read Sequence Data from Paired Samples. COPS is available at ftp://115.119.160.213 with username “cops” and password “cops”. [@Ref](http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0047812)
- **[multiSNV](https://bitbucket.org/joseph07/multisnv/wiki/Home)**: multiSNV is a tool for calling somatic single-nucleotide variants (SNVs) using NGS data from a normal and multiple tumour samples of the same patient. [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4482059/)
- **[SomaticSeq](http://bioinform.github.io/somaticseq/)**: SomaticSeq is a flexible post-somatic-mutation-calling workflow for improved accuracy. [@Ref](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574535/)
- **[Cosmos](http://seselab.org/cosmos/)**: COSMOS can detect somatic structural variations from whole genome short-read sequences. [@Ref](http://nar.oxfordjournals.org/content/44/8/e78)
- **[Platypus](http://www.well.ox.ac.uk/platypus)**: Platypus is a tool designed for efficient and accurate variant-detection in high-throughput sequencing data.

## Variant Filtering
- **[SnpSift](http://snpeff.sourceforge.net/SnpSift.html)**: SnpSift is a toolbox that allows you to filter and manipulate annotated files.
- **[Varapp](https://varapp-demo.vital-it.ch/#/login)**: Varapp is an open-source web application to filter variants from large sets of exome data stored in a relational database. 

## Variant Annotators
- **[ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)**: is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes (including human genome hg18, hg19, hg38, as well as mouse, worm, fly, yeast and many others). 
- **[SnpEff](http://snpeff.sourceforge.net/)**: Genetic variant annotation and effect prediction toolbox. [@Ref](http://www.tandfonline.com/doi/abs/10.4161/fly.19695)
- **[GEMINI](https://gemini.readthedocs.io/en/latest/)**: GEMINI (GEnome MINIng) is a flexible framework for exploring genetic variation in the context of the wealth of genome annotations available for the human genome. [@Ref](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003153)
- **[Variant Effect Predictor](http://asia.ensembl.org/Homo_sapiens/Tools/VEP?db=core)**: Analyse your own variants and predict the functional consequences of known and unknown variants via our Variant Effect Predictor (VEP) tool. [@Ref](http://bioinformatics.oxfordjournals.org/content/26/16/2069)
- **[VAT - Variant Annotation Tool](http://vat.gersteinlab.org/index.php)**: A computational framework to functionally annotate variants in personal genomes using a cloud-computing environment. [@Ref](http://bioinformatics.oxfordjournals.org/content/28/17/2267)
- **[SeattleSeq Variation Annotation](http://snp.gs.washington.edu/SeattleSeqAnnotation147/)**: The SeattleSeq Annotation server provides annotation of SNVs (single-nucleotide variations) and small indels, both known and novel. [@Ref](http://www.nature.com/nature/journal/v461/n7261/full/nature08250.html)
- **[Jannovar](https://github.com/charite/jannovar)**: A Java Library for Exome Annotation.
- **[Cellbase](https://github.com/opencb/cellbase)**: CellBase is a scalable and high-performance NoSQL database that integrates relevant biological information from well-known data sources such as Ensembl, Uniprot, IntAct or ClinVar among others. All this data can be queried through a comprehensive RESTful web services API or using the command line interface. Also, a built-in variant annotator has been developed and can be used to annotate files containing variants in Variant Call Format (VCF). [@Ref](http://nar.oxfordjournals.org/content/40/W1/W609.short)
- **[GenomeD3Plot](https://github.com/brinkmanlab/GenomeD3Plot/)**: GenomeD3Plot (formerly [Islandplot](https://github.com/lairdm/islandplot)) is an SVG based genome viewer written in javascript using D3.
- **[Variant Tools](http://varianttools.sourceforge.net/Main/HomePage)**: variant tools is a software tool for the manipulation, annotation, selection, simulation, and analysis of variants in the context of next-gen sequencing analysis. [@Ref](http://bioinformatics.oxfordjournals.org/content/28/3/421.abstract?sid=f64403e7-5050-4102-963c-e690efe003f7) 
- **[GWAVA](https://www.sanger.ac.uk/sanger/StatGen_Gwava)**: GWAVA is a tool which aims to predict the functional impact of non-coding genetic variants based on a wide range of annotations of non-coding elements (largely from ENCODE/GENCODE), along with genome-wide properties such as evolutionary conservation and GC-content.

## Variant prioritization
- **[Mutsig](http://archive.broadinstitute.org/cancer/cga/mutsig)**: MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes. [@Ref](http://www.nature.com/nature/journal/v499/n7457/full/nature12213.html)

## Variant Simulation
- **[SCNVSim](https://sourceforge.net/projects/scnvsim/)**: somatic copy number variation and structure variation simulator. [@Ref](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0502-7)

## Haplotype Estimation Tools
- **[PHASE](http://stephenslab.uchicago.edu/phase/download.html)**: A program for reconstructing haplotypes from population data. PHASE was limited by its speed and was not applicable to datasets from genome-wide association studies. [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1275651/) 
- **[fastPHASE](http://scheet.org/code/fastphase_doc_1.4.pdf)**: fastPHASE is software that implements methods for estimating missing genotypes and reconstructing haplotypes from unphased SNP genotype data of unrelated individuals.  [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1424677/) 
- **[IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)**: IMPUTE version 2 (also known as IMPUTE2) is a genotype imputation and haplotype phasing program based on ideas from [Howie et al. 2009](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000529)

## NGS Data/Variant/Genome Visulizers/Browsers/Diagrams
- **[UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks)**: The UCSC Genome Browser is an on-line genome browser hosted by the University of California, Santa Cruz (UCSC).
- **[JBrowse](http://jbrowse.org/)**: a JavaScript genome browser by the open-source Generic Model Organism Database project. [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2752129)
- **[Synthesis-View](http://visualization.ritchielab.psu.edu/synthesis_views/plot)**: A SNP visualization tool. [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3012023/)
- **[IGV - Integrative Genomics Viewer](http://www.broadinstitute.org/igv/)**: A high-performance visualization tool for interactive exploration of large, integrated genomic datasets. [@Ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3603213)
- **[pileup.js](http://www.hammerlab.org/2015/06/19/introducing-pileup-js-a-browser-based-genome-viewer/)**: a Browser-based Genome Viewer.
- **[Biodalliance](http://www.biodalliance.org/index.html)**: Biodalliance is a fast, interactive, genome visualization tool that's easy to embed in web pages and applications.
- **[TADkit](https://github.com/3DGenomes/TADkit)**: TADkit is a HTML5 and JavaScript-based 3D genome browser. It makes use of D3.js for rendering the 1D and 2D tracks and WebGl by Three.js for rendering the 3D track. 
- **[ngs.plot](https://github.com/shenlab-sinai/ngsplot)**: ngs.plot is a program that allows you to easily visualize your next-generation sequencing (NGS) samples at functional genomic regions. [@Ref](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-284)
- **[CHiCP](https://www.chicp.org/chicp/)**: a web-based tool for the integrative and interactive visualization of promoter capture Hi-C datasets. [@Ref](http://bioinformatics.oxfordjournals.org/content/early/2016/04/26/bioinformatics.btw173.full)
- **[WashU Epigenome Browser](http://epigenomegateway.wustl.edu/browser/)**: An Epigenome browser.

## NGS Data Analysis Pipeline/framework
- **[nextflow](https://www.nextflow.io/)**: Nextflow is a fluent DSL modelled around the UNIX pipe concept, that simplifies writing parallel and scalable pipelines in a portable manner.
- **[RUbioSeq+](http://rubioseq.bioinfo.cnio.es/)**: RUbioSeq+ is a stand-alone and multiplatform application for the integrated analysis of NGS data. More specifically, our software implements pipelines for the analysis of single nucleotide and copy-number variation, bisulfite-seq and ChIP-seq experiments using well-established tools to perform these common tasks.
- **[SpeedSeq](https://github.com/hall-lab/speedseq)**: A flexible framework for rapid genome analysis and interpretation. [@Ref](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3505.html)
- **[HICUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/)**: pipeline for mapping and processing Hi-C data. [@Ref](https://www.ncbi.nlm.nih.gov/pubmed/26835000)
- **[OpEx](http://www.icr.ac.uk/our-research/research-divisions/division-of-genetics-and-epidemiology/genetic-susceptibility/genetic-data-and-software-resources/the-opex-ngs-pipeline)**: Provides a fixed implementation of alignment, calling and annotation tools optimized for individual or multiple exome sequencing analysis in the research or clinical setting. [@Ref](http://www.nature.com/articles/srep31029)
