# Bioinformatics-cheatsheet
A cheat sheet for Bioinformatians.

## Biological Concepts
- **DNA**: Deoxyribonucleic acid is a molecule that carries the genetic instructions used in the growth, development, functioning and reproduction of all known living organisms and many viruses. [@Wiki](https://en.wikipedia.org/wiki/DNA), [@NIH](https://www.genome.gov/25520880/deoxyribonucleic-acid-dna-fact-sheet/)
- **RNA**: Ribonucleic acid (RNA) is a polymeric molecule implicated in various biological roles in coding, decoding, regulation, and expression of genes. [@Wiki](https://en.wikipedia.org/wiki/RNA)
- **Protein**: Proteins are large biomolecules, or macromolecules, consisting of one or more long chains of amino acid residues. [@Wiki](https://en.wikipedia.org/wiki/Protein)
- **Gene**: A gene is a locus (or region) of DNA which is made up of nucleotides and is the molecular unit of heredity. [@Wiki](https://en.wikipedia.org/wiki/Gene), [@NIH](https://ghr.nlm.nih.gov/primer/basics/gene)
- **SNP**: A single nucleotide polymorphism, often abbreviated to SNP (pronounced snip; plural snips), is a variation in a single nucleotide that occurs at a specific position in the genome, where each variation is present to some appreciable degree within a population (e.g. >1%). [@Wiki](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism), [@NIH](https://ghr.nlm.nih.gov/primer/genomicresearch/snp)
- **Pathway (Biological pathway)**: A biological pathway is a series of actions among molecules in a cell that leads to a certain product or a change in a cell. Such a pathway can trigger the assembly of new molecules, such as a fat or protein. Pathways can also turn genes on and off, or spur a cell to move. [@Wiki](https://en.wikipedia.org/wiki/Biological_pathway), [@NIH](https://www.genome.gov/27530687/biological-pathways-fact-sheet/)
- **CNV**: The chromosome now has two copies of this section of DNA, rather than one. Copy number variation (CNVs) is a relatively new field in genomics and it is defined as a phenomenon in which sections of the genome are repeated and the number of repeats in the genome varies between individuals in the human population. [@Wiki](https://en.wikipedia.org/wiki/Copy-number_variation)
- **Transcription**: Transcription is the first step of gene expression, in which a particular segment of DNA is copied into RNA (mRNA) by the enzyme RNA polymerase. [@Wiki](https://en.wikipedia.org/wiki/Transcription_(genetics))
- **Translation**: In molecular biology and genetics, translation is the process in which cellular ribosomes create proteins. In translation, messenger RNA (mRNA)—produced by transcription from DNA—is decoded by a ribosome to produce a specific amino acid chain, or polypeptide. [@Wiki](https://en.wikipedia.org/wiki/Translation_(biology))
- **Transcription factor**: In molecular biology and genetics, a transcription factor (sometimes called a sequence-specific DNA-binding factor) is a protein that binds to specific DNA sequences, thereby controlling the rate of transcription of genetic information from DNA to messenger RNA. [@Wiki](https://en.wikipedia.org/wiki/Transcription_factor), [@BroadInstitute](https://www.broadinstitute.org/education/glossary/transcription-factor), [@Scitable](http://www.nature.com/scitable/definition/general-transcription-factor-transcription-factor-167)
- **Expression (Gene expression)**: Gene expression is the process by which information from a gene is used in the synthesis of a functional gene product. These products are often proteins, but in non-protein coding genes such as transfer RNA (tRNA) or small nuclear RNA (snRNA) genes, the product is a functional RNA. [@Wiki](https://en.wikipedia.org/wiki/Gene_expression), [@Scitable](http://www.nature.com/scitable/topicpage/gene-expression-14121669)

## Experimental Techniques
### Next Generate Sequencing
- **DNA-seq**: DNA sequencing is the process of determining the precise order of nucleotides within a DNA molecule. [@Wiki](https://en.wikipedia.org/wiki/DNA_sequencing)
- **RNA-seq**: RNA-seq (RNA sequencing), also called whole transcriptome shotgun sequencing[1] (WTSS), uses next-generation sequencing (NGS) to reveal the presence and quantity of RNA in a biological sample at a given moment in time. [@Wiki](https://en.wikipedia.org/wiki/RNA-Seq)
- **ChIP-seq**: ChIP-sequencing, also known as ChIP-seq, is a method used to analyze protein interactions with DNA. ChIP-seq combines chromatin immunoprecipitation (ChIP) with massively parallel DNA sequencing to identify the binding sites of DNA-associated proteins. [@Wiki](https://en.wikipedia.org/wiki/ChIP-sequencing)
- **CLIP-Seq**: [@Wiki](https://en.wikipedia.org/wiki/CLIP)
- **FAIRE-seq**: [@Wiki](https://en.wikipedia.org/wiki/FAIRE-Seq)
- **DNase-Seq**: [@Wiki](https://en.wikipedia.org/wiki/DNase-Seq)
- **CAGE**: [@Wiki](https://en.wikipedia.org/wiki/Cap_analysis_gene_expression)
- **ChIA-PET**: [@Wiki](https://en.wikipedia.org/wiki/Chia_Pet)
- **5C/Hi-C**: [@Wiki](https://en.wikipedia.org/wiki/Chromosome_conformation_capture) 
### Uncategorized
- **ChIP**: Chromatin Immunoprecipitation (ChIP) is a type of immunoprecipitation experimental technique used to investigate the interaction between proteins and DNA in the cell. [@Wiki](https://en.wikipedia.org/wiki/Chromatin_immunoprecipitation)
- **ChIP-chip**: ChIP-chip (also known as ChIP-on-chip) is a technology that combines chromatin immunoprecipitation ('ChIP') with DNA microarray ("chip"). Like regular ChIP, ChIP-on-chip is used to investigate interactions between proteins and DNA in vivo. [@Wiki](https://en.wikipedia.org/wiki/ChIP-on-chip)

## NGS data repositories
- **[GEO](http://www.ncbi.nlm.nih.gov/geo/)**: GEO is a public functional genomics data repository supporting MIAME-compliant data submissions. Array- and sequence-based data are accepted.
- **[ENCODE](https://www.encodeproject.org/)**: The goal of ENCODE is to build a comprehensive parts list of functional elements in the human genome, including elements that act at the protein and RNA levels, and regulatory elements that control cells and circumstances in which a gene is active.
- **[TCGA](http://cancergenome.nih.gov/)**: The Cancer Genome Atlas (TCGA) is a collaboration between the National Cancer Institute (NCI) and the National Human Genome Research Institute (NHGRI) that has generated comprehensive, multi-dimensional maps of the key genomic changes in 33 types of cancer. 
- **[1000 Genomes Project](http://www.1000genomes.org/)**: The 1000 Genomes Project ran between 2008 and 2015, creating the largest public catalogue of human variation and genotype data.
- **[Array Express](http://www.ebi.ac.uk/arrayexpress/)**: an NIH-funded database at the European Molecular Biology Laboratory -European Bioinformatics Institute that collects and disseminates microarray-based gene-expression data.
- **[DDBJ](http://www.ddbj.nig.ac.jp/)**: DNA Data Bank of Japan (DDBJ) is a data bank organized by the National Institute of Genetics in Japan that collects sequence data.

## NGS data analysis
### Variants calling
#### SNP calling
- **[GATK](https://software.broadinstitute.org/gatk/)**: Genome Analysis Toolkit offers a wide variety of tools with a primary focus on variant discovery and genotyping.

## File formats
### Formats
- **[SAM](http://samtools.github.io/hts-specs/SAMv1.pdf)**: The SAM Format is a text format for storing sequence data in a series of tab delimited ASCII columns. [@Wiki](http://genome.sph.umich.edu/wiki/SAM)
- **BAM**: BAM is the compressed binary version of the Sequence Alignment/Map (SAM) format, a compact and index-able representation of nucleotide sequence alignments.
- **BED**: The BED format consists of one line per feature, each containing 3-12 columns of data, plus optional track definition lines. [@UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), [@bedtools](http://bedtools.readthedocs.io/en/latest/content/general-usage.html), [@Ensembl](http://useast.ensembl.org/info/website/upload/bed.html)
- **WIG**: The WIG (wiggle) format is designed for display of dense continuous data such as probability scores. Wiggle data elements must be equally sized; if you need to display continuous data that is sparse or contains elements of varying size, use the BedGraph format instead. [@UCSC](https://genome.ucsc.edu/goldenpath/help/wiggle.html), [@Ensembl](http://useast.ensembl.org/info/website/upload/wig.html)
- **BigWig**: The BigWig format is designed for dense, continuous data that is intended to be displayed as a graph. Files can be created from WIG or BedGraph files using the appropriate utility program. [@UCSC](http://genome.ucsc.edu/goldenPath/help/bigWig.html)
- **GFF**: The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines. [@UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format3), [@Ensembl](http://useast.ensembl.org/info/website/upload/gff.html), [GTF(GFFv2)@GMOD](http://gmod.org/wiki/GFF2), [@Wiki](https://en.wikipedia.org/wiki/General_feature_format)
- **VCF**: The Variant Call Format (VCF) specifies the format of a text file used in bioinformatics for storing gene sequence variations. [v4.0@1000genomes](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40), [@Wiki](https://en.wikipedia.org/wiki/Variant_Call_Format)
### Tools
- **[bedtools](https://github.com/arq5x/bedtools2)**: Collectively, the bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks.
- **[samtools](https://github.com/samtools/samtools)**: SAM Tools provide various utilities for manipulating alignments in the SAM/BAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- **[sambamba](https://github.com/lomereiter/sambamba)**: samtools functions with multi-threading support.
- **[BEDOPS](https://bedops.readthedocs.io/en/latest/)**: the fast, highly scalable and easily-parallelizable genome analysis toolkit.
- **[bwtool](https://github.com/CRG-Barcelona/bwtool)**: bwtool is a command-line utility for bigWig files.
- **[VCFtools](https://vcftools.github.io/index.html)**: A set of tools written in Perl and C++ for working with VCF files.


