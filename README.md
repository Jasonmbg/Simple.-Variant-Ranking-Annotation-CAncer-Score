# Scoring initial theme

### Efstathios-Iason Vlachavas
#### DKFZ-Department of Molecular Genome Analysis
##### Efstathios-Iason.Vlachavas@dkfz-heidelberg.de

## Description

The main goal of the developing a ranking scheme in translational cancer research, is to aid the biological interpretation of lists of annotated cancer variants at the single patient resolution. Briefly, taking into account different annotation resources, such as variant effect prediction, cancer evidence, expression and systems biology properties, an integrated scoring value is assigned to each variant as a *hollistic* score, ranking variants from a single patient mutations list. In addition, the first version of the scoring process is based on single nucleotide changes that occur in the protein-coding space; that is, somatic mutations that can lead to several possible changes to a protein. Overall, the main rationale is despite the fact that a significant number of variant annotation tools in cancer research are available for a robust exploitation of putative somatic variants, no direct score or assessment is available for a simple ranking or prioritization of the interrogated variants, based on the plethora of distinct evidence. 

The ranking score is based on the output of the OpenCRAVAT annotation platform (https://doi.org/10.1200/cci.19.00132) and can serve as an additional *plug-in*, aiding in the interpretation of the annotated variant calling results. The open-source OpenCRAVAT toolkit possesses significant and important features, by making it possible to integrate multiple sources of evidence for variant annotation and exploitation.

In addition, for the construction and the relative score cut-offs, various guidelines and publications were considered for estimating the oncogenicity of somatic mutations, based on big consortia, such as:

1. [Belgian ComPerMed Expert Panel](https://doi.org/10.3390/cancers11122030) 

2. [Variant Interpretation for Cancer Consortium (VICC)](https://cancervariants.org/research/standards/onc_path_sop/) 


## Implementation

1. The user has to initially run OpenCRAVAT web server (https://run.opencravat.org/) or install locally (https://open-cravat.readthedocs.io/en/latest/quickstart.html). The input can be a vcf file, or a txt with necessary columns (https://open-cravat.readthedocs.io/en/latest/File-Formats.html)

2. The following annotators should run: gnomAD, ClinVar, CIViC, Mutation Assessor, FATHMM XF Coding, VEST4, SpliceAI, COSMIC, CScape Coding, Cancer Gene Census,
Cancer Gene Landscape, Cancer Hotspots, SiPhy and Phast Cons (15 annotators if having hg19 as the reference genome, to also include **hg19 coordinates**).

3. Next, an RData file has to be created either from the download section of the web server, or locally using the installed version of OpenCRAVAT

```python
oc report example_input.sqlite -t rdata

```

For more details see [here](https://open-cravat.readthedocs.io/en/latest/Reporter.html#example)

4. Finally, after creating the necessary RData file including the variants from one patient/sample, the main function in R to run:

```r
ranked_snvs = predict_sysSVM2(rdata_dir,exp.genes=NULL)

```

## Dependencies

```r
install.packages(c("DT","tidyverse","jsonlite"))

```

## Reproducible example

Here we present a simple example using the mutations from a randomly selected colorectal cancer patient sample ("crc4") from published Reiter et al., 2018 study (https://doi.org/10.1126/science.aat7171), mainly utilizing the Kim et al., 2015 publication (https://clincancerres.aacrjournals.org/content/21/19/4461#:~:text=10.1158/1078-0432.CCR-14-2413). Then, the web version of OpenCRAVAT was used to perform integrative variant annotation using the 15 aforementioned annotators, and the relative RData file was created. Below, a snapshot of the created html file with the top 10 hits are depicted:

![Alt text](relative/path/to/img.jpg?raw=true "Top 10 ranked variants example")

## Utilization feedback

For any questions, suggestions or issues please directly use my email or the github issue page 

## Acknowledgements

Stefan Wiemann

Kym Pagel

Rick Kym

Olga Papadodima
