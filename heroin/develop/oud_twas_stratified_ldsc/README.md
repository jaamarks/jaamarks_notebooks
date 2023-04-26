**Author**: Jesse Marks (jmarks@rti.org)<br>
**NIH Project**: [Harnessing Knowledge of Gene Function in Brain Tissue for Discovering Biology Underlying Heroin Addiction](https://reporter.nih.gov/search/RC99reuHhEW0n_3WuFPU6g/project-details/10116351) <br>
**Charge Code**: 0218755.001.002.001<br>
**GitHub Issue**: [Opioid Use Disorder TWAS Meta-analysis (Uniform Processing) #183](https://github.com/RTIInternational/bioinformatics/issues/183)<br>


<br>


# Description
These notebooks detail the process used to test for heritability enrichment in and around differentially expressed genes from our Opioid Use Disorder TWAS meta-analysis. 

We perform Stratified LD Score Regression (S-LDSC) analyses using the LD score regression approach described in [Finucane et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5896795/) on specifically expressed genes (LDSC-SEG).
LDSC-SEG applies stratified LD score regression to test whether disease heritability is enriched in the regions surrounding genes with the highest specific expression in a given tissue.
This approach helps to interpret GWAS signal by leveraging gene expression data.

We used gene sets from a TWAS meta-analysis of 4 OUD case/control datasets:
- [Corradin et al. (2022) Molecular Psychiatry](https://doi.org/10.1038/s41380-022-01477-y)
- [Mendez et al. (2021) Molecular Psychiatry](https://doi.org/10.1038/s41380-021-01259-y)
- [Seney et al. (2021) Biological Psychiatry](https://doi.org/10.1016/j.biopsych.2021.06.007)
- [Sosnowski et al. (2022) Drug and Alcohol Dependence Reports](https://doi.org/10.1016/j.dadr.2022.100040)


We filtered the results to obtain two sets of significantly expressed genes:

- List of genes with a Benjamini-Hochberg FDR < 0.05
- List of genes with a Benjamini-Hochberg FDR < 0.10

The `wfisher_adj_pvalue` column was used to create these gene lists.
This column contains the [weighted Fisher's p-value](https://www.nature.com/articles/s41598-021-86465-y) that we applied the Benjamini-Hochberg FDR thresholds to.

Our aim was to determine whether genes (and their proximal genomic regions) that show differential expression by opioid use disorder are also enriched for genetic signal associated with opioid addiction and related phenotypes.


___

<br><br>

<details>
    <summary>phenotype list</summary>
    
    
___    
    
* Age of Initiation  (Liu et al., 2019 Nat Genet [30643251](https://pubmed.ncbi.nlm.nih.gov/30643251/))
* Alcohol Dependence (Walters et al., 2018 Nat Neurosci [30482948](https://pubmed.ncbi.nlm.nih.gov/30482948))
* Alcohol Drinks per Week (DPW) (Liu et al., 2019 Nat Genet [30643251]())
* Alzheimer's Disease (Lambert et al., 2013 Nat Genet [24162737](https://pubmed.ncbi.nlm.nih.gov/24162737))
* Amyotrophic Lateral Sclerosis (Rheenen et al., 2016 Nat Genet [27455348](https://pubmed.ncbi.nlm.nih.gov/27455348))
* Anorexia Nervosa (Watson et al., 2019 Nat Genet [31308545](https://pubmed.ncbi.nlm.nih.gov/31308545))
* Attention Deficit Hyperactivity Disorder (Demontis et al., 2019 Nat Genet [30478444]())
* Autism Spectrum Disorders (Grove et al., 2019 Nat Genet [30804558](https://pubmed.ncbi.nlm.nih.gov/30804558))
* Bipolar Disorder (Stahl et al., 2019 Nat Genet [31043756](https://pubmed.ncbi.nlm.nih.gov/31043756))
* Cannabis Use Disorder (CUD) (Demontis et al., 2019 Nat Neurosci [31209380](https://pubmed.ncbi.nlm.nih.gov/31209380))
* Childhood IQ (Benyamin et al., 2014 Mol Psychiatry [23358156](https://pubmed.ncbi.nlm.nih.gov/23358156))
* Cigarettes Per Day (Liu et al., 2019 Nat Genet [30643251](https://pubmed.ncbi.nlm.nih.gov/30643251/))
* College Completion (Rietveld et al., 2013 Science [23722424](https://pubmed.ncbi.nlm.nih.gov/23722424))
* Cotinine Levels (Ware et al., 2016 Sci Rep [26833182](https://pubmed.ncbi.nlm.nih.gov/26833182/))
* Fagerstrom Test for Nicotine Dependence (FTND) (Quach et al., 2020 Nat Commun [33144568](https://pubmed.ncbi.nlm.nih.gov/33144568/))
* Heaviness of Smoking Index (HSI) (Quach et al., 2020 Nat Commun [33144568](https://pubmed.ncbi.nlm.nih.gov/33144568/))
* Insomnia (Jansen et al., 2019 Nat Genet [30804565](https://pubmed.ncbi.nlm.nih.gov/30804565/))
* Insomnia (Lane et al., 2019 Nat Genet [30804566](https://pubmed.ncbi.nlm.nih.gov/30804566/))
* Intelligence (Sniekers et al., 2017 Nat Genet [28530673](https://pubmed.ncbi.nlm.nih.gov/28530673))
* Lifetime Cannabis Use (Ever vs. Never) (Pasman et al., 2018 Nat Neurosci [30150663](https://pubmed.ncbi.nlm.nih.gov/30150663))
* LongSleepDur (Dashti et al., 2019 Nat Commun [30846698](https://pubmed.ncbi.nlm.nih.gov/30846698/))
* Major Depressive Disorder (Howard et al., 2018 Nat Commun [29662059](https://pubmed.ncbi.nlm.nih.gov/29662059))
* Mean Accumbens Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Mean Caudate Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Mean Hippocampus Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Mean Pallidum Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Mean Putamen Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Mean Thalamus Volume (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Neo-conscientiousness (de Moor et al., 2012 Mol Psychiatry [21173776](https://pubmed.ncbi.nlm.nih.gov/21173776))
* Neo-openness to Experience (de Moor et al., 2012 Mol Psychiatry [21173776](https://pubmed.ncbi.nlm.nih.gov/21173776))
* Neuroticism (Okbay et al., 2016 Nat Genet [27089181]())
* Opioid Addiction: GENOA GWAS meta-analysis
* Opioid Addiction: gSEM OA GWAS meta-analysis (i.e., GENOA, MVP-SAGE-YP, PGC-SUD, and Partners Health)
* Parkinson's Disease (Sanchez et al., 2009 Nat Genet [19915575](https://pubmed.ncbi.nlm.nih.gov/19915575))
* Post-traumatic Stress Disorder (Nievergelt et al., 2019 Nat Commun [31594949](https://pubmed.ncbi.nlm.nih.gov/31594949))
* Psychiatric Genetics Consortium Cross-disorder GWAS (Schizophrenia, Bipolar Disorder, MDD, ASD and ADHD) (Cross-Disorder Group of the Psychiatric Genomics Consortium, 2013 Lancet [23453885](https://pubmed.ncbi.nlm.nih.gov/23453885))
* Schizophrenia (Ripke et al., 2014 Nature [25056061](https://pubmed.ncbi.nlm.nih.gov/25056061))
* ShortSleepDur (Dashti et al., 2019 Nat Commun [30846698](https://pubmed.ncbi.nlm.nih.gov/30846698/))
* sleepDuration (Dashti et al., 2019 Nat Commun [30846698](https://pubmed.ncbi.nlm.nih.gov/30846698/))
* Sleepdur (Jansen et al., 2019 Nat Genet [30804565](https://pubmed.ncbi.nlm.nih.gov/30804565/))
* Smoking Cessation (Liu et al., 2019 Nat Genet [30643251](https://pubmed.ncbi.nlm.nih.gov/30643251/))
* Smoking Initiation (Liu et al., 2019 Nat Genet [30643251](https://pubmed.ncbi.nlm.nih.gov/30643251/))
* Subjective Well Being (Okbay et al., 2016 Nat Genet [27089181](https://pubmed.ncbi.nlm.nih.gov/27089181))
* Total Intracranial Volume (ICV) (Hibar et al., 2015 Nature [25607358](https://pubmed.ncbi.nlm.nih.gov/25607358/))
* Years of Schooling (Okbay et al., 2016 Nature [27225129](https://pubmed.ncbi.nlm.nih.gov/27225129))
</details><br><br>
    
    


## `20230425_oud_stratified_ldsc_on_twas_meta.ipynb`
- 47 traits (including 6 additional sleep studies)
- meta-analysis significant genes at FDR < 0.05 and < 0.10 from the results file `meta_analysis_sumstats_no_singletons_20220727.tsv.gz.`
    - See [this GitHub comment](https://github.com/RTIInternational/bioinformatics/issues/183#issuecomment-1198311459) we used the results with singletons removed: genes that only appear in one study.



<br><br>

## `20230426_oud_stratified_ldsc_on_twas_meta.ipynb`
Rerun the s-LDSC using a new gene-set list. In particular, a new TWAS meta-analysis was performed thus we will obtain a new set of significantly differentially expressed genes.
- List of genes with a Benjamini-Hochberg FDR < 0.05
