##This is an automated pipeline for preparing the genotype data for a GWAS using the ProbABEL tool. The steps in this pipeline are:
##
##1. Merge chrX
##    * The chrX imputed data from the Michigan Imputation Server is split between males and females. We need to merge these data for the GWAS because the other data is not split by sex. We control for that with a sex covariate.
##2. Convert imputed data format
##    * The ProbABEL software requires the data to be formatted differently than the VCF format that is the output from MIS. In particular, we need to convert the VCF files (dosage) to MaCH format.
##3. Prune imputed data
##    * Prune the imputed data to only the subjects with phenotype data
##4. Reorder phenotype file 
##    * order the subjects in phenotype file to be the same order as the genotype data (takes less memory to reorder the phenotype file rather than all of the genotype files)
##5. Create legend file
##6. Reformat info file
##
##<br>
##<br>
##
##**Note1:** that the genomic data need to be inflated after they have been download from the Michigan Imputation Server.
##
##**Note2:** this pipeline will perform the actual GWAS (1df and 2df) as well. For the 2df analysis, you might have to adjust the `interaction` parameter in the ProbABEL command. The number you provide should be the number of columns to the right of the phenotype variable that the interaction variable is in your phenotype file. For example, the interaction variable is currently coded as 2 in this pipeline for the 2df GWAS. It is coded as a 2 because we are wanting to perform a GxSex analysis for the 2df GWAS and the sex variable `gender` is two columns to the right of the phenotype variable `hiv` in the phenotype file.
##
##```
##iid hiv age gender PC10 PC9 PC2 PC6
##109@1064714572_109@1064714572 0 26 1 0.0011 0.0075 0.0039 -3e-04
##202@1064714531_202@1064714531 0 27 0 0.0093 -0.0112 0.0023 0.0017
##312@1064714548_312@1064714548 0 34 1 -0.0012 -0.0012 4e-04 0.0173
##378@1064610814_378@1064610814 0 30 1 0.0113 0.0016 0.0026 0.0043
##```
##
##**Note3:** pay close attention to the variables above the double octothorpe line. Be sure to edit them accordingly, like for example the `ancestry, covars` and the `phenotype_file` variables. The `s3_upload` variable should be changed per study as well.
##
##**Note4:** we recommend running the first set of functions first, and then when they finish go on to set number two. This is because the first set of functions do not require as much memory as the second set of functions. For the first set of functions, use a smaller compute node (e.g. the m5.large) and then use a memory intensive compute node (e.g. x1e.2xlarge).
##
##**Note5:** compute the chromosomes in batches. For ProbABEL to run, one needs to decompress the chromosomes. This can create storage issues; chromosomes 2 decompressed for the UHS1234 AA group is 92GB! It is therefore advised to process just a few chromosomes at a time so as not to run into storage issues. 

####################################################################################################
####################################################################################################
## How to run pipeline
**Note on performing the GWAS:**
```
The ProbABEL software requires the genotype data to be decompressed so that it requires a lot of storage. For this reason, it is best to run only portions of the genome at a time. An advisable way to run the pipeline above is as follows:

1. Convert all of the genetic data to the correct format by running the functions:
    -merge_chrx  
    -convert_chrx 
    -convert_auto 
    -create_legend 
    -format_info
   
        Comment out the functions:
            -prune_geno
            -gwas_1df          
            -gwas_2df
            -status_check 
            
2. Prune the genotype data to match the phenotype file with:
    -prune_geno
    
        Comment out the functions:
            -merge_chrx  
            -convert_chrx 
            -convert_auto 
            -create_legend 
            -format_info
            -gwas_1df          
            -gwas_2df
            -status_check 
            
3. Perform the 1df & 2df GWAS and then gzip the genotype data afterwards:
    -gwas_1df          
    -gwas_2df
    -status_check 
            
        Comment out the functions:
            -merge_chrx  
            -convert_chrx 
            -convert_auto 
            -create_legend 
            -format_info
            -prune_geno
```

Perform step 1 for all of the chromosomes. Then perform step 2 & 3 with just a couple or three chromosomes at a time. For step 1 and 2, use a smaller (and less expensive) compute node like the m5.large. When performing (3)—the ProbABEL 1df & 2df GWAS—one will have to use a compute node with very large memory (e.g. x1e.2xlarge). See below for example of how much memory just one chromosome is using for a 1df GWAS.

```
HOSTNAME                ARCH         NCPU NSOC NCOR NTHR  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
----------------------------------------------------------------------------------------------
ip-172-31-22-183        lx-amd64        8    1    4    8  1.24  240.1G  211.4G     0.0     0.0
```

####################################################################################################
####################################################################################################
## bash script ##

chr=$1  # command line argument (chromosome number)

ancestry="ea"
covars="sex,age,PC1,PC9,PC10"  # comma separated
study="uhs1234"
version="001"  # gwas attempt/version (e.g. 001)
probabel="palogist" # palogist (logistic) or palinear (linear)

baseD=/shared/jmarks/hiv/uhs1234/acquisition_gwas
imputeD=$baseD/genotype/imputed/$ancestry
phenD=$baseD/phenotype/final
procD=$baseD/gwas/1df/$ancestry/$version # output of 1df test and log files
procD2=$baseD/gwas/2df/$ancestry/$version # output for 2df test

phenotype_file=$phenD/uhs1234_ea_hiv_age_gender_PC1+PC9+PC10.txt
pheno_order=$phenD/phenotype_ids_${ancestry} # list of ids extracted from the phenotype file 

# path to the order of the subject IDs in the pruned genotype files 
# (output from reording genotype files with python script)
id_order=$baseD/genotype/imputed/$ancestry/mach/chr$chr/chr$chr.genotype.id.order  

# where to upload the genotype data after GWAS are complete 
s3_upload=s3://rti-hiv/hiv_uhs1234/data/genotype/imputed/$ancestry/output/gwas_processed/chr$chr/

if [ $ancestry == "aa" ]; then
    group=afr
elif [ $ancestry == "ea" ]; then
    group=eur
fi

# prepare directory structure
for dir in 1df 2df; do
    mkdir -p $procD/final
    mkdir -p $procD/processing/chr{1..23}
    mkdir -p $procD2/final
    mkdir -p $procD2/processing/chr{1..23}

done


# Don't edit below this line
####################################################################################################
####################################################################################################
main() {
    echo -e "\n\nHello, initializing the ProbABEL GWAS pipeline."
    
    ## SET1: These functions can be run with smaller compute nodes (m5.large)
    merge_chrx 
    convert_chrx 
    convert_auto 
    create_legend 
    format_info
    prune_geno $chr $imputeD/mach/chr$chr $pheno_order

    ## SET2: Run the functions below with larger memory compute nodes (x1e.2xlarge)
    #gwas_1df          
    #gwas_2df
    #status_check "yes"  # "yes" (zip data and upload to S3) or "no" (don't zip/upload data)
}
####################################################################################################
####################################################################################################

# merge chrx male & female data
merge_chrx() {
    echo -e "\n\n\n################################################################################"
    echo -e "START: merge_chrx() function for chr$chr"
    echo -e "################################################################################"

    if [ $chr == 23 ]; then
        echo "Merging chrx for males & females..."

        /shared/bioinformatics/software/third_party/bcftools-1.6/bcftools merge \
            $imputeD/chrX.no.auto_male.dose.vcf.gz $imputeD/chrX.no.auto_female.dose.vcf.gz \
            -O z -o $imputeD/chrX.no.auto.dose.vcf.gz
    else
        echo "chr$chr is an autosome. Exiting merge_chrx() function."
    fi
    
    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: merge_chrx() function for chr$chr\n\n\n\n"
}

# convert chrx data to ProbABEL acceptable format
convert_chrx() {

    echo -e "\n\n\n################################################################################"
    echo "START: convert_chrx() function for chr$chr"
    echo -e "################################################################################\n"
    echo "Converting chrx data to ProbABEL acceptable format...\n"

    if [ $chr == 23 ]; then
        chr_location=$imputeD/mach/chr${chr}
        mkdir -p $chr_location
        /shared/bioinformatics/software/third_party/dosage_converter_v1.0.4/bin/DosageConvertor \
            --vcfDose $imputeD/chrX.no.auto.dose.vcf.gz \
            --prefix ${chr_location}/chr${chr} \
            --type mach \
            --format 1 # contains the expected alternate allele count (one value per sample per marker).
    else
        echo "chr$chr is an autosome. Exiting convert_chrx() function."
    fi

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: convert_chrx() function for chr$chr\n\n\n\n"
}



# convert autosomal data to ProbABEL acceptable format
convert_auto() {
    echo -e "\n\n\n################################################################################"
    echo "START: convert_auto() function for chr$chr"
    echo -e "################################################################################\n"

    chr_location=$imputeD/mach/chr$chr
    mkdir -p $chr_location
    /shared/bioinformatics/software/third_party/dosage_converter_v1.0.4/bin/DosageConvertor \
        --vcfDose $imputeD/chr$chr.dose.vcf.gz \
        --info $imputeD/chr$chr.info.gz \
        --prefix $chr_location/chr$chr \
        --type mach \
        --format 1 # contains the expected alternate allele count (one value per sample per marker).
    
    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: convert_auto() function for chr$chr\n\n\n\n"
}


# create legend file for ProbABEL
create_legend() {
    echo -e "\n\n\n################################################################################"
    echo "START: create_legend() function for chr$chr"
    echo -e "################################################################################\n"

    # HapMap "legend" file format
    # note that our data is not in rsID format - we have chr:position for the rsID instead
    echo "Creating the legend file for ProbABEL input."
    echo "id position 0 1" > $imputeD/mach/chr$chr/map.chr$chr.legend
    # grab the SNP, position, allele1 and allele2
    tail -n +2 $imputeD/mach/chr$chr/*info |\
    awk '{pos = $1; gsub(/^.+:/, "", pos); print $1,pos,$2,$3}' >>\
    $imputeD/mach/chr$chr/map.chr$chr.legend

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: create_legend() function for chr$chr\n\n\n\n"
}



# format the info file for ProbABEL
format_info() {
    echo -e "\n\n\n################################################################################"
    echo "START: format_info() function for chr$chr"
    echo -e "################################################################################\n"
    echo "Pruning down the info file to be compliant with ProbABEL."
    cut -f 1-7 $imputeD/mach/chr${chr}/chr${chr}.mach.info >\
        $imputeD/mach/chr$chr/chr$chr.mach.info.pruned  \

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: format_info() function  for chr$chr\n\n\n\n"
}



# prune the genotype data
prune_geno() {
    echo -e "\n\n\n################################################################################"
    echo "START: prune_geno() function for chr$chr"
    echo -e "################################################################################\n"

    baseD=$2
    idL=$3
    time python ~/bin/prune_genotype_files.py $chr $baseD $idL

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: prune_geno() function for chr$chr\n\n\n\n"
}


## perform 1df GWAS (logistic)
gwas_1df() {
    echo -e "\n\n\n################################################################################"
    echo "START: gwas_1df() function for chr$chr"
    echo -e "################################################################################\n"
    echo "First, reorder the phenotype file so that the subject IDs are in the same order as the genotype data..."

    if [ $chr == 23 ]; then
        phenoF=${phenotype_file::-4}_ordered_chrx.txt
    else
        phenoF=${phenotype_file::-4}_ordered.txt
    fi

    # check to see if file has already been created
    if [ -e $phenoF ]; then
        echo -e "\t...The phenotype file: $phenoF already exists. Now on to the 1df GWAS.\n"
    else
        head -1 $phenotype_file > $phenoF
        awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' \
            $phenotype_file $id_order >> $phenoF
        echo -e "Phenotype file: $phenoF has now been created. On to the 1df GWAS.\n"
    fi

    /shared/bioinformatics/software/third_party/probabel-0.5.0/bin/$probabel \
        --pheno $phenoF \
        --dose  $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned \
        --info  $imputeD/mach/chr$chr/chr$chr.mach.info.pruned \
        --map   $imputeD/mach/chr$chr/map.chr$chr.legend \
        --chrom $chr \
        --out   $procD/processing/chr$chr/chr$chr.$probabel.results

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: gwas_1df() function for chr$chr\n\n\n\n"
}


## perform 2df GWAS (logistic)
gwas_2df() {
    echo -e "\n\n\n################################################################################"
    echo "START gwas_2df() function for chr$chr"
    echo -e "################################################################################\n"
    echo "First, reorder the phenotype file so that the subject IDs are in the same order as the genotype data."

    if [ $chr == 23 ]; then
        phenoF=${phenotype_file::-4}_ordered_chrx.txt
    else
        phenoF=${phenotype_file::-4}_ordered.txt
    fi

    # check to see if file has already been created
    if [ -e $phenoF ]; then
        echo -e "The phenotype file: $phenoF already exists. Now on to the 2df GWAS. \n"
    else
        head -1 $phenotype_file > $phenoF
        awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' \
            $phenotype_file $id_order >> $phenoF
        echo -e "Phenotype file: $phenoF has now been created. On to the 2df GWAS.\n"
    fi

    /shared/bioinformatics/software/third_party/probabel-0.5.0/bin/$probabel \
        --pheno $phenoF \
        --dose  $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned \
        --info  $imputeD/mach/chr$chr/chr$chr.mach.info.pruned \
        --map   $imputeD/mach/chr$chr/map.chr$chr.legend \
        --chrom $chr \
        --interaction 2 \
        --out   $procD2/processing/chr$chr/chr$chr.$probabel.results 

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: gwas_2df() function for chr$chr\n\n\n\n"
}


## Determine if GWAS was completed and gzip imputed genetic data ##
status_check() {
    echo -e "\n\n\n################################################################################"
    echo "START: status_check() function for chr$chr"
    echo -e "################################################################################\n"


    myzip=$1 # should be yes or no
    
    # if there are two 100.00% then both GWAS finished
    stat="$(grep  100.00% $procD/processing/chr$chr/*log | wc -l)"
    
    if [ $stat == "2" ]; then
        echo "GWAS COMPLETE!!"
        echo "Both 1df and 2df GWAS are complete for chr$chr."

        if [ $myzip == "yes" ]; then
            echo -e "Attempting to gzip genotype file..."

            if [ -e $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned ]; then
                echo "Zipping genotype file..."
                gzip $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned
                echo -e "Gzip was successful.\n\n"
                echo -e "Uploading $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned.gz $to $s3_upload...\n"
                aws s3 sync $imputeD/mach/chr$chr $s3_upload 
                echo -e "\nUpload complete.\n"
            else
                echo "The file $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned does not exist."
                if [ -e $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned.gz ]; then
                    echo "This file has already been gzipped."
                    echo -e "Uploading $imputeD/mach/chr$chr/chr$chr.mach.dose.pruned.gz $to $s3_upload...\n"
                    aws s3 sync $imputeD/mach/chr$chr $s3_upload 
                    echo -e "\nUpload complete."
                fi
            fi
        else
           echo "You have chosen not to zip & upload. Exiting now. Thank you, have a nice day! "    
        fi
    else
        echo "WARNING: The 1df and/or 2df GWAS may not have finished for chr$chr."
        echo "Please check log file: $procD/processing/chr$chr/*log"
    fi

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: status_check() function for chr$chr\n"
}

####################################################################################################
####################################################################################################
time main $1  # $1 is chromosome, 


