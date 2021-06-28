### bash script ###
chr=$1

ancestry="aa"
study="uhs1234"
phenotype="hiv_acq"
version=001
degree_free="2df"
cmax=23 # number of chromosomes (22 or 23)
probabel=palogist


baseD=/shared/jmarks/hiv/uhs1234/acquisition_gwas
procD=$baseD/gwas/$degree_free/$ancestry/$version/processing
phenD=$baseD/phenotype/final   # phenotype directory

# path to chrX info file for males and females (including file name)
maleInfo=$baseD/genotype/imputed/$ancestry/chrX.no.auto_male.info.gz 
femaleInfo=$baseD/genotype/imputed/$ancestry/chrX.no.auto_female.info.gz 

if [ $ancestry == "aa" ]; then
    group=afr
    phenoXF=uhs1234_aa_hiv_age_gender_PC10+PC9+PC2+PC6_ordered_chrx.txt  # name of chromosome X phenotype file

    if [ $degree_free == "1df" ];then
        MODEL=HIV_ACQ~SNP+AGE+SEX+PC10+PC9+PC2+PC6  # for 1df
    else
        MODEL=HIV_ACQ~SNP+AGE+SEX+SNPbySEX+PC10+PC9+PC2+PC6  # for 2df
    fi
        
elif [ $ancestry == "ea" ]; then
    group=eur
    phenoXF=uhs1234_ea_hiv_age_gender_PC1+PC9+PC10_ordered_chrx.txt  # name of chromosome X phenotype file

    if [ $degree_free == "1df" ];then
        MODEL=HIV_ACQ~SNP+AGE+SEX+PC1+PC9+PC10  # for 1df
    else
        MODEL=HIV_ACQ~SNP+AGE+SEX+SNPbySEX+PC1+PC9+PC10  # for 2df
    fi
fi
####################################################################################################
####################################################################################################
main() {
    echo -e "Results Processing\n\n"

    calc_pvalue
    chrx_maf
    autosomes_1kg 
    chrx_1kg
    maf_study
    maf_1kg
    rsq
    #plot_results
    #pvalue_filter
}
####################################################################################################
####################################################################################################

## Calculate chi^2, Pvalue, and odds ratio ## 
calc_pvalue() {
    echo -e "\n################################################################################"
    echo  "START: calc_pvalue() function"
    echo -e "################################################################################\n\n"

    echo "$degree_free p-value calculation for $study $ancestry chr$chr..."
    if [ $degree_free == "1df" ]; then
        Rscript ~/bin/calculate_stats_for_probabel_results.R \
            --in_file $procD/chr$chr/chr$chr.$probabel.results_add.out.txt \
            --out_file $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats \
            --complete
    else

        Rscript ~/bin/calculate_2df_stats_from_probabel_results.R \
            --in $procD/chr$chr/chr$chr.$probabel.results_add.out.txt \
            --out $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats \
            --interaction_covar gender

        echo -e "Wrote to: $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats"
    fi
        
    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: calc_pvalue() function\n\n"
}

chrx_maf() {
    
    echo -e "\n\n################################################################################"
    echo -e "START: chrx_maf()"
    echo -e "################################################################################\n"

    if [ $chr == 23 ]; then
    
        ## calculate chrX MAF ##
        chrxD=$procD/chr23

        # count number of males and females (need this for new MAF calculation)
        femaleN=$(awk '$4==1' $phenD/$phenoXF | wc -l)
        maleN=$(awk '$4==0' $phenD/$phenoXF | wc -l)

        echo "Number of females: $femaleN"
        echo -e "Number of males: $maleN\n"

        # calculate new MAF by recalculating the combined sex specific stats (first get overlapping SNPs)
        cut -f1 <(zcat $femaleInfo | tail -n +2) | sort >\
           $chrxD/female_snps.txt 
        cut -f1 <(zcat $maleInfo | tail -n +2) | sort >\
           $chrxD/male_snps.txt 
        comm -12 $chrxD/male_snps.txt $chrxD/female_snps.txt >\
            $chrxD/male_female_snp_intersection.txt

        echo -e "Number of female SNPs: $(wc -l $chrxD/female_snps.txt)"
        echo -e "Number of male SNPs: $(wc -l $chrxD/male_snps.txt)\n"

        # filter info file down to the intersection set of SNPs (should be a large intersection)
        zcat $femaleInfo | head -1 > $chrxD/female_intersection_info.txt
        zcat $maleInfo| head -1  > $chrxD/male_intersection_info.txt

        # new female info file
        awk 'NR==FNR{map[$1]=$1; next} { if ($1 in map) {print $0} }' \
            $chrxD/male_female_snp_intersection.txt <(zcat $femaleInfo) >>\
            $chrxD/female_intersection_info.txt
        # new male info file
        awk 'NR==FNR{map[$1]=$1; next} { if ($1 in map) {print $0} }' \
            $chrxD/male_female_snp_intersection.txt <(zcat $maleInfo) >>\
            $chrxD/male_intersection_info.txt

        # Calculate new merged MAF
        awk -v fN=$femaleN -v mN=$maleN 'NR==FNR{female[$1] = $5; next} 
            FNR>=2{$5 = ( female[$1]*fN + $5*mN ) / (fN + mN)}
            {print $0}' $chrxD/female_intersection_info.txt\
            $chrxD/male_intersection_info.txt >\
            $chrxD/merged_male_female.info

        # insert recalculated MAF into stats file and filter to intersection
        head -1 $chrxD/$study.$ancestry.$degree_free.1000G_p3.chr23.$MODEL.stats >\
              $chrxD/$study.$ancestry.$degree_free.1000G_p3.chr23.$MODEL.recalc_maf.stats
              
        awk 'FNR==NR{map[$1]=$5; next} 
            { if ($1 in map) 
            {$7= map[$1]; print $0 } }'  $chrxD/merged_male_female.info \
            $chrxD/$study.$ancestry.$degree_free.1000G_p3.chr23.$MODEL.stats >>\
            $chrxD/$study.$ancestry.$degree_free.1000G_p3.chr23.$MODEL.recalc_maf.stats

        rm $chrxD/male_snps.txt
        rm $chrxD/female_snps.txt
        rm $chrxD/male_female_snp_intersection.txt
        rm $chrxD/female_intersection_info.txt
        rm $chrxD/male_intersection_info.txt
        rm $chrxD/merged_male_female.info
    else
        echo -e "WARNING: This is an autosome (not chrX). Exiting now.\n"
    fi


    echo -e "--------------------------------------------------------------------------------"
    echo -e "END: chrx_maf() function\n\n"
}

              
## convert SNP names to 1KG_p3
autosomes_1kg() {
    echo -e "\n\n################################################################################"
    echo -e "START: autosomes_1kg() function"
    echo -e "################################################################################\n"
    echo -e "Converting the SNP names to 1000 Genomes Phase 3 format (rsID:position:A1:A2).\n\n"

    if [ $chr -lt 23 ]; then
        if [ -e $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats ]; then
            perl /shared/bioinformatics/software/perl/id_conversion/convert_to_1000g_p3_ids.pl \
                --file_in $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats \
                --file_out $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats.converted \
                --legend /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.gz \
                --file_in_header 1 \
                --file_in_id_col  0 \
                --file_in_chr_col  1 \
                --file_in_pos_col  2 \
                --file_in_a1_col  3 \
                --file_in_a2_col  4 \
                --chr $chr
        else
            echo "File $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats does not exist."
        fi
    else
        echo "WARNING: chr$chr is not an autosome. Exiting"
    fi

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: autosomes_1kg() function.\n\n"
}

chrx_1kg() {
    
    echo -e "\n\n################################################################################"
    echo -e "START: chrx_1kg() function"
    echo -e "################################################################################\n"
    echo -e "Converting the SNP names to 1000 Genomes Phase 3 format (rsID:position:A1:A2).\n\n"

    if [ $chr == 23 ]; then
        if [ -e $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.recalc_maf.stats ]; then
            perl /shared/bioinformatics/software/perl/id_conversion/convert_to_1000g_p3_ids.pl \
                --file_in $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.recalc_maf.stats \
                --file_out $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats.converted \
                --legend /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chrX_NONPAR.legend.gz \
                --file_in_header 1 \
                --file_in_id_col  0 \
                --file_in_chr_col  1 \
                --file_in_pos_col  2 \
                --file_in_a1_col  3 \
                --file_in_a2_col  4 \
                --chr $chr 
        else
            echo "The file $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.recalc_maf.stats does not exist. Exiting."
        fi
    else
        echo -e "WARNING: chr$chr is not chrX. Exiting"
    fi
    
    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: chrx_1kg() function.\n\n"
}

# Remove variants with MAF <= 0.01 in the study data
maf_study() {
        
        echo -e "\n\n################################################################################"
        echo -e "START: maf_study() function."
        echo -e "################################################################################\n"
        echo -e "Removing variants with MAF <= 0.01 in the study data."

        echo "Processing chr${chr}_${ancestry} for MAF (study) filtering."
        head -n 1 $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats.converted > \
            $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject.stats

        # note column 7 corresponds to the MAF column 
        awk ' NR>=2 {if ($7 >= 0.01) {print $0}}' \
            $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.stats.converted \
            >> $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject.stats

        echo -e "\n--------------------------------------------------------------------------------"
        echo -e "END: maf_study() function.\n\n"
}
              
              
# Remove variants with MAF <= 0.01 in the 1KG data
maf_1kg() {

    echo -e "\n\n################################################################################"
    echo -e "START: maf_1kg() function."
    echo -e "################################################################################\n"
    echo -e "Removing variants with MAF <= 0.01 in the 1000 Genomes data."

    # creating a list of SNPs based off of 1000G population 
    # filter the variants to ones whose MAF <= 1%
    if [ $chr -lt 23 ]; then 
       awk ' { if ($9 >= 0.01) {print $1}}' \
       <(zcat /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.gz) >\
        /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.unique_ids.maf_gt_0.01_${group}

    # chrX
    elif [ $chr -eq 23 ]; then
        awk ' { if ($9 >= 0.01) {print $1}}' \
            <(zcat /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chrX_NONPAR.legend.gz) >\
            /shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.unique_ids.maf_gt_0.01_${group}
    else
        echo "ERROR: chr$chr is an invalid entry. Exiting"

    fi

    idList=/shared/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.unique_ids.maf_gt_0.01_${group} 
    echo "Processing chr${chr}_${ancestry}"      
    /shared/bioinformatics/software/perl/utilities/extract_rows.pl \
        --source $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject.stats \
        --id_list $idList \
        --out $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.stats \
        --header 1 \
        --id_column 0 

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: maf_1kg() function.\n\n"
}
              
            
## Remove SNPs with poor imputation quality (r^2 <= 0.30)
rsq() {
    echo -e "\n\n################################################################################"
    echo -e "START: rsq() function"
    echo -e "################################################################################\n"
    echo -e "Removing variants with poor imputation quality (r^2 <= 0.30).\n"
    
    # autosomes
    if [ $chr -lt 23 ]; then
        echo "Processing $ancestry chr$chr..."
        awk 'FNR==1 {print $0} FNR>=2{ if ($9 > 0.3){ print $0 } }' \
            $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.stats >\
             $procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats
    
    # chrX
    elif [ $chr -eq 23 ]; then

    # chrX imputed data are split up by males & females 
    # perform the filtering on both male and female data and then merge results
        echo -e "Processing ${ancestry} chrX..."

        # male variants passing rsq filter
        tail -n +2 <(zcat $maleInfo) | \
            awk '{ if( $7 > 0.3){ print $1":"$2":"$3 } }' \
            > $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep.tmp

        # female variants passing rsq filter
        tail -n +2 <(zcat $femaleInfo) | \
            awk '{ if( $7 > 0.3){ print $1":"$2":"$3 } }' \
            >> $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep.tmp

        # keep only SNPs that passed filters for both males and females
        sort $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep.tmp |\
            uniq -d > $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep 

        # change X to 23 to be compliant with our plotting script  
        awk -F":" '$1=23 {print $0}' $procD/chr23/$study.chr${chr}_variants_rsq_gt_0.3.keep >\
            $procD/chr23/$study.chr${chr}_variants_rsq_gt_0.3.keep.split

        # filter results
        awk ' NR==FNR { map[$1":"$2":"$3":"$4] = 1; next}
            FNR==1 {print $0}
            FNR>=2 {if  (map[$2":"$3":"$4":"$5] == 1)
            { print $0} }' $procD/chr23/$study.chr${chr}_variants_rsq_gt_0.3.keep.split \
            $procD/chr$chr/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.stats >\
            $procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats

        rm $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep.tmp
        rm $procD/chr$chr/$study.chr${chr}_variants_rsq_gt_0.3.keep
        rm $procD/chr23/$study.chr${chr}_variants_rsq_gt_0.3.keep.split

    # exit function if a valid chromosome was not supplied as input.
    else
        echo "ERROR:  chr$chr is an invalid entry. Please enter a chromosome number {1..23}. Exiting"
    fi

    echo -e "Done"
    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: rsq() function.\n\n"

}
              
              
## Make Q-Q and manhattan plots ##
plot_results() {

    echo -e "\n\n################################################################################"
    echo -e "START: plot_results() function"
    echo -e "#################################################################################\n"
    echo -e "Creating the Manhattan and QQ plots.\n\n"

    if [ $degree_free == "1df" ];then
        skip_cols=10   # 1df
    else
        skip_cols=13   # 2df
    fi
    
    # create table of SNPs and corresponding pvalues, etc
    outfile=$procD/${study}.${ancestry}.1000G_p3.${phenotype}.maf_gt_0.01.rsq_gt_0.3.assoc.table
    echo -e "VARIANT_ID\tCHR\tPOSITION\tP\tTYPE" > $outfile

    for ((chr=1; chr<=$cmax; chr++));do
        infile=$procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats
        echo Processing $infile
        tail -n +2 $infile |
          perl -slne '/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+\S+){$myvar}\s+(\S+)/;
                      if (($4 eq "A" || $4 eq "C" || $4 eq "G" || $4 eq "T") && ($5 eq "A" || $5 eq "C" || $5 eq "G" || $5 eq "T")) {
                        print join("\t",$1,$2,$3,$6,"snp");
                      } else {
                        print join("\t",$1,$2,$3,$6,"indel");
                      }' -- -myvar=$skip_cols >> $outfile
    done 

    Rscript /shared/bioinformatics/software/R/generate_gwas_plots.R \
        --in $outfile \
        --in_chromosomes autosomal_nonPAR \
        --in_header \
        --out $procD/../final/$study.${ancestry}.$degree_free.1000G.${phenotype}.maf_gt_0.01.rsq_gt_0.3.assoc.plot.all_chr \
        --col_id VARIANT_ID \
        --col_chromosome CHR \
        --col_position POSITION \
        --col_p P \
        --col_variant_type TYPE \
        --generate_snp_indel_manhattan_plot \
        --manhattan_odd_chr_color red3 \
        --manhattan_even_chr_color dodgerblue3 \
        --manhattan_points_cex 1.5 \
        --generate_snp_indel_qq_plot \
        --qq_lines \
        --qq_points_bg black \
        --qq_lambda  > $procD/$study.${ancestry}.$degree_free.1000G.${phenotype}.maf_gt_0.01.rsq_gt_0.3.assoc.plot 2>&1

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: plot_results() function.\n\n"
}
              
              
    
# filter to top SNPs (p <= 0.0001)
pvalue_filter() {

    echo -e "\n\n################################################################################"
    echo -e "START: pvalue_filter() function."
    echo -e "################################################################################\n\n"
    echo -e "Filtering the results to only the most signifanct SNPs (pvalue <= 0.0001).\n\n"

    if [ $degree_free == "1df" ];then
        pcol=15   # 1df
    else
        pcol=18   # 2df
    fi

    outFile=$procD/../final/$study.$ancestry.$degree_free.1000G_p3.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.p_lte_0.001
    

    head -n 1 $procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr22.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats > $outFile
    for ((chr=1; chr<=$cmax; chr++));do
            echo "Processing $procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats"
            tail -n +2 $procD/../final/$study.$ancestry.$degree_free.1000G_p3.chr$chr.$MODEL.maf_gt_0.01_subject+${group}.rsq_gt_0.30.stats |\
            perl -slane 'if ($F[$pcol] <= 0.001) { print; }' -- -pcol=$pcol >>  $outFile
    done

    echo -e "\n--------------------------------------------------------------------------------"
    echo -e "END: pvalue_filter() function.\n"
}

##################################################################################
time main
