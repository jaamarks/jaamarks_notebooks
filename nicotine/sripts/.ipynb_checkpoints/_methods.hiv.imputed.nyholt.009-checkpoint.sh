# Copy list from \\rcdcollaboration01.rti.ns\johnson_share\1.HIV GWAS II\technical\Meta\Meta001_GoGeneCandidates_ListForNyholtCorrection_05.26.16.xlsx

inFile=/share/nas04/bioinformatics_group/data/studies/hiv/imputed/nyholt/009/variant_list.original
outFile=/share/nas04/bioinformatics_group/data/studies/hiv/imputed/nyholt/009/variant_list.final
# Copy variants with phase 3 IDs to final file
perl -lane 'if ($F[0] =~ /\:/) { print $F[0]; }' $inFile >> $outFile
# Convert non-phase 3 variant IDs to phase 3 IDs based on name
for variant in $(perl -lane 'if ($F[0] !~ /\:/) { print $F[0]; }' $inFile); do
    chr=$(grep -P "$variant\s" $inFile | perl -lane 'print $F[1];')
    gunzip -c /data/common/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.gz |
        grep "$variant:" |
        perl -lane 'print $F[0];'
done > $outFile
# Convert non-phase 3 variant IDs to phase 3 IDs based on position
for position in $(perl -lane 'if ($F[0] !~ /\:/) { print $F[2]; }' $inFile); do
    chr=$(grep -P "$position" $inFile | grep -v ":" | perl -lane 'print $F[1];')
    gunzip -c /data/common/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.gz |
        grep ":$position:"
done > /share/nas04/bioinformatics_group/data/studies/hiv/imputed/nyholt/009/variant_list.nonphase3.legend
perl -lane 'BEGIN { $lastPosition = 0 } if ($F[1] != $lastPosition) { print $F[0]; $lastPosition = $F[1]; }' /share/nas04/bioinformatics_group/data/studies/hiv/imputed/nyholt/009/variant_list.nonphase3.legend >> $outFile

    
### START Extract replication SNPs from 1000G ###

# Extract variants from 1000G panel
for (( chr=1; chr<23; chr++ )); do
for chr in 1 2; do
    /data/common/software/scripts/qsub_job.sh \
      --job_name ALL_1000G \
      --script_prefix /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants \
      --mem 80 \
      --priority 0 \
      --program /data/common/software/perl/convert_reference_panels.v2.pl \
      --impute2_hap /data/common/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.hap.gz \
      --impute2_legend /data/common/data/ref_panels/1000G/2014.10/1000GP_Phase3_chr$chr.legend.gz \
      --impute2_sample /data/common/data/ref_panels/1000G/2014.10/1000GP_Phase3.sample \
      --extract /data/common/data/studies/hiv/imputed/nyholt/009/variant_list.final \
      --out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants \
      --generate_plink_ped_file \
      --generate_plink_map_file \
      --chr $chr
done

### END Extract replication SNPs from 1000G ###



### START Generate correlation matrix ###

for (( chr=1; chr<23; chr++ )); do
    /data/common/software/plink \
        --noweb \
        --file /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants.plink \
        --r \
        --matrix \
        --out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants.r
done

# for chr 2, rs568413:21235475:T:C appears to be causing trouble so remove it and recalculate matrix
echo rs568413:21235475:T:C > /data/common/data/studies/hiv/imputed/nyholt/009/rs568413
/data/common/software/plink \
    --noweb \
    --file /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.plink \
    --exclude /data/common/data/studies/hiv/imputed/nyholt/009/rs568413 \
    --recode \
    --out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.minus_rs568413.plink
/data/common/software/plink \
    --noweb \
    --file /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.minus_rs568413.plink \
    --r \
    --matrix \
    --out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.minus_rs568413.r

 ### END Generate correlation matrix ###


### START Run matSpD analysis ###

cd /data/common/software/matSpD/
for chr in 1 2 4 6 8 10 11 12 17 19; do
    echo $chr
    cp /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants.r.ld \
        /data/common/software/matSpD/correlation.matrix
    R CMD BATCH matSpDlite.R
    mv /data/common/software/matSpD/matSpDlite.out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants.matspdlite
done

cp /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.minus_rs568413.r.ld \
    /data/common/software/matSpD/correlation.matrix
R CMD BATCH matSpDlite.R
mv /data/common/software/matSpD/matSpDlite.out /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr2.extracted_variants.matspdlite

veffLi=0
for chr in 1 2 4 6 8 10 11 12 17 19; do
    chrVeffLi=$(grep -A 2 Equation /data/common/data/studies/hiv/imputed/nyholt/009/1000G_ALL.chr$chr.extracted_variants.matspdlite | tail -n +3 | perl -pe 's/^\s+//; s/\s+$//;')
    echo $chr
    echo $chrVeffLi
    veffLi=`echo $veffLi + $chrVeffLi | bc`
done
echo $veffLi

### END Run matSpD analysis ###


