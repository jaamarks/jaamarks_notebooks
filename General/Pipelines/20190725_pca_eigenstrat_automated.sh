#!/bin/bash

## Set global variables ##
anlist="ea" # "ea", "aa", or "ea aa"
study="vidus"
case_control="false"  # true if (case control phenotype) false if (not)
procD="/shared/jmarks/hiv/vidus/troubleshooting"  # processing directory
eig="$procD/eig" # do not alter
mkdir -p $eig/results # do not alter


## Genotype Data should be in the directory $eig
## The genotype data names should be of the form: $study.$an.genotype.{bed,bim,fam}. For Example:
## vidus.ea.genotype.bim, vidus.ea.genotype.bed, vidus.ea.genotype.fam


# Do not edit below this line
################################################################################
################################################################################

## Remove high-LD region variants ##
for an in $anlist; do
    perl -lane 'if (($F[0]==5 && $F[3] >= 43964243 && $F[3] <= 51464243) || ($F[0]==6 && $F[3] >= 24892021 && $F[3] <= 33392022) || ($F[0]==8 && $F[3] >= 7962590 && $F[3] <= 11962591) || ($F[0]==11 && $F[3] >= 45043424 && $F[3] <= 57243424)) { print $F[1]."\n"; }' \
        $eig/$study.$an.genotypes.bim > $eig/$study.$an.genotypes.high_ld_regions.remove

    # Remove SNPs in known high-LD regions
    /shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
        --noweb \
        --bfile $eig/$study.$an.genotypes \
        --exclude $eig/$study.$an.genotypes.high_ld_regions.remove \
        --make-bed \
        --out $eig/$study.$an.genotypes.high_ld_regions_removed
done

## Linkage disequilibrium pruning ##
for an in $anlist;do
    for chr in {1..22}; do
        /shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
            --noweb \
            --memory 3000 \
            --bfile $eig/$study.$an.genotypes.high_ld_regions_removed \
            --indep-pairwise 1500 150 0.2 \
            --chr ${chr} \
            --out $eig/${an}_chr${chr}_ld_pruned
    done
done

## Merge *prune.in files ##
for an in $anlist; do
    cat $eig/${an}_chr*_ld_pruned.prune.in > $eig/${an}_chr_all_ld_pruned.prune.in

    # Create new PLINK filesets with only lD pruned variants
    /shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
        --noweb \
        --bfile $eig/$study.$an.genotypes.high_ld_regions_removed \
        --extract $eig/${an}_chr_all_ld_pruned.prune.in \
        --make-bed \
        --out $eig/$study.$an.genotypes.ld_pruned
done

# Clean up
rm $eig/*ld_pruned.{prune.in,prune.out,log}
rm $eig/*ld_prune*qsub*
rm $eig/*high_ld_regions*
rm $eig/*nosex


## Rename BIM/FAM file IDs ##
for an in $anlist;do
    # Rename FAM file IDs
    awk '{$1="ID_"NR; $2="ID_"NR; print $0}' $eig/$study.$an.genotypes.ld_pruned.fam \
        > $eig/$study.$an.genotypes.ld_pruned_renamed.fam

    ## Rename BIM file IDs ##
    awk '{$2="ID_"NR; print $0}' $eig/$study.$an.genotypes.ld_pruned.bim \
        > $eig/$study.$an.genotypes.ld_pruned_renamed.bim
done

## Assign Case/Control to FAM file if applicable ##
if [ $case_control == "true" ]; then
    for an in $anlist; do
        pfile=$procD/$study_${an}_phenotype_for_pca.txt

        ## extract the IDs of controls
        tail -n +2 $pfile |\
            awk '{ if ($2==0) {print $1}}' > $eig/$an.controls.ids

        # Modify FAM file to include case/control status
        awk 'NR==FNR{map[$1] = $1; next} {if($2 in map) {$6=1} else {$6=2}} {print $0} ' \
            $eig/$an.controls.ids $eig/$study.$an.genotypes.ld_pruned_renamed.fam >\
            $eig/$study.$an.genotypes.ld_pruned_renamed_case_control.fam
    done
fi

## run eigenstrat ##
for an in $anlist; do
    if [ $case_control == "true" ]; then
        famfile="$eig/$study.$an.genotypes.ld_pruned_renamed_case_control.fam"
    else
        famfile="$eig/$study.$an.genotypes.ld_pruned_renamed.fam"
    fi
    bimfile="$eig/$study.$an.genotypes.ld_pruned_renamed.bim"
    bedfile="$eig/$study.$an.genotypes.ld_pruned.bed"

    /shared/bioinformatics/software/third_party/EIG-6.1.4/bin/smartpca.perl \
        -i $bedfile \
        -a $bimfile \
        -b $famfile \
        -o $eig/results/${an}_ld_pruned.pca \
        -p $eig/results/${an}_ld_pruned.plot \
        -e $eig/results/${ancestry}_ld_pruned.eval \
        -l $eig/results/${ancestry}_ld_pruned.pca.log \
        -m 0
done

## Extract eigenvectors (top 10 PCs) ##
for an in $anlist; do
    echo "FID IID EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10" > $eig/results/${an}_ld_pruned_top10_eigenvecs.txt
    tail -n +2 $eig/results/${an}_ld_pruned.pca.evec | \
        perl -lne 's/:/ /; @F=split; print join(" ",$F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[11]);' \
        >> $eig/results/${an}_ld_pruned_top10_eigenvecs.txt
done


## Next produce PCA plots