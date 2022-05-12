### Python3 ###
def main():
    """
    Carryout all the functions necessary to 
    perform the fwGWAS naive method. Please there are simply for variables to
    """ ### Variables to alter
    ###############################################################################################

    # config file (full path)
    config="/shared/jmarks/projects/heroin/fwgwas/snpEff/snpEff.config"
    # jar file (full path)
    jar="/shared/jmarks/projects/heroin/fwgwas/snpEff/snpEff.jar"

    # path to directory containing the GWAS results (should end with a forward # slash)
    gwas_dir = '/shared/jmarks/projects/heroin/fwgwas/gwas_results/'

    # name of the GWAS results file in the directory above (gzip results, extension .gz)
    gwas_file= 'wbc_astle2016.txt.gz'

    # path to directory that the processed data should be placed in (should end with a forward # slash)
    out_dir = '/shared/jmarks/projects/heroin/fwgwas/output/004/v04/'

    snp = 'SNP' # name of variant id in the header
    chrom = 'CHR' # name of chromosome in header
    position = 'POS' # name of position in header
    ref = 'A1' # name of reference allele
    alt = 'A2' # name of alternative allele
    pvalue= 'P' # name of pvalue in header


    ### DO NOT alter below this line
    ###############################################################################################
    rtvf_in = gwas_dir + gwas_file
    rtvf_out = out_dir + gwas_file[:-3] + '.vcf'
    manhattan_table = out_dir + gwas_file[:-3] + '.plot_table'


    # Convert GWAS results to VCF format
    se_inF, plot_table = results_to_vcf_format(in_file=rtvf_in, 
            out_file=rtvf_out,
            manhattan_table=manhattan_table,
            pvalue=pvalue,
            chromosome=chrom,
            snp=snp,
            position=position,
            ref=ref,
            alt=alt)

    se_outF = se_inF[:-4] + '.annotated.vcf'

    ## Run SnpEff to obtain annotations
    ex_inF = snp_eff(base_dir=out_dir, 
                     in_file=se_inF, 
                     out_file=se_outF,
                     config=config,
                     jar=jar)

    ex_outF = ex_inF[:-4] + '-cleaned'

    ## Extract annotations, IDs, and P-values
    ga_inF = extract_ann(in_file=ex_inF,
                         out_file=ex_outF)
    ga_outF = ga_inF + '-grouped'
    

    ## Group variants into the four function annotation categories
    ff_inF = group_annotations(in_file=ga_inF,
                               out_file=ga_outF)

    #ff_inF = out_dir + 'test.txt-ann-cleaned-grouped'
    ff_outF = ff_inF + '-filtered-fw'
    

    ## Filter using the four fw-thresholds
    rs_inF1 = filter_fw(in_file=ff_inF,
                        out_file=ff_outF)
    
    fs_outF = ff_inF + '-filtered-std'
   

    ## Filter using the standard GWAS threshold (5e-8)
    rs_inF2 = filter_standard(in_file=ga_outF, 
                             out_file=fs_outF)

    rs_out = rtvf_out[:-4] + '.results-summary'
    
    ## Get results summary
    novel_vars, lost_vars, duplicated_vars = results_summary(in_file=rs_inF1, in_file2=rs_inF2, out_file=rs_out)


    ## plot for highlighting novel, list, and duplicated variants
    for count, highlightlist in enumerate([novel_vars, lost_vars, duplicated_vars]):
        if count == 0: 
            color = 'green'
        if count == 1: 
            color = 'red'
        if count == 2: 
            color = 'yellow'

        gwas_plot(out_dir=out_dir, out_file=gwas_file, table=plot_table, highlightlist=highlightlist, color=color)

    ## plot for base Manhattan without highlights
    gwas_plot(out_dir=out_dir, out_file=gwas_file, table=plot_table)

############################################################################################################################


# Function1
def results_to_vcf_format(in_file, out_file, manhattan_table, pvalue, chromosome, snp, position, ref, alt):
    """
    The SnpEff software expects, as input, a file in VCF format. 
    This function performs the conversion of the GWAS results 
    to VCF format so that SnpEff can obtain the annotations.
    Also, generate the table needed for the plotting script that 
    creates the Manhattan GWAS plot. 
    
    INPUT:
    in_file - GWAS results file. 
    out_file - Name (and path) of the VCF file to be created.
    
    OUTPUT:
    VCF file
    GWAS plotting table to generate the Manhattan plots.
    """
    import gzip


    try:
        with gzip.open(in_file, 'rt') as inF:
            with open(out_file, 'wt') as outF, open(manhattan_table, 'wt') as tableF:
                line = inF.readline()

                # for GWAS plotting table
                table_header = ('VARIANT_ID', 'CHR', 'POSITION', 'P', 'TYPE')
                table_header = '\t'.join(table_header)
                tableF.write(table_header + '\n')
    
                # for VCF file
                head_line = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
                outF.write('\t'.join(head_line) + '\n')
                split_line = line.split()
     
                pval_index = split_line.index(pvalue)
                chr_index = split_line.index(chromosome)
                rsid_index = split_line.index(snp)
                position_index = split_line.index(position)
                ref_allele_index = split_line.index(ref)
                alt_allele_index = split_line.index(alt)

                # skip the header now
                line = inF.readline()
                while(line):
                    split_line = line.split()
     
                    # VCF
                    f1 = split_line[chr_index]
                    f2 = split_line[position_index]
                    f3 = split_line[rsid_index]
                    f4 = split_line[ref_allele_index]
                    f5 = split_line[alt_allele_index]
                    f6 = '.'
                    f7 = '.'
                    f8 = split_line[pval_index]
     
                    vcf_list = [f1, f2, f3, f4, f5, f6, f7, f8]
                    outF.write('\t'.join(vcf_list) + '\n')

                    # GWAS table
                    unique_id = '{}:{}:{}:{}'.format(f1,f2,f4,f5)

                    acgt = ('A', 'C', 'G', 'T')
                    if (f4 and f5) in acgt:
                        var_type = 'snp'
                    else:
                        var_type = 'indel'

                    line_values = (unique_id, f1, f2, f8, var_type)
                    table_line = '\t'.join(line_values) 
                    tableF.write(table_line + '\n')

                    line = inF.readline()

    except (OSError):
        print("Please gzip your input file.")
    return (out_file, manhattan_table)


## Function2
def snp_eff(base_dir, in_file, out_file, config, jar):
    """
    This function executes the SnpEff software that annotates the
    sequence variants using the Genome Build 37 as the reference.'
    
    INPUT:
    base_dir - path were results should be saved.
    in_file - name (and path) of VCF file for input to SnpEff
    out_file - name of the output annotated VCF.
    
    OUTPUT:
    This function returns a character string of the name and path of out_file.
    """
    
    import subprocess

    ### DO NOT modify these variables
    ###########################################################################
    config_path = config 
    snpEff_path = jar
    snp_eff = base_dir + 'snp-eff.sh'
    # -t for multithreading implies -noStats (speeds process way up)
    command_list = ['java', '-Xmx8g', '-jar', snpEff_path, '-c', config_path, '-v',
                     'GRCh37.75', in_file, '>', out_file]
    command_string = ' '.join(command_list)
    ###########################################################################

    # save command as a bash script
    with open(snp_eff, 'w') as outF:
        message = '#!/usr/bin/bash\n'
        message += command_string
        outF.write(message)

    # execute bash script
    run_command = ['bash', snp_eff]
    subprocess.run(run_command)
    return out_file

## Function3
def extract_ann(in_file, out_file):
    """
    Extract the annotation information from the SnpEff
    results VCF file. The unique ID and GWAS P-value
    for each variant is extracted as well.
    
    INPUT: 
    in_file - The name of the file that was output from SnpEff.
              The file should be in vcf format and have in the 
              INFO field the pval+annotation for each variant.
    out_file - Name of the file for which to save the results 
               of this function to. This file will have the 
               following three fields:
          1. unique ID (CHR:POSTION:A1:A2)
          2. sequence variant annotation (e.g. stop-gain,)
          3. P-value
    OUTPUT: 
    This function returns a character string of the name and path of out_file.
    
    """
    with open(in_file, 'r') as inF:
        with open(out_file, 'w') as outF:
            line = inF.readline()

            while line[0] == '#':
                line = inF.readline()
     
            while(line):
                split_line = line.split()
                unique_id = split_line[0] + ':' + split_line[1] + \
                    ':' + split_line[3] + ':' + split_line[4]
                rsID = split_line[2]
                info_field = split_line[7]
                all_annotations = info_field.split(';')
                pval = all_annotations[0]
                functional_annotations = all_annotations[1].split(',')
                for item in functional_annotations:
                    output = rsID + '\t' + unique_id + '\t' + item.split('|')[1] + '\t' + pval
                    outF.write(output + '\n')
                line = inF.readline()  
    return out_file;

## Function4
def group_annotations(in_file, out_file):
    """
    Categorized the variants into groups based off
    of their annotations.
    
    INPUT:
    in_file - the file name (and path) of the file
              that contains the variant annotations.
    out_file - the file name (and path) of the output file
               that will essentially append the variant group
               to the in_put file.
    
    OUTPUT:
    This function returns a character string of the name and path of out_file.
    """
    with open(in_file, 'r') as inF:
        with open(out_file, 'w') as outF:
            line = inF.readline()
            group_dict = {}
            # add loss-of-function variants
            group_dict['loss-of-function variants'] = ['stop_gained', 
                                                      'stop_lost', 
                                                      'frameshift_variant', 
                                                      'splice_donor_variant', 
                                                      'splice_acceptor_variant', 
                                                      'initiator_codon_variant']
            # add moderate-impact variants
            group_dict['moderate-impact variants'] = ['missense_variant', 
                                                     'inframe_insertion', 
                                                     'splice_region_variant']
            # add low-impact variants
            group_dict['low-impact variants'] = ['synonymous_variant', 
                                                '3_prime_UTR_variant', 
                                                '5_prime_UTR_variant', 
                                                'upstream_gene_variant', 
                                                'downstream_gene_variant']
            ## Group variants based on their categorization.
            while(line):
                split_line = line.split()
                # Search for ann in dict.
                for item in group_dict:
                    if split_line[2] in group_dict[item]:
                        split_line.append(item)
                        break
                # If the variant was not in any of the three groups;
                # then it is categorized as an 'other' variant.
                if len(split_line) == 4:
                    split_line.append('other')
                outF.write('\t'.join(split_line) + '\n')
                line = inF.readline()
    return out_file;

## Function5
def filter_fw(in_file, out_file):
    """
    This function filters the variants by comparing the sequence
    variant P-values obtained from the GWAS results with the 
    functional-weighted threshold.
    
    INPUT:
    in_file - name (and path) of the input file that contains the 
              variant groupings.
    out_file - name (and path) of the output file that will be the
               fw-threshold filtered file.
    
    OUTPUT:
    A file containing functionally weighted significant variants.
    """
    with open(in_file) as inF:

        # keep track of duplicate variants
        fw_dict = {}

        with open(out_file, 'w') as outF:
            line = inF.readline()
            loss_func = 5.5e-7
            mod_impact = 1.1e-7
            low_impact = 1.0e-8
            other = 1.7e-9

            # map group to new threshold. also give each group a rank
            # so that lower priority duplicated SNPs will be filtered out.
            thresh_dict = {'loss-of-function variants': [loss_func, 4],
                           'moderate-impact variants': [mod_impact,3],
                           'low-impact variants': [low_impact,2],
                           'other': [other,1]}
            while(line):
                split_line = line.split('\t')
                unique_id = split_line[1]
                group = split_line[4].strip()
                pval = float(split_line[3])
                if pval <= thresh_dict[group][0]:
                    if unique_id not in fw_dict:
                        # store variant information and its rank
                        fw_dict[unique_id] = [line, thresh_dict[group][1]]
                    else:
                        # if the rank of this duplicated SNP is higher than the that
                        # of the previously recored rank then keep this new line
                        if fw_dict[unique_id][1] < thresh_dict[group][1]:
                            fw_dict[unique_id] = [line, thresh_dict[group][1]]

                line = inF.readline()

            for item in fw_dict:
                line = fw_dict[item][0]
                outF.write(line)
                
    return out_file;

## Function6
def filter_standard(in_file, out_file):
    """
    This function filters the variants by comparing the sequence
    variant pvalues obtained from the GWAS results with the 
    standard GWAS threshold of 5e-8.
         
    INPUT:
    in_file - name (and path) of the input file that contains the 
              variant groupings.
    out_file - name (and path) of the output file that will be the
               standard-threshold filtered file.
    
    OUTPUT:
    This function returns a character string of the name and path of out_file.
    """
    
    # keep track of duplicates
    std_dict = {}

    # give each group a rank so that lower priority duplicated SNPs will be filtered out.
    thresh_dict = {'loss-of-function variants': [ 4],
                   'moderate-impact variants': [3],
                   'low-impact variants': [2],
                   'other': [1]}

    with open(in_file) as inF:
        with open(out_file, 'w') as outF:
            line = inF.readline()
            bf_threshold = 5e-8
            
            while(line):
                split_line = line.split('\t')
                unique_id = split_line[1]
                group = split_line[4].strip()
                pval = float(split_line[3])
                if pval <= bf_threshold:
                    if unique_id not in std_dict:
                        # store variant information and its rank
                        std_dict[unique_id] = [line, thresh_dict[group][0]]
                    else:
                        # if the rank of this duplicated SNP is higher than the that
                        # of the previously recored rank then keep this new line
                        if std_dict[unique_id][1] < thresh_dict[group][0]:
                            std_dict[unique_id] = [line, thresh_dict[group][0]]

                line = inF.readline()

            for item in std_dict:
                line = std_dict[item][0]
                outF.write(line)
    return out_file;


## Function7
def results_summary(in_file, in_file2, out_file):
    """
    The function takes as input the two results files from the threshold filtering
    performed above and outputs the summary statistics. Specifically, the output file 
    will contain two counts-dictionaries. One dict will count the number of variants 
    in each functional group that was exclusive to the fw-thresholded sequence
    variants. The other will be for the variants exclusive to the standard WGS P-value
    thresholded results.
    
    INPUT:
    in_file -  name of the fw-thresholded file
    in_file2 - name of the standard thresholded file
    out_file     - Name of the file to which the summary statistics will be saved.
    
    OUTPUT:
    Nothing is returned from this function, however three files are output. One is
    the out_file provided as input, and then two additional files are created.
    
    additional01 - File assumes the name <out_file>-fw-variants and includes a list of 
                   the variants that were deemed significant only when the fw-thresholds
                   were applied as thresholds of significance.
    additional02 - File assumes the name <out_file>-std-variants and includes a list of
                   the variants that were deemed significant only when the standard
                   threshold (5e-8) was applied.
    """
    import subprocess
    ## compare the two filtered files and print the variants 
    #  Exclusive to the fw-thresholded variants
    bash_command = 'bash -c "comm -23 <(sort {0}) <(sort {1})"'.format(in_file, in_file2)
    bash_commanda = 'bash -c "comm -23 <(sort {0}) <(sort {1}) | cut -f2"'.format(in_file, in_file2)
    fw_exclusive = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout
    fw_exclusivea = subprocess.run(bash_commanda, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout

    # Exlusive to the standard thresholded variants (5e-8)
    bash_command2 = 'bash -c "comm -13 <(sort {0}) <(sort {1})"'.format(in_file, in_file2)
    bash_command2a = 'bash -c "comm -13 <(sort {0}) <(sort {1})| cut -f2"'.format(in_file, in_file2)
    standard_exclusive = subprocess.run(bash_command2, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout
    standard_exclusivea = subprocess.run(bash_command2a, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout

    # Variants deemed siginificant using both thresholding methods
    bash_command3 = 'bash -c "comm -12 <(sort {0}) <(sort {1})"'.format(in_file, in_file2)
    bash_command3a = 'bash -c "comm -12 <(sort {0}) <(sort {1})| cut -f2"'.format(in_file, in_file2)
    both_methods = subprocess.run(bash_command3, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout
    both_methodsa = subprocess.run(bash_command3a, shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout

    with open(out_file, 'w') as outF:
        for count,string in enumerate([fw_exclusive, standard_exclusive, both_methods]):
            counts_dict = {'loss-of-function variants':0,
                           'moderate-impact variants':0,
                           'low-impact variants':0,
                           'other':0}
            if count == 0:
                message = '####\n####Novel variants.\n\n'

            elif count == 1:
                message = '\n\n####\n####Variants not replicated.\n\n'

            else:
               message = '\n\n####\n####Variants that were replicated.\n\n'

            split_lines = string.splitlines()
            for line in split_lines:
                ann = line.split('\t')[4]
                counts_dict[ann] += 1
            outF.write(message)
            dict_sum = str(sum(counts_dict.values()))
            dict_sum = '\nTotal number of elements: {}'.format(dict_sum)
            outF.write(str(counts_dict) + dict_sum + '\n')

    novelfile = out_file+'-novel-variants'
    with open(novelfile, 'w') as fwF:
        fwF.write(fw_exclusivea)

    lostfile = out_file+'-lost-variants'
    with open(lostfile, 'w') as stdF:
        stdF.write(standard_exclusivea)

    dupfile = out_file+'-duplicated-variants'
    with open(dupfile, 'w') as bothF:
        bothF.write(both_methodsa)

    return novelfile, lostfile, dupfile



# Function 8 
def gwas_plot(out_dir, out_file, table, highlightlist="", color=""):
    import subprocess

    if color == 'green':
        plot_script = out_dir + 'plot_script_novel_variants.sh'
        out_file = out_file[:-3] + '.gwas_plots_novel_variants'
    elif color == 'red':
        plot_script = out_dir + 'plot_script_lost_variants.sh'
        out_file = out_file[:-3] + '.gwas_plots_lost_variants'
    elif color == 'yellow':
        plot_script = out_dir + 'plot_script_duplicated_variants.sh'
        out_file = out_file[:-3] + '.gwas_plots_duplicated_variants'
    else:
        plot_script = out_dir + 'plot_script.sh'
        out_file = out_file[:-3] + '.gwas_plots'
        


    # base plot without highlights
    command_list = ['Rscript', '/home/ec2-user/bin/generate_gwas_plots.R',
                     '--in', table,
                     '--in_chromosomes autosomal_nonPAR',
                     '--in_header',
                     '--out', out_file, 
                     '--col_id VARIANT_ID',
                     '--col_chromosome CHR',
                     '--col_position POSITION',
                     '--col_p P',
                     '--col_variant_type TYPE ',
                     '--generate_snp_indel_manhattan_plot',
                     '--manhattan_odd_chr_color gray74',
                     '--manhattan_even_chr_color gray87',
                     '--manhattan_points_cex 1.5']

    # if we need to highlight
    if color != "":
        command_list.extend(['--manhattan_highlight_color', color])
        command_list.extend(['--highlight_list', highlightlist])
                                          

    command_string = ' '.join(command_list)

    with open(plot_script, 'w') as outF:
        message = '#!/usr/bin/bash\n'
        message += command_string
        outF.write(message)

    # execute bash script
    run_command = ['bash', plot_script]
    subprocess.run(run_command)

# End of function definitions
################################################################################
if __name__ == "__main__":
    main()
