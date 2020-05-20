import os, re, itertools

## view head of dictionaries
#def glance(d, size):
#    return dict(itertools.islice(d.items(), size))

os.chdir("C:/Users/jmarks/OneDrive - Research Triangle Institute/Projects/heroin/ngc/uhs4/phenotype")

phenotype_file = "unprocessed/hiv_all_merged_with_uhs_all_phenotype_data_08282017.csv"
var_list = ["hivwave", "pmmn", "personid", "presid", "date_acq", "gwasdate_acq",
            "gwaspresid", "hivstat", "hiv_status", "wave.x", "wave.y", "ind_id"]

def gp_map(phenotype_file, var_list, hardcode = False):
    with open(phenotype_file) as inF:
        head = inF.readline()
        head = re.findall(r'"(.*?)"', head)

        # initialized dictionaries that capture the variables-of-interest information
        # for all subjects in phenotype file
        cell_dic = {}
        serum_dic = {}
        gwas_dic = {}

        # get indices for each of the 3 IDs to map to in the phenotype file
        serum_index = head.index("serum")
        cell_line_index = head.index("cell_line")
        gwas_index = head.index("gwasserum")

        line = inF.readline()
        count = 1
        while line: # loop through phenotype file

            sl = re.findall(r'"(.*?)"', line) # split up line by entries that are surrounded in quotes
            # get the entries for cell_line, serum, and gwasserum
            celi = sl[cell_line_index]
            seri = sl[serum_index]
            gwai = sl[gwas_index]

            if hardcode == True:
                if count == 140:
                    return_line = line
                    celi = sl[cell_line_index-1]

            # initialize list for key-value in dictionary
            cell_dic[celi] = []
            serum_dic[seri] = []
            gwas_dic[gwai] = []

            # append variable entries to each dictionary
            for var in var_list:
                cell_dic[celi].append(sl[head.index(var)])
                serum_dic[seri].append(sl[head.index(var)])
                gwas_dic[gwai].append(sl[head.index(var)])
            count += 1
            line = inF.readline()
        print("done Jesse")
        return cell_dic, serum_dic, gwas_dic, line

cell_d, serum_d, gwas_d, myline = gp_map(phenotype_file, var_list, hardcode=True)


genotype_file = "processing/troubleshoot/all_genotype_ids.n3469"

with open(genotype_file) as inF:
    line = inF.readline()
    keep_list = []
    while line:
        sl = line.split("_")
        asnum = sl[0]
        hhg = sl[2]

        if hhg in cell_d:
            tmptup = (line.strip(), cell_d[hhg])
            keep_list.append(tmptup)
        elif asnum in serum_d:
            tmptup = (line.strip(), serum_d[asnum])
            keep_list.append(tmptup)
        elif asnum in gwas_d:
            tmptup = (line.strip(), gwas_d[asnum])
            keep_list.append(tmptup)
        else:
            print("Didn't find the following ID:\n")
            print(line)
            print(hhg)
        line = inF.readline()

