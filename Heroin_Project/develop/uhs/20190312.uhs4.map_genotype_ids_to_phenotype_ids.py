### python ###
import itertools, os

#os.chdir("/Users/jmarks/OneDrive - Research Triangle Institute/Projects/heroin/ngc/uhs4/phenotype")
#print(os.getcwd())

base_dir = "/Users/jmarks/OneDrive - Research Triangle Institute/Projects/heroin/ngc/uhs4/phenotype"
date = "20190313"
ha_percent = "8"
match_list = ["viralload_cperml.y", "viralload_cperml.x", "hiv_status", "gwashiv", "hivstat", "hiv"]
for match_var in match_list:
    
    out_dir = "ha_data"
    gen = "{}/master.genotype.ids.n3469".format(base_dir)
    phen = "{}/hiv_all_merged_with_uhs_all_phenotype_data_08282017.csv".format(base_dir)
    out_file = "{}/{}/{}.genotype.to.phenotype.ancestry.{}.ha_{}percent.map".format(base_dir, out_dir, date, match_var, ha_percent)
    ha_ids = "ha.ids.{}".format(ha_percent)

    def glance(d):
        return dict(itertools.islice(d.items(), 3))

    def map_fun(gen, phen, match_var):
        with open(gen) as asF, open(phen) as pF:
            phead = pF.readline().split(",")
            serum_index = phead.index("serum")
            cell_line_index = phead.index("cell_line")
            gwas_index = phead.index("gwasserum")
            ancestry_index = phead.index("ancestry_selfreport")
            hiv_index = phead.index(match_var)

            cell_dic = {}
            serum_dic = {}
            gwas_dic = {}
            line = pF.readline()
            while line:
                sl = line.split(",")
                cell_dic[sl[cell_line_index]] = (phead[cell_line_index], sl[ancestry_index], sl[hiv_index])
                serum_dic[sl[serum_index]] = (phead[serum_index], sl[ancestry_index], sl[hiv_index])
                gwas_dic[sl[gwas_index]] = (phead[gwas_index], sl[ancestry_index], sl[hiv_index])
                line = pF.readline()
            print(glance(cell_dic))
            print(glance(gwas_dic))
        #
            keep_list = []
            sline = asF.readline()
            while sline:
                spl = sline.split()
                if spl[1] in cell_dic:
                    tmptup = (spl[2], cell_dic[spl[1]])
                    keep_list.append(tmptup)
                elif spl[0] in serum_dic:
                    tmptup = (spl[2], serum_dic[spl[0]])
                    keep_list.append(tmptup)
                elif spl[0] in gwas_dic:
                    tmptup = (spl[2], gwas_dic[spl[0]])
                    keep_list.append(tmptup)
                else:
                    print(spl[2])
                sline = asF.readline()


        print(len(keep_list))
        mytup = keep_list[1]
        mytup = (mytup[0],) + mytup[1]

        mapped_ids = [(x[0],) + x[1] for x in keep_list]
        print(mapped_ids[:5])


        out_head = "{}\t{}\t{}\t{}".format("genotype_id", "phenotype_column", "ancestry_selfreport", match_var)
        with open(out_file, 'w') as outF:
            outF.write(out_head + "\n")
            for x in mapped_ids:
                line = "\t".join(str(i) for i in x)
                outF.write(line + "\n")
            print("done")
    map_fun(gen, phen, match_var)
    
    def ha_filter(ha_ids, map_file):
        out_file2 = "{}.ha_only".format(map_file)
        with open(ha_ids) as inF, open(map_file) as mF, open(out_file2, "w") as outF:
            head = mF.readline()
            outF.write(head)
            data_dic = {}
            line = mF.readline()
            while line:
                sl = line.split()
                data_dic[sl[0]] = line
                line = mF.readline()

            line = inF.readline()
            while line:
                sl = line.strip()
                outF.write(data_dic[sl])
                line = inF.readline()

    ha_filter(ha_ids, out_file)
