"""
Generate a FUMA file from GWAS results.
"""

import sys, os, operator

# function parses input arguments

# required arguments: 
# --pval
# 

def parseArguments(args):
    while(len(args) > 0):
        
        # input file name (full path)
        if(args[0] == "--in_file"):
                in_file = args[1]
                args = args[2:]

        # Name of the chromosome column
        elif(args[0] == "--chromosome"):
                chrom = args[1]
                args = args[1:]

        # name of genetic position column
        elif(args[0] == "--position"):
                position = args[1]
                args = args[1:]

        # name of marker name (rsID) column
        elif(args[0] == "--marker_name"):
                marker_name = args[1]
                args = args[1:]

        # name of P-value column
        elif(args[0] == "--p_value"):
                pval = args[1]
                args = args[2:]

        # name of a1 (effect) allele column
        elif(args[0] == "--a1"):
                a1 = args[1]
                args = args[1:]

        # name of a2 (other) allele column
        elif(args[0] == "--a2"):
                a2 = args[1]
                args = args[1:]

        # name of effect (beta) column
        elif(args[0] == "--effect"):
                beta = args[1]
                args = args[1:]

        # name of standard error column
        elif(args[0] == "--std_err"):
                std_err = args[1]
                args = args[1:]

        # FUMA file name (full path)
        elif(args[0] == "--out"):
                outFName = args[1]
                args = args[2:]

        # other arguments unused
        else:
            sys.exit("Unused arguments: " + " ".join(args))

    # check if missing required arguments
    if("pval" not in locals() or "marker_name" not in locals() or "chrom" not in locals()):
    #if("pval" not in locals() or "marker_name" not in locals() or "chrom" not in locals() 
       #or "position" not in locals() or "in_file" not in locals() or "out_file" not in locals()):
        sys.exit("Missing at least one required argument. \
                 Please specify --in_file, --out_file, --p_value, --position, or chromosome")


####################################################################################################
####################################################################################################


