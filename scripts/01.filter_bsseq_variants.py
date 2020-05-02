# Filters a vcf file for bsseq
# Input:
# 1. Output file -- vcf/bed
import numpy as np
import pandas as pd
import argparse
import logging
import re

logging.basicConfig(format='%(levelname)s:%(asctime)s: %(name)s:  %(message)s', level=logging.DEBUG)
log = logging.getLogger(__name__)

inOptions = argparse.ArgumentParser(description='filter snps')
inOptions.add_argument("-i", dest="snp_vcf", help="input snpvcf file")
inOptions.add_argument("-o", dest="output_file", default = "filtered_snpvcf.bed", type=str, help="output file")

args = inOptions.parse_args()

log.info("reading input files")

input_files = pd.read_csv( args.snp_vcf, header = None, sep = "\t" )
input_files = input_files.rename({
    0: "CHROM",
    1: "POS", 
    2: "REF", 
    3: "ALT",
    4: "AD", 
    5: "GT"}, axis='columns')
log.info("done!")


filter_c_t = (input_files['REF'] == 'C') & (input_files['ALT'] == 'T') 
# all_filters = (input_files['REF'] == 'T') & (input_files['ALT'] == 'C') 
# all_filters = (input_files['REF'] == 'T') & (input_files['ALT'] == 'C') 
filter_g_a = (input_files['REF'] == 'G') & (input_files['ALT'] == 'A') 
filter_ad = input_files['AD'].str.split(",", expand = True)
filter_ad[filter_ad == "."] = np.nan
filter_ad[filter_ad == "None"] = np.nan
filter_ad = filter_ad.astype(float).sum(axis = 1) > 0
input_files = input_files.iloc[ np.where( ~filter_c_t & ~filter_g_a & filter_ad )[0], ]
log.info("writing output files")
# Chr1    21829   .       T       C
# Chr1    341196  T       C
input_files.to_csv( args.output_file, sep = "\t", header = None, index = None)
# import ipdb; ipdb.set_trace()
log.info("finished!")



