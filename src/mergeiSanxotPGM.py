# -*- coding: utf-8 -*-

#
# Import Modules
#

import os
import sys
import argparse
import logging
import pandas as pd


#
# Main
# 

def main(args):
    '''
    main
    '''
    
    # read report
    logging.info(f"Reading iSanXoT report: {args.is_file}...")
    rep = pd.read_csv(args.is_file, sep='\t', header=[0,1], low_memory=False)

    logging.info(f"Reading PDM/PGM table: {args.pdm_file}...")
    pdmt = pd.read_csv(args.pdm_file, sep='\t', low_memory=False)

    logging.info("Providing second header to PDM/PGM table...")
    lvCol = [i for i,j in rep.columns if j == 'LEVEL']
    pdmt.columns = pd.MultiIndex.from_tuples([(i,'LEVEL') if i in lvCol else (i, 'REL') for i in pdmt.columns])

    logging.info("Merging iSanXoT report and PDM/PGM table...")
    rep2 = pd.merge(
        rep,
        pdmt,
        on=[i for i in pdmt.columns if i[1]=='LEVEL'],
        how='left'
    )

    logging.info(f"Writing output report: {args.outfile}")
    rep2.to_csv(args.outfile, sep='\t', index=False)




if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(
        description='mergeiSanXoTPGM',
        epilog='''
        Example:
            python mergeiSanXoTPGM.py
        ''')

    parser.add_argument('-i', '--is-file', required=True, help='Path to iSanXoT report file')
    parser.add_argument('-p', '--pdm-file', required=True, help='Path to PDM/PGM table')
    parser.add_argument('-o', '--outfile', required=True, help='Output with the merged values')

    args = parser.parse_args()


    # prepare workspace
    outdir = os.path.dirname(args.outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=False)


    logging.basicConfig(level=logging.INFO,
                        format='mergeiSanXoTPGM.py - '+str(os.getpid())+' - %(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[logging.FileHandler(
                            os.path.join(outdir, 'mergeiSanXoTPGM.log')
                            ),
                            logging.StreamHandler()])

    logging.info('Start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info(f'End script')

