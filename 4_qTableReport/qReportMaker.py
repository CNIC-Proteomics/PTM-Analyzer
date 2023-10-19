# -*- coding: utf-8 -*-

#
# Import Modules
#

import argparse
import logging
import numpy as np
import os
import pandas as pd
import openpyxl
import sys
import yaml
import re

import itertools

from scipy.stats import hypergeom

import multiprocessing


idx = pd.IndexSlice
#
# Local Functions
#

#os.chdir(r"S:\U_Proteomica\UNIDAD\software\MacrosRafa\data\Proteomics\GroupTools\ReportAnalysis\qTableReport")
from utils.BinomialSiteListMaker import main as BSLM

def getColumnNames(config, contrast):
    return [tuple(i) for i in [
        config['pdmCol'], 
        config['qCol'], 
        config['qDescCol'], 
        config['pdmFreq'],
        config['qFreq'], 
        [f"{config['sign'][0]}_{contrast}", config['sign'][1]], 
        [f"{config['signNM'][0]}_{contrast}", config['signNM'][1]], 
        [f"{config['FDR_dNM'][0]}_{contrast}", config['FDR_dNM'][1]],
        [f"{config['FDR_NM'][0]}_{contrast}", config['FDR_NM'][1]],
        ]
    ]

def generateFreqTable(config, sign_i, fdr_i, rep, contrast):
    '''
    Parameters
    ----------
    sign_i : TYPE
        DESCRIPTION.
    fdr_i : TYPE
        DESCRIPTION.

    Returns
    -------
    ptm : TYPE
        DESCRIPTION.
    '''
    pdmCol, qCol, qdCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM = getColumnNames(config, contrast)
    
    boolean = np.logical_and.reduce([
        rep[FDRdNM] < fdr_i,
        rep[sign] > 0 if sign_i == 'up' else rep[sign] < 0    
        ])
    
    rep_i = rep[boolean]
    
    rep_i = rep_i[[pdmCol, pdmFreq]].droplevel(1, axis=1)
    
    # If no pdm is filtered return empty list
    if rep_i.shape[0] == 0:
        return []
    
    bi, biPivot = BSLM({
        'infile': rep_i,
        'outfile': None,
        'peptidoform_column': pdmCol[0],
        'x': config['x'],
        'peakorph_column': None,
        'scanfreq_column': pdmFreq[0],
        'binom': config['binom'],
        'q_thr': config['q_thr'],
        'values_pivot': config['values_pivot']
        })
    
    with pd.ExcelWriter(os.path.join(config['outfolder'], 'freqTables', contrast, f'freqTable_{sign_i}_{fdr_i}.xlsx')) as writer:
        bi.to_excel(writer, sheet_name='Raw', index=False)
        biPivot.to_excel(writer, sheet_name=f'PIVOT-{config["binom"]}-FDR-{config["q_thr"]}-{config["values_pivot"]}')
    
    ptm = bi[bi[config['binom']]<config['q_thr']]
    ptm = list(zip(ptm.a, ptm.d))
    return ptm



def qReportPivot(config, fdr_i, sign_i, rep, ptmCol, contrast):
    '''

    Parameters
    ----------
    fdr_i : TYPE
        DESCRIPTION.
    sign_i : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    pdmCol, qCol, qdCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM = getColumnNames(config, contrast)
    
    # Generate PTM Freq Table
    biptm = generateFreqTable(config, sign_i, fdr_i, rep, contrast)
    
    repD = rep[np.logical_and.reduce([
        rep[FDRdNM] < fdr_i,
        rep[sign] > 0 if sign_i == 'up' else rep[sign] < 0,
        np.isin(rep[ptmCol], pd.Series(biptm, dtype='object')),
        ])].sort_values([qCol])
    
    if repD.shape[0] == 0:
        logging.info(f'No PTMs found at FDR = {fdr_i} and Sign = {sign_i}')
        return None
    
    qTableD = pd.pivot_table(
        repD,
        index=[qCol, pdmCol],
        columns=[ptmCol],
        values=[pdmFreq])\
        .droplevel(1, axis=1).fillna(0)
        
    qTableD = {
        'PSMs': qTableD,
        'Peptides': qTableD.applymap(lambda x: 0 if x==0 else 1)
        }
    
    return qTableD
        

def qReportAddData(config, fdr_i, sign_i, quan, qTableD, repNM, rep, contrast):
    '''
    
    Parameters
    ----------
    fdr_i : TYPE
        DESCRIPTION.
    sign_i : TYPE
        DESCRIPTION.
    quan : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    pdmCol, qCol, qdCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM = getColumnNames(config, contrast)
    
    # Collapse PTMs of the same protein
    qTableD = qTableD.reset_index().drop(pdmCol, axis=1).groupby([qCol]).agg(sum)
    
    # Get all PTMs inside each protein 
    qTableD[(quan, 'PTMs')] = qTableD.sum(axis=1)
        
    # Add qFreq and qDesc
    qTableD = qTableD.join(
        rep[[qCol, qFreq, qdCol]]\
            .drop_duplicates().set_index(qCol),
        how='left'
        )
        
    # Add NM freq considering all pdm
    qTableD = qTableD.join(
        repNM[[qCol, pdmFreq]].groupby(qCol).agg(sum).fillna(0)\
            .rename(columns={'pgmFreq': quan, 'REL': 'NM'}),
        how='left'
        )
    
    # Add NM freq of significative pdm changing in the opposite direction
    qTableD = qTableD.join(
        repNM.loc[
            np.logical_and.reduce([
                repNM[signNM]<0 if sign_i=='up' else repNM[signNM]>0,
                repNM[FDRNM]<fdr_i
                ]),
            [qCol, pdmFreq]
            ].groupby(qCol).agg(sum)\
            .rename(columns={'pgmFreq': quan, 'REL': 'NMsig'}),
        how='left'
        ).fillna(0)
    
    # Add Hypergeometric test
    qTableD[('Hypergeom', 'NMbased')] = [
        1-hypergeom.cdf(
            x-1, # x
            qTableD[(quan,'NM')].sum(), # M overall population
            qTableD[(quan,'NMsig')].sum(), # n Defect population
            N # N test population
        )
        for x,N in zip(qTableD[(quan,'NMsig')], qTableD[(quan,'NM')])
        ]

    qTableD[('Hypergeom', 'PTMbased')] = [
        1-hypergeom.cdf(
            x-1, # x
            qTableD[qFreq].sum()-qTableD[(quan,'NM')].sum(), # M overall population
            qTableD[(quan, 'PTMs')].sum(), # n Defect population
            N # N test population
        ) 
        for x,N in zip(
                qTableD[(quan, 'PTMs')], 
                qTableD[qFreq]-qTableD[(quan,'NM')]
                )
        ]
    
    return qTableD


def qReportDesign(config, quan, qTableD, contrast):
    '''

    Parameters
    ----------
    config : TYPE
        DESCRIPTION.
    qTableD : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    pdmCol, qCol, qdCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM = getColumnNames(config, contrast)
    
    # Sort columns
    infoCols = [qdCol,qFreq,(quan,'NM'),(quan,'NMsig'),('Hypergeom','NMbased'),('Hypergeom','PTMbased'), (quan, 'PTMs')]
    i = qTableD.loc[:, infoCols]
    qTableD = i.join(qTableD.loc[:, [pdmFreq[0]]].replace(0, np.nan)).sort_values((quan, 'PTMs'), ascending=False)


    # Add summ row
    sumRow = []
    for column in qTableD.columns:
        if qTableD[column].dtype == 'object':
            sumRow.append('')
            continue
        elif column[0] == 'Hypergeom':
            sumRow.append(np.nan)
            continue
        else:
            sumRow.append(qTableD[column].sum())
            continue

    qTableD = pd.concat([
        pd.DataFrame([sumRow], columns=qTableD.columns, index=['Sum']),
        qTableD
        ])

    # Sort PTMs by total number
    qTableD = qTableD[infoCols].join(
        qTableD[    
            qTableD.loc['Sum', [pdmFreq[0]]].sort_values(ascending=False).index
        ]
    ).reset_index(names=[qCol])

    qTableD.columns = pd.MultiIndex.from_tuples([
        (quan,j[0], j[1]) if i==pdmFreq[0] else (i,j,'') 
        for i,j in qTableD.columns
        ])
    
    
    plotted_q = [os.path.splitext(i)[0] for i in os.listdir(config['plotFolder'])]
    qTableD[qTableD.columns[0]] = \
        [f"=HYPERLINK(\"{os.path.join(config['plotFolder'], i)}.html\", \"{i}\")" if i in plotted_q else i for i in qTableD.iloc[:, 0]]
    
    return qTableD


def qReportWrite(config, fdr_i, sign_i, quan, qTableD, contrast):
    
    if not os.path.exists(os.path.join(config['outfolder'], 'qReport', contrast, f'FDR-{fdr_i}')):
        os.makedirs(os.path.join(config['outfolder'], 'qReport', contrast, f'FDR-{fdr_i}'), exist_ok=True)
    
    qReportPath = os.path.join(config['outfolder'], 'qReport', contrast, f'FDR-{fdr_i}', f'qReport_FDR-{fdr_i}_{sign_i}_{quan}.xlsx')
    
    header = list(zip(*qTableD.columns.tolist()))
    qTableD.columns = np.arange(0, qTableD.shape[1])
    
    qTableD = pd.concat([pd.DataFrame(header), qTableD])
    
    qTableD.to_excel(
        qReportPath,
        header=False,
        index=False#, sep='\t'
        )
    
    toFormat = [n+1 for n,i in enumerate(qTableD.iloc[:, 0]) if 'HYPERLINK' in i]
    
    book = openpyxl.load_workbook(qReportPath)
    sheet = book['Sheet1']
    # sheet.delete_rows(4, 1)
    # sheet.delete_cols(1, 1)
    
    for i in toFormat:
        sheet[f'A{i}'].font = openpyxl.styles.Font(color='0000FF', underline='single')
    
    book.save(qReportPath)
    
    

def qReportContrast(rep, config, contrast):
    '''
    

    Parameters
    ----------
    config : TYPE
        DESCRIPTION.
    contrast : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    pdmCol, qCol, qdCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM = getColumnNames(config, contrast)
    
    # Get required report fraction
    rep = rep.loc[:, list(set([pdmCol, qCol, pdmFreq, qFreq, sign, signNM, FDRdNM, FDRNM, qdCol]))].drop_duplicates()
    
    if config['pdmColFormat']==2:
        logging.info('Formatting pdm from ; to []')
        mypgm = [i.split(';') for i in rep[pdmCol]]
        rep[pdmCol] = [
            f'{i[0]}_{i[1]}' if len(i)==2 or i[2]=='' else f'{i[0][:int(i[2])]}[{i[1]}]{i[0][int(i[2]):]}'
            for i in mypgm
        ]
    
    rep = rep[~rep[pdmCol].duplicated()]
    
    # Get PTM from input report
    logging.info('Get PTMs from input report')
    ptmCol = ('PTM', 'REL')
    myptm = [re.search(r'(.)\[([^]]+)\]', pdm_i) for pdm_i in rep[pdmCol]]
    rep[ptmCol] = [i.groups() if i else (None, None) for i in myptm]
    
    # Extract NM elements from report
    repNM = rep.loc[[i==(None, None) for i in rep[ptmCol]], [qCol, pdmCol, pdmFreq, FDRNM, signNM]]
    repNM = {
        'PSMs': repNM.copy(),
        'Peptides': repNM.copy()
        }
    repNM['Peptides'][pdmFreq] = [0 if i==0 else 1 for i in repNM['Peptides'][pdmFreq]]
    
    
    # Create folder with output files
    if not os.path.exists(os.path.join(config['outfolder'], 'FreqTables', contrast)):
        os.makedirs(os.path.join(config['outfolder'], 'FreqTables', contrast), exist_ok=True)
    
    
    # All combinations FDR x Sign
    fdrxsign = list(itertools.product(config['FDR'], ['up', 'down']))
    
    params = [(config, fdr_i, sign_i, rep, ptmCol, contrast) for fdr_i, sign_i in fdrxsign]
    qReportList = [(i[1], i[2], qReportPivot(*i)) for i in params]
    
    logging.info('Pivot Report to obtain pre-qReport')
    params = [
     (config, fdr_i, sign_i, quan, qTableD[quan], repNM[quan], rep, contrast) 
     for fdr_i, sign_i, qTableD in qReportList if qTableD 
     for quan in qTableD
     ]
    
    # Single core
    # qReportList = [(i[1], i[2], i[3], qReportAddData(*i)) for i in params]

    # Multicore
    pool = multiprocessing.Pool(processes=config['n_cpu'])
    qReportList = pool.starmap(qReportAddData, params)
    pool.close()
    pool.join()
    qReportList = [(i[1], i[2], i[3], j) for i,j in zip(params, qReportList)]
    
    
    logging.info('Adding data to qReport')
    params = [
     [(config, quan, qTableD, contrast), (fdr_i, sign_i, quan)]
     for fdr_i, sign_i, quan, qTableD in qReportList
     ]
    qReportList = [(fdr_i, sign_i, quan, qReportDesign(*i)) for i, (fdr_i, sign_i, quan) in params]
    
    logging.info('Adapting qReport format')
    params = [
     (config, fdr_i, sign_i, quan, qTableD, contrast)
     for fdr_i, sign_i, quan, qTableD in qReportList
     ]
    

    if config['outfolder']:
        logging.info('Writing output')
        _ = [qReportWrite(*i) for i in params]
    
    else:
        return qReportList
    


#
# Main
# 


def main(config, file=None):
    '''
    main
    '''
    
    # Get pandas report from file or read it from config
    if file:
        rep = file.copy()
    else:
        logging.info(f"Reading Report: {config['infile']}")
        rep = pd.read_csv(
            config['infile'], 
            sep='\t', 
            low_memory=False, 
            header=[0,1],
            keep_default_na=True,
            na_values='#DIV/0!'
            )
    
    
    return [qReportContrast(rep, config, contrast) for contrast in config['groups']]
    



if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(
        description='qTableMaker',
        epilog='''
        Example:
            python qTableMaker.py
        ''')

    parser.add_argument('-c', '--config', default=os.path.join(os.path.dirname(__file__), 'qReportMaker.yaml'), type=str, help='Path to config file')

    args = parser.parse_args()


    with open(args.config) as file:
        config = yaml.load(file, yaml.FullLoader)
        

    logging.basicConfig(level=logging.INFO,
                        format='qTableMaker.py - '+str(os.getpid())+' - %(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[logging.FileHandler(
                            os.path.splitext(config['infile'])[0]+'.log'
                            ),
                            logging.StreamHandler()])

    logging.info('Start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(config)
    logging.info('End script')

