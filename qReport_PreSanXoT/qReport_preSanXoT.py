import yaml
import os
import pandas as pd
import argparse
import logging
import sys


os.chdir(r'D:\ReportAnalysis\qReport_PreSanXoT')
from utils.BinomialSiteListMaker import main as BSLM


def getCols(config):
    return (
        config['pdmCol'], config['pCol'], config['freqCol'], config['gCol'], 
        config['aCol'], config['mCol'], config['qCol'], config['descCol'],
        )

configPath = r"D:\ReportAnalysis\qReport_PreSanXoT\qReport_preSanXot.yaml"
with open(configPath) as file:
    config = yaml.load(file, yaml.FullLoader)

def main(config):
    pass
    
    # main
    pdmCol, pCol, freqCol, gCol, aCol, mCol, qCol, descCol = getCols(config)
    
    df = pd.read_csv(config['infile'], sep='\t', low_memory=False)
    
    if not descCol:
        descCol = '_desc'
        df[descCol] = df[qCol]
    
    df2 = df[[qCol, descCol, pdmCol, pCol, gCol, aCol, mCol, freqCol]].drop_duplicates()
    
    bi, biPivot = BSLM({
        'infile': df2,#df[[pdmCol, pCol, gCol, aCol, mCol, freqCol]].drop_duplicates(),
        'outfile': None,
        'peptidoform_column': pdmCol,
        'peptide_column': pCol,
        'modifcation_column': gCol,
        'modified_residue_column': aCol,
        'modified_position_column': mCol,
        'show_unassigned': True,
        'x': config['x'],
        'peakorph_column': None,
        'scanfreq_column': freqCol,
        'binom': config['binom'],
        'q_thr': config['q_thr'],
        'values_pivot': config['values_pivot']
        })
    
    
    
    # Escribir freqTable y generar qTable...
    outFolder = os.path.join(config['outfolder'])
    if not os.path.exists(outFolder):
        os.makedirs(outFolder, exist_ok=True)
    
    with pd.ExcelWriter(os.path.join(outFolder, 'freqTable.xlsx')) as writer:
        bi.to_excel(writer, sheet_name='Raw', index=False)
        biPivot.to_excel(writer, sheet_name=f'PIVOT-{config["binom"]}-{config["q_thr"]}-{config["values_pivot"]}')
    
    ptm = bi[bi[config['binom']]<config['q_thr']]
    ptm = list(zip(ptm[aCol], ptm[gCol]))
    
    qrep = df2.copy()
    qrep.index = list(zip(qrep[aCol], qrep[gCol]))
    qrep = qrep.loc[ptm, :].reset_index()
    qrep = qrep[['index', descCol, qCol, freqCol]].groupby(['index', descCol, qCol]).agg('sum').reset_index()
    qrep = pd.pivot_table(qrep, freqCol, [qCol, descCol], 'index')
    qrep.columns =[(j[0], j[1], int(k)) for j,k in zip(qrep.columns, qrep.sum(axis=0).values)]
    
    
    qrep.index = [(i[0], i[1], j) for i, j in zip(qrep.index.tolist(), qrep.sum(axis=1).values)]
    
    tmp = df2[df2[mCol].isna()][[qCol, freqCol]].groupby(qCol).agg('sum').rename(columns={freqCol: 'Unassigned'}).join(
        df2[~df2[mCol].isna()][[qCol, freqCol]].groupby(qCol).agg('sum').rename(columns={freqCol: 'PTM_total'}),
        how='outer'
        )
    
    tmp = tmp.loc[list(list(zip(*qrep.index.tolist()))[1])]
    
    qrep.index = [(i[0],i[1], j,k, i[2]) for i, j ,k in zip(qrep.index.tolist(), tmp.Unassigned, tmp.PTM_total)]
    qrep = qrep.loc[sorted(qrep.index, key=lambda x: x[4], reverse=True)]
    
    
    idxs = pd.DataFrame(qrep.index.tolist())
    idxs.columns = pd.MultiIndex.from_tuples([
        ('PSMs', descCol, ''), 
        ('PSMs',qCol, ''), 
        ('PSMs','Unassigend', 'total'), 
        ('PSMs', 'PTMs_all', 'total'), 
        ('PSMs', 'PTMs_sig', 'total')
        ])
    
    body = qrep.copy()
    
    body = body[sorted(body.columns, key=lambda x: x[2], reverse=True)]
    body.columns = pd.MultiIndex.from_tuples(body.columns)
    body = body.reset_index(drop=True)
    
    qrep = idxs.join(body)
    header = qrep.columns.to_frame().T
    qrep.columns = list(range(qrep.columns.shape[0]))
    header.columns = list(range(qrep.columns.shape[0]))
    qrep = pd.concat([header, qrep])
    
    if descCol=='_desc':
        qrep = qrep.iloc[:, 1:]
    
    qrep.to_excel(
        os.path.join(config['outfolder'], 'qReport.xlsx'),
        header=False,
        index=False
        )
    


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