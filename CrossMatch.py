
from mast import *
import pandas

def fetch_short_cad(verbose=False):
    """Fetch list of TESS targets with 2 minute cadence data in sectors 1-7"""
    short_cad = pandas.read_csv('data/TESS_short_cad.csv')
    short_cad_TICs = set(short_cad['TIC'].tolist())
    if verbose: print('TESS short cadence targets:',len(short_cad_TICs))
    return short_cad_TICs

def fetch_K2_targets(target_type='planets',verbose=False):
    K2_targets = pandas.read_csv('data/overlap_'+target_type+'.csv')
    if verbose: print('K2 targets:',K2_targets.shape[0])
    return K2_targets

def cross_match(target_type='planets',short_cad_TICs=None,verbose=False,vals=[]):
    if short_cad_TICs == None:
        short_cad_TICs = fetch_short_cad(verbose=verbose)
    else:
        if verbose: print('Short cadence targets:',len(short_cad_TICs))
    K2_targets = fetch_K2_targets(target_type=target_type,verbose=verbose)
    rows_list = []
    for index, row in K2_targets.iterrows():
        RA = row['RA']
        Dec = row['DEC']
        epic = row['EPIC']
        sector = None
        for i in range(1,27):
            if row['S'+str(i)] == 1:
                sector = i
        try:
            tic = tic_from_coords((float(RA),float(Dec)))[0]
        except ValueError:
            print('Failed on RA:',RA,'Dec:',Dec)
            pass
        if tic in short_cad_TICs:
            val_dict = {'EPIC':epic,'TIC':tic,'RA':RA,'DEC':Dec,'sector':sector}
            for additional_val in vals:
                val_dict[additional_val] = row[additional_val]
            rows_list.append(val_dict)
    matches = pandas.DataFrame(rows_list)
    filename = 'data/crossmatch_'+target_type+'.csv'
    matches.to_csv(filename)
    if verbose: print(matches.shape[0],'matches saved to',filename)


