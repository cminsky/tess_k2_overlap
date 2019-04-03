"""
Thanks to Christina Hedges for the NEXSCI data retrieval function
"""

import pandas as pd
import numpy as np

def retrieve_online_data():
    ''' Obtain a dataframe from NExScI of all planet data and store it in the
    package data directory.
    '''
    NEXSCI_API = 'http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI'
    try:
        planets = pd.read_csv(NEXSCI_API + '?table=planets&select=pl_hostname,pl_letter,'
                'pl_disc,ra,dec,pl_trandep,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,'
                'pl_tranflag,pl_trandur,pl_pnum,pl_k2flag,pl_kepflag,pl_facility,'
                'pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,pl_orblper,st_mass,st_masserr1,'
                'st_masserr2,st_rad,st_raderr1,st_raderr2,st_teff,st_tefferr1,'
                'st_tefferr2,st_optmag,st_j,st_h,pl_imppar,pl_impparerr1,pl_impparerr2', comment='#')
        composite = pd.read_csv(NEXSCI_API + '?table=compositepars&select=fpl_hostname,'
                    'fpl_letter,fpl_smax,fpl_smaxerr1,fpl_smaxerr2,fpl_radj,fpl_radjerr1,'
                    'fpl_radjerr2,fpl_bmassj,fpl_bmassjerr1,fpl_bmassjerr2,fpl_eqt,'
                    'fpl_orbper,fpl_orbpererr1,fpl_orbpererr2,fpl_eccen,'
                    'fpl_eccenerr1,fpl_eccenerr2, ', comment='#')
        composite.columns = ['pl_hostname', 'pl_letter', 'pl_orbsmax', 'pl_orbsmaxerr1',
                             'pl_orbsmaxerr2', 'pl_radj', 'pl_radjerr1', 'pl_radjerr2',
                             'pl_bmassj', 'pl_bmassjerr1', 'pl_bmassjerr2', 'pl_eqt',
                             'pl_orbper', 'pl_orbpererr1', 'pl_orbpererr2', 'pl_eccen',
                             'pl_eccenerr1', 'pl_eccenerr2']
    except:
        raise OnlineRetrievalFailure("Couldn't obtain data from NExScI. Do you have an internet connection?")
    df = pd.merge(left=planets, right=composite, how='left', left_on=['pl_hostname', 'pl_letter'],
         right_on=['pl_hostname', 'pl_letter'])
    return df[df.pl_tranflag == 1]

def save_nexsci_output():    
    df = retrieve_online_data()
    df.to_csv('data/nexsci_output.csv')

def read_nexsci_output(verbose=True):
    try:
        df = pd.read_csv('data/nexsci_output.csv')
    except FileNotFoundError:
        if verbose: print('File not found; retrieving data...')
        save_nexsci_output()
        df = read_nexsci_output()
    if verbose: print('Success')
    return df

def get_number_of_planets(hostname):
    df = read_nexsci_output(verbose=False)
    n = np.asarray(df.loc[df.pl_hostname == hostname][['pl_pnum']])[0]
    return n

def get_planet_params(hostname, planet_letter='b'):
    df = read_nexsci_output(verbose=False)
    period, t0, dur, depth, impact_param, r_planet, r_star = np.asarray(df.loc[(df.pl_hostname == hostname) & 
                                               (df.pl_letter == planet_letter)][['pl_orbper', 'pl_tranmid',
                                                'pl_trandur', 'pl_trandep', 'pl_imppar', 'pl_radj', 'st_rad']])[0]
    return period, t0, dur, depth, impact_param, r_planet, r_star

def define_offsets():
    tess_offset = 2457000
    k2_offset = 2454833
    return tess_offset, k2_offset

def get_params(hostname):
    df = read_nexsci_output(verbose=False)
    planet_num = get_number_of_planets(hostname)[0]
    planet_letters = 'bcdefghijklmnopqrstuvwxyz'
    params_dict = {}
    for i in range(planet_num):
        period, t0, dur, depth, impact_param, r_planet, r_star = get_planet_params(hostname,planet_letter=planet_letters[i])
        params_dict[planet_letters[i]] = [period, t0, dur, depth, impact_param, r_planet, r_star]
    return params_dict

def generate_inputs(period, t0, dur, depth, impact_param, r_planet, r_star):
    r, r_min, r_max = generate_radius_ratio(depth, r_planet, r_star)
    b_test = generate_impact_parameter(impact_param)
    
    return period, t0, dur, depth, impact_param, r_planet, r_star

def generate_radius_ratio(depth, r_planet, r_star ):
    if pd.notna(depth):
        r_ratio = np.sqrt(depth)
        min_radius = r_ratio - .5*r_ratio
        max_radius = r_ratio + .5*r_ratio
    elif (pd.notna(r_planet)) & (pd.notna(r_star)):
        r_ratio = (r_planet*71492) / (r_star*695700)
        min_radius = r_ratio - .5*r_ratio
        max_radius = r_ratio + .5*r_ratio
    else:
        r_ratio = .01
        min_radius=0.001
        max_radius=0.1
        
    return r_ratio, min_radius, max_radius

def generate_impact_parameter(impact_param):
    if pd.notna(impact_param):
        b = impact_param
    else:
        b = 0.5
    
    return b

