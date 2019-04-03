from __future__ import print_function
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import pandas
from astropy.time import *
from parameters import *
print('lk version:',lk.__version__)

# Fetching planet information
crossmatch_planets = pandas.read_csv('data/crossmatch_planets.csv')
crossmatch_candidates = pandas.read_csv('data/crossmatch_candidates.csv')
combined = pandas.concat([crossmatch_planets,crossmatch_candidates],sort=True)

"""
Helper functions
"""

def get_star_name(tic=None,epic=None):
    """Return the name of a host star based on its TIC ID."""
    if tic:
        index = combined.index[combined['TIC'] == tic].tolist()
    if epic:
        index = combined.index[combined['EPIC'] == epic].tolist()
    if index ==[]:
        return None
    name = combined['pl_hostname'][index].tolist()[0]
    return name

def get_tic_from_epic(epic):
    index = combined.index[combined['EPIC'] == epic].tolist()
    tic = combined['TIC'][index].tolist()[0]
    return tic

def get_epic_from_tic(tic):
    index = combined.index[combined['TIC'] == tic].tolist()
    epic = combined['EPIC'][index].tolist()[0]
    return epic

def get_tic_from_name(hostname):
    index = combined.index[combined['pl_hostname'] == hostname].tolist()
    tic = combined['TIC'][index].tolist()[0]
    return tic

"""
Transit predictions
"""

def predict_transits(hostname,planet='b',period=None,t0=None,mission='tess',verbose=True):
    if period == None or t0 == None:
        if verbose: print('Fetching planet parameters...')
        period, t0, dur, depth, impact_param, r_planet, r_star = get_planet_params(hostname,planet)

    if verbose: print('Predicting transits...')
    now = Time.now().jd
    latest_transit = t0
    transits = []
    while latest_transit < now:
        transits.append(latest_transit)
        latest_transit += period 

    if verbose: 'Calculating time window...'
    if mission.lower() == 'tess':
        planet_column_name = 'sector'
        file_name = 'data/times_TESS.csv'
        time_column_name = 'sector'
    elif mission.lower() == 'k2':
        planet_column_name = 'k2_campaign_str'
        file_name = 'data/times_K2.csv'
        time_column_name = 'sector'

    index = crossmatch_planets.index[crossmatch_planets['pl_hostname'] == hostname].tolist()
    sector_campaign = crossmatch_planets[planet_column_name][index].tolist()[0]
    times = pandas.read_csv(file_name)
    time_index = times.index[times[time_column_name] == sector_campaign].tolist()
    start_jd = times['start_jd'][time_index].tolist()[0]
    end_jd = times['end_jd'][time_index].tolist()[0]

    overlap_transits = [t for t in transits if t > start_jd and t < end_jd]
    
    if verbose: print('Success')
    return overlap_transits

def predict_multi_planet_transits(hostname,mission='tess',verbose=True):
    params = get_params(hostname)
    transit_dict = {}
    for planet in params:
        p = params[planet][0]
        t0 = params[planet][1]
        transits = predict_transits(hostname=hostname,planet=planet,period=p,t0=t0,verbose=False)
        transit_dict[planet] = transits
    return transit_dict

"""
Light curve generation and plotting
"""

def make_lc(tic=None,epic=None,hostname=None,mission='tess',verbose=True):
    """Return a LightCurve object from TESS or K2 planet.

    Keyword arguments:
    Need one of:
        tic -- TIC ID # of the host star
        epic -- EPIC ID # of the host star
        hostname -- host star name
    mission -- which mission to get the data from (default 'tess')
    verbose -- whether to print status updates (default False)
    """
    if tic == None:
        if epic:
            tic = get_tic_from_epic(epic)
        elif hostname:
            tic = get_tic_from_name(hostname)

    if verbose: print('Fetching TPF...')
    if mission.lower() == 'tess':
        tpf = lk.search_targetpixelfile(tic,mission='tess').download()
    elif mission.lower() == 'k2':
        epic = get_epic_from_tic(tic)
        tpf = lk.search_targetpixelfile(epic,mission='k2').download()
    else:
        raise ValueError('Mission must be \'tess\' or \'k2\'.')

    if verbose: print('Detrending...')
    if mission.lower() == 'tess':
        corrector = lk.PLDCorrector(tpf)
        lc = corrector.correct()
    elif mission.lower() == 'k2':
        lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
        lc.remove_outliers(sigma=6).flatten()    


    if verbose: print('Success')    
    return lc

def plot_lc(lc=None,tic=None,epic=None,hostname=None,mission='tess',transits=None,verbose=True):
    """Plot the light curve of a host star with predicted transit overlaid.

    Keyword arguments:
    Must have one:
        lc -- LightCurve object
        tic -- TIC ID # of the host star
        epic -- EPIC ID # of the host star
        hostname -- host star name
    mission -- which mission to get the data from (default 'tess')
    transit_dict -- predicted transits, of the form {'b':[JD1,JD2]}
    verbose -- whether to print status updates (default False)
    """
    colors = ['dodgerblue','hotpink','darkmagenta','navy','mediumseagreen']

    if lc == None:
        lc = make_lc(tic=tic,epic=epic,hostname=hostname,mission=mission,verbose=verbose)

    if hostname:
        title = hostname
    else:
        if tic == None and epic == None:
            title = 'Unknown'
        elif tic:
            hostname = get_star_name(tic=tic)
        elif epic:
            hostname = get_star_name(epic=epic)

    title = hostname+' '+mission

    if mission.lower() == 'tess':
        offset = 2457000
    elif mission.lower() == 'k2':
        offset = 2454833

    ax = lc.scatter(s=0.1)
    
    # for hardcoded single planet lists
    if isinstance(transits,list):
        for time in transits:
            plt.axvline(x=time-offset, color=colors[0], alpha=0.5, linewidth = 3)
    
    # for multi-planet systems and auto-generating dicts
    else:
        transits = predict_multi_planet_transits(hostname=hostname,mission=mission,verbose=verbose)
    
    if isinstance(transits,dict):
        j = 0
        for p in transits:
            color_index = j % len(colors)
            for time in transits[p]:
                plt.axvline(x=time-offset, color=colors[color_index], alpha=0.5, linewidth = 3, label=p)
            j += 1
    
    plt.title(title) 
    plt.legend()
    plt.savefig('images/'+title+'.png')
    plt.show()


