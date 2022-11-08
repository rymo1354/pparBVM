from concurrent.futures import ProcessPoolExecutor
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN
from pparBVM.calculator import GIICalculator
from pparBVM.parameterizer import BVMParameterizer
from scipy.stats import pearsonr
import functools
import multiprocessing as mp
import numpy as np
from tqdm import tqdm
from copy import deepcopy
import argparse
import json

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--structure_energy_dictionary', help='.json file with structures and energies', type=str, required=True)
    args = parser.parse_args()

    return args

def get_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def get_neighbors(sed, nnf, cmpd):
    structures = []
    for structure in sed[cmpd]['structures']:
        neighbors = []
        s = Structure.from_dict(structure)
        for i in range(len(s)):
            nn_info = nnf.get_nn_info(s, i)
            site_neighbors = [nn_dict['site'] for nn_dict in nn_info]
            neighbors.append(site_neighbors)
        s.add_site_property('neighbors', neighbors)
        structures.append(s)
    return (structures, sed[cmpd]['energies'])

def pool_neighbors_map(cmpds, sed, nnf, nprocs):
    partial_neighbors = functools.partial(get_neighbors, sed, nnf)
    with ProcessPoolExecutor(max_workers=nprocs) as executor:
        data_dict = {cmpd: {'structures': sed_tup[0], 'energies': sed_tup[1]} for cmpd, sed_tup in zip(cmpds, tqdm(executor.map(partial_neighbors, cmpds)))}
    return data_dict

def get_pairs(sed, giic, cmpd):
    pairs = []
    for structure in sed[cmpd]['structures']:
        s = Structure.from_dict(structure)
        gii = giic.GII(s, use_sym=True) 
        params = deepcopy(giic.params_dict)
    return params

def merge_dcts(dcts):
    params_dict = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    pairs = []
    for dct in dcts:
        for i in range(len(dct['Cation'])):
            pair = [dct['Cation'][i], dct['Anion'][i]]
            if pair not in pairs:
                params_dict['Cation'].append(dct['Cation'][i])
                params_dict['Anion'].append(dct['Anion'][i])
                params_dict['R0'].append(dct['R0'][i])
                params_dict['B'].append(dct['B'][i])
                pairs.append(pair)
    return params_dict

def pool_pairs_map(cmpds, sed, giic, nprocs):
    partial_pairs = functools.partial(get_pairs, sed, giic)
    with ProcessPoolExecutor(max_workers=nprocs) as executor:
        dcts = list(tqdm(executor.map(partial_pairs, cmpds)))
    params_dict = merge_dcts(dcts)
    return params_dict

if __name__ == '__main__':
    args = argument_parser()

    ### Read-in the structures and energies dictionary ###
    print('Loading dict...')
    sed = get_data(args.structure_energy_dictionary)
    cmpds = list(sed.keys())[0:1]

    ### Get the structures + energies with neighbors as site properties ###
    print('Finding neighbors...')
    nnf = CrystalNN()
    pmg_sed = pool_neighbors_map(cmpds, sed, nnf, mp.cpu_count())

    ### Get starting parameterization dictionary  ###
    print('Choosing starting parameters...')
    giic = GIICalculator()
    params_dict = pool_pairs_map(cmpds, sed, giic, mp.cpu_count())
    print(params_dict)  
  
    ### Parameterize using starting dictionary ###
    print('Parameterizing...')
    bvmp = BVMParameterizer(pmg_sed, params_dict)
    new_params = bvmp.optimizer()
    print(new_params)
