from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from multiprocessing import get_context
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
import time

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

def get_values(sed, nnf, cmpd):
    giic = GIICalculator()
    structures = []
    for structure in sed[cmpd]['structures']:
        neighbors = []
        s = Structure.from_dict(structure)
        gii = giic.GII(s)
        for i in range(len(s)):
            nn_info = nnf.get_nn_info(s, i)
            site_neighbors = [nn_dict['site'] for nn_dict in nn_info]
            neighbors.append(site_neighbors)
        s.add_site_property('neighbors', neighbors)
        structures.append(s)
    return (structures, sed[cmpd]['energies']), giic.params_dict

def pool_map(cmpds, sed, nnf, nprocs):
    partial_values = functools.partial(get_values, sed, nnf)
    with get_context('spawn').Pool(processes=nprocs) as pool:
        data, params = zip(*pool.map(partial_values, cmpds))
    data_dict = {cmpd: {'structures': sed_tup[0], 'energies': sed_tup[1]} for cmpd, sed_tup in zip(cmpds, data)}
    params_dict = merge_dcts(list(params))    
    return data_dict, params_dict

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

if __name__ == '__main__':
    args = argument_parser()

    ### Read-in the structures and energies dictionary ###
    print('Loading dict...\n')
    sed = get_data(args.structure_energy_dictionary)
    cmpds = list(sed.keys())[0:1]

    ### Get the structures + energies with neighbors and choose starting dictionary ###
    print('Finding neighbors and choosing starting dictionary...\n')
    nnf = CrystalNN()
    pmg_sed, params_dict = pool_map(cmpds, sed, nnf, mp.cpu_count())
    print('Starting parameters:')
    print(params_dict)
    print()

    ### Parameterize using starting dictionary ###
    print('Parameterizing...\n')
    bvmp = BVMParameterizer(pmg_sed, params_dict)
    new_params = bvmp.optimizer()
    print('Optimized parameters:')
    print(new_params)
    print()
