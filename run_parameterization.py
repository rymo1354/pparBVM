#!/usr/bin/env python

from mpi4py import MPI
from pymatgen.core.structure import Structure
from pymatgen.core.structure import PeriodicSite, PeriodicNeighbor
from pymatgen.core.periodic_table import Specie
from pparBVM import GIICalculator
from pparBVM import BVMParameterizer
from scipy.stats import pearsonr
from copy import deepcopy
import argparse
import json

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-rse', '--read_structures_energies', help='path to .json file with structures and energies', type=str, required=True)
    parser.add_argument(
        '-rp', '--read_parameters', help='path to .json file with BVM parameters', type=str, required=True)
    parser.add_argument(
        '-algo', '--algorithm', help='pyOpt algorithm to use', type=str, required=True)
    parser.add_argument(
        '-kw', '--optimizer_kwargs', help='.json convertible str of pyOpt optimizer kwargs, form \'{"key": "value"}\'', type=json.loads, required=True)
    parser.add_argument(
        '-opt', '--optimizer_options', help='.json convertible str of pyOpt optimizer options, form \'{"key": "value"}\'', type=json.loads, required=True)
    parser.add_argument(
        '-wp', '--write_parameters', help='path to .json file of parameterized bond valence parameters', type=str, required=False)
    args = parser.parse_args()

    return args

def get_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def get_site_neighbors(j_structure):
    ### Get site neighbors from PeriodicNeighbor.as_dict() ###
    all_neighbors = []
    for site in j_structure['sites']:
        neighbors_list = []
        for neighbor in site['properties']['neighbors']:
            ps = PeriodicSite.from_dict(neighbor)
            neighbors_list.append(ps)
        all_neighbors.append(neighbors_list)
    return all_neighbors

def sanitize_structure(j_structure):
    ### Remove site properties so Structure object can be read ### 
    sj_structure = deepcopy(j_structure)
    for site in sj_structure['sites']:
        del site['properties']
    return sj_structure

def get_structures_energies(dct):
    use_se = {}
    for cmpd in list(dct.keys()):
        use_se[cmpd] = {}
        structures = []
        for j_structure in dct[cmpd]['structures']:
            neighbors = get_site_neighbors(j_structure)
            structure = Structure.from_dict(sanitize_structure(j_structure))
            structure.add_site_property('neighbors', neighbors) 
            structures.append(structure)
        use_se[cmpd]['structures'] = structures
        use_se[cmpd]['energies'] = dct[cmpd]['energies']
    return use_se

def count_structures(dct):
    count = 0
    for cmpd in list(dct.keys()):
        for structure in dct[cmpd]['structures']:
            count += 1
    return count

def get_parameters(params):
    use_params = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    use_params['Cation'] = [Specie.from_dict(c) for c in params['Cation']]
    use_params['Anion'] = [Specie.from_dict(a) for a in params['Anion']]
    use_params['R0'] = params['R0']
    use_params['B'] = params['B']
    return use_params

def write_data(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f)
    return

def params_to_json(params):
    json_params = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    json_params['Cation'] = [c.as_dict() for c in params['Cation']]
    json_params['Anion'] = [a.as_dict() for a in params['Anion']]
    json_params['R0'] = params['R0']
    json_params['B'] = params['B']
    return json_params
 
if __name__ == '__main__':
    args = argument_parser()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        ### Read-in the structures and energies dictionary ###
        print('Loading dictonaries...\n', flush=True)
    sed = get_data(args.read_structures_energies)
    osed = get_structures_energies(sed)
    scount = count_structures(osed)

    params = get_data(args.read_parameters)
    oparams = get_parameters(params)

    if rank == 0:
        print('Optimizing %s parameters over %s structures comprising %s compositions\n' % (len(oparams['Cation']), scount, len(osed)), flush=True)
        print('Starting parameters:', flush=True)
        print(oparams, flush=True)
        print(flush=True)

    if rank == 0:
        ### Parameterize using starting dictionaries ###
        print('Parameterizing...', flush=True)
    bvmp = BVMParameterizer(osed, oparams)
    new_params = bvmp.optimizer(algo=args.algorithm, kwargs=args.optimizer_kwargs, options=args.optimizer_options)
    json_params = params_to_json(new_params)
    if args.write_parameters is not None:
        write_data(json_params, args.write_parameters)

    if rank == 0:
        print('Optimized parameters:', flush=True)
        print(new_params, flush=True)
        print(flush=True)
