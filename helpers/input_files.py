#!/usr/bin/env python

from mpi4py import MPI
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Specie
from pparBVM import GIICalculator
from pparBVM import BVMParameterizer
from pymatgen.analysis.local_env import CrystalNN
import argparse
import json

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-rse', '--read_structures_energies', help='path to .json file with structures and energies', type=str, required=True)
    parser.add_argument(
        '-wse', '--write_structures_energies', help='path to .json file to write structure and energies', type=str, required=True)
    parser.add_argument(
        '-wp', '--write_parameters', help='path to .json file to write starting BVM parameters', type=str, required=False)
    args = parser.parse_args()

    return args

def get_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def write_data(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f)
    return

def se_from_json(dct):
    use_se = {}
    for cmpd in list(dct.keys()):
        use_se[cmpd] = {}
        structures = [Structure.from_dict(s) for s in dct[cmpd]['structures']]
        use_se[cmpd]['structures'] = structures
        use_se[cmpd]['energies'] = dct[cmpd]['energies']
    return use_se

def se_to_json(dct):
    json_se = {}
    for cmpd in list(dct.keys()):
        json_se[cmpd] = {}
        structures = [s.as_dict() for s in dct[cmpd]['structures']]
        json_se[cmpd]['structures'] = structures
        json_se[cmpd]['energies'] = dct[cmpd]['energies']
    return json_se

def params_from_json(params):
    use_params = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    use_params['Cation'] = [Specie.from_dict(c) for c in params['Cation']]
    use_params['Anion'] = [Specie.from_dict(a) for a in params['Anion']]
    use_params['R0'] = params['R0']
    use_params['B'] = params['B']
    return use_params

def params_to_json(params):
    json_params = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    json_params['Cation'] = [c.as_dict() for c in params['Cation']]
    json_params['Anion'] = [a.as_dict() for a in params['Anion']]
    json_params['R0'] = params['R0']
    json_params['B'] = params['B']
    return json_params

def get_values(sed, nnf, cmpd):
    giic = GIICalculator()
    structures = []
    for structure in sed[cmpd]['structures']:
        neighbors = []
        gii = giic.GII(structure)
        for i in range(len(structure)):
            nn_info = nnf.get_nn_info(structure, i)
            site_neighbors = [nn_dict['site'].as_dict() for nn_dict in nn_info] # Make json serializable
            neighbors.append(site_neighbors)
        structure.add_site_property('neighbors', neighbors)
        structures.append(structure)
    return {cmpd: {'structures': structures, 'energies': sed[cmpd]['energies']}},  giic.params_dict

def screen_dcts(lsts):
    new_d = {}
    for d in lsts:
        if d is not None:
            new_d.update(d)
    return new_d

def merge_dcts(lsts):
    params_dict = {'Cation': [], 'Anion': [], 'R0': [], 'B': []}
    pairs = []
    for node_lst in lsts: # Appended from each run
        if node_lst is not []: # Node wasn't needed to run
            for dct in node_lst: # Individual dictionary from each run
                for i in range(len(dct['Cation'])): # length of each dictionary
                    pair = [dct['Cation'][i], dct['Anion'][i]]
                    if pair not in pairs:
                        params_dict['Cation'].append(dct['Cation'][i])
                        params_dict['Anion'].append(dct['Anion'][i])
                        params_dict['R0'].append(dct['R0'][i])
                        params_dict['B'].append(dct['B'][i])
                        pairs.append(pair)
    return params_dict

def get_dct_params(sed, nnf, cmpds):
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        if len(cmpds) > nprocs:
            ave, res = divmod(len(cmpds), nprocs)
            counts = [ave + 1 if p < res else ave for p in range(nprocs)]
            starts = [sum(counts[:p]) for p in range(nprocs)]
            ends = [sum(counts[:p+1]) for p in range(nprocs)]
            split_cmpds = [(starts[p], ends[p]) for p in range(nprocs)]
        else:
            split_cmpds = [(i, i+1) if i < len(cmpds) else (len(cmpds), len(cmpds)) for i in range(nprocs)]
    else:
        split_cmpds = None

    cmpds_idxs = comm.scatter(split_cmpds, root=0)
    params_dcts = []
    se_dicts = {}
    for i in range(cmpds_idxs[0], cmpds_idxs[1]):
        if cmpds_idxs[0] != cmpds_idxs[1]:
            c_sed, c_par = get_values(sed, nnf, cmpds[i])
            params_dcts.append(c_par)
            se_dicts.update(c_sed)
        else:
            params_dcts.append({})
            se_dicts.append({})
    g_params_dcts = comm.gather(params_dcts, root=0)
    g_se_dicts = comm.gather(se_dicts, root=0)
    
    if rank == 0:
        final_params = merge_dcts(g_params_dcts)
        final_dict = screen_dcts(g_se_dicts)
        return final_dict, final_params 

if __name__ == "__main__":
    args = argument_parser()
    sed = get_data(args.read_structures_energies)
    o_sed = se_from_json(sed)
    nnf = CrystalNN()

    cmpds = list(o_sed.keys()) 
    vals = get_dct_params(o_sed, nnf, cmpds)
    
    if vals is not None: # rank 0 processor
        pmg_sed, params_dict = vals[0], vals[1]
        w_pmg_sed = se_to_json(pmg_sed)
        write_data(w_pmg_sed, args.write_structures_energies)
       
        if args.write_parameters is not None:
            w_params_dict = params_to_json(params_dict)
            write_data(w_params_dict, args.write_parameters)
