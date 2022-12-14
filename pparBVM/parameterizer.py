from mpi4py import MPI
from copy import deepcopy
from pyOpt import Optimization
from pparBVM.calculator import GIICalculator
import numpy as np
from scipy.stats import pearsonr
import importlib
import sys

class BVMParameterizer():
    def __init__(self, structures_and_energies, starting_parameters):
        self.structures_and_energies = structures_and_energies
        self.cmpds = list(self.structures_and_energies.keys())
        self.starting_parameters = starting_parameters
        self.gs_structures_and_energies = self.get_gs_structures_and_energies()
   
    def __evaluator__(self, val):
        try:
            return eval(val)
        except NameError: 
            return val

    def get_gs_structures_and_energies(self):
        gs_structures_and_energies = {}
        for cmpd in self.cmpds:
            structures = self.structures_and_energies[cmpd]['structures']
            energies = self.structures_and_energies[cmpd]['energies']
            min_ind = energies.index(min(energies))
            gs_structures_and_energies[cmpd] = {'structures': [structures[min_ind]], 'energies': [energies[min_ind]]}
        return gs_structures_and_energies

    def get_params_dict(self, x):
        # Only supports R0 parameterization currently
        params_dict = deepcopy(self.starting_parameters)
        params_dict['R0'] = x
        return params_dict

    def mean_GIIGS(self, gii_calculator):
        GS_GIIs = 0
        for cmpd in self.cmpds:
            gii = gii_calculator.GII(self.gs_structures_and_energies[cmpd]['structures'][0])
            GS_GIIs += gii
        return np.divide(GS_GIIs, len(self.cmpds))

    def mu_GIIGS(self, x):
        params = self.get_params_dict(x)
        GIIcalc = GIICalculator(params_dict=params)
        val = self.mean_GIIGS(GIIcalc)
        return val

    def mean_Pearson(self, gii_calculator):
        Pearsons = 0
        for cmpd in self.cmpds:
            giis = [gii_calculator.GII(s) for s in self.structures_and_energies[cmpd]['structures']]
            energies = self.structures_and_energies[cmpd]['energies']
            if len(giis) > 1 and len(energies) > 1: # Pearsons of compositions with > 1 structure 
                pearson = pearsonr(giis, energies)[0]
                Pearsons += pearson
        return np.divide(Pearsons, len(self.cmpds))

    def mu_Pearson(self, x, C=0.75):
        params = self.get_params_dict(x)
        GIIcalc = GIICalculator(params_dict=params)
        pearson = self.mean_Pearson(GIIcalc)
        val = C - pearson
        return val

    def optimization_function(self):
        '''
        General Formulation: 
        for structures i composing unique chemical compositions alpha:
            min. sum(GII_GS)
            s.t. mean(pearson(GII, energy))  <= C, where -1 <= C <= 1

        algo (str): algorithm to be used from the pyOpt package
        kwargs (dct): dictionary of kwargs for the optimizer
        options (dct): dictionary of options for the optimizer

        see http://www.pyopt.org/reference/optimizers.html for supported optimizers, kwargs and options
        '''
        def obj_func(x):
            f = self.mu_GIIGS(x)
            g = [self.mu_Pearson(x)]
            fail = 0
            return f, g, fail

        def getlastsolution(prob: Optimization):
            new_index = prob.firstavailableindex(prob.getSolSet())
            return prob.getSol(new_index - 1)

        ### Instantiate Optimization object ###
        opt_prob = Optimization('GII_GS with Pearson Constraint', obj_func)

        ### Add optimization variables ###
        for i in range(len(self.starting_parameters['Cation'])):
            opt_prob.addVar('x'+str(i+1), 'c', lower=1.0, upper=4.0, value=self.starting_parameters['R0'][i])

        ### Specify objective function and constraint(s) ###
        opt_prob.addObj('f')
        opt_prob.addCon('g1', 'i')
        
        return opt_prob

    def optimizer(self, algo, kwargs=None, options=None):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        ### Import and Initialize Optimizer ###
        try:
            'Note: pyOpt must be on PATH' 
            i = importlib.import_module("pyOpt.py{0}.py{0}".format(algo))
            oi = getattr(i, algo)
            try:
                o = oi(pll_type=kwargs['pll_type']) # To specify function parallelization
                del kwargs['pll_type']
            except (TypeError, KeyError) as e:
                o = oi()
        except ImportError:
            print('Invalid optimizer %s; exiting' % algo)
            sys.exit(1)
        
        ### Set Optimizer options ###
        if options is not None:
            for key in list(options.keys()):
                try:
                    o.setOption(key, self.__evaluator__(options[key])) # Get correct data type
                except OSError:
                    print('%s not valid optimizer option; exiting' % key)
                    sys.exit(1)
        
        ### Set Optimizer **kwargs and optimize ###
        opt_prob = self.optimization_function()
        if kwargs is not None:
            try: 
                kwargs_converted = {key: self.__evaluator__(kwargs[key]) for key in list(kwargs.keys())} # Get correct data type
                o(opt_prob, **kwargs_converted) 
            except TypeError:
                print('Invalid keyword argument for obj_func; exiting')
                sys.exit(1)
        else:
            o(opt_prob)
        
        ### Get Optimizer solution ###
        res = opt_prob.solution(0)
        if rank == 0:
            print(res)
        vs = res.getVarSet()
        x = [np.round(vs[key].value, 3) for key in vs]
        new_params = deepcopy(self.starting_parameters)
        new_params['R0'] = x

        return new_params
