from copy import deepcopy
from pyOpt import SLSQP, Optimization
from pparBVM.calculator import GIICalculator
import numpy as np
from scipy.stats import pearsonr

class BVMParameterizer():
    def __init__(self, structures_and_energies, starting_parameters):
        self.structures_and_energies = structures_and_energies
        self.cmpds = list(self.structures_and_energies.keys())
        self.starting_parameters = starting_parameters
        self.gs_structures_and_energies = self.get_gs_structures_and_energies()

        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            myrank = comm.Get_rank()
        except:
            raise ImportError('mpi4py is required for parallelization') 

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
        print('Mean GIIGS: %s' % str(np.round(float(val), 3)))
        return val

    def mean_Pearson(self, gii_calculator):
        Pearsons = 0
        for cmpd in self.cmpds:
            giis = [gii_calculator.GII(s) for s in self.structures_and_energies[cmpd]['structures']]
            energies = self.structures_and_energies[cmpd]['energies']
            pearson = pearsonr(giis, energies)[0]
            Pearsons += pearson
        return np.divide(Pearsons, len(self.cmpds))

    def mu_Pearson(self, x, C=0.75):
        params = self.get_params_dict(x)
        GIIcalc = GIICalculator(params_dict=params)
        pearson = self.mean_Pearson(GIIcalc)
        val = C - pearson
        print('Mean pearson: %s' % str(np.round(float(pearson), 3)))
        return val

    def optimizer(self, algo=SLSQP, sens_type='CS', sens_mode='pgc'):
        '''
        for structures i composing unique chemical compositions alpha:
            min. sum(GII_GS)
            s.t. mean(pearson(GII, energy))  <= C, where -1 <= C <= 1
        '''
        def obj_func(x):
            f = self.mu_GIIGS(x)
            g = [self.mu_Pearson(x)]
            fail = 0
            return f, g, fail
        
        def getlastsolution(prob: Optimization):
            new_index = prob.firstavailableindex(prob.getSolSet())
            return prob.getSol(new_index - 1)

        opt_prob = Optimization('GII_GS with Pearson Constraint', obj_func)
        
        for i in range(len(self.starting_parameters['Cation'])):
            opt_prob.addVar('x'+str(i+1), 'c', lower=0.0, upper=4.0, value=self.starting_parameters['R0'][i])
        
        opt_prob.addObj('f')
        opt_prob.addCon('g1', 'i')
        o = algo()
        o.setOption('IPRINT', 1)
        o(opt_prob, sens_type=sens_type) 
        
        res = opt_prob.solution(0)
        print(res)
        vs = res.getVarSet()
        x = [np.round(vs[key].value, 3) for key in vs]
        new_params = deepcopy(self.starting_parameters)
        new_params['R0'] = x
        return new_params
