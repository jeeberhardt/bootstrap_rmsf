#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jérôme Eberhardt 2018
# Bootstrap RMSF
# Author: Jérôme Eberhardt <qksonoe@gmail.com>
#
# License: MIT


import numpy as np
import pandas as pd


class Bootstrap_RMSF():
    def __init__(self, n_sample=1000, n_iteration=1000):
        self.n_sample = n_sample
        self.n_iteration = n_iteration

    def _rmsf(self, atomgroup, frames):
        """ Compute the RMSF on the selected frames"""
        means = np.zeros((len(atomgroup.atoms), 3))
        sumsq = np.zeros_like(means)

        for k, ts in enumerate(atomgroup.universe.trajectory[frames]):
            sumsq += (k/(k + 1.0)) * (atomgroup.positions - means)**2
            means[:] = (k*means + atomgroup.positions)/(k + 1.0)
        
        rmsf = np.sqrt(sumsq.sum(axis=1)/(k + 1.0))

        return rmsf

    def run(self, atomgroup):
        """ Perform the analysis"""
        rmsf = np.zeros(shape=(len(atomgroup.atoms), self.n_iteration))
        trj_size = len(atomgroup.universe.trajectory)

        # Do a bootstrap
        for i in range(self.n_iteration):
            frames = np.random.choice(trj_size, self.n_sample, replace=True)
            rmsf[:,i] = self._rmsf(atomgroup, frames)

        rmsf_avg = np.mean(rmsf, axis=1)
        rmsf_std = np.std(rmsf, axis=1)

        # Create dataframe
        data = {'resid': atomgroup.resids, 'resname': atomgroup.resnames, 
                'segid': atomgroup.segids, 'mean': rmsf_avg, 'std': rmsf_std}
        columns=['resid', 'resname', 'segid', 'mean', 'std']
        result = pd.DataFrame(data=data, columns=columns)

        # Groupby residue
        result = result.groupby(['resid', 'resname', 'segid']).mean().reset_index()
        # Sort by resid and segid
        result = result.sort_values(['segid', 'resid'])
        # Store it
        self.result = result

        return result

    def save(self, filename):
        """ Write results into a CSV file"""
        self.result.to_csv(filename, index=False)
