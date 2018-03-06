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
    def __init__(self, n_iteration=100, alpha=0.05, n_sample=-1):
        self.n_iteration = n_iteration
        self.alpha = alpha
        self.n_sample = n_sample

    def _rmsf(self, atomgroup, frames):
        """ Compute the RMSF on the selected frames"""
        means = np.zeros((len(atomgroup.atoms), 3))
        sumsq = np.zeros_like(means)

        for k, ts in enumerate(atomgroup.universe.trajectory[frames]):
            sumsq += (k/(k + 1.0)) * (atomgroup.positions - means)**2
            means[:] = (k*means + atomgroup.positions)/(k + 1.0)
        
        rmsf = np.sqrt(sumsq.sum(axis=1)/(k + 1.0))

        return rmsf

    def _confidence_interval(self, sample, bootstrap, alpha=0.05):
        """ Compute the pivot confidence interval with alpha"""
        s = np.sort((bootstrap.transpose() - sample).transpose())
        low = sample + np.percentile(s, 100 * (alpha), axis=1)
        high = sample + np.percentile(s, 100 * (1-alpha), axis=1)
        return low, high

    def run(self, atomgroup):
        """ Perform the analysis"""
        trj_size = len(atomgroup.universe.trajectory)

        if self.n_sample == -1:
            frame_sample = np.arange(trj_size)
        else:
            frame_sample = np.random.choice(trj_size, self.n_sample)

        # We compute the RMSF on the whole trajectory 
        # This is empirical distribution F* from the true distribution F
        rmsf_sample = self._rmsf(atomgroup, frame_sample)

        # Do the resampling
        rmsf_bootstrap = np.zeros(shape=(len(atomgroup.atoms), self.n_iteration))
        # Generate the empirical boostrap samples
        for i in range(self.n_iteration):
            frame_bootstrap = np.random.choice(frame_sample, frame_sample.shape[0], replace=True)
            rmsf_bootstrap[:,i] = self._rmsf(atomgroup, frame_bootstrap)

        # Get the confidence interval
        rmsf_low, rmsf_high = self._confidence_interval(rmsf_sample, rmsf_bootstrap, self.alpha)

        # Create dataframe
        data = {'resid': atomgroup.resids, 'resname': atomgroup.resnames, 
                'segid': atomgroup.segids, 'rmsf': rmsf_sample,
                'low': rmsf_low, 'high': rmsf_high}
        columns=['resid', 'resname', 'segid', 'rmsf', 'low', 'high']
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
