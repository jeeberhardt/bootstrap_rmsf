#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jérôme Eberhardt 2018
# Bootstrap RMSF
# Author: Jérôme Eberhardt <qksonoe@gmail.com>
#
# License: MIT

import matplotlib.pyplot as plt


def plot_rmsf(fig_name, df, dssp_file=None, ymax=None, start_resid=0):
    color = 'dodgerblue'

    if start_resid != 0:
        start_resid -= 1

    if ymax == None:
        ymax = df['rmsf'].max() + 1.

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(df['resid'] + start_resid, df['rmsf'], color=color, linewidth=0.5)
    plt.fill_between(df['resid'] + start_resid, df['low'], df['high'], color=color, alpha=0.5)

    ax.set_ylim(0, ymax)

    plt.savefig(fig_name, dpi=300, format='png', bbox_inches='tight')