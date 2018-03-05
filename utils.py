#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jérôme Eberhardt 2018
# Bootstrap RMSF
# Author: Jérôme Eberhardt <qksonoe@gmail.com>
#
# License: MIT

import matplotlib.pyplot as plt


def plot_rmsf(fig_name, df, dssp_file=None, ymax=None):
    color = 'dodgerblue'

    if ymax == None:
        ymax = df['mean'].max() + 1.

    std_neg = df['mean'] - df['std']
    std_neg[std_neg < 0] = 0
    std_pos = df['mean'] + df['std']

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(df['resid'], df['mean'], color=color, linewidth=1)
    plt.fill_between(df['resid'], std_neg, std_pos, color=color, alpha=0.5)

    ax.set_ylim(0, ymax)

    plt.savefig(fig_name, dpi=300, format='png', bbox_inches='tight')