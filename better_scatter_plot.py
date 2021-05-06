#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description of script
"""

import sys
import argparse
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from uproot.interpretation import library

plt.rcParams.update({'font.size': 18})

signal_files=["dh50", "dh90", "dh130", "dh150", "dh130_zp3500", "dh130_zp500"]
# signal_files=["zp1500_dh90_p4062", "zp1500_dh90_p4434", "zp1500_dh90_p4434_UFO"]
# signal_files=["zp1500_dh90_p4434", "zp1500_dh90_p4434_UFO"]
signal_labels={
    "dh50":"m(Z')=1500 GeV, m(DM)=200 GeV, m(S)=50 GeV",
    "dh90":"m(Z')=1500 GeV, m(DM)=200 GeV, m(S)=90 GeV",
    "dh130":"m(Z')=1500 GeV, m(DM)=200 GeV, m(S)=130 GeV",
    "dh150":"m(Z')=1500 GeV, m(DM)=200 GeV, m(S)=150 GeV",
    "dh130_zp500":"m(Z')=500 GeV, m(DM)=200 GeV, m(S)=130 GeV",
    "dh130_zp3500":"m(Z')=3500 GeV, m(DM)=200 GeV, m(S)=130 GeV",
    "zp1500_dh90_p4062":"zp1500_dh90_dm200, p4062 LCTopo",
    "zp1500_dh90_p4434":"zp1500_dh90_dm200, p4434 LCTopo",
    "zp1500_dh90_p4434_UFO":"zp1500_dh90_dm200, p4434 UFO"
}
signal_colors={
    "dh50":"r",
    "dh90":"g",
    "dh130":"b",
    "dh150":"m",
    "dh130_zp500":"k",
    "dh130_zp3500":"c",
    "zp1500_dh90_p4062":"r",
    "zp1500_dh90_p4434":"g",
    "zp1500_dh90_p4434_UFO":"b"
}
signal_style={
    "dh50":"-",
    "dh90":"-",
    "dh130":"--",
    "dh150":"-",
    "dh130_zp500":"--",
    "dh130_zp3500":"--",
    "zp1500_dh90_p4062":"-",
    "zp1500_dh90_p4434":"-",
    "zp1500_dh90_p4434_UFO":"-"
}
# fmt: off
def comandline_argument_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    return parser
# fmt: on

def get_data_frame(file_name):
    file = uproot.open(file_name)
    tree = file["MonoH_Nominal;1"]
    df = tree.arrays(["DeltaR_bb", "MetTST_met", "N_Jets10", "pt_J_corr","m_J_corr"], library="pd")
    df = df[df["N_Jets10"]>0]
    df["2mass_over_pt"]=df.apply(lambda row: 2*row["m_J_corr"]/row["pt_J_corr"], axis=1)
    df = df[df["MetTST_met"]>1.5e5]
    return df

def calc_average(df, df_interval_var="MetTST_met", df_mean_var="DeltaR_bb"):
    bin_width=(binnings[df_interval_var][-1]-binnings[df_interval_var][0])/n_bins
    x_intervals=[(binnings[df_interval_var][0]+bin_width * i, binnings[df_interval_var][0]+bin_width * (i+1)) for i in range(n_bins)]
    x = []
    y = []
    for interval in x_intervals:
        interval_mean=sum(interval)/2
        interval_value=df[(df[df_interval_var]>=interval[0]) & (df[df_interval_var]<interval[1])][df_mean_var].median()
        x.append(interval_mean)
        y.append(interval_value)
    return x,y

n_bins=20
binnings={
    "2mass_over_pt":(0,1.5),
    "DeltaR_bb":(0,4),
}

def main(args):
    signal_dfs={}
    x_var="2mass_over_pt"
    y_var="DeltaR_bb"

    ax=None

    for sf in signal_files:
        fig,ax = plt.subplots()
        signal_dfs[sf] = get_data_frame(f"samples/{sf}.root")
        if ax:
            signal_dfs[sf].plot(kind='scatter', x=x_var, y=y_var, color=signal_colors[sf], ax=ax, alpha=0.2)
        else:
            ax=signal_dfs[sf].plot(kind='scatter', x=x_var, y=y_var, color=signal_colors[sf], alpha=0.2)
        x,y=calc_average(signal_dfs[sf], df_interval_var=x_var, df_mean_var=y_var)
        ax.plot(x,y,color=signal_colors[sf], label=signal_labels[sf], linewidth=3, linestyle=signal_style[sf])

        ax.set_xlim([binnings[x_var][0], binnings[x_var][-1]])
        ax.set_xlabel(x_var)
        ax.set_ylim([binnings[y_var][0], binnings[y_var][-1]])
        ax.set_ylabel(y_var)
        ax.legend()
        ax.grid(True)
        fig.savefig("plot_{}.png".format(sf.replace(".root","")))




    


if __name__ == "__main__":
    parser = comandline_argument_parser()
    command_line_arguments = parser.parse_args()
    sys.exit(main(command_line_arguments))
