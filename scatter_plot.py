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
import logging
import matplotlib.pyplot as plt

from uproot.interpretation import library
LOGGER = logging.getLogger(__name__)

plt.rcParams.update({'font.size': 18})

#signal_files=["dh50", "dh90", "dh130", "dh150", "dh130_zp3500", "dh130_zp500"]
# signal_files=["zp1500_dh90_p4062", "zp1500_dh90_p4434", "zp1500_dh90_p4434_UFO"]
signal_files=["zp1500_dh90_p4434", "zp1500_dh90_p4434_UFO"]
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
    if __name__ == "__main__":
        parser.add_argument("--loggers", nargs="*", default=[__name__], help="Changes the logging level of all the given loggers. 'root' is the global logger and __name__  is logger of this script")
        parser.add_argument("--logging-level", default="warning", choices=["notset", "debug", "info", "warning", "error", "critical"], help="Logging level")
        parser.add_argument("--logging-file", help="Logging file name")
    return parser
# fmt: on

def get_data_frame(file_name):
    file = uproot.open(file_name)
    tree = file["MonoH_Nominal;1"]
    # df = tree.arrays(["DeltaR_bb", "MetTST_met", "N_Jets10", "pt_J_corr","m_J_corr"], library="pd")
    df = tree.arrays(["MetTST_met", "N_Jets10", "pt_J_corr","m_J_corr"], library="pd")
    #df["2mass_over_pt"]=df.apply(lambda row: 2*row["m_J_corr"]/row["pt_J_corr"], axis=1)
    # df = df[df["N_Jets10"]>0]
    #df = df[df["m_J_corr"]>5e4]
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

def calc_efficiency(df, df_interval_var="DeltaR_bb"):
    bin_width=(binnings[df_interval_var][-1]-binnings[df_interval_var][0])/n_bins
    x_intervals=[(binnings[df_interval_var][0]+bin_width * i, binnings[df_interval_var][0]+bin_width * (i+1)) for i in range(n_bins)]
    x = []
    y = []
    for interval in x_intervals:
        interval_mean=sum(interval)/2
        tmp_df=df[(df[df_interval_var]>=interval[0]) & (df[df_interval_var]<interval[1])]
        tot_events=tmp_df["N_Jets10"].count()
        with_selection_events=tmp_df[tmp_df["N_Jets10"]>0]["N_Jets10"].count()
        x.append(interval_mean)
        if tot_events:
            y.append(float(with_selection_events)/float(tot_events))
        else:
            y.append(0)
    return x,y


n_bins=20
binnings={
    # "MetTST_met":(0,1e6),
    "MetTST_met":(1e5,5e5),
    "pt_J_corr":(0,1e6),
    "m_J_corr":(5e4,14e4),
    "2mass_over_pt":(0,1.5),
    "DeltaR_bb":(0,4),
    "eff":(0,1.1)
}

def make_hist_plot(signal_dfs, x_var):
    fig,ax = plt.subplots()
    ax.set_xlim([binnings[x_var][0], binnings[x_var][-1]])
    ax.set_ylim(0, 300)
    for samp_name, df in signal_dfs.items():
        # df.plot.hist(bins=40, ax=ax, column="m_J_corr")
    #     ax.plot(x,y,color=signal_colors[sf], label=signal_labels[sf], linewidth=3, linestyle=signal_style[sf])
        df["m_J_corr"].hist(bins=200, ax=ax, alpha=0.5,color=signal_colors[samp_name], label=signal_labels[samp_name])
    ax.set_xlabel(x_var)
    ax.set_ylabel("# Events")
    ax.legend()
    ax.grid(True)
    plt.show()

def main(args):
    signal_dfs={}
    x_var="MetTST_met"
    # x_var="m_J_corr"
    #x_var="2mass_over_pt"
    # x_var="DeltaR_bb"
    # x_var="pt_J_corr"
    # y_var="DeltaR_bb"
    y_var="eff"
    #y_var="pt_J_corr"

    ax=None
    fig,ax = plt.subplots()

    for sf in signal_files:
        # signal_dfs[sf] = get_data_frame(f"samples/dh{sf}.root")
        signal_dfs[sf] = get_data_frame(f"samples/{sf}.root")


        # if ax:
        #     signal_dfs[sf].plot(kind='scatter', x=x_var, y=y_var, color=signal_colors[sf], ax=ax, alpha=0.2)
        # else:
        #     ax=signal_dfs[sf].plot(kind='scatter', x=x_var, y=y_var, color=signal_colors[sf], alpha=0.2)
        # x,y=calc_average(signal_dfs[sf], df_interval_var=x_var, df_mean_var=y_var)
        # ax.plot(x,y,color=signal_colors[sf], label=signal_labels[sf], linewidth=3, linestyle=signal_style[sf])


        x,y=calc_efficiency(signal_dfs[sf], df_interval_var=x_var)
        ax.plot(x,y,color=signal_colors[sf], label=signal_labels[sf], linewidth=3, linestyle=signal_style[sf])




    ax.set_xlim([binnings[x_var][0], binnings[x_var][-1]])
    ax.set_xlabel(x_var)
    ax.set_ylim([binnings[y_var][0], binnings[y_var][-1]])
    ax.set_ylabel(y_var)
    ax.legend()
    ax.grid(True)

    # plt.show()

    make_hist_plot(signal_dfs, x_var)


    


if __name__ == "__main__":
    parser = comandline_argument_parser()
    command_line_arguments = parser.parse_args()
    logging.basicConfig(
        filename=command_line_arguments.logging_file,
        format="%(levelname)s [%(filename)s:%(lineno)s - %(funcName)s() ]: %(message)s",
    )
    logLevel = getattr(logging, command_line_arguments.logging_level.upper())
    for logger in command_line_arguments.loggers:
        if logger == "root":
            tmpLOGGER = logging.getLogger()
        else:
            tmpLOGGER = logging.getLogger(logger)
        tmpLOGGER.setLevel(logLevel)
    LOGGER.info(command_line_arguments)
    sys.exit(main(command_line_arguments))
