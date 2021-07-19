#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contour plots from TTrees
"""

from functools import update_wrapper
import os
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from uproot.interpretation import library

# matplotlib.use('qt5agg')
plt.rcParams.update({"font.size": 18})
sample_root_dir="samples/v07UFO"
N_max=2000

class Sample():

    def __init__(self, name, sample_list, topology):
        self.topology=topology
        self.name=name+"_"+self.topology
        self.sample_list=sample_list

    def is_signal(self):
        return True if "mzp" in self.name else False

    def get_data_frame(self, file_name):
        file = uproot.open(file_name)
        tree = file[f"{file_name.split('/')[-1].replace('.conf.root','')}_Nominal;1"]
        if self.topology=="resolved":
            # weight = "Weights*GenWeight*GenWeightMCSampleMerging*EleWeight*MuoWeight*TauWeight*JetWeightJVT*JetWeightBTag*MET_TriggerSF*Znunu_Normalization*Znunu_Merging"
            weight = "GenWeightMCSampleMerging*EleWeight*MuoWeight*TauWeight*MET_TriggerSF*Znunu_Normalization*Znunu_Merging"
            pt_higgs="pt_jj_corr"
        elif self.topology=="merged":
            pt_higgs="pt_J_corr"
            # weight = "Weights*GenWeight*GenWeightMCSampleMerging*EleWeight*MuoWeight*TauWeight*TrackJetWeight*MET_TriggerSF*Znunu_Normalization*Znunu_Merging"
            weight = "GenWeightMCSampleMerging*EleWeight*MuoWeight*TauWeight*TrackJetWeight*MET_TriggerSF*Znunu_Normalization*Znunu_Merging"
        df = tree.arrays( ["Met",pt_higgs,weight], library="pd")
        if self.topology=="resolved":
            df = df[df["Met"]<500]
            # df = df[df["Met"]>500]
        elif self.topology=="merged":
            # df = df[df["Met"]<500]
            df = df[df["Met"]>500]
        df=df.rename(columns={pt_higgs:"pt_higgs", weight:"weight"})
        return df

    def read_sample_list(self):
        self.data=pd.DataFrame()
        for samp in self.sample_list:
            df=self.get_data_frame(samp)
            self.data=self.data.append(df)
        self.data["sample"]=self.name
        if self.is_signal():
            self.data["sample_type"]="signal"
        else:
            self.data["sample_type"]="background"
        self.data=self.data.sample(frac=1)
        self.data=self.data[:N_max]
        if len(self.data)>0:
            self.data["metpt_ratio"]=self.data.apply(lambda row: row["Met"]/row["pt_higgs"], axis=1)
            self.data["metpt_diff"]=self.data.apply(lambda row: row["Met"]-row["pt_higgs"], axis=1)
            # # band=60
            # band=180
            # if self.topology=="merged":
            #     # band=180
            #     band=60
            # # a=1.2
            # # b=40
            # # upper= a*self.data["Met"]+b
            # # lower= self.data["Met"]*(2-a)-b
            # lincomb_eff=len(self.data[(self.data["pt_higgs"]<band+self.data["Met"]) & (self.data["pt_higgs"]>self.data["Met"]-band)]) / len(self.data)
            # metpt_eff=len(self.data[(self.data["metpt"]<1.3) & (self.data["metpt"]>0.8)]) / len(self.data)
            # print(f"{self.name}: {lincomb_eff/metpt_eff:.3f} ({lincomb_eff=:.2f}/{metpt_eff=:.2f})")
            # print("linear comb cut efficiency: ",len(self.data[(self.data["pt_higgs"]<upper) & (self.data["pt_higgs"]>lower)]) / len(self.data))
        # print(f"Loaded {len(self.data)} events from {self.name}!")



def get_samples(base_path, topology):
    Diboson = Sample(
        name="Diboson",
        topology=topology,
        sample_list=[
        "{}/WW.conf.root".format(base_path),
        "{}/WZ.conf.root".format(base_path),
        "{}/ZZ.conf.root".format(base_path)
        ])


    VHbb = Sample(
        name="VHbb",
        topology=topology,
        sample_list=[
        "{}/VHbb_WmH.conf.root".format(base_path),
        "{}/VHbb_WpH.conf.root".format(base_path),
        "{}/VHbb_ZHllbb.conf.root".format(base_path),
        "{}/VHbb_ZHvvbb.conf.root".format(base_path),
        "{}/VHbb_ggZHllbb.conf.root".format(base_path),
        "{}/VHbb_ggZHvvbb.conf.root".format(base_path),
        ])

    Zjets = Sample(
        name="Zjets",
        topology=topology,
        sample_list=[
        "{}/Zee_cl.conf.root".format(base_path),
        "{}/Zee_hf_0_70.conf.root".format(base_path),
        "{}/Zee_hf_140_280.conf.root".format(base_path),
        "{}/Zee_hf_280_500.conf.root".format(base_path),
        "{}/Zee_hf_70_140.conf.root".format(base_path),
        "{}/Zee_hpt_1000p.conf.root".format(base_path),
        "{}/Zee_hpt_500_1000.conf.root".format(base_path),
        "{}/Zee_l.conf.root".format(base_path),
        "{}/Zmumu_cl.conf.root".format(base_path),
        "{}/Zmumu_hf_0_70.conf.root".format(base_path),
        "{}/Zmumu_hf_140_280.conf.root".format(base_path),
        "{}/Zmumu_hf_70_140.conf.root".format(base_path),
        "{}/Zmumu_hpt_1000p.conf.root".format(base_path),
        "{}/Zmumu_hpt_500_1000.conf.root".format(base_path),
        "{}/Zmumu_l.conf.root".format(base_path),
        "{}/Znunu_cl_maxhtptv.conf.root".format(base_path),
        "{}/Znunu_cl_ptv.conf.root".format(base_path),
        "{}/Znunu_hf_maxhtptv_0_70.conf.root".format(base_path),
        "{}/Znunu_hf_maxhtptv_140_280.conf.root".format(base_path),
        "{}/Znunu_hf_maxhtptv_280_500.conf.root".format(base_path),
        "{}/Znunu_hf_maxhtptv_70_140.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_100_140_mjj_0_500.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_100_140_mjj_1000p.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_100_140_mjj_500_1000.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_140_280_mjj_0_500.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_140_280_mjj_1000p.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_140_280_mjj_500_1000.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_280_500.conf.root".format(base_path),
        "{}/Znunu_hf_ptv_70_100.conf.root".format(base_path),
        "{}/Znunu_hpt_mxhtptv.conf.root".format(base_path),
        "{}/Znunu_hpt_ptv.conf.root".format(base_path),
        "{}/Znunu_l.conf.root".format(base_path),
        "{}/Ztautau_cl.conf.root".format(base_path),
        "{}/Ztautau_hf.conf.root".format(base_path),
        "{}/Ztautau_hpt.conf.root".format(base_path),
        "{}/Ztautau_l.conf.root".format(base_path),
        ])

    Wjets = Sample(
        name="Wjets",
        topology=topology,
        sample_list=[
        "{}/Wenu_cl.conf.root".format(base_path),
        "{}/Wenu_hf.conf.root".format(base_path),
        "{}/Wenu_hpt.conf.root".format(base_path),
        "{}/Wenu_l.conf.root".format(base_path),
        "{}/Wmunu_cl.conf.root".format(base_path),
        "{}/Wmunu_hf_0_70.conf.root".format(base_path),
        "{}/Wmunu_hf_140_280.conf.root".format(base_path),
        "{}/Wmunu_hf_280_500.conf.root".format(base_path),
        "{}/Wmunu_hf_70_140.conf.root".format(base_path),
        "{}/Wmunu_hpt.conf.root".format(base_path),
        "{}/Wmunu_l.conf.root".format(base_path),
        "{}/Wtaunu_cl.conf.root".format(base_path),
        "{}/Wtaunu_hf.conf.root".format(base_path),
        "{}/Wtaunu_hpt.conf.root".format(base_path),
        "{}/Wtaunu_l.conf.root".format(base_path),
        ])


    ttbar = Sample(
        name="ttbar",
        topology=topology,
        sample_list=[
        "{}/ttbar_MET_100_200.conf.root".format(base_path),
        "{}/ttbar_MET_200_300.conf.root".format(base_path),
        "{}/ttbar_MET_300_400.conf.root".format(base_path),
        "{}/ttbar_MET_400p.conf.root".format(base_path),
        "{}/ttbar_nonallhad.conf.root".format(base_path),
        ])


    stop = Sample(
        name="stop",
        topology=topology,
        sample_list=[
        "{}/stopWt.conf.root".format(base_path),
        "{}/stops.conf.root".format(base_path),
        "{}/stopt.conf.root".format(base_path),
        ])


    signal_mzp3500_dm200_dh150 = Sample(
        name="mzp3500_dh150",
        topology=topology,
        sample_list=[
        "{}/monoSbb_zp3500_dm200_dh150.conf.root".format(base_path),
        ])


    signal_mzp3500_dm200_dh50 = Sample(
        name="mzp3500_dh50",
        topology=topology,
        sample_list=[
        "{}/monoSbb_zp3500_dm200_dh50.conf.root".format(base_path),
        ])


    signal_mzp500_dm200_dh150 = Sample(
        name="mzp500_dh150",
        topology=topology,
        sample_list=[
        "{}/monoSbb_zp500_dm200_dh150.conf.root".format(base_path),
        ])

    signal_mzp500_dm200_dh50 = Sample(
        name="mzp500_dh50",
        topology=topology,
        sample_list=[
        "{}/monoSbb_zp500_dm200_dh50.conf.root".format(base_path),
        ])

    # return (VHbb, Diboson, stop, ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50, signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    return ( ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50, signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    # return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh50)
    # return (VHbb, Diboson, stop, ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, signal_mzp500_dm200_dh150)
    # return (Zjets, signal_mzp3500_dm200_dh150)
    # return (signal_mzp3500_dm200_dh150, Zjets)
    # return (VHbb,)
    # return (Zjets,)
    # return (signal_mzp3500_dm200_dh50, signal_mzp500_dm200_dh50)
    # return (signal_mzp3500_dm200_dh150,)
    # return (signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50, signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    # return (signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50)
    # return (signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    # return (signal_mzp500_dm200_dh50,)


def make_scatter_plots(data):
    sample_list=data["sample"].unique()
    straight_line=np.linspace(0,1500,10)
    band=100
    straight_lower=np.linspace(0-band,1500-band,10)
    straight_upper=np.linspace(0+band,1500+band,10)
    for n, samp_name in enumerate(sample_list):
        if "merged" in samp_name: continue
        fig, ax = plt.subplots()
        sns.scatterplot(data=data[data["sample"]==samp_name], x="Met", y="pt_higgs", hue="sample", ax=ax)
        sns.scatterplot(data=data[data["sample"]==samp_name.replace("resolved","merged")], x="Met", y="pt_higgs", hue="sample", palette="inferno", ax=ax)
        # sns.contour(data=data[data["sample"]==samp_name], x="Met", y="pt_higgs", hue="sample", ax=ax)
        # sns.contour(data=data[data["sample"]==samp_name.replace("resolved","merged")], x="Met", y="pt_higgs", weights=data["weights"], levels=[0.1], hue="sample", palette="inferno", common_norm=True, ax=ax)
        # sns.contour(data=data, x="Met", y="pt_higgs", weights=data["weights"], levels=[0.1], hue="sample", palette="inferno", common_norm=True, ax=ax)
        ax.plot(straight_line, straight_line, "r-")
        ax.plot(straight_line, straight_upper, "b-")
        ax.plot(straight_line, straight_lower, "b-")
        fig.savefig(f"plots/contour_plots/scatter_plot_{samp_name}.png")

def kde_plot(data, x_var):
    for topo in ["resolved" ,"merged"]:
        fig, ax = plt.subplots(figsize=(13,10))
        tmp_data=data[data["sample"].str.contains(topo)]
        tmp_data=tmp_data[(tmp_data["metpt_ratio"]>0.7) & (tmp_data["metpt_ratio"]<1.3)]
        signals=[]
        if topo=="merged":
            signals=[f"mzp500_dh50_{topo}", f"mzp500_dh150_{topo}"]
        elif topo=="resolved":
            signals=[f"mzp3500_dh50_{topo}", f"mzp3500_dh150_{topo}"]
        tmp_data=tmp_data[~(tmp_data["sample"].isin(signals))]
        # s_data=tmp_data[tmp_data["sample_type"]=="signal"]
        # b_data=tmp_data[tmp_data["sample_type"]=="background"]
        # sns.kdeplot(data=s_data, x=x_var, hue="sample", weights=s_data["weight"], alpha=0.5, ax=ax, palette="tab20", common_norm=False, fill=True)
        sns.kdeplot(data=tmp_data, x=x_var, hue="sample", weights=tmp_data["weight"], alpha=0.5, ax=ax, palette="tab20", common_norm=False, multiple="layer")
        # plt.vlines([0.7, 1.3], ymin=0, ymax=ax.get_ylim()[1], colors="red")
        plt.vlines([-60, 60], ymin=0, ymax=ax.get_ylim()[1], colors="red")
        plt.vlines([-180, 180], ymin=0, ymax=ax.get_ylim()[1], colors="blue")
        # ax.legend()
        fig.savefig(f"plots/contour_plots/kde_plot_{x_var}_{topo}.png")

def correlation_plot(data):
    fig, ax = plt.subplots()
    straight_line=np.linspace(0,700,10)
    band=100
    straight_lower=np.linspace(0-band,700-band,10)
    straight_upper=np.linspace(0+band,700+band,10)
    ax.plot(straight_line, straight_line, "r-")
    ax.plot(straight_line, straight_upper, "b-")
    ax.plot(straight_line, straight_lower, "b-")
    sns.kdeplot(data=data, x="Met", y="pt_higgs", levels=[0.1], hue="sample", weights=data["weight"], alpha=0.5, ax=ax, palette="tab20", common_norm=False)
    fig.savefig(f"plots/contour_plots/correlation_plot.png")



def main():
    total_data=pd.DataFrame()
    for topology in ["resolved", "merged"]:
        sample_set=get_samples(os.path.join(sample_root_dir, topology), topology=topology)
        for samp in sample_set:
            samp.read_sample_list()
            total_data=total_data.append(samp.data)
            del samp

    # plot_density(total_data)
    # make_scatter_plots(total_data)
    kde_plot(total_data, "metpt_ratio")
    kde_plot(total_data, "metpt_diff")
    # correlation_plot((total_data))
    # plt.show()


if __name__ == "__main__":
    main()
