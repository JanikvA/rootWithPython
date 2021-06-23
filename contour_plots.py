#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contour plots from TTrees
"""

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
N_max=400

# topology="resolved"
# topology="merged"

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
            pt_higgs="pt_jj_corr"
        elif self.topology=="merged":
            pt_higgs="pt_J_corr"
        df = tree.arrays( ["Met",pt_higgs], library="pd")
        if self.topology=="resolved":
            df = df[df["Met"]<500]
            # df = df[df["Met"]>500]
        elif self.topology=="merged":
            # df = df[df["Met"]<500]
            df = df[df["Met"]>500]
        df=df.rename(columns={pt_higgs:"pt_higgs"})
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
            self.data["metpt"]=self.data.apply(lambda row: row["Met"]/row["pt_higgs"], axis=1)
        print(f"Loaded {len(self.data)} events from {self.name}!")
        # print(self.data.corr())



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
    return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh50)
    # return (VHbb, Diboson, stop, ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, Wjets, signal_mzp3500_dm200_dh150)
    # return (ttbar, Zjets, signal_mzp500_dm200_dh150)
    # return (Zjets, signal_mzp3500_dm200_dh150)
    # return (signal_mzp3500_dm200_dh150, Zjets)
    # return (Zjets,)
    # return (signal_mzp3500_dm200_dh150,)
    # return (signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50, signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    # return (signal_mzp3500_dm200_dh150, signal_mzp3500_dm200_dh50)
    # return (signal_mzp500_dm200_dh150, signal_mzp500_dm200_dh50)
    # return (signal_mzp500_dm200_dh50,)


def plot_samples(data):
    # g = sns.pairplot(data, hue="sample", diag_kind="hist")
    # g.map_lower(sns.kdeplot, levels=[0.1,0.5])

    # sns.pairplot(data, kind="kde", hue="sample", diag_kind="hist", levels=[0.1,0.5])
    # sns.pairplot(data, kind="kde", hue="sample_type", diag_kind="hist")
    # sns.pairplot(data, hue="sample", diag_kind="hist")
    # sns.kdeplot(data=data, x="Met", y=pt_higgs, hue="sample_type", levels=[0.1,0.5,0.9])
    # sns.kdeplot(data=data, x="Met", y=pt_higgs, hue="sample", levels=[0.1])
    # sns.scatterplot(data=data, x="Met", y=pt_higgs, hue="sample")

    # f = plt.figure(figsize=(19, 15))
    # plt.matshow(data.corr(), fignum=f.number)

    # ax = sns.regplot(x="Met", y="met/pt", data=data)
    # data.plot.hist(bins=40)

    # test=pd.DataFrame()
    # for samp_name in data["sample"].unique():
    #     tmp=data[data["sample"]==samp_name]["metpt"].copy(deep=True)
    #     test[samp_name]=tmp.values
    #     test[samp_name]=tmp.values
    # # test.plot.density()
    # test.plot(kind="kde")
    # sns.kdeplot(data=[data[data["sample"]==samp for samp in data["sample"].unique()], x="metpt", hue="sample")
    is_signal=(data["sample"] == "mzp3500_dh50_merged") | (data["sample"] == "mzp3500_dh50_resolved")
    sns.kdeplot(data=data[is_signal], x="metpt", hue="sample", fill=True, common_norm=False)
    sns.kdeplot(data=data[~is_signal], x="metpt", hue="sample", palette="tab20", common_norm=False)
    plt.show()

def make_scatter_plots(data):
    sample_list=data["sample"].unique()
    straight_line=np.linspace(0,1500,10)
    for n, samp_name in enumerate(sample_list):
        if "merged" in samp_name: continue
        fig=plt.figure()
        sns.scatterplot(data=data[data["sample"]==samp_name], x="Met", y="pt_higgs", hue="sample")
        sns.scatterplot(data=data[data["sample"]==samp_name.replace("resolved","merged")], x="Met", y="pt_higgs", hue="sample", palette="inferno")
        plt.plot(straight_line, straight_line, "r-")
        fig.savefig(f"plots/contour_plots/recycled_{samp_name}.png")

def main():
    total_data=pd.DataFrame()
    for topology in ["resolved", "merged"]:
        sample_set=get_samples(os.path.join(sample_root_dir, topology), topology=topology)
        for samp in sample_set:
            samp.read_sample_list()
            total_data=total_data.append(samp.data)
            del samp

    plot_samples(total_data)
    # make_scatter_plots(total_data)


if __name__ == "__main__":
    main()
