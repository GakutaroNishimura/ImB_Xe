import pandas as pd
import matplotlib.pyplot as plt
import ROOT
import numpy as np

df = pd.read_excel("./Pulse_0310_mod1.xlsx", sheet_name=1, usecols=[0, 10], header=2, names=["energy", "counts"])
# print(df.energy[:5])

# Counts = [4.8e7, 1.2e7, 1.2e6]
# Energy = [0.4, 1e6, 1e7]
# plt.plot(df.energy, df.counts, marker=".",  linestyle="None")
# plt.plot(Energy, Counts, marker=".",  linestyle="None")
# plt.yscale("log")
# plt.xscale("log")
# # plt.xlim([1.9, 4.0])
# # plt.ylim([6e10, 1.2e11])
# plt.show()
# print(df[:5])

gr = ROOT.TGraph(len(df.energy), np.array(df.energy), np.array(df.counts, dtype="double"))
# grFit = ROOT.TF1("f", "pol3", 0.0, 0.4)
# gr.Fit(grFit, "QR")
gr.Draw("AP")
gr.SetMarkerStyle(7)
gr.SetMarkerSize(10)
gr.SetTitle(" ")
gr.GetXaxis().SetTitle("Energy [eV]")
gr.GetYaxis().SetTitle("neutron intensity [n/cm^2/s/sr/eV]")
gr.GetXaxis().SetRangeUser(1.0*10**(-4), 2.0*10**(4))
gr.GetYaxis().SetRangeUser(1e8, 5e14)
# ROOT.gStyle.SetOptLogy(1)
c = ROOT.gROOT.FindObject("c1")
c.SetLogy(1)
c.SetLogx(1)
c.Draw("same")