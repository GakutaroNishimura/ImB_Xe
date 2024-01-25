import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ROOT

# df = pd.read_table("./He_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")
df = pd.read_table("./Xe131_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")

gr = ROOT.TGraph(len(df.energy), np.array(df.energy), np.array(df.CS))
# gr_fit = ROOT.TF1("f", "[0]/sqrt(x) + [1]", 2.0, 20.0)
# gr_fit.SetParameters(900.0, 1.0)
# gr_fit.FixParameter(1, 0.0)
# gr.Fit(gr_fit, "QR")
# par = [gr_fit.GetParameter(k) for k in range(gr_fit.GetNpar())]
# print(par)
gr.SetMarkerStyle(7)
gr.SetMarkerSize(10)
gr.Draw("APL")
gr.GetXaxis().SetRangeUser(2.0, 20.0)
# gr.GetYaxis().SetRangeUser(1e2, 1e3)
c1 = ROOT.gROOT.FindObject("c1")
c1.SetLogx(1)
c1.SetLogy(1)
c1.SetLogy()