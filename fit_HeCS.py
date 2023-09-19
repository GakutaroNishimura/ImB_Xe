import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ROOT

df = pd.read_table("./He_totalCS.txt", names=["energy", "CS"], sep="    ", skiprows=1, engine="python")

gr = ROOT.TGraph(len(df.energy), np.array(df.energy), np.array(df.CS))
gr_fit = ROOT.TF1("f", "[0]/((x+211000)*(x+211000)+[1]*[1])", 5.0, 10.0)
gr_fit.SetParameters(1., 1.)
gr.Fit(gr_fit, "QR")
gr.SetMarkerStyle(7)
gr.SetMarkerSize(10)
gr.Draw("APL")