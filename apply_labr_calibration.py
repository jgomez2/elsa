#!/usr/bin/env python
import ROOT as r
import os.path
import sys
import elsa
import array
r.gROOT.SetBatch()

what_gain = 45000.0

ofile = r.TFile('./stage1/labr0_cal.root','recreate')
tree = r.TTree("friend","friend thing")
nph_cal = array.array('f',[0.])
pph_cal = array.array('f',[0.])
#tree.Branch("nph_cal",nph_cal,'nph_cal/F')
tree.Branch("pph_cal",pph_cal,'pph_cal/F')

ifname = './macros/startana.C'
r.gROOT.ProcessLine('.x %s'%(ifname))
r.gROOT.ProcessLine('.x ./hist_pph0.C')
#r.gROOT.ProcessLine('.x ./hist_nph.C')

for i,event in enumerate(r.ch):
  if i%1e5==0:
    print i
  #nph_cal[0] = event.nph*what_gain/r.h_nph_1d.Interpolate(i)
  pph_cal[0] = event.pph*what_gain/r.h_pph0_1d.Interpolate(i)
  tree.Fill()

ofile.Write()
ofile.Close()
  

