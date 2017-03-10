#!/usr/bin/env python
import ROOT as r
import os.path
import sys
import elsa
r.gROOT.SetBatch()

def calculate_energy(ifname,fit_threshold_lo,fit_threshold_hi,detector,hname):
    if not os.path.isfile(ifname):
        print 'failed to find file',ifname
        return None
    r.gROOT.ProcessLine('.x %s'%(ifname))
    r.gROOT.ProcessLine('.x ./photopeak.C')
    #r.gROOT.Execute(ifname)
    # histogram limits (x)
    nbins = 250
    nbinsy = 1000
    xlo = 0.0
    xhi = r.ch.GetEntries() + 1.0
    # histogram limits (y)
    ylo = 0.0
    yhi = 60000.0
    drawcmd = r'%s:Entry$>>%s(%i,%f,%f,%i,%f,%f)'%(detector,hname,nbins,xlo,xhi,nbinsy,ylo,yhi)
    print drawcmd
    r.ch.Draw(drawcmd,"detector_id1==0 && photopeak")
    htemp2d = r.gDirectory.Get(hname)
    hname_1d = hname + '_1d'
    htemp1d = r.TH1F(hname_1d,hname_1d,nbins,xlo,xhi)
    for bini in range(1,htemp2d.GetNbinsX()+1):
      htemp2d.GetXaxis().SetRange(bini,bini)
      hist = htemp2d.ProjectionY()
      maxbin = hist.GetMaximumBin()
      maxcontent = hist.GetMaximum()
      howmanybins = hist.GetNbinsX()
      fit_uplim = 0
      fit_lolim = 0
      for i in range(maxbin,howmanybins):
        if hist.GetBinContent(i) < fit_threshold_hi*maxcontent:
          fit_uplim = i
          break
      for i in range(maxbin,1,-1):
        if hist.GetBinContent(i) < fit_threshold_lo*maxcontent:
          fit_lolim = i
          break
      lolim = hist.GetXaxis().GetBinCenter(fit_lolim)
      uplim = hist.GetXaxis().GetBinCenter(fit_uplim)
      hist.Fit('gaus',"M0Q","",lolim,uplim)
      sig = 0.0
      mean = 0.0
      try:
        sig = hist.GetFunction('gaus').GetParameter(2)
        mean = hist.GetFunction('gaus').GetParameter(1)
      except TypeError:
        sig = 0.0
        mean = 0.0
      res = sig*2.35
      htemp1d.SetBinContent(bini,mean)
      print bini,mean
    return htemp1d
  


if __name__=='__main__':
  ifname = './macros/startana.C'
  fit_threshold_lo = 0.4
  fit_threshold_hi = 0.4
  detector = 'pph'
  hname = 'h_pph0'
  hist_pph = calculate_energy(ifname,fit_threshold_lo,fit_threshold_hi,detector,hname)
  print hist_pph
  hname = 'h_nph'
  detector = 'nph'
  #hist_nph = calculate_energy(ifname,fit_threshold_lo,fit_threshold_hi,detector,hname)
  #print hist_nph
  hist_pph.SaveAs('hist_pph0.C')
  #hist_nph.SaveAs('hist_nph.C')

