#!/usr/bin/env python
import ROOT as r
r.gROOT.SetBatch()
r.gROOT.Macro("./macros/startana.C")

def calculate_timing(fit_threshold_lo,fit_threshold_hi,detnum):
  hname = 'hist%d'%(detnum)
  howmanybins = 1000
  drawcmd = r'ntcorr>>%s(%i,0,200)'%(hname,howmanybins)
  gatecmd = r'pph>1500 && nph>1500 && detector_id2==%d'%(detnum)
  #print gatecmd
  r.ch.Draw(drawcmd,gatecmd)
  hist = r.gDirectory.Get(hname)
  hist.SetDirectory(0)
  maxbin = hist.GetMaximumBin()
  max = hist.GetMaximum()
  fit_uplim = 0
  fit_lolim = 0
  for i in range(maxbin,howmanybins):
    if hist.GetBinContent(i) < fit_threshold_hi*max:
      fit_uplim = i
      break
  for i in range(maxbin,1,-1):
    if hist.GetBinContent(i) < fit_threshold_lo*max:
      fit_lolim = i
      break
  lolim = hist.GetXaxis().GetBinCenter(fit_lolim)
  uplim = hist.GetXaxis().GetBinCenter(fit_uplim)
  hist.Fit('gaus',"M0Q","",lolim,uplim)
  sig = hist.GetFunction('gaus').GetParameter(2)
  mean = hist.GetFunction('gaus').GetParameter(1)
  res = sig*2.35
  return mean,res,hist

if __name__=='__main__':
    detrange = range(9,16+1) + range(25,32+1) + range(41,78+1)
    threshlo = 0.3
    threshhi = 0.3
    ofname = './time_calibrations/ligl_tshifts.txt'
    ofile = open(ofname,'w')
    for detnum in detrange:
        stuff = calculate_timing(threshlo,threshhi,detnum)
        print stuff
        offset = 1.33 - stuff[0]
        ostring = '%i %f\n'%(detnum,offset)
        ofile.write(ostring)
    ofile.close()

