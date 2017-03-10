#!/usr/bin/env python
import ROOT as r
import os.path
r.gROOT.SetBatch()

def calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,rnum,detlo,dethi):
  ifile = r.TFile(ifname)
  try:
    ifile.tree.SetAlias("tof","((100000.*(t0tof-0.0)%(100000*1788.819875776))/100000.)")
  except AttributeError:
    return 0.0,0.0,None
  hname = 'hist%05d_%d_%d'%(rnum,detlo,dethi)
  howmanybins = 1000
  drawcmd = r'tof>>%s(%i,0,1000)'%(hname,howmanybins)
  gatecmd = r'pph>2000 && detector_id1>=%d && detector_id1<=%d'%(detlo,dethi)
  #print gatecmd
  ifile.tree.Draw(drawcmd,gatecmd)
  hist = r.gDirectory.Get(hname)
  hist.SetDirectory(0)
  ifile.Close()
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
  sig = 0.0
  mean = 0.0
  try:
    sig = hist.GetFunction('gaus').GetParameter(2)
    mean = hist.GetFunction('gaus').GetParameter(1)
  except TypeError:
    sig = 0.0
    mean = 0.0
  res = sig*2.35
  ifile.Close()
  return mean,res,hist

if __name__=='__main__':
    """
    rnum = 19151
    fit_threshold_lo = 0.1
    fit_threshold_hi = 0.7
    detlo = 5
    dethi = 8
    ifname = './stage1/run%d_s1.root'%(rnum)
    stuff = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,rnum,detlo,dethi)
    print stuff
    """

    fit_threshold_lo = 0.1
    fit_threshold_hi = 0.7
    #runs = range(19150,19154+1) + range(19156,19198+1) + range(19201,19222+1)
    #runs = range(19194,19198+1) + range(19201,19222+1)
    #runs = range(19582,19597+1) + range(19600,19608+1) + range(19610,19638+1) #range(19582,19638+1) # second run period 2015
    #runs = range(19150,19193+1)
    #rnums = range(18373,18410+1) + range(18418,18447+1) + range(18455,18481+1) + range(18482,18562+1)# jan 2015 pu239 run
    runs = range(18460,18562+1)# jan 2015 pu239 run
    for run in runs:
        ifname = './stage1-chinu/run%d_s1.root'%(run)
	if not os.path.isfile(ifname):
          continue
        print 'doing file',ifname
        detlo = 5
        dethi = 8
        board0 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,run,detlo,dethi)
        detlo = 21
        dethi = 23
        board1 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,run,detlo,dethi)
        detlo = 37
        dethi = 38
        board2 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,run,detlo,dethi)
        b01 = board1[0]-board0[0]
        b02 = board2[0]-board0[0]
        print b01,b02
        b01_int = round((b01+21.0)/8.0,0)
        b02_int = round((b02+48.0)/8.0,0)
        print b01_int,b02_int
        ofname = './time_calibrations/run%05d_bcal.txt'%(run)
        ofile = open(ofname,'w')
        for i in range(1,16+1):
            ostring = '%d %f\n'%(i,0.0)
            ofile.write(ostring)
        for i in range(17,32+1):
            ostring = '%d %f\n'%(i,-8.0*b01_int)
            ofile.write(ostring)
        for i in range(33,48+1):
            ostring = '%d %f\n'%(i,-8.0*b02_int)
            ofile.write(ostring)
        ofile.close()


