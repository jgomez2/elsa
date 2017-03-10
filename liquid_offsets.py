#!/usr/bin/env python
import ROOT as r
import os.path
r.gROOT.SetBatch()
r.gROOT.Macro("./macros/startana.C")

def calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detnum):
    if not os.path.isfile(ifname):
        print 'failed to find file',ifname
        return None
    ifile = r.TFile(ifname)
    hist_2d = ifile.Get("basic/c2tof_b0gated")
    lolim = hist_2d.GetYaxis().FindBin(detnum)
    hilim = hist_2d.GetYaxis().FindBin(detnum)
    hist_2d.GetYaxis().SetRange(lolim,hilim)
    hist = hist_2d.ProjectionX()
    hist.SetDirectory(0)
    ifile.Close()
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
    ifile.Close()
    return mean,res


if __name__=='__main__':
    ifname = './summed_histos/pu239_liquid.root'
    liglrange = range(9,16+1) + range(25,32+1) + range(41,78+1)
    threshlo = 0.4
    threshhi = 0.4
    ofname = './time_calibrations/liquid_tshifts_pu239_2015.txt'
    ofile = open(ofname,'w')
    for detnum in liglrange:
        stuff = calculate_timing(ifname,threshlo,threshhi,detnum)
        print stuff
        offset = 3.33 - stuff[0]
        ostring = '%i %f\n'%(detnum,offset)
        ofile.write(ostring)
    ofile.close()
