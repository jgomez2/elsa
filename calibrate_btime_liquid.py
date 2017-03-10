#!/usr/bin/env python
import ROOT as r
import os.path
import sys
r.gROOT.SetBatch()

def calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi):
    if not os.path.isfile(ifname):
        print 'failed to find file',ifname
        return None
    ifile = r.TFile(ifname)
    hist_2d = ifile.Get("basic/c2tof_b0gated")
    lolim = hist_2d.GetYaxis().FindBin(detlo)
    hilim = hist_2d.GetYaxis().FindBin(dethi)
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
    #print 'found mean of',mean
    return mean,res

if __name__=='__main__':
    if len(sys.argv)<3:
        print 'usage: %s lo hi'%(sys.argv[0])
        exit(0)
    for run in range(int(sys.argv[1]),int(sys.argv[2])+1):
        print 'doing run number',run
        ifname = './stage1h_liquid/run%05i.root'%(run)
        fit_threshold_lo = 0.45
        fit_threshold_hi = 0.45
        detlo = 9
        dethi = 16
        board0 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi)
        detlo = 28
        dethi = 32
        board1 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi)
        detlo = 41
        dethi = 47
        board2 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi)
        detlo = 55
        dethi = 62
        board3 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi)
        detlo = 70
        dethi = 78
        board4 = calculate_timing(ifname,fit_threshold_lo,fit_threshold_hi,detlo,dethi)

        if not board0 or not board1 or not board2 or not board3 or not board4:
            print 'at least one fit failed - quitting'
            exit(0)

        b10 = board1[0]-board0[0]
        b20 = board2[0]-board0[0]
        b30 = board3[0]-board0[0]
        b40 = board4[0]-board0[0]
        b10_baseoffset = 27.0
        b20_baseoffset = 58.5
        b30_baseoffset = 86.2
        b40_baseoffset = 116.1

        b10_int = round((b10+b10_baseoffset)/8.0,0)
        b20_int = round((b20+b20_baseoffset)/8.0,0)
        b30_int = round((b30+b30_baseoffset)/8.0,0)
        b40_int = round((b40+b40_baseoffset)/8.0,0)

        print 'b10',b10,b10_int,'b20',b20,b20_int,'b30',b30,b30_int,'b40',b40,b40_int

        ofname = './time_calibrations/run%05d_bcal.txt'%(run)
        ofile = open(ofname,'w')
        for i in range(1,16+1):
            ostring = '%d %f\n'%(i,0.0)
            ofile.write(ostring)
        for i in range(17,32+1):
            ostring = '%d %f\n'%(i,8.0*b10_int)
            ofile.write(ostring)
        for i in range(33,48+1):
            ostring = '%d %f\n'%(i,8.0*b20_int)
            ofile.write(ostring)
        for i in range(49,64+1):
            ostring = '%d %f\n'%(i,8.0*b30_int)
            ofile.write(ostring)
        for i in range(65,78+1):
            ostring = '%d %f\n'%(i,8.0*b40_int)
            ofile.write(ostring)
        ofile.close()



