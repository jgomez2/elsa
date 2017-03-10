#!/usr/bin/env python
import elsa
import sys

if __name__=='__main__':
   

    runstart = int(sys.argv[1])
    runend = int(sys.argv[2])
    rnums = range(runstart,runend+1)

#STAGE0 STUFF
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        liquidrange = range(9,16+1) + range(26,32+1) + range(41,78+1)
        n_tfilter_params = 4
        for i in range(80):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.polarity.append(-1.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.0)
        for i in liquidrange:
            testclass.detclass[i]=1
        testclass.detclass[1] = 0     
        testclass.toffset[1] = -1000.0
        testclass.t0ds = 15
        temp_tfiltersettings = [2,0.25,2,150] # delay, atten, smooth, thresh
    

        # stage 0: reduce the data
        basepath = './data/'
        ifname = '%s/run%05d.mid.gz'%(basepath,rnum)
        ofname = './stage0/run%05d_s0.bin'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1000000000

        testclass.dtype = 0
        testclass.interp_slope = -1.0
        testclass.extra_option = 0
        testclass.ftoption = 2
        testclass.MidasEventPrint = 0
        testclass.MidasEventPrintThresh = 600000;
        testclass.upconvert = 0
        args.append(testclass)
    #elsa.run_multithread_elsa(args,2)
    elsa.run_elsa(testclass)

#STAGE 1 Stuff
    rnums = range(runstart,runend+1)
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        liquidrange = range(9,16+1) + range(26,32+1) + range(41,78+1)
        for i in range(80):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.0)
        #tcalfname = './time_calibrations/run%05d_bcal.txt'%rnum
        tcalfname = './time_calibrations/licf_toffsets.txt'
        #testclass.read_detinfo_fromfile(tcalfname,testclass.toffset)
        ppacshiftname = './time_calibrations/ppac_tshifts_ligl_pu239_2015.txt'
        #testclass.read_tshift_fromfile(ppacshiftname)
        testclass.walkname = './time_calibrations/ppac_walk.txt'
        liglshiftname = './time_calibrations/ligl_tshifts_pu239_2015.txt'
        #testclass.read_tshift_fromfile(liglshiftname)
        thrfname = './gates/thresholds.txt'
        #testclass.read_detinfo_fromfile(thrfname,testclass.threshold)
        satfname = './gates/saturate.txt'
        #testclass.read_detinfo_fromfile(satfname,testclass.saturate)
        testclass.detclass[1] = 0
        testclass.toffset[1] = -800.0
        testclass.t0ds = 15
        testclass.calibrate_time = 1
        for i in liquidrange:
            testclass.detclass[i] = 1 # call ligls class 2

        ifname = './stage0/run%05d_s0.bin'%(rnum)
        ofname = './stage1/run%05d_s1.root'%(rnum)
        ohfname = './stage1h/run%05d.root'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.ohfname = ohfname
        testclass.STOPatEVENT = 1000000000

        testclass.t0ebuild = 1
        testclass.ebuild = 4
        testclass.extra_option = 0
        testclass.dtype = 4
        testclass.dump1r = 1
        testclass.dump1h = 0
        args.append(testclass)
    elsa.run_multithread_elsa(args,3)

