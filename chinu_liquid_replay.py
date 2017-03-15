#!/usr/bin/env python
import elsa

if __name__=='__main__':

    """
    rnums = range(30454,30485+1) # 235U 1st run period 2015
    print len(rnums)
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(22,25+1) + range(37,38+1)
        liquidrange = range(9,16+1) + range(26,32+1) + range(41,78+1)
        for i in range(80):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.polarity.append(-1.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.85)
        testclass.detclass[1] = 0
        testclass.t0ds = 15
        testclass.calibrate_time = 1
        for i in ppacrange:
            testclass.detclass[i] = 1 # call ppacs class 1
        for i in liquidrange:
            testclass.detclass[i] = 2 # call liquids class 2

        #thrfname = './gates/thresholds.txt'
        #testclass.read_detinfo_fromfile(thrfname,testclass.threshold)
        #satfname = './gates/saturate.txt'
        #testclass.read_detinfo_fromfile(satfname,testclass.saturate)

        basepath = '/mnt/chi-nu-data/2/chinu'
        ifname = '%s/run%05d.mid.gz'%(basepath,rnum)
        ofname = './stage0/run%05d_s0.bin'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1e9

        testclass.dtype = 0
        testclass.interp_slope = -1.0
        testclass.extra_option = 0
        testclass.ftoption = 1
        testclass.MidasEventPrint = 0
        testclass.MidasEventPrintThresh = 600000;
        testclass.upconvert = 0
        args.append(testclass)
    elsa.run_multithread_elsa(args,10)

    """
    rnums = range(30454,30479+1) # 235U 1st run period 2015
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(21,24+1) + range(37,38+1)
        liquidrange = range(9,16+1) + range(25,32+1) + range(41,47+1) + range(41,78+1)
        for i in range(80):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.85)
        tcalfname = './time_calibrations/run%05d_bcal.txt'%rnum
        #testclass.read_detinfo_fromfile(tcalfname,testclass.toffset)
        ppacshiftname = './time_calibrations/ppac_tshifts_liquid_pu239_2015.txt'
        #testclass.read_tshift_fromfile(ppacshiftname)
        testclass.walkname = './time_calibrations/ppac_walk.txt'
        liquidshiftname = './time_calibrations/liquid_tshifts_pu239_2015.txt'
        #testclass.read_tshift_fromfile(liquidshiftname)
        thrfname = './gates/liquid_thresholds.txt'
        #testclass.read_detinfo_fromfile(thrfname,testclass.threshold)
        satfname = './gates/liquid_saturation.txt'
        #testclass.read_detinfo_fromfile(satfname,testclass.saturate)
        psdfname = './gates/liquid_psd1d.txt'
        #testclass.read_detinfo_fromfile(psdfname,testclass.psd1d)
        testclass.detclass[1] = 0
        testclass.toffset[1] = -800.0
        testclass.t0ds = 15
        testclass.calibrate_time = 1
        for i in ppacrange:
            testclass.detclass[i] = 1 # call ppacs class 1
        for i in liquidrange:
            testclass.detclass[i] = 2 # call ligls class 2

        ifname = './stage0/run%05d_s0.bin'%(rnum)
        ofname = './stage1/run%05d_s1.root'%(rnum)
        ohfname = './stage1h/run%05d.root'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.ohfname = ohfname
        testclass.STOPatEVENT = 1e9

        testclass.t0ebuild = 1
        testclass.ebuild = 1
        testclass.extra_option = 0
        testclass.dtype = 4
        testclass.dump1r = 1
        testclass.dump1h = 0
        args.append(testclass)
    elsa.run_multithread_elsa(args,13)
