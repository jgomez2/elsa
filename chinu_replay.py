#!/usr/bin/env python
import elsa
import sys

if __name__=='__main__':
    """
    # single threaded basic testing
    testclass = elsa.Elsa()

    # detector class setup
    ppacrange = range(5,8+1) + range(20,23+1) + range(37,40+1)
    liglrange = range(9,16+1) + range(24,32+1) + range(41,48+1)
    n_tfilter_params = 4
    for i in range(50):
        testclass.detclass.append(-1)
        testclass.tslope.append(2.0)
        testclass.toffset.append(0.0)
        testclass.eslope.append(1.0)
        testclass.eoffset.append(0.0)
        testclass.polarity.append(-1.0)
        for j in range(n_tfilter_params):
            testclass.time_filter_params.append(0.0)
    testclass.detclass[1] = 0
    testclass.toffset[1] = -1000.0
    testclass.t0ds = 15
    temp_tfiltersettings = [2,0.25,2,150] # delay, atten, smooth, thresh
    for j,stuff in enumerate(temp_tfiltersettings):
        testclass.time_filter_params[n_tfilter_params*1+j] = stuff
    testclass.calibrate_time = 1
    for i in ppacrange:
        testclass.detclass[i] = 1 # call ppacs class 1
        temp_tfiltersettings = [9,0.75,4,150] # delay, atten, smooth, thresh
        for j,stuff in enumerate(temp_tfiltersettings):
            testclass.time_filter_params[n_tfilter_params*i+j] = stuff
    for i in liglrange:
        testclass.detclass[i] = 2 # call ligls class 2
        temp_tfiltersettings = [5,0.5,2,100] # delay, atten, smooth, thresh
        for j,stuff in enumerate(temp_tfiltersettings):
            testclass.time_filter_params[n_tfilter_params*i+j] = stuff

    # stage 0: reduce the data
    ifname = './test_files/run18375.mid.gz'
    ofname = './test_files/run18375_s0.bin'
    testclass.ifname = ifname
    testclass.ofname = ofname
    testclass.STOPatEVENT = 1e9

    testclass.dtype = 0
    testclass.extra_option = 1
    testclass.ftoption = 2
    testclass.MidasEventPrint = 0
    testclass.upconvert = 1
    #elsa.run_elsa(testclass)


    testclass.ifname = './test_files/run18375_s0.bin' 
    testclass.ofname = './test_files/run18375_s1.root' 
    testclass.t0ebuild = 1
    testclass.ebuild = 1
    testclass.dtype = 4
    #elsa.run_elsa(testclass)
    """

    # multithreaded "production-like" testing for stage 0
    #rnums = range(19150,19222+1) # 235U 1st run period 2015
    #rnums = range(19582,19638+1) # 235U 2nd run period 2015
    #rnums = range(18373,18410+1) + range(18418,18447+1) + range(18455,18481+1) + range(18482,18562+1)# jan 2015 pu239 run
    #rnums = range(20557,20698+1) + range(20482,20555+1) + range(20369,20481+1) # dec 2015 pu239 run
    #rnums = range(19126,20223+1) # complete low energy U5 2015
    """
    rnums = elsa.read_runlist('./runs/2015_ligl_u235.txt')
"""
    runstart = int(sys.argv[1])
    runend = int(sys.argv[2])
    rnums = range(runstart,runend+1)

    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(22,25+1) + range(37,38+1)
        liglrange = range(9,16+1) + range(26,32+1) + range(41,47+1)
        for i in range(50):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.polarity.append(-1.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.0)
        testclass.detclass[1] = 0
        testclass.t0ds = 15
        testclass.calibrate_time = 1
        for i in ppacrange:
            testclass.detclass[i] = 1 # call ppacs class 1
        for i in liglrange:
            testclass.detclass[i] = 2 # call ligls class 2

        #thrfname = './gates/thresholds.txt'
        #testclass.read_detinfo_fromfile(thrfname,testclass.threshold)
        #satfname = './gates/saturate.txt'
        #testclass.read_detinfo_fromfile(satfname,testclass.saturate)

        basepath = './data/'
#  if rnum > 18000 and rnum < 19100:
#      basepath = '/data/1/chinu'
#  elif rnum >19099 and rnum < 20200:
#      basepath = '/data/2/chinu'
#  elif rnum >20199 and rnum < 21842:
#      basepath = '/data/3/chinu'
#  elif rnum >21841 and rnum < 22292:
#      basepath = '/data/4/chinu'
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
    elsa.run_multithread_elsa(args,9)
   

    # multithreaded "production-like" testing for stage 1
    #rnums = range(19150,19222+1) # 1st run period 
    #rnums = range(19194,19222+1)
    #rnums = range(19582,19638+1) # 235U 2nd run period 2015
    #rnums = range(18373,18410+1) + range(18418,18447+1) + range(18455,18481+1) + range(18482,18562+1)# jan 2015 pu239 run
    #rnums = range(20557,20668+1) # dec 2015 pu239 run
    #rnums = range(20557,20564+1) # dec 2015 pu239 run
    #rnums = elsa.read_runlist('./runs/2015_ligl_u235.txt')
    runstart = int(sys.argv[1])
    runend = int(sys.argv[2])
    rnums = range(runstart,runend+1)
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(21,24+1) + range(37,38+1)
        liglrange = range(9,16+1) + range(25,32+1) + range(41,47+1)
        for i in range(50):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.0)
        tcalfname = './time_calibrations/run%05d_bcal.txt'%rnum
        testclass.read_detinfo_fromfile(tcalfname,testclass.toffset)
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
        testclass.toffset[1] = -900.0
        testclass.t0ds = 15
        testclass.calibrate_time = 1
        for i in ppacrange:
            testclass.detclass[i] = 1 # call ppacs class 1
        for i in liglrange:
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
    elsa.run_multithread_elsa(args,9)
