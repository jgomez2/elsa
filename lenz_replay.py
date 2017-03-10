#!/usr/bin/env python
import elsa

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
    #rnums = elsa.read_runlist('./runs/2015_ligl_u235.txt')
    #rnums = range(13540,13553+1)
    #rnums = range(12623,12623+1)
    #rnums = range(12623,12649+1)
    rnums = range(13679,13679+1)
    #rnums = range(13679,13827+1) + range(13890,14154+1);
    rnums = elsa.read_runlist('./runs/lenz_ssd_2015.txt')
    print len(rnums)


    """
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        siliconrange = range(16,48+1)
        n_tfilter_params = 4
        n_efilter_params = 4
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
            for j in range(n_tfilter_params):
                testclass.time_filter_params.append(0.0)
            for j in range(n_efilter_params):
                testclass.energy_filter_params.append(0.0)
        testclass.calibrate_time = 1

        # SET UP THE T0
        testclass.detclass[0] = 0
        temp_tfiltersettings = [2,0.0,1,1000] # delay, atten, smooth, thresh
        for j,stuff in enumerate(temp_tfiltersettings):
            testclass.time_filter_params[n_tfilter_params*0+j] = stuff
        temp_efiltersettings = [1.,20.,-5.0,5.0] # baseline lo, baseline hi, ci lo, ci hi
        for j,stuff in enumerate(temp_efiltersettings):
            testclass.energy_filter_params[n_efilter_params*0+j] = stuff
        testclass.t0ds = 0
        print testclass.time_filter_params

        # SET UP THE FLUX MONITOR
        testclass.detclass[2] = 1 # call the flux monitor class 1
        temp_tfiltersettings = [5,0.0,2,1000] # delay, atten, smooth, thresh
        for j,stuff in enumerate(temp_tfiltersettings):
            testclass.time_filter_params[n_tfilter_params*2+j] = stuff
        temp_efiltersettings = [1.,20.,20.,80.] # baseline lo, baseline hi, ci lo, ci hi
        for j,stuff in enumerate(temp_efiltersettings):
            testclass.energy_filter_params[n_efilter_params*2+j] = stuff

        # SET UP THE SILICONS
        temp_tfiltersettings = [7,0.0,3,200] # delay, atten, smooth, thresh
        for i in siliconrange:
            testclass.detclass[i] = 2 # call silicons class 2
            testclass.polarity[i] = 1.0
            for j,stuff in enumerate(temp_tfiltersettings):
                testclass.time_filter_params[n_tfilter_params*i+j] = stuff
            temp_efiltersettings = [1.,140.,10.,610.] # baseline lo, baseline hi, ci lo, ci hi
            for j,stuff in enumerate(temp_efiltersettings):
                testclass.energy_filter_params[n_efilter_params*i+j] = stuff

        #thrfname = './gates/thresholds.txt'
        #testclass.read_detinfo_fromfile(thrfname,testclass.threshold)
        #satfname = './gates/saturate.txt'
        #testclass.read_detinfo_fromfile(satfname,testclass.saturate)

        basepath = '/data/6/chinu/lenz/data2015/'
        ifname = '%s/run%05d.bin.gz'%(basepath,rnum)
        ofname = './stage0_lenz/run%05d_s0.bin'%(rnum)
        ohfname = './stage0h_lenz/run%05d.root'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.ohfname = ohfname
        testclass.STOPatEVENT = 1e9

        testclass.dtype = 7 # devt stupid code binary w/ 1 signal trace
        testclass.extra_option = 2 # over-write long integral with pulse height
        testclass.ftoption = 3 # calculate fine timne from double derivative
        testclass.MidasEventPrint = 0
        testclass.MidasEventPrintThresh = 600000;
        testclass.upconvert = 0
        testclass.dump0h = 0
        args.append(testclass)
    elsa.run_multithread_elsa(args,11)
    """

    # multithreaded "production-like" testing for stage 1
    #rnums = range(13679,13827+1) + range(13890,14154+1);
    #rnums = range(13679,13688+1)
    rnums = elsa.read_runlist('./runs/lenz_ssd_2015.txt')
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        siliconrange = range(16,48+1)
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
        testclass.detclass[0] = 0
        testclass.toffset[0] = -900.0
        testclass.detclass[1] = 1
        testclass.t0ds = 0
        testclass.calibrate_time = 1
        for i in siliconrange:
            testclass.detclass[i] = 2 # call silicons class 2

        ifname = './stage0_lenz/run%05d_s0.bin'%(rnum)
        ofname = './stage1_lenz/run%05d_s1.root'%(rnum)
        ohfname = './stage1h_lenz/run%05d.root'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.ohfname = ohfname
        testclass.STOPatEVENT = 1e9

        testclass.t0ebuild = 1
        testclass.ebuild = 5 # build lenz style events
        testclass.extra_option = 2
        testclass.dtype = 4
        testclass.dump1r = 0
        testclass.dump1h = 1
        testclass.t0ds = 6
        args.append(testclass)
    elsa.run_multithread_elsa(args,10)
