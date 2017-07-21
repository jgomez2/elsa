#!/usr/bin/env python
import elsa
import sys

if __name__=='__main__':
 
    runstart = int(sys.argv[1])
    runend = int(sys.argv[2])
    rnums = range(runstart,runend+1)
 
#    rnums = range(26150,26150+1) 
    args = []
    for rnum in rnums:
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
            testclass.threshold.append(0.0)
            testclass.saturate.append(70000.0)
            testclass.psd1d.append(0.0)
	    for j in range(n_tfilter_params):
		testclass.time_filter_params.append(0.0)
        testclass.detclass[1] = 0
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

        ifname = 'data/run%05d.mid.gz'%(rnum)
        ofname = './stage0/run%05d_qils.mid'%(rnum)

        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1e9

        testclass.dtype = 51
        testclass.dump0b = 1
        testclass.interp_slope = -1.0
        testclass.extra_option = 0
        testclass.ftoption = 1
        testclass.MidasEventPrint = 0
        testclass.MidasEventPrintThresh = 500000;
        testclass.upconvert = 0
        args.append(testclass)
    elsa.run_multithread_elsa(args,1)
