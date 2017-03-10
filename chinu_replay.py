#!/usr/bin/env python
import elsa
import sys

if __name__=='__main__':
   
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
    ifname = './data/run29559.mid.gz'
    ofname = './stage0/run29559_s0.bin'
    testclass.ifname = ifname
    testclass.ofname = ofname
    testclass.STOPatEVENT = 1e9

    testclass.dtype = 0
    testclass.extra_option = 1
    testclass.ftoption = 2
    testclass.MidasEventPrint = 0
    testclass.upconvert = 1
    elsa.run_elsa(testclass)


    testclass.ifname = './stage0/run29559_s0.bin' 
    testclass.ofname = './stage1/run29559_s1.root' 
    testclass.t0ebuild = 1
    testclass.ebuild = 1
    testclass.dtype = 4
    #elsa.run_elsa(testclass)
