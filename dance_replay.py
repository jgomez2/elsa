#!/usr/bin/env python
import elsa
import sys

if __name__=='__main__':
    if len(sys.argv)!=4:
        print 'usage: dance_replay path lo hi'
        exit(0)

    args = []
    for rnum in range(int(sys.argv[2]),(int(sys.argv[3])+1)):
        print 'here'
        # single threaded basic testing
        testclass = elsa.Elsa()
        testclass.MidasEventPrint = 0
        testclass.MidasEventPrintthresh = 0
        testclass.dump0b = 1

        # detector class setup
        for i in range(400):
            testclass.detclass.append(1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
        testclass.detclass[176]=0
        testclass.lsgate_ds.append(-1)
        testclass.t0ds = 0
        testclass.extra_option = 0
        testclass.calibrate_time = 0

        # stage 0: reduce the data
        ifname = '%s/run%06i.mid.gz'%(sys.argv[1],rnum)
        ofname = './stage0-dance/run%06i_s0.bin'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1e9

        testclass.dtype = 1
        #testclass.lsgate = './dance_gates/baf2_alphas.C'
        elsa.run_elsa(testclass)
        #args.append(testclass)


        testclass.ifname = './stage0-dance/run%06i_s0.bin'%(rnum)
        testclass.ofname = './stage1-dance/run%06i.root'%(rnum)
        testclass.ebuild = 2
        testclass.dtype = 4
        testclass.t0ebuild = 1
        elsa.run_elsa(testclass)
        #args.append(testclass)
    #elsa.run_multithread_elsa(args,4)

    """
    # multithreaded "production-like" testing for stage 0
    #rnums = range(19150,19170+1)
    #rnums = range(19171,19195+1)
    #rnums = range(19196,19222+1)
    rnums = range(19150,19222+1)
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(20,23+1) + range(37,40+1)
        liglrange = range(9,16+1) + range(24,32+1) + range(41,48+1)
        for i in range(50):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
        testclass.detclass[1] = 0
        testclass.t0ds = 15
        for i in ppacrange:
            testclass.detclass[i] = 1 # call ppacs class 1
        for i in liglrange:
            testclass.detclass[i] = 2 # call ligls class 2

        ifname = '/mnt/hygd-data/7/tke/chinu-mcnp/2015_data/run%05d.mid.gz'%(rnum)
        ofname = './stage0/run%05d_s0.bin'%(rnum)
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1e9

        testclass.dtype = 0
        testclass.interp_slope = -1.0
        args.append(testclass)
    elsa.run_multithread_elsa(args,2)
    """

    """
    # multithreaded "production-like" testing for stage 1
    rnums = range(19150,19222+1)
    #rnums = range(19194,19222+1)
    args = []
    for rnum in rnums:
        testclass = elsa.Elsa()
        # detector class setup
        ppacrange = range(5,8+1) + range(20,23+1) + range(37,40+1)
        liglrange = range(9,16+1) + range(24,32+1) + range(41,48+1)
        for i in range(50):
            testclass.detclass.append(-1)
            testclass.tslope.append(2.0)
            testclass.toffset.append(0.0)
            testclass.eslope.append(1.0)
            testclass.eoffset.append(0.0)
        tcalfname = './time_calibrations/run%05d_bcal.txt'%rnum
        testclass.read_tcal_fromfile(tcalfname)
        ppacshiftname = './time_calibrations/ppac_tshifts.txt'
        testclass.read_tshift_fromfile(ppacshiftname)
        liglshiftname = './time_calibrations/ligl_tshifts.txt'
        testclass.read_tshift_fromfile(liglshiftname)
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
        testclass.ifname = ifname
        testclass.ofname = ofname
        testclass.STOPatEVENT = 1e9

        testclass.t0ebuild = 1
        testclass.ebuild = 1
        testclass.dtype = 4
        args.append(testclass)
    elsa.run_multithread_elsa(args,2)
    """
