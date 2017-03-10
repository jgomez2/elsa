import ROOT as r
r.gSystem.Load("elsa_C.so")
from multiprocessing import Pool
import array
import os.path

def read_runlist(fname):
    runlist = []
    for line in open(fname):
        thing = line.split()
        runlist += range(int(thing[0]),int(thing[1])+1)
    return runlist

class Elsa:
    def __init__(self):
        #self.eroot = r.elsa_root()

        # file read options
        self.rnum = 0 # run number to process
        self.ifname = '' # path to input file
        self.wfrm = 1 # remove waveforms from data stream
        self.dtype = 0 # data stream type
                       # 0 = UAC QILS banks
                       # 1 = DANCE UAC banks
                       # 2 = ANNA binary output
                       # 3 = OLAF binary output
                       # 4 = ELSA stage 0 binary output
                       # 5 = UAC PILS banks (for archival data)
                       # 6 = stupid code output (devt) w/ no trace
                       # 7 = stupid code output (devt) w/ one trace
                       # 8 = stupid code output (devt) w/ dual trace
        self.interp_slope = 1.0 # deal with UAC bug where interpolated time goes the wrong way

        self.MidasEventPrint = 0
        self.MidasEventPrintThresh = 0
        self.STOPatEVENT = 0

        # file write options
        self.dump0b = 1 # dump stage 0 binary output
        self.dump0r = 0 # dump stage 0 root output
        self.dump0h = 0 # dump stage 0 histograms
        self.dump1r = 0 # dump stage 1 root output
        self.dump1h = 0 # dump stage 1 histograms
        self.ofname = '' # path to output data file
        self.ohfname = '' # path to output histogram file
        self.extra_option = 0 # extra word option for output file
                          # 0 = off
                          # 1 = write pulse height

        # T0 handling
        self.t0ds = 10 # t0 downscale factor
        self.t0to = 1e6 # t0 timeout in ns - resets the downscaler counter

        # detector description
        self.ndets = 48 # number of detectors in the system
        self.nconvention = 0 # detector numbering convention
                             # 0 = UAC (starts at 1)
                             # 1 = ANNA / OLAF (starts at 0)
        self.detid_shift = 0 # whether / how far to shift detector id
        self.detclass = []
        # timing options
        self.upconvert = 0 # whether to look for clock rollovers
        self.tcalibrate = 0 # use timing offsets
        self.ftoption = 1 # fine time option
                          # 0 = no fine time
                          # 1 = firmware time
                          # 2 = digital CFD - requires waveform processing
                          # 3 = digital double derivative - requires waveform processing
                          # 4 = DANCE timing algorithm
                          # 5 = asymmetrically smoothed CFD (experimental)
        self.ftopionind = [] # over-ride the global ftoption for individual channels

        self.time_filter_params = [] # timing filter params for waveform time extractions
        self.energy_filter_params = [] # energy filter params for waveform energy extractions
        self.polarity = [] # timing filter params for waveform time extractions

        # calibration options
        self.calibrate_energy = 0 # use energy slopes / offsets
        self.calibrate_time = 0 # use energy slopes / offsets
        self.tslope = []
        self.toffset = []
        self.eslope = []
        self.eoffset = []
        self.polarity = []
        self.walkname = ''

        # gating options
        self.threshold = []
        self.saturate = []
        self.psd1d = []

        # singles filtering options
        # things that are IN the gate get excluded here
        self.lsgate_ds = [] # long/short gate downscaling channelwise. Options:
                            # 0 = no gate at all
                            # -1 = get rid of all
                            # 1 -> whatever = downscale counter
        self.lsint_sub = 0 # whether to change long integral to be long - short
        self.lsgate = '' # string path to the gate


        # event build controls
        self.ebuild = 0 # whether to build events at all
                        # 0 = don't build events
                        # 1 = build chi-nu style events
                        # 2 = build dance style events
                        # 3 = build mona style events
                        # 4 = build rpi style events
                        # 5 = build lenz style events
        self.t0ebuild = 0 # whether to consider t0s as part of event build
                          # 0 = no
                          # 1 = yes

        # auxiliary mona stuff
        self.detid2side = []
        self.detid2bar = []


    def convert_args(self):
        # the class we're really hooking into
        self.eroot = r.elsa_root()

        self.eroot.ipath = r.TString(self.ifname)
        self.eroot.opath = r.TString(self.ofname)
        self.eroot.ohpath = r.TString(self.ohfname)
        print "hist path",self.eroot.ohpath
        self.eroot.ftoption = self.ftoption
        temp = array.array('i',self.detclass)
        self.eroot.detclass = temp

        self.eroot.tslope = array.array('d',self.tslope)
        self.eroot.toffset = array.array('d',self.toffset)
        self.eroot.eslope = array.array('d',self.eslope)
        self.eroot.eoffset = array.array('d',self.eoffset)

        # gating options
        self.eroot.threshold = array.array('d',self.threshold)
        self.eroot.saturate = array.array('d',self.saturate)
        self.eroot.psd1d = array.array('d',self.psd1d)

        # t0 downscaling
        self.eroot.t0ds = self.t0ds
        self.eroot.t0to = self.t0to
        
        # event printing
        self.eroot.MidasEventPrint = self.MidasEventPrint
        self.eroot.MidasEventPrintThresh = self.MidasEventPrintThresh
        self.eroot.STOPatEVENT = self.STOPatEVENT

        # file writing
        self.eroot.dump0b = self.dump0b # dump stage 0 binary output
        self.eroot.dump0r = self.dump0r # dump stage 0 root output
        self.eroot.dump0h = self.dump0h # dump stage 0 histogram output
        self.eroot.dump1r = self.dump1r # dump stage 1 root output
        self.eroot.dump1h = self.dump1h # dump stage 1 histogram output
        self.eroot.extra_option = self.extra_option

        # waveform processing
        if len(self.polarity)>0:
            self.eroot.polarity = array.array('f',self.polarity)
        if len(self.time_filter_params)>0:
            self.eroot.time_filter_params = array.array('f',self.time_filter_params)
        if len(self.energy_filter_params)>0:
            self.eroot.energy_filter_params = array.array('f',self.energy_filter_params)

        # timing options
        self.eroot.upconvert = self.upconvert

        # lsgate downscaling options
        if len(self.lsgate_ds)>0:
            self.eroot.lsgate_ds = array.array('i',self.lsgate_ds)
        if len(self.lsgate)>0:
            r.gROOT.Macro(self.lsgate)
            self.eroot.lsgate = r.cutg
        self.eroot.lsint_sub = self.lsint_sub

        if len(self.walkname)>0:
            self.eroot.walk = r.TGraph(self.walkname)

        # auxiliary mona stuffs
        if len(self.detid2side)>0:
            self.eroot.detid2side = array.array('i',self.detid2side)
        if len(self.detid2bar)>0:
            self.eroot.detid2bar = array.array('i',self.detid2bar)

        # detector id options
        self.eroot.detid_shift = self.detid_shift
        print 'out of convert args'

        """
        # individual ftoption choices
        if len(self.ftoptionind)>0:
            self.eroot.ftoptionind = array.array('i',self.ftoptionind)
        else:
	    self.eroot.ftoptionind = array.array("i",)
        """
        return

    def read_detinfo_fromfile(self,fname,temparray):
        if os.path.isfile(fname):
            for line in open(fname):
                temp = line.split()
                temparray[int(temp[0])] = float(temp[1])
        else:
            print "can't find file",fname,"so use tcals of 0"

    def read_tshift_fromfile(self,fname):
        if os.path.isfile(fname):
            for line in open(fname):
                temp = line.split()
                self.toffset[int(temp[0])] += float(temp[1])
        else:
            print "can't find file",fname,"so use existing tcals"



def run_elsa(args):
    args.convert_args()
    print 'looking for file',args.ifname
    if not os.path.isfile(args.ifname):
      return
    if args.dtype == 0:
        print 'running in qils',args.ifname
        args.eroot.execute_uac_qils(args.interp_slope)
    elif args.dtype == 1:
        print 'running in cevt',args.ifname
        args.eroot.execute_uac_cevt(args.interp_slope)
    elif args.dtype == 4:
        print 'running in s0bin'
        args.eroot.execute_s0bin(args.calibrate_time,args.calibrate_energy)
    elif args.dtype == 5:
        print 'running in pils',args.ifname
        args.eroot.execute_uac_pils(args.interp_slope)
    elif args.dtype >5 and args.dtype<9:
        print 'running in devt',args.ifname
        args.eroot.execute_devt(args.dtype)

    if args.ebuild == 1:
        print 'running in chinu'
        args.eroot.build_events_chinu(args.t0ebuild)
    elif args.ebuild == 2:
        args.eroot.build_events_dance(args.t0ebuild)
    elif args.ebuild == 3:
        args.eroot.build_events_mona(args.t0ebuild)
    elif args.ebuild == 4:
        args.eroot.build_events_rpi()
    elif args.ebuild == 5:
        print 'running in lenz'
        args.eroot.build_events_lenz()


def run_multithread_elsa(arglist,ncpu):
    p = Pool(int(ncpu))
    print arglist
    p.map(run_elsa,arglist)
