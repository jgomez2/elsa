// system includes
//#include "RTypes.h"
#include <cstdlib>
#include "stdint.h"
// ROOT includes
#include "TH1.h"
#include "TF1.h"
#include <cmath>

// My defines
#define WF_SIZE 200
#define WF_SIZE_BIG 15000
#define MAXNBOARDS 3
#define MAXNCHANNELS 16
#define NBITSCLOCK 31
// Random Variables Class
using namespace std;

  typedef struct {
    uint64_t timestamp;  // timestamp
    uint32_t Ns;         // number of samples in waveform
    uint16_t sgate;      // short gate
    uint16_t lgate;      // long gate
    uint16_t baseline;   // baseline
    uint8_t board;       // board number
    uint8_t channel;     // channel number
  } DEVT_BANK;

bool reset_bank (DEVT_BANK *bank)
{
    bank->timestamp = 0;  // timestamp
    bank->sgate = 0;      // short gate
    bank->lgate = 0;      // long gate
    bank->baseline = 0;   // baseline
    bank->board = 256;       // board number
    bank->channel = 256;     // channel number
    bank->Ns = 0;         // number of samples in waveform
}

float dance_timing(UShort_t *in_array)
{
  int howmany = 96;
  int bl_start = 1;
  int bl_end = 10;
  int cfd_start = 75;
  int cfd_end = 95;

  float num = 0.0;
  float denom = 0.0;
  float baseline = 0.0;
  for (int i=bl_start;i<bl_end;++i)
  {
    num += in_array[i];
    denom += 1.0;
  }
  baseline = num/denom;
  num = 0.0;
  denom = 0.0;
  float cfd_level = 0.0;
  for (int i=cfd_start;i<cfd_end;++i)
  {
    num += in_array[i]-baseline;
    denom += 1.0;
  }
  cfd_level = 0.9*num/denom;
  //cout << "cfd leve " << cfd_level << " lbaseline " << baseline << endl;
  
  for (int i=1;i<howmany;++i)
  {
    //cout << "sample " << in_array[i]-baseline << endl;
    if (in_array[i]-baseline < cfd_level)
    {
      float interp = 0.0;
      float delta_dat = in_array[i]-in_array[i-1];
      float delta_cfd = cfd_level - (in_array[i-1]-baseline);
      //cout << "delta data " << delta_dat << " delta cfd " << delta_cfd << endl;
      float frac = delta_cfd/delta_dat;
      interp = frac+i-1;
      //cout << " interp " << interp << endl << endl;
      return interp;
    }
  }
  return -1.0;
}

bool IsValid(double par)
{
  if (par == -1)
  {
   return false;
  } else
  {
    return true;
  }
}
class TKESettings
{
    public:
      Int_t sample_delta[MAXNBOARDS][MAXNCHANNELS];// how far away to sample for each point of the derivative
      Int_t derive_point_delta[MAXNBOARDS][MAXNCHANNELS]; // how far away to move to calculate the derivative points
      Int_t npts_wf[MAXNBOARDS][MAXNCHANNELS]; // length of waveform 
      Int_t polarity[MAXNBOARDS][MAXNCHANNELS]; // 0 positive, 1 negative
      double zc_thresh[MAXNBOARDS][MAXNCHANNELS]; // zero crossing threshold (for double derivative) UNUSED RIGHT NOW
      Int_t offset_baseline[MAXNBOARDS][MAXNCHANNELS]; // offset to go into the baseline
      Int_t npts_baseline[MAXNBOARDS][MAXNCHANNELS]; // number of points to go in the baseline
      Int_t ph_method[MAXNBOARDS][MAXNCHANNELS]; // method to use to get "peak height": 0=charge integrate, 1=exp fit
      Int_t npts_ph[MAXNBOARDS][MAXNCHANNELS]; // number of points to go into calculatin the peak height
      Int_t offset_ph[MAXNBOARDS][MAXNCHANNELS]; // number of points to go into calculatin the peak height
      float ph_guess[MAXNBOARDS][MAXNCHANNELS];  // initial values for exponential fitting - peak height
      float decay_guess[MAXNBOARDS][MAXNCHANNELS]; // initial values for exponential fitting - decay constant

    public:
      void reset();
};
void TKESettings::reset()
{
    for (int i=0;i<MAXNBOARDS;++i)
    {
      for (int j=0;j<MAXNCHANNELS;++j)
      {
        sample_delta[i][j] = 4;
        derive_point_delta[i][j] = 10;
        npts_wf[i][j] = 4096;
        offset_baseline[i][j] = 1;
        npts_baseline[i][j] = 10;
      }
    }
}
class MidasEvent
{
    public:
        UChar_t bnum;
        UChar_t chnum;
        ULong64_t ts;
        UShort_t sgate; // 32 bit unsigned integer
        UShort_t lgate; // 32 bit unsigned integer
        UShort_t wf[WF_SIZE_BIG];   // raw waveform
        UShort_t wf2[WF_SIZE_BIG];   // raw waveform
        Byte_t dp0[WF_SIZE_BIG];   // digital probe 0
        Byte_t dp1[WF_SIZE_BIG];   // digital probe 0
        float flt0[WF_SIZE_BIG]; // filtered waveform 0 - mostly 1st derivative
        float flt1[WF_SIZE_BIG]; // filtered waveform 1 - mostly 2nd derivative

    public:
        void reset();
        Int_t process_wf(TKESettings tvar,double *t,double *peak,double* baseline);

        void double_derivative_filter(Int_t npts,UShort_t *inarray,float *outarray,Int_t sample_delta,Int_t derive_point_delta);
        void derivative_filter(Int_t npts,UShort_t *inarray,float *outarray,Int_t sample_delta,Int_t derive_point_delta);
        float locate_zerocrossing_neg(Int_t npts, float *inarray, float armthresh);
        float locate_zerocrossing_pos(Int_t npts, float *inarray, float armthresh);
        float locate_peak_neg(Int_t npts, float *inarray);
        float locate_peak_pos(Int_t npts, float *inarray);
        float calculate_baseline(Int_t npts,Int_t offset, UShort_t *inarray);
        float calculate_peakheight_ci(Int_t npts,Int_t offset,float baseline, UShort_t *inarray); // charge integration
        float calculate_peakheight_ci_grid(Int_t npts,Int_t offset,float baseline, UShort_t *inarray); // charge integration for grid
        float calculate_peakheight_exp(Int_t wfpts,Int_t npts,Int_t offset,float baseline, UShort_t *inarray,float ph_guess, float decay_guess); // exponential fitting
};

void MidasEvent::reset()
{
    bnum = 255;
    chnum = 255;
    ts = 0;
    for (int i=0;i<WF_SIZE_BIG;++i)
    {
        wf[i]=0;
        flt0[i]=0;
        flt1[i]=0;
    }
}


class ProcessedRootEvent
{
    public:
        UChar_t bnum;
        UChar_t chnum;
        ULong64_t ts;
        double t;
        double peak;
        UInt_t sgate;
        UInt_t lgate;
        double baseline;
        UShort_t wf[128];   // raw waveform

    public:
        void reset();
};

void ProcessedRootEvent::reset()
{
    bnum = 255;
    chnum = 255;
    ts = 0;
    t = -1;
    peak = -1;
    baseline = -1;
    sgate = 0;
    lgate = 0;
}

class GammaSingles
{
    public:
        UChar_t bnum;
        UChar_t chnum;
        double t0t; // t0
        double dt; // uncalibrated time of flight
        double ge;
        double baseline;
    public:
        void reset();
};

void GammaSingles::reset()
{
    bnum = 255;
    chnum = 255;
    t0t = -1;
    dt = -1;
    ge = -1;
    baseline = -1;
}

class DANCEBall
{
    public:
        double esum;
        Int_t mcrys;
        double t; // time of event
        double dte; // width of event
        double energies[40];
        double times[40];
        Short_t ids[40];
    public:
        void reset();
};

void DANCEBall::reset()
{
    esum = 0;
    mcrys = 0;
    t = -1;
    dte = -1;
    for (int eye=0;eye<40;++eye)
    {
        energies[eye]=-1;
        times[eye]=-1;
        ids[eye]=-1;
    }
}

class EventInfo
{
    public:
        double t0t; // t0
        double ct; // cathode time
        double at; // anode time rel. to t0
        double gt; // grid time rel. to t0
        double ft; // flux time rel. to t0
        double dt; // uncalibrated time of flight
        double TOF; // calibrated time of flight	
        double cph; // cathode pulse height
        double aph; // anode pulse height
        double gph; // grid pulse height
        double fph; // flux pulse height
        bool ccoinc; // have all charged particle channels fired
        double Etof[32]; // DSSD time
        double Ecoinc[32]; // Cathode & DSSD coinc time
        double Epeak[32]; //DSSD peak
        double Elgate[32]; //DSSD integral

  //double gt[2]; // uncalibrated time of flight for gamma detectors
        double ge[2]; // "pulse height" for gamma detectors
        double gtl[16]; // extra big list for wide coinc gates - time
        double gel[16]; // extra big list for wide coinc gates - energy
        double gchl[16]; // extra big list for wide coinc gates - channel
        double singles_timing[10];
        double singles_energy[10];
        double gwf2[128]; // extra big list for wide coinc gates - channel
        double gwf3[128]; // extra big list for wide coinc gates - channel
        double gwf4[128]; // extra big list for wide coinc gates - channel
        double gwf5[128]; // extra big list for wide coinc gates - channel
        double gwf6[128]; // extra big list for wide coinc gates - channel
        double gwf7[128]; // extra big list for wide coinc gates - channel
        int gmult;    // coincidence level for gamma detectors

    public:
        void reset();
};

void EventInfo::reset()
{
        t0t=-1; // t0
        ct=-1; // cathode time
	at=-1; // cathode time
        gt=-1; // cathode time
        ft=-1; // cathode time
        dt=-1; // uncalibrated time of flight
	TOF=-1;
        cph=-1; // cathode pulse height
	aph=-1;  // anode pulse height
	gph=-1; // grid pulse height
	fph=-1; // grid pulse height
        ccoinc = false;
        for (int i=0;i<2;++i)
        {
	  //aph[i]=-1; // anode pulse height
	  // gph[i]=-1; // grid pulse height
	  // gt[i] = -1;
            ge[i] = -1;
        }
        gmult = 0;
        for (int i=0;i<32;++i)
        {
	  Etof[i]=-1e6;
	  Ecoinc[i]=-2000;
	  Epeak[i]=-2000;
	  Elgate[i]=-2000;
        }
}

float MidasEvent::calculate_peakheight_exp(Int_t wfpts, Int_t npts,Int_t offset,float baseline, UShort_t *inarray,float ph_guess, float decay_guess)
{
    TH1F *hist = new TH1F("thist","thist",wfpts,0,wfpts);
    for (int i=0;i<wfpts;++i)
    {
        hist->SetBinContent(i+1,inarray[i]);
        //cout << "wf 4 fit " << i << " " << inarray[i] << endl;
    }

    TF1 *expf = new TF1("expf","[0]*exp(-(x-[3])/[1])+[2]",offset,offset+npts);
    expf->SetParameter(0,ph_guess);
    expf->FixParameter(1,decay_guess);
    expf->FixParameter(2,baseline);
    expf->FixParameter(3,offset);
    //cout << "ph guess " << ph_guess << " decay " << decay_guess << " baseline " << baseline << " offset " << offset << " end " << offset+npts << endl;

 //   hist->Fit("expf","RMQNO");
    //hist->Fit("expf");
    float peakheight = expf->GetParameter(0);

    hist->Delete();
    expf->Delete();
    //cout << "peakheight " << peakheight << endl;
    return peakheight;
}

float MidasEvent::calculate_peakheight_ci(Int_t npts,Int_t offset,float baseline, UShort_t *inarray)
{
    // this function performs a simple charge integration
    float integral = 0;
    for (int i=offset;i<npts+offset;++i)
    {
        integral += (inarray[i]-baseline);
        //cout << "amp " << inarray[i] << " baseline " << baseline << " integral " << integral << endl;
    }
    //cout << "final integral " << integral << endl;
    return integral;
}

float MidasEvent::calculate_peakheight_ci_grid(Int_t npts,Int_t offset,float baseline, UShort_t *inarray)
{
    // this function performs a simple charge integration
    float integral = 0;
    //cout << "CI:input npts " << npts << ",  offset " << offset << endl;
    for (int i=offset;i<npts+offset;++i)
      {
        if(inarray[i]<4850)
          {
        integral += (inarray[i]-baseline);
        //cout << "i= "<<i <<", " << inarray[i] << " baseline " << baseline << " integral " << integral << endl;
          }
    }
    integral = -integral;
    // cout << "final integral " << integral << endl;
    return integral;
}

float MidasEvent::calculate_baseline(Int_t npts,Int_t offset, UShort_t *inarray)
{
    float num = 0;
    float denom = 0;
    //cout << "input npts " << npts << "offset " << offset << endl;
    for (int i=offset;i<offset+npts;++i)
    {
        //cout << inarray[i] << endl;
        num += inarray[i];
        denom += 1.;
    }
    //cout << "num " << num << "denom " << denom << endl;
    return num/denom;
}

float MidasEvent::locate_peak_neg(Int_t npts, float *inarray)
{
    
    float mag = 0.;
    int index = 0;
    for (int i=0;i<npts;++i)
    {
        //std::cout << "wfa " << inarray[i] << std::endl;
        if (inarray[i]<mag)
        {
            index = i;
            mag = inarray[i];
        }
    }
    return mag;
}
float MidasEvent::locate_peak_pos(Int_t npts, float *inarray)
{
    
    float mag = 0.;
    int index = 0;
    for (int i=0;i<npts;++i)
    {
        if (inarray[i]>mag)
        {
            index = i;
            mag = inarray[i];
        }
    }
    return mag;
}
float MidasEvent::locate_zerocrossing_neg(Int_t npts, float *inarray, float armthresh)
{
    bool armed = false;
    for (int i=0;i<npts;++i)
    {
        if (!armed)
        {
            if (inarray[i]<armthresh)
            {
                armed = true;
                continue;
            }
        } else
        {
            if (inarray[i]>0.)
            {
                float y1 = inarray[i-1];
                float y2 = inarray[i];
                float deltay = y2-y1;
                float frac = abs(y1/deltay);
                float indexinterp = (float)i+frac;
                return indexinterp;
            }
        }
    }
    return 0;
}
float MidasEvent::locate_zerocrossing_pos(Int_t npts, float *inarray, float armthresh)
{
    bool armed = false;
    for (int i=0;i<npts;++i)
    {
        if (!armed)
        {
            if (inarray[i]>armthresh)
            {
                armed = true;
                continue;
            }
        } else
        {
            if (inarray[i]<0.)
            {
                float y1 = inarray[i-1];
                float y2 = inarray[i];
                float deltay = y2-y1;
                float frac = abs(y1/deltay);
                float indexinterp = (float)i+frac;
                return indexinterp;
            }
        }
    }
    return 0;
}

void MidasEvent::double_derivative_filter(Int_t npts,UShort_t *inarray,float *outarray,Int_t sample_delta,Int_t derive_point_delta)
{
    // global loop
    int gllim = (2.*sample_delta+1)+derive_point_delta;
    int gulim = npts - ((2.*sample_delta+1)+derive_point_delta);
    for (int i=gllim;i<gulim;++i)
    {
        float loavg = 0.;
        // calculate low average
        int lolim = i-1-2.*sample_delta-derive_point_delta;
        int uplim = i-1-derive_point_delta;
        for (int j=lolim;j<=uplim;++j)
        {
            loavg += inarray[j];
        }
        loavg /= 2.*sample_delta+1.;
        float hiavg = 0.;
        // calculate hi average
        lolim = i-1-2*sample_delta;
        uplim = i-1;
        for (int j=lolim;j<=uplim;++j)
        {
            hiavg += inarray[j];
        }
        hiavg /= 2.*sample_delta+1.;
        float derivative = (hiavg-loavg)/(2.*derive_point_delta+1.);
        // other derivative
        loavg = 0.;
        // calculate low average
        lolim = i+1;
        uplim = i+1+2*sample_delta;
        for (int j=lolim;j<=uplim;++j)
        {
            loavg += inarray[j];
        }
        loavg /= 2.*sample_delta+1.;
        hiavg = 0.;
        // calculate hi average
        lolim = i+1+derive_point_delta;
        uplim = i+1+2.*sample_delta+derive_point_delta;
        for (int j=lolim;j<=uplim;++j)
        {
            hiavg += inarray[j];
        }
        hiavg /= 2.*sample_delta+1.;
        float derivative2 = (hiavg-loavg)/(2.*derive_point_delta+1.);
        float second_derivative = (derivative2-derivative)/derive_point_delta;
        outarray[i] = second_derivative;
        //cout<< "i " << i << " original waveform " << inarray[i] << " flt " << outarray[i] << "deriv 1 " << derivative << "deriv 2 " << derivative2 << endl;
    }
};
void MidasEvent::derivative_filter(Int_t npts,UShort_t *inarray,float *outarray,Int_t sample_delta,Int_t derive_point_delta)
{
    // global loop
    for (int i=sample_delta+derive_point_delta;i<npts-sample_delta-derive_point_delta;++i)
    {
        float loavg = 0.;
        // calculate low average
        int lolim = (i-sample_delta-derive_point_delta);
        int uplim = i+sample_delta-derive_point_delta;
        for (int j=lolim;j<=uplim;++j)
        {
            //cout << "filtering in array " << inarray[j] << endl;
            loavg += inarray[j];
        }
        loavg /= 2.*sample_delta+1.;
        float hiavg = 0.;
        // calculate hi average
        lolim = i-sample_delta+derive_point_delta;
        uplim = i+sample_delta+derive_point_delta;
        for (int j=lolim;j<=uplim;++j)
        {
            hiavg += inarray[j];
        }
        hiavg /= 2.*sample_delta+1.;
        float derivative = (hiavg-loavg)/(2.*derive_point_delta+1.);
        outarray[i] = derivative;
    }
};


// this is a critical function - it defines  how all the wavefoms get processed
Int_t MidasEvent::process_wf(TKESettings tvar,double *t,double *peak, double *baseline)
{
  // filter the waveform
  //printf("about to derivative filter\n");
  derivative_filter(tvar.npts_wf[bnum][chnum],this->wf,this->flt0,tvar.sample_delta[bnum][chnum],tvar.derive_point_delta[bnum][chnum]);
  //printf("about to double derivative filter\n");
  double_derivative_filter(tvar.npts_wf[bnum][chnum],this->wf,this->flt1,tvar.sample_delta[bnum][chnum],tvar.derive_point_delta[bnum][chnum]);
  //printf("done with filters\n");
  // locate zero crossing (pretty good timing and definition of stuff for signal height calculations)
  double zcross = 0;
  //printf("bnum %i chnum %i\n",bnum,chnum);
  if (tvar.polarity[bnum][chnum] == 0)
  {
    //printf("zero crossing positive find\n");
    float bdelta = locate_peak_pos(tvar.npts_wf[bnum][chnum],this->flt1);
    zcross = locate_zerocrossing_pos(tvar.npts_wf[bnum][chnum],this->flt1,bdelta/2.);
    //std::cout << "zero crossing pos" << zcross << std::endl;
    //cout <<"b= "<<bnum<<", ch= "<<chnum <<", time= " << zcross << endl;
  }
  else if (tvar.polarity[bnum][chnum] == 1)
  {
    //std::cout << "points in waveform " << tvar.npts_wf[bnum][chnum] << std::endl;
    float bdelta = locate_peak_neg(tvar.npts_wf[bnum][chnum],this->flt1);
    //std::cout << "bdelta" << bdelta << std::endl;
    
    zcross = locate_zerocrossing_neg(tvar.npts_wf[bnum][chnum],this->flt1,bdelta/2.);
    //std::cout << "zero crossing neg" << zcross << std::endl;
  }
  *t = zcross;

  // find the baseline
  float bl = calculate_baseline(tvar.npts_baseline[bnum][chnum],tvar.offset_baseline[bnum][chnum],this->wf);
  //cout << "baseline " << bl << endl;
  *baseline = bl;
  
  //calculate the peak height
  float ph = 0;
  if (tvar.ph_method[bnum][chnum] == 0)
  {
    //Int_t offset_temp = ((Int_t)zcross)+tvar.offset_ph[bnum][chnum];
    Int_t offset_temp = tvar.offset_ph[bnum][chnum];
    ph = calculate_peakheight_ci(tvar.npts_ph[bnum][chnum],offset_temp,bl,this->wf);
    *peak = ph;
  }
  else if (tvar.ph_method[bnum][chnum] == 1)
  {
    //cout << "number of points in waveform " << tvar.npts_wf[bnum][chnum] << endl;
    ph = calculate_peakheight_exp(tvar.npts_wf[bnum][chnum],tvar.npts_ph[bnum][chnum],(Int_t)zcross+tvar.offset_ph[bnum][chnum],bl,this->wf,tvar.ph_guess[bnum][chnum],tvar.decay_guess[bnum][chnum]);
    *peak = ph;
    //cout << "baseline " << bl << "peak height " << *peak << endl;
  }
  else if (tvar.ph_method[bnum][chnum] == 2) // grid signal for lenz: for charge integral
  {
    //Int_t offset_temp = ((Int_t)zcross)+tvar.offset_ph[bnum][chnum];
    ph = calculate_peakheight_ci_grid(tvar.npts_ph[bnum][chnum],tvar.offset_ph[bnum][chnum],bl,this->wf);
    *peak = ph;
  }
  
  return 0;
}
