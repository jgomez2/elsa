#include "TFile.h"
#include "TString.h"
#include "TCutG.h"
#include "TH1.h"
#include "TH2.h"

#include <zlib.h> 
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Structures.cp"
#include "TTree.h"

#define MAXNDETCLASS 4
#define MAXNDETS 1024
#define MAXNCHANNELS 16
#define NBITSCLOCK 31

// try out a notion of how to define hists
#include "RootStuff/Hist_S0LENZ_Define.C"

class elsa_root {
    public:
        // variables
        TString *orfilename;
        TFile *orfile;
        TString *ipath;
        TString *opath;
        TString *ohpath;
        bool MidasEventPrint;
        Int_t MidasEventPrintThresh;
        Int_t STOPatEVENT;
        Int_t detid_shift;
        bool wfrm;
        bool dump0b;
        bool dump0r;
        bool dump0h;
        bool dump1r;
        bool dump1h;
        bool lsint_sub;
        Int_t ftoption;
        Int_t ftoptionind[MAXNDETS];
        Int_t t0ds;
        Int_t t0ch;
        Double_t t0to;
        Int_t extra_option; // extra word option in elsa stage 0 output
                            // 0 = off
                            // 1 = write pulse height
                            // 2 = write pulse height into charge integral
        bool upconvert; // whether to look for clock rollovers

        

        // data holders
        std::vector<ELSA_BANK> ebanks[MAXNDETCLASS];
        std::vector<short> extras[MAXNDETCLASS]; // optional extra word

        // detector descriptions
        int detclass[MAXNDETS];
        Double_t tslope[MAXNDETS];
        Double_t toffset[MAXNDETS];
        Double_t eslope[MAXNDETS];
        Double_t eoffset[MAXNDETS];
        float polarity[MAXNDETS]; // polarity +-1
        float time_filter_params[MAXNDETS][4]; // timing filter parameters
                                               // CFD: 0 = delay in samples
                                               //      1 = attenuation factor
                                               //      2 = number of samples to smooth
                                               //      3 = arm threshold
        float energy_filter_params[MAXNDETS][4]; // energy filter parameters - depends on other options
                                                 // CI: 0 = baseline low limit
                                                 //     1 = baseline upper limit
                                                 //     2 = integral low limit
                                                 //     3 = integral upper limit
        // detector thresholds / max
        Double_t threshold[MAXNDETS];
        Double_t saturate[MAXNDETS];
        Double_t psd1d[MAXNDETS];

        int lsgate_ds[MAXNDETS];
        TCutG lsgate;
        TGraph *walk;

        // functions
        Int_t execute_uac_qils(double interp_slope);
        Int_t execute_uac_pils(double interp_slope);
        Int_t execute_uac_cevt(double interp_slope);
        Int_t execute_s0bin(Int_t calibrate_time, Int_t calibrate_energy);
        Int_t execute_devt(Int_t dtype);
        Int_t execute_anna();

        Int_t build_events_chinu(Int_t t0ebuild);
        Int_t build_events_rpi();
        Int_t build_events_dance(Int_t t0ebuild);
        Int_t build_events_mona(Int_t t0ebuild);
        Int_t build_events_lenz();
        Int_t build_events_ndse(Int_t t0ebuild);

        Int_t locate_good_indices(int *i, int *j, Double_t coincwin, int class1, int class2);
        Int_t locate_t0i(int *t0i,int *i,int t0class,int class1);

        float digital_cfd(short *wf, float *fwf, uint32_t wflen, uint16_t detid);
        float digital_double_derivative(short *wf, float *fwf, uint32_t wflen, uint16_t detid);
        short find_peakheight(short *wf, uint32_t wflen, uint16_t detid);
        short integral_peakheight(short *wf, uint32_t wflen, uint16_t detid, short position);

        // auxiliary functions for MoNA
        Int_t create_barhits(MONA_EVENT* mevt);
        Int_t insert_pmtstrike(MONABAR_HIT *barhit, ELSA_BANK *pmthit, bool isfirst);
        // auxiliary variables for MoNA
        Int_t detid2side[MAXNDETS];
        Int_t detid2bar[MAXNDETS];



    private:
};

short elsa_root::find_peakheight(short *wf, uint32_t wflen, uint16_t detid)
{
    short baseline = wf[0];
    short pol = (short)polarity[detid];
    short ph = 0;
    //printf("in the peak height finder\n");
    for (size_t i=1;i<wflen;++i)
    {
        short temp_ph = pol*(wf[i] - baseline);
        if (temp_ph > ph)
        {
            ph = temp_ph;
        }

    }
    //printf("the peak height is %d\n",ph);
    return ph;
}

short elsa_root::integral_peakheight(short *wf, uint32_t wflen, uint16_t detid, short position)
{
    if (detid<16)
    {
        return 0;
    }
    //printf("\n\n\n\nwflen %i id %i position %i\n",wflen,detid,position);
    float scale = 4; // scaling factor for the pulse height
    float baseline = 0;
    float ci = 0;
    //printf("baseline limits %f %f\n",energy_filter_params[detid][0],energy_filter_params[detid][1]);
    for (short i=energy_filter_params[detid][0];i<energy_filter_params[detid][1]+1;++i)
    {
        baseline += (float)wf[i];
    }
    //printf("baseline %f\n",baseline);
    baseline /= (float)(energy_filter_params[detid][1]-energy_filter_params[detid][0]);

    //printf("position %i ci lo %f ci hi %f wflen %i\n",position,energy_filter_params[detid][2]+position,energy_filter_params[detid][3]+position,wflen);

    if (energy_filter_params[detid][2]+position<0 || energy_filter_params[detid][3]+position>=wflen)
    {
        return 0;
    }
    for (size_t i=energy_filter_params[detid][2]+position;i<energy_filter_params[detid][3]+position+1;++i)
    {
        ci += (float)wf[i];
    }
    //printf("ci %f\n",ci);
    ci /= (float)(energy_filter_params[detid][3]-energy_filter_params[detid][2]);

    float totint = polarity[detid]*scale*(ci-baseline);
    return (short)totint;
}

float elsa_root::digital_cfd(short *wf, float *fwf, uint32_t wflen, uint16_t detid)
{
    int delay = (int)time_filter_params[detid][0];
    float atten = time_filter_params[detid][1];
    int smooth = (int)time_filter_params[detid][2];
    float thresh = time_filter_params[detid][3];

    float norm = 1./smooth;

    float losum = 0;
    float hisum = 0;

    // initialization
    for (int i=0;i<smooth;++i)
    {
        losum += wf[i];
        hisum += wf[i+delay];
    }
    float offset = -1.0*norm*(atten*hisum - losum);

    bool armed = false;
    for (size_t i=smooth;i<wflen-delay-smooth;++i)
    {
        losum += wf[i];
        losum -= wf[i-smooth];
        hisum += wf[i+delay];
        hisum -= wf[i+delay-smooth];
        float filt = norm*(atten*hisum - losum) + offset;
        fwf[i] = polarity[detid]*filt;
        if (!armed)
        {
            if (fwf[i]>thresh)
            {
                armed = true;
            }
        } else
        {
            if (fwf[i]<0.0)
            {
                float ydelta = fwf[i-1] - fwf[i];
                float interp = (fwf[i-1]/ydelta) + i - 1;
                return interp;
            }
        }
    }
    return -1;
}

float elsa_root::digital_double_derivative(short *wf, float *fwf, uint32_t wflen, uint16_t detid)
{
    int delay = (int)time_filter_params[detid][0];
    int smooth = (int)time_filter_params[detid][2];
    float thresh = time_filter_params[detid][3];
    //printf("id %i delay %i smooth %i thresh %f\n",detid,delay,smooth,thresh);

    float norm = 1./smooth;

    float losum = 0;
    float hisum = 0;
    float losum2 = 0;
    float hisum2 = 0;

    // initialization
    for (int i=0;i<smooth;++i)
    {
        losum += wf[i];
        hisum += wf[i+delay];
        losum2 += wf[i+delay+smooth];
        hisum2 += wf[i+delay+delay+smooth];

        if (dump0h)
        {
            h_wf[detid]->Fill(i,wf[i]);
        }
    }

    for (int i=0;i<smooth;++i)
    {
        float filt = norm*((hisum2-losum2) - (hisum - losum));
        fwf[i] = polarity[detid]*filt;

        if (dump0h)
        {
            h_fwf[detid]->Fill(i,fwf[i]);
        }
    }

    bool armed = false;
    bool triggered = false;
    float maxdev = 0.0;
    float interp;
    for (size_t i=smooth;i<wflen-2*delay-2*smooth;++i)
    {
        losum += wf[i];
        losum -= wf[i-smooth];
        hisum += wf[i+delay];
        hisum -= wf[i+delay-smooth];

        losum2 += wf[i +delay+smooth];
        losum2 -= wf[i-smooth +delay+smooth];
        hisum2 += wf[i+delay +delay+smooth];
        hisum2 -= wf[i+delay-smooth +delay+smooth];


        float filt = norm*((hisum2-losum2) - (hisum - losum));

        fwf[i] = polarity[detid]*filt;

        if (fwf[i]>maxdev)
        {
            maxdev = fwf[i];
        }

        if (dump0h)
        {
            h_wf[detid]->Fill(i,wf[i]);
            h_fwf[detid]->Fill(i,fwf[i]);
        }

        if (!armed)
        {
            if (fwf[i]>thresh)
            {
                armed = true;
            }
        } else
        {
            if (fwf[i]<0.0 && !triggered)
            {
                float ydelta = fwf[i-1] - fwf[i];
                interp = (fwf[i-1]/ydelta) + i - 1;
                triggered = true;
                if (dump0h)
                {
                    h_mdev->Fill(detid,maxdev);
                    h_ppos->Fill(i,detid);
                }
                else
                {
                    return interp;
                }
            }
        }
    }
    return interp;
}

Int_t elsa_root::locate_good_indices(int *i, int *j, Double_t coincwin, int class1, int class2)
{
    //printf("i am in locate_good_indices %d %d\n",*i,*j);
    //*i += 1;
    if (*i>(int)ebanks[class1].size()-1)
    {
        //printf("exit 0\n");
        return 0;
    }
    while (true)
    {
        if (*j<0)
        {
            *j = 0;
            //printf("exit 1\n");
            return 1;
        }
        Double_t deltat = ebanks[class2][*j].time - ebanks[class1][*i].time;
        if (deltat < coincwin)
        {
            //printf("exit 2\n");
            return 2;
        } else
        {
            *j -= 1;
        }
    }
    //printf("exit 3\n");
}

Int_t elsa_root::locate_t0i(int *t0i,int *i,int t0class,int class1)
{
    while (true)
    {
        // at the end of the t0 list... can't go farther
        if ((*t0i)+1>(int)ebanks[t0class].size()-1)
        {
            return 0;
        }
        if (ebanks[class1][*i].time<ebanks[t0class][*t0i].time)
        {
            *t0i -= 1;
            if (*t0i < 0)
            {
                *t0i = 0;
                return 1;
            }
        } else if (ebanks[class1][*i].time>=ebanks[t0class][*t0i].time && ebanks[class1][*i].time<ebanks[t0class][(*t0i)+1].time)
        {
            return 2;
        } else
        {
            *t0i += 1;
        }

    }
}

bool tsort(ELSA_BANK i, ELSA_BANK j) {return (i.time<j.time); }

bool mona_barsort(MONABAR_HIT i, MONABAR_HIT j) {return (i.tmean<j.tmean); }

Int_t reset_barhit(MONABAR_HIT *mhit)
{
  mhit->tmean = -1;
  mhit->tleft = -1;
  mhit->tright = -1;
  mhit->qleft = -1;
  mhit->qright = -1;
  mhit->pleft = -1;
  mhit->pright = -1;
  mhit->barnum = -1;
  mhit->iscomplete = false;
  
  return 0;
}

Int_t reset_mona_event(MONA_EVENT *mevt)
{
  for (int i=0;i<NMONABARS;++i)
  {
    mevt->tmean[i]=-1;
    mevt->qmean[i]=-1;
    mevt->qleft[i]=-1;
    mevt->qright[i]=-1;
    mevt->pleft[i]=-1;
    mevt->pright[i]=-1;
    mevt->tdiff[i]=-1;
    mevt->bar[i]=-1;
    mevt->x[i]=-1;
    mevt->y[i]=-1;
    mevt->z[i]=-1;
  }
  mevt->coincevts.clear();
  mevt->mult = 0;;

  return 0;
}

Int_t elsa_root::insert_pmtstrike(MONABAR_HIT *barhit, ELSA_BANK *pmthit, bool isfirst)
{
  Int_t leftright = detid2side[pmthit->detector_id];
  Int_t barnum = detid2bar[pmthit->detector_id];
  if (isfirst)
  {
    if (leftright==0) // left is 0
    {
      barhit->tleft = pmthit->time;
      barhit->qleft = pmthit->integral[1];
      barhit->pleft = pmthit->integral[0];
    }
    else
    {
      barhit->tright = pmthit->time;
      barhit->qright = pmthit->integral[1];
      barhit->pright = pmthit->integral[0];
    }
    barhit->barnum = barnum;
  } else
  {
    // check that this pmt matches this bar
    if (barnum == barhit->barnum)
    {
      //printf("found common bar %d\n",barnum);
      if (leftright==0)
      {
        barhit->tleft = pmthit->time;
        barhit->qleft = pmthit->integral[1];
        barhit->pleft = pmthit->integral[0];
      } else
      {
        barhit->tright = pmthit->time;
        barhit->qright = pmthit->integral[1];
        barhit->pright = pmthit->integral[0];
      }
      barhit->tmean = (barhit->tright + barhit->tleft)/2.;
    } else
    {
      return 0;
    }

  }
  
  return 1;
}

Int_t elsa_root::build_events_mona(Int_t t0ebuild)
{
    orfile = new TFile(opath->Data(),"recreate");
    TTree *otree = new TTree("tree","tree");
    MONA_EVENT mevt;
    Double_t t0time; // time of the t0
    Double_t evtime; // time of first light of event
    Double_t t0tof; // time of flight from t0 to first light in mona
    int t0i = 0;
    char buf[1024];
    otree->Branch("t0t",&t0time);
    otree->Branch("et",&evtime);
    otree->Branch("mult",&mevt.mult);
    sprintf(buf,"tmean[%i]",NMONABARS);
    otree->Branch(buf,mevt.tmean);
    sprintf(buf,"qmean[%i]",NMONABARS);
    otree->Branch(buf,mevt.qmean);
    sprintf(buf,"qleft[%i]",NMONABARS);
    otree->Branch(buf,mevt.qleft);
    sprintf(buf,"qright[%i]",NMONABARS);
    otree->Branch(buf,mevt.qright);
    sprintf(buf,"pleft[%i]",NMONABARS);
    otree->Branch(buf,mevt.pleft);
    sprintf(buf,"pright[%i]",NMONABARS);
    otree->Branch(buf,mevt.pright);
    sprintf(buf,"tdiff[%i]",NMONABARS);
    otree->Branch(buf,mevt.tdiff);
    sprintf(buf,"bar[%i]",NMONABARS);
    otree->Branch(buf,mevt.bar);
    sprintf(buf,"x[%i]",NMONABARS);
    otree->Branch(buf,mevt.x);
    sprintf(buf,"y[%i]",NMONABARS);
    otree->Branch(buf,mevt.y);
    sprintf(buf,"z[%i]",NMONABARS);
    otree->Branch(buf,mevt.z);
    
    // sort everyone first
    std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    std::sort(ebanks[2].begin(),ebanks[2].end(),tsort);
    printf("all events sorted\n");

    Double_t coincwin = 100.0;
    Int_t t0class = 0;
    Int_t class1 = 1;

    double tstart = ebanks[1][0].time;
    Int_t istart = 0;
    mevt.coincevts.push_back(ebanks[1][0]);
    for (int i=1;i<ebanks[1].size();++i)
    {
      //printf("this event time %f start time %f\n",ebanks[1][i].time,tstart);
      if (ebanks[1][i].time>tstart + coincwin) // out of coinc
      {
        //printf("creating barhits from size %d\n",mevt.coincevts.size());
        Int_t stat = create_barhits(&mevt);
        if (mevt.mult>0)
        {
          //printf("found a multiplicity of %d\n",mevt.mult);
        }
        //printf("bar hits created status %d\n",stat);
        if (!stat)
        {
          // looking for a good t0
          if (t0ebuild)
          {
              locate_t0i(&t0i,&istart,t0class,class1);
              t0time = ebanks[t0class][t0i].time;
              //printf("event index %d of %d\n",i,ebanks[class1].size());
              evtime = ebanks[class1][i].time;
              t0tof = evtime-t0time;
          }
          //printf("about to fill\n");
          otree->Fill();
          //printf("filled\n");
        }
        reset_mona_event(&mevt);
        tstart = ebanks[1][i].time;
        istart = i;
        mevt.coincevts.push_back(ebanks[1][i]);
      }
      else // inside coinc
      {
        mevt.coincevts.push_back(ebanks[1][i]);
      }
    }
    otree->Write();
    orfile->Close();
    return 0;
}

Int_t elsa_root::build_events_dance(Int_t t0ebuild)
{
    printf("i am in the dance event builder\n");
    orfile = new TFile(opath->Data(),"recreate");
    TTree *otree = new TTree("tree","tree");
    Int_t detid1;
    Int_t detid2;
    Double_t t0tof;
    Double_t tof;
    Double_t d1ph;
    Double_t d2ph;
    Double_t t0time;
    Double_t ntof1;
    Double_t ntof2;
    otree->Branch("detector_id1",&detid1);
    otree->Branch("detector_id2",&detid2);
    otree->Branch("t0tof",&t0tof);
    otree->Branch("ntof",&tof);
    otree->Branch("ntof1",&ntof1);
    otree->Branch("ntof2",&ntof2);
    otree->Branch("t0t",&t0time);
    otree->Branch("d1ph",&d1ph);
    otree->Branch("d2ph",&d2ph);
    
    // sort everyone first
    //std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    //std::sort(ebanks[2].begin(),ebanks[2].end(),tsort);
    printf("all events sorted\n");

    int t0i = 0;
    //int i = 0;
    int j = 0;
    Double_t coincwin = 1000.0;
    int t0class = 0;
    int class1 = 1;
    int class2 = 2;
    Double_t dancetime;

    int interesting_counter = 0;
    printf("about to try coinc\n");
    
    for (int i=0;i<ebanks[1].size()-2;++i)
    {
        if (t0ebuild)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            dancetime = ebanks[class1][i].time;
            t0tof = dancetime-t0time;
        }
        Double_t deltat = ebanks[1][i].time - ebanks[1][i+1].time;
        if (deltat > -1.0*coincwin && deltat < coincwin || true)
        {
            //printf("in coinc\n");
            detid1 = ebanks[1][i].detector_id;
            detid2 = ebanks[1][i+1].detector_id;
            d1ph = ebanks[1][i].integral[1];
            d2ph = ebanks[1][i+1].integral[1];
            tof = deltat;
            ntof1 = ebanks[1][i].time;
            ntof2 = ebanks[1][i+1].time;
            otree->Fill();
        }
        
    }

    /*
    while (true)
    {
        if (i>(int)ebanks[class1].size()-1)
        {
            break;
        }
        // looking for a good t0
        if (t0ebuild)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            ppactime = ebanks[class1][i].time;
            t0tof = ppactime-t0time;
        }
        // loop for neutron detectors
        while (true)
        {
            if (j>(int)ebanks[class2].size()-1)
            {
                break;
            }
            Double_t deltat = ebanks[class2][j].time - ebanks[class1][i].time;
            if (deltat < -1.0*coincwin) // outside low bound of coinc win
            {
                j += 1;
                continue;
            } else if (deltat > -1.0*coincwin && deltat < coincwin) // inside the coincidence window
            {
                detid1 = ebanks[class1][i].detector_id;
                detid2 = ebanks[class2][j].detector_id;
                tof = deltat;
                pph = ebanks[class1][i].integral[1];
                nph = ebanks[class2][j].integral[1];
                otree->Fill();
                j += 1;
            } else // outside of upper boundary of coincidence window
            {
                locate_good_indices(&i,&j,coincwin,class1,class2);
                break;
            }
        }
        i+=1;
        if ((j>(int)ebanks[class2].size()-1) && (i>(int)ebanks[class1].size()-1))
        {
            break;
        }
        //printf("end of i loop %d of %d\n",i,ebanks[class1].size());
    }
    */
    otree->Write();
    orfile->Close();
    return 0;
}

Int_t elsa_root::create_barhits(MONA_EVENT* mevt)
{
   if (mevt->coincevts.size()<2)
   {
     //printf("coinc events too small\n");
     return 1;
   }
   
   std::vector<MONABAR_HIT> bar_hits;
   MONABAR_HIT temphit;
   reset_barhit(&temphit);
   bar_hits.push_back(temphit);
   insert_pmtstrike(&bar_hits[0],&(mevt->coincevts[0]),true);
   
   for (size_t i=1;i<mevt->coincevts.size();++i)
   {
     Int_t success = 0;
     for (size_t j=0;j<bar_hits.size();++j)
     {
       Int_t success = insert_pmtstrike(&bar_hits[j], &(mevt->coincevts[i]),false);
       if (success)
       {
         bar_hits[j].iscomplete = true;
         break;
       }
     }
     if (!success)
     {
         reset_barhit(&temphit);
         bar_hits.push_back(temphit);
	 insert_pmtstrike(&bar_hits[bar_hits.size()-1],&(mevt->coincevts[i]),true);
         success = 0;
     }
   }

   std::sort(bar_hits.begin(),bar_hits.end(),mona_barsort);

   mevt->mult = 0;
   for (size_t i=0;i<bar_hits.size();++i)
   {
     if (i>NMONABARS)
     {
       printf("have more coincidences than mona bars %zu of %i from %i pmts\n",bar_hits.size(),NMONABARS,mevt->coincevts.size());
       break;
     }
     if (bar_hits[i].iscomplete)
     {
       mevt->qleft[mevt->mult] = bar_hits[i].qleft;
       mevt->qright[mevt->mult] = bar_hits[i].qright;
       mevt->pleft[mevt->mult] = bar_hits[i].pleft;
       mevt->pright[mevt->mult] = bar_hits[i].pright;
       mevt->tmean[mevt->mult] = bar_hits[i].tmean;
       mevt->qmean[mevt->mult] = (bar_hits[i].qleft + bar_hits[i].qright)/2.;
       mevt->tdiff[mevt->mult] = bar_hits[i].tleft - bar_hits[i].tright;
       mevt->bar[mevt->mult] = bar_hits[i].barnum;
       if (mevt->bar[mevt->mult]==0) {
		mevt->x[mevt->mult]=7.52*mevt->tdiff[mevt->mult] + -580.5;
		mevt->y[mevt->mult]=-10.254;
		mevt->z[mevt->mult]=105.12;
       } else if (mevt->bar[mevt->mult]==1) {
		mevt->x[mevt->mult]=7.43*mevt->tdiff[mevt->mult] + -364.1;
		mevt->y[mevt->mult]=0.;
		mevt->z[mevt->mult]=105.12;
       } else if (mevt->bar[mevt->mult]==2) {
		mevt->x[mevt->mult]=7.51*mevt->tdiff[mevt->mult] + -362.6;
		mevt->y[mevt->mult]=10.254;
		mevt->z[mevt->mult]=105.12;
       } else if (mevt->bar[mevt->mult]==3) {
		mevt->x[mevt->mult]=7.51*mevt->tdiff[mevt->mult] + -376.1;
		mevt->y[mevt->mult]=-10.254;
		mevt->z[mevt->mult]=115.374;
       } else if (mevt->bar[mevt->mult]==4) {
		mevt->x[mevt->mult]=7.51*mevt->tdiff[mevt->mult] + -369.3;
		mevt->y[mevt->mult]=0;
		mevt->z[mevt->mult]=115.374;
       } else if (mevt->bar[mevt->mult]==5) {
		mevt->x[mevt->mult]=7.43*mevt->tdiff[mevt->mult] + -351.9;
		mevt->y[mevt->mult]=10.254;
		mevt->z[mevt->mult]=115.374;
       } else if (mevt->bar[mevt->mult]==6) {
		mevt->x[mevt->mult]=7.43*mevt->tdiff[mevt->mult] + -361.9;
		mevt->y[mevt->mult]=-20.508;
		mevt->z[mevt->mult]=125.628;
       } else if (mevt->bar[mevt->mult]==7) {
		mevt->x[mevt->mult]=7.38*mevt->tdiff[mevt->mult] + -353.7;
		mevt->y[mevt->mult]=-10.254;
		mevt->z[mevt->mult]=125.628;
       } else if (mevt->bar[mevt->mult]==8) {
		mevt->x[mevt->mult]=7.72*mevt->tdiff[mevt->mult] + -367.7;
		mevt->y[mevt->mult]=0.;
		mevt->z[mevt->mult]=125.628;
       } else if (mevt->bar[mevt->mult]==9) {
		mevt->x[mevt->mult]=7.43*mevt->tdiff[mevt->mult] + -366.3;
		mevt->y[mevt->mult]=10.254;
		mevt->z[mevt->mult]=125.628;
       } else if (mevt->bar[mevt->mult]==10) {
		mevt->x[mevt->mult]=7.47*mevt->tdiff[mevt->mult] + -355.4;
		mevt->y[mevt->mult]=20.508;
		mevt->z[mevt->mult]=125.628;
       } else if (mevt->bar[mevt->mult]==11) {
		mevt->x[mevt->mult]=7.51*mevt->tdiff[mevt->mult] + -375.0;
		mevt->y[mevt->mult]=-20.628;
		mevt->z[mevt->mult]=135.882;
       } else if (mevt->bar[mevt->mult]==12) {
		mevt->x[mevt->mult]=7.42*mevt->tdiff[mevt->mult] + -356.2;
		mevt->y[mevt->mult]=-10.254;
		mevt->z[mevt->mult]=135.882;
       } else if (mevt->bar[mevt->mult]==13) {
		mevt->x[mevt->mult]=7.42*mevt->tdiff[mevt->mult] + -351.7;
		mevt->y[mevt->mult]=0;
		mevt->z[mevt->mult]=135.882;
       } else if (mevt->bar[mevt->mult]==14) {
		mevt->x[mevt->mult]=7.55*mevt->tdiff[mevt->mult] + -363.9;
		mevt->y[mevt->mult]=10.254;
		mevt->z[mevt->mult]=135.882;
       } else if (mevt->bar[mevt->mult]==15) {
		mevt->x[mevt->mult]=7.39*mevt->tdiff[mevt->mult] + -359.3;
		mevt->y[mevt->mult]=20.508;
		mevt->z[mevt->mult]=135.882;
       }		
       mevt->mult += 1;
     }
   }

   return 0;
}


Int_t elsa_root::build_events_rpi()
{
    orfile = new TFile(opath->Data(),"recreate");
    TTree *otree = new TTree("tree","tree");
    Int_t detid1;
    Double_t tof;
    Double_t t0ph;
    Double_t sgate;
    Double_t lgate;
    Double_t t0time;
    Double_t ntime;;
    otree->Branch("detid",&detid1);
    otree->Branch("tof",&tof);
    otree->Branch("t0t",&t0time);
    otree->Branch("nt",&ntime);
    otree->Branch("t0ph",&t0ph);
    otree->Branch("sgate",&sgate);
    otree->Branch("lgate",&lgate);
    
    // sort everyone first
    std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    printf("all rpi events sorted\n");

    int t0i = 0;
    int i = 0;
    int t0class = 0;
    int class1 = 1;

    int interesting_counter = 0;

    while (true)
    {
        if (i>(int)ebanks[class1].size()-1)
        {
            //printf("on index %d of %d\n",i,ebanks[class1].size());
            break;
        }
        // looking for a good t0
        if (1)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            ntime = ebanks[class1][i].time;
            tof = ntime-t0time;
            detid1 = ebanks[class1][i].detector_id;
            sgate = ebanks[class1][i].integral[0];
            lgate = ebanks[class1][i].integral[1];
            t0ph = ebanks[t0class][t0i].integral[1];
            otree->Fill();
        }
        i+=1;
    }
    otree->Write();
    orfile->Close();
    return 0;
}

Int_t elsa_root::build_events_lenz()
{
    if (dump1r)
    {
        orfile = new TFile(opath->Data(),"recreate");
    }
    TFile *ohfile;
    TTree *otree = new TTree("tree","tree");
    Int_t detid1;
    Double_t tof;
    Double_t tof_wrapped;
    Double_t t0ph;
    Double_t sgate;
    Double_t lgate;
    Double_t t0time;
    Double_t ntime;;
    if (dump1r)
    {
        otree->Branch("detid",&detid1);
        otree->Branch("tof",&tof);
        otree->Branch("t0t",&t0time);
        otree->Branch("nt",&ntime);
        otree->Branch("t0ph",&t0ph);
        otree->Branch("sgate",&sgate);
        otree->Branch("lgate",&lgate);
    }
    #include "RootStuff/Hist_S1LENZ_Define.C"
    if (dump1h)
    {
        printf("opening root histo file %s\n",ohpath->Data());
        ohfile = new TFile(ohpath->Data(),"recreate");
        #include "RootStuff/Hist_S1LENZ_Create.C"
    }
    
    // sort everyone first
    std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    std::sort(ebanks[2].begin(),ebanks[2].end(),tsort);
    printf("all lenz events sorted\n");

    int t0i = 0;
    int i = 0;
    int t0class = 0;
    int class1 = 2;

    int interesting_counter = 0;

    while (true)
    {
        if (i>(int)ebanks[class1].size()-1)
        {
            //printf("on index %d of %d\n",i,ebanks[class1].size());
            break;
        }
        // looking for a good t0
        if (1)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            ntime = ebanks[class1][i].time;
            tof = ntime-t0time;
            tof_wrapped = fmod(tof,1788.819875776);
            detid1 = ebanks[class1][i].detector_id;
            sgate = ebanks[class1][i].integral[0];
            lgate = ebanks[class1][i].integral[1];
            t0ph = ebanks[t0class][t0i].integral[1];
            if (dump1r)
            {
                otree->Fill();
            }
            if (dump1h)
            {
                #include "RootStuff/Hist_S1LENZ_Fill.C"
            }
        }
        i+=1;
    }
    if (dump1r)
    {
        otree->Write();
    }
    if (dump1h)
    {
        #include "RootStuff/Hist_S1LENZ_Write.C"
        ohfile->Close();
    }
    if (dump1r)
    {
        orfile->Close();
    }
    return 0;
}

Int_t elsa_root::build_events_chinu(Int_t t0ebuild)
{
    Int_t detid1;
    Int_t detid2;
    Double_t t0tof;
    Double_t t0tof_wrapped;
    Double_t tof;
    Double_t pph;
    Double_t nph;
    Double_t nsg;
    Double_t t0time;
    Double_t ppactime;;
    TTree *otree;
    TFile *ohfile;
    Double_t coincwin_lo = -2000000.0;
    Double_t coincwin_hi = 2000000.0;


    printf("dump1r %i dump1h %i\n",dump1r,dump1h);
    if (dump1r)
    {
        orfile = new TFile(opath->Data(),"recreate");
        otree = new TTree("tree","tree");
        otree->Branch("detector_id1",&detid1);
        otree->Branch("detector_id2",&detid2);
        otree->Branch("t0tof",&t0tof);
        otree->Branch("ntof",&tof);
        otree->Branch("t0t",&t0time);
        otree->Branch("ppt",&ppactime);
        otree->Branch("pph",&pph);
        otree->Branch("nsg",&nsg);
        otree->Branch("nph",&nph);
    }
    #include "RootStuff/Hist_S1CHINU_Define.C"
    if (dump1h)
      {
        printf("opening root histo file %s\n",ohpath->Data());
        ohfile = new TFile(ohpath->Data(),"recreate");
        #include "RootStuff/Hist_S1CHINU_Create.C"
      }
    
    // sort everyone first
    std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    std::sort(ebanks[2].begin(),ebanks[2].end(),tsort);
    printf("all events sorted\n");

    int t0i = 0;
    int i = 0;
    int j = 0;
    int t0class = 0;
    int class1 = 1;
    int class2 = 2;
        
    if (t0ebuild && ebanks[t0class].size()==0)
    {
      return -1;
    }
    if (ebanks[class1].size()==0 || ebanks[class2].size()==0)
    {
      return -2;
    }

    while (true)
    {
        if (i>(int)ebanks[class1].size()-1)
        {
            break;
        }
        // looking for a good t0
        if (t0ebuild)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            ppactime = ebanks[class1][i].time;
            t0tof = ppactime-t0time;
            //t0tof_wrapped = ((100000.*(t0tof-0.0)%(100000*1788.819875776))/100000.);
            t0tof_wrapped = fmod(t0tof,1788.819875776);
        }

	// loop for neutron detectors
        while (true)
        {
            if (j>(int)ebanks[class2].size()-1)
            {
                break;
            }
	    
	    Double_t deltat = ebanks[class2][j].time - ebanks[class1][i].time;
            if (deltat < coincwin_lo) // outside low bound of coinc win
            {
                j += 1;
                continue;
            } else if (deltat > coincwin_lo && deltat < coincwin_hi) // inside the coincidence window
            {
                detid1 = ebanks[class1][i].detector_id;
                detid2 = ebanks[class2][j].detector_id;
                tof = deltat;
                pph = ebanks[class1][i].integral[1];
                nsg = ebanks[class2][j].integral[0];
                nph = ebanks[class2][j].integral[1];
                if (dump1r)
                {
                    otree->Fill();
                }
                if (dump1h)
                {
                    if (pph>threshold[detid1]&&pph<saturate[detid1]&&nph>threshold[detid2]&&nph<saturate[detid2])
                    {
                      if (t0tof<((Double_t)t0ds)*1.3*1800.0)
                      {
                        #include "RootStuff/Hist_S1CHINU_Fill.C"
                      }
                    }
                }
                j += 1;
            } else // outside of upper boundary of coincidence window
            {
                locate_good_indices(&i,&j,coincwin_lo,class1,class2);
                break;
            }
        }
        i+=1;
        if ((j>(int)ebanks[class2].size()-1) && (i>(int)ebanks[class1].size()-1))
        {
            break;
        }
        //printf("end of i loop %d of %d\n",i,ebanks[class1].size());
    }

    
    if (dump1r)
    {
      otree->Write();
      orfile->Close();
    }
    if (dump1h)
    {
     #include "RootStuff/Hist_S1CHINU_Write.C"
        ohfile->Close();
    }

    return 0;
}


Int_t elsa_root::build_events_ndse(Int_t t0ebuild)
{
    Int_t detid1;
    Int_t detid2;
    Double_t t0tof;
    Double_t t0tof_wrapped;
    Double_t tof;
    Double_t pph;
    Double_t nph;
    Double_t nsg;
    Double_t t0time;
    Double_t ppactime;;
    TTree *otree;
    TFile *ohfile;
    Double_t coincwin_lo = -2000000.0;
    Double_t coincwin_hi = 2000000.0;

    //HISTOS FOR NDSE
    TH1D *liqdeadtime[10];
    TH1D *bgodeadtime[1];
    double lasttime[10];
    double blasttime[1];
    double singlelasttime;
    int lastdet;
    TH1D* liqbackground[10];//THIS LINE WILL BE USEFUL LATER
    TH1D* bgobackground[10];//THIS LINE WILL ALSO BE USEFUL LATER
    TH1D* rossialpha_liquid[10];
    TH1D* psdliquid2 = new TH1D("psd","psd",1000,0,.3);
    TH1D* shortraliquid;
    double ralasttime[10];
    //END OF NDSE PACIFIC STUFF


    printf("dump1r %i dump1h %i\n",dump1r,dump1h);
    if (dump1r)
    {
        orfile = new TFile(opath->Data(),"recreate");
        otree = new TTree("tree","tree");
        otree->Branch("detector_id1",&detid1);
        otree->Branch("detector_id2",&detid2);
        otree->Branch("t0tof",&t0tof);
        otree->Branch("ntof",&tof);
        otree->Branch("t0t",&t0time);
        otree->Branch("ppt",&ppactime);
        otree->Branch("pph",&pph);
        otree->Branch("nsg",&nsg);
        otree->Branch("nph",&nph);
	for (int iter=0;iter<10;iter++)
	  {
	    if (iter<1)
	      {
		char btitle[128];
		sprintf(btitle,"DeadTime_BGO%02i",iter);
		bgodeadtime[iter]=new TH1D(btitle,btitle,300000,0,coincwin_hi);
		
		char bgotitle[128];
		sprintf(bgotitle,"td1_05_back_sngl_ppac%02i",iter+1);
		bgobackground[iter]=new TH1D(bgotitle,bgotitle,150000,0,coincwin_hi);
	      }
	    if(iter>=2 && iter<9)
	      {
		char title[128];
		sprintf(title,"DeadTime_Liquid%02i",iter);
		liqdeadtime[iter]=new TH1D(title,title,300000,0,coincwin_hi);
		
		//background histograms must have same bin width as the final histograms
		char bgtitle[128];
		sprintf(bgtitle,"td1_05_back_sngl_ligl%02i",iter);
		liqbackground[iter]=new TH1D(bgtitle,bgtitle,450000,coincwin_lo,2*coincwin_hi);

		//rossi alpha histograms
		char ratitle[128];
		sprintf(ratitle,"RossiAlpha_liquid%02i",iter);
		rossialpha_liquid[iter]=new TH1D(ratitle,ratitle,2000050,-50,coincwin_hi);
	      }
	  }
	shortraliquid=new TH1D("ShortTimeRossiAlpha_AllLiquids","ShortTimeRossiAlpha_AllLiquids",2000050,-50,coincwin_hi);
    }
    #include "RootStuff/Hist_S1CHINU_Define.C"
    if (dump1h)
    {
        printf("opening root histo file %s\n",ohpath->Data());
        ohfile = new TFile(ohpath->Data(),"recreate");
        #include "RootStuff/Hist_S1CHINU_Create.C"
    }
    
    // sort everyone first
    std::sort(ebanks[0].begin(),ebanks[0].end(),tsort);
    std::sort(ebanks[1].begin(),ebanks[1].end(),tsort);
    std::sort(ebanks[2].begin(),ebanks[2].end(),tsort);
    printf("all events sorted\n");

    int t0i = 0;
    int i = 0;
    int j = 0;
    int t0class = 0;
    int class1 = 1;
    int class2 = 2;
    for(int iter=0;iter<10;iter++)
      {
	ralasttime[iter]=-100.0;
	if(iter>=2) lasttime[iter]=-100.0;
	if(iter<1) blasttime[iter]=-100.0;
	singlelasttime=-100;
	lastdet=-100;
      }
    
    if (t0ebuild && ebanks[t0class].size()==0)
    {
      return -1;
    }
    if (ebanks[class1].size()==0 || ebanks[class2].size()==0)
    {
      return -2;
    }
    /*    
    while (true)
    {
        if (i>(int)ebanks[class1].size()-1)
        {
            break;
        }
        // looking for a good t0
        if (t0ebuild)
        {
            locate_t0i(&t0i,&i,t0class,class1);
            t0time = ebanks[t0class][t0i].time;
            ppactime = ebanks[class1][i].time;
            t0tof = ppactime-t0time;
            //t0tof_wrapped = ((100000.*(t0tof-0.0)%(100000*1788.819875776))/100000.);
            t0tof_wrapped = fmod(t0tof,1788.819875776);
        }
	
        


	// loop for neutron detectors
        while (true)
        {
            if (j>(int)ebanks[class2].size()-1)
            {
                break;
            }
	    
	    Double_t deltat = ebanks[class2][j].time - ebanks[class1][i].time;
            if (deltat < coincwin_lo) // outside low bound of coinc win
            {
                j += 1;
                continue;
            } else if (deltat > coincwin_lo && deltat < coincwin_hi) // inside the coincidence window
            {
                detid1 = ebanks[class1][i].detector_id;
                detid2 = ebanks[class2][j].detector_id;
                tof = deltat;
                pph = ebanks[class1][i].integral[1];
                nsg = ebanks[class2][j].integral[0];
                nph = ebanks[class2][j].integral[1];
                if (dump1r)
                {
                    otree->Fill();
                }
                if (dump1h)
                {
                    if (pph>threshold[detid1]&&pph<saturate[detid1]&&nph>threshold[detid2]&&nph<saturate[detid2])
                    {
                      if (t0tof<((Double_t)t0ds)*1.3*1800.0)
                      {
                        #include "RootStuff/Hist_S1CHINU_Fill.C"
                      }
                    }
                }
                j += 1;
            } else // outside of upper boundary of coincidence window
            {
                locate_good_indices(&i,&j,coincwin_lo,class1,class2);
                break;
            }
        }
        i+=1;
        if ((j>(int)ebanks[class2].size()-1) && (i>(int)ebanks[class1].size()-1))
        {
            break;
        }
        //printf("end of i loop %d of %d\n",i,ebanks[class1].size());
    }
    */
    
    // //NEW FOR LOOP OVER EVENTS(liquids)
    for (int queue=0;queue<=ebanks[class2].size();queue++)
      {
	//////////////////////JAIMES COMPARE TEST//////////////////////
	Int_t mydid=ebanks[class2][queue].detector_id;
	Double_t current_time=ebanks[class2][queue].time;
	//if(mydid==5)cout<<current_time<<endl;
	Double_t current_time_wrapped=fmod(current_time,3*coincwin_hi);//usually 2, gonna try 3
	if (ebanks[class2][queue].integral[1]==0) continue;
	Double_t signal = ebanks[class2][queue].integral[1];
	Double_t background = ebanks[class2][queue].integral[0];
	//std::cout<<(signal-background)/signal<<std::endl;
	Double_t psd=(signal-background)/signal;

	if (lasttime[mydid]>0.0 &&(current_time!=lasttime[mydid]))
	  {
	    liqdeadtime[mydid]->Fill(current_time-lasttime[mydid]);
	    //if(mydid==2 && psd>0.21 && psd<0.42)//neutron cut
	    if(mydid==2 && psd<0.2)//gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		//if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector2
		  }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid ) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
		psdliquid2->Fill(psd);
				 
	      }//END OF DET2
	    //else if(mydid==3 && psd>0.15 && psd<0.3)//neutron cut
	    else if(mydid==3 && psd<0.13)//gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector3
		  }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
	      }//END OF DET3
	    //else if(mydid==4 && psd>0.12 && psd<0.3) //neutron cut
	    else if(mydid==4 && psd<0.1) //gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector4
		      }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
	      }//END OF DET4
	    //else if(mydid==5 && psd<0.33 && psd>0.15) //neutron cut
	    else if(mydid==5 && psd<0.13) //gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector5
		  }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid ) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
	      }//END OF DET5
	      //else if(mydid==6 && psd<0.33 && psd>0.16) //neutron cut
	    else if(mydid==6 && psd<0.14) //gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector6
		  }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
	      }//END OF DET6
	      //else if(mydid==7 && psd<0.33 && psd>0.16) //neutron cut
	    else if(mydid==7 && psd<0.14) //gamma cut
	      {
		if (lastdet<0) lastdet=mydid;
		if (ralasttime[mydid]<0.0) ralasttime[mydid]=current_time;//This gets us past the first event
		if (singlelasttime<0.0) singlelasttime=current_time;
		if ( ((current_time-ralasttime[mydid]) <= coincwin_hi) && (ralasttime[mydid]!=current_time))
		  {
		    rossialpha_liquid[mydid]->Fill( (current_time-ralasttime[mydid]) ); //detector7
		  }
		if ( (current_time-singlelasttime <= coincwin_hi) && lastdet!=mydid) shortraliquid->Fill(current_time-singlelasttime);
		liqbackground[mydid]->Fill(current_time_wrapped);
	      }//END OF DET7
	  }//END OF MAKING SURE LAST TIME ISNT THIS TIME
	lasttime[mydid]=ebanks[class2][queue].time;
	if( (current_time-ralasttime[mydid]) > coincwin_hi)  ralasttime[mydid]=ebanks[class2][queue].time; //If the event is greater than the coinc window away, that event becomes the new "start" clock
	if((current_time-singlelasttime)>coincwin_hi) 
	  {
	    singlelasttime=ebanks[class2][queue].time;
	    lastdet=mydid;
	  }
      }//end of new loop over liquids
    
    for (int zee=0;zee<=ebanks[class1].size();zee++)
      {
	Int_t did=ebanks[class1][zee].detector_id;
	Double_t now_time=ebanks[class1][zee].time;
	if (blasttime[did]>0.0 &&(now_time!=blasttime[did]))
	  {
	    bgodeadtime[did]->Fill(now_time-blasttime[did]);
	    if (ebanks[class1][zee].integral[1]>16000 && ebanks[class1][zee].integral[1]<18500)
	      {
		bgobackground[did]->Fill(fmod(now_time,3*coincwin_hi));//Multiply by 3 to maybe work?
	      }
	  }
	blasttime[did]=ebanks[class1][zee].time;
      }



    if (dump1r)
    {
      for (int q=0;q<10;q++)
	{
	  if(q>=2 && q<8)
	    {
	      liqbackground[q]->Write();
	      liqdeadtime[q]->Write();
	      rossialpha_liquid[q]->Write();
	    }
	  if(q<1)
	    {
	      bgodeadtime[q]->Write();
	      bgobackground[q]->Write();
	    }
	}
      psdliquid2->Write();
      shortraliquid->Write();
      otree->Write();
      orfile->Close();
    }
    if (dump1h)
    {
     #include "RootStuff/Hist_S1CHINU_Write.C"
        ohfile->Close();
    }

    return 0;
}


Int_t elsa_root::execute_s0bin(Int_t calibrate_time, Int_t calibrate_energy)
{
    gzFile in;
    ELSA_BANK elsa_event;
    short extra_word;
    int readval = 0;
    bool run = true;
    double prev_time[MAXNDETS];
    for (size_t i=0;i<MAXNDETS;++i)
    {
        prev_time[i] = 0;
    }

    printf("open file %s\n",ipath->Data());
	in=gzopen(ipath->Data(),"rb");

    do
    {
        readval = gzread(in,&elsa_event,sizeof(ELSA_BANK));
        //printf("readval %d %d\n",readval,sizeof(ELSA_BANK));
        if (readval != sizeof(ELSA_BANK))
        {
            run = false;
            break;
        }
        if (extra_option!=0)
        {
            //printf("readval %d %d\n",readval,sizeof(short));
            readval = gzread(in,&extra_word,sizeof(short));
	    if (readval != sizeof(short))
	    {
	      run = false;
	      break;
	    }
        }
        if (0)
        {
            if (prev_time[elsa_event.detector_id] > elsa_event.time)
            {
                printf("time discontinuous channel %d %f %f %f\n",elsa_event.detector_id,prev_time[elsa_event.detector_id],elsa_event.time,prev_time[elsa_event.detector_id]-elsa_event.time);
            }
            prev_time[elsa_event.detector_id] = elsa_event.time;
        }
        //printf("detector id %d\n",elsa_event.detector_id);
        if (calibrate_time)
        {
            elsa_event.time = tslope[elsa_event.detector_id]*elsa_event.time + toffset[elsa_event.detector_id];
        }
        if (calibrate_energy)
        {
            elsa_event.integral[0] = eslope[elsa_event.detector_id]*elsa_event.integral[0] + eoffset[elsa_event.detector_id];
            elsa_event.integral[1] = eslope[elsa_event.detector_id]*elsa_event.integral[1] + eoffset[elsa_event.detector_id];
        }
        int type = detclass[elsa_event.detector_id];
        if (type > -1)
        {
            ebanks[type].push_back(elsa_event);
            if (extra_option != 0)
            {
                extras[type].push_back(extra_word);
                if (extra_option==2)
                {
                    elsa_event.integral[0] = extra_word;
                }
            }
        }
    } while (run);
    printf("number events type 0 %d\n",(int)ebanks[0].size());
    printf("number events type 1 %d\n",(int)ebanks[1].size());
    printf("number events type 2 %d\n",(int)ebanks[2].size());
    return 0;
}

Int_t elsa_root::execute_uac_qils(double interp_slope)
{
    EventHeader_t head;	// Midas EventHeader
    EventHeader_t endrun; // end run buffer
    BankHeader_t bhead;		// Midas bank header
    Bank32_t bank32;		// Midas 32bit bank
    QILS_BANK *evinfo = new QILS_BANK();
    QILS_BANK *evinfo_proc = new QILS_BANK(); // caen event info
    ELSA_BANK *elsa_event = new ELSA_BANK();
    test_struct_qils *evaggr = new test_struct_qils();
    FILE *obfile = fopen(opath->Data(),"w");

    UInt_t last_ts[MAXNDETS];
    ULong64_t ts_base[MAXNDETS];
    ULong64_t ts_full = 0;
    for (int i=1;i<=MAXNDETS;++i)
    {
        last_ts[i]=0;
        ts_base[i]=0;
    }
    

    int TotalDataSize;
    int TotalBankSize; //=head.fDataSize;
    int EventBankSize; //=head.fDataSize;
    short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
    short waveform[20000];
    int nevt = 0;
    int n_read_evts = 0;
    int neweventId=0;
    long devt_padding = 0; // padding for the devt struct read
    Int_t t0_count = 0;
    Double_t t0_lasttime = -2e6;
    gzFile in;

	in=gzopen(ipath->Data(),"rb");
    bool run = true;
    do
    {
        if (nevt%100000==0)
        {
            printf("on event %d\n",nevt);
        }
        gzread(in,&head,sizeof(EventHeader_t));
        TotalDataSize=head.fDataSize;
        if(MidasEventPrint && nevt >= MidasEventPrintThresh){
            cout << "Event_HEADER " << endl;
            cout << hex << head.fEventId << endl;
            cout << dec << head.fTriggerMask << endl;
            cout << dec << head.fSerialNumber << endl;
            cout << dec << head.fTimeStamp << endl;
            cout << dec << head.fDataSize << endl;
        }	
        if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ){
            if(head.fEventId==0x8001)
            {
                endrun = head;
                break;
            }
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        else if(head.fEventId==1){
            //printf("This is event data\n");
            // this is event data
            gzread(in,&bhead,sizeof(BankHeader_t));	
            if(MidasEventPrint && nevt > MidasEventPrintThresh){
                cout << "Bank_HEADER " << endl;
                cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
                cout << dec << bhead.fFlags << endl;
            }

            TotalBankSize = bhead.fDataSize;
            int insidecounter = 0;
            while(TotalBankSize>0){
                insidecounter += 1;
                            
                gzread(in,&bank32,sizeof(Bank32_t));
                TotalBankSize-=sizeof(Bank32_t);
                if(MidasEventPrint && nevt > MidasEventPrintThresh){
                    cout << "BANK " << endl;
                    cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                    cout << dec << bank32.fType << endl;
                    cout << dec << bank32.fDataSize << endl;
                }
                EventBankSize = bank32.fDataSize;
                if (bank32.fName[0]=='Q' && bank32.fName[1]=='I')
                {
                    // TODO: this is a COMPLETE disaster but just happens to work. Fix it.
                    evaggr->N = 0; // reset how many events we've processed this event
                    int number_pils_events = bank32.fDataSize/sizeof(QILS_BANK);
                    //printf("reading %d events\n",number_pils_events);
                    for (int eye = 0; eye < number_pils_events; ++eye)
                    {
                        gzread(in,evinfo,sizeof(QILS_BANK));
                        //printf("the size of a pils is %d\n",sizeof(QILS_BANK));
                        gzseek(in,devt_padding,SEEK_CUR);
                        //printf("channel number %i\n",evinfo->channel);
                        TotalBankSize-=sizeof(QILS_BANK)+devt_padding;
                        EventBankSize-=sizeof(QILS_BANK)+devt_padding;
                        evaggr->P[evaggr->N] = *evinfo;
                        evaggr->N += 1;
                        n_read_evts += 1;
                    }

                                    
                    // snag the trig bank
                    gzread(in,&bank32,sizeof(Bank32_t));
                    TotalBankSize-=sizeof(Bank32_t);
                    EventBankSize = bank32.fDataSize;

                    char *fData;
                    fData=(char*)malloc(bank32.fDataSize);
                    //printf("data size %d\n",bank32.fDataSize);
                    gzread(in,fData,bank32.fDataSize);
                    TotalBankSize -= EventBankSize;
                    free (fData);	
                                
                    // begin funny place between peaks and cpu
                    while (true)
                    {
                        // the peaks bank should be here
                        gzread(in,&bank32,sizeof(Bank32_t));
                        TotalBankSize-=sizeof(Bank32_t);
                        EventBankSize = bank32.fDataSize;
                        if(bank32.fName[0]=='p')
                        {
                            int whichpeak = atoi(&bank32.fName[1]);
                            //printf("about to read peaks size %d\n",bank32.fDataSize);
                            gzread(in,imported_peaks[whichpeak],bank32.fDataSize);
                            //gzread(in,waveform,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                        } else
                        {
                            // get the cpu bank information
                            fData=(char*)malloc(bank32.fDataSize);
                            gzread(in,fData,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                            free (fData);	
                            break; // you break here because the cpu comes last
                        }

                    }
                    //printf("total bank left now is %i\n",TotalBankSize);
                    // and here is where we actually tie thangs back together
                    int last_detnum = evaggr->P[0].detector_id;
                    int where_in_peakbank = 0;
                    //printf("about to stitch together %d events\n",evaggr->N);
                    for (size_t evtnum=0;evtnum<evaggr->N;++evtnum)
                    {
                        if (evaggr->N > MaxHitsPerT0)
                        {
                          printf("on event %d of %d\n",evtnum,evaggr->N);
                          break;
                          //return 0;
                        }
                        INT current_detnum = evaggr->P[evtnum].detector_id;
                        if (current_detnum != last_detnum)
                        {
                            where_in_peakbank = 0;
                        }
                        else
                        {
                        }
                        uint32_t wflen = evaggr->P[evtnum].wavelet_stop - evaggr->P[evtnum].wavelet_start;
                        for (size_t wfindex=where_in_peakbank;wfindex<where_in_peakbank+wflen;++wfindex)
                        {
                            //printf("on index %d\n",wfindex);
                            evaggr->wavelets[evtnum][wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                            waveform[wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                        }
                        where_in_peakbank += wflen;
                        last_detnum = current_detnum;

                        // finally write to output
                        *evinfo_proc = evaggr->P[evtnum];
                        if (dump0b)
                        {
                            elsa_event->integral[0] = evinfo_proc->integral[1];
                            elsa_event->integral[1] = evinfo_proc->integral[3];
                            elsa_event->detector_id = evinfo_proc->detector_id;
                            if (ftoption == 0)
                            {
                                if (upconvert)
                                {
                                   elsa_event->time = (double)ts_full;
                                } else
                                {
                                   elsa_event->time = (double)evinfo_proc->position;
                                }
                            }
                            else if (ftoption == 1)
                            {
                              if (upconvert)
                              {
                                elsa_event->time = (double)ts_full;
                              } else
                              {
                                elsa_event->time = (double)evinfo_proc->position;
                              }
                              elsa_event->time += interp_slope*((double)evinfo_proc->interpolation);
                              //printf("fine time %f\n",evinfo_proc->interpolation);
                            } else if (ftoption == 2)
                            {
                                float ftime = digital_cfd(evaggr->wavelets[evtnum],evaggr->filtered_wavelets[evtnum],wflen,evinfo_proc->detector_id);
                                if (upconvert)
                                {
                                  elsa_event->time = (double)ts_full;
                                } else
                                {
                                  elsa_event->time = (double)evinfo_proc->position;
                                }
                                elsa_event->time += (double)ftime;
                            }
                            //printf("time in seconds %f\n",elsa_event->time/0.5e9);
                            // check for t0 downscaling situation
                            if (t0ds > 0 && (detclass[elsa_event->detector_id]==0))
                            {
                                if ((elsa_event->time - t0_lasttime)>t0to)
                                {
                                    t0_lasttime = elsa_event->time;
                                    t0_count = 0;
                                } else
                                {
                                    if ((t0_count+1) < t0ds)
                                    {
                                        t0_count += 1;
                                        continue;
                                    } else
                                    {
                                        t0_count = 0;
                                    }
                                }
                            }
                            //printf("writing to file\n");
                            if (elsa_event->integral[1]>threshold[elsa_event->detector_id]&&elsa_event->integral[1]<saturate[elsa_event->detector_id])
                            {
                              fwrite(elsa_event,sizeof(ELSA_BANK),1,obfile);
                                if (extra_option == 1)
                                {
                                    short ph = find_peakheight(evaggr->wavelets[evtnum],wflen,evinfo_proc->detector_id);
                                    fwrite(&ph,sizeof(short),1,obfile);
                                }
                            }
                        }

                    }
                }
                break;
            }
                
                    //printf("BROKEN OUT OF EVENT READ LOOP\n");
        }
        else {
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        nevt += 1;
        if (nevt>STOPatEVENT)
        {
            printf("stopped at event %d of %d\n",nevt,STOPatEVENT);
            run = false;
        }
    } while (run);

    fclose(obfile);
    delete evinfo;
    delete evinfo_proc;
    delete elsa_event;
    delete evaggr;
    return 0;
}


Int_t elsa_root::execute_uac_pils(double interp_slope)
{
    int readval = 0;
    EventHeader_t head;	// Midas EventHeader
    EventHeader_t endrun; // end run buffer
    BankHeader_t bhead;		// Midas bank header
    Bank32_t bank32;		// Midas 32bit bank
    PILS_BANK *evinfo = new PILS_BANK();
    PILS_BANK *evinfo_proc = new PILS_BANK(); // caen event info
    ELSA_BANK *elsa_event = new ELSA_BANK();
    test_struct_pils *evaggr = new test_struct_pils();
    FILE *obfile = fopen(opath->Data(),"w");

    UInt_t last_ts[MAXNDETS];
    ULong64_t ts_base[MAXNDETS];
    ULong64_t ts_full = 0;
    for (int i=1;i<=MAXNDETS;++i)
    {
        last_ts[i]=0;
        ts_base[i]=0;
    }
    

    int TotalDataSize;
    int TotalBankSize; //=head.fDataSize;
    int EventBankSize; //=head.fDataSize;
    short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
    short waveform[20000];
    int nevt = 0;
    int n_read_evts = 0;
    int neweventId=0;
    long devt_padding = 0; // padding for the devt struct read
    Int_t t0_count = 0;
    Double_t t0_lasttime = -2e6;
    gzFile in;

	in=gzopen(ipath->Data(),"rb");
    bool run = true;
    do
    {
        if (nevt%100000==0)
        {
            printf("on event %d\n",nevt);
        }
        gzread(in,&head,sizeof(EventHeader_t));
        TotalDataSize=head.fDataSize;
        if(MidasEventPrint && nevt >= MidasEventPrintThresh){
            cout << "Event_HEADER " << endl;
            cout << hex << head.fEventId << endl;
            cout << dec << head.fTriggerMask << endl;
            cout << dec << head.fSerialNumber << endl;
            cout << dec << head.fTimeStamp << endl;
            cout << dec << head.fDataSize << endl;
        }	
        if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ){
            if(head.fEventId==0x8001)
            {
                endrun = head;
                break;
            }
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        else if(head.fEventId==1){
            //printf("This is event data\n");
            // this is event data
            gzread(in,&bhead,sizeof(BankHeader_t));	
            if(MidasEventPrint && nevt > MidasEventPrintThresh){
                cout << "Bank_HEADER " << endl;
                cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
                cout << dec << bhead.fFlags << endl;
            }

            TotalBankSize = bhead.fDataSize;
            int insidecounter = 0;
            while(TotalBankSize>0){
                insidecounter += 1;
                            
                gzread(in,&bank32,sizeof(Bank32_t));
                TotalBankSize-=sizeof(Bank32_t);
                if(MidasEventPrint && nevt > MidasEventPrintThresh){
                    cout << "BANK " << endl;
                    cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                    cout << dec << bank32.fType << endl;
                    cout << dec << bank32.fDataSize << endl;
                }
                EventBankSize = bank32.fDataSize;
                if (bank32.fName[0]=='P' && bank32.fName[1]=='I')
                {
                    // TODO: this is a COMPLETE disaster but just happens to work. Fix it.
                    evaggr->N = 0; // reset how many events we've processed this event
                    int number_pils_events = bank32.fDataSize/sizeof(PILS_BANK);
                    //printf("reading %d events\n",number_pils_events);
                    for (int eye = 0; eye < number_pils_events; ++eye)
                    {
                        gzread(in,evinfo,sizeof(PILS_BANK));
                        //printf("the size of a pils is %d\n",sizeof(QILS_BANK));
                        gzseek(in,devt_padding,SEEK_CUR);
                        //printf("channel number %i\n",evinfo->channel);
                        TotalBankSize-=sizeof(PILS_BANK)+devt_padding;
                        EventBankSize-=sizeof(PILS_BANK)+devt_padding;
                        evaggr->P[evaggr->N] = *evinfo;
                        evaggr->N += 1;
                        n_read_evts += 1;
                    }

                                    
                    // snag the trig bank
                    gzread(in,&bank32,sizeof(Bank32_t));
                    TotalBankSize-=sizeof(Bank32_t);
                    EventBankSize = bank32.fDataSize;

                    char *fData;
                    fData=(char*)malloc(bank32.fDataSize);
                    //printf("data size %d\n",bank32.fDataSize);
                    gzread(in,fData,bank32.fDataSize);
                    TotalBankSize -= EventBankSize;
                    free (fData);	
                                
                    // begin funny place between peaks and cpu
                    while (true)
                    {
                        // the peaks bank should be here
                        readval = gzread(in,&bank32,sizeof(Bank32_t));
			if (readval != sizeof(Bank32_t))
			{
                            printf("Found end of file in weird place... did DAQ crash?\n");
			    run = false;
			    break;
			}
                        TotalBankSize-=sizeof(Bank32_t);
                        EventBankSize = bank32.fDataSize;
                        if(MidasEventPrint && nevt > MidasEventPrintThresh){
                          printf("Looking for cpu / peak banks\n");
                          cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                        }
                        if(bank32.fName[0]=='p')
                        {
                            int whichpeak = atoi(&bank32.fName[1]);
                            //printf("about to read peaks size %d\n",bank32.fDataSize);
                            gzread(in,imported_peaks[whichpeak],bank32.fDataSize);
                            //gzread(in,waveform,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                        } else
                        {
                            // get the cpu bank information
                            fData=(char*)malloc(bank32.fDataSize);
                            gzread(in,fData,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                            free (fData);	
                            break; // you break here because the cpu comes last
                        }

                    }
                    if (run == false)
                    {
                      break;
                    }
                    //printf("total bank left now is %i\n",TotalBankSize);
                    // and here is where we actually tie thangs back together
                    int last_detnum = evaggr->P[0].detector_id;
                    int where_in_peakbank = 0;
                    //printf("about to stitch together %d events\n",evaggr->N);
                    for (size_t evtnum=0;evtnum<evaggr->N;++evtnum)
                    {
                        if (evaggr->N > MaxHitsPerT0)
                        {
                          printf("on event %d of %d\n",evtnum,evaggr->N);
                          break;
                          //return 0;
                        }
                        INT current_detnum = evaggr->P[evtnum].detector_id;
                        if (current_detnum != last_detnum)
                        {
                            where_in_peakbank = 0;
                        }
                        else
                        {
                        }
                        uint32_t wflen = evaggr->P[evtnum].wavelet_stop - evaggr->P[evtnum].wavelet_start;
                        for (size_t wfindex=where_in_peakbank;wfindex<where_in_peakbank+wflen;++wfindex)
                        {
                            //printf("on index %d\n",wfindex);
                            evaggr->wavelets[evtnum][wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                            waveform[wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                        }
                        where_in_peakbank += wflen;
                        last_detnum = current_detnum;
                        //printf("made it here\n");

                        // finally write to output
                        *evinfo_proc = evaggr->P[evtnum];
                        if (upconvert) // whether to look for clock rollovers
                        {
                            if (evaggr->P[evtnum].position<last_ts[current_detnum])
                            {
                                //printf("hey i found rollover on channel %i\n",current_detnum);
                                ts_base[current_detnum] =ts_base[current_detnum] + pow(2,NBITSCLOCK);
                                //std::cout << "My ts base is now " << ts_base[current_detnum] << std::endl;
                            
                            }
                            ts_full = ts_base[current_detnum] + evaggr->P[evtnum].position;
                            last_ts[current_detnum] = evaggr->P[evtnum].position;
                        }
                        if (dump0b)
                        {
                            if(MidasEventPrint && nevt > MidasEventPrintThresh){
                              printf("in the binary dump stage\n");
                            }
                            elsa_event->integral[0] = evinfo_proc->integral[1];
                            elsa_event->integral[1] = evinfo_proc->integral[3];
                            elsa_event->detector_id = evinfo_proc->detector_id;
                            if (ftoption == 0)
                            {
                                if (upconvert)
                                {
                                   elsa_event->time = (double)ts_full;
                                } else
                                {
                                   elsa_event->time = (double)evinfo_proc->position;
                                }
                            }
                            else if (ftoption == 1)
                            {
                              if (upconvert)
                              {
                                elsa_event->time = (double)ts_full;
                              } else
                              {
                                elsa_event->time = (double)evinfo_proc->position;
                              }
                              elsa_event->time += interp_slope*((double)evinfo_proc->interpolation);
                              //printf("fine time %f\n",evinfo_proc->interpolation);
                            } else if (ftoption == 2)
                            {
                                float ftime = digital_cfd(evaggr->wavelets[evtnum],evaggr->filtered_wavelets[evtnum],wflen,evinfo_proc->detector_id);
                                if (upconvert)
                                {
                                  elsa_event->time = (double)ts_full;
                                } else
                                {
                                  elsa_event->time = (double)evinfo_proc->position;
                                }
                                elsa_event->time += (double)ftime;
                            }
                            //printf("time in seconds %f\n",elsa_event->time/0.5e9);
                            // check for t0 downscaling situation
                            if (t0ds > 0 && (detclass[elsa_event->detector_id]==0))
                            {
                                if ((elsa_event->time - t0_lasttime)>t0to)
                                {
                                    t0_lasttime = elsa_event->time;
                                    t0_count = 0;
                                } else
                                {
                                    if ((t0_count+1) < t0ds)
                                    {
                                        t0_count += 1;
                                        continue;
                                    } else
                                    {
                                        t0_count = 0;
                                    }
                                }
                            }
                            //printf("writing to file\n");
                            fwrite(elsa_event,sizeof(ELSA_BANK),1,obfile);
                            if (extra_option == 1)
                            {
                                short ph = find_peakheight(evaggr->wavelets[evtnum],wflen,evinfo_proc->detector_id);
                                fwrite(&ph,sizeof(short),1,obfile);
                            }
                        }

                    }
                }
                if(MidasEventPrint && nevt > MidasEventPrintThresh){
                  printf("done with PILS bank\n");
                }
                break;
            }
                
                    //printf("BROKEN OUT OF EVENT READ LOOP\n");
        }
        else {
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        nevt += 1;
        if (nevt>STOPatEVENT)
        {
            printf("stopped at event %d of %d\n",nevt,STOPatEVENT);
            run = false;
        }
    } while (run);

    fclose(obfile);
    delete evinfo;
    delete evinfo_proc;
    delete elsa_event;
    delete evaggr;
    return 0;
}

Int_t elsa_root::execute_anna()
{
    return 0;
}

Int_t elsa_root::execute_uac_cevt(double interp_slope)
{
    EventHeader_t head;	// Midas EventHeader
    EventHeader_t endrun; // end run buffer
    BankHeader_t bhead;		// Midas bank header
    Bank32_t bank32;		// Midas 32bit bank
    CEVT_BANK *evinfo = new CEVT_BANK();
    CEVT_BANK *evinfo_proc = new CEVT_BANK(); // caen event info
    ELSA_BANK *elsa_event = new ELSA_BANK();
    test_struct_cevt *evaggr = new test_struct_cevt();
    FILE *obfile;
    // Histogram stuff
    TFile *ohfile;
    //#include "RootStuff/Hist_S0DANCE_Define.C"

    if (dump0b)
    {
        obfile = fopen(opath->Data(),"w");
    }
    if (dump0h)
    {
        printf("opening root histo file %s\n",ohpath->Data());
        ohfile = new TFile(ohpath->Data(),"recreate");
        //#include "RootStuff/Hist_S0DANCE_Create.C"
    }
    printf("start cevt read\n");
    

    int TotalDataSize;
    int TotalBankSize; //=head.fDataSize;
    int EventBankSize; //=head.fDataSize;
    short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
    short waveform[20000];
    int nevt = 0;
    int n_read_evts = 0;
    int neweventId=0;
    long devt_padding = 0; // padding for the devt struct read
    Int_t t0_count = 0;
    Double_t t0_lasttime = -2e6;
    gzFile in;
    Int_t detector_count[MAXNDETS];

	in=gzopen(ipath->Data(),"rb");
    printf("opening file %s\n",ipath->Data());
    bool run = true;
    do
    {
        gzread(in,&head,sizeof(EventHeader_t));
        TotalDataSize=head.fDataSize;
        if(MidasEventPrint && nevt >= MidasEventPrintThresh){
            cout << "Event_HEADER " << endl;
            cout << hex << head.fEventId << endl;
            cout << dec << head.fTriggerMask << endl;
            cout << dec << head.fSerialNumber << endl;
            cout << dec << head.fTimeStamp << endl;
            cout << dec << head.fDataSize << endl;
        }	
        if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ){
            if(head.fEventId==0x8001)
            {
                endrun = head;
                break;
            }
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        else if(head.fEventId==1){
            //printf("This is event data\n");
            // this is event data
            gzread(in,&bhead,sizeof(BankHeader_t));	
            if(MidasEventPrint && nevt > MidasEventPrintThresh){
                cout << "Bank_HEADER " << endl;
                cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
                cout << dec << bhead.fFlags << endl;
            }

            TotalBankSize = bhead.fDataSize;
            int insidecounter = 0;
            while(TotalBankSize>0){
                insidecounter += 1;
                            
                gzread(in,&bank32,sizeof(Bank32_t));
                TotalBankSize-=sizeof(Bank32_t);
                if(MidasEventPrint && nevt > MidasEventPrintThresh){
                    cout << "BANK " << endl;
                    cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                    cout << dec << bank32.fType << endl;
                    cout << dec << bank32.fDataSize << endl;
                }
                EventBankSize = bank32.fDataSize;
                if (bank32.fName[0]=='C' && bank32.fName[1]=='E')
                {
                    // TODO: this is a COMPLETE disaster but just happens to work. Fix it.
                    evaggr->N = 0; // reset how many events we've processed this event
                    int number_pils_events = bank32.fDataSize/sizeof(CEVT_BANK);
                    for (int eye = 0; eye < number_pils_events; ++eye)
                    {
                        gzread(in,evinfo,sizeof(CEVT_BANK));
                        //printf("the size of a pils is %d\n",sizeof(QILS_BANK));
                        gzseek(in,devt_padding,SEEK_CUR);
                        //printf("channel number %i\n",evinfo->channel);
                        TotalBankSize-=sizeof(CEVT_BANK)+devt_padding;
                        EventBankSize-=sizeof(CEVT_BANK)+devt_padding;
                        evaggr->P[evaggr->N] = *evinfo;
                        evaggr->N += 1;
                        n_read_evts += 1;
                    }
                    if (MidasEventPrint)
                    {
                      printf("read the firmware events\n");
                    }

                                    
                    // snag the trig bank
                    gzread(in,&bank32,sizeof(Bank32_t));
                    TotalBankSize-=sizeof(Bank32_t);
                    EventBankSize = bank32.fDataSize;

                    char *fData;
                    fData=(char*)malloc(bank32.fDataSize);
                    gzread(in,fData,bank32.fDataSize);
                    TotalBankSize -= EventBankSize;
                    free (fData);	
                                
                    // begin funny place between peaks and cpu
                    while (true)
                    {
                        // the peaks bank should be here
                        gzread(in,&bank32,sizeof(Bank32_t));
                        TotalBankSize-=sizeof(Bank32_t);
                        EventBankSize = bank32.fDataSize;
                        if(bank32.fName[0]=='p')
                        {
                            int whichpeak = atoi(&bank32.fName[1]);
                            gzread(in,imported_peaks[whichpeak],bank32.fDataSize);
                            //gzread(in,waveform,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                        } else
                        {
                            // get the cpu bank information
                            fData=(char*)malloc(bank32.fDataSize);
                            gzread(in,fData,bank32.fDataSize);
                            TotalBankSize -= EventBankSize;
                            free (fData);	
                            break; // you break here because the cpu comes last
                        }

                    }
                    if (MidasEventPrint)
                    {
                      printf("total bank left now is %i\n",TotalBankSize);
                    }
                    // and here is where we actually tie thangs back together
                    int last_detnum = evaggr->P[0].detector_id;
                    int where_in_peakbank = 0;
                    for (size_t evtnum=0;evtnum<evaggr->N;++evtnum)
                    {
                        if (MidasEventPrint)
                        {
                            
                            if (evaggr->N > MaxHitsPerT0)
                            {
                              printf("on event %d of %d\n\n\n\n\n\n\n",evtnum,evaggr->N);
                              break;
                              //return 0;
                            }
                           
                        }
                        INT current_detnum = evaggr->P[evtnum].detector_id;
                        if (current_detnum != last_detnum)
                        {
                            where_in_peakbank = 0;
                        }
                        else
                        {
                        }
                        uint32_t wflen = evaggr->P[evtnum].width;
                        for (size_t wfindex=where_in_peakbank;wfindex<where_in_peakbank+wflen;++wfindex)
                        {
                            evaggr->wavelets[evtnum][wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                            waveform[wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                        }
                        where_in_peakbank += wflen;
                        last_detnum = current_detnum;

                        // finally write to output
                        *evinfo_proc = evaggr->P[evtnum];
                        if (1)
                        {
                            elsa_event->integral[0] = evinfo_proc->integral[0];
                            if (lsint_sub)
                            {
                                elsa_event->integral[1] = evinfo_proc->integral[1] - evinfo_proc->integral[0];
                            } else
                            {
                                elsa_event->integral[1] = evinfo_proc->integral[1];
                            }
                            elsa_event->detector_id = evinfo_proc->detector_id + detid_shift;
                            /*
                            if (lsgate_ds[elsa_event->detector_id]!=0) // lsgate related downscaling option
                            {
                              if (lsgate_ds[elsa_event->detector_id]==-1) // kill everything passing the gate
                              {
                                if (lsgate.IsInside(elsa_event->integral[1],elsa_event->integral[0]))
                                {
                                  break;
                                }
                              } else
                              {
                                if (lsgate.IsInside(elsa_event->integral[1],elsa_event->integral[0]))
                                {
                                  detector_count[elsa_event->detector_id] += 1;
                                  if (detector_count[elsa_event->detector_id] < lsgate_ds[elsa_event->detector_id])
                                  {
                                    break;
                                  } else
                                  {
                                    detector_count[elsa_event->detector_id] = 0;
                                  }
                                }
                              }
                            }
                            */
                            if (ftoption == 0)
                            {
                              //elsa_event->time = evinfo_proc->position;
                              elsa_event->time = (double)(evinfo_proc->position & 0x7FFFFFFFFFFF);
                            }
                            else if (ftoption == 1)
                            {
                              //elsa_event->time = (double)evinfo_proc->position + interp_slope*((double)evinfo_proc->interpolation);
                              elsa_event->time = (double)(evinfo_proc->position & 0x7FFFFFFFFFFF);
                              //printf("fine time %f\n",evinfo_proc->interpolation);
                            }
                            // check for t0 downscaling situation
                            if (t0ds > 0 && (detclass[elsa_event->detector_id]==0))
                            {
                                if ((elsa_event->time - t0_lasttime)>t0to)
                                {
                                    t0_lasttime = elsa_event->time;
                                    t0_count = 0;
                                } else
                                {
                                    if ((t0_count+1) < t0ds)
                                    {
                                        t0_count += 1;
                                        continue;
                                    } else
                                    {
                                        t0_count = 0;
                                    }
                                }
                            }
                            if (dump0b)
                            {
                                fwrite(elsa_event,sizeof(ELSA_BANK),1,obfile);
                            }
                            if (dump0h)
                            {
                                //#include "RootStuff/Hist_S0DANCE_Fill.C"
                            }
                        }

                    }
                }
                break;
            }
                
                    //printf("BROKEN OUT OF EVENT READ LOOP\n");
        }
        else {
            char *fData;
            fData=(char*)malloc(head.fDataSize);
            gzread(in,fData,head.fDataSize);
            if(MidasEventPrint && neweventId > MidasEventPrintThresh){
                for(size_t i=0;i<head.fDataSize;i++){
                    cout << fData[i];			
                }
            }	
            free (fData);	
        }
        nevt += 1;
        if (nevt>STOPatEVENT)
        {
            printf("stopped at event %d of %d",nevt,STOPatEVENT);
            run = false;
        }
    } while (run);

    if (dump0b)
    {
        fclose(obfile);
    }
    if (dump0h)
    {
        //#include "RootStuff/Hist_S0DANCE_Write.C"
      ohfile->Close();
    }
    return 0;
}

Int_t elsa_root::execute_devt(Int_t dtype)
{
    /*
    for (int i=0;i<32;++i)
    {
        for (int j=0;j<4;++j)
        {
            printf("%i %f %f %f %f\n",i,time_filter_params[i][j],time_filter_params[i][j],time_filter_params[i][j],time_filter_params[i][j]);
        }
    }

    return 0;
    */
    gzFile in;
    DEVT_BANK devt_event;
    ELSA_BANK elsa_event;
    short wf[34000];
    short wf2[34000];
    float fwf[34000];
    Byte_t dp0[34000];
    Byte_t dp1[34000];
    int readval = 0;
    bool run = true;
    float ftime;
    double prev_time[MAXNDETS];
    for (size_t i=0;i<MAXNDETS;++i)
    {
        prev_time[i] = 0;
    }
    // Histogram stuff
    TFile *ohfile;

    printf("open file %s\n",ipath->Data());
	in=gzopen(ipath->Data(),"rb");

    FILE *obfile = fopen(opath->Data(),"w");

    if (dump0h)
    {
        printf("opening root histo file %s\n",ohpath->Data());
        ohfile = new TFile(ohpath->Data(),"recreate");
        #include "RootStuff/Hist_S0LENZ_Create.C"
    }

    int counter = 0;
    do
    {
        if (counter>STOPatEVENT)
        {
            printf("stopping at event %i\n",counter);
            break;
        }
        counter += 1;
        ftime = 0;
        readval = gzread(in,&devt_event,sizeof(DEVT_BANK));
        if (dtype == 7 || dtype ==8)
        {
          gzread(in,wf,sizeof(UShort_t)*devt_event.Ns);
          if (dtype==8)
          {
            gzread(in,wf2,sizeof(UShort_t)*devt_event.Ns);
            gzread(in,dp0,sizeof(Byte_t)*devt_event.Ns);
            gzread(in,dp1,sizeof(Byte_t)*devt_event.Ns);
          }
        }
        //printf("readval %d %d\n",readval,sizeof(ELSA_BANK));
        if (readval != sizeof(DEVT_BANK))
        {
            run = false;
            break;
        }
        elsa_event.integral[0] = devt_event.sgate;
        elsa_event.integral[1] = devt_event.lgate;
        elsa_event.detector_id = devt_event.board*MAXNCHANNELS + devt_event.channel;
        if (ftoption==1)
        {
          elsa_event.time = (double)devt_event.timestamp + ((double)devt_event.baseline)/1023.;
        } else if (ftoption==3)
        {
            ftime = digital_double_derivative(wf,fwf,devt_event.Ns,elsa_event.detector_id);
            if (extra_option==2)
            {
                //printf("wf len %i \n",devt_event.Ns);
                elsa_event.integral[1] = integral_peakheight(wf, devt_event.Ns, elsa_event.detector_id, (short)ftime);
                /*
                if (elsa_event.detector_id>16)
                {
                    printf("detectof %i ftime %f ph %f\n",elsa_event.detector_id,ftime,elsa_event.integral[1]);
                }
                */
            }
            elsa_event.time = (double)devt_event.timestamp + (double)ftime;
        } else
        {
          elsa_event.time = (double)devt_event.timestamp;
        }
        fwrite(&elsa_event,sizeof(ELSA_BANK),1,obfile);
        if (ftoption==3 && extra_option==2)
        {
            UShort_t position = (UShort_t)(65536*(ftime/(float)devt_event.Ns));
            fwrite(&position,sizeof(UShort_t),1,obfile);
        }
        if (dump0h)
        {
            #include "RootStuff/Hist_S0LENZ_Fill.C"
        }
    } while (run);
    fclose(obfile);
    if (dump0h)
    {
        #include "RootStuff/Hist_S0LENZ_Write.C"
        ohfile->Close();
    }
    return 0;
}
