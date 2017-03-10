const int howmany = 48;
int npts_wf = 1024;

char buf[1024];
char buf2[1024];

for (int eye=0;eye<howmany;++eye)
{
    sprintf(buf,"wf_%03i",eye);
    h_wf[eye] = new TH2F(buf,buf,npts_wf,0,npts_wf,1600,0,16000);

    sprintf(buf,"fwf_%03i",eye);
    h_fwf[eye] = new TH2F(buf,buf,npts_wf,0,npts_wf,1600,-16000,16000);
}

h_mdev = new TH2F("mdev","Max deviation of filtered wf",howmany,0,howmany,1600,-16000,16000);

h_ppos = new TH2F("ppos","Position of Peak",npts_wf,0,npts_wf,howmany,0,howmany);

h_ph = new TH2F("ph","Pulse height",howmany,0,howmany,1600,-4000,16000);

h_walk = new TH2F("walk","Walk",npts_wf,0,npts_wf,1600,-4000,16000);
