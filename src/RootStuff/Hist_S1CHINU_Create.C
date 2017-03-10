const int howmany = 80;
// 2d detid/tof
h_c1tof_b0gated = new TH2F("c1tof_b0gated","Class 1 tof, board 0 gated",1000,-250,250,howmany,0,howmany);
h_c2tof_b0gated = new TH2F("c2tof_b0gated","Class 2 tof, board 0 gated",1000,-250,250,howmany,0,howmany);
h_c2tof_b0gated_psdl = new TH2F("c2tof_b0gated_psdl","Class 2 tof, board 0 gated psd lesser",500,-50,150,howmany,0,howmany);
h_c2tof_b0gated_psdg = new TH2F("c2tof_b0gated_psdg","Class 2 tof, board 0 gated psd greater",500,-50,150,howmany,0,howmany);

// 2d detector hit matrix'
h_c1_c2 = new TH2F("c1_c2","Hit Matrix",howmany,0,howmany,howmany,0,howmany);

// charge integral spectra
h_ci_id = new TH2F("ci_id","Charge Integral vs. Det ID",howmany,0,howmany,600,0,30000);

// detector id vs. t0 tof
h_id_t0tof = new TH2F("id_t0tof","Det ID vs. beam tof",4000,-100,1900,howmany,0,howmany);
h_id_t0tof_gam = new TH2F("id_t0tof_gam","Det ID vs. beam tof, PFG gated",4000,-100,1900,howmany,0,howmany);

// ppac charge integral vs. t0 to ppac timing
h_ppci_t0tof = new TH2F("ppci_t0tof","PPAC Charge Integral vs. time of flight",600,0,300,300,0,30000);

// ligl charge integral / time of flight
h_lci_ntof = new TH2F("lci_ntof","LiGl Charge Integral vs. Time of Flight",400,-50,150,300,0,15000);
h_lci_ntof_gated = new TH2F("lci_ntof_gated","LiGl Charge Integral vs. Time of Flight (Gated)",400,-50,150,300,0,15000);

h_psd_1d = new TH2F("psd1d","1d psd projections",howmany,0,howmany,125,0,0.25);

// psd
for (int eye=0;eye<howmany;++eye)
{
    char buf[1024];
    char buf2[1024];
    sprintf(buf,"psd_%02i",eye);
    sprintf(buf2,"psd for det %02i",eye);
    h_psd[eye] = new TH2F(buf,buf2,800,0,40000,125,0,0.25);
}
