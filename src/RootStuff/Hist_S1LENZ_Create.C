int howmany = 50;

// basics
h_tof = new TH2F("tof","Unwrapped Time of Flight",10000,0,1e5,howmany,0,howmany);
h_tof_wrapped = new TH2F("tof_wrapped","Wrapped Time of Flight",4000,-200,1800,howmany,0,howmany);
h_ph = new TH2F("ph","Pulse height",howmany,0,howmany,1600,-4000,32000);

// cleanup hists
h_pos_tof = new TH2F("pos_tof","Peak Pos vs. Time of Flight",4000,-200,1800,500,0,66000);
h_ph_pos = new TH2F("ph_pos","Pulse height vs. Peak position",500,0,66000,500,-4000,32000);

// physics
h_ph_tof = new TH2F("ph_tof","Pulse height vs. time of flight",4000,-200,1800,500,-4000,32000);
