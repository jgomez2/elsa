const int howmany = 48;

TH2F* h_wf[howmany]; // waveform plots for individual detectors
TH2F* h_fwf[howmany]; // filtered waveform plots for individual detectors

TH2F* h_mdev; // maximum deviation of the filtered waveform

TH2F* h_ppos; // position of the peak located on the filtered waveform

TH2F* h_ph; // pulse height for each detector

TH2F* h_walk; // pulse height vs. position in waveform
