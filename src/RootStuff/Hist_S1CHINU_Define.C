// 2d detid/tof
TH2F* h_c1tof_b0gated; // class 1 detids vs. ppac - neutron tof, gated on board 0 neutron detectors
TH2F* h_c2tof_b0gated; // class 2 detids vs. ppac - neutron tof, gated on board 0 ppacs
// special versions of that histogram
TH2F* h_c2tof_b0gated_psdl; // class 2 detids vs. ppac - neutron tof, gated on board 0 ppacs
TH2F* h_c2tof_b0gated_psdg; // class 2 detids vs. ppac - neutron tof, gated on board 0 ppacs

// 2d detector hit matrix plot
TH2F* h_c1_c2;

// charge integral spectra
TH2F* h_ci_id;

// detector id vs. t0 time of flight
TH2F* h_id_t0tof;
TH2F* h_id_t0tof_gam;

// ppac charge integral vs. t0 to ppac timing
TH2F* h_ppci_t0tof;

// ligl charge integral / time of flight
TH2F* h_lci_ntof;
TH2F* h_lci_ntof_gated;

// PSD stuffs
TH2F* h_psd[80];
TH2F* h_psd_1d;
