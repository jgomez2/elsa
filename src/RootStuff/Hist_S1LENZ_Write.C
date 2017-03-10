printf("writing stage 1 histos\n");

TDirectoryFile *basics=new TDirectoryFile("basics","basics");	
ohfile->Add(basics);
basics->cd();

h_tof->Write();
h_tof->Delete();

h_tof_wrapped->Write();
h_tof_wrapped->Delete();

h_ph->Write();
h_ph->Delete();

h_pos_tof->Write();
h_pos_tof->Delete();

h_ph_pos->Write();
h_ph_pos->Delete();

h_ph_tof->Write();
h_ph_tof->Delete();

printf("done writing histos\n");
