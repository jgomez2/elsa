printf("writing standard histos\n");
TDirectoryFile *basic_dir=new TDirectoryFile("basic","basic");	
ohfile->Add(basic_dir);
basic_dir->cd();

h_c1_c2->Write();
h_c1_c2->Delete();

h_c1tof_b0gated->Write();
h_c1tof_b0gated->Delete();

h_c2tof_b0gated->Write();
h_c2tof_b0gated->Delete();

h_c2tof_b0gated_psdl->Write();
h_c2tof_b0gated_psdl->Delete();

h_c2tof_b0gated_psdg->Write();
h_c2tof_b0gated_psdg->Delete();

h_ci_id->Write();
h_ci_id->Delete();

h_id_t0tof->Write();
h_id_t0tof->Delete();

h_id_t0tof_gam->Write();
h_id_t0tof_gam->Delete();

h_ppci_t0tof->Write();
h_ppci_t0tof->Delete();

h_lci_ntof->Write();
h_lci_ntof->Delete();

h_lci_ntof_gated->Write();
h_lci_ntof_gated->Delete();

ohfile->cd();
TDirectoryFile *psd_dir=new TDirectoryFile("psd","psd");	
ohfile->Add(psd_dir);
psd_dir->cd();
for (int eye=0;eye<80;++eye)
{
    h_psd[eye]->Write();
    h_psd[eye]->Delete();
}
h_psd_1d->Write();
h_psd_1d->Delete();
