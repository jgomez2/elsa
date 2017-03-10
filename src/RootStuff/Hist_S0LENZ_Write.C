printf("writing stage 0 histos\n");

TDirectoryFile *wf_dir=new TDirectoryFile("waveforms","waveforms");	
ohfile->Add(wf_dir);
wf_dir->cd();

const int howmany = 48;
for (int eye=0;eye<howmany;++eye)
{
    h_wf[eye]->Write();
    h_wf[eye]->Delete();

    h_fwf[eye]->Write();
    h_fwf[eye]->Delete();
}

h_mdev->Write();
h_mdev->Delete();
h_ppos->Write();
h_ppos->Delete();

h_ph->Write();
h_ph->Delete();

h_walk->Write();
h_walk->Delete();
