D{
  //set up psd gates for NE213 detectors
  ch->SetAlias("psd","(nph-nsg)/nph");
  TCut ambe="pph>16000 && pph<18500";
  TCut notambe="pph<12000";
  TCut g1="psd<0.2 && detector_id2==2";
  TCut n1="psd>0.21 && psd<0.42 && detector_id2==2";
  TCut g2="psd<0.13 && detector_id2==3";
  TCut n2="psd>0.15 && psd<0.3 && detector_id2==3";
  TCut g3="psd<0.1 && detector_id2==4";
  TCut n3="psd>0.12 && psd<0.3 && detector_id2==4";
  TCut g4="psd<0.13 && detector_id2==5";
  TCut n4="psd>0.15 && psd<0.33 && detector_id2==5";
  TCut g5="psd<0.14 && detector_id2==6";
  TCut n5="psd>0.16 && psd<0.33 && detector_id2==6";
  TCut g6="psd<0.14 && detector_id2==7";
  TCut n6="psd>0.16 && psd<0.33 && detector_id2==7";
  TCut gall=g1 || g2 || g3 || g4 || g5 || g6;
  TCut nall=n1 || n2 || n3 || n4 || n5 || n6;
  // define tof histograms with various gates
  // Modified for AmBe to have extra gate on BGO PH around 4.44 MeV gamma
  ch->Draw("ntof>>tof_gf(30000,-15000,15000)",gall && ambe);
  ch->Draw("ntof>>tof_gb(30000,-15000,15000)",gall && notambe);
  ch->Draw("ntof>>tof_nf(30000,-15000,15000)",nall && ambe);
  ch->Draw("ntof>>tof_nb(30000,-15000,15000)",nall && notambe);
  //Write results to text files
  // TOF for gammas foreground
  TH1F *tof_g_pointer=tof_gf;
  FILE *pFile;
  pFile=fopen("gamma_fgr_TOF_AmBe.txt","w");
  int i;
  int tofBins;
  float time,counts;
  tofBins=tof_g_pointer->GetNbinsX();
  for (i=1;i<=tofBins;++i) {
    time=tof_g_pointer->GetBinCenter(i);
    counts=tof_g_pointer->GetBinContent(i);
    fprintf(pFile,"%e %e \n",time,counts);
  }
  fclose(pFile);
  // TOF for gammas background
  TH1F *tof_g_pointer=tof_gb;
  FILE *pFile;
  pFile=fopen("gamma_bkg_TOF_AmBe.txt","w");
  int i;
  int tofBins;
  float time,counts;
  tofBins=tof_g_pointer->GetNbinsX();
  for (i=1;i<=tofBins;++i) {
    time=tof_g_pointer->GetBinCenter(i);
    counts=tof_g_pointer->GetBinContent(i);
    fprintf(pFile,"%e %e \n",time,counts);
  }
  fclose(pFile);
  // TOF for neutrons foreground
  TH1F *tof_n_pointer=tof_nf;
  FILE *pFile;
  pFile=fopen("neutron_fgr_TOF_AmBe.txt","w");
  tofBins=tof_n_pointer->GetNbinsX();
  for (i=1;i<=tofBins;++i) {
    time=tof_n_pointer->GetBinCenter(i);
    counts=tof_n_pointer->GetBinContent(i);
    fprintf(pFile,"%e %e \n",time,counts);
  }
  fclose(pFile);
}
  // TOF for neutrons background
  TH1F *tof_n_pointer=tof_nb;
  FILE *pFile;
  pFile=fopen("neutron_bkg_TOF_AmBe.txt","w");
  tofBins=tof_n_pointer->GetNbinsX();
  for (i=1;i<=tofBins;++i) {
    time=tof_n_pointer->GetBinCenter(i);
    counts=tof_n_pointer->GetBinContent(i);
    fprintf(pFile,"%e %e \n",time,counts);
  }
  fclose(pFile);
}
