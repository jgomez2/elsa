{
    gROOT->ProcessLine(".L macros/loopOnTree.C");
    //TChain *ch = makeChain_pk_DAF_AmBe_LiquidRotation();
    TChain *ch = makeChain_pk_DAF_AmBe_LiquidRotation_LongWindow();
    //ch->AddFriend("friend","./stage1/labr0_cal.root");
    ch->SetAlias("psd","(nph-nsg)/nph");
    /*
    TChain *ch0 = makeChain_pu239_jan2015();
    TChain *ch1 = makeChain_pu239_dec2015();
    */
    /*
    ch->SetAlias("tof","((100000.*(t0tof-0.0)%(100000*1788.819875776))/100000.)");
    ch->SetAlias("tcorr","(tof + pph*0.00283551 - pph*pph*2.48784e-07)");
    ch->SetAlias("ntcorr","(ntof - pph*0.00283551 + pph*pph*2.48784e-07)");
    */
}
