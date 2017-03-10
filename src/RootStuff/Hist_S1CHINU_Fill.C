// always fill the hit matrix
h_c1_c2->Fill(detid2,detid1);

// calculate a psd parameter
//Double_t psd = (nph-nsg)/nph;

if (detid2<17)
{
    h_c1tof_b0gated->Fill(tof,detid1);
}
if (detid1<17)
{
    h_c2tof_b0gated->Fill(tof,detid2);
    /*
    if (psd > psd1d[detid2])
    {
        h_c2tof_b0gated_psdg->Fill(tof,detid2);
    } else if (psd < psd1d[detid2])
    {
        h_c2tof_b0gated_psdl->Fill(tof,detid2);
    }
    */
}

h_ci_id->Fill(detid1,pph);
h_ci_id->Fill(detid2,nph);

h_id_t0tof->Fill(t0tof_wrapped,detid1);
h_id_t0tof->Fill(t0tof_wrapped+deltat,detid2);
/*
if (psd > psd1d[detid2] && tof>25.0 && tof<65.0)
{
    h_id_t0tof_gam->Fill(t0tof_wrapped,detid1);
    h_id_t0tof_gam->Fill(t0tof_wrapped+deltat,detid2);
}
*/

h_ppci_t0tof->Fill(t0tof_wrapped,pph);
//h_ppci_t0tof->Fill(tof,pph);

h_lci_ntof->Fill(tof,nph);
if ((t0tof_wrapped<1590.&&t0tof_wrapped>200.)&&(detid2!=32&&detid2!=13&&detid2!=11))
{
    h_lci_ntof_gated->Fill(tof,nph);
}


//h_psd[detid2]->Fill(nph,psd);
//h_psd_1d->Fill(detid2,psd);
