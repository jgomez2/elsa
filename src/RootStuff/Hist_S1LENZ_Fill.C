h_tof->Fill(tof,detid1);
h_ph->Fill(detid1,lgate);

if (tof<((Double_t)t0ds)*1.3*1800.0)
{
    h_tof_wrapped->Fill(tof_wrapped,detid1);
    h_pos_tof->Fill(tof_wrapped,sgate);
    h_ph_pos->Fill(sgate,lgate);

    // physics
    h_ph_tof->Fill(tof_wrapped,lgate);
}


