h_ph->Fill(elsa_event.detector_id,elsa_event.integral[1]);
if (elsa_event.detector_id>15)
{
  h_walk->Fill(ftime,elsa_event.integral[1]);
}
