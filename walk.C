{
//========= Macro generated from object: walk/Graph
//========= by ROOT version5.26/00
   
   TGraph *graph = new TGraph(8);
   graph->SetName("walk");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(2);
   graph->SetLineWidth(3);
   graph->SetPoint(0,58.3544,30322.9);
   graph->SetPoint(1,58.9713,6703.01);
   graph->SetPoint(2,59.8525,4595.14);
   graph->SetPoint(3,60.3812,3507.93);
   graph->SetPoint(4,61.2624,2864.48);
   graph->SetPoint(5,62.3198,2354.15);
   graph->SetPoint(6,68.5764,445.98);
   graph->SetPoint(7,69.2813,401.604);
   
   TH1F *Graph1 = new TH1F("Graph1","Graph",100,56.7242,69.9422);
   Graph1->SetMinimum(0);
   Graph1->SetMaximum(11054.1);
   Graph1->SetDirectory(0);
   Graph1->SetStats(0);
   Graph1->GetXaxis()->SetLabelFont(42);
   Graph1->GetXaxis()->SetTitleFont(42);
   Graph1->GetYaxis()->SetLabelFont(42);
   Graph1->GetYaxis()->SetTitleFont(42);
   Graph1->GetZaxis()->SetLabelFont(42);
   Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph1);
   
   graph->Draw("");
}
