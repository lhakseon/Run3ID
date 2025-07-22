{
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0001);

//=========Macro for Trigger Efficiency in MonoPhoton Analysis
//may the force be with you==========================================>
        TFile* datafile1 = TFile::Open("73_Barrel_Plots_updat_SR_TruePVID.root");
        TFile* datafile2 = TFile::Open("82_Barrel_Plots_updat_SR_TruePVID.root");
        TFile* datafile3 = TFile::Open("89_Barrel_Plots_updat_SR_TruePVID.root");
        TFile* datafile5 = TFile::Open("12_Barrel_Plots_updat_SR_TruePVID.root");
        TFile* datafile6 = TFile::Open("11_Barrel_Plots_updat_SR_TruePVID.root");
        TFile* datafile7 = TFile::Open("10_Barrel_Plots_updat_SR_TruePVID.root");


   //   TH1F *h_PhoEtVsPhotonBaselineTriggerEvents = dynamic_cast<TH1F*>(datafile->Get("h_PhoEtVsPhotonBaselineTriggerEvents"));
   //   TH1F *h_PhoEtVsPhotonSingleTriggerEvents = dynamic_cast<TH1F*>(datafile->Get("h_PhoEtVsPhotonSingleTriggerEvents"));

//    TH1F *h_1 = dynamic_cast<TH1F*>(datafile->Get("Lptb"));
    TH1F *h_3 = dynamic_cast<TH1F*>(datafile1->Get("Mpts"));
    TH1F *h_4 = dynamic_cast<TH1F*>(datafile2->Get("Mpts"));
    TH1F *h_1 = dynamic_cast<TH1F*>(datafile3->Get("Mpts"));
    TH1F *h_5 = dynamic_cast<TH1F*>(datafile5->Get("Mpts"));
    TH1F *h_6 = dynamic_cast<TH1F*>(datafile6->Get("Mpts"));
    TH1F *h_7 = dynamic_cast<TH1F*>(datafile7->Get("Mpts"));

h_1->SetStats(0);
h_3->SetStats(0);
h_4->SetStats(0);
h_5->SetStats(0);
h_6->SetStats(0);
h_7->SetStats(0);

    //    TH1F *h_3 = dynamic_cast<TH1F*>(datafile1->Get("Lpts"));
//    TH1F *h_4 = dynamic_cast<TH1F*>(datafile1->Get("Mpts"));
//    TH1F *h_1 = dynamic_cast<TH1F*>(datafile1->Get("Tpts"));
    //h1->GetYaxis()->SetRangeUser(0, 0.01);
   //TCanvas c1("c1","Background Efficiency");
    TCanvas *c1 = new TCanvas("c1","Pt Eff",400,400);

h_3->SetTitle("Signal Efficiency");

    h_3->GetYaxis()->SetRangeUser(0, 1.0);  
    h_3->SetLineColor(kCyan);
    h_3->SetMarkerColor(kCyan);
    h_3->SetMarkerStyle(3);
    h_4->SetLineColor(kOrange );
    h_4->SetMarkerColor(kOrange );
    h_1->SetLineColor(kBlue );
    h_1->SetMarkerColor(kBlue );
    h_5->SetLineColor(kGreen );
    h_5->SetMarkerColor(kGreen );
    h_6->SetLineColor(kRed );
    h_6->SetMarkerColor(kRed );
    h_7->SetLineColor(kViolet );
    h_7->SetMarkerColor(kViolet );



h_3->Draw();

//h_2->Draw("same");
    h_4->Draw("same");
    h_1->Draw("same");

    h_5->Draw("same");
    h_6->Draw("same");
    h_7->Draw("same");

auto mylegend = new TLegend(0.50, 0.15, 0.80, 0.35);
       mylegend->SetBorderSize(0);
       mylegend->SetFillColor(0);
       mylegend->AddEntry(h_1,"High pt Loose ID","lp");
       mylegend->AddEntry(h_4,"High pt Medium ID","lp");
       mylegend->AddEntry(h_3,"High pt Tight ID","lp");
       mylegend->AddEntry(h_7,"Egamma Loose ID","lp");
       mylegend->AddEntry(h_6,"Egamma Medium ID","lp");
       mylegend->AddEntry(h_5,"Egamma Tight ID","lp");
//       mylegend->AddEntry(h_1,"Loose","lp");
//       mylegend->AddEntry(h_4,"Medium","lp");
//       mylegend->AddEntry(h_3,"Tight","lp");
//       mylegend->AddEntry(h_7,"Loose (official)","lp");
//       mylegend->AddEntry(h_6,"Medium (official)","lp");
//       mylegend->AddEntry(h_5,"Tight (official)","lp");
mylegend->SetTextFont(42); 
mylegend->SetTextSize(0.035);
mylegend->SetBorderSize(0);

mylegend->SetFillColor(0);

mylegend->SetFillStyle(0);
mylegend->SetMargin(0.3);




       mylegend->Draw("same");     
  
     c1->SaveAs("comp_egamma_PfworstID_plot_eff.C");
}


