#ifndef HIST_UTILS_H
#define HIST_UTILS_H

//Utility function that takes as argument 2 histogram objects, and the xlabel, ylabel (assume same range) and plots their ratio
void hist_ratio(TH1F *hdata, TH1F *hsimc, TString xlabel="", TString ylabel="", TString title="")
{
  int font_type = 132;
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.06);

  //Create Canvas
  TCanvas *c = new TCanvas(title, "c", 700, 700);
  c->Divide(1,2, 0., 0.0);

  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  double dataI;  //data integral
  double simcI;  //simc integral
  
  //Set Histo Aesthetics
  hsimc->SetLineWidth(2);
  hsimc->SetLineColor(kRed);

  hdata->SetFillColorAlpha(kBlue, 0.35);
  hdata->SetFillStyle(3004);
  hdata->SetLineWidth(2);
  //Create new histo to store ratio
  TH1F *hratio = (TH1F*) hdata->Clone();

  //Set Histos Axis Label Size
  hdata->GetYaxis()->SetLabelSize(0.06);
  hratio->GetYaxis()->SetLabelSize(0.05);
  hratio->GetXaxis()->SetLabelSize(0.06);
  hdata->SetTitleSize(0.06, "XY");
  hratio->SetTitleSize(0.06, "XY");

  hratio->SetLineWidth(2);
  hratio->SetLineColor(kBlack);

  
  c->cd(1);
  //auto leg = new TLegend(0.1,0.8,0.28,0.9); 
  auto leg = new TLegend(0.14,0.84,0.25,0.64); 

  TPad* pad1 = (TPad*)c->GetPad(1);
  pad1->SetBottomMargin(0.02);
  pad1->SetRightMargin(0.05);
  pad1->SetLeftMargin(0.1);
  pad1->SetTopMargin(0.13);
  pad1->SetFrameLineWidth(2);
  hdata->Draw("samehistE0");
  hsimc->Draw("samesE0");
  hdata->SetTitle(title);

  hdata->GetXaxis()->SetLabelSize(0);
  hdata->GetYaxis()->SetTitle(ylabel);
  hdata->GetYaxis()->CenterTitle();
  hdata->GetYaxis()->SetRangeUser(0,hdata->GetMaximum()+0.6*hdata->GetMaximum());
  hdata->GetYaxis()->SetTitleOffset(0.6);
  hdata->SetLabelFont(font_type, "XY");
  hdata->SetTitleFont(font_type, "XY");

  double dataI_err, simcI_err;
  double nbins = hdata->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = hdata->IntegralAndError(1, nbins, dataI_err);
  simcI = hsimc->IntegralAndError(1, nbins, simcI_err);
  double R = (float)dataI / simcI;
  double R_err = R * sqrt(pow(dataI_err/dataI,2) + pow(simcI_err/simcI,2));
  
  leg->AddEntry(hdata,Form("Data | Integral: %.3f", dataI),"f");
  leg->AddEntry(hsimc,Form("SIMC | Integral: %.3f", simcI));
  leg->AddEntry((TObject*)0, Form("Ratio: %.3f #pm %.3f", R, R_err), "");
  leg->Draw();

  c->cd(2);
  TPad* pad = (TPad*)c->GetPad(2);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.2);
  pad->SetRightMargin(0.05);
  pad->SetLeftMargin(0.1);
  pad->SetFrameLineWidth(2);

  hratio->Divide(hdata, hsimc);
  hratio->GetYaxis()->SetRangeUser(0.,2.1);
  hratio->SetTitle("");

  hratio->GetXaxis()->SetTitle(xlabel);
  hratio->GetYaxis()->SetTitle("Y_{data} / Y_{SIMC}");
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();
  hratio->GetYaxis()->SetTitleOffset(0.6);
  hratio->SetLabelFont(font_type, "XY");
  hratio->SetTitleFont(font_type, "XY");

  hratio->Draw();
  pad->Update();
  //Draw lines at +/- 10 % of 1.
  TLine* lmin_10 = new TLine(pad->GetUxmin(),0.9,pad->GetUxmax(),0.9);
  TLine* lmax_10 = new TLine(pad->GetUxmin(),1.1,pad->GetUxmax(),1.1);

  lmin_10->SetLineColor(kBlack);
  lmin_10->SetLineStyle(2);
  lmin_10->Draw();
  lmax_10->SetLineColor(kBlack);
  lmax_10->SetLineStyle(2);
  lmax_10->Draw();
  
  TLine* lmin_20 = new TLine(pad->GetUxmin(),0.8,pad->GetUxmax(),0.8);
  TLine* lmax_20 = new TLine(pad->GetUxmin(),1.2,pad->GetUxmax(),1.2);

  lmin_20->SetLineColor(kBlue);
  lmin_20->SetLineStyle(2);
  lmin_20->Draw();
  lmax_20->SetLineColor(kBlue);
  lmax_20->SetLineStyle(2);
  lmax_20->Draw();
}

void plot_hist(TH1F *hist, TString xlabel="", TString ylabel="", TString title="", TString set_logy="")
{
  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);

  //Create Canvas
  TCanvas *c = new TCanvas(title, "c", 900, 700);
  c->SetFrameLineWidth(2);
  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  double dataI;  //data, simc integral
  double dataI_err;


  hist->SetFillColorAlpha(kBlue, 0.35);
  hist->SetFillStyle(3004);

  //Set Histos Axis Label Size
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->SetTitleSize(0.04, "XY");

  TLegend *leg = new TLegend(0.14,0.8,0.4,0.7); 

  if(set_logy=="logy"){
    c->SetLogy();
  }

  hist->Draw("samehistE0");

  
  hist->SetTitle(title);
  
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->SetTitle(xlabel);

  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->CenterTitle();

  
  //if(set_logy==""){hist->GetYaxis()->SetRangeUser(0,hist->GetMaximum()+0.6*hist->GetMaximum());}

  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleOffset(1.3);
  
  hist->SetLabelFont(font_type, "XY");
  hist->SetTitleFont(font_type, "XY");

  //Get DATA, SIMC histogram names to determine which range to integrate
  TString H_name = hist->GetName();

  double xmin, xmax;
  if(H_name=="H_pcal_etotTrkNorm_sys")
    { xmin = hist->FindBin(0.7), xmax = hist->FindBin(5.0);}
  else if(H_name=="H_ctime_sys")
    { xmin = hist->FindBin(10.5), xmax = hist->FindBin(14.5);}

  
  double nbins = hist->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = hist->IntegralAndError(xmin, xmax, dataI_err);
  
  leg->AddEntry(hist,Form("Data | Integral: %.3f #pm %.3f", dataI, dataI_err),"f");
  leg->Draw();

}



//---------------------

void compare_hist(TH1F *hdata, TH1F *hsimc, TString xlabel="", TString ylabel="", TString title="", TString set_logy="", Int_t norm=0)
{

  int font_type = 132;
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);

  //Create Canvas
  TCanvas *c = new TCanvas(title, "c", 900, 700);
  c->SetFrameLineWidth(2);

  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  double dataI;  //data integral
  double simcI;  //simc integral
  
  //Set Histo Aesthetics
  hsimc->SetLineWidth(2);
  hsimc->SetLineColor(kRed);

  hdata->SetFillColorAlpha(kBlue, 0.35);
  hdata->SetFillStyle(3004);

  //Set Histos Axis Label Size
  hdata->GetYaxis()->SetLabelSize(0.04);
  hdata->GetXaxis()->SetLabelSize(0.04);
  hdata->SetTitleSize(0.04, "XY");
  
  hdata->SetTitle(title);

  hdata->GetYaxis()->SetTitle(ylabel);
  hdata->GetXaxis()->SetTitle(xlabel);
  hdata->GetYaxis()->CenterTitle();
  hdata->GetXaxis()->CenterTitle();
  hdata->GetYaxis()->SetRangeUser(0.001,hdata->GetMaximum()+0.6*hdata->GetMaximum());
  hdata->GetXaxis()->SetTitleOffset(1.);  
  hdata->GetYaxis()->SetTitleOffset(1.);
  hdata->SetLabelFont(font_type, "XY");
  hdata->SetTitleFont(font_type, "XY");
  
  //auto leg = new TLegend(0.1,0.8,0.28,0.9); 
  TLegend *leg = new TLegend(0.14,0.88,0.25,0.73); 
  if(set_logy=="logy"){
    c->SetLogy();
  }

  if(norm==0)
    {
      hdata->Draw("samehistE0");
      hsimc->Draw("samesE0");
    }
  
  else if(norm==1)
    {
      hdata->DrawNormalized("samehistE0");
      hsimc->DrawNormalized("samesE0");
    }



  //Get DATA, SIMC histogram names to determine which range to integrate
  TString H_data_name = hdata->GetName();
  TString H_simc_name = hsimc->GetName();

  double xmin, xmax;
  double xmin_simc, xmax_simc;
  if(H_data_name=="H_Em_nuc_sys")
    { xmin = 14., xmax = 36.;}
  else if(H_data_name=="H_Q2_sys")
    { xmin = hdata->FindBin(4.0), xmax = (hdata->FindBin(5.0))-1.0;}
  else if(H_data_name=="H_hdelta_sys")
    { xmin = hdata->FindBin(-8.0), xmax = hdata->FindBin(7.99);}
  else if(H_data_name=="H_edelta_sys")
    { xmin = hdata->FindBin(-10.0), xmax = hdata->FindBin(21.99);}
  else if(H_data_name=="H_ztar_diff_sys")
    { xmin = 23., xmax = 35.,  xmin_simc = 24., xmax_simc = 36;}

  
  double dataI_err, simcI_err;
  double nbins = hdata->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = 22.653; // hdata->IntegralAndError(xmin, xmax, dataI_err);
  simcI = 23.908; //hsimc->IntegralAndError(xmin, xmax, simcI_err);
  dataI_err = 0.449;
  simcI_err = 0.428;
  
  double R = (float)dataI / simcI;
  double R_err = R * sqrt(pow(dataI_err/dataI,2) + pow(simcI_err/simcI,2));
  
  leg->AddEntry(hdata,Form("Data | Integral: %.3f #pm %.3f", dataI, dataI_err),"f");
  leg->AddEntry(hsimc,Form("SIMC | Integral: %.3f #pm %.3f", simcI, simcI_err));
  leg->AddEntry((TObject*)0, Form("Ratio: %.3f #pm %.3f", R, R_err), "");
  leg->Draw();


}



//-----------------------------------------------
void combine_sets(TH1F *hist80, TH1F *hist_580set1, TH1F *hist_580set2, TH1F *hist_750set1, TH1F *hist_750set2, TH1F *hist_750set3, TString xlabel="", TString ylabel="", TString title="", Double_t scale_factor=1, Bool_t combine=0)
{

  //This function adds histograms of multiple data sets from the same kinematic setting (i.e., 580_set1 + 580_set2)
 

  //Clone histograms and rename them total
  TH1F *hist_580tot = (TH1F*)hist_580set2->Clone("hist_580tot");
  TH1F *hist_750tot = (TH1F*)hist_750set1->Clone("hist_750tot");

    //Add remaining histograms to total
    //hist_580tot->Add(hist_580set2);

    //hist_750tot->Add(hist_750set2);
    //hist_750tot->Add(hist_750set3);
  
  
  int font_type = 132;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleFont(font_type, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(font_type);
  gStyle->SetLegendTextSize(0.03);

  //Create Canvas
  TCanvas *c = new TCanvas(title, "c", 900, 700);
  c->SetFrameLineWidth(2);
  c->cd();

  //Set Frame Margin spacing (smaller spacing->tighter frame)
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.12);
  
  hist80->SetLineColor(kBlue);
  hist80->SetLineWidth(2);
  hist80->SetFillColorAlpha(kBlue, 0.6);
  hist80->SetFillStyle(3004);
  hist80->SetMarkerStyle(21);
  hist80->SetMarkerColor(kBlue);

  hist_580tot->SetLineColor(kRed);
  hist_580tot->SetLineWidth(2);
  hist_580tot->SetFillColorAlpha(kRed, 0.6);
  hist_580tot->SetFillStyle(3004);
  hist_580tot->SetMarkerStyle(21);
  hist_580tot->SetMarkerColor(kRed);

  hist_750tot->SetLineColor(kMagenta);
  hist_750tot->SetLineWidth(2);
  hist_750tot->SetFillColorAlpha(kMagenta, 0.6);
  hist_750tot->SetFillStyle(3005);
  hist_750tot->SetMarkerStyle(21);
  hist_750tot->SetMarkerColor(kMagenta);
 
  //Set Histos Axis Label Size
  hist80->GetYaxis()->SetLabelSize(0.05);
  hist80->SetTitleSize(0.05, "XY");
  
  //Draw combined histos
  hist80->Scale(1./scale_factor);  //scale down the 80MeV setting so it may be comaprable to other settings
  hist80->Draw("samehistE0");
  hist_580tot->Draw("samehistE0");
  hist_750tot->Draw("samehistE0");

  hist80->SetTitle(title);
  
  hist80->GetXaxis()->SetLabelSize(0.05);
  hist80->GetYaxis()->SetLabelSize(0.05);
  
  hist80->GetYaxis()->SetTitle(ylabel);
  hist80->GetXaxis()->SetTitle(xlabel);

  hist80->GetYaxis()->CenterTitle();
  hist80->GetXaxis()->CenterTitle();
  
  hist80->GetYaxis()->SetTitleOffset(1.35);
  hist80->GetXaxis()->SetTitleOffset(1.);
  
  hist80->SetLabelFont(font_type, "XY");
  hist80->SetTitleFont(font_type, "XY");


}


#endif
