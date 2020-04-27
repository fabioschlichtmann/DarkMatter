//#include "/home/schlichtmann/ALICE/Utils/functions.h"
//#include <string.h>
#include <boost/algorithm/string.hpp>
#include <string>



static TF1* func_Gauss_fit;
static TF1* func_Gauss_fit2;
static TF1* func_polynom;
static TF1* func_polynom2;

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,
                      Float_t size=0.06,Int_t color=1,Float_t angle=0.0,
                      Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color
    // align: 1 left aligned, 32, right aligned

    // align = 10*HorizontalAlign + VerticalAlign
    // For horizontal alignment the following convention applies:
    // 1=left adjusted, 2=centered, 3=right adjusted
    // For vertical alignment the following convention applies:
    // 1=bottom adjusted, 2=centered, 3=top adjusted

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}

void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}

Double_t GaussFitFunc_plus2ndorderpolynom(Double_t* x_val, Double_t* par)
{
    Double_t x, y, A, mean, sigma, par3,par4,par5;
    A  = fabs(par[0]);
    mean  = par[1];
    sigma  = fabs(par[2]);
    par3  = par[3];
    par4  = par[4];
    par5 = par[5];
    x = x_val[0];
    y = A*TMath::Gaus(x,mean,sigma,0) + par3*x*x+par4*x+par5;
    return y;
}

Double_t polynom(Double_t* x_val, Double_t* par)
{
    double x,par0,par1,par2;
    x=x_val[0];
    par0=par[0];
    par1=par[1];
    par2=par[2];
    return par0*x*x+par1*x+par2;
}

void Gaus_for_K0(TH1D* blubber)
{
    cout<<"Gaus started"<<endl;
    func_Gauss_fit             = new TF1("func_Gauss_fit",GaussFitFunc_plus2ndorderpolynom,0,0.6,6);
    int maxbin=blubber->GetMaximumBin();
    double amplitude=blubber->GetBinContent(maxbin);
    double mean=blubber->GetBinCenter(maxbin);
    double sigma=0.005;
    double par3,par4,par5;


    //Parameter auf 0 setzen
    for(Int_t i = 0; i < 6; i++)
       {
           func_Gauss_fit ->ReleaseParameter(i);
           func_Gauss_fit ->SetParameter(i,0.0);
           func_Gauss_fit ->SetParError(i,0.0);
       }

       func_Gauss_fit ->SetParameter(0,amplitude);
       func_Gauss_fit ->SetParameter(1,mean);
       func_Gauss_fit ->SetParameter(2,sigma);
       func_Gauss_fit ->SetParameter(3,0.0);
       func_Gauss_fit ->SetParameter(5,15000.);

       //Fit
       blubber ->Fit("func_Gauss_fit","QMN","",0.4,0.6);

       //Parameter auslesen
       amplitude = func_Gauss_fit ->GetParameter(0);
       mean      = func_Gauss_fit ->GetParameter(1);
       sigma     = fabs(func_Gauss_fit ->GetParameter(2));
       par3   =  func_Gauss_fit ->GetParameter(3);
       par4   =  func_Gauss_fit ->GetParameter(4);
       par5   =  func_Gauss_fit ->GetParameter(5);

       cout<<amplitude<<endl;
       cout<<mean<<endl;
       cout<<sigma<<endl;

       //erneut fitten
       //Parameter auf 0 setzen
        for(Int_t i = 0; i < 6; i++)
       {
           func_Gauss_fit ->ReleaseParameter(i);
           func_Gauss_fit ->SetParameter(i,0.0);
           func_Gauss_fit ->SetParError(i,0.0);
       }

       //Parameter auf vorherige Fitwerte
       func_Gauss_fit ->SetParameter(0,amplitude);
       func_Gauss_fit ->SetParameter(1,mean);
       func_Gauss_fit ->SetParameter(2,sigma);
       func_Gauss_fit ->SetParameter(3,par3);
       func_Gauss_fit ->SetParameter(4,par4);
       func_Gauss_fit ->SetParameter(5,par5);
       //func_Gauss_fit ->SetParameter(3,0.0);

       //fitten
       blubber ->Fit("func_Gauss_fit","QMN","",0.4,0.6);

       //Parameter auslesen
       amplitude = func_Gauss_fit ->GetParameter(0);
       mean      = func_Gauss_fit ->GetParameter(1);
       sigma     = fabs(func_Gauss_fit ->GetParameter(2));

       cout<<""<<endl;
       cout<<amplitude<<endl;
       cout<<mean<<endl;
       cout<<sigma<<endl;
       cout<<par3<<endl;
       cout<<par4<<endl;
       cout<<par5<<endl;
       //etwas bessere Fitwerte

       func_Gauss_fit ->SetLineColor(kRed);
       func_Gauss_fit ->SetLineStyle(1);
       func_Gauss_fit ->SetRange(mean-4.0*sigma,mean+4.0*sigma);
       func_Gauss_fit ->DrawClone("same");

}


void Gaus_for_lambda(TH1D* blubber)
{
    cout<<"Gaus started"<<endl;
    func_Gauss_fit2             = new TF1("func_Gauss_fit2",GaussFitFunc_plus2ndorderpolynom,0,0.6,6);
    int maxbin=blubber->GetMaximumBin();
    double amplitude=blubber->GetBinContent(maxbin);
    double mean=blubber->GetBinCenter(maxbin);
    double sigma=0.005;
    double par3,par4,par5;


    //Parameter auf 0 setzen
    for(Int_t i = 0; i < 6; i++)
       {
           func_Gauss_fit2 ->ReleaseParameter(i);
           func_Gauss_fit2 ->SetParameter(i,0.0);
           func_Gauss_fit2 ->SetParError(i,0.0);
       }

       func_Gauss_fit2 ->SetParameter(0,amplitude);
       func_Gauss_fit2 ->SetParameter(1,mean);
       func_Gauss_fit2 ->SetParameter(2,sigma);
       func_Gauss_fit2 ->SetParameter(3,0.0);
       func_Gauss_fit2 ->SetParameter(5,600.);

       //Fit
       blubber ->Fit("func_Gauss_fit2","QMN","",mean-3.0*sigma,mean+3.0*sigma);

       //Parameter auslesen
       amplitude = func_Gauss_fit2 ->GetParameter(0);
       mean      = func_Gauss_fit2 ->GetParameter(1);
       sigma     = fabs(func_Gauss_fit2 ->GetParameter(2));
       par3   =  func_Gauss_fit2 ->GetParameter(3);
       par4   =  func_Gauss_fit2 ->GetParameter(4);
       par5   =  func_Gauss_fit2 ->GetParameter(5);

       cout<<amplitude<<endl;
       cout<<mean<<endl;
       cout<<sigma<<endl;

       //erneut fitten
       //Parameter auf 0 setzen
        for(Int_t i = 0; i < 6; i++)
       {
           func_Gauss_fit2 ->ReleaseParameter(i);
           func_Gauss_fit2 ->SetParameter(i,0.0);
           func_Gauss_fit2 ->SetParError(i,0.0);
       }

       //Parameter auf vorherige Fitwerte
       func_Gauss_fit2 ->SetParameter(0,amplitude);
       func_Gauss_fit2 ->SetParameter(1,mean);
       func_Gauss_fit2 ->SetParameter(2,sigma);
       func_Gauss_fit2 ->SetParameter(3,par3);
       func_Gauss_fit2 ->SetParameter(4,par4);
       func_Gauss_fit2 ->SetParameter(5,par5);
       //func_Gauss_fit2 ->SetParameter(3,0.0);

       //fitten
       blubber ->Fit("func_Gauss_fit2","QMN","",mean-4.0*sigma,mean+4.0*sigma);

       //Parameter auslesen
       amplitude = func_Gauss_fit2 ->GetParameter(0);
       mean      = func_Gauss_fit2 ->GetParameter(1);
       sigma     = fabs(func_Gauss_fit2 ->GetParameter(2));

       cout<<""<<endl;
       cout<<amplitude<<endl;
       cout<<mean<<endl;
       cout<<sigma<<endl;
       cout<<par3<<endl;
       cout<<par4<<endl;
       cout<<par5<<endl;
       //etwas bessere Fitwerte

       func_Gauss_fit2 ->SetLineColor(kRed);
       func_Gauss_fit2 ->SetLineStyle(1);
       //func_Gauss_fit2 ->SetRange(mean-4.0*sigma,mean+4.0*sigma);
       func_Gauss_fit2 ->SetRange(mean-4.0*sigma,mean+4.0*sigma);
       func_Gauss_fit2 ->DrawClone("same");

}


void Drawnicehisto()
{
    TFile* inputfile = TFile::Open("./Histos.root");
    TH1D* histo_invariantmass_lambda = (TH1D*)inputfile->Get("histo inv mass lambda");
    TH1D* histo_invariantmass_K0 = (TH1D*)inputfile->Get("histo inv mass K0");

    TCanvas* a = new TCanvas();
    TCanvas* b = new TCanvas();

    a->cd();
    SetRootGraphicStyle();
    a->SetTicks(1,1);
    gStyle->SetLegendBorderSize(0);
  
        
    histo_invariantmass_lambda->GetXaxis()->SetTitle("invariant mass (p,#pi^{-} and antip,#pi^{+}) [GeV/c^{2}]");
    //histo_invariantmass_lambda->SetTitle("K0 peak");

    histo_invariantmass_lambda->SetMarkerColor(kBlack);
    histo_invariantmass_lambda->SetMarkerSize(0.5);
    histo_invariantmass_lambda->SetMarkerStyle(20);
    histo_invariantmass_lambda->SetTitle("");
    histo_invariantmass_lambda->GetYaxis()->SetRangeUser(0,2000);
    histo_invariantmass_lambda->GetYaxis()->SetTitle("counts");
    histo_invariantmass_lambda->GetXaxis()->SetTitleSize(0.050);
    histo_invariantmass_lambda->GetYaxis()->SetMaxDigits(3);
    histo_invariantmass_lambda->GetYaxis()->SetTitleSize(0.050);
    histo_invariantmass_lambda->GetXaxis()->SetTitleOffset(0.9);
    histo_invariantmass_lambda->GetYaxis()->SetTitleOffset(0.7);
    histo_invariantmass_lambda->GetXaxis()->CenterTitle();
    histo_invariantmass_lambda->GetYaxis()->CenterTitle();
    histo_invariantmass_lambda->SetLineWidth(1);
    histo_invariantmass_lambda->SetLineColor(kBlack);
  
    gStyle->SetOptStat(0);



    histo_invariantmass_lambda->Draw("P E1");

    Gaus_for_lambda( histo_invariantmass_lambda);

    func_polynom2          = new TF1("func_polynom",polynom,0,1.2,3);
    func_polynom2->SetParameter(0,func_Gauss_fit2 ->GetParameter(3));
    func_polynom2->SetParameter(1,func_Gauss_fit2 ->GetParameter(4));
    func_polynom2->SetParameter(2,func_Gauss_fit2 ->GetParameter(5));
    func_polynom2 ->SetLineColor(kBlue);
    func_polynom2 ->SetLineStyle(2);
    func_polynom2 ->DrawClone("same");

    TString Inhalt,Inhalt1, Inhalt2, Inhalt3, Inhalt4;
    Inhalt="#Lambda -> p + #pi^{-}";
    Inhalt1="#Lambda -> antip + #pi^{+}";
    Inhalt2="p + Pb";
    Inhalt3="#mu: ";
    double mean     = fabs(func_Gauss_fit2 ->GetParameter(1));
    double sigma     = fabs(func_Gauss_fit2 ->GetParameter(2));
    cout<<mean<<endl;
    Inhalt3+=to_string(mean);
    Inhalt3+=" [GeV/c^{2}]";
    Inhalt4="#sigma: ";
    Inhalt4+=to_string(sigma);
    Inhalt4+=" [GeV/c^{2}]";
    plotTopLegend((char*)Inhalt.Data(),0.2,0.7+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt1.Data(),0.2,0.62+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt2.Data(),0.2,0.54+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt3.Data(),0.2,0.46+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt4.Data(),0.2,0.38+0.1,0.040,kBlack,0.0,42,1,11);

    auto legend_lambda = new TLegend(0.6,0.6,0.85,0.85);
    //legend_lambda->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    //legend_lambda->AddEntry(histo_invariantmass_K0,"Histogram filled with random numbers","lep");
    legend_lambda->AddEntry("func_Gauss_fit2","Gaussian fit","l");
    legend_lambda->AddEntry("func_polynom","background fit","l");
    //legend_lambda->AddEntry("gr","Graph with error bars","lep");
    legend_lambda->Draw();
    //histo_invariantmass_K0->SetTitle("K0 peak");

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //histo_invariantmass_lambda->Draw("P");



    b->cd();
    SetRootGraphicStyle();
    b->SetTicks(1,1);
    b->SetBottomMargin(3);
    gStyle->SetPadBottomMargin(0.5);
  
        
    histo_invariantmass_K0->GetXaxis()->SetTitle("invariant mass (#pi^{+},#pi^{-}) [GeV/c^{2}]");

  

    histo_invariantmass_K0->SetMarkerColor(kBlack);
    histo_invariantmass_K0->SetMarkerSize(0.5);
    histo_invariantmass_K0->SetMarkerStyle(20);
    histo_invariantmass_K0->SetTitle("");
    histo_invariantmass_K0->GetYaxis()->SetRangeUser(0,23000);
    histo_invariantmass_K0->GetYaxis()->SetTitle("counts");
    histo_invariantmass_K0->GetXaxis()->SetTitleSize(0.050);
    histo_invariantmass_K0->GetYaxis()->SetMaxDigits(3);
    histo_invariantmass_K0->GetYaxis()->SetTitleSize(0.050);
    histo_invariantmass_K0->GetXaxis()->SetTitleOffset(0.9);
    histo_invariantmass_K0->GetYaxis()->SetTitleOffset(0.7);
    histo_invariantmass_K0->GetXaxis()->CenterTitle();
    histo_invariantmass_K0->GetYaxis()->CenterTitle();
    histo_invariantmass_K0->SetLineWidth(1);
    histo_invariantmass_K0->SetLineColor(kBlack);
  
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    

    histo_invariantmass_K0->Draw("P E1");
    Gaus_for_K0(histo_invariantmass_K0);
    Inhalt="K^{0}_{s} -> #pi^{+} + #pi^{-}";
    Inhalt2="p + Pb";
    Inhalt3="#mu: ";
     mean     = fabs(func_Gauss_fit ->GetParameter(1));
     sigma     = fabs(func_Gauss_fit ->GetParameter(2));
    cout<<mean<<endl;
    Inhalt3+=to_string(mean);
    Inhalt3+=" [GeV/c^{2}]";
    Inhalt4="#sigma: ";
    Inhalt4+=to_string(sigma);
    Inhalt4+=" [GeV/c^{2}]";
    plotTopLegend((char*)Inhalt.Data(),0.2,0.7+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt2.Data(),0.2,0.62+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt3.Data(),0.2,0.54+0.1,0.040,kBlack,0.0,42,1,11);
    plotTopLegend((char*)Inhalt4.Data(),0.2,0.46+0.1,0.040,kBlack,0.0,42,1,11);

  
  

    func_polynom          = new TF1("func_polynom",polynom,0,0.6,3);
    func_polynom->SetParameter(0,func_Gauss_fit ->GetParameter(3));
    func_polynom->SetParameter(1,func_Gauss_fit ->GetParameter(4));
    func_polynom->SetParameter(2,func_Gauss_fit ->GetParameter(5));
    func_polynom ->SetLineColor(kBlue);
    func_polynom ->SetLineStyle(2);
    func_polynom ->DrawClone("same");
    auto legend = new TLegend(0.6,0.6,0.85,0.85);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    //legend->AddEntry(histo_invariantmass_K0,"Histogram filled with random numbers","lep");
    legend->AddEntry("func_Gauss_fit","Gaussian fit","l");
    legend->AddEntry("func_polynom","background fit","l");
    //legend->AddEntry("gr","Graph with error bars","lep");
    legend->Draw();

    

}