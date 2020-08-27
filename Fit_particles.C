Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, A, mean, sigma, par3,par4,par5;
    A  = fabs(par[0]);
    mean  = par[1];
    sigma  = fabs(par[2]);
    par3  = par[3];
    par4  = par[4];
    par5 = par[5];
    x = x_val[0];
    y = A*TMath::Gaus(x,mean,sigma,0);
    return y;
}

double tdistribution_pdf(double* x_val,double* par) {

    double x = x_val[0];
    double r = par[0];
    double x0 = par[1];
    double A = par[2];
    double b = par[3];

    x = x*b;

    return A*( (std::exp (ROOT::Math::lgamma((r + 1.0)/2.0) - ROOT::Math::lgamma(r/2.0)) / std::sqrt (M_PI * r))
              * std::pow ((1.0 + (x-x0)*(x-x0)/r), -(r + 1.0)/2.0) );


}

double double_fit_tdistribution_pdf(double* x_val,double* par) {

    double x = x_val[0];
    double r_1 = par[0];
    double x0_1 = par[1];
    double A_1 = par[2];
    double b_1 = par[3];

    double r_2 = par[4];
    double x0_2 = par[5];
    double A_2 = par[6];
    double b_2 = par[7];

   // x = x*b;

    return A_1*( (std::exp (ROOT::Math::lgamma((r_1 + 1.0)/2.0) - ROOT::Math::lgamma(r_1/2.0)) / std::sqrt (M_PI * r_1))
                * std::pow ((1.0 + (x*b_1-x0_1)*(x*b_1-x0_1)/r_1), -(r_1 + 1.0)/2.0) )
        +  A_2*( (std::exp (ROOT::Math::lgamma((r_2 + 1.0)/2.0) - ROOT::Math::lgamma(r_2/2.0)) / std::sqrt (M_PI * r_2))
              * std::pow ((1.0 + (x*b_2-x0_2)*(x*b_2-x0_2)/r_2), -(r_2 + 1.0)/2.0) );
 
}

double triple_fit_tdistribution_pdf(double* x_val,double* par) {

    double x = x_val[0];

    double r_1 = par[0];
    double x0_1 = par[1];
    double A_1 = par[2];
    double b_1 = par[3];

    double r_2 = par[4];
    double x0_2 = par[5];
    double A_2 = par[6];
    double b_2 = par[7];

    double r_3 = par[8];
    double x0_3 = par[9];
    double A_3 = par[10];
    double b_3 = par[11];

    // x = x*b;

    double res1 = A_1*( (std::exp (ROOT::Math::lgamma((r_1 + 1.0)/2.0) - ROOT::Math::lgamma(r_1/2.0)) / std::sqrt (M_PI * r_1))
                         * std::pow ((1.0 + (x*b_1-x0_1)*(x*b_1-x0_1)/r_1), -(r_1 + 1.0)/2.0) );
    double res2 = A_2*( (std::exp (ROOT::Math::lgamma((r_2 + 1.0)/2.0) - ROOT::Math::lgamma(r_2/2.0)) / std::sqrt (M_PI * r_2))
                * std::pow ((1.0 + (x*b_2-x0_2)*(x*b_2-x0_2)/r_2), -(r_2 + 1.0)/2.0) );
    double res3 = A_3*( (std::exp (ROOT::Math::lgamma((r_3 + 1.0)/2.0) - ROOT::Math::lgamma(r_3/2.0)) / std::sqrt (M_PI * r_3))
                       * std::pow ((1.0 + (x*b_3-x0_3)*(x*b_3-x0_3)/r_3), -(r_3 + 1.0)/2.0) ) ;

    //if(isnan(res1)){res1 =0; }
    if(isnan(res2)){res2= 0; }
    if(isnan(res3)){res3 = 0; }

    //cout<<"results: "<<endl;
    //cout<<res1<<endl;
    //cout<<res2<<endl;
    //cout<<res3<<endl;

    double result = res1 + res2 + res3;
    

    return result;

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
    /*gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    */
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


TF1* Gauss(TH1D* blubber, double& sigma, double& mean,bool draw,double x_min,double x_max,int mode=0,int save=0)
{
    //cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    //cout<<"Gaus started"<<endl;

    TF1* func_Gauss_fit3;
    func_Gauss_fit3             = new TF1("func_Gauss_fit3",GaussFitFunc,-0.4,1.4,6);

    int maxbin;
    double amplitude;
    if(mode==0)
    {
        maxbin=blubber->GetMaximumBin();
        amplitude=blubber->GetBinContent(maxbin);
        //mean=blubber->GetBinCenter(maxbin);
    }
    if(mode==1)
    {
        mean = 1.1157;
        maxbin = blubber->FindBin(mean);
        amplitude=blubber->GetBinContent(maxbin);
    }
    if(mode==2)
    {
        mean = 0.4981;
        maxbin = blubber->FindBin(mean);
        amplitude=blubber->GetBinContent(maxbin);
    }

    double par3,par4,par5;

    //Parameter auf 0 setzen
    for(Int_t i = 0; i < 6; i++)
    {
        func_Gauss_fit3 ->ReleaseParameter(i);
        func_Gauss_fit3 ->SetParameter(i,0.0);
        func_Gauss_fit3 ->SetParError(i,0.0);
    }

    func_Gauss_fit3 ->SetParameter(0,amplitude);
    func_Gauss_fit3 ->SetParameter(1,mean);
    func_Gauss_fit3 ->SetParameter(2,sigma);
    func_Gauss_fit3 ->SetParameter(3,0.);
    func_Gauss_fit3 ->SetParameter(4,0.);
    func_Gauss_fit3 ->SetParameter(5,0.);

    //Fit
    blubber ->Fit("func_Gauss_fit3","QMN","",mean-2*sigma,mean+2*sigma);

    //Parameter auslesen
    amplitude = func_Gauss_fit3 ->GetParameter(0);
    mean      = func_Gauss_fit3 ->GetParameter(1);
    sigma     = fabs(func_Gauss_fit3 ->GetParameter(2));
    par3   =  func_Gauss_fit3 ->GetParameter(3);
    par4   =  func_Gauss_fit3 ->GetParameter(4);
    par5   =  func_Gauss_fit3 ->GetParameter(5);

    /*
    cout<<amplitude<<endl;
    cout<<mean<<endl;
    cout<<sigma<<endl;
    */

    //erneut fitten
    //Parameter auf 0 setzen
    for(Int_t i = 0; i < 6; i++)
    {
        //if(i==3){continue;}
        func_Gauss_fit3 ->ReleaseParameter(i);
        func_Gauss_fit3 ->SetParameter(i,0.0);
        func_Gauss_fit3 ->SetParError(i,0.0);
    }

    //Parameter auf vorherige Fitwerte
    func_Gauss_fit3 ->SetParameter(0,amplitude);
    func_Gauss_fit3 ->SetParameter(1,mean);
    func_Gauss_fit3 ->SetParameter(2,sigma);
    func_Gauss_fit3 ->SetParameter(3,par3);
    func_Gauss_fit3 ->SetParameter(4,par4);
    func_Gauss_fit3 ->SetParameter(5,par5);
    //func_Gauss_fit3 ->SetParameter(3,0.0);

    //fitten
    //TFitResultPtr r = blubber ->Fit("func_Gauss_fit3","QMNS","",mean-7.0*sigma,mean+7.0*sigma);
    TFitResultPtr r = blubber ->Fit("func_Gauss_fit3","QMNS","",mean-2*sigma,mean+2*sigma);
    //TF1* fit = blubber->GetFunction("func_Gauss_fit3");
    //double chi2 = fit->GetChisquare();
    //cout<<"chi2: "<<chi2<<endl;
    //blubber ->Fit("func_Gauss_fit3","QMN","",1.105,1.127);
    Double_t chi2   = r->Chi2();
    cout<<"chi2: "<<chi2<<endl;

    double* error_array = new double[3];
    error_array[0]=-1;
    error_array[1]=-1;
    error_array[2]=-1;

    //if(chi2>200){return error_array;}

    //Parameter auslesen
    amplitude = func_Gauss_fit3 ->GetParameter(0);
    mean      = func_Gauss_fit3 ->GetParameter(1);
    sigma     = fabs(func_Gauss_fit3 ->GetParameter(2));
    par3   =  func_Gauss_fit3 ->GetParameter(3);
    par4   =  func_Gauss_fit3 ->GetParameter(4);
    par5   =  func_Gauss_fit3 ->GetParameter(5);

    cout<<"amplitude: "<<amplitude<<endl;
    cout<<"mean: "<<mean<<endl;
    cout<<"sigma: "<<sigma<<endl;
    /*
    cout<<par3<<endl;
    cout<<par4<<endl;
    cout<<par5<<endl;
    */
    cout<<"par3: "<<par3<<endl;
    cout<<"par4: "<<par4<<endl;
    cout<<"par5: "<<par5<<endl;
    //etwas bessere Fitwerte
    

    func_Gauss_fit3 ->SetLineColor(kRed);
    func_Gauss_fit3 ->SetLineStyle(1);
    func_Gauss_fit3 ->SetRange(-0.4,1.4);
    //func_Gauss_fit3 ->SetRange(1.105,1.127);
    if(draw) {func_Gauss_fit3 ->Draw("same");}

    
    //if(draw) {Gauss_only->Draw("same"); }

    par3   =  func_Gauss_fit3 ->GetParameter(3);
    par4   =  func_Gauss_fit3 ->GetParameter(4);
    par5   =  func_Gauss_fit3 ->GetParameter(5);

   
    if(draw)
    {
        TCanvas* a = new TCanvas();
        blubber->Draw("P E1");
        func_Gauss_fit3 ->Draw("same");
        //polynomial_only->Draw("same");

        if(save==1)
        {
            TString save_name = "K0_radius>10";
            save_name+=".png";
            a->SaveAs(save_name.Data());

        }
        if(save==2)
        {
            TString save_name = "Lambda_radius>10";
            save_name+=".png";
            a->SaveAs(save_name.Data());

        }
        if(save==3)
        {
            TString save_name = "K0_radius>20";
            save_name+=".png";
            a->SaveAs(save_name.Data());

        }
        if(save==4)
        {
            TString save_name = "Lambda_radius>20";
            save_name+=".png";
            a->SaveAs(save_name.Data());

        }
        if(save==4)
        {
            TString save_name = "mass_squared_kaons_get_purity";
            save_name+=".png";
            a->SaveAs(save_name.Data());
        }

    }
    return func_Gauss_fit3;
}

void tdistribution(TH1D* histo_2,int mode=0)
{
    TH1D* histo_original = (TH1D*)histo_2->Clone();
    TH1D* histo = (TH1D*)histo_2->Clone();

    TF1* func_tdist;
    func_tdist             = new TF1("func_tdist",tdistribution_pdf,-0.4,1.4,4);

    func_tdist->SetParameter(0,100);
    func_tdist->SetParameter(1,0);
    if(mode==0){func_tdist->SetParameter(2,300);}
    if(mode==1)func_tdist->SetParameter(2,200e3);

    histo ->Fit("func_tdist","QMNS","",-0.2,0.1);

    double par0,par1,par2,par3;
    par0 = func_tdist->GetParameter(0);
    par1 = func_tdist->GetParameter(1);
    par2 = func_tdist->GetParameter(2);
    par3 = func_tdist->GetParameter(3);

    func_tdist->SetParameter(0,par0);
    func_tdist->SetParameter(1,par1);
    func_tdist->SetParameter(2,par2);
    func_tdist->SetParameter(3,par3);


    histo ->Fit("func_tdist","QMNS","",-0.2,0.1);

    cout<<"tdistribution"<<endl;
    cout<<par0<<endl;
    cout<<par1<<endl;
    cout<<par2<<endl;
    cout<<par3<<endl;


    TCanvas* a = new TCanvas();
    histo_original->Draw("P E1");
    func_tdist->Draw("same");

    histo->Add(func_tdist,-1);
    histo->Rebin(2);
    TCanvas* b = new TCanvas();
    histo->Draw("P E1");
    func_tdist->Draw("same");

    double sigma = 0.1;
    double mean = 0.243;
    TF1* func_Gauss_fit3 = Gauss(histo, sigma,mean,1,0 ,0.6,0,0);

    TGraph* gr = new TGraph();
    int i_point=0;

    for(double x=-0.2;x<0.6;x+=0.01)
    {
        double t = func_tdist->Eval(x);
        double g = func_Gauss_fit3->Eval(x);
        double ratio = g / t;
        gr->SetPoint(i_point,x,ratio);
        i_point++;
    }
    TCanvas* c = new TCanvas();
    gr->Draw();

    TString save_name = "mass_squared_kaons_get_purity";
    save_name+=".png";
    //b->SaveAs(save_name.Data());

    save_name = "purity_distr_for_kaons";
    save_name+=".png";
    //c->SaveAs(save_name.Data());


}

vector<double> tdistribution_double_fit(TH1D* histo_2,int mode=0)
{
    TH1D* histo_original = (TH1D*)histo_2->Clone();
    TH1D* histo = (TH1D*)histo_2->Clone();

    TF1* func_tdist;
    func_tdist             = new TF1("func_tdist_double",double_fit_tdistribution_pdf,-0.4,1.4,8);

    func_tdist->SetParameter(0,1);
    func_tdist->SetParameter(1,0.25);
    func_tdist->SetParameter(2,10e3);

    func_tdist->SetParameter(4,1);
    func_tdist->SetParameter(5,0.9);
    func_tdist->SetParameter(6,5e3);

    histo ->Fit("func_tdist_double","QMNS","",0.2,1.4);

    for(int i=0;i<8;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }

    histo ->Fit("func_tdist_double","QMNS","",0.2,1.4);

    for(int i=0;i<8;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }

    histo ->Fit("func_tdist_double","QMNS","",0.2,1.4);

    TCanvas* a = new TCanvas();
    histo_original->Draw("P E1");
    func_tdist->Draw("same");

    vector<double> vec;
    cout<<"double fit: "<<endl;
    for(int i=0;i<8;i++)
    {
        cout<<func_tdist->GetParameter(i)<<endl;
        vec.push_back(func_tdist->GetParameter(i));
    }
    return vec;
}

vector<double> tdistribution_triple_fit(TH1D* histo,vector<double> vec,TCanvas* a,TF1* func_tdist)
{
    TH1D* histo_original = (TH1D*)histo->Clone();

    //TF1* func_tdist;
    //func_tdist             = new TF1("func_tdist_triple",triple_fit_tdistribution_pdf,-0.4,1.4,12);

    for(int i=0;i<8;i++)
    {
        func_tdist->SetParameter(i+4,vec[i]);
    }

    //first fix 2nd and 3rd peak to zero
    func_tdist->FixParameter(6,0);
    func_tdist->FixParameter(10,0);
    
    func_tdist->SetParameter(0,149);
    func_tdist->SetParameter(1,0);
    func_tdist->SetParameter(2,300);
    //Fit first peak
    histo ->Fit("func_tdist_triple","QMNS","",-0.2,0.1);
    for(int i=0;i<12;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }
    histo ->Fit("func_tdist_triple","QMNS","",-0.2,0.1);
    for(int i=0;i<12;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }

    //------------------------------------------------------------

    func_tdist->ReleaseParameter(6);
    func_tdist->ReleaseParameter(10);
   
    for(int i=0;i<4;i++)
    {
        func_tdist->FixParameter(i,func_tdist->GetParameter(i));
    }

    for(int i=0;i<8;i++)
    {
        func_tdist->SetParameter(i+4,vec[i]);
    }
    

    //func_tdist->SetParameter(0+4,1);
   // func_tdist->SetParameter(1+4,-5);
    //func_tdist->SetParameter(2+4,10700);

    histo ->Fit("func_tdist_triple","QMNS","",0.2,1.4);
    for(int i=0;i<12;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }
    histo ->Fit("func_tdist_triple","QMNS","",0.2,1.4);


    //histo ->Fit("func_tdist_triple","QMNS","",-0.2,1.4);
    for(int i=0;i<4;i++)
    {
        func_tdist->ReleaseParameter(i);
    }
    histo ->Fit("func_tdist_triple","QMNS","",-0.4,1.4);
    for(int i=0;i<12;i++)
    {
        func_tdist->SetParameter(i,func_tdist->GetParameter(i));
    }
    histo ->Fit("func_tdist_triple","QMNS","",-0.4,1.4);

    vector<double> res;

    for(int i=0;i<12;i++)
    {
        cout<<"param: "<<func_tdist->GetParameter(i)<<endl;
        res.push_back(func_tdist->GetParameter(i));
    }

    //style
    //TString title = "m^{2} spectra of Kaons";
    TString title = "";
    TString xlabel = "m^{2} [(GeV/c^{2})^{2}]";
    TString ylabel = "counts";

    histo_original->SetTitle(title.Data());
    histo_original->GetYaxis()->SetTitle(ylabel.Data());
    histo_original->GetXaxis()->SetTitle(xlabel.Data());
    histo_original->GetXaxis()->SetTitleSize(0.050);
    //histo_original->GetYaxis()->SetMaxDigits(3);
    histo_original->GetYaxis()->SetTitleSize(0.050);
    histo_original->GetXaxis()->SetTitleOffset(0.9);
    histo_original->GetYaxis()->SetTitleOffset(0.9);
    histo_original->GetXaxis()->CenterTitle();
    histo_original->GetYaxis()->CenterTitle();
    histo_original->GetXaxis()->SetRangeUser(-0.38,1.38);
    histo_original->GetYaxis()->SetRangeUser(400,3e5);
    histo_original->SetMarkerSize(10);

    gStyle->SetLegendBorderSize(0);
    //SetRootGraphicStyle();
    gStyle->SetOptStat(0);


    //TCanvas* a = new TCanvas();
    a->SetTicks(1,1);
    a->SetLogy();
    histo_original->Draw("P E1");
    func_tdist->SetLineColor(kRed);
    func_tdist ->SetFillStyle(3001);
    func_tdist ->SetFillColor(kRed);
    func_tdist->Draw("same fc");
    func_tdist->Draw("same fc");

    return res;

}



void Fit_particles()
{
    TFile* file3 = TFile::Open("/misc/alidata120/alice_u/schlichtmann/out/Merge_S_search_V17_high_stat.root");

    TString name ="mass_squared_antip_type3";

    TH1D* histo_antip_type3 = (TH1D*)file3->Get(name.Data());

    TString xaxislabel = "mass squared [GeV^{2}]";
    TString yaxislabel = "counts";

    double sigma = 0.2;
    double mean = 0.;
    Gauss(histo_antip_type3, sigma,mean,1,-0.2,0.2);

    name ="mass_squared_kaons_type3";
    TH1D* histo_kaonen_type3= (TH1D*)file3->Get(name.Data());

    tdistribution(histo_kaonen_type3,1);


    //--------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------

    TFile* file4 = TFile::Open("/misc/alidata120/alice_u/schlichtmann/out/Merge_S_search_V18_high_stat.root");
    name ="mass_squared_all_kaons_type2";
    TH1D* histo_kaonen_type2 = (TH1D*)file4->Get(name.Data());
    tdistribution(histo_kaonen_type2);

    vector<double> vec;
    vec.resize(8);
    vec = tdistribution_double_fit(histo_kaonen_type2);

    vector<double> res_3_fit;
    res_3_fit.resize(12);
    TCanvas* can = new TCanvas();

    TF1* triple_func;
    triple_func             = new TF1("func_tdist_triple",triple_fit_tdistribution_pdf,-0.4,1.4,12);
    triple_func->SetNpx(300);

    res_3_fit = tdistribution_triple_fit(histo_kaonen_type2,vec,can,triple_func);

    TF1* func_tdist_1;
    func_tdist_1             = new TF1("func_tdist_1",tdistribution_pdf,-0.4,1.4,4);
    func_tdist_1->SetNpx(300);
    for(int i=0;i<4;i++)
    {
        func_tdist_1->SetParameter(i,res_3_fit[i]);
    }
    func_tdist_1->SetLineColor(kBlue);
    func_tdist_1 ->SetFillStyle(3001);
    func_tdist_1 ->SetFillColor(kBlue);
    func_tdist_1->Draw("same fc");

    TF1* func_tdist_2;
    func_tdist_2             = new TF1("func_tdist_2",tdistribution_pdf,-0.4,1.4,4);
    func_tdist_2->SetNpx(300);
    for(int i=0;i<4;i++)
    {
        func_tdist_2->SetParameter(i,res_3_fit[i+8]);
    }
    func_tdist_2->SetLineColor(kGreen);
    func_tdist_2 ->SetFillStyle(3001);
    func_tdist_2 ->SetFillColor(kGreen);
    func_tdist_2->Draw("same fc");
    func_tdist_2->Draw("fcsame");


    TF1* func_tdist_3;
    func_tdist_3             = new TF1("func_tdist_3",tdistribution_pdf,-0.4,1.4,4);
    func_tdist_3->SetNpx(300);
    for(int i=0;i<4;i++)
    {
        func_tdist_3->SetParameter(i,res_3_fit[i+4]);
    }
    func_tdist_3->SetLineColor(kMagenta);
    func_tdist_3 ->SetFillStyle(3001);
    func_tdist_3 ->SetFillColor(kMagenta);
    func_tdist_3->Draw("same fc");

    auto legend = new TLegend(0.6,0.5,0.85,0.85);
    //legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetHeader("dE/dx selected Kaons",""); // option "C" allows to center the header
    legend->AddEntry("func_tdist_triple","all","l");
    legend->AddEntry("func_tdist_1","#pi","l");
    legend->AddEntry("func_tdist_2","K^{+}","l");
    legend->AddEntry("func_tdist_3","p","l");
    //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
    //legend->AddEntry("gr","Graph with error bars","l");
    legend->Draw();

    can->SaveAs("m_squared_Kaons_triple_fit.png");


    int i_point=0;
    TGraph* gr = new TGraph();
    for(double x=-0.2;x<0.6;x+=0.01)
    {
        double signal = func_tdist_2->Eval(x);
        double back1 = func_tdist_1->Eval(x);
        double back2 = func_tdist_3->Eval(x);
        double ratio = signal / (back1+back2);
        gr->SetPoint(i_point,x,ratio);
        i_point++;
    }
    TCanvas* c = new TCanvas();
    c->SetTicks(1,1);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle("signal/background");
    gr->GetXaxis()->SetTitle("m^{2} [(GeV/c^{2})^{2}]");
    gr->GetXaxis()->SetTitleSize(0.050);
    //gr->GetYaxis()->SetMaxDigits(3);
    gr->GetYaxis()->SetTitleSize(0.050);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetTitleOffset(0.9);
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
    gr->SetMarkerSize(10);

    gr->Draw();

    

    


    //set logy
    //schattieren
    //histogrammpunkte groeser
    //in Legende griechische Buchstaben
    //set npx
    //bis 1.38
    
}

