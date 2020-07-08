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

Double_t polynom(Double_t* x_val, Double_t* par)
{
    double x,par0,par1,par2;
    x=x_val[0];
    par0=par[0];
    par1=par[1];
    par2=par[2];
    return par0*x*x+par1*x+par2;
}

void print_vec(vector<double> vec)
{
    for(int i=0;i<vec.size();i++)
    {
        cout<<vec[i]<<endl;
    }
}

int find_max_element_position(vector<double> vec)
{
    double max_element=-1;
    int max_position = -1;

    for(int i=0;i<vec.size();i++)
    {
        if(vec[i]>max_element)
        {
            max_element = vec[i];
            max_position = i;
        }
    }
    cout<<"max_postion: "<<max_position<<" max_element: "<<max_element<<endl;
    return max_position;
}

//void calculate_sig(TH1D* blubber, TF* func_gauss)
//{
//
//}

double Gauss(TH1D* blubber, double sigma,bool draw)
{
    //cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    //cout<<"Gaus started"<<endl;

    TF1* func_Gauss_fit3;
    func_Gauss_fit3             = new TF1("func_Gauss_fit3",GaussFitFunc_plus2ndorderpolynom,0,1.5,6);

    int maxbin=blubber->GetMaximumBin();
    double amplitude=blubber->GetBinContent(maxbin);
    double mean=blubber->GetBinCenter(maxbin);
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
    func_Gauss_fit3 ->SetParameter(3,0.0);
    func_Gauss_fit3 ->SetParameter(4,0.0);
    func_Gauss_fit3 ->SetParameter(5,600.);

    //Fit
    blubber ->Fit("func_Gauss_fit3","QMN","",mean-5.0*sigma,mean+5.0*sigma);

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
    blubber ->Fit("func_Gauss_fit3","QMN","",mean-5.0*sigma,mean+5.0*sigma);

    //Parameter auslesen
    amplitude = func_Gauss_fit3 ->GetParameter(0);
    mean      = func_Gauss_fit3 ->GetParameter(1);
    sigma     = fabs(func_Gauss_fit3 ->GetParameter(2));

    
    /*
    cout<<amplitude<<endl;
    cout<<mean<<endl;
    cout<<sigma<<endl;
    cout<<par3<<endl;
    cout<<par4<<endl;
    cout<<par5<<endl;
    */
    //etwas bessere Fitwerte
    if(draw)
    {
        TCanvas* a = new TCanvas();
        blubber->Draw("P E1");
    }

    func_Gauss_fit3 ->SetLineColor(kRed);
    func_Gauss_fit3 ->SetLineStyle(1);
    func_Gauss_fit3 ->SetRange(mean-4.0*sigma,mean+4.0*sigma);
    //func_Gauss_fit3 ->SetRange(1.3,1.36);
    if(draw) {func_Gauss_fit3 ->Draw("same");}

    //calculate Sig
    TF1* Gauss_only;
    Gauss_only  = new TF1("Gauss_only",GaussFitFunc,0,1.5,3);
    Gauss_only -> SetParameter(0,amplitude);
    Gauss_only -> SetParameter(1,mean);
    Gauss_only -> SetParameter(2,sigma);
    Gauss_only ->SetLineColor(kBlack);
    if(draw) {Gauss_only->Draw("same"); }

    par3   =  func_Gauss_fit3 ->GetParameter(3);
    par4   =  func_Gauss_fit3 ->GetParameter(4);
    par5   =  func_Gauss_fit3 ->GetParameter(5);

    TF1* polynomial_only;
    polynomial_only = new TF1("polynomial_only",polynom,0,1.5,3);
    polynomial_only->SetParameter(0,par3);
    polynomial_only->SetParameter(1,par4);
    polynomial_only->SetParameter(2,par5);
    polynomial_only ->SetLineColor(kBlue);
    if(draw){ polynomial_only->Draw("same");  }

    double signal = Gauss_only->Integral(mean-sigma,mean+sigma);
    double untergrund  = polynomial_only->Integral(mean-sigma,mean+sigma);
    double Sig = signal/sqrt(signal+untergrund);
    printf("signal: %f, untergrund: %f, Sig: %f \n",signal, untergrund, Sig);
    cout<<""<<endl;
    return Sig;
}


void Gauss_fit()
{
    TString inputdir = "/misc/alidata120/alice_u/schlichtmann/out/";
    inputdir += "Pb_V3_list.txt_out.root";
    TFile* file = TFile::Open(inputdir.Data());

    vector<double> vec_Sig;
    float arr_distance_prim_sec[4]{1.,2.5,4.,6.};
    float arr_distance_daughter_particles[4]{0.05,0.07,0.1,0.2};
    float arr_dca_daughter_prim[4]{0.3,0.5,0.7,1.0};
    //for(int i=0;i<64;i++)
    //{
    int i=0;
    double sigma = 0.001;
    for(int i=0;i<64;i++)
    {
        cout<<i<<endl;
        TString name = "radius >  ";
        name+=arr_distance_prim_sec[i/16];
        name+=" dcaV0 < ";
        name+=arr_distance_daughter_particles[(i%16)/4];
        name+=" dca_daughter_prim> ";
        name+=arr_dca_daughter_prim[i%4];
        TH1D* histo = (TH1D*)file->Get(name.Data());

        double Sig = Gauss(histo, sigma,0);
        vec_Sig.push_back(Sig);
       
    }

    int max_pos = find_max_element_position(vec_Sig);

    i=max_pos;
    TString name = "radius >  ";
    name+=arr_distance_prim_sec[i/16];
    name+=" dcaV0 < ";
    name+=arr_distance_daughter_particles[(i%16)/4];
    name+=" dca_daughter_prim> ";
    name+=arr_dca_daughter_prim[i%4];
    TH1D* histo = (TH1D*)file->Get(name.Data());

    double Sig = Gauss(histo, sigma,1);
    

    //sort(vec_Sig.begin(),vec_Sig.end());
    //print_vec(vec_Sig);

    

    



}