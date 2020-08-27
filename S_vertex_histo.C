
void S_vertex_histo()
{
    TFile* inputfile = TFile::Open("./Histos.root");
    TH1D* histo_s_vertex_radius = (TH1D*)inputfile->Get("histo S vertex radius");
    TH1D* histo_S_mass_r_larger_10 = (TH1D*)inputfile->Get("histo S mass r>10");
    TH1D* histo_S_mass_r_larger_20 = (TH1D*)inputfile->Get("histo S mass r>20");


    TH2D* histo_S_x_and_y = (TH2D*)inputfile->Get("histo S x and y");


    TCanvas* a = new TCanvas();
    TCanvas* b = new TCanvas();
    TCanvas* c = new TCanvas();
    TCanvas* d = new TCanvas();

    a->cd();
    histo_s_vertex_radius->Draw();


    b->cd();
    histo_S_x_and_y->Draw();

    c->cd();
    histo_S_mass_r_larger_10->Draw();

    d->cd();
    histo_S_mass_r_larger_20->Draw();

}
    