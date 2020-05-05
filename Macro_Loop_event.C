
#include "Dark_Matter_Read.h"


void Macro_Loop_event()
{
    printf("Macro_Loop_event started \n");
    TFile* outputfile = new TFile("Histos.root","RECREATE");
    

    Dark_Matter_Read* DM_Read = new Dark_Matter_Read();
    DM_Read ->Init_tree("file_list.txt");

    Long64_t numentries = DM_Read -> getnumberentries();
    printf("num entries: %i \n",numentries);

    TH1D* histo_invariantmass_lambda = new TH1D("histo inv mass lambda","histo inv mass lambda",50*2,1.1,1.13);
    TH1D* histo_invariantmass_K0 = new TH1D("histo inv mass K0","histo inv mass K0",50*3,0.4,0.6);

    TH1D* histo_lambda_vertex_radius = new TH1D("histo lambda vertex radius ","histo lambda vertex radius ",50,0,200);
    TH1D* histo_K0_vertex_radius = new TH1D("histo K0 vertex radius ","histo K0 vertex radius ",50,0,200);
    

    int binning2D = 200;
    TH2D* histo_lambda_x_and_y = new TH2D("histo lambda x and y ","histo lambda x and y ",binning2D,-200,200,binning2D,-200,200);
    TH2D* histo_K0_x_and_y = new TH2D("histo K0 x and y ","histo K0 x and y ",binning2D,-200,200,binning2D,-200,200);
   

    TH1D* histo_S_vertex_radius = new TH1D("histo S vertex radius ","histo S vertex radius ",200,0,200);
    
    TH2D* histo_S_x_and_y = new TH2D("histo S x and y ","histo S x and y ",binning2D,-200,200,binning2D,-200,200);

    TH1D* histo_S_mass_r_larger_10 = new TH1D("histo S mass r>10","histo S mass r>10",200,0,20);
    TH1D* histo_S_mass_r_larger_20 = new TH1D("histo S mass r>20","histo S mass r>20",200,0,20);

    TH1D* histo_num_tracks = new TH1D("histo number of tracks","histo number of tracks",200,0,1000);
    TH1D* histo_vertex_z_pos = new TH1D("histo vertex z position","histo vertex z position",100,-30,30);

    vector<TH1D*> histos_1D;
    vector<TH2D*> histos_2D;

    histos_1D.push_back(histo_invariantmass_lambda);
    histos_1D.push_back(histo_invariantmass_K0);
    histos_1D.push_back(histo_lambda_vertex_radius);
    histos_1D.push_back(histo_K0_vertex_radius);
    histos_1D.push_back(histo_S_vertex_radius);
    histos_1D.push_back(histo_S_mass_r_larger_10);
    histos_1D.push_back(histo_S_mass_r_larger_20);
    histos_1D.push_back(histo_num_tracks);      //7
    histos_1D.push_back(histo_vertex_z_pos);      //87

    histos_2D.push_back(histo_lambda_x_and_y);
    histos_2D.push_back(histo_K0_x_and_y);
    histos_2D.push_back(histo_S_x_and_y);
    double event_counter=0;


    //--------------------------------------------------------------------------------------------------
    //for(int i =0 ; i<numentries; i++)
    for(int i =0 ; i<1000; i++)
    {
        //printf("eventcounter: %f",event_counter);
        DM_Read ->Loop_event(i,histos_1D,histos_2D,event_counter);
    }



    printf("number of entries: %d, number of events: %f", numentries, event_counter);

   /* TCanvas* a = new TCanvas();
    TCanvas* b = new TCanvas();
    TCanvas* c = new TCanvas();
    TCanvas* d = new TCanvas();
    TCanvas* e = new TCanvas();
    TCanvas* f = new TCanvas();
    TCanvas* g = new TCanvas();
    TCanvas* h = new TCanvas();
    TCanvas* i = new TCanvas();
    TCanvas* j = new TCanvas(); */
    TCanvas* k = new TCanvas();
    TCanvas* l = new TCanvas();

   /*
    a->cd();
    //histo_invariantmass_lambda->Draw();
   // histo_invariantmass_K0->Draw();
    histo_lambda_vertex_radius->Draw();
    b->cd();
    histo_K0_vertex_radius->Draw();

    c->cd();
    histo_lambda_x_and_y->Draw("colz");
    d->cd();
    histo_K0_x_and_y->Draw("colz");
    e->cd();
    histo_S_vertex_radius->Draw();
    f->cd();
    histo_S_x_and_y->Draw("colz");
    g->cd();
    histo_invariantmass_lambda->Draw();
    h->cd();
    histo_invariantmass_K0->Draw();

    i->cd();
    histo_S_mass_r_larger_10->Draw();

    j->cd();
    histo_S_mass_r_larger_20->Draw();
    */
    k->cd();
    histo_num_tracks->Draw();

    l->cd();
    histo_vertex_z_pos->Draw();


    /*
    outputfile ->cd();
    histo_invariantmass_lambda   ->Write();
    histo_invariantmass_K0       ->Write();
    histo_lambda_vertex_radius   ->Write();
    histo_lambda_x_and_y         ->Write();
    histo_K0_x_and_y             ->Write();
    histo_S_vertex_radius        ->Write();
    histo_S_x_and_y              ->Write();
    histo_S_mass_r_larger_10     ->Write();
    histo_S_mass_r_larger_20     ->Write();
    histo_num_tracks             ->Write();
    histo_vertex_z_pos           ->Write();

    outputfile->Close();
    */

}