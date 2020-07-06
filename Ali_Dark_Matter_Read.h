
#ifndef __ALI_DARK_MATTER_READ_H__
#define __ALI_DARK_MATTER_READ_H__

using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TObject.h"


#include<TMath.h>
#include<cmath> 

//for generator
#include<vector>
#include<algorithm>
#include<iterator>
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TChain.h"
#include "TTree.h"
//#include "AliHelix.h"
#include "AliTRDpadPlane.h"


#include "Ali_AS_Event_old.h"
#include "Ali_AS_EventLinkDef.h"

ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_V0)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_Event)
ClassImp(Ali_AS_DM_particle);



//----------------------------------------------------------------------------------------
class Ali_Dark_Matter_Read
{
private:

    Ali_AS_Event*     AS_Event;
    Ali_AS_Track*     AS_Track;
    Ali_AS_V0*        AS_V0;
    Ali_AS_Tracklet*  AS_Tracklet;
    Ali_AS_offline_Tracklet*  AS_offline_Tracklet;
    Ali_AS_TRD_digit* AS_Digit;
    Ali_AS_DM_particle* AS_DM_particle;
    Ali_AS_Track*     AS_DM_Track;

    TTree       *Tree_AS_DM_particle;

    TGraph* gr = new TGraph();
    TGraph* gr2 = new TGraph();

    int counter_path_trackA = 0;
    int counter_path_trackB = 0;
    int counter_path_trackC = 0;
    int counter_path_trackA_negativ = 0;
    int counter_path_trackB_negativ = 0;
    int counter_path_trackC_negativ = 0;

    int counter_lambdas = 0;
    int counter_anti_lambdas = 0;

    TH2D* mass_squared_vs_charge_dot_momentum = new TH2D("mass_squared_vs_charge_dot_momentum","mass_squared_vs_charge_dot_momentum",
                                                         200,-8,8,200,-0.5,0.5);


    //histos from Macro
    TH1D* histo_invariantmass_lambda = new TH1D("histo inv mass lambda","histo inv mass lambda",50*2,1.1,1.13);
    TH1D* histo_invariantmass_anti_lambda = new TH1D("histo inv mass anti lambda","histo inv mass anti lambda",50*2,1.0,1.5);

    TH1D* histo_invariantmass_K0 = new TH1D("histo inv mass K0","histo inv mass K0",50*3,0.4,0.6);

    TH1D* histo_lambda_vertex_radius = new TH1D("histo lambda vertex radius ","histo lambda vertex radius ",50,0,200);
    TH1D* histo_K0_vertex_radius = new TH1D("histo K0 vertex radius ","histo K0 vertex radius ",50,0,200);

    TH1D* mass_squared_kaons = new TH1D("mass_squared_kaons", "mass_squared_kaons",100,-0.1,0.4);
    TH1D* mass_squared_kaons_and_background = new TH1D("mass_squared_kaons_and_background", "mass_squared_kaons_and_background",100,-0.1,0.4);

    TH1D* histo_counter = new TH1D("histo counter 1.5: S-vertices with two pions; 2.5: S-vertices that fulfill all cuts",
                                   "histo counter 1.5: S-vertices with two pions; 2.5: S-vertices that fulfill all cuts",10,0,10);

    int binning2D = 400;
    TH2D* histo_lambda_x_and_y = new TH2D("histo lambda x and y ","histo lambda x and y ",binning2D,-200,200,binning2D,-200,200);
    TH2D* histo_K0_x_and_y = new TH2D("histo K0 x and y ","histo K0 x and y ",binning2D,-200,200,binning2D,-200,200);

    TH2D* histo_V0_with_K_plus_x_and_y = new TH2D("histo V0 with K + x and y","histo V0 with K + x and y",binning2D,-100,100,binning2D,-100,100);
   

    TH1D* histo_S_vertex_radius = new TH1D("histo S vertex radius ","histo S vertex radius ",200,0,200);
    
    TH2D* histo_S_x_and_y = new TH2D("histo S x and y ","histo S x and y ",binning2D,-100,100,binning2D,-100,100);

    TH1D* histo_S_mass_r_larger_10 = new TH1D("histo S mass r>10","histo S mass r>10",200,0,20);
    TH1D* histo_S_mass_r_larger_20 = new TH1D("histo S mass r>20","histo S mass r>20",200,0,20);

    TH1D* histo_num_tracks = new TH1D("histo number of tracks","histo number of tracks",200,0,1000);
    TH1D* histo_vertex_z_pos = new TH1D("histo vertex z position","histo vertex z position",100,-30,30);

    TH1D* histo_V0_with_K_plus_radius = new TH1D("histo V0 with K plus radius","histo V0 with K plus radius",100,0,200);

    vector<TH1D*> histos_1D;
    vector<TH2D*> histos_2D;
    //--------------------------------------------------------------

    //counters----------------------------------------------------
    int event_counter=0;
    int counter_of_2pions_close_to_S_vertex=0;
    int counter_of_V0s_of_antiproton_and_K_plus_and_another_K=0;
    int counter_of_S_vertices_without_pions=0;
    int counter_1pion_close=0;
    int counter_V0s_antiproton_and_K_plus=0;
    int counter_correct_S_vertices=0;

    int counter_vertices_antip_K_plus_K_plus=0;
    int counter_vertices_antip_K_plus_K_plus_r_larger_5=0;
    int counter_vertices_antip_K_plus_K_plus_r_larger_5_and_dot_product=0;
    int counter_skipped_track=0;

    vector<int> counters;

    //-----------------------------------------------------------------


    TH2D* dEdx_vs_charge_dot_momentum;
    //TH2D* mass_squared_vs_charge_dot_momentum_kaons;

    TH1D* mass_squared_no_pions = new TH1D("mass squared no pions", "mass squared no pions",100,-0.1,0.1);

    TH1D* histo_reference_vertex_radius_3_pionen = new TH1D("histo_reference_vertex_radius_3_pionen","histo_reference_vertex_radius_3_pionen",50*5,0,200);
    TH1D* histo_reference_vertex_radius_4_pionen = new TH1D("histo_reference_vertex_radius_4_pionen","histo_reference_vertex_radius_4_pionen",50*5,0,200);

    TH2D* histo_reference_x_and_y_3_pionen = new TH2D("histo_reference_x_and_y_3_pionen ","histo_reference_x_and_y_3_pionen ",200,-100,100,200,-100,100);
    TH2D* histo_reference_x_and_y_4_pionen = new TH2D("histo_reference_x_and_y_4_pionen ","histo_reference_x_and_y_4_pionen ",200,-100,100,200,-100,100);

    TFile* outputfile;
    TFile* outputfile_trkl;

    TFile* outputfile_histos;
    TFile* output_histos;

    //TFile* ofile = new TFile("ntuple.root","RECREATE");
    //TFile* ofile2 = new TFile("mass_squared_and_dEdx.root","RECREATE");

    TFile* calibration_params;

    TH1D* histo_invariantmass_gamma = new TH1D("histo inv mass gamma","histo inv mass gamma",50*2*4,-0.1,1.);
    TH1D* histo_invariantmass_electron_plus = new TH1D("histo inv mass squared electron plus","histo inv mass squared electron plus",50*2*4,-0.2,1.3);
    TH1D* histo_invariantmass_electron_minus = new TH1D("histo inv mass electron minus","histo inv mass electron minus",50*2,2e-4,8e-4);

    TH1D* histo_invariantmass_xi_minus_baryon = new TH1D("invariant mass xi-","invariant mass xi-",200,1.2,1.4);
    TH1D* histo_invariantmass_xi_plus_baryon = new TH1D("invariant mass xi+","invariant mass xi+",200,1.2,1.4);

    TH1D* histo_invariantmass_omega_minus_baryon = new TH1D("invariant mass omega-","invariant mass omega-",200,1.5,1.8);
    TH1D* histo_invariantmass_omega_plus_baryon = new TH1D("invariant mass omega+","invariant mass omega+",200,1.5,1.8);

    vector<TH1D*> vec_histo_omega_minus;
    vector<TH1D*> vec_histo_omega_plus;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;


    Long64_t Event_active = 0;

    TChain* input_SE;
    TString JPsi_TREE   = "Tree_AS_Event";
    TString JPsi_BRANCH = "Tree_AS_Event_branch";
    Long64_t file_entries_total;

    Float_t digit_pos[3];


    TString HistName;
    char NoP[50];

    double counter =0;

    //TNtuple* tpl = new TNtuple ("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");


public:
    Ali_Dark_Matter_Read(){}
    Ali_Dark_Matter_Read(TString list);
    //Ali_Dark_Matter_Read(TString list);
    ~Ali_Dark_Matter_Read();
    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t event );
    Long64_t getnumberentries(){return file_entries_total;};
    void copy_track_params(Ali_AS_Track* track_in, Ali_AS_Track* track_out);
    void Save();

    ClassDef(Ali_Dark_Matter_Read, 1)
};
//----------------------------------------------------------------------------------------



#endif // __ALI_DARK_MATTER_READ_H__




