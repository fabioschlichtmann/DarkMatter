
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
#include "TProfile.h"

#include "TChain.h"
#include "TTree.h"
//#include "AliHelix.h"
#include "AliTRDpadPlane.h"


#include "Ali_AS_Event_V3.h"
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
    int counter_NUCLEV = 0;

    const Float_t mass_proton = 0.93827208816 ;  //in GeV?
    const Float_t mass_pion = 0.139657061 ;  //in GeV?
    const Float_t mass_electron = 0.510998950 * 1e-3 ;  //in GeV?
    const double  mass_K = 0.493677 ;  //in GeV?

    double mass_K0 = 0.493677;
    double mass_Lambda = 1.1155683;
    double mass_neutron = 0.939565;
    double S_mass = -1;

    TH2D* mass_squared_vs_charge_dot_momentum = new TH2D("mass_squared_vs_charge_dot_momentum","mass_squared_vs_charge_dot_momentum",
                                                         200,-8,8,200,-0.5,0.5);


    int counter_DM_saved = 0;

    //histos from Macro
    TH1D* histo_invariantmass_lambda = new TH1D("histo inv mass lambda","histo inv mass lambda",50*2,1.1,1.13);

    vector<TH1D*> vec_histo_invariantmass_Lambda;
    vector<TH1D*> vec_histo_invariantmass_K0;
    vector<TH1D*> vec_histo_radius_variation_invariantmass_K0;
    vector<TH1D*> vec_histo_radius_variation_invariantmass_Lambda;
    

    TH1D* histo_invariantmass_anti_lambda = new TH1D("histo inv mass anti lambda","histo inv mass anti lambda",50*2,1.1,1.13);

    TH1D* histo_invariantmass_K0 = new TH1D("histo inv mass K0","histo inv mass K0",50*3,0.4,0.6);

    //TH1D* histo_invariantmass_K0_type3 = new TH1D("histo_invariantmass_K0_type3","histo_invariantmass_K0_type3",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type3_with_cuts_on_antip_and_K =
        new TH1D("histo_invariantmass_K0_type3_with_cuts_on_antip_and_K","histo_invariantmass_K0_type3_with_cuts_on_antip_and_K",50*3,0.4,0.6);

    TH1D* histo_invariantmass_K0_type3_with_cut_on_antip =
        new TH1D("histo_invariantmass_K0_type3_with_cut_on_antip","histo_invariantmass_K0_type3_with_cut_on_antip",50*3,0.4,0.6);

    TH1D* histo_invariantmass_K0_type3_with_cuts_on_K =
        new TH1D("histo_invariantmass_K0_type3_with_cuts_on_K","histo_invariantmass_K0_type3_with_cuts_on_K",50*3,0.4,0.6);



    TH1D* histo_lambda_vertex_radius = new TH1D("histo lambda vertex radius ","histo lambda vertex radius ",50,0,200);
    TH1D* histo_K0_vertex_radius = new TH1D("histo K0 vertex radius ","histo K0 vertex radius ",50,0,200);

    TH1D* mass_squared_kaons = new TH1D("mass_squared_kaons", "mass_squared_kaons",100,-0.1,0.4);
    TH1D* mass_squared_kaons_type3 = new TH1D("mass_squared_kaons_type3", "mass_squared_kaons_type3",100,-1.1,1.4);
    TH1D* mass_squared_kaons_type3_with_cut_on_antip = new TH1D("mass_squared_kaons_type3_with_cut_on_antip", "mass_squared_kaons_type3_with_cut_on_antip",100,-0.4,1.4);
    TH1D* mass_squared_kaons_type3_with_cut_on_invmass_K0 = new TH1D("mass_squared_kaons_type3_with_cut_on_invmass_K0", "mass_squared_kaons_type3_with_cut_on_invmass_K0",100,-0.4,1.4);
    TH1D* mass_squared_antip_type3 = new TH1D("mass_squared_antip_type3", "mass_squared_antip_type3",100,-1.1,1.4);
    TH1D* mass_squared_antip_type3_with_cut_on_K = new TH1D("mass_squared_antip_type3_with_cut_on_K", "mass_squared_antip_type3_with_cut_on_K",100,-0.4,1.4);
    TH1D* mass_squared_antip_type3_with_cut_on_invmass_K0 = new TH1D("mass_squared_antip_type3_with_cut_on_invmass_K0", "mass_squared_antip_type3_with_cut_on_invmass_K0",100,-0.4,1.4);
    TH1D* mass_squared_kaons_and_background = new TH1D("mass_squared_kaons_and_background", "mass_squared_kaons_and_background",100,-0.1,0.4);

    //type2
    TH1D* mass_squared_all_kaons_type2 = new TH1D("mass_squared_all_kaons_type2", "mass_squared_all_kaons_type2",100,-0.4,1.4);
    TH1D* mass_squared_all_kaons_overlap_cuts_type2 = new TH1D("mass_squared_all_kaons_overlap_cuts_type2", "mass_squared_all_kaons_overlap_cuts_type2",100,-0.4,1.4);
    TH1D* histo_m_squared_kaon_no_other_cuts_type2 = new TH1D("histo_m_squared_kaon_no_other_cuts_type2", "histo_m_squared_kaon_no_other_cuts_type2",100,-0.1,0.4);
    TH1D* m_squared_anti_proton_type2= new TH1D("m_squared_anti_proton_type2", "m_squared_anti_proton_type2",100,-0.1,1.5);
    TH1D* histo_m_squared_kaon_with_all_other_cuts_type2 = new TH1D("histo_m_squared_kaon_with_all_other_cuts_type2", "histo_m_squared_kaon_with_all_other_cuts_type2",100,-0.1,0.4);
    TH1D* histo_m_squared_kaon_with_some_other_cuts_type2 = new TH1D("histo_m_squared_kaon_with_some_other_cuts_type2", "histo_m_squared_kaon_with_some_other_cuts_type2",100,-0.1,0.4);

    vector<TH1D*> vec_mass_squared_dEdx_selected;


    TH1D* histo_invariant_mass_kaon_no_other_cuts = new TH1D("histo_invariant_mass_kaon_no_other_cuts", "histo_invariant_mass_kaon_no_other_cuts",100,-0.1,0.4);
    TH1D* histo_invariant_mass_kaon_with_other_cuts = new TH1D("histo_invariant_mass_kaon_with_other_cuts", "histo_invariant_mass_kaon_with_other_cuts",100,-0.1,0.4);

    TH1D* histo_counter = new TH1D("histo counter 1.5: S-vertices with two pions; 2.5: S-vertices that fulfill all cuts",
                                   "histo counter 1.5: S-vertices with two pions; 2.5: S-vertices that fulfill all cuts",10,0,10);

    TH1D* histo_counter_overlap = new TH1D("histo_counter_overlap","histo_counter_overlap",10,0.5,10.5);

    TH1D* histo_counter_S = new TH1D("histo_counter_S","histo_counter_S",200,0.5,200.5);
    TH1D* histo_counter_type5 = new TH1D("histo_counter_type5","histo_counter_type5",40,0.5,40.5);

    TH1D* histo_numberDMs = new TH1D("histo_number_DMs","histo_number_DMs",10,0,10);

    int binning2D = 400;
    TH2D* histo_lambda_x_and_y = new TH2D("histo lambda x and y ","histo lambda x and y ",binning2D,-200,200,binning2D,-200,200);
    TH2D* histo_K0_x_and_y = new TH2D("histo K0 x and y ","histo K0 x and y ",binning2D,-200,200,binning2D,-200,200);

    TH2D* histo_V0_with_K_plus_x_and_y = new TH2D("histo V0 with K + x and y","histo V0 with K + x and y",binning2D,-100,100,binning2D,-100,100);
   

    TH1D* histo_S_vertex_radius = new TH1D("histo S vertex radius ","histo S vertex radius ",200,0,200);
    
    TH2D* histo_S_x_and_y = new TH2D("histo S x and y ","histo S x and y ",binning2D,-100,100,binning2D,-100,100);

    TH2D* histo_eta = new TH2D("histo_eta","histo_eta",100,-1,1,100,-1,1);
    TH2D* histo_eta_larger_range = new TH2D("histo_eta_-3_3","histo_eta-3_3",300,-3,3,300,-3,3);
    TH2D* histo_eta_z_cut = new TH2D("histo_eta_z_cut","histo_eta_z_cut",300,-3,3,300,-3,3);

    TH1D* histo_delta_eta = new TH1D("histo_delta_eta","histo_delta_eta",400,-4,4);
    TH1D* histo_delta_eta_z_cut = new TH1D("histo_delta_eta_z_cut","histo_delta_eta_z_cut",400,-4,4);

    TH1D* histo_sum_eta_z_cut = new TH1D("histo_sum_eta_z_cut","histo_sum_eta_z_cut",400,-4,4);

    TH1D* histo_S_mass_r_larger_10 = new TH1D("histo S mass r>10","histo S mass r>10",200,0,20);
    TH1D* histo_S_mass_r_larger_20 = new TH1D("histo S mass r>20","histo S mass r>20",200,0,20);

    TH1D* histo_num_tracks = new TH1D("histo number of tracks","histo number of tracks",200,0,1000);
    TH1D* histo_vertex_z_pos = new TH1D("histo vertex z position","histo vertex z position",100,-30,30);

    TH1D* histo_V0_with_K_plus_radius = new TH1D("histo V0 with K plus radius","histo V0 with K plus radius",100,0,200);

    vector<TH1D*> histos_1D;
    vector<TH2D*> histos_2D;

    vector<TH1D*> histo_delta;

    vector<TH1D*> histos_m_squared_type5;
    vector<TH1D*> histos_m_squared_type51;

    vector<TH1D*> m_squared_type5_K_r;
    vector<TH1D*> m_squared_type5_pi_minus_r;
    vector<TH1D*> m_squared_type5_pi_plus_r;
    vector<TH1D*> m_squared_type5_p_r;

    vector<TH1D*> m_squared_type51_K_r;
    vector<TH1D*> m_squared_type51_pi_minus_r;
    vector<TH1D*> m_squared_type51_pi_plus_r;
    vector<TH1D*> m_squared_type51_p_r;
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
    int counter_draw =0;

    int counter_arr[6]={0,0,0,0,0,0};

    vector<int> counters;

    //-----------------------------------------------------------------


    TH2D* dEdx_vs_charge_dot_momentum;
    TH2D* dEdx_vs_charge_dot_momentum_S_cand;

    TH2D* dEdx_type5_anticuts;

    vector<TH2D*> vec_dEdx_S;

    vector<TH2D*> vec_dEdx;
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
    bool draw=1;


    TString HistName;
    char NoP[50];

    double counter =0;

    int binning_nuclev=300;
    int range_min_nuclev=0;
    int range_max_nuclev=200;

    TH1D* histo_radius_nuclev_2_or_more = new TH1D("nuclev_radius_2_or_more_tracks ","nuclev_radius_2_or_more_tracks",binning_nuclev,range_min_nuclev,range_max_nuclev);
    TH1D* histo_radius_nuclev_4_or_more = new TH1D("nuclev_radius_4_or_more_tracks ","nuclev_radius_4_or_more_tracks",binning_nuclev,range_min_nuclev,range_max_nuclev);
    TH1D* histo_radius_nuclev_4_or_more_momentum = new TH1D("nuclev_radius_4_or_more_tracks_momentum_cut ","nuclev_radius_4_or_more_tracks_momentum_cut",binning_nuclev,range_min_nuclev,range_max_nuclev);
    TH1D* histo_radius_nuclev_5_or_more = new TH1D("nuclev_radius_5_or_more_tracks ","nuclev_radius_5_or_more_tracks",binning_nuclev,range_min_nuclev,range_max_nuclev);
    TH1D* histo_radius_nuclev_3_or_more = new TH1D("nuclev_radius_3_or_more_tracks ","nuclev_radius_3_or_more_tracks",binning_nuclev,range_min_nuclev,range_max_nuclev);

    vector<TH1D*> vec_histo_radius_nuclev;

    TH2D* x_and_y_nuclev_3 = new TH2D("x_and_y_nuclev_3","x_and_y_nuclev_3",200,-100,100,200,-100,100);
    TH2D* x_and_y_nuclev_4 = new TH2D("x_and_y_nuclev_4","x_and_y_nuclev_4",200,-100,100,200,-100,100);
    TH2D* x_and_y_nuclev_5 = new TH2D("x_and_y_nuclev_5","x_and_y_nuclev_5",200,-100,100,200,-100,100);
    TH2D* x_and_y_nuclev_6 = new TH2D("x_and_y_nuclev_6","x_and_y_nuclev_6",200,-100,100,200,-100,100);

    vector<TH2D*> vec_radius_ortho_vs_z;

    vector<TH2D*> vec_x_y_slices;
    vector<TH1D*> vec_radius_slices;
    vector<TH1D*> vec_z_slices_in_radius;
    
    TH1D* histo_path = new TH1D("histo_path ","histo_path",100,-10,10);
    TH1D* histo_nitscls = new TH1D("histo_nitscls ","histo_nitscls",80,0,80);
    TH1D* histo_path_has_its = new TH1D("histo_path_has_its","histo_path_has_its",50,-100,100);
    TH1D* histo_path_has_no_its = new TH1D("histo_path_has_no_its","histo_path_has_no_its",50,-100,100);
    TH1D* histo_path_has_atleast_3_its_hits = new TH1D("histo_path_has_atleast_3_its_hits","histo_path_has_atleast_3_its_hits",80,-100,100);

    //TH1D* histo_ = new TH1D("histo_path_has_atleast_3_its_hits","histo_path_has_atleast_3_its_hits",80,-100,100);
    TH1D* histo_radius_nuclev_4_or_more_sum = new TH1D("nuclev_radius_4_or_more_tracks_sum","nuclev_radius_4_or_more_tracks_sum",binning_nuclev,range_min_nuclev,range_max_nuclev);

    TH1D* histo_its = new TH1D("its number","its number",70,0,70);


    //TNtuple* tpl = new TNtuple ("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");

    TH1D* histo_type_of_S = new TH1D("S_type","S_type",60,0.5,60.5);
    TH1D* histo_angle = new TH1D("histo_angle","histo_angle",100,0,7);
    TH1D* histo_angle_z_0_10 = new TH1D("histo_angle_z_0_10","histo_angle_z_0_10",100,0,7);
    TH1D* histo_angle_z_50_60 = new TH1D("histo_angle_z_50_60","histo_angle_z_50_60",100,0,7);

    TH1D* S_radius_from_origin = new TH1D("S_radius_from_origin","S_radius_from_origin",200,0,200);
    vector<TH1D*> vec_S_radius_from_origin;

    TH2D* radius_ortho_vs_z_eta = new TH2D("R_vs_z_sum_and_diff_eta_larger0.2","R_vs_z_sum_and_diff_eta_larger0.2",500,-200,200,500,0,200);
    TH2D* radius_ortho_vs_z_eta2 = new TH2D("R_vs_z_sum_and_diff_eta_smaller0.2","R_vs_z_sum_and_diff_eta_smaller0.2",500,-200,200,500,0,200);
    TH2D* radius_ortho_vs_z_eta3 = new TH2D("R_vs_z_sum_and_diff_eta_both_smaller0.2","R_vs_z_sum_and_diff_eta_both_smaller0.2",500,-200,200,500,0,200);
    TH2D* radius_ortho_vs_z_eta4 = new TH2D("R_vs_z_sum_and_diff_eta_larger0.2andmom","R_vs_z_sum_and_diff_eta_larger0.2andmom",500,-200,200,500,0,200);
    TH2D* radius_ortho_vs_z_eta5 = new TH2D("R_vs_z_sum_and_diff_eta_larger0.2andangle","R_vs_z_sum_and_diff_eta_larger0.2andangle",500,-200,200,500,0,200);

    TH1D* histo_pt_S = new TH1D("histo_pt_S","histo_pt_S",100,0.,10.);

    TProfile* tprof = new TProfile("tprof","tprof",20,0.5,1.5);

    vector<TH1D*> vec_invmass_Lambda_type5;
    vector<TH1D*> vec_invmass_Lambda_type51;
    vector<TH1D*> vec_invmass_S_type5;
    vector<TH1D*> vec_invmass_S_type51;

    TH1D* histo_invmass_Lambda_type5 = new TH1D("histo_invmass_lambda_type5","histo_invmass_lambda_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_by_V0_type5 = new TH1D("histo_invmass_Lambda_by_V0_type5","histo_invmass_Lambda_by_V0_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_dcaprim_line_type5 = new TH1D("histo_invmass_Lambda_dcaprim_line_type5","histo_invmass_Lambda_dcaprim_line_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_dcaprim1_type5 = new TH1D("histo_invmass_Lambda_dcaprim1_type5","histo_invmass_Lambda_dcaprim1_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_antip_cut_type5 = new TH1D("histo_invmass_Lambda_antip_cut_type5","histo_invmass_Lambda_antip_cut_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_antip_cut_byV0_type51 = new TH1D("histo_invmass_Lambda_antip_cut_byV0_type51","histo_invmass_Lambda_antip_cut_byV0_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_radius_cut_type5 = new TH1D("histo_invmass_Lambda_radius_cut_type5","histo_invmass_Lambda_radius_cut_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_anti_cut_type5 = new TH1D("histo_invmass_Lambda_anti_cut_type5","histo_invmass_Lambda_anti_cut_type5",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_final_cut_type5 = new TH1D("histo_invmass_Lambda_final_cut_type5","histo_invmass_Lambda_final_cut_type5",50*2,1.1,1.13);

    TH1D* histo_invmass_Lambda_type51 = new TH1D("histo_invmass_lambda_type51","histo_invmass_lambda_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_dcaprim_type51 = new TH1D("histo_invmass_Lambda_dcaprim_type51","histo_invmass_Lambda_dcaprim_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_antip_cut_type51 = new TH1D("histo_invmass_Lambda_antip_cut_type51","histo_invmass_Lambda_antip_cut_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_radius_cut_type51 = new TH1D("histo_invmass_Lambda_radius_cut_type51","histo_invmass_Lambda_radius_cut_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_anti_cut_type51 = new TH1D("histo_invmass_Lambda_anti_cut_type51","histo_invmass_Lambda_anti_cut_type51",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_final_cut_type51 = new TH1D("histo_invmass_Lambda_final_cut_type51","histo_invmass_Lambda_final_cut_type51",50*2,1.1,1.13);

    TH1D* histo_S_mass_type5 = new TH1D("histo_S_mass_type5","histo_S_mass_type5",200,0,10);
    TH1D* histo_S_mass_type5_cut_on_antip = new TH1D("histo_S_mass_type5_cut_on_antip","histo_S_mass_type5_cut_on_antip",200,0,10);
    TH1D* histo_S_mass_type5_anticuts = new TH1D("histo_S_mass_type5_anticuts","histo_S_mass_type5_anticuts",200,0,10);

    TH1D* histo_S_mass_type5_anticuts_angle_90 = new TH1D("histo_S_mass_type5_anticuts_angle_90","histo_S_mass_type5_anticuts_angle_90",200,0,10);
    TH1D* histo_S_mass_type5_anticuts_angle_120 = new TH1D("histo_S_mass_type5_anticuts_angle_120","histo_S_mass_type5_anticuts_angle_120",200,0,10);
    TH1D* histo_S_mass_type5_anticuts_angle_150 = new TH1D("histo_S_mass_type5_anticuts_angle_150","histo_S_mass_type5_anticuts_angle_150",200,0,10);

    TH1D* histo_S_mass_type5_dcaprim = new TH1D("histo_S_mass_type5_dcaprim","histo_S_mass_type5_dcaprim",200,0,10);

    TH1D* histo_S_mass_type5_angle_lambda = new TH1D("histo_S_mass_type5_angle_lambda","histo_S_mass_type5_angle_lambda",200,0,10);

    TH1D* histo_S_mass_type5_angle_lambda_and_dcaprim =
        new TH1D("histo_S_mass_type5_angle_lambda_and_dcaprim","histo_S_mass_type5_angle_lambda_and_dcaprim",200,0,10);

    TH1D* histo_S_mass_type51_angle_lambda_and_dcaprim =
        new TH1D("histo_S_mass_type51_angle_lambda_and_dcaprim","histo_S_mass_type51_angle_lambda_and_dcaprim",200,0,10);


    TH1D* histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10
        = new TH1D("histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10","histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10",200,0,10);

    TH1D* histo_eta_type_5 = new TH1D("histo_eta_type_5","histo_eta_type_5",200,-5,5);
    TH1D* histo_phi_type_5 = new TH1D("histo_phi_type_5","histo_phi_type_5",200,-7,2*3.5);

    TH1D* histo_S_mass_mixed_event_no_dcaprim =
        new TH1D("histo_S_mass_mixed_event_no_dcaprim","histo_S_mass_mixed_event_no_dcaprim",200,0,10);
    TH1D* histo_S_mass_mixed_event_with_dcaprim =
        new TH1D("histo_S_mass_mixed_event_with_dcaprim","histo_S_mass_mixed_event_with_dcaprim",200,0,10);


    TH1D* histo_invmass_Lambda_type1 = new TH1D("histo_invmass_Lambda_type1","histo_invmass_Lambda_type1",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_type1_by_track = new TH1D("histo_invmass_Lambda_type1_by_track","histo_invmass_Lambda_type1_by_track",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_type11 = new TH1D("histo_invmass_Lambda_type11","histo_invmass_Lambda_type11",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_type11_by_track = new TH1D("histo_invmass_Lambda_type11_by_track","histo_invmass_Lambda_type11_by_track",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_type1_dcaprim = new TH1D("histo_invmass_Lambda_type1_dcaprim","histo_invmass_Lambda_type1_dcaprim",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_type1_dcaprim_1 = new TH1D("histo_invmass_Lambda_type1_dcaprim_1","histo_invmass_Lambda_type1_dcaprim_1",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_overlap_cut_type1 = new TH1D("histo_invmass_Lambda_overlap_cut_type1","histo_invmass_Lambda_overlap_cut_type1",50*2,1.1,1.13);
    TH1D* histo_invmass_Lambda_overlap_cut_dcaprim_type1 = new TH1D("histo_invmass_Lambda_overlap_cut_dcaprim_type1","histo_invmass_Lambda_overlap_cut_dcaprim_type1",50*2,1.1,1.13);
    TH1D* m_squared_proton_type1= new TH1D("m_squared_proton_type1", "m_squared_proton_type1",100,-0.1,1.5);
    TH1D* m_squared_proton_type11= new TH1D("m_squared_proton_type11", "m_squared_proton_type11",100,-0.1,1.5);

    TH1D* histo_diff_inv_mass_type1 = new TH1D("histo_diff_inv_mass_type1","histo_diff_inv_mass_type1",50*2,-0.04,0.04);

    TH1D* histo_invariantmass_K0_type1 = new TH1D("histo_invariantmass_K0_type1","histo_invariantmass_K0_type1",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type1_dcaprim = new TH1D("histo_invariantmass_K0_type1_dcaprim","histo_invariantmass_K0_type1_dcaprim",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type1_overlap = new TH1D("histo_invariantmass_K0_type1_overlap","histo_invariantmass_K0_type1_overlap",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type1_cut_on_L = new TH1D("histo_invariantmass_K0_type1_cut_on_L","histo_invariantmass_K0_type1_cut_on_L",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type11 = new TH1D("histo_invariantmass_K0_type11","histo_invariantmass_K0_type11",50*3,0.4,0.6);

    TH1D* histo_invariantmass_K0_type3 = new TH1D("histo_invariantmass_K0_type3","histo_invariantmass_K0_type3",50*3,0.4,0.6);
    TH1D* histo_invariantmass_K0_type31 = new TH1D("histo_invariantmass_K0_type31","histo_invariantmass_K0_type31",50*3,0.4,0.6);

    vector<TH1D*> vec_invmass_K0_type3;
    vector<TH1D*> vec_invmass_K0_type31;


    vector<vector<TH1D*>> vec_vec_invmass_K0_type3;
    vector<vector<TH1D*>> vec_vec_invmass_K0_type31;

    vector<vector<TH1D*>> vec_vec_invmass_L_type5;
    vector<vector<TH1D*>> vec_vec_invmass_L_type51;

    vector<TH1D*> vec_m_squared_type_3_K;
    vector<TH1D*> vec_m_squared_type_3_p;
    vector<TH1D*> vec_m_squared_type_3_add_pi;
    vector<TH1D*> vec_m_squared_type_3_K01;
    vector<TH1D*> vec_m_squared_type_3_K02;

    vector<TH1D*> vec_m_squared_type_31_K;
    vector<TH1D*> vec_m_squared_type_31_p;
    vector<TH1D*> vec_m_squared_type_31_add_pi;
    vector<TH1D*> vec_m_squared_type_31_K01;
    vector<TH1D*> vec_m_squared_type_31_K02;

    vector<TH1D*> vec_m_sq_p_type_3_for_dcas;
    vector<TH1D*> vec_m_sq_p_type_31_for_dcas;

    vector<TH1D*> vec_m_squared_type_2_p;
    vector<TH1D*> vec_m_squared_type_2_K_plus1;
    vector<TH1D*> vec_m_squared_type_2_K_plus2;

    vector<TH1D*> vec_m_squared_type_21_p;
    vector<TH1D*> vec_m_squared_type_21_K_plus1;
    vector<TH1D*> vec_m_squared_type_21_K_plus2;

    vector<TH1D*> vec_S_mass_ch3;
    vector<TH1D*> vec_S_mass_ch31;

    TH1D* histo_angle_ch3 = new TH1D("histo_angle_ch3","histo_angle_ch3",360,-180,180);

    TH1D* TOF_eff_K_ch3 = new TH1D("TOF_eff_K_ch3","TOF_eff_K_ch3",2,0.5,2.5);

    TH1D* m_squared_p_ch3 = new TH1D("m_squared_p_ch3","m_squared_p_ch3",400,-0.4,1.4);

    vector<TH1D*> vec_S_mass_ch3_normalband;
    vector<TH1D*> vec_S_mass_ch3_sideband;

    vector<TH1D*> vec_S_mass_ch5_normalband;
    vector<TH1D*> vec_S_mass_ch5_sideband;

    TH1D* hist_radiusK0 = new TH1D("hist_radiusK0","hist_radiusK0",1000,0.,100.);

    TH1D* histo_n_tracks = new TH1D("histo_n_tracks","histo_n_tracks",15000,0.,15000.);



public:
    Ali_Dark_Matter_Read(){}
    Ali_Dark_Matter_Read(TString list);
    //Ali_Dark_Matter_Read(TString list);
    ~Ali_Dark_Matter_Read();
    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t event );
    Int_t Loop_event_S_search(Long64_t event );
    Long64_t getnumberentries(){return file_entries_total;};
    void copy_track_params(Ali_AS_Track* track_in, Ali_AS_Track* track_out);
    void copy_dm_params(Ali_AS_DM_particle* dm_in, Ali_AS_DM_particle* dm_out);
    void Save();
    int DM_Analysis_type5(Ali_AS_DM_particle* DM,int mode);
    void Analyse_Mixed_Events();
    vector<Ali_AS_DM_particle> mix_event(vector<vector<vector<Ali_AS_DM_particle>>> &vec,int eta_cat,int phi_cat);

    float arr_distance_prim_sec[15]={};
    float arr_distance_daughter_particles[1]={0.1};
    float arr_dca_AB[15]={};
    float arr_dca_daughter_prim[15]={};
    float arr_radius_variation[45]={};
    vector<vector<vector<Ali_AS_DM_particle>>> mixed_event_vec;
    vector<Ali_AS_DM_particle> all_mixed_events;

    ClassDef(Ali_Dark_Matter_Read, 1)


};
//----------------------------------------------------------------------------------------



#endif // __ALI_DARK_MATTER_READ_H__




