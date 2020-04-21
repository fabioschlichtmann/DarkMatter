#ifndef __TBASE_TRD_CALIB_H__
#define __TBASE_TRD_CALIB_H__

using namespace std;
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TChain.h"

#include "TObject.h"

#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"
//#include "Ana_Digits_functions.h"

ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_V0)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_Event)



//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Dark_Matter_Read
{
private:
    Ali_AS_Event*     AS_Event;
    Ali_AS_Track*     AS_Track;
    Ali_AS_V0*        AS_V0;
    Ali_AS_Tracklet*  AS_Tracklet;
    Ali_AS_offline_Tracklet*  AS_offline_Tracklet;
    Ali_AS_TRD_digit* AS_Digit;


    TFile* outputfile;
    TFile* outputfile_trkl;

    TFile* calibration_params;


    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;


    Long64_t Event_active = 0;
    

    TChain* input_SE;
    TString JPsi_TREE   = "Tree_AS_Event";
    TString JPsi_BRANCH = "Tree_AS_Event_branch";
    Long64_t file_entries_total;

    Float_t digit_pos[3];

    AliHelix aliHelix;
   

    TString HistName;
    char NoP[50];


public:
    Dark_Matter_Read();
    ~Dark_Matter_Read();
    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t event);
   

    ClassDef(Dark_Matter_Read, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Dark_Matter_Read::Dark_Matter_Read()
{
    //outputfile = new TFile("./TRD_Calib.root","RECREATE");
    outputfile = new TFile("./Results_Dark_Matter.root","RECREATE");
   // outputfile_trkl = new TFile("./TRD_Calib_on_trkl.root","RECREATE");



}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Dark_Matter_Read::Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;
    // TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_offline_calib/Data/";
    TString pinputdir = "/home/ceres/schlichtmann/ESD_Analysis/";
    //TString pinputdir = "/home/ceres/berdnikova/TRD-Run3-Calibration/";

    AS_Event = new Ali_AS_Event();
    AS_V0    = new Ali_AS_V0();
    AS_Track = new Ali_AS_Track();
    AS_Tracklet = new Ali_AS_Tracklet();
    AS_Digit = new Ali_AS_TRD_digit();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( JPsi_TREE.Data(), JPsi_TREE.Data() );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_SE ->AddFile(addfile.Data(),-1, JPsi_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( JPsi_BRANCH, &AS_Event );
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
        }
    }

    file_entries_total = input_SE->GetEntries();
    N_Events = file_entries_total;
    cout << "Total number of events in tree: " << file_entries_total << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Dark_Matter_Read::Loop_event(Long64_t event)
{
    printf("Loop event number: %lld \n",event);

    Event_active = event;

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;

    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Int_t    N_TRD_tracklets_online = AS_Event ->getNumTracklets(); // online tracklet
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    UShort_t NumV0s = AS_Event ->getN_V0s();

    //loop over V0s
    for(int V0counter=0; v0counter<NumV0s ; NumV0s++)
    {
        AS_V0 = AS_Event -> getV0(V0counter);

        //get position of V0
        Float_t x;
        Float_t y;
        Float_t z;
        AS_V0 -> getx(x);
        AS_V0 -> gety(y);
        AS_V0 -> getz(z);

        //get two tracks for each V0
        as_trackP = AS_V0 -> getTrack(0);
        as_trackN = AS_V0 -> getTrack(1);



        //getters for track P
        /*
         Float_t nsig_e_TPC[2] = {as_trackP->getnsigma_e_TPC(),as_trackN->getnsigma_e_TPC()};
	Float_t getnsigma_e_TOF();
	Float_t getnsigma_pi_TPC();
	Float_t getnsigma_pi_TOF();
	Float_t getnsigma_K_TPC();
	Float_t getnsigma_K_TOF();
	Float_t getnsigma_p_TPC();
	Float_t getnsigma_p_TOF();
	Float_t getTRDSignal();
	Float_t getTRDsumADC();
	Float_t getdca();
	TLorentzVector get_TLV_part();
	UShort_t getNTPCcls();
	UShort_t getNTRDcls();
	UShort_t getNITScls();
	UShort_t getStatus();
	Float_t  getTPCchi2();
        ULong64_t getTRD_layer(Int_t i_layer);
        Float_t   getimpact_angle_on_TRD();
	Float_t   getTPCdEdx();
	Float_t   getTOFsignal();
        Float_t   getTrack_length();
        Float_t   getHelix_param(Int_t i_param);
        */

        // Lambda0 -> proton + pi-
        // anti-Lambda0 -> anti-proton + pi+
        // if(fabs(nsig_pi_TPC[0]) < 2.5)

        // TLorentzVector* tlv_pos = new TLorentzVector();
        // TLorentzVector* tlv_neg;

        // Double_t Energy = mass
        // tlv_pos ->SetPxPyPzE(Px,Py,Pz,Energy);

        // TLorentzVector* tlv_Lambda = tlv_pos + tlv_neg;
        // Double_t M_inv = tlv_Lambda.M();
        // fill 1D hist with M_inv
    }




    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------

    //printf("N_TRD_tracklets_online: %d \n",N_TRD_tracklets_online);

    N_Tracks = NumTracks;

   /* vec_track_info.clear();
    vec_digit_track_info.clear();
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_Tracklet_pos.clear();
        vec_TV3_Tracklet_dir.clear();
    } */

   /* // Loop over all tracklets
    Double_t scale_factor_length = 2.0;
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_online; ++i_tracklet) // loop over all tracklets of the actual event
    {
        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset(); // online tracklets
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();    // online tracklets
        Short_t  i_det_tracklet = AS_Tracklet ->get_detector();

        vec_TV3_Tracklet_pos[i_det_tracklet].push_back(TV3_offset);
        vec_TV3_Tracklet_dir[i_det_tracklet].push_back(TV3_dir);

        Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
        if(impact_angle > TMath::Pi()*0.5) impact_angle -= TMath::Pi();

        //if(!(i_det_tracklet >= 114 && i_det_tracklet <= 119)) continue;
        //if(!(i_det_tracklet >= 378 && i_det_tracklet <= 382)) continue;

        Float_t points[6] =
        {
            (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
            (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
        };
        //printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_online,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);

        vec_TPL3D_tracklets.push_back(new TEveLine());
        vec_TPL3D_tracklets[(Int_t)vec_TPL3D_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]));
        vec_TPL3D_tracklets[(Int_t)vec_TPL3D_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2]));
    }   */

    // Loop over all tracks
   /* for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
    {
        //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = AS_Track ->getTRDSignal();
        Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
        Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = AS_Track ->get_TLV_part();
        UShort_t NTPCcls        = AS_Track ->getNTPCcls();
        UShort_t NTRDcls        = AS_Track ->getNTRDcls();
        UShort_t NITScls        = AS_Track ->getNITScls();
        Float_t TPCchi2         = AS_Track ->getTPCchi2();
        Float_t TPCdEdx         = AS_Track ->getTPCdEdx();
        Float_t TOFsignal       = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = AS_Track ->getTrack_length();

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        vec_track_single_info[0]  = dca;
        vec_track_single_info[1]  = TPCdEdx;
        vec_track_single_info[2]  = momentum;
        vec_track_single_info[3]  = eta_track;
        vec_track_single_info[4]  = pT_track;
        vec_track_single_info[5]  = TOFsignal;
        vec_track_single_info[6]  = Track_length;
        vec_track_single_info[7]  = TRDsumADC;
        vec_track_single_info[8]  = TRD_signal;
        vec_track_single_info[9]  = nsigma_TPC_e;
        vec_track_single_info[10] = nsigma_TPC_pi;
        vec_track_single_info[11] = nsigma_TPC_p;
        vec_track_info.push_back(vec_track_single_info);

        //----------------------------------------------
        // TRD digit information
        UShort_t  fNumTRDdigits        = AS_Track ->getNumTRD_digits();
        UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
        //--------------------------


        //--------------------------
        // Offline tracklet loop
        for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
        {
            AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
            TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
            TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
            Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();

            printf("offline, i_tracklet: %d, offset: {%4.3f, %4.3f, %4.3f}, dir: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,(Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),(Float_t)TV3_dir[0],(Float_t)TV3_dir[1],(Float_t)TV3_dir[2]);

            vec_TPL3D_offline_tracklets.push_back(new TEveLine());
            vec_TPL3D_offline_tracklets[(Int_t)vec_TPL3D_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]));
            vec_TPL3D_offline_tracklets[(Int_t)vec_TPL3D_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2]));
        }
        //--------------------------


        printf("i_track: %d, pT: %4.3f, phi: %4.3f, eta: %4.3f, N_off_trkl: %d \n",i_track,pT_track,phi_track,eta_track,fNumOfflineTracklets);

        //printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

        vec_digit_info.clear();
        for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
        {
            //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
            AS_Digit              = AS_Track ->getTRD_digit(i_digits);
            Int_t    layer        = AS_Digit ->get_layer();
            Int_t    sector       = AS_Digit ->get_sector();
            Int_t    column       = AS_Digit ->get_column();
            Int_t    stack        = AS_Digit ->get_stack();
            Int_t    row          = AS_Digit ->get_row();
            Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);
            Float_t  dca_to_track = AS_Digit ->getdca_to_track();
            Float_t  dca_x        = AS_Digit ->getdca_x();
            Float_t  dca_y        = AS_Digit ->getdca_y();
            Float_t  dca_z        = AS_Digit ->getdca_z();
            Float_t  ImpactAngle  = AS_Digit ->getImpactAngle();

            for(Int_t i_time = 0; i_time < 24; i_time++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz);
                }
                //printf("track: %d/%d, digit: %d/%d \n",i_track,NumTracks,i_digits,fNumTRDdigits);
                //printf("pos: {%4.3f, %4.3f, %4.3f} \n,",digit_pos[0],digit_pos[1],digit_pos[2]);
                vec_TPM3D_digits[layer] ->SetNextPoint(digit_pos[0],digit_pos[1],digit_pos[2]);

                Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time);


                // x,y,z,time,ADC,sector,stack,layer,row,column,dca
                vec_digit_single_info[0]  = digit_pos[0];
                vec_digit_single_info[1]  = digit_pos[1];
                vec_digit_single_info[2]  = digit_pos[2];
                vec_digit_single_info[3]  = i_time;
                vec_digit_single_info[4]  = ADC;
                vec_digit_single_info[5]  = sector;
                vec_digit_single_info[6]  = stack;
                vec_digit_single_info[7]  = layer;
                vec_digit_single_info[8]  = row;
                vec_digit_single_info[9]  = column;
                vec_digit_single_info[10] = dca_to_track;
                vec_digit_single_info[11] = dca_x;
                vec_digit_single_info[12] = dca_y;
                vec_digit_single_info[13] = dca_z;

                //printf("dca_full: %4.3f, dca: {%4.3f, %4.3f, %4.3f} \n",dca_to_track,dca_x,dca_y,dca_z);
                vec_digit_info.push_back(vec_digit_single_info);

                N_Digits ++;
            }

        } // end of digits loop

        vec_digit_track_info.push_back(vec_digit_info);

    } // end of track loop */

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------




#endif // __TBASE_TRD_CALIB_H__