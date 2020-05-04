

#ifndef __TBASE_TRD_CALIB_H__
#define __TBASE_TRD_CALIB_H__

using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TChain.h"

#include "TObject.h"

#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"
#include "vertex_modified.h"
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
    Int_t Loop_event(Long64_t event, vector<TH1D*> histos_1D, vector<TH2D*> histos_2D, double& number_event_counter);
    Long64_t getnumberentries(){return file_entries_total;};
   

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
    //TString pinputdir = "/home/ceres/schlichtmann/ESD_Analysis/";
    TString pinputdir = "/misc/alidata120/alice_u/schlichtmann/dark_matter/";
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
Int_t Dark_Matter_Read::Loop_event(Long64_t event, vector<TH1D*> histos_1D,vector<TH2D*> histos_2D,double& number_event_counter)
{
    printf("Loop event number: %lld \n",event);
    number_event_counter++;
    cout<<""<<endl;

    

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

    histos_1D[7]->Fill(NumTracks);
    histos_1D[8]->Fill(EventVertexZ);

    //printf("Event vertex: %f", EventVertexX);
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    UShort_t NumV0s = AS_Event ->getN_V0s();
    Float_t* pos = new Float_t[3];
    Float_t* momP = new Float_t[3];
    Float_t* momN = new Float_t[3];
    TLorentzVector* tlv_pos = new TLorentzVector();
    TLorentzVector* tlv_neg = new TLorentzVector();
    TLorentzVector* tlv_Lambda = new TLorentzVector();
    TLorentzVector* tlv_Kaon = new TLorentzVector();
    Ali_AS_Track* as_trackP = new Ali_AS_Track;
    Ali_AS_Track* as_trackN = new Ali_AS_Track;
    Float_t* sigma_proton_TPC = new Float_t[2];
    Float_t* sigma_pion_TPC = new Float_t[2];
    Float_t radius;
    TVector3 position_SV2;
    TVector3 position_SV3;
    TVector3 direction_SV2;
    TVector3 direction_SV3;


    vector<TVector3> vec_position_SV2;
    vector<TVector3> vec_position_SV3;
    vector<TVector3> vec_direction_SV2;
    vector<TVector3> vec_direction_SV3;

    vector<vector<vector<TVector3>>> cat_direction_SV2;
    vector<vector<vector<TVector3>>> cat_position_SV2;

    vector<vector<vector<TVector3>>> cat_direction_SV3;
    vector<vector<vector<TVector3>>> cat_position_SV3;

    cat_direction_SV2.resize(4);
    cat_direction_SV2[0].resize(4);
    cat_direction_SV2[1].resize(4);
    cat_direction_SV2[2].resize(4);
    cat_direction_SV2[3].resize(4);

    cat_direction_SV3.resize(4);
    cat_direction_SV3[0].resize(4);
    cat_direction_SV3[1].resize(4);
    cat_direction_SV3[2].resize(4);
    cat_direction_SV3[3].resize(4);

    cat_position_SV2.resize(4);
    cat_position_SV2[0].resize(4);
    cat_position_SV2[1].resize(4);
    cat_position_SV2[2].resize(4);
    cat_position_SV2[3].resize(4);

    cat_position_SV3.resize(4);
    cat_position_SV3[0].resize(4);
    cat_position_SV3[1].resize(4);
    cat_position_SV3[2].resize(4);
    cat_position_SV3[3].resize(4);
   

    vector<double> limits_num_tracks;
    limits_num_tracks.push_back(0);
    limits_num_tracks.push_back(47);
    limits_num_tracks.push_back(67);
    limits_num_tracks.push_back(87);
    limits_num_tracks.push_back(1000);

    vector<double> limits_z_vertex;
    limits_z_vertex.push_back(-30);
    limits_z_vertex.push_back(-3.3);
    limits_z_vertex.push_back(1.5);
    limits_z_vertex.push_back(5.7);
    limits_z_vertex.push_back(30);

    const Float_t mass_proton = 0.938 ;  //in GeV?
    const Float_t mass_pion = 0.1396 ;  //in GeV?

    //loop over V0s
    for(int V0counter=0; V0counter<NumV0s ; V0counter++)
    {
        AS_V0 = AS_Event -> getV0(V0counter);

        //get position of V0
        pos = AS_V0 -> getxyz();
        //cout<<"posx: "<<pos[0]<<endl;
        //position.SetXYZ(pos[0],pos[1],pos[2]);
        radius = sqrt( pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2] );

        //printf("x %f,y %f, z %f \n",pos[0],pos[1],pos[2]);

        //get momentum of posititve particle
        momP = AS_V0 -> getPpxpypz();

        //printf("momentum of positive particle px: %f,py: %f,pz: %f \n",momP[0],momP[1],momP[2]);

        //get momentum for negative particle
      
        momN = AS_V0 -> getNpxpypz();

        //printf("momentum of negative particle px: %f,py: %f,pz: %f \n",momN[0],momN[1],momN[2]);

        //--------------------------------------------------------------------------------
        //get two tracks for each V0

        as_trackP = AS_V0 -> getTrack(0);
        as_trackN = AS_V0 -> getTrack(1);

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        Double_t invariantmass = -1.;
        //_____________________________________________________________________----

        sigma_proton_TPC[0] = as_trackP -> getnsigma_p_TPC();
        sigma_proton_TPC[1] = as_trackN -> getnsigma_p_TPC();
        //printf("sigmas proton %f  %f", sigma_proton_TPC[0], sigma_proton_TPC[1]) ;



        sigma_pion_TPC[0]  = as_trackP -> getnsigma_pi_TPC();
        sigma_pion_TPC[1]  = as_trackN -> getnsigma_pi_TPC();


        //Lambda0 -> proton + pi-
        //check if positive particle is proton and if negative particle is pion-
        if(fabs(sigma_proton_TPC[0])<2.5 && fabs(sigma_pion_TPC[1])<2.5)
        {
           // printf("particles are proton and pion- \n");

            Float_t energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            Float_t energy_pion = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();
            //cout<<invariantmass<<endl;
            histos_1D[0]->Fill(invariantmass);

            histos_1D[2]->Fill(radius);

            histos_2D[0]->Fill(pos[0],pos[1]);


            //cut on mass
            if(invariantmass<1.1157+0.001495*2 && invariantmass > 1.1157-0.001495*2)
            {

                position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                for(int cat_track=0;cat_track<4;cat_track++)
                {
                    //cout<<cat_track<<endl;
                    for(int cat_z=0;cat_z<4;cat_z++)
                    {
                        //cout<<cat_z<<endl;
                        if(NumTracks>=limits_num_tracks[cat_track] && NumTracks<limits_num_tracks[cat_track+1])
                        {
                            if(EventVertexZ>=limits_z_vertex[cat_z] && EventVertexZ<limits_z_vertex[cat_z+1])
                            {
                               cat_direction_SV3[cat_track][cat_z].push_back(direction_SV3);
                               cat_position_SV3[cat_track][cat_z].push_back(position_SV3);

                            }

                        }

                    }

                }


            }
        }

        //check if negative particle is antiproton and if positive particle is pion+
        if(fabs(sigma_proton_TPC[1])<2.5 && fabs(sigma_pion_TPC[0])<2.5)
        {
           // printf("particles are antiproton and pion+ \n");
            Float_t energy_antiproton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            Float_t energy_pion = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();
           // cout<<invariantmass<<endl;
            histos_1D[0]->Fill(invariantmass);

            histos_1D[2]->Fill(radius);

            histos_2D[0]->Fill(pos[0],pos[1]);

            //cut on mass
            if(invariantmass<1.1157+0.001495*2 && invariantmass > 1.1157-0.001495*2)
            {

                position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                 for(int cat_track=0;cat_track<4;cat_track++)
                {
                    for(int cat_z=0;cat_z<4;cat_z++)
                    {
                        if(NumTracks>=limits_num_tracks[cat_track] && NumTracks<limits_num_tracks[cat_track+1])
                        {
                            if(EventVertexZ>=limits_z_vertex[cat_z] && EventVertexZ<limits_z_vertex[cat_z+1])
                            {
                               cat_direction_SV3[cat_track][cat_z].push_back(direction_SV3);
                               cat_position_SV3[cat_track][cat_z].push_back(position_SV3);

                            }

                        }

                    }
                }
            }
        }

        //K0 - > pi+ and pi-
        //check if pion+ and pion-
        if(fabs(sigma_pion_TPC[0])<2.5 && fabs(sigma_pion_TPC[1])<2.5)
        {
           //printf("particles are pion+ and pion- \n");
            Float_t energy_pion_plus = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            Float_t energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();
            //cout<<invariantmass<<endl;
            histos_1D[1]->Fill(invariantmass);

            histos_1D[3]->Fill(radius);

            histos_2D[1]->Fill(pos[0],pos[1]);

            //cut on mass
            if(invariantmass<0.4981+0.0042*2 && invariantmass > 0.4981-0.0042*2)
            {
                printf("particles are pion+ and pion- \n");

                position_SV2.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV2.push_back(position_SV2);

                direction_SV2.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                vec_direction_SV2.push_back(direction_SV2);

                printf("position of vertex:  %f %f %f \n",pos[0],pos[1],pos[2]);
                //cout<<""<<endl;
                printf("momentum of pion + : %f %f %f \n",momP[0],momP[1],momP[2]);
                printf("momentum of pion - : %f %f %f \n",momN[0],momN[1],momN[2]);
                cout<<""<<endl;
                //cout<<"direction SV2: "<<direction_SV2[0]<<endl;
                //cout<<"positionSV2: "<<position_SV2[0]<<endl;
                //cout<<"position SV2: "<<position_SV2[0]<<endl;

                 for(int cat_track=0;cat_track<4;cat_track++)
                {
                    for(int cat_z=0;cat_z<4;cat_z++)
                    {
                        if(NumTracks>=limits_num_tracks[cat_track] && NumTracks<limits_num_tracks[cat_track+1])
                        {
                            if(EventVertexZ>=limits_z_vertex[cat_z] && EventVertexZ<limits_z_vertex[cat_z+1])
                            {
                               cat_direction_SV2[cat_track][cat_z].push_back(direction_SV2);
                               cat_position_SV2[cat_track][cat_z].push_back(position_SV2);

                            }

                        }

                    }
                }
            }
        }

        
    }     //end of V0 loop



    float radiusS;
    TLorentzVector* tlv_SV2 = new TLorentzVector();  //kaon
    TLorentzVector* tlv_SV3 = new TLorentzVector();  //lambda
    TLorentzVector* tlv_neutron = new TLorentzVector();  //neutron

    TLorentzVector* tlv_SV1 = new TLorentzVector();  //tlv_SV2+tlv_SV3-tlv_neutron

    double mass_K0 = 0.493677;
    double mass_Lambda = 1.1155683;
    double mass_neutron = 0.939565;
    double S_mass = -1;

    tlv_neutron->SetPxPyPzE(0.,0.,0.,mass_neutron);

    TVector3 vec_primary_vertex_to_SV1;
    TVector3 unit_vec_primary_vertex_to_SV1;

    TVector3 momentum_SV1;
    TVector3 unit_momentum_SV1;


    for (int vector_loop_1=0; vector_loop_1<vec_position_SV3.size();vector_loop_1++)
    {
       for( int vector_loop_2=0; vector_loop_2<vec_position_SV2.size();vector_loop_2++)
           if(vec_position_SV3.size()>0 && vec_position_SV2.size()>0 && vec_direction_SV3.size()>0 &&  vec_direction_SV2.size()>0)
           {
               if(calculateMinimumDistance(vec_position_SV3[vector_loop_1],vec_direction_SV3[vector_loop_1],vec_position_SV2[vector_loop_2],vec_direction_SV2[vector_loop_2])>1.){continue;}

               TVector3 S_vertex_pos = calcVertexAnalytical(vec_position_SV3[vector_loop_1],vec_direction_SV3[vector_loop_1],vec_position_SV2[vector_loop_2],vec_direction_SV2[vector_loop_2]);

               //build lorentz vector of SV2 and SV3
               //SV2
               double px,py,pz,E;
               px = vec_direction_SV2[vector_loop_2][0];
               py = vec_direction_SV2[vector_loop_2][1];
               pz = vec_direction_SV2[vector_loop_2][2];
               E = (sqrt(px*px+py*py+pz*pz)+(mass_K0*mass_K0));
               tlv_SV2 -> SetPxPyPzE(px,py,pz,E);

               //SV3
               px = vec_direction_SV3[vector_loop_1][0];
               py = vec_direction_SV3[vector_loop_1][1];
               pz = vec_direction_SV3[vector_loop_1][2];
               E = (sqrt(px*px+py*py+pz*pz)+(mass_Lambda*mass_Lambda));
               tlv_SV3 -> SetPxPyPzE(px,py,pz,E);

               //calculate SV1
               *tlv_SV1 = *tlv_SV2 + *tlv_SV3 - *tlv_neutron;
               S_mass = tlv_SV1->M();


               //check if SV1 vector is parallel to vertex from primary vertex to SV1 (s-vertex)

               // vector from primary vertex to SV1 (S-vertex)
               vec_primary_vertex_to_SV1.SetXYZ(S_vertex_pos[0]-EventVertexX ,S_vertex_pos[1]-EventVertexY , S_vertex_pos[2]-EventVertexZ);
               unit_vec_primary_vertex_to_SV1 = vec_primary_vertex_to_SV1.Unit();

               //build unit vector of momentum from tlv_SV1
               momentum_SV1.SetXYZ(tlv_SV1->Px(),tlv_SV1->Py(),tlv_SV1->Pz());
               unit_momentum_SV1 = momentum_SV1.Unit();

               double dot_product = unit_vec_primary_vertex_to_SV1.Dot(unit_momentum_SV1);
              // printf("dot product: %f \n", dot_product);

               //check if dot product is larger than 0.8
               if(dot_product<0.8){continue;}
               
               radiusS = sqrt ( S_vertex_pos[0]*S_vertex_pos[0]+S_vertex_pos[1]*S_vertex_pos[1]+S_vertex_pos[2]*S_vertex_pos[2] ) ;
               if( fabs(radiusS) < 200 )
               {
                   //cout<<S_vertex_pos[0]<<endl;
                   histos_1D[4]->Fill(radiusS);

                   histos_2D[2]->Fill(S_vertex_pos[0],S_vertex_pos[1]);

                   if(fabs(radiusS)>10)
                   {
                       histos_1D[5]->Fill(S_mass);
                   }

                   if(fabs(radiusS)>20)
                   {
                       histos_1D[6]->Fill(S_mass);
                   }
               }


               //printf("S mass: %f", S_mass);

              


           }


    }
    for (int i = 0;i<3;i++)
    {
        //cout<<i<<endl;
        if(cat_direction_SV2[0][0].size()==0 ){continue;}
        //printf("2 direction SV2: %f \n",cat_direction_SV2[0][0][i][0]);
    }

    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------


    N_Tracks = NumTracks;


    return 1;

    
}





//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------




#endif // __TBASE_TRD_CALIB_H__