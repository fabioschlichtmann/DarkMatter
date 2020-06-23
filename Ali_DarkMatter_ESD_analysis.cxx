#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

//------------------------
#include "AliHelix.h"
#include "TLorentzVector.h"
#include "TSystem.h"
//------------------------

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliKalmanTrack.h"

//#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
//#include "AliTRDseedV1.h"
#include "AliESDfriend.h"

//#include "AliTRDdigitsManager.h"
//#include "AliTRDarrayADC.h"

#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliESDtrackCuts.h"

#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDRun.h"

#include "AliMultSelection.h"

#include "AliCDBEntry.h"
#include "TClonesArray.h"
#include "TGeoMatrix.h"
#include "AliAlignObjParams.h"

//#include "AliTRDdigitsParam.h"
#include "AliRunTag.h"
#include "TObjString.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
//#include "AliTRDCalPad.h"
//#include "AliTRDCalDet.h"
//#include "AliTRDCalOnlineGainTable.h"
//#include "AliTRDCalROC.h"
#include "TPolyMarker.h"

#include "AliTRDCommonParam.h"

#include "AliESDTrdTracklet.h"
#include "TProfile.h"
#include "TGraph.h"

//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
//------------------------

#include "Ali_DarkMatter_ESD_analysis.h"


#include <iostream>
#include <iomanip>
using namespace std;

static Int_t flag_plot_event = 0;
static TString HistName;

static TFile* dfile;
static const char *pathdatabase="alien://folder=/alice/data/2016/OCDB"; // for pPb
//static const char *pathdatabase="alien://folder=/alice/data/2015/OCDB"; // for PbPb

static AliTRDCommonParam* fParam;
static AliESDfriend *esdFr = NULL;
static AliESDInputHandler *esdH = NULL;
static TString esdFriendTreeFName;

ClassImp(Ali_AS_Event)
ClassImp(Ali_AS_V0)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_DarkMatter_ESD_analysis)

    //________________________________________________________________________
    Ali_DarkMatter_ESD_analysis::Ali_DarkMatter_ESD_analysis(const char *name)
    : AliAnalysisTaskSE(name),
    AS_Event(0),AS_V0(0),AS_Track(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0),
    fListOfHistos(0x0),fTree(0x0),h_dca(0x0),h_dca_xyz(0x0), h2D_TPC_dEdx_vs_momentum(0x0), fPIDResponse(0), EsdTrackCuts(0)
{
    // Constructor

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());

    // Output slot #0 id reserved by the base class for AOD
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());

}
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//new class definition


void print_int_vector(vector<Int_t> vec)
{
    for(Int_t i=0;i<vec.size();i++)
    {
        cout<<"Vektor  "<<i<<": "<<vec[i]<<endl;

    }
}

//_______________________________________________________________________

Bool_t Ali_DarkMatter_ESD_analysis::UserNotify()
{
    cout << "" << endl;
    cout << "In UserNotify" << endl;

    fParam = AliTRDCommonParam::Instance();


    cout << "Digits file pointers deleted" << endl;

    cout << "All pointers deleted" << endl;

    //AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
    //    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if ( ! esdH ) return kFALSE;
    if ( ! esdH->GetTree() ) return kFALSE;
    if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;


    //-----------------------------------
    TList* list = esdH->GetUserInfo();
    

    TString fname = esdH->GetTree()->GetCurrentFile()->GetName();
    TString Tree_name = esdH->GetTree()->GetName();
    FileStat_t file_stat;
    Int_t PathInfo = gSystem->GetPathInfo(fname.Data(),file_stat);
    cout << "PathInfo: " << PathInfo << ", fname: " << fname << endl;
    //TFile* file = TFile::Open(fname.Data());
    //cout << "Zombie: " << file->IsZombie() << ", header size: " << file->Sizeof() << ", FileBytesRead: " << file->GetFileBytesRead() << endl;

    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!inputHandler)
    {
	printf("WARNING: Inputhandler not available \n");
    }
    else
    {
	printf("Inputhandler available \n");

	fPIDResponse = inputHandler->GetPIDResponse();

        cout << "Got PID response" << endl;
    }



    fEventNoInFile = -1;
    N_good_events  = 0;

#if 0
    // Connect the friends
    //fESD ->SetESDfriend(esdFr);

    esdFriendTreeFName = fname;
    TString basename = gSystem->BaseName(esdFriendTreeFName);
    Int_t index = basename.Index("#")+1;
    basename.Remove(index);
    basename += "AliESDfriends.root";
    TString dirname = gSystem->DirName(esdFriendTreeFName);
    dirname += "/";
    esdFriendTreeFName = dirname + basename;
    cout << "Friend name: " << esdFriendTreeFName.Data() << endl;
#endif

    EsdTrackCuts = new AliESDtrackCuts();

    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    EsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.52);
    EsdTrackCuts->AliESDtrackCuts::SetMinNClustersTPC(50); // 60, Automatically requires TPC refitted tracks?
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXY(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetPtRange(0.15,200.0); // 0.15, 200.0
    EsdTrackCuts->AliESDtrackCuts::SetEtaRange(-1.0,1.0); // 0.85

    // create the digits manager
    cout << "" << endl;
    cout << "________________________________________________________________________" << endl;
    cout << "Created AliTRDdigitsManager" << endl;
    cout << "" << endl;

   
    return kTRUE;
}



//________________________________________________________________________
void Ali_DarkMatter_ESD_analysis::UserCreateOutputObjects()
{
    cout << "" << endl;
    cout << "In UserCreateOutputObjects" << endl;


    OpenFile(1);
    cout << "File opened" << endl;

    fListOfHistos = new TList();
    fListOfHistos ->SetOwner();

   
    OpenFile(2);
    cout << "File opened" << endl;

    AS_Event       = new Ali_AS_Event();
    AS_Track       = new Ali_AS_Track();
   
    Tree_AS_Event  = NULL;
    Tree_AS_Event  = new TTree("Tree_AS_Event" , "AS_Events" );
    Tree_AS_Event  ->Branch("Tree_AS_Event_branch"  , "AS_Event", AS_Event );

    PostData(1,fListOfHistos);
    PostData(2,Tree_AS_Event);

    cout << "PostData called" << endl;

}



//________________________________________________________________________
Bool_t Ali_DarkMatter_ESD_analysis::NextEvent(Bool_t preload)
{
    fEventNoInFile++;
    //cout << "fEventNoInFile: " << fEventNoInFile << endl;
    //fDigitsLoadedFlag = kFALSE;


    if(preload)
    {
	//cout << "Preload"  << endl;
        //return ReadDigits();
        return kTRUE;
    }
    else
    {
	//cout << "No preload"  << endl;
	return kTRUE;
    }
}

//________________________________________________________________________
void Ali_DarkMatter_ESD_analysis::UserExec(Option_t *)
{
    //cout << "" << endl;
    //cout << "Analysis started" << endl;
    //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;


    //-----------------------------------------------------------------
    // IMPORTANT: call NextEvent() for book-keeping
    NextEvent();
    //if(fEventNoInFile > 50) return;
    //-----------------------------------------------------------------

    counter_events++;
    cout<<"counter events: "<<counter_events<<endl;

    //if(counter_events<250){return;}

    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    // prepare event data structures
    AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD)
    {
	printf("ERROR: fESD not available\n");
	return;
    }

    //-----------------------------------------------------------------
    // Check if TRD digits (raw data) are available for this ESD event
    //if(!ReadDigits()) return;
    //-----------------------------------------------------------------


#if 0
    //-----------------------------------------------------------------
    // Connect friends
    printf("Connect friends \n");
    fESD->SetESDfriend(esdFr);

    TTree* cTree = esdH->GetTree();
    cTree->AddFriend("esdFriendTree", esdFriendTreeFName.Data());
    cTree->SetBranchStatus("ESDfriend.", 1);
    esdFr = (AliESDfriend*)(fESD->FindListObject("AliESDfriend"));
    if (esdFr) cTree->SetBranchAddress("ESDfriend.", &esdFr);
    //-----------------------------------------------------------------
#endif

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man)
    {
        //Int_t run_id = man->GetRunFromPath(); // doesn't work
	//cout << "Got AliAnalysisManager, run_id: " << run_id << endl;
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	if(inputHandler)
	{
	    //cout << "Got AliInputEventHandler" << endl;
	    fPIDResponse = inputHandler->GetPIDResponse();
	}
    }
    //cout << "cent: " << fPIDResponse->GetCurrentCentrality() << endl;


    Int_t          N_tracks         = fESD ->GetNumberOfTracks();
    Int_t          N_TRD_tracks     = fESD ->GetNumberOfTrdTracks();
    Float_t        magF             = fESD ->GetMagneticField();
    const AliESDVertex* PrimVertex  = fESD ->GetPrimaryVertex();
    Int_t          RunNum           = fESD ->GetRunNumber();
    Double_t       T0zVertex        = fESD ->GetT0zVertex();
    AliCentrality* Centrality       = fESD ->GetCentrality();
    Double_t       MeanBeamIntAA    = fESD ->GetESDRun()->GetMeanIntensity(0,0);


    Int_t ncascades =  fESD-> GetNumberOfCascades();
    Int_t numberV0  =  fESD ->GetNumberOfV0s () ;


    cout<<"numberV0: "<<numberV0<<endl;
    /*
    if(numberV0>60000)
    {
        //counter_many_V0s++:
       // return;
    }
    */
    
    //printf("RunNum: %d, ncascades: %d , numberV0: %d /n ",RunNum,ncascades,numberV0);

    Double_t Sign_magnetic_field = (magF/fabs(magF));
    //cout << "Trigger: " <<  fESD->GetFiredTriggerClasses() << endl;
 


    // Fill event information
    AS_Event ->clearTrackList();
    AS_Event ->clearV0List();
    AS_Event ->setTriggerWord(fESD->GetFiredTriggerClasses());
    AS_Event ->setx(PrimVertex->GetX());
    AS_Event ->sety(PrimVertex->GetY());
    AS_Event ->setz(PrimVertex->GetZ());
    AS_Event ->setid(RunNum);
    AS_Event ->setN_tracks(N_tracks);
    AS_Event ->setBeamIntAA(MeanBeamIntAA);
    AS_Event ->setT0zVertex(T0zVertex);
    AS_Event ->setN_V0s(numberV0);
   

    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if(MultSelection)
    {
	// V0MEq, V0AEq, V0CEq, SPDTracklets

	AS_Event ->setcent_class_ZNA(MultSelection->GetMultiplicityPercentile("ZNA"));
	AS_Event ->setcent_class_ZNC(MultSelection->GetMultiplicityPercentile("ZNC"));
	AS_Event ->setcent_class_V0A(MultSelection->GetMultiplicityPercentile("V0A"));
	AS_Event ->setcent_class_V0C(MultSelection->GetMultiplicityPercentile("V0C"));
	AS_Event ->setcent_class_V0M(MultSelection->GetMultiplicityPercentile("V0M"));
	AS_Event ->setcent_class_CL0(MultSelection->GetMultiplicityPercentile("CL0"));
	AS_Event ->setcent_class_CL1(MultSelection->GetMultiplicityPercentile("CL1"));
	AS_Event ->setcent_class_SPD(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	AS_Event ->setcent_class_V0MEq(MultSelection->GetMultiplicityPercentile("V0MEq"));
	AS_Event ->setcent_class_V0AEq(MultSelection->GetMultiplicityPercentile("V0AEq"));
        AS_Event ->setcent_class_V0CEq(MultSelection->GetMultiplicityPercentile("V0CEq"));

    }

     Double_t *pN;
     Double_t *pP;
     pN = new Double_t[3];
     pP = new Double_t[3];


    for (Int_t V0_counter=0; V0_counter<numberV0; V0_counter++)
    {
        //get position of V0-----------------------
        Double_t x=0;
        Double_t y=0;
        Double_t z=0;
        AliESDv0 *V0=fESD->GetV0(V0_counter);
        V0->AliESDv0::GetXYZ(x,y,z);

        //cout<<""<<endl;
        //printf("V0 number: %d \n",V0_counter);
        //printf("x: %f,y: %f, z: %f \n",x,y,z);
        
        //-------------------------------

        //get impulse of particle N and P by AliESDv0 class
        Double_t pxN,pyN,pzN;
        V0->GetNPxPyPz(pxN,pyN,pzN);
        Double_t pxP,pyP,pzP;
        V0->GetPPxPyPz(pxP,pyP,pzP);
       //-------------------------------------


        //get track (by using index)-----------------------------
        Int_t indexN = V0->GetNindex();
        Int_t indexP = V0->GetPindex();

        //all_positive_track_ids.push_back(indexP);
        //all_negative_track_ids.push_back(indexN);
        //cout<<indexN<<endl;
        //cout<<indexP<<endl;
        AliESDtrack* trackN = fESD->AliESDEvent::GetTrack(indexN);
        AliESDtrack* trackP = fESD->AliESDEvent::GetTrack(indexP);


        Double_t momentumP = sqrt(pxP*pxP + pyP*pyP + pzP*pzP);
        Double_t momentumN = sqrt(pxN*pxN + pyN*pyN + pzN*pzN);

        //printf("trackidP %d, trackidN %d \n",indexP,indexN);
        //printf("momentum of positive particle: %f   momentum of negative particle: %f \n",momentumP,momentumN);
        //----------------------------------------------------------------

        //get impulse vector by using AliESDtrack class-------------------------
       
        trackN->GetInnerPxPyPz(pN);
        trackP->GetInnerPxPyPz(pP);

        //cout<<pxN<<pN[0]<<endl;
        //cout<<pyN<<pN[1]<<endl;
        //-------------------------------------------------


        //create element of Ali_AS_V0 class

        AS_V0  = AS_Event ->createV0();
        /*
        //use set functions
        AS_V0 -> setxyz(x,y,z);
        AS_V0 -> setNpxpypz(pxN,pyN,pzN);
        AS_V0 -> setPpxpypz(pxP,pyP,pzP);

        AS_V0 -> setdcaV0( V0->GetDcaV0Daughters() );


        //create tracks for positive and negative particle
        Ali_AS_Track* as_trackP = AS_V0->createTrack();
        Ali_AS_Track* as_trackN = AS_V0->createTrack();


        //fill tracks------------------------------------------------------------------
        //as_trackP  ->clearTRD_digit_list();
        //as_trackP  ->clearOfflineTrackletList();

        TLorentzVector TL_vec;
        Double_t Track_pT     = trackP ->Pt();
        Double_t Track_eta    = trackP ->Eta();
        Double_t Track_phi    = trackP ->Phi();
        TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,0.1349766);       //what to choose for M ?

        as_trackP  ->set_TLV_part(TL_vec);

        Float_t track_xy_impact,track_z_impact;
	trackP ->GetImpactParameters(track_xy_impact,track_z_impact);
        Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);
        Int_t    charge       = trackP ->Charge();
        as_trackP  ->setdca(((Double_t)charge)*track_total_impact);

        as_trackP -> setTRDSignal(trackP->GetTRDsignal());
        //as_trackP -> setTRDsumADC(-1);        //not found ok?
        as_trackP  ->setStatus(trackP->GetStatus());
        as_trackP  ->setNITScls(trackP->GetITSNcls());   
        as_trackP -> setTPCchi2(trackP->GetTPCchi2());
        as_trackP -> setTrack_length(trackP->GetIntegratedLength());
        as_trackP -> setNTPCcls(trackP ->GetTPCNcls());
        as_trackP -> setTOFsignal(trackP ->GetTOFsignal());
        as_trackP -> setTPCdEdx(trackP ->GetTPCsignal());

        FillHelix(trackP,magF);
        as_trackP ->setHelix(aliHelix.fHelix[0],aliHelix.fHelix[1],aliHelix.fHelix[2],aliHelix.fHelix[3],aliHelix.fHelix[4],aliHelix.fHelix[5],aliHelix.fHelix[6],aliHelix.fHelix[7],aliHelix.fHelix[8]);

        //track_PID
        // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID_P[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        // nSigma TPC
	Track_PID_P[0] = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kElectron);
	Track_PID_P[1] = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kMuon);
	Track_PID_P[2] = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kPion);
	Track_PID_P[3] = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kKaon);
	Track_PID_P[4] = fPIDResponse->NumberOfSigmasTPC(trackP,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID_P[5] = fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kElectron);
	Track_PID_P[6] = fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kMuon);
	Track_PID_P[7] = fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kPion);
	Track_PID_P[8] = fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kKaon);
        Track_PID_P[9] = fPIDResponse->NumberOfSigmasTOF(trackP,AliPID::kProton);


        as_trackP  ->setnsigma_e_TPC(Track_PID_P[0]);
	as_trackP  ->setnsigma_e_TOF(Track_PID_P[5]);
	as_trackP  ->setnsigma_pi_TPC(Track_PID_P[2]);
	as_trackP  ->setnsigma_pi_TOF(Track_PID_P[7]);
	as_trackP  ->setnsigma_K_TPC(Track_PID_P[3]);
	as_trackP  ->setnsigma_K_TOF(Track_PID_P[8]);
	as_trackP  ->setnsigma_p_TPC(Track_PID_P[4]);
	as_trackP  ->setnsigma_p_TOF(Track_PID_P[9]);

        as_trackP ->settrackid(indexP);
        as_trackN ->settrackid(indexN);
        //end of filling for P particle track-------------------------------------------------------------------

        //do the same for N particle track-----------------------------------------------------------------------
        //as_trackN  ->clearTRD_digit_list();
        //as_trackN  ->clearOfflineTrackletList();

        Track_pT     = trackN ->Pt();
        Track_eta    = trackN ->Eta();
        Track_phi    = trackN ->Phi();
        TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,0.1349766);       //what to choose for M ?

        as_trackN  ->set_TLV_part(TL_vec);

	trackN ->GetImpactParameters(track_xy_impact,track_z_impact);
        track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);
        charge       = trackN ->Charge();
        as_trackN  ->setdca(((Double_t)charge)*track_total_impact);

        as_trackN -> setTRDSignal(trackN->GetTRDsignal());
        //as_trackN -> setTRDsumADC(-1);        //not found ok?
        as_trackN  ->setStatus(trackN->GetStatus());
        as_trackN  ->setNITScls(trackN->GetITSNcls());
        as_trackN -> setTPCchi2(trackN->GetTPCchi2());
        as_trackN -> setTrack_length(trackN->GetIntegratedLength());
        as_trackN -> setNTPCcls(trackN ->GetTPCNcls());
        as_trackN -> setTOFsignal(trackN ->GetTOFsignal());
        as_trackN -> setTPCdEdx(trackN ->GetTPCsignal());

        FillHelix(trackN,magF);
        as_trackN ->setHelix(aliHelix.fHelix[0],aliHelix.fHelix[1],aliHelix.fHelix[2],aliHelix.fHelix[3],aliHelix.fHelix[4],aliHelix.fHelix[5],aliHelix.fHelix[6],aliHelix.fHelix[7],aliHelix.fHelix[8]);


        //track_PID
        // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID_N[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        // nSigma TPC
	Track_PID_N[0] = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kElectron);
	Track_PID_N[1] = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kMuon);
	Track_PID_N[2] = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kPion);
	Track_PID_N[3] = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kKaon);
	Track_PID_N[4] = fPIDResponse->NumberOfSigmasTPC(trackN,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID_N[5] = fPIDResponse->NumberOfSigmasTOF(trackN,AliPID::kElectron);
	Track_PID_N[6] = fPIDResponse->NumberOfSigmasTOF(trackN,AliPID::kMuon);
	Track_PID_N[7] = fPIDResponse->NumberOfSigmasTOF(trackN,AliPID::kPion);
	Track_PID_N[8] = fPIDResponse->NumberOfSigmasTOF(trackN,AliPID::kKaon);
        Track_PID_N[9] = fPIDResponse->NumberOfSigmasTOF(trackN,AliPID::kProton);


        as_trackN  ->setnsigma_e_TPC(Track_PID_N[0]);
	as_trackN  ->setnsigma_e_TOF(Track_PID_N[5]);
	as_trackN  ->setnsigma_pi_TPC(Track_PID_N[2]);
	as_trackN  ->setnsigma_pi_TOF(Track_PID_N[7]);
	as_trackN  ->setnsigma_K_TPC(Track_PID_N[3]);
	as_trackN  ->setnsigma_K_TOF(Track_PID_N[8]);
	as_trackN  ->setnsigma_p_TPC(Track_PID_N[4]);
        as_trackN  ->setnsigma_p_TOF(Track_PID_N[9]);

        //end of filling for N particle track-------------------------------------------------------------------
        */
    }
    cout<<"end of V0 loop"<<endl;
    cout<<""<<endl;
    //end of V0 loop

    //----------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------

    //cout << "Event number: " << fEventNoInFile << ", event number with TRD digits: " << N_good_events << endl;
    //printf("Event number: %d, N_tracks: %d, cent(V0M): %f , cent(CL0): %f \n",fEventNoInFile,N_tracks,MultSelection->GetMultiplicityPercentile("SPDTracklets"));
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    Double_t TPC_TRD_matching_window = 10.0;  // Matching window between TRD digits (pad position) and TPC track in cm
    Double_t TPC_min_radius_plot     = 114.0/2.0; // Minimum radius for TPC track to be plotted
    Double_t TPC_max_radius_plot     = 368.0; // Maximum radius for TPC track to be plotted
    Double_t TPC_radius_scan         = 368.0 - 30.0; // Radius for TPC track to be used for first match with TRD digits

    if(N_tracks == 0)
    {
	// Skip empty event
	return;
    }

    //printf("There are %d tracks in this event\n", N_tracks);
    //printf("There are %d TRD tracks in this event\n", N_TRD_tracks);
    //printf("There are %d TRD tracklets in this event\n", N_TRD_tracklets);

    Int_t N_good_tracks = 0;

    //-----------------------------------------------------------------
    // Track loop
    //cout << "" << endl;
    //cout << "-----------------------------------------------------------------" << endl;
    //cout << "Start matching " << N_tracks << " TPC tracks with " << TV3_TRD_hits_middle.size() << " TRD pads" << endl;
    N_good_tracks = 0;


    for(Int_t iTracks = 0; iTracks < N_tracks; iTracks++)
    {
	//---------------------------------------------------------------
	// Gather track information

	// We always want the ESD track
	AliESDtrack* track = fESD->GetTrack(iTracks);
	if(!track)
	{
	    printf("ERROR: Could not receive track %d\n", iTracks);
	    continue;
	}


        if(!EsdTrackCuts->AcceptTrack(track)) continue;

        TBits tbits_fit  = track->GetTPCFitMap();
        //Int_t countbits = tbit.CountBits();
        //cout<<"bits1 fit: : "<<countbits<<endl;
        //Int_t nbits = tbit.GetNbits();
        //cout<<"nbits fit: : "<<nbits<<endl;

        TBits tbits_shared = track->GetTPCSharedMap();
        //Int_t countbits_sh = shared.CountBits();
        //cout<<"bits1 shared:: "<<countbits_sh<<endl;
        //Int_t nbits_sh = shared.GetNbits();
        //cout<<"nbits shared:: "<<nbits_sh<<endl;

        Double_t TRD_signal   = track ->GetTRDsignal(); // truncated mean signal?
        Double_t Track_pT     = track ->Pt();
        Double_t Track_p      = track ->P();
        Double_t p_vec[3];
        track->GetPxPyPz(p_vec);
        Int_t    charge       = track ->Charge();
        Double_t Track_phi    = track ->Phi();
	Double_t Track_theta  = track ->Theta();
	Double_t Track_eta    = track ->Eta();
	Double_t TPC_chi2     = track ->GetTPCchi2();
	Double_t TPC_signal   = track ->GetTPCsignal(); // dE/dx?
	Double_t TOF_signal   = track ->GetTOFsignal(); // time-of-flight?
        Double_t Track_length = track ->GetIntegratedLength();
        UShort_t N_TPC_cls    = track ->GetTPCNcls();
        Int_t   trackid       = track ->GetID();

        

        Double_t r[3];
        Double_t p[3];
        track -> GetInnerXYZ(r);
        track -> GetInnerPxPyPz(p);
        Double_t momentum = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

        //printf("track number: %d, charge:  %d, momentum: %f \n",iTracks, charge, momentum);

	TLorentzVector TLV_p_vec;
	Double_t p_vec_energy = TMath::Sqrt(p_vec[0]*p_vec[0] + p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2] + 0.938*0.938);
	TLV_p_vec.SetPxPyPzE(p_vec[0],p_vec[1],p_vec[2],p_vec_energy);
	//cout << "TLV_p_vec.P: " << TLV_p_vec.P() << ", P: " << Track_p << ", TLV_p_vec.Theta: " << TLV_p_vec.Theta() << ", Theta: " << Track_theta
	//<< ", TLV_p_vec.Phi: " << TLV_p_vec.Phi() << ", phi: " << Track_phi  << endl;


	ULong_t status = track->GetStatus();
	Int_t ITS_refit = 0;
        Int_t TPC_refit = 0;
        Int_t track_status = 0;
	if(((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit))
	{
	    ITS_refit = 1;
	    track_status |= 1 << 0; // setting bit 0 to 1
	}
	if(((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit))
	{
	    TPC_refit = 1;
	    track_status |= 1 << 1; // setting bit 1 to 1
	}


          // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        // nSigma TPC
	Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
	Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kMuon);
	Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
	Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
	Track_PID[4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
	Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kMuon);
	Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
	Track_PID[8] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
	Track_PID[9] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);


	Float_t track_xy_impact,track_z_impact;
	track->GetImpactParameters(track_xy_impact,track_z_impact);
	Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);


	Double_t TRD_ADC_bin_width = 100.0;

	//-------------------
	Int_t N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(track ->HasPointOnITSLayer(i_ITS_layer))
	    {
		N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
	    }
	}
	//-------------------

	TLorentzVector TL_vec;
        TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,0.1349766);


	AS_Track  = AS_Event ->createTrack();


        AS_Track  ->set_TLV_part(TL_vec);
	AS_Track  ->setdca(((Double_t)charge)*track_total_impact);
	AS_Track  ->setnsigma_e_TPC(Track_PID[0]);
	AS_Track  ->setnsigma_e_TOF(Track_PID[5]);
	AS_Track  ->setnsigma_pi_TPC(Track_PID[2]);
	AS_Track  ->setnsigma_pi_TOF(Track_PID[7]);
	AS_Track  ->setnsigma_K_TPC(Track_PID[3]);
	AS_Track  ->setnsigma_K_TOF(Track_PID[8]);
	AS_Track  ->setnsigma_p_TPC(Track_PID[4]);
	AS_Track  ->setnsigma_p_TOF(Track_PID[9]);
	AS_Track  ->setTRDSignal(TRD_signal);
	AS_Track  ->setNTPCcls(N_TPC_cls);
	AS_Track  ->setNITScls(N_ITS_cls);
	AS_Track  ->setStatus(track_status);
	AS_Track  ->setTPCchi2(TPC_chi2);
	AS_Track  ->setTPCdEdx(TPC_signal);
	AS_Track  ->setTOFsignal(TOF_signal);
        AS_Track  ->setTrack_length(Track_length);
        AS_Track  ->settrackid(trackid);

        AS_Track -> settbitsfit(tbits_fit);
        AS_Track -> settbitsshared(tbits_shared);

#if 1
        //--------------------------------------
        //printf("Loop over friend tracks \n");
        const AliESDfriendTrack *trkFr = track->GetFriendTrack();
        if(!trkFr)
        {
            continue;
        }

        Int_t ESDtrackID       = trkFr ->GetESDtrackID();
        Float_t one_over_p     = trkFr ->Get1P();
        //Int_t N_MaxTPCclusters = trkF ->GetMaxTPCcluster();
        //printf("Friend track available, ESDtrackID: %d, one_over_p: %4.3f \n",ESDtrackID,one_over_p);
        const AliTRDtrackV1 *trdTrack = 0;
        const TObject *calibObject = 0;

        for(Int_t idx = 0; (calibObject = trkFr->GetCalibObject(idx)); ++idx)
        {
            //printf("idx: %d \n",idx);
            if(calibObject->IsA() != AliTRDtrackV1::Class())
            {
                continue;
            }
            trdTrack = (AliTRDtrackV1*) calibObject;
        }

        //printf("Calib object available \n");

       
        //--------------------------------------
#endif



#if 0
        //---------------------
        //---------------------
#endif

	//Helix
        FillHelix(track,magF);

        AS_Track ->setHelix(aliHelix.fHelix[0],aliHelix.fHelix[1],aliHelix.fHelix[2],aliHelix.fHelix[3],aliHelix.fHelix[4],aliHelix.fHelix[5],aliHelix.fHelix[6],aliHelix.fHelix[7],aliHelix.fHelix[8]);

      

    } // End of TPC track loop
    cout << "Tracks matched" << endl;
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    Tree_AS_Event ->Fill();

    //Ali_AS_V0* AS_V0;

    for (Int_t V0_counter=0; V0_counter<numberV0; V0_counter++)
    {
        //if(counter_events>280) {cout<<"V0 counter2: "<<V0_counter<<endl;}
        AS_V0 =  AS_Event->getV0(V0_counter);
        AS_V0 -> clearTrackList();
        //AS_V0 -> clearV0List();
        //AS_V0->~Ali_AS_V0();
    }


    int Ntracks = AS_Event->getN_tracks();
    //cout<<"Ntracks: "<<Ntracks<<endl;

    /*
    for(Int_t iTracks = 0; iTracks < Ntracks; iTracks++)
    {
        //cout<<iTracks<<endl;
        Ali_AS_Track* AS_Track = AS_Event->getTrack(iTracks);
        //AS_Track->~Ali_AS_Track();
        //if(AS_Track==NULL){continue;}
        delete AS_Track;
    }
    */

    //AS_Event->clearTrackList();
   // AS_Event ->clearV0List();

    N_good_events++;

}



//________________________________________________________________________
void Ali_DarkMatter_ESD_analysis::Terminate(Option_t *)
{
    cout << "In terminate" << endl;
}


//________________________________________________________________________
void Ali_DarkMatter_ESD_analysis::FillHelix(AliESDtrack* track_in, Double_t magF_in)
{
    //-------------------
    // Get helix
    // Track parametrization:
    // https://www.physi.uni-heidelberg.de/~sma/alice/LukasLayer_bachelor.pdf
    Double_t alpha, alpha_alt, x_param, p_param[5]; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
    track_in->GetExternalParameters(x_param,p_param); // Those are the parameters used for vertexing
    alpha_alt=track_in->GetAlpha();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "A i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //track_in->GetOuterExternalParameters(alpha,x_param,p_param); //
    //alpha_alt = alpha;
    Double_t fX     = track_in->GetX();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "B i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //cout << "x_param: " << x_param << endl;

    //const AliExternalTrackParam* TPC_track_param = track_in->GetOuterParam();
    //const Double_t* p_param_b = TPC_track_param->GetParameter();
    //const Double_t* p_param_b = track_in->GetOuterParam()->GetParameter();
    //const Double_t* p_param_b = track_in->GetInnerParam()->GetParameter();
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    p_param[i_param] = p_param_b[i_param];
    //}
    //x_param = fX;

    //-------------------
    // Correct way of filling aliHelix from
    // http://personalpages.to.infn.it/~puccio/htmldoc/src/AliHelix.cxx.html#PalT1E
    // line 52

    // CONCLUSION: AliTracker::GetBz() is not identical to GetC(magF_in), magnetic field is slightly different but it doesn't matter...

    Double_t fHelix_alt[9];
    Double_t x_alt,cs_alt,sn_alt;
    //track_in->GetExternalParameters(x_alt,fHelix_alt); // Those are the parameters used for vertexing
    x_alt = x_param;

    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
	fHelix_alt[i_param] = p_param[i_param];
    }

    //cout << "alpha: " << alpha << ", alpha_alt: " << alpha_alt << endl;

    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...

    // kB2C=-0.299792458e-3; // from /AliRoot/STEER/STEERBase/AliVParticle.h
    //Double_t kB2C_test =-0.299792458e-3;
    //Double_t GetC(Double_t b) const
    //{return fP[4]*b*kB2C;}

    //Double_t par4test = fHelix_alt[4]*magF_in*kB2C_test;
    fHelix_alt[4] = track_in->GetC(magF_in);
    //fHelix_alt[4] = par4test; // take the one with the magnetic field directly from the ESD file

    cs_alt = TMath::Cos(alpha_alt);
    sn_alt = TMath::Sin(alpha_alt);

    Double_t xc_alt, yc_alt, rc_alt;
    rc_alt  =  1/fHelix_alt[4];
    xc_alt  =  x_alt-fHelix_alt[2]*rc_alt;
    Double_t dummy = 1-(x_alt-xc_alt)*(x_alt-xc_alt)*fHelix_alt[4]*fHelix_alt[4];
    yc_alt  =  fHelix_alt[0]+TMath::Sqrt(dummy)/fHelix_alt[4];

    fHelix_alt[6] = xc_alt*cs_alt - yc_alt*sn_alt;
    fHelix_alt[7] = xc_alt*sn_alt + yc_alt*cs_alt;
    fHelix_alt[8] =  TMath::Abs(rc_alt);
    //
    //
    fHelix_alt[5]=x_alt*cs_alt - fHelix_alt[0]*sn_alt;            // x0
    fHelix_alt[0]=x_alt*sn_alt + fHelix_alt[0]*cs_alt;            // y0
    fHelix_alt[2]=TMath::ATan2(-(fHelix_alt[5]-fHelix_alt[6]),fHelix_alt[0]-fHelix_alt[7]); // phi0
    if (fHelix_alt[4]>0) fHelix_alt[2]-=TMath::Pi();
    fHelix_alt[5]   = fHelix_alt[6];
    fHelix_alt[0]   = fHelix_alt[7];

    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
	aliHelix.fHelix[i_param] = fHelix_alt[i_param];
    }
    //-------------------
}



//________________________________________________________________________
void Ali_DarkMatter_ESD_analysis::FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB)
{
    // V1.0
    Float_t pA[2] = {path_initA,path_initB}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    for(Int_t r = 0; r < 2; r++)
    {
	Double_t helix_point[3];
	helixA.Evaluate(pA[r],helix_point);
	testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[r]
	distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 30.0;
    while(fabs(scale_length) > 0.1 && loopcounter < 100) // stops when the length is too small
    {
	//cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
	//    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
	//    << ", d[1] = " << distarray[1] << ", flip = " << flip
	//    << ", scale_length = " << scale_length << endl;
	if(distarray[0] > distarray[1])
	{
	    if(loopcounter != 0)
	    {
		if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
		else scale = 0.7; // go on in this direction but only by the way * 0.7
	    }
	    scale_length = (pA[1]-pA[0])*scale; // the next length interval
	    pA[0]     = pA[1] + scale_length; // the new path

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[0],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[0]
	    distarray[0] = (testA-space_vec).Mag(); // new dca
	    flip = 1.0;
	}
	else
	{
	    if(loopcounter != 0)
	    {
		if(flip == -1.0) scale = 0.4;
		else scale = 0.7;
	    }
	    scale_length = (pA[0]-pA[1])*scale;
	    pA[1]     = pA[0] + scale_length;

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[1],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[1]
	    distarray[1] = (testA-space_vec).Mag();
	    flip = -1.0;
	}
	loopcounter++;
    }

    if(loopcounter >= 100) cout << "WARNING: FindDCAHelixPoint exceeded maximum of 100 loops" << endl;

    if(distarray[0] < distarray[1])
    {
	pathA = pA[0];
	dcaAB = distarray[0];
    }
    else
    {
	pathA = pA[1];
	dcaAB = distarray[1];
    }
}


//----------------------------------------------------------------------------------------

void Ali_DarkMatter_ESD_analysis::func_tail_cancellation(Short_t *arr, Int_t nexp)
{
    // Tail cancellation by deconvolution for PASA v4 TRF
    //

    Int_t fBaseline = 10;

    Float_t rates[2];
    Float_t coefficients[2];

    // Initialization (coefficient = alpha, rates = lambda)
    Float_t r1 = 1.0;
    Float_t r2 = 1.0;
    Float_t c1 = 0.5;
    Float_t c2 = 0.5;

    r1 = 1.156;
    r2 = 0.130;
    c1 = 0.066;
    c2 = 0.000;

    if (nexp == 1) {   // 1 Exponentials
        r1 = 1.156;
        r2 = 0.130;
        c1 = 0.066;
        c2 = 0.000;
    }
    if (nexp == 2) {   // 2 Exponentials
        //Double_t par[4];
        //fReconstructor->GetRecoParam()->GetTCParams(par);
        //r1 = par[0];//1.156;
        //r2 = par[1];//0.130;
        //c1 = par[2];//0.114;
        //c2 = par[3];//0.624;

        r1 = 1.156;
        r2 = 0.130;
        c1 = 0.114;
        c2 = 0.624;
    }

    coefficients[0] = c1;
    coefficients[1] = c2;

    Double_t dt = 0.1;

    rates[0] = TMath::Exp(-dt/(r1));
    rates[1] = (nexp == 1) ? .0 : TMath::Exp(-dt/(r2));

    Float_t reminder[2] = { .0, .0 };
    Float_t correction = 0.0;
    Float_t result     = 0.0;

    for (int i = 0; i < 24; i++) {

        result = arr[i] - correction - fBaseline;    // No rescaling
        arr[i] = (Short_t)(result + fBaseline + 0.5f);

        correction = 0.0;
        for (int k = 0; k < 2; k++) {
            correction += reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
        }
    }
}
//----------------------------------------------------------------------------------------

