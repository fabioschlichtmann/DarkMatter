
#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"



//________________________________________________________________________
void FindDCAHelixPoint(TVector3 space_vec, Ali_AS_Track* helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB)
{
    // V1.0
    Float_t pA[2] = {path_initA,path_initB}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    for(Int_t r = 0; r < 2; r++)
    {
	Double_t helix_point[3];
	helixA->Evaluate(pA[r],helix_point);
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
	    helixA->Evaluate(pA[0],helix_point);
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
	    helixA->Evaluate(pA[1],helix_point);
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



//----------------------------------------------------------------------------------------
class Ali_Draw_SEvent
{
private:
    TChain* input_SE;
    TString TRD_DM_TREE   = "Tree_AS_DM_particle";
    TString TRD_DM_BRANCH = "Tree_AS_DM_branch";
    Long64_t file_entries_total;
    Long64_t N_Events;
    Ali_AS_Track*           AS_Track;
    Ali_AS_DM_particle*     AS_SParticle;

    TString HistName;

    TEveLine*     TEveLine_beam_axis = NULL;
    TEveLine*     TPL3D_helix        = NULL;
    TEvePointSet* TEvePoint_vertex   = NULL;
    vector<TEveLine*>     vec_TPL3D_helix;
    vector<TEveLine*>     vec_TPL3D_helix_inner;
    vector<TEveLine*>     vec_TPL3D_helix_hull;
    vector<TEvePointSet*> vec_TEvePoint_vertices;


public:
    Ali_Draw_SEvent();
    ~Ali_Draw_SEvent();
    void Init_tree(TString SEList);
    void Draw_event(Long64_t i_event);
    int Cutevent(Long64_t i_event);
    Long64_t getfileentries();

    ClassDef(Ali_Draw_SEvent, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Ali_Draw_SEvent::Ali_Draw_SEvent()
{
    // constructor
    TEveManager::Create();

    TPL3D_helix = new TEveLine();
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-30.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,30.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainColor(kGray+2);
    gEve->AddElement(TEveLine_beam_axis);

    TEvePoint_vertex = new TEvePointSet();
    vec_TPL3D_helix.resize(6);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_Draw_SEvent::Init_tree(TString SEList)
{
    printf("Ali_Draw_SEvent::Init_tree \n");
    //TString pinputdir = "./";
    TString pinputdir = "/misc/alidata120/alice_u/schlichtmann/out/";

    AS_Track     = new Ali_AS_Track();
    AS_SParticle = new Ali_AS_DM_particle();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( TRD_DM_TREE.Data(), TRD_DM_TREE.Data() );
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
                    input_SE ->AddFile(addfile.Data(),-1, TRD_DM_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( TRD_DM_BRANCH, &AS_SParticle );
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
double calculate_m_squared_by_TOF(Ali_AS_Track* track_in_func)
{

        Float_t TPCdEdx   = track_in_func->getTPCdEdx();
        Float_t tofsignal = track_in_func->getTOFsignal();
        Float_t dca       = track_in_func->getdca();
        Float_t tracklength = track_in_func->getTrack_length();
        int charge;

        if(dca>0){charge = 1;}
        else {charge = -1;}

        TLorentzVector tlv = track_in_func->get_TLV_part();
        double momentum = tlv.P();

        //printf("dEdx: %f, momentum: %f, dca: %f, charge: %d, tofsignal: %f \n"
          //     ,TPCdEdx, momentum,dca, charge, tofsignal);
        if(tofsignal>99990){return -1;}

        double velocity = tracklength/tofsignal;

        //printf("velocity: %f \n", velocity);

        velocity = velocity * 1e10;

         //printf("velocity: %f \n", velocity);

        double speed_of_light_SI = 299792458.;

        velocity = velocity /speed_of_light_SI;  //now in units of c

        // printf("velocity: %f \n", velocity);

        double gamma_squared = 1. / (1-velocity*velocity) ;
        //printf("momentum: %f, gamma: %f, velocity: %f \n",momentum,gamma,velocity);

        //m^2 =  ( p/ (gamma * v) )^2
        double m_squared = ( momentum / velocity)  *  ( momentum /velocity) * 1./gamma_squared ;
        return m_squared;
}


int Ali_Draw_SEvent::Cutevent(Long64_t i_event)
{
    if (!input_SE->GetEntry( i_event )) return -1;
    AS_Track = AS_SParticle ->getTrack(0);
    double mass1 = calculate_m_squared_by_TOF(AS_Track);
    //cout<<mass1<<endl;
    //return 0;

    AS_Track = AS_SParticle ->getTrack(2);
    double mass2 = calculate_m_squared_by_TOF(AS_Track);
    //cout<<mass2<<endl;
    //printf("event number: %lld, mass first K: %f, mass second K: %f \n",i_event,mass1,mass2);
    if(mass1>0.2 && mass1<0.35){cout<<i_event<<endl;}
    if(mass2>0.2 && mass2<0.35){cout<<i_event<<endl;}
    if(mass1>0.2 && mass1<0.35 && mass2>0.2 && mass2<0.35 ){cout<<"Jackpot"<<endl;}
    return 0;


}

Long64_t Ali_Draw_SEvent::getfileentries()
{
    return file_entries_total;
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void Ali_Draw_SEvent::Draw_event(Long64_t i_event)
{
    printf("Ali_Draw_SEvent::Draw_event \n");

    //cout<<"1"<<endl;
    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event
    //cout<<"2"<<endl;


    TVector3 TV3_primVertex = AS_SParticle ->get_primVertex();
    TVector3 TV3_S1Vertex   = AS_SParticle ->get_S1Vertex();
    TVector3 TV3_S2Vertex   = AS_SParticle ->get_S2Vertex();
    TVector3 TV3_S3Vertex   = AS_SParticle ->get_S3Vertex();
    TVector3 TV3_DirSV1     = AS_SParticle ->get_DirSV1();
    TVector3 TV3_DirSV2     = AS_SParticle ->get_DirSV2();
    TVector3 TV3_DirSV3     = AS_SParticle ->get_DirSV3();
    Int_t    N_V0s          = AS_SParticle ->getN_V0s();
    UShort_t N_tracks       = AS_SParticle ->getNumTracks();

    //track 0: pion1
    //track 1: pion2
    //track 2: SV2 track
    //track 3: SV2 track
    //track 4: SV3 track
    //track 5: SV3 track

    //printf("N_tracks: %d \n",N_tracks);

    Int_t SV1_color = kRed;
    Int_t SV2_color = kGreen;
    Int_t SV3_color = kBlue;


    // Primary vertex
    TEvePoint_vertex ->SetPoint(0,TV3_primVertex.X(),TV3_primVertex.Y(),TV3_primVertex.Z());
    vec_TEvePoint_vertices.push_back((TEvePointSet*)TEvePoint_vertex->Clone());
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerSize(3);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerStyle(20);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerColor(kGray+1);
    gEve->AddElement(vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1]);

    // S1 vertex
    TEvePoint_vertex ->SetPoint(0,TV3_S1Vertex.X(),TV3_S1Vertex.Y(),TV3_S1Vertex.Z());
    vec_TEvePoint_vertices.push_back((TEvePointSet*)TEvePoint_vertex->Clone());
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerSize(3);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerStyle(20);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerColor(SV1_color);
    gEve->AddElement(vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1]);


    // S2 vertex
    TEvePoint_vertex ->SetPoint(0,TV3_S2Vertex.X(),TV3_S2Vertex.Y(),TV3_S2Vertex.Z());
    vec_TEvePoint_vertices.push_back((TEvePointSet*)TEvePoint_vertex->Clone());
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerSize(3);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerStyle(20);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerColor(SV2_color);
    gEve->AddElement(vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1]);

    // S3 vertex
    TEvePoint_vertex ->SetPoint(0,TV3_S3Vertex.X(),TV3_S3Vertex.Y(),TV3_S3Vertex.Z());
    vec_TEvePoint_vertices.push_back((TEvePointSet*)TEvePoint_vertex->Clone());
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerSize(3);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerStyle(20);
    vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1] ->SetMarkerColor(SV3_color);
    gEve->AddElement(vec_TEvePoint_vertices[vec_TEvePoint_vertices.size()-1]);



    TEveLine* TEL_conn_prim_S1 = new TEveLine();
    TEL_conn_prim_S1 ->SetNextPoint(TV3_primVertex.X(),TV3_primVertex.Y(),TV3_primVertex.Z());
    TEL_conn_prim_S1 ->SetNextPoint(TV3_S1Vertex.X(),TV3_S1Vertex.Y(),TV3_S1Vertex.Z());
    TEL_conn_prim_S1 ->SetName("prim to SV1");
    TEL_conn_prim_S1 ->SetLineStyle(9);
    TEL_conn_prim_S1 ->SetLineWidth(3);
    TEL_conn_prim_S1 ->SetMainColor(kGray);
    gEve->AddElement(TEL_conn_prim_S1);


    TEveLine* TEL_conn_S1_S2 = new TEveLine();
    TEL_conn_S1_S2 ->SetNextPoint(TV3_S1Vertex.X(),TV3_S1Vertex.Y(),TV3_S1Vertex.Z());
    TEL_conn_S1_S2 ->SetNextPoint(TV3_S2Vertex.X(),TV3_S2Vertex.Y(),TV3_S2Vertex.Z());
    TEL_conn_S1_S2 ->SetName("SV1 to SV2");
    TEL_conn_S1_S2 ->SetLineStyle(9);
    TEL_conn_S1_S2 ->SetLineWidth(3);
    TEL_conn_S1_S2 ->SetMainColor(kGray);
    gEve->AddElement(TEL_conn_S1_S2);

    TEveLine* TEL_conn_S1_S3 = new TEveLine();
    TEL_conn_S1_S3 ->SetNextPoint(TV3_S1Vertex.X(),TV3_S1Vertex.Y(),TV3_S1Vertex.Z());
    TEL_conn_S1_S3 ->SetNextPoint(TV3_S3Vertex.X(),TV3_S3Vertex.Y(),TV3_S3Vertex.Z());
    TEL_conn_S1_S3 ->SetName("S1 to S3");
    TEL_conn_S1_S3 ->SetLineStyle(9);
    TEL_conn_S1_S3 ->SetLineWidth(3);
    TEL_conn_S1_S3 ->SetMainColor(kGray);
    gEve->AddElement(TEL_conn_S1_S3);

    Double_t scale_dir = 20.0;


    TEveArrow* TEA_dir_SV2 = new TEveArrow();
    TEA_dir_SV2 ->SetOrigin(TV3_S2Vertex.X(),TV3_S2Vertex.Y(),TV3_S2Vertex.Z());
    TEA_dir_SV2 ->SetVector(scale_dir*TV3_DirSV2.X(),scale_dir*TV3_DirSV2.Y(),scale_dir*TV3_DirSV2.Z());
    TEA_dir_SV2 ->SetMainColor(SV2_color);
    gEve->AddElement(TEA_dir_SV2);

    TEveArrow* TEA_dir_SV3 = new TEveArrow();
    TEA_dir_SV3 ->SetOrigin(TV3_S3Vertex.X(),TV3_S3Vertex.Y(),TV3_S3Vertex.Z());
    TEA_dir_SV3 ->SetVector(scale_dir*TV3_DirSV3.X(),scale_dir*TV3_DirSV3.Y(),scale_dir*TV3_DirSV3.Z());
    TEA_dir_SV3 ->SetMainColor(SV3_color);
    gEve->AddElement(TEA_dir_SV3);


    AS_Track = AS_SParticle ->getTrack(2);
    TLorentzVector tlv = AS_Track->get_TLV_part();

    //cout<<tlv[0]<<endl;

    /*
    TEveArrow* TEA_dir_track0 = new TEveArrow();
    TEA_dir_track0 ->SetOrigin(TV3_S1Vertex.X(),TV3_S1Vertex.Y(),TV3_S1Vertex.Z());
    TEA_dir_track0 ->SetVector(scale_dir*tlv[0],scale_dir*tlv[1],scale_dir*tlv[2]);
    TEA_dir_track0 ->SetMainColor(kBlue);
    gEve->AddElement(TEA_dir_track0);
    */
    Double_t track_pos[3];
    Double_t radius_helix;

    Int_t arr_line_style[6] = {1,9,1,9,1,9};
    Int_t arr_line_color[6] = {SV1_color,SV1_color,SV2_color,SV2_color,SV3_color,SV3_color};


    for(Int_t i_track = 0; i_track < 6; i_track++)
    //for(Int_t i_track = 0; i_track < 3; i_track++)
    //for(Int_t i_track = 0; i_track < 2; i_track++)
    {
        AS_Track = AS_SParticle ->getTrack(i_track);
        Int_t track_id = AS_Track ->gettrackid();

        printf("i_track: %d, track_id: %d \n",i_track,track_id);

        TLorentzVector tlv = AS_Track->get_TLV_part();
        cout<<"tlv: "<<tlv[0]<<" "<<tlv[1]<<" "<<tlv[2]<<endl;

        double dedx = AS_Track->getTPCdEdx();
        cout<<"dEdx: "<<dedx<<endl;


        Float_t path_closest_to_point = 0;
        Float_t dca_closest_to_point  = 0;
        Float_t path_initA = 0.0;
        Float_t path_initB = 30.0;

        //calculate dca from vertex to particle track
        if(i_track == 0 || i_track == 1) FindDCAHelixPoint(TV3_S1Vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
        if(i_track == 2 || i_track == 3) FindDCAHelixPoint(TV3_S2Vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
        if(i_track == 4 || i_track == 5) FindDCAHelixPoint(TV3_S3Vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);


        vec_TPL3D_helix[i_track] = new TEveLine();
        for(Double_t track_path = path_closest_to_point; track_path < 1000; track_path += 1.0)
        //for(Double_t track_path = 0; track_path < 1000; track_path += 1.0)
        {
            //cout<<"track_path"<<track_path<<endl;
            AS_Track ->Evaluate(track_path,track_pos);
            radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
            if(radius_helix > 200.0) break;
            if(fabs(track_pos[2]) > 200.0) break;
            vec_TPL3D_helix[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
        }

        HistName = "track ";
        HistName += i_track;
        vec_TPL3D_helix[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix[i_track]    ->SetLineStyle(arr_line_style[i_track]);
        vec_TPL3D_helix[i_track]    ->SetLineWidth(3);
        vec_TPL3D_helix[i_track]    ->SetMainColor(arr_line_color[i_track]);
        vec_TPL3D_helix[i_track]    ->SetMainAlpha(1.0);

        gEve->AddElement(vec_TPL3D_helix[i_track]);
    }


    gEve->Redraw3D(kTRUE);

}
//----------------------------------------------------------------------------------------



