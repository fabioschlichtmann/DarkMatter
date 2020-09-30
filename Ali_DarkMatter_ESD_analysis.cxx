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
#include "TMath.h"
//------------------------

#include "Ali_DarkMatter_ESD_analysis.h"

#include "vertex_modified.h"


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
ClassImp(Ali_AS_NUCLEV)
ClassImp(Ali_AS_DM_particle)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_DarkMatter_ESD_analysis)

    //________________________________________________________________________
    Ali_DarkMatter_ESD_analysis::Ali_DarkMatter_ESD_analysis(const char *name)
    : AliAnalysisTaskSE(name),
    AS_Event(0),AS_V0(0),AS_Track(0),AS_NUCLEV(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0),
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

bool check_if_int_is_in_vector(int a, vector<int> vec)
{
    for(int i=0;i<vec.size();i++)
    {
        if(a == vec[i]) {return 1;}

    }
    return 0;

}

bool check_if_two_vectors_have_same_element(vector<int> vec1, vector<int> vec2)
{
    for(int i=0;i<vec1.size();i++)
    {
        for(int j=0;j<vec2.size();j++)
        {
            if(vec1[i]==vec2[j]){return 1;}
        }

    }
    return 0;

}

bool check_if_value_is_doppelt_in_vector(vector<int> vec)
{
    for(int i=0;i<(Int_t)vec.size();i++)
    {
        int counter = 0;
        for(int j=0;j<(Int_t)vec.size();j++)
        {
            if(vec[i]==vec[j]){counter++;}
        }
        if(counter>1){return true;}
    }
    return false;
}

void copy_track_params(Ali_AS_Track* track_in, Ali_AS_Track* track_out)
{
    track_out ->setnsigma_e_TPC(track_in ->getnsigma_e_TPC());
    track_out ->setnsigma_e_TOF(track_in ->getnsigma_e_TOF());
    track_out ->setnsigma_pi_TPC(track_in ->getnsigma_pi_TPC());
    track_out ->setnsigma_pi_TOF(track_in ->getnsigma_pi_TOF());
    track_out ->setnsigma_K_TPC(track_in ->getnsigma_K_TPC());
    track_out ->setnsigma_K_TOF(track_in ->getnsigma_K_TOF());
    track_out ->setnsigma_p_TPC(track_in ->getnsigma_p_TPC());
    track_out ->setnsigma_p_TOF(track_in ->getnsigma_p_TOF());
    track_out ->setTRDSignal(track_in ->getTRDSignal());
    track_out ->setTRDsumADC(track_in ->getTRDsumADC());
    track_out ->setdca(track_in ->getdca());
    track_out ->set_TLV_part(track_in ->get_TLV_part());
    track_out ->setNTPCcls(track_in ->getNTPCcls());
    track_out ->setNTRDcls(track_in ->getNTRDcls());
    track_out ->setNITScls(track_in ->getNITScls());
    track_out ->setStatus(track_in ->getStatus());
    track_out ->setTPCchi2(track_in ->getTPCchi2());
    for(Int_t i_layer = 0; i_layer < 6; i_layer++){track_out ->setTRD_layer(i_layer, track_in ->getTRD_layer(i_layer));}
    track_out ->setimpact_angle_on_TRD(track_in ->getimpact_angle_on_TRD());
    track_out ->setTPCdEdx(track_in ->getTPCdEdx());
    track_out ->setTOFsignal(track_in ->getTOFsignal());
    track_out ->setTrack_length(track_in ->getTrack_length());
    track_out ->setHelix(track_in->getHelix_param(0),track_in->getHelix_param(1),track_in->getHelix_param(2),track_in->getHelix_param(3),track_in->getHelix_param(4),track_in->getHelix_param(5),track_in->getHelix_param(6),track_in->getHelix_param(7),track_in->getHelix_param(8));
    track_out ->settrackid(track_in->gettrackid());
}

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

float calc_momentum_squared(float* mom)
{
     return  mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2] ;

}

TLorentzVector get_tlv_by_tracks(Ali_AS_Track* trackP, Ali_AS_Track* trackN, double massP, double massN, float* momP, float* momN )
{
    double energyP = sqrt(massP*massP+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
    double energyN = sqrt(massN*massN+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));

    TLorentzVector tlvP;
    tlvP.SetPxPyPzE(momP[0],momP[1],momP[2],energyP);
    TLorentzVector tlvN;
    tlvN.SetPxPyPzE(momN[0],momN[1],momN[2],energyN);

    TLorentzVector tlv = tlvP + tlvN;

    return tlv;

}


void fHelixAtoPointdca(TVector3 space_vec, Ali_AS_Track* helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {100.0,-100.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixA ->Evaluate(pA[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.01 && loopcounter < 100) // stops when the length is too small
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

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

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
            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-space_vec).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
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
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}

Int_t fCircle_Interception(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                                               Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
        Double_t dif_x = -(x1 - x2);
	Double_t dif_y = -(y1 - y2);
	Double_t dist=TMath::Sqrt( (dif_x*dif_x) + (dif_y*dif_y) );
        if(dist==0) // circles overlap -> most likely identical helices
        {
            x1_c = 0.0;
            y1_c = 0.0;
            x2_c = 0.0;
            y2_c = 0.0;
            return 0;
        };
        Double_t dist_inv =1/dist;
	
        if((dist<(r1+r2))&& (dist> TMath::Abs(r1-r2)))
        {
	//2 intersections
//		8<i<i<i<iii<i<iiii>ii


		Double_t x_loc 		=((dist*dist +r1*r1 -r2*r2)*(0.5*dist_inv));
		Double_t y_loc		=TMath::Sqrt(r1*r1 -x_loc*x_loc);
		
                x1_c	=(x1*dist+x_loc*dif_x -y_loc*dif_y)*dist_inv;
                y1_c	=(y1*dist+x_loc*dif_y +y_loc*dif_x)*dist_inv;
		
                x2_c	=(x1*dist+x_loc*dif_x +y_loc*dif_y)*dist_inv;
                y2_c	=(y1*dist+x_loc*dif_y -y_loc*dif_x)*dist_inv;

                return 1;
		}
		else if (dist<= TMath::Abs(r1-r2))
		{
			Double_t normrad	=(r1-dist+r2)*dist_inv*0.5;
            
            //Double_t normrad	= (r1 + 0.5*(dist - r1 - r2))*dist_inv;
			if(r2>r1)
			{	
                x1_c		=x1 - normrad*dif_x;
                y1_c		=y1 - normrad*dif_y;
			}
			else
			{
				x1_c		=x2 + normrad*dif_x;
                y1_c		=y2 + normrad*dif_y;
			
			}	
                x2_c		=x1_c;
                y2_c		=y1_c;
        
            	return 2;
		}	
        else
        {
	//no intersection (maybe 1)
            Double_t normrad	=(r1+dist-r2)*dist_inv*0.5;
            
            //Double_t normrad	= (r1 + 0.5*(dist - r1 - r2))*dist_inv;
           	x1_c		=x1 + normrad*dif_x;
       	    y1_c		=y1 + normrad*dif_y;
		
  	        x2_c		=x1_c;
            y2_c		=y1_c;
        
         	return 2;
		}	
	
	return 3;
        
}    

pair<Double_t,Double_t> fpathLength(Double_t r,Ali_AS_Track* helixA)
{
	//taken from https://www.star.bnl.gov/webdata/dox/html/StHelix_8cc_source.html
	
	pair<Double_t,Double_t> value;
	pair<Double_t,Double_t> VALUE(999999999.,999999999.);
	Double_t curvature		= fabs(helixA->getHelix_param(4));
	Double_t radius			= fabs(1/curvature);
	Double_t x0				= helixA->getHelix_param(5) +radius*TMath::Sin(helixA->getHelix_param(2));
	Double_t y0				= helixA->getHelix_param(0) -radius*TMath::Cos(helixA->getHelix_param(2));
	Double_t z0				= helixA->getHelix_param(1);

	Double_t phase			= helixA->getHelix_param(2) -TMath::Pi()/2;
	//Double_t phase=atan2f(-(x0-helixA->getHelix_param(5)),(y0-helixA->getHelix_param(0)));
	
	while(phase<-TMath::Pi()) phase	+= 2.*TMath::Pi() ;
	while(phase>TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	
	Double_t dipangle		= TMath::ATan(helixA->getHelix_param(3));
	Double_t cosdipangle	= TMath::Cos(dipangle);
	Double_t sindipangle	= TMath::Sin(dipangle);
	Double_t h				= TMath::Sign(1.,(Float_t)helixA->getHelix_param(4));
	if(phase> TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	Double_t cosphase		= TMath::Cos(phase);
	Double_t sinphase		= TMath::Sin(phase);

        //printf("h: %4.3f, phase: %4.3f, dipangle: %4.3f \n",h,phase,dipangle);

	Double_t t1 	= y0*curvature;
	Double_t t2 	= sinphase;
  	Double_t t3 	= curvature*curvature;
  	Double_t t4 	= y0*t2;
	Double_t t5 	= cosphase;
	Double_t t6 	= x0*t5;
	Double_t t8 	= x0*x0;
	Double_t t11	= y0*y0;
	Double_t t14 	= r*r;
	Double_t t15 	= t14*curvature;
	Double_t t17 	= t8*t8;
	Double_t t19 	= t11*t11;
	Double_t t21 	= t11*t3;
	Double_t t23 	= t5*t5;
	Double_t t32 	= t14*t14;
	Double_t t35 	= t14*t3;
	Double_t t38 	= 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*curvature*t6 +
				4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
				4.0*t8*x0*curvature*t5 - 4.0*t11*t23 -
				4.0*t11*y0*curvature*t2 + 4.0*t11 - 4.0*t14 +
				t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
	//cout<<t38<<endl;
	Double_t t40 	= (-t3*t38);
	if (t40<0.) return VALUE;
        t40 = ::sqrt(t40);

        if (h<0) phase += TMath::Pi();
	
	if(phase>TMath::Pi()) phase	-= 2.*TMath::Pi();

	Double_t t43 	= x0*curvature;
	Double_t t45 	= 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
	Double_t t46 	= h*cosdipangle*curvature;

	value.first 	= (-phase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
	value.second 	= -(phase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

	//
	//   Solution can be off by +/- one period, select smallest
	//
	Double_t p 		= fabs(2*TMath::Pi()/(h*curvature*cosdipangle));
	if (! std::isnan(value.first)) {
		if (fabs(value.first-p) < fabs(value.first)) value.first = value.first-p;
	   	else if (fabs(value.first+p) < fabs(value.first)) value.first = value.first+p;
	}
	if (! std::isnan(value.second)) {
	   	if (fabs(value.second-p) < fabs(value.second)) value.second = value.second-p;
	   	else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
	}
	
	if (value.first > value.second)
	swap(value.first,value.second);
	return(value);

	
}

Int_t fDCA_Helix_Estimate(Ali_AS_Track* helixA, Ali_AS_Track* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{

    // Calculates the 2D crossing point, calculates the corresponding 3D point and returns pathA and pathB
    //cout<<"a"<<endl;
    Double_t helix_point[3];

    Double_t x1 = helixA->getHelix_param(5);
    Double_t y1 = helixA->getHelix_param(0);
    Double_t x2 = helixB->getHelix_param(5);
    Double_t y2 = helixB->getHelix_param(0);
    Double_t c1 = helixA->getHelix_param(4);
    Double_t c2 = helixB->getHelix_param(4);
    Double_t r1 = 0.0;
    Double_t r2 = 0.0;
    if(c1 != 0 && c2 != 0)
    {
        r1 = fabs(1.0/c1);
        r2 = fabs(1.0/c2);
    } else return 0;
    //cout<<"B"<<endl;
    Double_t x1_c = 0.0;
    Double_t y1_c = 0.0;
    Double_t x2_c = 0.0;
    Double_t y2_c = 0.0;

    //Int_t bCross_points = fCross_points_Circles(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points (Alex): {%4.3f, %4.3f}, {%4.3f, %4.3f} \n",x1_c,y1_c,x2_c,y2_c);

    pathA = 0.0;
    pathB = 0.0;
    dcaAB = 0.0;

    //cout<<"C"<<endl;
    Int_t cCross_points = fCircle_Interception(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points (Sven): {%4.3f, %4.3f}, {%4.3f, %4.3f}, return: %d \n",x1_c,y1_c,x2_c,y2_c,cCross_points);
    //cout<<"C1"<<endl;
    Double_t radiusA = sqrt(x1_c*x1_c+y1_c*y1_c);
    Double_t radiusB = sqrt(x2_c*x2_c+y2_c*y2_c);
    //printf("radiusA: %4.3f, radiusB: %4.3f \n",radiusA,radiusB);
    //cout<<"C2"<<endl;
    ////cout << "bCross_points = " << bCross_points << ", xyr(1) = {" << x1 << ", " << y1 << ", " << r1
    //    << "}, xyr(2) = {"  << x2 << ", " << y2 << ", " << r2 << "}, p1 = {" << x1_c << ", " << y1_c << "}, p2 = {" << x2_c << ", " << y2_c << "}" << endl;

    //if(bCross_points == 0) return 0;
    //cout<<"D"<<endl;
    TVector3 pointA,pointB,pointA1,pointB1,pointA2,pointB2;

    Double_t path_lengthA_c1,path_lengthA_c2,path_lengthB_c1,path_lengthB_c2;

    // first crossing point for helix A
    pair< double, double > path_lengthA = fpathLength(radiusA,helixA);
    Double_t path_lengthA1 = path_lengthA.first;
    Double_t path_lengthA2 = path_lengthA.second;

    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    //cout<<"E"<<endl;
    if( ((x1_c-pointA1.x())*(x1_c-pointA1.x()) + (y1_c-pointA1.y())*(y1_c-pointA1.y())) <
       ((x1_c-pointA2.x())*(x1_c-pointA2.x()) + (y1_c-pointA2.y())*(y1_c-pointA2.y())))
    {
        path_lengthA_c1 = path_lengthA1;
    }
    else
    {
        path_lengthA_c1 = path_lengthA2;
    }
    //cout<<"f"<<endl;
    // second crossing point for helix A
    path_lengthA = fpathLength(radiusB,helixA);
    path_lengthA1 = path_lengthA.first;
    path_lengthA2 = path_lengthA.second;
    //cout<<"G"<<endl;
    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x2_c-pointA1.x())*(x2_c-pointA1.x()) + (y2_c-pointA1.y())*(y2_c-pointA1.y())) <
       ((x2_c-pointA2.x())*(x2_c-pointA2.x()) + (y2_c-pointA2.y())*(y2_c-pointA2.y())))
    {
        path_lengthA_c2 = path_lengthA1;
    }
    else
    {
        path_lengthA_c2 = path_lengthA2;
    }
    //cout<<"H"<<endl;
    // first crossing point for helix B
    pair< double, double > path_lengthB = fpathLength(radiusA,helixB);
    Double_t path_lengthB1 = path_lengthB.first;
    Double_t path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    //cout<<"I"<<endl;
    if( ((x1_c-pointB1.x())*(x1_c-pointB1.x()) + (y1_c-pointB1.y())*(y1_c-pointB1.y())) <
       ((x1_c-pointB2.x())*(x1_c-pointB2.x()) + (y1_c-pointB2.y())*(y1_c-pointB2.y())))
    {
        path_lengthB_c1 = path_lengthB1;
    }
    else
    {
        path_lengthB_c1 = path_lengthB2;
    }
    //cout<<"J"<<endl;
    // second crossing point for helix B
    path_lengthB = fpathLength(radiusB,helixB);
    path_lengthB1 = path_lengthB.first;
    path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    //cout<<"K"<<endl;
    if( ((x2_c-pointB1.x())*(x2_c-pointB1.x()) + (y2_c-pointB1.y())*(y2_c-pointB1.y())) <
       ((x2_c-pointB2.x())*(x2_c-pointB2.x()) + (y2_c-pointB2.y())*(y2_c-pointB2.y())))
    {
        path_lengthB_c2 = path_lengthB1;
    }
    else
    {
        path_lengthB_c2 = path_lengthB2;
    }
    //cout<<"L"<<endl;
    helixA ->Evaluate(path_lengthA_c1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA_c2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    helixB ->Evaluate(path_lengthB_c1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB_c2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    //cout<<"m"<<endl;

    Double_t dcaAB1 = (pointA1-pointB1).Mag();
    Double_t dcaAB2 = (pointA2-pointB2).Mag();
    Double_t dcaAB3 = (pointA1-pointB2).Mag();
    Double_t dcaAB4 = (pointA2-pointB1).Mag();

#if 0
    printf("pointA1: {%4.3f, %4.3f, %4.3f} \n",pointA1.X(),pointA1.Y(),pointA1.Z());
    printf("pointA2: {%4.3f, %4.3f, %4.3f} \n",pointA2.X(),pointA2.Y(),pointA2.Z());
    printf("pointB1: {%4.3f, %4.3f, %4.3f} \n",pointB1.X(),pointB1.Y(),pointB1.Z());
    printf("pointB2: {%4.3f, %4.3f, %4.3f} \n",pointB2.X(),pointB2.Y(),pointB2.Z());
    printf("dcaAB1: %4.3f, dcaAB2: %4.3f, dcaAB3: %4.3f, dcaAB4: %4.3f \n",dcaAB1,dcaAB2,dcaAB3,dcaAB4);
#endif
    //cout<<"N"<<endl;
    if(dcaAB1 < dcaAB2)
    {
        pathA = path_lengthA_c1;
        pathB = path_lengthB_c1;
        dcaAB = dcaAB1;
    }
    else
    {
        pathA = path_lengthA_c2;
        pathB = path_lengthB_c2;
        dcaAB = dcaAB2;
    }
    //cout<<"O"<<endl;

    return 1;

}



void fHelixABdca(Ali_AS_Track* helixA, Ali_AS_Track* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB,
                                    Float_t pathA_in, Float_t pathB_in)
{
    //cout << "Standard fHelixABdca called..." << endl;
    Float_t pA[2];
    Float_t pB[2]; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    if(pathA_in < -9990.0 && pathB_in < -9990.0)
    {
        pA[0] = 0.0;
        pA[1] = 0.0;
        pB[0] = 0.0;
        pB[1] = -70.0;
    }
    else
    {
        pA[0] = pathA_in+5.0;
        pA[1] = pathA_in-5.0;
        pB[0] = pathB_in+5.0;
        pB[1] = pathB_in-5.0;
    }
    Float_t distarray[2];
    TVector3 testA, testB;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixB ->Evaluate(pB[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation

        helixA ->Evaluate(pathA_dca,helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-testB).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            helixB ->Evaluate(pB[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[0] = (testA-testB).Mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            helixB ->Evaluate(pB[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-testB).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}


TVector3 find_sec_vertex(Ali_AS_Track* tracka, Ali_AS_Track* trackb )
{
    Float_t pathA_est, pathB_est, dcaAB_est;
    //printf("i_track_A: %d, i_track_B: %d, i_comb: %d \n",i_track_A,i_track_B,i_comb);
    //printf("3D cross point: {%4.3f, %4.3f, %4.3f} \n",vertex_point[0],vertex_point[1],vertex_point[2]);
    Int_t est_return = fDCA_Helix_Estimate(tracka,trackb,pathA_est,pathB_est,dcaAB_est);
    Float_t pathA, pathB, dcaAB;
    fHelixABdca(tracka,trackb,pathA,pathB,dcaAB,pathA_est,pathB_est);
    double vec_a[3];
    double vec_b[3];
    tracka->Evaluate(pathA,vec_a);
    trackb->Evaluate(pathB,vec_b);
    TVector3 sec_vertex;
    sec_vertex.SetXYZ( (vec_a[0]+vec_b[0])*0.5, (vec_a[1]+vec_b[1])*0.5 , (vec_a[2]+vec_b[2])*0.5 );

    return sec_vertex;

}

double get_dca_between_track(Ali_AS_Track* tracka, Ali_AS_Track* trackb )
{
    Float_t pathA_est, pathB_est, dcaAB_est;
    //printf("i_track_A: %d, i_track_B: %d, i_comb: %d \n",i_track_A,i_track_B,i_comb);
    //printf("3D cross point: {%4.3f, %4.3f, %4.3f} \n",vertex_point[0],vertex_point[1],vertex_point[2]);
    Int_t est_return = fDCA_Helix_Estimate(tracka,trackb,pathA_est,pathB_est,dcaAB_est);
    //cout<<"2"<<endl;
    Float_t pathA, pathB, dcaAB;
    fHelixABdca(tracka,trackb,pathA,pathB,dcaAB,pathA_est,pathB_est);

    return dcaAB;

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


    histo_delta.resize(3);
    delta_dca_vs_delta.resize(3);

    histo_delta[0] = new TH1D("deltax","deltax",120,-20,20);
    fListOfHistos->Add(histo_delta[0]);

    histo_delta[1] = new TH1D("deltay","deltay",120,-20,20);
    fListOfHistos->Add(histo_delta[1]);

    histo_delta[2] = new TH1D("deltaz","deltaz",120,-20,20);
    fListOfHistos->Add(histo_delta[2]);

    vec_histo_counter.resize(2);
    vec_histo_counter[0] = new TH1D("histo_counter","histo_counter",80,0,80);
    vec_histo_counter[1] = new TH1D("histo_pT_pions","histo_pT_pions",100,0,3);
    vec_histo_counter[2] = new TH1D("histo_pT_protons","histo_pT_protons",100,0,3);
    fListOfHistos->Add(vec_histo_counter[0]);
    fListOfHistos->Add(vec_histo_counter[1]);
    fListOfHistos->Add(vec_histo_counter[2]);


    TString arrn[3]={"x","y","z"};

    for(int i=0;i<3;i++)
    {
        TString n = "delta_dca_vs_delta_";
        n+=arrn[i];
        delta_dca_vs_delta[i]= new TH2D(n.Data(),n.Data(),120,-20,20,120,-20,20);
        fListOfHistos->Add(delta_dca_vs_delta[i]);


    }

    vec_t_prof.resize(2);
    vec_t_prof[0]= new TProfile("tprof","tprof",20,0,20);
    vec_t_prof[1]= new TProfile("tprof_range_0_2","",20,-0.05,1.95);
    vec_t_prof[1]->GetXaxis()->SetTitle("dca lower limit [cm]");
    vec_t_prof[1]->GetYaxis()->SetTitle("efficiency of TOF-detector");

    fListOfHistos->Add(vec_t_prof[0]);
    fListOfHistos->Add(vec_t_prof[1]);


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
    //Tree_AS_Event->SetAutoSave( -4900000 );

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

    //if(counter_events>25){return;}

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

    if(N_tracks == 0)
    {
	// Skip empty event
	return;
    }


    Int_t ncascades =  fESD-> GetNumberOfCascades();
    Int_t numberV0  =  fESD ->GetNumberOfV0s () ;


    //cout<<"numberV0: "<<numberV0<<endl;
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
    AS_Event ->clearNUCLEVList();
    AS_Event ->clearDMparticleList();
    AS_Event ->clearV0List();
    AS_Event ->setTriggerWord(fESD->GetFiredTriggerClasses());
    AS_Event ->setx(PrimVertex->GetX());
    AS_Event ->sety(PrimVertex->GetY());
    AS_Event ->setz(PrimVertex->GetZ());
    AS_Event ->setid(RunNum);
    AS_Event ->setN_tracks(N_tracks);
    AS_Event ->setBeamIntAA(MeanBeamIntAA);
    AS_Event ->setT0zVertex(T0zVertex);
    //AS_Event ->setN_V0s(numberV0);
   

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

     double xprim, yprim, zprim;
     xprim = PrimVertex->GetX();
     yprim = PrimVertex->GetY();
     zprim = PrimVertex->GetZ();

     TVector3 pos_primary_vertex;
     pos_primary_vertex.SetXYZ(xprim,yprim,zprim);

     vector<int> all_used_positive_track_ids_for_V0s;
     vector<int> all_used_negative_track_ids_for_V0s;

     TLorentzVector* tlv_pos = new TLorentzVector();
     TLorentzVector* tlv_neg = new TLorentzVector();
     TLorentzVector* tlv_Lambda = new TLorentzVector();
     TLorentzVector* tlv_Kaon = new TLorentzVector();
     TLorentzVector* tlv_gamma = new TLorentzVector();
     TLorentzVector  tlv_anti_p_and_K_plus;

     const Float_t mass_proton = 0.93827208816 ;  //in GeV?
     const Float_t mass_pion = 0.139657061 ;  //in GeV?
     const Float_t mass_electron = 0.510998950 * 1e-3 ;  //in GeV?
     const double  mass_K = 0.493677 ;  //in GeV?
     const double mass_K0 = 0.497614;
     const double mass_Lambda = 1.1155683;
     const double mass_neutron = 0.939565;
     double S_mass = -1;

     TVector3 position_SV2;
     TVector3 position_SV3;
     TVector3 direction_SV2;
     TVector3 direction_SV3;

     vector<TVector3> vec_position_SV2;
     vector<TVector3> vec_position_SV3;
     vector<TVector3> vec_direction_SV2;
     vector<TVector3> vec_direction_SV3;

     vector<int> vec_SV2_number;
     vector<int> vec_SV3_number;

     vector<Ali_AS_Track> vec_SV2_tracks;
     vector<Ali_AS_Track> vec_SV3_tracks;

     vector<int> vec_SV2_track_ids;
     vector<int> vec_SV3_track_ids;

     
     Ali_AS_V0* AS_V0 = new Ali_AS_V0;
     Ali_AS_Track* as_trackP = new Ali_AS_Track;
     Ali_AS_Track* as_trackP_save = new Ali_AS_Track;
     Ali_AS_Track* as_trackN = new Ali_AS_Track;
     Ali_AS_Track* as_trackN_save = new Ali_AS_Track;
     Ali_AS_Track* tracka  =  new Ali_AS_Track;
     Ali_AS_Track* trackb  =  new Ali_AS_Track;

     Float_t* pos = new Float_t[3];
     Float_t* momP = new Float_t[3];
     Float_t* momN = new Float_t[3];

     double dcaP=-10000;
     double dcaN=-10000;
     Int_t trackidP,trackidN;
     Float_t energy_proton,energy_pion,energy_antiproton,energy_pion_plus,energy_pion_minus, energy_K_plus,energy_anti_proton;
     Float_t energy_electron_plus,energy_electron_minus;
     Double_t invariantmass = -1.;
     double dcaV0;
     Double_t x=0;
     Double_t y=0;
     Double_t z=0;

     double sigma_proton_TPC;
     double sigma_antiproton_TPC;
     double sigma_pion_plus_TPC ;
     double sigma_pion_minus_TPC;

     double sigma_K_plus_TPC;

     Float_t path_closest_to_point = 0;
     Float_t path_closest_to_point2 = 0;
     Float_t dca_closest_to_point  = 0;
     Float_t path_initA = 0.0;
     Float_t path_initB = 30.0;

     float radiuscuts[4]{1,3,5,7};
     float dcaprimcuts[4]{0.5,1,2,5};
     Double_t momentumP;
     Double_t momentumN;
     double radius;
     TVector3 vec_primtoV0;

     vector<Ali_AS_Track> vec_antip_and_K_plus_V0s;
     vector<Ali_AS_Track> vec_antip_and_pi_plus;
     vector<TVector3> vec_S_vertex_positions;
     vector<TLorentzVector> vec_tlv_antip_and_K;
     vector<TLorentzVector> vec_tlv_type4;
     vector<TVector3> vec_ch4_S_vertex_positions;
     vector<Ali_AS_Track> vec_K0_V0s;
     vector<TLorentzVector> vec_K0_tlvs;
     vector<TVector3> vec_K0_V0_pos;

     vector<int> brute_force_type1;
     vector<int> brute_force_type2;
     vector<int> brute_force_type3;
     vector<int> brute_force_type4;


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
                //cout<<"nitscls: "<<N_ITS_cls<<endl;
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

        float path_closest_to_point,dca;
        FindDCAHelixPoint2(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca);
        printf("dca: %f, momentum: %f, phi: %f, eta: %f, pathlength: %f, TOF: %f \n", dca ,Track_p,Track_phi,Track_eta,Track_length,TOF_signal);

        Double_t TOF_raw = track -> GetTOFsignalRaw();
        Double_t TOF_chi2 = track -> GetTOFchi2();
        double TOFToT =  track->GetTOFsignalToT();
        int tofnclus = track-> GetNTOFclusters ();
        double tofsigdx = track->    GetTOFsignalDx();

        cout<<"TOF raw: "<<TOF_raw<<" TOF_chi2: "<<TOF_chi2<<" TOFTot: "<<TOFToT<<"tofnclus:" <<tofnclus
            <<" tofsigdx: "<<tofsigdx<<endl;
        cout<<""<<endl;
        //count tracks
        vec_histo_counter[0]->Fill(6.5);

        //tracks with TOF hit
        if(TOF_signal<99990){vec_histo_counter[0]->Fill(7.5);}

        if(fabs( Track_pT ) >0.3 && fabs(Track_eta) < 0.8 )
        {
            vec_histo_counter[0]->Fill(8.5);
            if(TOF_signal<99990){vec_histo_counter[0]->Fill(9.5);}
        }

        if(dca>5)
        {
            vec_histo_counter[0]->Fill(10.5);
            if(TOF_signal<99990){vec_histo_counter[0]->Fill(11.5);}

        }

        int step=0;
        for(double i=0;i<2;i+=0.3)
        {
            if(dca>i && fabs( Track_pT ) >0.4 && fabs(Track_eta) < 0.8)
            {
                vec_histo_counter[0]->Fill(20.5+step*2);
                if(TOF_signal<99990){vec_histo_counter[0]->Fill(20.5+step*2+1);}

            }

            if(dca>i && fabs( Track_pT ) >0.4 && fabs(Track_eta) < 0.8 && fabs(Track_PID[3])<2.5)
            {
                vec_histo_counter[0]->Fill(40.5+step*2);
                if(TOF_signal<99990){vec_histo_counter[0]->Fill(40.5+step*2+1);}
            }

            
            step++;
        }

        if(fabs( Track_pT ) >0.4 && fabs(Track_eta) < 0.8)
        {
            for(double i=0;i<20 ; i++)
            {
                if(dca>i)
                {
                    if(TOF_signal<99990) {vec_t_prof[0]->Fill(i+0.5,1);  }
                    else {vec_t_prof[0]->Fill(i+0.5,0); }
                }
            }

            for(double i=0;i<2;i+=0.1)
            {
                if(dca>i)
                {
                    if(TOF_signal<99990) {vec_t_prof[1]->Fill(i,1);  }
                    else {vec_t_prof[1]->Fill(i,0); }
                }
            }
        }

        double m_squared = calculate_m_squared_by_TOF(AS_Track);


        //pion
        if(fabs(Track_PID[2])<2.5 )
        {
            vec_histo_counter[1]->Fill(Track_pT);
        }

        //protons

        if(fabs(Track_PID[4])<2.5 && m_squared>0.6 && m_squared < 1.2)
        {
            vec_histo_counter[2]->Fill(Track_pT);
        }




    } // End of TPC track loop

    //-----------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event

    //-------------------------------------------------------------------------------------------------
    //find similiar tracks
    vector<int> trackids_bits_shared;
    vector<int> trackids_similiar_tracks;
    vector<int> trackids_similiar_used;


    for(Int_t atrack = 0; atrack < NumTracks; atrack++)
    {
        tracka =  AS_Event -> getTrack(atrack);
        TLorentzVector tlva = tracka->get_TLV_part();
        int trackida = tracka->gettrackid();
        //if ( check_if_int_is_in_vector(trackida,trackids_similiar_used)){continue;}

        for(int btrack=atrack+1; btrack<NumTracks; btrack++)
        {
            if(btrack==atrack){continue;}

            trackb = AS_Event -> getTrack(btrack);
            int trackidb = trackb->gettrackid();
            //if ( check_if_int_is_in_vector(trackidb,trackids_similiar_used)){continue;}

            TLorentzVector tlvb = trackb->get_TLV_part();

            TBits tbitsshared = tracka->getbitsshared();
            int numbitsshared = tbitsshared.CountBits();


            if( fabs( 1-tlva[0]/tlvb[0]) < 0.05 && fabs( 1-tlva[1]/tlvb[1]) < 0.05 && fabs( 1-tlva[2]/tlvb[2]) < 0.05)
            {
                
                printf("tracka: %d, trackb %d \n",trackida,trackidb);
                cout<<"numbitsshared: "<<numbitsshared<<endl;
                cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

                double  pos[3];
                float dca=-1;

                tracka->Evaluate(0,pos);
                TVector3 position;
                position[0]=pos[0];
                position[1]=pos[1];
                position[2]=pos[2];
                FindDCAHelixPoint2(position,trackb,path_initA,path_initB,path_closest_to_point,dca);
                cout<<"dca: "<<dca<<endl;
                if(dca<0.5)
                {
                    cout<<"pushed back"<<endl;
                    trackids_similiar_tracks.push_back(trackidb);

                    trackids_similiar_used.push_back(trackida);
                    trackids_similiar_used.push_back(trackidb);
                }


            }


        }

    }   //end of finding similiar tracks



    for (Int_t V0_counter=0; V0_counter<numberV0; V0_counter++)
    {
        //get position of V0-----------------------
        AliESDv0 *V0=fESD->GetV0(V0_counter);
        x=0;
        y=0;
        z=0;
        V0->AliESDv0::GetXYZ(x,y,z);

        TVector3 V0_position;
        V0_position.SetXYZ(x,y,z);


        //cout<<""<<endl;
       // printf("V0 number: %d \n",V0_counter);
       // printf("x: %f,y: %f, z: %f \n",x,y,z);
        
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


        momentumP = sqrt(pxP*pxP + pyP*pyP + pzP*pzP);
        momentumN = sqrt(pxN*pxN + pyN*pyN + pzN*pzN);

        //printf("trackidP %d, trackidN %d \n",indexP,indexN);
        //printf("momentum of positive particle: %f   momentum of negative particle: %f \n",momentumP,momentumN);
        //----------------------------------------------------------------

        //get impulse vector by using AliESDtrack class-------------------------
       
        trackN->GetInnerPxPyPz(pN);
        trackP->GetInnerPxPyPz(pP);

        //cout<<pxN<<pN[0]<<endl;
        //cout<<pyN<<pN[1]<<endl;
        //-------------------------------------------------

        radius = sqrt ( (x-xprim)*(x-xprim) + (y-yprim)*(y-yprim) + (z-zprim)*(z-zprim)  );
        //radius cut???


        //create element of Ali_AS_V0 class

        //AS_V0  = AS_Event ->createV0();

        //AS_V0->~Ali_AS_V0();

        //use set functions
        AS_V0 -> setxyz(x,y,z);
        AS_V0 -> setNpxpypz(pxN,pyN,pzN);
        AS_V0 -> setPpxpypz(pxP,pyP,pzP);

        AS_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

        TVector3 vec_primtoV0;
        vec_primtoV0.SetXYZ((x-xprim),(y-yprim),(z-zprim));
        TVector3 unit_prim_to_V0;
        unit_prim_to_V0 = vec_primtoV0.Unit();


        //create tracks for positive and negative particle
        //Ali_AS_Track* as_trackP = AS_V0->createTrack();
        //Ali_AS_Track* as_trackN = AS_V0->createTrack();


        //fill tracks------------------------------------------------------------------
        //as_trackP  ->clearTRD_digit_list();
        //as_trackP  ->clearOfflineTrackletList();

        Int_t N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(trackP ->HasPointOnITSLayer(i_ITS_layer))
	    {
                N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
                //cout<<"nitscls: "<<N_ITS_cls<<endl;
	    }
	}

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
        as_trackP  ->setNITScls(N_ITS_cls);   
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

        N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(trackN ->HasPointOnITSLayer(i_ITS_layer))
	    {
                N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
                //cout<<"nitscls: "<<N_ITS_cls<<endl;
	    }
	}

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
        as_trackN -> setNTPCcls(N_ITS_cls);
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

        //AS_V0->clearTrackList();
        //end of filling for N particle track-------------------------------------------------------------------
        //FindDCAHelixPoint2()

        //cout<<"a"<<endl;
        /*
        TVector3 sec_vertex_found = find_sec_vertex(as_trackP,as_trackN);

        double delta_x = sec_vertex_found[0]-V0_position[0];
        double delta_y = sec_vertex_found[1]-V0_position[1];
        double delta_z = sec_vertex_found[2]-V0_position[2];

        histo_delta[0]->Fill(delta_x);
        histo_delta[1]->Fill(delta_y);
        histo_delta[2]->Fill(delta_z);

        double dca_wir = get_dca_between_track(as_trackP,as_trackN);
        double dca_ESD = AS_V0 -> getdcaV0();

        delta_dca_vs_delta[0]->Fill(delta_x,dca_wir-dca_ESD);
        delta_dca_vs_delta[1]->Fill(delta_y,dca_wir-dca_ESD);
        delta_dca_vs_delta[2]->Fill(delta_z,dca_wir-dca_ESD);
        */
        //cout<<"b"<<endl;

        //start of analysis-----------------------------------------------------------
        //variables:
        dcaP=-10000;
        dcaN=-10000;
        
        //double momentumP,momentumN;
        
        
        //------------------------------------------------------------------------------
        pos = AS_V0 -> getxyz();
        //cout<<"posx: "<<pos[0]<<endl;
        //position.SetXYZ(pos[0],pos[1],pos[2]);
        radius = sqrt( (pos[0]-xprim) *(pos[0]-xprim)+(pos[1]-yprim)*(pos[1]-yprim)+(pos[2]-zprim)*(pos[2]-zprim) );

        
        vec_primtoV0.SetXYZ((pos[0]-xprim),(pos[1]-yprim),(pos[2]-zprim));
        //printf("x %f,y %f, z %f \n",pos[0],pos[1],pos[2]);


        //get momentum of posititve particle
        //momP = AS_V0 -> getPpxpypz();
        momP[0]=pxP;
        momP[1]=pyP;
        momP[2]=pzP;
        //get momentum for negative particle

        //momN = AS_V0 -> getNpxpypz();
        momN[0]=pxN;
        momN[1]=pyN;
        momN[2]=pzN;


        dcaV0 = AS_V0 -> getdcaV0();
        dcaP = as_trackP->getdca();
        dcaN = as_trackN->getdca();

        float dcaV01,dcaV02;
        float dcaAprim,dcaBprim;

        FindDCAHelixPoint2(pos_primary_vertex,as_trackP,path_initA,path_initB,path_closest_to_point,dcaAprim);
        FindDCAHelixPoint2(pos_primary_vertex,as_trackN,path_initA,path_initB,path_closest_to_point,dcaBprim);

        
        //FindDCAHelixPoint2(pos,as_trackP,path_initA,path_initB,path_closest_to_point,dcaV01);
        //FindDCAHelixPoint2(pos,as_trackN,path_initA,path_initB,path_closest_to_point,dcaV02);

        //printf("dcaV0(esd): %f, dcaV01: %f, dcaV02: %f \n",dcaV0,dcaV01,dcaV02);
        //printf("dcaprim1(esd): %f, dcaprim2: %f, dcaprim1(calc): %f dcaprim2: %f \n",dcaP,dcaN,dcaAprim,dcaBprim);

        if(dcaP < 0){continue;}
        if(dcaN > 0){continue;}

        trackidP = as_trackP->gettrackid();
        trackidN = as_trackN->gettrackid();

        if ( check_if_int_is_in_vector(trackidP,trackids_similiar_tracks)){continue;}
        if ( check_if_int_is_in_vector(trackidN,trackids_similiar_tracks)){continue;}

        

        //momentumP = sqrt( momP[0]*momP[0] + momP[1]*momP[1]+ momP[2]*momP[2] );
        //momentumN = sqrt( momN[0]*momN[0] + momN[1]*momN[1]+ momN[2]*momN[2] );
        //all_used_positive_track_ids_for_V0s.push_back(trackidP);
        //all_used_negative_track_ids_for_V0s.push_back(trackidN);

        sigma_proton_TPC = as_trackP -> getnsigma_p_TPC();
        sigma_antiproton_TPC = as_trackN -> getnsigma_p_TPC();
        sigma_pion_plus_TPC = as_trackP -> getnsigma_pi_TPC();
        sigma_pion_minus_TPC = as_trackN -> getnsigma_pi_TPC();

        sigma_K_plus_TPC = as_trackP -> getnsigma_K_TPC();

        double sigma_e_P  = as_trackP->getnsigma_e_TPC();
        double sigma_e_N  = as_trackN->getnsigma_e_TPC();
        double sigma_pi_P = as_trackP-> getnsigma_pi_TPC();
        double sigma_pi_N = as_trackN-> getnsigma_pi_TPC();
        double sigma_K_P  = as_trackP-> getnsigma_K_TPC();
        double sigma_K_N  = as_trackN-> getnsigma_K_TPC();
        double sigma_p_P  = as_trackP-> getnsigma_p_TPC();
        double sigma_p_N  = as_trackN-> getnsigma_p_TPC();

        path_closest_to_point = 0;
        path_closest_to_point2 = 0;
        dca_closest_to_point  = 0;
        path_initA = 0.0;
        path_initB = 30.0;

        //****************************************************************************************+
        //****************************************************************************************+
        //-----------------------------------------------------------------------------------
        //S-particle search
        //----------------------------------------------------------------------------

        //----------------------------------------------------------------------------
        //channel1: antiS + n -> antiLambda + K0 + pi- + pi+


        //antilambda
        if(fabs(sigma_antiproton_TPC) < 2.0 && fabs(sigma_pion_plus_TPC) < 2.0)
        {
            energy_antiproton = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            energy_pion       = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1] ,momN[2],energy_antiproton);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();

            //cut on mass
            if(invariantmass < (1.1157+0.001495*6) && invariantmass > (1.1157-0.001495*6))
            {
                //create V0
                /*
                Ali_AS_V0* save_V0 = AS_Event->createV0();


                //Setters
                save_V0 -> setxyz(x,y,z);
                save_V0 -> setNpxpypz(pxN,pyN,pzN);
                save_V0 -> setPpxpypz(pxP,pyP,pzP);
                save_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

                //createtracks
                as_trackP_save = save_V0->createTrack();
                as_trackN_save = save_V0->createTrack();

                copy_track_params(as_trackP , as_trackP_save);
                copy_track_params(as_trackN , as_trackN_save);
                */

                position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                vec_SV3_tracks.push_back(*as_trackP);
                vec_SV3_tracks.push_back(*as_trackN);

                vec_SV3_track_ids.push_back(trackidP);
                vec_SV3_track_ids.push_back(trackidN);

                vec_SV3_number.push_back(V0_counter);
                //cout<<"SV3"<<endl;
                //cout<<"invariantmass: "<<invariantmass<<endl;
                SV3_counter++;

            }

        }


        //K0 - > pi+ and pi-
        //check if pion+ and pion-
        if(fabs(sigma_pion_plus_TPC) < 2.0 && fabs(sigma_pion_minus_TPC) < 2.0)
        {
            energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();

            //cut on mass
            if(invariantmass< (0.4981+0.0042*6) && invariantmass > (0.4981-0.0042*6))
            {
                /*
                Ali_AS_V0* save_V0 = AS_Event->createV0();

                //Setters
                save_V0 -> setxyz(x,y,z);
                save_V0 -> setNpxpypz(pxN,pyN,pzN);
                save_V0 -> setPpxpypz(pxP,pyP,pzP);
                save_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

                //createtracks
                as_trackP_save = save_V0->createTrack();
                as_trackN_save = save_V0->createTrack();

                copy_track_params(as_trackP , as_trackP_save);
                copy_track_params(as_trackN , as_trackN_save);
                */

                position_SV2.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV2.push_back(position_SV2);

                direction_SV2.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                vec_direction_SV2.push_back(direction_SV2);

                vec_SV2_tracks.push_back(*as_trackP);
                vec_SV2_tracks.push_back(*as_trackN);

                vec_SV2_track_ids.push_back(trackidP);
                vec_SV2_track_ids.push_back(trackidN);

                vec_SV2_number.push_back(V0_counter);
                //cout<<"SV2"<<endl;


            }
        }


        //**********************************************************************************+
        //**********************************************************************************+

        //channel2: antiS + p -> antiproton + K+ + K+ + (pi0)

        //search for V0 coming from anti-proton and K+
        if (fabs(sigma_antiproton_TPC) < 2.0 && fabs(sigma_K_plus_TPC < 2.0))
        {
            //  cout<<"V0 from anti-proton and K+"<<endl;
            //energy_K_plus      = sqrt(mass_K * mass_K + (momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_K_plus      = sqrt(mass_K * mass_K + calc_momentum_squared(momP) );
            energy_anti_proton = sqrt(mass_proton * mass_proton + calc_momentum_squared(momN)) ;

            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_K_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_anti_proton);

            tlv_anti_p_and_K_plus = *tlv_pos + *tlv_neg;

            invariantmass = tlv_anti_p_and_K_plus.M();


            //counters[5]++;

            //assume position of V0 is also position of potential S-vertex
            //loop over all other tracks in order to find K+ that comes from this position
            vector<int> save_track_ids;
            vector<int> save_track_nums;

            for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
            {

                AS_Track = AS_Event->getTrack(i_track_A);
                Int_t Trackid  = AS_Track->gettrackid();
                if ( check_if_int_is_in_vector(Trackid,trackids_similiar_tracks)){continue;}
                //cout<<Trackid<<endl;

                if(Trackid==trackidP){continue;}
                if(Trackid==trackidN){continue;}
                if ( check_if_int_is_in_vector(Trackid,brute_force_type2)){continue;}

                //if( check_if_int_is_in_vector(Trackid,used_track_ids_of_pions) == 1 ){continue;}

                double sigma = AS_Track -> getnsigma_K_TPC();

                // PID for Kaon
                if(fabs(sigma) > 2.0) {continue;}
                //initial parameters
                Float_t path_closest_to_point = 0;
                Float_t dca_closest_to_point  = 0;

                //calculate dca from vertex to particle track
                FindDCAHelixPoint2(pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                if(dca_closest_to_point > 1.) {continue;}
                //if(radius<5) {continue;}

                //check dca of additional K+ to primary vertex
                Float_t dca_to_primary_vertex = -1;
                path_closest_to_point=0;

                FindDCAHelixPoint2(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca_to_primary_vertex);
                if(dca_to_primary_vertex<3.){continue;}

                //cout<<"Found V0 of antiproton and K+ with another K+"<<endl;

                save_track_ids.push_back(Trackid);
                save_track_nums.push_back(i_track_A);



            }

            if(save_track_ids.size()==1)
            {
                //counter_vertices_antip_K_plus_K_plus++;
                //counters[2]++;
                //histos_1D[9]->Fill(radius);
                //histos_2D[3]->Fill(pos[0],pos[1]);

                AS_Track = AS_Event->getTrack(save_track_nums[0]);
                Int_t Trackid  = AS_Track->gettrackid();
                TLorentzVector tlv = AS_Track->get_TLV_part();

                TLorentzVector tlv_S = tlv_anti_p_and_K_plus + tlv;


                TVector3 dir;
                dir.SetXYZ(tlv_anti_p_and_K_plus[0]+tlv[0] , tlv_anti_p_and_K_plus[1]+tlv[1], tlv_anti_p_and_K_plus[2]+tlv[2]);

                TVector3 unit_dir;
                unit_dir = dir.Unit();

                double dot_product = unit_dir.Dot(unit_prim_to_V0);

                //TVector3 null;
                //null.SetXYZ(0.,0.,0.);


                if(radius>5)
                {
                    //counter_vertices_antip_K_plus_K_plus_r_larger_5++;

                    if(dot_product<0.6){continue;}

                    //counter_vertices_antip_K_plus_K_plus_r_larger_5_and_dot_product++;

                    //histo_counter->Fill(3.5);

                    double m_squared1 = calculate_m_squared_by_TOF(AS_Track);
                    double m_squared2 = calculate_m_squared_by_TOF(as_trackP);

                    //if(m_squared1!=-1.) {mass_squared_kaons_and_background->Fill(m_squared1);}
                    //if(m_squared2!=-1.) {mass_squared_kaons_and_background->Fill(m_squared2);}

                    //printf("mass squared Kaon 1: %f, Kaon 2: %f \n",m_squared1,m_squared2);
                    //printf("trackids: %d %d %d \n",trackidP,trackidN,Trackid);

                    //if(m_squared1>0.2 && m_squared1<0.35 && m_squared2>0.2 && m_squared2<0.35)
                   // {
                        //histo_counter->Fill(7.5);
                        ////cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
                        ////cout<<"hi"<<endl;
                        ////cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

                    if ( check_if_int_is_in_vector(trackidP,brute_force_type2)){continue;}
                    if ( check_if_int_is_in_vector(trackidN,brute_force_type2)){continue;}
                    if ( check_if_int_is_in_vector(Trackid,brute_force_type2)){continue;}


                    TVector3 S_vertex_pos;
                    //S_vertex_pos.SetXYZ(pos[0],pos[1],pos[2]);
                    S_vertex_pos = V0_position;

                    Ali_AS_DM_particle* DMparticle = AS_Event->createDMparticle();

                    DMparticle->settype(2);

                    DMparticle->set_primVertex(pos_primary_vertex);
                    DMparticle->set_S1Vertex(S_vertex_pos);
                    
                    DMparticle->set_DirSV1(dir);
                    DMparticle->setN_V0s(1);

                    DMparticle->set_tlv(tlv_S);

                    Ali_AS_Track* dmtrack0 = DMparticle->createTrack();
                    copy_track_params(as_trackP,dmtrack0);

                    Ali_AS_Track* dmtrack1 = DMparticle->createTrack();
                    copy_track_params(as_trackN,dmtrack1);

                    Ali_AS_Track* dmtrack2 = DMparticle->createTrack();
                    copy_track_params(AS_Track,dmtrack2);

                    printf("type2; position x: %f, y: %f, z: %f \n",S_vertex_pos[0],S_vertex_pos[1],S_vertex_pos[2]);


                    brute_force_type2.push_back(trackidP);
                    brute_force_type2.push_back(trackidN);
                    brute_force_type2.push_back(Trackid);
                    /*
                    AS_DM_particle ->set_primVertex(pos_primary_vertex);
                    AS_DM_particle ->set_S1Vertex(pos);
                    //AS_DM_particle ->set_S2Vertex(null);
                    //AS_DM_particle ->set_S3Vertex(null);
                    AS_DM_particle ->set_DirSV1(dir);
                    //AS_DM_particle ->set_DirSV2(null);
                    //AS_DM_particle ->set_DirSV3(null);
                    AS_DM_particle ->setN_V0s(1);
                    AS_DM_particle ->clearTrackList();

                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(as_trackP,AS_DM_Track);
                    ////cout<<ASTrack1->gettrackid()<<endl;

                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(as_trackN,AS_DM_Track);
                    ////cout<<ASTrack2->gettrackid()<<endl;

                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(AS_Track,AS_DM_Track);
                    ////cout<<vec_SV2_tracks[2*vector_loop_SV2]->gettrackid()<<endl;

                    Tree_AS_DM_particle ->Fill();
                    //cout<<"filled Tree"<<endl;
                    //cout<<""<<endl;

                    TBits tbitsshared = AS_Track->getbitsshared();
                    int numbitsshared = tbitsshared.CountBits();
                    //cout<<"numbitsshared: "<<numbitsshared<<endl;
                    //histo_counter->SetBinContent(7,numbitsshared);
                    */
               // }

                }

            }

            ////cout<<"invariantmass: "<<invariantmass<<endl;
            save_track_ids.clear();

        }

        //**********************************************************************************+
        //**********************************************************************************+


        //channel3:
        //AntiS + p  ->  antip + K+  + K0 + pi+

        //mass_squared mitnehmen zusaetzlich

        if( fabs(sigma_p_N)<2.0  && fabs(sigma_K_P)<2.0 )
        {
            vec_antip_and_K_plus_V0s.push_back(*as_trackP);
            vec_antip_and_K_plus_V0s.push_back(*as_trackN);

            double energy_kaon_plus = sqrt(mass_K*mass_K+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            double energy_proton_minus = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));

            tlv_pos-> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            TLorentzVector tlv_antip_and_K = *tlv_pos + *tlv_neg;
            vec_tlv_antip_and_K.push_back(tlv_antip_and_K);

            TVector3 S_vertex_position;
            S_vertex_position.SetXYZ(pos[0],pos[1],pos[2]);
            vec_S_vertex_positions.push_back(S_vertex_position);
        }


        //K0s
        if(fabs(sigma_pion_plus_TPC) < 2.0 && fabs(sigma_pion_minus_TPC) < 2.0)
        {
            energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();

            //cut on mass
            if(invariantmass< (0.4981+0.0042*6) && invariantmass > (0.4981-0.0042*6))
            {
                vec_K0_V0s.push_back(*as_trackP);
                vec_K0_V0s.push_back(*as_trackN);
                vec_K0_tlvs.push_back(*tlv_Kaon);
                TVector3 K0_V0_position;
                K0_V0_position.SetXYZ(pos[0],pos[1],pos[2]);
                vec_K0_V0_pos.push_back(K0_V0_position);
            }
        }


        //**********************************************************************************+
        //**********************************************************************************+
        //channel4

        //antiS  +  n  -> antip  +  2*K0  +  pi+

        if( fabs(sigma_p_N)<2.0  &&  fabs(sigma_pi_P)<2.0 )
        {
            vec_antip_and_pi_plus.push_back(*as_trackP);
            vec_antip_and_pi_plus.push_back(*as_trackN);
            TVector3 S_vertex_position;
            S_vertex_position.SetXYZ(pos[0],pos[1],pos[2]);
            vec_ch4_S_vertex_positions.push_back(S_vertex_position);

            TLorentzVector tlv_type4 = get_tlv_by_tracks(as_trackP,as_trackN, mass_pion, mass_proton, momP, momN );
            vec_tlv_type4.push_back(tlv_type4);

        }

        //use K0s from above!

      






        //----------------------------------------------------------------------------
        //NUCLEV search
        //---------------------------------------------------------------------------

        /*
        double sigma_e_P  = as_trackP->getnsigma_e_TPC();
        double sigma_e_N  = as_trackN->getnsigma_e_TPC();
        double sigma_pi_P = as_trackP-> getnsigma_pi_TPC();
        double sigma_pi_N = as_trackN-> getnsigma_pi_TPC();
        double sigma_K_P  = as_trackP-> getnsigma_K_TPC();
        double sigma_K_N  = as_trackN-> getnsigma_K_TPC();
        double sigma_p_P  = as_trackP-> getnsigma_p_TPC();
        double sigma_p_N  = as_trackN-> getnsigma_p_TPC();
        */
        Float_t TPCdEdx   = as_trackP->getTPCdEdx();
        Float_t tofsignal = as_trackP->getTOFsignal();
        Float_t dca       = as_trackP->getdca();
        Float_t tracklength = as_trackP->getTrack_length();
        Float_t energy_proton,energy_pion,energy_antiproton,energy_pion_plus,energy_pion_minus, energy_K_plus,energy_anti_proton;
        Float_t energy_electron_plus,energy_electron_minus;
        double invariantmass;

        TLorentzVector tlv_in_loop;
        double m_squared;

        if(dcaV0>3.){continue;}
        if(radius<2.){continue;}

        if(fabs(sigma_e_P)<2.5 )
        {
            energy_electron_plus  = sqrt(mass_electron*mass_electron+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_electron_plus);
            invariantmass = tlv_pos ->M();


            tlv_in_loop = as_trackP->get_TLV_part();
            double momentum = tlv_in_loop.P();

            //printf("dEdx: %f, momentum: %f, dca: %f, charge: %d, tofsignal: %f \n"
            //     ,TPCdEdx, momentum,dca, charge, tofsignal);


            if(tofsignal>99990){continue;}
            double velocity = tracklength/tofsignal;

            //printf("velocity: %f \n", velocity);

            velocity = velocity * 1e10;


            double speed_of_light_SI = 299792458.;

            velocity = velocity /speed_of_light_SI;  //now in units of c

            // printf("velocity: %f \n", velocity);

            double gamma_squared = 1. / (1-velocity*velocity) ;
            //printf("momentum: %f, gamma: %f, velocity: %f \n",momentum,gamma,velocity);

            //m^2 =  ( p/ (gamma * v) )^2
            m_squared = ( momentum / velocity)  *  ( momentum /velocity) * 1./gamma_squared ;

            //histo_invariantmass_electron_plus->Fill(m_squared);

            //if(m_squared<0){printf("mass squared: %f \n", m_squared);}
            //printf("mass squared: %f \n", m_squared);
        }


        //gamma  (two electrons)
        if(fabs(sigma_e_P)<2.5 && fabs(sigma_e_N)<2.5) 
        {
            
            energy_electron_plus  = sqrt(mass_electron*mass_electron+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_electron_minus = sqrt(mass_electron*mass_electron+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_electron_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_electron_minus);
            *tlv_gamma= *tlv_pos + *tlv_neg;
            invariantmass = tlv_gamma->M();

            //printf("invariante Masse gamma: %f \n",invariantmass);
            //histo_invariantmass_gamma->Fill(invariantmass);

            if( (m_squared<0.4 && tofsignal<99990) || tofsignal>99990)
            {
                if(invariantmass<0+0.05){continue;}
            }
        }

        //K0s (two pions)
        if(fabs(sigma_pi_P)<2.5 && fabs(sigma_pi_N)<2.5 )
        {
            energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            ////cout<<energy_proton<<endl;
            ////cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();

            if(invariantmass< (0.4981+0.0042*2) && invariantmass > (0.4981-0.0042*2))
            {
                continue;
            }
        }

        //Lambda -> proton + pi-
        if(fabs(sigma_p_P)<2.5 && fabs(sigma_pi_N)<2.5)
        {
            Float_t energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            Float_t energy_pion   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();
            ////cout<<invariantmass<<endl;
            if(invariantmass < (1.1157+0.001495*2) && invariantmass > (1.1157-0.001495*2))
            {
                continue;
            }

        }
         //cout<<"n"<<endl;
        //Anti-Lambda -> antiproton + p+
        if(fabs(sigma_p_N)<2.5 && fabs(sigma_pi_P)<2.5)
        {
            energy_antiproton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_plus   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();
            ////cout<<invariantmass<<endl;
            if(invariantmass < (1.1157+0.001495*2) && invariantmass > (1.1157-0.001495*2))
            {
                continue;
            }
        }

        vector<Ali_AS_Track*> vec_tracks;
        vector<double> vec_dca_to_V0;
        vector<double> vec_path;

        //loop over all tracks
        for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
        {

            AS_Track = AS_Event->getTrack(i_track_A);
            Int_t Trackid  = AS_Track->gettrackid();
            if ( check_if_int_is_in_vector(Trackid,trackids_similiar_tracks)){continue;}
            //cout<<Trackid<<endl;

            if(Trackid==trackidP){continue;}
            if(Trackid==trackidN){continue;}

            Float_t path_closest_to_point = 0;
            Float_t dca_closest_to_point  = 0;

            //calculate dca from vertex to particle track
            FindDCAHelixPoint2(pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

            if(dca_closest_to_point > 1.) {continue;}

            Float_t dca_to_prim = -1;
            //FindDCAHelixPoint2(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point2,dca_to_prim);
            dca_to_prim = fabs( AS_Track->getdca()  );
            if(dca_to_prim<0.5){continue;}

            vec_tracks.push_back(AS_Track);
            vec_dca_to_V0.push_back(dca_closest_to_point);
            vec_path.push_back(path_closest_to_point);
        }

        //if(vec_tracks.size()>=1)
        // {

            /*
            AS_NUCLEV = AS_Event->createNUCLEV();

            AS_NUCLEV ->set_secVertex(pos);
            AS_NUCLEV ->setdcaV0_ESD(dcaV0);
            Ali_AS_Track* track1 = AS_NUCLEV->createTrack();
            Ali_AS_Track* track2 = AS_NUCLEV->createTrack();
            copy_track_params(as_trackP , track1);
            copy_track_params(as_trackN , track2);

            FindDCAHelixPoint2(pos,track1,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            AS_NUCLEV->addpath(path_closest_to_point);
            AS_NUCLEV->adddca(dca_closest_to_point);

            FindDCAHelixPoint2(pos,track2,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            AS_NUCLEV->addpath(path_closest_to_point);
            AS_NUCLEV->adddca(dca_closest_to_point);

            for(int i=0;i<vec_tracks.size();i++)
            {
                AS_Track = AS_NUCLEV->createTrack();
                copy_track_params(vec_tracks[i] , AS_Track);
                AS_NUCLEV->addpath(vec_path[i]);
                AS_NUCLEV->adddca(vec_dca_to_V0[i]);

            }
            */

       // }




        /*
        //lambda
        if(fabs(sigma_proton_TPC) < 2.5 && fabs(sigma_pion_minus_TPC) < 2.5)
        {
            energy_pion  = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            energy_proton       = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();

            //cut on mass
            if(invariantmass < (1.1157+0.001495*8) && invariantmass > (1.1157-0.001495*8))
            {
                /*
                //create V0
                Ali_AS_V0* save_V0 = AS_Event->createV0();

                //Setters
                save_V0 -> setxyz(x,y,z);
                save_V0 -> setNpxpypz(pxN,pyN,pzN);
                save_V0 -> setPpxpypz(pxP,pyP,pzP);
                save_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

                //createtracks
                as_trackP_save = save_V0->createTrack();
                as_trackN_save = save_V0->createTrack();

                copy_track_params(as_trackP , as_trackP_save);
                copy_track_params(as_trackN , as_trackN_save);

            }
        }

        //antilambda
        if(fabs(sigma_antiproton_TPC) < 2.5 && fabs(sigma_pion_plus_TPC) < 2.5)
        {
            energy_antiproton = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            energy_pion       = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();

            //cut on mass
            if(invariantmass < (1.1157+0.001495*8) && invariantmass > (1.1157-0.001495*8))
            {
                //create V0
                /*
                Ali_AS_V0* save_V0 = AS_Event->createV0();


                //Setters
                save_V0 -> setxyz(x,y,z);
                save_V0 -> setNpxpypz(pxN,pyN,pzN);
                save_V0 -> setPpxpypz(pxP,pyP,pzP);
                save_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

                //createtracks
                as_trackP_save = save_V0->createTrack();
                as_trackN_save = save_V0->createTrack();

                copy_track_params(as_trackP , as_trackP_save);
                copy_track_params(as_trackN , as_trackN_save);

                /*
                position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                vec_SV3_tracks.push_back(as_trackP);
                vec_SV3_tracks.push_back(as_trackN);

                vec_SV3_track_ids.push_back(trackidP);
                vec_SV3_track_ids.push_back(trackidN);

                vec_SV3_number.push_back(V0_counter);
                //cout<<"SV3"<<endl;
                //cout<<"invariantmass: "<<invariantmass<<endl;
                SV3_counter++;

            }

        }

        //K0 - > pi+ and pi-
        //check if pion+ and pion-
        if(fabs(sigma_pion_plus_TPC) < 2.5 && fabs(sigma_pion_minus_TPC) < 2.5)
        {
            energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();

            //cut on mass
            if(invariantmass< (0.4981+0.0042*8) && invariantmass > (0.4981-0.0042*8))
            {
                /*
                Ali_AS_V0* save_V0 = AS_Event->createV0();

                //Setters
                save_V0 -> setxyz(x,y,z);
                save_V0 -> setNpxpypz(pxN,pyN,pzN);
                save_V0 -> setPpxpypz(pxP,pyP,pzP);
                save_V0 -> setdcaV0( V0->GetDcaV0Daughters() );

                //createtracks
                as_trackP_save = save_V0->createTrack();
                as_trackN_save = save_V0->createTrack();

                copy_track_params(as_trackP , as_trackP_save);
                copy_track_params(as_trackN , as_trackN_save);

                /*
                position_SV2.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV2.push_back(position_SV2);

                direction_SV2.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                vec_direction_SV2.push_back(direction_SV2);

                vec_SV2_tracks.push_back(as_trackP);
                vec_SV2_tracks.push_back(as_trackN);

                vec_SV2_track_ids.push_back(trackidP);
                vec_SV2_track_ids.push_back(trackidN);

                vec_SV2_number.push_back(V0_counter);
                //cout<<"SV2"<<endl;


            }
        }
        */


    }
    //cout<<"end of V0 loop"<<endl;
    //cout<<""<<endl;
    //end of V0 loop

    



    //----------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------

    //cout << "Event number: " << fEventNoInFile << ", event number with TRD digits: " << N_good_events << endl;
    //printf("Event number: %d, N_tracks: %d, cent(V0M): %f , cent(CL0): %f \n",fEventNoInFile,N_tracks,MultSelection->GetMultiplicityPercentile("SPDTracklets"));
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------


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

    /*
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
    */
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------

    TVector3 S_vertex_pos;
    TLorentzVector* tlv_SV1 = new TLorentzVector();  //tlv_SV2+tlv_SV3-tlv_neutron
    TLorentzVector* tlv_neutron = new TLorentzVector();  //neutron
    TLorentzVector* tlv_proton = new TLorentzVector();  //proton
    tlv_neutron -> SetPxPyPzE(0.,0.,0.,mass_neutron);
    tlv_proton ->SetPxPyPzE(0.,0.,0.,mass_proton);

    float radiusS;
    TLorentzVector* tlv_SV2 = new TLorentzVector();  //kaon
    TLorentzVector* tlv_SV3 = new TLorentzVector();  //lambda
    TVector3 vec_primary_vertex_to_SV1;
    TVector3 unit_vec_primary_vertex_to_SV1;

    TVector3 momentum_SV1;
    TVector3 unit_momentum_SV1;
    
    path_closest_to_point = 0;
    dca_closest_to_point  = 0;
    path_initA = 0.0;
    path_initB = 30.0;

    Ali_AS_Track* ASTrack1 = new Ali_AS_Track;
    Ali_AS_Track* ASTrack2 = new Ali_AS_Track;

    
    vector<int> createdV0s;
    vector<int> tracknumbers;
    vector<int> trackids;
    vector<int> alltrackids;
    vector<int> V0numbers;
    vector<int> trackidsSV2andSV3;

    int counter_V0s = 0;

    //------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------
    //find S-vertex:
    //channel1


    for(Int_t vector_loop_SV3 = 0; vector_loop_SV3 < vec_position_SV3.size(); vector_loop_SV3++)
    {
        for(Int_t vector_loop_SV2 = 0; vector_loop_SV2 < vec_position_SV2.size(); vector_loop_SV2++)
        {
            if(vec_position_SV3.size() > 0 && vec_position_SV2.size() > 0 && vec_direction_SV3.size() > 0 &&  vec_direction_SV2.size() > 0)
            {
                if(calculateMinimumDistance(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]) > 0.5){ continue;}
                S_vertex_pos = calcVertexAnalytical(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]);

                double px,py,pz,E;
                //SV2
                px = vec_direction_SV2[vector_loop_SV2][0];
                py = vec_direction_SV2[vector_loop_SV2][1];
                pz = vec_direction_SV2[vector_loop_SV2][2];
                E = (sqrt(px*px+py*py+pz*pz)+(mass_K0*mass_K0));
                tlv_SV2 -> SetPxPyPzE(px,py,pz,E);

                //SV3
                px = vec_direction_SV3[vector_loop_SV3][0];
                py = vec_direction_SV3[vector_loop_SV3][1];
                pz = vec_direction_SV3[vector_loop_SV3][2];
                E = (sqrt(px*px+py*py+pz*pz)+(mass_Lambda*mass_Lambda));
                tlv_SV3 -> SetPxPyPzE(px,py,pz,E);

               
                //calculate SV1
                *tlv_SV1 = *tlv_SV2 + *tlv_SV3 - *tlv_neutron;

                //neglecting pions -> wrong mass
                //S_mass = tlv_SV1->M();

                vec_primary_vertex_to_SV1.SetXYZ(S_vertex_pos[0]-xprim ,S_vertex_pos[1]-yprim , S_vertex_pos[2]-zprim);
                unit_vec_primary_vertex_to_SV1 = vec_primary_vertex_to_SV1.Unit();

                radiusS = vec_primary_vertex_to_SV1.Mag();

                if(radiusS<5){continue;}

                //decay length cut
                TVector3 vec_S_to_Lambda;
                TVector3 vec_S_to_K0;
                TVector3 pos_Lambda;
                TVector3 pos_K0;

                pos_Lambda = vec_position_SV3[vector_loop_SV3];
                pos_K0 = vec_position_SV2[vector_loop_SV2];

                vec_S_to_Lambda = pos_Lambda - S_vertex_pos;
                vec_S_to_K0 = pos_K0 - S_vertex_pos;

                double radius_Lambda = vec_S_to_Lambda.Mag();
                double radius_K0     = vec_S_to_K0.Mag();

                if(radius_Lambda < 1){continue;}
                if(radius_K0 < 0.5) {continue;}


                //check direction daughter particles not backwards:
                TVector3 unit_vec_S_to_Lambda = vec_S_to_Lambda.Unit();
                TVector3 unit_vec_S_to_K0 = vec_S_to_K0.Unit();
                TVector3 unit_tlv_Lambda = vec_direction_SV3[vector_loop_SV3].Unit();
                TVector3 unit_tlv_K0 = vec_direction_SV2[vector_loop_SV2].Unit();

                if(unit_vec_S_to_Lambda.Dot(unit_tlv_Lambda)< 0.6 ){continue;}
                if(unit_vec_S_to_K0.Dot(unit_tlv_K0) < 0.6 ){continue;}

                Int_t counter_pions_close_to_S_vertex = 0;

                //store tracks of all pions that come from S-vertex
                
                tracknumbers.clear();
                trackids.clear();
                trackidsSV2andSV3.clear();

                trackidsSV2andSV3.push_back(vec_SV3_track_ids[2*vector_loop_SV3]);
                trackidsSV2andSV3.push_back(vec_SV3_track_ids[2*vector_loop_SV3+1]);
                trackidsSV2andSV3.push_back(vec_SV2_track_ids[2*vector_loop_SV2]);
                trackidsSV2andSV3.push_back(vec_SV2_track_ids[2*vector_loop_SV2+1]);


                //for each S-vertex loop over all tracks of event
                for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
                {
                    AS_Track = AS_Event->getTrack(i_track_A);
                    int trackid2 = AS_Track->gettrackid();

                    if ( check_if_int_is_in_vector(trackid2,trackids_similiar_tracks)){continue;}
                    if( check_if_int_is_in_vector(trackid2,trackidsSV2andSV3) ){continue;}

                    double sigma = AS_Track -> getnsigma_pi_TPC();
                    // Do some PID here for pi+ and pi-
                    if(fabs(sigma)>2.0){continue;}

                    float dcaprim=0;
                    FindDCAHelixPoint2(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dcaprim);
                    if(dcaprim<1.){continue;}


                    //calculate dca from vertex to particle track
                    FindDCAHelixPoint2(S_vertex_pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                    //cut on dca and distance S vertex from primary vertex
                    //cut on r because otherwise many pions from primary vertex
                    if(dca_closest_to_point < 0.5) // && radiusS>5)
                    {
                        //printf("track: %d is pion and close to vertex \n",i_track_A);
                        counter_pions_close_to_S_vertex++;
                        //save track number for tracks that are close to S-vertex
                        tracknumbers.push_back(i_track_A);
                        trackids.push_back(trackid2);
                    }
                }
                if(trackids.size()==1)
                {
                    // counters[4]++;
                }

                if(trackids.size()==2)
                {

                    //check all tracks if one of them is already used:
                    int a = check_if_int_is_in_vector(trackids[0],brute_force_type1);
                    int b = check_if_int_is_in_vector(trackids[1],brute_force_type1);
                    int c = check_if_int_is_in_vector(vec_SV2_track_ids[2*vector_loop_SV2],brute_force_type1);
                    int d = check_if_int_is_in_vector(vec_SV2_track_ids[2*vector_loop_SV2+1],brute_force_type1);
                    int e = check_if_int_is_in_vector(vec_SV3_track_ids[2*vector_loop_SV3],brute_force_type1);
                    int f = check_if_int_is_in_vector(vec_SV3_track_ids[2*vector_loop_SV3+1],brute_force_type1);

                    int check = a+b+c+d+e+f;

                    if(check>=1){continue;}

                    Double_t r1[3];
                    Double_t r2[3];

                    ASTrack1 = AS_Event->getTrack(tracknumbers[0]);
                    ASTrack2 = AS_Event->getTrack(tracknumbers[1]);

                    double dca1 = ASTrack1->getdca();
                    double dca2 = ASTrack2->getdca();

                    //check if one is negative and one is positive
                    if (dca1 * dca2 >0) {continue;}

                    //do all for pion 1--------------------------------------------------------------
                    //---------------------------------------------------------------------------------
                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;

                    //calculate again path and dca
                    FindDCAHelixPoint2(S_vertex_pos,ASTrack1,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                    //find direction of pion 1
                    ASTrack1->Evaluate(path_closest_to_point,r1);
                    ASTrack1->Evaluate(path_closest_to_point+0.01,r2);

                    TVector3 mom_dir_pion1;
                    TVector3 unit_mom_dir_pion1;
                    TVector3 vec_momentum;

                    mom_dir_pion1.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);

                    //unit vector in direction of pion 1
                    unit_mom_dir_pion1 = mom_dir_pion1.Unit();

                    //cout<<r1[0]<<endl;
                    //cout<<r2[0]<<endl;

                    TLorentzVector tlv = ASTrack1->get_TLV_part();
                    double momentum = tlv.P();

                    vec_momentum.SetXYZ(unit_mom_dir_pion1[0]*momentum,unit_mom_dir_pion1[1]*momentum,unit_mom_dir_pion1[2]*momentum);

                    //cout<<"momentum vec: "<<vec_momentum[0]<<endl;

                    TLorentzVector tlv_pion1;
                    double energy_pion1 = sqrt(mass_pion * mass_pion + momentum*momentum) ;

                    tlv_pion1.SetPxPyPzE(vec_momentum[0],vec_momentum[1],vec_momentum[2],energy_pion1);

                    *tlv_SV1+=tlv_pion1;

                    //do the same for pion 2
                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;

                    //calculate again path and dca
                    FindDCAHelixPoint2(S_vertex_pos,ASTrack2,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                    //find direction of pion 2
                    ASTrack2->Evaluate(path_closest_to_point,r1);
                    ASTrack2->Evaluate(path_closest_to_point+0.01,r2);

                    TVector3 mom_dir_pion2;
                    TVector3 unit_mom_dir_pion2;
                    TVector3 vec_momentum_pion2;

                    mom_dir_pion2.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);

                    //unit vector in direction of pion 1
                    unit_mom_dir_pion2 = mom_dir_pion2.Unit();

                    TLorentzVector tlv2 = ASTrack2->get_TLV_part();
                    double momentum2 = tlv2.P();

                    vec_momentum_pion2.SetXYZ(unit_mom_dir_pion2[0]*momentum2,unit_mom_dir_pion2[1]*momentum2,unit_mom_dir_pion2[2]*momentum2);

                    TLorentzVector tlv_pion2;
                    double energy_pion2 = sqrt(mass_pion * mass_pion + momentum2*momentum2) ;
                    tlv_pion2.SetPxPyPzE(vec_momentum_pion2[0],vec_momentum_pion2[1],vec_momentum_pion2[2],energy_pion2);

                    *tlv_SV1+=tlv_pion2;

                    //build unit vector of momentum from tlv_SV1
                    momentum_SV1.SetXYZ(tlv_SV1->Px(),tlv_SV1->Py(),tlv_SV1->Pz());
                    unit_momentum_SV1 = momentum_SV1.Unit();

                    double dot_product = unit_vec_primary_vertex_to_SV1.Dot(unit_momentum_SV1);
                    // printf("dot product: %f \n", dot_product);

                    //check if dot product is larger than 0
                    if(dot_product<0.8){continue;}



                    //--------------------------------------------------------------------------------------------------------------

                    //calculate again S-mass----------------------------------------------
                    double S_mass_correct = tlv_SV1->M();
                    //cout<<"falsche S-mass: "<<S_mass<<endl;
                    //cout<<"korrekte S-mass: "<<S_mass_correct<<endl;



                    alltrackids.clear();
                    alltrackids.push_back(trackids[0]);   //pion1
                    alltrackids.push_back(trackids[1]);   //pion2
                    alltrackids.push_back(vec_SV2_track_ids[2*vector_loop_SV2]);   //tracks of SV2
                    alltrackids.push_back(vec_SV2_track_ids[2*vector_loop_SV2+1]);   //tracks of SV2
                    alltrackids.push_back(vec_SV3_track_ids[2*vector_loop_SV3]);   //tracks of SV3
                    alltrackids.push_back(vec_SV3_track_ids[2*vector_loop_SV3+1]);   //tracks of SV3

                    if ( check_if_value_is_doppelt_in_vector(alltrackids) ) {continue;}


                    cout<<"all trackids: "<<endl;
                    printf("%d %d %d %d %d %d \n",alltrackids[0],alltrackids[1],alltrackids[2],alltrackids[3],alltrackids[4],alltrackids[5]);
                    //vertices that fulfill all cuts

                    //after all cuts store all 6 tracks in brute force
                    brute_force_type1.push_back(trackids[0]);   //pion1
                    brute_force_type1.push_back(trackids[1]);   //pion2
                    brute_force_type1.push_back(vec_SV2_track_ids[2*vector_loop_SV2]);   //tracks of SV2
                    brute_force_type1.push_back(vec_SV2_track_ids[2*vector_loop_SV2+1]);   //tracks of SV2
                    brute_force_type1.push_back(vec_SV3_track_ids[2*vector_loop_SV3]);   //tracks of SV3
                    brute_force_type1.push_back(vec_SV3_track_ids[2*vector_loop_SV3+1]);   //tracks of SV3

                    //get position of V0-----------------------

                    V0numbers.clear();
                    V0numbers.push_back(vec_SV2_number[vector_loop_SV2]);
                    V0numbers.push_back(vec_SV3_number[vector_loop_SV3]);

                    Ali_AS_DM_particle* DMparticle = AS_Event->createDMparticle();

                    DMparticle->settype(1);

                    DMparticle->set_primVertex(pos_primary_vertex);
                    DMparticle->set_S1Vertex(S_vertex_pos);

                    printf("type1; position x: %f, y: %f, z: %f \n",S_vertex_pos[0],S_vertex_pos[1],S_vertex_pos[2]);

                    DMparticle->set_S2Vertex(vec_position_SV2[vector_loop_SV2]);
                    DMparticle->set_S3Vertex(vec_position_SV3[vector_loop_SV3]);
                    DMparticle->set_DirSV1(momentum_SV1);
                    DMparticle->set_DirSV2(vec_direction_SV2[vector_loop_SV2]);
                    DMparticle->set_DirSV3(vec_direction_SV3[vector_loop_SV3]);
                    DMparticle->setN_V0s(2);
                    DMparticle->set_tlv(*tlv_SV1);

                    Ali_AS_Track* dmtrack0 = DMparticle->createTrack();
                    copy_track_params(ASTrack1,dmtrack0);       //first pion

                    Ali_AS_Track* dmtrack1 = DMparticle->createTrack();
                    copy_track_params(ASTrack2,dmtrack1);   //second pion
 
                    Ali_AS_Track* dmtrack2 = DMparticle->createTrack();
                    Ali_AS_Track* ptr_track = &vec_SV2_tracks[2*vector_loop_SV2];
                    copy_track_params(ptr_track,dmtrack2);
 
                    Ali_AS_Track* dmtrack3 = DMparticle->createTrack();
                    ptr_track = &vec_SV2_tracks[2*vector_loop_SV2+1];
                    copy_track_params(ptr_track,dmtrack3);
 
                    Ali_AS_Track* dmtrack4 = DMparticle->createTrack();
                    ptr_track = &vec_SV3_tracks[2*vector_loop_SV3];
                    copy_track_params(ptr_track,dmtrack4);
 
                    Ali_AS_Track* dmtrack5 = DMparticle->createTrack();
                    ptr_track = &vec_SV3_tracks[2*vector_loop_SV3+1];
                    copy_track_params(ptr_track,dmtrack5);

                    /*
                    //cout<<"trackids: "<<endl;
                    cout<<ASTrack1->gettrackid()<<endl;
                    cout<<ASTrack2->gettrackid()<<endl;
                    cout<<vec_SV2_tracks[2*vector_loop_SV2].gettrackid()<<endl;
                    cout<<vec_SV2_tracks[2*vector_loop_SV2+1].gettrackid()<<endl;
                    cout<<vec_SV3_tracks[2*vector_loop_SV3].gettrackid()<<endl;
                    cout<<vec_SV3_tracks[2*vector_loop_SV3+1].gettrackid()<<endl;
                    */
                    /*
                    for(int i=0;i<6;i++)
                    {
                        Ali_AS_Track* track = DMparticle->getTrack(i);
                        int trackid = track->gettrackid();
                       // cout<<"trackid: "<<trackid<<endl;
 
                    }


                    //cout<<"all trackids: "<<endl;
                    printf("%d %d %d %d %d %d \n",trackids[0],trackids[1],vec_SV2_track_ids[2*vector_loop_SV2],vec_SV2_track_ids[2*vector_loop_SV2+1],
                           vec_SV3_track_ids[2*vector_loop_SV3],vec_SV3_track_ids[2*vector_loop_SV3+1]);
                    */

                       /*
                   for(int i=0;i<2;i++)
                   {
                       Double_t x=0;
                       Double_t y=0;
                       Double_t z=0;
                       AliESDv0 *V0=fESD->GetV0(V0numbers[i]);
                       V0->AliESDv0::GetXYZ(x,y,z);

                       //get impulse of particle N and P by AliESDv0 class
                       Double_t pxN,pyN,pzN;
                       V0->GetNPxPyPz(pxN,pyN,pzN);
                       Double_t pxP,pyP,pzP;
                       V0->GetPPxPyPz(pxP,pyP,pzP);
                       //-------------------------------------
                       //get track (by using index)-----------------------------
                       Int_t indexN = V0->GetNindex();
                       Int_t indexP = V0->GetPindex();

                       AliESDtrack* trackN = fESD->AliESDEvent::GetTrack(indexN);
                       AliESDtrack* trackP = fESD->AliESDEvent::GetTrack(indexP);

                       Double_t momentumP = sqrt(pxP*pxP + pyP*pyP + pzP*pzP);
                       Double_t momentumN = sqrt(pxN*pxN + pyN*pyN + pzN*pzN);

                       trackN->GetInnerPxPyPz(pN);
                       trackP->GetInnerPxPyPz(pP);
                       double radius = sqrt ( (x-xprim)*(x-xprim) + (y-yprim)*(y-yprim) + (z-zprim)*(z-zprim)  );
                       if(radius<5){continue;}
                       if(check_if_int_is_in_vector(V0numbers[i],createdV0s)){continue;}

                       cout<<"created V0" <<endl;
                       cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
                       cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;


                       AS_V0  = AS_Event ->createV0();
                       counter_V0s++;
                       //cout<<"NumV0s: "<<AS_Event->getN_V0s()<<endl;
                       
                       createdV0s.push_back(V0numbers[i]);
                       //use set functions
                       AS_V0 -> setxyz(x,y,z);
                       //printf("V0 number: %d \n",V0_counter);
                       printf("x: %f,y: %f, z: %f \n",x,y,z);
                       //cout<<"x: "<<x<<endl;
                       
                       momP[0]=pxP;
                       momP[1]=pyP;
                       momP[2]=pzP;
                       //get momentum for negative particle

                       //momN = AS_V0 -> getNpxpypz();
                       momN[0]=pxN;
                       momN[1]=pyN;
                       momN[2]=pzN;

                       printf("pxP: %f,pyP: %f, pzP: %f \n",momP[0],momP[1],momP[2]);
                       printf("pxN: %f,pyN: %f, pzN: %f \n",momN[0],momN[1],momN[2]);

                       AS_V0 -> setNpxpypz(pxN,pyN,pzN);
                       AS_V0 -> setPpxpypz(pxP,pyP,pzP);
                       
                       energy_antiproton = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                       energy_pion       = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));

                       cout<<"unter Wurzel: "<<mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
                       cout<<"masspion: "<<mass_pion<<endl;
                       cout<<"energy_antiproton: "<<energy_antiproton<<endl;
                       cout<<"energy_pion: "<<energy_pion<<endl;

                       tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
                       tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

                       *tlv_Lambda = *tlv_pos + *tlv_neg;
                       invariantmass = tlv_Lambda->M();

                       cout<<"invariantmass: "<<invariantmass<<endl;

                       AS_V0 -> setdcaV0( V0->GetDcaV0Daughters() );
                       //create tracks for positive and negative particle
                       as_trackP = AS_V0->createTrack();
                       as_trackN = AS_V0->createTrack();

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

                   }
                   */
                //-----------------------------------------------------------------------------------------------
                //-----------------------------------------------------------------------------------------------
                   //-----------------------
                   //------------------

                }


            }

        }
    }


    //********************************************************************************
    //********************************************************************************
    //channel3  find S-vertex position

    int sizeV1 =  vec_S_vertex_positions.size();
    int sizeV2 =  vec_K0_V0_pos.size();
    

    for(int antipKplus_loop=0; antipKplus_loop<sizeV1; antipKplus_loop++)
    {
        for(int K0_loop=0; K0_loop<sizeV2; K0_loop++)
        {
            vector<Ali_AS_Track> all_tracks_channel3;
            all_tracks_channel3.resize(5);

            vector<int> all_tracks_ids_channel3;
            all_tracks_ids_channel3.resize(5);

            TVector3 S_vertex_pos = vec_S_vertex_positions[antipKplus_loop];
            TVector3 K0_pos = vec_K0_V0_pos[K0_loop];

            TLorentzVector tlv_K0_V0 = vec_K0_tlvs[K0_loop];
            TVector3 dir_K0_V0;
            dir_K0_V0.SetXYZ(tlv_K0_V0.X(),tlv_K0_V0.Y(),tlv_K0_V0.Z());
            TVector3 pos_K0_V0 = vec_K0_V0_pos[K0_loop];

            TVector3 vec_prim_to_S_vertex;
            TVector3 unit_vec_prim_to_S_vertex;
            vec_prim_to_S_vertex.SetXYZ(S_vertex_pos[0]-xprim,S_vertex_pos[1]-yprim,S_vertex_pos[2]-zprim);
            double radiusS = vec_prim_to_S_vertex.Mag();
            unit_vec_prim_to_S_vertex = vec_prim_to_S_vertex.Unit();

            if(radiusS<5.){continue;}

            TVector3 vec_S_to_K0 = K0_pos-S_vertex_pos;
            TVector3 unit_vec_S_to_K0 = vec_S_to_K0.Unit();
            TVector3 unit_dir_K0_V0 = dir_K0_V0.Unit();

            if ( unit_vec_S_to_K0.Dot( unit_dir_K0_V0 ) < 0.8  ) {continue;}


            all_tracks_channel3[0] = vec_antip_and_K_plus_V0s[2*antipKplus_loop];
            all_tracks_channel3[1] = vec_antip_and_K_plus_V0s[2*antipKplus_loop+1];
            all_tracks_channel3[2] = vec_K0_V0s[2*K0_loop];
            all_tracks_channel3[3] = vec_K0_V0s[2*K0_loop+1];

            for(int i=0;i<4;i++)
            {
                int trackid = all_tracks_channel3[i].gettrackid();
                all_tracks_ids_channel3[i]= trackid;
            }


            double distance  =  calculateMinimumDistanceStraightToPoint(pos_K0_V0,dir_K0_V0,S_vertex_pos);

            if(distance>0.5){continue;}

            vector<Ali_AS_Track> vec_pion_track;
            for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
            {
                AS_Track = AS_Event->getTrack(i_track_A);
                int trackid = AS_Track->gettrackid();
                if ( check_if_int_is_in_vector(trackid,trackids_similiar_tracks)){continue;}
                if (check_if_int_is_in_vector(trackid , all_tracks_ids_channel3) ) {continue;}
                if (check_if_int_is_in_vector(trackid , brute_force_type3) ) {continue;}

                double sigma_pi = AS_Track ->getnsigma_pi_TPC();
                if(sigma_pi>2.0){continue;}

                path_closest_to_point = 0;
                dca_closest_to_point  = 0;
                FindDCAHelixPoint2(S_vertex_pos, AS_Track, path_initA, path_initB, path_closest_to_point,dca_closest_to_point);
                if(dca_closest_to_point>0.5){continue;}

                float dcaprim;
                FindDCAHelixPoint2(pos_primary_vertex, AS_Track, path_initA, path_initB, path_closest_to_point,dcaprim);
                if(dcaprim<1.){continue;}

                vec_pion_track.push_back(*AS_Track);

            }

            if(vec_pion_track.size()==1)
            {
                all_tracks_channel3[4]=vec_pion_track[0];
                TLorentzVector tlv_pion = vec_pion_track[0].get_TLV_part();

                if ( check_if_value_is_doppelt_in_vector(all_tracks_ids_channel3) ){continue;}
                if (check_if_two_vectors_have_same_element(all_tracks_ids_channel3,brute_force_type3) ){continue;}

                TLorentzVector tlv_type3 = vec_tlv_antip_and_K[antipKplus_loop];
                tlv_type3 += tlv_K0_V0 ;
                tlv_type3 += tlv_pion;
                tlv_type3 -= *tlv_proton;

                TVector3 dir_SV1;
                TVector3 unit_dir_SV1;
                dir_SV1.SetXYZ(tlv_type3[0],tlv_type3[1],tlv_type3[2]);
                unit_dir_SV1 = dir_SV1.Unit();

                if( unit_dir_SV1.Dot(unit_vec_prim_to_S_vertex) < 0.8 ){continue;}

                Ali_AS_DM_particle* DMparticle = AS_Event->createDMparticle();

                DMparticle->settype(3);

                DMparticle->set_primVertex(pos_primary_vertex);
                DMparticle->set_S1Vertex(S_vertex_pos);

                //printf("type3; position x: %f, y: %f, z: %f \n",S_vertex_pos[0],S_vertex_pos[1],S_vertex_pos[2]);

                DMparticle->set_S2Vertex(pos_K0_V0);

                DMparticle -> set_tlv(tlv_type3);
                //DMparticle->set_DirSV1(momentum_SV1);
                //DMparticle->set_DirSV2(vec_direction_SV2[vector_loop_SV2]);
                //DMparticle->set_DirSV3(vec_direction_SV3[vector_loop_SV3]);
                //DMparticle->setN_V0s(2);

                for(int i=0;i<5;i++)
                {
                    Ali_AS_Track* dmtrack = DMparticle->createTrack();
                    Ali_AS_Track* ptr_all_tracks = &all_tracks_channel3[i];
                    copy_track_params(ptr_all_tracks,dmtrack);

                    int trackid = all_tracks_channel3[i].gettrackid();
                    brute_force_type3.push_back(trackid);
                    //cout<<"trackid: "<<trackid<<endl;

                }

            }


        }
    }


    //********************************************************************************
    //********************************************************************************
    //channel4  do combinatoric

    //from vector of K0s choose 2
    //build lines from each vector
    //for each line calculate dca to S_vertex pos

    int size3 =  vec_ch4_S_vertex_positions.size();

    for(int K0_loop1=0; K0_loop1<sizeV2; K0_loop1++)
    {
        for(int K0_loop2=K0_loop1+1; K0_loop2<sizeV2; K0_loop2++)
        {
            for(int S_vertex_loop=0;S_vertex_loop<size3;S_vertex_loop++)
            {
                TVector3 S_vertex_pos =  vec_ch4_S_vertex_positions[S_vertex_loop];
                TVector3 pos_K0_V01 = vec_K0_V0_pos[K0_loop1];
                TVector3 pos_K0_V02 = vec_K0_V0_pos[K0_loop2];

                TVector3 vec_prim_to_S_vertex;
                TVector3 unit_vec_prim_to_S_vertex;
                vec_prim_to_S_vertex.SetXYZ(S_vertex_pos[0]-xprim,S_vertex_pos[1]-yprim,S_vertex_pos[2]-zprim);
                double radiusS = vec_prim_to_S_vertex.Mag();
                unit_vec_prim_to_S_vertex = vec_prim_to_S_vertex.Unit();

                if(radiusS<5.){continue;}

                TLorentzVector tlv_K0_V01 = vec_K0_tlvs[K0_loop1];
                TVector3 dir_K0_V01;
                dir_K0_V01.SetXYZ(tlv_K0_V01.X(),tlv_K0_V01.Y(),tlv_K0_V01.Z());
                TVector3 unit_dir_K0_V01 = dir_K0_V01.Unit();

                TLorentzVector tlv_K0_V02 = vec_K0_tlvs[K0_loop2];
                TVector3 dir_K0_V02;
                dir_K0_V02.SetXYZ(tlv_K0_V02.X(),tlv_K0_V02.Y(),tlv_K0_V02.Z());
                TVector3 unit_dir_K0_V02 = dir_K0_V02.Unit();

                TVector3 vec_S_to_K01 = pos_K0_V01 - S_vertex_pos;
                TVector3 vec_S_to_K02 = pos_K0_V02 - S_vertex_pos;

                if(vec_S_to_K01.Mag()<0.5){continue;}
                if(vec_S_to_K02.Mag()<0.5){continue;}

                TVector3 unit_vec_S_to_K01 = vec_S_to_K01.Unit();
                TVector3 unit_vec_S_to_K02 = vec_S_to_K02.Unit();

                if ( unit_vec_S_to_K01.Dot(unit_dir_K0_V01) < 0.8){continue;}
                if ( unit_vec_S_to_K02.Dot(unit_dir_K0_V02) < 0.8){continue;}

                double dist1 = calculateMinimumDistanceStraightToPoint(pos_K0_V01,dir_K0_V01,S_vertex_pos);
                if(dist1>0.5) {continue;}

                double dist2 = calculateMinimumDistanceStraightToPoint(pos_K0_V02,dir_K0_V02,S_vertex_pos);
                if(dist2>0.5) {continue;}

                vector<Ali_AS_Track> all_tracks;
                vector<int> all_track_ids;

                all_tracks.push_back( vec_antip_and_pi_plus[2*S_vertex_loop] );
                all_tracks.push_back( vec_antip_and_pi_plus[2*S_vertex_loop + 1] );
                all_tracks.push_back( vec_K0_V0s[2*K0_loop1] );
                all_tracks.push_back( vec_K0_V0s[2*K0_loop1 + 1] );
                all_tracks.push_back( vec_K0_V0s[2*K0_loop2] );
                all_tracks.push_back( vec_K0_V0s[2*K0_loop2 + 1] );

                for(int i=0;i<6;i++)
                {
                    int trackid = all_tracks[i].gettrackid();
                    all_track_ids.push_back(trackid);

                }

                if ( check_if_value_is_doppelt_in_vector(all_track_ids) ){continue;}
                if (check_if_two_vectors_have_same_element(all_track_ids,brute_force_type4) ){continue;}    

                //calculate tlv
                TLorentzVector tlv_type4 = vec_tlv_type4[S_vertex_loop];
                tlv_type4 += tlv_K0_V01;
                tlv_type4 += tlv_K0_V02;
                tlv_type4 -= *tlv_neutron;

                TVector3 dir_type4;
                TVector3 unit_dir_type4;
                dir_type4.SetXYZ(tlv_type4[0],tlv_type4[1],tlv_type4[2]);
                unit_dir_type4 = dir_type4.Unit();

                if ( unit_vec_prim_to_S_vertex.Dot(unit_dir_type4) < 0.8 ) {continue;}




                Ali_AS_DM_particle* DMparticle = AS_Event->createDMparticle();

                DMparticle->settype(4);

                DMparticle->set_primVertex(pos_primary_vertex);
                DMparticle->set_S1Vertex(S_vertex_pos);
                DMparticle->set_S2Vertex(pos_K0_V01);
                DMparticle->set_S3Vertex(pos_K0_V02);

                DMparticle -> set_tlv(tlv_type4);

                for(int i=0;i<6;i++)
                {
                    Ali_AS_Track* dmtrack = DMparticle->createTrack();
                    Ali_AS_Track* ptr_all_tracks = &all_tracks[i];
                    copy_track_params(ptr_all_tracks,dmtrack);
                    //int trackid = all_tracks_channel3[i].gettrackid();
                    //cout<<"trackid: "<<trackid<<endl;

                }

                for(int i=0;i<6;i++)
                {
                    brute_force_type4.push_back(all_track_ids[i]);
                }





                //printf("type4; position x: %f, y: %f, z: %f \n",S_vertex_pos[0],S_vertex_pos[1],S_vertex_pos[2]);

            }

        }

    }











    AS_Event->setN_V0s(counter_V0s);
    //cout<<"NumV0s: "<<AS_Event->getN_V0s<<endl;

    int numtracks = AS_Event->getNumTracks();
    int id_of_track;

    /*
    for(int i_track=0;i_track<numtracks;i_track++)
    {
       Ali_AS_Track* track = AS_Event->getTrack(i_track);
        id_of_track = track -> gettrackid();
        if (check_if_int_is_in_vector(id_of_track , brute_force_type1) ) {continue;}
        //delete track;

    }
        */

    AS_Event->clearTrackList();



    Tree_AS_Event ->Fill();
    //Ali_AS_V0* AS_V0;

    numberV0 = AS_Event->getNumV0s();
    //cout<<"numberV0: "<<numberV0<<endl;
    for (Int_t V0_counter=0; V0_counter<numberV0; V0_counter++)
    {
        //if(counter_events>280) {cout<<"V0 counter2: "<<V0_counter<<endl;}
        AS_V0 =  AS_Event->getV0(V0_counter);
        //int numtracks1= AS_V0->getNumTracks();
        //cout<<"numtracks1:"<<numtracks1<<endl;
        AS_V0 -> clearTrackList();
        //int numtracks2= AS_V0->getNumTracks();
        //cout<<"numtracks2:"<<numtracks2<<endl;
        //AS_V0 -> clearV0List();
        AS_V0->~Ali_AS_V0();
    }
    //cout<<"a1"<<endl;
    int Ntracks = AS_Event->getN_tracks();

    if(numberV0<1)
    {
        as_trackP->~Ali_AS_Track();
        as_trackN->~Ali_AS_Track();
    }

    ASTrack1->~Ali_AS_Track();
    ASTrack2->~Ali_AS_Track();

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
    AS_Event ->clearNUCLEVList();

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


void Ali_DarkMatter_ESD_analysis::FindDCAHelixPoint2(TVector3 space_vec, Ali_AS_Track* helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB)
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

