
#include "Ali_Dark_Matter_Read.h"
#include "Ali_Dark_Matter_Read_LinkDef.h"

ClassImp(Ali_Dark_Matter_Read)



Double_t calcDeterminant(TVector3& v1,TVector3& v2,TVector3& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix 
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    . 
  // |                       |      |   \\     \\   / \\  | /     /     . 
  // |                       |      |    \\     \\ /   \\ |/     /      . 
  // |                       |      |     \\     X     \\/     /       . 
  // |                       |      |      \\   / \\    /\\    /        .  
  // |                       |      |       \\ /   \\  / |\\  /         . 
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          . 
  // |                       |      |      /   \\  /  \\ |/  \\         . 
  // |                       |      |     /     \\/    \\/    \\        . 
  // |                       |      |    /      /\\    /\\     \\       . 
  // |                       |      |   /      /  \\  / |\\     \\      .  
  // |                       |      |  /      /    \\/  | \\     \\     . 
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) .  
  //                                 /      /     /  \\     \\     \\   .
  //                                                                
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2) 
	   + v2(0) * v3(1) * v1(2) 
	   + v3(0) * v1(1) * v2(2) 
	   - v3(0) * v2(1) * v1(2) 
	   - v1(0) * v3(1) * v2(2) 
	   - v2(0) * v1(1) * v3(2)); 
}


TVector3 calculatePointOfClosestApproach(TVector3 &base1, TVector3 &dir1,
								    TVector3 &base2, TVector3 &dir2)
{
  //  calculating point of closest approach
  //        
  //        from the equations of the straight lines of g and h 
  //        g: x1 = base1 + l * dir1 
  //        h: x2 = base2 + m * dir2 
  //        
  //        you can construct the following planes:
  //        
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //        
  //        now the intersection point of E1 with g2 = {P1} 
  //        and the intersection point of E2 with g1 = {P2}
  //        
  //        form the base points of the perpendicular to both straight lines.
  //        
  //        The point of closest approach is the middle point between P1 and P2: 
  //        
  //        vertex = (p2 - p1)/2
  // 
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2 
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)       
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2        
  //                                          (b)
  //           
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)] 
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //           
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  //  
  //            m  = Dm/D12
  //           
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1 
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)       
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1        
  //                                          (t)
  //           
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)] 
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //           
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //           
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //           
  //           
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                      

  TVector3 Cross = dir1.Cross(dir2); // Cross product: dir1 x dir2

  // straight lines are either skew or have a Cross point
	      
  TVector3 diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
		
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,Cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and Cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");
      return TVector3(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, Cross);
  Double_t Dl = -calcDeterminant(diff , dir2, Cross);

  TVector3 vertex;
  TVector3 dm;
  TVector3 dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return TVector3(vertex);
}



TVector3 calculateCrossPoint(TVector3 &base1, TVector3 &dir1,
							TVector3 &base2, TVector3 &dir2)
{ 
  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
  
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
  
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating Cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");
      return TVector3(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
  
  Double_t lambda = Dlambda/D;
  
  TVector3 vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  return TVector3(vertex);

}


Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}


Double_t calculateMinimumDistance(TVector3 &base1, TVector3 &dir1,
							  TVector3 &base2, TVector3 &dir2)
{
  // calculates the minimum distance of two tracks given as parametric straights x = base + n * dir

  TVector3 cross = dir1.Cross(dir2);

  TVector3 ab = base1 - base2;

  if ( !( fabs(cross.Mag())>0.)) // dir1 || dir2
    {
      return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);
    }
 
  return fabs(ab.Dot(cross)/cross.Mag());
}


TVector3 calcVertexAnalytical(TVector3 &base1, TVector3 &dir1,
							 TVector3 &base2, TVector3 &dir2)
{
  // Calculates the Vertex of two straight lines define by the vectors base and dir
  //
  //      g: x1 = base1 + l * dir1 
  //      h: x2 = base2 + m * dir2 , where l,m are real numbers 
  //
  // 1. are g and h
  //       parallel / identical, i.e. are dir1 and dir2 linear dependent?
  //       
  //                                        /-                               
  //                                        |
  //                                        |   = 0    linear dependent, no unique solution, returning dummy  
  //      => Cross product : dir1 x dir2 = -|  
  //                                        |  != 0    linear independent
  //                                        |
  //                                        \\-         
  //
  // 2. are g and h 
  //       skew or do they have a Crossing point, i.e are dir1, dir2 and (base1 - base2) linear dependent ?
  //
  //                                                    /-                               
  //                                                    |
  //                                                    |   = 0    linear dependent
  //                                                    |          g and h are intersecting
  //                                                    |          calculating vertex as point of intersection
  //                                                    |
  //    => determinant: det[ dir1, dir2, base1-base2]= -|
  //                                                    |  != 0    linear independent
  //                                                    |          g and h are skew
  //                                                    |          calulating vertex as point of closest approach
  //                                                    |
  //                                                    \\-         
  //  
  // 3.
  //    (a) calculating intersection point
  //    (b) calculating point of closest approach





  // 1. exists a unique solution ?

  if ((dir1.Cross(dir2)).Mag()> 0.) // dir1 and dir2 linear independent
  {
     
      // straight lines are either skew or have a Cross point

      TVector3 diff = base1;
      diff-=base2; // Difference of two base vectors base1 - base2

      // 2. skew or intersecting ?
	
      if (fabs(calcDeterminant(dir2, dir1 ,diff))>0.) 
	{
            // 3. (b) skew
            //modified by Fabio--------------------------------------------------------------------------------------
            //if(calculateMinimumDistance(base1,dir1, base2, dir2)>1.){return TVector3(-7000000000.,-7000000000.,-7000000000.);}
            return TVector3(calculatePointOfClosestApproach(base1, dir1, base2, dir2));
            //-------------------------------------------------------------------------------------------------------------
	}
      else
	{
	  // 3. (a) intersection 
	  return TVector3(calculateCrossPoint(base1 ,dir1, base2 ,dir2));
	}
    }
  else
    {
      // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
      return TVector3(-10000000.,-10000000.,-10000000.);
    }
  return TVector3(-10000000.,-10000000.,-10000000.);
}




bool check_if_int_is_in_vector(int a, vector<int> vec)
{
    for(int i=0;i<(Int_t)vec.size();i++)
    {
        if(a == vec[i]) {return 1;}

    }
    return 0;

}

void print_int_vector(vector<int> vec)
{
    for(int i=0;i<(Int_t)vec.size();i++)
    {
        cout<<"Vektor  "<<i<<": "<<vec[i]<<endl;

    }
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

void fHelixAtoLinedca(TVector3 dirB, TVector3 spaceB, Ali_AS_Track* helixA, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{
    Float_t pA[2] = {0.0,0.0};
    Float_t pB[2] = {0.0,-5.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    Float_t path_closest_to_point = 0;
    Float_t dca_closest_to_point  = 0;
    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;
    double arr_testA[3];

    TVector3 testA, testB;
    for(Int_t r = 0; r < 2; r++)
    {
        testB     = spaceB+pB[r]*dirB;  // 3D-vector of helixB point at path pB[r]
        FindDCAHelixPoint(testB, helixA, path_initA,path_initB, path_closest_to_point, dca_closest_to_point);
        helixA->Evaluate(path_closest_to_point,arr_testA); // 3D-vector of helixA point at dca to testB
        testA.SetXYZ(arr_testA[0],arr_testA[1],arr_testA[2]);
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
            testB     = spaceB+pB[0]*dirB;  // 3D-vector of helixB point at path pB[r]
            FindDCAHelixPoint(testB, helixA, path_initA,path_initB, path_closest_to_point, dca_closest_to_point);
            pA[0]     = path_closest_to_point; // pathA at dca to helixB
            helixA->Evaluate(pA[0],arr_testA); // new vector testA
            testA.SetXYZ(arr_testA[0],arr_testA[1],arr_testA[2]);
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
            testB     = spaceB+pB[1]*dirB;  // 3D-vector of helixB point at path pB[r]
            FindDCAHelixPoint(testB, helixA, path_initA,path_initB, path_closest_to_point, dca_closest_to_point);
            pA[1]     = path_closest_to_point; // pathA at dca to helixB
            helixA->Evaluate(pA[1],arr_testA); // pathA at dca to helixB
            testA.SetXYZ(arr_testA[0],arr_testA[1],arr_testA[2]);
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



//----------------------------------------------------------------------------------------
template<typename T = double>
class Logspace {
private:
    T curValue, base, step;

public:
    Logspace(T first, T last, int num, T base = 10.0) : curValue(first), base(base){
       step = (last - first)/(num-1);
    }

    T operator()() {
        T retval = pow(base, curValue);
        curValue += step;
        return retval;
    }
};





//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Ali_Dark_Matter_Read::Ali_Dark_Matter_Read(TString list)
{
    TString in_list_name = list;
    TString outfile_dir = "/misc/alidata120/alice_u/schlichtmann/out/";
    TString outfile_name = outfile_dir + in_list_name + "_out.root";
    //outputfile = new TFile("./TRD_Calib.root","RECREATE");
    //outputfile = new TFile(outfile_name.Data(),"RECREATE");
    output_histos = new TFile(outfile_name.Data(),"RECREATE");
    // outputfile_trkl = new TFile("./TRD_Calib_on_trkl.root","RECREATE");
    //outputfile_histos = new TFile("Histos.root","RECREATE");

    AS_DM_particle = new Ali_AS_DM_particle();
    AS_DM_Track    = new Ali_AS_Track();

    Tree_AS_DM_particle  = NULL;
    Tree_AS_DM_particle  = new TTree("Tree_AS_DM_particle" , "Ali_AS_DM_particles" );
    //Tree_AS_DM_particle  = new TTree();
    Tree_AS_DM_particle  ->Branch("Tree_AS_DM_branch"  , "Ali_AS_DM_particle", AS_DM_particle );


}


//----------------------------------------------------------------------------------------
Ali_Dark_Matter_Read::~Ali_Dark_Matter_Read()
{

}
//----------------------------------------------------------------------------------------

#if 1

//----------------------------------------------------------------------------------------
void Ali_Dark_Matter_Read::Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;
    //TString pinputdir = "/home/ceres/schlichtmann/ESD_Analysis/";
    TString inlistdir = "/home/ceres/schlichtmann/ESD_Analysis/Lists/";

    TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_V3/";

    //TString pinputdir = "/misc/alidata120/alice_u/schlichtmann/dark_matter_V11/";
    //
    //TString pinputdir = "/misc/alidata120/alice_u/schlichtmann/dark_matter_PbPb_V4/";
    //TString pinputdir = "/misc/alidata120/alice_u/schlichtmann/dark_matter/";
    //TString pinputdir = "/home/ceres/berdnikova/TRD-Run3-Calibration/";

    TString in_list_name = SEList;
    SEList = inlistdir + SEList;

    AS_Event = new Ali_AS_Event();
    AS_V0    = new Ali_AS_V0();
    AS_Track = new Ali_AS_Track();
    //AS_Tracklet = new Ali_AS_Tracklet();
    //AS_Digit = new Ali_AS_TRD_digit();

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


    double start = -1;
    double stop = TMath::Log10(10);
    int num = 200;
    std::vector<double> vals;
    std::generate_n(std::back_inserter(vals), num, Logspace<>(start,stop,num));
    //for(double num : vals) std::cout << num << '\n';
    double x_vals[num];
    for(int i=0;i<(Int_t)vals.size();i++)
    {
        x_vals[i] = vals[i];
        //cout<<x_vals[i]<<endl;
    }

    start = 1;
    stop = TMath::Log10(4000);
    num = 200;
    std::vector<double> valsy;
    std::generate_n(std::back_inserter(valsy), num, Logspace<>(start,stop,num));
    //for(double num : vals) std::cout << num << '\n';
    double y_vals[num];
    for(int i=0;i<(Int_t)valsy.size();i++)
    {
        y_vals[i] = valsy[i];
        //cout<<y_vals[i]<<endl;
    }




   
    dEdx_vs_charge_dot_momentum = new TH2D("dEdx_vs_charge_dot_momentum","dEdx_vs_charge_dot_momentum",
                                           199,x_vals,199,y_vals);

    histos_1D.push_back(histo_invariantmass_lambda);
    histos_1D.push_back(histo_invariantmass_K0);
    histos_1D.push_back(histo_lambda_vertex_radius); //2
    histos_1D.push_back(histo_K0_vertex_radius);     //3
    histos_1D.push_back(histo_S_vertex_radius);      //4
    histos_1D.push_back(histo_S_mass_r_larger_10);   //5
    histos_1D.push_back(histo_S_mass_r_larger_20);   //6
    histos_1D.push_back(histo_num_tracks);      //7
    histos_1D.push_back(histo_vertex_z_pos);      //8
    histos_1D.push_back(histo_V0_with_K_plus_radius);      //9

    histos_2D.push_back(histo_lambda_x_and_y);
    histos_2D.push_back(histo_K0_x_and_y);
    histos_2D.push_back(histo_S_x_and_y);
    histos_2D.push_back(histo_V0_with_K_plus_x_and_y); //3

    counters.push_back(event_counter);
    counters.push_back(counter_of_2pions_close_to_S_vertex);
    counters.push_back(counter_of_V0s_of_antiproton_and_K_plus_and_another_K);
    counters.push_back(counter_of_S_vertices_without_pions);     //3
    counters.push_back(counter_1pion_close);     //4
    counters.push_back(counter_V0s_antiproton_and_K_plus);     //5
    counters.push_back(counter_correct_S_vertices);     //6

   
    for(int i=0;i<16;i++)
    {
        TString name="omega minus ";
        name+=i;
        TH1D* histo_invariantmass_omega_minus_baryon = new TH1D(name.Data(),name.Data(),200,1.5,1.8);
        vec_histo_omega_minus.push_back(histo_invariantmass_omega_minus_baryon);

        name="omega plus ";
        name+=i;
        TH1D* histo_invariantmass_omega_plus_baryon = new TH1D(name.Data(),name.Data(),200,1.5,1.8);
        vec_histo_omega_plus.push_back(histo_invariantmass_omega_plus_baryon);

    }

    //float arr_distance_prim_sec[10]{0.5,0.7,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7};
    //float arr_distance_daughter_particles[10]{0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.1,1.3};
    //float arr_dca_daughter_prim[10]{0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.1,0.2};
    for(int i=0;i<15;i++)
    {
        arr_distance_prim_sec[i]=6./15*(i+1);
        arr_dca_daughter_prim[i] = 1.5/15*(i+1);
        arr_dca_AB[i] = 0.1*(i+1);
    }

    for(int i=0;i<45;i++)
    {
        arr_radius_variation[i]=1.*(i+1);
    }


    for(int i=0;i<15*15;i++)
    {
        TString name = "radius >  ";
        name+=arr_distance_prim_sec[i/15];
        name+=" dcaV0 < ";
        name+=arr_distance_daughter_particles[0];
        name+=" dca_daughter_prim> ";
        name+=arr_dca_daughter_prim[i%15];
        TH1D* histo_invariantmass_lambda = new TH1D(name.Data(),name.Data(),50*2,1.1,1.13);
        vec_histo_invariantmass_Lambda.push_back(histo_invariantmass_lambda);

        TString name2 = "K0: radius >  ";
        name2+=arr_distance_prim_sec[i/15];
        name2+=" dcaV0 < ";
        name2+=arr_distance_daughter_particles[0];
        name2+=" dca_daughter_prim> ";
        name2+=arr_dca_daughter_prim[i%15];
        TH1D* histo_invariantmass_K0 = new TH1D(name2.Data(),name2.Data(),50*3,0.4,0.6);
        vec_histo_invariantmass_K0.push_back(histo_invariantmass_K0);

        //vec_histo_invariantmass_Lambda[top1*16+top2*4+top3]->Fill(invariantmass);
    }

    for(int i=0;i<45;i++)
    {
        TString name = "K0: radius >  ";
        name+= arr_radius_variation[i];
        name+=" dcaV0 < ";
        name+=1.;
        name+=" dca_daughter_prim> ";
        name+=0.2;
        TH1D* histo_invariantmass_K0 = new TH1D(name.Data(),name.Data(),50*3,0.4,0.6);
        vec_histo_radius_variation_invariantmass_K0.push_back(histo_invariantmass_K0);
    }

    for(int i=0;i<45;i++)
    {
        TString name = "Lambda: radius >  ";
        name+= arr_radius_variation[i];
        name+=" dcaV0 < ";
        name+=1.;
        name+=" dca_daughter_prim> ";
        name+=0.2;
        TH1D* histo_invariantmass_lambda = new TH1D(name.Data(),name.Data(),50*2,1.1,1.13);
        vec_histo_radius_variation_invariantmass_Lambda.push_back(histo_invariantmass_lambda);
    }

    for(int j=0;j<15*15;j++)
    {
        TString name = "Lambda radius > 45 ";
        name+=" dcaV0 < ";
        name+=arr_dca_AB[j/15];
        name+=" dca_daughter_prim> ";
        name+=arr_dca_daughter_prim[j%15];
        TH1D* histo_invariantmass_lambda = new TH1D(name.Data(),name.Data(),50*2,1.1,1.13);
        vec_histo_invariantmass_Lambda.push_back(histo_invariantmass_lambda);

        TString name2 = "K0 radius > 45 ";
        name2+=" dcaV0 < ";
        name2+=arr_dca_AB[j/15];
        name2+=" dca_daughter_prim> ";
        name2+=arr_dca_daughter_prim[j%15];
        TH1D* histo_invariantmass_K0 = new TH1D(name2.Data(),name2.Data(),50*3,0.4,0.6);
        vec_histo_invariantmass_K0.push_back(histo_invariantmass_K0);


    }

    


    


    
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_Dark_Matter_Read::copy_track_params(Ali_AS_Track* track_in, Ali_AS_Track* track_out)
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
//----------------------------------------------------------------------------------------

float calc_momentum(float* mom)
{
     return sqrt( mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2] );
}

float calc_momentum_squared(float* mom)
{
     return  mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2] ;
}

//----------------------------------------------------------------------------------------
Int_t Ali_Dark_Matter_Read::Loop_event(Long64_t event)
{
    if(event%100==0)
    {
        printf("Loop event number: %lld \n",event);
    }
    printf("Loop event number: %lld \n",event);
    counters[0]++;
    //cout<<""<<endl;


    Event_active = event;

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;
    //cout<<"a"<<endl;
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

    TVector3 pos_primary_vertex;
    pos_primary_vertex.SetXYZ(EventVertexX,EventVertexY,EventVertexZ);

     //cout<<"b"<<endl;

    //printf("Event vertex: %f", EventVertexX);
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    Int_t NumV0s = AS_Event ->getNumV0s();

    cout<<"NumV0s: "<<NumV0s<<endl;

    Float_t* pos = new Float_t[3];
    Float_t* momP = new Float_t[3];
    Float_t* momN = new Float_t[3];
    TLorentzVector* tlv_pos = new TLorentzVector();
    TLorentzVector* tlv_neg = new TLorentzVector();
    TLorentzVector* tlv_Lambda = new TLorentzVector();
    TLorentzVector* tlv_Kaon = new TLorentzVector();
    TLorentzVector* tlv_gamma = new TLorentzVector();
    Ali_AS_Track* as_trackP = new Ali_AS_Track;
    Ali_AS_Track* as_trackN = new Ali_AS_Track;
    Ali_AS_Track* as_Track5 = new Ali_AS_Track;
    Ali_AS_Track* as_tracksave = new Ali_AS_Track;
    Ali_AS_Track* ASTrack1 = new Ali_AS_Track;
    Ali_AS_Track* ASTrack2 = new Ali_AS_Track;
    Ali_AS_Track* AS_TrackA = new Ali_AS_Track;
    Ali_AS_Track* AS_TrackB = new Ali_AS_Track;
    Ali_AS_Track* AS_TrackC = new Ali_AS_Track;
    Ali_AS_Track* AS_TrackD = new Ali_AS_Track;
    Ali_AS_Track* track_in_loop  =  new Ali_AS_Track;
    Ali_AS_Track* tracka  =  new Ali_AS_Track;
    Ali_AS_Track* trackb  =  new Ali_AS_Track;
    Float_t* sigma_proton_TPC = new Float_t[2];
    Float_t* sigma_pion_TPC = new Float_t[2];
    Float_t radius;
    TVector3 position_SV2;
    TVector3 position_SV3;
    TVector3 direction_SV2;
    TVector3 direction_SV3;

    TLorentzVector tlv_in_loop;

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

    const Float_t mass_proton = 0.93827208816 ;  //in GeV?
    const Float_t mass_pion = 0.139657061 ;  //in GeV?
    const Float_t mass_electron = 0.510998950 * 1e-3 ;  //in GeV?
    const double  mass_K = 0.493677 ;  //in GeV?

    vector<int> used_track_ids_of_pions;
    vector<int> used_track_ids_V0;
    vector<int> all_used_positive_track_ids_for_V0s;
    vector<int> all_used_negative_track_ids_for_V0s;

    double dcaP=-10000;
    double dcaN=-10000;

    Int_t trackidP,trackidN;

    double momentumP,momentumN;

    Double_t invariantmass = -1.;

    double sigma_K_plus_TPC;

    Float_t energy_proton,energy_pion,energy_antiproton,energy_pion_plus,energy_pion_minus, energy_K_plus,energy_anti_proton;
    Float_t energy_electron_plus,energy_electron_minus;

    TLorentzVector tlv_anti_p_and_K_plus;

    vector<Ali_AS_Track*> vec_SV2_tracks;
    vector<Ali_AS_Track*> vec_SV3_tracks;

    vector<int> vec_SV2_track_ids;
    vector<int> vec_SV3_track_ids;


    vector<int> brute_force;

    vector<int> trackids_bits_shared;
    vector<int> trackids_similiar_tracks;
    vector<int> trackids_similiar_used;

     //cout<<"c"<<endl;

    //find similiar tracks

    int counter_tracks=0;
    /*
    for(Int_t atrack = 0; atrack < NumTracks; atrack++)
    {
        tracka =  AS_Event -> getTrack(atrack);
        TLorentzVector tlva = tracka->get_TLV_part();
        int trackida = tracka->gettrackid();
        if ( check_if_int_is_in_vector(trackida,trackids_similiar_used)){continue;}

        for(int btrack=0; btrack<NumTracks; btrack++)
        {
            if(btrack==atrack){continue;}

            trackb = AS_Event -> getTrack(btrack);
            int trackidb = trackb->gettrackid();
            if ( check_if_int_is_in_vector(trackidb,trackids_similiar_used)){continue;}

            TLorentzVector tlvb = trackb->get_TLV_part();

             TBits tbitsshared = tracka->getbitsshared();
             int numbitsshared = tbitsshared.CountBits();


            if( fabs( 1-tlva[0]/tlvb[0]) < 0.005 && fabs( 1-tlva[1]/tlvb[1]) < 0.005 && fabs( 1-tlva[2]/tlvb[2]) < 0.005)
            {
                
                printf("tracka: %d, trackb %d \n",trackida,trackidb);
                cout<<"numbitsshared: "<<numbitsshared<<endl;
                cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

                Float_t path_closest_to_point = 0;
                Float_t dca_closest_to_point  = 0;
                Float_t path_initA = 0.0;
                Float_t path_initB = 30.0;
                double  pos[3];
                float dca=-1;

                tracka->Evaluate(0,pos);
                TVector3 position;
                position[0]=pos[0];
                position[1]=pos[1];
                position[2]=pos[2];
                FindDCAHelixPoint(position,trackb,path_initA,path_initB,path_closest_to_point,dca);
                cout<<"dca: "<<dca<<endl;
                if(dca<0.5)
                {
                    cout<<"pushed back"<<endl;
                    trackids_similiar_tracks.push_back(trackidb);

                    trackids_similiar_used.push_back(trackida);
                    trackids_similiar_used.push_back(trackidb);
                }

                AS_DM_particle ->set_primVertex(pos_primary_vertex);
                AS_DM_particle ->clearTrackList();
                AS_DM_Track = AS_DM_particle ->createTrack();
                copy_track_params(tracka,AS_DM_Track);
                AS_DM_Track = AS_DM_particle ->createTrack();
                copy_track_params(trackb,AS_DM_Track);
                //cout<<"a"<<endl;
                Tree_AS_DM_particle ->Fill();
                //cout<<"b"<<endl;
            }


            /*
            if(atrack==10 && btrack ==15 && counter_tracks<10)
            {
                //printf("tracka: %d, trackb %d \n",trackida,trackidb);
                //cout<<"numbitsshared: "<<numbitsshared<<endl;
                //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
                AS_DM_particle ->set_primVertex(pos_primary_vertex);
                AS_DM_particle ->clearTrackList();
                AS_DM_Track = AS_DM_particle ->createTrack();
                copy_track_params(tracka,AS_DM_Track);
                AS_DM_Track = AS_DM_particle ->createTrack();
                copy_track_params(trackb,AS_DM_Track);

                Tree_AS_DM_particle ->Fill();

                counter_tracks++;

            }




        }

    }
    */
     //cout<<"d"<<endl;

    //loop over all tracks
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        //histo_counter->Fill(7.5);

        track_in_loop = AS_Event -> getTrack(i_track);

        int trackid = track_in_loop->gettrackid();
        TBits tbitsshared = track_in_loop->getbitsshared();
        int numbitsshared = tbitsshared.CountBits();
        if(numbitsshared>10)
        {
            //histo_counter->Fill(8.5);
            trackids_bits_shared.push_back(trackid);
        }

    }
    //cout<<"a"<<endl;
    
    //-------------------------------------------------------------------------------------------------------------------------------
    //loop over all tracks of event in order to make Bethe Bloch plot: dEdx as function of charge*momentum
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        track_in_loop = AS_Event -> getTrack(i_track);
        int trackid = track_in_loop->gettrackid();
        if( check_if_int_is_in_vector(trackid,trackids_similiar_tracks) )
        {
            counter_skipped_track++;
            continue;
        }


        Float_t TPCdEdx   = track_in_loop->getTPCdEdx();
        Float_t tofsignal = track_in_loop->getTOFsignal();
        Float_t dca       = track_in_loop->getdca();
        Float_t tracklength = track_in_loop->getTrack_length();
        int charge;

        if(dca>0){charge = 1;}
        else {charge = -1;}

        tlv_in_loop = track_in_loop->get_TLV_part();
        double momentum = tlv_in_loop.P();

        //printf("dEdx: %f, momentum: %f, dca: %f, charge: %d, tofsignal: %f \n"
          //     ,TPCdEdx, momentum,dca, charge, tofsignal);

        int num_points = gr ->GetN();
        gr ->SetPoint(num_points,charge * momentum, TPCdEdx);
        dEdx_vs_charge_dot_momentum->Fill(charge * momentum, TPCdEdx);

        if(tofsignal>99990){continue;}
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

        //if(m_squared<0){printf("mass squared: %f \n", m_squared);}
        if(isnan(m_squared)) {continue;}

        int num_points2 = gr2 ->GetN();
        gr2 ->SetPoint(num_points,charge * momentum, m_squared);

        double m_squared2 =   calculate_m_squared_by_TOF(track_in_loop);

        mass_squared_vs_charge_dot_momentum->Fill( charge * momentum , m_squared);

        if ( fabs( track_in_loop -> getnsigma_K_TPC()) < 1.0 )
        {
            //mass_squared_vs_charge_dot_momentum_kaons->Fill(charge * momentum , m_squared);
        }

        //cut auf pi - und momentum<0.4

        
        if( momentum<0.4 && fabs( track_in_loop -> getnsigma_pi_TPC()) > 1.5 )
        {
          mass_squared_no_pions->Fill(m_squared);
        }

        if(fabs( track_in_loop-> getnsigma_K_TPC())<2.5 && tofsignal<99990)
        {
            mass_squared_kaons->Fill(m_squared2);
        }



    }

     //cout<<"f"<<endl;
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------

    bool V0_is_used = 0;

    vector<vector<int>> track_ids_all;
    track_ids_all.resize(NumV0s);

  
    //cout<<"b"<<endl;
    //loop over V0s
    for(Int_t V0counter = 0; V0counter < NumV0s; V0counter++)
    {
       // cout<<"c"<<endl;
        V0_is_used =0;

        AS_V0 = AS_Event -> getV0(V0counter);
        //cout<<"d"<<endl;
        //get position of V0
        pos = AS_V0 -> getxyz();
        //cout<<"d1"<<endl;
        ////cout<<"posx: "<<pos[0]<<endl;
        //position.SetXYZ(pos[0],pos[1],pos[2]);
        radius = sqrt( (pos[0]-EventVertexX) *(pos[0]-EventVertexX)+(pos[1]-EventVertexY)*(pos[1]-EventVertexY)+(pos[2]-EventVertexZ)*(pos[2]-EventVertexZ) );
        // cout<<"d2"<<endl;
        TVector3 vec_primtoV0;
        vec_primtoV0.SetXYZ((pos[0]-EventVertexX),(pos[1]-EventVertexY),(pos[2]-EventVertexZ));
        TVector3 unit_prim_to_V0;
        unit_prim_to_V0 = vec_primtoV0.Unit();
        //printf("x %f,y %f, z %f \n",pos[0],pos[1],pos[2]);
        //cout<<"d3"<<endl;

        //get momentum
        momP = AS_V0 -> getPpxpypz();
        momN = AS_V0 -> getNpxpypz();

        //printf("pxP: %f,pyP: %f, pzP: %f \n",momP[0],momP[1],momP[2]);
        //printf("pxN: %f,pyN: %f, pzN: %f \n",momN[0],momN[1],momN[2]);

        energy_antiproton = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
        energy_pion       = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
        //cout<<"hi: "<<mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
        //cout<<"masspion: "<<mass_pion<<endl;
        //cout<<"energy_antiproton: "<<energy_antiproton<<endl;
        //cout<<"energy_pion: "<<energy_pion<<endl;

        tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
        tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

        *tlv_Lambda = *tlv_pos + *tlv_neg;
        invariantmass = tlv_Lambda->M();

        //cout<<"invariantmass: "<<invariantmass<<endl;
        //cout<<"momentumN: "<<momN[0]<<endl;
        //cout<<"d4"<<endl;
        double dcaV0 = AS_V0 -> getdcaV0();
        //printf("momentum of negative particle px: %f,py: %f,pz: %f \n",momN[0],momN[1],momN[2]);

        //--------------------------------------------------------------------------------
        //get two tracks for each V0
        as_trackP = AS_V0 -> getTrack(0);
        as_trackN = AS_V0 -> getTrack(1);
        dcaP = as_trackP->getdca();
        dcaN = as_trackN->getdca();
        //printf("dcaP: %f, dcaN: %f \n",dcaP,dcaN);

        //cout<<"e"<<endl;
        // Double check that the charge is correct
        if(dcaP < 0){continue;}
        if(dcaN > 0){continue;}

        trackidP = as_trackP->gettrackid();
        trackidN = as_trackN->gettrackid();
        //printf("trackidN: %d \n",trackidN);

        momentumP = sqrt( momP[0]*momP[0] + momP[1]*momP[1]+ momP[2]*momP[2] );
        momentumN = sqrt( momN[0]*momN[0] + momN[1]*momN[1]+ momN[2]*momN[2] );

        //printf("trackidP %d, trackidN %d \n",trackidP,trackidN )  ;
        //printf("momentum of positive particle px: %f,py: %f,pz: %f   momentum of negative particle px: %f,py: %f,pz: %f \n",momP[0],momP[1],momP[2],momN[0],momN[1],momN[2]);
        //printf("momentum of positive particle: %f   momentum of negative particle: %f \n",momentumP,momentumN);

        all_used_positive_track_ids_for_V0s.push_back(trackidP);
        all_used_negative_track_ids_for_V0s.push_back(trackidN);
        ///if(event%100==0 && V0counter%10 ==0)
        // {
        // printf("trackidP %d, trackidN %d \n",trackidP,trackidN )  ;
        // }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        //_____________________________________________________________________----

        sigma_proton_TPC[0] = as_trackP -> getnsigma_p_TPC();
        sigma_proton_TPC[1] = as_trackN -> getnsigma_p_TPC();
        //printf("sigmas proton %f  %f", sigma_proton_TPC[0], sigma_proton_TPC[1]) ;



        sigma_pion_TPC[0]  = as_trackP -> getnsigma_pi_TPC();
        sigma_pion_TPC[1]  = as_trackN -> getnsigma_pi_TPC();

        sigma_K_plus_TPC   = as_trackP -> getnsigma_K_TPC();

        Float_t path_closest_to_point = 0;
        Float_t dca_closest_to_point  = 0;
        Float_t dca_closest_to_pointP  = 0;
        Float_t dca_closest_to_pointN  = 0;
        Float_t path_initA = 0.0;
        Float_t path_initB = 30.0;

        float radiuscuts[4]{1,3,5,7};
        float dcaprimcuts[4]{0.5,1,2,5};
        // cout<<"d"<<endl;

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //test topology cuts
      
        //float arr_distance_prim_sec[10]{0.3,0.5,0.7,0.9,1.1,1.2,1.3,1.4,1.5,1.6};
        //float arr_distance_daughter_particles[10]{0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.1,1.3};
        //float arr_dca_daughter_prim[10]{0.02,0.03,0.035,0.04,0.045,0.05,0.055,0.1,0.2,0.7};
        // float arr_distance_prim_sec[10]{1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85};
        //float arr_distance_daughter_particles[10]{0.35,0.4,0.45,0.47,0.5,0.55,0.6,0.65,0.7,0.75};
        //float arr_dca_daughter_prim[10]{0.04,0.042,0.044,0.046,0.048,0.05,0.052,0.054,0.056,0.058};
        for(int i=0;i<45;i++)
        {
            if(radius<arr_radius_variation[i]){continue;}
            if(dcaV0>1.){continue;}
            if(fabs(dcaP) < 0.2){continue;}
            if(fabs(dcaN) < 0.2){continue;}

            if(fabs(sigma_pion_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
            {
                energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
                tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

                *tlv_Kaon = *tlv_pos + *tlv_neg;
                TVector3 dir_K0;
                dir_K0.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                TVector3 unit_dir_K0 = dir_K0.Unit();
                double dot_product = unit_dir_K0.Dot(unit_prim_to_V0);
                if(dot_product<0.){continue;}

                invariantmass = tlv_Kaon->M();
                vec_histo_radius_variation_invariantmass_K0[i]->Fill(invariantmass);

            }

            if(fabs(sigma_proton_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
            {
                energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                energy_pion   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
                tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

                *tlv_Lambda = *tlv_pos + *tlv_neg;
                TVector3 dir_Lambda;
                dir_Lambda.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                TVector3 unit_dir_Lambda = dir_Lambda.Unit();
                double dot_product = unit_dir_Lambda.Dot(unit_prim_to_V0);
                if(dot_product<0.){continue;}

                invariantmass = tlv_Lambda->M();
                vec_histo_radius_variation_invariantmass_Lambda[i]->Fill(invariantmass);
               // cout<<"filled"<<endl;
            }

        }

        for(int dca_AB = 0; dca_AB<15; dca_AB++)
        {
            for(int dca_daughter_prim=0;dca_daughter_prim<15;dca_daughter_prim++)
            {
                if(radius<45){continue;}
                if(dcaV0>arr_dca_AB[dca_AB]){continue;}
                if(fabs(dcaP) < arr_dca_daughter_prim[dca_daughter_prim]){continue;}
                if(fabs(dcaN) < arr_dca_daughter_prim[dca_daughter_prim]){continue;}

                int i_val = 15*15+dca_AB*15+dca_daughter_prim;

                if(fabs(sigma_proton_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
                {
                    energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                    energy_pion   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                    tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
                    tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

                    *tlv_Lambda = *tlv_pos + *tlv_neg;
                    TVector3 dir_Lambda;
                    dir_Lambda.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                    TVector3 unit_dir_Lambda = dir_Lambda.Unit();
                    double dot_product = unit_dir_Lambda.Dot(unit_prim_to_V0);
                    if(dot_product<0.){continue;}

                    invariantmass = tlv_Lambda->M();
                    vec_histo_invariantmass_Lambda[i_val]->Fill(invariantmass);

                }

                //PID K0
                if(fabs(sigma_pion_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
                {
                    energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                    energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                    tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
                    tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

                    *tlv_Kaon = *tlv_pos + *tlv_neg;
                    TVector3 dir_K0;
                    dir_K0.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                    TVector3 unit_dir_K0 = dir_K0.Unit();
                    double dot_product = unit_dir_K0.Dot(unit_prim_to_V0);
                    if(dot_product<0.){continue;}

                    invariantmass = tlv_Kaon->M();
                    vec_histo_invariantmass_K0[i_val]->Fill(invariantmass);

                }


            }
        }


        for(int top1 = 0; top1<15; top1++)
        {
            for(int top2=0;top2<1;top2++)
            {
                for(int top3=0;top3<15;top3++)
                {
                    float distance_prim_sec = arr_distance_prim_sec[top1];
                    float distance_daughter_particles = arr_distance_daughter_particles[top2];
                    float dca_daughter_prim = arr_dca_daughter_prim[top3];

                    if(radius < distance_prim_sec){continue;}
                    if(dcaV0 > distance_daughter_particles){continue;}
                    if(fabs(dcaP) < dca_daughter_prim){continue;}
                    if(fabs(dcaN) < dca_daughter_prim){continue;}

                    int i_val = top1*15+top3;

                    //PID for Lambdas
                    if(fabs(sigma_proton_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
                    {
                        energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                        energy_pion   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                        tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
                        tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

                        *tlv_Lambda = *tlv_pos + *tlv_neg;
                        TVector3 dir_Lambda;
                        dir_Lambda.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                        TVector3 unit_dir_Lambda = dir_Lambda.Unit();
                        double dot_product = unit_dir_Lambda.Dot(unit_prim_to_V0);
                        if(dot_product<0.){continue;}

                        invariantmass = tlv_Lambda->M();
                        vec_histo_invariantmass_Lambda[i_val]->Fill(invariantmass);

                    }

                    //PID K0
                    if(fabs(sigma_pion_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
                    {
                        energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
                        energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
                        tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
                        tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

                        *tlv_Kaon = *tlv_pos + *tlv_neg;
                        TVector3 dir_K0;
                        dir_K0.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                        TVector3 unit_dir_K0 = dir_K0.Unit();
                        double dot_product = unit_dir_K0.Dot(unit_prim_to_V0);
                        if(dot_product<0.){continue;}

                        invariantmass = tlv_Kaon->M();
                        vec_histo_invariantmass_K0[i_val]->Fill(invariantmass);

                    }
                }


            }


        }


      

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        //Lambda0 -> proton + pi-
        //check if positive particle is proton and if negative particle is pion-
        if(fabs(sigma_proton_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
        {
            //printf("particles are proton and pion- \n");
           

            energy_proton = sqrt(mass_proton*mass_proton+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion   = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_proton);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();
            //cout<<invariantmass<<endl;
            //histos_1D[0]->Fill(invariantmass);
            if(radius>45)
            {
                histo_invariantmass_lambda->Fill(invariantmass);
            }

            histos_1D[2]->Fill(radius);

            histos_2D[0]->Fill(pos[0],pos[1]);



            //cut on mass
            if(invariantmass < (1.1157+0.001495*2) && invariantmass > (1.1157-0.001495*2))
            {
                counter_lambdas++;
                V0_is_used = 1;

                used_track_ids_of_pions.push_back(trackidN);
                //printf(" proton and pi- ; trackidP %u, trackidN %u \n",trackidP,trackidN )  ;
                //position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                //vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                //vec_SV3_tracks.push_back(as_trackP);
                //vec_SV3_tracks.push_back(as_trackN);

                //vec_SV3_track_ids.push_back(trackidP);
                //vec_SV3_track_ids.push_back(trackidN);


                //search for Xi-baryons xi- -> pi- + lambda ; lambda -> proton + pi-
                //search for Xi-baryons xi+ -> pi+ + lambda ; lambda -> proton + pi-
                //search for Omega-baryons omega- -> K- + lambda ; lambda -> proton + pi-
                //search for Omega-baryons omega+ -> K+ + lambda ; lambda -> proton + pi-
                for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
                {
                    AS_Track = AS_Event->getTrack(i_track_A);
                    Int_t trackid  = AS_Track->gettrackid();
                    float dca = AS_Track->getdca();
                    int charge;
                    if(dca>0) {charge=1;}
                    else {charge = -1;}

                    if(trackid==trackidP || trackid==trackidN){continue;}

                    double sigma_pi = AS_Track->getnsigma_pi_TPC();
                    double sigma_K = AS_Track->getnsigma_K_TPC();

                    Float_t pathA,pathB,dcaAB;

                    fHelixAtoLinedca(direction_SV3, position_SV3, AS_Track, pathA, pathB, dcaAB);
                    if(dcaAB>0.5){continue;}

                    float dcaprim;
                    FindDCAHelixPoint(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dcaprim);
                    if(dcaprim>0.5)
                    {
                        //xi- and x+
                        if(fabs(sigma_pi)<2.5)
                        {
                            if(radius<10){continue;}
                            TLorentzVector tlv;
                            TLorentzVector tlv_pion;
                            TLorentzVector tlv_xi;

                            tlv = AS_Track->get_TLV_part();
                            double mom_pion = tlv.P();
                            double energy_pion = sqrt( mass_pion*mass_pion +  mom_pion*mom_pion );
                            Double_t r1[3];
                            Double_t r2[3];

                            AS_Track->Evaluate(pathA,r1);
                            AS_Track->Evaluate(pathA+0.01,r2);

                            TVector3 mom_dir_pion;
                            TVector3 unit_mom_dir_pion;
                            TVector3 vec_momentum;

                            mom_dir_pion.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
                            unit_mom_dir_pion = mom_dir_pion.Unit();
                            vec_momentum.SetXYZ(unit_mom_dir_pion[0]*mom_pion,unit_mom_dir_pion[1]*mom_pion,unit_mom_dir_pion[2]*mom_pion);

                            tlv_pion.SetPxPyPzE(vec_momentum[0],vec_momentum[1],vec_momentum[2],energy_pion);

                            tlv_xi = *tlv_Lambda + tlv_pion;

                            double xi_mass = tlv_xi.M();

                            if(charge==-1) {histo_invariantmass_xi_minus_baryon->Fill(xi_mass);}
                            //if(charge==1) {histo_invariantmass_xi_plus_baryon->Fill(xi_mass);cout<<"filledminus"<<endl;}

                        }

                    }

                    for(int i =0;i<4;i++)
                    {
                        for(int j=0;j<4;j++)
                        {
                            if(dcaprim>dcaprimcuts[i])
                            {

                                if(fabs(sigma_K)<2.5)
                                {
                                    if(radius<radiuscuts[j]){continue;}
                                    TLorentzVector tlv;
                                    TLorentzVector tlv_kaon;
                                    TLorentzVector tlv_omega;

                                    tlv = AS_Track->get_TLV_part();
                                    double mom_kaon = tlv.P();
                                    double energy_kaon = sqrt( mass_K*mass_K +  mom_kaon*mom_kaon );
                                    Double_t r1[3];
                                    Double_t r2[3];

                                    AS_Track->Evaluate(pathA,r1);
                                    AS_Track->Evaluate(pathA+0.01,r2);

                                    TVector3 mom_dir_kaon;
                                    TVector3 unit_mom_dir_kaon;
                                    TVector3 vec_momentum_kaon;

                                    mom_dir_kaon.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
                                    unit_mom_dir_kaon = mom_dir_kaon.Unit();
                                    vec_momentum_kaon.SetXYZ(unit_mom_dir_kaon[0]*mom_kaon,unit_mom_dir_kaon[1]*mom_kaon,unit_mom_dir_kaon[2]*mom_kaon);

                                    tlv_kaon.SetPxPyPzE(vec_momentum_kaon[0],vec_momentum_kaon[1],vec_momentum_kaon[2],energy_kaon);

                                    tlv_omega = *tlv_Lambda + tlv_kaon;

                                    double omega_mass = tlv_omega.M();

                                    if(charge==-1) {vec_histo_omega_minus[i*4+j]->Fill(omega_mass);}
                                    //if(charge==1) {histo_invariantmass_omega_plus_baryon->Fill(omega_mass);}

                                }
                            }
                        }
                    }

                    //if(fabs(sigma_pi)>2.5){continue;}

                   
                    //cout<<"dca: "<<dcaAB<<endl;

                    
                    //find dca helix to straight
                    //helix is AS_Track
                    //straight: base = pos , direction=direction_SV3

                    //if(dca> 0.5) {continue;}


                    

                    //cout<<"xi mass"<<xi_mass<<endl;


                 }




            }
        }
         //cout<<"h"<<endl;
        //if (V0_is_used ==1 )  {continue;}



        //check if negative particle is antiproton and if positive particle is pion+
        //antilambda
        if(fabs(sigma_proton_TPC[1]) < 2.5 && fabs(sigma_pion_TPC[0]) < 2.5)
        {
            //printf("particles are antiproton and pion+ \n");
            energy_antiproton = sqrt(mass_proton*mass_proton+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            energy_pion       = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_antiproton);

            *tlv_Lambda = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Lambda->M();

            histo_invariantmass_anti_lambda->Fill(invariantmass);
            // cout<<invariantmass<<endl;
            //histos_1D[0]->Fill(invariantmass);

            histos_1D[2]->Fill(radius);

            histos_2D[0]->Fill(pos[0],pos[1]);


            //cut on mass
            if(invariantmass< (1.1157+0.001495*2) && invariantmass > (1.1157-0.001495*2))
            {
                counter_anti_lambdas++;
                V0_is_used = 1;
                used_track_ids_of_pions.push_back(trackidP);
                //printf(" antiproton and pi+ ; trackidP %u, trackidN %u \n",trackidP,trackidN )  ;
                position_SV3.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV3.push_back(position_SV3);

                direction_SV3.SetXYZ(tlv_Lambda->Px(),tlv_Lambda->Py(),tlv_Lambda->Pz());
                vec_direction_SV3.push_back(direction_SV3);

                vec_SV3_tracks.push_back(as_trackP);
                vec_SV3_tracks.push_back(as_trackN);

                vec_SV3_track_ids.push_back(trackidP);
                vec_SV3_track_ids.push_back(trackidN);

                //search for Xi-baryons xi- -> pi- + lambda ; lambda -> proton + pi-
                //search for Xi-baryons xi+ -> pi+ + lambda ; lambda -> proton + pi-
                //search for Omega-baryons omega- -> K- + lambda ; lambda -> proton + pi-
                //search for Omega-baryons omega+ -> K+ + lambda ; lambda -> proton + pi-
                for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
                {
                    AS_Track = AS_Event->getTrack(i_track_A);
                    Int_t trackid  = AS_Track->gettrackid();
                    float dca = AS_Track->getdca();
                    int charge;
                    if(dca>0) {charge=1;}
                    else {charge = -1;}

                    if(trackid==trackidP || trackid==trackidN){continue;}

                    Float_t pathA,pathB,dcaAB;

                    float dcaprim;
                    FindDCAHelixPoint(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dcaprim);
                   

                    fHelixAtoLinedca(direction_SV3, position_SV3, AS_Track, pathA, pathB, dcaAB);
                    if(dcaAB>0.5){continue;}

                    //if(charge>0){continue;}

                    double sigma_pi = AS_Track->getnsigma_pi_TPC();
                    double sigma_K = AS_Track->getnsigma_K_TPC();

                    if(dcaprim>0.5)
                    {
                        //xi- and x+
                        if(fabs(sigma_pi)<2.5)
                        {
                            if(radius<10){continue;}
                            TLorentzVector tlv;
                            TLorentzVector tlv_pion;
                            TLorentzVector tlv_xi;

                            tlv = AS_Track->get_TLV_part();
                            double mom_pion = tlv.P();
                            double energy_pion = sqrt( mass_pion*mass_pion +  mom_pion*mom_pion );
                            Double_t r1[3];
                            Double_t r2[3];

                            AS_Track->Evaluate(pathA,r1);
                            AS_Track->Evaluate(pathA+0.01,r2);

                            TVector3 mom_dir_pion;
                            TVector3 unit_mom_dir_pion;
                            TVector3 vec_momentum;

                            mom_dir_pion.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
                            unit_mom_dir_pion = mom_dir_pion.Unit();
                            vec_momentum.SetXYZ(unit_mom_dir_pion[0]*mom_pion,unit_mom_dir_pion[1]*mom_pion,unit_mom_dir_pion[2]*mom_pion);

                            tlv_pion.SetPxPyPzE(vec_momentum[0],vec_momentum[1],vec_momentum[2],energy_pion);

                            tlv_xi = *tlv_Lambda + tlv_pion;

                            double xi_mass = tlv_xi.M();

                            //if(charge==-1) {histo_invariantmass_xi_minus_baryon->Fill(xi_mass); cout<<"filledplus"<<endl;}
                            if(charge==1) {histo_invariantmass_xi_plus_baryon->Fill(xi_mass);}

                        }
                    }

                    for(int i =0;i<4;i++)
                    {
                        for(int j=0;j<4;j++)
                        {
                            if(dcaprim>dcaprimcuts[i])
                            {
                                //omega + and omega-
                                if(fabs(sigma_K)<2.5)
                                {
                                    if(radius<radiuscuts[j]){continue;}
                                    TLorentzVector tlv;
                                    TLorentzVector tlv_kaon;
                                    TLorentzVector tlv_omega;

                                    tlv = AS_Track->get_TLV_part();
                                    double mom_kaon = tlv.P();
                                    double energy_kaon = sqrt( mass_K*mass_K +  mom_kaon*mom_kaon );
                                    Double_t r1[3];
                                    Double_t r2[3];

                                    AS_Track->Evaluate(pathA,r1);
                                    AS_Track->Evaluate(pathA+0.01,r2);

                                    TVector3 mom_dir_kaon;
                                    TVector3 unit_mom_dir_kaon;
                                    TVector3 vec_momentum_kaon;

                                    mom_dir_kaon.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
                                    unit_mom_dir_kaon = mom_dir_kaon.Unit();
                                    vec_momentum_kaon.SetXYZ(unit_mom_dir_kaon[0]*mom_kaon,unit_mom_dir_kaon[1]*mom_kaon,unit_mom_dir_kaon[2]*mom_kaon);

                                    tlv_kaon.SetPxPyPzE(vec_momentum_kaon[0],vec_momentum_kaon[1],vec_momentum_kaon[2],energy_kaon);

                                    tlv_omega = *tlv_Lambda + tlv_kaon;

                                    double omega_mass = tlv_omega.M();

                                    //if(charge==-1) {histo_invariantmass_omega_minus_baryon->Fill(omega_mass);}
                                    if(charge==1) {vec_histo_omega_plus[i*4+j]->Fill(omega_mass);}
                                }
                            }
                        }
                    }
                }

                


            }
        }
         //cout<<"i"<<endl;
        //if(V0_is_used==1){continue;}

        //K0 - > pi+ and pi-
        //check if pion+ and pion-
        if(fabs(sigma_pion_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
        {
            V0_is_used=1;
            //printf("particles are pion+ and pion- \n");
            energy_pion_plus  = sqrt(mass_pion*mass_pion+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_pion_minus = sqrt(mass_pion*mass_pion+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            //cout<<energy_proton<<endl;
            //cout<<sqrt(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2])<<endl;
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_pion_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_pion_minus);

            *tlv_Kaon = *tlv_pos + *tlv_neg;
            invariantmass = tlv_Kaon->M();
            //cout<<invariantmass<<endl;
            if(radius>45)
            {
                histo_invariantmass_K0->Fill(invariantmass);
            }
            //histos_1D[1]->Fill(invariantmass);

            histos_1D[3]->Fill(radius);

            histos_2D[1]->Fill(pos[0],pos[1]);

            //cut on mass
            if(invariantmass< (0.4981+0.0042*2) && invariantmass > (0.4981-0.0042*2))
            {

                used_track_ids_of_pions.push_back(trackidP);
                used_track_ids_of_pions.push_back(trackidN);
                //printf("particles are pion+ and pion- \n");

                position_SV2.SetXYZ(pos[0],pos[1],pos[2]);
                vec_position_SV2.push_back(position_SV2);

                direction_SV2.SetXYZ(tlv_Kaon->Px(),tlv_Kaon->Py(),tlv_Kaon->Pz());
                vec_direction_SV2.push_back(direction_SV2);

                vec_SV2_tracks.push_back(as_trackP);
                vec_SV2_tracks.push_back(as_trackN);

                vec_SV2_track_ids.push_back(trackidP);
                vec_SV2_track_ids.push_back(trackidN);

                //printf("position of vertex:  %f %f %f \n",pos[0],pos[1],pos[2]);
                //cout<<""<<endl;
                //printf("momentum of pion + : %f %f %f \n",momP[0],momP[1],momP[2]);
                //printf("momentum of pion - : %f %f %f \n",momN[0],momN[1],momN[2]);
                //cout<<""<<endl;
                //cout<<"direction SV2: "<<direction_SV2[0]<<endl;
                //cout<<"positionSV2: "<<position_SV2[0]<<endl;
                //cout<<"position SV2: "<<position_SV2[0]<<endl;

            }
        }
        //if(V0_is_used==1){continue;}  */

        //-------------------------------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------------------------------
         //cout<<"j"<<endl;

        //search for V0 coming from anti-proton and K+
        if (fabs(sigma_proton_TPC[1]) < 2.5 && fabs(sigma_K_plus_TPC < 2.5))
        {
            //  cout<<"V0 from anti-proton and K+"<<endl;
            //energy_K_plus      = sqrt(mass_K * mass_K + (momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_K_plus      = sqrt(mass_K * mass_K + calc_momentum_squared(momP) );
            energy_anti_proton = sqrt(mass_proton * mass_proton + calc_momentum_squared(momN)) ;

            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_K_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_anti_proton);

            tlv_anti_p_and_K_plus = *tlv_pos + *tlv_neg;

            invariantmass = tlv_anti_p_and_K_plus.M();


            counters[5]++;

            //assume position of V0 is also position of potential S-vertex
            //loop over all other tracks in order to find K+ that comes from this position
            vector<int> save_track_ids;
            vector<int> save_track_nums;

            for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
            {

                AS_Track = AS_Event->getTrack(i_track_A);
                Int_t Trackid  = AS_Track->gettrackid();
                //cout<<Trackid<<endl;

                if(Trackid==trackidP){continue;}
                if(Trackid==trackidN){continue;}

                //if( check_if_int_is_in_vector(Trackid,used_track_ids_of_pions) == 1 ){continue;}

                double sigma = AS_Track -> getnsigma_K_TPC();

                // PID for Kaon
                if(fabs(sigma) > 2.5) {continue;}
                //initial parameters
                Float_t path_closest_to_point = 0;
                Float_t dca_closest_to_point  = 0;
                Float_t path_initA = 0.0;
                Float_t path_initB = 30.0;

                //calculate dca from vertex to particle track
                FindDCAHelixPoint(pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                if(dca_closest_to_point > 0.5) {continue;}
                //if(radius<5) {continue;}

                //check dca of additional K+ to primary vertex
                Float_t dca_to_primary_vertex = -1;
                FindDCAHelixPoint(pos_primary_vertex,AS_Track,path_initA,path_initB,path_closest_to_point,dca_to_primary_vertex);
                if(dca_to_primary_vertex<3.){continue;}

                //cout<<"Found V0 of antiproton and K+ with another K+"<<endl;

                save_track_ids.push_back(Trackid);
                save_track_nums.push_back(i_track_A);



            }

            if(save_track_ids.size()==1)
            {
                counter_vertices_antip_K_plus_K_plus++;
                counters[2]++;
                histos_1D[9]->Fill(radius);
                histos_2D[3]->Fill(pos[0],pos[1]);

                AS_Track = AS_Event->getTrack(save_track_nums[0]);
                Int_t Trackid  = AS_Track->gettrackid();
                TLorentzVector tlv = AS_Track->get_TLV_part();

                TVector3 dir;
                dir.SetXYZ(tlv_anti_p_and_K_plus[0]+tlv[0] , tlv_anti_p_and_K_plus[1]+tlv[1], tlv_anti_p_and_K_plus[2]+tlv[2]);

                TVector3 unit_dir;
                unit_dir = dir.Unit();


                double dot_product = unit_dir.Dot(unit_prim_to_V0);

                TVector3 null;
                null.SetXYZ(0.,0.,0.);


                if(radius>15)
                {
                    counter_vertices_antip_K_plus_K_plus_r_larger_5++;

                    if(dot_product<0.8){continue;}

                    counter_vertices_antip_K_plus_K_plus_r_larger_5_and_dot_product++;

                    //histo_counter->Fill(3.5);

                    double m_squared1 = calculate_m_squared_by_TOF(AS_Track);
                    double m_squared2 = calculate_m_squared_by_TOF(as_trackP);

                    if(m_squared1!=-1.) {mass_squared_kaons_and_background->Fill(m_squared1);}
                    if(m_squared2!=-1.) {mass_squared_kaons_and_background->Fill(m_squared2);}

                    //printf("mass squared Kaon 1: %f, Kaon 2: %f \n",m_squared1,m_squared2);
                    //printf("trackids: %d %d %d \n",trackidP,trackidN,Trackid);
                    if(m_squared1>0.2 && m_squared1<0.35 && m_squared2>0.2 && m_squared2<0.35)
                    {
                        histo_counter->Fill(7.5);
                        ////cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
                        ////cout<<"hi"<<endl;
                        ////cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

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
                    }

                }

            }

            ////cout<<"invariantmass: "<<invariantmass<<endl;
            save_track_ids.clear();

        }

        ///-----------------------------------------------------------------------------------
        //reference by looking at antineutron annhihilation
        //searching for V0 coming from pi+ and pi-  (K0-vertex)
        //then looking for 2 or 3 pions also from that vertex
        /*
        int num_of_pions =0; //total number of pions
        vector<int> tracks_of_V0;
        vector<int> all_track_ids_of_pions;
        vector<int> track_numbers_of_pions;
        Int_t trackid3=-1;

        if(fabs(sigma_pion_TPC[0]) < 2.5 && fabs(sigma_pion_TPC[1]) < 2.5)
        {
            tracks_of_V0.push_back(trackidP);
            tracks_of_V0.push_back(trackidN);

            all_track_ids_of_pions.push_back(trackidP);
            all_track_ids_of_pions.push_back(trackidN);

            counter++;
            num_of_pions = 2;

            for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
            {
                ////cout<<i_track_A<<endl;
                as_Track5 = AS_Event->getTrack(i_track_A);
                trackid3 = as_Track5->gettrackid();

                ////cout<<"trackids: "<<trackid3<<endl;

                if( check_if_int_is_in_vector(trackid3,tracks_of_V0) ){continue;}

                double sigma = as_Track5 -> getnsigma_pi_TPC();

                if(sigma>2.5){continue;}

                FindDCAHelixPoint(pos,as_Track5,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                if(dca_closest_to_point>0.5){continue;}

                num_of_pions++;

                all_track_ids_of_pions.push_back(trackid3);
                track_numbers_of_pions.push_back(i_track_A);
            }
            //printf("num of pions: %d \n",num_of_pions);

            if(num_of_pions == 3)
            {
                //if(radius<15){continue;}
                ////cout<<"filled"<<endl;
                histo_reference_vertex_radius_3_pionen->Fill(radius);
                histo_reference_x_and_y_3_pionen->Fill(pos[0],pos[1]);

            }

            if(num_of_pions == 4)
            {
                //if(radius<15){continue;}

                histo_reference_vertex_radius_4_pionen->Fill(radius);
                histo_reference_x_and_y_4_pionen->Fill(pos[0],pos[1]);

                //pion A is first V0 particle (positive)
                //pion B is second V0 particle (negative)
                //pion C is first extra pion
                //pion D is second extra pion

                //dca A is dca of particle A to V0


                AS_TrackA = as_trackP;
                AS_TrackB = as_trackN;
                AS_TrackC = AS_Event->getTrack(track_numbers_of_pions[0]);
                AS_TrackD = AS_Event->getTrack(track_numbers_of_pions[1]);

                //printf("track ids: %d %d %d %d \n",all_tracks_of_pions[0],all_tracks_of_pions[1],all_tracks_of_pions[2],all_tracks_of_pions[3]);

                if(!AS_TrackA || !AS_TrackB || !AS_TrackC || !AS_TrackD){continue;}

                float dcaA,dcaB,dcaC,dcaD;
                FindDCAHelixPoint(pos,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaA);
                FindDCAHelixPoint(pos,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaB);
                FindDCAHelixPoint(pos,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaC);
                FindDCAHelixPoint(pos,AS_TrackD,path_initA,path_initB,path_closest_to_point,dcaD);

                //printf("dcaA: %f, dcaB: %f, dcaC: %f, dcaD: %f \n",dcaA,dcaB,dcaC,dcaD);
                //pA is momentum of particle A
                float pA,pB,pC,pD;
                TLorentzVector tl_vec = AS_TrackA->get_TLV_part();
                pA = tl_vec.P();

                tl_vec = AS_TrackB->get_TLV_part();
                pB = tl_vec.P();

                tl_vec = AS_TrackC->get_TLV_part();
                pC= tl_vec.P();

                tl_vec = AS_TrackD->get_TLV_part();
                pD = tl_vec.P();

                //dcaAprim is dca of particle A to primary vertex
                float dcaAprim,dcaBprim,dcaCprim,dcaDprim;
                FindDCAHelixPoint(pos_primary_vertex,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaAprim);
                FindDCAHelixPoint(pos_primary_vertex,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaBprim);
                FindDCAHelixPoint(pos_primary_vertex,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaCprim);
                FindDCAHelixPoint(pos_primary_vertex,AS_TrackD,path_initA,path_initB,path_closest_to_point,dcaDprim);

                //("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");
                tpl->Fill(pos[0],pos[1],pos[2],dcaA,dcaB,dcaC,dcaD,pA,pB,pC,pD,dcaAprim,dcaBprim,dcaCprim,dcaDprim);
                ////cout<<"filled ntuple"<<endl;

            }




        }
        all_track_ids_of_pions.clear();
        tracks_of_V0.clear();


        */




        //-------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------
        //do it like in Master Thesis of Friederike Bock
         //cout<<"k"<<endl;
        if(dcaV0>0.5){continue;}
        //if(radius<30){continue;}

        AS_TrackA = as_trackP;
        AS_TrackB = as_trackN;

        float pA,pB;
        TLorentzVector tl_vec = AS_TrackA->get_TLV_part();
        pA = tl_vec.P();

        tl_vec = AS_TrackB->get_TLV_part();
        pB = tl_vec.P();

        //if(pA<0.5){continue;}
        //if(pB<0.5){continue;}
        //if(pA>2.){continue;}
        //if(pB>2.){continue;}

        


        UShort_t ntpcclsP = as_trackP->getNTPCcls();
        UShort_t nitsclsP = as_trackP->getNITScls();
        UShort_t statusP = as_trackP->getStatus();
        Float_t tpcchi2P = as_trackP->getTPCchi2();

        UShort_t ntpcclsN = as_trackN->getNTPCcls();
        UShort_t nitsclsN = as_trackN->getNITScls();
        UShort_t statusN = as_trackN->getStatus();
        Float_t tpcchi2N = as_trackN->getTPCchi2();

        //if(ntpcclsP<120 || ntpcclsN<120){continue;}

        /*//cout<<"tpccls: "<<ntpcclsP<<" "<<ntpcclsN<<endl;
        //cout<<"nitscls: "<<nitsclsP<<" "<<nitsclsN<<endl;
        //cout<<"status: "<<statusP<<" "<<statusN<<endl;
        //cout<<"tpcchi2: "<<tpcchi2P<<" "<<tpcchi2N<<endl;
        */
        
        //printf("track ids: %d %d %d %d \n",all_tracks_of_pions[0],all_tracks_of_pions[1],all_tracks_of_pions[2],all_tracks_of_pions[3]);
        //if(!AS_TrackA || !AS_TrackB ){continue;}

        float dcaA,dcaB;
        FindDCAHelixPoint(pos,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaA);

        if(tpcchi2P<150)
        {
            counter_path_trackA++;
            ////cout<<"pathA: "<<path_closest_to_point<<endl;
            if(path_closest_to_point<0){counter_path_trackA_negativ++;}
        }

        FindDCAHelixPoint(pos,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaB);

        if(tpcchi2N<150)
        {
            counter_path_trackB++;
            ////cout<<"pathB: "<<path_closest_to_point<<endl;
            ////cout<<""<<endl;
            if(path_closest_to_point<0){counter_path_trackB_negativ++;}
        }

        double sigma_e_P  = as_trackP->getnsigma_e_TPC();
        double sigma_e_N  = as_trackN->getnsigma_e_TPC();
        double sigma_pi_P = as_trackP-> getnsigma_pi_TPC();
        double sigma_pi_N = as_trackN-> getnsigma_pi_TPC();
        double sigma_K_P  = as_trackP-> getnsigma_K_TPC();
        double sigma_K_N  = as_trackN-> getnsigma_K_TPC();
        double sigma_p_P  = as_trackP-> getnsigma_p_TPC();
        double sigma_p_N  = as_trackN-> getnsigma_p_TPC();

        double m_squared;

        //do PID
        //if V0 particle is gamma or K0 or Lambda or Anti-Lambda then continue

        Float_t TPCdEdx   = as_trackP->getTPCdEdx();
        Float_t tofsignal = as_trackP->getTOFsignal();
        Float_t dca       = as_trackP->getdca();
        Float_t tracklength = as_trackP->getTrack_length();
        int charge;
        //invariant mass of e+
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

            histo_invariantmass_electron_plus->Fill(m_squared);

            //if(m_squared<0){printf("mass squared: %f \n", m_squared);}
            //printf("mass squared: %f \n", m_squared);
        }
         //cout<<"l"<<endl;
        //invariant mass of e-
        if(fabs(sigma_e_N)<2.5 )
        {
            energy_electron_minus  = sqrt(mass_electron*mass_electron+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            tlv_neg -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_electron_minus);
            invariantmass = tlv_neg ->M();
            histo_invariantmass_electron_minus->Fill(invariantmass);
        }

        //gamma  (two electrons)
        if(fabs(sigma_e_P)<2.5 && fabs(sigma_e_N)<2.5) //fabs(as_trackP->getnsigma_e_TOF())<2.5 && fabs(as_trackN->getnsigma_e_TOF())<2.5
        {
            
            energy_electron_plus  = sqrt(mass_electron*mass_electron+(momP[0]*momP[0]+momP[1]*momP[1]+momP[2]*momP[2]));
            energy_electron_minus = sqrt(mass_electron*mass_electron+(momN[0]*momN[0]+momN[1]*momN[1]+momN[2]*momN[2]));
            tlv_pos -> SetPxPyPzE(momP[0],momP[1],momP[2],energy_electron_plus);
            tlv_neg -> SetPxPyPzE(momN[0],momN[1],momN[2],energy_electron_minus);
            *tlv_gamma= *tlv_pos + *tlv_neg;
            invariantmass = tlv_gamma->M();

            //printf("invariante Masse gamma: %f \n",invariantmass);
            histo_invariantmass_gamma->Fill(invariantmass);

            if( (m_squared<0.4 && tofsignal<99990) || tofsignal>99990)
            {
                if(invariantmass<0+0.05){continue;}
            }
        }
         //cout<<"m"<<endl;
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


        vector<int> track_ids_pions;
        vector<int> track_numbers_pions;
        vector<int> track_ids_V0;


        

        track_ids_V0.push_back(trackidP);
        track_ids_V0.push_back(trackidN);
         //cout<<"o"<<endl;
        for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
        {

            AS_Track = AS_Event->getTrack(i_track_A);
            int trackid2 = AS_Track->gettrackid();

            if( check_if_int_is_in_vector(trackid2,track_ids_V0) == 1 ) {continue;};

            double sigma = AS_Track -> getnsigma_pi_TPC();

            // Do some PID here for pi+ and pi-
            if(fabs(sigma)>2.5){continue;}

            //calculate dca from vertex to particle track
            FindDCAHelixPoint(pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            if(dca_closest_to_point>1){continue;}

            track_ids_pions.push_back(trackid2);
            track_numbers_pions.push_back(i_track_A);

            //track_ids_all[V0counter].push_back(trackidP);
            //track_ids_all[V0counter].push_back(trackidN);
            //track_ids_all[V0counter].push_back(trackid2);



        }
         //cout<<"p"<<endl;
        if(track_ids_pions.size() == 0)
        {

            //pion A is first V0 particle (positive)
            //pion B is second V0 particle (negative)
            //pion C is first extra pion

            //dca A is dca of particle A to V0
            //if( check_if_int_is_in_vector(trackiP,brute_force) == 1 ) {continue;};
           // if( check_if_int_is_in_vector(trackiN,brute_force) == 1 ) {continue;};


            AS_TrackA = as_trackP;
            AS_TrackB = as_trackN;

            //printf("track ids: %d %d %d %d \n",all_tracks_of_pions[0],all_tracks_of_pions[1],all_tracks_of_pions[2],all_tracks_of_pions[3]);
            if(!AS_TrackA || !AS_TrackB ){continue;}

            float dcaA,dcaB;
            FindDCAHelixPoint(pos,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaA);
            FindDCAHelixPoint(pos,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaB);

            //pA is momentum of particle A
            float pA,pB;
            TLorentzVector tl_vec = AS_TrackA->get_TLV_part();
            pA = tl_vec.P();

            tl_vec = AS_TrackB->get_TLV_part();
            pB = tl_vec.P();

            //dcaAprim is dca of particle A to primary vertex
            float dcaAprim,dcaBprim;
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaAprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaBprim);

            //("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");
            // XHERE
            //tpl->Fill(pos[0],pos[1],pos[2],dcaA,dcaB,-999.,-999.,pA,pB,-999.,-999.,dcaAprim,dcaBprim,-999.,-999.);

            ////cout<<"filled ntuple"<<endl;
            //for(int i = 0;i < track_ids_pions.size();i++)
        }
         //cout<<"q"<<endl;

        if(track_ids_pions.size() == 1)
        {

            
            //pion A is first V0 particle (positive)
            //pion B is second V0 particle (negative)
            //pion C is first extra pion

            //dca A is dca of particle A to V0

            int a = check_if_int_is_in_vector(trackidP,brute_force);
            int b = check_if_int_is_in_vector(trackidN,brute_force);
            int c = check_if_int_is_in_vector(track_ids_pions[0],brute_force);

            int check = a+b+c;

            if(check>=1){continue;}


            AS_TrackA = as_trackP;
            AS_TrackB = as_trackN;
            AS_TrackC = AS_Event->getTrack(track_numbers_pions[0]);

            //printf("track ids: %d %d %d %d \n",all_tracks_of_pions[0],all_tracks_of_pions[1],all_tracks_of_pions[2],all_tracks_of_pions[3]);

            if(!AS_TrackA || !AS_TrackB || !AS_TrackC){continue;}

            float dcaA,dcaB,dcaC;
            FindDCAHelixPoint(pos,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaA);
            FindDCAHelixPoint(pos,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaB);
            FindDCAHelixPoint(pos,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaC);

            //printf("dcaA: %f, dcaB: %f, dcaC: %f, dcaD: %f \n",dcaA,dcaB,dcaC,dcaD);
            //pA is momentum of particle A
            float pA,pB,pC;
            TLorentzVector tl_vec = AS_TrackA->get_TLV_part();
            pA = tl_vec.P();

            tl_vec = AS_TrackB->get_TLV_part();
            pB = tl_vec.P();

            tl_vec = AS_TrackC->get_TLV_part();
            pC= tl_vec.P();

            //dcaAprim is dca of particle A to primary vertex
            float dcaAprim,dcaBprim,dcaCprim;
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaAprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaBprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaCprim);

            //("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");
            // XHERE
            //tpl->Fill(pos[0],pos[1],pos[2],dcaA,dcaB,dcaC,-999.,pA,pB,pC,-999.,dcaAprim,dcaBprim,dcaCprim,-999.);

            brute_force.push_back(track_ids_pions[0]);

            
            ////cout<<"filled ntuple"<<endl;
            //for(int i = 0;i < track_ids_pions.size();i++)
        }

         //cout<<"r"<<endl;
        if(track_ids_pions.size() == 2)
        {

            //pion A is first V0 particle (positive)
            //pion B is second V0 particle (negative)
            //pion C is first extra pion
            //pion D is second extra pion

            //dca A is dca of particle A to V0

            //check all tracks if one of them is already used:
            int a = check_if_int_is_in_vector(trackidP,brute_force);
            int b = check_if_int_is_in_vector(trackidN,brute_force);
            int c = check_if_int_is_in_vector(track_ids_pions[0],brute_force);
            int d = check_if_int_is_in_vector(track_ids_pions[1],brute_force);

            int check = a+b+c+d;

            if(check>=1){continue;}

            AS_TrackA = as_trackP;
            AS_TrackB = as_trackN;
            AS_TrackC = AS_Event->getTrack(track_numbers_pions[0]);
            AS_TrackD = AS_Event->getTrack(track_numbers_pions[1]);

            //printf("track ids: %d %d %d %d \n",all_tracks_of_pions[0],all_tracks_of_pions[1],all_tracks_of_pions[2],all_tracks_of_pions[3]);

            if(!AS_TrackA || !AS_TrackB || !AS_TrackC || !AS_TrackD){continue;}

            float dcaA,dcaB,dcaC,dcaD;
            FindDCAHelixPoint(pos,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaA);
            ////cout<<"pathA: "<<path_closest_to_point<<endl;
            //counter_path_trackA++;
           // if(path_closest_to_point<0){counter_path_trackA_negativ++;}

            FindDCAHelixPoint(pos,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaB);
            ////cout<<"pathB: "<<path_closest_to_point<<endl;

            FindDCAHelixPoint(pos,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaC);
            ////cout<<"pathC: "<<path_closest_to_point<<endl;
            counter_path_trackC++;
            if(path_closest_to_point<0){counter_path_trackC_negativ++;}

            FindDCAHelixPoint(pos,AS_TrackD,path_initA,path_initB,path_closest_to_point,dcaD);
            ////cout<<"pathD: "<<path_closest_to_point<<endl;

            //printf("dcaA: %f, dcaB: %f, dcaC: %f, dcaD: %f \n",dcaA,dcaB,dcaC,dcaD);
            //pA is momentum of particle A
            float pA,pB,pC,pD;
            TLorentzVector tl_vec = AS_TrackA->get_TLV_part();
            pA = tl_vec.P();

            tl_vec = AS_TrackB->get_TLV_part();
            pB = tl_vec.P();

            tl_vec = AS_TrackC->get_TLV_part();
            pC= tl_vec.P();

            tl_vec = AS_TrackD->get_TLV_part();
            pD = tl_vec.P();

            //dcaAprim is dca of particle A to primary vertex
            float dcaAprim,dcaBprim,dcaCprim,dcaDprim;
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackA,path_initA,path_initB,path_closest_to_point,dcaAprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackB,path_initA,path_initB,path_closest_to_point,dcaBprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackC,path_initA,path_initB,path_closest_to_point,dcaCprim);
            FindDCAHelixPoint(pos_primary_vertex,AS_TrackD,path_initA,path_initB,path_closest_to_point,dcaDprim);

            //("ntuple","ntuple","x:y:z:dcaA:dcaB:dcaC:dcaD:pA:pB:pC:pD:dcaAprim:dcaBprim:dcaCprim:dcaDprim");
            // XHERE
            //tpl->Fill(pos[0],pos[1],pos[2],dcaA,dcaB,dcaC,dcaD,pA,pB,pC,pD,dcaAprim,dcaBprim,dcaCprim,dcaDprim);

            brute_force.push_back(track_ids_pions[0]);
            brute_force.push_back(track_ids_pions[1]);
            ////cout<<"filled ntuple"<<endl;
            //for(int i = 0;i < track_ids_pions.size();i++)
        }
         //cout<<"s"<<endl;

        brute_force.push_back(trackidP);
        brute_force.push_back(trackidN);




        






    }     //end of V0 loop

    cout<<"end of V0 loop"<<endl;
    trackids_bits_shared.clear();
    brute_force.clear();
  

    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------------
    



    sort(all_used_positive_track_ids_for_V0s.begin(),all_used_positive_track_ids_for_V0s.end());
    sort(all_used_negative_track_ids_for_V0s.begin(),all_used_negative_track_ids_for_V0s.end());

    //print_int_vector(all_used_positive_track_ids_for_V0s);
    ////cout<<""<<endl;
    //print_int_vector(all_used_negative_track_ids_for_V0s);

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

    TVector3 S_vertex_pos;

    float path_closest_to_point = 0;
    float dca_closest_to_point  = 0;
    float path_initA = 0.0;
    float path_initB = 30.0;




    for(Int_t vector_loop_SV3 = 0; vector_loop_SV3 < (Int_t)vec_position_SV3.size(); vector_loop_SV3++)
    {
        for(Int_t vector_loop_SV2 = 0; vector_loop_SV2 < (Int_t)vec_position_SV2.size(); vector_loop_SV2++)
        {
            // printf("vector loop 1: %d, vector loop2: %d \n",vector_loop_SV3,vector_loop_SV2 ) ;
            if(vec_position_SV3.size() > 0 && vec_position_SV2.size() > 0 && vec_direction_SV3.size() > 0 &&  vec_direction_SV2.size() > 0)
            {
                if(calculateMinimumDistance(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]) > 0.5){ continue;}

                S_vertex_pos = calcVertexAnalytical(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]);
                counters[3]++;

                vector<int> trackidsV0s;
                trackidsV0s.push_back(vec_SV3_track_ids[vector_loop_SV3]);
                trackidsV0s.push_back(vec_SV3_track_ids[vector_loop_SV3+1]);
                trackidsV0s.push_back(vec_SV2_track_ids[vector_loop_SV2]);
                trackidsV0s.push_back(vec_SV2_track_ids[vector_loop_SV2+1]);
                //build lorentz vector of SV2 and SV3
                //SV2
                double px,py,pz,E;
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
                S_mass = tlv_SV1->M();

                //check if SV1 vector is parallel to vertex from primary vertex to SV1 (s-vertex)

                // vector from primary vertex to SV1 (S-vertex)
                vec_primary_vertex_to_SV1.SetXYZ(S_vertex_pos[0]-EventVertexX ,S_vertex_pos[1]-EventVertexY , S_vertex_pos[2]-EventVertexZ);
                unit_vec_primary_vertex_to_SV1 = vec_primary_vertex_to_SV1.Unit();



                radiusS = vec_primary_vertex_to_SV1.Mag();

                //------------------------------------------------------------------------------------
                //search for two pions coming from S-vertex------------------------------------------------------
                //------------------------------------------------------------------------------------

                //set to zero for each vertex
                Int_t counter_pions_close_to_S_vertex = 0;

                //store tracks of all pions that come from S-vertex
                vector<int> tracknumbers;
                vector<int> trackids;

                //------------------------------------------------------------
                //for each S-vertex loop over all tracks of event
                for(Int_t i_track_A = 0; i_track_A < NumTracks; i_track_A++)
                {
                    AS_Track = AS_Event->getTrack(i_track_A);
                    int trackid2 = AS_Track->gettrackid();

                    if( check_if_int_is_in_vector(trackid2,trackidsV0s) ) {continue;};


                    double sigma = AS_Track -> getnsigma_pi_TPC();

                    // Do some PID here for pi+ and pi-
                    if(fabs(sigma)>2.5){continue;}

                    //calculate dca from vertex to particle track
                    FindDCAHelixPoint(S_vertex_pos,AS_Track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                    //printf("i_track_A: %d, path_closest_to_point: %4.3f, dca_closest_to_point: %4.3f \n",i_track_A,path_closest_to_point,dca_closest_to_point);

                    // if dca is good then calculate Lorentzvectors at vertex position, add them to the S and calculate invariant mass

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

                ////cout<<"number pions close to S-vertex: "<<counter_pions_close_to_S_vertex<<endl;

                //if(tracks.size()<2){continue;}
                ////cout<<tracks[0]<<"  "<<tracks[1]<<endl;
                if(trackids.size()==1)
                {
                    counters[4]++;
                }

                //check if exactly 2 pions come from S-vertex
                if(trackids.size()==2)
                {

                    //check all tracks if one of them is already used:
                    int a = check_if_int_is_in_vector(trackids[0],brute_force);
                    int b = check_if_int_is_in_vector(trackids[1],brute_force);
                    int c = check_if_int_is_in_vector(vec_SV2_track_ids[2*vector_loop_SV2],brute_force);
                    int d = check_if_int_is_in_vector(vec_SV2_track_ids[2*vector_loop_SV2+1],brute_force);
                    int e = check_if_int_is_in_vector(vec_SV3_track_ids[2*vector_loop_SV3],brute_force);
                    int f = check_if_int_is_in_vector(vec_SV3_track_ids[2*vector_loop_SV3+1],brute_force);

                    int check = a+b+c+d+e+f;

                    if(check>=1){continue;}

                    vector<int> alltrackids;
                    alltrackids.push_back(trackids[0]);   //pion1
                    alltrackids.push_back(trackids[1]);   //pion2
                    alltrackids.push_back(vec_SV2_track_ids[2*vector_loop_SV2]);   //tracks of SV2
                    alltrackids.push_back(vec_SV2_track_ids[2*vector_loop_SV2+1]);   //tracks of SV2
                    alltrackids.push_back(vec_SV3_track_ids[2*vector_loop_SV3]);   //tracks of SV3
                    alltrackids.push_back(vec_SV3_track_ids[2*vector_loop_SV3+1]);   //tracks of SV3

                    if ( check_if_value_is_doppelt_in_vector(alltrackids) ) {continue;}



                    Double_t r1[3];
                    Double_t r2[3];

                    //count S-vertices with two 2 pions coming from there for all events
                    counters[1]++;
                    histo_counter->Fill(1.5);

                    //get tracks of 2 pions
                    //if(!ASTrack1 || !ASTrack2){continue;}
                    ASTrack1 = AS_Event->getTrack(tracknumbers[0]);
                    ASTrack2 = AS_Event->getTrack(tracknumbers[1]);

                    double dca1 = ASTrack1->getdca();
                    double dca2 = ASTrack2->getdca();

                    //check if one is negative and one is positive
                    if (dca1 * dca2 >0) {continue;}


                    //check dca to primary vertex:
                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;
                    FindDCAHelixPoint(pos_primary_vertex,ASTrack1,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
                    if(dca_closest_to_point<0.5){continue;}

                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;
                    FindDCAHelixPoint(pos_primary_vertex,ASTrack2,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
                    if(dca_closest_to_point<0.5){continue;}

                    //do all for pion 1--------------------------------------------------------------
                    //---------------------------------------------------------------------------------
                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;

                    //calculate again path and dca
                    FindDCAHelixPoint(S_vertex_pos,ASTrack1,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

                    //find direction of pion 1
                    ASTrack1->Evaluate(path_closest_to_point,r1);
                    ASTrack1->Evaluate(path_closest_to_point+0.01,r2);

                    TVector3 mom_dir_pion1;
                    TVector3 unit_mom_dir_pion1;
                    TVector3 vec_momentum;

                    mom_dir_pion1.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);

                    //unit vector in direction of pion 1
                    unit_mom_dir_pion1 = mom_dir_pion1.Unit();

                    ////cout<<r1[0]<<endl;
                    ////cout<<r2[0]<<endl;

                    TLorentzVector tlv = ASTrack1->get_TLV_part();
                    double momentum = tlv.P();

                    vec_momentum.SetXYZ(unit_mom_dir_pion1[0]*momentum,unit_mom_dir_pion1[1]*momentum,unit_mom_dir_pion1[2]*momentum);

                    ////cout<<"momentum vec: "<<vec_momentum[0]<<endl;


                    TLorentzVector tlv_pion1;
                    double energy_pion1 = sqrt(mass_pion * mass_pion + momentum*momentum) ;

                    tlv_pion1.SetPxPyPzE(vec_momentum[0],vec_momentum[1],vec_momentum[2],energy_pion1);

                    *tlv_SV1+=tlv_pion1;


                    //-----------------------------------------------------------------------------------------------
                    //-----------------------------------------------------------------------------------------------
                    //--------------------------------------------------------------------------------
                    //do the same for pion 2
                    path_closest_to_point = 0;
                    dca_closest_to_point  = 0;
                    path_initA = 0.0;
                    path_initB = 30.0;

                    //calculate again path and dca
                    FindDCAHelixPoint(S_vertex_pos,ASTrack2,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);

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
                    if(dot_product<0.){continue;}

                   

                    //--------------------------------------------------------------------------------------------------------------

                    //calculate again S-mass----------------------------------------------
                    double S_mass_correct = tlv_SV1->M();
                    ////cout<<"falsche S-mass: "<<S_mass<<endl;
                    ////cout<<"korrekte S-mass: "<<S_mass_correct<<endl;

                    if( fabs(radiusS) < 200 )
                    {
                        ////cout<<S_vertex_pos[0]<<endl;
                        histos_1D[4]->Fill(radiusS);

                        histos_2D[2]->Fill(S_vertex_pos[0],S_vertex_pos[1]);

                        if(fabs(radiusS)>10)
                        {
                            histos_1D[5]->Fill(S_mass_correct);
                        }

                        if(fabs(radiusS)>20)
                        {
                            histos_1D[6]->Fill(S_mass_correct);
                        }
                    }

                    

                    TVector3 TV3_S1(tlv_SV1->Px(),tlv_SV1->Py(),tlv_SV1->Pz());

                    histo_counter->Fill(2.5);
                    
                    if(radiusS<5){continue;}

                    
                    counters[6]++;

                    //vertices that fulfill all cuts
                    histo_counter->Fill(3.5);


                    //after all cuts store all 6 tracks in brute force
                    brute_force.push_back(trackids[0]);   //pion1
                    brute_force.push_back(trackids[1]);   //pion2
                    brute_force.push_back(vec_SV2_track_ids[2*vector_loop_SV2]);   //tracks of SV2
                    brute_force.push_back(vec_SV2_track_ids[2*vector_loop_SV2+1]);   //tracks of SV2
                    brute_force.push_back(vec_SV3_track_ids[2*vector_loop_SV3]);   //tracks of SV3
                    brute_force.push_back(vec_SV3_track_ids[2*vector_loop_SV3+1]);   //tracks of SV3

                    //-------------------------
                    // Store the S information in the output tree

                    AS_DM_particle ->set_primVertex(pos_primary_vertex);
                    AS_DM_particle ->set_S1Vertex(S_vertex_pos);
                    AS_DM_particle ->set_S2Vertex(vec_position_SV2[vector_loop_SV2]);
                    AS_DM_particle ->set_S3Vertex(vec_position_SV3[vector_loop_SV3]);
                    AS_DM_particle ->set_DirSV1(TV3_S1);
                    AS_DM_particle ->set_DirSV2(vec_direction_SV2[vector_loop_SV2]);
                    AS_DM_particle ->set_DirSV3(vec_direction_SV3[vector_loop_SV3]);
                    AS_DM_particle ->setN_V0s(3);
                    AS_DM_particle ->clearTrackList();
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(ASTrack1,AS_DM_Track);
                    ////cout<<ASTrack1->gettrackid()<<endl;
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(ASTrack2,AS_DM_Track);
                    ////cout<<ASTrack2->gettrackid()<<endl;
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(vec_SV2_tracks[2*vector_loop_SV2],AS_DM_Track);
                    ////cout<<vec_SV2_tracks[2*vector_loop_SV2]->gettrackid()<<endl;
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(vec_SV2_tracks[2*vector_loop_SV2+1],AS_DM_Track);
                    ////cout<<vec_SV2_tracks[2*vector_loop_SV2+1]->gettrackid()<<endl;
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(vec_SV3_tracks[2*vector_loop_SV3],AS_DM_Track);
                    ////cout<<vec_SV3_tracks[2*vector_loop_SV3]->gettrackid()<<endl;
                    AS_DM_Track = AS_DM_particle ->createTrack();
                    copy_track_params(vec_SV3_tracks[2*vector_loop_SV3+1],AS_DM_Track);
                    ////cout<<vec_SV3_tracks[2*vector_loop_SV3+1]->gettrackid()<<endl;
                    ////cout<<""<<endl;

                    Tree_AS_DM_particle ->Fill();
                    cout<<"filled Tree"<<endl;
                    //-------------------------

                }
                trackids.clear();
                tracknumbers.clear();

                //for reaction anti S
                //check if exactly 1 pion comes from S-vertex



                //------------------------------------------------------------



            }
        }


    }  //end of vector loop
    cout<<"end of vector loop"<<endl;

    /*for (int i = 0;i<3;i++)
     {
     //cout<<i<<endl;
     if(cat_direction_SV2[0][0].size()==0 ){continue;}
     //printf("2 direction SV2: %f \n",cat_direction_SV2[0][0][i][0]);
     }
     */

    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------


    N_Tracks = NumTracks;
    //Tree_AS_DM_particle->Write();

    return 1;


}

void Ali_Dark_Matter_Read::Save()
{
    /*
    outputfile->cd();
    Tree_AS_DM_particle->Write();
    histo_counter->Write();
    */
    
    TCanvas* can6 = new TCanvas;
    TCanvas* can7 = new TCanvas;
    TCanvas* can8 = new TCanvas;
    TCanvas* can9 = new TCanvas;
    TCanvas* can10 = new TCanvas;
    TCanvas* can11 = new TCanvas;
    TCanvas* can12 = new TCanvas;


    //printf("number of entries: %d, number of events over that we looped: %d \n", numentries, counters[0]);
    cout<<""<<endl;
    cout<<"considering reaction: anti-S + n -> K0s + anti-lambda + pi+ + pi- :"<<endl;
    printf("number of S-vertices without 2 pions: %d \n", counters[3]);
    printf("number of exactly 1 pion close to S-vertex: %d \n", counters[4]);
    printf("number of exactly 2 pions close to S-vertex: %d \n", counters[1]);
    printf("number of S vertices that fullfill all cuts: %d \n", counters[6]);
    
    cout<<""<<endl;
    cout<<"consindering reaction: anti-S + p -> anti-p + K+ + K+ " <<endl;
    printf("number of V0s of antiproton and K+: %d \n",counters[5]);
    printf("number of V0s of antiproton and K+ with another K+: %d \n", counter_vertices_antip_K_plus_K_plus);
    printf("number of V0s of antiproton and K+ with another K+ r>15: %d \n", counter_vertices_antip_K_plus_K_plus_r_larger_5);
    printf("number of V0s of antiproton and K+ with another K+ r>15 and dot product>0.8: %d \n", counter_vertices_antip_K_plus_K_plus_r_larger_5_and_dot_product);


    cout<<"counter: "<<counter<<endl;

    printf("counterA: %d \n",counter_path_trackA);
    printf("counterA neg: %d \n",counter_path_trackA_negativ);
    printf("counterB : %d \n",counter_path_trackB);
    printf("counterB neg: %d \n",counter_path_trackB_negativ);

    cout<<"counter skipped tracks: "<<counter_skipped_track<<endl;

   
    //printf("prozent A: %d \n",counter_path_trackA_negativ/counter_path_trackA);
    //printf("counterC: %d \n",counter_path_trackC);
    //printf("counterC neg: %d \n",counter_path_trackC_negativ);
    //printf("prozent c: %d \n",counter_path_trackC_negativ/counter_path_trackC);

    /*can->cd();
    mass_squared_vs_charge_dot_momentum->Draw("colz");

    can2->cd();
    dEdx_vs_charge_dot_momentum->Draw("colz");

    can3->cd();
    histo_invariantmass_gamma->Draw();

    can4->cd();
    histo_invariantmass_electron_plus->Draw();

    can5->cd();
    histo_invariantmass_electron_minus->Draw();
    */

    can6->cd();
    histo_invariantmass_xi_minus_baryon->Draw();

    can7->cd();
    histo_invariantmass_xi_plus_baryon->Draw();

    can8->cd();
    histo_invariantmass_omega_minus_baryon->Draw();

    can9->cd();
    histo_invariantmass_omega_plus_baryon->Draw();

    can10->cd();
    histo_invariantmass_anti_lambda->Draw();

    can11->cd();
    histo_invariantmass_lambda->Draw();

    can12->cd();
    histo_invariantmass_K0->Draw();

    can11->SaveAs("invariant_mass_lambda.png");
    can12->SaveAs("invariant_mass_K0.png");
    /*
    vector<TCanvas*> vec_can;

    for(int i=0;i<16;i++)
    {
        TCanvas* can = new TCanvas();
        can->cd();
        vec_histo_omega_minus[i]->Draw();

    }

    for(int i=0;i<16;i++)
    {
        TCanvas* can = new TCanvas();
        can->cd();
        vec_histo_omega_plus[i]->Draw();

    }
    */

    //histo_reference_vertex_radius_3_pionen->Draw();

    /*
    can2->cd();
    histo_reference_vertex_radius_4_pionen->Draw();

    can3->cd();
    histo_reference_x_and_y_3_pionen->Draw("colz");

    can4->cd();
    histo_reference_x_and_y_4_pionen->Draw("colz");
    */

    //ofile->cd();
    //tpl->Write();


    //ofile2->cd();
    //mass_squared_vs_charge_dot_momentum->Write();
    //dEdx_vs_charge_dot_momentum->Write();

    printf("antilambdas: %d, lambdas: %d \n",counter_anti_lambdas,counter_lambdas);

    //int total_size = arr_distance_prim_sec.size()*arr_distance_daughter_particles.size()*arr_dca_daughter_prim.size();
    int total_size = 15*15*2;

    output_histos->cd();
    histo_invariantmass_lambda->Write();
    histo_invariantmass_anti_lambda->Write();
    histo_invariantmass_K0->Write();
    histo_S_vertex_radius->Write();
    histo_invariantmass_xi_minus_baryon->Write();
    histo_invariantmass_xi_plus_baryon->Write();
    mass_squared_kaons_and_background->Write();
    mass_squared_kaons->Write();
    histo_counter->Write();
    for(int i=0;i<total_size;i++)
    {
        vec_histo_invariantmass_Lambda[i]->Write();
        vec_histo_invariantmass_K0[i]->Write();
    }

    for(int i=0;i<45;i++)
    {
        vec_histo_radius_variation_invariantmass_K0[i]->Write();
        vec_histo_radius_variation_invariantmass_Lambda[i]->Write();
        cout<<"wrote"<<endl;
    }

    Tree_AS_DM_particle->Write();
    /*
    for(int i=0;i<16;i++)
    {
        vec_histo_omega_minus[i]->Write();
        vec_histo_omega_plus[i]->Write();
    }
    */

}
#endif
