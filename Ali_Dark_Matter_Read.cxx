
#include "Ali_Dark_Matter_Read.h"
#include "Ali_Dark_Matter_Read_LinkDef.h"
//#include "TRD_ST_Analyze_tracklets.cxx"
#include "TEveLine.h"
#include "TEvePointSet.h"
#include <TEveManager.h>
//#incclude "additionalfunc.h"

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


/*
void Draw_TPC_track(Ali_AS_Track* TRD_ST_TPC_Track, Int_t color, Double_t line_width, int i_track,TVector3 pos_prim)
{
    Double_t track_pos[3];
    Double_t radius_helix;

   // Ali_AS_Track* TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);
    vector<TEveLine*> vec_TPL3D_helix;
    vector<TEveLine*> vec_TPL3D_helix_inner;
    vector<TEveLine*> vec_TPL3D_helix_hull;

    vec_TPL3D_helix.resize(i_track+1);
    vec_TPL3D_helix_hull.resize(i_track+1);
    vec_TPL3D_helix_inner.resize(i_track+1);
    vec_TPL3D_helix[i_track] = new TEveLine();
    vec_TPL3D_helix_hull[i_track] = new TEveLine();
    vec_TPL3D_helix_inner[i_track] = new TEveLine();
    //for(Int_t i_param=0;i_param<6;i_param++)
    //    //cout<<TRD_ST_TPC_Track ->getHelix_param(i_param)<<" ";
    ////cout<<endl;
    for(Double_t track_path = 0.0; track_path < 1000; track_path += 1.0)
    {
        TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
        radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
        if(radius_helix > 370.0) break;
        if(fabs(track_pos[2]) > 360.0) break;
        if(radius_helix > 80.0)
        {
            vec_TPL3D_helix[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
            vec_TPL3D_helix_hull[i_track]   ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
        }
        if(radius_helix < 80.0)
        {
            vec_TPL3D_helix_inner[i_track] ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
        }


        //if(i_track == 0) printf("track_path: %4.3f, pos: {%4.2f, %4.2f, %4.2f} \n",track_path,track_pos[0],track_pos[1],track_pos[2]);
    }

    TString HistName = "track ";
    HistName += i_track;
    vec_TPL3D_helix[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix[i_track]    ->SetLineWidth(line_width);
    vec_TPL3D_helix[i_track]    ->SetMainColor(color);
    vec_TPL3D_helix[i_track]    ->SetMainAlpha(1.0);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_hull[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_hull[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_hull[i_track]    ->SetLineWidth(8);
    vec_TPL3D_helix_hull[i_track]    ->SetMainColor(kWhite);
    vec_TPL3D_helix_hull[i_track]    ->SetMainAlpha(0.3);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_inner[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_inner[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_inner[i_track]    ->SetLineWidth(2);
    vec_TPL3D_helix_inner[i_track]    ->SetMainColor(kGray);
    vec_TPL3D_helix_inner[i_track]    ->SetMainAlpha(0.8);
    //if(i_track == 3)
    
    {
        gEve->AddElement(vec_TPL3D_helix[i_track]);
        gEve->AddElement(vec_TPL3D_helix_hull[i_track]);
        gEve->AddElement(vec_TPL3D_helix_inner[i_track]);
    }

    TEvePointSet* TEvePoint_vertex   = NULL;
    TEveLine*     TEveLine_beam_axis = NULL;

    TEvePoint_vertex = new TEvePointSet();
    TEveLine_beam_axis = new TEveLine();

    TEvePoint_vertex ->SetPoint(0,pos_prim[0],pos_prim[1],pos_prim[2]);
    TEvePoint_vertex-> SetMarkerSize(3);
    TEvePoint_vertex->  SetMarkerStyle(20);
    TEvePoint_vertex-> SetMarkerColor(kGray+1);
    gEve->AddElement(TEvePoint_vertex);

    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-30.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,30.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainColor(kGray+2);
    gEve->AddElement(TEveLine_beam_axis);

    //printf("%s Ali_TRD_ST_Analyze::Draw_TPC_track, i_track: %d %s \n",KRED,i_track,KNRM);
}

*/





bool check_if_int_is_in_vector(int a, vector<int> vec)
{
    for(int i=0;i<(Int_t)vec.size();i++)
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

void print_int_vector(vector<int> vec)
{
    for(int i=0;i<(Int_t)vec.size();i++)
    {
        //cout<<"Vektor  "<<i<<": "<<vec[i]<<endl;

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


vector<int> inttobinaryarray(int x)
{
    double quotient = x;
    double remainder;
    vector<int> res;


    while(quotient!=0)
    {
        quotient = x/2;
        remainder = x%2;
        x = quotient;
        ////cout<<"remainder: "<<remainder<<endl;
        res.insert(res.begin(),remainder);
    }
    while(res.size()<8)
    {
        res.insert(res.begin(),0);

    }
    return res;
}

int getnumones(vector<int> vec)
{
    int counter = 0;
    for(int i=0;i<vec.size();i++)
    {
        if(vec[i]==1){counter++;}
    }
    return counter;

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
        if(fabs(tofsignal)<1e-10){return -1;}
        if(tofsignal<0){return -1;}

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

bool overlap_PID(Ali_AS_Track* track, TString particle)
{
    double sigma_pi = fabs( track -> getnsigma_pi_TPC() );
    double sigma_K = fabs ( track -> getnsigma_K_TPC() );
    double sigma_p = fabs ( track -> getnsigma_p_TPC() );
    double sigma_e = fabs ( track -> getnsigma_e_TPC() );

    double m_squared = calculate_m_squared_by_TOF(track);

    if(particle=="K")
    {
        if(sigma_pi<2.5 || sigma_p<2.5 || sigma_e<2.5)
        {
            if(m_squared>0.2 && m_squared < 0.35) return true;
            else return false;
        }
        return true;
    }

    if(particle=="p")
    {
        if(sigma_pi<2.5 || sigma_K<2.5 || sigma_e<2.5)
        {
            if(m_squared>0.6 && m_squared < 1.2) return true;
            else return false;
        }
        return true;

    }


}

bool isoverlap(Ali_AS_Track* track, TString particle)
{
    double sigma_pi = fabs( track -> getnsigma_pi_TPC() );
    double sigma_K = fabs ( track -> getnsigma_K_TPC() );
    double sigma_p = fabs ( track -> getnsigma_p_TPC() );
    double sigma_e = fabs ( track -> getnsigma_e_TPC() );

    if(particle=="K")
    {
        if(sigma_pi<2.5 || sigma_p<2.5 || sigma_e<2.5)
        {
            return true;
        }
        return false;
    }

    if(particle=="p")
    {
        if(sigma_pi<2.5 || sigma_K<2.5 || sigma_e<2.5)
        {
            return true;
        }
        return false;

    }

}

bool m2_cut_if_available(Ali_AS_Track* track, TString particle)
{
    double m = calculate_m_squared_by_TOF(track);
    if(m == -1) return 1;

    if(particle=="pi")
    {
        if( m < -0.2 || m > 0.2) return 0;
    }

    if(particle=="K")
    {
        if( m < 0.2 || m > 0.35) return 0;
    }

    if(particle=="p")
    {
        if( m < 0.6 || m > 1.2) return 0;
    }

    return 1;

}

void printtv3(TVector3 v,TString n)
{
    cout<<n<<": "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;

}

void printtlv(TLorentzVector v,TString n)
{
    cout<<n<<": "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;

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

void print_TV3(TVector3 vec)
{
    cout<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<endl;
}

void print_dm_event(Ali_AS_DM_particle* dm_in)
{
    cout<<"prim: ";
    print_TV3(dm_in->get_primVertex());
    cout<<"S: ";
    print_TV3(dm_in->get_S1Vertex());
    cout<<"S2: ";
    print_TV3(dm_in->get_S2Vertex());
    cout<<"S3: ";
    print_TV3(dm_in->get_S3Vertex());
    cout<<"type: "<<
    dm_in->gettype()<<endl;

    int numtracks = dm_in->getNumTracks();
    cout<<"number of tracks: "<<numtracks<<endl;

    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;
    Float_t path_closest_to_point = 0;
    Float_t dca_closest_to_point  = 0;

    for(int i=0;i<numtracks;i++)
    {
        Ali_AS_Track* track = dm_in -> getTrack(i);
        int trackid = track->gettrackid();
        cout<<"trackid: "<<trackid<<endl;

        float dca;
        if(i==0 || i==1 || i==2){ FindDCAHelixPoint(dm_in->get_S1Vertex(),track,path_initA,path_initB,path_closest_to_point,dca);}
        if(i==3 || i==4){ FindDCAHelixPoint(dm_in->get_S2Vertex(),track,path_initA,path_initB,path_closest_to_point,dca);}
        cout<<"dca: "<<dca<<endl;
    }

    
    
}

TLorentzVector get_tlv(Ali_AS_Track* track, double mass, TVector3 vertex)
{
    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;

    TLorentzVector tlv = track->get_TLV_part();
    double momentum = tlv.P();
    double energy  = sqrt(mass*mass+momentum*momentum);

    Float_t path_closest_to_point,dca_closest_to_point;
    Double_t r1[3];
    Double_t r2[3];
    FindDCAHelixPoint(vertex,track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
    track->Evaluate(path_closest_to_point,r1);
    track->Evaluate(path_closest_to_point+0.01,r2);
    TVector3 dir;
    dir.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
    TVector3 unit_dir = dir.Unit();
    TVector3 momentum_vec = unit_dir * momentum;
    tlv.SetPxPyPzE(momentum_vec[0],momentum_vec[1],momentum_vec[2],energy);

    return tlv;


}

TLorentzVector get_tlv_by_V0(double massP, double massN, float* momP, float* momN )
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



double get_invariant_mass(Ali_AS_Track* trackK0P, Ali_AS_Track* trackK0N, double massP, double massN , TVector3 K0_vertex)
{
    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;

    TLorentzVector tlvP = trackK0P->get_TLV_part();
    TLorentzVector tlvN = trackK0N->get_TLV_part();

    double momentumP = tlvP.P();
    double momentumN = tlvN.P();

    double energy_pion_plus  = sqrt(massP*massP+momentumP*momentumP);
    double energy_pion_minus  = sqrt(massN*massN+momentumN*momentumN);

    Float_t path_closest_to_point,dca_closest_to_point;
    Double_t r1[3];
    Double_t r2[3];

    //track1
    FindDCAHelixPoint(K0_vertex,trackK0P,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
    trackK0P->Evaluate(path_closest_to_point,r1);
    trackK0P->Evaluate(path_closest_to_point+0.01,r2);
    TVector3 dir;
    dir.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
    TVector3 unit_dir = dir.Unit();
    TVector3 momentum_vec = unit_dir * momentumP;
    tlvP.SetPxPyPzE(momentum_vec[0],momentum_vec[1],momentum_vec[2],energy_pion_plus);

    //track2
    FindDCAHelixPoint(K0_vertex,trackK0N,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
    trackK0N->Evaluate(path_closest_to_point,r1);
    trackK0N->Evaluate(path_closest_to_point+0.01,r2);
    dir.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
    unit_dir = dir.Unit();
    momentum_vec = unit_dir * momentumN;
    tlvN.SetPxPyPzE(momentum_vec[0],momentum_vec[1],momentum_vec[2],energy_pion_minus);

    TLorentzVector tlvK0 = tlvP + tlvN;
    double invariantmass = tlvK0.M();

    return invariantmass;
}




bool checkifsimiliar(Ali_AS_Track* tracka,Ali_AS_Track* trackb)
{
    TLorentzVector tlva = tracka->get_TLV_part();
    TLorentzVector tlvb = trackb->get_TLV_part();

    if( fabs( 1-tlva[0]/tlvb[0]) < 0.05 && fabs( 1-tlva[1]/tlvb[1]) < 0.05 && fabs( 1-tlva[2]/tlvb[2]) < 0.05)
    {
        double  pos[3];
        float dca=-1;
        Float_t path_closest_to_point = 0;
        Float_t dca_closest_to_point  = 0;
        Float_t path_initA = 0.0;
        Float_t path_initB = 30.0;
        int counter=0;

        for(int i=0;i<100;i+=10)
        {
            tracka->Evaluate(i,pos);
            TVector3 position;
            position[0]=pos[0];
            position[1]=pos[1];
            position[2]=pos[2];
            FindDCAHelixPoint(position,trackb,path_initA,path_initB,path_closest_to_point,dca);
            if(dca<0.5)
            {
                counter++;
            }
            if(counter>2){return 1;}
        }

    }

    return 0;

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
    TString outfile_name = outfile_dir + in_list_name + "NUCLEV_out.root";
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


    //full stat new server alidata100
    TString pinputdir = "/misc/alidata100/alice_u/schlichtmann/Pb_Pb_3_channels_1/";


    //3 channels
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_3_channels_1/";


    //Pb-Pb ch3 and ch5
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_ch1_and_ch5_1/";



    //newest Pb-Pb
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_antichannels_3/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_S_Search_V2/";

    //newest p-Pb
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/p_Pb_S_search_V10_more_stat/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/p_Pb_S_search_V9_more_stat/";

    //nuclev
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/p_Pb_S_search_V3/";

    //Pb-Pb type5     
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_ch5_5/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_ch5_2/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_ch5_1/";



    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/p_Pb_S_search_V9_more_stat/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/p_Pb_NUCLEV_V4/";
    //TString pinputdir = "/misc/alidata121/alice_u/schlichtmann/Pb_Pb_V3/";

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
        ////cout<<x_vals[i]<<endl;
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
        ////cout<<y_vals[i]<<endl;
    }


    TString str_type_count[4]={"1","2","3","4"};
    TString str_particles[3]={"pion","kaon","proton"};

   
    dEdx_vs_charge_dot_momentum = new TH2D("dEdx_vs_charge_dot_momentum","dEdx_vs_charge_dot_momentum",
                                           199,x_vals,199,y_vals);

    dEdx_vs_charge_dot_momentum_S_cand = new TH2D("dEdx_vs_charge_dot_momentum_S_cand","dEdx_vs_charge_dot_momentum_S_cand",
                                                  199,x_vals,199,y_vals);

    for(int i=0;i<4;i++)
    {
        TString name= "dEdx_vs_charge_dot_momentum_S_cand_type";
        name+=str_type_count[i];
        TH2D* histo  = new TH2D(name.Data(),name.Data(),199,x_vals,199,y_vals);
        vec_dEdx.push_back(histo);

    }

    for(int i=0;i<3;i++)
    {
        TString name= "dEdx_vs_charge_dot_momentum_S_for_particle";
        name+=str_particles[i];
        TH2D* histo  = new TH2D(name.Data(),name.Data(),199,x_vals,199,y_vals);
        vec_dEdx_S.push_back(histo);

    }

    TString nom= "dEdx_type5_anticuts";
    dEdx_type5_anticuts = new TH2D(nom.Data(),nom.Data(),199,x_vals,199,y_vals);


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

    for(int i=2;i<6;i++)
    {
        TString name = "radius_ortho_vs_z_";
        name+=i;
        TH2D* radius_ortho_vs_z = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
        vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z);
    }

    for(int i=2;i<6;i++)
    {
        TString name = "radius_ortho_vs_z_";
        name+=i;
        name+=" dcaprim>2";
        TH2D* radius_ortho_vs_z = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
        vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z);
    }

    for(int i=2;i<6;i++)
    {
        TString name = "radius_ortho_vs_z_";
        name+=i;
        name+=" dcaprim>5";
        TH2D* radius_ortho_vs_z = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
        vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z);
    }

    TString name = "radius_ortho_vs_z_";
    name+=4;
    name+=" dcaprim>5 momentum cut";
    TH2D* radius_ortho_vs_z = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
    vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z);

    name = "radius_ortho_vs_z_";
    name+=4;
    name+=" dcaprim>7";
    TH2D* radius_ortho_vs_z2 = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
    vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z2);

    name = "radius_ortho_vs_z_";
    name+=4;
    name+=" dcaprim>7 momentum cut";
    TH2D* radius_ortho_vs_z3 = new TH2D(name.Data(),name.Data(),500,-200,200,500,0,200);
    vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z3);

    int binning_nuclev=300;
    int range_min_nuclev=0;
    int range_max_nuclev=200;

    TString names[] ={"nuclev_radius_4_or_less","nuclev_radius_4_or_less dcaprim>5","nuclev_radius_4_or_less dcaprim>5 and momentum",
    "nuclev_radius_4_or_less dcaprim>7","nuclev_radius_4_or_less dcaprim>7 and momentum","nuclev_radius_4_or_less dcaprim>5 and momentum and its","nuclev_radius_4_or_less dcaprim>5 and its"} ;
    TString names2[] ={"nuclev_radius_4_or_less_vs_z","nuclev_radius_4_or_less dcaprim>5_vs_z",
    "nuclev_radius_4_or_less dcaprim>5 and momentum_vs_Z","nuclev_radius_4_or_less dcaprim>7_vs_z","nuclev_radius_4_or_less dcaprim>7 and momentum_vs_z",
    "nuclev_radius_4_or_less dcaprim>5 and momentum and its_vs_z","nuclev_radius_4_or_less dcaprim>5 and its_vs_z"} ;

    for(int i=0;i<sizeof(names)/sizeof(names[0]);i++)
    {
        TH1D* histo_radius_nuclev_4_or_less = new TH1D(names[i].Data(),names[i].Data(),binning_nuclev,range_min_nuclev,range_max_nuclev);
        vec_histo_radius_nuclev.push_back(histo_radius_nuclev_4_or_less);

        TH2D* radius_ortho_vs_z3 = new TH2D(names2[i].Data(),names[i].Data(),500,-200,200,500,0,200);
        vec_radius_ortho_vs_z.push_back(radius_ortho_vs_z3);
        
    }


    for(int i=0;i<10;i++)
    {
        int range_min=-100+ i*20;
        int range_max=-80 + i*20;
        TString name ="x_y_slices_z";
        name+=to_string(range_min);
        name+=to_string(range_max);
        TH2D* histo_x_y_slices = new TH2D(name.Data(),name.Data(),500,-200,200,500,-200,200);
        vec_x_y_slices.push_back(histo_x_y_slices);

        name = "radius_slices_z";
        name+=to_string(range_min);
        name+=to_string(range_max);
        TH1D* histo_radius = new TH1D(name.Data(),name.Data(),binning_nuclev,range_min_nuclev,range_max_nuclev);
        vec_radius_slices.push_back(histo_radius);
    }

    for(int i=0;i<20;i++)
    {
        int range_min=0+ i*10;
        int range_max=10 + i*10;
        name = "z_for_slice_r_";
        name+=to_string(range_min);
        name+=to_string(range_max);
        TH1D* histo_radius = new TH1D(name.Data(),name.Data(),500,-200,200);
        vec_z_slices_in_radius.push_back(histo_radius);

    }

    TString nameshist[3] = {"deltax,deltay,deltaz"};

    for(int i=0;i<3;i++)
    {
        TH1D* hist = new TH1D(nameshist[i].Data(),nameshist[i].Data(),40,-20,20);
        histo_delta.push_back(hist);
    }

    TString str_type[6]={"type1","type2","type3","type4","type5","type51"};

    for(int i=0;i<6;i++)
    {
        TString n = "S_radius_from_origin_";
        n+=str_type[i];
        TH1D* histo = new TH1D(n.Data(),n.Data(),200,0,200);
        vec_S_radius_from_origin.push_back(histo);
    }

    vec_mass_squared_dEdx_selected.resize(3);

    TString particles[3]={"pions","kaons","protons"};
    for(int i=0;i<3;i++)
    {
        TString n = "m_squared_for_dEdx_selected_";
        n+=particles[i];
        TH1D* histo = new TH1D(n.Data(),n.Data(),100,-0.4,1.4);
        vec_mass_squared_dEdx_selected[i]=histo;
    }

    histos_m_squared_type5.resize(8);
    histos_m_squared_type51.resize(8);

    TString part[4]={"K+","pi-","pi+","antip"};
    for(int i=0;i<8;i++)
    {
        if(i<4)
        {
            TString n = "m_squared_ch5_for_";
            TString n2 = "m_squared_ch51_for_";
            n+=part[i];
            n2+=part[i];
            histos_m_squared_type5[i] = new TH1D(n.Data(),n.Data(),400,-0.4,1.4);
            histos_m_squared_type51[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
        }
        if(i>=4)
        {
            TString n = "m_squared_ch5_cutonK_for_";
            TString n2 = "m_squared_ch51_cutonK_for_";
            n+=part[i-4];
            n2+=part[i-4];
            histos_m_squared_type5[i] = new TH1D(n.Data(),n.Data(),400,-0.4,1.4);
            histos_m_squared_type51[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
        }
    }

    m_squared_type5_K_r.resize(34);
    m_squared_type5_pi_minus_r.resize(34);
    m_squared_type5_pi_plus_r.resize(34);
    m_squared_type5_p_r.resize(34);

    m_squared_type51_K_r.resize(34);
    m_squared_type51_pi_minus_r.resize(34);
    m_squared_type51_pi_plus_r.resize(34);
    m_squared_type51_p_r.resize(34);


    for(int i=0;i<34;i++)
    {
        if(i<17)
        {
            TString n1 = "m_squared_ch5_K_r_";
            TString n2 = "m_squared_ch5_pi_minus_r_";
            TString n3 = "m_squared_ch5_pi_plus_r_";
            TString n4 = "m_squared_ch5_p_r_";

            TString n5 = "m_squared_ch51_K_r_";
            TString n6 = "m_squared_ch51_pi_minus_r_";
            TString n7 = "m_squared_ch51_pi_plus_r_";
            TString n8 = "m_squared_ch51_p_r_";

            int j = 5+i*5;

            n1+=to_string(j);
            n2+=to_string(j);
            n3+=to_string(j);
            n4+=to_string(j);
            n5+=to_string(j);
            n6+=to_string(j);
            n7+=to_string(j);
            n8+=to_string(j);

            m_squared_type5_K_r[i] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
            m_squared_type5_pi_minus_r[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
            m_squared_type5_pi_plus_r[i] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);
            m_squared_type5_p_r[i] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);

            m_squared_type51_K_r[i] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
            m_squared_type51_pi_minus_r[i] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);
            m_squared_type51_pi_plus_r[i] = new TH1D(n7.Data(),n7.Data(),400,-0.4,1.4);
            m_squared_type51_p_r[i] = new TH1D(n8.Data(),n8.Data(),400,-0.4,1.4);
        }

        if(i>=17)
        {
            int j= 5+(i-17)*5;

            TString n1 = "m_squared_ch5_cut_on_K_K_r_";
            TString n2 = "m_squared_ch5_cut_on_K_pi_minus_r_";
            TString n3 = "m_squared_ch5_cut_on_K_pi_plus_r_";
            TString n4 = "m_squared_ch5_cut_on_K_p_r_";

            TString n5 = "m_squared_ch51_cut_on_K_K_r_";
            TString n6 = "m_squared_ch51_cut_on_K_pi_minus_r_";
            TString n7 = "m_squared_ch51_cut_on_K_pi_plus_r_";
            TString n8 = "m_squared_ch51_cut_on_K_p_r_";

            n1+=to_string(j);
            n2+=to_string(j);
            n3+=to_string(j);
            n4+=to_string(j);
            n5+=to_string(j);
            n6+=to_string(j);
            n7+=to_string(j);
            n8+=to_string(j);

            m_squared_type5_K_r[i] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
            m_squared_type5_pi_minus_r[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
            m_squared_type5_pi_plus_r[i] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);
            m_squared_type5_p_r[i] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);

            m_squared_type51_K_r[i] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
            m_squared_type51_pi_minus_r[i] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);
            m_squared_type51_pi_plus_r[i] = new TH1D(n7.Data(),n7.Data(),400,-0.4,1.4);
            m_squared_type51_p_r[i] = new TH1D(n8.Data(),n8.Data(),400,-0.4,1.4);
        }

    }

    vec_m_squared_type_3_K.resize(34);
    vec_m_squared_type_3_p.resize(34);
    vec_m_squared_type_3_add_pi.resize(34);
    vec_m_squared_type_3_K01.resize(34);
    vec_m_squared_type_3_K02.resize(34);

    vec_m_squared_type_31_K.resize(34);
    vec_m_squared_type_31_p.resize(34);
    vec_m_squared_type_31_add_pi.resize(34);
    vec_m_squared_type_31_K01.resize(34);
    vec_m_squared_type_31_K02.resize(34);


    for(int i=0;i<34;i++)
    {
        

        if(i<17)
        {
            int j = 5+i*5;
            TString n1 = "m_squared_ch3_K_r_";
            TString n2 = "m_squared_ch3_p_r_";
            TString n3 = "m_squared_ch3_add_pi_r_";
            TString n4 = "m_squared_ch3_K01_r_";
            TString n5 = "m_squared_ch3_K02_r_";

            TString n6 = "m_squared_ch31_K_r_";
            TString n7 = "m_squared_ch31_p_r_";
            TString n8 = "m_squared_ch31_add_pi_r_";
            TString n9 = "m_squared_ch31_K01_r_";
            TString n10 = "m_squared_ch31_K02_r_";


            n1+=to_string(j);
            n2+=to_string(j);
            n3+=to_string(j);
            n4+=to_string(j);
            n5+=to_string(j);
    
            n6+=to_string(j);
            n7+=to_string(j);
            n8+=to_string(j);
            n9+=to_string(j);
            n10+=to_string(j);
    
            vec_m_squared_type_3_K[i] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_p[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_add_pi[i] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_K01[i] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_K02[i] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
    
            vec_m_squared_type_31_K[i] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_p[i] = new TH1D(n7.Data(),n7.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_add_pi[i] = new TH1D(n8.Data(),n8.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_K01[i] = new TH1D(n9.Data(),n9.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_K02[i] = new TH1D(n10.Data(),n10.Data(),400,-0.4,1.4);
        }

        if(i>=17)
        {
            int j = 5+(i-17)*5;

            TString n1 = "m_squared_ch3_cut_on_K_and_p_K_r_";
            TString n2 = "m_squared_ch3_cut_on_K_and_p_p_r_";
            TString n3 = "m_squared_ch3_cut_on_K_and_p_add_pi_r_";
            TString n4 = "m_squared_ch3_cut_on_K_and_p_K01_r_";
            TString n5 = "m_squared_ch3_cut_on_K_and_p_K02_r_";

            TString n6 = "m_squared_ch31_cut_on_K_and_p_K_r_";
            TString n7 = "m_squared_ch31_cut_on_K_and_p_p_r_";
            TString n8 = "m_squared_ch31_cut_on_K_and_p_add_pi_r_";
            TString n9 = "m_squared_ch31_cut_on_K_and_p_K01_r_";
            TString n10 = "m_squared_ch31_cut_on_K_and_p_K02_r_";


            n1+=to_string(j);
            n2+=to_string(j);
            n3+=to_string(j);
            n4+=to_string(j);
            n5+=to_string(j);
    
            n6+=to_string(j);
            n7+=to_string(j);
            n8+=to_string(j);
            n9+=to_string(j);
            n10+=to_string(j);
    
            vec_m_squared_type_3_K[i] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_p[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_add_pi[i] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_K01[i] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);
            vec_m_squared_type_3_K02[i] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
    
            vec_m_squared_type_31_K[i] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_p[i] = new TH1D(n7.Data(),n7.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_add_pi[i] = new TH1D(n8.Data(),n8.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_K01[i] = new TH1D(n9.Data(),n9.Data(),400,-0.4,1.4);
            vec_m_squared_type_31_K02[i] = new TH1D(n10.Data(),n10.Data(),400,-0.4,1.4);
        }

    }


    vec_m_sq_p_type_3_for_dcas.resize(5);
    vec_m_sq_p_type_31_for_dcas.resize(5);
    TString n[5];
    TString n2[5];
    TString arr_dca[5] = {"0.5","1.0","1.5","2.0","3.0"};

    for(int i=0;i<5;i++)
    {
        n[i]="m_squared_p_type3_for_dca_";
        n2[i]="m_squared_p_type31_for_dca_";
        n[i] += arr_dca[i];
        n2[i] += arr_dca[i];
        vec_m_sq_p_type_3_for_dcas[i] =  new TH1D(n[i].Data(),n[i].Data(),400,-0.4,1.4);
        vec_m_sq_p_type_31_for_dcas[i] =  new TH1D(n2[i].Data(),n2[i].Data(),400,-0.4,1.4);
    }



    
    mixed_event_vec.resize(4);
    for(int i=0;i<4;i++)
    {
        mixed_event_vec[i].resize(10);
        for(int j=0;j<10;j++)
        {
            mixed_event_vec[i][j].reserve(5);
        }
    }

    TString arr_r[6]={"5","10","15","20","30","50"};
    for(int i=0;i<6;i++)
    {
        TString n ="histo_invmass_lambda_type5_r_";
        TString n51 ="histo_invmass_lambda_type51_r_";
        TString n2 ="histo_invmass_S_type5_r_";
        TString n3 ="histo_invmass_S_type51_r_";
        n+=arr_r[i];
        n2+=arr_r[i];
        n3+=arr_r[i];
        n51+=arr_r[i];
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n.Data(),n.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n51.Data(),n51.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_S_type5 = new TH1D(n2.Data(),n2.Data(),200,0.,10);
        TH1D* histo_invmass_S_type51 = new TH1D(n3.Data(),n3.Data(),200,0.,10);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
        vec_invmass_S_type5.push_back(histo_invmass_S_type5);
        vec_invmass_S_type51.push_back(histo_invmass_S_type51);

    }


    for(int i=0;i<6;i++)
    {
        TString n ="histo_invmass_lambda_type5_overlap_r_";
        TString n51 ="histo_invmass_lambda_type51_overlap_r_";
        n+=arr_r[i];
        n51+=arr_r[i];
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n.Data(),n.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n51.Data(),n51.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);

    }

    //TString arr_r2[5]={"10","20","25","30","50"};


    for(int i=0;i<17;i++)
    {
        TString n1="invmass_lambda_dca_cuts_type5_r_";
        TString n2="invmass_lambda_dca_cuts_type51_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_lambda_dca_cuts_and_Kaoncut_type5_r_";
        TString n2="invmass_lambda_dca_cuts_and_Kaoncut_type51_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
    }


    for(int i=0;i<17;i++)
    {
        TString n1="invmass_lambda_dca_cuts_type5_r_of_L_";
        TString n2="invmass_lambda_dca_cuts_type51_r_of_L_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_lambda_dca_cuts_and_Kaoncut_type5_r_of_L_";
        TString n2="invmass_lambda_dca_cuts_and_Kaoncut_type51_r_of_L_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_lambda_dca_cuts_and_Kaoncut_distL_S4_type5_r_";
        TString n2="invmass_lambda_dca_cuts_and_Kaoncut_distL_S4_type51_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_Lambda_type5 = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
        TH1D* histo_invmass_Lambda_type51 = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        vec_invmass_Lambda_type5.push_back(histo_invmass_Lambda_type5);
        vec_invmass_Lambda_type51.push_back(histo_invmass_Lambda_type51);
    }


    for(int i=0;i<17;i++)
    {
        TString n1="S_mass_ch5_normalband_r_";
        TString n2="S_mass_ch5_sideband_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_S_normal = new TH1D(n1.Data(),n1.Data(),500,0,10);
        TH1D* histo_S_side  = new TH1D(n2.Data(),n2.Data(),500,0,10);
        vec_S_mass_ch5_normalband.push_back(histo_S_normal);
        vec_S_mass_ch5_sideband.push_back(histo_S_side);
    }

    //-----------------------------------------------------------------------

    int binningK0 = 50*6;

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_type3_r_";
        TString n2="invmass_K0_type31_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }


    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_dca_cuts_type3_r_";
        TString n2="invmass_K0_dca_cuts_type31_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_dca_cuts_K_and_antipcut_type3_r_";
        TString n2="invmass_K0_dca_cuts_K_and_antipcut_type31_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_type3_r_of_K0";
        TString n2="invmass_K0_type31_r_of_K0";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_dca_cutscombi2_type3_r_";
        TString n2="invmass_K0_dca_cutscombi2_type31_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }

    for(int i=0;i<17;i++)
    {
        TString n1="invmass_K0_dca_cutscombi2_K_and_antipcut_type3_r_";
        TString n2="invmass_K0_dca_cutscombi2_K_and_antipcut_type31_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_invmass_K0_type3 = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);
        TH1D* histo_invmass_K0_type31 = new TH1D(n2.Data(),n2.Data(),binningK0,0.4,0.6);
        vec_invmass_K0_type3.push_back(histo_invmass_K0_type3);
        vec_invmass_K0_type31.push_back(histo_invmass_K0_type31);
    }

    vec_vec_invmass_K0_type3.resize(17*2);
    vec_vec_invmass_K0_type31.resize(17*2);

    vec_vec_invmass_L_type5.resize(17*2);
    vec_vec_invmass_L_type51.resize(17*2);

    for(int i=0;i<17*2;i++)
    {
        vec_vec_invmass_K0_type3[i].resize(6);
        vec_vec_invmass_K0_type31[i].resize(6);
        
    }

    for(int i=0;i<2*17;i++)
    {
        vec_vec_invmass_L_type5[i].resize(6);
        vec_vec_invmass_L_type51[i].resize(6);
    }

    float dcaprim_arr[6]={0.5,1,2,3,4,5};

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_K0_dca_cuts_m2_cuts_type3_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimK0_";
            n1+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_K0_type3[i][j] = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);

        }
    }

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_K0_dca_cutscombi2_m2_cuts_type3_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimK0_";
            n1+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_K0_type3[i+17][j] = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);

        }
    }

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_K0_dca_cuts_m2_cuts_type31_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimK0_";
            n1+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_K0_type31[i][j] = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);

        }
    }

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_K0_dca_cutscombi2_m2_cuts_type31_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimK0_";
            n1+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_K0_type31[i+17][j] = new TH1D(n1.Data(),n1.Data(),binningK0,0.4,0.6);

        }
    }

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_L_dca_cuts_Kcut_type5_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimL_";
            n1+=to_string(dcaprim_arr[j]);
            TString n2 = "invmass_L_dca_cuts_Kcut_type51_r_";
            n2+=to_string(5+i*5);
            n2+="_dcaprimL_";
            n2+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_L_type5[i][j] = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
            vec_vec_invmass_L_type51[i][j] = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        }
    }

    for(int i=0;i<17;i++)
    {
        for(int j=0;j<6;j++)
        {
            TString n1 = "invmass_L_dca_cuts_Kcut_distLtoS4_type5_r_";
            n1+=to_string(5+i*5);
            n1+="_dcaprimL_";
            n1+=to_string(dcaprim_arr[j]);
            TString n2 = "invmass_L_dca_cuts_Kcut_distLtoS4_type51_r_";
            n2+=to_string(5+i*5);
            n2+="_dcaprimL_";
            n2+=to_string(dcaprim_arr[j]);
            vec_vec_invmass_L_type5[i+17][j] = new TH1D(n1.Data(),n1.Data(),50*2,1.1,1.13);
            vec_vec_invmass_L_type51[i+17][j] = new TH1D(n2.Data(),n2.Data(),50*2,1.1,1.13);
        }
    }


    for(int i=0;i<17;i++)
    {
        TString n1="S_mass_ch3_normalband_r_";
        TString n2="S_mass_ch3_sideband_r_";
        n1+=to_string(5+i*5);
        n2+=to_string(5+i*5);
        TH1D* histo_S_normal = new TH1D(n1.Data(),n1.Data(),500,0,10);
        TH1D* histo_S_side  = new TH1D(n2.Data(),n2.Data(),500,0,10);
        vec_S_mass_ch3_normalband.push_back(histo_S_normal);
        vec_S_mass_ch3_sideband.push_back(histo_S_side);
    }



    vec_m_squared_type_2_p.resize(17*3);
    vec_m_squared_type_2_K_plus1.resize(17*3);
    vec_m_squared_type_2_K_plus2.resize(17*3);

    vec_m_squared_type_21_p.resize(17*3);
    vec_m_squared_type_21_K_plus1.resize(17*3);
    vec_m_squared_type_21_K_plus2.resize(17*3);


    for(int i=0;i<17;i++)
    {
        int j = 5+i*5;
        TString n1 = "m_squared_ch2_p_r_";
        TString n2 = "m_squared_ch2_K_plus1_r_";
        TString n3 = "m_squared_ch2_K_plus2_r_";

        TString n4 = "m_squared_ch21_p_r_";
        TString n5 = "m_squared_ch21_K_plus1_r_";
        TString n6 = "m_squared_ch21_K_plus2_r_";

        n1+=to_string(j);
        n2+=to_string(j);
        n3+=to_string(j);
        n4+=to_string(j);
        n5+=to_string(j);
        n6+=to_string(j);
    
        vec_m_squared_type_2_p[i] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
        vec_m_squared_type_2_K_plus1[i] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
        vec_m_squared_type_2_K_plus2[i] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);

        vec_m_squared_type_21_p[i] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);
        vec_m_squared_type_21_K_plus1[i] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
        vec_m_squared_type_21_K_plus2[i] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);


    }

    for(int i=0;i<17;i++)
    {
        int j = 5+i*5;
        TString n1 = "m_squared_ch2_dcacuts_p_r_";
        TString n2 = "m_squared_ch2_dcacuts_K_plus1_r_";
        TString n3 = "m_squared_ch2_dcacuts_K_plus2_r_";

        TString n4 = "m_squared_ch21_dcacuts_p_r_";
        TString n5 = "m_squared_ch21_dcacuts_K_plus1_r_";
        TString n6 = "m_squared_ch21_dcacuts_K_plus2_r_";

        n1+=to_string(j);
        n2+=to_string(j);
        n3+=to_string(j);
        n4+=to_string(j);
        n5+=to_string(j);
        n6+=to_string(j);
    
        vec_m_squared_type_2_p[i+17] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);
        vec_m_squared_type_2_K_plus1[i+17] = new TH1D(n2.Data(),n2.Data(),400,-0.4,1.4);
        vec_m_squared_type_2_K_plus2[i+17] = new TH1D(n3.Data(),n3.Data(),400,-0.4,1.4);

        vec_m_squared_type_21_p[i+17] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);
        vec_m_squared_type_21_K_plus1[i+17] = new TH1D(n5.Data(),n5.Data(),400,-0.4,1.4);
        vec_m_squared_type_21_K_plus2[i+17] = new TH1D(n6.Data(),n6.Data(),400,-0.4,1.4);


    }

    for(int i=0;i<17;i++)
    {
        int j = 5+i*5;
        TString n1 = "m_squared_ch2_dcacuts_and_Kaoncut_p_r_";

        TString n4 = "m_squared_ch21_dcacuts_and_Kaoncut_p_r_";

        n1+=to_string(j);
        n4+=to_string(j);
    
        vec_m_squared_type_2_p[i+17*2] = new TH1D(n1.Data(),n1.Data(),400,-0.4,1.4);

        vec_m_squared_type_21_p[i+17*2] = new TH1D(n4.Data(),n4.Data(),400,-0.4,1.4);

    }

    

    TString str1 = "S_mass_ch3_normalband_r_40";
    TString str2 = "S_mass_ch3_sideband_r_40";
    TString str3 = "S_mass_ch31_normalband_r_40";
    TString str4 = "S_mass_ch31_sideband_r_40";
    TH1D* histo_S_mass_n = new TH1D(str1.Data(),str1.Data(),200,0,10);
    TH1D* histo_S_mass_s = new TH1D(str2.Data(),str2.Data(),200,0,10);

    TH1D* histo_S_mass_n_2 = new TH1D(str3.Data(),str3.Data(),200,0,10);
    TH1D* histo_S_mass_s_2 = new TH1D(str4.Data(),str4.Data(),200,0,10);
    
    vec_S_mass_ch3.push_back(histo_S_mass_n);
    vec_S_mass_ch3.push_back(histo_S_mass_s);

    vec_S_mass_ch31.push_back(histo_S_mass_n_2);
    vec_S_mass_ch31.push_back(histo_S_mass_s_2);

    //dcacombi2
    str1 = "S_mass_ch3_dcacombi2_normalband_r_40";
    str2 = "S_mass_ch3_dcacombi2_sideband_r_40";
    str3 = "S_mass_ch31_dcacombi2_normalband_r_40";
    str4 = "S_mass_ch31_dcacombi2_sideband_r_40";
    TH1D* histo_S_mass_n_3 = new TH1D(str1.Data(),str1.Data(),200,0,10);
    TH1D* histo_S_mass_s_3 = new TH1D(str2.Data(),str2.Data(),200,0,10);

    TH1D* histo_S_mass_n_4 = new TH1D(str3.Data(),str3.Data(),200,0,10);
    TH1D* histo_S_mass_s_4 = new TH1D(str4.Data(),str4.Data(),200,0,10);
    
    vec_S_mass_ch3.push_back(histo_S_mass_n_3);
    vec_S_mass_ch3.push_back(histo_S_mass_s_3);

    vec_S_mass_ch31.push_back(histo_S_mass_n_4);
    vec_S_mass_ch31.push_back(histo_S_mass_s_4);




   cout<<"end of initialize"<<endl;
   

    
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
void copyV0params(Ali_AS_V0* V0_in,Ali_AS_V0* V0_out)
{
    float* xyz = V0_in->getxyz();
    float* pN = V0_in->getNpxpypz();
    float* pP = V0_in -> getPpxpypz();

    V0_out -> setxyz (xyz[0],xyz[1],xyz[2] ) ;
    V0_out -> setNpxpypz (pN[0],pN[1],pN[2] ) ;
    V0_out -> setPpxpypz (pP[0],pP[1],pP[2] ) ;

    V0_out ->setdcaV0(V0_in->getdcaV0());

}

void Ali_Dark_Matter_Read::copy_dm_params(Ali_AS_DM_particle* dm_in, Ali_AS_DM_particle* dm_out)
{
    dm_out->  set_primVertex(dm_in->get_primVertex());
    dm_out->  set_S1Vertex(dm_in->get_S1Vertex());
    dm_out->  set_S2Vertex(dm_in->get_S2Vertex()) ;
    dm_out->  set_S3Vertex(dm_in->get_S3Vertex()) ;
    dm_out->  set_DirSV1(dm_in->get_DirSV1());
    dm_out->  set_DirSV2(dm_in->get_DirSV2());
    dm_out->  set_DirSV3(dm_in->get_DirSV3());
    dm_out->  setN_V0s(dm_in->getN_V0s());
    dm_out->  settype(dm_in->gettype());
    dm_out->  set_tlv(dm_in->get_tlv());

    int type= dm_in->gettype();

    int numtracks = dm_in->getNumTracks();

    for(int i=0;i<numtracks;i++)
    {
        Ali_AS_Track* track_in = dm_in -> getTrack(i);
        Ali_AS_Track* track_out = dm_out->createTrack();
        copy_track_params(track_in,track_out);
    }

    for(int i=0;i<5;i++)
    {
        dm_out->settrack(i,dm_in->gettrack(i));
    }

    int num_V0s=0;
    if(type==2 || type==21) num_V0s = 1;
    if(type==3 || type==31) num_V0s = 2;
    if(type==5 || type==51) num_V0s = 2;

    for(int i=0;i<num_V0s;i++)
    {
        Ali_AS_V0* V0_in = dm_in->getV0(i);
        Ali_AS_V0* V0_out = dm_out->createV0();
        copyV0params(V0_in,V0_out);

    }

}


//


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
    //printf("Loop event number: %lld \n",event);
    counters[0]++;
    ////cout<<""<<endl;
    histo_counter_S->Fill(1);


    Event_active = event;

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;
    ////cout<<"a"<<endl;
    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Int_t    N_TRD_tracklets_online = AS_Event ->getNumTracklets(); // online tracklet
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    histo_n_tracks->Fill(N_tracks_event);

    histos_1D[7]->Fill(NumTracks);
    histos_1D[8]->Fill(EventVertexZ);

    TVector3 pos_primary_vertex;
    pos_primary_vertex.SetXYZ(EventVertexX,EventVertexY,EventVertexZ);

     ////cout<<"b"<<endl;

    //printf("Event vertex: %f", EventVertexX);
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    Int_t NumV0s = AS_Event ->getNumV0s();
    int numDMs = AS_Event->getNumDMs();
    //if(numDMs==0){return 0;}

    ////cout<<"NumV0s: "<<NumV0s<<endl;

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

    /*
    mixed_event_vec.resize(4);
    for(int i=0;i<4;i++)
    {
        mixed_event_vec[i].resize(10);
        for(int j=0;j<10;j++)
        {
            //mixed_event_vec[i][j].resize(5);
        }
    }
    */

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

    /*O
    const Float_t mass_proton = 0.93827208816 ;  //in GeV?
    const Float_t mass_pion = 0.139657061 ;  //in GeV?
    const Float_t mass_electron = 0.510998950 * 1e-3 ;  //in GeV?
    const double  mass_K = 0.493677 ;  //in GeV?

    double mass_K0 = 0.493677;
    double mass_Lambda = 1.1155683;
    double mass_neutron = 0.939565;
    double S_mass = -1;
    */

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
    vector<int> brute_force_type5;
    vector<int> brute_force_type51;

    vector<int> trackids_bits_shared;
    vector<int> trackids_similiar_tracks;
    vector<int> trackids_similiar_used;

    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;

     ////cout<<"c"<<endl;

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
                //cout<<"numbitsshared: "<<numbitsshared<<endl;
                //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

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
                //cout<<"dca: "<<dca<<endl;
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

    /*
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
    */
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

    Float_t path_closest_to_point = 0;
    Float_t dca_closest_to_point  = 0;
    //Float_t path_initA = 0.0;
    //Float_t path_initB = 30.0;

    vector<vector<int>> track_ids_all;
    track_ids_all.resize(NumV0s);


    int numNUCLEV = AS_Event->getNumNUCLEVs();
    vector<int> vec_trackids;
    
    vector<int> vec_skipped_nuclevs;

    //nuclear event loop
    //find double used tracks
    for(int i=0;i<numNUCLEV;i++)
    {
        //cout<<"nuclev: "<<i<<endl;
        bool skip=0;

        Ali_AS_NUCLEV* nuc = AS_Event-> getNUCLEV(i);
        TVector3 sec_vertex = nuc->get_secVertex();
        int numtracks = nuc ->getNumTracks();
        vector<int> vec_trackids_in_loop;

        //cout<<"numtracks: "<<numtracks<<endl;
        for(int track=0;track<numtracks;track++)
        {
           Ali_AS_Track* as_track = nuc->getTrack(track);
           int trackid = as_track->gettrackid();
           if ( check_if_int_is_in_vector(trackid,vec_trackids_in_loop) )
           {
                skip=1;
                break;
           }
           vec_trackids_in_loop.push_back(trackid);
        }

        if(skip)
        {
            vec_skipped_nuclevs.push_back(i);
            continue;
        }
        for(int track=0;track<numtracks;track++)
        {
            Ali_AS_Track* as_track = nuc->getTrack(track);
            int trackid = as_track->gettrackid();
            if ( check_if_int_is_in_vector(trackid,vec_trackids) )
            {
                skip=1;
                break;
            }
        }
        if(skip)
        {
            vec_skipped_nuclevs.push_back(i);
            continue;
        }

        for(int track=0;track<numtracks;track++)
        {
            Ali_AS_Track* as_track = nuc->getTrack(track);
            int trackid = as_track->gettrackid();
            vec_trackids.push_back(trackid);
        }


        //FindDCAHelixPoint(sec_vertex,as_track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
        //   if(path_closest_to_point<0)
        //   {
               //skip=1;
               //break;
        //   }
        counter_NUCLEV++;
        //cout<<""<<endl;
    }



    
    //TEveManager::Create();

    //numNUCLEV = -1;

    //cout<<"a"<<endl;
    for(int i=0;i<numNUCLEV;i++)
    {
        if(check_if_int_is_in_vector(i,vec_skipped_nuclevs)){continue;}

        Ali_AS_NUCLEV* nuc = AS_Event-> getNUCLEV(i);
        int numtracks = nuc ->getNumTracks();
        TVector3 sec_vertex = nuc->get_secVertex();
        double radius = sqrt ( (sec_vertex[0]-EventVertexX)*(sec_vertex[0]-EventVertexX)+ (sec_vertex[1]-EventVertexY)*(sec_vertex[1]-EventVertexY ) + (sec_vertex[2]-EventVertexZ)*(sec_vertex[2]-EventVertexZ ) );
        //double radius_ortho = sqrt ((sec_vertex[0]-EventVertexX)*(sec_vertex[0]-EventVertexX)+ (sec_vertex[1]-EventVertexY)*(sec_vertex[1]-EventVertexY ));
        double radius_ortho = sqrt (sec_vertex[0]*sec_vertex[0]+ sec_vertex[1]*sec_vertex[1]);

        TVector3 vec_primtosec;
        vec_primtosec.SetXYZ((sec_vertex[0]-EventVertexX),(sec_vertex[1]-EventVertexY),(sec_vertex[2]-EventVertexZ));
        TVector3 unit_prim_to_sec;
        unit_prim_to_sec = vec_primtosec.Unit();

        Ali_AS_Track* track0 = nuc->getTrack(0);
        Ali_AS_Track* track1 = nuc->getTrack(1);

        double dcaAprim = track0->getdca();
        double dcaBprim = track1->getdca();

        int itscls0 = track0->getNITScls();
        int itscls1 = track1->getNITScls();

        histo_its->Fill(itscls0);
        histo_its->Fill(itscls1);

        vector<int> binary_array0 = inttobinaryarray(itscls0);
        vector<int> binary_array1 = inttobinaryarray(itscls1);

        int firstlayer0 = -1;
        int firstlayer1 = -1;

        float arr_its[6]={4.,7.,15.,24.,39.,44.};

        /*
        TVector3 sec_vertex_found = find_sec_vertex(track0,track1);
        double delta_x = sec_vertex_found[0]-sec_vertex[0];
        double delta_y = sec_vertex_found[1]-sec_vertex[1];
        double delta_z = sec_vertex_found[2]-sec_vertex[2];

        histo_delta[0]->Fill(delta_x);
        histo_delta[1]->Fill(delta_y);
        histo_delta[2]->Fill(delta_z);
        */

        for(int j=7;j>=0;j--)
        {
            if(binary_array0[j]==1)
            {
                firstlayer0 = 7-j;
                break;
            }
        }

        for(int j=7;j>=0;j--)
        {
            if(binary_array1[j]==1)
            {
                firstlayer1 = 7-j;
                break;
            }
        }

        bool itscut=1;

        if(firstlayer0!=-1)
        {
            if(radius>arr_its[firstlayer0]){itscut=0;}
        }

        if(firstlayer1!=-1)
        {
            if(radius>arr_its[firstlayer1]){itscut=0;}
        }



        if(i==0 && firstlayer0!=-1)
        {
            printf("itscls0: %d \n",itscls0);
            print_int_vector(binary_array0);
            //cout<<"firstlayer0: "<<firstlayer0<<endl;
        }


        TLorentzVector tlv0 = track0->get_TLV_part();
        TLorentzVector tlv1 = track1->get_TLV_part();

        TLorentzVector tlv_ges = tlv0 + tlv1;
        TVector3 dir;
        TVector3 unit_dir;
        dir.SetXYZ(tlv_ges.Px(),tlv_ges.Py(),tlv_ges.Pz());
        unit_dir = dir.Unit();

        double dot_product = unit_dir.Dot(unit_prim_to_sec);

        if(dot_product<0.8){continue;}


        UShort_t ntpcclsP = track0->getNTPCcls();
        UShort_t nitsclsP = track0->getNITScls();
        UShort_t statusP = track0->getStatus();
        Float_t tpcchi2P = track0->getTPCchi2();

        UShort_t ntpcclsN = track1->getNTPCcls();
        UShort_t nitsclsN = track1->getNITScls();
        UShort_t statusN = track1->getStatus();
        Float_t tpcchi2N = track1->getTPCchi2();

        //if(ntpcclsP<120 || ntpcclsN<120){continue;}

        bool dcaprim2check=1;
        bool dcaprim5check=1;
        bool dcaprim7check=1;
        bool sigmacheck=0;
        bool momentumcheck=1;
        bool anglecheck=1;

        if(radius<2){continue;}

        TLorentzVector tlv = track0->get_TLV_part();
        double momentum = tlv.P();
        //if(momentum<0.1){continue;}

        tlv = track1->get_TLV_part();
        momentum = tlv.P();
        //if(momentum<0.1){continue;}

        float dcaA,dcaB;
        FindDCAHelixPoint(sec_vertex,track0,path_initA,path_initB,path_closest_to_point,dcaA);
        //cout<<"dca: "<<dca_closest_to_point<<endl;
        if(dcaA>0.5){continue;}

        FindDCAHelixPoint(sec_vertex,track1,path_initA,path_initB,path_closest_to_point,dcaB);
        if(dcaB>0.5){continue;}

        double dcaV0_ESD = nuc->getdcaV0_ESD();
        if(dcaV0_ESD>0.5){continue;}

        vector<TVector3> vec_unit_dir;
        Double_t r1[3];
        Double_t r2[3];

        for(int track=0;track<2;track++)
        //for(int track=0;track<numtracks;track++)
        {
            Ali_AS_Track* as_track = nuc->getTrack(track);
            double dcaprim = fabs ( as_track->getdca() );

            float dca_closest_to_point;
            FindDCAHelixPoint(pos_primary_vertex,as_track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            dcaprim = dca_closest_to_point;

            if(dcaprim<2.){dcaprim2check=0;}
            if(dcaprim<5.){dcaprim5check=0;}
            if(dcaprim<7.){dcaprim7check=0;}

            TLorentzVector tlv = as_track->get_TLV_part();
            double momentum = tlv.P();
            if(0.5>momentum || momentum>2){momentumcheck=0;}

            FindDCAHelixPoint(sec_vertex,as_track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            as_track->Evaluate(path_closest_to_point,r1);
            as_track->Evaluate(path_closest_to_point+0.01,r2);

            TVector3 dir;
            dir.SetXYZ(r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]);
            TVector3 unit_dir = dir.Unit();
            vec_unit_dir.push_back(unit_dir);

        }

        double angle = vec_unit_dir[0].Angle(vec_unit_dir[1]);
        histo_angle->Fill(angle);

        double z = sec_vertex[2];

        if(z>0 && z<10){ histo_angle_z_0_10 ->Fill(angle); }
        if(z>50 && z<60){ histo_angle_z_50_60 ->Fill(angle); }


        for(int i=2;i<numtracks;i++)
        {
            Ali_AS_Track* as_track = nuc->getTrack(i);
            double sigma = fabs ( as_track-> getnsigma_pi_TPC() );
            if(sigma>2.5){sigmacheck=1;}

        }
        //if(dcacheck){continue;}

        if(angle>3.){anglecheck=0;}

        if(sigmacheck){continue;}

        if(numtracks==2)
        {
            histo_radius_nuclev_2_or_more-> Fill(radius_ortho);
            x_and_y_nuclev_3->Fill(sec_vertex[0],sec_vertex[1]);
            vec_radius_ortho_vs_z[0]->Fill(sec_vertex[2],radius_ortho);
            if(dcaprim2check){ vec_radius_ortho_vs_z[0+4]->Fill(sec_vertex[2],radius_ortho); }
            if(dcaprim5check){ vec_radius_ortho_vs_z[0+8]->Fill(sec_vertex[2],radius_ortho); }


            Ali_AS_Track* as_track1 = nuc->getTrack(0);
            Ali_AS_Track* as_track2 = nuc->getTrack(1);

            TLorentzVector tlv1  = as_track1 -> get_TLV_part();
            TLorentzVector tlv2  = as_track2 -> get_TLV_part();

            double eta1 = tlv1.Eta();
            double eta2 = tlv2.Eta();

            if(15<radius_ortho && radius_ortho<25 && -5<sec_vertex[2] && sec_vertex[2] < 5)
            {
                histo_eta -> Fill(eta1,eta2);
                histo_eta_larger_range -> Fill(eta1,eta2);
                histo_delta_eta->Fill(eta2-eta1);

            }

            if(fabs(sec_vertex[2])>20)
            {
                histo_eta_z_cut -> Fill(eta1,eta2);
                histo_delta_eta_z_cut->Fill(eta2-eta1);

                histo_sum_eta_z_cut->Fill(eta1+eta2);
            }

            if( fabs(eta1+eta2) > 0.2 && fabs(eta1-eta2) > 0.2 && dcaprim5check )
            {
                radius_ortho_vs_z_eta->Fill(sec_vertex[2],radius_ortho);
                if(momentumcheck){radius_ortho_vs_z_eta4->Fill(z,radius_ortho);}
                if(anglecheck){radius_ortho_vs_z_eta5->Fill(z,radius_ortho);}
            }

            if( ( fabs(eta1+eta2) < 0.2 || fabs(eta1-eta2) < 0.2 ) && dcaprim5check )
            {
                radius_ortho_vs_z_eta2->Fill(sec_vertex[2],radius_ortho);
            }

            if( ( fabs(eta1+eta2) < 0.2 && fabs(eta1-eta2) < 0.2 ) && dcaprim5check )
            {
                radius_ortho_vs_z_eta3->Fill(sec_vertex[2],radius_ortho);
            }





        }

        if(numtracks==3)
        {
            histo_radius_nuclev_3_or_more-> Fill(radius_ortho);
            x_and_y_nuclev_4->Fill(sec_vertex[0],sec_vertex[1]);
            vec_radius_ortho_vs_z[1]->Fill(sec_vertex[2],radius_ortho);
            if(dcaprim2check){ vec_radius_ortho_vs_z[1+4]->Fill(sec_vertex[2],radius_ortho); }
            if(dcaprim5check && momentumcheck){ vec_radius_ortho_vs_z[1+8]->Fill(sec_vertex[2],radius_ortho); }
        }

        if(numtracks<=4)
        {
            if(dcaprim5check){histo_radius_nuclev_4_or_more-> Fill(radius_ortho);}
            if(dcaprim5check && momentumcheck){histo_radius_nuclev_4_or_more_momentum-> Fill(radius_ortho);}
            x_and_y_nuclev_5->Fill(sec_vertex[0],sec_vertex[1]);
            vec_radius_ortho_vs_z[2]->Fill(sec_vertex[2],radius_ortho);
            vec_histo_radius_nuclev[0]->Fill(radius_ortho);
            vec_radius_ortho_vs_z[15]->Fill(sec_vertex[2],radius_ortho);
            if(dcaprim2check){ vec_radius_ortho_vs_z[2+4]->Fill(sec_vertex[2],radius_ortho);}

            if(dcaprim5check)
            {
                for(Int_t i=0;i<10;i++)
                {
                    int range_min=-100+ i*20;
                    int range_max=-80 + i*20;
                    if(range_min<sec_vertex[2] && range_max>sec_vertex[2])
                    {
                        vec_x_y_slices[i]->Fill(sec_vertex[0],sec_vertex[1]);
                        vec_radius_slices[i]->Fill(radius_ortho);
                    }
                }
                for(Int_t i=0;i<20;i++)
                {
                    int range_min=0+ i*10;
                    int range_max=10 + i*10;

                    if(range_min<radius_ortho && radius_ortho<range_max)
                    {
                        vec_z_slices_in_radius[i]->Fill(sec_vertex[2]);
                    }
                    
                }


                vec_radius_ortho_vs_z[2+8]->Fill(sec_vertex[2],radius_ortho);
                vec_histo_radius_nuclev[1]->Fill(radius_ortho);
                vec_radius_ortho_vs_z[16]->Fill(sec_vertex[2],radius_ortho);

                if(itscut)
                {
                    vec_histo_radius_nuclev[6]->Fill(radius_ortho);
                    vec_radius_ortho_vs_z[21]->Fill(sec_vertex[2],radius_ortho);
                }

                if(momentumcheck)
                {
                    vec_radius_ortho_vs_z[12]->Fill(sec_vertex[2],radius_ortho);
                    vec_histo_radius_nuclev[2]->Fill(radius_ortho);
                    vec_radius_ortho_vs_z[17]->Fill(sec_vertex[2],radius_ortho);
                    
                    if(itscut)
                    {
                        vec_histo_radius_nuclev[5]->Fill(radius_ortho);
                        vec_radius_ortho_vs_z[20]->Fill(sec_vertex[2],radius_ortho);
                    }
                }
            }

            if(dcaprim7check)
            {
                vec_radius_ortho_vs_z[13]->Fill(sec_vertex[2],radius_ortho);
                vec_histo_radius_nuclev[3]->Fill(radius_ortho);
                vec_radius_ortho_vs_z[18]->Fill(sec_vertex[2],radius_ortho);
                if(momentumcheck)
                {
                    vec_radius_ortho_vs_z[14]->Fill(sec_vertex[2],radius_ortho);
                    vec_histo_radius_nuclev[4]->Fill(radius_ortho);
                    vec_radius_ortho_vs_z[19]->Fill(sec_vertex[2],radius_ortho);

                }
            }
        }

        if(numtracks>=5)
        {
            histo_radius_nuclev_5_or_more-> Fill(radius_ortho);
            x_and_y_nuclev_3->Fill(sec_vertex[0],sec_vertex[1]);
            vec_radius_ortho_vs_z[3]->Fill(sec_vertex[2],radius_ortho);
            if(dcaprim2check){ vec_radius_ortho_vs_z[3+4]->Fill(sec_vertex[2],radius_ortho); }
            if(dcaprim5check && momentumcheck){ vec_radius_ortho_vs_z[3+8]->Fill(sec_vertex[2],radius_ortho); }
        }

        if(dcaprim5check && radius_ortho<1.)
        {
            for(int track=0;track<numtracks;track++)
            {
                //cout<<"track: "<<track<<endl;
                Ali_AS_Track* as_track = nuc->getTrack(track);

                FindDCAHelixPoint(pos_primary_vertex,as_track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
                double dcaprim = fabs ( as_track->getdca() );

                //cout<<"dcaprim1: "<<dca_closest_to_point<<endl;
                //cout<<"dcaprim2: "<<dcaprim<<endl;
                //if(counter_draw==-1) { Draw_TPC_track(as_track, kGreen, 2.,track,pos_primary_vertex); }
            }
            draw=0;
            counter_draw++;
        }

        for(int track=0;track<numtracks;track++)
        {
            Ali_AS_Track* as_track = nuc->getTrack(track);

            FindDCAHelixPoint(sec_vertex,as_track,path_initA,path_initB,path_closest_to_point,dca_closest_to_point);
            //cout<<"path: "<<path_closest_to_point<<endl;
            //histo_path->Fill(path_closest_to_point);

            double path = path_closest_to_point;
            //double path = nuc->getpath(track);
            histo_path->Fill(path);

            int itscls = as_track->getNITScls();
            //cout<<"nitscls: "<<itscls<<endl;
            histo_nitscls->Fill(itscls);

            if(radius>30)
            {
                if(itscls>0){histo_path_has_its -> Fill(path);}
                if(itscls==0){histo_path_has_no_its -> Fill(path);}
            }

            vector<int> binary_array = inttobinaryarray(itscls);
            if(getnumones(binary_array)>2 && radius>30)
            {
                histo_path_has_atleast_3_its_hits->Fill(path);

            }


        }

        


    }
    

    //end of nuclear event loop

    //DM loop
    //int numDMs = AS_Event->getNumDMs();
    histo_numberDMs->Fill(numDMs);
    
    cout<<"dmloop"<<endl;
    for(int dm_loop = 0; dm_loop<numDMs; dm_loop++)
    {
        histo_counter_S->Fill(2);

        Ali_AS_DM_particle* DM  = AS_Event->getDMparticle(dm_loop);
        int type = DM->gettype();
        cout<<"type: "<<type<<endl;
        histo_type_of_S->Fill(type);

        int numtracks = DM->getNumTracks();

        TVector3 S_vertex_pos = DM->get_S1Vertex();
        TVector3 vec_primtosec;
        vec_primtosec.SetXYZ((S_vertex_pos[0]-EventVertexX),(S_vertex_pos[1]-EventVertexY),(S_vertex_pos[2]-EventVertexZ));
        double radius = vec_primtosec.Mag();
        double radius_from_origin = S_vertex_pos.Mag();
        bool dcaprimcheck=1;
        bool dcaprimcheck_1=1;

        S_radius_from_origin->Fill(radius_from_origin);
        
        for(int i=0;i<5;i++)
        {
            if(type==i+1){vec_S_radius_from_origin[i]->Fill(radius_from_origin);}
        }

        if(type==5) vec_S_radius_from_origin[4]->Fill(radius_from_origin);
        if(type==51) vec_S_radius_from_origin[5]->Fill(radius_from_origin);


        cout<<"numtracks: "<<numtracks<<endl;
        //dE/dx plots
        for(int i=0;i<numtracks;i++)
        {
            //if(type==5 || type==51){continue;}
            //if(type>5){continue;}
            Ali_AS_Track* track = DM->getTrack(i);
            int trackid = track->gettrackid();

            //dEdx
            Float_t TPCdEdx   = track->getTPCdEdx();
            Float_t tofsignal = track->getTOFsignal();
            Float_t dca       = track->getdca();
            Float_t tracklength = track->getTrack_length();
            int charge;

            if(dca>0){charge = 1;}
            else {charge = -1;}
            TLorentzVector tlv = track->get_TLV_part();
            double momentum = tlv.P();

            dEdx_vs_charge_dot_momentum_S_cand->Fill(charge*momentum,TPCdEdx);

            //vec_dEdx[type-1]->Fill(charge*momentum,TPCdEdx);

            histo_pt_S->Fill(tlv.Pt());

            if(tofsignal<99990) {tprof->Fill(1,1);}
            else {tprof->Fill(1,0); }

            double sigma_pi = track -> getnsigma_pi_TPC();
            double sigma_K = track -> getnsigma_K_TPC();
            double sigma_p = track -> getnsigma_p_TPC();

            if(fabs(sigma_pi)<2.5){vec_dEdx_S[0]->Fill(charge*momentum,TPCdEdx);};
            if(fabs(sigma_K)<2.5){vec_dEdx_S[1]->Fill(charge*momentum,TPCdEdx);};
            if(fabs(sigma_p)<2.5){vec_dEdx_S[2]->Fill(charge*momentum,TPCdEdx);};

            if(fabs(sigma_pi<2.))
            {
                double m_squared_pi =calculate_m_squared_by_TOF(track);
                vec_mass_squared_dEdx_selected[0]->Fill(m_squared_pi);
            }

            if(fabs(sigma_K<2.))
            {
                double m_squared_K =calculate_m_squared_by_TOF(track);
                vec_mass_squared_dEdx_selected[1]->Fill(m_squared_K);
            }

            if(fabs(sigma_p<2.))
            {
                double m_squared_p =calculate_m_squared_by_TOF(track);
                vec_mass_squared_dEdx_selected[2]->Fill(m_squared_p);

            }


            float dcaprim;
            FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
            if(dcaprim<0.5){dcaprimcheck=0;}
            if(dcaprim<1.){dcaprimcheck_1=0;}
        }

        /*
        //for all types fill 10 for 3D view  //type 2,21,3,31,5,51
        int typarr[6]={2,21,3,31,5,51};
        for(int t = 0; t<6 ; t++)
        {
            if(type==typarr[t])
            {
                if(counter_arr[t]>=10){continue;}
                AS_DM_particle->clearTrackList();
                AS_DM_particle->clearV0List();
                copy_dm_params(DM,AS_DM_particle);
                Tree_AS_DM_particle ->Fill();
                counter_arr[t]++;
            }
        }
        */

        //type 1--------------------------------------------------------------------------------------
        if(type==1 || type==11)
        {
            if(type==1) histo_counter_S->Fill(101);
            if(type==11) histo_counter_S->Fill(111);

            TVector3 S_vertex_pos = DM -> get_S1Vertex();
            TVector3 K0_pos = DM->get_S2Vertex();
            TVector3 L_pos = DM->get_S3Vertex();

            Ali_AS_Track* track_S_pion1  = DM -> getTrack(0);
            Ali_AS_Track* track_S_pion2  = DM -> getTrack(1);
            Ali_AS_Track* track1_K0 = DM -> getTrack(2);
            Ali_AS_Track* track2_K0 = DM -> getTrack(3);

            Ali_AS_Track* track_L_pi;
            Ali_AS_Track* track_L_p;

            if(type==1)
            {
                track_L_pi = DM -> getTrack(4);
                track_L_p = DM -> getTrack(5);
            }
            if(type==11)
            {
                track_L_pi = DM -> getTrack(5);
                track_L_p = DM -> getTrack(4);
            }

            Ali_AS_V0* S_V0 = DM->getV0(0);
            Ali_AS_V0* K0_V0 = DM->getV0(1);
            Ali_AS_V0* L_V0 = DM->getV0(2);

            double m_squared_S_pion1 = calculate_m_squared_by_TOF(track_S_pion1);
            double m_squared_S_pion2 = calculate_m_squared_by_TOF(track_S_pion2);
            double m_squared_p = calculate_m_squared_by_TOF(track_L_p);

            if(type==1) m_squared_proton_type1->Fill(m_squared_p);
            if(type==11) m_squared_proton_type11->Fill(m_squared_p);


            bool dcaprimcheck_type1=1;
            for(int i=0;i<2;i++)  //check dcaprim only for pions!
            {
                Ali_AS_Track* track = DM->getTrack(i);
                float dcaprim;
                FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
                if(dcaprim<1.){dcaprimcheck_type1=0;}
            }

            double invmass_antilambda =  get_invariant_mass(track_L_pi, track_L_p , mass_pion,mass_proton , S_vertex_pos);
            //double invmass_K0 =  get_invariant_mass(track1_K0, track2_K0, mass_pion, mass_pion , S_vertex_pos) ;

            TLorentzVector tlv_L;
            if(type==1) tlv_L = get_tlv_by_V0(mass_pion, mass_proton, L_V0->getPpxpypz(), L_V0->getNpxpypz() );
            if(type==11) tlv_L = get_tlv_by_V0(mass_proton, mass_pion, L_V0->getPpxpypz(), L_V0->getNpxpypz() );
            double invmass_antilambda_1 = tlv_L.M();

            TLorentzVector tlv_K0;
            if(type==1 || type==11) tlv_K0 = get_tlv_by_V0(mass_pion, mass_pion, K0_V0->getPpxpypz(), K0_V0->getNpxpypz() );
            double invmass_K0 = tlv_K0.M();

            TVector3 dir_L;
            dir_L.SetXYZ(tlv_L.Px(),tlv_L.Py(),tlv_L.Pz());
            double dist = calculateMinimumDistanceStraightToPoint(L_pos,dir_L,pos_primary_vertex);

            if(dist>1.5 && type==1)
            {


            }

            if(type==1)
            {
                histo_invmass_Lambda_type1->Fill(invmass_antilambda_1);
                histo_invmass_Lambda_type1_by_track->Fill(invmass_antilambda);
                if(dcaprimcheck ) histo_invmass_Lambda_type1_dcaprim->Fill(invmass_antilambda_1);
                if(dcaprimcheck_1 ) histo_invmass_Lambda_type1_dcaprim_1->Fill(invmass_antilambda_1);
                histo_invariantmass_K0_type1->Fill(invmass_K0);
                if(dcaprimcheck) histo_invariantmass_K0_type1_dcaprim->Fill(invmass_K0);
            }
            if(type==11)
            {
                histo_invmass_Lambda_type11_by_track->Fill(invmass_antilambda);
                histo_invmass_Lambda_type11->Fill(invmass_antilambda_1);
                histo_invariantmass_K0_type11->Fill(invmass_K0);

            }

            cout<<"invmass: "<< invmass_antilambda <<endl;
            cout<<"invmass2: "<< invmass_antilambda_1 <<endl;

            histo_diff_inv_mass_type1 ->Fill(invmass_antilambda_1-invmass_antilambda);

            //overlap PID:
            if( overlap_PID(track_L_p, "p") && type==1)
            {
                histo_invmass_Lambda_overlap_cut_type1-> Fill(invmass_antilambda_1);
                if(dcaprimcheck) histo_invmass_Lambda_overlap_cut_dcaprim_type1-> Fill(invmass_antilambda_1);
                histo_invariantmass_K0_type1_overlap->Fill(invmass_K0);
            }

            double invmass_L = invmass_antilambda_1;
            //double invmass_K0;

            if(invmass_L < (1.1157+0.001495*2) && invmass_L > (1.1157-0.001495*2) )
            {
                histo_invariantmass_K0_type1_cut_on_L->Fill(invmass_K0);
            }

            //normal analysis  both 2-sigma
            if(invmass_L < (1.1157+0.001495*2) && invmass_L > (1.1157-0.001495*2)
               && invmass_K0 < (0.4981+0.0042*2) && invmass_K0 > (0.4981-0.0042*2))
            {
                if(type==1)histo_counter_S->Fill(102);
                if(type==11) histo_counter_S->Fill(112);

                if( overlap_PID(track_L_p, "p") )
                {
                    if(type==1)
                    {
                        histo_counter_S->Fill(103);
                        if(dcaprimcheck) histo_counter_S->Fill(104);
                        if(dcaprimcheck && radius > 10.) histo_counter_S->Fill(105);
                    }
                    if(type==11)
                    {
                        histo_counter_S->Fill(113);
                        if(dcaprimcheck) histo_counter_S->Fill(114);
                        if(dcaprimcheck && radius > 10.) histo_counter_S->Fill(115);
                    }
                }

                if(radius>10)
                {
                    if(type==1) histo_counter_S->Fill(11);

                    if(dcaprimcheck_type1)
                    {
                        /*
                        AS_DM_particle->clearTrackList();
                        copy_dm_params(DM,AS_DM_particle);
                        Tree_AS_DM_particle ->Fill();
                        histo_counter_S->Fill(21);
                        */
                    }
                }


            }

            //sideband analysis   K0 2-sigma ; Lambda sideband
            if( ( (invmass_L < (1.1157+0.001495*4) && invmass_L > (1.1157+0.001495*2))
                 || ( invmass_L > (1.1157-0.001495*4) && invmass_L < (1.1157-0.001495*2)) )
               &&  ( invmass_K0 < (0.4981+0.0042*2) && invmass_K0 > (0.4981-0.0042*2) ) )
            {

            }

            //sideband analysis   K0 sideband ; Lambda 2-sigma
            if( (invmass_L < (1.1157+0.001495*2) && invmass_L > (1.1157-0.001495*2))
               && (invmass_K0 < (0.4981+0.0042*4) && invmass_K0 > (0.4981+0.0042*2))
                  || ( invmass_K0 > (0.4981-0.0042*4) && invmass_K0 < (0.4981-0.0042*2) ) )
            {

            }


            if(radius>20)
            {
                histo_counter_S->Fill(3);
                
            }

        }


        //-------------------------------------------------------------------------------------------
        //typ 2-------------------------------------------------------------------------------------------
        if(type==2 || type==21)
        {
            if(type==2) histo_counter_S->Fill(121);
            if(type==21) histo_counter_S->Fill(131);

            TVector3 S_vertex;
            S_vertex = DM->get_S1Vertex();
            //printf("S vertex: %f %f %f \n",S_vertex[0],S_vertex[1],S_vertex[2]);

            Ali_AS_Track* track0;  //K+
            Ali_AS_Track* track1;  //antip
            Ali_AS_Track* track2; //add K+

            if(type==2)
            {
                track0 = DM->getTrack(0); //K+
                track1 = DM->getTrack(1); //antip
                track2 = DM->getTrack(2); //add K+
            }
            if(type==21)
            {
                track0 = DM->getTrack(1); //K-
                track1 = DM->getTrack(0); //p
                track2 = DM->getTrack(2); //add K-
            }

            //do m2 cut if available
            if ( !m2_cut_if_available(track0, "K")) continue;
            if ( !m2_cut_if_available(track1, "p")) continue;
            if ( !m2_cut_if_available(track2, "K")) continue;


            Ali_AS_V0* S_V0 =  DM->getV0(0);

            //if(checkifsimiliar(track0,track1)){continue;}
            //if(checkifsimiliar(track0,track2)){continue;}
            //if(checkifsimiliar(track1,track2)){continue;}

            Float_t dcaprim_p,dcaprim_K1,dcaprim_K2;
            FindDCAHelixPoint(pos_primary_vertex,track0,path_initA,path_initB,path_closest_to_point,dcaprim_K1);
            FindDCAHelixPoint(pos_primary_vertex,track1,path_initA,path_initB,path_closest_to_point,dcaprim_p);
            FindDCAHelixPoint(pos_primary_vertex,track2,path_initA,path_initB,path_closest_to_point,dcaprim_K2);


            TLorentzVector tlv_add_K = track2 -> get_TLV_part();

            double m_squared_add_Kaon = calculate_m_squared_by_TOF(track2);
            histo_m_squared_kaon_no_other_cuts_type2 ->Fill(m_squared_add_Kaon);

            double m_squared_Kaon =calculate_m_squared_by_TOF(track0);
            double m_squared_antip =calculate_m_squared_by_TOF(track1);

            m_squared_anti_proton_type2->Fill(m_squared_antip);
            mass_squared_all_kaons_type2->Fill(m_squared_add_Kaon);
            mass_squared_all_kaons_type2->Fill(m_squared_Kaon);

            //radius dependend m_squared
            for(int i=0;i<17;i++)
            {
                float r = 5+i*5;
                if(radius>r)
                {
                    //no dca cuts
                    if(type==2)
                    {
                        vec_m_squared_type_2_p[i]->Fill(m_squared_antip);
                        vec_m_squared_type_2_K_plus1[i]->Fill(m_squared_Kaon);
                        vec_m_squared_type_2_K_plus2[i]->Fill(m_squared_add_Kaon);

                    }

                    if(type==21)
                    {
                        vec_m_squared_type_21_p[i]->Fill(m_squared_antip);
                        vec_m_squared_type_21_K_plus1[i]->Fill(m_squared_Kaon);
                        vec_m_squared_type_21_K_plus2[i]->Fill(m_squared_add_Kaon);
                    }

                    //dcacuts
                    if(dcaprim_p>2. && dcaprim_K1 > 2. && dcaprim_K2 > 2.)
                    {
                        if(type==2)
                        {
                            vec_m_squared_type_2_p[i+17]->Fill(m_squared_antip);
                            vec_m_squared_type_2_K_plus1[i+17]->Fill(m_squared_Kaon);
                            vec_m_squared_type_2_K_plus2[i+17]->Fill(m_squared_add_Kaon);

                        }

                        if(type==21)
                        {
                            vec_m_squared_type_21_p[i+17]->Fill(m_squared_antip);
                            vec_m_squared_type_21_K_plus1[i+17]->Fill(m_squared_Kaon);
                            vec_m_squared_type_21_K_plus2[i+17]->Fill(m_squared_add_Kaon);
                        }

                        //additional overlapcuts
                        if( overlap_PID(track0, "K") &&  overlap_PID(track2, "K") )
                        {
                            if(type==2) vec_m_squared_type_2_p[i+17*2]->Fill(m_squared_antip);
                            if(type==21) vec_m_squared_type_21_p[i+17*2]->Fill(m_squared_antip);
                        }


                    }

                }
            }



            if( overlap_PID(track0, "K") &&  overlap_PID(track1, "p") &&  overlap_PID(track2, "K") )
            {
                mass_squared_all_kaons_overlap_cuts_type2->Fill(m_squared_add_Kaon);
                mass_squared_all_kaons_overlap_cuts_type2->Fill(m_squared_Kaon);
                if(type==2) histo_counter_S->Fill(122);
                if(type==2 && dcaprimcheck) histo_counter_S->Fill(123);
                if(type==2 && dcaprimcheck && radius > 10) histo_counter_S->Fill(124);

                if(type==21) histo_counter_S->Fill(132);
                if(type==21 && dcaprimcheck) histo_counter_S->Fill(133);
                if(type==21 && dcaprimcheck && radius > 10) histo_counter_S->Fill(134);

            }



            if( m_squared_Kaon> 0.2 && m_squared_Kaon < 0.35  && m_squared_antip > 0.6 &&  m_squared_antip< 1.2
               && m_squared_add_Kaon>0.2 && m_squared_add_Kaon<0.35)
            {
                histo_m_squared_kaon_with_all_other_cuts_type2->Fill(m_squared_add_Kaon);
                //histo_invariant_mass_kaon_with_other_cuts->Fill(invariantmass_add_Kaon);
                histo_counter_S->Fill(13);
                if(radius>10)
                {
                    
                    if(dcaprimcheck==1)
                    {
                        /*
                        AS_DM_particle->clearTrackList();
                        copy_dm_params(DM,AS_DM_particle);
                        Tree_AS_DM_particle ->Fill();
                        histo_counter_S->Fill(22);
                        */
                    }
                }

            }

            if(m_squared_antip > 0.8  &&  m_squared_antip< 1.1  )
            {
                histo_m_squared_kaon_with_some_other_cuts_type2->Fill(m_squared_add_Kaon);
                histo_counter_S->Fill(12);
            }

                                
            if(radius>20)
            {
                histo_counter_S->Fill(4);
                
            }
        }

        
        if(type==3 || type==31)
        {
            if(type==3) histo_counter_S->Fill(141);
            if(type==31) histo_counter_S->Fill(151);

            Ali_AS_Track* track_K;
            Ali_AS_Track* track_p;

            if(type==3)
            {
                track_K = DM->getTrack(0);
                track_p = DM->getTrack(1);
            }

            if(type==31)
            {
                track_K = DM->getTrack(1);
                track_p = DM->getTrack(0);
            }


            Ali_AS_Track* track_add_pi  = DM->getTrack(2);

            Ali_AS_Track* track_K0_pi1 = DM->getTrack(3);        //jetzt richtig! //noch falsch!!!! da vorher falsch auf Grid
            Ali_AS_Track* track_K0_pi2 = DM->getTrack(4);

            double m2_p = calculate_m_squared_by_TOF(track_p);

            if( fabs(track_p->getnsigma_p_TPC()) <2.5 )
            {
                if(type==3) m_squared_p_ch3->Fill(m2_p);
            }


            //tlv
            TLorentzVector tlv_S = DM->get_tlv();
            TLorentzVector tlv_proton;
            tlv_proton.SetPxPyPzE(0.,0.,0.,mass_proton);

            tlv_S = tlv_S-tlv_proton;

            double S_mass = tlv_S.M();

            //double TOF efficiency for Kaon
            double m_squared_K_1 = calculate_m_squared_by_TOF(track_K);
            TOF_eff_K_ch3->Fill(1);
            if(m_squared_K_1!=-1) TOF_eff_K_ch3->Fill(2);

            //------------------------------------------------------------------------------------
            //msquared
            double m_squared_K = calculate_m_squared_by_TOF(track_K);
            double m_squared_p = calculate_m_squared_by_TOF(track_p);
            
            double m_squared_K01 = calculate_m_squared_by_TOF(track_K0_pi1);
            double m_squared_K02 = calculate_m_squared_by_TOF(track_K0_pi2);

            double m_squared_add_pi = calculate_m_squared_by_TOF(track_add_pi);

            float dca_arr[5]={0.5,1.,1.5,2.,3.};

            float dcaprim_p;
            FindDCAHelixPoint(pos_primary_vertex,track_p,path_initA,path_initB,path_closest_to_point,dcaprim_p);

            if(radius>25)
            {
                for(int i=0;i<5;i++)
                {
                    if(dcaprim_p<dca_arr[i]){continue;}
                    if(type==3) vec_m_sq_p_type_3_for_dcas[i]->Fill(m_squared_p);
                    if(type==31) vec_m_sq_p_type_31_for_dcas[i]->Fill(m_squared_p);
                }
            }

            for(int i=0;i<17;i++)
            {
                float r = 5+i*5;
                if(radius>r)
                {
                    if(type==3)
                    {
                        vec_m_squared_type_3_K[i]->Fill(m_squared_K);
                        vec_m_squared_type_3_p[i]->Fill(m_squared_p);
                        vec_m_squared_type_3_add_pi[i]->Fill(m_squared_add_pi);
                        vec_m_squared_type_3_K01[i]->Fill(m_squared_K01);
                        vec_m_squared_type_3_K02[i]->Fill(m_squared_K02);

                    }

                    if(type==31)
                    {
                        vec_m_squared_type_31_K[i]->Fill(m_squared_K);
                        vec_m_squared_type_31_p[i]->Fill(m_squared_p);
                        vec_m_squared_type_31_add_pi[i]->Fill(m_squared_add_pi);
                        vec_m_squared_type_31_K01[i]->Fill(m_squared_K01);
                        vec_m_squared_type_31_K02[i]->Fill(m_squared_K02);

                    }

                    if(overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                    {
                        int j=i;
                        j+=17;

                        if(type==3)
                        {
                            vec_m_squared_type_3_K[j]->Fill(m_squared_K);
                            vec_m_squared_type_3_p[j]->Fill(m_squared_p);
                            vec_m_squared_type_3_add_pi[j]->Fill(m_squared_add_pi);
                            vec_m_squared_type_3_K01[j]->Fill(m_squared_K01);
                            vec_m_squared_type_3_K02[j]->Fill(m_squared_K02);

                        }

                        if(type==31)
                        {
                            vec_m_squared_type_31_K[j]->Fill(m_squared_K);
                            vec_m_squared_type_31_p[j]->Fill(m_squared_p);
                            vec_m_squared_type_31_add_pi[j]->Fill(m_squared_add_pi);
                            vec_m_squared_type_31_K01[j]->Fill(m_squared_K01);
                            vec_m_squared_type_31_K02[j]->Fill(m_squared_K02);

                        }

                    }
                }

               

            }
            //------------------------------------------------------------------------------


            //do for all tracks m2 cut if available
            if ( !m2_cut_if_available(track_K, "K")) continue;
            if ( !m2_cut_if_available(track_p, "p")) continue;
            if ( !m2_cut_if_available(track_add_pi, "pi")) continue;
            if ( !m2_cut_if_available(track_K0_pi1, "pi")) continue;
            if ( !m2_cut_if_available(track_K0_pi2, "pi")) continue;


            TVector3 pos_K0 = DM->get_S2Vertex();
            TVector3 S_to_K0 = pos_K0-S_vertex_pos;
            TVector3 prim_to_K0 = pos_K0-pos_primary_vertex;

            double dist_K0_to_S = S_to_K0.Mag();
            double radiusK0 = prim_to_K0.Mag();

            hist_radiusK0->Fill(radiusK0);

            printtv3(vec_primtosec,"primtosec");
            printtv3(S_to_K0,"S_to_K0");

            double angle = vec_primtosec.Angle(S_to_K0);
            angle = angle/3.14159*180;

            histo_angle_ch3->Fill(angle);

            Ali_AS_V0* S_V0 = DM->getV0(0);
            Ali_AS_V0* K0_V0 = DM->getV0(1);

            double invmass_K0_by_tracks =  get_invariant_mass(track_K0_pi1, track_K0_pi2 , mass_pion,mass_pion , S_vertex_pos);

            TLorentzVector tlv_K0 = get_tlv_by_V0(mass_pion, mass_pion, K0_V0->getPpxpypz(), K0_V0->getNpxpypz() );
            printtlv(tlv_K0,"tlvK0");
            cout<<"angle: "<<angle<<endl;
            double invmass_K0 = tlv_K0.M();
            TVector3 dir_K0;
            dir_K0.SetXYZ(tlv_K0.Px(),tlv_K0.Py(),tlv_K0.Pz());

            //cut on m2 if information is available:
            

            if(type==3) histo_invariantmass_K0_type3->Fill(invmass_K0);
            if(type==31) histo_invariantmass_K0_type31->Fill(invmass_K0);

            //double m_squared_p = calculate_m_squared_by_TOF(track_p);
            //double m_squared_K = calculate_m_squared_by_TOF(track_K);

            //fill r-dependent histos-------------------------------------------------------

            for(int i=0;i<17;i++)
            {
                float min_rad = 5+i*5;

                if(radius>min_rad)
                {
                    if(type==3) vec_invmass_K0_type3[i]->Fill(invmass_K0);
                    if(type==31) vec_invmass_K0_type31[i]->Fill(invmass_K0);

                    
                }
            }

            //---------------------------------------------------------------------

            float dcaprim_pi1, dcaprim_pi2, dca_to_S_pi1, dca_to_S_pi2, dcaprim_add_pi, dcaprim_K, dcaprim_pro;

            FindDCAHelixPoint(pos_primary_vertex,track_K0_pi1,path_initA,path_initB,path_closest_to_point,dcaprim_pi1);
            FindDCAHelixPoint(pos_primary_vertex,track_K0_pi2,path_initA,path_initB,path_closest_to_point,dcaprim_pi2);
            FindDCAHelixPoint(pos_primary_vertex,track_add_pi,path_initA,path_initB,path_closest_to_point,dcaprim_add_pi);
            FindDCAHelixPoint(pos_primary_vertex,track_K,path_initA,path_initB,path_closest_to_point,dcaprim_K);
            FindDCAHelixPoint(pos_primary_vertex,track_p,path_initA,path_initB,path_closest_to_point,dcaprim_pro);

            FindDCAHelixPoint(S_vertex_pos,track_K0_pi1,path_initA,path_initB,path_closest_to_point,dca_to_S_pi1);
            FindDCAHelixPoint(S_vertex_pos,track_K0_pi1,path_initA,path_initB,path_closest_to_point,dca_to_S_pi2);

            //double dcaprim_K0 = calculateMinimumDistanceStraightToPoint(pos_K0,dir_K0,pos_primary_vertex);

            double dca_K0_prim  =  calculateMinimumDistanceStraightToPoint(pos_K0,dir_K0,pos_primary_vertex);

            float dcaprim_arr[6]={0.5,1,2,3,4,5};

            //loose cuts and save all with r>25
            if(dist_K0_to_S>0.5 && dcaprim_pi1>1. && dcaprim_pi2 > 1. && dca_to_S_pi1>0.1 && dca_to_S_pi2 >0.1)
            {
                //find out how often overlap
                if( dca_K0_prim > 0.5 && radius > 25)
                {
                    histo_counter_overlap->Fill(1.);
                    if( isoverlap(track_p, "p") ) histo_counter_overlap->Fill(2.);
                    if( isoverlap(track_K, "K") ) histo_counter_overlap->Fill(3.);

                }

                if( dca_K0_prim > 0.5 && radius > 25 && overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                {
                    AS_DM_particle->clearTrackList();
                    AS_DM_particle->clearV0List();
                    copy_dm_params(DM,AS_DM_particle);
                    AS_DM_particle->setN_V0s(7);
                    Tree_AS_DM_particle ->Fill();
                    counter_DM_saved++;

                }

            }



           // if(dist_K0_to_S>0.5 && dist_K0_to_S<20.  && dcaprim_pi1>2. && dcaprim_pi2 > 2. && dca_to_S_pi1>0.5 && dca_to_S_pi2 >0.5)
            if(dist_K0_to_S>0.5 && dcaprim_pi1>2. && dcaprim_pi2 > 2. && dca_to_S_pi1>0.2 && dca_to_S_pi2 >0.2)
             // && dcaprim_add_pi>0.5 && dcaprim_K > 0.5 && dcaprim_pro  > 0.5)
            {
                 for(int i=0;i<17;i++)
                 {
                     float min_rad = 5+i*5;
                     //double min_dist = 0.075*min_rad+0.125;
                     double min_dist = 0.075*radiusK0+0.125;

                     //vary dcaprim cut
                     if(radius>min_rad && overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                     {
                         for(int j=0;j<6;j++)
                         {
                             if(dca_K0_prim>dcaprim_arr[j])
                             {
                                 if(type==3) vec_vec_invmass_K0_type3[i][j]->Fill(invmass_K0);
                                 if(type==31) vec_vec_invmass_K0_type31[i][j]->Fill(invmass_K0);
                             }
                         }
                         
                     }

                     //radius-dependent backtracking cut
                     if(dca_K0_prim<min_dist){continue;}

                     if(radius>min_rad)
                     {
                         if(type==3) vec_invmass_K0_type3[i+17]->Fill(invmass_K0);
                         if(type==31) vec_invmass_K0_type31[i+17]->Fill(invmass_K0);

                         if(overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                         {
                             if(type==3) vec_invmass_K0_type3[i+17*2]->Fill(invmass_K0);
                             if(type==31) vec_invmass_K0_type31[i+17*2]->Fill(invmass_K0);

                             //3D for r_min=40
                             if(min_rad==40 && invmass_K0 > 0.488 && invmass_K0 < 0.515 && counter_DM_saved < 20000)
                             {
                                 /*
                                 AS_DM_particle->clearTrackList();
                                 AS_DM_particle->clearV0List();
                                 copy_dm_params(DM,AS_DM_particle);
                                 AS_DM_particle->setN_V0s(1);
                                 Tree_AS_DM_particle ->Fill();
                                 counter_DM_saved++;
                                 */
                             }


                             //----------------------------------------------------------------------
                             //normalband
                             if(min_rad==40 && invmass_K0 > 0.488 && invmass_K0 < 0.515 )
                             {
                                 if(type==3) vec_S_mass_ch3[0]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[0] ->Fill(S_mass);

                             }

                             //sideband1
                             if(min_rad==40 && invmass_K0 > 0.461 && invmass_K0 < 0.488 )
                             {
                                 if(type==3) vec_S_mass_ch3[1]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[1] ->Fill(S_mass);
                             }

                             //sideband2
                             if(min_rad==40 && invmass_K0 > 0.515 && invmass_K0 < 0.542)
                             {
                                 if(type==3) vec_S_mass_ch3[1]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[1] ->Fill(S_mass);
                             }

                         }
                     }

                     if(radiusK0>min_rad)
                     {
                         if(type==3) vec_invmass_K0_type3[i+17*3]->Fill(invmass_K0);
                         if(type==31) vec_invmass_K0_type31[i+17*3]->Fill(invmass_K0);

                     }


                 }

                
            }


            //dca combi 2-------------------------------------------------------------
           // if(dist_K0_to_S>1.5 && dist_K0_to_S<20.  && dcaprim_pi1>1.5 && dcaprim_pi2 > 1.5 &&                  //dcaprim auf 1.5, winkel 120 statt 90
             //  dcaprim_add_pi>1.5 && dca_to_S_pi1>0.5 && dca_to_S_pi2 >0.5
            //  && dcaprim_K > 1.5 && dcaprim_pro > 1.5 && fabs(angle)<120)

            if( dist_K0_to_S>0.5 && dist_K0_to_S<20.  && fabs(angle)<120
               && dcaprim_pi1>2. && dcaprim_pi2 > 2.  && dca_to_S_pi1>0.2 && dca_to_S_pi2 >0.2
               && dcaprim_add_pi>0.5 && dcaprim_K > 0.5 && dcaprim_pro  > 0.5 )
            {
                 for(int i=0;i<17;i++)
                 {
                     float min_rad = 5+i*5;
                     //double min_dist = 0.075*min_rad+0.125;
                     double min_dist = 0.075*radiusK0+0.125;
    
    
                     //vary dcaprim cut
                     if(radius>min_rad && overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                     {
                         for(int j=0;j<6;j++)
                         {
                             if(dca_K0_prim>dcaprim_arr[j])
                             {
                                 if(type==3) vec_vec_invmass_K0_type3[i+17][j]->Fill(invmass_K0);
                                 if(type==31) vec_vec_invmass_K0_type31[i+17][j]->Fill(invmass_K0);
                             }
                         }
    
                     }
    
                     //radius-dependent backtracking cut
                     if(dca_K0_prim<min_dist){continue;}
    
                     if(radius>min_rad)
                     {
                         if(type==3) vec_invmass_K0_type3[i+17*4]->Fill(invmass_K0);
                         if(type==31) vec_invmass_K0_type31[i+17*4]->Fill(invmass_K0);
    
                         if(overlap_PID(track_p, "p") && overlap_PID(track_K, "K"))
                         {
                             if(type==3) vec_invmass_K0_type3[i+17*5]->Fill(invmass_K0);
                             if(type==31) vec_invmass_K0_type31[i+17*5]->Fill(invmass_K0);


                             //-----------------------------------------------------------------------------
                             //complete radial dependend sidebandanalysis
                             double mean = 0.498227;
                             double sig = 0.004206;

                             //normalband
                             if(invmass_K0 > mean-2*sig && invmass_K0 < mean+2*sig)
                             {
                                 if(type==3) vec_S_mass_ch3_normalband[i]->Fill(S_mass);

                             }

                             //sideband
                             if
                             (
                              ( mean-4*sig < invmass_K0 && invmass_K0 < mean-2*sig ) ||
                              ( mean+2*sig < invmass_K0 && invmass_K0 < mean + 4*sig)
                             )
                             {
                                 if(type==3) vec_S_mass_ch3_sideband[i]->Fill(S_mass);
                             }



                             //-----------------------------------------------------------------------------

                             //3D for r_min=40
                             if(min_rad==40 && invmass_K0 > 0.488 && invmass_K0 < 0.515 && counter_DM_saved < 20000)
                             {
                                 /*
                                 AS_DM_particle->clearTrackList();
                                 AS_DM_particle->clearV0List();
                                 copy_dm_params(DM,AS_DM_particle);
                                 AS_DM_particle->setN_V0s(2);
                                 Tree_AS_DM_particle ->Fill();
                                 counter_DM_saved++;
                                 */
                             }
    
                             //normalband
                             if(min_rad==40 && invmass_K0 > 0.488 && invmass_K0 < 0.515 )
                             {
                                 if(type==3) vec_S_mass_ch3[0+2]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[0+2] ->Fill(S_mass);
    
                             }
    
                             //sideband1
                             if(min_rad==40 && invmass_K0 > 0.461 && invmass_K0 < 0.488 )
                             {
                                 if(type==3) vec_S_mass_ch3[1+2]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[1+2] ->Fill(S_mass);
                             }
    
                             //sideband2
                             if(min_rad==40 && invmass_K0 > 0.515 && invmass_K0 < 0.542)
                             {
                                 if(type==3) vec_S_mass_ch3[1+2]->Fill(S_mass);
                                 if(type==31) vec_S_mass_ch31[1+2] ->Fill(S_mass);
                             }
    
                         }
                     }
                 }


            }



            //----------------------------------------------------------------------------

            

            //-----------------------------------------------------------------------


            if(invmass_K0 < (0.4981+0.0042*2) && invmass_K0 > (0.4981-0.0042*2))
            {
                if(type==3) histo_counter_S->Fill(142);
                if(type==31) histo_counter_S->Fill(152);
                if( overlap_PID(track_K, "K") &&  overlap_PID(track_p, "p") )
                {
                    if(type==3)
                    {
                        histo_counter_S->Fill(143);
                        if(dcaprimcheck) histo_counter_S->Fill(144);
                        if(dcaprimcheck && radius>10) histo_counter_S->Fill(145);
                    }

                    if(type==31)
                    {
                        histo_counter_S->Fill(153);
                        if(dcaprimcheck) histo_counter_S->Fill(154);
                        if(dcaprimcheck && radius>10) histo_counter_S->Fill(155);
                    }
                }
            }

            if(invmass_K0 < (0.4981+0.0042*2) && invmass_K0 > (0.4981-0.0042*2))
            {
                histo_counter_S->Fill(14);

                if( overlap_PID(track_K, "K") &&  overlap_PID(track_p, "p") )
                {
                    histo_counter_S->Fill(16);
                    //if(m_squared_K>0.2 && m_squared_K < 0.35)
                    //{
                        histo_counter_S->Fill(17);

                        //dca prim just for add pion
                        bool dcaprimcheck_type3=1;
                        Ali_AS_Track* track = DM->getTrack(4);
                        float dcaprim;
                        FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
                        if(dcaprim<1.){dcaprimcheck_type3=0;}

                        if(dcaprimcheck_type3)
                        {
                            /*
                            AS_DM_particle->clearTrackList();
                            copy_dm_params(DM,AS_DM_particle);
                            Tree_AS_DM_particle ->Fill();
                            histo_counter_S->Fill(23);
                            */
                        }
                   // }
                }

                if(radius>20)
                {
                    histo_counter_S->Fill(15);
                }
            }

            if(radius>20)
            {
                histo_counter_S->Fill(5);
            }

        }

        if(type==3)
        {

            //dca prim >1
            //radius > 10

            Ali_AS_V0* S_V0 = DM->getV0(0);
            Ali_AS_V0* K0_V0 = DM->getV0(1);

            TVector3 K0_vertex = DM->get_S2Vertex();

            Ali_AS_Track* trackK0P = DM->getTrack(2);
            Ali_AS_Track* trackK0N = DM->getTrack(3);

            Ali_AS_Track* track_K_plus  = DM->getTrack(0);
            Ali_AS_Track* track_antip  = DM->getTrack(1);
            Ali_AS_Track* track_pion  = DM->getTrack(4);

            Float_t dca_prim_K,dca_prim_antip,dca_prim_pion;
            FindDCAHelixPoint(pos_primary_vertex,track_K_plus,path_initA,path_initB,path_closest_to_point,dca_prim_K);
            FindDCAHelixPoint(pos_primary_vertex,track_antip,path_initA,path_initB,path_closest_to_point,dca_prim_antip);
            FindDCAHelixPoint(pos_primary_vertex,track_pion,path_initA,path_initB,path_closest_to_point,dca_prim_pion);

            if(dca_prim_K >1. && dca_prim_antip >1. && dca_prim_pion >1.)
            {
                double invariantmass = get_invariant_mass( trackK0P, trackK0N, mass_pion,mass_pion,  K0_vertex);
                //histo_invariantmass_K0_type3 ->Fill(invariantmass);
    
                double m_squared_K_plus = calculate_m_squared_by_TOF(track_K_plus);
                double m_squared_antip = calculate_m_squared_by_TOF(track_antip);
                double m_squared_pion = calculate_m_squared_by_TOF(track_pion);
    
                mass_squared_kaons_type3->Fill(m_squared_K_plus);
                mass_squared_antip_type3->Fill(m_squared_antip);
    
                if(m_squared_K_plus>0.2 && m_squared_K_plus < 0.35 && m_squared_antip>0.6
                   && m_squared_antip <1.2)
                {
                    histo_invariantmass_K0_type3_with_cuts_on_antip_and_K->Fill(invariantmass);
    
                }
    
                if(m_squared_antip>0.6 && m_squared_antip <1.2)
                {
                    histo_invariantmass_K0_type3_with_cut_on_antip->Fill(invariantmass);
                    mass_squared_kaons_type3_with_cut_on_antip->Fill(m_squared_K_plus);
                }

                if(m_squared_K_plus>0.2 && m_squared_K_plus<0.35)
                {
                    mass_squared_antip_type3_with_cut_on_K->Fill(m_squared_antip);
                    histo_invariantmass_K0_type3_with_cuts_on_K->Fill(invariantmass);
                }

                if(invariantmass< (0.4981+0.0042*2) && invariantmass > (0.4981-0.0042*2))
                {
                    mass_squared_kaons_type3_with_cut_on_invmass_K0->Fill(m_squared_K_plus);
                    mass_squared_antip_type3_with_cut_on_invmass_K0->Fill(m_squared_antip);

                }

               // AS_DM_particle->clearTrackList();
                //copy_dm_params(DM,AS_DM_particle);
                //Tree_AS_DM_particle ->Fill();

            }







        }


        if(type==4 || type==41)
        {

            if(type==4) histo_counter_S->Fill(161);
            if(type==41) histo_counter_S->Fill(171);

            if(radius>20)
            {
                

            }
            bool dot_product_check=1;


            TVector3 prim_vertex = DM->get_primVertex();
            TVector3 pos_S_vertex = DM -> get_S1Vertex();
            TVector3 pos_K0_1 = DM -> get_S2Vertex();
            TVector3 pos_K0_2 = DM -> get_S3Vertex();

            TVector3 vec_prim_to_S;
            vec_prim_to_S =  pos_S_vertex- prim_vertex;
            TVector3 unit_vec_prim_to_S;
            unit_vec_prim_to_S = vec_prim_to_S.Unit();

            TLorentzVector tlv = DM->get_tlv();

            TVector3 dir_type4;
            TVector3 unit_dir_type4;
            dir_type4.SetXYZ(tlv[0],tlv[1],tlv[2]);
            unit_dir_type4 = dir_type4.Unit();

            if(unit_dir_type4.Dot(unit_vec_prim_to_S)<0.8){continue;}

            Ali_AS_Track* trackK0_1_1 = DM->getTrack(2);
            Ali_AS_Track* trackK0_1_2 = DM->getTrack(3);
            Ali_AS_Track* trackK0_2_1 = DM->getTrack(4);
            Ali_AS_Track* trackK0_2_2 = DM->getTrack(5);

            Ali_AS_Track* track_pi;
            Ali_AS_Track* track_p;
            if(type==4)
            {
                track_pi = DM->getTrack(0);
                track_p= DM->getTrack(1);
            }
            if(type==41)
            {
                track_pi = DM->getTrack(1);
                track_p= DM->getTrack(0);
            }
            double m_squared_pi_plus = calculate_m_squared_by_TOF(track_pi);
            double m_squared_antip = calculate_m_squared_by_TOF(track_p);
            double invm_K01  = get_invariant_mass( trackK0_1_1, trackK0_1_2, mass_pion,mass_pion,  pos_K0_1);
            double invm_K02  = get_invariant_mass( trackK0_2_1, trackK0_2_2, mass_pion,mass_pion,  pos_K0_2);


            if(invm_K01 < (0.4981+0.0042*2) &&  invm_K01 > (0.4981-0.0042*2)
               && invm_K02 < (0.4981+0.0042*2) &&  invm_K02 > (0.4981-0.0042*2))
            {
                 if(type==4) histo_counter_S->Fill(162);
                 if(type==41) histo_counter_S->Fill(172);

                 if( overlap_PID(track_p, "p")  )
                 {
                     if(type==4)
                     {
                         histo_counter_S->Fill(163);
                         if(dcaprimcheck) histo_counter_S->Fill(164);
                         if(dcaprimcheck && radius > 10) histo_counter_S->Fill(165);
                     }
                     if(type==41)
                     {
                         histo_counter_S->Fill(173);
                         if(dcaprimcheck) histo_counter_S->Fill(174);
                         if(dcaprimcheck && radius > 10) histo_counter_S->Fill(175);
                     }
                 }
                histo_counter_S->Fill(18);
                if(dcaprimcheck==1){histo_counter_S->Fill(24);}

                if(m_squared_antip>0.6 && m_squared_antip<1.2)
                {
                    /*
                    histo_counter_S->Fill(25);
                    AS_DM_particle->clearTrackList();
                    copy_dm_params(DM,AS_DM_particle);
                    Tree_AS_DM_particle ->Fill();
                    */
                }
            }




        }

        if(type==51)
        {
            histo_counter_type5->Fill(1+10);
            vector<Ali_AS_Track*> all_tracks_ch51;
            vector<int> all_track_ids_ch51;
            all_tracks_ch51.resize(5);
            all_track_ids_ch51.resize(5);

            for(int i=0;i<5;i++)
            {
                all_tracks_ch51[i]= DM->getTrack(i);
                all_track_ids_ch51[i] = all_tracks_ch51[i]->gettrackid();
            }


            if(check_if_two_vectors_have_same_element(all_track_ids_ch51,brute_force_type51)){continue;}

            for(int i=0;i<5;i++)
            {
                brute_force_type51.push_back(all_track_ids_ch51[i]);
            }


            int egal = DM_Analysis_type5(DM,2);

        }

        

        if(type==5)
        {
            
            cout<<"type5" <<endl;
            histo_counter_type5->Fill(1);

            //check for same track ids twice-------------------------------------------------
            vector<Ali_AS_Track*> all_tracks_ch5;
            vector<int> all_track_ids_ch5;
            all_tracks_ch5.resize(5);
            all_track_ids_ch5.resize(5);

            for(int i=0;i<5;i++)
            {
                all_tracks_ch5[i]= DM->getTrack(i);
                all_track_ids_ch5[i] = all_tracks_ch5[i]->gettrackid();
            }


            if(check_if_two_vectors_have_same_element(all_track_ids_ch5,brute_force_type5)){continue;}

           for(int i=0;i<5;i++)
           {
               brute_force_type5.push_back(all_track_ids_ch5[i]);
           }
            //------------------------------------------------------------------

            TVector3 S_vertex_pos = DM->get_S1Vertex();
            TVector3 Lambda_vertex_pos = DM->get_S2Vertex();
            TVector3 prim_to_S = S_vertex_pos-pos_primary_vertex;
            TVector3 prim_to_Lambda = Lambda_vertex_pos - pos_primary_vertex;

            

            //event mixing start------------------------------------------------------
            //---------------------------------------------------------------------
            double eta = prim_to_S.Eta();
            double phi = prim_to_S.Phi();

            histo_eta_type_5->Fill(eta);
            histo_phi_type_5->Fill(phi);

            int eta_cat = -1;
            int phi_cat = -1;

            //divide eta in 4 categories
            for(int i=0;i<4;i++)
            {
                if(eta>=-0.8+i*0.4 && eta < -0.4+i*0.4){eta_cat =i;}
            }
            printf("eta: %f, eta_cat: %d \n",eta,eta_cat);

            //divide phi in 10 categories
            double pi = 3.141592654;
            for(int i=0;i<10;i++)
            {
                if( phi>=-pi + i*2*pi/10. && phi< 0+(i+1)*2*pi/10. ){phi_cat = i;}
            }
            printf("phi: %f, phi_cat: %d \n",phi,phi_cat);

            Ali_AS_DM_particle DM_new;
            copy_dm_params(DM,&DM_new);
            UShort_t numtracks = DM->getNumTracks();
            for(int i=0;i<numtracks;i++)
            {
                Ali_AS_Track* track = DM->getTrack(i);
                DM_new.settrack(i,*track);
            }

            if(eta_cat!=-1 && phi_cat!=-1)
            //if(eta_cat==2 && phi_cat==5)
            {

                int size = mixed_event_vec[eta_cat][phi_cat].size();
                //cout<<"size: "<<size<<endl;
                if(size<5){mixed_event_vec[eta_cat][phi_cat].push_back(DM_new);}
                size = mixed_event_vec[eta_cat][phi_cat].size();
                cout<<"size: "<<size<<endl;
                //print_dm_event(DM);
                //cout<<""<<endl;
                //print_dm_event(&mixed_event_vec[eta_cat][phi_cat][size-1]);

                if(size==5)
                {
                    cout<<"before mixing"<<endl;
                    print_dm_event(DM);
                    cout<<""<<endl;
                    print_dm_event(&mixed_event_vec[eta_cat][phi_cat][size-2]);
                    cout<<""<<endl;
                    print_dm_event(&mixed_event_vec[eta_cat][phi_cat][size-1]);
                    cout<<""<<endl;

                    vector<Ali_AS_DM_particle> mixed_event = mix_event(mixed_event_vec,eta_cat,phi_cat);
                    for(int i=0;i<mixed_event.size();i++)
                    {
                        all_mixed_events.push_back(mixed_event[i]);
                    }

                    //mixed_event_vec[eta_cat][phi_cat].erase(mixed_event_vec[eta_cat][phi_cat].begin(),mixed_event_vec[eta_cat][phi_cat].end());
                    mixed_event_vec[eta_cat][phi_cat].clear();
                    mixed_event_vec[eta_cat][phi_cat].reserve(5);
                }
            }
            //end of mixing

            //cout<<""<<endl;
            //-------------------------------------------------------------------------
            //do analysis by function***********************************************
            int egal = DM_Analysis_type5(DM,0);
            //***********************************************************************+


            /*
            Ali_AS_Track* track_kaon = DM->getTrack(0);
            Ali_AS_Track* track_pi_minus = DM->getTrack(1);
            Ali_AS_Track* track_pi_plus = DM->getTrack(2);
            Ali_AS_Track* track_lambda_pi = DM->getTrack(3);
            Ali_AS_Track* track_lambda_antip = DM->getTrack(4);

            TLorentzVector tlv_lambda_pi = get_tlv(track_lambda_pi,mass_pion,Lambda_vertex_pos);
            TLorentzVector tlv_lambda_antip = get_tlv(track_lambda_antip,mass_proton,Lambda_vertex_pos);
            TLorentzVector tlv_lambda = tlv_lambda_pi + tlv_lambda_antip;
            TVector3 p_lambda = tlv_lambda.Vect();
            TVector3 p_lambda2 = DM->get_DirSV2();

            printf("x: %f, y: %f, z: %f \n",p_lambda[0],p_lambda[1],p_lambda[2]);
            printf("x: %f, y: %f, z: %f \n",p_lambda2[0],p_lambda2[1],p_lambda2[2]);
                                         
            //check dca to vertices
            bool dca_check=1;
            for(int i=0;i<5;i++)
            {
                Ali_AS_Track* track = all_tracks_ch5[i];
                float dca;
 
                if(i==0 || i==1 || i==2){ FindDCAHelixPoint(S_vertex_pos,track,path_initA,path_initB,path_closest_to_point,dca);}
                if(i==3 || i==4){ FindDCAHelixPoint(Lambda_vertex_pos,track,path_initA,path_initB,path_closest_to_point,dca);}
                if(dca>0.5){dca_check=0;}
            }
 
            if(!dca_check){continue;}

            histo_counter_type5->Fill(2);
 
            
            TVector3 S_to_L = Lambda_vertex_pos-S_vertex_pos;

            double dist_S_to_L = S_to_L.Mag();
            if(dist_S_to_L<1.){continue;}

            histo_counter_type5->Fill(3);
 
 
            double invmass_Lambda =  get_invariant_mass(track_lambda_pi, track_lambda_antip , mass_pion, mass_proton , Lambda_vertex_pos);
 
            double m_squared_kaon = calculate_m_squared_by_TOF(track_kaon);
            double m_squared_pi_minus = calculate_m_squared_by_TOF(track_pi_minus);
            double m_squared_pi_plus = calculate_m_squared_by_TOF(track_pi_plus);
            double m_squared_lambda_antip = calculate_m_squared_by_TOF(track_lambda_antip);
 
            histos_m_squared_type5[0]->Fill(m_squared_kaon);
            histos_m_squared_type5[1]->Fill(m_squared_pi_minus);
            histos_m_squared_type5[2]->Fill(m_squared_pi_plus);
            histos_m_squared_type5[3]->Fill(m_squared_lambda_antip);
 
            histo_invmass_Lambda_type5 -> Fill(invmass_Lambda);
 
            TLorentzVector tlv_proton;
            tlv_proton.SetPxPyPzE(0.,0.,0.,mass_proton);
 
            TLorentzVector tlv_type5 = DM->get_tlv();
 
            //cout<<"tlv: "<<tlv_type5[0]<<" "<<tlv_type5[1]<<" "<<tlv_type5[2]<<endl;
 
            tlv_type5-=tlv_proton;
 
            double S_mass_type5 = tlv_type5.M();
            histo_S_mass_type5->Fill(S_mass_type5);
 
            /*
            if(S_mass_type5<0.9)
            {
                cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                cout<<"tlv: "<<tlv_type5[0]<<" "<<tlv_type5[1]<<" "<<tlv_type5[2]<<"  "<<tlv_type5[3]<<endl;
            }

 
            if(0.6<m_squared_lambda_antip && m_squared_lambda_antip<1.2)
            {
                histo_S_mass_type5_cut_on_antip->Fill(S_mass_type5);
 
            }
 
           
            double sigma_kaon = track_kaon -> getnsigma_K_TPC();
            double sigma_pi_minus = track_pi_minus -> getnsigma_pi_TPC();
            double sigma_pi_plus = track_pi_plus -> getnsigma_pi_TPC();
            double sigma_lambda_antip = track_lambda_antip -> getnsigma_p_TPC();
            double sigma_lambda_pi = track_lambda_pi -> getnsigma_p_TPC();

            //show that if no anticuts apply, then no overlap
            Float_t TPCdEdx   = track_kaon->getTPCdEdx();
            Float_t tofsignal = track_kaon->getTOFsignal();
            Float_t dca       = track_kaon->getdca();
            int charge;
            if(dca>0){charge = 1;}
            else {charge = -1;}
            TLorentzVector tlv = track_kaon->get_TLV_part();
            double momentum = tlv.P();

            if(fabs(track_kaon->getnsigma_pi_TPC() )>2.5 &&  fabs(track_kaon->getnsigma_p_TPC() )>2.5 && fabs(track_kaon->getnsigma_e_TPC() )>2.5 )
            {
                dEdx_type5_anticuts->Fill(charge*momentum,TPCdEdx);
            }
            TPCdEdx   = track_lambda_antip->getTPCdEdx();
            tofsignal = track_lambda_antip->getTOFsignal();
            dca       = track_lambda_antip->getdca();
            if(dca>0){charge = 1;}
            else {charge = -1;}
            tlv = track_lambda_antip->get_TLV_part();
            momentum = tlv.P();
            if(fabs(track_lambda_antip->getnsigma_pi_TPC() )>2.5 && fabs(track_lambda_antip->getnsigma_K_TPC() )>2.5 && fabs(track_lambda_antip->getnsigma_e_TPC() )>2.5)
            {
                dEdx_type5_anticuts->Fill(charge*momentum*-1,TPCdEdx);
            }
            //---------------------------------------------------------------------------------------------
            //anticuts decide if TOF-hit is needed  //********************************anticuts******************************************
            if(fabs(track_kaon->getnsigma_pi_TPC() )<2.5 || fabs(track_kaon->getnsigma_p_TPC() )<2.5 || fabs(track_kaon->getnsigma_e_TPC() )<2.5 )
            {
                if(m_squared_kaon<0.2 || m_squared_kaon>0.35 ){continue;}
            }
 
            if(fabs(track_lambda_antip->getnsigma_pi_TPC() )<2.5 || fabs(track_lambda_antip->getnsigma_K_TPC() )<2.5 || fabs(track_lambda_antip->getnsigma_e_TPC() )<2.5)
            {
                if(m_squared_lambda_antip<0.6 || m_squared_lambda_antip>1.2 ){continue;}
 
            }
            histo_counter_type5->Fill(4);
 
            histo_S_mass_type5_anticuts->Fill(S_mass_type5);

            double angle = S_to_L.Angle(prim_to_S);
            if(angle<90*3.14159/180.){histo_S_mass_type5_anticuts_angle_90->Fill(S_mass_type5);}
            if(angle<120*3.14159/180.){histo_S_mass_type5_anticuts_angle_120->Fill(S_mass_type5);}
            if(angle<150*3.14159/180.){histo_S_mass_type5_anticuts_angle_150->Fill(S_mass_type5);}
 
            bool check_dcaprim=1;
            for(int i=0;i<5;i++)
            {
                Ali_AS_Track* track = all_tracks_ch5[i];
                float dcaprim;
                FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
                if(dcaprim<0.5){check_dcaprim=0;}
 
            }

            if(check_dcaprim)
            {
                histo_counter_type5->Fill(5);
                histo_S_mass_type5_dcaprim->Fill(S_mass_type5);
            }
            histo_counter_type5->Fill(6);

            if(p_lambda.Angle(prim_to_Lambda)<10 /180.*3.141){continue;}
            histo_counter_type5->Fill(7);
            histo_S_mass_type5_angle_lambda->Fill(S_mass_type5);
            if(check_dcaprim)histo_S_mass_type5_angle_lambda_and_dcaprim->Fill(S_mass_type5);

            //fill AS_DM_particle
            if(check_dcaprim && angle<120*3.14159/180. && radius > 10)
            {
                histo_counter_type5->Fill(8);
                histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10->Fill(S_mass_type5);
                AS_DM_particle->clearTrackList();
                copy_dm_params(DM,AS_DM_particle);
                Tree_AS_DM_particle ->Fill();
            }
 
 
            //FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
 
            */
 
            
 
 
 
 
 




        }

        




    }

    //




    //cout<<"b"<<endl;
    //loop over V0s
    /*
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
        //if(V0_is_used==1){continue;}

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

        //cout<<"tpccls: "<<ntpcclsP<<" "<<ntpcclsN<<endl;
        //cout<<"nitscls: "<<nitsclsP<<" "<<nitsclsN<<endl;
        //cout<<"status: "<<statusP<<" "<<statusN<<endl;
        //cout<<"tpcchi2: "<<tpcchi2P<<" "<<tpcchi2N<<endl;

        
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
*/

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

    

    tlv_neutron->SetPxPyPzE(0.,0.,0.,mass_neutron);

    TVector3 vec_primary_vertex_to_SV1;
    TVector3 unit_vec_primary_vertex_to_SV1;

    TVector3 momentum_SV1;
    TVector3 unit_momentum_SV1;

    TVector3 S_vertex_pos;
    /*
    float path_closest_to_point = 0;
    float dca_closest_to_point  = 0;
    float path_initA = 0.0;
    float path_initB = 30.0;
     */


    /*
    for(Int_t vector_loop_SV3 = 0; vector_loop_SV3 < (Int_t)vec_position_SV3.size(); vector_loop_SV3++)
    {
        for(Int_t vector_loop_SV2 = 0; vector_loop_SV2 < (Int_t)vec_position_SV2.size(); vector_loop_SV2++)
        {
            // printf("vector loop 1: %d, vector loop2: %d \n",vector_loop_SV3,vector_loop_SV2 ) ;
            if(vec_position_SV3.size() > 0 && vec_position_SV2.size() > 0 && vec_direction_SV3.size() > 0 &&  vec_direction_SV2.size() > 0)
            {
                if(calculateMinimumDistance(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]) > 0.5){ continue;}

                S_vertex_pos = calcVertexgAnalytical(vec_position_SV3[vector_loop_SV3],vec_direction_SV3[vector_loop_SV3],vec_position_SV2[vector_loop_SV2],vec_direction_SV2[vector_loop_SV2]);
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
    */

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

int Ali_Dark_Matter_Read::DM_Analysis_type5 (Ali_AS_DM_particle* DM,int mode)
{
    //mode 0: type5
    //mode 1: mixed event
    //mode 2: type 51

    int sideband=0;

    //sideband 0: analyse mean-2*sigma, mean+2*sigma for Lambda invariant mass
    //sideband 1: analyse mean-4*sigma, mean-2*sigma and mean+2*sigma, mean+4sigma for Lambda invariant mass

    if(mode==0) histo_counter_S->Fill(181);
    if(mode==2) histo_counter_S->Fill(191);

    
    //get all tracks and track ids------------------------------------------------------
    vector<Ali_AS_Track*> all_tracks_ch5;
    vector<Ali_AS_Track> all_tracks_ch5_no_ptr;
    vector<int> all_track_ids_ch5;
    all_tracks_ch5.resize(5);
    all_tracks_ch5_no_ptr.resize(5);
    all_track_ids_ch5.resize(5);

    if(mode==0 || mode ==2)
    {
        for(int i=0;i<5;i++)
        {
            all_tracks_ch5[i]= DM->getTrack(i);
            all_track_ids_ch5[i] = all_tracks_ch5[i]->gettrackid();
        }
    }

    if(mode==1)
    {
        for(int i=0;i<5;i++)
        {
            all_tracks_ch5_no_ptr[i]= DM->gettrack(i);
            all_tracks_ch5[i]=&all_tracks_ch5_no_ptr[i];
            all_track_ids_ch5[i] = all_tracks_ch5[i]->gettrackid();
        }
    }

    Ali_AS_Track* track_kaon;
    Ali_AS_Track* track_pi_minus;
    Ali_AS_Track* track_pi_plus;
    Ali_AS_Track* track_lambda_pi;
    Ali_AS_Track* track_lambda_antip;

    Ali_AS_Track track0;
    Ali_AS_Track track1;
    Ali_AS_Track track2;
    Ali_AS_Track track3;
    Ali_AS_Track track4;

    if(mode==0)
    {
        track_kaon = DM->getTrack(0);
        track_pi_minus = DM->getTrack(1);
        track_pi_plus = DM->getTrack(2);
        track_lambda_pi = DM->getTrack(3);
        track_lambda_antip = DM->getTrack(4);
    }

    if(mode==2)
    {
        track_kaon = DM->getTrack(1);
        track_pi_minus = DM->getTrack(2);
        track_pi_plus = DM->getTrack(0);
        track_lambda_pi = DM->getTrack(4);
        track_lambda_antip = DM->getTrack(3);

    }

    if(mode==1)
    {
        track0 = DM->gettrack(0);
        track1 = DM->gettrack(1);
        track2 = DM->gettrack(2);
        track3 = DM->gettrack(3);
        track4 = DM->gettrack(4);

        track_kaon=&track0;
        track_pi_minus=&track1;
        track_pi_plus=&track2;
        track_lambda_pi=&track3;
        track_lambda_antip=&track4;
    }
    //------------------------------------------------------------------


    //do m2 cut if available
    if ( !m2_cut_if_available(track_kaon, "K")) return 0;
    if ( !m2_cut_if_available(track_pi_minus, "pi")) return 0;
    if ( !m2_cut_if_available(track_pi_plus, "pi")) return 0;
    if ( !m2_cut_if_available(track_lambda_pi, "pi")) return 0;
    if ( !m2_cut_if_available(track_lambda_antip, "p")) return 0;


    //get TVector3s
    TVector3 pos_primary_vertex = DM->get_primVertex();
    TVector3 S_vertex_pos = DM->get_S1Vertex();
    TVector3 Lambda_vertex_pos = DM->get_S2Vertex();
    TVector3 original_Lambda_vertex_pos = DM->get_S3Vertex();
    TVector3 prim_to_S = S_vertex_pos-pos_primary_vertex;
    TVector3 prim_to_Lambda = Lambda_vertex_pos - pos_primary_vertex;
    TVector3 L_to_S = S_vertex_pos-Lambda_vertex_pos;
    TVector3 S_to_L = Lambda_vertex_pos-S_vertex_pos;
    double radius_L = prim_to_Lambda.Mag();
    TLorentzVector tlv_type5 = DM->get_tlv();
    double radius = prim_to_S.Mag();
    double dist_L_to_S = L_to_S.Mag();
    double angle_1 = prim_to_S.Angle(S_to_L);
    angle_1 = angle_1/3.14159*180;

    //cut on r!!!!!!!!!!!!!!!!!
    if(radius<5){return 0;}
    //----------------------------------------------------------------

    // get tlv by tracks
    TLorentzVector tlv_lambda_pi = get_tlv(track_lambda_pi,mass_pion,Lambda_vertex_pos);
    TLorentzVector tlv_lambda_antip = get_tlv(track_lambda_antip,mass_proton,Lambda_vertex_pos);
    TLorentzVector tlv_lambda = tlv_lambda_pi + tlv_lambda_antip;
    //------------------------------------------------

    //calc LorentzVector by V0 and invariantmass
    double invmass_antilambda_1;
    TLorentzVector tlv_L;
    if(mode==0 || mode==2)
    {
        Ali_AS_V0* L_V0 = DM->getV0(1);
        if(mode==0) tlv_L= get_tlv_by_V0(mass_pion, mass_proton, L_V0->getPpxpypz(), L_V0->getNpxpypz() );
        if(mode==2) tlv_L= get_tlv_by_V0(mass_proton, mass_pion, L_V0->getPpxpypz(), L_V0->getNpxpypz() );
        invmass_antilambda_1 = tlv_L.M();
    }

    if(mode==1)
    {
        invmass_antilambda_1 = tlv_lambda.M();
        tlv_L = tlv_lambda;
    }

    TVector3 p_lambda = tlv_L.Vect();
    TVector3 dir_L;
    dir_L.SetXYZ(tlv_L.Px(),tlv_L.Py(),tlv_L.Pz());

    TVector3 p_lambda2 = DM->get_DirSV2();
    //-------------------------------------------------------------------------------

    float radius_array[6] ={5,10,15,20,30,50};

    if(mode==0)
    {
        for(int i=0;i<6;i++)
        {
            if(radius>radius_array[i]) vec_invmass_Lambda_type5[i]->Fill(invmass_antilambda_1);
        }

        if( overlap_PID(track_lambda_antip, "p") )
        {
            for(int i=6;i<12;i++)
            {
                if(radius>radius_array[i-6]) vec_invmass_Lambda_type5[i]->Fill(invmass_antilambda_1);
            }

        }
    }

    if(mode==2)
    {
        for(int i=0;i<6;i++)
        {
            if(radius>radius_array[i]) vec_invmass_Lambda_type51[i]->Fill(invmass_antilambda_1);
        }

        if( overlap_PID(track_lambda_antip, "p") )
        {
            for(int i=6;i<12;i++)
            {
                if(radius>radius_array[i-6]) vec_invmass_Lambda_type51[i]->Fill(invmass_antilambda_1);
            }

        }
    }



    //look for lambdas and antilambdas at high radius
    Float_t path_initA = 0.0;
    Float_t path_initB = 30.0;
    Float_t path_closest_to_point = 0;
    Float_t dca_closest_to_point  = 0;
    float dcaprim_antip;
    float dcaprim_pi;
    float dcaprim_pi1;
    float dcaprim_pi2;
    float dcaprim_K;
    float dca_to_S_antip;
    float dca_to_S_pi;
    FindDCAHelixPoint(pos_primary_vertex,track_lambda_antip,path_initA,path_initB,path_closest_to_point,dcaprim_antip);
    FindDCAHelixPoint(pos_primary_vertex,track_lambda_pi,path_initA,path_initB,path_closest_to_point,dcaprim_pi);
    FindDCAHelixPoint(pos_primary_vertex,track_pi_plus,path_initA,path_initB,path_closest_to_point,dcaprim_pi1);
    FindDCAHelixPoint(pos_primary_vertex,track_pi_minus,path_initA,path_initB,path_closest_to_point,dcaprim_pi2);
    FindDCAHelixPoint(pos_primary_vertex,track_kaon,path_initA,path_initB,path_closest_to_point,dcaprim_K);

    FindDCAHelixPoint(S_vertex_pos,track_lambda_antip,path_initA,path_initB,path_closest_to_point,dca_to_S_antip);
    FindDCAHelixPoint(S_vertex_pos,track_lambda_pi,path_initA,path_initB,path_closest_to_point,dca_to_S_pi);

    double dca_L_prim  =  calculateMinimumDistanceStraightToPoint(Lambda_vertex_pos,dir_L,pos_primary_vertex);

    float dcaprim_arr[6]={0.5,1,2,3,4,5};

    TLorentzVector tlv_proton;
    tlv_proton.SetPxPyPzE(0.,0.,0.,mass_proton);

    //cout<<"tlv: "<<tlv_type5[0]<<" "<<tlv_type5[1]<<" "<<tlv_type5[2]<<endl;
    tlv_type5-=tlv_proton;

    double S_mass_type5 = tlv_type5.M();
    histo_S_mass_type5->Fill(S_mass_type5);



    //save all with loose cuts and r>25
    if(dist_L_to_S>1. && dcaprim_antip>1 && dcaprim_pi > 1. && dca_to_S_antip>0.1 && dca_to_S_pi >0.1 && overlap_PID(track_lambda_antip, "p"))
    {
        if(dca_L_prim>0.5  && radius>25)
        {
             AS_DM_particle->clearTrackList();
             AS_DM_particle->clearV0List();
             copy_dm_params(DM,AS_DM_particle);
             AS_DM_particle->setN_V0s(8);
             Tree_AS_DM_particle ->Fill();
        }
    }


    //if(dist_L_to_S>1. && dist_L_to_S<35. && dcaprim_antip>1. && dcaprim_pi > 1. && dca_to_S_antip>0.4 && dca_to_S_pi >0.4 && overlap_PID(track_lambda_antip, "p")  && fabs(angle_1)<120
    if(dist_L_to_S>1. && dist_L_to_S<35. && dcaprim_antip>1. && dcaprim_pi > 1. && dca_to_S_antip>0.4 && dca_to_S_pi >0.4   && fabs(angle_1)<120
      && dcaprim_pi1>0.5 && dcaprim_pi2 >0.5 && dcaprim_K>0.5 )
    {
       // float arr_r2[5]={15,20,25,30,50};
        for(int i=0;i<17;i++)
        {
            double r = 5+i*5;
            double min_dist = 0.075*radius_L+0.125;

            for(int j=0;j<6;j++)
            {
                if(radius<r){continue;}
                if(dca_L_prim<dcaprim_arr[j]){continue;}

                if(overlap_PID(track_kaon, "K") && overlap_PID(track_lambda_antip, "p") )
                {
                    if(mode==0) vec_vec_invmass_L_type5[i][j]->Fill(invmass_antilambda_1);
                    if(mode==2) vec_vec_invmass_L_type51[i][j]->Fill(invmass_antilambda_1);

                    if(mode==0 && dist_L_to_S>4.) vec_vec_invmass_L_type5[i+17][j]->Fill(invmass_antilambda_1);
                    if(mode==2 && dist_L_to_S>4.) vec_vec_invmass_L_type51[i+17][j]->Fill(invmass_antilambda_1);
                }


            }


            if(dca_L_prim<min_dist){continue;}

            if(radius>r && mode==0)  vec_invmass_Lambda_type5[12+i]->Fill(invmass_antilambda_1);
            if(radius>r && mode==2)  vec_invmass_Lambda_type51[12+i]->Fill(invmass_antilambda_1);

            if(radius_L>r && mode==0)  vec_invmass_Lambda_type5[12+17*2+i]->Fill(invmass_antilambda_1);
            if(radius_L>r && mode==2)  vec_invmass_Lambda_type51[12+17*2+i]->Fill(invmass_antilambda_1);

            

            if(overlap_PID(track_kaon, "K") && overlap_PID(track_lambda_antip, "p") )
            {
                 if(r==40 && counter_DM_saved<20000 && radius>r)
                 {
                     /*
                     AS_DM_particle->clearTrackList();
                     AS_DM_particle->clearV0List();
                     copy_dm_params(DM,AS_DM_particle);
                     AS_DM_particle->setN_V0s(1);
                     Tree_AS_DM_particle ->Fill();
                     counter_DM_saved++;
                     */
                 }

                 if(radius>r && mode==0)  vec_invmass_Lambda_type5[12+17+i]->Fill(invmass_antilambda_1);
                 if(radius>r && mode==2)  vec_invmass_Lambda_type51[12+17+i]->Fill(invmass_antilambda_1);

                 if(radius_L>r && mode==0)  vec_invmass_Lambda_type5[12+17*3+i]->Fill(invmass_antilambda_1);
                 if(radius_L>r && mode==2)  vec_invmass_Lambda_type51[12+17*3+i]->Fill(invmass_antilambda_1);

                 //----------------------------------------------------------------------------------------------
                 //do sideband analysis with given dca cuts


                 //normalband
                 if( invmass_antilambda_1 > (1.115718-0.001484*2) && invmass_antilambda_1 < (1.115718+0.001484*2))
                 {
                     if(radius>r && mode==0)
                     {
                         vec_S_mass_ch5_normalband[i]->Fill(S_mass_type5);

                     }

                 }

                 //sideband
                 if(
                    ( invmass_antilambda_1 > (1.115718-0.001484*4) && invmass_antilambda_1 < (1.115718-0.001484*2) )
                    ||  ( invmass_antilambda_1 > (1.115718+0.001484*2) && invmass_antilambda_1 < (1.115718+0.001484*4) )
                   )

                 {
                     if(radius>r && mode==0)
                     {
                         vec_S_mass_ch5_sideband[i]->Fill(S_mass_type5);
                     }

                 }


                 //----------------------------------------------------------------------------------------------



                 if(radius>r && dist_L_to_S>4. && dcaprim_pi1>1.5 && dcaprim_pi2>1.5 && dcaprim_K>1.5 && fabs(angle_1)<120)
                 {
                     if(mode==0) vec_invmass_Lambda_type5[12+17*4+i]->Fill(invmass_antilambda_1);
                     if(mode==2) vec_invmass_Lambda_type51[12+17*4+i]->Fill(invmass_antilambda_1);
                     if(r==40 && counter_DM_saved<20000)
                     {
                         /*
                         AS_DM_particle->clearTrackList();
                         AS_DM_particle->clearV0List();
                         copy_dm_params(DM,AS_DM_particle);
                         AS_DM_particle->setN_V0s(2);
                         Tree_AS_DM_particle ->Fill();
                         counter_DM_saved++;
                         */
                     }
                 }

            }
        }



        

        

    }

    //important sideband or not!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(sideband==0)
    {
        //2-sigma cut----------------------------------------------------------------------
        if( invmass_antilambda_1 > (1.1157+0.001495*2) || invmass_antilambda_1 < (1.1157-0.001495*2))
        {
            return 0;
        }
    }

    if(sideband==1)
    {
        //sideband cut----------------------------------------------------------------------
        if( invmass_antilambda_1< (1.1157-0.001495*4)){return 0;}
        if( invmass_antilambda_1> (1.1157+0.001495*4)){return 0;}
        if( invmass_antilambda_1> (1.1157-0.001495*2) && invmass_antilambda_1 < (1.1157+0.001495*2) ){return 0;}
        
    }


    if(mode==0) histo_counter_S->Fill(182);
    if(mode==2) histo_counter_S->Fill(192);
    //------------------------------------------------------------------------------------------


    bool check_dcaprim=1;
    bool dca_check=1;
    
    for(int i=0;i<5;i++)
    {
        Ali_AS_Track* track = all_tracks_ch5[i];
        float dcaprim;
        FindDCAHelixPoint(pos_primary_vertex,track,path_initA,path_initB,path_closest_to_point,dcaprim);
        if(dcaprim<0.5){check_dcaprim=0;}

    }


    //********************************************************************************

    //do topology cuts from grid!! because not done for mixed
    //TVector3 dir_L = DM->get_DirSV2();
    //TVector3 dir_L = p_lambda;
    double dist = calculateMinimumDistanceStraightToPoint(Lambda_vertex_pos,dir_L,S_vertex_pos);
    if(dist>0.5){return 0;}

    TVector3 vec_S_to_L = Lambda_vertex_pos-S_vertex_pos;
    TVector3 unit_vec_S_to_L;
    unit_vec_S_to_L = vec_S_to_L.Unit();
    TVector3 unit_dir_L = dir_L.Unit();
    if(unit_vec_S_to_L.Dot(unit_dir_L)<0.8){return 0;}

    TVector3 dir_type5;
    TVector3 unit_dir_type5;
    dir_type5.SetXYZ(tlv_type5[0],tlv_type5[1],tlv_type5[2]);
    unit_dir_type5 = dir_type5.Unit();
    TVector3 unit_vec_prim_to_S;
    unit_vec_prim_to_S = prim_to_S.Unit();
    if(unit_vec_prim_to_S.Dot(unit_dir_type5)<0.8){return 0;}
    //----------------------------------------------------------------------------------------
    if(mode==0) histo_counter_S->Fill(183);
    if(mode==2) histo_counter_S->Fill(193);
                                 
    //check dca to vertices
    for(int i=0;i<5;i++)
    {
        Ali_AS_Track* track = all_tracks_ch5[i];
        float dca;
        if(i==0 || i==1 || i==2){ FindDCAHelixPoint(S_vertex_pos,track,path_initA,path_initB,path_closest_to_point,dca);}
        if(mode==0 || mode==2){if(i==3 || i==4){ FindDCAHelixPoint(Lambda_vertex_pos,track,path_initA,path_initB,path_closest_to_point,dca);}}
        if(mode==1){if(i==3 || i==4){ FindDCAHelixPoint(original_Lambda_vertex_pos,track,path_initA,path_initB,path_closest_to_point,dca);}}
        //cout<<"dca: "<<dca<<endl;
        if(dca>0.5){dca_check=0;}
    }

    if(!dca_check){return 0;}

    if(mode==0) histo_counter_S->Fill(184);
    if(mode==2) histo_counter_S->Fill(194);

    if(mode == 0 ) histo_counter_type5->Fill(2);
    if(mode == 2 ) histo_counter_type5->Fill(2+10);

    
    //TVector3 S_to_L = Lambda_vertex_pos-S_vertex_pos;

    //double dist_S_to_L = S_to_L.Mag();
    //if(dist_S_to_L<1.){return 0;;}

    //histo_counter_type5->Fill(3);

    //fill various---------------------------------------------------------------------------------------------------------
    double invmass_Lambda =  get_invariant_mass(track_lambda_pi, track_lambda_antip , mass_pion, mass_proton , Lambda_vertex_pos);
                   
    double m_squared_kaon = calculate_m_squared_by_TOF(track_kaon);
    double m_squared_pi_minus = calculate_m_squared_by_TOF(track_pi_minus);
    double m_squared_pi_plus = calculate_m_squared_by_TOF(track_pi_plus);
    double m_squared_lambda_antip = calculate_m_squared_by_TOF(track_lambda_antip);

    if(mode==0)
    {
        histos_m_squared_type5[0]->Fill(m_squared_kaon);
        histos_m_squared_type5[1]->Fill(m_squared_pi_minus);
        histos_m_squared_type5[2]->Fill(m_squared_pi_plus);
        histos_m_squared_type5[3]->Fill(m_squared_lambda_antip);
        if( overlap_PID(track_kaon, "K"))
        {
            histos_m_squared_type5[0+4]->Fill(m_squared_kaon);
            histos_m_squared_type5[1+4]->Fill(m_squared_pi_minus);
            histos_m_squared_type5[2+4]->Fill(m_squared_pi_plus);
            histos_m_squared_type5[3+4]->Fill(m_squared_lambda_antip);
        }

        for(int i=0;i<17;i++)
        {
            float r = 5+i*5;
            if(radius>r)
            {
                m_squared_type5_K_r[i]->Fill(m_squared_kaon);
                m_squared_type5_pi_minus_r[i]->Fill(m_squared_pi_minus);
                m_squared_type5_pi_plus_r[i]->Fill(m_squared_pi_plus);
                m_squared_type5_p_r[i]->Fill(m_squared_lambda_antip);

                if( overlap_PID(track_kaon, "K"))
                {
                    m_squared_type5_K_r[i+17]->Fill(m_squared_kaon);
                    m_squared_type5_pi_minus_r[i+17]->Fill(m_squared_pi_minus);
                    m_squared_type5_pi_plus_r[i+17]->Fill(m_squared_pi_plus);
                    m_squared_type5_p_r[i+17]->Fill(m_squared_lambda_antip);
                }
            }
        }
    }

    if(mode==2)
    {
        histos_m_squared_type51[0]->Fill(m_squared_kaon);
        histos_m_squared_type51[1]->Fill(m_squared_pi_minus);
        histos_m_squared_type51[2]->Fill(m_squared_pi_plus);
        histos_m_squared_type51[3]->Fill(m_squared_lambda_antip);
        if( overlap_PID(track_kaon, "K"))
        {
            histos_m_squared_type51[0+4]->Fill(m_squared_kaon);
            histos_m_squared_type51[1+4]->Fill(m_squared_pi_minus);
            histos_m_squared_type51[2+4]->Fill(m_squared_pi_plus);
            histos_m_squared_type51[3+4]->Fill(m_squared_lambda_antip);
        }

        for(int i=0;i<17;i++)
        {
            float r = 5+i*5;
            if(radius>r)
            {
                m_squared_type51_K_r[i]->Fill(m_squared_kaon);
                m_squared_type51_pi_minus_r[i]->Fill(m_squared_pi_minus);
                m_squared_type51_pi_plus_r[i]->Fill(m_squared_pi_plus);
                m_squared_type51_p_r[i]->Fill(m_squared_lambda_antip);

                if( overlap_PID(track_kaon, "K"))
                {
                    m_squared_type51_K_r[i+17]->Fill(m_squared_kaon);
                    m_squared_type51_pi_minus_r[i+17]->Fill(m_squared_pi_minus);
                    m_squared_type51_pi_plus_r[i+17]->Fill(m_squared_pi_plus);
                    m_squared_type51_p_r[i+17]->Fill(m_squared_lambda_antip);
                }
            }
        }

    }
    /*
    if(mode==1)
    {
        histos_m_squared_type5[0+4]->Fill(m_squared_kaon);
        histos_m_squared_type5[1+4]->Fill(m_squared_pi_minus);
        histos_m_squared_type5[2+4]->Fill(m_squared_pi_plus);
        histos_m_squared_type5[3+4]->Fill(m_squared_lambda_antip);
    }
    */
    if(mode == 0) histo_invmass_Lambda_type5 -> Fill(invmass_Lambda);
    if(mode == 0) histo_invmass_Lambda_by_V0_type5 -> Fill(invmass_antilambda_1);
    if(mode == 2) histo_invmass_Lambda_type51 -> Fill(invmass_Lambda);

    double dist_L_prim = calculateMinimumDistanceStraightToPoint(Lambda_vertex_pos,dir_L,pos_primary_vertex);
    if(dist_L_prim>1.)
    {
        if(mode == 0) histo_invmass_Lambda_dcaprim_line_type5 -> Fill(invmass_antilambda_1);
        if(mode == 2) histo_invmass_Lambda_dcaprim_type51 -> Fill(invmass_antilambda_1);
    }

    if(mode==0 && check_dcaprim) histo_invmass_Lambda_dcaprim1_type5->Fill(invmass_antilambda_1);

    if(0.6<m_squared_lambda_antip && m_squared_lambda_antip<1.2)
    {
        if(mode == 0) histo_invmass_Lambda_antip_cut_type5 -> Fill(invmass_antilambda_1);
        if(mode == 2) histo_invmass_Lambda_antip_cut_type51 -> Fill(invmass_antilambda_1);
        if(mode == 2) histo_invmass_Lambda_antip_cut_byV0_type51 -> Fill(invmass_antilambda_1);
    }

    if(radius>10)
    {
        if(mode == 0) histo_invmass_Lambda_radius_cut_type5 -> Fill(invmass_antilambda_1);
        if(mode == 2) histo_invmass_Lambda_radius_cut_type51 -> Fill(invmass_antilambda_1);
    }

    //--------------------------------------------------------------------------------------------------------------------------

    


    if(0.6<m_squared_lambda_antip && m_squared_lambda_antip<1.2)
    {
        histo_S_mass_type5_cut_on_antip->Fill(S_mass_type5);

    }

   
    double sigma_kaon = track_kaon -> getnsigma_K_TPC();
    double sigma_pi_minus = track_pi_minus -> getnsigma_pi_TPC();
    double sigma_pi_plus = track_pi_plus -> getnsigma_pi_TPC();
    double sigma_lambda_antip = track_lambda_antip -> getnsigma_p_TPC();
    double sigma_lambda_pi = track_lambda_pi -> getnsigma_p_TPC();

    /*
    //show that if no anticuts apply, then no overlap
    Float_t TPCdEdx   = track_kaon->getTPCdEdx();
    Float_t tofsignal = track_kaon->getTOFsignal();
    Float_t dca       = track_kaon->getdca();
    int charge;
    if(dca>0){charge = 1;}
    else {charge = -1;}
    TLorentzVector tlv = track_kaon->get_TLV_part();
    double momentum = tlv.P();

    if(fabs(track_kaon->getnsigma_pi_TPC() )>2.5 &&  fabs(track_kaon->getnsigma_p_TPC() )>2.5 && fabs(track_kaon->getnsigma_e_TPC() )>2.5 )
    {
        dEdx_type5_anticuts->Fill(charge*momentum,TPCdEdx);
    }
    TPCdEdx   = track_lambda_antip->getTPCdEdx();
    tofsignal = track_lambda_antip->getTOFsignal();
    dca       = track_lambda_antip->getdca();
    if(dca>0){charge = 1;}
    else {charge = -1;}
    tlv = track_lambda_antip->get_TLV_part();
    momentum = tlv.P();
    if(fabs(track_lambda_antip->getnsigma_pi_TPC() )>2.5 && fabs(track_lambda_antip->getnsigma_K_TPC() )>2.5 && fabs(track_lambda_antip->getnsigma_e_TPC() )>2.5)
    {
        dEdx_type5_anticuts->Fill(charge*momentum*-1,TPCdEdx);
    }
        */
    //---------------------------------------------------------------------------------------------
    //anticuts decide if TOF-hit is needed  //********************************anticuts******************************************
    if(fabs(track_kaon->getnsigma_pi_TPC() )<2.5 || fabs(track_kaon->getnsigma_p_TPC() )<2.5 || fabs(track_kaon->getnsigma_e_TPC() )<2.5 )
    {
        if(m_squared_kaon<0.2 || m_squared_kaon>0.35 ){return 0;}
    }

    if(fabs(track_lambda_antip->getnsigma_pi_TPC() )<2.5 || fabs(track_lambda_antip->getnsigma_K_TPC() )<2.5 || fabs(track_lambda_antip->getnsigma_e_TPC() )<2.5)
    {
        if(m_squared_lambda_antip<0.6 || m_squared_lambda_antip>1.2 ){return 0;}

    }
    //anticuts/overlap done-------------------------------------------------------------------


    if(mode==0) histo_counter_S->Fill(185);
    if(mode==2) histo_counter_S->Fill(195);
    if(mode==0 && check_dcaprim) histo_counter_S->Fill(186);
    if(mode==2 && check_dcaprim) histo_counter_S->Fill(196);
    if(mode==0 && check_dcaprim && radius > 10) histo_counter_S->Fill(187);
    if(mode==2 && check_dcaprim && radius > 10) histo_counter_S->Fill(197);

    if(mode == 0) histo_invmass_Lambda_anti_cut_type5 -> Fill(invmass_antilambda_1);
    if(mode == 2) histo_invmass_Lambda_anti_cut_type51 -> Fill(invmass_antilambda_1);


    if(mode==0) histo_counter_type5->Fill(4);
    if(mode==2) histo_counter_type5->Fill(4+10);

    double angle_prim_L = p_lambda.Angle(prim_to_Lambda); //in radiant

    if(dist_L_prim>1. && angle_prim_L>10 /180.*3.141)
    {
        for(int i=0;i<6;i++)
        {
            if(mode==0) if(radius>radius_array[i]) vec_invmass_S_type5[i]->Fill(S_mass_type5);
            if(mode==2) if(radius>radius_array[i]) vec_invmass_S_type51[i]->Fill(S_mass_type5);
        }

    }

    histo_S_mass_type5_anticuts->Fill(S_mass_type5);

    double angle = S_to_L.Angle(prim_to_S);
    if(mode==0)
    {
        if(angle<90*3.14159/180.){histo_S_mass_type5_anticuts_angle_90->Fill(S_mass_type5);}
        if(angle<120*3.14159/180.){histo_S_mass_type5_anticuts_angle_120->Fill(S_mass_type5);}
        if(angle<150*3.14159/180.){histo_S_mass_type5_anticuts_angle_150->Fill(S_mass_type5);}
    }

    if(check_dcaprim)
    {
        if(mode==0) histo_counter_type5->Fill(5);
        if(mode==2) histo_counter_type5->Fill(5+10);
        histo_S_mass_type5_dcaprim->Fill(S_mass_type5);
    }

    if(check_dcaprim)
    {
        if(mode==0)
        {

        }

    }

    



    //histo_counter_type5->Fill(6);
    

    if(p_lambda.Angle(prim_to_Lambda)<10 /180.*3.141){return 0;}
    if(mode==0) histo_counter_type5->Fill(7);
    if(mode==2) histo_counter_type5->Fill(7+10);

    if(mode==0) histo_S_mass_type5_angle_lambda->Fill(S_mass_type5);
    if(mode==1) histo_S_mass_mixed_event_no_dcaprim->Fill(S_mass_type5);

    if(check_dcaprim)
    {
        if(mode==0) histo_S_mass_type5_angle_lambda_and_dcaprim->Fill(S_mass_type5);
        if(mode==1) histo_S_mass_mixed_event_with_dcaprim->Fill(S_mass_type5);
        if(mode==2) histo_S_mass_type51_angle_lambda_and_dcaprim->Fill(S_mass_type5);
    }

    //fill AS_DM_particle

    if(check_dcaprim && angle<120*3.14159/180. && radius > 7)
    {
        histo_counter_type5->Fill(8);
        histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10->Fill(S_mass_type5);
        
    }

    if(check_dcaprim && radius > 7)
    {
        if(mode == 0) histo_invmass_Lambda_final_cut_type5 -> Fill(invmass_antilambda_1);
        if(mode == 2) histo_invmass_Lambda_final_cut_type51 -> Fill(invmass_antilambda_1);
    }

    //AS_DM_particle->clearTrackList();
    //copy_dm_params(DM,AS_DM_particle);
    //Tree_AS_DM_particle ->Fill();
    //Tree_AS_DM_particle->Write();


    //fill brute force!   do it here or before??
    for(int i=0;i<5;i++)
    {
       // if(mode==0) brute_force_type5.push_back(all_track_ids_ch5[i]);
        //if(mode==2) brute_force_type51.push_back(all_track_ids_ch5[i]);
    }

    return 0;
}





vector<Ali_AS_DM_particle> Ali_Dark_Matter_Read::mix_event(vector<vector<vector<Ali_AS_DM_particle>>> &vec,int eta_cat,int phi_cat)
{
    cout<<"********************************************"<<endl;
    cout<<"mix started"<<endl;
    cout<<"eta: "<<eta_cat<<endl;
    cout<<"phi: "<<phi_cat<<endl;
    vector<Ali_AS_Track*> vec_tracks;
    Ali_AS_DM_particle* DMa = new Ali_AS_DM_particle();
    Ali_AS_DM_particle* DMb = new Ali_AS_DM_particle();

    Ali_AS_DM_particle* DM_mixed = new Ali_AS_DM_particle();

    vector<Ali_AS_DM_particle> mixed_events;

    for(int i=0;i<5;i++)
    {
      print_dm_event(&vec[eta_cat][phi_cat][i]);
    }
    cout<<"end of print"<<endl;


    for(int i=0;i<5;i++)
    {
        DMa = &vec[eta_cat][phi_cat][i];
        TVector3 prima = DMa->get_primVertex();
        TVector3 S_vertexa = DMa->get_S1Vertex();
        TVector3 prim_to_Sa = S_vertexa-prima;
        
        for(int j=0;j<5;j++)
        {
            if(i==j){continue;}
            DMb = &vec[eta_cat][phi_cat][j];
            TVector3 primb = DMb->get_primVertex();
            TVector3 S_vertexb = DMb->get_S1Vertex();
            TVector3 Lambda_b = DMb->get_S2Vertex();
            TVector3 S_a_to_S_b = S_vertexb - S_vertexa;

            
            DM_mixed -> settype(6);
            DM_mixed -> set_primVertex(prima);
            DM_mixed -> set_S1Vertex(S_vertexa);
            //verschiebe Lambda Vertex
            DM_mixed -> set_S2Vertex(Lambda_b - S_a_to_S_b);

            //save original Lambda vertex
            DM_mixed -> set_S3Vertex(Lambda_b);

            int numtracks = DM_mixed->getNumTracks();

            for(int i_track =0;i_track<3;i_track++)
            {
                //cout<<"i_track: "<<i_track<<endl;
                //cout<<"numtracks: "<<DMa->getNumTracks()<<endl;
                Ali_AS_Track tracka = DMa->gettrack(i_track);
                Ali_AS_Track* track_mixed = new Ali_AS_Track();
                copy_track_params(&tracka,track_mixed);
                DM_mixed->settrack(i_track,*track_mixed);
            }

            //numtracks = DM_mixed->getNumTracks();

            //helix noch nicht verschoben!
            for(int i_track =3;i_track<5;i_track++)
            {
                Ali_AS_Track trackb = DMb->gettrack(i_track);
                Ali_AS_Track* track_mixed = new Ali_AS_Track();
                copy_track_params(&trackb,track_mixed);
                DM_mixed->settrack(i_track,*track_mixed);
            }

           // numtracks = DM_mixed->getNumTracks();

            double masses[5] = {mass_K,mass_pion,mass_pion,mass_pion,mass_proton};
            TLorentzVector tlv_sum;
            tlv_sum.SetPxPyPzE(0.,0.,0,0.);
            for(int i_track =0;i_track<5;i_track++)
            {
                Ali_AS_Track track;
                TLorentzVector tlv;
                if(i_track<3)
                {
                    track= DMa->gettrack(i_track);
                    tlv = get_tlv(&track, masses[i], S_vertexa);
                }
                else
                {
                    track=DMb->gettrack(i_track);
                    tlv = get_tlv(&track, masses[i], Lambda_b - S_a_to_S_b);
                }

                tlv_sum+=tlv;
            }
            
            DM_mixed->set_tlv(tlv_sum);
            DM_mixed->clearTrackList();
            mixed_events.push_back(*DM_mixed);
            cout<<"DM mixed: "<<endl;
            print_dm_event(DM_mixed);
            

        }

    }
    return mixed_events;
    cout<<"end of mix"<<endl;



}





void Ali_Dark_Matter_Read::Analyse_Mixed_Events()
{
    int size = all_mixed_events.size();
    cout<<"size: "<<size<<endl;
    for(int i=0;i<size;i++)
    {
        DM_Analysis_type5(&all_mixed_events[i],1);
    }
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

    /*
    TCanvas* can10 = new TCanvas;
    TCanvas* can11 = new TCanvas;
    TCanvas* can12 = new TCanvas;
    TCanvas* can13 = new TCanvas;
    TCanvas* can14 = new TCanvas;
    TCanvas* can15 = new TCanvas;
    TCanvas* can16 = new TCanvas;
    TCanvas* can17 = new TCanvas;
    TCanvas* can18 = new TCanvas;
    */

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
    /*
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
    /*
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
    printf("nuclear events: %d \n",counter_NUCLEV);

    /*
    can6->cd();
    histo_radius_nuclev_3_or_more->Draw();
    can7->cd();
    histo_radius_nuclev_4_or_more->Draw();
    can8->cd();
    histo_radius_nuclev_5_or_more->Draw();
    can9->cd();
    histo_radius_nuclev_2_or_more->Draw();
    */

    //sum up one histo
    //TAxis* axis = histo_radius_nuclev_4_or_more->GetXaxis();
    int nbins = histo_radius_nuclev_4_or_more->GetNbinsX();
    cout<<"nbins: "<<nbins<<endl;

    for(int bin = 1;bin<nbins+1;bin++)
    {
        double integral = histo_radius_nuclev_4_or_more->Integral(1,bin);
        //cout<<"integral: "<<integral<<endl;
        histo_radius_nuclev_4_or_more_sum->SetBinContent(bin,integral);
    }

    can10->cd();
    histo_radius_nuclev_4_or_more_sum->Draw();
    

    /*
    can10->cd();
    x_and_y_nuclev_3->Draw("colz");
    can11->cd();
    x_and_y_nuclev_4->Draw("colz");
    can12->cd();
    x_and_y_nuclev_5->Draw("colz");
    can13->cd();
    x_and_y_nuclev_6->Draw("colz");
    */

    //for(int i=0;i<12;i++)
    //{
    int i =10;
    TCanvas* cani = new TCanvas();
    cani->cd();
    vec_radius_ortho_vs_z[i]->Draw("colz");
    TString name ="/home/ceres/schlichtmann/ESD_Analysis/Plots/";
    name+=vec_radius_ortho_vs_z[i]->GetTitle();
    name+=".png";
    cani->SaveAs(name.Data());

    TCanvas* cani2 = new TCanvas();
    cani2->cd();
    histo_radius_nuclev_4_or_more_momentum->Draw();

    //************************************************************************************
    //Write outputfile--------------------------------------------------------------------------------
    output_histos->cd();

    Tree_AS_DM_particle->Write();

    if (histo_numberDMs->GetEntries()>0) {histo_numberDMs->Write();}
    histo_counter_S->Write();
    histo_type_of_S->Write();
    TOF_eff_K_ch3->Write();
    histo_counter_overlap->Write();
    histo_n_tracks->Write();
    

    //nuclev----------------------------------------------
    output_histos ->mkdir("nuclev");
    output_histos ->cd("nuclev");
    histo_angle->Write();
    histo_angle_z_0_10->Write();
    histo_angle_z_50_60->Write();

    vec_radius_ortho_vs_z[10]->Write();
    for(int i=15;i<=21;i++)
    {
        vec_radius_ortho_vs_z[i]->Write();
    }

    vec_radius_ortho_vs_z[1]->Write();
    vec_radius_ortho_vs_z[5]->Write();
    vec_radius_ortho_vs_z[9]->Write();

    for(int i=0;i<7;i++)
    {
        vec_histo_radius_nuclev[i]->Write();
    }
    for(int i=0;i<10;i++)
    {
        vec_x_y_slices[i]->Write();
        vec_radius_slices[i]->Write();
    }

    for(int i=0;i<20;i++){vec_z_slices_in_radius[i]->Write();}

    histo_eta->Write();
    histo_eta_larger_range->Write();
    histo_delta_eta->Write();
    histo_eta_z_cut ->Write();
    histo_delta_eta_z_cut ->Write();
    histo_sum_eta_z_cut->Write();
    radius_ortho_vs_z_eta->Write();
    radius_ortho_vs_z_eta2->Write();
    radius_ortho_vs_z_eta3->Write();
    radius_ortho_vs_z_eta4->Write();
    radius_ortho_vs_z_eta5->Write();

    //-------------------------------------------------------------

    output_histos ->cd();
    output_histos ->mkdir("S_radius");
    output_histos ->cd("S_radius");
    S_radius_from_origin->Write();

    for(int i=0;i<6;i++)
    {
        vec_S_radius_from_origin[i]->Write();
    }
    //------------------------------------------------------------------

    output_histos ->cd();
    output_histos ->mkdir("type1");
    output_histos ->cd("type1");
    histo_invmass_Lambda_type1 ->Write();
    histo_invmass_Lambda_type1_by_track->Write();
    histo_invmass_Lambda_type11->Write();
    histo_invmass_Lambda_type11_by_track->Write();
    histo_invmass_Lambda_type1_dcaprim ->Write();
    histo_invmass_Lambda_type1_dcaprim_1 ->Write();
    histo_invmass_Lambda_overlap_cut_type1 ->Write();
    histo_invmass_Lambda_overlap_cut_dcaprim_type1 ->Write();
    m_squared_proton_type1->Write();
    m_squared_proton_type11->Write();
    histo_diff_inv_mass_type1->Write();

    histo_invariantmass_K0_type1->Write();
    histo_invariantmass_K0_type11->Write();
    histo_invariantmass_K0_type1_dcaprim ->Write();
    histo_invariantmass_K0_type1_overlap ->Write();
    histo_invariantmass_K0_type1_cut_on_L ->Write();

    //---------------------------------------------------

    //type2
    output_histos ->cd();
    output_histos ->mkdir("type2");
    output_histos ->cd("type2");
    histo_m_squared_kaon_no_other_cuts_type2->Write();
    histo_m_squared_kaon_with_all_other_cuts_type2->Write();
    histo_m_squared_kaon_with_some_other_cuts_type2->Write();
    mass_squared_all_kaons_type2->Write();
    mass_squared_all_kaons_overlap_cuts_type2->Write();
    histo_invariant_mass_kaon_no_other_cuts->Write();
    histo_invariant_mass_kaon_with_other_cuts->Write();
    m_squared_anti_proton_type2->Write();

    output_histos ->cd();
    output_histos ->mkdir("type2_msquared");
    output_histos ->cd("type2_msquared");

    for(int i=0;i<17*2;i++)
    {
        vec_m_squared_type_2_p[i]->Write();
        vec_m_squared_type_2_K_plus1[i]->Write();
        vec_m_squared_type_2_K_plus2[i]->Write();

    }

    for(int i=0;i<17*2;i++)
    {
        vec_m_squared_type_21_p[i]->Write();
        vec_m_squared_type_21_K_plus1[i]->Write();
        vec_m_squared_type_21_K_plus2[i]->Write();
    }

    for(int i=17*2;i<17*3;i++)
    {
        vec_m_squared_type_2_p[i]->Write();
    }

    for(int i=17*2;i<17*3;i++)
    {
        vec_m_squared_type_21_p[i]->Write();
    }

    //--------------------------------------------------------------
    //type3
    output_histos ->cd();
    output_histos ->mkdir("type3");
    output_histos ->cd("type3");
    mass_squared_kaons_type3->Write();
    mass_squared_kaons_type3_with_cut_on_antip->Write();
    mass_squared_kaons_type3_with_cut_on_invmass_K0->Write();

    mass_squared_antip_type3->Write();
    mass_squared_antip_type3_with_cut_on_K->Write();
    mass_squared_antip_type3_with_cut_on_invmass_K0->Write();

    //histo_invariantmass_K0_type3_with_cuts->Write();
    //histo_invariantmass_K0_type3->Write();
    histo_invariantmass_K0_type3_with_cut_on_antip->Write();
    histo_invariantmass_K0_type3_with_cuts_on_antip_and_K->Write();
    histo_invariantmass_K0_type3_with_cuts_on_K->Write();

    histo_angle_ch3->Write();
    histo_invariantmass_K0_type3->Write();

    vec_S_mass_ch3[0]->Write();
    vec_S_mass_ch3[1]->Write();
    vec_S_mass_ch3[2]->Write();
    vec_S_mass_ch3[3]->Write();
    vec_S_mass_ch31[0]->Write();
    vec_S_mass_ch31[1]->Write();
    vec_S_mass_ch31[2]->Write();
    vec_S_mass_ch31[3]->Write();

    hist_radiusK0->Write();

    output_histos ->cd();
    output_histos ->mkdir("invmass_K0_type3");
    output_histos ->cd("invmass_K0_type3");
    for(int i=0;i<6*17;i++)
    {
        vec_invmass_K0_type3[i]->Write();
    }

    for(int i=0;i<6*17;i++)
    {
        vec_invmass_K0_type31[i]->Write();
    }

    for(int i=0;i<17*2;i++)
    {
        for(int j=0;j<6;j++)
        {
            vec_vec_invmass_K0_type3[i][j]->Write();
        }
    }

    for(int i=0;i<17*2;i++)
    {
        for(int j=0;j<6;j++)
        {
            vec_vec_invmass_K0_type31[i][j]->Write();
        }
    }


    output_histos ->cd();
    output_histos ->mkdir("m_squared_type3");
    output_histos ->cd("m_squared_type3");

    for(int i=0;i<5;i++)
    {
        vec_m_sq_p_type_3_for_dcas[i]->Write();
        vec_m_sq_p_type_31_for_dcas[i]->Write();
    }

    for(int i=0;i<34;i++)
    {
        vec_m_squared_type_3_K[i]->Write();
        vec_m_squared_type_3_p[i]->Write();
        vec_m_squared_type_3_add_pi[i]->Write();
        vec_m_squared_type_3_K01[i]->Write();
        vec_m_squared_type_3_K02[i]->Write();

        vec_m_squared_type_31_K[i]->Write();
        vec_m_squared_type_31_p[i]->Write();
        vec_m_squared_type_31_add_pi[i]->Write();
        vec_m_squared_type_31_K01[i]->Write();
        vec_m_squared_type_31_K02[i]->Write();
    }


    output_histos ->cd();
    output_histos ->mkdir("S_mass_ch3");
    output_histos ->cd("S_mass_ch3");

    for(int i=0;i<vec_S_mass_ch3_normalband.size();i++)
    {
        vec_S_mass_ch3_normalband[i]->Write();
        vec_S_mass_ch3_sideband[i]->Write();
    }


    //-----------------------------------------------------------------

    /*
    histo_delta[0] ->Write();
    histo_delta[1] ->Write();
    histo_delta[2] ->Write();
    */
    //dEdx
    output_histos ->cd();
    output_histos ->mkdir("dEdx");
    output_histos ->cd("dEdx");
    
    dEdx_vs_charge_dot_momentum_S_cand->Write();
    vec_dEdx[0]->Write();
    vec_dEdx[1]->Write();
    vec_dEdx[2]->Write();
    vec_dEdx[3]->Write();

    vec_dEdx_S[0]->Write();
    vec_dEdx_S[1]->Write();
    vec_dEdx_S[2]->Write();

    vec_mass_squared_dEdx_selected[0]->Write();
    vec_mass_squared_dEdx_selected[1]->Write();
    vec_mass_squared_dEdx_selected[2]->Write();

    histo_pt_S->Write();
    tprof->Write();

    //--------------------------------------------------------
    

    //ch5------------------------------------------------------
    output_histos ->mkdir("m_squared_ch5_and_51");
    output_histos ->cd("m_squared_ch5_and_51");

    

    for(int i=0;i<8;i++)
    {
        histos_m_squared_type5[i]->Write();
        histos_m_squared_type51[i]->Write();
    }
    for(int i=0;i<17;i++)
    {
        m_squared_type5_K_r[i]->Write();
        m_squared_type5_K_r[i+17]->Write();
        m_squared_type5_pi_minus_r[i]->Write();
        m_squared_type5_pi_minus_r[i+17]->Write();
        m_squared_type5_pi_plus_r[i]->Write();
        m_squared_type5_pi_plus_r[i+17]->Write();
        m_squared_type5_p_r[i]->Write();
        m_squared_type5_p_r[i+17]->Write();

        m_squared_type51_K_r[i]->Write();
        m_squared_type51_K_r[i+17]->Write();
        m_squared_type51_pi_minus_r[i]->Write();
        m_squared_type51_pi_minus_r[i+17]->Write();
        m_squared_type51_pi_plus_r[i]->Write();
        m_squared_type51_pi_plus_r[i+17]->Write();
        m_squared_type51_p_r[i]->Write();
        m_squared_type51_p_r[i+17]->Write();


    }
    output_histos ->cd();
    output_histos ->mkdir("invmass_L_and_AntiL_ch5");
    output_histos ->cd("invmass_L_and_AntiL_ch5");

    for(int i=0;i<12+5*17;i++)
    {
        vec_invmass_Lambda_type5[i]->Write();
    }
    for(int i=0;i<12+5*17;i++)
    {
        vec_invmass_Lambda_type51[i]->Write();
    }

    for(int i=0;i<2*17;i++)
    {
        for(int j=0;j<6;j++)
        {
            vec_vec_invmass_L_type5[i][j]->Write();
        }
    }

    for(int i=0;i<17*2;i++)
    {
        for(int j=0;j<6;j++)
        {
            vec_vec_invmass_L_type51[i][j]->Write();
        }
    }

    output_histos ->cd();
    output_histos ->mkdir("S_mass_ch5");
    output_histos ->cd("S_mass_ch5");

    for(int i=0;i<vec_S_mass_ch5_normalband.size();i++)
    {
        vec_S_mass_ch5_normalband[i]->Write();
        vec_S_mass_ch5_sideband[i]->Write();
    }



    output_histos ->cd();
    output_histos ->mkdir("other_ch5");
    output_histos ->cd("other_ch5");


    histo_invmass_Lambda_type5->Write();
    histo_invmass_Lambda_by_V0_type5 ->Write();
    histo_invmass_Lambda_type51->Write();
    histo_invmass_Lambda_dcaprim_line_type5->Write();
    histo_invmass_Lambda_dcaprim_type51->Write();
    histo_invmass_Lambda_dcaprim1_type5->Write();
    histo_invmass_Lambda_antip_cut_type5 ->Write();
    histo_invmass_Lambda_antip_cut_type51 ->Write();
    histo_invmass_Lambda_antip_cut_byV0_type51->Write();
    histo_invmass_Lambda_radius_cut_type5 ->Write();
    histo_invmass_Lambda_radius_cut_type51 ->Write();
    histo_invmass_Lambda_anti_cut_type5->Write();
    histo_invmass_Lambda_anti_cut_type51->Write();
    histo_invmass_Lambda_final_cut_type5 ->Write();
    histo_invmass_Lambda_final_cut_type51->Write();

    histo_S_mass_type5->Write();
    histo_S_mass_type5_cut_on_antip->Write();
    histo_S_mass_type5_anticuts->Write();
    histo_S_mass_type5_dcaprim->Write();

    histo_S_mass_type5_anticuts_angle_90->Write();
    histo_S_mass_type5_anticuts_angle_120->Write();
    histo_S_mass_type5_anticuts_angle_150->Write();
    dEdx_type5_anticuts->Write();
    histo_S_mass_type5_anticuts_dcaprim_angle_120_and_radius_10->Write();
    histo_S_mass_type5_angle_lambda->Write();
    histo_counter_type5->Write();
    histo_S_mass_type5_angle_lambda_and_dcaprim->Write();
    histo_S_mass_type51_angle_lambda_and_dcaprim->Write();
    histo_eta_type_5->Write();
    histo_phi_type_5->Write();
    histo_S_mass_mixed_event_no_dcaprim->Write();
    histo_S_mass_mixed_event_with_dcaprim->Write();

    for(int i=0;i<6;i++)
    {
        vec_invmass_S_type5[i]->Write();
    }
    for(int i=0;i<6;i++)
    {
        vec_invmass_S_type51[i]->Write();
    }

    //--------------------------------------------------------

    

    

    

    
    


    //-----------------------------------------------------

    TCanvas* mycan = new TCanvas();
    mycan->cd();
    histo_its->Draw();

    //}
    /*
    can14->cd();
    histo_path->Draw();

    can15->cd();
    histo_nitscls->Draw();

    can16->cd();
    histo_path_has_no_its->Draw();

    can17->cd();
    histo_path_has_its->Draw();

    can18->cd();
    histo_path_has_atleast_3_its_hits->Draw();
    */

    cout<<"counterdraw: "<<counter_draw<<endl;

}
#endif
