
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
  // calculating Cross point 
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of 
  // base1 + l * dir2 =  base1 + m * dir2 
  //
  // set of lin. equations:
  //  
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0) 
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //     
  // to solve multiply the equation by the transposed Matrix of M from the left: M 
  //     
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                    
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .
  // M * M = |                                                                                                                |
  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                        
  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  // 
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  // 
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  // 
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |              
  //
  // l = Dl/D = l_Cross
  // 
  // vertex = base1 + l_Cross * dir1
  //

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

  //cout << "Cross point calculated" << endl;
  return TVector3(vertex);

 // return TVector3(-20000.,-20000.,-20000.);
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

