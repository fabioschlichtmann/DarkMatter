#ifndef __ALI_AS_EVENT_V2_H__
#define __ALI_AS_EVENT_V2_H__

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//----------------------------------------------------------------------------------------
class Ali_AS_TRD_digit : public TObject
{
private:
    // Digit data
    UShort_t hit_ids[2]; // contains the full information of the TRD hit position (sector, stack, layer, row, column), definition see below
    Short_t ADC_time_values[24]; // raw ADC values for all 24 time bins of a single pad
    Short_t arr_pos[3][24];
    //Short_t arr_pos_uncalib[3][24]
    //Short_t ADC_time_values_corrected_tc[24]; // raw ADC values for all 24 time bins of a single pad
    Short_t  dca_to_track; // distance of closest approach of digit TRD hit to TPC track
    Short_t  dca_x;
    Short_t  dca_y;
    Short_t  dca_z;
    Short_t  ImpactAngle; // impact angle of TPC track to this TRD digit

public:
    Ali_AS_TRD_digit() :
	//hit_ids(), ADC_time_values(), ADC_time_values_corrected(), ADC_time_values_corrected_tc(), dca_to_track(-1), dca_x(-1), dca_y(-1), dca_z(-1), ImpactAngle(-1)
        hit_ids(), ADC_time_values(), arr_pos(), dca_to_track(-1), dca_x(-1), dca_y(-1), dca_z(-1), ImpactAngle(-1)
    {
    }
	~Ali_AS_TRD_digit() {}

	// setters
	void sethit_ids(UShort_t x, UShort_t y)                     { hit_ids[0] = x; hit_ids[1] = y; }
	void setADC_time_value(Int_t time_bin, Short_t ADC_value)  { ADC_time_values[time_bin] = ADC_value; }
        void set_pos(Int_t time_bin, Float_t x_pos, Float_t y_pos, Float_t z_pos) //{arr_pos[0][time_bin] = (Short_t)(100.0*x_pos); arr_pos[1][time_bin] = (Short_t)(100.0*y_pos); arr_pos[2][time_bin] = (Short_t)(100.0*z_pos);}
        {
            // Subtract a vector of length 200 cm from digit positions to reduce the range. Usually ~300-360 cm, now 100-160 cm.
            // This is important to not get out of Short_t range: 8 byte: 2^16/2 = 32768. 360*100 > 32768
            TVector3 vec_200cm(x_pos,y_pos,z_pos);
            //printf("pos: {%4.3f, %4.3f, %4.3f} \n",x_pos,y_pos,z_pos);
            if(vec_200cm.Mag() > 0.0)
            {
                vec_200cm *= 200.0/vec_200cm.Mag();
                //printf("vector length: %4.3f \n",vec_200cm.Mag());
                arr_pos[0][time_bin] = (Short_t)(100.0*(x_pos-vec_200cm.X()));
                arr_pos[1][time_bin] = (Short_t)(100.0*(y_pos-vec_200cm.Y()));
                arr_pos[2][time_bin] = (Short_t)(100.0*(z_pos-vec_200cm.Z()));
            }
            else
            {
                arr_pos[0][time_bin] = 0.0;
                arr_pos[1][time_bin] = 0.0;
                arr_pos[2][time_bin] = 0.0;
            }

            //cout << "x_pos: " << x_pos << endl;
            //cout << "y_pos: " << y_pos << endl;
            //cout << "z_pos: " << z_pos << endl;
            //cout << "x_pos-200: " << arr_pos[0][time_bin] << endl;
            //cout << "y_pos-200: " << arr_pos[1][time_bin] << endl;
            //cout << "z_pos-200: " << arr_pos[2][time_bin] << endl;
        }

        //void set_pos_uncalib(Int_t time_bin, Float_t x_pos_uncalib, Float_t y_pos_uncalib, Float_t z_pos_uncalib)  { arr_pos_uncalib[0][time_bin] = (Short_t)(100.0*x_pos_uncalib); arr_pos[1][time_bin] = (Short_t)(100.0*y_pos); arr_pos_uncalib[2][time_bin] = (Short_t)(100.0*z_pos_uncalib); }


        //void setADC_time_value_corrected_tc(Int_t time_bin, Short_t ADC_value)  { ADC_time_values_corrected_tc[time_bin] = ADC_value; }
        void setdca_to_track(Float_t f, Float_t f_x, Float_t f_y, Float_t f_z)
        {
            dca_to_track = (Short_t)(100.0*f);
            dca_x = (Short_t)(100.0*f_x);
            dca_y = (Short_t)(100.0*f_y);
            dca_z = (Short_t)(100.0*f_z);
        }
	void setImpactAngle(Float_t f)                              { ImpactAngle = (Short_t)(100.0*f); }

	// getters
        UInt_t   gethit_ids(Int_t i)                const           { return hit_ids[i]; }
	Short_t  getADC_time_value(Int_t time_bin)  const           { return ADC_time_values[time_bin]; }
        Float_t  get_pos(Int_t time_bin, Int_t index)  const //{ return ((Float_t)arr_pos[index][time_bin])/100.0; }
        {
            TVector3 vec_200cm((Float_t)arr_pos[0][time_bin],(Float_t)arr_pos[1][time_bin],(Float_t)arr_pos[2][time_bin]);

            // now vec_200cm is again in cm
            if(vec_200cm.Mag() > 0.0)
            {
                vec_200cm *= 200.0/vec_200cm.Mag();     //when we get positions XYZ-200cm
                return (Float_t)(arr_pos[index][time_bin]/100.0 + vec_200cm(index));
            }
            else return 0.0;
        }

        //Float_t  get_pos_uncalib(Int_t time_bin, Int_t index)  const { return ((Float_t)arr_pos_uncalib[index][time_bin])/100.0; }
        //Short_t  getADC_time_value_corrected_tc(Int_t time_bin)  const { return ADC_time_values_corrected_tc[time_bin]; }
	Float_t  getdca_to_track()                  const           { return ((Float_t)dca_to_track)/100.0; }
	Float_t  getdca_x()                         const           { return ((Float_t)dca_x)/100.0; }
	Float_t  getdca_y()                         const           { return ((Float_t)dca_y)/100.0; }
        Float_t  getdca_z()                         const           { return ((Float_t)dca_z)/100.0; }
	Float_t  getImpactAngle()                   const           { return ((Float_t)ImpactAngle)/100.0; }

	// x = TRD_col + TRD_sec*144;
        // y = TRD_row + TRD_stack*16 + TRD_lay*16*5;
	Int_t    get_sector()                       const           { return (Int_t)((Float_t)hit_ids[0]/144.0); }
	Int_t    get_column()                       const           { return (hit_ids[0]%144); }
	Int_t    get_layer()                        const           { return (Int_t)((Float_t)hit_ids[1]/(16.0*5.0)); }
	Int_t    get_stack()                        const           { return ((Int_t)((Float_t)hit_ids[1]/16.0))%5; }
	Int_t    get_row()                          const           { return (hit_ids[1]%16); }

        Int_t    get_detector(Int_t layer, Int_t stack, Int_t sector) const { return (layer + stack * 6 + sector * 6 * 5);}

	ClassDef(Ali_AS_TRD_digit,1);  //
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Ali_AS_offline_Tracklet : public TObject
{
private:
    // Tracklet properties
    Short_t detector;
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Float_t  chi2;
    Float_t  refY;
    Float_t  refZ;
    Float_t  locY;
    Float_t  locZ;
    Float_t  refdYdx;
    Float_t  refdZdx;
    Float_t  locdYdx;
    Float_t  locdZdx;

public:
    Ali_AS_offline_Tracklet() :
        detector(0), TV3_offset(), TV3_dir(), chi2(0), refY(0), refZ(0), locY(0), locZ(0), refdYdx(0),
        refdZdx(0), locdYdx(0), locdZdx(0)
    {
    }
        ~Ali_AS_offline_Tracklet()
        {
        }

	// setters
        void set_detector(Short_t s)                     { detector = s;         }
        void set_TV3_offset(TVector3 tv3)                { TV3_offset = tv3;     }
        void set_TV3_dir(TVector3 tv3)                   { TV3_dir = tv3;        }
        void set_chi2(Float_t chi2in)                    { chi2 = chi2in;        }
        void set_refYZ(Float_t refYin, Float_t refZin)   { refY = refYin; refZ = refZin;}
        void set_refdYdZdx(Float_t refdYdxin, Float_t refdZdxin)   { refdYdx = refdYdxin; refdZdx = refdZdxin;}
        void set_locYZ(Float_t locYin, Float_t locZin)   { locY = locYin; locZ = locZin;}
        void set_locdYdZdx(Float_t locdYdxin, Float_t locdZdxin)   { locdYdx = locdYdxin; locdZdx = locdZdxin;}

	// getters
	Short_t get_detector() const                     { return detector;         }
        TVector3 get_TV3_offset() const                  { return TV3_offset;       }
        TVector3 get_TV3_dir() const                     { return TV3_dir;          }
        Float_t get_chi2() const                         { return chi2;             }
        Float_t get_refZ() const                         { return refZ;}
        Float_t get_refY() const                         { return refY;}
        Float_t get_refdZdx() const                      { return refdZdx;}
        Float_t get_refdYdx() const                      { return refdYdx;}
        Float_t get_locZ() const                         { return locZ;}
        Float_t get_locY() const                         { return locY;}
        Float_t get_locdZdx() const                      { return locdZdx;}
        Float_t get_locdYdx() const                      { return locdYdx;}

        ClassDef(Ali_AS_offline_Tracklet,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Ali_AS_Track : public TObject
{
private:
    // Track properties
    Float_t        nsigma_e_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_e_TOF; // nsigma TOF of particle
    Float_t        nsigma_pi_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_pi_TOF; // nsigma TOF of particle
    Float_t        nsigma_K_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_K_TOF; // nsigma TOF of particle
    Float_t        nsigma_p_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_p_TOF; // nsigma TOF of particle
    Float_t        TRD_signal; // TRD PID signal
    Float_t        TRDsumADC; // TRD electron PID probability
    Float_t        dca; // distance of closest approach of particle A
    TLorentzVector TLV_part; // Lorentz vector properties of this particle
    UShort_t       NTPCcls; // Number of TPC clusters
    UShort_t       NTRDcls; // Number of TRD clusters
    UShort_t       NITScls; // Number of TRD clusters
    UShort_t       status; // status of track: bit 0: ITS refit, bit1: TPC refit
    Float_t        TPCchi2; // TPC chi2
    ULong64_t      TRD_ADC_time_layer[6];
    Float_t        impact_angle_on_TRD; // Track impact angle on TRD
    Float_t        TPCdEdx; // Energy loss information of TPC
    Float_t        TOFsignal; // Time-of-flight
    Float_t        Track_length; // length of track
    Float_t        aliHelix_params[9];

    //UShort_t      fNumTRDdigits; // number of TRD digits for this track
    //UShort_t      fNumOfflineTracklets; // number of offline tracklets

    Int_t track_id;

    //TClonesArray* fTRD_digits;          //->
    //TClonesArray* fOfflineTracklets;    //->

    TBits tbits_fit;
    TBits tbits_shared;

public:
    Ali_AS_Track() :
	nsigma_e_TPC(-1),nsigma_e_TOF(-1),nsigma_pi_TPC(-1),nsigma_pi_TOF(-1),nsigma_K_TPC(-1),nsigma_K_TOF(-1),nsigma_p_TPC(-1),nsigma_p_TOF(-1),TRD_signal(-1),
        TRDsumADC(-1),dca(-1),TLV_part(),NTPCcls(-1),NTRDcls(-1),NITScls(-1),status(-1),TPCchi2(-1),TRD_ADC_time_layer(),
        impact_angle_on_TRD(-1),TPCdEdx(-1),TOFsignal(-1),Track_length(-1),aliHelix_params(),track_id(0)
    {
        //fTRD_digits       = new TClonesArray( "Ali_AS_TRD_digit", 10 );
        //fOfflineTracklets = new TClonesArray( "Ali_AS_offline_Tracklet", 10 );
    }
	~Ali_AS_Track()
	{
            //delete fTRD_digits;
            //delete fOfflineTracklets;
            //fTRD_digits       = NULL;
            //fOfflineTracklets = NULL;
	}

	// setters
	void setnsigma_e_TPC(Float_t f)                     { nsigma_e_TPC = f;         }
	void setnsigma_e_TOF(Float_t f)                     { nsigma_e_TOF = f;         }
	void setnsigma_pi_TPC(Float_t f)                     { nsigma_pi_TPC = f;         }
	void setnsigma_pi_TOF(Float_t f)                     { nsigma_pi_TOF = f;         }
	void setnsigma_K_TPC(Float_t f)                     { nsigma_K_TPC = f;         }
	void setnsigma_K_TOF(Float_t f)                     { nsigma_K_TOF = f;         }
	void setnsigma_p_TPC(Float_t f)                     { nsigma_p_TPC = f;         }
	void setnsigma_p_TOF(Float_t f)                     { nsigma_p_TOF = f;         }
	void setTRDSignal(Float_t f)                     { TRD_signal = f;         }
	void setTRDsumADC(Float_t f)                     { TRDsumADC = f;         }
	void setdca(Float_t f)                    { dca = f;        }
	void set_TLV_part(TLorentzVector tlv)     { TLV_part = tlv; }
	void setNTPCcls(UShort_t s)               { NTPCcls = s;}
	void setNTRDcls(UShort_t s)               { NTRDcls = s;}
	void setNITScls(UShort_t s)               { NITScls = s;}
	void setStatus(UShort_t s)                { status = s;}
	void setTPCchi2(Float_t f)                { TPCchi2 = f;}
        void setTRD_layer(Int_t i_layer, ULong64_t l)  { TRD_ADC_time_layer[i_layer] = l;}
        void setimpact_angle_on_TRD(Float_t f)           {impact_angle_on_TRD = f;}
	void setTPCdEdx(Float_t f)                       {TPCdEdx = f;}
	void setTOFsignal(Float_t f)                     {TOFsignal = f;}
        void setTrack_length(Float_t f)                  {Track_length = f;}
        void setHelix(Float_t a, Float_t b,Float_t c,Float_t d,Float_t e,Float_t f,Float_t g,Float_t h,Float_t i)
        {
            aliHelix_params[0] = a;
            aliHelix_params[1] = b;
            aliHelix_params[2] = c;
            aliHelix_params[3] = d;
            aliHelix_params[4] = e;
            aliHelix_params[5] = f;
            aliHelix_params[6] = g;
            aliHelix_params[7] = h;
            aliHelix_params[8] = i;
        }

        void settrackid(Int_t id)  {track_id = id;}

        void settbitsfit(TBits bits) {tbits_fit=bits;}
        void settbitsshared(TBits bits) {tbits_shared=bits;}

	// getters
	Float_t getnsigma_e_TPC() const                     { return nsigma_e_TPC;         }
	Float_t getnsigma_e_TOF() const                     { return nsigma_e_TOF;         }
	Float_t getnsigma_pi_TPC() const                     { return nsigma_pi_TPC;         }
	Float_t getnsigma_pi_TOF() const                     { return nsigma_pi_TOF;         }
	Float_t getnsigma_K_TPC() const                     { return nsigma_K_TPC;         }
	Float_t getnsigma_K_TOF() const                     { return nsigma_K_TOF;         }
	Float_t getnsigma_p_TPC() const                     { return nsigma_p_TPC;         }
	Float_t getnsigma_p_TOF() const                     { return nsigma_p_TOF;         }
	Float_t getTRDSignal() const                     { return TRD_signal;         }
	Float_t getTRDsumADC() const                     { return TRDsumADC;         }
	Float_t getdca() const                    { return dca;        }
	TLorentzVector get_TLV_part() const       { return TLV_part;   }
	UShort_t getNTPCcls() const               { return NTPCcls;    }
	UShort_t getNTRDcls() const               { return NTRDcls;    }
	UShort_t getNITScls() const               { return NITScls;    }
	UShort_t getStatus() const               { return status;    }
	Float_t  getTPCchi2() const              { return TPCchi2; }
        ULong64_t getTRD_layer(Int_t i_layer) const   { return TRD_ADC_time_layer[i_layer]; }
        Float_t   getimpact_angle_on_TRD() const    { return impact_angle_on_TRD; }
	Float_t   getTPCdEdx() const                { return TPCdEdx; }
	Float_t   getTOFsignal() const              { return TOFsignal; }
        Float_t   getTrack_length() const           { return Track_length; }
        Float_t   getHelix_param(Int_t i_param) const              {return aliHelix_params[i_param]; }

        Int_t gettrackid()  {return track_id;}
        void Evaluate(Double_t t, // helix evaluation, taken from AliHelix
                      Double_t r[3]);  //radius vector

        TBits getbitsfit() {return tbits_fit;}
        TBits getbitsshared() {return tbits_shared;}


	Float_t   getTRD_ADC(Int_t i_layer, Int_t i_time_bin) const
	{
	    if(i_layer < 0 || i_layer > 5 || i_time_bin < 0 || i_time_bin > 7) return -1; // out of range
	    ULong64_t TRD_value = 0;
	    for(Int_t i_bit = 0; i_bit < 8; i_bit++) // One single time slice 8 bit = 256
	    {
		Int_t bitcheck = i_bit + 8*i_time_bin; // range: 0..63 = 64 bit = 8 byte = Long64_t
		Int_t bit_status = (TRD_ADC_time_layer[i_layer] >> bitcheck) & 1; // check bit bitcheck
		if(bit_status) TRD_value |= (ULong64_t)1 << i_bit; // setting bit i_bit to 1
	    }
	    Float_t TRD_value_decode = (Float_t)TRD_value * 100.0; // * TRD_ADC_bin_width
	    return TRD_value_decode;
	}

	Int_t     HasITShit_on_layer(Int_t ilayer) { return ((NITScls >> ilayer) & 1);}  // ITShit -> LOL


        //----------------------------
	



	ClassDef(Ali_AS_Track,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------


//________________________________________________________________________
void Ali_AS_Track::Evaluate(Double_t t,Double_t r[3])  //radius vector
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  float phase=aliHelix_params[4]*t+aliHelix_params[2];
  Double_t sn=sinf(phase), cs=cosf(phase);
  //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = aliHelix_params[5] + sn/aliHelix_params[4];
  r[1] = aliHelix_params[0] - cs/aliHelix_params[4];
  r[2] = aliHelix_params[1] + aliHelix_params[3]*t;
}
//________________________________________________________________________




//----------------------------------------------------------------------------------------
class Ali_AS_Tracklet : public TObject
{
private:
    // Tracklet properties
    Short_t detector;
    TVector3 TV3_offset;
    TVector3 TV3_dir;

public:
    Ali_AS_Tracklet() :
	detector(0), TV3_offset(), TV3_dir()
    {
    }
        ~Ali_AS_Tracklet()
        {
        }

	// setters
        void set_detector(Short_t s)                     { detector = s;         }
        void set_TV3_offset(TVector3 tv3)                { TV3_offset = tv3;     }
        void set_TV3_dir(TVector3 tv3)                   { TV3_dir = tv3;        }

	// getters
	Short_t get_detector() const                     { return detector;         }
        TVector3 get_TV3_offset() const                  { return TV3_offset;       }
        TVector3 get_TV3_dir() const                     { return TV3_dir;          }

        ClassDef(Ali_AS_Tracklet,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------
class Ali_AS_V0 : public TObject
{
private:
    Float_t x; // Event vertex x
    Float_t y; // Event vertex y
    Float_t z; // Event vertex z
    Float_t Npx;   //momentum of negative charged particle
    Float_t Npy;
    Float_t Npz;

    Float_t Ppx;    //momentum of positive charged particle
    Float_t Ppy;
    Float_t Ppz;

    Float_t dca_V0;


    Float_t pos[3];
    Float_t momN[3];
    Float_t momP[3];
    //Float_t* mom1 = new Float_t[3];
    
    

    Int_t   N_tracks; // total number of tracks
   

    Int_t      fNumTracks; // number of tracks in event

    TClonesArray* fTracks;      //->

public:
    Ali_AS_V0() :
        x(-1),y(-1),z(-1),Npx(-1),Npy(-1),Npz(-1),Ppx(-1),Ppy(-1),Ppz(-1),dca_V0(0),pos(),momN(),momP(),N_tracks(0),fNumTracks(0)

        //brauchen wir auch nicht?
        /*         
        ,N_TRD_tracklets(0),
	cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),TriggerWord(),fNumTracks(0),fNumTracklets(0),
        ADC_sum_det()  */
    {
        fTracks         = new TClonesArray( "Ali_AS_Track", 2 );
       // fTracklets      = new TClonesArray( "Ali_AS_Tracklet", 10 );
    }
	~Ali_AS_V0()
	{
	    delete fTracks;
            fTracks = NULL;
            //delete fTracklets;
	    //fTracklets = NULL;
	}

        //set and get functions ------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------
	void       setxyz(Float_t r, Float_t s, Float_t t)                    { x = r; y = s; z = t;}

        Float_t*    getxyz()
        {
            pos[0] = x;
            pos[1] = y;
            pos[2] = z;

            return pos;
        }

        

        Float_t getx() const
        {
            return x;

        }

        void setNpxpypz(Float_t r, Float_t s, Float_t t) {Npx=r;Npy=s;Npz=t;}
        void setPpxpypz(Float_t r, Float_t s, Float_t t) {Ppx=r;Ppy=s;Ppz=t;}


        Float_t* getNpxpypz()
        {
            momN[0] = Npx;
            momN[1] = Npy;
            momN[2] = Npz;

            return momN;

        }
        

        Float_t* getPpxpypz()
        {
            momP[0] = Ppx;
            momP[1] = Ppy;
            momP[2] = Ppz;

            return momP;

        }

        void setdcaV0(Float_t dca) {dca_V0 = dca;}

        Float_t getdcaV0() {return dca_V0;}

        //-----------------------------------------------------------------------------
        //--------------------------------------------------------------------------------

	//void       setid(Int_t  r)                    { id = r;                        }
       // Int_t      getid() const                      { return id;                     }

	void       setN_tracks(Int_t r)                 { N_tracks = r;                    }
	Int_t      getN_tracks() const                    { return N_tracks;                 }

       


        //----------------------------

        Ali_AS_Track* createTrack()
	{
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 1 );
	    if (fNumTracks >= 100000)
	    {
		Fatal( "Ali_AS_V0::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

	    new((*fTracks)[fNumTracks++]) Ali_AS_Track;
	    return (Ali_AS_Track*)((*fTracks)[fNumTracks - 1]);
        }

        

	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
        }



	UShort_t getNumTracks() const
	{
	    return fNumTracks;
        }


	Ali_AS_Track* getTrack(Int_t i) const
	{
	    return i < fNumTracks ? (Ali_AS_Track*)((*fTracks)[i]) : NULL;
        }

        //----------------------------


        //----------------------------
       
        //----------------------------

ClassDef(Ali_AS_V0,1);  // A simple event compiled of tracks
};

class Ali_AS_NUCLEV : public TObject
{
private:
    TVector3 secVertex;
    double dcaV0_ESD;

    vector<double> paths;
    vector<double> dcas_to_V0;

    Int_t      fNumTracks;

    TClonesArray* fTracks;   //->

public:
    Ali_AS_NUCLEV() :
	secVertex(),fNumTracks(0),dcaV0_ESD(0)
    {
        fTracks         = new TClonesArray( "Ali_AS_Track", 10 );
    }
        ~Ali_AS_NUCLEV()
        {
	    delete fTracks;
            fTracks = NULL;
        }


        // setters
        void set_secVertex (TVector3 tv3)                { secVertex   = tv3;     }
        void addpath(double path) {paths.push_back(path); }
        void adddca(double dca) {dcas_to_V0.push_back(dca); }
        void setdcaV0_ESD(double dca) {dcaV0_ESD=dca;};

	// getters
        TVector3 get_secVertex()   const                 { return secVertex;      }
        double getpath(int i)
        {
            return paths[i];
        }
        double getdcaV0_ESD(){return dcaV0_ESD;}

        //----------------------------
	Ali_AS_Track* createTrack()
        {
            //cout<<"as_sec_vertex created Track"<<endl;
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 10 );
	    if (fNumTracks >= 100000)
	    {
		Fatal( "Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

	    new((*fTracks)[fNumTracks++]) Ali_AS_Track;
	    return (Ali_AS_Track*)((*fTracks)[fNumTracks - 1]);
        }

	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
        }

	UShort_t getNumTracks() const
	{
	    return fNumTracks;
        }

	Ali_AS_Track* getTrack(Int_t i) const
	{
	    return i < fNumTracks ? (Ali_AS_Track*)((*fTracks)[i]) : NULL;
        }
        //----------------------------

ClassDef(Ali_AS_NUCLEV,1);  // A simple event compiled of tracks
};


class Ali_AS_DM_particle : public TObject
{
private:
    TVector3 primVertex;
    TVector3 S1Vertex;
    TVector3 S2Vertex;
    TVector3 S3Vertex;
    TVector3 DirSV1;
    TVector3 DirSV2;
    TVector3 DirSV3;
    Int_t    N_V0s;

    int type_of_reaction_channel;

    UShort_t      fNumTracks; // number of tracks in event
    UShort_t      fNumV0s; //

    TClonesArray* fTracks;   //->

public:
    Ali_AS_DM_particle() :
	primVertex(),S1Vertex(),S2Vertex(),S3Vertex(),DirSV1(),DirSV2(),DirSV3(),N_V0s(0),fNumTracks(0),fNumV0s(0)
    {
        fTracks         = new TClonesArray( "Ali_AS_Track", 10 );
    }
	~Ali_AS_DM_particle()
	{
	    delete fTracks;
            fTracks = NULL;
        }


        // setters
        void set_primVertex(TVector3 tv3)                { primVertex = tv3;     }
        void set_S1Vertex(TVector3 tv3)                  { S1Vertex   = tv3;     }
        void set_S2Vertex(TVector3 tv3)                  { S2Vertex   = tv3;     }
        void set_S3Vertex(TVector3 tv3)                  { S3Vertex   = tv3;     }
        void set_DirSV1(TVector3 tv3)                    { DirSV1     = tv3;     }
        void set_DirSV2(TVector3 tv3)                    { DirSV2     = tv3;     }
        void set_DirSV3(TVector3 tv3)                    { DirSV3     = tv3;     }
        void setN_V0s(Int_t r)                           { N_V0s = r;            }

        //type 1: antiS + n -> antiLambda + K0 + pi- + pi+

        void settype(Int_t r)                           { type_of_reaction_channel = r;            }

	// getters
        TVector3 get_primVertex() const                  { return primVertex;    }
        TVector3 get_S1Vertex()   const                  { return S1Vertex;      }
        TVector3 get_S2Vertex()   const                  { return S2Vertex;      }
        TVector3 get_S3Vertex()   const                  { return S3Vertex;      }
        TVector3 get_DirSV1()     const                  { return DirSV1;        }
        TVector3 get_DirSV2()     const                  { return DirSV2;        }
        TVector3 get_DirSV3()     const                  { return DirSV3;        }
	Int_t    getN_V0s()       const                  { return N_V0s;         }
        int gettype(){return type_of_reaction_channel;}

        //----------------------------
	Ali_AS_Track* createTrack()
        {
            std::cout<<"as_dm created Track"<<endl;
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 10 );
	    if (fNumTracks >= 100000)
	    {
		Fatal( "Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

	    new((*fTracks)[fNumTracks++]) Ali_AS_Track;
	    return (Ali_AS_Track*)((*fTracks)[fNumTracks - 1]);
        }

	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
        }

	UShort_t getNumTracks() const
	{
	    return fNumTracks;
        }

	Ali_AS_Track* getTrack(Int_t i) const
	{
	    return i < fNumTracks ? (Ali_AS_Track*)((*fTracks)[i]) : NULL;
        }

        
        //----------------------------

ClassDef(Ali_AS_DM_particle,1);  // A simple event compiled of tracks
};


//end of new class definition
//_______________________________________________________________________
//_______________________________________________________________________
//_______________________________________________________________________
//_______________________________________________________________________



//----------------------------------------------------------------------------------------
class Ali_AS_Event : public TObject
{
private:
    Float_t x; // Event vertex x
    Float_t y; // Event vertex y
    Float_t z; // Event vertex z
    Int_t   id; // Run id
    Int_t   N_tracks; // total number of tracks
    Int_t   N_V0s; // total number of tracks
    Int_t   N_TRD_tracklets; // total number of TRD tracklets
    Float_t   cent_class_ZNA; // ZDC neutral A
    Float_t   cent_class_ZNC; // ZDC neutral C
    Float_t   cent_class_V0A; // V0 A
    Float_t   cent_class_V0C; // V0 C
    Float_t   cent_class_V0M; // V0 average
    Float_t   cent_class_CL0; // clusters in layer 0
    Float_t   cent_class_CL1; // clusters in layer 1
    Float_t   cent_class_SPD; // SPD
    Float_t   cent_class_V0MEq; //
    Float_t   cent_class_V0AEq; //
    Float_t   cent_class_V0CEq; //
    Float_t   ADC_sum_det[540];


    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    TString TriggerWord; // Trigger word

    UShort_t      fNumTracks; // number of tracks in event
    Int_t      fNumNUCLEVs; // number of nuclear events in event
    Int_t      fNumDMparticles; // number of DMparticles in event
    Int_t         fNumV0s; // number of tracks in event
    UShort_t      fNumTracklets; // number of tracks in event

    TClonesArray* fTracks;      //->
    TClonesArray* fNUCLEVs;      //->
    TClonesArray* fDMparticles;      //->
    TClonesArray* fV0s;      //->
    TClonesArray* fTracklets;      //->

public:
    Ali_AS_Event() :
	x(-1),y(-1),z(-1),id(-1),N_tracks(0),N_V0s(0),N_TRD_tracklets(0),
	cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),TriggerWord(),fNumTracks(0),fNumNUCLEVs(0),fNumDMparticles(0),fNumV0s(0),fNumTracklets(0),
        ADC_sum_det()
    {
        fTracks         = new TClonesArray( "Ali_AS_Track", 10 );
        fNUCLEVs         = new TClonesArray( "Ali_AS_NUCLEV", 10 );
        fDMparticles         = new TClonesArray( "Ali_AS_DM_particle", 10 );
        fV0s            = new TClonesArray( "Ali_AS_V0", 10 );
        fTracklets      = new TClonesArray( "Ali_AS_Tracklet", 10 );
    }
	~Ali_AS_Event()
	{
	    delete fTracks;
            fTracks = NULL;
            delete fNUCLEVs;
            fNUCLEVs = NULL;
            delete fDMparticles;
            fDMparticles= NULL;
            delete fV0s;
            fV0s = NULL;
            delete fTracklets;
	    fTracklets = NULL;
	}

	void       setx(Float_t r)                    { x = r;                         }
	Float_t    getx() const                       { return x;                      }

	void       sety(Float_t r)                    { y = r;                         }
	Float_t    gety() const                       { return y;                      }

	void       setz(Float_t r)                    { z = r;                         }
	Float_t    getz() const                       { return z;                      }

	void       setid(Int_t  r)                    { id = r;                        }
	Int_t      getid() const                      { return id;                     }

	void       setN_tracks(Int_t r)                 { N_tracks = r;                    }
        Int_t      getN_tracks() const                    { return N_tracks;                 }

        void       setN_V0s(Int_t r)                 { N_V0s = r;                    }
	Int_t      getN_V0s() const                    { return N_V0s;                 }

	void       setN_TRD_tracklets(Int_t r)                 { N_TRD_tracklets = r;                    }
	Int_t      getN_TRD_tracklets() const                    { return N_TRD_tracklets;                 }

	void       setcent_class_ZNA(Float_t r)             { cent_class_ZNA = r;                }
	Float_t      getcent_class_ZNA() const              { return cent_class_ZNA;             }

	void       setcent_class_ZNC(Float_t r)             { cent_class_ZNC = r;                }
	Float_t      getcent_class_ZNC() const              { return cent_class_ZNC;             }

	void       setcent_class_V0A(Float_t r)             { cent_class_V0A = r;                }
	Float_t      getcent_class_V0A() const              { return cent_class_V0A;             }

	void       setcent_class_V0C(Float_t r)             { cent_class_V0C = r;                }
	Float_t      getcent_class_V0C() const              { return cent_class_V0C;             }

	void       setcent_class_V0M(Float_t r)             { cent_class_V0M = r;                }
	Float_t      getcent_class_V0M() const              { return cent_class_V0M;             }

	void       setcent_class_CL0(Float_t r)             { cent_class_CL0 = r;                }
	Float_t      getcent_class_CL0() const              { return cent_class_CL0;             }

	void       setcent_class_CL1(Float_t r)             { cent_class_CL1 = r;                }
	Float_t      getcent_class_CL1() const              { return cent_class_CL1;             }

	void       setcent_class_SPD(Float_t r)             { cent_class_SPD = r;                }
	Float_t      getcent_class_SPD() const              { return cent_class_SPD;             }

	void       setcent_class_V0MEq(Float_t r)             { cent_class_V0MEq = r;                }
	Float_t      getcent_class_V0MEq() const              { return cent_class_V0MEq;             }

	void       setcent_class_V0AEq(Float_t r)             { cent_class_V0AEq = r;                }
	Float_t      getcent_class_V0AEq() const              { return cent_class_V0AEq;             }

	void       setcent_class_V0CEq(Float_t r)             { cent_class_V0CEq = r;                }
	Float_t      getcent_class_V0CEq() const              { return cent_class_V0CEq;             }

	void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
	Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

	void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
	Float_t    getT0zVertex() const               { return T0zVertex;                  }

	void       setTriggerWord(TString s)          { TriggerWord = s;}
        TString    getTriggerWord() const             { return TriggerWord; }

        void       setADC_sum_det(Int_t i_det, Float_t r) { ADC_sum_det[i_det] = r;}
        ULong64_t  getADC_sum_det(Int_t i_det) const   { return ADC_sum_det[i_det]; }


        //----------------------------
	Ali_AS_Track* createTrack()
	{
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 10 );
	    if (fNumTracks >= 100000)
	    {
		Fatal( "Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

            new((*fTracks)[fNumTracks++]) Ali_AS_Track;
            //printf("fNumTracks: %d \n",fNumTracks);
	    return (Ali_AS_Track*)((*fTracks)[fNumTracks - 1]);
        }

        Ali_AS_V0* createV0()
        {
            //cout<<"size: "<< fV0s->GetSize()<<endl;
	    if (fNumV0s == fV0s->GetSize())
            {
                fV0s->Expand( fNumV0s + 10 );
            }
	    if (fNumV0s >= 500000)
	    {
		Fatal( "Ali_AS_Event::createV0()", "ERROR: Too many V0s (>500000)!" );
		exit( 2 );
	    }

	    new((*fV0s)[fNumV0s++]) Ali_AS_V0;
	    return (Ali_AS_V0*)((*fV0s)[fNumV0s - 1]);
        }

        Ali_AS_NUCLEV* createNUCLEV()
        {
            if (fNumNUCLEVs == fNUCLEVs->GetSize())
		fNUCLEVs->Expand( fNumNUCLEVs + 10 );
	    if (fNumNUCLEVs >= 100000)
	    {
		Fatal( "Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

            new((*fNUCLEVs)[fNumNUCLEVs++]) Ali_AS_NUCLEV;
            //printf("fNumNUCLEVs: %d \n",fNumNUCLEVs);
	    return (Ali_AS_NUCLEV*)((*fNUCLEVs)[fNumNUCLEVs - 1]);

        }

        Ali_AS_DM_particle* createDMparticle()
        {
            if (fNumDMparticles == fDMparticles->GetSize())
		fDMparticles->Expand( fNumDMparticles + 10 );
	    if (fNumDMparticles >= 100000)
	    {
		Fatal( "Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>100000)!" );
		exit( 2 );
	    }

            new((*fDMparticles)[fNumDMparticles++]) Ali_AS_DM_particle;
            //printf("fNumDMparticles: %d \n",fNumDMparticles);
	    return (Ali_AS_DM_particle*)((*fDMparticles)[fNumDMparticles - 1]);

        }

       // Ali_AS_NUCLEV* createNUCLEV()
       // {

       // }

	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
        }

        void clearV0List()
	{
	    fNumV0s   = 0;
	    fV0s      ->Clear();
        }

        void clearNUCLEVList()
	{
	    fNumNUCLEVs   = 0;
	    fNUCLEVs      ->Clear();
        }

        void clearDMparticleList()
	{
	    fNumDMparticles   = 0;
	    fDMparticles      ->Clear();
	}

	UShort_t getNumTracks() const
	{
	    return fNumTracks;
        }

        Int_t getNumNUCLEVs() const
	{
	    return fNumNUCLEVs;
        }

        UShort_t getNumV0s() const
	{
	    return fNumV0s;
        }

        int getNumDMs()
        {
            return fNumDMparticles;
        }

	Ali_AS_Track* getTrack(Int_t i) const
	{
	    return i < fNumTracks ? (Ali_AS_Track*)((*fTracks)[i]) : NULL;
        }
        //----------------------------
        Ali_AS_V0* getV0(Int_t i) const
	{
	    return i < fNumV0s ? (Ali_AS_V0*)((*fV0s)[i]) : NULL;
        }

        Ali_AS_NUCLEV* getNUCLEV(Int_t i) const
	{
	    return i < fNumNUCLEVs ? (Ali_AS_NUCLEV*)((*fNUCLEVs)[i]) : NULL;
        }

        Ali_AS_DM_particle* getDMparticle(int i)
        {
            return i < fNumDMparticles ? (Ali_AS_DM_particle*)((*fDMparticles)[i]) : NULL;
        }


        //----------------------------
        Ali_AS_Tracklet* createTracklet() // online tracklet
	{
	    if (fNumTracklets == fTracklets->GetSize())
		fTracklets->Expand( fNumTracklets + 10 );
	    if (fNumTracklets >= 10000)
	    {
		Fatal( "Ali_AS_Event::createTracklet()", "ERROR: Too many tracklets (>100000)!" );
		exit( 2 );
	    }

	    new((*fTracklets)[fNumTracklets++]) Ali_AS_Tracklet;
	    return (Ali_AS_Tracklet*)((*fTracklets)[fNumTracklets - 1]);
	}
	void clearTrackletList()
	{
	    fNumTracklets   = 0;
	    fTracklets      ->Clear();
        }

        
	UShort_t getNumTracklets() const
	{
	    return fNumTracklets;
        }

	Ali_AS_Tracklet* getTracklet(UShort_t i) const
	{
	    return i < fNumTracklets ? (Ali_AS_Tracklet*)((*fTracklets)[i]) : NULL;
        }
        //----------------------------

ClassDef(Ali_AS_Event,1);  // A simple event compiled of tracks
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------





#endif // __ALI_AS_EVENT_V2_H__
