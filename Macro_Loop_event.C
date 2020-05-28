
#include "Dark_Matter_Read.h"


void Macro_Loop_event()
{
    printf("Macro_Loop_event started \n");
    //TFile* outputfile = new TFile("Histos.root","RECREATE");
    
    Dark_Matter_Read* DM_Read = new Dark_Matter_Read();
    DM_Read ->Init_tree("file_list.txt");

    Long64_t numentries = DM_Read -> getnumberentries();
    printf("num entries: %i \n",numentries);

    
    //--------------------------------------------------------------------------------------------------
    //for(int i =0 ; i<numentries; i++)
    for(int i =0 ; i<1000; i++)
    {
        //printf("eventcounter: %f",event_counter);
        DM_Read ->Loop_event(i);
    }

    DM_Read->Save();

    //----------------------------------------------------------------------------------------------------


}