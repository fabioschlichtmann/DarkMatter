
#include "Dark_Matter_Read.h"

void Macro_Loop_event()
{
    printf("Macro_Loop_event started \n");

    Dark_Matter_Read* DM_Read = new Dark_Matter_Read();
    DM_Read ->Init_tree("file_list.txt");

    Long64_t numentries = DM_Read -> getnumberentries();
    printf("num entries: %i \n",numentries);

    TH1D* histo_invariantmass = new TH1D("histo inv mass","histo inv mass",10,1,2);

    for(int i =0 ; i<numentries; i++)
    {
        DM_Read ->Loop_event(i,histo_invariantmass);
    }
    histo_invariantmass->Draw();
}