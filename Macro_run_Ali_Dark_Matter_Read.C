
R__LOAD_LIBRARY(Ali_Dark_Matter_Read_cxx.so);

void Macro_run_Ali_Dark_Matter_Read(TString in_list)
{
    gSystem ->Load("Ali_Dark_Matter_Read_cxx.so");

    Ali_Dark_Matter_Read* DM_Read = new Ali_Dark_Matter_Read(in_list.Data());

    DM_Read ->Init_tree(in_list.Data());

    Long64_t numentries = DM_Read -> getnumberentries();
    printf("num entries: %lld \n",numentries);

    
    //--------------------------------------------------------------------------------------------------
    for(int i_event =0 ; i_event < numentries; i_event++)
    //for(Int_t i_event = 0 ; i_event < 100; i_event++)
    {
        //printf("eventcounter: %f",event_counter);
        DM_Read ->Loop_event(i_event);
    }

    DM_Read->Save();

}
