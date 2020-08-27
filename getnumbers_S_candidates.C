void getnumbers_S_candidates()
{
    TFile* file = TFile::Open("/misc/alidata120/alice_u/schlichtmann/out/Merge_S_search_V14_high_stat.root");
    TH1D* histo_counter = (TH1D*)file->Get("histo_counter_S");
    TH1D* histo_type = (TH1D*)file->Get("S_type");

    double number_events = histo_counter->GetBinContent(1);
    double total_number_S_candidates = histo_counter->GetBinContent(2);

    double number_type1 = histo_type ->GetBinContent(1);
    double number_type2 = histo_type ->GetBinContent(2);
    double number_type3 = histo_type ->GetBinContent(3);
    double number_type4 = histo_type ->GetBinContent(4);

    double number_type1_radius_larger20 = histo_counter ->GetBinContent(3);
    double number_type2_radius_larger20 = histo_counter ->GetBinContent(4);
    double number_type3_radius_larger20 = histo_counter ->GetBinContent(5);
    double number_type4_radius_larger20 = histo_counter ->GetBinContent(6);

    cout<<"number events: "<<number_events<<endl;
    cout<<"total number S candidates: "<<total_number_S_candidates<<endl;
    cout<<"type 1: "<<number_type1<<endl;
    cout<<"type 2: "<<number_type2<<endl;
    cout<<"type 3: "<<number_type3<<endl;
    cout<<"type 4: "<<number_type4<<endl;

    cout<<""<<endl;

    cout<<"type 1 radius cut > 20: "<<number_type1_radius_larger20<<endl;
    cout<<"type 2 radius cut > 20: "<<number_type2_radius_larger20<<endl;
    cout<<"type 3 radius cut > 20: "<<number_type3_radius_larger20<<endl;
    cout<<"type 4 radius cut > 20: "<<number_type4_radius_larger20<<endl;


}