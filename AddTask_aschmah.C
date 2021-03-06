AliAnalysisTask *AddTask_aschmah(Int_t trigger = 0)
{
    cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << " in add macro " << endl;
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
	Error("AddTask_TRDMatching", "No analysis manager found.");
	return 0;
    }

    //========= Add task to the ANALYSIS manager =====
    Ali_DarkMatter_ESD_analysis *task = new Ali_DarkMatter_ESD_analysis("DigitsTask");
    mgr->AddTask(task); // <-
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":TRD_analysis_hists";
    cout << "outputFileName: " << outputFileName << endl;

    TString outputFileName_B = AliAnalysisManager::GetCommonFileName();
    outputFileName_B += ":EventsAndTracks";
    cout << "outputFileName_B: " << outputFileName_B << endl;

    //----------------------------------
    TString containerName = "TRD_Digits_output";

    TString results = "AnalysisResults.root";

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,results.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,results.Data());

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1); // TList
    mgr->ConnectOutput(task, 2, coutput2); // TTree
    //----------------------------------

    return task;
}
