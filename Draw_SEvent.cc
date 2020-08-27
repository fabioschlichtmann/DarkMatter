
#include "Draw_SEvent.h"

void Draw_SEvent()
{
    printf("Draw_SEvent \n");
    Ali_Draw_SEvent* CDraw_SEvent = new Ali_Draw_SEvent();
    CDraw_SEvent ->Init_tree("List_data.txt");
    Long64_t event = 15;  //216
    //for(Long64_t event=0;event<7;event++)
   // {
        CDraw_SEvent ->Draw_event(event);
   // }
    /*
    for(int i=0;i<CDraw_SEvent->getfileentries();i++)
    {
        CDraw_SEvent->Cutevent(i);
    }
        */
}