
#include "Dark_Matter_Read.h"

void Macro_Loop_event()
{
    printf("Macro_Loop_event started \n");

    Dark_Matter_Read* DM_Read = new Dark_Matter_Read();
    DM_Read ->Init_tree("file_list.txt");
    DM_Read ->Loop_event(0);

}