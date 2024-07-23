#include <string.h>
#include "helpers.h"
#include "bam.h"
#include "binary.h"

BinaryType binary_type_from_string(char *binary_type_str)
{
    if (strcmp(binary_type_str, "BNS") == 0)
    {
        return BNS;
    }
    if (strcmp(binary_type_str, "BBH") == 0)
    {
        return BBH;
    }
    if (strcmp(binary_type_str, "BHNS") == 0)
    {
        return BHNS;
    }
    return -1;
}

void set_binary_parameters(BinaryInfo *bi)
{
    Setd("mass1", bi->mass1);
    Setd("mass2", bi->mass2);
    Setd("px1", bi->position_x1);
    Setd("px2", bi->position_x2);

    log_printf("set binary parameters: mass1=%lf, mass2=%lf, px1=%lf, px2=%lf",
               bi->mass1, bi->mass2, bi->position_x1, bi->position_x2);
}