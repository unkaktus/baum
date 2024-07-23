/* FUKA ID reader */
/* Ivan Markin, 05/09/2023 */

#include "fuka_reader.h"
#include "helpers.h"

void bam_fuka_reader()
{
    // Bail if we are not needed
    if (!Getv("physics", "fuka_reader"))
    {
        log_printf("skipped fuka_reader");
        return;
    }
    printf("Adding fuka_reader\n");

    // Hooks
    AddFun(PRE_GRID, fuka_reader_pre_grid, "Set configuration parameters before making the grid");
    AddFun(INITIALDATA_SET, fuka_reader_set_fields, "Set the field values on the grid");
    AddFun(INITIALDATA_FINISH, fuka_reader_finish, "Mark the interpolated data as ready");

    // Parameters
    AddPar("fuka_reader_binary_type", "", "Type of the binary [BBH, BNS, BHNS]");
    AddPar("fuka_reader_info_filename", "", "Path to the FUKA .info file for the computed ID");
    AddPar("fuka_reader_interpolation_offset", "0.0", "Interpolation offset (in units of r_AH)");
    AddPar("fuka_reader_interpolation_order", "8", "Interpolation order for smooth junk");
    AddPar("fuka_reader_relative_dr_spacing", "0.3", "Relative dr spacing for the interpolation polynomial");
    // To create objects on the grid
    AddPar("mass1", "1.0", "");
    AddPar("mass2", "1.0", "");
}