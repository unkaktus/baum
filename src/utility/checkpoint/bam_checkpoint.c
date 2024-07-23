/* bam_checkpoint.c */
/* Bernd Bruegmann 01/2006, Wolfgang Tichy 11/2006 */

#include "bam.h"
#include "checkpoint.h"



void bam_checkpoint()
{
  printf("Adding checkpoint\n");

  /* functions */
  AddFun(PRE_GRID, checkpoint_init, "function to init checkpoint");

  /* parameters */
  AddPar("checkpoint", "no", "whether to checkpoint [no,yes]");
  if (Getv("checkpoint", "no")) return;

  AddPar("checkpoint_di", "1", "how often to checkpoint");
  AddPar("checkpoint_dt", "0", "how often to checkpoint");

  AddPar("checkpoint_dt_hours", "0", "how often to checkpoint");
  AddPar("checkpoint_dt_hours_quit", "0", "when to quit after checkpoint");



  AddPar("checkpoint_DeltaT", "0",
	 "run for that much beyond previous checkpoint");

  AddPar("checkpoint_variables", "auto",
	 "which variables to checkpoint [auto,all]");
  /* additional variables which also need to be saved,
     besides the ones implied by auto */
  AddPar("checkpoint_additional_variables", "",
	 "additional variables we need to save [any varnames]");
  AddPar("checkpoint_previous", "no",
	 "whether to keep previous directory [no]");

  /* special for checkpointing:
     check whether we have a complete set of checkpoint files in 
     the "outdir_previous" directory
  */
  if (checkpoint_checkforfiles("_previous")) {
    /* go into restart mode */
    Sets("checkpoint", "restart");
  }

  AddPar("checkpoint_backup", "no", "whether we want to save the checkpoints somewhere");
  if (Getv("checkpoint_backup", "yes")) {
    AddPar("checkpoint_backup_dir", "", "dir where to save the checkpoints");
  }

  AddPar("checkpoint_resub", "no", "whether we want to resub the simulation when walltime is over");
  if (Getv("checkpoint_resub", "yes")) {
    AddPar("checkpoint_resub_command", "", "command to resub simulation");
  }

}
