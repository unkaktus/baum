/* bam_boundary.h */
/* Bernd Bruegmann 01/00 */

void set_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void set_boundary_elliptic(tL *level, tVarList *u);
void set_boundary_elliptic_1d(tL *level, tVarList *u);
extern void (*set_boundary_elliptic_ptr)(tL *level, tVarList *u);
void set_boundary_excision(tL *level, int unew, int ucur);
void copy_onto_excision_boundary(tL *level, tVarList *vlu);
void set_boundary_radiative(tL *level, 
			    int unew, int upre, double c, int ucur,
			    double var0, double v);
void set_boundary_cartoon_radiative(tL *level, 
                                    int unew, int upre, double c, int ucur,
                                    double var0, double v);
void set_boundary_rot_rad(tL *level, 
	                  int ivarnew, int ivarpre, double dt, int ivarcur,
	                  double var0, double v,
		  	  int ibetax_cur, int ibetax_far, int rot, 
		  	  int betaAdv);
void set_boundary_radiative_shells(tL *level, 
                                   int unew, int upre, double c, int ucur,
                                   double var0, double v);
void set_boundary_radiative_shells_centered(tL *level, 
                                   int unew, int upre, double c, int ucur,
                                   double var0, double v);
void set_boundary_radiative_shibata(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void set_boundary_extrapolate(tL *level, int unew);
void extrapolate_to_flagged_boundary(tL *level, tVarList *vlu, int Bflag);

void set_boundary_symmetry(tL *level, tVarList *varlist);
void set_boundary_inversion(tL *level, tVarList *varlist);
void set_boundary_periodic(tL *level, tVarList *varlist);
void set_boundary_plane(tL *level, tVarList *varlist);
void set_boundary_tube(tL *level, tVarList *varlist);
void set_boundary_cartoon_1d(tL *level, tVarList *varlist);
void set_boundary_cartoon_2d(tL *level, tVarList *varlist);

void set_boundary_flags_periodic(tL *level);
void set_boundary_flags_reflection(tL *level);
void set_boundary_flags_symmetry(tL *level);

void setmask(tG *g); 
void setmasklevel(tL *level, int type, 
		  double r0, double x0, double y0, double z0);
 
int set_mask_bhdata(tL *level);
int set_mask_manual(tL *level);
int set_mask_lapseneg(tL *level);

int set_boundary_flags_excision(tL *level);
void find_mask_boundary(tL *level);
void find_mask_normal(tL *level);

void find_robin_normal(tL *level);
void set_boundary_robin(tL *level, tVarList *ul);

void epolmasklevel(tL *level, double **var, int nvars);

void add_rotant_points_to_buffer(tL *level, int vi,
  int *ptr_nbuffer, double **ptr_buffer,
  double *px0, double *py0, double *pz0, int *pnx, int *pny, int *pnz);

/* make visible for external provider */
extern void (*Gauge_add_betarot)(tVarList *vl, double sign);
void boundary_register(void (*f)(tVarList *unew, tVarList *upre, double c, tVarList *ucur));

/* normals.c */
int give_cube_normal (int i, int j, int k, 
		      int imax, int jmax, int kmax, int o,
		      double *shat1, double *shat2, double *shat3);


