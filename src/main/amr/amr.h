/* amr.h */
/* Bernd Bruegmann, 12/99 */


/* box.c */
void level_add_box_withdata(tL *level, tL *level0);
tL *make_level_nested_boxes(tL *level, int pr);
tL *make_level_bboxes(tL *level, int nboxes, tSBox *box, int pr);

/* flag.c */
void set_boundary_flags(tL *level);
void setRPflags(tL *lc, tL *lnew);

/* grid.c */
tG *make_grid_box_local(int l, int m, int n, int o, double x0, double y0,
                        double z0, double dx, double dy, double dz);
tG *make_grid_box(int l, int m, int n, int o, double xmin, double ymin,
                  double zmin, double dx, double dy, double dz, int pr);
tL *make_level_box(int l, int m, int n, int o, double xmin, double ymin,
                   double zmin, double dx, double dy, double dz, int pr);
tL *make_level_bbox(tL *level, double *bbox, int *ibbox, int pr);
tL *make_level_flagregrid(tL *level, int pr);
tL *make_level_box_coarsealigned(int l, int m, int n, int o, 
				 double xmin, double ymin, double zmin, 
				 double dx, double dy, double dz, int pr);
tL *make_level_box_local(int l, int m, int n, int o,
			 double xmin, double ymin, double zmin, 
			 double dx, double dy, double dz);

/* move.c */
void box_join(double *a, double *b);
void find_nested_bboxes(tL *level, int *pnboxes, tSBox *box);

/* regrid.c */
int add_level_while_evolving(tL *level);

/* shells.c */
void restrict_prolong_shells(tG* g, int lfine, int lcoarse, tVarList *uc, tVarList *uf);
void set_boundary_shells_symmetry(tL *level, tVarList *varlist);

/* storage.c */
tL *alloc_level(tG *grid, int nlevel, int nnodes); 















