/* AHmod_output.h */
/* mth 12/07 */
/* NL 02/09 */

void AHmod_output(const int ntheta, const double dtheta, 
    const int nphi, const double dphi, const double xc, const double yc, 
    const double zc, double **rr, const char *outdir, 
    const int number, const double time, const double ahf_time);

void AHmod_output_sphercoeff(const int LMAX1,
      const double xc, const double yc, const double zc, 
      double *a0,  double **ac, double **as, const char *outdir, 
      const int number, const double time, const double ahf_time);


void AHmod_output_xyt_vtk(const int ntheta, const double dtheta, 
      const int nphi, const double dphi, 
      const double xc, const double yc, const double zc, double **rr, 
      const char *outdir, const int number, const double time);

void AHmod_Appendxyt(char *rawname,
      const int ntheta, const double dtheta, const int nphi, const double dphi, 
      const double xc, const double yc, const double zc, 
      double **rr, const double time);
      
void AHmod_ConvertRaw2VTK(char *rawname, char *vtkname,
      const int number, const int nphi);
      
