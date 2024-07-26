/* TOV equations solver
   sbernuz 12/2011  */

#include "bam.h"
#include "BBIDdataReader.h"

#define PRNT 0

/* rest-mass density calculation */
double computeMb(double *r, double *m, double *rho,
                 int n)
{
  int i;
  double sum = 0.;
  double dr, tmp;
  for (i = 1; i < n; i++)
  {
    dr = r[i] - r[i - 1];
    tmp = 1. / sqrt(1. - 2. * m[i] / r[i]);
    sum += r[i] * r[i] * rho[i] * tmp * dr;
  }
  return 4. * PI * sum;
}

/* Lambda calculation */

double computeLambda(double M, double R, double G, double H)
{

  double C, y, tmp1, tmp2, tmp3, Lambda;

  C = M / R;
  y = R * G / H;

  tmp1 = 16. / 15. * (1. - 2. * C) * (1. - 2. * C) * (2. + 2. * C * (y - 1.) - y);
  tmp2 = 2. * C * (6. - 3. * y + 3. * C * (5. * y - 8.)) + 4. * C * C * C * (13. - 11. * y + C * (3. * y - 2.) + 2. * C * C * (1. + y));
  tmp3 = 3. * (1. - 2. * C) * (1. - 2. * C) * (2. - y + 2. * C * (y - 1.)) * log(1. - 2. * C);

  Lambda = tmp1 / (tmp2 + tmp3);

  return Lambda;
}

/* Invariant radius computation */

double computeRinv(double *r, double *m, int n)
{

  int i;
  double sum = 0.;
  double tmp, dr;

  for (i = 1; i < n; i++)
  {
    dr = r[i] - r[i - 1];
    tmp = 1. / sqrt(1. - 2. * m[i] / r[i]);
    sum += dr * tmp;
  }

  return sum;
}

/* Method 1
   x: h coord
   Ref: Lindblom 1992ApJ...398..569L */

#define nvar 4
#define forv for (v = 0; v < nvar; v++)

/* eos wrapper for total energy density */
void eos(double h, double *e, double *p)
{
  double rho, epsl;
  EOS.comp("H", "", "", "pre", "", "", h, p, &rho, &epsl);
  *e = rho * (1. + epsl);
}

/* tov h rhs */
void tov_h_rhs(double dh, double *u, double *k)
{

  double h = u[0];
  double r = u[1];
  double m = u[2];
  double I = u[3];

  double e, p;
  eos(h, &e, &p);

  double r3 = r * r * r;
  double tmp = -(r - 2. * m) / (m + 4 * PI * r3 * p);
  double drdh = r * tmp;
  double dmdh = 4 * PI * e * r3 * tmp;
  double f = sqrt(1. - 2. * m / r);
  double dIdh = (1. - f) / (r * f) * drdh;

  k[0] = 0;
  k[1] = drdh;
  k[2] = dmdh;
  k[3] = dIdh;
}

/* tov h solver */
int tov_h(double hc, int npts,
          double *p_r, double *p_m, double *p_h,
          double *p_rho, double *p_pre, double *p_phi,
          double *p_riso)
{
  int n, v;
  double h, m, r, IR, C, phi, phiR, hR, Mb;
  double e, p, ec, pc, decdh, rho, epsl;

  double **u = (double **)malloc((npts) * sizeof(double *));
  for (n = 0; n < npts; n++)
    u[n] = (double *)malloc((nvar) * sizeof(double));
  double u1[nvar], u2[nvar], u3[nvar], k[nvar];
  double fact = 1. / 6.;

  // spacing
  double dh = hc / npts;

  // central values from EoS
  eos(hc, &ec, &pc);

  // initial data
  h = hc - dh;

  eos(h, &e, &p);
  decdh = (e - ec) / (h - hc);

  r = sqrt(3. * dh / (2. * PI * (ec + 3. * pc))) *
      (1. - 0.25 * (ec - 3. * pc - 3. / 5. * decdh) * (dh / (ec + 3. * pc)));
  m = 4. / 3. * PI * ec * r * r * r * (1. - 3. / 5. * decdh * dh / ec);

  // get into
  u[0][0] = h;
  u[0][1] = r;
  u[0][2] = m;
  u[0][3] = 0;

  printf("tov_h:  solve TOV star (only once):\n");
  printf("    dh   = %.16e npts = %d\n", dh, npts);
  printf("    hc   = %.16e\n", hc);
  printf("    ec   = %.16e\n", ec);
  printf("    pc   = %.16e\n", pc);

  if (PRNT)
    printf("id: h  = %.16e r    = %.16e m  = %.16e\n", h, r, m);
  if (PRNT)
    printf("it: h r m\n");

  double stp = -dh;

  for (n = 0; n < npts - 1; n++)
  {

    tov_h_rhs(stp, u[n], k);

    // u_1 = u + dt/2 k
    forv u1[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_1)
    tov_h_rhs(stp, u1, k);

    // u_2 = u + dt/2 k
    forv u2[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_2)
    tov_h_rhs(stp, u2, k);

    // u_3 = u + dt k
    forv u3[v] = u[n][v] + stp * k[v];

    // r = rhs(u_3)
    tov_h_rhs(stp, u3, k);

    // u = 1/6 ( -2 u + 2 u_1 + 4 u_2 + 2 u_3 + dt k )
    forv u[n + 1][v] = fact * (2. * (-u[n][v] + u1[v] + u3[v]) + 4. * u2[v] + stp * k[v]);

    u[n + 1][0] += stp;

    if (PRNT)
      printf("%d %.16e %.16e %.16e %.16e\n", n, u[n][0], u[n][1], u[n][2], u[n][3]);
  }

  hR = u[npts - 1][0];
  r = u[npts - 1][1];
  m = u[npts - 1][2];
  IR = u[npts - 1][3];
  phiR = 0.5 * log(1. - 2. * m / r);

  for (n = 0; n < npts; n++)
  {

    p_h[n] = u[n][0];
    p_r[n] = u[n][1];
    p_m[n] = u[n][2];

    // grav potential
    phi = phiR + hR - p_h[n];

    // eos
    EOS.comp("H", "", "", "pre", "", "", p_h[n], &p, &rho, &epsl);

    p_rho[n] = rho;
    p_pre[n] = p;
    p_phi[n] = phi;

    C = 1 / (2 * r) * (sqrt(r * r - 2 * m * r) + r - m) * exp(-IR);
    p_riso[n] = p_r[n] * C * exp(u[n][3]);

    if (PRNT)
      printf("%d %.16e %.16e %.16e %.16e %.16e | %.16e\n", n, p_r[n], p_m[n], p_phi[n], p_rho[n], p_pre[n], p_riso[n]);
  }

  Mb = computeMb(p_r, p_m, p_rho, npts - 1);

  printf("    R    = %.16e\n", r);
  printf("    M    = %.16e\n", m);
  printf("    Mb   = %.16e\n", Mb);
  if (Getv("BBID_solve", "yes") && !(checkpoint_checkforfiles("_previous")))
  {
    Setd("BBID_rb", r);
  }

  return 1;
}

/* Method 2
   x: r coord */

#undef nvar
#undef forvar
#define nvar 5
#define forv for (v = 0; v < nvar; v++)

/* tov r rhs */
int tov_r_rhs(double dr, double *u, double *k)
{
  double r = u[0];
  double rho = u[1];
  double m = u[2];
  double phi = u[3];
  double I = u[4];

  double p, eps, dpdrho, T;
  if (EOS.comp("r", "", "", "pe", "r", "", rho, &p, &eps, &dpdrho))
  {
    if (PRNT)
      printf("problem in EoS:  %e %e %e %e\n", rho, p, eps, dpdrho);
    return 1;
  }
  double e = rho * (1. + eps);

  if (r == 0)
    r = 1e-10;
  double tmp1 = m + 4. * PI * r * r * r * p;
  double tmp2 = r * r * (1. - 2. * m / r);
  double tmp = (r == 0.) ? 0. : tmp1 / tmp2;

  double drhodr = -(e + p) * tmp / dpdrho;
  double dmdr = 4. * PI * r * r * e;
  double dphidr = tmp;
  double f = sqrt(1. - 2. * m / r);
  double dIdr = (1. - f) / (r * f);

  k[0] = 0;
  k[1] = drhodr;
  k[2] = dmdr;
  k[3] = dphidr;
  k[4] = dIdr;

  return 0;
}

/* tov r solver */
int tov_r(double rhoc, double R, int *npts,
          double **p_r, double **p_m,
          double **p_rho, double **p_pre, double **p_phi,
          double **p_riso)
{

  int n, v, i;

  double **u = (double **)malloc((*npts) * sizeof(double *));
  for (n = 0; n < *npts; n++)
    u[n] = (double *)malloc((nvar) * sizeof(double));
  double u1[nvar], u2[nvar], u3[nvar], k[nvar];
  double fact = 1. / 6.;
  double Tc;

  double stp = R / (*npts);

  double pc, epslc;
  if (EOS.comp("r", "", "", "pe", "", "", rhoc, &pc, &epslc))
    return 1;
  double ec = rhoc * (1. + epslc);

  // get into
  u[0][0] = 0;
  u[0][1] = rhoc;
  u[0][2] = 0;
  u[0][3] = 0;
  u[0][4] = 0;

  printf("tov_r: solve TOV star (only once):\n");
  printf("    drho = %.16e npts = %d\n", stp, *npts);
  printf("    rhoc = %.16e\n", rhoc);
  printf("    ec   = %.16e\n", ec);
  printf("    pc   = %.16e\n", pc);

  double rhoo = u[0][1];
  int stop = 0;
  n = 0;

  double rhomin;

  if (EOS.type == TAB1D || EOS.type == TAB1DHOT)
  {
    rhomin = EOS.rhomin;
  }
  else if (EOS.type == TAB3D)
    rhomin = EOS.rhomin1D;
  else
    rhomin = 0.;

  while (u[n][1] > rhomin && u[n][1] <= 1.01 * rhoo && stop == 0)
  {

    stop += tov_r_rhs(stp, u[n], k);

    // u_1 = u + dt/2 k
    forv u1[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_1)
    stop += tov_r_rhs(stp, u1, k);

    // u_2 = u + dt/2 k
    forv u2[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_2)
    stop += tov_r_rhs(stp, u2, k);

    // u_3 = u + dt k
    forv u3[v] = u[n][v] + stp * k[v];

    // r = rhs(u_3)
    stop += tov_r_rhs(stp, u3, k);

    // u = 1/6 ( -2 u + 2 u_1 + 4 u_2 + 2 u_3 + dt k )
    forv u[n + 1][v] = fact * (2. * (-u[n][v] + u1[v] + u3[v]) + 4. * u2[v] + stp * k[v]);

    u[n + 1][0] += stp;

    if (PRNT)
      printf("%d %.16e %.16e %.16e %.16e %.16e    %e %d\n", n, u[n + 1][0], u[n + 1][1], u[n + 1][2], u[n + 1][3], u[n + 1][4], rhoo, nvar);

    if (n >= (*npts) - 5)
    {
      u = (double **)realloc(u, (*npts * 2) * sizeof(double *));
      for (i = *npts; i < *npts * 2; i++)
        u[i] = (double *)malloc((nvar) * sizeof(double));
      *npts = (*npts) * 2;
      // printf("  expand grid\n");
      // errorexit("");
    }
    rhoo = u[n][1];

    n++;
  }

  double p, dpdrho, phi, C, M, Mb, IR, phiR, phiRa;

  *npts = n;
  R = u[*npts - 1][0];
  M = u[*npts - 1][2];
  phiR = u[*npts - 1][3];
  IR = u[*npts - 1][4];
  phiRa = 0.5 * log(1. - 2. * M / R);
  C = 1 / (2 * R) * (sqrt(R * R - 2 * M * R) + R - M) * exp(-IR);

  *p_r = (double *)malloc(*npts * sizeof(double));
  *p_m = (double *)malloc(*npts * sizeof(double));
  *p_rho = (double *)malloc(*npts * sizeof(double));
  *p_pre = (double *)malloc(*npts * sizeof(double));
  *p_phi = (double *)malloc(*npts * sizeof(double));
  *p_riso = (double *)malloc(*npts * sizeof(double));

  for (n = 0; n < *npts; n++)
  {

    (*p_r)[n] = u[n][0];
    (*p_rho)[n] = u[n][1];
    (*p_m)[n] = u[n][2];
    (*p_phi)[n] = (u[n][3] - phiR + phiRa);
    // EOS.comp("r","","","p","r","", (*p_rho)[n], &p,&dpdrho);
    //(*p_pre)[n] = p;
    (*p_riso)[n] = (*p_r)[n] * C * exp(u[n][4]);
  }

  Mb = computeMb(*p_r, *p_m, *p_rho, *npts - 1);

  printf("    R    = %.16e   (%.16e)\n", R, (*p_riso)[*npts - 1]);
  printf("    M    = %.16e\n", M);
  printf("    Mb   = %.16e\n", Mb);

  // free memory
  for (i = 0; i < *npts; i++)
    free(u[i]);
  free(u);

  return 0;
}

/* hg: TOV solver using pressure, tov_p_rhs and tov_p with Lambda computation*/

#undef nvar
#undef forv
#define nvar 7
#define forv for (v = 0; v < nvar; v++)

/* tov p rhs */
int tov_p_rhs(double dr, double *u, double *k)
{
  double r = u[0];
  double p = u[1];
  double m = u[2];
  double phi = u[3];
  double I = u[4];
  double H = u[5];
  double G = u[6];

  double e, rho, dedp;
  if (EOS.comp("p", "", "", "re", "e", "", p, &rho, &e, &dedp))
  {
    if (PRNT)
      printf("problem in EoS:  %e %e %e %e\n", rho, p, e, dedp);
    return 1;
  }

  if (r == 0)
  {
    r = 1e-10;
    H = u[5] = r * r;
    G = u[6] = 2. * r;
  }
  double tmp1 = m + 4. * PI * r * r * r * p;
  double tmp2 = r * r * (1. - 2. * m / r);
  double tmp = (r == 0.) ? 0. : tmp1 / tmp2;

  double dpdr = -(e + p) * tmp;
  double dmdr = 4. * PI * r * r * e;
  double dphidr = tmp;
  double f = sqrt(1. - 2. * m / r);
  double dIdr = (1. - f) / (r * f);
  double dnudr = -2. * dpdr / (e + p);

  double a = 2. / r + r * r / tmp2 * (2. * m / (r * r) + 4. * M_PI * r * (p - e));
  double b = r * r / tmp2 * (-6. / (r * r) + 4. * M_PI * (5. * e + 9. * p) + 4. * M_PI * dedp * (p + e)) - dnudr * dnudr;

  double dGdr = -a * G - b * H;
  double dHdr = G;

  k[0] = 0;
  k[1] = dpdr;
  k[2] = dmdr;
  k[3] = dphidr;
  k[4] = dIdr;
  k[5] = dHdr;
  k[6] = dGdr;

  return 0;
}

/* tov p solver */
int tov_p(double pc, double R, int *npts,
          double **p_r, double **p_m,
          double **p_rho, double **p_pre, double **p_phi,
          double **p_riso, double *Lamb)
{

  int n, v, i;

  double **u = (double **)malloc((*npts) * sizeof(double *));
  for (n = 0; n < *npts; n++)
    u[n] = (double *)malloc((nvar) * sizeof(double));
  double u1[nvar], u2[nvar], u3[nvar], k[nvar];
  double fact = 1. / 6.;

  double stp = R / (*npts);

  double rhoc, ec, dummy;
  if (EOS.comp("p", "", "", "re", "", "", pc, &rhoc, &ec))
    return 1;

  // get into
  u[0][0] = 0.;
  u[0][1] = pc;
  u[0][2] = 0.;
  u[0][3] = 0.;
  u[0][4] = 0.;
  u[0][5] = 0.;
  u[0][6] = 0.;

  printf("tov_p: solve TOV star (only once):\n");
  printf("    dp = %.16e npts = %d\n", stp, *npts);
  printf("    pc = %.16e\n", pc);
  printf("    ec   = %.16e\n", ec);
  printf("    rhoc   = %.16e\n", rhoc);

  double poo = u[0][1];
  int stop = 0;
  n = 0;

  double pmin;

  if (EOS.type == TAB1D || EOS.type == TAB1DHOT)
    EOS.comp("r", "", "", "pe", "", "", EOS.rhomin, &pmin, &dummy);

  while (u[n][1] > pmin && u[n][1] <= 1.01 * poo && stop == 0)
  {

    stop += tov_p_rhs(stp, u[n], k);

    // u_1 = u + dt/2 k
    forv u1[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_1)
    stop += tov_p_rhs(stp, u1, k);

    // u_2 = u + dt/2 k
    forv u2[v] = u[n][v] + 0.5 * stp * k[v];

    // r = rhs(u_2)
    stop += tov_p_rhs(stp, u2, k);

    // u_3 = u + dt k
    forv u3[v] = u[n][v] + stp * k[v];

    // r = rhs(u_3)
    stop += tov_p_rhs(stp, u3, k);

    // u = 1/6 ( -2 u + 2 u_1 + 4 u_2 + 2 u_3 + dt k )
    forv u[n + 1][v] = fact * (2. * (-u[n][v] + u1[v] + u3[v]) + 4. * u2[v] + stp * k[v]);

    u[n + 1][0] += stp;

    if (PRNT)
      printf("%d %.16e %.16e %.16e %.16e %.16e    %e %d\n", n, u[n + 1][0], u[n + 1][1], u[n + 1][2], u[n + 1][3], u[n + 1][4], poo, nvar);

    if (n >= (*npts) - 5)
    {
      u = (double **)realloc(u, (*npts * 2) * sizeof(double *));
      for (i = *npts; i < *npts * 2; i++)
        u[i] = (double *)malloc((nvar) * sizeof(double));
      *npts = (*npts) * 2;
      // printf("  expand grid\n");
      // errorexit("");
    }
    poo = u[n][1];

    n++;
  }

  double rho, dpdrho, phi, C, M, Mb, IR, phiR, phiRa, H, G;

  *npts = n;
  R = u[*npts - 1][0];
  M = u[*npts - 1][2];
  phiR = u[*npts - 1][3];
  IR = u[*npts - 1][4];
  phiRa = 0.5 * log(1. - 2. * M / R);
  C = 1 / (2 * R) * (sqrt(R * R - 2 * M * R) + R - M) * exp(-IR);
  H = u[*npts - 1][5];
  G = u[*npts - 1][6];

  *p_r = (double *)malloc(*npts * sizeof(double));
  *p_m = (double *)malloc(*npts * sizeof(double));
  *p_rho = (double *)malloc(*npts * sizeof(double));
  *p_pre = (double *)malloc(*npts * sizeof(double));
  *p_phi = (double *)malloc(*npts * sizeof(double));
  *p_riso = (double *)malloc(*npts * sizeof(double));

  for (n = 0; n < *npts; n++)
  {

    (*p_r)[n] = u[n][0];
    (*p_pre)[n] = u[n][1];
    (*p_m)[n] = u[n][2];
    (*p_phi)[n] = (u[n][3] - phiR + phiRa);
    EOS.comp("p", "", "", "re", "", "", (*p_pre)[n], &rho, &dummy);
    (*p_rho)[n] = rho;
    (*p_riso)[n] = (*p_r)[n] * C * exp(u[n][4]);
  }

  *Lamb = computeLambda(M, R, G, H);
  Mb = computeMb(*p_r, *p_m, *p_rho, *npts - 1);

  printf("    R    = %.16e   (%.16e)\n", R, (*p_riso)[*npts - 1]);
  printf("    M    = %.16e\n", M);
  printf("    Mb   = %.16e\n", Mb);
  printf("    Lamb = %.16e\n", *Lamb);

  // free memory
  for (i = 0; i < *npts; i++)
    free(u[i]);
  free(u);

  return 0;
}

/* test function
   use tov_r ! */
int test_tov_MvsR(tL *level)
{

  static double *p_r, *p_m, *p_h, *p_rho, *p_pre, *p_phi, *p_riso;

  double rho0 = Getd("BBID_tov_rho0");
  double drho = Getd("BBID_tov_drho");
  double R0 = Getd("BBID_tov_R0");
  int Npts = Geti("BBID_tov_npts");
  int N = Geti("BBID_tov_N");

  int n = Npts;

  if (EOS.type == TAB1D || EOS.type == TAB1DHOT)
  {
    drho = pow(EOS.rhomax / rho0, (1. / (double)N));
  }


  // Table:
  //   0      1    2     3    4   5    6
  // rho_c, R_iso, M, status, Mb, R, R_inv
  double data[7][N];
  for (int i = 0; i < N; i++)
  {

    data[0][i] = rho0 * pow(drho, i);

    data[3][i] = tov_r(data[0][i], R0, &n,
                       &p_r, &p_m, &p_rho,
                       &p_pre, &p_phi,
                       &p_riso);

    data[1][i] = p_riso[n - 1];
    data[2][i] = p_m[n - 1];

    // rest mass
    data[4][i] = computeMb(p_r, p_m, p_rho, n - 1);

    // radius
    data[5][i] = p_r[n - 1];

    // invariant radius
    data[6][i] = computeRinv(p_r, p_m, n - 1);

    free(p_r);
    free(p_m);
    free(p_rho);
    free(p_pre);
    free(p_phi);
    free(p_riso);

    n = Npts;
  }

  char file[10000];
  sprintf(file, "%s/%s", Gets("outdir"), Gets("BBID_tov_file"));
  if (processor0)
  {
    FILE *ofp = fopen(file, "w");
    if (!ofp)
      errorexit("cannot write file");
    fprintf(ofp, "# MvsR for %s %s\n", Gets("eos"), Gets("eos_tab_file"));
    if (Getv("eos", "poly"))
    {
      fprintf(ofp, "#   Gamma = %e\n", EOS.GAMMA);
      fprintf(ofp, "#   K     = %e\n", EOS.K);
    }
    // Write field description header
    fprintf(ofp, "# rho_c, R_iso, M, Mb, R, ec, R_inv\n");
    for (int i = 0; i < N; i++)
    {
      int status = data[3][i];
      if (status != 0)
      {
        continue;
      }
      double rho_c = data[0][i];
      double R_iso = data[1][i];
      double M = data[2][i];
      double Mb = data[4][i];
      double R = data[5][i];
      double R_inv = data[6][i];
      double epsl_c;

      EOS.comp("r", "", "", "e", "", "", rho_c, &epsl_c);
      double e_c = rho_c * (epsl_c + 1.);

      fprintf(ofp, " %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          rho_c, R_iso, M, Mb, R, e_c, R_inv);
    }
    fclose(ofp);
  }

  errorexit("stop here, output is ready");
}

int test_tov_MvsR_p(tL *level)
{

  static double *p_r, *p_m, *p_h, *p_rho, *p_pre, *p_phi, *p_riso;

  double p0 = Getd("BBID_tov_p0");
  double dp = Getd("BBID_tov_dp");
  double R0 = Getd("BBID_tov_R0");
  int Npts = Geti("BBID_tov_npts");
  int N = Geti("BBID_tov_N");
  double pmax, dummy;

  int n = Npts;

  double p2, p3;

  if (EOS.type == TAB1D || EOS.type == TAB1DHOT)
  {
    EOS.comp("r", "", "", "p", "", "", EOS.rhomax, &pmax);
    dp = pow(pmax / p0, (1. / (double)N));
  }
  else
    dp = pow(5e-03 / p0, (1. / (double)N));

  // Table:
  //  0     1    2     3    4   5    6       7
  // p_c, R_iso, M, status, Mb, R, Lambda, R_inv
  double data[8][N];
  double Lamb;
  for (int i = 0; i < N; i++)
  {

    data[0][i] = p0 * pow(dp, i);

    data[3][i] = tov_p(data[0][i], R0, &n,
                       &p_r, &p_m, &p_rho,
                       &p_pre, &p_phi,
                       &p_riso, &Lamb);

    data[1][i] = p_riso[n - 1];
    data[2][i] = p_m[n - 1];

    // rest mass
    data[4][i] = computeMb(p_r, p_m, p_rho, n - 1);

    // radius
    data[5][i] = p_r[n - 1];

    // Lambda
    data[6][i] = Lamb;

    // invariant radius
    data[7][i] = computeRinv(p_r, p_m, n - 1);

    free(p_r);
    free(p_m);
    free(p_rho);
    free(p_pre);
    free(p_phi);
    free(p_riso);

    n = Npts;
  }

  char file[10000];
  sprintf(file, "%s/%s", Gets("outdir"), Gets("BBID_tov_file"));
  if (processor0)
  {
    FILE *ofp = fopen(file, "w");
    if (!ofp)
      errorexit("cannot write file");
    fprintf(ofp, "# MvsR for %s %s\n", Gets("eos"), Gets("eos_tab_file"));
    if (Getv("eos", "poly"))
    {
      fprintf(ofp, "#   Gamma = %e\n", EOS.GAMMA);
      fprintf(ofp, "#   K     = %e\n", EOS.K);
    }
    // Write field description header
    fprintf(ofp, "# rho_c, R_iso, M, Mb, R, Lambda, e_c, R_inv\n");

    for (int i = 0; i < N; i++)
    {
      if (data[3][i] == 0)
      {
        double rho_c, p_c, e_c;
        p_c = data[0][i];
        EOS.comp("p", "", "", "re", "", "", p_c, &rho_c, &e_c);
        fprintf(ofp, " %e %e %e %e %e %e %e %e\n", rho_c, data[1][i], data[2][i], data[4][i], data[5][i], data[6][i], e_c, data[7][i]);
      }
    }
    fclose(ofp);
  }

  errorexit("stop here, output is ready");
}
