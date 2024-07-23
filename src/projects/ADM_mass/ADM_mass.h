/* Adm_mass.h */
/* Jose Gonzalez, 9/04 */

int ADM_mass(tL *level);

void adm_mass_integrand(tL *level, int i_x, int i_g, int i_integrand);
void adm_mass_integrand_chi(tL *level, int i_x, int i_chi, int i_integrand);
void adm_pj_integrand(tL *level, int i_x, int i_g, int i_K, int i_TrK);





