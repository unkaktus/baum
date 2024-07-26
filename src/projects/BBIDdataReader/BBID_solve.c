
/* BBID_solve.c */
/* nmol 10/13 */

#include "bam.h"
#include "BBIDdataReader.h"
#include <stdio.h>


#define debug 1
#define PR 1

int solve_apsi=0;
int solve_omega=0;
int solve_beta=0;

double n = -1 ; //Rescaling parameter: n = -1 <==> no rescaling

#define VY (xp[ijk]>=0?vy1:vy2)

void DPflatlinear(tL *level, tVarList *v, tVarList *u);

void load_data(tL *level){

	double *variable;


	char cmd[1000],file[1000],line[1024],prev[1000],eos[1000];
	static char prog[1000],data[1000],path[1000],finish[1000];

	FILE *ifp;

	int *index =  malloc (5*sizeof(int));
	index[0] = Ind("BBID_u_omega");
	index[1] = Ind("BBID_u_apsi");
	index[2] = Ind("BBID_u_betax");
	index[3] = Ind("BBID_u_betay");
	index[4] = Ind("BBID_u_betaz");

	char boxstring[2] = {0, 0};

	int i, proci, vari, ijk;
	/* Take care of box a and b */
	for (i = 0; i < level->nboxes; i++) {
		if (level->nboxes == 1) boxstring[0] = 0;
		else boxstring[0] = 'a' + i;
	}

	double tmp;
	tL *level0 = level;
	for (i = 0; i < level0->nboxes; i++) {
		level = one_box_level(level0, i);
		if (level0->nboxes == 1) boxstring[0] = 0;
		else boxstring[0] = 'a' + i;
		
		for(proci = 0; proci < bampi_size(); proci++){
			
			if(proci == bampi_rank()){


				for (vari=0; vari<5; vari++){
				
					sprintf(path,   "%s/%sL%d%s.%d", Gets("BBID_solve_loadDir"),VarName(index[vari]),level->l,boxstring,bampi_rank());
					printf("Reading: %s\n", path);
					if (!system_isfile(path)) {
						printf("can not find saved BBID_solve_ID output ( %s )",path);
// 						errorexit("EXIT");}
					}
					// printf("Reading variable %s\n",path);
					ifp = fopen(path, "r");
					if (!ifp) errorexit("Cannot open file");



					variable	= Ptr(level, VarName(index[vari]));
					ijk=0;
					while(!feof(ifp)){
						fscanf (ifp, "%lf", &tmp); //value
						if (ijk< (level->box[0]->m)*(level->box[0]->n)*(level->box[0]->o))// errorexit("Too many points in saved initial data file");
						variable[ijk]=tmp;
						//printf("%e\n",f);
						ijk+=1;
					}
					
					fclose(ifp);
				}
			}
			bampi_barrier();
			
		}
	}
	

}

void save_data(tL *level){

	double *variable;
	int *index =  malloc (5*sizeof(int));
	index[0] = Ind("BBID_u_omega");
	index[1] = Ind("BBID_u_apsi");
	index[2] = Ind("BBID_u_betax");
	index[3] = Ind("BBID_u_betay");
	index[4] = Ind("BBID_u_betaz");
		
	FILE *ofp;
	char outpath[1000],filename[1000];

	sprintf(outpath,   "%s/BBID_solve_data", Gets("outdir"));

	if(level->l==level->grid->lmin)
		system_mkdir(outpath);
	bampi_barrier();

	int vari, proci, i, ijk;

	char boxstring[2] = {0, 0};

	/* Take care of box a and b */
	tL *level0 = level;
	for (i = 0; i < level0->nboxes; i++) {
		level = one_box_level(level0, i);
		if (level0->nboxes == 1) boxstring[0] = 0;
		else boxstring[0] = 'a' + i;
			
		for(proci = 0; proci < bampi_size(); proci++){


			if(proci == bampi_rank()){
				
				printf("Saving level %d from processor %d \n",level->l,bampi_rank());
				for (vari=0; vari<5; vari++){
					variable	= Ptr(level, VarName(index[vari]));
					sprintf(filename,   "%s/BBID_solve_data/%sL%d%s.%d", Gets("outdir"),VarName(index[vari]),level->l,boxstring,bampi_rank());
					/* open file for writing */
					ofp = fopen(filename, "wb");
						
					ijk=0;
					
					if (!ofp) errorexits("failed opening %s", filename);
					
					
					while(ijk<(level->box[0]->m)*(level->box[0]->n)*(level->box[0]->o)){
						fprintf(ofp,"%.15f\n",variable[ijk]);
						ijk+=1;
					}
					fclose(ofp);
				}
			}
			bampi_barrier();
		}
	}
}

double massintegral(tL *level, int star){

	double local_sum=0, integral=0;
	
	double *eta 	= Ptr(level, "BBID_eta");
	double *psiq 	= Ptr(level, "BBID_psi");
	double *alpha 	= Ptr(level, "alpha");
	double *ut 	= Ptr(level, "BBID_ut");
	double *xp	= Ptr(level, "x");
	//double eta0  	= exp(Getd("BBID_tov_hc"))-1.; 						
	double K	= EOS.K; 
	double polytrope= 1./EOS.GAMMAMO;
	
	double fppp,fppc,fpcp,fpcc,fcpp,fcpc,fccp,fccc;
	
	double rho0;
	
	
	double dx = level->dx;
	double dy = level->dy;
	double dz = level->dz;

	
	
	forinner27(level){
		if((star==0 && xp[ccc]>0) || (star==1 && xp[ccc]<0)){
			fppp = pow(eta[ppp]/(K*(1.+polytrope)),polytrope)*ut[ppp]*alpha[ppp]*pow(psiq[ppp]+1,6);
			fppc = pow(eta[ppc]/(K*(1.+polytrope)),polytrope)*ut[ppc]*alpha[ppc]*pow(psiq[ppc]+1,6);
			fpcp = pow(eta[pcp]/(K*(1.+polytrope)),polytrope)*ut[pcp]*alpha[pcp]*pow(psiq[pcp]+1,6);
			fpcc = pow(eta[pcc]/(K*(1.+polytrope)),polytrope)*ut[pcc]*alpha[pcc]*pow(psiq[pcc]+1,6);
			fcpp = pow(eta[cpp]/(K*(1.+polytrope)),polytrope)*ut[cpp]*alpha[cpp]*pow(psiq[cpp]+1,6);
			fcpc = pow(eta[cpc]/(K*(1.+polytrope)),polytrope)*ut[cpc]*alpha[cpc]*pow(psiq[cpc]+1,6);
			fccp = pow(eta[ccp]/(K*(1.+polytrope)),polytrope)*ut[ccp]*alpha[ccp]*pow(psiq[ccp]+1,6);
			fccc = pow(eta[ccc]/(K*(1.+polytrope)),polytrope)*ut[ccc]*alpha[ccc]*pow(psiq[ccc]+1,6);
			
			local_sum+=fppp+fppc+fpcp+fpcc+fcpp+fcpc+fccp+fccc;
		}
		
	}endfor;
	
	bampi_allreduce_sum_vector(&local_sum, &integral, 1);
	integral*=dx*dy*dz/8;
	

	return(integral);
	
}

/*
 * Computes the constants Epsilon, Omega and r_B out of the first Integral of Euler's eq
 * evaluated at three points (2 intersections of star and x-axis + 1 center of star)
 * 
 */
void find_constants(tL *level, double Omega, double Econst_0, double Econst_1){ 	//Current solution xn

//	printf("║\n║ Computing the constants with Newton-Raphson on level %d\n",level->l);
	//Define Delta_xn which updates xn
	double del_rb, del_Omega, del_Econst_0, del_Econst_1;
	double etaC0	= exp(Getd("BBID_tov_hc1"))-1; 	//Defined density
	double etaC1	= exp(Getd("BBID_tov_hc2"))-1; 	//Defined density
	
	double C0, C1;
	double e = Getd("BBID_solve_ecc");		//eccentricity
	double xcom = Getd("BBID_solve_com");		//center of mass
	
	double *psiq	= Ptr(level, "BBID_psi");
	double *alpha	= Ptr(level, "alpha");
	double pnew[3], p0old[3], p1old[3] = {0,0,0};
	double V   =  - Getd("my1");
	int ii;
	double betaxC0, betayC0, betazC0, alphaC0, psiC0, beta2C0, utC0;
	double betaxC1, betayC1, betazC1, alphaC1, psiC1, beta2B1, utC1;
	double DalphaC0, DpsiC0, DbetaxC0, DbetayC0, DbetazC0, DutC0;
	double DalphaC1, DpsiC1, DbetaxC1, DbetayC1, DbetazC1, DutC1;
	double dx = level->dx;

	double m01 = Getd("BBID_solve_M_central"); 	
	
	//define functions f=0 to solve
	double fA, fB, fC;

	p0old[0] = level->grid->puncpos[0][0];
	p1old[0] = level->grid->puncpos[1][0];
	
	
	
	

	//track_moving_puncture_extremum(level, 0, p0old, pnew);
	C0 = p0old[0];
		
	//track_moving_puncture_extremum(level, 0, p1old, pnew);
	C1 = p1old[0];
	
	
	
	//printf("╠═══ Intersection and central points on x-axis C %f\n",C0);		

	
	bampi_synchronize(level, Ind("BBID_u_betax"));
	bampi_synchronize(level, Ind("BBID_psi"));
	bampi_synchronize(level, Ind("alpha"));
	bampi_synchronize(level, Ind("BBID_ut"));
	
	
	betaxC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE); 
	betayC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE); 
	betazC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE); 
	alphaC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE); 
	psiC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE)+1; 	
	betaxC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE); 
	betayC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE); 
	betazC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE); 
	alphaC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE); 
	psiC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE)+1; 	
	
	utC0	= interpolate_xyz_scalar(level, C0, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE); 
	utC1	= interpolate_xyz_scalar(level, C1, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE); 
	
	DalphaC0 = 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE));
	DpsiC0 	= 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE));

	DalphaC1 = 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("alpha"), Geti("order_RP"),LAGRANGE));
	DpsiC1 	= 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("BBID_psi"), Geti("order_RP"),LAGRANGE));

	
	DutC0 	= 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE));
	DutC1 	= 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("BBID_ut"), Geti("order_RP"),LAGRANGE));

	DbetaxC0 = 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE));
	DbetayC0 = 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE));
	DbetazC0 = 0.5/dx*(interpolate_xyz_scalar(level, C0+dx, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C0-dx, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE));

	DbetaxC1 = 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("BBID_u_betax"), Geti("order_RP"),LAGRANGE));
	DbetayC1 = 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("BBID_u_betay"), Geti("order_RP"),LAGRANGE));
	DbetazC1 = 0.5/dx*(interpolate_xyz_scalar(level, C1+dx, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE) - interpolate_xyz_scalar(level, C1-dx, 0, 0, Ind("BBID_u_betaz"), Geti("order_RP"),LAGRANGE));

	

	
	double xC0=C0;

	double xC1=C1;

		
	
	if(!Getv("BBID_modus","NSNS")){
		del_Omega =0; Omega =0;
		
		del_Econst_0 = - Econst_0 + (1+etaC0)*sqrt(alphaC0*alphaC0);
		del_Econst_1 = - Econst_1 + (1+etaC1)*sqrt(alphaC1*alphaC1);
	}else if(Getv("BBID_solve_spin","corotational")){
		
		///EQUAL MASS
// 		printf("╠═══ Using corotational configuration\n");					
// 		del_Omega = - Omega + (-(betayC0*pow(psiC0,4) + 4*betayC0*DpsiC0*pow(psiC0,3)*C0 + DbetayC0*pow(psiC0,4)*C0 - 
// 			sqrt(-4*pow(psiC0,3)*(-(alphaC0*DalphaC0) + pow(psiC0,3)*(2*pow(betaxC0,2)*DpsiC0 + 2*pow(betayC0,2)*DpsiC0 + 2*pow(betazC0,2)*DpsiC0 + betaxC0*DbetaxC0*psiC0 + 
// 			betayC0*DbetayC0*psiC0 + betazC0*DbetazC0*psiC0))*C0*(psiC0 + 2*DpsiC0*C0) + pow(psiC0,6)*pow(DbetayC0*psiC0*C0 + betayC0*(psiC0 + 4*DpsiC0*C0),2)))/
// 			(2.*pow(psiC0,3)*C0*(psiC0 + 2*DpsiC0*C0)));
// 			
// 			
// 		Omega += del_Omega;
// 		
// 		if(Getv("BBID_solve_fix_at_center","M")){
// 	
// 			double *eta 	= Ptr(level, "BBID_eta"); 
// 			printf("\nMassintegral:before %f - Central Density: %f\n", massintegral(level), eta0);
// 			double fak = m01/massintegral(level);
// 			eta0 *= fak;
// 			forallpoints_ijk(level) { 
// 				eta[ijk] *=fak;
// 			} endfor_ijk;
// 			
// 			Setd("BBID_tov_hc",log(eta0+1));
// 		}
// 		
// 		del_Econst_0 = - Econst_0 + (1+eta0)*sqrt(alphaC0*alphaC0-pow(psiC0,4)*pow(Omega*C0+betayC0,2));
// 		del_Econst_1 = - Econst_1 + (1+eta0)*sqrt(alphaC1*alphaC1-pow(psiC1,4)*pow(Omega*C1+betayC1,2));
		//printf("╠═══ Using corotational configuration\n");					
		int NR_loop;
				
		double fx, fy;
		//define coefficients of Jacobian of f
		double Jxx, Jxy, Jyx, Jyy;
		//define coefficients of inverse Jacobian of f
		double Ixx, Ixy, Iyx, Iyy;
		for(NR_loop =0; NR_loop<15; NR_loop++){
			
			fx = ((2*alphaC0*DalphaC0 - 2*pow(psiC0,4)*(betaxC0*DbetaxC0 + betazC0*DbetazC0 + (DbetayC0 + Omega)*(betayC0 + Omega*(xC0 - xcom))) - 
				4*DpsiC0*pow(psiC0,3)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)))*(1+etaC0))/
				(2.*sqrt(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2))));
			fy = ((2*alphaC1*DalphaC1 - 2*pow(psiC1,4)*(betaxC1*DbetaxC1 + betazC1*DbetazC1 + (DbetayC1 + Omega)*(betayC1 + Omega*(xC1 - xcom))) - 
				4*DpsiC1*pow(psiC1,3)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)))*(1+etaC1))/
				(2.*sqrt(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2))));
			
			Jxx = ((-2*pow(psiC0,4)*(betayC0 + Omega*(xC0 - xcom) + (DbetayC0 + Omega)*(xC0 - xcom)) - 8*DpsiC0*pow(psiC0,3)*(betayC0 + Omega*(xC0 - xcom))*(xC0 - xcom))*(1+etaC0))/
				(2.*sqrt(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)))) + 
				(pow(psiC0,4)*(2*alphaC0*DalphaC0 - 2*pow(psiC0,4)*(betaxC0*DbetaxC0 + betazC0*DbetazC0 + (DbetayC0 + Omega)*(betayC0 + Omega*(xC0 - xcom))) - 
				4*DpsiC0*pow(psiC0,3)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)))*(betayC0 + Omega*(xC0 - xcom))*(xC0 - xcom)*(1+etaC0))/
				(2.*pow(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)),1.5));
			
			Jxy = ((2*Omega*(DbetayC0 + Omega)*pow(psiC0,4) + 8*DpsiC0*Omega*pow(psiC0,3)*(betayC0 + Omega*(xC0 - xcom)))*(1+etaC0))/
				(2.*sqrt(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)))) - 
				(Omega*pow(psiC0,4)*(2*alphaC0*DalphaC0 - 2*pow(psiC0,4)*(betaxC0*DbetaxC0 + betazC0*DbetazC0 + (DbetayC0 + Omega)*(betayC0 + Omega*(xC0 - xcom))) - 
				4*DpsiC0*pow(psiC0,3)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)))*(betayC0 + Omega*(xC0 - xcom))*(1+etaC0))/
				(2.*pow(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betazC0,2) + pow(betayC0 + Omega*(xC0 - xcom),2)),1.5));
			
			Jyx =((-2*pow(psiC1,4)*(betayC1 + Omega*(xC1 - xcom) + (DbetayC1 + Omega)*(xC1 - xcom)) - 8*DpsiC1*pow(psiC1,3)*(betayC1 + Omega*(xC1 - xcom))*(xC1 - xcom))*(1+etaC1))/
				(2.*sqrt(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)))) + 
				(pow(psiC1,4)*(2*alphaC1*DalphaC1 - 2*pow(psiC1,4)*(betaxC1*DbetaxC1 + betazC1*DbetazC1 + (DbetayC1 + Omega)*(betayC1 + Omega*(xC1 - xcom))) - 
				4*DpsiC1*pow(psiC1,3)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)))*(betayC1 + Omega*(xC1 - xcom))*(xC1 - xcom)*(1+etaC1))/
				(2.*pow(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)),1.5));
			Jyy = ((2*Omega*(DbetayC1 + Omega)*pow(psiC1,4) + 8*DpsiC1*Omega*pow(psiC1,3)*(betayC1 + Omega*(xC1 - xcom)))*(1+etaC1))/
				(2.*sqrt(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)))) - 
				(Omega*pow(psiC1,4)*(2*alphaC1*DalphaC1 - 2*pow(psiC1,4)*(betaxC1*DbetaxC1 + betazC1*DbetazC1 + (DbetayC1 + Omega)*(betayC1 + Omega*(xC1 - xcom))) - 
				4*DpsiC1*pow(psiC1,3)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)))*(betayC1 + Omega*(xC1 - xcom))*(1+etaC1))/
				(2.*pow(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betazC1,2) + pow(betayC1 + Omega*(xC1 - xcom),2)),1.5));
				
			Ixx= 1/(Jxx*Jyy-Jxy*Jyx)*Jyy;
			Ixy= -1/(Jxx*Jyy-Jxy*Jyx)*Jxy;
			Iyx= -1/(Jxx*Jyy-Jxy*Jyx)*Jyx;
			Iyy= 1/(Jxx*Jyy-Jxy*Jyx)*Jxx;
				
			del_Omega=-(Ixx*fx+Ixy*fy);
			Omega+=del_Omega;
			xcom -=Iyx*fx+Iyy*fy;
// 			printf("Adjusting Omega (delta %9.19f) and Center of Mass (%f) \n",Omega,xcom);
			if(del_Omega < 1e-12){
// 				printf("Convergence reached in Newton-Raphson\n");
				break;
			}
			
		} if(NR_loop >= 14) printf("WARNING: Convergenge NOT reached in Newton-Raphson. This is probably going to fail!\n");
		del_Econst_0 = - Econst_0 + (1+etaC0)*sqrt(alphaC0*alphaC0-pow(psiC0,4)*pow(Omega*(xC0-xcom)+betayC0,2)); 
		del_Econst_1 = - Econst_1 + (1+etaC1)*sqrt(alphaC1*alphaC1-pow(psiC1,4)*pow(Omega*(xC1-xcom)+betayC1,2));
		
		
	}else if(Getv("BBID_solve_spin","irrotational")){
	
		//printf("╠═══ Using irrotational configuration\n");
		del_Omega = - Omega  + (-(betayC0*pow(psiC0,4)*utC0) - betayC0*DutC0*pow(psiC0,4)*C0 - betayC0*DutC0*sqrt(1 - e)*pow(psiC0,4)*C0 - 4*betayC0*DpsiC0*pow(psiC0,3)*utC0*C0 - 
			4*betayC0*DpsiC0*sqrt(1 - e)*pow(psiC0,3)*utC0*C0 - DbetayC0*pow(psiC0,4)*utC0*C0 - DbetayC0*sqrt(1 - e)*pow(psiC0,4)*utC0*C0 + 
			sqrt(-4*sqrt(1 - e)*pow(psiC0,3)*(-(pow(alphaC0,2)*DutC0) - 2*alphaC0*DalphaC0*utC0 + 
			pow(psiC0,3)*(2*betaxC0*DbetaxC0*psiC0*utC0 + 2*betayC0*DbetayC0*psiC0*utC0 + pow(betaxC0,2)*(DutC0*psiC0 + 4*DpsiC0*utC0) + pow(betayC0,2)*(DutC0*psiC0 + 4*DpsiC0*utC0) + 
			betazC0*(betazC0*DutC0*psiC0 + 4*betazC0*DpsiC0*utC0 + 2*DbetazC0*psiC0*utC0)))*C0*(DutC0*psiC0*C0 + utC0*(psiC0 + 4*DpsiC0*C0)) + 
			pow(psiC0,6)*pow(betayC0*DutC0*(1 + sqrt(1 - e))*psiC0*C0 + utC0*(DbetayC0*(1 + sqrt(1 - e))*psiC0*C0 + betayC0*(psiC0 + 4*DpsiC0*(1 + sqrt(1 - e))*C0)),2)))/
			(2.*sqrt(1 - e)*pow(psiC0,3)*C0*(DutC0*psiC0*C0 + utC0*(psiC0 + 4*DpsiC0*C0)));
		Omega += del_Omega;
		
		del_Econst_0 = - Econst_0 + ((1 + etaC0)*utC0*(pow(alphaC0,2) - pow(psiC0,4)*(pow(betaxC0,2) + pow(betayC0,2) + pow(betazC0,2) + betayC0*sqrt(1 - e)*Omega*C0 + Omega*C0*(betayC0 + sqrt(1 - e)*Omega*C0))));
		del_Econst_1 = - Econst_1 + ((1 + etaC1)*utC1*(pow(alphaC1,2) - pow(psiC1,4)*(pow(betaxC1,2) + pow(betayC1,2) + pow(betazC1,2) + betayC1*sqrt(1 - e)*Omega*C1 + Omega*C1*(betayC1 + sqrt(1 - e)*Omega*C1))));
		
	}else errorexit("Specify co- or irrotational stars!");
	
	
	Econst_0 += del_Econst_0;
	Econst_1 += del_Econst_1;
	
// 	if(debug) printf("Star 0: alpha %f,  psi %f,  betax %.10f,  betay %f,  betaz %f,  utC %f\n" , alphaC0, psiC0, betaxC0, betayC0, betazC0, utC0);
// 	if(debug) printf("        Dalpha %f, Dpsi %f, Dbetax %.10f, Dbetay %f, Dbetaz %f, DutC %f\n" , DalphaC0, DpsiC0, DbetaxC0, DbetayC0, DbetazC0, DutC0);
// 	if(debug) printf("Star 1: alpha %f,  psi %f,  betax %.10f,  betay %f,  betaz %f,  utC %f\n" , alphaC1, psiC1, betaxC1, betayC1, betazC1, utC1);
// 	if(debug) printf("        Dalpha %f, Dpsi %f, Dbetax %.10f, Dbetay %f, Dbetaz %f, DutC %f\n" , DalphaC1, DpsiC1, DbetaxC1, DbetayC1, DbetazC1, DutC1);
// 	printf("Quantitites %f, %f", Econst_0/sqrt(alphaC0*alphaC0-pow(psiC0,4)*pow((Omega*C0+betayC0),2))-1 ,(1+etaC0)*sqrt(alphaC0*alphaC0-pow(psiC0,4)*pow((Omega*C0+betayC0),2)));
// 	
// 	printf("\n╠═══ Massintegral: %f / %f - Central Density: %f / %f - Center of Mass: %f\n", massintegral(level,0),massintegral(level,1), etaC0, etaC1, xcom);
// 
// 
// 	printf("╠═══ Omega=%.6e, Econst0=%.6e, Econst1=%.6e \n", Omega, Econst_0, Econst_1);
// 	
	Setd("BBID_solve_com",xcom);
	Setd("BBID_solve_Omega",Omega);
	Setd("BBID_Econst_0", Econst_0);
	Setd("BBID_Econst_1", Econst_1);
	//printf("╠═══ Change in Omega: %.5f %%  and in Econst: %.5f %%\n", 100*del_Omega/Omega , 50*(del_Econst_0+del_Econst_1)/Econst_0);
}

void compute_eta(tL *level) //also computes ut
{	//printf("║\n║ Updating density distribution on level %d\n",level->l);


	double *psiq	= Ptr(level, "BBID_psi");
	double *omegaq	= Ptr(level, "BBID_u_omega");
	double *apsiq	= Ptr(level, "BBID_u_apsi");
	double *betax	= Ptr(level, "BBID_u_betax");
	double *betay	= Ptr(level, "BBID_u_betay");
	double *betaz	= Ptr(level, "BBID_u_betaz");
	double *alpha 	= Ptr(level, "alpha");
	double *ut	= Ptr(level, "BBID_ut");
	double *eta	= Ptr(level, "BBID_eta");
	
	double vy1   =  - Getd("my1");
	double vy2   =  - Getd("my2");
	int tmp, beta2;
	double EconstC0, EconstC1;
	
	double pos[2][3]; 		//pos[i][j] position of star i at coordinate j
	
	double e = Getd("BBID_solve_ecc");
	double xcom = Getd("BBID_solve_com");
	
	double Omega	= 0;
	if (Getd("BBID_solve_Omega")!=0) Omega=Getd("BBID_solve_Omega");
		else if(level->grid->puncpos[0][0] != 0 ) Omega=vy1/(pos[1][0]); //setting the orbital angular frequency
	
		
	pos[0][0] = level->grid->puncpos[0][0];
	pos[1][0] = level->grid->puncpos[1][0];
	
	double R;
	
	double *xp	= Ptr(level, "x");
	double *yp	= Ptr(level, "y");	
	double *zp	= Ptr(level, "z");	
	double rb 	= Getd("BBID_rb");

	
	double V;

	EconstC0 = Getd("BBID_Econst_0"); 
	EconstC1 = Getd("BBID_Econst_1"); 	
	
		
	double soft = Getd("BBID_solve_soft_eta"); 
	
	
	if(Getv("BBID_solve_spin","corotational")){
		forallpoints_ijk(level) {     
			beta2 = betax[ijk]*betax[ijk]+betay[ijk]*betay[ijk]+betaz[ijk]*betaz[ijk]; //beta^i*beta^i    
			//alpha[ijk]= (1+apsiq[ijk])/pow(1+psiq[ijk],1);
			
			eta[ijk] = soft*((xp[ijk]>=0 ? EconstC0:EconstC1)/sqrt(pow(alpha[ijk],2)-pow(1+psiq[ijk],4)*(pow(-Omega*yp[ijk]+betax[ijk],2)+pow(Omega*(xp[ijk]-xcom)+betay[ijk],2)+pow(betaz[ijk],2)))-1)+(1-soft)*eta[ijk];
			
			if(eta[ijk]<0 || xp[ijk] > (rb+2) || xp[ijk]< - (rb+2) || yp[ijk] > (rb+2) || yp[ijk] < -(rb+2) || zp[ijk] > (rb+2) || zp[ijk] < -(rb+2)) eta[ijk]=0;
			ut[ijk] = (1+eta[ijk])/(xp[ijk]>=0 ? EconstC0:EconstC1);
			
			if(debug) if (!finite(eta[ijk])) {printf("HERE %f,%f,%f,%f",eta[ijk],ut[ijk],psiq[ijk],apsiq[ijk]);errorexit("something bad happend with eta");}
		} endfor_ijk;
	}
	else if(Getv("BBID_solve_spin","irrotational")){
		
		forallpoints_ijk(level) {     
			if (xp[ijk]>0) {
				R =level->grid->puncpos[0][0];
			}else {
				R =level->grid->puncpos[1][0];
			}
			
			eta[ijk] = soft*((xp[ijk]>=0 ? EconstC0:EconstC1)/(ut[ijk]*(pow(alpha[ijk],2) - pow(psiq[ijk]+1,4)*(pow(betax[ijk],2) + pow(betay[ijk],2) + pow(betaz[ijk],2) + sqrt(1 - e)*pow(Omega,2)*R*xp[ijk] + 
				Omega*betay[ijk]*(sqrt(1 - e)*R + xp[ijk]) - Omega*betax[ijk]*yp[ijk])))-1)+(1-soft)*eta[ijk];
			
			
			//if(sqrt(pow(xp[ijk]-R,2)+pow(yp[ijk],2)+pow(zp[ijk],2))>rb-abs(R)+2) eta[ijk]=0;
			if(eta[ijk]<0 || xp[ijk] > (rb+2) || xp[ijk]< - (rb+2) || yp[ijk] > (rb+2) || yp[ijk] < -(rb+2) || zp[ijk] > (rb+2) || zp[ijk] < -(rb+2)) eta[ijk]=0;
			//if(eta[ijk]<0 || xp[ijk] > (rb+2) || xp[ijk]< - (rb+2) || yp[ijk] > (rb+2) || yp[ijk] < -(rb+2) || zp[ijk] > (rb+2) || zp[ijk] < -(rb+2)) eta[ijk]=0;
			if(debug) if (!finite(eta[ijk])) {printf("HERE %f,%f,%f,%f",eta[ijk],ut[ijk],psiq[ijk],apsiq[ijk]);errorexit("something bad happend with eta");}
// 			ut[ijk] = (xp[ijk]>=0 ? E_const0a:E_const0b)/((1+eta[ijk])*(pow(alpha[ijk],2) - pow(psiq[ijk]+1,4)*(pow(betax[ijk],2) + pow(betay[ijk],2) + pow(betaz[ijk],2) + sqrt(1 - e)*pow(Omega,2)*R*xp[ijk] + 
// 				Omega*betay[ijk]*(sqrt(1 - e)*R + xp[ijk]) - Omega*betax[ijk]*yp[ijk])));
			
 			ut[ijk] = 1./sqrt((alpha[ijk]*alpha[ijk]-pow(1+psiq[ijk],4)*(pow(betax[ijk],2) + pow(betay[ijk],2) + pow(betaz[ijk],2)+2*Omega*R*sqrt(1-e)*betay[ijk] +R*R*Omega*Omega*(1-e) )));
			V = Omega* R;
// 			fC = -(xp[ijk]>=0 ? E_const0a:E_const0b)- (1 + eta[ijk])*ut[ijk]*(-pow(alpha[ijk],2) + pow(1+psiq[ijk],4)*(betax[ijk]*betax[ijk]+betay[ijk]*betay[ijk]+betaz[ijk]*betaz[ijk] + V*betay[ijk] + Omega*xp[ijk]*(betay[ijk] + V)-Omega*yp[ijk]*betax[ijk]));
// 			if(eta[ijk]>0) printf("fC =%f\n",fC  );
			
			
		} endfor_ijk;
		
	}


	bampi_synchronize(level, Ind("BBID_ut"));
	bampi_synchronize(level, Ind("BBID_eta"));	
	//printf("BBID_Econst %f compared to %f\n",EconstC0, EconstC1);
	/*
	*///printf("╠═══ Ratio of central eta to specified etaC: Star 1... %f, Star 2... %f\n",interpolate_xyz_scalar(level, pos[0][0], 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)/etaC0,interpolate_xyz_scalar(level, pos[1][0], 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)/etaC1);
// 	  printf("╠═══ Derivative at center of Star 1... %f, Star 2... %f\n",	0.5/level->dx*(interpolate_xyz_scalar(level, pos[0][0]+level->dx, 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)-interpolate_xyz_scalar(level, pos[0][0]-level->dx, 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)),
// 										0.5/level->dx*(interpolate_xyz_scalar(level, pos[1][0]+level->dx, 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)-interpolate_xyz_scalar(level, pos[1][0]-level->dx, 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)));

}  /* function */


void compute_sources(tL *level/*, tVarList *vlu, tVarList *vlut, tVarList *vleta, tVarList *vlrhoH, tVarList *vlSi, tVarList *vlS, tVarList *vlAA*/){
	

	
	
	double *psiq	= Ptr(level, "BBID_psi");
	double *omegaq 	= Ptr(level, "BBID_u_omega");
	double *apsiq	= Ptr(level, "BBID_u_apsi");
	double *betax	= Ptr(level, "BBID_u_betax");
	double *betay	= Ptr(level, "BBID_u_betay");
	double *betaz	= Ptr(level, "BBID_u_betaz");
	
	

	double *ut	= Ptr(level, "BBID_ut");
	double *eta	= Ptr(level, "BBID_eta");
	double *S	= Ptr(level, "BBID_S");
	double *Sx	= Ptr(level, "BBID_Six");
	double *Sy	= Ptr(level, "BBID_Siy");
	double *Sz	= Ptr(level, "BBID_Siz");
	double *AA	= Ptr(level, "BBID_AA");
	double *rhoH	= Ptr(level, "BBID_rhoH");
	
	double p, epsilon, alpha;
	
	double eta0  		= exp(Getd("BBID_tov_hc"))-1.; 						
	double K		= EOS.K; 
	double polytrope	= 1./EOS.GAMMAMO;
	double e 		= Getd("BBID_solve_ecc");
	double Omega = Getd("BBID_solve_Omega");
	//double Omega = 0;
	//if(level->grid->puncpos[0][0] != 0 ) Omega=vy1/(level->grid->puncpos[0][0]);

	double rho0=pow(eta0/(K*(1.+polytrope)),polytrope);
	
	double dpsix, dpsiy, dpsiz; //Derivatives of psi
	double dapsix, dapsiy, dapsiz; //Derivatives of apsi
	double dbetaxx,dbetaxy,dbetaxz,dbetayx,dbetayy,dbetayz,dbetazx,dbetazy,dbetazz;
	
	double cx = 1/(level->dx*level->dx);			//for second derivatives
	double cy = 1/(level->dy*level->dy);
	double cz = 1/(level->dz*level->dz);
	double cc = -2.*(cx + cy + cz);

	double cx_1 = 0.5/(level->dx);			    //for first derivatives
	double cy_1=  0.5/(level->dy);
	double cz_1=  0.5/(level->dz);
	double Axx,Ayy, Azz, Axy, Axz, Ayz;
	double *xp = Ptr(level, "x");
	double *yp = Ptr(level, "y");
	double *zp = Ptr(level, "z");
	
	forallpoints_ijk(level) {
		psiq[ijk] = pow(1+omegaq[ijk],-1./n)-1;
	} endfor_ijk;
	
	
	
	bampi_synchronize(level, Ind("BBID_psi"));
	bampi_synchronize(level, Ind("BBID_u_betax"));
	bampi_synchronize(level, Ind("BBID_ut"));
	
	forinner19(level) {
		
		alpha	= (1+apsiq[ccc])/(1+psiq[ccc]);
		
		
		//calculate our matter quantities 

		p       = rho0*pow(eta[ijk]/eta0,polytrope+1)*eta0/(1+polytrope);
		epsilon = rho0*pow(eta[ijk]/eta0,polytrope)*(1+polytrope/(1+polytrope)*eta[ijk]);

		if(debug) if (!finite(p)) {printf("BLUB %f",eta[ijk]); errorexit("something bad happend with p nonlin");}
		if(debug) if (!finite(epsilon)) errorexit(" something bad happend with epsilon");


		// ddxVx = cx*(u_vx[pcc]-2.*u_vx[ccc]+u_vx[mcc]);
		// dxdyVx = cx_1*cy_1*(u_vx[ppc]-u_vx[pmc]-u_vx[mpc]+u_vx[mmc])
		// dzVx = cz_1*(u_vx[ccp]-u_vx[ccm])
			

		dbetaxx = cx_1*(betax[pcc]-betax[mcc]);
		dbetaxy = cy_1*(betax[cpc]-betax[cmc]);
		dbetaxz = cz_1*(betax[ccp]-betax[ccm]);
		dbetayx = cx_1*(betay[pcc]-betay[mcc]);
		dbetayy = cy_1*(betay[cpc]-betay[cmc]);
		dbetayz = cz_1*(betay[ccp]-betay[ccm]);
		dbetazx = cx_1*(betaz[pcc]-betaz[mcc]);
		dbetazy = cy_1*(betaz[cpc]-betaz[cmc]);
		dbetazz = cz_1*(betaz[ccp]-betaz[ccm]);

		dpsix  = cx_1*(psiq[pcc]-psiq[mcc]);
		dpsiy  = cy_1*(psiq[cpc]-psiq[cmc]);
		dpsiz  = cz_1*(psiq[ccp]-psiq[ccm]);
		dapsix = cx_1*(apsiq[pcc]-apsiq[mcc]);
		dapsiy = cy_1*(apsiq[cpc]-apsiq[cmc]);
		dapsiz = cz_1*(apsiq[ccp]-apsiq[ccm]);

		Axx = 0.5*pow(1+psiq[ccc],6)/alpha*(4./3.*dbetaxx-2./3.*(dbetayy+dbetazz));
		Ayy = 0.5*pow(1+psiq[ccc],6)/alpha*(4./3.*dbetayy-2./3.*(dbetaxx+dbetazz));
		Azz = 0.5*pow(1+psiq[ccc],6)/alpha*(4./3.*dbetazz-2./3.*(dbetayy+dbetaxx));
		Axy = 0.5*pow(1+psiq[ccc],6)/alpha*(dbetayx+dbetaxy);
		Axz = 0.5*pow(1+psiq[ccc],6)/alpha*(dbetazx+dbetaxz);
		Ayz = 0.5*pow(1+psiq[ccc],6)/alpha*(dbetazy+dbetayz);
		
		AA[ccc]  = (Axx*Axx+Ayy*Ayy+Azz*Azz+2.*(Axy*Axy+Axz*Axz+Ayz*Ayz));

		if(Getv("BBID_solve_spin","corotational")){
			Sx[ccc] = 2.*pow(1+psiq[ccc],-7)*(Axx*(dapsix-7.*alpha*dpsix)+Axy*(dapsiy-7.*alpha*dpsiy)+Axz*(dapsiz-7*alpha*dpsiz)) +
				16.*PI*(epsilon+p)*alpha*alpha*ut[ccc]*ut[ccc]*pow(1+psiq[ccc],4)*(betax[ccc]-Omega*yp[ijk]);
			Sy[ccc] = 2.*pow(1+psiq[ccc],-7)*(Axy*(dapsix-7.*alpha*dpsix)+Ayy*(dapsiy-7.*alpha*dpsiy)+Ayz*(dapsiz-7*alpha*dpsiz)) +
				16.*PI*(epsilon+p)*alpha*alpha*ut[ccc]*ut[ccc]*pow(1+psiq[ccc],4)*(betay[ccc]+Omega*xp[ijk]);
		}
		else if(Getv("BBID_solve_spin","irrotational")){
			Sx[ccc] = 2.*pow(1+psiq[ccc],-7)*(Axx*(dapsix-7.*alpha*dpsix)+Axy*(dapsiy-7.*alpha*dpsiy)+Axz*(dapsiz-7*alpha*dpsiz)) +
				16.*PI*(epsilon+p)*alpha*alpha*ut[ccc]*ut[ccc]*pow(1+psiq[ccc],4)*(betax[ccc]);
			Sy[ccc] = 2.*pow(1+psiq[ccc],-7)*(Axy*(dapsix-7.*alpha*dpsix)+Ayy*(dapsiy-7.*alpha*dpsiy)+Ayz*(dapsiz-7*alpha*dpsiz)) +
				16.*PI*(epsilon+p)*alpha*alpha*ut[ccc]*ut[ccc]*pow(1+psiq[ccc],4)*(betay[ccc]+Omega*sqrt(1-e)*(xp[ijk]>0?level->grid->puncpos[0][0]:level->grid->puncpos[1][0] ));
		}
		Sz[ccc] = 2.*pow(1+psiq[ccc],-7)*(Axz*(dapsix-7.*alpha*dpsix)+Ayz*(dapsiy-7.*alpha*dpsiy)+Azz*(dapsiz-7*alpha*dpsiz)) +
			16.*PI*(epsilon+p)*alpha*alpha*ut[ccc]*ut[ccc]*pow(1+psiq[ccc],4)*(betaz[ccc]);

		S[ccc]     = ((epsilon+p)*(alpha*alpha*ut[ccc]*ut[ccc]-1)+ 3.*p);
		rhoH[ccc]  = (epsilon + (epsilon+p)*(pow(alpha*ut[ccc],2)-1));

		if(debug) if (!finite(AA[ccc]))  {printf("HERE %f,%f,%f,%f",psiq[ccc],ut[ccc],apsiq[ccc],alpha);errorexit("something bad happend with AA");}
	} endfor;

	
	bampi_synchronize(level, Ind("BBID_AA"));
	bampi_synchronize(level, Ind("BBID_Six"));
	bampi_synchronize(level, Ind("BBID_S"));
	bampi_synchronize(level, Ind("BBID_rhoH"));	

}


void calculate_densities(tL *level){
	
	double *alpha, *apsiq, *psiq, *omegaq;
	double Omega, del_Omega=0;
	
	alpha = Ptr(level, "alpha");
	apsiq = Ptr(level, "BBID_u_apsi");
	psiq  = Ptr(level, "BBID_psi");
	omegaq  = Ptr(level, "BBID_u_omega");
	
	double etaC0	= exp(Getd("BBID_tov_hc1"))-1; 	//Defined density
	double etaC1	= exp(Getd("BBID_tov_hc2"))-1; 	//Defined density
	
	forallpoints_ijk(level) { //making sure that the correct values are there
	
		psiq[ijk] = pow(1+omegaq[ijk],-1/n)-1;
		alpha[ijk]= (apsiq[ijk]+1)/(psiq[ijk]+1);
	} endfor_ijk;
	
	int tmp;
	for(tmp = 0 ; tmp <30; tmp++){
		bampi_barrier();
		Omega = Getd("BBID_solve_Omega");
		
		/* Only use finest level to compute constants */
		if(level->l == level->grid->lmax){ 	
			if(tmp==0){ 		//prevent a lot of unecessary output
				printf("║\n║ Finding constants Omega, xcom and Econst0/1 with Newton-Raphson on level %d\n",level->l);
				printf("╠═══ Using %s configuration\n",Gets("BBID_solve_spin"));
				
			}
			find_constants(level,Omega, Getd("BBID_Econst_0"),Getd("BBID_Econst_1"));
			del_Omega = Omega-Getd("BBID_solve_Omega");
			if(tmp==0) printf("╠═══ First change in Omega: %f%%\n",del_Omega/Omega*100);
		}
		
		/* Now use constants to compute new eta. This also changes u^t, which will change the constants again. Therefore iterate afterwards"*/
		if(tmp==0) printf("║\n║ Updating density distribution on level %d\n",level->l);
		compute_eta(level);
		
		if(level->l == level->grid->lmax){ 	
			if(tmp==0) printf("║\n║ Now iterate until constants don't change: ");
			printf("■ ");
		}
		if(fabs(del_Omega/Omega)<1e-12)  break;
	}
	if(tmp>=29)
		printf("WARNING: Constants kept changing during iteration. Something is going to explode!\n");
	if(level->l ==level->grid->lmax) 
		printf("\n╠═══ Found Constants: Omega=%f, xcom=%f, Econst0=%f, Econst1=%f\n",Getd("BBID_solve_Omega"), Getd("BBID_solve_com"), Getd("BBID_Econst_0"),Getd("BBID_Econst_1"));
	printf("╠═══ Ratio of central eta to specified etaC: Star 1... %f, Star 2... %f\n",interpolate_xyz_scalar(level, level->grid->puncpos[0][0], 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)/etaC0,interpolate_xyz_scalar(level, level->grid->puncpos[0][0], 0, 0, Ind("BBID_eta"), Geti("order_RP"),LAGRANGE)/etaC1);
	printf("║ Compute AA and new source terms rho_H, S^i, S on level %d\n",level->l);
	compute_sources(level);
}



/* non-linear Gauss-Seidel:  v = u + (f - Lu)/Lii */
/* 
*/


void LBBID_GS(tL *level, tVarList *vlv, tVarList *vlu)
{	

	//if (PR) printf("solve   linear system\n");
	
	double cx = 1/(level->dx*level->dx);			//for second derivatives
	double cy = 1/(level->dy*level->dy);
	double cz = 1/(level->dz*level->dz);
	double cc = -2.*(cx + cy + cz);

	double cx_1 = 0.5/(level->dx);			    //for first derivatives
	double cy_1=  0.5/(level->dy);
	double cz_1=  0.5/(level->dz);

	double *f_psi	= Ptr(level, "BBID_f_omega");
	double *f_apsi	= Ptr(level, "BBID_f_apsi");
	double *f_betax	= Ptr(level, "BBID_f_betax");
	double *f_betay	= Ptr(level, "BBID_f_betay");
	double *f_betaz	= Ptr(level, "BBID_f_betaz");
	double *v_omegaq, *u_omegaq;					//psiq=psi-1
	double *v_apsiq, *u_apsiq;
	double *v_betax, *u_betax;
	double *v_betay, *u_betay;
	double *v_betaz, *u_betaz;
	double lu, lii;
	double *xp, *yp, *zp;
	
 	
	//double *ut 	= Ptr(level, "BBID_ut");
	double *eta 	= Ptr(level, "BBID_eta");
	double *S	= Ptr(level, "BBID_S");
	double *Sx	= Ptr(level, "BBID_Six");
	double *Sy	= Ptr(level, "BBID_Siy");
	double *Sz	= Ptr(level, "BBID_Siz");
	double *rhoH	= Ptr(level, "BBID_rhoH");
	double *AA	= Ptr(level, "BBID_AA");
	
// 	compute_sources(level);
	
	u_omegaq = VLPtr(vlu, 0);
	v_omegaq = VLPtr(vlv, 0);

	u_apsiq = VLPtr(vlu, 1);
	v_apsiq = VLPtr(vlv, 1);

	if (solve_beta == 1){
	u_betax = VLPtr(vlu, 2);
	v_betax = VLPtr(vlv, 2);

	u_betay = VLPtr(vlu, 3);
	v_betay = VLPtr(vlv, 3);

	u_betaz = VLPtr(vlu, 4);
	v_betaz = VLPtr(vlv, 4);
	}



	//printf("Ea=%f...Eb=%f - level: %d \n",E_const0a, E_const0b, level->l);
	/* interior */
	forinner19(level) {

		// Laplace Psi
		if (solve_omega==1){
			lu = cc * u_omegaq[ccc] +
				cx * (u_omegaq[mcc] + u_omegaq[pcc]) +
				cy * (u_omegaq[cmc] + u_omegaq[cpc]) +
				cz * (u_omegaq[ccm] + u_omegaq[ccp]) 
				- (n+1)/n*(pow(cx_1*(u_omegaq[pcc]-u_omegaq[mcc]),2)+pow(cy_1*(u_omegaq[cpc]-u_omegaq[cmc]),2)+pow(cz_1*(u_omegaq[ccp]-u_omegaq[ccm]),2))*pow(1+u_omegaq[ccc],-1)
				- 2*PI*n*rhoH[ccc]*pow(1+u_omegaq[ccc],-4./n+1) 
				- n/8. * pow(1+u_omegaq[ccc],8./n+1)*AA[ccc];
				
			lii = cc 	+ (n+1)/n*(pow(cx_1*(u_omegaq[pcc]-u_omegaq[mcc]),2)+pow(cy_1*(u_omegaq[cpc]-u_omegaq[cmc]),2)+pow(cz_1*(u_omegaq[ccp]-u_omegaq[ccm]),2))*pow(1+u_omegaq[ccc],-2)
				- 2*PI*(-4+n)*rhoH[ccc]*pow(1+u_omegaq[ccc],-4./n) 
				- (8+n)/8. * pow(1+u_omegaq[ccc],8./n)*AA[ccc];
		
			v_omegaq[ccc] = u_omegaq[ccc] + (f_psi[ccc] - lu)/lii;
		} else {
			v_omegaq[ccc] = u_omegaq[ccc];
		}

		// Laplace apsi
		if (solve_apsi==1){
			lu = cc * u_apsiq[ccc] +
				cx * (u_apsiq[mcc] + u_apsiq[pcc]) +
				cy * (u_apsiq[cmc] + u_apsiq[cpc]) +
				cz * (u_apsiq[ccm] + u_apsiq[ccp]) -
				(1+u_apsiq[ccc])*(0.875*pow(1+u_omegaq[ccc],8./n)*AA[ccc]+2*PI*pow(1+u_omegaq[ccc],-4./n)*(rhoH[ccc]+2*S[ccc]));
			lii = cc - (0.875*pow(1+u_omegaq[ccc],8./n)*AA[ccc]+2*PI*pow(1+u_omegaq[ccc],-4./n)*(rhoH[ccc]+2*S[ccc]));


			v_apsiq[ccc] = u_apsiq[ccc] + (f_apsi[ccc] - lu)/lii;
		}  else {
			v_apsiq[ccc] = u_apsiq[ccc];
		}
			
		if (solve_beta == 1){
			// Laplace Beta^x
					
			lu = cc * u_betax[ccc] +
					cx * (u_betax[mcc] + u_betax[pcc]) +
					cy * (u_betax[cmc] + u_betax[cpc]) +
					cz * (u_betax[ccm] + u_betax[ccp]) +
					1./3. * ( cx*(u_betax[pcc]-2.*u_betax[ccc]+u_betax[mcc])
						+  cx_1*cy_1*(u_betay[ppc]-u_betay[pmc]-u_betay[mpc]+u_betay[mmc])
						+  cx_1*cz_1*(u_betaz[pcp]-u_betaz[pcm]-u_betaz[mcp]+u_betaz[mcm]))
					-Sx[ccc];
			lii = cc -2./3.*cx;

			v_betax[ccc] =  u_betax[ccc] + (f_betax[ccc] - lu)/lii;

			// Laplace Beta^y

			lu =  cc * u_betay[ccc] +
					cx * (u_betay[mcc] + u_betay[pcc]) +
					cy * (u_betay[cmc] + u_betay[cpc]) +
					cz * (u_betay[ccm] + u_betay[ccp]) +
					1./3. * ( cy*(u_betay[cpc]-2.*u_betay[ccc]+u_betay[cmc])
						+  cx_1*cy_1*(u_betax[ppc]-u_betax[pmc]-u_betax[mpc]+u_betax[mmc])
						+  cy_1*cz_1*(u_betaz[cpp]-u_betaz[cpm]-u_betaz[cmp]+u_betaz[cmm]))
					-Sy[ccc];
			lii = cc -2./3.*cx;

			v_betay[ccc] = u_betay[ccc] + (f_betay[ccc] - lu)/lii;


			// Laplace Beta^z

			lu =  cc * u_betaz[ccc] +
					cx * (u_betaz[mcc] + u_betaz[pcc]) +
					cy * (u_betaz[cmc] + u_betaz[cpc]) +
					cz * (u_betaz[ccm] + u_betaz[ccp]) +
					1./3. * ( cz*(u_betaz[ccp]-2.*u_betaz[ccc]+u_betaz[ccm])
						+  cx_1*cz_1*(u_betax[pcp]-u_betax[pcm]-u_betax[mcp]+u_betax[mcm])
						+  cy_1*cz_1*(u_betay[cpp]-u_betay[cpm]-u_betay[cmp]+u_betay[cmm]))
					-Sz[ccc];
			lii = cc -2./3.*cx;

			v_betaz[ccc] = u_betaz[ccc] + (f_betaz[ccc] - lu)/lii;
		}  else {
			v_betax[ccc] = u_betax[ccc];
			v_betay[ccc] = u_betay[ccc];
			v_betaz[ccc] = u_betaz[ccc];
		}


	} endfor;



	/* which is actually needed ? */
	bampi_vlsynchronize(vlv);
	//printf("\nu_apsi[211]=\%f %f %f %f\n",u_apsi[211], v_apsi[211], u_psi[211], v_psi[211]);


	set_boundary_elliptic(level, vlv);

}





/* apply non-linear elliptic operator:  lu = L(u) */
/*
*/
void LBBID_nonlin(tL *level, tVarList *vllu, tVarList *vlu)
{  
	

	double cx = 1/(level->dx*level->dx);			//for second derivatives
	double cy = 1/(level->dy*level->dy);
	double cz = 1/(level->dz*level->dz);
	double cc = -2.*(cx + cy + cz);

	double cx_1 = 0.5/(level->dx);			    //for first derivatives
	double cy_1=  0.5/(level->dy);
	double cz_1=  0.5/(level->dz);


	double *lu_omegaq, *lu_apsi;
	double *lu_betax,*lu_betay,*lu_betaz;

	double *u_omegaq, *u_apsiq,*u_betax, *u_betay, *u_betaz;
	double Axx,Ayy, Azz, Axy, Axz, Ayz;

	tVarList *vlut = VLPtrEnable1(level, "BBID_ut");
	tVarList *vleta = VLPtrEnable1(level, "BBID_eta");
	tVarList *vlrhoH = VLPtrEnable1(level, "BBID_rhoH");
	tVarList *vlS = VLPtrEnable1(level, "BBID_S");
	tVarList *vlSi  = VLPtrEnable1(level, "BBID_Six");
	tVarList *vlAA = VLPtrEnable1(level, "BBID_AA");
	
	double *ut = Ptr(level, "BBID_ut");
	double *eta = Ptr(level, "BBID_eta");
	double *S	= Ptr(level, "BBID_S");
	double *Sx	= Ptr(level, "BBID_Six");
	double *Sy	= Ptr(level, "BBID_Siy");
	double *Sz	= Ptr(level, "BBID_Siz");
	double *rhoH	= Ptr(level, "BBID_rhoH");
	double *AA	= Ptr(level, "BBID_AA");
	
// 	compute_sources(level);
	
	u_omegaq   = VLPtr(vlu, 0);
	u_apsiq    = VLPtr(vlu, 1);
	if (solve_beta == 1){
	u_betax    = VLPtr(vlu, 2);
	u_betay    = VLPtr(vlu, 3);
	u_betaz    = VLPtr(vlu, 4);
	}

	lu_omegaq  = VLPtr(vllu, 0);
	lu_apsi    = VLPtr(vllu, 1);
	if (solve_beta == 1){lu_betax   = VLPtr(vllu, 2);
	lu_betay   = VLPtr(vllu, 3);
	lu_betaz   = VLPtr(vllu, 4);}

	

	if (1) bampi_vlsynchronize(vlu);

	/* apply boundary conditions */
	set_boundary_elliptic(level, vlu);



	
	//printf("Nonlin: Ea=%f...Eb=%f - level: %d \n",E_const0a, E_const0b, level->l);
	/* interior */
	forinner19(level) {


		//  Laplace Psi
		if (solve_omega==1){
			lu_omegaq[ccc] = cc * u_omegaq[ccc] +
				cx * (u_omegaq[mcc] + u_omegaq[pcc]) +
				cy * (u_omegaq[cmc] + u_omegaq[cpc]) +
				cz * (u_omegaq[ccm] + u_omegaq[ccp]) +
					- (n+1)/n*(pow(cx_1*(u_omegaq[pcc]-u_omegaq[mcc]),2)+pow(cy_1*(u_omegaq[cpc]-u_omegaq[cmc]),2)+pow(cz_1*(u_omegaq[ccp]-u_omegaq[ccm]),2))*pow(1+u_omegaq[ccc],-1)
					- 2*PI*n*rhoH[ccc]*pow(1+u_omegaq[ccc],-4./n+1) 
					- n/8. * pow(1+u_omegaq[ccc],8./n+1)*AA[ccc];
		} else {
			lu_omegaq[ccc] = 0.;
		}

		// Laplace apsi
		if (solve_apsi==1){
			lu_apsi[ccc] = cc * u_apsiq[ccc] +
				cx * (u_apsiq[mcc] + u_apsiq[pcc]) +
				cy * (u_apsiq[cmc] + u_apsiq[cpc]) +
				cz * (u_apsiq[ccm] + u_apsiq[ccp]) -
				(1+u_apsiq[ccc])*(0.875*pow(1+u_omegaq[ccc],8./n)*AA[ccc]+2*PI*pow(1+u_omegaq[ccc],-4./n)*(rhoH[ccc]+2*S[ccc]));
		} else {
			lu_apsi[ccc] = 0;
		}

		if(solve_beta == 1){
			// Laplace Beta^x
			lu_betax[ccc] = cc * u_betax[ccc] +
					cx * (u_betax[mcc] + u_betax[pcc]) +
					cy * (u_betax[cmc] + u_betax[cpc]) +
					cz * (u_betax[ccm] + u_betax[ccp]) +
					1./3. * ( cx*(u_betax[pcc]-2.*u_betax[ccc]+u_betax[mcc])
						+  cx_1*cy_1*(u_betay[ppc]-u_betay[pmc]-u_betay[mpc]+u_betay[mmc])
						+  cx_1*cz_1*(u_betaz[pcp]-u_betaz[pcm]-u_betaz[mcp]+u_betaz[mcm]))
					-Sx[ccc];

			// Laplace Beta^y
			lu_betay[ccc] = cc * u_betay[ccc] +
					cx * (u_betay[mcc] + u_betay[pcc]) +
					cy * (u_betay[cmc] + u_betay[cpc]) +
					cz * (u_betay[ccm] + u_betay[ccp]) +
					1./3. * ( cy*(u_betay[cpc]-2.*u_betay[ccc]+u_betay[cmc])
						+  cx_1*cy_1*(u_betax[ppc]-u_betax[pmc]-u_betax[mpc]+u_betax[mmc])
						+  cy_1*cz_1*(u_betaz[cpp]-u_betaz[cpm]-u_betaz[cmp]+u_betaz[cmm]))
					-Sy[ccc];

			// Laplace Beta^z
			lu_betaz[ccc] = cc * u_betaz[ccc] +
					cx * (u_betaz[mcc] + u_betaz[pcc]) +
					cy * (u_betaz[cmc] + u_betaz[cpc]) +
					cz * (u_betaz[ccm] + u_betaz[ccp]) +
					1./3. * ( cz*(u_betaz[ccp]-2.*u_betaz[ccc]+u_betaz[ccm])
						+  cx_1*cz_1*(u_betax[pcp]-u_betax[pcm]-u_betax[mcp]+u_betax[mcm])
						+  cy_1*cz_1*(u_betay[cpp]-u_betay[cpm]-u_betay[cmp]+u_betay[cmm]))
					-Sz[ccc];
		}/*else{
			lu_betax[ccc] = 0;
			lu_betay[ccc] = 0;
			lu_betaz[ccc] = 0;
		}*/
			

		if(debug) if (!finite(lu_omegaq[ccc]))  {printf("HERE %f,%f,%f,%f",eta[ccc],ut[ccc],AA[ccc],rhoH[ccc]);errorexit("something bad happend with lu_psi");}
		if(debug) if (!finite(lu_apsi[ccc])) errorexit(" something bad happend with lu_apsi");


	} endfor;
	//errorexit("asd");
	bampi_vlsynchronize(vllu);

}


int BBID_solve(tL* level)
{


	if (level->l!=0)
	return 0;

	if (level->shells)
	errorexit("hahaha, nice joke, nope, this is not implemented for shells");

	tG *g = level->grid;

	tVarList *vlu =vlalloc(level);
	tVarList *vlv =vlalloc(level);
	tVarList *vlf =vlalloc(level);
	tVarList *vlr =vlalloc(level);
	tVarList *vlc = vlalloc(level);
	
	
	
	int itmax = Geti("BBID_solve_itmax");
	double tol = Getd("BBID_solve_tolerance");
	double normres, normresnonlin, u;
	int i, inewton, l;

	/* initialize+set computation variables */
	for (l = g->lmax; l >= level->l; l--) {
		enablevar(g->level[l], Ind("BBID_ut"));
		enablevar(g->level[l], Ind("BBID_eta"));
		enablevar(g->level[l], Ind("BBID_energyconst"));
		enablevar(g->level[l], Ind("BBID_S"));
		enablevar(g->level[l], Ind("BBID_Six"));	
		enablevar(g->level[l], Ind("BBID_rhoH"));
		enablevar(g->level[l], Ind("BBID_AA"));
		enablevar(g->level[l], Ind("BBID_psi"));
		enablevar(g->level[l], Ind("BBID_omega_old"));
		enablevar(g->level[l], Ind("BBID_apsi_old"));
		enablevar(g->level[l], Ind("BBID_beta_oldx"));
	}

	/* varname, farlimit, falloff, prolongation speed */
	VarNameSetBoundaryInfo("BBID_u_omega", 0, 1, .0);
	VarNameSetBoundaryInfo("BBID_u_apsi",	0, 1, .0);
	VarNameSetBoundaryInfo("BBID_u_betax",	0, 2, .0);
	VarNameSetBoundaryInfo("BBID_u_betay",	0, 2, .0);
	VarNameSetBoundaryInfo("BBID_u_betaz",	0, 4, .0);
	Appends("boundary", "robin");

	/* only coarsest box is allowed to have physical boundary for now */
	find_robin_normal(g->level[g->lmin]);
	//find_robin_normal(level);

	/* make variable list for coefficients */ 
	vlpush( vlc, Ind("BBID_ut"));
	vlpush( vlc, Ind("BBID_eta"));
	vlpush( vlc, Ind("BBID_S"));
	vlpush( vlc, Ind("BBID_Six"));
	vlpush( vlc, Ind("BBID_rhoH"));
	vlpush( vlc, Ind("BBID_AA"));
	vlpush( vlc, Ind("BBID_psi"));
	//vlpush( vlc, Ind("alpha"));
	//vlpush( vlc, Ind("BBID_energyconst"));


	vlpush( vlu, Ind("BBID_u_omega"));
	vlpush( vlu, Ind("BBID_u_apsi"));
	vlpush( vlu, Ind("BBID_u_betax"));

	vlpush( vlf, Ind("BBID_f_omega"));
	vlpush( vlf, Ind("BBID_f_apsi"));
	vlpush( vlf, Ind("BBID_f_betax"));


	vlpush( vlr, Ind("BBID_r_omega"));
	vlpush( vlr, Ind("BBID_r_apsi"));
	vlpush( vlr, Ind("BBID_r_betax"));

	vlpush( vlv, Ind("BBID_v_omega"));
	vlpush( vlv, Ind("BBID_v_apsi"));
	vlpush( vlv, Ind("BBID_v_betax"));
	
	

	for (l = g->lmax; l >= level->l; l--) {
		vlenablelevel(g->level[l],vlu);
		vlenablelevel(g->level[l],vlf);
		vlenablelevel(g->level[l],vlv);
		vlenablelevel(g->level[l],vlr);
		vlenablelevel(g->level[l],vlc);
		
	}


	/* Initializing Values */
	double psi4, *psiq, *apsiq, *betaxq,*betayq,*betazq, *omegaq, *omega_old, *apsi_old, *beta_oldx, *beta_oldy, *beta_oldz;
	double *alpha, *betax,*betay,*betaz;
	double *gxx,*gxy,*gxz,*gyy,*gyz,*gzz;
	double *Kxx,*Kxy,*Kxz,*Kyy,*Kyz,*Kzz;
	double *xp,*rho,*epsl,*vx,*vy,*vz, *ut, *eta, *ADMrho;
	double *yp,*zp, *E_const;
	double vy1   =  - Getd("my1");
	double vy2   =  - Getd("my2");

	double eta0     = exp(Getd("BBID_tov_hc"))-1;
	double polytrope= 1./EOS.GAMMAMO;
	double K        = EOS.K; 
	


	
	for (l = g->lmin; l <= g->lmax; l++) {
		// BAM VARS
		alpha = Ptr(g->level[l], "alpha");
		betax = Ptr(g->level[l], "betax");
		betay = Ptr(g->level[l], "betay");
		betaz = Ptr(g->level[l], "betaz");
		gxx   = Ptr(g->level[l], "adm_gxx");
		gxy   = Ptr(g->level[l], "adm_gxy");
		gxz   = Ptr(g->level[l], "adm_gxz");
		gyy   = Ptr(g->level[l], "adm_gyy");
		gyz   = Ptr(g->level[l], "adm_gyz");
		gzz   = Ptr(g->level[l], "adm_gzz");
		Kxx   = Ptr(g->level[l], "adm_Kxx");
		Kxy   = Ptr(g->level[l], "adm_Kxy");
		Kxz   = Ptr(g->level[l], "adm_Kxz");
		Kyy   = Ptr(g->level[l], "adm_Kyy");
		Kyz   = Ptr(g->level[l], "adm_Kyz");
		Kzz   = Ptr(g->level[l], "adm_Kzz");
		rho   = Ptr(g->level[l], "grhd_rho");
		epsl  = Ptr(g->level[l], "grhd_epsl");
		vx    = Ptr(g->level[l], "grhd_vx");
		vy    = Ptr(g->level[l], "grhd_vy");
		vz    = Ptr(g->level[l], "grhd_vz");
		xp     = Ptr(g->level[l], "x");
		// ELLIPTIC VARS
		omegaq  = Ptr(g->level[l], "BBID_u_omega");
		apsiq	= Ptr(g->level[l], "BBID_u_apsi");
		betaxq	= Ptr(g->level[l], "BBID_u_betax");
		betayq	= Ptr(g->level[l], "BBID_u_betay");
		betazq	= Ptr(g->level[l], "BBID_u_betaz");
		ut     	= Ptr(g->level[l], "BBID_ut");
		eta     = Ptr(g->level[l], "BBID_eta");	
		psiq    = Ptr(g->level[l], "BBID_psi");
		
		forallpoints_ijk(g->level[l]) {
			//assigning values 
			psi4 = gxx[ijk];
			omegaq[ijk]	= pow(psi4,-0.25*n)-1;
			apsiq[ijk]	= alpha[ijk]*pow(omegaq[ijk]+1,-1./n)-1;
			betaxq[ijk]	= betax[ijk];
			betayq[ijk]	= -betay[ijk];
			betazq[ijk]	= betaz[ijk];
			eta[ijk] 	= K*(1+polytrope)*pow(rho[ijk],1/polytrope);
			ut[ijk]		= 1./alpha[ijk];
// 			ut[ijk]		= 1./sqrt(alpha[ijk]*alpha[ijk]-pow(1+psiq[ijk],4)*(betax[ijk]*betax[ijk]+betay[ijk]*betay[ijk]+betaz[ijk]*betaz[ijk]+2*betay[ijk]*VY+VY*VY));
		} endfor_ijk;
	}

	double p1old[3] =	{g->level[g->lmax]->grid->puncpos[0][0],g->level[g->lmax]->grid->puncpos[0][1],g->level[g->lmax]->grid->puncpos[0][2]};
	double p2old[3] =	{g->level[g->lmax]->grid->puncpos[1][0],g->level[g->lmax]->grid->puncpos[1][1],g->level[g->lmax]->grid->puncpos[1][2]};
	double p1new[3], p2new[3];
	
	double e = Getd("BBID_solve_ecc");
	double Omega = 0;				//Bad Initial Guess for Omega and xcom
	if(level->grid->puncpos[0][0] != 0 ) Omega=(vy1>0?vy1:0.1)/(level->grid->puncpos[0][0]);
	Setd("BBID_solve_Omega",Omega);
	
	double xcom =(massintegral(g->level[g->lmax],0)*p1old[0]+massintegral(g->level[g->lmax],1)*p2old[0])/(massintegral(g->level[g->lmax],0)+massintegral(g->level[g->lmax],1));
	printf("Initial guess for Center of mass: %f\n",xcom);
	Setd("BBID_solve_com", xcom);

	/* iterate over all */
	int iter;
	double E_const0a, E_const0b, beta2, ham_old, *ham_new;
	tVarList *vl_constr = vlalloc(level);
	vlpush(vl_constr, Ind("BBID_omega_old"));
	ham_old=1;

	int tmp;
	double rb = Getd("BBID_rb")+p1old[0];
	Setd("BBID_rb",rb);
	
	printf("Initial resaling parameter  r_B: %.16e\n", rb);
	double rho0=pow(eta0/(K*(1+polytrope)),polytrope);


	
	/* Load data */
	
	if(Getv("BBID_solve_load","yes")){
		printf("\nRead binary neutron star provided by BBID_solve\n Iterations are set to 0. \n\n");
		for (l = g->lmin; l <= g->lmax; l++){
			load_data(g->level[l]);
		}
		itmax = 0;
		Setd("BBID_solve_soft_ell",1.0);
		Setd("BBID_solve_soft_eta",1.0);
	}
	
	for (iter=0; iter <=itmax; iter++) {
		

		/* compute Source terms */
		
		
		for (l = g->lmax; l >= g->lmin; l--){
			
			calculate_densities(g->level[l]);
		}
		
		
		
		
		track_moving_puncture_extremum(level, 0, p1old, p1new);
		track_moving_puncture_extremum(level, 1, p2old, p2new);
		
		printf("║\n║\n║ Position of star 1 (on lowest level): %f,%f,%f\n", p1new[0],p1new[1],p1new[2]);
		printf("║ Position of star 2 (on lowest level): %f,%f,%f\n", p2new[0],p2new[1],p2new[2]);
		
			
		for (l = g->lmax; l >= g->lmin; l--){
			
			omega_old  	= Ptr(g->level[l], "BBID_omega_old");
			omegaq  	= Ptr(g->level[l], "BBID_u_omega");
			apsi_old  	= Ptr(g->level[l], "BBID_apsi_old");
			apsiq  		= Ptr(g->level[l], "BBID_u_apsi");
			beta_oldx 	= Ptr(g->level[l], "BBID_beta_oldx");
			beta_oldy  	= Ptr(g->level[l], "BBID_beta_oldz");
			beta_oldz  	= Ptr(g->level[l], "BBID_beta_oldy");
			betaxq  	= Ptr(g->level[l], "BBID_u_betax"); 
			betayq  	= Ptr(g->level[l], "BBID_u_betay"); 
			betazq  	= Ptr(g->level[l], "BBID_u_betaz"); 
			forallpoints_ijk(g->level[l]) { //assigning values to ut	
				omega_old[ijk] = omegaq[ijk];	
				apsi_old[ijk]  = apsiq[ijk];
				beta_oldx[ijk] = betaxq[ijk];
				beta_oldy[ijk] = betayq[ijk];
				beta_oldz[ijk] = betazq[ijk];
				
			} endfor_ijk;
	
		}
		
		if(iter< itmax) { //Iteration over elliptic equations, last just enforces \varepsilon = const
			printf("║\n╠════════════════════════════════╦═══════════════╗");
			printf( "\n║Solving for elliptic equations: ║ Iteration %0d   ║", iter);
			printf( "\n╚════════════════════════════════╩═══════════════╝\n\n");

			solve_omega = 1;	//Better don't use
			solve_apsi  = 1;
			solve_beta  = 1;
			

// 			/* fill in overlap (disable for debugging) */
			restrict_prolong_grid(g, vlu);
						
			multigrid(level, vlu, vlf, vlv, vlc,
			itmax, tol, &normres,
			LBBID_nonlin, LBBID_GS);

			/* set boundary for final result */
			for (l = g->lmin; l <= g->lmax; l++)
				set_boundary_elliptic(g->level[l], vlu);
			
			/* fill in overlap (disable for debugging) */
			restrict_prolong_grid(g, vlu);
			
			
			
			
		}

		
		double soft=Getd("BBID_solve_soft_ell");
		
		for (l = g->lmax; l >= g->lmin; l--){ //AddDuplicateEnable may be more appropriate 
			
			omega_old 	= Ptr(g->level[l], "BBID_omega_old");
			omegaq  	= Ptr(g->level[l], "BBID_u_omega");
			apsi_old  	= Ptr(g->level[l], "BBID_apsi_old");
			apsiq  		= Ptr(g->level[l], "BBID_u_apsi");
			alpha  		= Ptr(g->level[l], "alpha");
			beta_oldx 	= Ptr(g->level[l], "BBID_beta_oldx");
			beta_oldy  	= Ptr(g->level[l], "BBID_beta_oldz");
			beta_oldz  	= Ptr(g->level[l], "BBID_beta_oldy");
			betaxq  	= Ptr(g->level[l], "BBID_u_betax"); 
			betayq  	= Ptr(g->level[l], "BBID_u_betay"); 
			betazq  	= Ptr(g->level[l], "BBID_u_betaz"); 
			xp     = Ptr(g->level[l], "x");
			yp     = Ptr(g->level[l], "y");
			zp     = Ptr(g->level[l], "z");
			forallpoints_ijk(g->level[l]) { //assigning values to ut	
				omegaq[ijk] = soft*omegaq[ijk]+(1-soft)*omega_old[ijk];	
				apsiq[ijk]  = soft*apsiq[ijk]+(1-soft)*apsi_old[ijk];
 				betaxq[ijk] = soft*betaxq[ijk]+(1-soft)*beta_oldx[ijk];	
 				betayq[ijk] = soft*betayq[ijk]+(1-soft)*beta_oldy[ijk];
 				betazq[ijk] = soft*betazq[ijk]+(1-soft)*beta_oldz[ijk];
				//alpha[ijk] = (1+apsiq[ijk])/pow(1+omegaq[ijk],-n);
				omega_old[ijk] -= omegaq[ijk];
				
			} endfor_ijk;
			
	
		}
		bampi_allreduce_norm2(vl_constr, &ham_new);
		printf("║ Change in Psi : %.5e\n", ham_new[0]);
		if(ham_new[0]<10e-11) break;


	} //END OF ITER

	/* do some final boundary things */
	for (l = g->lmin; l <= g->lmax; l++) { 
		set_boundary_elliptic(g->level[l], vlu);
		bampi_synchronize(g->level[l],Ind("BBID_ut"));

	}

	
	restrict_prolong_grid(g, vlu);

	


	
	/* Save elliptic variables for each processor and level */
	printf("Saving initial data for future runs\n");
	for (l = g->lmin; l <= g->lmax; l++){
		save_data(g->level[l]);
	}
	

	
	/* copy solved data */

	double dbetaxx,dbetaxy,dbetaxz,dbetayx,dbetayy,dbetayz,dbetazx,dbetazy,dbetazz;
	for (l = g->lmin; l <= g->lmax; l++) {
		
		bampi_synchronize(g->level[l],Ind("BBID_u_betax"));
		set_boundary_extrapolate(g->level[l], Ind("BBID_u_betax"));
		set_boundary_extrapolate(g->level[l], Ind("BBID_u_betay"));
		set_boundary_extrapolate(g->level[l], Ind("BBID_u_betaz"));

		double cx_1 = 0.5/(g->level[l]->dx);			//for first derivatives
		double cy_1=  0.5/(g->level[l]->dy);
		double cz_1=  0.5/(g->level[l]->dz);
		double xcom = Getd("BBID_solve_com");			//center of mass
		Omega = Getd("BBID_solve_Omega");
		
		// BAM VARS
		alpha = Ptr(g->level[l], "alpha");
		betax = Ptr(g->level[l], "betax");
		betay = Ptr(g->level[l], "betay");
		betaz = Ptr(g->level[l], "betaz");
		gxx   = Ptr(g->level[l], "adm_gxx");
		gxy   = Ptr(g->level[l], "adm_gxy");
		gxz   = Ptr(g->level[l], "adm_gxz");
		gyy   = Ptr(g->level[l], "adm_gyy");
		gyz   = Ptr(g->level[l], "adm_gyz");
		gzz   = Ptr(g->level[l], "adm_gzz");
		Kxx   = Ptr(g->level[l], "adm_Kxx");
		Kxy   = Ptr(g->level[l], "adm_Kxy");
		Kxz   = Ptr(g->level[l], "adm_Kxz");
		Kyy   = Ptr(g->level[l], "adm_Kyy");
		Kyz   = Ptr(g->level[l], "adm_Kyz");
		Kzz   = Ptr(g->level[l], "adm_Kzz");
		rho   = Ptr(g->level[l], "grhd_rho");
		epsl  = Ptr(g->level[l], "grhd_epsl");
		vx    = Ptr(g->level[l], "grhd_vx");
		vy    = Ptr(g->level[l], "grhd_vy");
		vz    = Ptr(g->level[l], "grhd_vz");
		// ELLIPTIC VARS
		ut      = Ptr(g->level[l], "BBID_ut");	
		eta	= Ptr(g->level[l], "BBID_eta");
		omegaq	= Ptr(g->level[l], "BBID_u_omega");	
		apsiq   = Ptr(g->level[l], "BBID_u_apsi");	
		psiq	= Ptr(g->level[l], "BBID_psi");	
		// COORDINATES
		xp     = Ptr(g->level[l], "x");
		yp     = Ptr(g->level[l], "y");
		zp     = Ptr(g->level[l], "z");


		set_boundary_extrapolate(g->level[l], Ind("BBID_ut"));
		

		forallpoints_ijk(g->level[l]){

			//ut[ijk] = (1+psiq[ijk])/(1+apsiq[ijk]);

			alpha[ijk] = (1+apsiq[ijk])/(pow(1+omegaq[ijk],-1./n));
			
			gxx[ijk] = gyy[ijk] = gzz[ijk] = pow(1.+omegaq[ijk],-4./n);
			gxy[ijk] = gxz[ijk] = gyz[ijk] = 0;

			
			rho[ijk] 	= rho0*pow(eta[ijk]/eta0,polytrope);
			epsl[ijk] 	= eta[ijk]*polytrope/(1+polytrope);		 //total energy density, GRHD_epsl

			betax[ijk] = Ptr(g->level[l], "BBID_u_betax")[ijk]; 
			betay[ijk] = Ptr(g->level[l], "BBID_u_betay")[ijk];
			betaz[ijk] = Ptr(g->level[l], "BBID_u_betaz")[ijk];

		} endfor_ijk;


		bampi_synchronize(g->level[l],Ind("betax"));
		set_boundary_extrapolate(g->level[l], Ind("betax"));

		forinner19(g->level[l]) {
			dbetaxx = cx_1*(betax[pcc]-betax[mcc]);
			dbetaxy = cy_1*(betax[cpc]-betax[cmc]);
			dbetaxz = cz_1*(betax[ccp]-betax[ccm]);
			dbetayx = cx_1*(betay[pcc]-betay[mcc]);
			dbetayy = cy_1*(betay[cpc]-betay[cmc]);
			dbetayz = cz_1*(betay[ccp]-betay[ccm]);
			dbetazx = cx_1*(betaz[pcc]-betaz[mcc]);
			dbetazy = cy_1*(betaz[cpc]-betaz[cmc]);
			dbetazz = cz_1*(betaz[ccp]-betaz[ccm]);
			Kxx[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(4./3.*dbetaxx-2./3.*(dbetayy+dbetazz));
			Kyy[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(4./3.*dbetayy-2./3.*(dbetaxx+dbetazz));
			Kzz[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(4./3.*dbetazz-2./3.*(dbetayy+dbetaxx));
			Kxy[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(dbetayx+dbetaxy);
			Kxz[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(dbetazx+dbetaxz);
			Kyz[ijk]	= 0.5*pow(1+(Ptr(g->level[l], "BBID_psi")[ijk]),4)/alpha[ijk]*(dbetazy+dbetayz);
		}endfor;
		
		set_boundary_extrapolate(g->level[l], Ind("adm_Kxy"));
		set_boundary_extrapolate(g->level[l], Ind("adm_Kxz"));
		set_boundary_extrapolate(g->level[l], Ind("adm_Kyz"));
		set_boundary_extrapolate(g->level[l], Ind("adm_Kxx"));
		set_boundary_extrapolate(g->level[l], Ind("adm_Kyy"));
		set_boundary_extrapolate(g->level[l], Ind("adm_Kzz"));			
		bampi_synchronize(g->level[l],Ind("adm_Kxx"));
		
		if(Getv("BBID_solve_spin","corotational")){			//TODO: This distinction is really annoying. Should use a new variable Killingvector instead.
			forallpoints_ijk(g->level[l]){
				vx[ijk] =  1./alpha[ijk]*(betax[ijk]-Omega*yp[ijk]);
				vy[ijk] =  1./alpha[ijk]*(betay[ijk]+Omega*(xp[ijk]-xcom)); 
				vz[ijk] =  1./alpha[ijk]*(betaz[ijk]);
			} endfor_ijk;
		}
		else if(Getv("BBID_solve_spin","irrotational")){
			forallpoints_ijk(g->level[l]){
				vx[ijk] =  1./alpha[ijk]*(betax[ijk]);
				vy[ijk] =  1./alpha[ijk]*(betay[ijk]+Omega*(xp[ijk]>0?level->grid->puncpos[0][0]:level->grid->puncpos[1][0] )*sqrt(1-e)); 
				vz[ijk] =  1./alpha[ijk]*(betaz[ijk]);
			} endfor_ijk;
		}
		
		
	}
	
	printf("Final Omega: %f\n", Omega);

	/* free unused data */
	for (l = level->l; l <= g->lmax; l++) {
		vlf->level = g->level[l];
		vlr->level = g->level[l];
		vlc->level = g->level[l];
		vlv->level = g->level[l];
		vlu->level = g->level[l];
		vldisable(vlf);
		vldisable(vlr);
		vldisable(vlc);
		vldisable(vlv);
		vldisable(vlu);
		
	}

	vlfree(vlf);
	vlfree(vlr);
	vlfree(vlc);
	vlfree(vlv);
	vlfree(vlu);
	
	return 0;
}



