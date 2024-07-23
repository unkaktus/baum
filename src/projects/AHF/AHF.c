/* AHF.c */
/* Jose Gonzalez 03/07 */
/* mth 06/09 */
/* George Reifenberger 04/2010 */
/* dtim 04/13 */

#include "bam.h"
#include "AHF.h"




int AHF(tL *level) {

    int LMAX1 = Geti("ahf_LMAX");

    static int search0,search1,search2;   // check: consistent with checkpointing?
    static int find_common=0;             // check: consistent with checkpointing?
    // static int ncall=0;   

    int number,l,m;
    int ntheta,nphi;
    int nhorizons;
    double *a0,**ac,**as;
    double pi=acos(-1.0);
    double r_init,r_max;
    double xc=0,yc=0,zc=0;
    double xmn,xmx,ymn,ymx,zmn,zmx;
    double dx,dy,dz,idx,idy,idz;
    tB *fbox;
    //tG *grid=level->grid;
    double **rr,**P,**Y0,**Yc,**Ys; 
    double dtheta,dphi,theta,phi,facl,fac2l;
    int i,j,k,l1,i1,nv,cont;

    double Kxx,Kxy,Kxz,Kyy,Kyz,Kzz;
    double gxx,gxy,gxz,gyy,gyz,gzz;
    double dgxxdx,dgxydx,dgxzdx,dgyydx,dgyzdx,dgzzdx;
    double dgxxdy,dgxydy,dgxzdy,dgyydy,dgyzdy,dgzzdy;
    double dgxxdz,dgxydz,dgxzdz,dgyydz,dgyzdz,dgzzdz;
    double ginvdetxx,ginvdetxy,ginvdetxz,ginvdetyy,ginvdetyz,ginvdetzz;
    double ginvxx,ginvxy,ginvxz,ginvyx,ginvyy,ginvyz,ginvzx,ginvzy,ginvzz,detg;

    double **dPdtheta,**dY0dtheta,**dYcdtheta,**dYsdtheta,**dYcdphi,**dYsdphi;
    double x,y,z,xp,yp,zp,rp,rhop,**u;
    double **dFdx,**dFdy,**dFdz,**sigma;
    double drdx,drdy,drdz,dthetadx,dthetady,dthetadz,dphidx,dphidy,dphidz;
    double drdxdx,drdxdy,drdxdz,drdydy,drdydz,drdzdz;
    double dthetadxdx,dthetadxdy,dthetadxdz,dthetadydy,dthetadydz,dthetadzdz;
    double dphidxdx,dphidxdy,dphidxdz,dphidydy,dphidydz,dphidzdz;
    double **dPdthetadtheta,**dY0dthetadtheta,**dYcdthetadtheta,**dYcdthetadphi,**dYcdphidphi,**dYsdthetadtheta;
    double **dYsdthetadphi,**dYsdphidphi;
    double dFdxdx,dFdxdy,dFdxdz,dFdydx,dFdydy,dFdydz,dFdzdx,dFdzdy,dFdzdz;
    double d2F,**H,**rho,dFdadFdbKab,dFdadFdbFdadb,dFupdx,dFupdy,dFupdz,K;

    double nnFxx,nnFxy,nnFxz,nnFyy,nnFyz,nnFzz;
    double **int1_spec0,*spec0,*spec0_l,**int1_spec_c,*spec_c,*spec_cl,**int1_spec_s,*spec_s,*spec_sl;
    double A,B,alpha=1.0,beta=0.5;
    int flow_iter;

    int half;
    double rtp1,rtm1,rpp1,rpm1,h11,h12,h22,**deth;
    double drdt,drdp,*int1_area,area_l,area,mass,mass_old,mass_tol,*int1_coarea,coarea_l,coarea;
    double int1_hrms[100],hrms_l,hrms;
    double int1_hmean[100],hmean_l,hmean,hmean_tol;
    double *localp,*globalp;
    double *dxdt,*dxdp,*dydt,*dydp,*dzdt,*dzdp;

    int order = Geti("ahf_interpolation_order");
    int rank = bampi_rank();

    double *gxxg,*gxyg,*gxzg,*gyyg,*gyzg,*gzzg;
    double *dgxxdxg,*dgxydxg,*dgxzdxg,*dgyydxg,*dgyzdxg,*dgzzdxg;
    double *dgxxdyg,*dgxydyg,*dgxzdyg,*dgyydyg,*dgyzdyg,*dgzzdyg;
    double *dgxxdzg,*dgxydzg,*dgxzdzg,*dgyydzg,*dgyzdzg,*dgzzdzg;
    tVarList *vl,*wl;
    double *vinterp;
    int flag,dtime,u0;
    double time=level->time;
    double ahf_time=Getd("ahf_time");
    double Rx,Ry,Rz,**intSx,**intSy,**intSz;
    double *integral_Sx,*integral_Sy,*integral_Sz;
    double Sx_l,Sy_l,Sz_l,Sx,Sy,Sz,S;
    double phix_x,phix_y,phix_z;
    double phiy_x,phiy_y,phiy_z;
    double phiz_x,phiz_y,phiz_z;
    double v_l[7], v[7];

    double sum;

    char ext[101],name[201],*outdir = Gets("outdir");
    FILE *fp1;
    int common=0;
    int pr = !Getv("ahf_verbose", "no");

    //4/2010: Added variables for symmetry 
    double xc_sym,yc_sym,zc_sym;
    int fxc_sym,fyc_sym,fzc_sym;

    flow_iter = Geti("ahf_flow_iter");

    dtime = (time+dequaleps)/ahf_time;
    if(!dequal(time-dtime*ahf_time, 0)) return 0;
    timer_start(0, "AHF");

   // I commented this to get everything compiling, I do not know the influence! 
  //  interpolate_setorder(order);

//    if (Getv("grid", "bitant")) errorexit("AHF: not implemented in bitant");
//    if(Getv("grid", "quadrant")) errorexit("AHF: not implemented in quadrant");

    A = alpha/(LMAX1*(LMAX1+1))+beta; //0.01
    B = beta/alpha; //0.5

    mass_tol = Getd("ahf_mass_tol");
    hmean_tol = Getd("ahf_hmean_tol");

    nhorizons = Geti("ahf_nhorizons");
    if((nhorizons<1)||(nhorizons>2)) errorexit("AHF: ahf_nhorizons should be 1 or 2");

    if(time>=Getd("ahf_common_time")&&Getv("ahf_common", "yes")&&(level->nboxes==1)) common=1;

    nhorizons+=common;





    if ((level->l==Geti("amr_lmax")) || (level->l==Geti("amr_lmax2"))) {
      if (Geti("amr_lmax2")!=-1) 
        if ((level->l>Geti("amr_lmax")) || (level->l>Geti("amr_lmax2"))) return 0;
      search0 = 1;
      search1 = 1;
      search2 = 1;
    }
    if((nhorizons==1)&&(search0==0)) return(0);
    if((nhorizons==2)&&(search0==0)&&(search1==0)) return(0);
    if((nhorizons==3)&&(search0==0)&&(search1==0)&&(search2==0)) return(0);

    //Angular parameters for the surface
    
    ntheta = Geti("ahf_ntheta");
    nphi   = Geti("ahf_nphi");
	
    dtheta = pi/ntheta;
    dphi   = 2.0*pi/nphi;

    //Allocate

    a0  = dvector(0,LMAX1);
    ac  = dmatrix(0,LMAX1,0,LMAX1);
    as  = dmatrix(0,LMAX1,0,LMAX1);
    rr = dmatrix(0,ntheta-1,0,nphi-1);
    u  = dmatrix(0,ntheta-1,0,nphi-1);
    H  = dmatrix(0,ntheta-1,0,nphi-1);
    rho = dmatrix(0,ntheta-1,0,nphi-1);
    P  = dmatrix(0,LMAX1,0,LMAX1);
    dPdtheta = dmatrix(0,LMAX1,0,LMAX1);
    dPdthetadtheta = dmatrix(0,LMAX1,0,LMAX1);
    Y0 = dmatrix(0,ntheta*nphi,0,LMAX1);
    Yc = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    Ys = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dY0dtheta = dmatrix(0,ntheta*nphi,0,LMAX1);
    dYcdtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYsdtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYcdphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYsdphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dY0dthetadtheta = dmatrix(0,ntheta*nphi,0,LMAX1);
    dYcdthetadtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYcdthetadphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYcdphidphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYsdthetadtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYsdthetadphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dYsdphidphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    dFdx = dmatrix(0,ntheta-1,0,nphi-1);
    dFdy = dmatrix(0,ntheta-1,0,nphi-1);
    dFdz = dmatrix(0,ntheta-1,0,nphi-1);
    sigma = dmatrix(0,ntheta-1,0,nphi-1);
    int1_spec0 = dmatrix(0,ntheta,0,LMAX1);
    spec0 = dvector(0,LMAX1);
    spec0_l = dvector(0,LMAX1);
    int1_spec_c = dmatrix(0,ntheta,0,(LMAX1+1)*(LMAX1+1));
    spec_c = dvector(0,(LMAX1+1)*(LMAX1+1));
    spec_cl = dvector(0,(LMAX1+1)*(LMAX1+1));
    int1_spec_s = dmatrix(0,ntheta,0,(LMAX1+1)*(LMAX1+1));
    spec_s = dvector(0,(LMAX1+1)*(LMAX1+1));
    spec_sl = dvector(0,(LMAX1+1)*(LMAX1+1));
    deth = dmatrix(0,ntheta-1,0,nphi-1);
    int1_area = dvector(0,ntheta-1);
    int1_coarea = dvector(0,ntheta-1);
    localp = dvector(0,ntheta*nphi);
    globalp = dvector(0,ntheta*nphi);
    dxdt = dvector(0,ntheta*nphi);
    dxdp = dvector(0,ntheta*nphi);
    dydt = dvector(0,ntheta*nphi);
    dydp = dvector(0,ntheta*nphi);
    dzdt = dvector(0,ntheta*nphi);
    dzdp = dvector(0,ntheta*nphi);
    intSx = dmatrix(0,ntheta-1,0,nphi-1);
    intSy = dmatrix(0,ntheta-1,0,nphi-1);
    intSz = dmatrix(0,ntheta-1,0,nphi-1);
    integral_Sx = dvector(0,ntheta-1);
    integral_Sy = dvector(0,ntheta-1);
    integral_Sz = dvector(0,ntheta-1);


    //Compute legendre and spherical harmonics only once, 
    //because the spherical grid is the same for all the horizons.
    //Loop over the surface

    for(i=0;i<ntheta;i++){
      theta = dtheta*(1.0/2.0 + i);
      for(j=0;j<nphi;j++){
	
	phi = dphi*(1.0/2.0 + j);

	i1=i*nphi+j;

	//set P dPdtheta dPdthtetadtheta to 0.0

	for(l=0;l<=LMAX1;l++){
	  for(m=0;m<LMAX1;m++){
	    P[l][m] = 0.0;
	    dPdtheta[l][m] = 0.0;
	    dPdthetadtheta[l][m] = 0.0;
	  }	
	}

	//Compute the Legendre functions

	facl  = 1.0;
	fac2l = 1.0;

	for(l=0;l<=LMAX1;l++){
	  if(l>0){
	    facl  *= l; 
	    fac2l *= 2*l*(2*l-1); 
	  } 
	  P[l][l] = sqrt((2*l+1)*fac2l/(4.0*pi))/(pow(2,l)*facl)*pow((-sin(theta)),l);
	}

	P[1][0] = sqrt(3)*cos(theta)*P[0][0];

	for(l=2;l<=LMAX1;l++){
	  for(m=0;m<l;m++){
            if (m<l-1) // if added NL 01/08
  	      P[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*cos(theta)*P[l-1][m]
		  			       - sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*P[l-2][m]);
            else   
              P[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*cos(theta)*P[l-1][m]);
	  }
	}

	//Compute first derivatives of the Legendre functions

	facl  = 1.0;
	fac2l = 1.0;
	  
	for(l=0;l<=LMAX1;l++){
	  if(l>0){
	    facl  *= l; 
	    fac2l *= 2*l*(2*l-1); 
	  } 
	  dPdtheta[l][l] = sqrt((2*l+1)*fac2l/(4.0*pi))/(pow(2,l)*facl)*l*pow((-sin(theta)),l-1)*(-cos(theta));
	}
	  
	dPdtheta[1][0] = sqrt(3)*(-sin(theta)*P[0][0]+cos(theta)*dPdtheta[0][0]);

	for(l=2;l<=LMAX1;l++){
	  for(m=0;m<l;m++){
	    dPdtheta[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-sin(theta)*P[l-1][m]
			   + cos(theta)*dPdtheta[l-1][m])-sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*dPdtheta[l-2][m]);
	  }
	}

	//Compute second derivatives of the Legendre functions

	facl  = 1.0;
	fac2l = 1.0;
	
	for(l=0;l<=LMAX1;l++){
	  if(l>0){
	    facl  *= l; 
	    fac2l *= 2*l*(2*l-1); 
	  } 
	  dPdthetadtheta[l][l] = sqrt((2*l+1)*fac2l/(4.0*pi))/(pow(2,l)*facl)*l
	                       *((l-1)*pow((-sin(theta)),l-2)*cos(theta)*cos(theta) 
			       + pow((-sin(theta)),l-1)*sin(theta));
	}

	dPdthetadtheta[1][0] = sqrt(3)*(-cos(theta)*P[0][0]-2.0*sin(theta)*dPdtheta[0][0]
         		     + cos(theta)*dPdthetadtheta[0][0]);

	for(l=2;l<=LMAX1;l++){
	  for(m=0;m<l;m++){
	    dPdthetadtheta[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-cos(theta)*P[l-1][m]
				 - 2.0*sin(theta)*dPdtheta[l-1][m] + cos(theta)*dPdthetadtheta[l-1][m])
 			         - sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*dPdthetadtheta[l-2][m]);
	  }
	}

	//Compute the spherical harmonics
	
	for(l=0;l<=LMAX1;l++){
	  Y0[i1][l] = P[l][0];
	}

	for(l=1;l<=LMAX1;l++){
	  for(m=1;m<=l;m++){
	    l1=l*(LMAX1+1)+m;
	    Yc[i1][l1] = sqrt(2.0)*P[l][m]*cos(m*phi);
	    Ys[i1][l1] = sqrt(2.0)*P[l][m]*sin(m*phi);
	  }
	}

	//Compute first derivatives of the spherical harmonics

	for(l=0;l<=LMAX1;l++){
	  dY0dtheta[i1][l] = dPdtheta[l][0];
	}

	for(l=1;l<=LMAX1;l++){
	  for(m=1;m<=l;m++){
	    l1=l*(LMAX1+1)+m;
	    dYcdtheta[i1][l1] = sqrt(2.0)*dPdtheta[l][m]*cos(m*phi);
	    dYsdtheta[i1][l1] = sqrt(2.0)*dPdtheta[l][m]*sin(m*phi);
	    dYcdphi[i1][l1] = -sqrt(2.0)*P[l][m]*m*sin(m*phi);
	    dYsdphi[i1][l1] =  sqrt(2.0)*P[l][m]*m*cos(m*phi);
	  }
	}
	
	//Compute second derivatives of the spherical harmonics
	
	for(l=0;l<=LMAX1;l++){
	  dY0dthetadtheta[i1][l] = dPdthetadtheta[l][0];
	}
	  
	for(l=1;l<=LMAX1;l++){
	  for(m=1;m<=l;m++){
	    l1=l*(LMAX1+1)+m;
	    dYcdthetadtheta[i1][l1] = sqrt(2.0)*dPdthetadtheta[l][m]*cos(m*phi);
	    dYcdthetadphi[i1][l1] = -sqrt(2.0)*dPdtheta[l][m]*m*sin(m*phi);
	    dYsdthetadtheta[i1][l1] = sqrt(2.0)*dPdthetadtheta[l][m]*sin(m*phi);
	    dYsdthetadphi[i1][l1] = sqrt(2.0)*dPdtheta[l][m]*m*cos(m*phi);
	    dYcdphidphi[i1][l1] = -sqrt(2.0)*P[l][m]*m*m*cos(m*phi);
	    dYsdphidphi[i1][l1] = -sqrt(2.0)*P[l][m]*m*m*sin(m*phi);
	  }
	}
      }
    }

    //dx dy dz

    dx = level->dx;
    dy = level->dy;
    dz = level->dz;

    idx = 1.0/(2.0*dx);
    idy = 1.0/(2.0*dy);
    idz = 1.0/(2.0*dz);

    //Compute the derivatives of the metric in the level
    //and store it in the variables dgxxdxg,....
    vl = vlalloc(level);
    vlpush(vl,Ind("ahf_dgdxxx"));

    enablevarlist(vl);

    gxxg = Ptr(level, "adm_gxx");
    gxyg = Ptr(level, "adm_gxy");
    gxzg = Ptr(level, "adm_gxz");
    gyyg = Ptr(level, "adm_gyy");
    gyzg = Ptr(level, "adm_gyz");
    gzzg = Ptr(level, "adm_gzz");
    dgxxdxg = vldataptr(vl,0);
    dgxydxg = vldataptr(vl,1);
    dgxzdxg = vldataptr(vl,2);
    dgyydxg = vldataptr(vl,3);
    dgyzdxg = vldataptr(vl,4);
    dgzzdxg = vldataptr(vl,5);
    dgxxdyg = vldataptr(vl,6);
    dgxydyg = vldataptr(vl,7);
    dgxzdyg = vldataptr(vl,8);
    dgyydyg = vldataptr(vl,9);
    dgyzdyg = vldataptr(vl,10);
    dgzzdyg = vldataptr(vl,11);
    dgxxdzg = vldataptr(vl,12);
    dgxydzg = vldataptr(vl,13);
    dgxzdzg = vldataptr(vl,14);
    dgyydzg = vldataptr(vl,15);
    dgyzdzg = vldataptr(vl,16);
    dgzzdzg = vldataptr(vl,17);

    wl = vlalloc(level);
    vlpush(wl,Ind("adm_gxx"));
    vlpush(wl,Ind("adm_Kxx"));
    vlpushvl(wl,vl);
    
    vinterp = dmalloc(wl->n);

    forinnerpoints_ijk(level){

      dgxxdxg[ijk] =  idx*(gxxg[ijk+di] - gxxg[ijk-di]);
      dgxydxg[ijk] =  idx*(gxyg[ijk+di] - gxyg[ijk-di]);
      dgxzdxg[ijk] =  idx*(gxzg[ijk+di] - gxzg[ijk-di]);
      dgyydxg[ijk] =  idx*(gyyg[ijk+di] - gyyg[ijk-di]);
      dgyzdxg[ijk] =  idx*(gyzg[ijk+di] - gyzg[ijk-di]);
      dgzzdxg[ijk] =  idx*(gzzg[ijk+di] - gzzg[ijk-di]);
      dgxxdyg[ijk] =  idy*(gxxg[ijk+dj] - gxxg[ijk-dj]);
      dgxydyg[ijk] =  idy*(gxyg[ijk+dj] - gxyg[ijk-dj]);
      dgxzdyg[ijk] =  idy*(gxzg[ijk+dj] - gxzg[ijk-dj]);
      dgyydyg[ijk] =  idy*(gyyg[ijk+dj] - gyyg[ijk-dj]);
      dgyzdyg[ijk] =  idy*(gyzg[ijk+dj] - gyzg[ijk-dj]);
      dgzzdyg[ijk] =  idy*(gzzg[ijk+dj] - gzzg[ijk-dj]);
      dgxxdzg[ijk] =  idz*(gxxg[ijk+dk] - gxxg[ijk-dk]);
      dgxydzg[ijk] =  idz*(gxyg[ijk+dk] - gxyg[ijk-dk]);
      dgxzdzg[ijk] =  idz*(gxzg[ijk+dk] - gxzg[ijk-dk]);
      dgyydzg[ijk] =  idz*(gyyg[ijk+dk] - gyyg[ijk-dk]);
      dgyzdzg[ijk] =  idz*(gyzg[ijk+dk] - gyzg[ijk-dk]);
      dgzzdzg[ijk] =  idz*(gzzg[ijk+dk] - gzzg[ijk-dk]);

    }endfor_ijk;

    /* set boundaries */
    set_boundary_symmetry(level, vl);
    
    /* synchronize */
    bampi_vlsynchronize(vl);

    if(find_common){
      search0 = 0;
      if(nhorizons==3) search1=0;
    }

    //Start loop over number of horizons
    for(number=0;number<nhorizons;number++){
      if((number==0&&search0==1)||(number==1&&search1==1)||(number==2&&search2==1)){

	mass=0.0;
	hmean=0.0;

	if(rank==0) printf("Searching for horizon %d in level %d at time %f\n",number,level->l,time);

	//Initialize a_lm
	
	for(l=0;l<=LMAX1;l++){
	  a0[l] = 0.0; 
	}
	
	for(l=0;l<=LMAX1;l++){
	  for(m=0;m<=l;m++){
	    ac[l][m] = 0.0; 
	    as[l][m] = 0.0; 
	  }
	}

	//Position of the BH

	if(number==0){
	  xc = Getd("moving_puncture1_x");
	  yc = Getd("moving_puncture1_y");
	  zc = Getd("moving_puncture1_z");
	}else if(number==1){
	  if(nhorizons==3){
	    xc = Getd("moving_puncture2_x");
	    yc = Getd("moving_puncture2_y");
	    zc = Getd("moving_puncture2_z");
	  }else if(nhorizons==2 &&common==0){
	    xc = Getd("moving_puncture2_x");
	    yc = Getd("moving_puncture2_y");
	    zc = Getd("moving_puncture2_z");
	  }else{
	    xc = 0.0;
	    yc = 0.0;
	    zc = 0.0;
	  }
	}else if(number==2){
	  xc = 0.0;
	  yc = 0.0;
	  zc = 0.0;
	}

	if(rank==0) printf("centered in:  x = %f  y = %f  z = %f\n",xc,yc,zc);

	//4/2010: Implemented to ensure BH is in symmetry-produced box (Ignored if "grid=full")
	
	xc_sym = xc;
	yc_sym = yc;
	zc_sym = zc;
	
	map_xyz_withsym(level,&xc_sym,&yc_sym,&zc_sym,&fxc_sym,&fyc_sym,&fzc_sym);

	//Limits of the box containing the BH.
	//4/2010: Edited fbox for "xc/yc/zc"_sym 

	fbox = box_containing_xyz(level, xc_sym, yc_sym, zc_sym);

	xmn = fbox->bbox[0];
	xmx = fbox->bbox[1];
	ymn = fbox->bbox[2];
	ymx = fbox->bbox[3];
	zmn = fbox->bbox[4];
	zmx = fbox->bbox[5];

	//Initial guess for a_lm

	r_max = min2(xc-xmn,xmx-xc);  //for rotant only take the maximum x
	//r_max = min2(r_max,yc-ymn);
	//r_max = min2(r_max,ymx-yc);
	//r_max = min2(r_max,zc-zmn);
	//r_max = min2(r_max,zmx-zc);

	r_init = 0.5*r_max;
	
	a0[0] = sqrt(4.0*pi)*r_init;

	//printf("\n");
	//printf("INITIAL GUESS\n");
	//printf("\n");
	//printf("radius = %22.16e\n",r_init);
	//printf("\n");

	//Flow loop      
	for(k=0;k<flow_iter;k++){

	  if (pr) {
	    printf("\n");
	    printf("ITERATION %d\n",k);
	    printf("\n");
	  }

	  //Initialize
	  for(i=0;i<ntheta;i++){
	    for(j=0;j<nphi;j++){
	      rr[i][j] = 0.0; 
	      u[i][j]  = 0.0; 
	      H[i][j]  = 0.0; 
	      rho[i][j] = 0.0;
	      deth[i][j] = 0.0;
	    }
	  }

	  //Compute the radius of the surface r=a_lm Y_lm
	  //loop over the surface
	  cont = 1;
	  for(i=0;i<ntheta;i++){
	    theta = dtheta*(1.0/2.0 + i);
	    for(j=0;j<nphi;j++){
	      phi = dphi*(1.0/2.0 + j);
	      i1=i*nphi+j;
	      for(l=0;l<=LMAX1;l++){
		rr[i][j] += a0[l]*Y0[i1][l];
	      }

	      for(l=1;l<=LMAX1;l++){
		for(m=1;m<=l;m++){
		  rr[i][j] += Yc[i1][l*(LMAX1+1)+m]*ac[l][m] + Ys[i1][l*(LMAX1+1)+m]*as[l][m];
		  if(rr[i][j]>r_max) cont=0;
		  if(rr[i][j]<0) cont=0;
		}
	      }
	    }
	  }

          // is this safe for parallelization? if not, try 
          cont = bampi_allreduce_min_int(cont);
	  if(!cont){
	    if(rank==0){
	      printf("AHF failed!\n");
	      printf("r>rmax or r<0\n");
	      printf("\n");
	    }
	    break;
	  }	 

	  //Compute derivatives of (x,y,z) with respect to theta and phi
	  for(i=0;i<ntheta;i++){
	    theta = dtheta*(1.0/2.0 + i);
	    for(j=0;j<nphi;j++){
	      phi = dphi*(1.0/2.0 + j);
	    
	      i1=i*nphi+j;
	    
	      if ((nphi+1)%2==0) errorexit("AHF: odd number of points");
	      //if (ntheta%2==0) errorexit("AHF: even number of circles");

	      half = (int)(nphi/2);
	      if(i==0){
		rtp1 = rr[1][j];
		if(j<half)  rtm1 = rr[0][j+half];
		if(j>=half) rtm1 = rr[0][j-half];
	      }else if(i==ntheta-1){
		if(j<half)  rtp1 = rr[ntheta-1][j+half];
		if(j>=half) rtp1 = rr[ntheta-1][j-half];
		rtm1 = rr[ntheta-2][j];
	      }else{
		rtp1 = rr[i+1][j];
		rtm1 = rr[i-1][j];
	      }
	  
	      if(j==0){
		rpp1 = rr[i][1];
		rpm1 = rr[i][nphi-1];
	      }else if(j==nphi-1){
		rpp1 = rr[i][0];
		rpm1 = rr[i][nphi-2];
	      }else{
		rpp1 = rr[i][j+1];
		rpm1 = rr[i][j-1];
	      }
		    
	      drdt = (rtp1-rtm1)/(2.0*dtheta);
	      drdp = (rpp1-rpm1)/(2.0*dphi);

	      //Derivatives of (x,y,z) with respect to theta
	      dxdt[i1] = (drdt*sin(theta) + rr[i][j]*cos(theta))*cos(phi);
	      dydt[i1] = (drdt*sin(theta) + rr[i][j]*cos(theta))*sin(phi);
	      dzdt[i1] =  drdt*cos(theta) - rr[i][j]*sin(theta);
		    
	      //Derivatives of (x,y,z) with respect to phi
	      dxdp[i1] = (drdp*cos(phi) - rr[i][j]*sin(phi))*sin(theta);
	      dydp[i1] = (drdp*sin(phi) + rr[i][j]*cos(phi))*sin(theta);
	      dzdp[i1] =  drdp*cos(theta);
	    }
	  }

	  timer_start(0, "AHF_local");
	
          
          
          
          //loop over the surface
	  for(i=0;i<ntheta;i++){
	    theta = dtheta*(1.0/2.0 + i);
	    for(j=0;j<nphi;j++){
	      phi = dphi*(1.0/2.0 + j);
	     
	      i1=i*nphi+j;
	      
	      //global coordinates of the surface

	      x = xc + rr[i][j]*sin(theta)*cos(phi);
	      y = yc + rr[i][j]*sin(theta)*sin(phi);
	      z = zc + rr[i][j]*cos(theta);
	       
             
	      //Interpolate the extrinsic curvature, the 3 metric and the
	      //derivatives of the 3 metric over the surface
	      //The interpolation is only done on the points in this processor

	      flag = check_interpolation_cube_local_withsym(level,x,y,z,order);

	      if(flag){

		interpolate_xyz_local_minimal_withsym(level,x,y,z,wl->n,wl->index,vinterp,order,LAGRANGE);

		nv = 0;

		gxx = vinterp[nv++];
		gxy = vinterp[nv++];
		gxz = vinterp[nv++];
		gyy = vinterp[nv++];
		gyz = vinterp[nv++];
		gzz = vinterp[nv++];

		Kxx = vinterp[nv++];
		Kxy = vinterp[nv++];
		Kxz = vinterp[nv++];
		Kyy = vinterp[nv++];
		Kyz = vinterp[nv++];
		Kzz = vinterp[nv++];
	
		dgxxdx = vinterp[nv++];
		dgxydx = vinterp[nv++];
		dgxzdx = vinterp[nv++];
		dgyydx = vinterp[nv++];
		dgyzdx = vinterp[nv++];
		dgzzdx = vinterp[nv++];
		dgxxdy = vinterp[nv++];
		dgxydy = vinterp[nv++];
		dgxzdy = vinterp[nv++];
		dgyydy = vinterp[nv++];
		dgyzdy = vinterp[nv++];
		dgzzdy = vinterp[nv++];
		dgxxdz = vinterp[nv++];
		dgxydz = vinterp[nv++];
		dgxzdz = vinterp[nv++];
		dgyydz = vinterp[nv++];
		dgyzdz = vinterp[nv++];
		dgzzdz = vinterp[nv++];

		//Inverse 3 metric.

		ginvdetxx = gyy*gzz - gyz*gyz;
		ginvdetxy = gxz*gyz - gxy*gzz;
		ginvdetxz = gxy*gyz - gxz*gyy;
		ginvdetyy = gxx*gzz - gxz*gxz;
		ginvdetyz = gxy*gxz - gxx*gyz;
		ginvdetzz = gxx*gyy - gxy*gxy;
	      
		detg = gxx*ginvdetxx + gxy*ginvdetxy + gxz*ginvdetxz;
	      
		ginvxx = ginvdetxx/detg;
		ginvxy = ginvdetxy/detg;
		ginvxz = ginvdetxz/detg;
		ginvyx = ginvdetxy/detg;
		ginvyy = ginvdetyy/detg;
		ginvyz = ginvdetyz/detg;
		ginvzx = ginvdetxz/detg;
		ginvzy = ginvdetyz/detg;
		ginvzz = ginvdetzz/detg;

		//Trace of K

		K = ginvxx*Kxx + ginvxy*Kxy + ginvxz*Kxz
		  + ginvyx*Kxy + ginvyy*Kyy + ginvyz*Kyz
		  + ginvzx*Kxz + ginvzy*Kyz + ginvzz*Kzz;
	      
		//local coordinates of the surface

		xp   = x-xc;
		yp   = y-yc;
		zp   = z-zc;
		rp   = sqrt(xp*xp + yp*yp + zp*zp);
		rhop = sqrt(xp*xp + yp*yp);
	      
		//first derivatives of (r,theta,phi) with respect to (x,y,z)

		drdx = xp/rp;
		drdy = yp/rp;
		drdz = zp/rp;

		dthetadx = zp*xp/(rp*rp*rhop);
		dthetady = zp*yp/(rp*rp*rhop);
		dthetadz = -rhop/(rp*rp);
	      
		dphidx = -yp/(rhop*rhop);
		dphidy = xp/(rhop*rhop);
		dphidz = 0.0;

		//second derivatives of (r,theta,phi) with respect to (x,y,z)

		drdxdx = 1.0/rp - xp*xp/(rp*rp*rp);
		drdxdy = - xp*yp/(rp*rp*rp);
		drdxdz = - xp*zp/(rp*rp*rp);
		drdydy = 1.0/rp - yp*yp/(rp*rp*rp);
		drdydz = - yp*zp/(rp*rp*rp);
		drdzdz = 1.0/rp - zp*zp/(rp*rp*rp);
	    
		dthetadxdx = zp*(-2.0*xp*xp*xp*xp-xp*xp*yp*yp+yp*yp*yp*yp+zp*zp*yp*yp)/(rp*rp*rp*rp*rhop*rhop*rhop);
		dthetadxdy = - xp*yp*zp*(3.0*xp*xp+3.0*yp*yp+zp*zp)/(rp*rp*rp*rp*rhop*rhop*rhop);
		dthetadxdz = xp*(xp*xp+yp*yp-zp*zp)/(rp*rp*rp*rp*rhop);
		dthetadydy = zp*(-2.0*yp*yp*yp*yp-yp*yp*xp*xp+xp*xp*xp*xp+zp*zp*xp*xp)/(rp*rp*rp*rp*rhop*rhop*rhop);
		dthetadydz = yp*(xp*xp+yp*yp-zp*zp)/(rp*rp*rp*rp*rhop);
		dthetadzdz = 2.0*zp*rhop/(rp*rp*rp*rp);
	      
		dphidxdx = 2.0*yp*xp/(rhop*rhop*rhop*rhop);
		dphidxdy = (yp*yp-xp*xp)/(rhop*rhop*rhop*rhop);
		dphidxdz = 0.0;
		dphidydy = - 2.0*yp*xp/(rhop*rhop*rhop*rhop);
		dphidydz = 0.0;
		dphidzdz = 0.0;
	    
		//Compute first derivatives of F: dFdi

		dFdx[i][j] = drdx;
		dFdy[i][j] = drdy;
		dFdz[i][j] = drdz;

		for(l=0;l<=LMAX1;l++){
		  dFdx[i][j] -= a0[l]*dthetadx*dY0dtheta[i1][l]; 
		  dFdy[i][j] -= a0[l]*dthetady*dY0dtheta[i1][l]; 
		  dFdz[i][j] -= a0[l]*dthetadz*dY0dtheta[i1][l]; 
		}
		
		for(l=1;l<=LMAX1;l++){
		  for(m=1;m<=l;m++){
		    l1=l*(LMAX1+1)+m;
		    dFdx[i][j] -= ac[l][m]*(dthetadx*dYcdtheta[i1][l1] + dphidx*dYcdphi[i1][l1]) 
		                + as[l][m]*(dthetadx*dYsdtheta[i1][l1] + dphidx*dYsdphi[i1][l1]); 
		    dFdy[i][j] -= ac[l][m]*(dthetady*dYcdtheta[i1][l1] + dphidy*dYcdphi[i1][l1]) 
		                + as[l][m]*(dthetady*dYsdtheta[i1][l1] + dphidy*dYsdphi[i1][l1]); 
		    dFdz[i][j] -= ac[l][m]*(dthetadz*dYcdtheta[i1][l1] + dphidz*dYcdphi[i1][l1]) 
		                + as[l][m]*(dthetadz*dYsdtheta[i1][l1] + dphidz*dYsdphi[i1][l1]); 
		  }
		}

		//Compute second derivatives of F: dFdidj

		dFdxdx = drdxdx;
		dFdxdy = drdxdy;
		dFdxdz = drdxdz;
		dFdydy = drdydy;
		dFdydz = drdydz;
		dFdzdz = drdzdz;
	      
		for(l=0;l<=LMAX1;l++){
		  dFdxdx -= a0[l]*(dthetadxdx*dY0dtheta[i1][l] + dthetadx*dthetadx*dY0dthetadtheta[i1][l]); 
		  dFdxdy -= a0[l]*(dthetadxdy*dY0dtheta[i1][l] + dthetadx*dthetady*dY0dthetadtheta[i1][l]); 
		  dFdxdz -= a0[l]*(dthetadxdz*dY0dtheta[i1][l] + dthetadx*dthetadz*dY0dthetadtheta[i1][l]); 
		  dFdydx -= a0[l]*(dthetadxdy*dY0dtheta[i1][l] + dthetady*dthetadx*dY0dthetadtheta[i1][l]); 
		  dFdydy -= a0[l]*(dthetadydy*dY0dtheta[i1][l] + dthetady*dthetady*dY0dthetadtheta[i1][l]); 
		  dFdydz -= a0[l]*(dthetadydz*dY0dtheta[i1][l] + dthetady*dthetadz*dY0dthetadtheta[i1][l]);
		  dFdzdx -= a0[l]*(dthetadxdz*dY0dtheta[i1][l] + dthetadz*dthetadx*dY0dthetadtheta[i1][l]); 
		  dFdzdy -= a0[l]*(dthetadydz*dY0dtheta[i1][l] + dthetadz*dthetady*dY0dthetadtheta[i1][l]);  
		  dFdzdz -= a0[l]*(dthetadzdz*dY0dtheta[i1][l] + dthetadz*dthetadz*dY0dthetadtheta[i1][l]); 
		}
	    
		for(l=1;l<=LMAX1;l++){
		  for(m=1;m<=l;m++){
		    l1=l*(LMAX1+1)+m;
		    dFdxdx -= ac[l][m]*(dthetadxdx*dYcdtheta[i1][l1] + dthetadx*(dthetadx*dYcdthetadtheta[i1][l1] 
		            + dphidx*dYcdthetadphi[i1][l1]) + dphidxdx*dYcdphi[i1][l1] 
			    + dphidx*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1])) 
			    + as[l][m]*(dthetadxdx*dYsdtheta[i1][l1] + dthetadx*(dthetadx*dYsdthetadtheta[i1][l1] 
                            + dphidx*dYsdthetadphi[i1][l1]) + dphidxdx*dYsdphi[i1][l1] 
                            + dphidx*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));
		    dFdxdy -= ac[l][m]*(dthetadxdy*dYcdtheta[i1][l1] + dthetadx*(dthetady*dYcdthetadtheta[i1][l1] 
                            + dphidy*dYcdthetadphi[i1][l1]) + dphidxdy*dYcdphi[i1][l1] 
                            + dphidx*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1])) 
		  	    + as[l][m]*(dthetadxdy*dYsdtheta[i1][l1] + dthetadx*(dthetady*dYsdthetadtheta[i1][l1] 
                            + dphidy*dYsdthetadphi[i1][l1]) + dphidxdy*dYsdphi[i1][l1] 
                            + dphidx*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));
		    dFdxdz -= ac[l][m]*(dthetadxdz*dYcdtheta[i1][l1] + dthetadx*(dthetadz*dYcdthetadtheta[i1][l1] 
                            + dphidz*dYcdthetadphi[i1][l1]) + dphidxdz*dYcdphi[i1][l1] 
                            + dphidx*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1])) 
                            + as[l][m]*(dthetadxdz*dYsdtheta[i1][l1] + dthetadx*(dthetadz*dYsdthetadtheta[i1][l1] 
                            + dphidz*dYsdthetadphi[i1][l1]) + dphidxdz*dYsdphi[i1][l1] 
                            + dphidx*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));
		    dFdydx -= ac[l][m]*(dthetadxdy*dYcdtheta[i1][l1] + dthetady*(dthetadx*dYcdthetadtheta[i1][l1] 
 		            + dphidx*dYcdthetadphi[i1][l1]) + dphidxdy*dYcdphi[i1][l1] 
			    + dphidy*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1])) 
			    + as[l][m]*(dthetadxdy*dYsdtheta[i1][l1] + dthetady*(dthetadx*dYsdthetadtheta[i1][l1] 
                            + dphidx*dYsdthetadphi[i1][l1]) + dphidxdy*dYsdphi[i1][l1] 
                            + dphidy*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));
		    dFdydy -= ac[l][m]*(dthetadydy*dYcdtheta[i1][l1] + dthetady*(dthetady*dYcdthetadtheta[i1][l1] 
                            + dphidy*dYcdthetadphi[i1][l1]) + dphidydy*dYcdphi[i1][l1] 
                            + dphidy*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1])) 
			    + as[l][m]*(dthetadydy*dYsdtheta[i1][l1] + dthetady*(dthetady*dYsdthetadtheta[i1][l1] 
                            + dphidy*dYsdthetadphi[i1][l1]) + dphidydy*dYsdphi[i1][l1] 
                            + dphidy*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));
		    dFdydz -= ac[l][m]*(dthetadydz*dYcdtheta[i1][l1] + dthetady*(dthetadz*dYcdthetadtheta[i1][l1] 
                            + dphidz*dYcdthetadphi[i1][l1]) + dphidydz*dYcdphi[i1][l1] 
                            + dphidy*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1])) 
                            + as[l][m]*(dthetadydz*dYsdtheta[i1][l1] + dthetady*(dthetadz*dYsdthetadtheta[i1][l1] 
                            + dphidz*dYsdthetadphi[i1][l1]) + dphidydz*dYsdphi[i1][l1] 
                            + dphidy*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));
		    dFdzdx -= ac[l][m]*(dthetadxdz*dYcdtheta[i1][l1] + dthetadz*(dthetadx*dYcdthetadtheta[i1][l1] 
		            + dphidx*dYcdthetadphi[i1][l1]) + dphidxdz*dYcdphi[i1][l1] 
                            + dphidz*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1])) 
			    + as[l][m]*(dthetadxdz*dYsdtheta[i1][l1] + dthetadz*(dthetadx*dYsdthetadtheta[i1][l1] 
                            + dphidx*dYsdthetadphi[i1][l1]) + dphidxdz*dYsdphi[i1][l1] 
                            + dphidz*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));
		    dFdzdy -= ac[l][m]*(dthetadydz*dYcdtheta[i1][l1] + dthetadz*(dthetady*dYcdthetadtheta[i1][l1] 
                            + dphidy*dYcdthetadphi[i1][l1]) + dphidydz*dYcdphi[i1][l1] 
                            + dphidz*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1])) 
			    + as[l][m]*(dthetadydz*dYsdtheta[i1][l1] + dthetadz*(dthetady*dYsdthetadtheta[i1][l1] 
                            + dphidy*dYsdthetadphi[i1][l1]) + dphidydz*dYsdphi[i1][l1] 
                            + dphidz*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));
		    dFdzdz -= ac[l][m]*(dthetadzdz*dYcdtheta[i1][l1] + dthetadz*(dthetadz*dYcdthetadtheta[i1][l1] 
                            + dphidz*dYcdthetadphi[i1][l1]) + dphidzdz*dYcdphi[i1][l1] 
                            + dphidz*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1])) 
                            + as[l][m]*(dthetadzdz*dYsdtheta[i1][l1] + dthetadz*(dthetadz*dYsdthetadtheta[i1][l1] 
                            + dphidz*dYsdthetadphi[i1][l1]) + dphidzdz*dYsdphi[i1][l1] 
                            + dphidz*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));
		  }
		}

		//Compute dFdi with the index up
		
		dFupdx = ginvxx*dFdx[i][j] + ginvxy*dFdy[i][j] + ginvxz*dFdz[i][j];
		dFupdy = ginvyx*dFdx[i][j] + ginvyy*dFdy[i][j] + ginvyz*dFdz[i][j];
		dFupdz = ginvzx*dFdx[i][j] + ginvzy*dFdy[i][j] + ginvzz*dFdz[i][j];

		//Compute norm of dFdi

		sum = dFupdx*dFdx[i][j] + dFupdy*dFdy[i][j] + dFupdz*dFdz[i][j];

		if(sum>=0){
		  u[i][j] = sqrt(dFupdx*dFdx[i][j] + dFupdy*dFdy[i][j] + dFupdz*dFdz[i][j]);
		  u0=0;
		}else{
		  u[i][j] = 0.0;
		  u0=1;
		}

		//Compute nabla_a nabla_b F
	    
		nnFxx = dFdxdx - 1.0/2.0*(dFupdx*dgxxdx + dFupdy*(2.0*dgxydx-dgxxdy) + dFupdz*(2.0*dgxzdx-dgxxdz));
		nnFxy = dFdxdy - 1.0/2.0*(dFupdx*dgxxdy + dFupdy*dgyydx + dFupdz*(dgyzdx+dgxzdy-dgxydz));
		nnFxz = dFdxdz - 1.0/2.0*(dFupdx*dgxxdz + dFupdy*(dgyzdx + dgxydz - dgxzdy) + dFupdz*dgzzdx);
		nnFyy = dFdydy - 1.0/2.0*(dFupdy*dgyydy + dFupdx*(2.0*dgxydy - dgyydx) + dFupdz*(2.0*dgyzdy-dgyydz));
		nnFyz = dFdydz - 1.0/2.0*(dFupdy*dgyydz + dFupdx*(dgxzdy + dgxydz - dgyzdx) + dFupdz*dgzzdy);
		nnFzz = dFdzdz - 1.0/2.0*(dFupdz*dgzzdz + dFupdx*(2.0*dgxzdz - dgzzdx) + dFupdy*(2.0*dgyzdz-dgzzdy));

		//Compute d2F = g^{ab} nabla_a nabla_b F 

		d2F = ginvxx*nnFxx + ginvxy*nnFxy + ginvxz*nnFxz + ginvyx*nnFxy + ginvyy*nnFyy + ginvyz*nnFyz
  	            + ginvzx*nnFxz + ginvzy*nnFyz + ginvzz*nnFzz; 

		//Compute dFd^a dFd^b Kab
		
		dFdadFdbKab = dFupdx*dFupdx*Kxx + 2.0*dFupdx*dFupdy*Kxy + 2.0*dFupdx*dFupdz*Kxz
		            + dFupdy*dFupdy*Kyy + 2.0*dFupdy*dFupdz*Kyz + dFupdz*dFupdz*Kzz;

		//Compute dFd^a dFd^b nabla_a nabla_b F
		
		dFdadFdbFdadb = dFupdx*dFupdx*nnFxx + 2.0*dFupdx*dFupdy*nnFxy + 2.0*dFupdx*dFupdz*nnFxz 
                              + dFupdy*dFupdy*nnFyy + 2.0*dFupdy*dFupdz*nnFyz + dFupdz*dFupdz*nnFzz;

		//Compute H

		if(u0==0){
		  H[i][j] = d2F/u[i][j] + dFdadFdbKab/(u[i][j]*u[i][j]) - dFdadFdbFdadb/(u[i][j]*u[i][j]*u[i][j]) - K;
		}else{
		  H[i][j] = 0.0;
		}

		//Compute sigma
	      
		sigma[i][j] = 1.0;

		//Compute rho = H*u*sigma;

		rho[i][j] = H[i][j]*u[i][j]*sigma[i][j];

		//Induced metric on the horizon.
		h11 = dxdt[i1]*(dxdt[i1]*gxx+dydt[i1]*gxy+dzdt[i1]*gxz) 
                    + dydt[i1]*(dxdt[i1]*gxy+dydt[i1]*gyy+dzdt[i1]*gyz) 
                    + dzdt[i1]*(dxdt[i1]*gxz+dydt[i1]*gyz+dzdt[i1]*gzz);
		h12 = dxdt[i1]*(dxdp[i1]*gxx+dydp[i1]*gxy+dzdp[i1]*gxz) 
                    + dydt[i1]*(dxdp[i1]*gxy+dydp[i1]*gyy+dzdp[i1]*gyz) 
                    + dzdt[i1]*(dxdp[i1]*gxz+dydp[i1]*gyz+dzdp[i1]*gzz);
		h22 = dxdp[i1]*(dxdp[i1]*gxx+dydp[i1]*gxy+dzdp[i1]*gxz) 
                    + dydp[i1]*(dxdp[i1]*gxy+dydp[i1]*gyy+dzdp[i1]*gyz) 
                    + dzdp[i1]*(dxdp[i1]*gxz+dydp[i1]*gyz+dzdp[i1]*gzz);
	      
		//Determinant of the induced metric.
		deth[i][j] = h11*h22-h12*h12;
		if(deth[i][j]<0) deth[i][j] = 0.0;

		//Flat space coordinate rotational killing vectors.

		phix_x =  0;
		phix_y = -(z-zc);
		phix_z =  (y-yc);
		phiy_x =  (z-zc);
		phiy_y =  0;
		phiy_z = -(x-xc);
		phiz_x = -(y-yc);
		phiz_y =  (x-xc);
		phiz_z =  0;
	
		//Normal vector

		Rx = dFupdx/u[i][j];
		Ry = dFupdy/u[i][j];
		Rz = dFupdz/u[i][j];

		//Integrands of S.
	      
		intSx[i][j] = phix_x*Rx*Kxx + phix_x*Ry*Kxy + phix_x*Rz*Kxz
		            + phix_y*Rx*Kxy + phix_y*Ry*Kyy + phix_y*Rz*Kyz
		            + phix_z*Rx*Kxz + phix_z*Ry*Kyz + phix_z*Rz*Kzz;
		intSy[i][j] = phiy_x*Rx*Kxx + phiy_x*Ry*Kxy + phiy_x*Rz*Kxz
		            + phiy_y*Rx*Kxy + phiy_y*Ry*Kyy + phiy_y*Rz*Kyz
  			    + phiy_z*Rx*Kxz + phiy_z*Ry*Kyz + phiy_z*Rz*Kzz;
		intSz[i][j] = phiz_x*Rx*Kxx + phiz_x*Ry*Kxy + phiz_x*Rz*Kxz
                            + phiz_y*Rx*Kxy + phiz_y*Ry*Kyy + phiz_y*Rz*Kyz
  			    + phiz_z*Rx*Kxz + phiz_z*Ry*Kyz + phiz_z*Rz*Kzz;

		intSx[i][j] = intSx[i][j]*sqrt(deth[i][j]);
		intSy[i][j] = intSy[i][j]*sqrt(deth[i][j]);
		intSz[i][j] = intSz[i][j]*sqrt(deth[i][j]);

		localp[i1] = rank;

	      }else{
		localp[i1] = -1;
	      }

	    }
	  } //end of loop over the sphere

	  timer_stop(0, "AHF_local");

	  ////////////////////////////////////////////////////////////////////

	  //determine the processor with the highest rank for each point
	  bampi_allreduce_max_vector(localp, globalp, i1+1); 

	  //Integrals over the surface

	  //Integrate the area element, the mean of the
	  //expansion and the spin components

	  for(i=0;i<ntheta;i++){
	    int1_area[i]=0.0;
      	    int1_coarea[i]=0.0;
	    int1_hrms[i]=0.0;
	    int1_hmean[i]=0.0;
	    integral_Sx[i]=0.0;
	    integral_Sy[i]=0.0;
	    integral_Sz[i]=0.0;
	  }

	  cont=1;
	  for(i=0;i<ntheta;i++){ 
	    for(j=0;j<nphi;j++){
	      i1=i*nphi+j; 
	      if (globalp[i1] == rank){
		  theta = dtheta*(1.0/2.0 + i);
		  int1_area[i]+=dphi*sqrt(deth[i][j]);
		  // coordinate area, no metric involved:
     		  int1_coarea[i]+=rr[i][j]*rr[i][j]*sin(theta)*dphi;
		  int1_hrms[i]+=dphi*H[i][j]*H[i][j]*sqrt(deth[i][j]);
		  int1_hmean[i]+=dphi*H[i][j]*sqrt(deth[i][j]);
		  integral_Sx[i]+=dphi*intSx[i][j];
		  integral_Sy[i]+=dphi*intSy[i][j];
		  integral_Sz[i]+=dphi*intSz[i][j];
	      }else if (globalp[i1] == -1){
		  cont=0;
	      }	    
	    } 
	  }

	  if(!cont){
	      printf("AHF failed!\n");
	      printf("point not found in the processors!\n");
	      break;
	  }

	  area_l = 0.0;
    	  coarea_l = 0.0;
	  hrms_l = 0.0;
	  hmean_l = 0.0;
	  Sx_l = 0.0;
	  Sy_l = 0.0;
	  Sz_l = 0.0;
	  
	  for(i=0;i<ntheta;i++){
	    area_l+=dtheta*int1_area[i];
      	    coarea_l+=dtheta*int1_coarea[i];
	    hrms_l+=dtheta*int1_hrms[i];
	    hmean_l+=dtheta*int1_hmean[i];
	    Sx_l+=dtheta*integral_Sx[i];
	    Sy_l+=dtheta*integral_Sy[i];
	    Sz_l+=dtheta*integral_Sz[i];
	  }

	  //the integral is the sum over all local sums

	  v_l[0] = area_l; 
	  v_l[1] = hmean_l; 
	  v_l[2] = hrms_l; 
	  v_l[3] = Sx_l; 
	  v_l[4] = Sy_l; 
	  v_l[5] = Sz_l;
          v_l[6] = coarea_l;
	    
	  bampi_allreduce_sum_vector(v_l, v, 7);

	  area  = v[0];
	  hmean = v[1];
	  hrms  = v[2];
	  Sx    = v[3];
	  Sy    = v[4];
	  Sz    = v[5];
          coarea = v[6];
	  mass_old = mass;
	  mass = sqrt(area/(16.0*pi));
	  hrms/= area;

	  Sx = Sx/(8.0*acos(-1.0));
	  Sy = Sy/(8.0*acos(-1.0));
	  Sz = Sz/(8.0*acos(-1.0));
	  S = sqrt(Sx*Sx + Sy*Sy + Sz*Sz);

	  if (pr) {
	    printf("\n");
	    printf("mass = %22.16e\n",mass);
	    printf("radius = %22.16e\n",a0[0]/sqrt(4.0*pi));
	    printf("hrms = %22.16e\n",hrms);
	    printf("hmean = %22.16e\n",hmean);
	    printf("\n");
	  }

	  //exit(0);

          if(fabs(mass_old-mass)<mass_tol){
            SetAHF(number, time,xc,yc,zc, a0[0]/sqrt(4.0*pi),mass);
            if(rank==0 || rank==1){ //also write into stdout.1 
	      printf("\n");
	      printf("Horizon found\n");
	      printf("\n");
	      printf("mass = %22.16e\n",mass);
	      printf("radius = %22.16e\n",a0[0]/sqrt(4.0*pi));
	      printf("hrms = %22.16e\n",hrms);
	      printf("hmean = %22.16e\n",hmean);
	      printf("Sx = %22.16e\n",Sx);
	      printf("Sy = %22.16e\n",Sy);
	      printf("Sz = %22.16e\n",Sz);
	      printf("S = %22.16e\n",S);
	      printf("\n");

	     //5/2010: Edit for rank=0 Only (Original: rank=0 || rank =1)
                if (rank==0) {
		sprintf(ext, "_%d", number);
		sprintf(name,"%s/horizon%s",outdir, ext);

/*              // first version 4/07, remove when new version has been tested
#if 0
		if (ncall==0) {
		  fp1 = fopen(name, "wb");
		  fprintf(fp1, "%s\n", "#    time              x               y               z              mass             Sx              Sy              Sz              S");
		} else {
		  fp1 = fopen(name, "ab");
		}
		fprintf(fp1, "%14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",time,xc,yc,zc,mass,Sx,Sy,Sz,S);
		fclose(fp1);
		if((ncall==0)&&(number==nhorizons-1)) ncall=1;
#endif
*/
		// check whether file exists, use standard read for portability
		fp1 = fopen(name, "r");

		// if it does not exist, start new file and write header line
		if (!fp1) {
		  fp1 = fopen(name, "wb");
		  if (!fp1) errorexits("failed opening %s", name);
                  fprintf(fp1, "%s\n", "#    time              x               y               z              mass             Sx              Sy              Sz              S              coord. area");
		}
		
		// if it does exist, reopen for append
		else {
		  fclose(fp1);
		  fp1 = fopen(name, "ab");
		  if (!fp1) errorexits("failed opening %s", name);
		}

               
		// write data
		fprintf(fp1, "%14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e %14.6e  %14.6e %14.6e %14.6e\n",time,xc,yc,zc,mass,Sx,Sy,Sz,S,coarea);
		fclose(fp1);
                
                
                if (Getv("ahf_output", "yes")) AHF_output(ntheta, dtheta, nphi, dphi, xc, yc, zc, &rr, &outdir, number, time, ahf_time);		
		if (Getv("ahf_output_xyt", "yes")) AHF_output_xyt(ntheta, dtheta, nphi, dphi, xc, yc, zc, &rr, &outdir, number, time, ahf_time);	
	      } 
	    }
	    if(number==0) search0 = 0;
	    if(number==1) search1 = 0;
	    if(number==2) search2 = 0;

	    if(number==1){
	      if(nhorizons==2 &&common==1) find_common=1;
	    }
	    if(number==2) find_common=1;
	    break;
	  }

	  if(fabs(hmean)>hmean_tol){
	    if(rank==0){
	      printf("AHF failed!\n");
	      printf("hmean > hmean_tol\n");
	      printf("\n");
	    }
	    break;
	  }	 

	  if(k==flow_iter-1){
	    if(rank==0){
	      printf("AHF failed!\n");
	      printf("k = flow_iter\n");
	      printf("\n");
	    }
	  }
	  ///////////////////////////////////////////////////////////////////

	  //Find spectral components
	
	  for(l=0;l<=LMAX1;l++){
	    
	    for(i=0;i<ntheta;i++){
	      int1_spec0[i][l]=0.0;
	    }

	    //sum all the results that this processor is responsible for
	    //and check at the same time whether any points were left over

	    cont=1;
	    for(i=0;i<ntheta;i++){
	      for(j=0;j<nphi;j++){ 
		i1=i*nphi+j;
		if (globalp[i1] == rank){
		  int1_spec0[i][l]+=dphi*Y0[i1][l]*rho[i][j];
		}else if (globalp[i1] == -1){
		    cont=0;
		}	   
	      } 
	    }

	    spec0_l[l] = 0.0;
	  
	    for(i=0;i<ntheta;i++){
	      theta = dtheta*(1.0/2.0 + i);
	      spec0_l[l]+=dtheta*int1_spec0[i][l]*sin(theta);
	    }

	    //the integral is the sum over all local sums
	    bampi_allreduce_sum_vector(spec0_l+l, spec0+l, 1);

	    a0[l] -= A/(1.0+B*l*(l+1))*spec0[l];
	      
	    //printf("a0[%d] = %22.16e\n",l,a0[l]);

	    for(m=0;m<=l;m++){
	    
	      l1=l*(LMAX1+1)+m;
	      
	      for(i=0;i<ntheta;i++){
		int1_spec_c[i][l1]=0.0;
		int1_spec_s[i][l1]=0.0;
	      }

	      cont=1;
	      for(i=0;i<ntheta;i++){ 
		for(j=0;j<nphi;j++){ 
		  i1=i*nphi+j;
		  if (globalp[i1] == rank){
		    int1_spec_c[i][l1]+=dphi*Yc[i1][l1]*rho[i][j];
		    int1_spec_s[i][l1]+=dphi*Ys[i1][l1]*rho[i][j];
		  }else if (globalp[i1] == -1){
		      cont=0;
		  }
		} 
	      }

	      spec_cl[l1] = 0.0;
	      spec_sl[l1] = 0.0;
	  
	      for(i=0;i<ntheta;i++){
                theta = dtheta*(1.0/2.0 + i);
		spec_cl[l1]+=dtheta*int1_spec_c[i][l1]*sin(theta);
		spec_sl[l1]+=dtheta*int1_spec_s[i][l1]*sin(theta);
	      }

	      bampi_allreduce_sum_vector(spec_cl+l1, spec_c+l1, 1);
	      bampi_allreduce_sum_vector(spec_sl+l1, spec_s+l1, 1);

	      ac[l][m] -= A/(1.0+B*l*(l+1))*spec_c[l1];
	      as[l][m] -= A/(1.0+B*l*(l+1))*spec_s[l1];
	    
	      //printf("ac[%d][%d] = %22.16e\n",l,m,ac[l][m]);
	      //printf("as[%d][%d] = %22.16e\n",l,m,as[l][m]);
	    }
	  }

	  if(!cont){
	      printf("AHF failed!\n");
	      printf("point not found in the processors!\n");
	      break;
	  }

	}//end flow's loop 

      }//end if search==1 

    }//end loop over number of horizons

    /* free variable lists for 3d data */
    disablevarlist(vl);
    vlfree(vl);
    vlfree(wl);
    free(vinterp);
    
    // I commented this to get everything compiling, I do not know the influence! 
    //interpolate_setorder(Geti("order_RP"));

    /* free */
    free_dvector(a0,0,LMAX1);
    free_dmatrix(ac,0,LMAX1,0,LMAX1);
    free_dmatrix(as,0,LMAX1,0,LMAX1);
    free_dmatrix(rr,0,ntheta-1,0,nphi-1);
    free_dmatrix(u,0,ntheta-1,0,nphi-1);
    free_dmatrix(H,0,ntheta-1,0,nphi-1);
    free_dmatrix(rho,0,ntheta-1,0,nphi-1);
    free_dmatrix(P,0,LMAX1,0,LMAX1);
    free_dmatrix(dPdtheta,0,LMAX1,0,LMAX1);
    free_dmatrix(dPdthetadtheta,0,LMAX1,0,LMAX1);
    free_dmatrix(Y0,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(Yc,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(Ys,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dY0dtheta,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(dYcdtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dY0dthetadtheta,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(dYcdthetadtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdthetadphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdphidphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdthetadtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdthetadphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdphidphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dFdx,0,ntheta-1,0,nphi-1);
    free_dmatrix(dFdy,0,ntheta-1,0,nphi-1);
    free_dmatrix(dFdz,0,ntheta-1,0,nphi-1);
    free_dmatrix(int1_spec0,0,ntheta,0,LMAX1);
    free_dvector(spec0,0,LMAX1);
    free_dvector(spec0_l,0,LMAX1);
    free_dmatrix(int1_spec_c,0,ntheta,0,(LMAX1+1)*(LMAX1+1));
    free_dvector(spec_c,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(int1_spec_s,0,ntheta,0,(LMAX1+1)*(LMAX1+1));
    free_dvector(spec_s,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(deth,0,ntheta-1,0,nphi-1);
    free_dvector(int1_area,0,ntheta-1);
    free_dvector(int1_coarea,0,ntheta-1);
    free_dvector(localp,0,ntheta*nphi);
    free_dvector(globalp,0,ntheta*nphi);
    free_dvector(dxdt,0,ntheta*nphi);
    free_dvector(dxdp,0,ntheta*nphi);
    free_dvector(dydt,0,ntheta*nphi);
    free_dvector(dydp,0,ntheta*nphi);
    free_dvector(dzdt,0,ntheta*nphi);
    free_dvector(dzdp,0,ntheta*nphi);
    free_dmatrix(intSx,0,ntheta-1,0,nphi-1);
    free_dmatrix(intSy,0,ntheta-1,0,nphi-1);
    free_dmatrix(intSz,0,ntheta-1,0,nphi-1);
    free_dvector(integral_Sx,0,ntheta-1);
    free_dvector(integral_Sy,0,ntheta-1);
    free_dvector(integral_Sz,0,ntheta-1);
    free_dmatrix(sigma,0,ntheta-1,0,nphi-1);
    free_dvector(spec_sl,0,(LMAX1+1)*(LMAX1+1));
    free_dvector(spec_cl,0,(LMAX1+1)*(LMAX1+1));
    
    timer_stop(0, "AHF");

    return 0;

}












void SetAHF(int N, double t,double x,double y,double z, double r, double mass)
{
    if (N==0) {
        Setd("ahf0_x",x);
        Setd("ahf0_y",y);
        Setd("ahf0_z",z);
        Setd("ahf0_m",mass);
        Setd("ahf0_r",r);
    } else if (N==1) {
        Setd("ahf1_x",x);
        Setd("ahf1_y",y);
        Setd("ahf1_z",z);
        Setd("ahf1_m",mass);
        Setd("ahf1_r",r);
    } else if (N==2) {
        Setd("ahf2_x",x);
        Setd("ahf2_y",y);
        Setd("ahf2_z",z);
        Setd("ahf2_m",mass);
        Setd("ahf2_r",r);
    } else
        printf("only 3 horizons are implemented");
}



















