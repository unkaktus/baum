// compute Dvs: symmetrized covariant derivative of vlow *STILL JUST PARTIAL DERIVATIVE*
void grhd_compute_Dvs(double* vx, double* vy, double* vz, 
		 double* Dvs11, double* Dvs12, double* Dvs13, double* Dvs22, double* Dvs23, double* Dvs33, 
		 tL* level){
	const int order_centered = Geti("order_centered");
	
	bampi_openmp_start
	
	double dx = level->dx;
	double dy = level->dy;
	double dz = level->dz;
	double oodx = 1/dx, oo2dy = 1/dy, oo2dz = 1/dz;
	
	forinnerpoints_ijk_openmp(level) {
	
	if (order_centered == 2 || boundaryNaway(1)) { 
		Dvs11[ijk] =  oodx*(-vx[-di + ijk] + vx[di + ijk]);						//11
		Dvs12[ijk] =  oodx/2*(-vy[-di + ijk] + vy[di + ijk]) + oody/2*(-vx[-dj + ijk] + vx[dj + ijk]);	//12
		Dvs13[ijk] =  oodx/2*(-vz[-di + ijk] + vz[di + ijk]) + oodz/2*(-vx[-dk + ijk] + vx[dk + ijk]);	//13
		Dvs22[ijk] =  oo2dy*(-vy[-dj + ijk] + vy[dj + ijk]);						//22
		Dvs23[ijk] =  oody/2*(-vz[-dj + ijk] + vz[dj + ijk]) + oodz/2*(-vy[-dk + ijk] + vy[dk + ijk]);	//23
		Dvs33[ijk] =  oo2dz*(-vz[-dk + ijk] + vz[dk + ijk]);						//33
		
	} else if (order_centered == 4 || boundaryNaway(2)) { 
	
		Dvs11[ijk] =  oodx/6*(vx[-2*di + ijk] - 8*vx[-di + ijk] + 8*vx[di + ijk] - vx[2*di + ijk]);	//11
		
		Dvs12[ijk] =  oodx/12*(vy[-2*di + ijk] - 8*vy[-di + ijk] + 8*vy[di + ijk] - vy[2*di + ijk]) +
			       oody/12*(vx[-2*dj + ijk] - 8*vx[-dj + ijk] + 8*vx[dj + ijk] - vx[2*dj + ijk]);	//12
			       
		Dvs13[ijk] =  oodx/12*(vz[-2*di + ijk] - 8*vz[-di + ijk] + 8*vz[di + ijk] - vz[2*di + ijk]) + 
			       oodz/12*(vx[-2*dk + ijk] - 8*vx[-dk + ijk] + 8*vx[dk + ijk] - vx[2*dk + ijk]);	//13
			       
		Dvs22[ijk] =  oody/6*(vy[-2*dj + ijk] - 8*vy[-dj + ijk] + 8*vy[dj + ijk] - vy[2*dj + ijk]);	//22
		
		Dvs23[ijk] =  oody/12*(vz[-2*dj + ijk] - 8*vz[-dj + ijk] + 8*vz[dj + ijk] - vz[2*dj + ijk]) +
			       oodz/12*(vy[-2*dk + ijk] - 8*vy[-dk + ijk] + 8*vy[dk + ijk] - vy[2*dk + ijk]);	//23
		
		Dvs33[ijk] =  oo2dz/6*(vz[-2*dk + ijk] - 8*vz[-dk + ijk] + 8*vz[dk + ijk] - vz[2*dk + ijk]);	//33
		
	} else if (order_centered == 6 || boundaryNaway(3)) { 
		Dvs11[ijk] =  oodx/15*(-vx[-3*di + ijk] + 9*vx[-2*di + ijk] - 45*vx[-di + ijk] + 45*vx[di + ijk] - 9*vx[2*di + ijk] + vx[3*di + ijk]);	//11
		
		Dvs12[ijk] =  oodx/30*(-vy[-3*di + ijk] + 9*vy[-2*di + ijk] - 45*vy[-di + ijk] + 45*vy[di + ijk] - 9*vy[2*di + ijk] + vy[3*di + ijk]) +
			       oody/30*(-vx[-3*dj + ijk] + 9*vx[-2*dj + ijk] - 45*vx[-dj + ijk] + 45*vx[dj + ijk] - 9*vx[2*dj + ijk] + vx[3*dj + ijk]);	//12
			       
		Dvs13[ijk] =  oodx/30*(-vz[-3*di + ijk] + 9*vz[-2*di + ijk] - 45*vz[-di + ijk] + 45*vz[di + ijk] - 9*vz[2*di + ijk] + vz[3*di + ijk]) + 
			       oodz/30*(-vx[-3*dk + ijk] + 9*vx[-2*dk + ijk] - 45*vx[-dk + ijk] + 45*vx[dk + ijk] - 9*vx[2*dk + ijk] + vx[3*dk + ijk]);	//13
			       
		Dvs22[ijk] =  oody/15*(-vy[-3*dj + ijk] + 9*vy[-2*dj + ijk] - 45*vy[-dj + ijk] + 45*vy[dj + ijk] - 9*vy[2*dj + ijk] + vy[3*dj + ijk]);	//22
		
		Dvs23[ijk] =  oody/30*(-vz[-3*dj + ijk] + 9*vz[-2*dj + ijk] - 45*vz[-dj + ijk] + 45*vz[dj + ijk] - 9*vz[2*dj + ijk] + vz[3*dj + ijk]) +
			       oodz/30*(-vy[-3*dk + ijk] + 9*vy[-2*dk + ijk] - 45*vy[-dk + ijk] + 45*vy[dk + ijk] - 9*vy[2*dk + ijk] + vy[3*dk + ijk]);	//23
		
		Dvs33[ijk] =  oodz/15*(-vz[-3*dk + ijk] + 9*vz[-2*dk + ijk] - 45*vz[-dk + ijk] + 45*vz[dk + ijk] - 9*vz[2*dk + ijk] + vz[3*dk + ijk]);	//33
	}
	
	} endfor_ijk_openmp;

	bampi_openmp_stop
	}
