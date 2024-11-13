#include <string.h>
#include "syspara.h"

void input_para(FILE *fpin)
{
	int i,ii;

	fscanf(fpin,"%d",&simtype);
	fscanf(fpin,"%d",&cluster);
	fscanf(fpin,"%d",&Nx1);
	fscanf(fpin,"%d",&Nx2);
	fscanf(fpin,"%d",&Nx3);
	fscanf(fpin,"%d",&Nx4);
	fscanf(fpin,"%d",&Ny1);
	fscanf(fpin,"%d",&Ny2);
	fscanf(fpin,"%d",&Ny3);
	fscanf(fpin,"%d",&Ny4);
	fscanf(fpin,"%lf",&Diff_l);
	fscanf(fpin,"%lf",&Diff_t);
	fscanf(fpin,"%lf",&rateD);
	fscanf(fpin,"%lf",&BCL);
	fscanf(fpin,"%lf",&kr_rate0);
	fscanf(fpin,"%lf",&kr_rate1);
	fscanf(fpin,"%lf",&ks_rate);
	fscanf(fpin,"%lf",&rategca);
	fscanf(fpin,"%lf",&ratef);
	fscanf(fpin,"%lf",&Istim_base);
	fscanf(fpin,"%d",&stim_sw);
	fscanf(fpin,"%d",&datas);
	for (ii = 0; ii < datas; ii++){
		for (i=0;i<NN;i++){
			fscanf(fpin,"%lf",&x0[ii][i]);
		}
	}
	fscanf(fpin,"%d",&div_num);
	fscanf(fpin,"%d",&write);
	fscanf(fpin,"%d",&reduce);
	fscanf(fpin,"%d",&outdata);
	fscanf(fpin,"%d",&outdata_plus);
	fscanf(fpin,"%d",&debug);

}

void input_init(FILE *fpin1, FILE *fpin2)
{

	int i,j;
	int check1=0,check2=0;
	double init1[NN],init2[NN];

	printf("input init, middle point=%d\n",MEDIA_SITE/2);

	for (i=0;i<NN;i++){
		fscanf(fpin1,"%lf",&init1[i]);
	}

	for (i=0;i<NN;i++){
		fscanf(fpin2,"%lf",&init2[i]);
	}

	for (i=0;i<MEDIA_SITE;i++){
		if(utype[i]==1){
			if(check1==0) printf("EAD unit initial values\n");
			for (j=0;j<NN;j++){
				//var.u[j][i]=init2[j];
				vmf[j][i]=init2[j];
				//if(check1==0) printf("x[%d][%d] = %lf\n",j,i,var.u[j][i]);
				if(check1==0) printf("x[%d][%d] = %lf\n",j,i,vmf[j][i]);
			} check1++;
			printf("%d, ",i);
		} else if(utype[i]==0){
			if(check2==0) printf("noEAD unit initial values\n");
			for (j=0;j<NN;j++){
				//var.u[j][i]=init1[j];
				vmf[j][i]=init1[j];
				//if(check2==0) printf("x[%d][%d] = %lf\n",j,i,var.u[j][i]);
				if(check2==0) printf("x[%d][%d] = %lf\n",j,i,vmf[j][i]);
			} check2++;
		}
	}
	printf("\n");

}

void val_consts(FILE *fp1)
{
	int i,w,c;
	double v_old,dvdt,dvdt_new;

	// Cell Geometry */
		radius = 0.0011;
		length = 0.01;
		//del_x = length;
		del_x = 2.0*length;
		//CapHV = 153.4;  // Membrane Capacitance (pF) : K05 model
		CapHV = 175.0;  // Membrane Capacitance (pF)
		//vcell = 38.0*0.000000000001; // Cell Volume:3.801e-5 (uL)
		//vi = vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		//vup = vcell*0.0552;    // NSR volume (uL)
		//vrel = vcell*0.0048;    // JSR volume (uL)
		vmyo = 16.404; 		// Cytoplasmic volume (pico L)
		vsr  =  1.094;		// Sarcoplasmic volume (pico L)
		vss  =  0.05468;	// Subspace volume (pico L)

		vr1 = CapHV/F/(vmyo + vss);
		vr2 = CapHV/zca/F/vmyo;
		vr6 = CapHV/zca/F/vss;
		vr3 = vsr/vmyo;
		vr4 = vsr/vss;
		vr5 = vmyo/vss;

	// invariant constant
		RTonF = R*T/F;
		RTon2F = R*T/zca/F;
		RTon4F = R*T/zca/zca/F/F;

	// Extracellular Concentrations
		Nao = 140.0;    // (mM) 
		Ko = 5.4;      // (mM)
		Cao = 2.0;     // (mM)

		//prnak = 0.01833;     // Na/K Permiability Ratio
		prnak = 0.03;     // Na/K Permiability Ratio
		pkna = 0.12;     // Na/K Permiability Ratio

	// cell types

		if(simtype==0){ // only Endo cell unit myofiber
			for(i=0;i<MEDIA_SITE;i++) celltype[i]=0;	
		} else if(simtype==1){ // only Epi cell unit myofiber
			for(i=0;i<MEDIA_SITE;i++) celltype[i]=1;	
		} else if(simtype==2){ // only Mid cell unit myofiber
			for(i=0;i<MEDIA_SITE;i++) celltype[i]=2;	
		} else if(simtype==5){ // Transmural case (#100 units)
			for(i=0;i<15;i++) celltype[i]=0;	
			for(i=15;i<70;i++) celltype[i]=2;	
			for(i=70;i<100;i++) celltype[i]=1;	
		}
	
	// parameter setup information
	printf("Nx1=%d,Ny1=%d,Nx2=%d,Ny2=%d,Nx3=%d,Ny3=%d,Nx4=%d,Ny4=%d\n",Nx1,Ny1,Nx2,Ny2,Nx3,Ny3,Nx4,Ny4);

	for(i=0;i<MEDIA_SITE;i++) utype[i]=0;

	switch(cluster){
		case 2: // 2 cluster (20path)
			for(i=0 ; i < MEDIA_SITE; i++){
				if( ((i >= Ny1*DX && i < Ny2*DX) && (i%DX >= Nx1 && i%DX < Nx2 )) || 
				  	((i >= Ny3*DX && i < Ny4*DX) && (i%DX >= Nx3 && i%DX < Nx4 )) 
				) {utype[i]=1;	c++;}
			}
		break;
	}

	// Fast sodium current
		for(i=0;i<MEDIA_SITE;i++){
			//Gnaf[i] = 7.8;	// (mS/uF). : K05 model
			Gnaf[i] = 14.838;	// (mS/uF). : TP06 original
		}
		
	// Transient outward current
		//pnato = 0.043;
		
		//dVgq0= -12.0; //gating variables (q/r) % Endo (from Bernus 2002 model)
		//dVgq1 = 0.0; //gating variables (q/r) % Epi (default for Kurata 2005 model)
		//dVgq2 = -4.0; //gating variables (q/r) % Mid (from Bernus 2002 model)

		for(i=0;i<MEDIA_SITE;i++){
			if(celltype[i]==0){
				//Gto[i] = 0.13;	// Endo (from Bernus 2002 model).
				Gto[i] = 0.073;		// Endo (nS/pF).
			} else if(celltype[i]==1){
				//Gto[i] = 0.4;		// Epi (default for Kurata 2005 model).
				Gto[i] = 0.294;		// Epi (nS/pF).
			} else if(celltype[i]==2){
				//Gto[i] = 0.35;	// Mid (from Bernus 2002 model).
				Gto[i] = 0.294;		// Mid (nS/pF).
			}
		}

	// Rapid delayed rectifier potassium current (Ikr)
		
		for(i=0;i<MEDIA_SITE;i++){
			Gkr[i] = 1.5*0.153*sqrt(Ko/5.4);	// (nS/pF)
			if(utype[i]==1){ // EAD case
				gkr_rate[i]=kr_rate1;
			} else {		// EAD case
				gkr_rate[i]=kr_rate0;
			}
		}
		printf("gkr_rate[%d]=%lf\n",0,gkr_rate[0]);

	// Slow delayed rectifier potassium current (Iks)
		for(i=0;i<MEDIA_SITE;i++){
			if(celltype[i] == 0){
				Gks[i] = 0.4*0.392;	// Endo (nS/pF)
				gks_rate[i]=ks_rate;
			} else if(celltype[i] == 1){
				Gks[i] = 0.4*0.392;	// Epi (nS/pF)
				gks_rate[i]=ks_rate;
			} else if(celltype[i] == 2){
				Gks[i] = 0.4*0.098;	// Mid (nS/pF)
				gks_rate[i]=ks_rate;
			}
		}
		printf("gks_rate[%d]=%lf\n",0,gks_rate[0]);
	
	// L-type calcium current
		//kmfca = 0.00035;
		//Ecal = 52.8; // (mV)
		for(i=0;i<MEDIA_SITE;i++){
			Gca[i] = 2.0*3.980E-5;
		}

	// Inward rectifier K current: Ik1
		for(i=0;i<MEDIA_SITE;i++){
			Gk1[i] = 5.405;	// (nS/pF)
		}
		rategk1 = sqrt(Ko/5.4);

	// Sodium-Calcium Exchanger V-S
		kmnaex = 87.5;
		kmcaex = 1.38;
		rncx = 0.35;
		ksat = 0.1;

		for(i=0;i<MEDIA_SITE;i++){
			if(celltype[i] == 0){
				kncx[i] = 1000.0*(1.0/1.1);	// Endo (based on O'Hara-Rudy 2011 model);
			} else {
				kncx[i] = 1000.0;			// Epi (default for Kurata 2005 model) or Mid;
			}
		}
		c1 = 1.0/(pow(kmnaex,3)+pow(Nao,3))/(kmcaex + Cao);

	// Sodium-Potassium Pump
		//kmnap = 10.0;
		//kmkp = 1.5;
		//nna = 1.5;
		kmna = 40.0;
		kmk  = 1.0;

		for(i=0;i<MEDIA_SITE;i++){
			Gnak[i] = 2.724;	// (pA/pF)
		}
		//rhonak = (exp(Nao/67.3)-1.0)/7.0;

	// Sarcolemmal Ca Pump
		Gpca = 0.1238;		// Max. Ca current through sarcolemmal Ca pump (mS/uF)
		Km_pca = 0.0005;	// Half-saturation concentration of sarcolemmal Ca pump (mM)

	// plateau K current 
		Gkp = 0.0146;	// (nS/pF)

	// Ca Background Current 
		Gcab = 0.000592;	// (nS/pF)

	// Na Background Current 
		Gnab = 0.000290;	// (nS/pF)

	// calcium uptake via SERCA pump (Jup)
		kup = 0.00025;
		for(i=0;i<MEDIA_SITE;i++){
			Gup[i] = 0.6*0.006375;	// maximal uptake conductance (mM)
		}

	// SR calcium release flux, via RyR (Jrel)
	/*	nrel = 3.0;
		for(i=0;i<MEDIA_SITE;i++){
			if(celltype[i]!=2){
				Grel[i] = 50.0;	// Epi (default for Kurata 2005 model) or Endo
			} else {
				Grel[i] = 50.0*1.7;	// Mid (based on O'Hara-Rudy 2011 model)
			}
		}
		tau_dr = 1.0*4.0;
		tau_df = 4.0*4.0;
	*/
		for(i=0;i<MEDIA_SITE;i++){
			Grel[i] = 0.102;	// (mM/ms)
		}
		//vrel = 0.102; // (mM/ms)
		k1_prime = 0.15;   // R to O and RI to I Irel transition rate (1/mM^2/ms)
		k2_prime = 0.045;  // O to I and R to RI Irel transition rate (1/mM/ms)
		relk3 = 0.060;                // O to R and I to RI Irel transition rate (1/ms)
		relk4 = 0.005;                // I to O and RI to I Irel transition rate (1/ms)
		EC = 1.5;  // Ca_sr half-saturation constant of kcasr (mM)
		maxsr = 2.5;
		minsr = 1.0;

	// Translocation of Ca Ions from cytoplasm to subspace
		//tau_tr = 180.0;	// Time constant of Ca transfer from NSR to JSR (ms)
		vxfer = 0.0038;   	// Maximal Ixfer conductance (mM/ms)


	// Ca leakage from cytoplasm to subspace 
		//pleak = 0.00026;	// K05 model case
		pleak = 0.00036;	// Maximal Ileak conductance (mM/ms)

	//Troponin Ca Buffering (in Myoplasm)
	//	concTc = 0.07;
	//	rftn = 10.0;
	//	kftc = 40.0*rftn;
	//	kbtc = 0.02*rftn;
		
	//Calmodulin Ca Buffering (in Myoplasm) --- Rapid Buffering Approximation
	/*	for(i=0;i<MEDIA_SITE;i++){
			if(celltype[i] != 1){
				conccm[i] = 0.05*(1.0/1.3);  //% Mid or Endo (based on O'Hara-Rudy 2011 model)
			} else {
				conccm[i] = 0.05;   //% Epi (default for Kurata 2005 model)
			}
		}
		kdcm = 0.00238;	// % Kd for Ca-binding to Calmoduline 
	*/
	//Calsequestrin Ca Buffering (in SR) --- Rapid Buffering Approximation	
	/*
		conccq = 10.0;     // Max. [Ca] buffered in CSQN (mM)
		kdcq = 0.8;     // % Kd for Ca-binding to Calsequestrin
	*/
	// Myoplasmic Ca Ion Concentration Changes 
		bufc = 0.20;	// Max. [Ca] buffer concentration (mM)
		kcai = 0.001;	// Cai half-saturation constant for cytoplasmic buffer (mM)

		bufsr = 10.0;	// Total sarcoplasmic buffer concentration (mM)
		kcasr = 0.30;	// Ca_sr half-saturation constant for sarcoplasmic buffer (mM)

		bufss = 0.40000;	// Total subspace buffer concentration (mM)
		kcass = 0.00025;	// Ca_ss half-saturation constant for subspace buffer (mM)

}

void utype_data(FILE *fp17)
{

	int i,j;

	fprintf(fp17,"%lf ",0.0);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp17,"%d ",utype[i]);
	}
	fprintf(fp17,"\n");

} 

void vm_data(FILE *fp5, double time)
{

	int i,j;

	fprintf(fp5,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		//fprintf(fp5,"%lf ",var.u[0][i]);
		fprintf(fp5,"%lf ",vmf[0][i]);
	}
	fprintf(fp5,"\n");

} 

void out_data(FILE *fp7, FILE *fp8, FILE *fp9, FILE *fp10, FILE *fp11, FILE *fp12, double time)
{
	int i;

// ikr
	fprintf(fp7,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp7,"%lf ",ikr[i]);
	}
	fprintf(fp7,"\n");

// iks
	fprintf(fp8,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp8,"%lf ",iks[i]);
	}
	fprintf(fp8,"\n");

// ical
	fprintf(fp9,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp9,"%lf ",ical[i]);
	}
	fprintf(fp9,"\n");

// incx
	fprintf(fp10,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp10,"%lf ",jncx[i]);
	}
	fprintf(fp10,"\n");

// inak
	fprintf(fp11,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp11,"%lf ",inak[i]);
	}
	fprintf(fp11,"\n");

// Cai
	fprintf(fp12,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		//fprintf(fp12,"%lf ",var.u[12][i]);
		fprintf(fp12,"%lf ",vmf[14][i]);
	}
	fprintf(fp12,"\n");

}

//void out_data_plus(FILE *fp13, FILE *fp14, FILE *fp15, FILE *fp16, double time)
void out_data_plus(FILE *fp13, FILE *fp14, double time)
{
	int i;

// detailed Vm data
	fprintf(fp13,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		//fprintf(fp13,"%lf ",var.u[0][i]);
		fprintf(fp13,"%lf ",vmf[0][i]);
	}
	fprintf(fp13,"\n");
// ito
	fprintf(fp14,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp14,"%lf ",jrel[i]);
	}
	fprintf(fp14,"\n");

/*	
// ik1
	fprintf(fp15,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp15,"%lf ",ik1.ik[i]);
	}
	fprintf(fp15,"\n");

// ipca
	fprintf(fp16,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		fprintf(fp16,"%lf ",ipca.ca[i]);
	}
	fprintf(fp16,"\n");
*/ 
}

void ECG_data(FILE *fp18, double time, double ECG)
{

	fprintf(fp18,"%lf %lf\n",time,ECG);

}
