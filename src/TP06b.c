#include <string.h>
#include "syspara.h"

void eular(int n, double h, double t)
{
    int i,j;  
    double dev,diff_c;
    double vm;

	//#pragma vector aligned
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_Erev(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ina(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ik1(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ito(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ikr(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_iks(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ical(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_inak(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ikp(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_icab(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_inab(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_incx(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_ipca(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_jup(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_jleak(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_jrel(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_jxfer(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_itot(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_concentration(i);
	#pragma ivdep
	for(i=0;i<MEDIA_SITE;i++) comp_diff_current(i);
	//#pragma ivdep
	#pragma omp parallel for simd
	for(i=0;i<MEDIA_SITE;i++) function(i,t);
	//for(i=MEDIA_SITE-1;i>=0;i--) function(i,t);

	#pragma ivdep
	for(j=0;j<NN;j++){
		for(i=0;i<MEDIA_SITE;i++) vmf[j][i] += h*dvmf[j][i];
		//for(i=MEDIA_SITE-1;i>=0;i--) vmf[j][i] += h*dvmf[j][i];
	}
}

void function(int s, double t){

	dvmf[0][s] = -itotal[s] + imo[s]; // Vm
	dvmf[1][s] = (mss[s] - vmf[1][s])/taum[s]; // m
	dvmf[2][s] = (hss[s] - vmf[2][s])/tauh[s]; // h
	dvmf[3][s] = (jss[s] - vmf[3][s])/tauj[s]; // j
	dvmf[4][s] = (rss[s] - vmf[4][s])/taur[s]; // r
	dvmf[5][s] = (sss[s] - vmf[5][s])/taus[s]; // s 
	dvmf[6][s] = (xr1ss[s] - vmf[6][s])/tauxr1[s]; // xr1
	dvmf[7][s] = (xr2ss[s] - vmf[7][s])/tauxr2[s]; // xr2
	dvmf[8][s] = (xsss[s] - vmf[8][s])/tauxs[s]; // xs
	dvmf[9][s] = (dss[s] - vmf[9][s])/taud[s]; // d 
	dvmf[10][s] = (fss[s] - vmf[10][s])/tauf[s]; // f 
	dvmf[11][s] = (f2ss[s] - vmf[11][s])/tauf2[s]; // f2
	dvmf[12][s] = (fcass[s] - vmf[12][s])/taufca[s]; // fca
	dvmf[13][s] = -relk2*vmf[16][s]*vmf[13][s] + relk4*(1.0 - vmf[13][s]); // R

	dvmf[14][s] = bufcai[s]*(-jcainet[s]*vr2 + (jleak[s] - jup[s])*vr3 + jxfer[s]); // Cai
	dvmf[15][s] = bufcasr[s]*(jup[s] - jleak[s] - jrel[s]);	// Casr
	dvmf[16][s] = bufcass[s]*(-jcasnet[s]*vr6 + jrel[s]*vr4 - jxfer[s]*vr5); // Cass
	dvmf[17][s] = -jnanet[s]*vr1;	// Nai
	dvmf[18][s] = -jknet[s]*vr1;		// Ki

}

void comp_diff_current(int site)
{
    int i;  
	double temp,dx,rD,lD,tD;

	dx = del_x;
	rD = rateD;
	lD = Diff_l;
	tD = Diff_t;

	if(site==0){ // site=0
		temp = -rD/dx/dx*(tD*(vmf[0][site] - vmf[0][site+1]) + lD*(vmf[0][site] - vmf[0][site+DX]));
	} else if (site==MEDIA_SITE-1){ //site = 8
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - vmf[0][site]) + lD*(vmf[0][site-DX] - vmf[0][site]));
	} else if (site==DX-1){ // site == 2
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - vmf[0][site]) - lD*(vmf[0][site] - vmf[0][site+DX]));
	} else if (site==DX*(DY-1)){ // site == 6
		temp = -rD/dx/dx*(tD*(vmf[0][site] - vmf[0][site+1]) - lD*(vmf[0][site-DX] - vmf[0][site]));
	} else if (site>0 && site<DX-1){ // site == 1
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - 2.0*vmf[0][site] + vmf[0][site+1]) - lD*(vmf[0][site] - vmf[0][site+DX]));
	} else if (site>DX*(DY-1) && site<DX*DY-1){ // site == 7
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - 2.0*vmf[0][site] + vmf[0][site+1]) + lD*(vmf[0][site-DX] - vmf[0][site]));
	} else if (site > 0 && site%DX == 0 && site < DX*(DY-1)){ // site == 3
		temp = -rD/dx/dx*(tD*(vmf[0][site] - vmf[0][site+1]) - lD*(vmf[0][site-DX] - 2.0*vmf[0][site] + vmf[0][site+DX]));
	} else if (site>DX-1 && site%DX==DX-1 && site<DX*DY-1){ // site == 5
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - vmf[0][site]) + lD*(vmf[0][site-DX] - 2.0*vmf[0][site] + vmf[0][site+DX]));
	} else { // site == 4
		temp =  rD/dx/dx*(tD*(vmf[0][site-1] - 2.0*vmf[0][site] + vmf[0][site+1]) + lD*(vmf[0][site-DX] - 2.0*vmf[0][site] + vmf[0][site+DX]));
	}
	imo[site] = temp;

}

void comp_itot(int s){

	int i;

	// Total Membrane Current (pA/pF)
	if(stim_sw==0){
		if (s >= 0 && s < DX ){ // row stimulation
		//if(s%100==0) 
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s] + Istim;
		} else {
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s];
		}
	} else if(stim_sw==1){ // column stimulation
		if (s%DX==0){
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s] + Istim;
		} else {
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s];
		}
	} else if(stim_sw==2){	// corner stimulation
		if ((s >= 0 && s < (int)(0.05*DX))||(s%DX<(int)(0.05*DX) && s< (int)(0.05*DY*DX)) ){
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s] + Istim;
		} else {
			itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s];
		}
	} else if(stim_sw > 2){ // i.e., uncoupled system with all units stimulation
		itotal[s] = inaf[s] + ical[s] + ito[s] + ikr[s] + iks[s] + ik1[s]
						+ inab[s] + icab[s] + ikp[s] + inak[s] + jncx[s] + ipca[s] + Istim;
		Diff_l = 0.0; Diff_t = 0.0;
	}

	// Net Ion Fluxes (umol/cell)

	jnanet[s] = inaf[s] + inab[s] + 3.0*inak[s] + 3.0*jncx[s];
	jknet[s] = ito[s] + ikr[s] + iks[s] + ik1[s] + ikp[s] - 2.0*inak[s] + Istim;
	jcainet[s] = icab[s] - 2.0*jncx[s] + ipca[s];
	jcasnet[s] = ical[s];

}

void comp_ina(int site)
{
	double vm,m,h,j,gna,ena;
	int iV=0;
	double V1,V2,d1,d2;
	vm = vmf[0][site];
	m = vmf[1][site];
	h = vmf[2][site];
	j = vmf[3][site];
	gna = Gnaf[site];
	ena = Enaf[site];

	V1 = ( vm + Emax )*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	
	mss[site] = Tmss[iV]*d2 + Tmss[iV+1]*d1;
	taum[site] = Ttaum[iV]*d2 + Ttaum[iV+1]*d1;
	hss[site] = Thss[iV]*d2 + Thss[iV+1]*d1;
	tauh[site] = Ttauh[iV]*d2 + Ttauh[iV+1]*d1;
	jss[site] = Tjss[iV]*d2 + Tjss[iV+1]*d1;
	tauj[site] = Ttauj[iV]*d2 + Ttauj[iV+1]*d1;
	
	inaf[site] = gna*(vm - ena)*m*m*m*h*j;
	
	//printf("inaf[%ld]=%e\n",site,inaf[site]);
}

// Ito Transient Outward Current
void comp_ito(int site)
{
	double vm,r,s,Eko,gto;
	vm = vmf[0][site];
	r = vmf[4][site];
	s = vmf[5][site];
	Eko = Ek[site];
	gto = Gto[site];

	int iV=0;
	double V1,V2,d1,d2;
	V1 = (vm + Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	rss[site] = Trss[iV]*d2 + Trss[iV+1]*d1;
	taur[site] = Ttaur[iV]*d2 + Ttaur[iV+1]*d1;

	sss[site] = Tsss[iV]*d2 + Tsss[iV+1]*d1;
	taus[site] = Ttaus[iV]*d2 + Ttaus[iV+1]*d1;

	ito[site] = gto*(vm - Eko)*r*s;

//	printf("ito[%ld]=%e\n",site,ito[site]);
}

// Rapidly Activating Potassium Current 
void comp_ikr(int site)
{
	double vm,xr1,xr2,gkr,Ekr;
	double grate;

	int iV=0;	
	double V1,V2,d1,d2;
	vm = vmf[0][site];
	xr1 = vmf[6][site];
	xr2 = vmf[7][site];
	gkr = Gkr[site];
	Ekr = Ek[site];
	grate = gkr_rate[site];

	V1 = (vm+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	 xr1ss[site] = Txr1ss[iV]*d2 + Txr1ss[iV+1]*d1;
	tauxr1[site] = Ttauxr1[iV]*d2 + Ttauxr1[iV+1]*d1;

	 xr2ss[site] = Txr2ss[iV]*d2 + Txr2ss[iV+1]*d1;
	tauxr2[site] = Ttauxr2[iV]*d2 + Ttauxr2[iV+1]*d1;

	ikr[site] = grate*gkr*xr1*xr2*(vm - Ekr);

//	printf("ikr[%ld]=%e\n",site,ikr[site]);
}

// Slowly Activating Potassium Current 
void comp_iks(int site)
{
	double vm,xs,gks,eks;
	double grate;

	int iV=0;
	double V1,V2,d1,d2;
	vm = vmf[0][site];
	xs = vmf[8][site];
	gks = Gks[site];
	eks = Eks[site];
	grate = gks_rate[site];

	V1 = (vm+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	 xsss[site] = Txsss[iV]*d2 + Txsss[iV+1]*d1;
	tauxs[site] = Ttauxs[iV]*d2 + Ttauxs[iV+1]*d1;
	
	iks[site] = grate*gks*(vm - eks)*xs*xs;

//	printf("iks[%ld]=%e\n",site,iks[site]);

}

// L-type calcium current
void comp_ical(int site)
{

	double vm,d,f,f2,fca,rCal,Eca;
	double gca,Cai,Cass,km;
	double tmp1,tmp2;

	int iV=0;
	double V1,V2,d1,d2;
	vm    = vmf[0][site];
	d     = vmf[9][site];
	f     = vmf[10][site];
	f2    = vmf[11][site];
	fca   = vmf[12][site];
	Cai   = vmf[14][site];
	Cass   = vmf[16][site];
	//fca = gfCa_inf;
	Eca = Ecal;
	gca = Gca[site];
	//km = kmfca;
	//fca = km/(km + Cai);
	
	V1 = ( vm + Emax )*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// VDA
	dss[site] = Tdss[iV]*d2 + Tdss[iV+1]*d1;
	taud[site] = Ttaud[iV]*d2 + Ttaud[iV+1]*d1;
	// VDI 
	fss[site] = Tfss[iV]*d2 + Tfss[iV+1]*d1;
	tauf[site] = Ttauf[iV]*d2 + Ttauf[iV+1]*d1;
	// VDI2 
	f2ss[site] = Tf2ss[iV]*d2 + Tf2ss[iV+1]*d1;
	tauf2[site] = Ttauf2[iV]*d2 + Ttauf2[iV+1]*d1;
	// CDI 
	fcass[site] = 0.4+0.6/(1.0+(Cass/0.05)*(Cass/0.05));
	taufca[site] = 2.0+80.0/(1.0+(Cass/0.05)*(Cass/0.05));
	//fca = km/(km + Cai);
	//gfCa_inf = fca;
	rCal = Trcal[iV]*d2 + Trcal[iV+1]*d1;

  //ical[site] =rategca*gca*d*f*f2*fca*(vm - 15.0)/RTon4F*(0.25*Cass*exp((vm-15.0)/RTon2F)-Cao)/(exp((vm-15.0)/RTon2F)-1.0);
	ical[site] =rategca*gca*d*f*f2*fca*(vm - 15.0)/RTon4F*(0.25*Cass*rCal - Cao)/(rCal-1.0);

//	printf("ical[%ld]=%e\n",site,ical[site]);

}

// Inward rectifier potassium current (Ik1)
void comp_ik1(int site)
{
	double vm,gk1,ek;
	double ak1,bk1,k1inf;

	vm = vmf[0][site];
	ek = Ek[site];
	gk1 = Gk1[site];


	ak1 = 0.1/(1.0+exp(0.06*(vm - ek -200.0)));
	bk1 = (3.0*exp(0.0002*(vm - ek + 100.0))+exp(0.1*(vm - ek -10.0)))/(1.0+exp(-0.5*(vm - ek)));
	k1inf = ak1/(ak1+bk1);
		
	ik1[site] = gk1*rategk1*k1inf*(vm - ek);

//	printf("ik1[%ld]=%e\n",site,ik1[site]);
}

// Sodium-Calcium Exchanger V-S

void comp_incx(int site)
{
	double c2,c3,exp0,exp1,Gncx;
	double vm,Nai,Cai;
	vm = vmf[0][site];
	Nai = vmf[17][site];
	Cai = vmf[14][site];
	Gncx = kncx[site];

	int iV=0;
	double V1,V2,d1,d2;
	V1 = (vm + Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp0 = Tncx0[iV]*d2 + Tncx0[iV+1]*d1;
	exp1 = Tncx1[iV]*d2 + Tncx1[iV+1]*d1;
	c2 = Cao*Nai*Nai*Nai;
	c3 = Nao*Nao*Nao*Cai*2.5;	// for TNNP06model:Nao*Nao*Nao*Cai*ncx.a;
	//c3 = Nao*Nao*Nao*Cai;	// for K05model
	
	//h3 = kmna*kmna*kmna + Nao*Nao*Nao;
	//h4 = kmca + Cao;
	// j = Gncx*(exp0*c2-exp1*c3)/(h3*h4*(1.0+ksat*exp1));	

	jncx[site] = Gncx*c1*(c2*exp0 - c3*exp1)/(1.0 + ksat*exp1);

//	printf("jncx[%ld]=%e\n",site,jncx[site]);
}

// Sodium-Potassium Pump
void comp_inak(int site)
{
	double exp0,exp1,fnak;
	double vm,gnak,Nai;
	//double k,rho;

	vm  = vmf[0][site];
	Nai = vmf[17][site];
	gnak = Gnak[site];
	//Kmkp = kmk;
	//Kmnap = kmna;
	//k = nna;
	//rho = rhonak;

	int iV=0;
	double V1,V2,d1,d2;
	V1 = (vm + Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp0 = Tnak0[iV]*d2 + Tnak0[iV+1]*d1;
	exp1 = Tnak1[iV]*d2 + Tnak1[iV+1]*d1;
	
	fnak = 1.0/(1.0 + exp0 + exp1);
	
	inak[site] = gnak*fnak*(Ko*Nai/((Ko + kmk)*(Nai + kmna)));
	// for K05 model
	// inak[site] = gnak*(Ko/(Ko+Kmkp))/(1.0+pow((Kmnap/Nai),nna))/(1.0+0.1245*exp0+0.0365*rho*exp1);

//	printf("inak[%ld]=%e\n",site,inak[site]);
}

// Sarcolemmal Ca Pump 
void comp_ipca(int site)
{
	double Cai,w1;
	double gpca = Gpca;
	double kmpca = Km_pca;

	Cai = vmf[14][site];
	ipca[site] = gpca*Cai/(Cai + kmpca);

//	printf("ipca[%d]=%e\n",site,ipca[site]);
}

// plateau K current 
void comp_ikp(int site)
{
	double gkp = Gkp;
	double uu,vm,ek;
	vm = vmf[0][site];
	ek = Ek[site];

	int iV=0;
	double V1,V2,d1,d2;
	V1 = (vm + Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	uu = Tkpu[iV]*d2 + Tkpu[iV+1]*d1;
	ikp[site] = gkp*uu*(vm - ek);

//	printf("ikp[%d]=%e\n",site,ikp[site]);
}

// Ca Background Current 

void comp_icab(int site)
{
	double ECa,vm,w1;
	double gcab = Gcab;

	vm = vmf[0][site];
	ECa = Eca[site];
	icab[site] = gcab*(vm - ECa);

}

// Na Background Current 

void comp_inab(int site)
{
	double ENa,vm,w1;
	double gnab = Gnab;
	
	vm = vmf[0][site];
	ENa = Ena[site];
	inab[site] = gnab*(vm - ENa);

}

void comp_jrel(int site)
{
	double Rp,Casr,Cass;
	double O,kcasr,grel;

	Rp   = vmf[13][site];
	Casr = vmf[15][site];
	Cass = vmf[16][site];
	grel = Grel[site];

	kcasr = maxsr - (maxsr - minsr)/(1.0+(EC/Casr)*(EC/Casr));
	relk1    = k1_prime/kcasr;
	relk2    = k2_prime*kcasr;

	O = relk1*Cass*Cass*Rp/(relk3 + relk1*Cass*Cass);
	jrel[site] = grel*O*(Casr -Cass);

}

void comp_jup(int site)
{
	int i;
	double gup,Kup,Cai;
	Cai = vmf[14][site];
	Kup = kup;
	gup = Gup[site];

	jup[site] = gup*Cai*Cai/(Cai*Cai + Kup*Kup);
}

void comp_jleak(int site)
{
	int i;
	double Gleak = pleak;
	double Cai,Casr;

	Cai  = vmf[14][site];
	Casr = vmf[15][site];
	jleak[site] = Gleak*(Casr - Cai);

}

//void comp_jtr(int site)
void comp_jxfer(int site)
{

	double Cai,Cass;
	Cai  = vmf[14][site];
	Cass = vmf[16][site];
	jxfer[site] = vxfer*(Cass - Cai);
}

void comp_concentration(int site)
{
	double Cai,Casr,Cass;
	Cai  = vmf[14][site];
	Casr = vmf[15][site];
	Cass = vmf[16][site];

	bufcai[site] =  1.0/(1.0 + bufc*kcai/((kcai + Cai)*(kcai + Cai)));
	bufcasr[site] = 1.0/(1.0 + bufsr*kcasr/((kcasr + Casr)*(kcasr + Casr)));
	bufcass[site] = 1.0/(1.0 + bufss*kcass/((kcass + Cass)*(kcass + Cass)));

}


// Reversal potentials */

void comp_Erev(int site)
{
	double Nai,Ki,Cai;
	double x1 = pkna;
	double x2 = prnak;
	
	Nai = vmf[17][site];
	Ki  = vmf[18][site];
	Cai = vmf[14][site];
	
	Ena[site] = RTonF *log(Nao/Nai);
	Ek[site]  = RTonF *log(Ko/Ki);
	Eks[site] = RTonF *log((x2*Nao+Ko)/(x2*Nai+Ki));
	Eca[site] = RTon2F*log(Cao/Cai);

}

