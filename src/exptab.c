#include <string.h>
#include "syspara.h"

void make_ExpTable()
{

	int vindex;
	double v;
	double ah,bh,ad,bd,aj,bj,ae,be,s;
	double h_tau,d_tau,j_tau,e_tau;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-(double)Emax;
		
	/** for ina **/
		/** m-gate for ina **/
		Tmss[vindex] = 1.0/((1.0+exp((-v-56.86)/9.03))*(1.0+exp((-v-56.86)/9.03)));
		Ttaum[vindex] = (1.0/(1.0+exp((-60.0-v)/5)))*(0.1/(1.0+exp((v+35.0)/5.0))+0.1/(1.0+exp((v-50)/200)));

		s = 1.0/(1.0+exp(v+40.0));
		ah = 0.057*exp(-(v+80.0)/6.8);
		ad = 0.0;
		bh = 2.7*exp(0.079*v)+3.1E+5*exp(0.3485*v);
		bd = 0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
		h_tau = 1.0/(ah+bh);
		d_tau = 1.0/(ad+bd);

		/** h-gate for ina **/
		Thss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
		Ttauh[vindex] = h_tau*s + d_tau*(1.0-s);
	
		/** j-gate for ina **/
		Tjss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));

		aj = (v+37.78)*(-2.5428E+4*exp(0.2444*v)-6.948E-6*exp(-0.04391*v))/(1.0+exp(0.311*(v+79.23)));
		ae = 0.0;

		bj = 0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
		be = 0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
		j_tau = 1.0/(aj+bj);
		e_tau = 1.0/(ae+be);

		Ttauj[vindex] = j_tau*s + e_tau*(1.0-s);

		// for ito
		Trss[vindex] = 1.0/(1.0+exp(-(v-20.0)/6.0));
		Ttaur[vindex] = 9.5*exp(-(v+40.0)*(v+40.0)/1800.0)+0.8;
		if(simtype == 0){ // simtype == 0: Endo
			Tsss[vindex] = 1.0/(1.0+exp((v+28.0)/5.0));
			Ttaus[vindex] = 1000.0*exp(-(v+67.0)*(v+67.0)/1000.0)+8.0;
		} else if(simtype == 1 || simtype == 2){ // simtype == 1: Epi, simtype == 2: Mid
			Tsss[vindex] = 1.0/(1.0+exp((v+20.0)/5.0));
			Ttaus[vindex] = 85.0*exp(-(v+45.0)*(v+45.0)/320.0)+5.0/(1.0+exp((v-20.0)/5.0))+3.0;
		}

		// for ikr 
		Txr1ss[vindex] = 1.0/(1.0+exp(-(v+26.0)/7.0));
		//for K05model = 1.0/(1.0+exp(-(v+14.0)/7.7));
		Ttauxr1[vindex] = (450/(1.0+exp((-45.0-v)/10.0)))*(6.0/(1.0+exp((v+30.0)/11.5)));

		Txr2ss[vindex] = 1.0/(1.0+exp((v+88.0)/24.0));
		Ttauxr2[vindex] = (3.0/(1.0+exp((-60-v)/20.0)))*(1.12/(1.0+exp((v-60.0)/20.0)));

		// for iks 
		 Txsss[vindex] = 1.0/(1.0+exp((-5.0-v)/14.0));
		//for K05model = 1.0/sqrt(1.0+exp(-(v-9.4)/11.8));
		Ttauxs[vindex] = (1400.0/sqrt(1.0+exp((5.0-v)/6.0)))*(1.0/(1.0+exp((v-35.0)/15.0)))+80.0;

		// for ical
		Tdss[vindex] = 1.0/(1.0+exp((-8.0-v)/7.5));
		Ttaud[vindex] = (0.25+1.4/(1.0+exp((-35.0-v)/13.0)))*(1.4/(1.0+exp((v+5.0)/5.0)))+(1.0/(1.0+exp((50.0-v)/20.0)));

		Tfss[vindex] = 1.0/(1.0+exp((v+20.0)/7.0));
		Ttauf[vindex] = ratef*(1102.5*exp(-(v+27.0)*(v+27.0)/225.0)+200.0/(1.0+exp((13.0-v)/10.0))+180.0/(1.0+exp((v+30.0)/10.0))+20.0);

		Tf2ss[vindex] = 0.67/(1.0+exp((v+35.0)/7.0))+0.33;
		Ttauf2[vindex] = ratef*(600.0*exp(-(v+25.0)*(v+25.0)/170.0)+31.0/(1.0+exp((25.0-v)/10.0))+16.0/(1.0+exp((v+30.0)/10.0)));

		Trcal[vindex] = exp((v-15.0)/RTon2F);

		// incx
		Tncx0[vindex] = exp(rncx*v/RTonF);
		Tncx1[vindex] = exp((rncx-1.0)*v/RTonF);

		// inak 
		Tnak0[vindex] = 0.1245*exp((-0.1*v)/RTonF);
		Tnak1[vindex] = 0.0353*exp((-1.0*v)/RTonF);
		
		// ikp
		Tkpu[vindex] = 1.0/(1.0+exp((25.0-v)/5.98));

	}
}

