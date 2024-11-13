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

		// ito
		Tar[vindex] = 0.5266*exp(-0.0166*(v-42.2912))/(1.0+exp(-0.0943*(v-42.2912)));
		Tbr[vindex] = (0.5149*exp(-0.1344*(v-5.0027))+0.00005186*v)/(1.0+exp(-0.1348*(v-0.00005186)));
		// for Endo case
		Taq10[vindex] = (0.0721*exp(-0.173*(v+34.2531-dVgq0))+0.00005612*(v-dVgq0))/(1.0+exp(-0.1732*(v+34.2531-dVgq0)));
		Tbq10[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235-dVgq0))+0.0001215*(v-dVgq0))/(1.0+exp(-0.1604*(v+34.0235-dVgq0)));
		// for Epi case
		Taq11[vindex] = (0.0721*exp(-0.173*(v+34.2531-dVgq1))+0.00005612*(v-dVgq1))/(1.0+exp(-0.1732*(v+34.2531-dVgq1)));
		Tbq11[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235-dVgq1))+0.0001215*(v-dVgq1))/(1.0+exp(-0.1604*(v+34.0235-dVgq1)));
		// for Mid case
		Taq12[vindex] = (0.0721*exp(-0.173*(v+34.2531-dVgq2))+0.00005612*(v-dVgq2))/(1.0+exp(-0.1732*(v+34.2531-dVgq2)));
		Tbq12[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235-dVgq2))+0.0001215*(v-dVgq2))/(1.0+exp(-0.1604*(v+34.0235-dVgq2)));

		Taq2[vindex] = (0.0721*exp(-0.173*(v+34.2531))+0.00005612*v)/(1.0+exp(-0.1732*(v+34.2531)));
		Tbq2[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235))+0.0001215*v)/(1.0+exp(-0.1604*(v+34.0235)));

		// for ical
		Tdinf[vindex] = 1.0/(1.0+exp(-(v+7.64)/6.32));
		Ttaud[vindex] = (1.4/(1.0+exp(-(v+35.0)/13.0))+0.25)*1.4/(1.0+exp((v+5.0)/5.0))+1.0/(1.0+exp(-(v-50.0)/20.0));

		Tfinf[vindex] = 1.0/(1.0+exp((v+24.6)/6.9));
		Ttauf[vindex] = 23.9*0.75/(0.1389*exp(-(0.0358*(v-10.9))*(0.0358*(v-10.9)))+0.0519);

		// for ikr 
		Tpinf[vindex] = 1.0/(1.0+exp(-(v+14.0)/7.7));
		if(fabs(v+14.0)<1E-3){
			Tap[vindex] = 0.0015;
		} else {
			Tap[vindex] = 0.0003*(v+14.0)/(1.0-exp(-(v+14.0)/5.0));
		}
		if(fabs(v-3.4328)<1E-3){
			Tbp[vindex] = 11.0/29052.0;
		} else {
			Tbp[vindex] = 0.000073898*(v-3.4328)/(exp((v-3.4328)/5.1237)-1.0);
		}
		Tpiinf[vindex] = 1.0/(1.0+exp((v+15.0)/22.4));

		// for iks 
		Tninf[vindex] = 1.0/sqrt(1.0+exp(-(v-9.4)/11.8));
		Ttaun[vindex] = 555.0/(1.0+exp(-(v+22.0)/11.3))+129.0;

		// ik1 
		//ik1.Tak1[vindex] = 0.1/(1.0+exp(0.06*(v-var.Ek-200.0)));
		//ik1.Tbk1[vindex] = (3.0*exp(0.0002*(v-var.Ek+100.0))+exp(0.1*(v-var.Ek-10.0)))/(1.0+exp(-0.5*(v-var.Ek)));

		// inaca
		Tncx0[vindex] = exp(rncx*FRT*v);
		Tncx1[vindex] = exp((rncx-1.0)*FRT*v);
		// inak 
		//Tnak0[vindex] = exp(-0.1*FRT*v);
		//Tnak1[vindex] = exp(-1.0*FRT*v);

	}
}

