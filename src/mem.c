#include <string.h>
#include "syspara.h"

typedef double Number;
typedef long long Lint;

void initial_mem()
{
	int i,k,m,l,z;
	int site;
// initialized memory for state variables and arguments for Eular method 
	vmf=(double**)malloc2d(sizeof(double*),NN,MEDIA_SITE);
	dvmf=(double**)malloc2d(sizeof(double*),NN,MEDIA_SITE);

	//var.u=(double**)malloc2d(sizeof(double*),NN,MEDIA_SITE);
	//var.f=(double**)malloc2d(sizeof(double*),NN,MEDIA_SITE);
	
	imo=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if(imo==NULL) exit(1);
	itotal=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if(itotal==NULL) exit(1);
	celltype=(int *)calloc(MEDIA_SITE,sizeof(int));
	if(celltype==NULL) exit(1);
	utype=(int *)calloc(MEDIA_SITE,sizeof(int));
	if(utype==NULL) exit(1);

	Enaf=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Ena=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Ek=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Eto=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Eks=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Eca=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Enaf==NULL || Ena== NULL || Ek==NULL || Eto ==NULL || Eks==NULL || Eca==NULL) exit(1);

	jnanet=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jknet=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jcainet=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jcasnet=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if(jnanet==NULL || jknet==NULL || jcainet==NULL || jcasnet==NULL) exit(1);

// initialized tablization memorys for Exp functions
	// ina_fast
	Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	Thss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauh=(Number *)calloc(VNMAX,sizeof(Number));
	Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauj=(Number *)calloc(VNMAX,sizeof(Number));

	Gnaf=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	inaf=(Number *)calloc(MEDIA_SITE,sizeof(Number));

	mss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	hss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	taum=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauh=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauj=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Tmss==NULL || Ttaum==NULL || Thss==NULL || Ttauh==NULL || Tjss==NULL || Ttauj == NULL
		|| Gnaf==NULL || inaf==NULL 
		|| mss==NULL || hss==NULL || jss==NULL || taum ==NULL || tauh == NULL || tauj == NULL ) exit(1);
	// ito
	Trss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttaur=(Number *)calloc(VNMAX,sizeof(Number));
	Tsss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttaus=(Number *)calloc(VNMAX,sizeof(Number));
	ito=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Gto=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	rss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	sss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	taur=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	taus=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Trss==NULL || Ttaur==NULL || Tsss==NULL || Ttaus==NULL
		|| ito==NULL || Gto==NULL || rss==NULL || sss==NULL || taur==NULL || taus==NULL) exit(1);
	// ikr
	Txr1ss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauxr1=(Number *)calloc(VNMAX,sizeof(Number));
	Txr2ss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauxr2=(Number *)calloc(VNMAX,sizeof(Number));
	ikr=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Gkr=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	xr1ss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauxr1=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	xr2ss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauxr2=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	gkr_rate=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Txr1ss==NULL || Ttauxr1==NULL || Txr2ss==NULL || Ttauxr2==NULL || 
		ikr==NULL || Gkr==NULL || xr1ss==NULL || tauxr1==NULL || xr2ss==NULL || tauxr2==NULL || gkr_rate==NULL) exit(1);
	// iks
	Txsss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauxs=(Number *)calloc(VNMAX,sizeof(Number));
	iks=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Gks=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	xsss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauxs=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	gks_rate=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Txsss==NULL || Ttauxs==NULL || iks==NULL || Gks==NULL || xsss==NULL || tauxs==NULL || gks_rate==NULL) exit(1);
	// ical
	Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauf=(Number *)calloc(VNMAX,sizeof(Number));
	Tf2ss=(Number *)calloc(VNMAX,sizeof(Number));
	Ttauf2=(Number *)calloc(VNMAX,sizeof(Number));
	Trcal=(Number *)calloc(VNMAX,sizeof(Number));
	ical=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Gca=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	dss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	fss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	f2ss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	fcass=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	taud=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauf=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	tauf2=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	taufca=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Tdss==NULL || Ttaud==NULL || Tfss==NULL || Ttauf==NULL || Tf2ss==NULL || Ttauf2==NULL || Trcal==NULL ||
		ical==NULL || Gca==NULL || dss==NULL || fss==NULL || f2ss==NULL || fcass==NULL ||
		taud==NULL || tauf==NULL || tauf2==NULL || taufca==NULL) exit(1);
	// ik1
	ik1=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Gk1=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( ik1 == NULL || Gk1==NULL ) exit(1);
	// NCX
	Tncx0=(Number *)calloc(VNMAX,sizeof(Number));
	Tncx1=(Number *)calloc(VNMAX,sizeof(Number));
	kncx=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jncx=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Tncx0==NULL || Tncx1==NULL || kncx==NULL || jncx==NULL ) exit(1);
	// NaK 
	Tnak0=(Number *)calloc(VNMAX,sizeof(Number));
	Tnak1=(Number *)calloc(VNMAX,sizeof(Number));
	Gnak=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	inak=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Tnak0==NULL || Tnak1==NULL || Gnak==NULL || inak==NULL ) exit(1);
	// Ikp
	Tkpu=(Number *)calloc(VNMAX,sizeof(Number));
	ikp=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// IpCa
	ipca=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// INaB
	inab=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// ICaB
	icab=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( Tkpu == NULL || ikp == NULL || ipca==NULL || inab==NULL || icab==NULL ) exit(1);
	// Jrel
	//drss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	//dfss=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	Grel=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jrel=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// Jup
	Gup=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	jup=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// Jleak
	jleak=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	// Jtr
	jxfer=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( jrel==NULL || Grel==NULL || jup==NULL || Gup==NULL || jleak==NULL || jxfer==NULL) exit(1);
	// Ca Buffer, Calmodulin
	bufcai=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	bufcasr=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	bufcass=(Number *)calloc(MEDIA_SITE,sizeof(Number));
	if( bufcai==NULL || bufcasr==NULL || bufcass==NULL ) exit(1);

}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	int i,j;

		free(vmf);free(dvmf);
		//free(var.u);//free(var.f);
		free(celltype);free(utype);
		
		free(Ena);free(Ek);free(Eca);
		free(Enaf); free(Eto); free(Eks);
		free(jnanet);free(jknet);free(jcainet);free(jcasnet);
		free(imo);free(itotal);

		// INa
		free(Tmss);free(Ttaum);free(Thss);free(Ttauh);free(Tjss);free(Ttauj);
		free(Gnaf);free(inaf);
		free(mss);free(hss);free(jss);free(taum);free(tauh);free(tauj);
		// Ito
		free(Trss);free(Ttaur);free(Tsss);free(Ttaus);
		free(ito);free(Gto);free(rss);free(taur);free(sss);free(taus);
		// IKr
		free(Txr1ss);free(Ttauxr1);free(Txr2ss);free(Ttauxr2);free(ikr);free(Gkr);free(xr1ss);free(tauxr1);free(xr2ss);free(tauxr2);free(gkr_rate);
		// IKs
		free(Txsss);free(Ttauxs);free(iks);free(Gks);free(xsss);free(tauxs);free(gks_rate);
		// ICaL
		free(Tdss);free(Ttaud);free(Tfss);free(Ttauf);free(Tf2ss);free(Ttauf2);free(Trcal);
		free(ical);free(Gca);
		free(dss);free(fss);free(f2ss);free(taud);free(tauf);free(tauf2),free(fcass),free(taufca);
		// IK1
		free(ik1); free(Gk1);
		// Incx
		free(Tncx0); free(Tncx1);free(kncx);free(jncx);
		// INaK
		free(Tnak0);free(Tnak1);free(Gnak);free(inak);
		// other currents
		free(ipca);free(Tkpu);free(ikp);free(inab);free(icab);
		// Ca handling
		//free(drss);free(dfss);
		free(jrel);free(Grel);free(jup);free(Gup);free(jleak);free(jxfer);
		free(bufcai);free(bufcasr);free(bufcass);

}

// メモリ領域が連続な2次元配列
void *malloc2d(size_t size, int row, int col)
{

	char **a, *b;
	int t = size*col;
	int i;

	// インデックスと要素を一気に確保
	a = (char**)malloc((sizeof(*a) + t) * row);

	if (a) {
		// [インデックス, インデックス, ..., 要素, 要素, 要素, ...]と整列させるため要素の開始位置をずらす
		b = (char*)(a + row);

		// 各行の先頭アドレスを与える
		for (i = 0; i < row; i++) {
			a[i] = b;
			b += t; // 要素のサイズ×列の長さの分だけずらす
		}
		return a;
	}

	return NULL;
}

// メモリ領域が連続な3次元配列
void *malloc3d(size_t size, int i, int j, int k)
{

	char ***a, **b, *c;
	int t = size*k;
	int idx1,idx2;

	// インデックスと要素を一気に確保
	a = (char***)malloc((sizeof(*a) + sizeof(**a) * j + t * j) * i);

	if (a) {
		b = (char**)(a + i);
		c = (char *)(b + i * j);

		for (idx1 = 0; idx1 < i; idx1++) {
			a[idx1] = b;
			for (idx2 = 0; idx2 < j; idx2++) {
				b[idx2] = c;
				c += t;
			}
			b += j;
		}
		return a;
	}

	return NULL;
}

