#include "Rcpp.hpp"
#include <R_ext/Rdynload.h>
 

RcppExport SEXP fabic(SEXP xS, SEXP PsiS,SEXP LS,SEXP laplaS,SEXP cycS, SEXP alphaS,SEXP epsS,SEXP eps1S,SEXP splS,SEXP spzS,SEXP scaleS,SEXP lapS) {

    RcppMatrix<double> xR(xS);

    int n = xR.getDim1();
    int nn = xR.getDim2();

    double **x=xR.cMatrix();


    RcppMatrix<double> LR(LS);
    double **L=        LR.cMatrix();
    if (n != LR.getDim1())
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 L %d\n",n,LR.getDim1());
	return NULL;
    }
    int K = LR.getDim2();

    RcppMatrix<double> laplaR(laplaS);
    double **lapla=    laplaR.cMatrix();

    if (nn != laplaR.getDim1())
    {
	Rprintf("STOP: Dim2 X: %d BUT Dim1 lapla %d\n",nn,laplaR.getDim1());
	return NULL;
    }
    if (K != laplaR.getDim2())
    {
	Rprintf("STOP: Dim2 L: %d BUT Dim2 lapla %d\n",K,laplaR.getDim2());
	return NULL;
    }
 
    RcppVector<double> PsiR(PsiS);
    double *Psi=     PsiR.cVector();

    if (n != PsiR.size())
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 Psi %d\n",n,PsiR.size());
	return NULL;

    }




    RcppVector<double> alphaR(alphaS);
    RcppVector<double> epsR(epsS);
    RcppVector<double> eps1R(eps1S);
    RcppVector<double> splR(splS);
    RcppVector<double> spzR(spzS);
    RcppVector<double> scaleR(scaleS);
    RcppVector<double> lapR(lapS);


    double alpha = alphaR(0);
    double eps = epsR(0);
    double eps1 = eps1R(0);
    double spl = splR(0);
    double spz = spzR(0);
    double scale = scaleR(0);
    double lap = lapR(0);



    RcppVector<int> cycR(cycS);

    int cyc = cycR(0);


    RcppVector<double> XXR(n);
    RcppVector<double> ipsiR(n);
    RcppVector<double> e_sx_nR(n);
    RcppVector<double> e_ssxx_nR(K);

    double *XX=      XXR.cVector();
    double *ipsi=    ipsiR.cVector();
    double *e_sx_n=  e_sx_nR.cVector();
    double *e_ssxx_n=e_ssxx_nR.cVector();

    RcppMatrix<double> sum1R(n,K);
    RcppMatrix<double> LPsiR(n,K);
    RcppMatrix<double> icholLPsiR(n,K);
    RcppMatrix<double> E_SX_nR(K,nn);
    RcppMatrix<double> LPsiLR(K,K);
    RcppMatrix<double> tLPsiLR(K,K);
    RcppMatrix<double> sum2R(K,K);
    RcppMatrix<double> icholR(K,K);

    double **sum1=     sum1R.cMatrix();
    double **LPsi=     LPsiR.cMatrix();
    double **icholLPsi= icholLPsiR.cMatrix();
    double **E_SX_n=   E_SX_nR.cMatrix();
    double **LPsiL=    LPsiLR.cMatrix();
    double **tLPsiL=   tLPsiLR.cMatrix();
    double **sum2=     sum2R.cMatrix();
    double **ichol=     icholR.cMatrix();


    int i,j,i1,i2,i3,zz;

    double in,s,sgn,t;



    SEXP rl;

    RcppResultSet rs;

    
    spl = -spl;
    spz = -spz;
    in = 1.0/nn;

    if (lap<eps)
    {
	lap = eps;
    }

    for (i1=0;i1<n;i1++)
    {
	s = 0.0;
	for (i2 = 0; i2 < nn; i2++)
	    s += x[i1][i2] * x[i1][i2];
	XX[i1] = s*in;
	if (XX[i1]<eps) XX[i1]=eps;

    }



    for (i=0;i<cyc;i++) {


	for (i1=0;i1<n;i1++)
	{
	    ipsi[i1]= 1.0/Psi[i1];
	}

	for (i1=0;i1<n;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		LPsi[i1][i2] =  ipsi[i1]*L[i1][i2];
	}

	for (i1=0;i1<K;i1++)
	{
	    for (i3=0;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = s;
	    }
	}


	for (i1=0;i1<n;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		sum1[i1][i2] = 0.0;
	}

	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		sum2[i1][i2] = (i1==i2 ? eps : 0.0);
	}

	for (j=0;j<nn;j++)
	{

	    for (i1=0;i1<K;i1++)
	    {
		tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
	    }


	    for (i1=0;i1<K;i1++) {
		for (i2=i1;i2<K;i2++) {
		    for (s=tLPsiL[i1][i2],i3=i1-1;i3>=0;i3--) s -= tLPsiL[i1][i3]*tLPsiL[i2][i3];
		    ichol[i1][i2]=ichol[i2][i1]=0.0;
		    if (i1 == i2) {
			tLPsiL[i1][i1]=sqrt(s);
		    } else tLPsiL[i2][i1]=s/tLPsiL[i1][i1];
		}
	    }


	    for (i1=0;i1<K;i1++) 
		for (i2=0;i2<=i1;i2++){
		    s = (i1==i2 ? 1.0 : 0.0);
		    for (i3=i1-1;i3>=i2;i3--) s -= tLPsiL[i1][i3]*ichol[i2][i3];
		    ichol[i2][i1]= s/tLPsiL[i1][i1];
		}


	    for (i1=K-1;i1>=0;i1--) 
		for (i2=0;i2<=i1;i2++){
		    s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		    for (i3=i1+1;i3<K;i3++) s -= tLPsiL[i3][i1]*ichol[i2][i3];
		    ichol[i1][i2] = ichol[i2][i1] = s/tLPsiL[i1][i1];
		}

	    for (i1=0;i1<n;i1++)
	    {
		for (i3 = 0; i3 < K; i3++){
		    s=0.0;
		    for (i2 = 0; i2 < K; i2++)
			s +=  ichol[i3][i2]*LPsi[i1][i2];
		    icholLPsi[i1][i3] = s;
		}
	    }


	    for (i2 = 0; i2 < K; i2++) {
		s = 0.0;
		for (i1=0;i1<n;i1++)
		{
		    s +=  icholLPsi[i1][i2]*x[i1][j];
		}
		e_sx_n[i2] = s;
	    }


	    for (i1=0;i1<n;i1++)
	    {
		for (i2 = 0; i2 < K; i2++)
		    sum1[i1][i2] += x[i1][j]*e_sx_n[i2];
	    }


	    for (i1=0;i1<K;i1++)
	    {
		for (i2 = 0; i2 < K; i2++)
		{
		    s = ichol[i1][i2] + e_sx_n[i1]*e_sx_n[i2];
		    sum2[i1][i2] += s;
		    if (i1==i2) e_ssxx_n[i1] = s;
		}
	    }

	    for (i1=0;i1<K;i1++)
	    {
		s = pow((eps+e_ssxx_n[i1]),spz);
		if (s<lap)
		{
		    lapla[j][i1] = lap;
		} else {
		    lapla[j][i1] = s; 
		}
	    }
	}


	for (i1=0;i1<K;i1++) {
	    for (i2=i1;i2<K;i2++) {
		for (s=sum2[i1][i2],i3=i1-1;i3>=0;i3--) s -= sum2[i1][i3]*sum2[i2][i3];
		ichol[i1][i2]=ichol[i2][i1]=0.0;
		if (i1 == i2) {
		    sum2[i1][i1]=sqrt(s);
		} else sum2[i2][i1]=s/sum2[i1][i1];
	    }
	}

	for (i1=0;i1<K;i1++) 
	    for (i2=0;i2<=i1;i2++){
		s = (i1==i2 ? 1.0 : 0.0);
		for (i3=i1-1;i3>=i2;i3--) s -= sum2[i1][i3]*ichol[i2][i3];
		ichol[i2][i1]= s/sum2[i1][i1];
	    }

	for (i1=K-1;i1>=0;i1--) 
	    for (i2=0;i2<=i1;i2++){
		s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		for (i3=i1+1;i3<K;i3++) s -= sum2[i3][i1]*ichol[i2][i3];
		ichol[i1][i2] = ichol[i2][i1] = s/sum2[i1][i1];
	    }


	for (i1=0;i1<n;i1++)
	{
	    for (i3 = 0; i3 < K; i3++) {
		s = 0.0;
		for (i2 = 0; i2 < K; i2++) {
		    s += sum1[i1][i2]*ichol[i2][i3];
		}
		L[i1][i3] = s;
	    }
	}
	for (i1=0;i1<n;i1++)
	{
	    for (i3 = 0; i3 < K; i3++) {
		s = 0.0;
		for (i2 = 0; i2 < K; i2++) {
		    s += pow((eps1+fabs(L[i1][i3])),spl)*ichol[i2][i3];
		}
		LPsi[i1][i3] = s;
	    }
	}


	zz=0;
	for (i1=0;i1<n;i1++)
	{
	    for (i2 = 0; i2 < K; i2++) {
	      s=L[i1][i2];
	      t = fabs(Psi[i1]*alpha*LPsi[i1][i2]);
	      if (fabs(s)>t){
		sgn = (s>0.0) ? 1.0 : ((s == 0.0) ? 0.0 : -1.0);
		L[i1][i2] = s - sgn* t;
		zz=1;
	      }
	      else {
		L[i1][i2] = 0;
	      }
	    }
	}

	if (zz==0)
	{

	    for (i1=0;i1<n;i1++)
	    {
		for (i2 = 0; i2 < K; i2++)
		{
		    L[i1][i2] = 0.1;
		}
	    }


	}


	for (i1=0;i1<n;i1++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < K; i2++) {
		s += L[i1][i2]*sum1[i1][i2];
	    }
	    ipsi[i1] = s;
	}

	for (i1=0;i1<n;i1++)
	{
	    Psi[i1] = XX[i1] - in*ipsi[i1];
	    if (Psi[i1]<eps)
	    {
		Psi[i1] = eps;
	    }
	}

	if (scale>0)
	{
	    for (i2 = 0; i2 < K; i2++)
	    {
		
	    
		s = 0.0;
		for (i1=0;i1<n;i1++)
		{
		    s+=L[i1][i2]*L[i1][i2];
		}
		s*=in;
		s=sqrt(s);
		s=scale/s;
		for (i1=0;i1<n;i1++)
		{
		    L[i1][i2]*=s;
		}
		s*=s;
		s = pow(s,spz);
		for (j=0;j<nn;j++)
		{
		    lapla[j][i2]*=s;
		}
	    
	    }
	}		

	if (i%20==0) {
	    Rprintf("Cycle: %d\n", i);
	    R_CheckUserInterrupt();
	}

    }



    for (i1=0;i1<n;i1++)
    {
	ipsi[i1]= 1.0/Psi[i1];
    }

    for (i1=0;i1<n;i1++)
    {
	for (i2 = 0; i2 < K; i2++)
	    LPsi[i1][i2] =  ipsi[i1]*L[i1][i2];
    }


    for (i1=0;i1<K;i1++)
    {
	for (i3=0;i3<K;i3++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < n; i2++)
		s += LPsi[i2][i1]*L[i2][i3];
	    LPsiL[i1][i3] = s;
	}
    }


    for (j=0;j<nn;j++)
    {

	for (i1=0;i1<K;i1++)
	{
	    tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
	}

	for (i1=0;i1<K;i1++) {
	    for (i2=i1;i2<K;i2++) {
		for (s=tLPsiL[i1][i2],i3=i1-1;i3>=0;i3--) s -= tLPsiL[i1][i3]*tLPsiL[i2][i3];
		ichol[i1][i2]=ichol[i2][i1]=0.0;
		if (i1 == i2) {
		    tLPsiL[i1][i1]=sqrt(s);
		} else tLPsiL[i2][i1]=s/tLPsiL[i1][i1];
	    }
	}

	for (i1=0;i1<K;i1++) for (i2=0;i2<=i1;i2++){
		s = (i1==i2 ? 1.0 : 0.0);
		for (i3=i1-1;i3>=i2;i3--) s -= tLPsiL[i1][i3]*ichol[i2][i3];
		ichol[i2][i1]= s/tLPsiL[i1][i1];
	    }

	for (i1=K-1;i1>=0;i1--) 
	    for (i2=0;i2<=i1;i2++){
		s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		for (i3=i1+1;i3<K;i3++) s -= tLPsiL[i3][i1]*ichol[i2][i3];
		ichol[i1][i2] = ichol[i2][i1] = s/tLPsiL[i1][i1];
	    }


	for (i1=0;i1<n;i1++)
	{
	    for (i3 = 0; i3 < K; i3++){
		s=0.0;
		for (i2 = 0; i2 < K; i2++)
		    s +=  ichol[i3][i2]*LPsi[i1][i2];
		icholLPsi[i1][i3] = s;
	    }
	}


	for (i2 = 0; i2 < K; i2++) {
	    s = 0.0;
	    for (i1=0;i1<n;i1++)
	    {
		s +=  icholLPsi[i1][i2]*x[i1][j];
	    }
	    E_SX_n[i2][j] = s;
	}

    }



    rs.add("L", L,n,K);
    rs.add("E_SX_n", E_SX_n,K,nn);
    rs.add("Psi",Psi,n);
    rs.add("lapla", lapla,nn,K);




    rl = rs.getReturnList();




    return rl;

}







RcppExport SEXP fabics(SEXP xS, SEXP PsiS,SEXP LS,SEXP laplaS,SEXP cycS, SEXP alphaS,SEXP epsS,SEXP spzS,SEXP lapS) {



    RcppMatrix<double> xR(xS);

    int n = xR.getDim1();
    int nn = xR.getDim2();

    double **x=xR.cMatrix();


    RcppMatrix<double> LR(LS);
    double **L=        LR.cMatrix();
    if (n != LR.getDim1())
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 L %d\n",n,LR.getDim1());
	return NULL;
    }
    int K = LR.getDim2();

    RcppMatrix<double> laplaR(laplaS);
    double **lapla=    laplaR.cMatrix();

    if (nn != laplaR.getDim1())
    {
	Rprintf("STOP: Dim2 X: %d BUT Dim1 lapla %d\n",nn,laplaR.getDim1());
	return NULL;
    }
    if (K != laplaR.getDim2())
    {
	Rprintf("STOP: Dim2 L: %d BUT Dim2 lapla %d\n",K,laplaR.getDim2());
	return NULL;
    }
 
    RcppVector<double> PsiR(PsiS);
    double *Psi=     PsiR.cVector();

    if (n != PsiR.size())
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 Psi %d\n",n,PsiR.size());
	return NULL;

    }


    RcppVector<double> alphaR(alphaS);
    RcppVector<double> epsR(epsS);
    RcppVector<double> spzR(spzS);
    RcppVector<double> lapR(lapS);

 
    double alpha = alphaR(0);
    double eps = epsR(0);
    double spz = spzR(0);
    double lap = lapR(0);



    RcppVector<int> cycR(cycS);

    int cyc = cycR(0);


    RcppVector<double> XXR(n);
    RcppVector<double> vR(n);
    RcppVector<double> wR(n);
    RcppVector<double> ipsiR(n);
    RcppVector<double> e_sx_nR(n);
    RcppVector<double> e_ssxx_nR(K);

    double *XX=      XXR.cVector();
    double *w=       wR.cVector();
    double *v=       vR.cVector();
    double *ipsi=    ipsiR.cVector();
    double *e_sx_n=  e_sx_nR.cVector();
    double *e_ssxx_n=e_ssxx_nR.cVector();

    RcppMatrix<double> sum1R(n,K);
    RcppMatrix<double> LPsiR(n,K);
    RcppMatrix<double> icholLPsiR(n,K);
    RcppMatrix<double> E_SX_nR(K,nn);
    RcppMatrix<double> LPsiLR(K,K);
    RcppMatrix<double> tLPsiLR(K,K);
    RcppMatrix<double> sum2R(K,K);
    RcppMatrix<double> icholR(K,K);

    double **sum1=     sum1R.cMatrix();
    double **LPsi=     LPsiR.cMatrix();
    double **icholLPsi= icholLPsiR.cMatrix();
    double **E_SX_n=   E_SX_nR.cMatrix();
    double **LPsiL=    LPsiLR.cMatrix();
    double **tLPsiL=   tLPsiLR.cMatrix();
    double **sum2=     sum2R.cMatrix();
    double **ichol=     icholR.cMatrix();

    RcppVector<int> zerosR(n);

    int *zeros=      zerosR.cVector();

    int i,j,i1,i2,i3,zz,ende,h1;

    double in,s,sgn,a,b,c,t,alphap,k2,seps;



    SEXP rl;

    RcppResultSet rs;

    seps=10e-10;
    alpha = sqrt((1.0*n))-(sqrt((1.0*n))-1.0)*alpha;
    in = 1.0/nn;
    spz = -spz;

    if (lap<eps)
    {
	lap = eps;
    }

    for (i1=0;i1<n;i1++)
    {
	s = 0.0;
	for (i2 = 0; i2 < nn; i2++)
	    s += x[i1][i2] * x[i1][i2];
	XX[i1] = s*in;
	if (XX[i1]<eps) XX[i1]=eps;

    }


 
    for (i=0;i<cyc;i++) {


	for (i1=0;i1<n;i1++)
	{
	    ipsi[i1]= 1.0/Psi[i1];
	}

	for (i1=0;i1<n;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		LPsi[i1][i2] =  ipsi[i1]*L[i1][i2];
	}

	for (i1=0;i1<K;i1++)
	{
	    for (i3=0;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = s;
	    }
	}


	for (i1=0;i1<n;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		sum1[i1][i2] = 0.0;
	}

	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < K; i2++)
		sum2[i1][i2] = (i1==i2 ? eps : 0.0);
	}

	for (j=0;j<nn;j++)
	{

	    for (i1=0;i1<K;i1++)
	    {
		tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
	    }


	    for (i1=0;i1<K;i1++) {
		for (i2=i1;i2<K;i2++) {
		    for (s=tLPsiL[i1][i2],i3=i1-1;i3>=0;i3--) s -= tLPsiL[i1][i3]*tLPsiL[i2][i3];
		    ichol[i1][i2]=ichol[i2][i1]=0.0;
		    if (i1 == i2) {
			tLPsiL[i1][i1]=sqrt(s);
		    } else tLPsiL[i2][i1]=s/tLPsiL[i1][i1];
		}
	    }


	    for (i1=0;i1<K;i1++) for (i2=0;i2<=i1;i2++){
		    s = (i1==i2 ? 1.0 : 0.0);
		    for (i3=i1-1;i3>=i2;i3--) s -= tLPsiL[i1][i3]*ichol[i2][i3];
		    ichol[i2][i1]= s/tLPsiL[i1][i1];
		}


	    for (i1=K-1;i1>=0;i1--) for (i2=0;i2<=i1;i2++){
		    s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		    for (i3=i1+1;i3<K;i3++) s -= tLPsiL[i3][i1]*ichol[i2][i3];
		    ichol[i1][i2] = ichol[i2][i1] = s/tLPsiL[i1][i1];
		}


	    for (i1=0;i1<n;i1++)
	    {
		for (i3 = 0; i3 < K; i3++){
		    s=0.0;
		    for (i2 = 0; i2 < K; i2++)
			s +=  ichol[i3][i2]*LPsi[i1][i2];
		    icholLPsi[i1][i3] = s;
		}
	    }


	    for (i2 = 0; i2 < K; i2++) {
		s = 0.0;
		for (i1=0;i1<n;i1++)
		{
		    s +=  icholLPsi[i1][i2]*x[i1][j];
		}
		e_sx_n[i2] = s;
	    }


	    for (i1=0;i1<n;i1++)
	    {
		for (i2 = 0; i2 < K; i2++)
		    sum1[i1][i2] += x[i1][j]*e_sx_n[i2];
	    }


	    for (i1=0;i1<K;i1++)
	    {
		for (i2 = 0; i2 < K; i2++)
		{
		    s = ichol[i1][i2] + e_sx_n[i1]*e_sx_n[i2];
		    sum2[i1][i2] += s;
		    if (i1==i2) e_ssxx_n[i1] = s;
		}
	    }

	    for (i1=0;i1<K;i1++)
	    {
		s = pow((eps+e_ssxx_n[i1]),spz);
		if (s<lap)
		{
		    lapla[j][i1] = lap;
		} else {
		    lapla[j][i1] = s; 
		}
	    }
	}


	for (i1=0;i1<K;i1++) {
	    for (i2=i1;i2<K;i2++) {
		for (s=sum2[i1][i2],i3=i1-1;i3>=0;i3--) s -= sum2[i1][i3]*sum2[i2][i3];
		ichol[i1][i2]=ichol[i2][i1]=0.0;
		if (i1 == i2) {
		    sum2[i1][i1]=sqrt(s);
		} else sum2[i2][i1]=s/sum2[i1][i1];
	    }
	}

	for (i1=0;i1<K;i1++) for (i2=0;i2<=i1;i2++){
		s = (i1==i2 ? 1.0 : 0.0);
		for (i3=i1-1;i3>=i2;i3--) s -= sum2[i1][i3]*ichol[i2][i3];
		ichol[i2][i1]= s/sum2[i1][i1];
	    }

	for (i1=K-1;i1>=0;i1--) for (i2=0;i2<=i1;i2++){
		s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		for (i3=i1+1;i3<K;i3++) s -= sum2[i3][i1]*ichol[i2][i3];
		ichol[i1][i2] = ichol[i2][i1] = s/sum2[i1][i1];
	    }



	for (i1=0;i1<n;i1++)
	{
	    for (i3 = 0; i3 < K; i3++) {
		s = 0.0;
		for (i2 = 0; i2 < K; i2++) {
		    s += sum1[i1][i2]*ichol[i2][i3];
		}
		L[i1][i3] = s;
	    }
	}









	for (i1=0;i1<n;i1++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < K; i2++) {
		s += L[i1][i2]*sum1[i1][i2];
	    }
	    ipsi[i1] = s;
	}

	for (i1=0;i1<n;i1++)
	{
	    Psi[i1] = XX[i1] - in*ipsi[i1];
	    if (Psi[i1]<eps)
	    {
		Psi[i1] = eps;
	    }
	}


//----------

        for (i3 = 0; i3 < K; i3++) {

	    s = 0.0;
	    for (i1=0;i1<n;i1++)
	    {
		e_sx_n[i1] = (L[i1][i3]<0.0 ? -1.0 : 1.0);
		sgn=fabs(L[i1][i3]);
		ipsi[i1] = sgn;
		s += sgn;
	    }


	    k2=1.0;

	    zz=0;
	    s=(alpha-s)/(1.0*n);
	    for (i1=0;i1<n;i1++)
	    {
		sgn = ipsi[i1] + s;
		if (sgn<=0.0) {
		    zeros[i1] = 2;
		    zz++;
		    v[i1] = 0.0;
		}
		else {
		    zeros[i1] = 0;
		    v[i1] = sgn;
		}

	    }


	    j=0;
	    ende=0;
	    while (ende<1) {

		a=0.0;
		b=0.0;
		c=0.0;
		h1 = n-zz;
		s=0.0;
		if (h1>0) s=alpha/(1.0*h1);
		for (i1=0;i1<n;i1++)
		{
		    if (zeros[i1]<1) {
			sgn = v[i1]-s;

			a+=sgn*sgn;
			b+=sgn*v[i1];
			c+=v[i1]*v[i1];
			w[i1]=sgn;
		    }
		}
		c -= k2;
		b*=2.0;
		t=b*b-4.0*a*c;
		if (t<0.0) t=0.0;
		if (a<seps) a = seps;
		alphap = (-b+sqrt(t))/(2.0*a);

		s=0.0;
		ende=2;
		for (i1=0;i1<n;i1++)
		{
		    if (zeros[i1]<1) {
			v[i1] += alphap*w[i1];
			if (v[i1]<=0.0) {
			    ende = 0;
			    zeros[i1] = 2;
			    zz++;
			    v[i1]=0.0;
			}
			else {
			    s+=v[i1];
			}
		    }

		}
		if (j<2) ende = 0;
		if (ende<1) {
		    j++;
		    h1 = n-zz;
		    if (h1>0)
		    {
			s= (alpha-s)/(1.0*h1);
			for (i1=0;i1<n;i1++)
			{
			    if (zeros[i1]<1) v[i1] += s;
			}

		    }
		    else
		    {
			ende=2;
		    }
		}

	    }

	    for (i1=0;i1<n;i1++)
	    {
		L[i1][i3]= e_sx_n[i1]*v[i1];
//		Rprintf("%d %lf %lf\n",i1,L[i1][i3],e_sx_n[i1]*v[i1]);
	    }
	}

//----------






	if (i%20==0) {
	    Rprintf("Cycle: %d\n", i);
	    R_CheckUserInterrupt();
	}
    }



    for (i1=0;i1<n;i1++)
    {
	ipsi[i1]= 1.0/Psi[i1];
    }

    for (i1=0;i1<n;i1++)
    {
	for (i2 = 0; i2 < K; i2++)
	    LPsi[i1][i2] =  ipsi[i1]*L[i1][i2];
    }


    for (i1=0;i1<K;i1++)
    {
	for (i3=0;i3<K;i3++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < n; i2++)
		s += LPsi[i2][i1]*L[i2][i3];
	    LPsiL[i1][i3] = s;
	}
    }


    for (j=0;j<nn;j++)
    {

	for (i1=0;i1<K;i1++)
	{
	    tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
	}

	for (i1=0;i1<K;i1++) {
	    for (i2=i1;i2<K;i2++) {
		for (s=tLPsiL[i1][i2],i3=i1-1;i3>=0;i3--) s -= tLPsiL[i1][i3]*tLPsiL[i2][i3];
		ichol[i1][i2]=ichol[i2][i1]=0.0;
		if (i1 == i2) {
		    tLPsiL[i1][i1]=sqrt(s);
		} else tLPsiL[i2][i1]=s/tLPsiL[i1][i1];
	    }
	}

	for (i1=0;i1<K;i1++) for (i2=0;i2<=i1;i2++){
		s = (i1==i2 ? 1.0 : 0.0);
		for (i3=i1-1;i3>=i2;i3--) s -= tLPsiL[i1][i3]*ichol[i2][i3];
		ichol[i2][i1]= s/tLPsiL[i1][i1];
	    }

	for (i1=K-1;i1>=0;i1--) for (i2=0;i2<=i1;i2++){
		s = (i1<i2 ? 0.0 : ichol[i2][i1]);
		for (i3=i1+1;i3<K;i3++) s -= tLPsiL[i3][i1]*ichol[i2][i3];
		ichol[i1][i2] = ichol[i2][i1] = s/tLPsiL[i1][i1];
	    }


	for (i1=0;i1<n;i1++)
	{
	    for (i3 = 0; i3 < K; i3++){
		s=0.0;
		for (i2 = 0; i2 < K; i2++)
		    s +=  ichol[i3][i2]*LPsi[i1][i2];
		icholLPsi[i1][i3] = s;
	    }
	}


	for (i2 = 0; i2 < K; i2++) {
	    s = 0.0;
	    for (i1=0;i1<n;i1++)
	    {
		s +=  icholLPsi[i1][i2]*x[i1][j];
	    }
	    E_SX_n[i2][j] = s;
	}

    }



    rs.add("L", L,n,K);
    rs.add("E_SX_n", E_SX_n,K,nn);
    rs.add("Psi",Psi,n);
    rs.add("lapla", lapla,nn,K);




    rl = rs.getReturnList();




    return rl;

}

 R_CallMethodDef callMethods[]  = {
       {"fabic", (DL_FUNC) &fabic, 12},
       {"fabics", (DL_FUNC) &fabics, 9},
       {NULL, NULL, 0}
     };



 void  R_init_myLib(DllInfo *info)
     {
	 R_registerRoutines(info, NULL, callMethods, NULL, NULL);
     }


int main() {

}


