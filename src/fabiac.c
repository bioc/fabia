#include <stdlib.h>
#include <stdio.h>  
#include <limits.h> 
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


SEXP fabic(SEXP xS, SEXP PsiS,SEXP LS,SEXP laplaS,SEXP cycS, SEXP alphaS,SEXP epsS,SEXP eps1S,SEXP splS,SEXP spzS,SEXP scaleS,SEXP lapS,SEXP nLS,SEXP lLS,SEXP bLS,SEXP non_negativeS) {

    int i,j,i1,i2,i3;

    double in,s,sgn,t;


    if(!isNumeric(xS) || !isMatrix(xS)) {
	Rprintf("STOP: X not a numeric matrix\n");
	return R_NilValue;
    }
    
    SEXP dimAttr;

    dimAttr=getAttrib(xS, R_DimSymbol);
    int n=INTEGER(dimAttr)[0];
    int nn=INTEGER(dimAttr)[1];
    
    double *x=REAL(xS);

    if(!isNumeric(LS) || !isMatrix(LS)) {
	Rprintf("STOP: L not a numeric matrix\n");
	return R_NilValue;
    }

    dimAttr = getAttrib(LS, R_DimSymbol);
    
    if (n != ((int) INTEGER(dimAttr)[0]) )
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 L %d\n",n,((int) INTEGER(dimAttr)[0]));
	return R_NilValue;
    }

    int K = INTEGER(dimAttr)[1];

    double **L = R_Calloc(n, double *);
    L[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	L[i] = L[0]+ i*K;
	for(j=0; j < K; j++)
	{
	    L[i][j] = (double)(REAL(LS)[i+n*j]);
	}
    }

    if(!isNumeric(laplaS) || !isMatrix(laplaS)) {
	Rprintf("STOP: lapla not a numeric matrix\n");
	return R_NilValue;
    }
    
    dimAttr = getAttrib(laplaS, R_DimSymbol);
    
    if (nn !=  ((int) INTEGER(dimAttr)[0]))
    {
	Rprintf("STOP: Dim2 X: %d BUT Dim1 lapla %d\n",nn,((int) INTEGER(dimAttr)[0]));
	return R_NilValue;
    }
    
    if (K != ((int) INTEGER(dimAttr)[1]))
    {
	Rprintf("STOP: Dim2 L: %d BUT Dim2 lapla %d\n",K,((int) INTEGER(dimAttr)[1]));
	return R_NilValue;
    }


    double **lapla = R_Calloc(nn, double *);
    lapla[0] = R_Calloc((long) nn*K, double);
    for(i=0; i < nn; i++)
    {
	lapla[i] = lapla[0] + i*K;
	for(j=0; j < K; j++)
	{
	    lapla[i][j] = (double)(REAL(laplaS)[i+nn*j]);
	}
    }


 
    if(!isNumeric(PsiS) || isMatrix(PsiS) || isLogical(PsiS))
    {
	Rprintf("STOP: Psi not a numeric vector\n");
	return R_NilValue;
    }

    if (n != length(PsiS))
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 Psi %d\n",n,length(PsiS));
	return R_NilValue;

    }

    double *Psi = R_Calloc(n, double); 

    for(i=0; i < n; i++)
    {
	Psi[i] = (double)(REAL(PsiS)[i]);  
    }



    double alpha = (double)(REAL(alphaS)[0]);
    double eps = (double)(REAL(epsS)[0]);
    double eps1 = (double)(REAL(eps1S)[0]);
    double spl = (double)(REAL(splS)[0]);
    double spz = (double)(REAL(spzS)[0]);
    double scale = (double)(REAL(scaleS)[0]);
    double lap = (double)(REAL(lapS)[0]);


    int non_negative = (int)(INTEGER(non_negativeS)[0]);
    int cyc = (int)(INTEGER(cycS)[0]);
    int nL = (int)(INTEGER(nLS)[0]);
    int lL = (int)(INTEGER(lLS)[0]);
    int bL = (int)(INTEGER(bLS)[0]);


    double *XX = R_Calloc(n, double); 
    double *e_sx_n = R_Calloc(K, double); 
    double *e_ssxx_n = R_Calloc(K, double); 



    double **ichol = R_Calloc(K, double *);
    ichol[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	ichol[i] = ichol[0] + i*K;
    }

    double **sum2 = R_Calloc(K, double *);
    sum2[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	sum2[i] = sum2[0]+i*K;
    }

    double **tLPsiL = R_Calloc(K, double *);
    tLPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	tLPsiL[i] = tLPsiL[0]+i*K;
    }

    double **LPsiL = R_Calloc(K, double *);
    LPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	LPsiL[i] =  LPsiL[0]+i*K;
    }

    double **LPsi = R_Calloc(n, double *);
    LPsi[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	LPsi[i] =  LPsi[0]+i*K;
    }

    double **sum1 = R_Calloc(n, double *);
    sum1[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	sum1[i] =  sum1[0]+ i*K;
    }

  


    
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
	    s += x[i1+n*i2] * x[i1+n*i2];
// x[i1][i2]=x[i1+n*i2]
	XX[i1] = s*in;
	if (XX[i1]<eps) XX[i1]=eps;

    }



    for (i=0;i<cyc;i++) {


	for (i1 = 0; i1 < K; i1++)
	{
	    for (i2=0;i2<n;i2++)
		LPsi[i2][i1] =  L[i2][i1]/Psi[i2];
	}

	for (i1=0;i1<K;i1++)
	{
	  s = 0.0;
	  for (i2 = 0; i2 < n; i2++)
	    s += LPsi[i2][i1]*L[i2][i1];
	  LPsiL[i1][i1] = s;

	}

	for (i1=0;i1<K-1;i1++)
	{
	    for (i3=i1+1;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = LPsiL[i3][i1] = s;
	    }
	}


	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < n; i2++)
		sum1[i2][i1] = 0.0;
	    for (i2 = 0; i2 < K; i2++)
		sum2[i1][i2] = (i1==i2 ? eps : 0.0);
	}

	for (j=0;j<nn;j++)
	{

	    for (i1=0;i1<K;i1++)
	    {
		for (i2=0;i2<K;i2++)
		{
		    
		    if (i1==i2) {
			tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		    } else {
			tLPsiL[i1][i2] = LPsiL[i1][i2];
		    }
		}
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


	    for (i3 = 0; i3 < K; i3++) {
	      t=0.0;
	      for (i2 = 0; i2 < n; i2++){
		t+= LPsi[i2][i3]*x[i2+j*n];// x[i2][j]=x[i2+j*n]
	      }
	      e_ssxx_n[i3]=t;
	    }



	    for (i1=0;i1<K;i1++)
	    {
		t=0.0;
		for (i3 = 0; i3 < K; i3++)
		  t +=  ichol[i1][i3]*e_ssxx_n[i3];
		if ((t<eps1) && (non_negative>0))
		{
		    t=0.0;
		    e_sx_n[i1] = 0.0;
		} else {
		  e_sx_n[i1] = t;
		  for (i2 = 0; i2 < n; i2++)
		    sum1[i2][i1] += x[i2+j*n]*t;
		}
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


	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < n; i2++){
		s=0.0;
		for (i3 = 0; i3 < K; i3++)
		    s +=  sum1[i2][i3]*ichol[i3][i1];

		sgn = (s>0.0) ? 1.0 : ((s == 0.0) ? 0.0 : -1.0);
		if ((sgn>0)||(non_negative<=0))
		{
		    t = fabs(Psi[i2]*alpha*pow((eps1+fabs(s)),spl));

		    if (fabs(s)>t){
			L[i2][i1] = s - sgn* t;
		    }
		    else {
			L[i2][i1] = 0;
		    }
		} else {
		    L[i2][i1] = 0;
		}
	    }
	}
    


	if ((nL>0)&&(nL<K)&&(cyc>bL))
	{

	    for (i1=0;i1<n;i1++)
	    {
		for (i2 = 0; i2 < nL; i2++)
		{
		    e_ssxx_n[i2] = -1.0;
		}
		for (i2 = 0; i2 < K; i2++)
		{
		    s=fabs(L[i1][i2]);
		    if (s>e_ssxx_n[nL-1])
		    {
			i3=nL-1;
			while ((i3>0)&&(s>e_ssxx_n[i3-1])) {
			    e_ssxx_n[i3]=e_ssxx_n[i3-1];
			    i3--;
			}
			e_ssxx_n[i3]=s;
		    }
		}
		s=e_ssxx_n[nL-1];
		for (i2 = 0; i2 < K; i2++)
		{
		    if(s>fabs(L[i1][i2]))
		    {
			L[i1][i2]=0.0;
		    }
		}
	

	    }
	}

	if ((lL>0)&&(lL<n)&&(cyc>bL))
	{

	  for (i2 = 0; i2 < K; i2++)
	    {
		for (i1=0;i1<lL;i1++)
		{
		    Psi[i1] = -1.0;
		}
		for (i1=0;i1<n;i1++)
		{
		    s=fabs(L[i1][i2]);
		    if (s>Psi[lL-1])
		    {
			i3=lL-1;
			while ((i3>0)&&(s>Psi[i3-1])) {
			    Psi[i3]=Psi[i3-1];
			    i3--;
			}
			Psi[i3]=s;
		    }
		}
		s=Psi[lL-1];
		for (i1=0;i1<n;i1++)
		{
		    if(s>fabs(L[i1][i2]))
		    {
			L[i1][i2]=0.0;
		    }
		}
	

	    }
	}


	for (i1=0;i1<n;i1++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < K; i2++) {
		s += L[i1][i2]*sum1[i1][i2];
	    }
	    Psi[i1] = XX[i1] - in*s;
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
		s=sqrt(s)+eps1;
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
	    Rprintf("Cycle: %d\r", i);
	    R_CheckUserInterrupt();
	}

    }
    
    Rprintf("Cycle: %d\n", (cyc-1));

    R_Free (sum1[0]);
    R_Free (sum1 );

    R_Free (sum2[0]);
    R_Free (sum2 );

    R_Free (XX );
    R_Free (e_sx_n );


    for (i1=0;i1<K;i1++)
    {
	for (i2 = 0; i2 < n; i2++)
	    LPsi[i2][i1] =  L[i2][i1]/Psi[i2];
    }


	for (i1=0;i1<K;i1++)
	{
	  s = 0.0;
	  for (i2 = 0; i2 < n; i2++)
	    s += LPsi[i2][i1]*L[i2][i1];
	  LPsiL[i1][i1] = s;

	}

	for (i1=0;i1<K-1;i1++)
	{
	    for (i3=i1+1;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = LPsiL[i3][i1] = s;
	    }
	}




    SEXP E_SX_n;
    PROTECT(E_SX_n = allocMatrix(REALSXP, K, nn));

    for (j=0;j<nn;j++)
    {

	for (i1=0;i1<K;i1++)
	{
	    for (i2=0;i2<K;i2++)
	    {
		
		if (i1==i2) {
		    tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		} else {
		    tLPsiL[i1][i2] = LPsiL[i1][i2];
		}
	    }
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


	    for (i3 = 0; i3 < K; i3++) {
	      t=0.0;
	      for (i2 = 0; i2 < n; i2++){
		t+= LPsi[i2][i3]*x[i2+j*n];// x[i2][j]=x[i2+j*n]
	      }
	      e_ssxx_n[i3]=t;
	    }

	    for (i1=0;i1<K;i1++)
	    {
		t=0.0;
		for (i3 = 0; i3 < K; i3++)
		  t +=  ichol[i1][i3]*e_ssxx_n[i3];
		if ((t<eps1) && (non_negative>0))
		{
		    t=0.0;
		}
		REAL(E_SX_n)[i1 + K*j] = (double) t;
	    }



    }


    R_Free (e_ssxx_n );

    R_Free (ichol[0]);
    R_Free (ichol );

    R_Free (tLPsiL[0]);
    R_Free (tLPsiL );

    R_Free (LPsiL[0]);
    R_Free (LPsiL );

    R_Free (LPsi[0]);
    R_Free (LPsi);


    SEXP L_n;
    PROTECT(L_n = allocMatrix(REALSXP, n, K));


    for(i = 0; i < n; i++)
	for(j = 0; j < K; j++)
	    REAL(L_n)[i + n*j] = (double) L[i][j];


    R_Free (L[0]);
    R_Free( L );


    SEXP Psi_n;
    PROTECT(Psi_n = allocVector(REALSXP, n));

    for(i = 0; i < n; i++)
	REAL(Psi_n)[i] = (double) Psi[i];

    R_Free ( Psi );

    SEXP lapla_n;
    PROTECT(lapla_n = allocMatrix(REALSXP, nn,K));
    for(i = 0; i < nn; i++)
	for(j = 0; j < K; j++)
	    REAL(lapla_n)[i + nn*j] = (double) lapla[i][j];



    R_Free (lapla[0]);
    R_Free( lapla );



    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesRET, 0, mkChar("L"));
    SET_STRING_ELT(namesRET, 1, mkChar("E_SX_n"));
    SET_STRING_ELT(namesRET, 2, mkChar("Psi"));
    SET_STRING_ELT(namesRET, 3, mkChar("lapla"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(RET, 0, L_n);
    SET_VECTOR_ELT(RET, 1, E_SX_n);
    SET_VECTOR_ELT(RET, 2, Psi_n);
    SET_VECTOR_ELT(RET, 3, lapla_n);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(6);
    return(RET);

}







SEXP fabics(SEXP xS, SEXP PsiS,SEXP LS,SEXP laplaS,SEXP cycS, SEXP alphaS,SEXP epsS,SEXP spzS,SEXP lapS,SEXP nLS,SEXP lLS,SEXP bLS,SEXP non_negativeS) {

    int i,j,i1,i2,i3,zz,ende,h1;

    double in,s,sgn,a,b,c,t,alphap,k2,seps,eps1=1e-10;


    if(!isNumeric(xS) || !isMatrix(xS)) {
	Rprintf("STOP: X not a numeric matrix\n");
	return R_NilValue;
    }

    SEXP dimAttr;

    dimAttr=getAttrib(xS, R_DimSymbol);
    int n=INTEGER(dimAttr)[0];
    int nn=INTEGER(dimAttr)[1];
    
    
    double *x=REAL(xS);


    if(!isNumeric(LS) || !isMatrix(LS)) {
	Rprintf("STOP: L not a numeric matrix\n");
	return R_NilValue;
    }

    dimAttr = getAttrib(LS, R_DimSymbol);
    
    if (n != ((int) INTEGER(dimAttr)[0]) )
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 L %d\n",n,((int) INTEGER(dimAttr)[0]));
	return R_NilValue;
    }

    int K = INTEGER(dimAttr)[1];

    double **L = R_Calloc(n, double *);
    L[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	L[i] =  L[0] + i*K;
	for(j=0; j < K; j++)
	{
	    L[i][j] = (double)(REAL(LS)[i+n*j]);
	}
    }


    if(!isNumeric(laplaS) || !isMatrix(laplaS)) {
	Rprintf("STOP: lapla not a numeric matrix\n");
	return R_NilValue;
    }

    dimAttr = getAttrib(laplaS, R_DimSymbol);
    
    if (nn !=  ((int) INTEGER(dimAttr)[0]))
    {
	Rprintf("STOP: Dim2 X: %d BUT Dim1 lapla %d\n",nn,((int) INTEGER(dimAttr)[0]));
	return R_NilValue;
    }
    
    if (K != ((int) INTEGER(dimAttr)[1]))
    {
	Rprintf("STOP: Dim2 L: %d BUT Dim2 lapla %d\n",K,((int) INTEGER(dimAttr)[1]));
	return R_NilValue;
    }


    double **lapla = R_Calloc(nn, double *);
    lapla[0] = R_Calloc((long) nn*K, double);
    for(i=0; i < nn; i++)
    {
	lapla[i] = lapla[0] + i*K;
	for(j=0; j < K; j++)
	{
	    lapla[i][j] = (double)(REAL(laplaS)[i+nn*j]);
	}
    }

 
 
    if(!isNumeric(PsiS) || isMatrix(PsiS) || isLogical(PsiS))
    {
	Rprintf("STOP: Psi not a numeric vector\n");
	return R_NilValue;
    }

    if (n != length(PsiS))
    {
	Rprintf("STOP: Dim1 X: %d BUT Dim1 Psi %d\n",n,length(PsiS));
	return R_NilValue;

    }

    double *Psi = R_Calloc(n, double); 
    
    for(i=0; i < n; i++)
    {
	Psi[i] = (double)(REAL(PsiS)[i]);  
    }
    

    double alpha = (double)(REAL(alphaS)[0]);
    double eps = (double)(REAL(epsS)[0]);
    double spz = (double)(REAL(spzS)[0]);
    double lap = (double)(REAL(lapS)[0]);

    int non_negative = (int)(INTEGER(non_negativeS)[0]);
    int cyc = (int)(INTEGER(cycS)[0]);
    int nL = (int)(INTEGER(nLS)[0]);
    int lL = (int)(INTEGER(lLS)[0]);
    int bL = (int)(INTEGER(bLS)[0]);

    double *XX = R_Calloc(n, double); 
    double *e_sx_n = R_Calloc(K, double); 
    double *e_ssxx_n = R_Calloc(K, double); 
    double *v = R_Calloc(n, double); 
    double *w = R_Calloc(n, double); 
    double *r_es_x = R_Calloc(n, double); 

    int *zeros = R_Calloc(n, int); 


    double **ichol = R_Calloc(K, double *);
    ichol[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	ichol[i] =  ichol[0] + i*K;
    }

    double **sum2 = R_Calloc(K, double *);
    sum2[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	sum2[i] = sum2[0] + i*K;
    }

    double **tLPsiL = R_Calloc(K, double *);
    tLPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	tLPsiL[i] =  tLPsiL[0] + i*K;
    }

    double **LPsiL = R_Calloc(K, double *);
    LPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	LPsiL[i] =  LPsiL[0] + i*K;
    }

    double **LPsi = R_Calloc(n, double *);
    LPsi[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	LPsi[i] = LPsi[0] + i*K;
    }

    double **sum1 = R_Calloc(n, double *);
    sum1[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	sum1[i] = sum1[0] + i*K;
    }

  



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
	    s += x[i1+n*i2] * x[i1+n*i2];
	XX[i1] = s*in;
	if (XX[i1]<eps) XX[i1]=eps;

    }


 
    for (i=0;i<cyc;i++) {


	for (i1 = 0; i1 < K; i1++)
	{
	    for (i2=0;i2<n;i2++)
		LPsi[i2][i1] =  L[i2][i1]/Psi[i2];
	}

	for (i1=0;i1<K;i1++)
	{
	  s = 0.0;
	  for (i2 = 0; i2 < n; i2++)
	    s += LPsi[i2][i1]*L[i2][i1];
	  LPsiL[i1][i1] = s;

	}

	for (i1=0;i1<K-1;i1++)
	{
	    for (i3=i1+1;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = LPsiL[i3][i1] = s;
	    }
	}



	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < n; i2++)
		sum1[i2][i1] = 0.0;
	    for (i2 = 0; i2 < K; i2++)
		sum2[i1][i2] = (i1==i2 ? eps : 0.0);
	}


	for (j=0;j<nn;j++)
	{

	    for (i1=0;i1<K;i1++)
	    {
		for (i2=0;i2<K;i2++)
		{
		    
		    if (i1==i2) {
			tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		    } else {
			tLPsiL[i1][i2] = LPsiL[i1][i2];
		    }
		}
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



	    for (i3 = 0; i3 < K; i3++) {
	      t=0.0;
	      for (i2 = 0; i2 < n; i2++){
		t+= LPsi[i2][i3]*x[i2+j*n];// x[i2][j]=x[i2+j*n]
	      }
	      e_ssxx_n[i3]=t;
	    }

	    for (i1=0;i1<K;i1++)
	    {
		t=0.0;
		for (i3 = 0; i3 < K; i3++)
		  t +=  ichol[i1][i3]*e_ssxx_n[i3];
		if ((t<eps1) && (non_negative>0))
		{
		    t=0.0;
		    e_sx_n[i1] = 0.0;
		} else {
		  e_sx_n[i1] = t;
		  for (i2 = 0; i2 < n; i2++)
		    sum1[i2][i1] += x[i2+j*n]*t;
		}
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



	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < n; i2++){
		s=0.0;
		for (i3 = 0; i3 < K; i3++)
		    s +=  sum1[i2][i3]*ichol[i3][i1];
		L[i2][i1] = s;
	    }
	}
	
	

//----------

        for (i3 = 0; i3 < K; i3++) {

	    s = 0.0;
	    for (i1=0;i1<n;i1++)
	    {
		r_es_x[i1] = (L[i1][i3]<0.0 ? -1.0 : 1.0);
		sgn=fabs(L[i1][i3]);
		w[i1] = sgn;
		s += sgn;
	    }


	    k2=1.0;

	    zz=0;
	    s=(alpha-s)/(1.0*n);
	    for (i1=0;i1<n;i1++)
	    {
		sgn = w[i1] + s;
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
		if (non_negative<=0)
		{
		    L[i1][i3] = r_es_x[i1]*v[i1];
		} else {
		    L[i1][i3] = v[i1];
		}
//		Rprintf("%d %lf %lf\n",i1,L[i1][i3],r_es_x[i1]*v[i1]);
	    }
	}

//----------


	if ((nL>0)&&(nL<K)&&(cyc>bL))
	{

	    for (i1=0;i1<n;i1++)
	    {
		for (i2 = 0; i2 < nL; i2++)
		{
		    e_ssxx_n[i2] = 0.0;
		}
		for (i2 = 0; i2 < K; i2++)
		{
		    s=fabs(L[i1][i2]);
		    if (s>e_ssxx_n[nL-1])
		    {
			i3=nL-1;
			while ((i3>0)&&(s>e_ssxx_n[i3-1])) {
			    e_ssxx_n[i3]=e_ssxx_n[i3-1];
			    i3--;
			}
			e_ssxx_n[i3]=s;
		    }
		}
		s=e_ssxx_n[nL-1];
		for (i2 = 0; i2 < K; i2++)
		{
		    if(s>fabs(L[i1][i2]))
		    {
			L[i1][i2]=0.0;
		    }
		}
	

	    }
	}




	if ((lL>0)&&(lL<n)&&(cyc>bL))
	{

	  for (i2 = 0; i2 < K; i2++)
	    {
		for (i1=0;i1<lL;i1++)
		{
		    Psi[i1] = -1.0;
		}
		for (i1=0;i1<n;i1++)
		{
		    s=fabs(L[i1][i2]);
		    if (s>Psi[lL-1])
		    {
			i3=lL-1;
			while ((i3>0)&&(s>Psi[i3-1])) {
			    Psi[i3]=Psi[i3-1];
			    i3--;
			}
			Psi[i3]=s;
		    }
		}
		s=Psi[lL-1];
		for (i1=0;i1<n;i1++)
		{
		    if(s>fabs(L[i1][i2]))
		    {
			L[i1][i2]=0.0;
		    }
		}
	

	    }
	}




	for (i1=0;i1<n;i1++)
	{
	    s = 0.0;
	    for (i2 = 0; i2 < K; i2++) {
		s += L[i1][i2]*sum1[i1][i2];
	    }
	    Psi[i1] = XX[i1] - in*s;
	    if (Psi[i1]<eps)
	    {
		Psi[i1] = eps;
	    }
	}





	if (i%20==0) {
	    Rprintf("Cycle: %d\r", i);
	    R_CheckUserInterrupt();
	}
    }

    Rprintf("Cycle: %d\n", (cyc-1));

    R_Free (sum1[0]);
    R_Free (sum1 );

    R_Free (sum2[0]);
    R_Free (sum2 );

    R_Free (XX );
    R_Free (e_sx_n );
    R_Free (v );
    R_Free (w );
    R_Free (r_es_x );
    R_Free (zeros );

    for (i1=0;i1<K;i1++)
    {
	for (i2 = 0; i2 < n; i2++)
	    LPsi[i2][i1] =  L[i2][i1]/Psi[i2];
    }


	for (i1=0;i1<K;i1++)
	{
	  s = 0.0;
	  for (i2 = 0; i2 < n; i2++)
	    s += LPsi[i2][i1]*L[i2][i1];
	  LPsiL[i1][i1] = s;

	}

	for (i1=0;i1<K-1;i1++)
	{
	    for (i3=i1+1;i3<K;i3++)
	    {
		s = 0.0;
		for (i2 = 0; i2 < n; i2++)
		    s += LPsi[i2][i1]*L[i2][i3];
		LPsiL[i1][i3] = LPsiL[i3][i1] = s;
	    }
	}




    SEXP E_SX_n;
    PROTECT(E_SX_n = allocMatrix(REALSXP, K, nn));

    for (j=0;j<nn;j++)
    {

	for (i1=0;i1<K;i1++)
	{
	    for (i2=0;i2<K;i2++)
	    {
		
		if (i1==i2) {
		    tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		} else {
		    tLPsiL[i1][i2] = LPsiL[i1][i2];
		}
	    }
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




	    for (i3 = 0; i3 < K; i3++) {
	      t=0.0;
	      for (i2 = 0; i2 < n; i2++){
		t+= LPsi[i2][i3]*x[i2+j*n];// x[i2][j]=x[i2+j*n]
	      }
	      e_ssxx_n[i3]=t;
	    }

	    for (i1=0;i1<K;i1++)
	    {
		t=0.0;
		for (i3 = 0; i3 < K; i3++)
		  t +=  ichol[i1][i3]*e_ssxx_n[i3];
		if ((t<eps1) && (non_negative>0))
		{
		    t=0.0;
		}
		REAL(E_SX_n)[i1 + K*j] = (double) t;
	    }


    }

    R_Free (e_ssxx_n );

    R_Free (ichol[0]);
    R_Free (ichol );

    R_Free (tLPsiL[0]);
    R_Free (tLPsiL );

    R_Free (LPsiL[0]);
    R_Free (LPsiL );

    R_Free (LPsi[0]);
    R_Free (LPsi);



    SEXP L_n;
    PROTECT(L_n = allocMatrix(REALSXP, n, K));


    for(i = 0; i < n; i++)
	for(j = 0; j < K; j++)
	    REAL(L_n)[i + n*j] = (double) L[i][j];


    R_Free (L[0]);
    R_Free( L );


    SEXP Psi_n;
    PROTECT(Psi_n = allocVector(REALSXP, n));

    for(i = 0; i < n; i++)
	REAL(Psi_n)[i] = (double) Psi[i];

    R_Free ( Psi );

    SEXP lapla_n;
    PROTECT(lapla_n = allocMatrix(REALSXP, nn,K));
    for(i = 0; i < nn; i++)
	for(j = 0; j < K; j++)
	    REAL(lapla_n)[i + nn*j] = (double) lapla[i][j];



    R_Free (lapla[0]);
    R_Free( lapla );






    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesRET, 0, mkChar("L"));
    SET_STRING_ELT(namesRET, 1, mkChar("E_SX_n"));
    SET_STRING_ELT(namesRET, 2, mkChar("Psi"));
    SET_STRING_ELT(namesRET, 3, mkChar("lapla"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(RET, 0, L_n);
    SET_VECTOR_ELT(RET, 1, E_SX_n);
    SET_VECTOR_ELT(RET, 2, Psi_n);
    SET_VECTOR_ELT(RET, 3, lapla_n);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(6);
    return(RET);

}





SEXP spfabic(SEXP file_nameS, SEXP KS, SEXP alphaS, SEXP cycS, SEXP splS,SEXP spzS, SEXP non_negativeS,SEXP randomS, SEXP write_fileS, SEXP init_psiS, SEXP init_laplaS, SEXP normS,SEXP scaleS,SEXP lapS,SEXP nLS, SEXP lLS,SEXP bLS,SEXP epsS,SEXP eps1S,SEXP samplesS,SEXP initLS,SEXP iterS,SEXP quantS,SEXP lowerBS, SEXP upperBS) {


    FILE *pFile;


    char sst[200]; 
    char iterc [4];

    int hpp,ig,jg,samp,inLL;

    int  i,j,i1,i2,i3,i4,n,nn,ite,nquant,la1,la2;
    
    double fs;

    double in,s,sgn,t=0.0;




    const char *file_name=CHAR(STRING_ELT(file_nameS,0));

  
    double init_lapla = (double)(REAL(init_laplaS)[0]);
    double init_psi = (double)(REAL(init_psiS)[0]);
    double random = (double)(REAL(randomS)[0]);
    double alpha = (double)(REAL(alphaS)[0]);
    double eps = (double)(REAL(epsS)[0]);
    double eps1 = (double)(REAL(eps1S)[0]);
    double spl = (double)(REAL(splS)[0]);
    double spz = (double)(REAL(spzS)[0]);
    double scale = (double)(REAL(scaleS)[0]);
    double lap = (double)(REAL(lapS)[0]);
    double quant = (double)(REAL(quantS)[0]);
    double lowerB = (double)(REAL(lowerBS)[0]);
    double upperB = (double)(REAL(upperBS)[0]);


    int write_file =  (int)(INTEGER(write_fileS)[0]);
    int non_negative =  (int)(INTEGER(non_negativeS)[0]);
    int norm =  (int)(INTEGER(normS)[0]);
    int cyc =  (int)(INTEGER(cycS)[0]);
    int K =  (int)(INTEGER(KS)[0]);
    int nL =  (int)(INTEGER(nLS)[0]);
    int lL = (int)(INTEGER(lLS)[0]);
    int bL = (int)(INTEGER(bLS)[0]);
    int iter = (int)(INTEGER(iterS)[0]);

    int *xa;
    int **xind;
    double **xval;

    int nsamp = length(samplesS);

    int *samples =  INTEGER(samplesS);
    if (samples[0]>0) {
      samp = 0;
      
    } else {
      samp=-1;
    }

    GetRNGstate();

    n=0;
    nn=0;

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,".txt");

    pFile = fopen (sst,"r");

    if (!(pFile>0)) {
	Rprintf("File >%s< not found! Stop.\n", sst);
	return R_NilValue;
    }

    fscanf(pFile,"%d\n",&nn);  

    if (!(nn>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    fscanf(pFile,"%d\n",&n);  

    if (!(n>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    if (samp<0) {

      xa = (int *) R_Calloc(nn, int); 
      xind = (int **) R_Calloc(nn, int *);
      xval = (double **) R_Calloc(nn, double *);
    
      for(i = 0; i < nn; i ++)
	{
	  fscanf(pFile,"%d\n",&ig); 
	  xa[i]=ig;
	  xind[i] = R_Calloc((ig+1), int);
	  xval[i] = R_Calloc((ig+1), double);
	  for(j = 0; j <  ig; j ++) {
	    fscanf(pFile,"%d",&hpp);
	    xind[i][j]=hpp;
	  }
	  fscanf(pFile,"\n");
	  for(j = 0; j < ig; j ++) {
	    fscanf(pFile,"%lf",&fs);
	    xval[i][j] = fs;
	  }
	  fscanf(pFile,"\n");
	}

    } else {
      xa = (int *) R_Calloc(nsamp, int); 
      xind = (int **) R_Calloc(nsamp, int *);
      xval = (double **) R_Calloc(nsamp, double *);

      for(i = 0; i < nn; i ++)
	{
	  if ((samples[samp]-1)>i)
	    {
	      fscanf(pFile,"%d\n",&ig);
	      for(j = 0; j <  ig; j ++) {
		fscanf(pFile,"%d",&hpp);
	      }
	      fscanf(pFile,"\n");
	      for(j = 0; j < ig; j ++) {
		fscanf(pFile,"%lf",&fs);
	      }
	      fscanf(pFile,"\n");
	      
	    } else {
	    fscanf(pFile,"%d\n",&ig); 
	    xa[samp]=ig;
	    xind[samp] = R_Calloc((ig+1), int);
	    xval[samp] = R_Calloc((ig+1), double);
	    for(j = 0; j <  ig; j ++) {
	      fscanf(pFile,"%d",&hpp);
	      xind[samp][j]=hpp;
	    }
	    fscanf(pFile,"\n");
	    for(j = 0; j < ig; j ++) {
	      fscanf(pFile,"%lf",&fs);
	      xval[samp][j] = fs;
	    }
	    fscanf(pFile,"\n");
	    samp++;
	    if (samp == nsamp) break;
	  }

	}
      if (samp!=nsamp)
	{
	  Rprintf("Only %d of %d samples found! Some sample numbers are too large. Continue.\n", samp,nsamp);
	}
      nn=samp;
      Rprintf("Using %d samples!\n",samp);
    }
    fclose (pFile);

    int ninitL = length(initLS);
    if (ninitL > K) ninitL = K;
    int *initL =  INTEGER(initLS);
    if ((initL[0]>0)&&(initL[0]<=nn)) {
      inLL = 0;
      
    } else {
      if (initL[0]==-2) {
	inLL=-2;
      } else {
	inLL=-1;
      }
    }

    int *La = R_Calloc(K, int); 
    int **Lind = R_Calloc(K, int *);
    Lind[0] = R_Calloc((long) K*n, int);
    for(i=0; i < K; i++)
    {
	Lind[i] = Lind[0] + i*n;
    }
    double **Lval = R_Calloc(K, double *);
    Lval[0] = R_Calloc((long) K*n, double);
    for(i=0; i < K; i++)
    {
	Lval[i] =  Lval[0] + i*n;
    }






    double *LLval = R_Calloc(K, double); 
    int *LLind = R_Calloc(K, int); 
    int *ig_vec = R_Calloc(K, int); 

    double *Psi = R_Calloc(n, double); 
    double *XX = R_Calloc(n, double); 
    double *e_sx_n = R_Calloc(K, double); 
    double *e_ssxx_n = R_Calloc(K, double); 



    int *LPsia = R_Calloc(K, int); 
    int **LPsiind = R_Calloc(K, int *);
    LPsiind[0] = R_Calloc((long) K*n, int);
    for(i=0; i < K; i++)
    {
	LPsiind[i] = LPsiind[0] + i*n;
    }

    double **LPsival = R_Calloc(K, double *);
    LPsival[0] = R_Calloc((long) K*n, double);
    for(i=0; i < K; i++)
    {
	LPsival[i] = LPsival[0] + i*n;
    }

    double **lapla = R_Calloc(nn, double *);
    lapla[0] = R_Calloc((long) nn*K, double);
    for(i=0; i < nn; i++)
    {
	lapla[i] = lapla[0] + i*K;
    }


    double **ichol = R_Calloc(K, double *);
    ichol[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	ichol[i] = ichol[0] + i*K;
    }

    double **sum2 = R_Calloc(K, double *);
    sum2[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	sum2[i] =  sum2[0]+ i*K;
    }

    double **tLPsiL = R_Calloc(K, double *);
    tLPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	tLPsiL[i] =  tLPsiL[0] + i*K;
    }

    double **LPsiL = R_Calloc(K, double *);
    LPsiL[0] = R_Calloc(K*K, double);
    for(i=0; i <K ; i++)
    {
	LPsiL[i] = LPsiL[0] + i*K;
    }

    double **LPsi = R_Calloc(n, double *);
    LPsi[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	LPsi[i] = LPsi[0]+i*K;
    }

    double **sum1 = R_Calloc(n, double *);
    sum1[0] = R_Calloc((long) n*K, double);
    for(i=0; i < n; i++)
    {
	sum1[i] = sum1[0] + i*K;
    }

    double **E_SX_tmp = R_Calloc(nn, double *);
    E_SX_tmp[0] = R_Calloc((long) nn*K*iter, double);
    for(i=0; i < nn; i++)
    {
	E_SX_tmp[i] = E_SX_tmp[0] + i*K*iter;
    }


    nquant = (int) floor(quant*n);
    spl = -spl;
    spz = -spz;
    in = 1.0/nn;

    if (lap<eps)
    {
	lap = eps;
    }

    if (iter<1) iter = 1;

    SEXP L_n;
    PROTECT(L_n = allocMatrix(REALSXP, n, iter*K));

    SEXP Psi_n;
    PROTECT(Psi_n = allocVector(REALSXP, (long) iter*n));

    SEXP lapla_n;
    PROTECT(lapla_n = allocMatrix(REALSXP, nn,iter*K));

    SEXP E_SX_n;
    PROTECT(E_SX_n = allocMatrix(REALSXP, iter*K, nn));


 
    for (i1=0;i1<n;i1++)
      {
	XX[i1] = 0.0;
      }	
    
   for (i2 = 0; i2 < nn; i2++) {
	for (ig=0; ig< xa[i2];ig++) {
	    XX[xind[i2][ig]] += xval[i2][ig];
	}
    }


      for (i2 = 0; i2 < n; i2++) {
	if ((XX[i2]>lowerB)&&(XX[i2]<upperB)) {
	  LPsiind[0][i2] = 0;
	} else {
	  LPsiind[0][i2] = 1;
	}
      }


      for (i2 = 0; i2 < nn; i2++) {
	for (ig=0,jg=0; (ig+jg) < xa[i2];) {
	  if ( LPsiind[0][xind[i2][ig+jg]]==0 )
	    {
	      if (jg>0) {
		xval[i2][ig] = xval[i2][ig+jg];
		xind[i2][ig] = xind[i2][ig+jg];
	      }
	      ig++;
	    } else 
	    {
	      jg++;
	    }
	}
	xa[i2]-=jg;
      }





   for (ite=0;ite<iter;ite++) {



  
    for (i1=0;i1<n;i1++)
      {
	XX[i1] = 0.0;
      }	
    
   for (i2 = 0; i2 < nn; i2++) {
	for (ig=0; ig< xa[i2];ig++) {
	    s=xval[i2][ig];
	    XX[xind[i2][ig]] += s*s;
	}
    }
  
    if (inLL==-2) {
      t = 0.0;
      for (i1=0;i1<n;i1++)
	{
	  if (XX[i1] > t) {
	    t=XX[i1];
	  }
	}	
      t=sqrt(t)+eps;
    }





   if (inLL < 0) {
     if (non_negative>0) {
       
       for(i = 0; i < K; i ++)
	 {
	   for(i1=0,j = 0; j < n; j ++) {
	     if (XX[j]>eps) {
	       
	       Lind[i][i1]= j;
	       // Lval[i][j] = random*(rand()%100001)/100000.0;
	       if (inLL==-2) {
		 Lval[i][i1] = (sqrt(XX[j])/t)*random*fabs(norm_rand());
	       } else {
		 Lval[i][i1] = random*fabs(norm_rand());
	       }
	       i1++;
	     }
	   }
	   La[i]=i1;
	   
	 }
     } else {
       for(i = 0; i < K; i ++)
	 {
	   for(i1=0,j = 0; j < n; j ++) {
	     if (XX[j]>eps) {
	       Lind[i][i1]= j;
//		if (rand()%2>0) {
//		    s = 1.0;
//		}  else {
//		    s= -1.0;
//		}
//		Lval[i][j] = random*s*(rand()%100001)/100000.0;
	       if (inLL==-2) {
		 Lval[i][i1] = (sqrt(XX[j])/t)*random*norm_rand();
	       } else {
		 Lval[i][i1] = random*norm_rand();
	       }
	       i1++;
	     }
	   }
	   La[i]=i1;
	   
	 }
     }
   } else {
     
     if (non_negative>0) {
       
       for(i = 0; i < ninitL; i ++)
	 {
	   i2= initL[i]-1;
	   if ((i2>nn)||(i2<0))
	     {
	       i2=i;
	     }
	   La[i]=xa[i2];
	   for (ig=0; ig< xa[i2];ig++) {
	     Lval[i][ig]=xval[i2][ig];
	     Lind[i][ig]=xind[i2][ig];
	   }
	 }
       for(i = ninitL; i < K; i ++)
	 {
	   La[i]=n;
	   for(j = 0; j < n; j ++) {
	     Lind[i][j]= j;
	     // Lval[i][j] = random*(rand()%100001)/100000.0;
	     Lval[i][j] = random*fabs(norm_rand());
	   }
	 }
       
     } else {
       for(i = 0; i < ninitL; i ++)
	 {
	   i2= initL[i]-1;
	   if ((i2>nn)||(i2<0))
	     {
	       i2=i;
	     }
	   La[i]=xa[i2];
	   for (ig=0; ig< xa[i2];ig++) {
	     Lval[i][ig]=xval[i2][ig];
	     Lind[i][ig]=xind[i2][ig];
	   }
	   
	 }
       for(i = ninitL; i < K; i ++)
	 {
	   La[i]=n;
	   for(j = 0; j < n; j ++) {
	     Lind[i][j]= j;
	     //		if (rand()%2>0) {
	     //		    s = 1.0;
	     //		}  else {
	     //		    s= -1.0;
	     //		}
	     //		Lval[i][j] = random*s*(rand()%100001)/100000.0;
	     Lval[i][j] = random*norm_rand();
	   }
	  }
       
     }
     
   }
   



 
  
    for (i1=0;i1<n;i1++) {
	s = XX[i1] * in;
	if (s<eps) s=eps;
	Psi[i1] = sqrt(s);
	if (norm>0) {
	 XX[i1] = 1.0;
	} else {
	 XX[i1] = s;
	}
    }




   if (norm>0) {
	
	for (i2 = 0; i2 < nn; i2++) {
	    for (ig=0; ig< xa[i2];ig++) {
		xval[i2][ig]/=Psi[xind[i2][ig]];
	    }
	}
    }


    for (i1=0;i1<n;i1++) {
	
	if (norm>0) {
	    Psi[i1]= init_psi;
	} else {
	    Psi[i1]= init_psi*XX[i1];
	}

    }

    for (i1=0;i1<K;i1++) {
	for (j=0;j<nn;j++)
	    lapla[j][i1] = init_lapla;
    }


 
    for (i=0;i<cyc;i++) {


	for (i2 = 0; i2 < K; i2++) {
	    LPsia[i2]=La[i2];
	    for (ig=0; ig< La[i2];ig++) {
		LPsiind[i2][ig] = Lind[i2][ig];
		LPsival[i2][ig] = Lval[i2][ig]/Psi[Lind[i2][ig]];
	    }
	}



 
	for (i2=0;i2<K;i2++)
	{
	  s = 0.0;  
	  for (ig=0;ig < La[i2];ig++) {
	    s += Lval[i2][ig] * LPsival[i2][ig];
	  }
	  LPsiL[i2][i2] = s;
	}

	for (i2=0;i2<K-1;i2++)
	{
	    for (i3=i2+1;i3<K;i3++)
	    {
	      s = 0.0; 
	      ig = 0; 
	      jg = 0;
	      la1= La[i3];
	      la2= LPsia[i2];
	      while ((ig < la1) && (jg < la2)) {
		
		if (Lind[i3][ig] < LPsiind[i2][jg]) {
		  ig++; 
		  while ((ig < la1) && (Lind[i3][ig] < LPsiind[i2][jg])) {
		    ig++; 
		  }
		  if (ig >= la1) break;
		}
		
		if (Lind[i3][ig] > LPsiind[i2][jg]) {
		  jg++; 
		  while ((jg < la2) && (Lind[i3][ig] > LPsiind[i2][jg])) {
		    jg++; 
		  }
		  if (jg >= la2) break;
		  if (Lind[i3][ig] == LPsiind[i2][jg]) {
		    s += Lval[i3][ig] * LPsival[i2][jg]; 
		    ig++;
		    jg++;
		  }
		} else {
		  s += Lval[i3][ig] * LPsival[i2][jg]; 
		  ig++;
		  jg++;
		}


	      } 

	      LPsiL[i2][i3] = LPsiL[i3][i2] = s;

	    }
	}
	    


	for (i1=0;i1<K;i1++)
	{
	    for (i2 = 0; i2 < n; i2++)
		sum1[i2][i1] = 0.0;
 	    for (i2 = 0; i2 < K; i2++)
		sum2[i1][i2] = (i1==i2 ? eps : 0.0);
	}




	for (j=0;j<nn;j++)
	{
	  la1=xa[j];

	    for (i1=0;i1<K;i1++)
	    {
		for (i2=0;i2<K;i2++)
		{
		    
		    if (i1==i2) {
			tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		    } else {
			tLPsiL[i1][i2] = LPsiL[i1][i2];
		    }
		}
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


 
	    for (i3 = 0; i3 < K; i3++)
	      {
		s=0.0;
		ig = 0; 
		jg = 0;
		la2= LPsia[i3];
		while ((ig < la1) && (jg < la2)) {

		  if (xind[j][ig] < LPsiind[i3][jg]) {
		    ig++;
		    while ((ig < la1) && (xind[j][ig] < LPsiind[i3][jg])) {
		      ig++; 
		    }
		    if (ig >= la1) break;
		  }
		  if (xind[j][ig] > LPsiind[i3][jg]) {
		    jg++;
		    while ((jg < la2) && (xind[j][ig] > LPsiind[i3][jg])) {
		      jg++;
		    }
		    if (jg >= la2) break;
		    if (xind[j][ig] == LPsiind[i3][jg]) {
		      s += xval[j][ig] * LPsival[i3][jg]; 
		      ig++;
		      jg++;
		    } 
		  } else  {
		    s += xval[j][ig] * LPsival[i3][jg]; 
		    ig++;
		    jg++;
		  } 
		}
		e_ssxx_n[i3]= s;
	      }


	    for (i1=0;i1<K;i1++)
	    {
		
		t=0.0;

		for (i3 = 0; i3 < K; i3++)
		  {
		    t +=  e_ssxx_n[i3]*ichol[i1][i3];
		  }
		if ((t<eps1) && (non_negative>0))
		{
		    t=0.0;
		    e_sx_n[i1] = 0.0;
		} else {
		  e_sx_n[i1] = t;
		  for (ig=0; ig< xa[j];ig++) {
		    sum1[xind[j][ig]][i1] += xval[j][ig]*t;
		  }
		}
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


 	if ((nL<=0)||(nL>=K)||(cyc<=bL))
	{
	    for (i1=0;i1<K;i1++)
	    {
		ig=0;
		for (i2 = 0; i2 < n; i2++){
		    s=0.0;
		    for (i3 = 0; i3 < K; i3++)
			s +=  sum1[i2][i3]*ichol[i3][i1];
		    
		    sgn = (s>0.0) ? 1.0 : ((s == 0.0) ? 0.0 : -1.0);

		    if ((sgn>0)||(random<0))
		    {
			t = fabs(Psi[i2]*alpha*pow((eps1+fabs(s)),spl));
		    
			if (fabs(s)>t){
			    Lind[i1][ig] = i2;
			    Lval[i1][ig] = s - sgn* t;
			    ig++;
			}
		    }
		}
		La[i1] = ig;
	    }
	    

	} else {


	    for (i1=0;i1<K;i1++)
	    {
		ig_vec[i1]=0;
	    }
    
	    for (i2 = 0; i2 < n; i2++){
		i4=0;
		for (i1=0;i1<K;i1++)
		{
		    s=0.0;
		    for (i3 = 0; i3 < K; i3++)
			s +=  sum1[i2][i3]*ichol[i3][i1];
		    
		    sgn = (s>0.0) ? 1.0 : ((s == 0.0) ? 0.0 : -1.0);

		    if ((sgn>0)||(non_negative<=0))
		    {
			t = fabs(Psi[i2]*alpha*pow((eps1+fabs(s)),spl));
			
			if (fabs(s)>t){
			    LLind[i4] = i1;
			    LLval[i4] = s - sgn* t;
			    i4++;
			}
		    }
		}
	
		if (i4>nL-1)
		{
		    for (i3 = 0; i3 < nL; i3++)
		    {
			e_ssxx_n[i3] = -1.0;
		    }
		    for (i3 = 0; i3 < i4; i3++)
		    {
			s=fabs(LLval[i3]);
			if (s>e_ssxx_n[nL-1])
			{
			    jg=nL-1;
			    while ((jg>0)&&(s>e_ssxx_n[jg-1])) {
				e_ssxx_n[jg]=e_ssxx_n[jg-1];
				jg--;
			    }
			    e_ssxx_n[jg]=s;
			}
		    }
		    s=e_ssxx_n[nL-1];
		    
		    for (i3 = 0; i3 < i4; i3++)
		    {
			if (s<=fabs(LLval[i3]))
			{
			    Lind[LLind[i3]][ig_vec[i1]] = i2;
			    Lval[LLind[i3]][ig_vec[i1]] = LLval[i3];
			    ig_vec[LLind[i3]]++;
		    
			}
		    }
		    
		} else {
		    for (i3 = 0; i3 < i4; i3++)
		    {
			Lind[LLind[i3]][ig_vec[i1]] = i2;
			Lval[LLind[i3]][ig_vec[i1]] = LLval[i3];
			ig_vec[LLind[i3]]++;
		    }
		    
		}

	
	    }
    
    
	    for (i1=0;i1<K;i1++)
	    {
		La[i1] = ig_vec[i1];
	    }


	}



	if ((lL>0)&&(lL<n)&&(cyc>bL))
	{

	  for (i2 = 0; i2 < K; i2++)
	    {
	      if (La[i2]>lL) {
		for (i1=0;i1<lL;i1++)
		  {
		    Psi[i1] = -1.0;
		  }
		for (ig=0; ig< La[i2];ig++) {
		  s=fabs(Lval[i2][ig]);
		  if (s>Psi[lL-1])
		    {
		      i3=lL-1;
		      while ((i3>0)&&(s>Psi[i3-1])) {
			Psi[i3]=Psi[i3-1];
			i3--;
		      }
		      Psi[i3]=s;
		    }
		}
		s=Psi[lL-1];
		i3=0;
		for (ig=0; ig< La[i2];ig++) {
		  if(s<fabs(Lval[i2][ig]))
		    {
		      Lval[i2][i3]= Lval[i2][ig];
		      Lind[i2][i3] = Lind[i2][ig];
		      i3++;
		    }
		}
		La[i2]=i3;
	      }
	    }
	}




  	for (i1=0;i1<n;i1++)
	    Psi[i1] = 0.0;
	
	for (i2 = 0; i2 < K; i2++) {
	    for (ig=0; ig< La[i2];ig++) {
		Psi[Lind[i2][ig]]+= Lval[i2][ig]*sum1[Lind[i2][ig]][i2];
	    }
	}


	for (i1=0;i1<n;i1++)
	{
	    
	    Psi[i1] = XX[i1] - in*Psi[i1];
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
		
		for (ig=0; ig< La[i2];ig++) {
		    s += Lval[i2][ig]*Lval[i2][ig];
		}
		
		s*=in;
		s=sqrt(s)+eps1;
		s=scale/s;
		for (ig=0; ig< La[i2];ig++) {
		    Lval[i2][ig] *=s;
		}
		s*=s;
		s = pow(s,spz);
		for (j=0;j<nn;j++)
		{
		    lapla[j][i2]*=s;
		}
		
	    }
	}		




//	if (i%20==0) {
	if (iter>1) {
	  Rprintf("Iteration: %d || Cycle: %d\r", (ite+1), i);
	  R_CheckUserInterrupt();
	} else {
	  Rprintf("Cycle: %d\r", i);
	  R_CheckUserInterrupt();
	}
//	}

    }

    if (iter>1) {
      Rprintf("Iteration: %d || Cycle: %d\n", (ite+1), (cyc-1));
    } else {
      Rprintf("Cycle: %d\n", (cyc-1));
    }


    for (i2 = 0; i2 < K; i2++) {
	LPsia[i2]=La[i2];
	for (ig=0; ig< La[i2];ig++) {
	    LPsiind[i2][ig] = Lind[i2][ig];
	    LPsival[i2][ig] = Lval[i2][ig]/Psi[Lind[i2][ig]];
	}
    }



 
    for (i2=0;i2<K;i2++)
      {
	s = 0.0;  
	for (ig=0;ig < La[i2];ig++) {
	  s += Lval[i2][ig] * LPsival[i2][ig];
	}
	LPsiL[i2][i2] = s;
      }
    
	for (i2=0;i2<K-1;i2++)
	{
	    for (i3=i2+1;i3<K;i3++)
	    {
	      s = 0.0; 
	      ig = 0; 
	      jg = 0;
	      la1= La[i3];
	      la2= LPsia[i2];
	      while ((ig < la1) && (jg < la2)) {
		
		if (Lind[i3][ig] < LPsiind[i2][jg]) {
		  ig++; 
		  while ((ig < la1) && (Lind[i3][ig] < LPsiind[i2][jg])) {
		    ig++; 
		  }
		  if (ig >= la1) break;
		}
		
		if (Lind[i3][ig] > LPsiind[i2][jg]) {
		  jg++; 
		  while ((jg < la2) && (Lind[i3][ig] > LPsiind[i2][jg])) {
		    jg++; 
		  }
		  if (jg >= la2) break;
		  if (Lind[i3][ig] == LPsiind[i2][jg]) {
		    s += Lval[i3][ig] * LPsival[i2][jg]; 
		    ig++;
		    jg++;
		  }
		} else {
		  s += Lval[i3][ig] * LPsival[i2][jg]; 
		  ig++;
		  jg++;
		}


	      } 

	      LPsiL[i2][i3] = LPsiL[i3][i2] = s;

	    }
	}
	    



   for (j=0;j<nn;j++)
    {

      la1=xa[j];

	for (i1=0;i1<K;i1++)
	{
	    for (i2=0;i2<K;i2++)
	    {
		
		if (i1==i2) {
		    tLPsiL[i1][i1] = LPsiL[i1][i1] + lapla[j][i1];
		} else {
		    tLPsiL[i1][i2] = LPsiL[i1][i2];
		}
	    }
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




	for (i3 = 0; i3 < K; i3++)
	  {
	    s=0.0;
	    ig = 0; 
	    jg = 0;
	    la2= LPsia[i3];
	    while ((ig < la1) && (jg < la2)) {
	      
	      if (xind[j][ig] < LPsiind[i3][jg]) {
		ig++;
		while ((ig < la1) && (xind[j][ig] < LPsiind[i3][jg])) {
		  ig++; 
		}
		if (ig >= la1) break;
	      }
	      if (xind[j][ig] > LPsiind[i3][jg]) {
		jg++;
		while ((jg < la2) && (xind[j][ig] > LPsiind[i3][jg])) {
		  jg++;
		}
		if (jg >= la2) break;
		if (xind[j][ig] == LPsiind[i3][jg]) {
		  s += xval[j][ig] * LPsival[i3][jg]; 
		  ig++;
		  jg++;
		} 
	      } else  {
		s += xval[j][ig] * LPsival[i3][jg]; 
		ig++;
		jg++;
	      } 
	    }
	    e_ssxx_n[i3]= s;
	  }
	

	for (i1=0;i1<K;i1++)
	{
	    
	    t=0.0;
	    for (i3 = 0; i3 < K; i3++)
	    {
		t += e_ssxx_n[i3]*ichol[i1][i3];
	    }
	    if ((t<eps1) && (non_negative>0))
	      {
		t=0.0;
	      }

	    E_SX_tmp[j][ite*K+i1] = (double) t;
	}


    }



    if (write_file>0)
    {

    if (iter == 1)
    {

	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_res_L.txt");
	pFile = fopen (sst,"w");
	fprintf(pFile,"%d\n",K); 
	fprintf(pFile,"%d\n",n); 

	for(i = 0; i < K; i ++) {

	    fprintf(pFile,"%d\n",La[i]); 
	    for(j = 0; j < La[i]; j ++)
		fprintf(pFile,"%d ",Lind[i][j]);
	    fprintf(pFile,"\n");
	    for(j = 0; j < La[i]; j ++)
		fprintf(pFile,"%.8f ",Lval[i][j]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);

	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_res_Z_full.txt");
	pFile = fopen (sst,"w");
	for(i = 0; i < K; i ++) {
	    for(j = 0; j < nn; j ++)
	      fprintf(pFile,"%.8f ",E_SX_tmp[j][i]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);

	{
	    int ind[nn];
	    double val[nn];
	    sst[0]=0;
	    strcat(sst,file_name);
	    strcat(sst,"_res_Z.txt");
	    pFile = fopen (sst,"w");
	    fprintf(pFile,"%d\n",K); 
	    fprintf(pFile,"%d\n",nn);
	    for(i = 0; i < K; i ++) {
		ig=0;
		for(j = 0; j < nn; j ++) {
		    s= E_SX_tmp[j][i];
		    if (fabs(s)>0.00001){
			ind[ig] = j;
			val[ig] = s;
			ig++;
		    }
		}
		fprintf(pFile,"%d\n",ig); 
		for(j = 0; j < ig; j ++)
		    fprintf(pFile,"%d ",ind[j]);
		fprintf(pFile,"\n");
		for(j = 0; j < ig; j ++)
		    fprintf(pFile,"%.8f ",val[j]);
		fprintf(pFile,"\n");
	    }
	    fclose (pFile);
	}

	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_res_Psi.txt");
	pFile = fopen (sst,"w");
	for(j = 0; j < nn; j ++)
	    fprintf(pFile,"%.8f ",Psi[j]);
	fprintf(pFile,"\n");
	fclose (pFile);

	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_res_lapla.txt");
	pFile = fopen (sst,"w");
	for(i = 0; i < K; i ++) {
	    for(j = 0; j < nn; j ++)
		fprintf(pFile,"%.8f ",lapla[j][i]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);
    } else {


	iterc[0]=0;
        sprintf(iterc,"%d",(ite+1)); 

	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_");
	strcat(sst,iterc);
	strcat(sst,"_res_L.txt");
	pFile = fopen (sst,"w");
	fprintf(pFile,"%d\n",K); 
	fprintf(pFile,"%d\n",n); 

	for(i = 0; i < K; i ++) {

	    fprintf(pFile,"%d\n",La[i]); 
	    for(j = 0; j < La[i]; j ++)
		fprintf(pFile,"%d ",Lind[i][j]);
	    fprintf(pFile,"\n");
	    for(j = 0; j < La[i]; j ++)
		fprintf(pFile,"%.8f ",Lval[i][j]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);


	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_");
	strcat(sst,iterc);
	strcat(sst,"_res_Z_full.txt");
	pFile = fopen (sst,"w");
	for(i = 0; i < K; i ++) {
	    for(j = 0; j < nn; j ++)
	      fprintf(pFile,"%.8f ", E_SX_tmp[j][ite*K+i]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);


	{
	    int ind[nn];
	    double val[nn];
	    sst[0]=0;
	    strcat(sst,file_name);
	    strcat(sst,"_");
	    strcat(sst,iterc);
	    strcat(sst,"_res_Z.txt");
	    pFile = fopen (sst,"w");
	    fprintf(pFile,"%d\n",K); 
	    fprintf(pFile,"%d\n",nn);
	    for(i = 0; i < K; i ++) {
	        ig=0;
		for(j = 0; j < nn; j ++) {
		    s= E_SX_tmp[j][ite*K+i];
		    if (fabs(s)>0.00001){
			ind[ig] = j;
			val[ig] = s;
			ig++;
		    }
		}
		fprintf(pFile,"%d\n",ig); 
		for(j = 0; j < ig; j ++)
		    fprintf(pFile,"%d ",ind[j]);
		fprintf(pFile,"\n");
		for(j = 0; j < ig; j ++)
		    fprintf(pFile,"%.8f ",val[j]);
		fprintf(pFile,"\n");
	    }
	    fclose (pFile);
	}

 	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_");
	strcat(sst,iterc);
	strcat(sst,"_res_Psi.txt");
	pFile = fopen (sst,"w");
	for(j = 0; j < nn; j ++)
	    fprintf(pFile,"%.8f ",Psi[j]);
	fprintf(pFile,"\n");
	fclose (pFile);

 
	sst[0]=0;
	strcat(sst,file_name);
	strcat(sst,"_");
	strcat(sst,iterc);
	strcat(sst,"_res_lapla.txt");
	pFile = fopen (sst,"w");
	for(i = 0; i < K; i ++) {
	    for(j = 0; j < nn; j ++)
		fprintf(pFile,"%.8f ",lapla[j][i]);
	    fprintf(pFile,"\n");
	}
	fclose (pFile);
 

    }

    }

 

    for (i1=0;i1<K;i1++)
    {
	for (i2 = 0; i2 < n; i2++)
	  REAL(L_n)[(long) (ite*K*n+i2 + n*i1)] = 0.0;
    }
    for (i2 = 0; i2 < K; i2++) {
	for (ig=0; ig< La[i2];ig++) {
	  REAL(L_n)[(long) (ite*K*n+Lind[i2][ig]+n*i2)]= Lval[i2][ig];
	}
    }


    for(i = 0; i < n; i++)
      REAL(Psi_n)[(long) (ite*n+i)] = (double) Psi[i];


    for(i = 0; i < nn; i++)
	for(j = 0; j < K; j++)
	  REAL(lapla_n)[(long) (ite*K*nn+i + nn*j)] = (double) lapla[i][j];


    if ((iter>1)&&(ite<(iter-1))) {

      for (i2 = 0; i2 < K; i2++) {
	for (ig=0; ig< La[i2];ig++) {
	  LPsival[0][ig] = - fabs(Lval[i2][ig]);
	}
	R_rsort(LPsival[0], La[i2]);
	LLval[i2] = -LPsival[0][nquant];
      }
      



      for (i2 = 0; i2 < n; i2++) {
	LPsiind[0][i2] = 0;
      }

      for (i2 = 0; i2 < K; i2++) {
	for (ig=0; ig< La[i2];ig++) {
	  if (fabs(Lval[i2][ig])>=LLval[i2]) {
	    LPsiind[0][Lind[i2][ig]]=1;
	  }
	}
      }

      for (i2 = 0; i2 < nn; i2++) {
	for (ig=0,jg=0; (ig+jg) < xa[i2];) {
	  if ( LPsiind[0][xind[i2][ig+jg]]==0 )
	    {
	      if (jg>0) {
		xval[i2][ig] = xval[i2][ig+jg];
		xind[i2][ig] = xind[i2][ig+jg];
	      }
	      ig++;
	    } else 
	    {
	      jg++;
	    }
	}
	xa[i2]-=jg;
      }
  
    }

    } // iter


    for(j = 0; j < nn; j++)
      for(i = 0; i < K*iter; i++) {
	REAL(E_SX_n)[(long) (i + iter*K*j)] = E_SX_tmp[j][i];
      }



    R_Free (xa );
    R_Free (xind[0]);
    R_Free (xind );
    R_Free (xval[0]);
    R_Free (xval );

    R_Free (E_SX_tmp[0]);
    R_Free (E_SX_tmp);
 
    R_Free (sum1[0]);
    R_Free (sum1 );

    R_Free (sum2[0]);
    R_Free (sum2 );

    R_Free (e_sx_n );

    R_Free (e_ssxx_n );



    R_Free (LLval );

    R_Free (LLind );

    R_Free (ig_vec );


    R_Free (XX );
    
    R_Free (LPsia );
    R_Free (LPsiind[0]);
    R_Free (LPsiind );
    R_Free (LPsival[0]);
    R_Free (LPsival );


    R_Free (ichol[0]);
    R_Free (ichol);

    R_Free (tLPsiL[0]);
    R_Free (tLPsiL );

    R_Free (LPsiL[0]);
    R_Free (LPsiL);
 
    R_Free (LPsi[0]);
    R_Free (LPsi );

    R_Free (La );
    R_Free (Lind[0]);
    R_Free (Lind );
    R_Free (Lval[0]);
    R_Free (Lval );





    R_Free (Psi );





    R_Free (lapla[0]);
    R_Free (lapla );

    PutRNGstate();

    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesRET, 0, mkChar("L"));
    SET_STRING_ELT(namesRET, 1, mkChar("E_SX_n"));
    SET_STRING_ELT(namesRET, 2, mkChar("Psi"));
    SET_STRING_ELT(namesRET, 3, mkChar("lapla"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(RET, 0, L_n);
    SET_VECTOR_ELT(RET, 1, E_SX_n);
    SET_VECTOR_ELT(RET, 2, Psi_n);
    SET_VECTOR_ELT(RET, 3, lapla_n);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(6);
    return(RET);

}


SEXP samplesPerFeature(SEXP file_nameS, SEXP samplesS,SEXP lowerBS, SEXP upperBS) {


    FILE *pFile;


    char sst[200]; 


    int hpp,ig,samp;



    int  i,j,i1,i2,n,nn;
    
    double fs;

    const char *file_name=CHAR(STRING_ELT(file_nameS,0));


    int *xa;
    int **xind;
    double **xval;

    double lowerB = (double)(REAL(lowerBS)[0]);
    double upperB = (double)(REAL(upperBS)[0]);

    int nsamp = length(samplesS);

    int *samples =  INTEGER(samplesS);
    if (samples[0]>0) {
      samp = 0;
      
    } else {
      samp=-1;
    }

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,".txt");

    pFile = fopen (sst,"r");

    if (!(pFile>0)) {
	Rprintf("File >%s< not found! Stop.\n", sst);
	return R_NilValue;
    }

    fscanf(pFile,"%d\n",&nn);  

    if (!(nn>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    fscanf(pFile,"%d\n",&n);  

    if (!(n>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    SEXP xLL,PhiA;

    PROTECT(xLL = allocVector(VECSXP, n)); 
    PROTECT(PhiA = allocVector(INTSXP, n));
    int *Phi = INTEGER(PhiA);
    int *Psi = R_Calloc(n, int); 
    double *XX = R_Calloc(n, double); 


    if (samp<0) {

      xa = (int *) R_Calloc(nn, int); 
      xind = (int **) R_Calloc(nn, int *);
      xind[0] = R_Calloc((long) nn*n, int);
      for(i=0; i < nn ; i++)
	{
	  xind[i] =  xind[0] + i*n;
	}
      xval = (double **) R_Calloc(nn, double *);
      xval[0] = R_Calloc((long) nn*n, double);
      for(i=0; i < nn; i++)
	{
	  xval[i] = xval[0] + i*n;
	}

    
      for(i = 0; i < nn; i ++)
	{
	  fscanf(pFile,"%d\n",&ig); 
	  xa[i]=ig;
	  for(j = 0; j <  ig; j ++) {
	    fscanf(pFile,"%d",&hpp);
	    xind[i][j]=hpp;
	  }
	  fscanf(pFile,"\n");
	  for(j = 0; j < ig; j ++) {
	    fscanf(pFile,"%lf",&fs);
	    xval[i][j] = fs;
	  }
	  fscanf(pFile,"\n");
	}

    } else {
      xa = (int *) R_Calloc(nsamp, int); 
      xind = (int **) R_Calloc(nsamp, int *);
      xind[0] = R_Calloc((long) nsamp*n, int);
      for(i=0; i < nsamp ; i++)
	{
	  xind[i] =  xind[0] + i*n;
	}
      xval = (double **) R_Calloc(nsamp, double *);
      xval[0] = R_Calloc((long) nsamp*n, double);
      for(i=0; i < nsamp; i++)
	{
	  xval[i] = xval[0] + i*n;
	}

      for(i = 0; i < nn; i ++)
	{
	  if ((samples[samp]-1)>i)
	    {
	      fscanf(pFile,"%d\n",&ig);
	      for(j = 0; j <  ig; j ++) {
		fscanf(pFile,"%d",&hpp);
	      }
	      fscanf(pFile,"\n");
	      for(j = 0; j < ig; j ++) {
		fscanf(pFile,"%lf",&fs);
	      }
	      fscanf(pFile,"\n");
	      
	    } else {
	      fscanf(pFile,"%d\n",&ig); 
	      xa[samp]=ig;
	      for(j = 0; j <  ig; j ++) {
	        fscanf(pFile,"%d",&hpp);
	        xind[samp][j]=hpp;
	      }
	      fscanf(pFile,"\n");
	      for(j = 0; j < ig; j ++) {
	        fscanf(pFile,"%lf",&fs);
	        xval[samp][j] = fs;
	      }
	      fscanf(pFile,"\n");
	      samp++;
	      if (samp == nsamp) break;
	    }

	}
      if (samp!=nsamp)
	{
	  Rprintf("Only %d of %d samples found! Some sample numbers are too large. Continue.\n", samp,nsamp);
	}
      nn=samp;
      Rprintf("Using %d samples!\n",samp);
    }
    fclose (pFile);



    for (i1=0;i1<n;i1++)
      {
	XX[i1] = 0.0;
	Psi[i1] = 0;
	Phi[i1] = 0;
      }	

    
   for (i2 = 0; i2 < nn; i2++) {
	for (ig=0; ig< xa[i2];ig++) {
	    XX[xind[i2][ig]] += xval[i2][ig];
	    Psi[xind[i2][ig]]++;
	}
    }



      for (i2 = 0; i2 < n; i2++) {
	if ((XX[i2]>lowerB)&&(XX[i2]<upperB)&&(Psi[i2]>0)) {
	  SET_VECTOR_ELT(xLL, i2, allocVector(INTSXP, Psi[i2])); 
	} else {
	  Psi[i2] = 0;
	  SET_VECTOR_ELT(xLL, i2, allocVector(INTSXP, 1)); 
          INTEGER(VECTOR_ELT(xLL, i2))[0] = 0; 

	}
      }



   for (i2 = 0; i2 < nn; i2++) {
     for (ig=0; ig< xa[i2];ig++) {
       if (Psi[xind[i2][ig]]>0) {
	 INTEGER(VECTOR_ELT(xLL, xind[i2][ig]))[Phi[xind[i2][ig]]] = i2+1;
	 Phi[xind[i2][ig]]++;
       }
     }
   }


    R_Free (xa );
    R_Free (xind[0]);
    R_Free (xind );
    R_Free (xval[0]);
    R_Free (xval );

    R_Free (XX);
    R_Free (Psi);


    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 2));
    SET_STRING_ELT(namesRET, 0, mkChar("sL"));
    SET_STRING_ELT(namesRET, 1, mkChar("nsL"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(RET, 0, xLL);
    SET_VECTOR_ELT(RET, 1, PhiA);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(4);
    return(RET);


 }



SEXP readSamplesSpfabic(SEXP file_nameS, SEXP samplesS,SEXP lowerBS, SEXP upperBS) {


    FILE *pFile;


    char sst[200]; 


    int hpp,ig,jg,samp;

    int  i,j,i1,i2,n,nn;
    
    double fs;





    const char *file_name=CHAR(STRING_ELT(file_nameS,0));

  
    int *xa;
    int **xind;
    double **xval;

    double lowerB = (double)(REAL(lowerBS)[0]);
    double upperB = (double)(REAL(upperBS)[0]);

    int nsamp = length(samplesS);

    int *samples =  INTEGER(samplesS);
    if (samples[0]>0) {
      samp = 0;
      
    } else {
      samp=-1;
    }

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,".txt");

    pFile = fopen (sst,"r");

    if (!(pFile>0)) {
	Rprintf("File >%s< not found! Stop.\n", sst);
	return R_NilValue;
    }

    fscanf(pFile,"%d\n",&nn);  

    if (!(nn>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    fscanf(pFile,"%d\n",&n);  

    if (!(n>0)) {
      fclose (pFile);
      Rprintf("Wrong file format (sparse file format required)! Stop.\n");
      return R_NilValue;
    }

    if (samp<0) {

      xa = (int *) R_Calloc(nn, int); 
      xind = (int **) R_Calloc(nn, int *);
      xind[0] = R_Calloc((long) nn*n, int);
      for(i=0; i < nn ; i++)
	{
	  xind[i] =  xind[0] + i*n;
	}
      xval = (double **) R_Calloc(nn, double *);
      xval[0] = R_Calloc((long) nn*n, double);
      for(i=0; i < nn; i++)
	{
	  xval[i] = xval[0] + i*n;
	}

    
      for(i = 0; i < nn; i ++)
	{
	  fscanf(pFile,"%d\n",&ig); 
	  xa[i]=ig;
	  for(j = 0; j <  ig; j ++) {
	    fscanf(pFile,"%d",&hpp);
	    xind[i][j]=hpp;
	  }
	  fscanf(pFile,"\n");
	  for(j = 0; j < ig; j ++) {
	    fscanf(pFile,"%lf",&fs);
	    xval[i][j] = fs;
	  }
	  fscanf(pFile,"\n");
	}

    } else {
      xa = (int *) R_Calloc(nsamp, int); 
      xind = (int **) R_Calloc(nsamp, int *);
      xind[0] = R_Calloc((long) nsamp*n, int);
      for(i=0; i < nsamp ; i++)
	{
	  xind[i] =  xind[0] + i*n;
	}
      xval = (double **) R_Calloc(nsamp, double *);
      xval[0] = R_Calloc((long) nsamp*n, double);
      for(i=0; i < nsamp; i++)
	{
	  xval[i] = xval[0] + i*n;
	}

      for(i = 0; i < nn; i ++)
	{
	  if ((samples[samp]-1)>i)
	    {
	      fscanf(pFile,"%d\n",&ig);
	      for(j = 0; j <  ig; j ++) {
		fscanf(pFile,"%d",&hpp);
	      }
	      fscanf(pFile,"\n");
	      for(j = 0; j < ig; j ++) {
		fscanf(pFile,"%lf",&fs);
	      }
	      fscanf(pFile,"\n");
	      
	    } else {
	      fscanf(pFile,"%d\n",&ig); 
	      xa[samp]=ig;
	      for(j = 0; j <  ig; j ++) {
	        fscanf(pFile,"%d",&hpp);
	        xind[samp][j]=hpp;
	      }
	      fscanf(pFile,"\n");
	      for(j = 0; j < ig; j ++) {
	        fscanf(pFile,"%lf",&fs);
	        xval[samp][j] = fs;
	      }
	      fscanf(pFile,"\n");
	      samp++;
	      if (samp == nsamp) break;
	    }

	}
      if (samp!=nsamp)
	{
	  Rprintf("Only %d of %d samples found! Some sample numbers are too large. Continue.\n", samp,nsamp);
	}
      nn=samp;
      Rprintf("Using %d samples!\n",samp);
    }
    fclose (pFile);

    int *Psi = R_Calloc(n, int); 
    double *XX = R_Calloc(n, double); 

    for (i1=0;i1<n;i1++)
      {
	XX[i1] = 0.0;
      }	
    
   for (i2 = 0; i2 < nn; i2++) {
	for (ig=0; ig< xa[i2];ig++) {
	    XX[xind[i2][ig]] += xval[i2][ig];
	}
    }


      for (i2 = 0; i2 < n; i2++) {
	if ((XX[i2]>lowerB)&&(XX[i2]<upperB)) {
	  Psi[i2] = 0;
	} else {
	  Psi[i2] = 1;
	}
      }


      for (i2 = 0; i2 < nn; i2++) {
	for (ig=0,jg=0; (ig+jg) < xa[i2];) {
	  if ( Psi[xind[i2][ig+jg]]==0 )
	    {
	      if (jg>0) {
		xval[i2][ig] = xval[i2][ig+jg];
		xind[i2][ig] = xind[i2][ig+jg];
	      }
	      ig++;
	    } else 
	    {
	      jg++;
	    }
	}
	xa[i2]-=jg;
      }




 
 
    SEXP X_n;
    PROTECT(X_n = allocMatrix(REALSXP, n, nn));

    for (i1=0;i1<nn;i1++)
    {
	for (i2 = 0; i2 < n; i2++)
	    REAL(X_n)[i2 + n*i1] = 0.0;
    }
    for (i2 = 0; i2 < nn; i2++) {
	for (ig=0; ig< xa[i2];ig++) {
	     REAL(X_n)[xind[i2][ig]+n*i2]= xval[i2][ig];
	}
    }

    R_Free (xa );
    R_Free (xind[0]);
    R_Free (xind );
    R_Free (xval[0]);
    R_Free (xval );

    R_Free (XX);
    R_Free (Psi);


    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 1));
    SET_STRING_ELT(namesRET, 0, mkChar("X"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(RET, 0, X_n);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(3);
    return(RET);

}




SEXP readSpfabicResult(SEXP file_nameS) {


    FILE *pFile;


    char sst[200]; 


    int  i,j,K,n,nn,ini,i1,i2,ig;
    
    double inf;



    const char *file_name=CHAR(STRING_ELT(file_nameS,0));



    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,"_res_L.txt");
    pFile = fopen (sst,"r");

     if (!(pFile>0)) {
       fclose (pFile);
       Rprintf("File >%s< not found! Stop.\n", sst);
       return R_NilValue;
    }

    fscanf(pFile,"%d\n",&K); 
    if (!(K>0)) {
       fclose (pFile);
       Rprintf("Wrong file format  >%s< (K) (sparse file format required)! Stop.\n", sst);
       return R_NilValue;
    }


    fscanf(pFile,"%d\n",&n); 
     if (!(n>0)) {
       fclose (pFile);
       Rprintf("Wrong file format  >%s< (n) (sparse file format required)! Stop.\n", sst);
       return R_NilValue;
    }
   
    int *La = R_Calloc(K, int); 
    int **Lind = R_Calloc(K, int *);
    Lind[0] = R_Calloc((long) K*n, int);
    for(i=0; i < K; i++)
    {
	Lind[i] = Lind[0] + i*n;
    }
    double **Lval = R_Calloc(K, double *);
    Lval[0] = R_Calloc((long) K*n, double);
    for(i=0; i < K; i++)
    {
	Lval[i] = Lval[0] + i*n;
    }

    for(i = 0; i < K; i ++) {
	
	fscanf(pFile,"%d\n",&ini);
	La[i]=ini;
	for(j = 0; j < La[i]; j ++)
	  {
	    fscanf(pFile,"%d ",&ini);
	    Lind[i][j] = ini;
	  }
	fscanf(pFile,"\n");
	for(j = 0; j < La[i]; j ++)
	  {
	    fscanf(pFile,"%lf ",&inf);
	    Lval[i][j] = inf;
	  }
	fscanf(pFile,"\n");
    }
    fclose (pFile);
    
    SEXP L_n;
    PROTECT(L_n = allocMatrix(REALSXP, n, K));

    for (i1=0;i1<K;i1++)
    {
	for (i2 = 0; i2 < n; i2++)
	    REAL(L_n)[i2 + n*i1] = 0.0;
    }
    for (i2 = 0; i2 < K; i2++) {
	for (ig=0; ig< La[i2];ig++) {
	     REAL(L_n)[Lind[i2][ig]+n*i2]= Lval[i2][ig];
	}
    }


    R_Free (La );
    R_Free (Lind[0]);
    R_Free (Lind );
    R_Free (Lval[0]);
    R_Free (Lval );

    

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,"_res_Z.txt");
    pFile = fopen (sst,"r");
    fscanf(pFile,"%d\n",&ini); 
    if (ini!=K) {
      fclose (pFile);
      Rprintf("K in >%s< is %d whereas K in >%s_res_L.txt< is %d! Stop.\n", sst,ini,file_name,K);
      return R_NilValue;
    }
    fscanf(pFile,"%d\n",&nn);
    if (!(nn>0)) {
      fclose (pFile);
      Rprintf("Wrong file format  >%s< (sparse file format required)! Stop.\n", sst);
      return R_NilValue;
    }
    fclose (pFile);



    SEXP E_SX_n;
    PROTECT(E_SX_n = allocMatrix(REALSXP, K, nn));


    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,"_res_Z_full.txt");
    pFile = fopen (sst,"r");
    for(i = 0; i < K; i ++) {
	for(j = 0; j < nn; j ++)
	{
	    fscanf(pFile,"%lf ",&inf);
	    REAL(E_SX_n)[i+K*j] = (double) inf;
	}
	fscanf(pFile,"\n");
    }
    fclose (pFile);
    


    SEXP Psi_n;
    PROTECT(Psi_n = allocVector(REALSXP, n));

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,"_res_Psi.txt");
    pFile = fopen (sst,"r");
    for(j = 0; j < nn; j ++)
    {
	fscanf(pFile,"%lf ",&inf);
	REAL(Psi_n)[j] = (double) inf;
    }
    fclose (pFile);
    

    SEXP lapla_n;
    PROTECT(lapla_n = allocMatrix(REALSXP, nn,K));

    sst[0]=0;
    strcat(sst,file_name);
    strcat(sst,"_res_lapla.txt");
    pFile = fopen (sst,"r");
    for(i = 0; i < K; i ++) {
	for(j = 0; j < nn; j ++)
	{
	    fscanf(pFile,"%lf ",&inf);
	    REAL(lapla_n)[i*nn + j] = (double) inf;
	}
	fscanf(pFile,"\n");
    }
    fclose (pFile);
    





    SEXP namesRET;
    PROTECT(namesRET = allocVector(STRSXP, 4));
    SET_STRING_ELT(namesRET, 0, mkChar("L"));
    SET_STRING_ELT(namesRET, 1, mkChar("E_SX_n"));
    SET_STRING_ELT(namesRET, 2, mkChar("Psi"));
    SET_STRING_ELT(namesRET, 3, mkChar("lapla"));
    
    SEXP RET;
    PROTECT(RET = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(RET, 0, L_n);
    SET_VECTOR_ELT(RET, 1, E_SX_n);
    SET_VECTOR_ELT(RET, 2, Psi_n);
    SET_VECTOR_ELT(RET, 3, lapla_n);
    setAttrib(RET, R_NamesSymbol, namesRET);
    UNPROTECT(6);
    return(RET);

}



 R_CallMethodDef callMethods[]  = {
       {"fabic", (DL_FUNC) &fabic, 16},
       {"fabics", (DL_FUNC) &fabics, 13},
       {"spfabic", (DL_FUNC) &spfabic, 25},
       {"readSamplesSpfabic", (DL_FUNC) &readSamplesSpfabic, 4},
       {"samplesPerFeature", (DL_FUNC) &samplesPerFeature, 4},
       {"readSpfabicResult", (DL_FUNC) &readSpfabicResult, 1},
       {NULL, NULL, 0}
     };



 void  R_init_myLib(DllInfo *info)
     {
	 R_registerRoutines(info, NULL, callMethods, NULL, NULL);
     }


int main() {

  return(1);
}


