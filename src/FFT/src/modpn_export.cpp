// -*- C++ -*- 
// modpn_export.cilk

//@author Yuzhen Xie

#include "../../../include/FFT/src/modpn.h"
#include "../../../include/FFT/src/modpn_export.h"

namespace MODPN {

  /**
   * Mont_dft_OPT2_AS_GENE:
   * @n: the FFT size
   * @power: log2 of n
   * @rootsPtr: powers of the n-th primitive root
   * @tmpVecPtr: (output) work space (int array of size n)
   * @degA: degree of input polynomial
   * @APtr: coefficient vector of input polynomial
   * @pPtr: prime number structure
   * 
   * 
   * 
   * 
   * Return value: DFT of the input polynomial
   **/
  sfixn * 
  Mont_dft_OPT2_AS_GENE ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr,  sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr )
	
  {
	signed int i;
	register sfixn p=pPtr->P;
	sfixn t1,t2;
	register sfixn tmpw;
	sfixn m, halfn=n>>1, s, t, Tno,  start=0, tmp1, tmp2;
	sfixn  * lPtr1, * rPtr1, *endPtr;
	

	if(power<0) return tmpVecPtr;
	
	
	m=halfn;
	
	cleanVec(n-1, tmpVecPtr);
	// DFTs begin=========================================>
	
	// compute the top level. (level-1 to level-2).
	if((degA+1)>halfn){
	  tmp1=degA-halfn+1;
	  for (i=0; i<tmp1; i++) {
		tmpVecPtr[i]=AddMod(APtr[i],APtr[i+halfn], p);
		tmpVecPtr[i+halfn]=SubMod(APtr[i],APtr[i+halfn], p);}
	  
	  for (i=tmp1;i<halfn;i++) tmpVecPtr[i]=tmpVecPtr[i+halfn]=APtr[i];
	}
	else
	  {for(i=0; i<=degA;i++) tmpVecPtr[i]=tmpVecPtr[i+halfn]=APtr[i];}
	
	
	// compute the rest levels.
	m=halfn>>1;
	for (s=2; s<power; s++){
	  // # of teams with "full butterflies";
	  Tno=n/(m<<1);
	  // # of "half bufflies" in the last team;
	  start=0;
	  
	  // middle-loop begin
	  for (t=0; t<Tno; t++){
		
		i=t<<1;
		i=partialBitRev(i,s);
		i*=m;
		tmpw=rootsPtr[i];
		//tmpw<<=pPtr->Base_Rpow;
		lPtr1=tmpVecPtr+start; rPtr1=lPtr1+m;
		//t1=MontMulMod_OPT2_AS_GENE(rPtr1[0],tmpw, pPtr);
		//t2=MontMulMod_OPT2_AS_GENE(rPtr1[1],tmpw, pPtr);
		
		t1=MontMulMod_OPT2_AS_Double_GENE(&t2, rPtr1[0], tmpw, rPtr1[1],tmpw, pPtr);
		// inner-loop begin
		endPtr=lPtr1+m-2;
		
		for(; lPtr1<endPtr; lPtr1+=2, rPtr1+=2){
		  rPtr1[0]=SubMod(lPtr1[0],t1,p);
		  lPtr1[0]=AddMod(lPtr1[0],t1,p);
		  
		  rPtr1[1]=SubMod(lPtr1[1],t2,p);
		  lPtr1[1]=AddMod(lPtr1[1],t2,p);
		  
		  t1=MontMulMod_OPT2_AS_Double_GENE(&t2, rPtr1[2],tmpw, rPtr1[3], tmpw, pPtr);
		}
		// inner-loop end.
		
		tmp1=lPtr1[0];
		tmp2=lPtr1[1];
		rPtr1[0]=SubMod(tmp1,t1,p);
		rPtr1[1]=SubMod(tmp2,t2,p);
		lPtr1[0]=AddMod(tmp1,t1,p);
		lPtr1[1]=AddMod(tmp2,t2,p);
		start+=(m<<1); }
	  
	  // middle-loop end.
	  
	  m>>=1;
	}
	
	if(power>1){
	  for (i=0; i<n; i+=2) {
		tmpw=rootsPtr[partialBitRev(i,power)];
		//tmpw<<=pPtr->Base_Rpow;
		t1=MontMulMod_OPT2_AS_GENE(tmpVecPtr[i+1],tmpw, pPtr);     
		tmp1=tmpVecPtr[i];
		tmpVecPtr[i]=AddMod(tmp1,t1, p);
		tmpVecPtr[i+1]=SubMod(tmp1,t1,p);
	  }
	}
	
	return tmpVecPtr; 
	
  }

  /**
   * Mont_dft_OPT2_AS_GENE_SPE:
   * @n: the FFT size
   * @power: log2 of n
   * @rootsPtr: powers of the n-th primitive root
   * @tmpVecPtr: (output) work space (int array of size n)
   * @degA: degree of input polynomial
   * @APtr: coefficient vector of input polynomial
   * @pPtr: prime number structure
   * 
   * Assume p-1 is a power of 2
   * 
   * 
   * Return value: DFT of the input polynomial
   **/
  sfixn * 
  Mont_dft_OPT2_AS_GENE_SPE ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr,  sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr )
  {
	
	signed int i;
	register sfixn  p=pPtr->P;
	sfixn t1,t2;
		
	sfixn m, halfn=n>>1, s, t, Tno,  start=0, tmp1, tmp2, tmpw;
	sfixn  * lPtr1, * rPtr1, *endPtr;
	
	if(power<0) return tmpVecPtr;
	
	
	m=halfn;
	cleanVec(n-1, tmpVecPtr);
	
	// DFTs begin=========================================>
	
	// compute the top level. (level-1 to level-2).
	if((degA+1)>halfn){
	  tmp1=degA-halfn+1;
	  for (i=0; i<tmp1; i++) {
		tmpVecPtr[i]=AddMod(APtr[i],APtr[i+halfn], p);
		tmpVecPtr[i+halfn]=SubMod(APtr[i],APtr[i+halfn], p);}
	  
	  for (i=tmp1;i<halfn;i++) tmpVecPtr[i]=tmpVecPtr[i+halfn]=APtr[i];
	}
	else
	  {for(i=0; i<=degA;i++) tmpVecPtr[i]=tmpVecPtr[i+halfn]=APtr[i];}
	
	// compute the rest levels.
	m=halfn>>1;
	for (s=2; s<power; s++){
	  // # of teams with "full butterflies";
	  Tno=n/(m<<1);
	  // # of "half bufflies" in the last team;
	  start=0;
	  
	  // middle-loop begin
	  for (t=0; t<Tno; t++){
		
		tmpw=rootsPtr[partialBitRev(t<<1,s)*m];
		//tmpw<<=pPtr->Base_Rpow;
		lPtr1=tmpVecPtr+start; rPtr1=lPtr1+m;
		//t1=MontMulMod_OPT2_AS_GENE(rPtr1[0],tmpw, pPtr);
		//t2=MontMulMod_OPT2_AS_GENE(rPtr1[1],tmpw, pPtr);
		
		t1=MontMulMod_OPT2_AS_Double_GENE_SPE(&t2, rPtr1[0], tmpw, rPtr1[1],tmpw, pPtr);
		// inner-loop begin
		endPtr=lPtr1+m-2;
		
		for(; lPtr1<endPtr; lPtr1+=2, rPtr1+=2){
		  rPtr1[0]=SubMod(lPtr1[0],t1,p);
		  lPtr1[0]=AddMod(lPtr1[0],t1,p);
		  
		  rPtr1[1]=SubMod(lPtr1[1],t2,p);
		  lPtr1[1]=AddMod(lPtr1[1],t2,p);
		  
		  t1=MontMulMod_OPT2_AS_Double_GENE_SPE(&t2, rPtr1[2],tmpw, rPtr1[3], tmpw, pPtr);
		}
		// inner-loop end.
		
		tmp1=lPtr1[0];
		tmp2=lPtr1[1];
		rPtr1[0]=SubMod(tmp1,t1,p);
		rPtr1[1]=SubMod(tmp2,t2,p);
		lPtr1[0]=AddMod(tmp1,t1,p);
		lPtr1[1]=AddMod(tmp2,t2,p);
		start+=(m<<1); }
	  
	  // middle-loop end.
	  
	  m>>=1;
	}
	
	if(power>1){
	  for (i=0; i<n; i+=2) {
		tmpw=rootsPtr[partialBitRev(i,power)];
		//tmpw<<=pPtr->Base_Rpow;
		t1=MontMulMod_OPT2_AS_GENE_SPE(tmpVecPtr[i+1],tmpw, pPtr);     
		tmp1=tmpVecPtr[i];
		tmpVecPtr[i]=AddMod(tmp1,t1, p);
		tmpVecPtr[i+1]=SubMod(tmp1,t1,p);
	  }
	  
	}
	return tmpVecPtr; 
  }
  
  /**
   * Mont_invdft_OPT2_AS_SPE_R:
   * @n: the FFT size
   * @power: log2 of n
   * @rootsPtr: powers of the n-th primitive root
   * @tmpVecPtr:  work space (int array of size n)
   * @degRes: degree of result
   * @ResPtr: (output) coefficient vector of result
   * @pPtr: prime number structure
   * 
   * Inverse DFT.  Case where p-1 is a power of 2.
   * 
   * 
   * Return value: Inverse DFT.
   **/
  sfixn * 
  Mont_invdft_OPT2_AS_GENE_SPE_R ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr)
	
  {
	register sfixn i, p=pPtr->P;
	// l is the size of result poly. 
	sfixn m, halfn=n>>1, s, u, t, invn, start=0, tmp1,tmp2, tmp, BRsft=pPtr->Base_Rpow;
	sfixn * wPtr, * lPtr, * rPtr;
	sfixn R=1L<<(pPtr->Rpow); R=R%pPtr->P;
	
	if (power==0) return tmpVecPtr;
	wPtr=(sfixn *)my_calloc(halfn, sizeof(sfixn));
	// lsPtr[0] not in use.
	
	// lsPtr[i] keeps the number of points needed to compute at level i for poly1/2 DFT tree.
	// the bottom level has one extra to compute, but leave it, won't hurt.
	m=halfn;
	
	for(i=0; i<n; i+=2){
      tmp=tmpVecPtr[i];
      tmpVecPtr[i]=AddMod(tmpVecPtr[i], tmpVecPtr[i+1], p);
      tmpVecPtr[i+1]=SubMod(tmp, tmpVecPtr[i+1], p);}
	m=2;
	wPtr[0]=R<<pPtr->Base_Rpow;
	
	
	for(s=power-1; s>1; s--){
	  for(i=1; i<m; i++) wPtr[i]=(rootsPtr[n-(i<<(s-1))]); //<<pPtr->Base_Rpow;
	  start=0;   
	  for(t=0; t<(n/(m<<1)); t++){
		
		lPtr=tmpVecPtr+start;
		rPtr=lPtr+m;
		for(u=0; u<m; u+=2){
		  tmp1=MontMulMod_OPT2_AS_Double_GENE_SPE(&tmp2,rPtr[u],wPtr[u],rPtr[u+1],wPtr[u+1], pPtr);
		  rPtr[u]=SubMod(lPtr[u],tmp1,p);
		  lPtr[u]=AddMod(lPtr[u],tmp1,p);
		  rPtr[u+1]=SubMod(lPtr[u+1],tmp2,p);
		  lPtr[u+1]=AddMod(lPtr[u+1],tmp2,p);}
		start+=m<<1;}
	  m<<=1;}
	

	
	if(power>1){
	  
	  lPtr=tmpVecPtr;
	  rPtr=lPtr+halfn;
	  
	  tmp1=MontMulMod_OPT2_AS_Double_GENE_SPE(&tmp2,rPtr[0],rootsPtr[0]/*<<pPtr->Base_Rpow*/,rPtr[1],rootsPtr[n-1]/*<<pPtr->Base_Rpow*/, pPtr);
	  rPtr[0]=SubMod(lPtr[0],tmp1,p);
	  lPtr[0]=AddMod(lPtr[0],tmp1,p);
	  rPtr[1]=SubMod(lPtr[1],tmp2,p);
	  lPtr[1]=AddMod(lPtr[1],tmp2,p);
	  
	  
	  for(u=2; u<halfn; u+=2){
		tmp1=MontMulMod_OPT2_AS_Double_GENE_SPE(&tmp2,rPtr[u],rootsPtr[n-u]/*<<pPtr->Base_Rpow*/,rPtr[u+1],rootsPtr[n-u-1]/*<<pPtr->Base_Rpow*/, pPtr);
		//tmp=MontMulMod_OPT2_AS_GENE(rPtr[u],wPtr[u], pPtr);
		rPtr[u]=SubMod(lPtr[u],tmp1,p);
		lPtr[u]=AddMod(lPtr[u],tmp1,p);
        //tmp=MontMulMod_OPT2_AS_GENE(rPtr[u+1],wPtr[u+1], pPtr);
		rPtr[u+1]=SubMod(lPtr[u+1],tmp2,p);
		lPtr[u+1]=AddMod(lPtr[u+1],tmp2,p);
      }
	  

	}
	
	invn=inverseMod(n,p);
	invn=MulMod(R,invn,p);
	invn=MulMod(R,invn,p);
	invn<<=BRsft;
	for(i=0; i<=degRes; i++) ResPtr[i]=MontMulMod_OPT2_AS_GENE_SPE(tmpVecPtr[i], invn, pPtr);
	my_free(wPtr);
	return(tmpVecPtr);
	
  }

  /**
   * Mont_invdft_OPT2_AS_GENE_R:
   * @n: the FFT size
   * @power: log2 of n
   * @rootsPtr: powers of the n-th primitive root
   * @tmpVecPtr:  work space (int array of size n)
   * @degRes: degree of result
   * @ResPtr: (output) coefficient vector of result
   * @pPtr: prime number structure
   * 
   * Inverse DFT. 
   * 
   * 
   * Return value: Inverse DFT  of the input polynomial.
   **/
  sfixn * 
  Mont_invdft_OPT2_AS_GENE_R ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr, sfixn degRes, sfixn * ResPtr, MONTP_OPT2_AS_GENE * pPtr)
	
  {	
	register sfixn i, p=pPtr->P;
	// l is the size of result poly. 
	sfixn m, halfn=n>>1, s, u, t, invn, start=0, tmp1,tmp2, tmp, BRsft=pPtr->Base_Rpow;
	sfixn * wPtr, * lPtr, * rPtr;
	sfixn R=1L<<(pPtr->Rpow); R=R%pPtr->P;
	
	if (power==0) return tmpVecPtr;
	wPtr=(sfixn *)my_calloc(halfn, sizeof(sfixn));
	// lsPtr[0] not in use.
	
	// lsPtr[i] keeps the number of points needed to compute at level i for poly1/2 DFT tree.
	// the bottom level has one extra to compute, but leave it, won't hurt.
	m=halfn;
		
	for(i=0; i<n; i+=2){
      tmp=tmpVecPtr[i];
      tmpVecPtr[i]=AddMod(tmpVecPtr[i], tmpVecPtr[i+1], p);
      tmpVecPtr[i+1]=SubMod(tmp, tmpVecPtr[i+1], p);}
	m=2;
	wPtr[0]=R<<pPtr->Base_Rpow;
	
	for(s=power-1; s>1; s--){
	  for(i=1; i<m; i++) wPtr[i]=(rootsPtr[n-(i<<(s-1))]); //<<pPtr->Base_Rpow;
	  start=0;   
	  for(t=0; t<(n/(m<<1)); t++){
		
		lPtr=tmpVecPtr+start;
		rPtr=lPtr+m;
		for(u=0; u<m; u+=2){
		  tmp1=MontMulMod_OPT2_AS_Double_GENE(&tmp2,rPtr[u],wPtr[u],rPtr[u+1],wPtr[u+1], pPtr);
		  rPtr[u]=SubMod(lPtr[u],tmp1,p);
		  lPtr[u]=AddMod(lPtr[u],tmp1,p);
		  rPtr[u+1]=SubMod(lPtr[u+1],tmp2,p);
		  lPtr[u+1]=AddMod(lPtr[u+1],tmp2,p);
		}
		start+=m<<1;}
	  m<<=1;}
	
	
	if(power>1){
	  lPtr=tmpVecPtr;
	  rPtr=lPtr+halfn;
	  tmp1=MontMulMod_OPT2_AS_Double_GENE(&tmp2,rPtr[0],rootsPtr[0]/*<<pPtr->Base_Rpow*/,rPtr[1],rootsPtr[n-1]/*<<pPtr->Base_Rpow*/, pPtr);
	  rPtr[0]=SubMod(lPtr[0],tmp1,p);
	  lPtr[0]=AddMod(lPtr[0],tmp1,p);
	  rPtr[1]=SubMod(lPtr[1],tmp2,p);
	  lPtr[1]=AddMod(lPtr[1],tmp2,p);
	  for(u=2; u<halfn; u+=2){
		tmp1=MontMulMod_OPT2_AS_Double_GENE(&tmp2,rPtr[u],rootsPtr[n-u]/*<<pPtr->Base_Rpow*/,rPtr[u+1],rootsPtr[n-u-1]/*<<pPtr->Base_Rpow*/, pPtr);
		//tmp=MontMulMod_OPT2_AS_GENE(rPtr[u],wPtr[u], pPtr);
		rPtr[u]=SubMod(lPtr[u],tmp1,p);
		lPtr[u]=AddMod(lPtr[u],tmp1,p);
        //tmp=MontMulMod_OPT2_AS_GENE(rPtr[u+1],wPtr[u+1], pPtr);
		rPtr[u+1]=SubMod(lPtr[u+1],tmp2,p);
		lPtr[u+1]=AddMod(lPtr[u+1],tmp2,p);
      }
	  
	}
		
	invn=inverseMod(n,p);
	invn=MulMod(R,invn,p);
	invn=MulMod(R,invn,p);
	invn<<=BRsft;
	for(i=0; i<=degRes; i++) ResPtr[i]=MontMulMod_OPT2_AS_GENE(tmpVecPtr[i], invn, pPtr);
	my_free(wPtr);
	return(tmpVecPtr);
	
  }

  //===================================================
  // To check the degrees in a given triangular set.
  // -1 means bounds of TS is NOT ok.
  //  0 means ok.
  //===================================================
  int checkDgsOfST(sfixn N, TriSet * ts){
	register int i;
	for(i=1; i<=N; i++) if(BDSI(ts, i)<0) return -1;
	return 0;
  }
  
}//end of MODPN
