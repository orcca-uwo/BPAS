#include <bpas.h>

sfixn * oldMont_dft_OPT2_AS_GENE ( sfixn n, sfixn power, sfixn * rootsPtr, sfixn * tmpVecPtr,  sfixn degA, sfixn * APtr, MONTP_OPT2_AS_GENE * pPtr )

{
 int32 i;
 register sfixn p=pPtr->P;
 sfixn t1,t2;
 register sfixn tmpw;
 sfixn m, halfn=n>>1, s, t, Tno,  start=0, tmp1, tmp2;
 sfixn  * lPtr1, * rPtr1, *endPtr;

  if(power > (pPtr->Npow)){
    //printf("The given Fourier prime number %ld can't hanlde this FFT multiplication in terms of the FFT-size. Please choose another Fourier prime with bigger FFT-size.\n", pPtr->P);
    // fflush(stdout);
   return tmpVecPtr;
 }

  if(n<2){
    //printf("n<2 won't work for FFT\n");
    //fflush(stdout);
   return tmpVecPtr;
 }

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
