#include "../../include/RationalNumberPolynomial/AlgebricFActorization.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support_Recursive-AA.h"
//#include "../../include/RationalNumberPolynomial/SMQP_Support-AA.h"

// #include <fstream>      // std::fstream, to write to file

// using namespace std;
// using namespace NTL;

AltArr_t* Factorization_AA (AltArr_t* a, AltArr_t* b, int *shift, mpq_t cont , AltArr_t**  HornerRet, AltArr_t **OriginalPassPoly,AltArr_t**contentPoly,AltArr_t**ZFac)
{
   /*  AltArr_t *minimalPoly = makePolynomial_AA(3, 2);// 3 term and 2 variable
     mpq_init(minimalPoly->elems[0].coef);
    mpq_set_si(minimalPoly->elems[0].coef, -3, 1); */

   // partialDegreeTerm_AA(minimalPoly, 1, 1); return second degree of the second term
    

   //printAA(minimalPoly);
   
 /*  const degree_t degsList[2] = {1, 1};
 
setDegrees_AA_inp(minimalPoly, 0,degsList, 2);// first term ,nvar=2
 */
    //    void Factorization_AA ( minimalPoly,  minimalPoly,  factors, exponet );*/
 
/* minimalPoly->size=1;
    const char* syms[]={"x","y"};
printPoly_AA(stderr, minimalPoly,syms,2);
 const char* syms[]={"z","a"};
fprintf(stderr,"\na=");
printPoly_AA(stderr,a,syms,2); */

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////above is the test for the poly
/////////////////////////////////////////////////////////////////////////

//RecArr_t* convertedpoly= convertToRecursiveArrayAtIdx (b, 1);



long *S;
long k=0;
S=&k;
mpq_t ret;
mpq_init(ret);
 AltArr_t* res_a = NULL;
 AltArr_t* res_r = NULL;
(*ZFac)=NULL;

/* const char* syms[]={"a"};

  fprintf(stderr,"\n\n a=");
  printPoly_AA(stderr,a,syms,1); */


/*   if(partialDegreeTerm_AA(b,b->size,0)!=0){// set Zfactor=z^e;
    AltArr_t* Zfactor;
    Zfactor=makePolynomial_AA(1,2);
    (Zfactor->size)=1;
    mpq_init(Zfactor->elems[0].coef);
    mpq_set_si(Zfactor->elems[0].coef,1,1);
    setPartialDegreeTerm_AA(Zfactor,0,0, partialDegreeTerm_AA(b,(b->size)-1,0));
    setPartialDegreeTerm_AA(Zfactor,0,1, 0);
         
    dividePolynomials_AA(b, Zfactor, &res_a,  &res_r, 2);
    freePolynomial_AA(b);
    b=deepCopyPolynomial_AA(res_a);                
    freePolynomial_AA(res_a);
    freePolynomial_AA(res_r);
    (*ZFac)=deepCopyPolynomial_AA(Zfactor);
    freePolynomial_AA(Zfactor);     
       
  }  */
    const char* tmpvars[]={"x","y"};
          integralContent_AA(b,ret);//factor a // return mainleading var
          gmp_fprintf(stderr,"\n ret before =%Qd\n",ret);
          printPoly_AA(stderr, b, tmpvars ,2);  
         AltArr_t* primitive=primitivePartAndContent_AA(b, ret);
         gmp_fprintf(stderr,"\n after ret**=%Qd\n",ret);  
                       

         AltArr_t* MainLeadingCoefPrimitive=mainLeadingCoefficient_AA(primitive);
         fprintf(stderr, "\n\n MainleadingCoef-.nvar=%d", MainLeadingCoefPrimitive->nvar);
         //////////////////////////////////////////content
           mpq_set_si(cont,1,1);
            mpq_div(cont,cont,ret);
            mpq_set(cont,ret);
            //////////////////////////////
           
          
         multiplyByRational_AA_inp(MainLeadingCoefPrimitive,cont);  // =content * contentPoly
        (*contentPoly)=deepCopyPolynomial_AA(MainLeadingCoefPrimitive);
        AltArr_t*InvLeadCoefPri= algInv(MainLeadingCoefPrimitive,a);
     
         fprintf(stderr, "\n\n Our  gmp_ret ===============================" );
            gmp_fprintf(stderr,"\n ret=%Qd",ret);  
             fprintf(stderr, "\n\n Our  primitivePoly===============================" );
                        const char* symp[]={"z","a"};            
             printPoly_AA(stderr,primitive,symp,2);
                 fprintf(stderr, "\n\n Our  contentPoly===============================" );
                  const char* symm[]={"a"};            
             printPoly_AA(stderr,MainLeadingCoefPrimitive,symm,1);
                     fprintf(stderr, "\n\n Our  contentPoly===============================" ); 
                        fprintf(stderr, "\n\n Our  contentPoly===============================" );
                        const char* syms[]={"a"};
        printPoly_AA(stderr,(*contentPoly),syms,1);
         fprintf(stderr, "\n\n Our  contentPoly===============================\n " );
        expandNumVarsLeft_AA(InvLeadCoefPri, 2);
         AltArr_t* primitivePoly= multiplyPolynomials_AA(InvLeadCoefPri,primitive,2);
         
        expandNumVarsLeft_AA(a, 2);

        dividePolynomials_AA(primitivePoly, a, &res_a,  &res_r, 2);
        (*OriginalPassPoly)=deepCopyPolynomial_AA(res_r);
           shrinkNumVarsAtIdx_AA(a,0);
         
           AltArr_t* squareFree= SquareFree(res_r, a, S,HornerRet);

           fprintf(stderr, "\n\n S=%ld", *S);
         
           gmp_fprintf(stderr,"\n\n cont=%Qd",ret);     

            *shift=*S; 
            
           gmp_fprintf(stderr,"\n\n ret=%Qd",cont);  
          //freePolynomial_AA(squareFree);
          freePolynomial_AA(res_a);
          freePolynomial_AA(res_r);
          freePolynomial_AA(primitive);
          freePolynomial_AA(MainLeadingCoefPrimitive);
          freePolynomial_AA(InvLeadCoefPri);
          freePolynomial_AA(primitivePoly);
          mpq_clear(ret);          


    return squareFree ;
}

AltArr_t * Hornershift(AltArr_t* rawPoly, AltArr_t* minimalPoly, int S)
{

  if(S==0)
  {
    return rawPoly;
  }
 
    degree_t tempDeg_levelOne=0;      
    AltArr_t* bivariatePoly;  //for levelone + bivarpoly
    AltArr_t* MinimalBivar;
    bivariatePoly=makePolynomial_AA(rawPoly->size ,2);
 
    MinimalBivar=makePolynomial_AA(minimalPoly->size ,2);

    for(int i=0;i<minimalPoly->size;i++){  // creat a BIVARIATE poly from minimalpoly
       MinimalBivar->size=(MinimalBivar->size)+1;
       mpq_init(MinimalBivar->elems[i].coef);
      //fprintf(stderr, "\n \n j===%d",i);
      mpq_set(MinimalBivar->elems[i].coef, minimalPoly->elems[i].coef);
      // degree_t a_deg=partialDegreeTerm_AA(levelOnePoly,i,0);// degree of a
       //fprintf(stderr, "\n a_deg=%d",a_deg);


      degree_t newDegs[2];
       newDegs[1] = partialDegreeTerm_AA(minimalPoly,i, 0);              
       newDegs[0] =0 ; 
       setDegrees_AA_inp(MinimalBivar, i, newDegs, 2);
               
      //setPartialDegreeTerm_AA(oneToTwoPoly,i,1, a_deg);
      // setPartialDegreeTerm_AA(oneToTwoPoly,i,0,0);
                
    }


    
    for (int i=0; i<rawPoly->size;i++){
      int j=i;
      // find the end of x^i
      AltArr_t *levelOnePoly =makePolynomial_AA(rawPoly->size, 1);  // 3 term and 2 variable for this while
        
        int countOne=0;
        while(partialDegreeTerm_AA(rawPoly,i,0)==partialDegreeTerm_AA(rawPoly,j,0)){
           /* fprintf(stderr,"\n\n i=%d, j=%d \n\n",i,j);
           fprintf(stderr, " first partialDegree_AA_rawpoly(%d)=%d, partialDegree_rowpoly(%d)=%d ",i,partialDegreeTerm_AA(rawPoly,i,0),j,partialDegreeTerm_AA(rawPoly,j,0));
           fprintf(stderr, "secondpartialDegree_AA_rawpoly(%d)=%d, partialDegree_rowpoly(%d)=%d ",i,partialDegreeTerm_AA(rawPoly,i,1),j,partialDegreeTerm_AA(rawPoly,j,1)); */

            countOne=countOne +1;
           
            tempDeg_levelOne=partialDegreeTerm_AA(rawPoly,j,1);  // degree of every same degree
        
            levelOnePoly->size=countOne;
            /* fprintf (stderr,"\n\n counter=%d",countOne);
            fprintf(stderr, "\n \n size=%d  partialdegree=%d \n",levelOnePoly->size,tempDeg_levelOne ); */
            
            mpq_init(levelOnePoly->elems[countOne-1].coef);
           // gmp_fprintf(stderr, "rawPoly->elems[%d].coef=%Qd ",j,rawPoly->elems[j].coef);
            mpq_set(levelOnePoly->elems[countOne-1].coef, rawPoly->elems[j].coef); 
            setPartialDegreeTerm_AA(levelOnePoly, countOne-1 , 0,tempDeg_levelOne );
                       
           //creating the level one poly 
        
            
             /* const char* syms1[]={"a"};
            fprintf(stderr,"\n\n onelevelpoly=");
            printPoly_AA(stderr,levelOnePoly,syms1,1);  */

          
           // tempDegMain_levelfirst=partialDegreeTerm(rawPoly,j-1,1); 
           
             j=j+1;
             if (j>=rawPoly->size){
             break;
            }


        } //end while
         
           degree_t diffdeg=0;

           if (rawPoly->size!=j){
            diffdeg=partialDegreeTerm_AA(rawPoly,i,0)-partialDegreeTerm_AA(rawPoly,j,0);
          
          }
          else
          {
            diffdeg=partialDegreeTerm_AA(rawPoly,i,0);
          }  
          
           fprintf(stderr,"diff=%ld",diffdeg);
          if(diffdeg==0){
             //AltArr_t *levelOnePoly_bivar =makePolynomial_AA(levelOnePoly->size, 1);
            AltArr_t * levelOnePoly_bivar=deepCopyPolynomial_AA(levelOnePoly); 
            expandNumVarsLeft_AA( levelOnePoly_bivar, 2);
             const char* symsbiva[]={"z","a"};
               fprintf(stderr, "\n\n levelonebivar=");
               printPoly_AA(stderr  ,levelOnePoly_bivar,symsbiva, 2);
               
               fprintf(stderr, "\n ");
            bivariatePoly=addPolynomials_AA_inp(bivariatePoly,levelOnePoly_bivar,2);
            freePolynomial_AA(levelOnePoly_bivar);
             const char* symscomb[]={"z","a"};
               fprintf(stderr, "\n\n?????????????????????????????????? finalhorner=");
               printPoly_AA(stderr  ,bivariatePoly,symscomb, 2);
               
               fprintf(stderr, "\n ");
          }
             for(int l=0;l<diffdeg;l++){
               AltArr_t* tempPoly;

                const char* syms1[]={"a"};
                fprintf(stderr,"\n\n LEVELONEPOLY=");
                printPoly_AA(stderr,levelOnePoly,syms1,1);



                   tempPoly=combPoly( levelOnePoly,  bivariatePoly , S);
                  const char* symste[]={"z","a"};
                               
                 tempPoly=sortPolynomial_AA(tempPoly);
                fprintf(stderr,"\n\n tempoly=");
                printPoly_AA(stderr,tempPoly,symste,2);

                   bivariatePoly=deepCopyPolynomial_AA(tempPoly);
                  // freePolynomial_AA(tempPoly);
                   
                 const char* syms[]={"z","a"};
                               
               
                fprintf(stderr,"\n\n BIVARPOLy=");
                printPoly_AA(stderr,bivariatePoly,syms,2);
                   // tempPoly=changePoly_AA(bivariatePoly,minimalPoly);
                       AltArr_t* res_a = NULL;
                       AltArr_t* res_r = NULL; 
               
                   dividePolynomials_AA(bivariatePoly, MinimalBivar,  &res_a,  &res_r, 2);

                   // tempPoly=deepCopyPolynomial_AA(bivariatePoly);// in other case should be commented
                      
                     bivariatePoly=deepCopyPolynomial_AA(res_r);
                    fprintf(stderr,"\n\n BIVARPOLy_rs=");
                    printPoly_AA(stderr,bivariatePoly,syms,2);
                      
                    freePolynomial_AA(tempPoly);
                    freePolynomial_AA(levelOnePoly);
                    freePolynomial_AA(res_a);
                    freePolynomial_AA(res_r);

                    levelOnePoly=NULL;
                       
              }  
             
            freePolynomial_AA(levelOnePoly);
            //freePolynomial_AA(MinimalBivar);
        i=j-1;  
    }

 
/*    const char* syms[]={"z","a"};

  fprintf(stderr,"\n\n bivariatePoly=");
  printPoly_AA(stderr,bivariatePoly,syms,2);  */
  freePolynomial_AA(MinimalBivar);
 
      return bivariatePoly;
}
////////////////////////////////////////////////  changePoly

AltArr_t * changePoly_AA(AltArr_t* rawPoly, AltArr_t* minimalPoly)
{
  /*const char* syms1[]={"z","a"};
   fprintf(stderr,"\n\nRRRRRRRRrawpoly=");
  printPoly_AA(stderr,rawPoly,syms1,2); */
    
  degree_t tempDeg_levelOne=0;
  AltArr_t* compeletPoly=makePolynomial_AA(rawPoly->size,2 );
   for (int i=0; i<rawPoly->size;i++){
      int j=i;
      // find the end of x^i
      //AltArr_t* tempLevelFirstPoly=makePolynomial_AA(1,2);
         AltArr_t *levelOnePoly =makePolynomial_AA(rawPoly->size, 1);  // 3 term and 2 variable for this while

        int countOne=0;
        while(partialDegreeTerm_AA(rawPoly,i,0)==partialDegreeTerm_AA(rawPoly,j,0)){
        /*  fprintf(stderr,"\n\n i=%d, j=%d \n",i,j);
            fprintf(stderr, " first partialDegree_AA_rawpoly(%d)=%d, partialDegree_rowpoly(%d)=%d ",i,partialDegreeTerm_AA(rawPoly,i,0),j,partialDegreeTerm_AA(rawPoly,j,0));
           fprintf(stderr, "secondpartialDegree_AA_rawpoly(%d)=%d, partialDegree_rowpoly(%d)=%d ",i,partialDegreeTerm_AA(rawPoly,i,1),j,partialDegreeTerm_AA(rawPoly,j,1)); */
 
            countOne=countOne +1;
           
            tempDeg_levelOne=partialDegreeTerm_AA(rawPoly,j,1);  // degree of every same degree
          
            levelOnePoly->size=countOne;
          /*   fprintf (stderr,"\n\n counter=%d",countOne);
            fprintf(stderr, "\n \n size=%d  partialdegree=%d \n",levelOnePoly->size,tempDeg_levelOne );
         */         
            mpq_init(levelOnePoly->elems[countOne-1].coef);
            // gmp_fprintf(stderr, "rawPoly->elems[%d].coef=%Qd ",j,rawPoly->elems[j].coef);
            mpq_set(levelOnePoly->elems[countOne-1].coef, rawPoly->elems[j].coef); 
            setPartialDegreeTerm_AA(levelOnePoly, countOne-1 , 0,  tempDeg_levelOne );
                       
           //creating the level one poly 
                
           j=j+1;
           if (j>=rawPoly->size){
             break;
            }
            

        } //end while
 /*            const char* syms1[]={"a"};
            fprintf(stderr,"\n\nLLLLL levenOnepoly=");
            printPoly_AA(stderr,levelOnePoly,syms1,1);

            const char* syms3[]={"a"};
            fprintf(stderr,"\n\nLLLLL MinimalPoly=");
            printPoly_AA(stderr,minimalPoly,syms3,1); */

          AltArr_t* res_a = NULL;
       
          AltArr_t* res_r = NULL; 
               

          dividePolynomials_AA(levelOnePoly, minimalPoly,  &res_a,  &res_r, 1);
               /* const char* syms2[]={"a"};
              fprintf(stderr,"\n\n res_r=");
              printPoly_AA(stderr,res_r,syms2,1);
              fprintf(stderr,"\n");  */
////////////////////////////////////////////////////////////////////////
              for(int k=0;k<res_r->size;k++){
                //fprintf(stderr," \n\nk=%d \n ",k);
   
              compeletPoly->size=(compeletPoly->size)+1;
             // fprintf(stderr, "\n compeletpoly->size=%d",compeletPoly->size);

             //fprintf(stderr,"\n\n k=%d , res_r->size=%d \n ", k,res_r->size);
              mpq_init(compeletPoly->elems[compeletPoly->size-1].coef);
              
                mpq_set(compeletPoly->elems[compeletPoly->size-1].coef, res_r->elems[k].coef);
             
                // fprintf(stderr, "\n\n z_deg =%d  \n",z_deg);
                // setPartialDegreeTerm_AA(compeletPoly,(compeletPoly->size)-1,0,partialDegreeTerm_AA(rawPoly,i,0)); 
                // fprintf(stderr, "\n\n i=%d, z_deg :new=%d   \n",i,partialDegreeTerm_AA(compeletPoly,(compeletPoly->size)-1,0));

                degree_t newDegs[compeletPoly->nvar];
                newDegs[0] = partialDegreeTerm_AA(rawPoly, i, 0);
                newDegs[1] = partialDegreeTerm_AA(res_r, k, 0); 
                setDegrees_AA_inp(compeletPoly, compeletPoly->size-1, newDegs, compeletPoly->nvar);
                             
                /* const char* symstot[]={"z","a"};
                fprintf(stderr,"\n compeletpoly=");

                 printPoly_AA(stderr  ,compeletPoly,symstot, 2);
                 fprintf(stderr,"\n end"); */
                               
              }
              compeletPoly=sortPolynomial_AA(compeletPoly);
                           
            freePolynomial_AA(levelOnePoly);
            freePolynomial_AA(res_a);
            freePolynomial_AA(res_r); 
        
         
        i=j-1;
    }     
    
      return  compeletPoly ;
}
/////////////////////////////// combPoly       //////////  level one poly +bivariatr*(z-a)
AltArr_t * combPoly(AltArr_t* levelOnePoly, AltArr_t* bivariatePoly, int S ){
            fprintf(stderr,"\n\n");
            const char* symi1[]={"a"};
            fprintf(stderr, "levelonepoly=");
            printPoly_AA(stderr, levelOnePoly,symi1,1);
            fprintf(stderr,"\n");
            const char* symi2[]={"z","a"};

            fprintf(stderr, "bivariatepoly=");
            printPoly_AA(stderr, bivariatePoly,symi2,2);
            fprintf(stderr,"\n");


  
            if(levelOnePoly==NULL ||(levelOnePoly->size)==0){
            levelOnePoly=makePolynomial_AA(1,2);
            levelOnePoly->size=1;
            mpq_init(levelOnePoly->elems[0].coef);
            mpq_set_d(levelOnePoly->elems[0].coef,0);
            setPartialDegreeTerm_AA(levelOnePoly,0,0, 0);

            } 
           
               AltArr_t* oneToTwoPoly=makePolynomial_AA((levelOnePoly->size)+(bivariatePoly->size),2);
               AltArr_t * combinePoly =makePolynomial_AA((levelOnePoly->size)+(bivariatePoly->size), 2);
        
         
            for(int i=0;i<levelOnePoly->size;i++){  // creat a BIVARIATE poly from levelonepoly
                oneToTwoPoly->size=(oneToTwoPoly->size)+1;
                mpq_init(oneToTwoPoly->elems[i].coef);
                //fprintf(stderr, "\n \n j===%d",i);
                mpq_set(oneToTwoPoly->elems[i].coef, levelOnePoly->elems[i].coef);
               // degree_t a_deg=partialDegreeTerm_AA(levelOnePoly,i,0);// degree of a
                //fprintf(stderr, "\n a_deg=%d",a_deg);


                degree_t newDegs[2];
                newDegs[1] = partialDegreeTerm_AA(levelOnePoly,i, 0);              
                newDegs[0] =0 ; 
                setDegrees_AA_inp(oneToTwoPoly, i, newDegs, 2);

                
                //setPartialDegreeTerm_AA(oneToTwoPoly,i,1, a_deg);
               // setPartialDegreeTerm_AA(oneToTwoPoly,i,0,0);
                
            }
            /* const char* symi1[]={"z","a"};
            fprintf(stderr, "oneToTwoPoly=");
            printPoly_AA(stderr, oneToTwoPoly,symi1,2);
            fprintf(stderr,"\n"); */
            oneToTwoPoly= addPolynomials_AA_inp(oneToTwoPoly, bivariatePoly, 2);//add bivariate with new one
            oneToTwoPoly=sortPolynomial_AA(oneToTwoPoly);

            /*  const char* symi[]={"z","a"};
            fprintf(stderr, "main_comb poly=");
            printPoly_AA(stderr, bivariatePoly,symi,2);
            fprintf(stderr,"\n");

            fprintf(stderr, "added constant_comb poly =");
            printPoly_AA(stderr, oneToTwoPoly,symi,2);
            fprintf(stderr,"\n");  */

            
            combinePoly=deepCopyPolynomial_AA(oneToTwoPoly);
             
            AltArr_t* combineShadowPoly= deepCopyPolynomial_AA(combinePoly);
           // creat -a*f=totalPoly


               
            for(int f=0;f<combinePoly->size;f++){//deg of  z+1 , a+1
              setPartialDegreeTerm_AA(combinePoly,f,0,partialDegreeTerm_AA(combinePoly, f,0)+1);
              //fprintf(stderr,"\n %d=   %d=",f, )
              setPartialDegreeTerm_AA(combineShadowPoly,f,1,partialDegreeTerm_AA(combineShadowPoly, f,1)+1);
              

            }
             /* fprintf(stderr, "first part =");
            printPoly_AA(stderr, combinePoly,symi,2);
            fprintf(stderr,"\n"); 

             fprintf(stderr, "second part =");
            printPoly_AA(stderr, combineShadowPoly,symi,2);
            fprintf(stderr,"\n");
            fprintf(stderr,"\n\n S=%d \n",S);  */


            mpq_t SS;
            mpq_init(SS);
            mpq_set_d(SS,S);
            multiplyByRational_AA_inp(combineShadowPoly,SS);

            combinePoly =subPolynomials_AA_inp(combinePoly,combineShadowPoly,2 );
            freePolynomial_AA(combineShadowPoly);
            freePolynomial_AA(oneToTwoPoly);
            mpq_clear(SS);

              fprintf(stderr, "\n\ncombinePoly=");
              const char* symscomb[]={"z","a"};
               printPoly_AA(stderr  ,combinePoly,symscomb, 2);
               fprintf(stderr, "\n "); 
               fprintf(stderr, "\n @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

              return combinePoly;
           }
////////////////////////////////////////////////////// Square Free Polynomial
AltArr_t* SquareFree(AltArr_t* Poly, AltArr_t* minimalPoly, long *S, AltArr_t**  HornerRet){
  
  AltArr_t* oneToTwoPoly=deepCopyPolynomial_AA(minimalPoly);
  expandNumVarsLeft_AA(oneToTwoPoly,2);
  AltArr_t* originalPoly=deepCopyPolynomial_AA(Poly);



            int  varMap[2]={1,0};
            int varMapsize=2;
              fprintf(stderr, "\n\n onetot=");
                const char* symsc[]={"a","z"};
                printPoly_AA(stderr  ,oneToTwoPoly,symsc, 2);
                fprintf( stderr,"\nnvar_oneto=%d\n", oneToTwoPoly->nvar);


             reorderVars_AA( oneToTwoPoly,  varMap,  varMapsize);
             reorderVars_AA( Poly,  varMap,  varMapsize);
  
               fprintf(stderr, "\n\n minimal=");
                //const char* symsc[]={"a","z"};
                printPoly_AA(stderr  ,oneToTwoPoly,symsc, 2);
                
                  fprintf(stderr, "\n\n poly=");
               
                printPoly_AA(stderr  ,Poly,symsc, 2);
                 fprintf(stderr, "\n\n ");


  const char* symscomb1[]={"a","z"};
  AltArr_t* Norm= DucosResultant(oneToTwoPoly, Poly);
   fprintf(stderr, "\n\n      $$$$$$$$$$$$$$$$$$$$$$$$$$444=");
   fprintf(stderr, "\n\n Norm=");
  printPoly_AA(stderr  ,Norm,symscomb1, 2);
   fprintf(stderr, "\n ");  
   // AltArr_t * tempNorm=makePolynomial_AA((Norm->size),1);
   AltArr_t * tempNorm=deepCopyPolynomial_AA(Norm);
         
  //shrinkNumVarsAtIdx_AA(tempNorm, 0);// reduce var to the second one
 
/*   fprintf(stderr, "\n\n tempnorm=");
     const char* symsc[]={"z"};

  printPoly_AA(stderr  ,tempNorm,symsc, 1);
   fprintf(stderr, "\n "); */
   

    AltArr_t* DerNorm=derivative_AA( tempNorm, 0, 1);
     fprintf(stderr, "\n\n DerNorm=");
    printPoly_AA(stderr  ,DerNorm,symsc, 1);
    fprintf(stderr, "\n "); 

    AltArr_t* GCD=gcd_AA (DerNorm, tempNorm);
   fprintf(stderr, "\n\n GCD=");
    printPoly_AA(stderr  ,GCD,symscomb1, 2);
    fprintf(stderr, "\n "); 
    int deg=partialDegreeTerm_AA(GCD,0, 0);

   
    freePolynomial_AA(DerNorm);
    freePolynomial_AA(GCD);
    freePolynomial_AA(tempNorm);

  AltArr_t* finalNorm;
  if(deg==0){
      
      finalNorm=deepCopyPolynomial_AA(Norm);
      *HornerRet=deepCopyPolynomial_AA(originalPoly);     
      free(originalPoly);
    }
    freePolynomial_AA(Norm);
    
 // int S=0;
        AltArr_t* tempPoly;
     // S=S+1;
/*       fprintf(stderr, "\n \n Poly=");
      printPoly_AA(stderr, Poly,symscomb,2);
      fprintf(stderr, "end");

      fprintf(stderr, "\n \n MinimalPoly=");
      printPoly_AA(stderr, minimalPoly,symscomb,2);
      fprintf(stderr, "end"); */


   while (deg!=0)
  {
    
      *S=*S+1;
      int value=*S;

      tempPoly=Hornershift(originalPoly,minimalPoly,value);
     /*fprintf(stderr, "\n \n tempPoly=");
      printPoly_AA(stderr, tempPoly,symscomb,2); */
      int varMap[2]={1,0};
      int varMapsize=2;

      reorderVars_AA( tempPoly,  varMap,  varMapsize);
       

      Norm=DucosResultant(oneToTwoPoly, tempPoly);

      tempNorm=deepCopyPolynomial_AA(Norm);         
      shrinkNumVarsAtIdx_AA(tempNorm, 0);
        
      DerNorm=derivative_AA(tempNorm,0, 1);
/*       fprintf(stderr, "\n\n DerNorm_new=");
    printPoly_AA(stderr  ,DerNorm,symsc, 1);
    fprintf(stderr, "\n ");  */

      GCD=gcd_AA (DerNorm, tempNorm);
  /*     fprintf(stderr, "\n\n GCD_new=");
    printPoly_AA(stderr  ,GCD,symscomb1, 2);
    fprintf(stderr, "\n "); */
    deg=partialDegreeTerm_AA(GCD,0, 0);
         
    freePolynomial_AA(DerNorm);
    freePolynomial_AA(GCD);
    freePolynomial_AA(tempNorm);
    
    if(deg==0){
      
      finalNorm=deepCopyPolynomial_AA(Norm);
      int varMap[2]={1,0};
      int varMapsize=2;

      reorderVars_AA( tempPoly,  varMap,  varMapsize);
         *HornerRet=tempPoly;
     
    } 
    freePolynomial_AA(Norm);
   //freePolynomial_AA(tempPoly);
   //break;
  }
  freePolynomial_AA(oneToTwoPoly);
  //freePolynomial_AA(originalPoly);
  

  return finalNorm;
}

void  Combination(AltArr_t* rawpoly , AltArr_t*** factors,long NumFac , AltArr_t** HornerRet, AltArr_t** OriginalPassPoly,int* shift){


  if (NumFac==1){
    (*factors)[0]=deepCopyPolynomial_AA((*OriginalPassPoly));
    
        fprintf(stderr, "\n\n jAVAB_numfac=1  =++++++++++++++++++++++++++>");
        const char* symsja[]={"z","a"};
        
        printPoly_AA(stderr, (*factors)[0],symsja,2);
    return;
  }
  *shift=(*shift )*( -1);
   AltArr_t*  tempFacPoly;

   AltArr_t*  tempMinimal;
   AltArr_t* res_a = NULL;  
   AltArr_t* res_r = NULL;
    AltArr_t*ans;
    AltArr_t* MainLeadingCoef;
    AltArr_t* Inv;
     tempMinimal=makePolynomial_AA(rawpoly->size,2);

     for(int k=0;k<rawpoly->size;k++){  // creat a BIVARIATE poly from minimalpoly
      
        tempMinimal->size+=1;
        mpq_init(tempMinimal->elems[k].coef);
        mpq_set(tempMinimal->elems[k].coef, rawpoly->elems[k].coef);
                      
        degree_t Degs[2];
       
        Degs[0] = 0;              
        Degs[1] =partialDegreeTerm_AA(rawpoly,k, 0) ; 
        setDegrees_AA_inp(tempMinimal, k, Degs, 2);

      }
      
    
   for(int i=0; i<NumFac;i++){
     
      tempFacPoly=makePolynomial_AA((*factors)[i]->size,2);
      for(int j=0;j<(*factors)[i]->size;j++){  // creat a BIVARIATE poly from factors
      
        tempFacPoly->size+=1;
        mpq_init(tempFacPoly->elems[j].coef);
        mpq_set(tempFacPoly->elems[j].coef, (*factors)[i]->elems[j].coef);
                      
        degree_t newDegs[2];
       
        newDegs[0] = partialDegreeTerm_AA((*factors)[i],j, 0);              
        newDegs[1] =0 ; 
        setDegrees_AA_inp(tempFacPoly, j, newDegs, 2);

      }
             /*    fprintf(stderr, "\n\n mainPoly=");
                const char* symsc[]={"z","a"};
                printPoly_AA(stderr  ,tempFacPoly,symsc, 2); */

 /* 
                 int  varMap[2]={1,0};
                 int varMapsize=2;

             
             reorderVars_AA( (*HornerRet),  varMap,  varMapsize);  */ 

             fprintf(stderr, "\n\n Hornerpart=");
            const char* symsh[]={"z","a"};
             printPoly_AA(stderr, (*HornerRet),symsh,2);
            fprintf(stderr, "\n");
    
                     
               
        ans= modLastNonZeroChain_AA (tempFacPoly, (*HornerRet),tempMinimal,0);
        
          
        
        freePolynomial_AA((*factors)[i]);
       

           
        ans= Hornershift(ans,  rawpoly,  *shift);

        MainLeadingCoef=mainLeadingCoefficient_AA(ans);
        fprintf(stderr, "\n\n MainLeadingCoef=");
        const char* syms[]={"a"};
        printPoly_AA(stderr, MainLeadingCoef,syms,1);
        fprintf(stderr, "\n");


         Inv= algInv(MainLeadingCoef,rawpoly);
        /* fprintf(stderr, "\n\n Inv=");
        const char* symsin[]={"a"};
        printPoly_AA(stderr, Inv,symsin,1);
        fprintf(stderr, "\n");  */
        expandNumVarsLeft_AA( Inv, 2);
        ans=multiplyPolynomials_AA_inp(ans,Inv,2);
   
       dividePolynomials_AA(ans,rawpoly,&res_a,&res_r,2);
        // AltArr_t* primitive=primitivePart_AA( res_r);
        (*factors)[i]=deepCopyPolynomial_AA(res_r);
        fprintf(stderr, "\n\n jAVAB=++++++++++++++++++++++++++>");
        const char* symsj[]={"z","a"};
        
        printPoly_AA(stderr, (*factors)[i],symsj,2);
        freePolynomial_AA(ans);        
        freePolynomial_AA(tempFacPoly);
        freePolynomial_AA(Inv);
        freePolynomial_AA(MainLeadingCoef);
        freePolynomial_AA(res_a);
        freePolynomial_AA(res_r);
                   
      }    
      
  freePolynomial_AA(tempMinimal);
}

AltArr_t* algInv(AltArr_t* leadingPoly,AltArr_t * minimalPoly){
  AltArr_t* r0;
  AltArr_t* r1;
  AltArr_t* r2;
  AltArr_t* t0;
  AltArr_t* t1;
  AltArr_t* t2;
  AltArr_t* s0;
  AltArr_t* s1;
  AltArr_t* s2;
  AltArr_t* q1;
  
  r0=deepCopyPolynomial_AA(minimalPoly) ;
  r1=deepCopyPolynomial_AA(leadingPoly);
  //s0=1
  s0=makePolynomial_AA((leadingPoly->size )+( minimalPoly->size),1);
  s0->size=1;
  mpq_init(s0->elems[0].coef);
  mpq_set_si(s0->elems[0].coef,1,1);
  setPartialDegreeTerm_AA(s0,0,0,0);
  

//s1=0
  //s1=deepCopyPolynomial_AA(s0);
  s1=subPolynomials_AA(s0,s0,1);

  //t0=0
   t0=deepCopyPolynomial_AA(s1);

  

//t1=1
  t1= deepCopyPolynomial_AA(s0);

fprintf(stderr, "\n\n      #####################################  leadingPloy=############################## ");
            const char* syms[]={"a"};
             printPoly_AA(stderr, leadingPoly,syms,1);
            fprintf(stderr, "\n\n");



  AltArr_t* tempq;
  AltArr_t* res_a = NULL;
  AltArr_t* res_r = NULL;  
  while(isZero_AA(r1) ==0){
   
 
    // q1=r0 quo r1           
    dividePolynomials_AA(r0, r1,  &res_a, &res_r, 1);   
    q1=deepCopyPolynomial_AA(res_a);      
    //r2=r0-q1r1 
    res_a= multiplyPolynomials_AA_inp( res_a, r1, 1);

    r2=subPolynomials_AA(r0,res_a,1);
    freePolynomial_AA(res_a);
    res_a=deepCopyPolynomial_AA(q1);

     
    //s2=s0-q1s1 
    res_a= multiplyPolynomials_AA_inp( res_a, s1, 1);
    s2=subPolynomials_AA(s0,res_a,1);
    freePolynomial_AA(res_a);
    res_a=deepCopyPolynomial_AA(q1);

     
    //t2=t0-q1t1 

    res_a= multiplyPolynomials_AA_inp( res_a, t1, 1);
    t2=subPolynomials_AA(t0,res_a,1);
    freePolynomial_AA(res_a);
    freePolynomial_AA(res_r);
    freePolynomial_AA(q1); 
   
      freePolynomial_AA(r0); 
      r0=deepCopyPolynomial_AA(r1);
     
      freePolynomial_AA(r1); 
      r1=deepCopyPolynomial_AA(r2);
     
      freePolynomial_AA(s0); 
      s0=deepCopyPolynomial_AA(s1);
      freePolynomial_AA(s1); 
      s1=deepCopyPolynomial_AA(s2);
      freePolynomial_AA(s2); 
      freePolynomial_AA(t0); 
      t0=deepCopyPolynomial_AA(t1);
      freePolynomial_AA(t1); 
      t1=deepCopyPolynomial_AA(t2); 
      freePolynomial_AA(t2); 
      freePolynomial_AA(r2); 
      if(isZero_AA(r1) !=0){
         dividePolynomials_AA(t0,r0,&res_a,&res_r,1);
         freePolynomial_AA(res_r);
        }

     
  }
  
   // dividePolynomials_AA(s0, minimalPoly,  &res_a, &res_r, 1);
     
    freePolynomial_AA(r1);
    freePolynomial_AA(s1);
    freePolynomial_AA(t1);
    freePolynomial_AA(r0);
    freePolynomial_AA(s0);
    freePolynomial_AA(t0);
  
    return res_a;

}
