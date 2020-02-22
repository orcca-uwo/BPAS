#include "../../include/RationalNumberPolynomial/mrpolynomial.h"
#include "../../include/RationalNumberPolynomial/SMQP_Support-AA.h"
#include "../../include/DataStructures/Factors.hpp"
#include "bpas.h"
#include <fstream>  
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;
/*
The C++ Function for Alhgebric Factorization ove extention  field 
*/
Factors<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::AlgebricFactorization 
(SparseMultivariateRationalPolynomial& minimal,SparseMultivariateRationalPolynomial* polyTotalcontent){
  mpq_t firstContent;
  mpq_init(firstContent);
  int sameObject;
 integralContent_AA(this->poly,firstContent);
  Factors<SparseMultivariateRationalPolynomial> Mainfacts=this->squareFree(); 


  for (size_t p = 0; p < Mainfacts.size(); p++)
  {
   // fprintf(stderr,"\nMainfacts[p].mainvarible=%s\n",Mainfacts[p].first.leadingVariable());

    std::cerr<<"\n main nvar print="<< Mainfacts[p].first.poly->nvar<<"   \n Mainfacts-befor="<<Mainfacts[p]<<"\n";
  }
  

  Factors<SparseMultivariateRationalPolynomial> Factors, subfact; 
  std::vector<Symbol> vec={'a'};
  AltArr_t* AlgebricContent;
  AltArr_t* totalAlContent;

/////////////////////////////////////////// z=z?
  AltArr_t* onez=makePolynomial_AA(1,1);
  onez->size=1;
  mpq_init(onez->elems[0].coef);
  mpq_set_si(onez->elems[0].coef,1,1);
  setPartialDegreeTerm_AA(onez,0,0,1);
  expandNumVars_AA(onez , 2);
//////////////////////////////////////
  totalAlContent=makePolynomial_AA(this->poly->alloc,1); //totalAlContent=1;
  totalAlContent->size=1;
  mpq_init(totalAlContent->elems[0].coef);
  mpq_set_si(totalAlContent->elems[0].coef,1,1);
  setPartialDegreeTerm_AA(totalAlContent,0,0,0);
  const char* syms[]={"a"};
  fprintf(stderr, "\n totalAlContent=");
  printPoly_AA(stderr,totalAlContent,syms,1);
////////////////////////////////////////////////////////////////


  for(int i=0;i<Mainfacts.size();i++){
    
  //  std::cerr << "\nlv[0] := " << Mainfacts[0].first<< std::endl;
   // std::cerr << "\nlv[1] := " << Mainfacts[1].first<< std::endl;
    std::cerr << "Mainfacts.size() = " << Mainfacts.size() << std::endl;

    std::cerr << "\t\n\n\n\n\n\n iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiIIIIIIIIIIIIIIIIIII = " << i << "\n\n\n\n\n\n\n";


    if(Mainfacts[i].first.poly->nvar==1){
     
      fprintf(stderr,"\n    vorooooooooooood a_:   \n ");
       std::cerr << "Mainfacts[i].first.poly = " << Mainfacts[i].first << std::endl;
      AlgebricContent= exponentiatePoly_AA( Mainfacts[i].first.poly,Mainfacts[i].second,1);
           

      totalAlContent=multiplyPolynomials_AA_inp(totalAlContent,AlgebricContent,1);
          freePolynomial_AA(AlgebricContent);
      fprintf(stderr, "\n totalAlContent=");
      printPoly_AA(stderr,totalAlContent,syms,1);
      //Factors.addFactor(Mainfacts[i].first ,Mainfacts[i].second);



    }else if(isExactlyEqual_AA(Mainfacts[i].first.poly, onez))// compare z with z
    {  fprintf(stderr,"\n    vorooooooooooood  z:  \n ");
      Factors.addFactor(Mainfacts[i].first ,Mainfacts[i].second);
    }else
    {  fprintf(stderr,"\n    vorooooooooooood  els:  \n ");
      
      subfact=Mainfacts[i].first.Factoring (minimal);
      fprintf(stderr,"\nsubfact.size=%d\n",subfact.size());
      for(int j=0;j<subfact.size();j++){
        sameObject=0; // it controls the same factors

       /*  std::vector<Symbol> rvars = subfact[j].first.ringVariables();
        for (auto rsym : rvars) {
          std::cerr << "rvar: " << rsym << std::endl;
        } */

             const char* syms[]={"z","a"};
             fprintf(stderr,"\n subfact[%d].first",j);
             printPoly_AA(stderr, subfact[j].first.poly,syms,2);

               std::cerr<<"\n ---------subfact.first--------->="<<"j="<< j<<"   "<<subfact[j].first<<"\n";
               
          if(subfact[j].first.poly->nvar==1){
              fprintf(stderr,"\n j=%d   vorooooooooooooooooood for content \n\n\n\n\n nahae   \n ",j);
            //std::cerr<<"\n ------------------>="<<"j="<< j<<"i="<<i<<"sub="<<subfact[j].first<<"mainfac="<<Mainfacts[i].first<<"\n";
            
            fprintf(stderr, "\n Algebriccont=");
          // printPoly_AA(stderr,AlgebricContent,syms,1);
            fprintf(stderr, "\n");  
            totalAlContent=multiplyPolynomials_AA_inp(totalAlContent,subfact[j].first.poly,1);
          }else
          {
            for(int k=0;k<Factors.size();k++){  // Is factors equal? yes?: so change the exponent
              if(subfact[j].first==Factors[k].first && sameObject==0){
                sameObject=1;
                Factors[k].second=((subfact[j].second)*(Mainfacts[i].second))+(Factors[k].second);
              }
                          
            }
            if(sameObject==0)
            {
              Factors.addFactor((subfact[j].first),(subfact[j].second)*(Mainfacts[i].second));            
            }


          }
        
        
      }
         

    }
           
          
         //    printPoly_AA(stderr, polyTotalcontent.poly,symp,1);
    
         
       
    }
/* SparseMultivariateRationalPolynomial PP=Mainfacts[0].first;
std::cerr<<"\n PPPPPPPPPPPPPPPPPP="<<PP<<"\n";
std::cerr<<"\n QQQQQQQQQQQQQQQQQQ="<<Mainfacts[1].first<<"\n";

std::cerr<<"\n OOOOOOOOOOOOOOOOOO="<<Mainfacts[2].first<<"\n";
std::cerr<< "\n nvar="<<Mainfacts[2].first.nvar<<"\n";
int d=Mainfacts[0].second;
fprintf(stderr,"dddd=%d\n ",d);
int l=Mainfacts.size();
std::cerr<<"\n l="<<l<<"\n";
 */const char* symp[]={"a"};
   fprintf(stderr,"\n totalAlContent=");
             printPoly_AA(stderr, totalAlContent,symp,1);
    Symbol symone[2] = {'1','a'};
      //  expandNumVarsLeft_AA(contentPoly,2);
        multiplyByRational_AA_inp(totalAlContent,firstContent);
        mpq_clear(firstContent);
       //*polyTotalcontent =SparseMultivariateRationalPolynomial((totalAlContent), 1,symone);

        AltArr_t* res_a = NULL;
        AltArr_t* res_r = NULL; 
               
        dividePolynomials_AA(totalAlContent, minimal.poly,  &res_a,  &res_r, 1);

       Factors.setRingElement(SparseMultivariateRationalPolynomial((res_r), 1,symone));

       
     std::cerr<<"\n polyTotalcontent-befor="<<Factors.ringElement()<<"\n";
   
     freePolynomial_AA(onez);

     freePolynomial_AA(res_a);
     free(totalAlContent);
  
    for(int h=0;h<Factors.size();h++){
      std::cerr<<"\n TotalFacts = "<<Factors[h].first<<"\n"<< "exponent="<<Factors[h].second<<"\n";
    }
   
   return Factors;
}
Factors<SparseMultivariateRationalPolynomial> SparseMultivariateRationalPolynomial::Factoring (SparseMultivariateRationalPolynomial& a)  {
/* AltArr_t* firstpoly = this->poly;
AltArr_t* secondpoly=a.poly; */
///////////////////////////////////////////////////////////////////////////
/* AltArr_t *minimalPoly = makePolynomial_AA(3, 2);// 3 term and 2 variable
    mpq_init(minimalPoly->elems[0].coef);
    mpq_set_si(minimalPoly->elems[0].coef, -3, 1);
const degree_t degsList[2] = {1, 1};
setDegrees_AA_inp(minimalPoly, 0,degsList, 2);// first term ,nvar=2 

 
minimalPoly->size=1;
const char* syms[]={"x","y"};
printPoly_AA(stderr, minimalPoly,syms,2);

SparseMultivariateRationalPolynomial polyfirst;
polyfirst=SparseMultivariateRationalPolynomial();*/

///////////////////////////////////////////////////////////////////////////

//SMQP name=SMQP("2*x^2");
//cout << "name="<<name;

 int *shift;
 int k=0;
 shift=&k;
 mpq_t cont;
 mpq_init (cont);

  mpq_t cont_NTL;
  mpq_init (cont_NTL);
 
   AltArr_t** factors;
   AltArr_t* contentPoly;
   AltArr_t* ZFac;
 long *exponet;
 AltArr_t* HornerRet = NULL;
 AltArr_t* OriginalPassPoly=NULL;
 AltArr_t* minimalPoly;
 Factors<SparseMultivariateRationalPolynomial> facts;
 
 // takes poly and then 
 AltArr_t* Norm= Factorization_AA (a.poly, this->poly, shift,  cont , &HornerRet,&OriginalPassPoly, &contentPoly,&ZFac);

 fprintf(stderr, "\n\n Shift=%d", *shift);
  const char* symcontent[]={"z","a"};
  fprintf(stderr, "\n contentPoly===========================>>>>%d",contentPoly->nvar);
             //printPoly_AA(stderr, (contentPoly),symcontent,2);
             fprintf(stderr,"\n");

           fprintf(stderr, "\n Zfac===========================>>>>");
             printPoly_AA(stderr, (ZFac),symcontent,2);
             fprintf(stderr,"\n");
  

/* const char* syms1[]={"z","a"};
fprintf(stderr,"\na=");

printPoly_AA(stderr,this->poly,syms1,2); */

//Symbol
//polyfirst =SparseMultivariateRationalPolynomial(minimalPoly, 2, syms) ;
 long NumFac = NtlFactor ( Norm, shift, cont_NTL ,  &factors, &exponet );
  mpq_mul(cont,cont_NTL,cont);
  mpq_clear(cont_NTL);
    
 const char* symi2[]={"z","a"};
 
           
             fprintf(stderr, "\n HornRet=");
             printPoly_AA(stderr, (HornerRet),symi2,2);
             fprintf(stderr,"\n");


   Combination( a.poly , &factors, NumFac , & HornerRet,&OriginalPassPoly,shift);
   freePolynomial_AA(HornerRet);
   // we pass the factors to SMQP class
   Symbol syms[3] = {'1','z','a'};
   for(int i=0;i<NumFac;i++){
     
      facts.addFactor(SparseMultivariateRationalPolynomial((factors[i]), 2,syms),exponet[i]);
       
    }   
    gmp_fprintf(stderr,"\n\n   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ cont=%Qd\n",cont);
    ///////////////////////////////////////////////////////////////////////////////////// cont add to factlist
 // multiplyByRational_AA_inp(contentPoly,cont);
  shrinkNumVarsAtIdx_AA(contentPoly,0);
  //!isConstant_AA ( contentPoly)
       Symbol symone[2] = {'1','a'};
      //  expandNumVarsLeft_AA(contentPoly,2);
        fprintf(stderr, "\n contentpoly-new=");
             fprintf(stderr, "\n");
             printPoly_AA(stderr, (contentPoly),symcontent,1);
             fprintf(stderr,"\n");      
        facts.addFactor(SparseMultivariateRationalPolynomial((contentPoly) , 1,symone),1);
      std::cerr<<"\n \n \n\nfacts[facts.size()].first=******----------------------------------->"
      <<facts[facts.size()-1].first<<"facts.size()="<<facts.size()<<"\n";
    
  
   ///////////////////////////////////////////////////////////////////////////////////////////
   /* if (!isConstant_AA ( ZFac)){

        expandNumVars_AA(ZFac,2);
        degree_t deg=partialDegreeTerm_AA(ZFac,1,0);        
        setPartialDegreeTerm_AA(ZFac,0,0, 1);

         fprintf(stderr, "\n Zfac_new=");
             fprintf(stderr, "\n");
             printPoly_AA(stderr, (ZFac),symcontent,2);
             fprintf(stderr,"\n");
        facts.addFactor(SparseMultivariateRationalPolynomial((ZFac), 2,syms),deg);
            
   } */  freePolynomial_AA(ZFac);

  mpq_clear(cont);
  free(exponet);
  freePolynomial_AA(OriginalPassPoly);
 return facts;
}

/*
The C Function for Alhgebric Factorization ove extention  field 
*/
long NtlFactor (AltArr_t* Norm,int *shift, mpq_t cont , AltArr_t*** factors, long **exponet ){

  ZZX f;

  NTL::Vec< NTL::Pair< NTL::ZZX, long >> facts;
  ZZ c;
  //shrinkNumVarsAtIdx_AA(Norm, 0);

	// std::string sc;
  NTL::ZZ nc;
  char* cc = NULL; 
  long count=-1;
   mpq_t zero;
  mpq_init(zero);
  // mpq_t mpc;
  for (int i = 0;i< Norm->size; i++) {// parsing to intl
  
  count+=1;
		cc = mpq_get_str(cc, 10, Norm->elems[(Norm->size)-1-i].coef); // base = 10
	  // const char* cc = sc.c_str();
		
		conv (nc, cc);
		SetCoeff (f, count, nc);	
    free(cc);
    cc = NULL;
   
    if(((Norm->size-1)-i-1)>=0){
      long def= partialDegreeTerm_AA(Norm,(Norm->size-1)-i-1,1)
        -partialDegreeTerm_AA(Norm,(Norm->size-1)-i,1)-1;
   
      
      if(partialDegreeTerm_AA(Norm,(Norm->size-1)-i-1,1)
        -partialDegreeTerm_AA(Norm,(Norm->size-1)-i,1)>1 ){
      
        
          
        for(long l=1;l<=def;l++){
          

         
          cc = mpq_get_str(cc, 10, zero); 
          conv (nc, cc);
          count+=1;
          SetCoeff (f, count, nc);
        }     
      }
    }
  }
   mpq_clear(zero);
  factor(c, facts, f);
  std ::cout<<"\n f="<<f<<"\n";
  std::cout << "\n c="<<c << "\n";
  std::cout << "\n facts="<<facts << "\n";
  std::cout << "\n lengh="<<facts.length() << "\n";


  ///////////////////////////////////////// convert content to mpq
  std::stringstream ssa;
  ssa << c;
  mpq_set_str( cont, ssa.str().c_str(),10);
////////////////////////////////////////  example
// ZZ ncco;
     /*  GetCoeff (ncco, facts[i].a, j);
      std::stringstream ssc;
      ssc << ncco;
	    mpq_t mpc;
      mpq_init(*factors[i]->elems[j].coef);
      mpq_set_str(*factors[i]->elems[j].coef, ssc.str().c_str(), 10);
      gmp_fprintf(stderr,"\n\n mpc=%Qd",mpc); */
////////////////////////////////////////////////

   *factors = (AltArr_t **) malloc(sizeof(AltArr_t*)*facts.length());
   *exponet=(long *)calloc(facts.length(), sizeof(long)); 
  // mpq_t no;
  // mpq_init(no);
  for (long i = 0; i < facts.length(); i++) {
  

     long degree=deg(facts[i].a);
    (*factors)[i] = makePolynomial_AA (degree+2,1);


  ////////////////////// UNCOMMENT ////////////////////
    std::stringstream ssc;
    ssc << facts[i].b;

    long k=(long)atoi( ssc.str().c_str());
  //  fprintf(stderr, "\n\n *********************************** k=%ld ,i=%d", k,i);
    (*exponet)[i]=k;
   fprintf(stderr,"\n\n ##################exponent[%ld]=%ld",i,k/* (*exponet)[i] */);

   int count =0;
   
   long lTemp(0);
   for(long j=0;j<=degree;j++){
       //pointer[i] === *(pointer + i)
       //**pointer[i] === **(*(pointer + i))
    conv(lTemp,facts[i].a[j]);
    
    if (lTemp == 0) {
      continue;
    }

    
    (*factors)[i]->size += 1;
    mpq_init ((*factors)[i]->elems[count].coef);
    mpq_set_si ((*factors)[i]->elems[count].coef,lTemp, 1ul);

    /* gmp_fprintf(stderr, "\ncount(Ind) =%d nvar=%d size=%d alloc=%d deg=%d lTemp=%ld coef=%Qd\n", 
                count, (*factors)[i]->nvar, (*factors)[i]->size, (*factors)[i]->alloc ,j ,
                lTemp, (*factors)[i]->elems[count].coef); */
  
    setPartialDegreeTerm_AA ((*factors)[i],count,0,j);
    count += 1; 
   // fprintf(stderr, "for i = %ld and j=%ld ... PASS\n",i,j);

    /*GetCoeff (ncco, facts[i].a, j);

    std::stringstream ssc;
    ssc << ncco;

    mpq_init(mpc);
    mpq_set_str(mpc, ssc.str().c_str(), 10); */


   /*  gmp_fprintf(stderr,"\n \n%Qd=", mpc);
    fprintf(stderr,"\n \n i=%d j=%d",i,  j); */
    //  gmp_fprintf(stderr,"\n\n cont=%Qd",no);


  //   if(mpq_cmp(mpc,no)!=0){
    
  //     (*factors)[i]->size=((*factors)[i]->size)+1;
  //    mpq_init((*factors)[i]->elems[count].coef);
  //   // int p=atoi( ssc.str().c_str());
  //    //fprintf(stderr,"count=%d",count);
  //    mpq_set((*factors)[i]->elems[count].coef,mpc);
  //    //gmp_fprintf(stderr, "coef[%d]=%Qd",count,(*factors)[i]->elems[count].coef );
  //    gmp_fprintf(stderr, "\ncount(Ind) =%d nvar=%d size=%d alloc=%d deg=%d coef=%Qd\n",count, (*factors)[i]->nvar, (*factors)[i]->size, (*factors)[i]->alloc ,degree-j ,(*factors)[i]->elems[count].coef);
  //     setPartialDegreeTerm_AA((*factors)[i],count,0,degree-j);
  //   //88 sortPolynomial_AA((*factors)[i]); 
  //    count =count +1;
  //    gmp_fprintf(stderr, "   mpc=%Qd\n\n",mpc);
  //  } 
     //mpq_clear(mpc); 
     
 }  
  // TODO: size could be set here with count!

  mergeSortPolynomial_AA((*factors)[i]);
  
  const char* syms1[]={"z"};
 fprintf(stderr,"\nfactors[%ld]=",i);
 printPoly_AA(stderr,(*factors)[i],syms1,1);
   
}
 freePolynomial_AA(Norm);
/* fprintf(stderr,"\n\n final  exponent[0]=%ld",(*exponet)[0]);
fprintf(stderr,"\n exponent[1]=%ld\n",(*exponet)[1]);
gmp_fprintf(stderr, "\n\n cont=%Qd",cont ); */
// mpq_clear(no);
return facts.length();
}

