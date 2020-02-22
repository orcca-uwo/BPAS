#include "RationalNumberPolynomial/urpolynomial-altarr.hpp"
/** 
 * Sparse Multivariate Rational Polynomial.
 * A multivariate polynomial with rational number coefficients using a 
 * sparse representation. 
 */
           
SparseUnivariateRationalPolynomial::SparseUnivariateRationalPolynomial() : 
    poly(NULL)
    //nvar(0)
{
	name = "x";
    // names = new Symbol[1];
    // names[0] = "1";
}
/**
 * Construct with a variable name such that f(x) = x;
 *
 * @param x: The variable name
 **/
SparseUnivariateRationalPolynomial::SparseUnivariateRationalPolynomial (const Symbol& x) :
    poly(NULL) 
  //  nvar(1)
{



    name = x;

    mpq_t coef;
    mpq_init(coef);
    mpq_set_ui(coef, 1ul, 1ul);
    poly = makeConstPolynomial_AAU(1, coef);
    poly->elems->deg = 1;
    mpq_clear(coef);

}

/**
 * Copy Constructor.
 * 
 * Does not reuse underlying memory allocated by b. 
 *
 * @param b: A sparse multivariate polynomial
 **/
SparseUnivariateRationalPolynomial::SparseUnivariateRationalPolynomial(const SparseUnivariateRationalPolynomial& b) 

{

    poly = deepCopyPolynomial_AAU(b.poly);
    name = b.name;

           
       
   //slp = std::vector<SLPRepresentation>(b.slp);
}

/**
 * Move Constructor.
 *
 * @params b: The r-value reference polynomial.
 */
SparseUnivariateRationalPolynomial::SparseUnivariateRationalPolynomial(SparseUnivariateRationalPolynomial&& b) {
    //nvar = b.nvar;//??
    poly = b.poly;
    name = b.name;


    // names = new Symbol[nvar+1];
    // std::copy(b.names, b.names+nvar+1, names);

   // slp = b.slp;

    b.poly = NULL;
  //  b.slp.clear();
}
//////////////////////////////////////////////////////////////////////////////////////

SparseUnivariateRationalPolynomial::  SparseUnivariateRationalPolynomial(AltArrU_t * aa,Symbol s){
 
    this->poly =aa;
    this-> name=s;
}

///////////////////////////////////////////////////////////////////////// ExpressionTree convertToExpressionTree()



/**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */
ExprTreeNode* SparseUnivariateRationalPolynomial::exprTreeNodeFromAAElem( AAElemU_t* n,  Symbol& sym) const{
    if (n == NULL) {
        return new ExprTreeNode(0l);
    }

    degree_t deg = n->deg;
    ExprTreeNode* t = new ExprTreeNode(mpq_class(n->coef));
    
    
      //  degree_t deg = n->deg;

        if (deg > 1) {
            //TODO pass symbol directly to expression tree
            ExprTreeNode* var = new ExprTreeNode(sym);
            ExprTreeNode* num = new ExprTreeNode(deg);
            ExprTreeNode* exp = ExprTreeNode::combineExprTreeNodes(var, num, EXPR_EXP);
            t = ExprTreeNode::combineExprTreeNodes(t, exp, EXPR_MULT);
        } else if (deg == 1) {
            ExprTreeNode* var = new ExprTreeNode(sym);
            t = ExprTreeNode::combineExprTreeNodes(t, var, EXPR_MULT);
        }
    

    return t;
}

//////////////////////////////////////////////////////////////////////////////////
ExpressionTree SparseUnivariateRationalPolynomial::convertToExpressionTree() const {
  
    if (poly==NULL||poly->size==0) { //define zero
        ExprTreeNode* r = new ExprTreeNode(0l);
        ExpressionTree t(r);
        return t; 
    }
   
    Symbol sym = this->name;
    ExprTreeNode* prev = exprTreeNodeFromAAElem(poly->elems, sym);
    

    for (int i = 1; i < poly->size; ++i) {
       
        ExprTreeNode* thisNode = exprTreeNodeFromAAElem(&(poly->elems[i]), sym);
        

        prev = ExprTreeNode::combineExprTreeNodes(prev, thisNode, EXPR_ADD);
    }

   
    return ExpressionTree(prev);
}

/**
 * Output the string representation of *this to the input ostream.
 */
void SparseUnivariateRationalPolynomial::print(std::ostream& os) const {
    os << this->polyToString();
}

////////////////////////////////////////////////// Polytostring
 std::string  SparseUnivariateRationalPolynomial::polyToString() const {
 
  if (poly == NULL) {
    return "0";
  }
 int n_size= poly->size;

  std::stringstream ss;

  bool first = true;
  bool needsMult = false;
  bool isConst = true;
  mpq_t coef;
   int counter=0;
   mpq_t i;
     mpq_init(i);
     mpq_set_d(i,1);
  while (/*node != NULL*/ counter < n_size ) {
   mpq_init(coef);


    mpq_set(coef,poly->elems[counter].coef);

    //coef = ratNum_class(n->elems[counter].coef);

    // coef = node->coef;
    isConst = true;
    if (mpq_sgn(coef)< 0) {
      mpq_neg(coef,coef);
      //coef *= -1; 
      ss << " - ";
    } else if (!first) {
      ss << " + ";
    }

     
    if (/*coef != 1*/  (mpq_cmp(coef,i))) {
      ss << coef;
      needsMult = true;
    }
    degree_t degs = poly->elems[counter].deg;
    
      if (degs == 0) {
      counter =counter +1;

        continue;
      }
      isConst = false;
      if (needsMult) {
        ss << "*";
      }
      ss << name;//"x";
      if (degs > 1) {
        ss << "^" << degs;
      }
      needsMult = true;
    

    /*node = node->next;*/ 
    counter =counter +1;
    first = false;
    needsMult = false;
  }
  
  if (isConst && !mpq_cmp(coef,i)/*coef == i*/) {
    ss << coef;
  }

  return ss.str();
}
/////////////////////////////////////////////////////////////////////////////////////////
////////// BPASRing /////////////////////////////////

bool SparseUnivariateRationalPolynomial::isZero() const {


    if (poly == NULL || poly->size == 0 || 
         ((poly->elems[0].deg)==0 && mpq_sgn(poly->elems[0].coef) == 0)) {
        return 1;
    }

    return 0;
}

void SparseUnivariateRationalPolynomial::zero() {
    freePolynomial_AAU(poly);
    poly = NULL;
    }

bool SparseUnivariateRationalPolynomial::isOne() const {
    if (poly != NULL) {
        if ((poly->elems[0].deg==0) && mpq_cmp_ui(poly->elems[0].coef, 1ul, 1ul) == 0) { 
            return 1;
        }
    }
    return 0;
}

void SparseUnivariateRationalPolynomial::one() {
    freePolynomial_AAU(poly);
    RationalNumber r(1);
    poly = makeConstPolynomial_AAU(1, r.get_mpq_t());
}
/////////////////////////////////////////////////////////////////////////////////////////  isconstant

int SparseUnivariateRationalPolynomial::isConstant() const {
    if (poly == NULL) {
        return 1;
    }
    if ((poly->elems->deg)==0) {
        if (mpq_sgn(poly->elems->coef) >= 0) {
            return 1;
        } else {
            return -1;
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
inline SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::unitCanonical(SparseUnivariateRationalPolynomial* u , SparseUnivariateRationalPolynomial* v ) const {
		RationalNumber lead = leadingCoefficient();
        	RationalNumber leadInv = lead.inverse();
        	if (u != NULL) {
        		*u = lead;
        	}
        	if (v != NULL) {
        		*v = leadInv;
        	}
            
        	return (*this * leadInv);
}

SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator= (const SparseUnivariateRationalPolynomial&b) {
if (this != &b) {
        freePolynomial_AAU(poly);
        poly = deepCopyPolynomial_AAU(b.poly);
        this->name=b.name;
        /*nvar = b.nvar;
        delete[] names; 
        names = new Symbol[nvar+1];
        std::copy(b.names, b.names+nvar+1, names);
        slp = b.slp;*/
    }


    return *this;
}

SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator+ (const SparseUnivariateRationalPolynomial&b) const {
 if (b.isZero()) {
        return *this;
    }   
    if (isZero()) {
        return b;
    }
    if (this->isConstant() != 0) {
        return (b + this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this + b.poly->elems->coef);
    }


 
        AltArrU_t* sum = addPolynomials_AAU(poly, b.poly);
        SparseUnivariateRationalPolynomial ret;
        ret.poly = sum;
        ret.name = this->name;
        return ret;
    
}

/**
 * Negate all the coefficients of *this. Note, that due to the 
 * sharing nature of underling nodes, this may alter the Nodes of
 * other SMQP.
 */
void SparseUnivariateRationalPolynomial::negate() {
    if (isZero()) {
        return;
    }
    negatePolynomial_AAU(poly);
}


  SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator+= (const SparseUnivariateRationalPolynomial&b) {
 *this = (*this + b);
    return *this;
  }  

    /**
     * Negation.
     */
      SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator- () const {
      	SparseUnivariateRationalPolynomial temp = *this;
    temp.negate();
    return temp;
      }


SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator- (const SparseUnivariateRationalPolynomial&b) const {

     if (b.isZero()) {
        return *this;
    } 
    if (isZero()) {
        return -b;
    }    
    if (this->isConstant() != 0) {
        SparseUnivariateRationalPolynomial negB = -b;
        return (negB + this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this - b.poly->elems->coef);
    }

    SparseUnivariateRationalPolynomial negB = -b;
    return (*this + negB);
}
    /**
     * Subtraction assignment.
     */
      SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator-= (const SparseUnivariateRationalPolynomial&b) {
      	  *this = (*this - b);
    return *this;
      }


    /**
     * Multiplication.
     */
      SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator* (const SparseUnivariateRationalPolynomial&b) const {
      	
   if (b.isZero() || this->isZero()) {
        return SparseUnivariateRationalPolynomial(NULL, name);
    }
    if (this->isConstant() != 0) {
        return (b * this->poly->elems->coef);
    }
    if (b.isConstant() != 0) {
        return (*this * b.poly->elems->coef);
    }

   
         AltArrU_t* Multy = multiplyPolynomials_AAU(poly, b.poly);
        SparseUnivariateRationalPolynomial ret;
        ret.poly = Multy;
        ret.name = this->name;
        return ret;
}
    /**
     * Multiplication assignment.
     */
      SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator*= (const SparseUnivariateRationalPolynomial&b) {
      	 *this = (*this * b);
         return *this;
      }
    
/**
 * Exponentiate *this by the input exponent integer.
 * Treats negative exponents as positive.
 */ 
      SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator^ (long long int e) const {
     
       if (e == 0) {
        SparseUnivariateRationalPolynomial ret(NULL, name);
        ret.one();
        return ret;
    } 

    

    if (isZero()) {
        return SparseUnivariateRationalPolynomial(NULL, name);
    }

    if (e == 1) {
        return SparseUnivariateRationalPolynomial(*this);
    }

    e = (e < 0) ? -e : e;
    AltArrU_t* retPoly = exponentiatePoly_AAU(poly, e);
    
    SparseUnivariateRationalPolynomial ret(retPoly, name);
    return ret;
      }
    
    /**
 * Exponentiate *this by the input exponent integer.
 * Treats negative exponents as positive.
 */ 
      SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator^= (long long int e) {
     *this = (*this ^ e);
      return *this;
      }

    /**
     * Equality test,
     *
     * returns true iff equal
     */
      bool SparseUnivariateRationalPolynomial::operator== (const SparseUnivariateRationalPolynomial&b) const {
    return this->isEqual(b);
      }

    /**
     * Inequality test,
     *
     * returns true iff not equal.
     */
     bool SparseUnivariateRationalPolynomial::operator!= (const SparseUnivariateRationalPolynomial&b) const {
      return !this->isEqual(b);
        }

///////////////////////////////////////////////////////////////////////////////////polynomial 48-32

   
   SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator= (const RationalNumber&r) {
           	 if (poly != NULL) {
        freePolynomial_AAU(poly);
        poly = NULL;
    }
    
    //slp.clear();
    poly = makeConstPolynomial_AAU(1, r.get_mpq_t());

    return *this;

           }

	 SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator+ (const RationalNumber&c) const {
		    SparseUnivariateRationalPolynomial ret = *this;

    if (isZero()) {
        if (ret.poly != NULL) {
            freePolynomial_AAU(ret.poly);
        }
        ret.poly = makeConstPolynomial_AAU(1,  c.get_mpq_t());
        return ret;
    }

    if ((ret.poly->elems[ret.poly->size-1].deg)==0) {
        mpq_add(ret.poly->elems[ret.poly->size-1].coef, ret.poly->elems[ret.poly->size-1].coef, c.get_mpq_t());
        return ret;
    } 

    
    if (ret.poly->size >= ret.poly->alloc) {

        resizePolynomial_AAU(ret.poly, ret.poly->alloc+10);
    }

   

    
    mpq_init(ret.poly->elems[ret.poly->size].coef);
    mpq_set(ret.poly->elems[ret.poly->size].coef, c.get_mpq_t());
   
    ret.poly->elems[ret.poly->size].deg = 0;
    ++(ret.poly->size);

    return ret;
		   };

    SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator+= (const RationalNumber&c) {
               //slp.clear();
		   // SparseUnivariateRationalPolynomial ret = *this;

    if (isZero()) {
        if (poly != NULL) {
            freePolynomial_AAU(poly);
        }
        poly = makeConstPolynomial_AAU(1, c.get_mpq_t());
        return *this;
    }

    if ((poly->elems[poly->size-1].deg)==0) {
        mpq_add(poly->elems[poly->size-1].coef, poly->elems[poly->size-1].coef, c.get_mpq_t());
        return *this;
    } 

    if (poly->size >= poly->alloc) {
        resizePolynomial_AAU(poly, poly->alloc+10);
    }


    mpq_init(poly->elems[poly->size].coef);
    mpq_set(poly->elems[poly->size].coef, c.get_mpq_t());
    poly->elems[poly->size].deg = 0;
    ++(poly->size);
    

    return *this;
	  };
		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator- (const RationalNumber&c) const {
           	 RationalNumber negC = -c;
              return *this + negC;

           }

		   

		   SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator-= (const RationalNumber&c) {
           	RationalNumber negC = -c;
              return *this += negC;

           }

		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator* (const RationalNumber&c) const {
            if (c == 0 || isZero()) {
                return SparseUnivariateRationalPolynomial(NULL,name);
            }

            SparseUnivariateRationalPolynomial ret = *this;

            if (c == 1) {
                return ret;
            }

            mpq_t mult;
            mpq_init(mult);
            mpq_set(mult, c.get_mpq_t());
            for(int i = 0; i < ret.poly->size; ++i) {
                mpq_mul(ret.poly->elems[i].coef, ret.poly->elems[i].coef, mult);
            }
            mpq_clear(mult);

            return ret;
        }


    SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator*= (const RationalNumber& c) {
           // slp.clear();
                if (isZero()) {
                    return *this;
                }
                if (c == 0) {
                    freePolynomial_AAU(poly);
                    poly = NULL;
                    return *this;
                }
                if (c == 1) {
                    return *this;
                }
               

                mpq_t mult;
                mpq_init(mult);
                mpq_set(mult, c.get_mpq_t());
                for(int i = 0; i < poly->size; ++i) {
                    mpq_mul(poly->elems[i].coef, poly->elems[i].coef, mult);
                    
                }
                mpq_clear(mult);

                return *this;
            }
            SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator/ (const RationalNumber&c) const {
                if (c == 0) {
                    std::cout << "BPAS: error, dividend is zero from SUQP." << std::endl;
                    exit(1);
                }

                if (isZero()) {
                    return SparseUnivariateRationalPolynomial(NULL,  name);
                }

                if (c == 1) {
                    return SparseUnivariateRationalPolynomial(*this);
                }
                
                ratNum_t rInv;
                mpq_init(rInv);
                mpq_inv(rInv, c.get_mpq_t());
                SparseUnivariateRationalPolynomial ret = (*this * rInv);
                mpq_clear(rInv);
                return ret;

          
        }
 SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator/= (const RationalNumber& c)  {
                if (c == 0) {
                    std::cout << "BPAS: error, dividend is zero from SUQP." << std::endl;
                    exit(1);
                }

                if (isZero()) {
                    return *this;
                }

                if (c == 1) {
                    return *this;
                }
                
                ratNum_t rInv;
                mpq_init(rInv);
                mpq_inv(rInv, c.get_mpq_t());
                *this = (*this * rInv);
                mpq_clear(rInv);
                return *this;

          
        }
           //////////////////////////////////////////////////////////////////////////////////////////////////// 
             // total degree
           Integer SparseUnivariateRationalPolynomial::degree() const {
           	 if (isConstant()) {
            return 0;
              }
         
              degree_t totalMax = 0;
              totalMax=poly->elems[0].deg;
              return totalMax;
      }  

		   RationalNumber SparseUnivariateRationalPolynomial::leadingCoefficient() const {
		   	  if (isZero()) {
               return RationalNumber(0);
             }

            return RationalNumber(poly->elems->coef);
		   } 

		   RationalNumber SparseUnivariateRationalPolynomial::trailingCoefficient() const {
		      if (isZero()) {
             return RationalNumber(0);
              }
              return RationalNumber(poly->elems[poly->size-1].coef);
		   }



		     bool SparseUnivariateRationalPolynomial::isConstantTermZero() const {
		    if (isZero()) {
                return 1;
              }
         
                 return !isZeroExponentVector(poly->elems[poly->size-1].deg);
               }
		   

		   Integer SparseUnivariateRationalPolynomial::numberOfTerms() const {
                                if (poly != NULL) {
                                return poly->size;
                            }
                            return 0;
                        }

          RationalNumber SparseUnivariateRationalPolynomial::content() const {
                if (isZero()) {
                    return RationalNumber(0);
                }
                if (isConstant()) {
                    RationalNumber ret(poly->elems->coef);
                    return ret;
                }

                mpq_t ret;
                mpq_init(ret);
                integralContent_AAU(poly, ret);
                RationalNumber rn(ret);
                mpq_clear(ret);
                return rn;
		   }

		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::primitivePart() const {
                    if (isZero()) {
                    return SparseUnivariateRationalPolynomial(NULL, name);
                }

                AltArrU_t* pp = primitivePart_AAU(this->poly);
                return SparseUnivariateRationalPolynomial(pp, name);
            }


///////////////////////////////////////////////////////////////////////////////////////////////// some mathematical functions
		   void SparseUnivariateRationalPolynomial::differentiate() {
                if (isZero() || isConstant()) {
                    this->zero();
                }

               AltArrU_t* temp = derivative_AAU(poly, 1);


                  SparseUnivariateRationalPolynomial ret;
                   ret.poly = temp;
                   ret.name = this->name;
                    *this = ret;

		   }; // p = dp/dx
		   void SparseUnivariateRationalPolynomial::differentiate(int k) {
                 
                if (k <= 0) {
                  return;
                 }

                if (isZero() || isConstant()) {
                    this->zero();
                }

               AltArrU_t* temp = derivative_AAU(poly, k);


                  SparseUnivariateRationalPolynomial ret;
                   ret.poly = temp;
                   ret.name = this->name;
                    *this = ret;
                  
		   }


		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::derivative() const {
 	         
                           if (isZero() || isConstant()) {
                    return SparseUnivariateRationalPolynomial(NULL, name);
                }

               AltArrU_t* temp = derivative_AAU(poly, 1);
            return SparseUnivariateRationalPolynomial(temp, name);
		   } // q = dp/dx



		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::derivative(int k) const {
                 if (k <= 0) {
                  return *this;
                 }

                if (isZero() || isConstant()) {
                    return SparseUnivariateRationalPolynomial(NULL, name);
                }

               AltArrU_t* temp = derivative_AAU(poly, k);
            return SparseUnivariateRationalPolynomial(temp, name);
		   }

            SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::Integral(int k) const {
                 if (k <= 0) {
                  return *this;
                 }

                if (isZero() ) {
                    return SparseUnivariateRationalPolynomial(NULL, name);

                }
           



               AltArrU_t* temp = integral_AAU(poly, k);
             
            
                                
               return SparseUnivariateRationalPolynomial(temp, name);
		   }
                   /**
		 * Convert current object to its k-th integral
		 *
		 * @param s: Symbol to differentiate with respect to
		 * @param k: Order of the derivative, k > 0
		 **/ 
    	 void SparseUnivariateRationalPolynomial::integrate( int k) {
    		*this = this->Integral( k);
    	}

       

 RationalNumber SparseUnivariateRationalPolynomial::evaluate(const RationalNumber& r) const {
                       
                    if (isZero() ) {
                        return 0;
                    }
                    if (isConstant()){
                        return (RationalNumber)poly->elems[0].coef;
                    }
                AltArrU_t *aa=this->poly;
                mpq_t val;
                mpq_init(val);

                evaluatePoly_AAU(aa, r.get_mpq_t(),val);

              RationalNumber ret(val) ;

            //RationalNumber ret(mpq_class(val)) ;
               mpq_clear(val);
                
                return ret;
		   	 	
		   }



		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::monicDivide(const SparseUnivariateRationalPolynomial&b) {
               if (poly==NULL)
               {
                 SparseUnivariateRationalPolynomial q;
                
             
               q.poly=NULL;


              q.name=b.name;
		    	return q;
               }


                if (b.poly==NULL )
                {
                    std::cout << "BPAS: error, the devisor is NULL." << std::endl;
                exit(1);
                }
                if (b.isZero()) {
                std::cout << "BPAS: error, dividend is zero from SUQP." << std::endl;
                exit(1);
            }
            if (name != b.name) {
			std::cout << "BPAS: error, trying to monic divide between Ring[" << name << "] and Ring[" << b.name << "]." << std::endl;
			exit(1);
		}


           AltArrU_t *aa=NULL,*cc=NULL;

            SparseUnivariateRationalPolynomial q;
                dividePolynomials_AAU(poly,b.poly,&aa,&cc);
             poly=cc;
             q.poly=aa;

             q.name=b.name; 
		   	return q;
		   }

		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::monicDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem) const {
		   std::cout << "*this " <<  *this << std::endl;
		   *rem = *this; std::cout<< "salam"<<std::endl;
	    	return rem->monicDivide(b);
           	
		   } 
		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::lazyPseudoDivide(const SparseUnivariateRationalPolynomial&b, RationalNumber* k, RationalNumber*mult) {
                         if (k==NULL)
                         {
                             AltArrU_t* qAA = NULL, *rAA = NULL;
                            int e = 0;
                            bool lazy=1;
                            univariatePseudoDividePolynomials_AAU(poly, b.poly, &qAA, &rAA, &e, lazy);
                        
      

                            if (mult != NULL) {
                                RationalNumber h(b.poly->elems->coef);
                            h ^= e;
                                *mult = h;                
                        }

                        poly=rAA;
                        SparseUnivariateRationalPolynomial ret;
                        ret.name=b.name;
                        ret.poly=qAA;
                        return ret;
                                 
                         }
                         
                         *k=1;
                        
                        AltArrU_t* qAA = NULL, *rAA = NULL;
                            int e = 0;
                            bool lazy=1;
                            univariatePseudoDividePolynomials_AAU(poly, b.poly, &qAA, &rAA, &e, lazy);
                        
      

                            if (mult != NULL) {
                                RationalNumber h(b.poly->elems->coef);
                            h ^= e;
                                *mult = h;                
                        }

                        poly=rAA;
                        SparseUnivariateRationalPolynomial ret;
                        ret.name=b.name;
                        ret.poly=qAA;
                        return ret;
                    
                            
                

		   }
		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::lazyPseudoDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem, RationalNumber*k, RationalNumber*mult) const {
             
             if (k=NULL){
                 AltArrU_t* qAA = NULL, *rAA = NULL;
                                            int e = 0;
                                            bool lazy=1;
                                            univariatePseudoDividePolynomials_AAU(poly, b.poly, &qAA, &rAA, &e, lazy);
                                        
                    

                                            if (mult != NULL) {
                                                RationalNumber h(b.poly->elems->coef);
                                            h ^= e;
                                                *mult = h;                
                                        }

                                        rem->poly=rAA;
                                        rem->name=b.name;
                                        SparseUnivariateRationalPolynomial ret;
                                        ret.name=b.name;
                                        ret.poly=qAA;
                                        return ret;

             }
                                  *k=1;
                                        
                                        AltArrU_t* qAA = NULL, *rAA = NULL;
                                            int e = 0;
                                            bool lazy=1;
                                            univariatePseudoDividePolynomials_AAU(poly, b.poly, &qAA, &rAA, &e, lazy);
                                        
                    

                                            if (mult != NULL) {
                                                RationalNumber h(b.poly->elems->coef);
                                            h ^= e;
                                                *mult = h;                
                                        }

                                        rem->poly=rAA;
                                        rem->name=b.name;
                                        SparseUnivariateRationalPolynomial ret;
                                        ret.name=b.name;
                                        ret.poly=qAA;
                                        return ret;
                                    
                            
		   }
		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::pseudoDivide(const SparseUnivariateRationalPolynomial&b, RationalNumber*mult) {
		                  AltArrU_t* qAA = NULL, *rAA = NULL;
                            int e = 0;
                            bool lazy=0;
                            univariatePseudoDividePolynomials_AAU(this->poly, b.poly, &qAA, &rAA, &e, lazy);
                        
                        

                            if (mult != NULL) {
                                RationalNumber h(b.poly->elems->coef);
                            h ^= e;
                                *mult = h;                
                        }
                        poly=rAA;
                        SparseUnivariateRationalPolynomial ret;
                        ret.name=b.name;
                        ret.poly=qAA;
                        return ret;
                                      
                    
		   }

		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::pseudoDivide(const SparseUnivariateRationalPolynomial&b, SparseUnivariateRationalPolynomial*rem, RationalNumber*mult) const {
		   	AltArrU_t* qAA = NULL, *rAA = NULL;
                            int e = 0;
                            bool lazy=0;
                            univariatePseudoDividePolynomials_AAU(this->poly, b.poly, &qAA, &rAA, &e, lazy);
                        
                        

                            if (mult != NULL) {
                                RationalNumber h(b.poly->elems->coef);
                            h ^= e;
                                *mult = h;                
                        }
                         rem->poly=rAA;
                         rem->name=b.name;
                                        SparseUnivariateRationalPolynomial ret;
                                        ret.name=b.name;
                                        ret.poly=qAA;
                                        return ret;
                       
		   }
           ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /**
         * Get a coefficient, given the exponent of  variable in d.
                */
		   RationalNumber SparseUnivariateRationalPolynomial::coefficient(int d) const {
            degree_t deg=poly->elems[0].deg;
                         
            if (isZero()) {
                return 0;
            }
            RationalNumber  ret;
            for (int i=0;  i<=deg;i++){
                if(poly->elems[i].deg==d){
               
                 RationalNumber ret(poly->elems[i].coef) ;
                    return ret;

                }

            }
           return ret;
           
		   }
           //////////////////////////////////////// set coeficient ////////////////////////////////////////////////////////////
		   void SparseUnivariateRationalPolynomial::setCoefficient(int d , const RationalNumber&r) {
               if (poly==NULL)
               { 
                   poly=makePolynomial_AAU(1);
                   mpq_init(poly->elems[0].coef);
                   mpq_set(poly->elems[0].coef,r.get_mpq_t());
                   poly->elems[0].deg=d;
                   poly->size=1;
                  
                     return;
               }
               mpq_t temp_coef,temp_coef2;
               degree_t temp_deg ,temp_deg2;
               mpq_init(temp_coef);
               mpq_init(temp_coef2);
               
              int j=0;

              if (poly->elems[0].deg<d){


                   mpq_set(temp_coef,poly->elems[j].coef);// saving in temp   
                        temp_deg=poly->elems[j].deg;

                        mpq_set(poly->elems[j].coef,r.get_mpq_t());
                        poly->elems[j].deg=d;
                   j=j+1;
                    for (j;j<(poly->size);j++)
                    {                        
                        mpq_set(temp_coef2,poly->elems[j].coef);// saving in temp
                        temp_deg2=poly->elems[j].deg;
                   

                         mpq_set(poly->elems[j].coef,temp_coef);
                         poly->elems[j].deg=  temp_deg;

                          mpq_set(temp_coef,temp_coef2);// saving in temp
                        temp_deg=temp_deg2;
                     
                        
                    }
                         mpq_set(poly->elems[j].coef,temp_coef);
                         poly->elems[j].deg=  temp_deg;

                         
                        poly->size=(poly->size)+1;
                        
                  return;              
              }

               for (j=0; poly->elems[j].deg>=d && j<poly->size;j++)
              {
                 // fprintf(stderr,)
                    if (d==poly->elems[j].deg){ 
                        mpq_set(poly->elems[j].coef,r.get_mpq_t());                                              
                                                              
                        return;
                            
                    }
                              
                  
              }     
                      
                  
                  if (poly->size==j){
                      mpq_init(poly->elems[j].coef);
                      mpq_set(poly->elems[j].coef,r.get_mpq_t());// saving in temp
                     poly->elems[j].deg=d;
                     poly->size=(poly->size) +1;
                     

                     return;
                                   
                  }
                  
                        mpq_set(temp_coef,poly->elems[j].coef);// saving in temp
                        temp_deg=poly->elems[j].deg;

                        mpq_set(poly->elems[j].coef,r.get_mpq_t());
                        poly->elems[j].deg=d;
                          j=j+1;

                    for (j;j<(poly->size);j++)
                    {
                        mpq_set(temp_coef2,poly->elems[j].coef);// saving in temp
                        temp_deg2=poly->elems[j].deg;
                   

                         mpq_set(poly->elems[j].coef,temp_coef);
                         poly->elems[j].deg=  temp_deg;

                          mpq_set(temp_coef,temp_coef2);// saving in temp
                        temp_deg=temp_deg2;
                        
                    }
                    mpq_init(poly->elems[j].coef);
                         mpq_set(poly->elems[j].coef,temp_coef);
                         poly->elems[j].deg=  temp_deg;
                  poly->size=(poly->size)+1;

                  return;

             
	    		}

//////////////////////////////////////////////////////////////// set variable ////////////////////////////////////////////////////////////



		   void SparseUnivariateRationalPolynomial::setVariableName (const Symbol&sym) {
               this->name=sym;
             return;
		   }

           //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		   Symbol SparseUnivariateRationalPolynomial::variable() const  {
		   	return name;
		   }
           ////////////////////////////////////////////////////////////////////////////////////////////////////////////// shifting right i times



		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator<< (int i) const {
               if (poly==NULL||isZero()==1)
               {
                   return *this;
               }
               SparseUnivariateRationalPolynomial  ret ;
               ret=* this;
               
                for (int j=0;j<ret.poly->size;j++){
                

                
                                    
                   ret.poly->elems[j].deg=(ret.poly->elems[j].deg)+i;

               // fprintf(stderr, "\n%s %d", "poly->elems[].deg=\n",poly->elems[j].deg);

               }
		   	return ret;
		   } // q = p * (x^i);


//////////////////////////////////////////////////////////////////////////// shifting i time right 

		   /////////////////////////////////////////////////////////////////////////// shifting i times right 

		   SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator<<= (int i) {
                 if (poly==NULL)
               {
                   return *this;
               }
                          for (int j=0;j<this->poly->size;j++){
                         this->poly->elems[j].deg+=i;
                          
                    }
               
		   	return *this;



		   } // p = p *(x^i)




           
		   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator>> (int i) const {
		   	if (poly==NULL)
               {
                   return *this;
               }
               SparseUnivariateRationalPolynomial  ret ;
               ret=* this;
              int j= (ret.poly->size)-1;
                                                                  

                for (j; j>=0; j--){
               
               
                          if   (ret.poly->elems[j].deg -i<0){
                                                     mpq_clear(ret.poly->elems[j].coef);
                                                     ret.poly->elems[j].deg=0;
                                              ret.poly->size -=1;
                                               
                          }
                          else{
                            ret.poly->elems[j].deg=(ret.poly->elems[j].deg)-i;
                          }
                   

                
               }
		   	return ret;
		   } // q = p / (x^i);
		   SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator>>= (int i) {
		   	  	if (poly==NULL)
               {
                   return *this;
               }
               
              int j= (poly->size)-1;

                for (j; j>=0; j--){
                             
                          if   (poly->elems[j].deg -i<0){
                                                     mpq_clear(poly->elems[j].coef);
                                                     poly->elems[j].deg=0;
                                                   poly->size -=1;
                          }
                          else{
                            poly->elems[j].deg=(poly->elems[j].deg)-i;
                          }
                                   
               }
		   	return *this;
		   }
//////////////////////////////////////////////////////////////////////////////////////////////// operator /                //////////////////////////////

 		SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::operator/ (const SparseUnivariateRationalPolynomial& b) const {
        SparseUnivariateRationalPolynomial q,r;
                if (this->isZero()) {
                    SparseUnivariateRationalPolynomial ret;
                    ret.zero();
                    return ret;
                }

            if (b.isZero()) {
                std::cout << "BPAS: error, dividend is zero from SUQP." << std::endl;
                exit(1);
            }

            if (b.isConstant()) {
                q = *this / b.poly->elems->coef;
                r = SparseUnivariateRationalPolynomial(NULL, name);
                return q;
            }
            if (isConstant()) {
                q = SparseUnivariateRationalPolynomial(NULL,  name);
                r = *this;
                
                mpq_t a;
                 mpq_init(a);
               mpq_set_ui(a,0,1);
                mpq_set(r.poly->elems[0].coef,a);
                return r;
            }

            dividePolynomials_AAU(poly, b.poly, &q.poly, &r.poly);

                  if (r.poly != NULL && r.poly->size != 0){
            
                   //std::cerr << "q: ==" << q<< std::endl;
                   std::cout << "BPAS: error, there is a reminder." << std::endl;
                    exit(1);         
                }    
             //  std::cerr << "r: ==" << r<< std::endl;
                SparseUnivariateRationalPolynomial ret;
                ret.poly = q.poly;
                ret.name = this->name;
                return ret;
		} 

        
	   SparseUnivariateRationalPolynomial& SparseUnivariateRationalPolynomial::operator/= (const SparseUnivariateRationalPolynomial& b) {
            *this = (*this / b);
            return *this;
	   } 
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	    Factors<SparseUnivariateRationalPolynomial> SparseUnivariateRationalPolynomial::squareFree() const {
	    	Factors<SparseUnivariateRationalPolynomial> f;

            if (isZero())
            {
                   std::cout << "BPAS: error, the input is Null or Zero" << std::endl;

                exit(1) ;
            }
            else{
                int nFactors; 
                AltArrU_t**  Fact=squareFree_AAU(poly, &nFactors);
                 
                int arraysize=nFactors;
              

                SparseUnivariateRationalPolynomial ret;
                ret.poly = Fact[nFactors];
                  ret.name = this->name;
                f.setRingElement(ret); 

                for(int i=0;i<arraysize;i++){
                   SparseUnivariateRationalPolynomial ret;
                   
                     ret.poly = Fact[i];
                      ret.name = this->name;
                

                    f.addFactor(ret,i+1);
                                     
                }
                         
            }
              
         return f;

	    }
   
   SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::gcd(const SparseUnivariateRationalPolynomial&b) const {
    if (this->isZero() || b.isZero())
   {
         std::cout << "BPAS: error, the input is Null or Zero" << std::endl;
        exit(1) ;
    }
    if(this ->name !=b.name)
    {
                 std::cout << "BPAS:  names of two polynomials are different " << std::endl;
        exit(1) ;
    }
  SparseUnivariateRationalPolynomial ret;

    ret.poly= univariateGCD_AAU(this->poly,b.poly ); 
    ret.name=this->name;
       
   	return ret;
   }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Modular GCD


SparseUnivariateRationalPolynomial SparseUnivariateRationalPolynomial::ModularGCD_AAU(const SparseUnivariateRationalPolynomial&b) const {
    if (this->isZero() || b.isZero())
   {
         std::cout << "BPAS: error, the input is Null or Zero" << std::endl;
        exit(1) ;
    }
    if(this ->name !=b.name)
    {
                 std::cout << "BPAS:  names of two polynomials are different " << std::endl;
        exit(1) ;
    }
  SparseUnivariateRationalPolynomial ret;

    ret.poly= modularGCD(this->poly,b.poly ); 
    ret.name=this->name;
       
   	return ret;
   }
////////////////////////////////////////////////////////////////// SUQP-Specific ////////////////////////////////////

bool SparseUnivariateRationalPolynomial::isEqual(const SparseUnivariateRationalPolynomial &b) const {
    if (this->isZero()) {
        if (b.isZero()) {
            return 1;
        } else {
            return 0;
        }
    }
    if (b.isZero()) {
        return 0;
    }
    if (isConstant()) {
        if (b.isConstant()) {
            mpq_class tempa(poly->elems->coef);
            mpq_class tempb(b.poly->elems->coef);
            return tempa == tempb;
        } else {
            return 0;
        }
    }
    if (poly->size != b.poly->size) {
        return 0;
    }
    //   order comparison
    bool isOrdered;
    if  (b.name==this ->name){
          isOrdered =  1;   //comparing symbols
    }
    else{
        isOrdered=0;
    }

    //comparing to polyomial with degrees and queficients 
        int ret = 1;
    degree_t adeg, bdeg;
    for (int i = 0; i < poly->size; ++ i){
        if (mpq_cmp(poly->elems[i].coef, b.poly->elems[i].coef) != 0) {
            ret = 0;
            break;
        }
        
         if (poly->elems[i].deg!=b.poly->elems[i].deg) {
            ret = 0;
            break;
        }
               
    }
    
    return ret;
}
/////////////////////////////////////////////////////////////////////// Random Poly
void SparseUnivariateRationalPolynomial::randomPolynomial( time_t seed,int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
 
    *this = SparseUnivariateRationalPolynomial();
   
        AltArrU_t* ret = buildRandomPolyU( seed,nterms, coefBound, sparsity, includeNeg);
 
        this->poly=ret;
        
    
} 
 void SparseUnivariateRationalPolynomial::randomPolynomialU( int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
 time_t seed=time(NULL);

  
     this->   randomPolynomial( seed, nterms,  coefBound, sparsity,  includeNeg);
    
} 
/////////////////////////////////////////// random max poly
void SparseUnivariateRationalPolynomial::buildRandomPolyFromMax_seeded(time_t seed,  const int maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {
 
 
    *this = SparseUnivariateRationalPolynomial();
   
        AltArrU_t* ret =  buildRandomPolyFromMax_seededgenU(seed,  maxDegs, coefBound,  sparsity,  includeNeg) ;
       this->poly=ret;
        
} 
void SparseUnivariateRationalPolynomial::buildRandomPolyFromMax_seededU( const int maxDegs, unsigned long int coefBound, float sparsity, int includeNeg) {
 time_t seed=time(NULL);
 fprintf(stderr, "\n\n seed=%ld",seed);
     this-> buildRandomPolyFromMax_seeded( seed, maxDegs, coefBound, sparsity, includeNeg);
    
} 
