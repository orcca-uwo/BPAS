#include "FiniteFields/SmallPrimeField_Support.h"


//helper function
//compute R^2 mod prime
//R = 2^64
long long int getRsquare(long long int prime){
  
  long long int res;
   __asm__(
      "movq %%rax,%%rsi\n\t"
     // "push %%rax\n\t"
      //"push %%rdx\n\t"
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "divq %1\n\t"
      "movq %%rdx,%%rax\n\t"
      "mulq %%rdx\n\t"
      "divq %1\n\t"
      //"movq %%rdx,%0\n\t"
      //"pop %%rdx\n\t"
      //"pop %%rax\n\t"
      "movq %%rsi,%%rax\n\t"
      : "=&d" (res)
      : "rm"(prime)
      :"rsi","rax");
  return res;
}

//RR' - pp' = 1
//compute p' using egcd
unsigned long long int getPp (long long int prime){
  //printf("getPp\n");
  long long int t, A, B, C, D, u, v,q;
  __asm__(                          //compute the first step of egcd where R = q*prime + v
     // "movq %%rax,%%rsi\n\t"
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "divq %2\n\t"
     // "movq %%rsi,%%rax\n\t"
     // "movq %%rdx,%0\n\t"
     // "movq %%rax,%1\n\t"
      : "=&d" (v),"=&a" (q) 
      : "b"(prime)
      :"rsi","rdi");
  A = 1;
  B = 0;
  C = 0;
  D = 1;
 // q = R/prime;
  u = prime;
  //v = R - R mod prime;
  t = A;
  A = B;
  B = t - q * B;
  t = C;
  C = D;
  D = t - q * D;

  while (v != 0){
    q = u / v;
    t = u;
    u = v;
    v = t - q * v;
    t = A;
    A = B;
    B = t - q * B;
    t = C;
    C = D;
    D = t - q * D;
 //   printf("%llu,%llu\n",A,C);
  }
   unsigned long long int res;
  if(C < 0){
    C = 0-C;
    return (unsigned long long int)C;
  }
  else{
    __asm__(                          //compute -C mod R;
     // "movq %%rax,%%rsi\n\t"
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "sub %1,%%rax\n\t"
      "sbb $0,%%rdx\n\t"
     // "movq %%rsi,%%rax\n\t"
     // "movq %%rdx,%0\n\t"
      "movq %%rax,%%rdx\n\t"
      : "=&d" (res)
      : "b"((unsigned long long int)C)
      :"rsi","rdi");
    return res;

  }
}

//return a pointer that contains the constants needed for the prime field
Prime_ptr* smallprimefield_get_prime_constants(long long int prime){
  Prime_ptr* Pptr = (Prime_ptr*)malloc(sizeof(Prime_ptr));
  Pptr->prime = prime;
  Pptr->rsquare = getRsquare(prime);
  Pptr->prime_inv = getPp(prime);
  return Pptr;
}


long long int smallprimefield_add(long long int a, long long int b, const Prime_ptr* Pptr){
  a = a + b;
  a -= Pptr->prime;
  a += (a >> 63) & Pptr->prime;
  return a;
}

long long int smallprimefield_sub(long long int a, long long int b, const Prime_ptr* Pptr){
   a = a - b;
   a += (a >> 63) & Pptr->prime;
   return a;
 }

//Montgomery multiplication based on Syvt's optimization
 //compute a* b / R mod prime
long long int smallprimefield_mul(long long int a, long long int b, const Prime_ptr* Pptr){
  __asm__(
      "mulq %2\n\t"
      "movq %%rax,%%rsi\n\t"
      "movq %%rdx,%%rdi\n\t"
      "imulq %3,%%rax\n\t"
      "mulq %4\n\t"
      "add %%rsi,%%rax\n\t"
      "adc %%rdi,%%rdx\n\t"
      "subq %4,%%rdx\n\t"
      "mov %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %4,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (a)
      : "a"(a),"rm"(b),"b"((unsigned long long int)Pptr->prime_inv),"c"(Pptr->prime)
      :"rsi","rdi");
  return a;
}

long long int smallprimefield_exp(long long int a, long long int e, const Prime_ptr* Pptr){
    long long int one = smallprimefield_convert_in(1LL,Pptr);
   // printf("one:%lld\n",one);
    if(a == 0 || a == one || e == 1)
      return a;
    if(e == 0)
      return one;
    int neg_flag = 0;
    if(e < 0)
      neg_flag = 1;
    while (e != 0) {
        if (e % 2)
          one = smallprimefield_mul(one,a,Pptr);
        a = smallprimefield_mul(a,a,Pptr);
     //   printf("e: %lld  a: %lld   one:%lld\n",e,a,one);
        e >>= 1;
    }
    if(neg_flag ==1)
      one = smallprimefield_inv(one,Pptr);
    return one;
}



//first step of Montgomery inversion
long long int* pinverse(long long int a, const Prime_ptr* Pptr){
    long long int u = a;
    long long int v = Pptr->prime;
    long long int x2 = 0;
    long long int* r = (long long int *)malloc(2*sizeof(long long int));
    r[0] = 1;
    r[1] = 0;
    while(v > 0){
	if (v%2 == 0){
	    v = v >> 1;
	    r[0] = 2 * r[0];
	}
	else if (u%2 ==0){
	    u = u >> 1;
	    x2 = 2 * x2;
	}
	else if(v>=u){
	    v = (v - u) >> 1;
	    x2 = x2 + r[0];
	    r[0] = 2 * r[0];
	}
	else{
	    u = (u - v) >> 1;
	    r[0] = x2 + r[0];
	    x2 = 2 * x2;
	}
	r[1]++;
    }
    if(u != 1){
		// printf("Not invertible!\n");
		return NULL;
    }
    if(r[0] > Pptr->prime){
		r[0] -= Pptr->prime;    
    }
    return r;
}


//Montgomery inversion 
long long int smallprimefield_inv(long long int a, const Prime_ptr* Pptr){
  if (a == 0)
  {
    // printf("Not invertible!\n");
    return 0;
  }
   // long long int rsquare =  getRsquare(prime);
    long long int* result = pinverse(a, Pptr);
    if(result[1] < 64){
		  result[0] = smallprimefield_mul(result[0], Pptr->rsquare,Pptr);
		  result[1] += 64;
    }
    result[0] = smallprimefield_mul(result[0], Pptr->rsquare,Pptr);
    long long int tmp;
    if(result[1] != 64){
	    tmp = 1L << (128 - result[1]);
	    result[0] = smallprimefield_mul(result[0], tmp,Pptr);
    }
    tmp = result[0];
    free(result);
    return tmp;
}

//simple version of division using multiplication and inversion
//TODO:optimize division in Montgomery representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
long long int smallprimefield_div(long long int a, long long int b, const Prime_ptr* Pptr){
    return smallprimefield_mul(a, smallprimefield_inv(b, Pptr),Pptr);
}

//compute aR mod p  R = 2^64
//a can be both positive and negative
long long int smallprimefield_convert_in(long long int a, const Prime_ptr* Pptr){
       __asm__(
      "xor %%rax,%%rax\n\t"
      "movq $1,%%rdx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "imulq %%rbx\n\t"
      "idivq %2\n\t"
      "movq %%rdx,%%rax\n\t"
      "sar $63,%%rax\n\t"
      "andq %2,%%rax\n\t"
      "addq %%rax,%%rdx\n\t"
      : "=&d" (a)
      : "b"(a),"c"(Pptr->prime));
  return a;
}

//compute a/R mod p using Montgomery multiplication (a,1)
long long int smallprimefield_convert_out(long long int a, const Prime_ptr* Pptr){
    return smallprimefield_mul(a, (long long int) 1,Pptr);
}


long long int smallprimefield_PrimitiveRootofUnity(long long int n, const Prime_ptr* Pptr){
    if ( (Pptr->prime - 1) % n != 0){
      printf("Error: %lld does not divide %lld - 1!\nCannot compute %lld-th primitive root of unity of %lld.\n",n,Pptr->prime,n,Pptr->prime);
      return 0;
    }
    long long int  q = (Pptr->prime - 1) / n;
    long long int c;
    long long int p1 = smallprimefield_convert_in((Pptr->prime-1),Pptr);
    int i = 0;
    long long int test = q * n / 2;
    srand (time(NULL));
    while(i < 20){
      c =  rand();
      c = smallprimefield_convert_in(c,Pptr);
      if (smallprimefield_exp(c,test,Pptr) == p1) {
          return smallprimefield_convert_out(smallprimefield_exp(c,q,Pptr),Pptr);
      }
      i ++;
    }      
    printf("No primitive root found!\nPlease try again.\n");
    return 0;
  
    // If no primitive root found
  }


///////////////////////////////////////////////////////////////////////////
//same functions as above for primes that is less than 32bits

int getRsquare32(int prime){
  
  int res;
   __asm__(
      "movl %%eax,%%esi\n\t"
      "xor %%eax,%%eax\n\t"
      "movl $1,%%edx\n\t"
      "divl %1\n\t"
      "movl %%edx,%%eax\n\t"
      "mull %%edx\n\t"
      "divl %1\n\t"
      "movl %%esi,%%eax\n\t"
      : "=&d" (res)
      : "rm"(prime)
      :"esi","edi");
  return res;
}

int getPp32 (int prime){
  //printf("getPp\n");
  int t, A, B, C, D, u, v,q;
  __asm__(
     // "movq %%rax,%%rsi\n\t"
      "xor %%eax,%%eax\n\t"
      "movl $1,%%edx\n\t"
      "divl %2\n\t"
     // "movq %%rsi,%%rax\n\t"
     // "movq %%rdx,%0\n\t"
     // "movq %%rax,%1\n\t"
      : "=&d" (v),"=&a" (q) 
      : "b"(prime)
      :"esi","edi");
  A = 1;
  B = 0;
  C = 0;
  D = 1;
 // q = R/prime;
  u = prime;
  //v = R - R mod prime;
  t = A;
  A = B;
  B = t - q * B;
  t = C;
  C = D;
  D = t - q * D;

  while (v != 0){
    q = u / v;
    t = u;
    u = v;
    v = t - q * v;
    t = A;
    A = B;
    B = t - q * B;
    t = C;
    C = D;
    D = t - q * D;
 //   printf("%llu,%llu\n",A,C);
  }
  if(C < 0){
    C = 0-C;
  }
  return C;
}
Prime_ptr32* smallprimefield_get_prime_constants32(int prime){
  Prime_ptr32* Pptr = (Prime_ptr32*)malloc(sizeof(Prime_ptr32));
  Pptr->prime = prime;
  Pptr->rsquare = getRsquare32(prime);
  Pptr->prime_inv = getPp32(prime);
  return Pptr;
}


int smallprimefield_add32(int a, int b, const Prime_ptr32* Pptr){
  a = a + b;
  a -= Pptr->prime;
  a += (a >> 31) & Pptr->prime;
  return a;
}

int smallprimefield_sub32(int a,int b, const Prime_ptr32* Pptr){
   a = a - b;
   a += (a >> 31) & Pptr->prime;
   return a;
 }


int smallprimefield_mul32(int a, int b, const Prime_ptr32* Pptr){
  __asm__(
      "mull %2\n\t"
      "movl %%eax,%%esi\n\t"
      "movl %%edx,%%edi\n\t"
      "imull %3,%%eax\n\t"
      "mull %4\n\t"
      "add %%esi,%%eax\n\t"
      "adc %%edi,%%edx\n\t"
      "subl %4,%%edx\n\t"
      "mov %%edx,%%eax\n\t"
      "sar $31,%%eax\n\t"
      "andl %4,%%eax\n\t"
      "addl %%eax,%%edx\n\t"
      : "=&d" (a)
      : "a"(a),"rm"(b),"b"(Pptr->prime_inv),"c"(Pptr->prime)
      :"esi","edi");
  return a;
}


int* pinverse32(int a, const Prime_ptr32* Pptr){
    int u = a;
    int v = Pptr->prime;
    int x2 = 0;
    int* r = (int *)malloc(2*sizeof(int));
    r[0] = 1;
    r[1] = 0;
    while(v > 0){
      if (v%2 == 0){
        v = v >> 1;
        r[0] = 2 * r[0];
      }
      else if (u%2 ==0){
        u = u >> 1;
        x2 = 2 * x2;
      }
      else if(v>=u){
        v = (v - u) >> 1;
        x2 = x2 + r[0];
        r[0] = 2 * r[0];
      }
      else{
        u = (u - v) >> 1;
        r[0] = x2 + r[0];
        x2 = 2 * x2;
      }
      r[1]++;
    }
    if(u != 1){
      // printf("Not invertible!\n");
      return NULL;
    }
    if(r[0] > Pptr->prime){
        r[0] -= Pptr->prime;    
    }
    return r;
}


int smallprimefield_inv32(int a, const Prime_ptr32* Pptr){
  if (a == 0)
  {
    // printf("Not invertible!\n");
    return 0;
  }
   // long long int rsquare =  getRsquare(prime);
    int* result = pinverse32(a, Pptr);
    if(result[1] < 32){
      result[0] = smallprimefield_mul32(result[0], Pptr->rsquare,Pptr);
      result[1] += 32;
    }
    result[0] = smallprimefield_mul32(result[0], Pptr->rsquare,Pptr);
    int tmp;
    if(result[1] != 32)
    {
      tmp = 1L << (64 - result[1]);
      result[0] = smallprimefield_mul32(result[0], tmp,Pptr);
    }
    tmp = result[0];
    free(result);
    return tmp;
}

int smallprimefield_div32(int a,int b, const Prime_ptr32* Pptr){
    return smallprimefield_mul32(a, smallprimefield_inv32(b, Pptr),Pptr);
}

int smallprimefield_convert_in32(int a, const Prime_ptr32* Pptr){
   __asm__(
      "xor %%eax,%%eax\n\t"
      "movl $1,%%edx\n\t"
      "idivl %2\n\t"
      "movl %%edx,%%eax\n\t"
      "imull %%ebx\n\t"
      "idivl %2\n\t"
      "mov %%edx,%%eax\n\t"
      "sar $31,%%eax\n\t"
      "andl %2,%%eax\n\t"
      "addl %%eax,%%edx\n\t"
      : "=&d" (a)
      : "b"(a),"c"(Pptr->prime));
  return a;
}

int smallprimefield_convert_out32(int a, const Prime_ptr32* Pptr){
    return smallprimefield_mul32(a, (int) 1,Pptr);
}

int smallprimefield_exp32(int a,int e, const Prime_ptr32* Pptr){
    int one = smallprimefield_convert_in32(1LL,Pptr);
   // printf("one:%lld\n",one);
    if(a == 0 || a == one || e == 1)
      return a;
    if(e == 0)
      return one;
    int neg_flag = 0;
    if(e < 0)
      neg_flag = 1;
    while (e != 0) {
        if (e % 2)
          one = smallprimefield_mul32(one,a,Pptr);
        a = smallprimefield_mul32(a,a,Pptr);
     //   printf("e: %lld  a: %lld   one:%lld\n",e,a,one);
        e >>= 1;
    }
    if(neg_flag ==1)
      one = smallprimefield_inv32(one,Pptr);
    return one;
}

int smallprimefield_PrimitiveRootofUnity32(int n, const Prime_ptr32* Pptr){
    if ( (Pptr->prime - 1) % n != 0){
      printf("Error: %lld does not divide %lld - 1!\nCannot compute %lld-th primitive root of unity of %lld.\n",n,Pptr->prime,n,Pptr->prime);
      return 0;
    }
    int  q = (Pptr->prime - 1) / n;
    int c;
    int p1 = smallprimefield_convert_in32((Pptr->prime-1),Pptr);
    int i = 0;
    int test = q * n / 2;
    srand (time(NULL));
    while(i < 20){
      c =  rand();
      c = smallprimefield_convert_in32(c,Pptr);
      if (smallprimefield_exp32(c,test,Pptr) == p1) {
          return smallprimefield_convert_out32(smallprimefield_exp32(c,q,Pptr),Pptr);
      }
      i ++;
    }      
    printf("No primitive root found!\nPlease try again.\n");
    return 0;
  
    // If no primitive root found
  }