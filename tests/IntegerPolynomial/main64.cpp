/**
 *Test file for Univariate Integer Dense Polynomial Multiplication.
 */

#include <bpas.h>
#include "../../include/global.h"

using namespace std;

#define	NAIVE_SERIAL				0
#define BIG_INTEGER_SERIAL              	1
#define RECURSIVE_BIG_INTEGER_PARALLEL  	2
#define TOOM_4					3
#define DnC					4
#define TOOM_8					5
#define SSA					6
#define FINAL_POINTER                           8
#define FINAL_POINTER_RETURN			9
#define FINAL					10

const static char* A_FILE_NAME = "a.dat";
const static char* B_FILE_NAME = "b.dat";
const static char* C_FILE_NAME = "c.dat";
const static char* C_NAIVE_FILE_NAME = "c-naive.dat";

int degree = 1 << 10;
int bitCount = degree;
int type = 6;
int check = 0;

const static int MAX_DEGREE = 1 << 20;

MulNaive naive;
MulDnC dnc;
MulToom4 toom4;
MulToom8 toom8;
MulSSA ssa;
MulKS ks;
Multiplication m;


/**
 * Generate random coefficients and store them in files.
 */
void generateData(int d, int bCount){

	ifstream ifile(A_FILE_NAME);

	if(!ifile){
        	UnivariateIntegerPolynomial a(d, bCount);
	        UnivariateIntegerPolynomial b(d, bCount);
	        a.generateRandom(); sleep(1); b.generateRandom();
	
	        a.write(A_FILE_NAME); b.write(B_FILE_NAME);
	}
}

void putDegreeToFile(int d){
	std::ofstream ofile("degree.dat");
	std::ofstream ofs("degree.dat", std::ofstream::out);
        ofs << d;
        ofs.close();
}


int main(int argc, char *argv[]){
  cout << "#workers:\t" << __cilkrts_get_nworkers() << endl;
  /*
  unsigned long int i1 = 7881299347898369;  //min=-9223372036854775807
  unsigned long int i2 = 7881299347898369; //max=9223372036854775807
  unsigned __int128 i3 = i1*i2;
  unsigned __int128 i4 = i3<<64; 
  unsigned long int i5 = (unsigned long int)(i4 >> 64); //15762598695796737
  unsigned long int i6 = (unsigned long int)(i3);
  
  cout<<"(unsigned long int) i5: "<<i5<<endl;
  cout<<"(unsigned long int) i6: "<<i6<<endl;
  */

  /*small ToMPZ=====================
  mpz_t A;           
  mpz_init(A);                                 
     
  mpz_set_str(A, "37779450224004514418863840514288674912571", 10);
  //mpz_set_str(A, "16774051254337737770216744589354326591067",10);
  cout<<"---input\n";
  for (int i=0; i<A->_mp_size; ++i)
    cout<<(A->_mp_d)[i]<<", ";
  cout<<"\n\n";

  int K1 = 8;
  int M1 = 17;
  signed long limb_bits = 128;
  sfixn q[8] = {112955, (7881299347898369-82235), 97331, 0, 7881299347898369-3453, 80621, 0, 7881299347898369-56844};
  //18445071950640596677, 432345395153131471, 111
  //-37779318093080006086716886581262620575429
  //sfixn q[8] = {112955, (7881299347898369-82235), 97331, 22749, 7881299347898369-3453, 80621, 26057, 7881299347898369-56844};
  //4112366136533993157, 425182901531891661, 111
  //-37779185968393245498771601188919404938949
  //sfixn q[8] = {112955, 82235, 97331, 22749, 3453, 80621, 26057, 56844};
  //sfixn q[8] = {34395, 94566, 39992, 51068, 115289, 23609, 103178, 25238};
  
  cout<<"---qqq\n";
  for (int i=0; i<8; ++i)
    cout<<q[i]<<", ";
  cout<<"\n--result--\n";
  mpz_t R;
  ssa.ToMPZ(q, limb_bits, R, K1, M1); 
  cout<<"R->_mp_alloc,size " <<R->_mp_alloc <<", "<<(R->_mp_size)<<endl;;
  for (int i=0; i<abs(R->_mp_size); ++i)
    cout<<(R->_mp_d)[i]<<", ";
  cout<<"\n";
  end small ToMPZ==============================*/

  /*ToMPZ-------------------------
  //int K1 = 2;
  //int M1 = 5;
  //int num_limb = 2;
  //sfixn q[30] = {728, 54, 744, 352, 1120, 528, 804, 560, 1252, 576, 2054, 618, 1370, 1092, 2292, 1008, 1138, 932, 1528, 740, 1027, 518, 716, 638, 1127, 332, 32, 324, 468, 176};
  //2456, 12008, 18016, 18724, 19684, 21830, 36314, 34548, 30962, 25208, 17603, 21132, 11751, 10400, 6100
  //sfixn q[30] = {730, 54, 768, 352, 1216, 528, 960, 560, 1408, 576, 2206, 618, 1536, 1092, 2560, 1008, 1444, 932, 1732, 740, 1167, 518, 820, 638, 1257, 332, 144, 324, 500, 176};
  //2458, 12032, 18112, 18880, 19840, 21982, 36480, 34816, 31268, 25412, 17743, 21236, 11881, 10512, 6132

  //2457 12020 18064 18802 19762 21906 36397 34682 31115 25310 17673 21184 11816 10456 6116 //wrong
  //3481 24308 67216 98674 99634 99730 121389 171898 187787 129758 89353 74432 78376 67800 22500
  int K1 = 16;
  int M1 = 17;
  int limb_bits = 6*64;
  sfixn q[48] = {7881283497105298, 7881277190893621, 7881293476456361, 7881293728322103, 7881277537255091, 7881287600908287, 7881295741051664, 23115027548, 26443805266, 29448972606, 44692415442, 32344295388, 49936432213, 59718311056, 58251264327, 44537730058, 73185862148, 68178273538, 71427890344, 32552527156, 44243947722, 23305371160, 11367098588, 7881283796533189, 7881250795760051, 7881260374269833, 7881220540142669, 7881199207772975, 7881195504630801, 7881171020609007, 7881189024715861, 7881184696092711, 7881233976799334, 7881256025190073, 7881245434300482, 7881265283113387, 7881267337337323, 7881276868277993, 7881287916129952, 7881293006967177, 15966513420, 20625596260, 30224072345, 52180597746, 56290872505, 63041343604, 66075772793, 60438692852};

  cout<<"---qqq\n";
  for (int i=0; i<48; ++i)
    cout<<q[i]<<", ";
  cout<<"\n--result--\n";
  for (int i=0; i<3; ++i){
    mpz_t R;
    ssa.ToMPZ(&q[i*16], limb_bits, R, K1, M1); 
    cout<<"R->_mp_alloc,size " <<R->_mp_alloc <<", "<<(R->_mp_size)<<endl;;
    for (int i=0; i<abs(R->_mp_size); ++i)
      cout<<(R->_mp_d)[i]<<", ";
   }
  cout<<"\n\n";
  
  ToMPZ end ==============================*/

  /*---------------------------------------------
  mpz_t A;
  mpz_init(A);

  mpz_set_str(A, "17877077239923462799432033330941128169972487760789857836964748591629018706764729302402499064222537903586709144248648067805005177082595426595244865341445460", 10); //89478484, 44739242, 91625968640
  int M = 27;
  int K = 19;

  //mpz_set_str(A, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096", 10); //67108864*2^(27*18);
  //int M = 27;
  //int K = 19;

  //mpz_set_str(A, "115792089237316195423570985008687907853269984665640564039457584007913129639936", 10); //16* 
  //int M = 6;
  //int K = 43;
  
  //int d1 = 8192; //0.024s
  //int M = 18;
  //int K = 512;

  //int d1 = 16384; //0.093s
  //int M = 18;
  //int K = 1024; 
  
  //int d1 = 131073; //5.743
  //int M = 18;
  //int K = 8192;

  ---------------------------------------------*/
  /*---------------------------------------------
  //int d1 = 128;
  //int K = 8;
  //int M = 17;

  //B:=2^17;
  //C:=[112955, 82235, 97331, 22749, 3453, 80621, 26057, 12];
  //X:=eval(C[1]+C[2]*B+C[3]*B^2+C[4]*B^3+C[5]*B^4+C[6]*B^5+C[7]*B^6+C[8]*B^7);
  //v: 8107495780344974429259089362403440955

  //int d1 = 8192; //1-core 0.024s, 2-core 0.015s
  //int K = 512;
  //int M = 18;

  int d1 = 16384; //1-core 0.093s, 2-core 0.052s, 3-core 0.04s
  int M = 18;
  int K = 1024; 

  //int d1 = 131072; // /8 = 16384 //1-core 5.825s, 2-core 3.01s, 3-core 2.095s, 4-core 1.57s, 8-core 0.881s, 
  //int M = 18;
  //int K = 8192; 
  
  //int d1 = 6;
  //int K = 2;
  //int M = 5;

  sfixn prime = 7881299347898369;

  gmp_randclass rr(gmp_randinit_default);
  rr.seed(1);

  mpz_class v = -rr.get_z_bits(d1);
  //cout<<"v: "<<v<<endl;

  sfixn *X = (sfixn *)my_calloc(K, sizeof(sfixn));

  cilkview_data_t start_data, end_data;
  unsigned long long start = __cilkview_getticks();
  __cilkview_query(start_data);
  //setup grain size, need to tune
  //#pragma cilk_grainsize = 8192;
  cilk_for (int i=0; i<d1; ++i)
  //mpzToPoly(A, M, X);
    ssa.mpzToPoly(v.get_mpz_t(), M, prime, X);

  __cilkview_query(end_data);
  unsigned long long end = __cilkview_getticks();

  cout << (end - start) /1.f << " (ms)"<<endl;

  //for (int i=0; i<K; ++i)
  //   cout<<X[i]<<", ";
  //cout<<"\n\n";

  my_free(X);

  return 0;
  ---------------------------------------------*/

  ///*main----------------------------------------------
  
 
  // First generate the coefficients of the source polynomial.
  //generateData(MAX_DEGREE, bitCount);
  //generateData(degree, bitCount);
  // Get the arguments
  if(argc > 1) degree = atoi(argv[1]);
  if(argc > 2) bitCount = atoi(argv[2]);
  if(argc > 3) type = atoi(argv[3]);
  if(argc > 4) check = atoi(argv[4]);
  
  UnivariateIntegerPolynomial a(degree, bitCount);
  UnivariateIntegerPolynomial b(degree, bitCount); //check d1+d2-1 is a power of 2?
  UnivariateIntegerPolynomial c;
  
  a.generateRandom(); sleep(1); b.generateRandom();
  //a.print();
  #if TDEBUG
	std::cout<<"a and b generated"<<std::endl;
  #endif
//if(check==1){
  putDegreeToFile(degree);
  a.write(A_FILE_NAME);
  b.write(B_FILE_NAME);
//}


  // cilkview_data_t start_data, end_data;
  unsigned long long start;
  startTimer(&start);
  // __cilkview_query(start_data);
  
  switch(type){
  case NAIVE_SERIAL:
    c = naive.multiply(&a, &b);
    break;
  case BIG_INTEGER_SERIAL:
    c = ks.multiply(&a, &b);
    break;
  case RECURSIVE_BIG_INTEGER_PARALLEL:
    c = dnc.multiply(&a, &b);
    break;
  case TOOM_4:
    c = toom4.multiply(&a, &b);
    break;
  case DnC:
    c = dnc.multiply(&a, &b);
    break;
  case TOOM_8:
    c = toom8.multiply(&a, &b);
    break;
  case SSA:
    c = ssa.multiply(&a, &b);
    break;
  case FINAL_POINTER:
    m.multiply(c.getCoefficients(), a.getCoefficients(), b.getCoefficients(), a.getSize(), b.getSize());
    break;
  case FINAL_POINTER_RETURN:
    c.setCoefficients(m.multiply(a.getCoefficients(), b.getCoefficients(), a.getSize(), b.getSize()));
    break;
  case FINAL:
    c = m.multiply(&a, &b);
    break;
  default:
    c = m.multiply(&a, &b);
    break;
	}
  
#if BPASDEBUG
std::cout<<"multiply done"<<std::endl;
#endif
  // __cilkview_query(end_data);
  float elapsed;
  stopTimer(&start, &elapsed);

  cout << degree <<", \t" <<bitCount<< ",\t" << elapsed << " (sec)"<<endl;

  //__cilkview_do_report(&start_data, &end_data, (char *)"mul-naive-serial", CV_REPORT_WRITE_TO_LOG | CV_REPORT_WRITE_TO_RESULTS);
  
  if(check == 1){
    c.write(C_FILE_NAME);
    UnivariateIntegerPolynomial cNaive;
    //cNaive = naive.multiply(&a, &b);
    cNaive = ks.multiply(&a, &b);
    //cNaive = toom8.multiply(&a, &b);
    cNaive.write(C_NAIVE_FILE_NAME);
    cNaive.freeHeap();
    if (system("diff -b c.dat c-naive.dat") == 0) {
	cout << "Pass." << endl;
    }
  }
  
#if BPASDEBUG
std::cout<<"free a and b"<<std::endl;
#endif
  a.freeHeap(); b.freeHeap();
#if BPASDEBUG
std::cout<<"free c"<<std::endl;
#endif
c.freeHeap();
  
  return 0;
  //--------------------------------------*/
}
