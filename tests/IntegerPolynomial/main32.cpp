/**
 *Test file for Univariate Integer Dense Polynomial Multiplication.
 */

#include <bpas.h>
#include <cilktools/cilkview.h>

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

int degree = 1 << 2;
int bitCount = degree;
int type = 10;
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
 * Test the dense integer univariate polynomial arithmetic using the PolyExposed interface...
 */
void testExposedInterface(){
	DenseUnivariateIntegerPolynomial t1(4);
	DenseUnivariateIntegerPolynomial t2(4);

	t1.setCoefficient(0,0);
	t1.setCoefficient(1,0);
	t1.setCoefficient(2,3);
	t1.setCoefficient(3,1);

	t2.setCoefficient(0,8);
	t2.setCoefficient(1,8);
	t2.setCoefficient(2,3);
	t2.setCoefficient(3,1);

	cout << t1 << endl; 
	cout << t2 << endl;
	
	DenseUnivariateIntegerPolynomial t3 = t2 / t1;
	cout << t3 << endl;
}

/**
 * Generate random coefficients and store them in files.
 */
void generateData(int d, int bCount){

	std::ifstream ifile(A_FILE_NAME);

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
  std::cout << "#workers:\t" << __cilkrts_get_nworkers() << std::endl;
 
  // First generate the coefficients of the source polynomial.
  //generateData(MAX_DEGREE, bitCount);
  //generateData(degree, bitCount);
  // Get the arguments
  if(argc > 1) degree = atoi(argv[1]);
  if(argc > 2) bitCount = atoi(argv[2]);
  if(argc > 3) type = atoi(argv[3]);
  if(argc > 4) check = atoi(argv[4]);

//testExposedInterface();

  UnivariateIntegerPolynomial a(degree, bitCount);
  UnivariateIntegerPolynomial b(degree, bitCount); //check d1+d2-1 is a power of 2?
  UnivariateIntegerPolynomial c;

  a.generateRandom(); sleep(1); b.generateRandom();

  putDegreeToFile(degree);
  a.write(A_FILE_NAME);
  b.write(B_FILE_NAME);

  cilkview_data_t start_data, end_data;
  unsigned long long start = __cilkview_getticks();
  __cilkview_query(start_data);
  
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
    c = dnc.mulDnC(&a, &b);
    break;
  case TOOM_8:
    c = toom8.multiply(&a, &b);
    break;
  case SSA:
    c = ssa.multiply(&a, &b);
    break;
  case FINAL_POINTER:
    c.setSize(a.getSize() + b.getSize() - 1);
    m.multiply(c.getCoefficients(), a.getCoefficients(), b.getCoefficients(), a.getSize(), b.getSize());
    break;
  case FINAL_POINTER_RETURN:
    c.setSize(a.getSize() + b.getSize() - 1);
    c.setCoefficients(m.multiply(a.getCoefficients(), b.getCoefficients(), a.getSize(), b.getSize()));
    break;
  case FINAL:
    c = m.multiply(&a, &b);
    break;
  default:
    c = m.multiply(&a, &b);
    break;
}

  __cilkview_query(end_data);
  unsigned long long end = __cilkview_getticks();
  
  std::cout << degree <<", \t" <<bitCount<< ",\t" << (end - start) /1000.f << " (sec)"<< std::endl;

  //__cilkview_do_report(&start_data, &end_data, (char *)"mul-naive-serial", CV_REPORT_WRITE_TO_LOG | CV_REPORT_WRITE_TO_RESULTS);

  if(check == 1){
    c.write(C_FILE_NAME);
    UnivariateIntegerPolynomial cNaive;
    //cNaive = naive.multiply(&a, &b);
    cNaive = ks.multiply(&a, &b);
    //cNaive = toom8.multiply(&a, &b);
    cNaive.write(C_NAIVE_FILE_NAME);
    cNaive.freeHeap();
    std::cout << "Check Result: " << (system("diff -b c.dat c-naive.dat") == 0) << std::endl;
  }
  
  a.freeHeap(); b.freeHeap(); c.freeHeap();
  
  return 0;
}
