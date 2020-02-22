//#define MontSmallPrimeField 
#include <bpas.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
//#include "../../include/FiniteFields/SmallPrimeField_Support.h"

using namespace std;
//#if ( MontSmallPrimeField )
//R = 2^32;
//long long int SmallPrimeField::R = 4294967296;
//Pp = -p^-1 mod R;
//long long int SmallPrimeField::Pp = 883949567;
//N = 2^30;
//long long int SmallPrimeField::N = 1073741824;
//#endif
extern DistributedDenseMultivariateModularPolynomial<RationalNumber> randomDDMMP(int, sfixn);
extern SFDDMMP randomSFDDMMP(int, sfixn, int);
extern void testDDMMP(int, sfixn);
extern void testSFDDMMP(int, sfixn);
extern void smallprimefieldtests();
extern void smallprimefieldsupporttests();
extern void bigprimefieldtests();
extern void generalizedfermatprimefieldtests();
extern void GFPF(int k);




int main(int argc, char *argv[]) {
	sfixn p = 469762049; // 4179340454199820289; 
	int var = 2;
	if (argc > 1) { var = atoi(argv[1]); }

	srand (time(NULL));
	

	testDDMMP(var, p);
	testSFDDMMP(var, p);
	smallprimefieldtests();
	smallprimefieldsupporttests();
	bigprimefieldtests();
	generalizedfermatprimefieldtests();
 	GFPF(64);
	
	return 0;
}

void testDDMMP(int var, sfixn p) {
	//string names[] = {"z", "y", "x"};
	vector<Symbol> names;
	names.push_back(Symbol("z"));
	names.push_back(Symbol("y"));
	//names.push_back("x");

	Field<RationalNumber>::DDMMP f = randomDDMMP(var, p);
	f.setRingVariables(names);
	//cout << "f := " << f << endl;
	DistributedDenseMultivariateModularPolynomial<RationalNumber> g = randomDDMMP(var, p);
	g.setRingVariables(names);
	//cout << "g := " << g << endl;

	DistributedDenseMultivariateModularPolynomial<RationalNumber> r = f + g;
	//cout << "f + g: " << r << endl;
	cout << "addition pass." << endl;

	if (r == f) { cout << "comparison pass." << endl; }

	r = f - g;
	//cout << "f - g: " << r << endl;
	cout << "subtraction pass." << endl;

	//cout << "y^" << r.degree("y") << std::endl;

	r = f * g;
	//cout << "f * g: " << r << endl;
	cout << "multiplication pass." << endl;

	ofstream inputfile("input.dat", ofstream::out);
	inputfile << "f := " << f << ";\n";
	inputfile << "g := " << g << ";\n";
	inputfile << "r := " << r << ";\n";
	// inputfile.close();
}
//Main: make dev
void testSFDDMMP(int var, sfixn p) {
	vector<Symbol> names;
	names.push_back(Symbol("z"));
	names.push_back(Symbol("y"));
	//names.push_back("x");
	//To test DMPMul, make d1+d2+1 >=1024 and a power of 2
	int d1 = 16, d2 = d1 - 1;
	//When d1 > 16 FFT multiplication will fail on gcc-6 but works for gcc-4.8
  SFDDMMP f(p), g(p), r(p); 
	for (int i = 0; i < 3; i++) {
		f = randomSFDDMMP(var, p, d1);
		f.setRingVariables(names);
		//cout << "f := " << f << endl;
  	g = randomSFDDMMP(var, p, d2);
		g.setRingVariables(names);
		//cout << "g := " << g << endl;
		if (!i) {
			r = f + g;
			//cout << "f + g: " << r << endl;
			cout << "addition pass." << endl;
			r = f - g;
			//cout << "f - g: " << r << endl;
			cout << "subtraction pass." << endl;
		}
		r = f * g;
		if (!i) { cout << "FFT-based (two-convolution) "; }
		else if (i == 1) { cout << "FFT-based (KS) "; }
		else if (i == 2) { cout << "Plain "; }
		//cout << "f * g: " << r << endl;
		cout << "multiplication pass." << endl;

		d1 >>= i + 1;
		d2 = d1 - 1;

		if (!i) {
		        ofstream inputfile("input.dat", ofstream::out);
			inputfile << "f := " << f << ";\n";
			inputfile << "g := " << g << ";\n";
			inputfile << "r := " << r << ";\n";
			inputfile.close();
		}
	}
}

DistributedDenseMultivariateModularPolynomial<RationalNumber> randomDDMMP(int var, sfixn p) {
	int* ds = new int[var];
	for (int i = 0; i < var; ++i) {
	  ds[i] = rand() % var + 1;
	}

	DistributedDenseMultivariateModularPolynomial<RationalNumber> f(var, ds, p);
	for (int i = 0; i < f.size(); ++i)
		f.setCoefficient(i, (rand()%p));
	delete [] ds;
	return f;
}

SFDDMMP randomSFDDMMP(int var, sfixn p, int d) {
	int* ds = new int[var];
	for (int i = 0; i < var; ++i)
	  ds[i] = d; //rand() % var + 1;

	SFDDMMP f(var, ds, p);
	for (int i = 0; i < f.size(); ++i)
		f.setCoefficient(i, 1);
		//f.setCoefficient(i, (rand()%11));
	delete [] ds;
	return f;
}

// void mul(uint32_t x,uint32_t y,uint32_t *hi, uint32_t *lo){
//     uint64_t res = (uint64_t) x * (uint64_t) y;
//     *lo = (uint32_t)res;
//     *hi = (uint32_t)(res >> 32);
// }

// long long int multi(long long int a, long long int b, long long int prime, long long int R, long long int Pp){
//   uint32_t lo, hi;
//     mul(a, b, &hi, &lo);
//     printf("a;%d,b:%d,hi:%lu,lo:%lu\n",a,b,hi,lo);
//     //  long long int x = (b * a);
//       long long int w = lo*Pp;
//       printf("w:%d\n",w);
//      uint32_t x = w&(R-1);
//      printf("x:%d\n",x);
//       uint32_t lo1, hi1;
//      mul(prime, x, &hi1, &lo1);
//      printf("prime;%d,x:%d,hi:%lu,lo:%lu\n",prime,x,hi1,lo1);
//      // printf("w = %lld\n",w);
//       //long long int y = ;
//     //  printf("y = %lld\n",y);
//       long long int z =  hi + hi1;
//       if(lo1 >= R - lo1)
//         z++;
//      // printf("z = %lld\n",z);
//       if(z >= prime){
//         z = z - prime;
//       }
//       return z;
// }

void smallprimefieldtests() {
  cout << "***************************************************************************" <<endl;
  cout<<"SmallPrimeField test"<<endl;
  clock_t start,end;
  start = clock();
 SmallPrimeField::setPrime(4179340454199820289);

  int n = 234456778;
  int i = 6850867;
  int j = 40191135;
  int k = -988;
  int g = -935;
 
  SmallPrimeField a(n);
   SmallPrimeField b(i);
   SmallPrimeField c(j);
   SmallPrimeField d(a);
   SmallPrimeField zero(0);

//    cout << a.number() << " " << b.number() << endl;

// //  b*= c;
// //  cout << "b:" << b << endl;
//  // cout << a.getint() << endl;

// //  cout << a.Mont(b)<<endl;

 long long int p = a.Prime();
   cout << "Prime: " << p << endl;
   cout <<"Root: "<< SmallPrimeField::findPrimitiveRootofUnity(32) << endl;
// cout << n << " " << i << endl;
//    cout << a.number() << " " << b.number() << endl;
   a=a+b;
//   cout << a.number() << endl;
   c+=d;
//   cout << c.number() << endl;
   d=d+k;
//   cout << d.number() << endl;
  if(a==(n + i)%p  && c == (j + n)%p && d == (n + k)%p)
    cout << "Addition pass." << endl;
  else{
    cout << "Addition error." << endl;
    exit(0);
  }

//cout << n << " " << i << endl;
  a = n;
  b = i;
  //cout << a.number() << " " << b.number() << endl;
  a = a-b;
  
  if ( a.number() == (n - i)%p){

    cout << "Subtraction pass." << endl;
  }
  else{
    cout << a.number() << " " << (n-i)%p<<endl;
    cout << "Subtraction error." << endl;
    exit(0);
  }
//   //  cout<< a <<" "<< b <<" "<< c <<" "<< d << endl;

  a = n;
  b = 7;
  a = a * b;

  if ( a.number() == (n * 7)%p)
    cout << "Multiplication pass." << endl;
  else{
    cout << "Multiplication error." << endl;
    exit(0);
  }

   a = b;
   b = a.inverse();
   c = a.inverse2();

   b = b - c;

  if(b.isZero()){
    cout << "Inversion pass." << endl;
  }
  else{
    cout << "Inversion error." << endl;
  }

//   //cout<< a <<" "<< b <<" "<< c <<" "<< d << endl;
  b = i;
  c = 540191135;
 // cout << b << " " << c << endl;
 // c = c.inverse();
  b /= c;
 // cout << b << endl;
  c = c.inverse();
  //cout<< a <<" "<< b <<" "<< c <<" "<< d << endl;
  if(b == (c * i)){
    cout << "Division pass." << endl;
  }
  else{
    cout << "Division error." << endl;
    exit(1);
  }

  
  a = n;
  b = 1;
  for(int index = 0; index < 5; index ++){
    b = (b * a);
  }


  a = a^5;

  if (a == b)
    cout << "Exponentiation pass."<< endl;
  else
    cout <<"Exponentiation error." << endl;
 // cout << c << endl;



  for (int in = 0; in < 100000000; in ++){
    SmallPrimeField test(988);
    a+=test;
  }

  end = clock();

  float ti = (float(end - start))/CLOCKS_PER_SEC;
  cout << ti << endl;

}


void smallprimefieldsupporttests(){
  printf("************************************\n");
  printf("Small Prime Field C support test\n");
  long long int prime = 4179340454199820289;
 // long long int R = 4294967296;
  Prime_ptr* ptr;
  ptr = smallprimefield_get_prime_constants(prime);
  long long int a,b,c,d,e;
  a = c = 40191135;
  b = d = 23768;
  printf("Rsquare: %lld\n",ptr->rsquare);
  //long long int Pp = 5000824070590443489;
  //smallprimefield_getPp(prime,R);
  //long long int Pp = 883949567;
  //cout << Pp << endl;
  // d = 1;
  // a =  smallprimefield_convert_in(a, ptr);
  // printf("%lld\n",a );
  // a = smallprimefield_exp(a,19,ptr);

  // printf("%lld\n",a );
  // a = smallprimefield_convert_out(a,ptr);

  // printf("%lld\n",a );
  // // for(int i = 0; i < 19; i ++)
  // //   d = (d * c)%prime;
  //  //printf("%lld\n",b );

  // // if(b == d)
  // //   printf("Exponentiation pass!\n");
  // // else
  // //   printf("Exponentiation Error!\n");
  a = smallprimefield_PrimitiveRootofUnity(8,ptr);
  printf("%lld\n",a);
  return;



  a = c = 40191135;
  b = d = 23768;
  a =  smallprimefield_convert_in(a, ptr);
//  cout << a << endl;
  b =  smallprimefield_convert_in(b,ptr);
  a = smallprimefield_add(a,b,ptr);
  //cout << a << endl;
  c = (c + d)%prime;
  //cout << c << endl;  
  e = smallprimefield_convert_out(a,ptr);
  //cout << e << endl;
  if( e == c)
    printf("Addition pass!\n");
  else
    printf("Addition error!\n");

  a = smallprimefield_sub(a,b,ptr);
  c = (c - d)%prime;  
  //cout << a << endl;
  //cout << c << endl;
  if(smallprimefield_convert_out(a,ptr) == c)
    printf("Subtraction pass!\n");
  else
    printf("Subtraction error!\n");

  a = smallprimefield_mul(a,b,ptr);
  c = (c * d)%prime;  
  //cout << a << endl;
  //cout << c << endl;
  if(smallprimefield_convert_out(a,ptr) == c)
    printf("Multiplication pass!\n");
  else
    printf("Multiplication error!\n");
  int j = 0;
  for(long long int i = 1; i <=10; i ++){
  //int i = 50196;
      a = i;
     c = i;
      a = smallprimefield_convert_in(a, ptr);
      c = smallprimefield_convert_in(c, ptr);
     a = smallprimefield_inv(a,ptr);
  //   a  = smallprimefield_convert_out(a,ptr);
     c = smallprimefield_convert_out(smallprimefield_mul(a,c,ptr),ptr); 
  //cout << a << endl;
  //cout << c << endl;
     if(c != 1){
        printf("c=%llu\n",c);
        j++;
        printf("a = %llu Inverse error!\n",i);
        //exit(0);
     }
  }
  printf("Inverse pass!\n");
  a = c = 40191135;
  b = d = 23768;
  a = smallprimefield_convert_in(a, ptr);
  b = smallprimefield_convert_in(b, ptr);
  c = smallprimefield_convert_in(c, ptr);
  d = smallprimefield_convert_in(d, ptr);
  a = smallprimefield_div(a,b,ptr);
  c = smallprimefield_mul(c,smallprimefield_inv(d,ptr),ptr);
  if(a != c ){
    printf("%llu %llu\n",a,c );
    printf("Division error!\n");
    exit(0);
  }
  else{
    printf("Division pass!\n");
  }

}




//mpz_class BigPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
void bigprimefieldtests() {
	cout << "***************************************************************************" <<endl;
  cout<<"BigPrimeField test"<<endl;
   mpz_class p("1091938819714247251320658718168450543670093214387642108906879549792740419068019517845067610696781282694278022873069240300680386317836693183540970830186786985787741723867145505541444129201020462084275294019814069366712462188994638523911468743310660808964250718198710337237688439001846059037619156172139078604354610041219492061864257864170553024248804863147877084098575702138952132003490519810700821020037345977895788329970443359476762062483048470722575265401075008734152928478394567658083649017991680414112738374046589906005272839526139298949026987916016608888246319627688473880291658694657",10);
   // mpz_class p("257",10);
  BigPrimeField::setPrime(p);
  mpz_class n("123543654879877656453342",10);
  mpz_class i("31249817516293874612350493845678909238475269347652938746529874562938746523475023794875239746928350827308497748258680430287450832704587203845720384572084357203294857203489572571203984563798417",10);
  mpz_class j("24345678987656",10);
  mpz_class k("-9889897478965634588",10);
  mpz_class g("-932457658764523434665",10);
  Integer zz(j);
 // cout << zz << endl;
	BigPrimeField a(3);
	BigPrimeField b(j);
  BigPrimeField c(k);
  BigPrimeField d(n);
  BigPrimeField e(zz);
 // mpz_class p = a.Prime();
 // cout << e << endl;
//  cout <<"Root:"<< BigPrimeField::findPrimitiveRootofUnity(8) << endl;
  //	cout<< a <<" "<< b<<" "<<c<<" "<<d<< endl;

  a=a+b;
 // d=d+k;
  if(a == (3+j)%p)
      cout << "Addition pass." << endl;
  else{
    cout << "Addition error." << endl;
    exit(0);
  }
 // cout<< a <<" "<< b<<" "<<c<<" "<<d<< endl;
  a=a-b;
  c-=d;
  if(a == 3 && c == (k-n)%p)
      cout << "Subtraction pass." << endl;
  else{
    cout << "Subtraction error." << endl;
    exit(0);
  }


  a = 3;
  b = j;
  c = k;
  d = n;

  a=a*b;
  c*=d;
  d*=g;
  if(a != (3*j)%p){
      cout << "Multiplication error." << endl;
      exit(0);
    }
  else if(c != (k*n)%p){
      cout << "Multiplication error." << endl;
      exit(0);
    }
  else if(d != (n*g)%p){
      cout << "Multiplication error." << endl;
      exit(0);
    }
  else
    cout << "Multiplication pass." << endl;

  b = n;
  b = b.inverse();
  mpz_t n_mpz;
  mpz_init (n_mpz);
  mpz_set(n_mpz,n.get_mpz_t());
  mpz_t p_mpz;
  mpz_init (p_mpz);
  mpz_set(p_mpz,p.get_mpz_t());
  mpz_invert (n_mpz, n_mpz, p_mpz);

  mpz_class n_inv(n_mpz);

  if(b == n_inv)
      cout << "Inversion pass." << endl;
  else{
    cout << "Inversion error." << endl;
    exit(0);
  }



  b = n;
  c = i;

  c = c / b;
  if(c == (i*n_inv)%p)
      cout << "Division pass." << endl;
  else{
    cout << "Division error." << endl;
    exit(0);  
  }

  d = k;
  i = 1;
  for(int index = 0; index < 7; index ++){
    i = (i * k)%p;
  }


  d=d^7;
 // cout << d << endl;
  if (d == i)
      cout << "Exponentiation pass." << endl;
  else{
    cout << "Exponentiation error." << endl;
    exit(0);
  }
//  cout<< a <<" "<< b <<" "<< c <<" "<< d << endl;
// BigPrimeField a1(i);
//   clock_t start,end;
//   start = clock();
// cout << a1 << endl;
//   for (int in = 0; in < 1000001; in ++){
//     a1*=a1;
//   }

//   end = clock();

//   float ti = (float(end - start))/CLOCKS_PER_SEC;
//   cout << a1 << endl;
//   cout << ti << endl;

  //   start = clock();
  // for (int in = 0; in < 100001; in ++){
  //   a=a+a;
  //  // a*=a;
  // }

  // end = clock();

  //  ti = (float(end - start))/CLOCKS_PER_SEC;
  // cout << a << endl;
  // cout << ti << endl;


}

//mpz_class GeneralizedFermatPrimeField::prime ("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
//unsigned long long int GeneralizedFermatPrimeField::r = 9223372054034644992ULL;
//int GeneralizedFermatPrimeField::k = 8;
//1000000000000000000000000000010000000000000000000000000000000000
void generalizedfermatprimefieldtests() {
  cout << "***************************************************************************" <<endl;
  cout<<"GeneralizedFermatPrimeField test"<<endl;
 // mpz_class p("52374250506775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
  //prime3
  //mpz_class p("41855814947416230160905824077102044107723669195743478941314613188961602997492974631771080284897031989884853955219823218284507222180203198418365663867454557929637857661786690253443410218509127538830031125593957763469901695303351030684228602570766096943274989297544663448958866408079360000000000000001",10);
  //prime4
 // mpz_class p("12194330274682935332610851851342563535752111249347079401277908877457504387594031661749056332752689482214279944317440318635455211497264750002177",10);
  //P8 for convolution implementation

 // mpz_class p("2269007733883464950444098216005816624661703879071974046166685917801868624984623568198674748106441557286217931676692278525983582602185970000241984176998213779951672720166897400277217166581802610820158765053097766381266070359787232937484689528129877885327166019946625216859944779777",10);
 //P16 for convolution implementation
  mpz_class p("279095111887779954448924290726583779438089639280034152600808727329141991929077956242494036568581164294946162439815234343375631006424021366822125241124173924805843708228196553364537992087680696670599673528986371229408625593683913121298725206630323629795781183113001619569754604853844077242993446043054578793865535056795673454357640364516367328989624305428281581356163174496495735111702305793866662480254936855935014275821099766294397933267790991891349993522918798914291606459279972211454820578951014801309903134572795051967533444810152083457",10);
//P32 for convolution implementation

  // long long int r = 4611686087146864640ULL; //r for prime 4
 // long long int r = 9223372054034644992ULL; // r for prime 3
 // long long int r = 4683743612465315840ULL; // r for prime 5
  //long long int r = 576460752303489024ULL; //R8
  //long long int r = 288230376151712768ULL; //R16
 long long int r = 72057594040025088ULL; //R32

 // int K = 8; // pirme 3
  int K = 32; // prime 5
 // int K = 16; // prime 4
  GeneralizedFermatPrimeField::setPrime(p,r,K);
  mpz_class n("123",10);
  mpz_class i("26187125253387706293540091008842954506639669630097560675975923979393277866127545231347033330913504906656138429677493633359612409895490708092711084228608",10);
  mpz_class j("24567791472635936451027341765192693615",10);
  mpz_class k("25941838401945782284278425092375",10);
  mpz_class g("7237005631252155698924907988729572412395441007692723318311900648330647044096",10);
  mpz_class l("6775412587080182017685909013279339260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457218",10);
  mpz_class m("12345678912345678912345678912345678939260195121351951847958786555732255090462694066661827009813312276859354987266719224819790981416185422168457217",10);
  mpz_t r_mpz;
  mpz_init (r_mpz);
  mpz_import(r_mpz, 1, -1, sizeof r, 0, 0, &r);
  mpz_class r2(r_mpz);

  GeneralizedFermatPrimeField a(3);
 // cout << a.findPrimitiveRootofUnity_plain(128) << endl;
  cout << GeneralizedFermatPrimeField::findPrimitiveRootofUnity(64) << endl;
  GeneralizedFermatPrimeField b(i);
  GeneralizedFermatPrimeField c(k);
  GeneralizedFermatPrimeField d(-1);
  GeneralizedFermatPrimeField e(l);
  GeneralizedFermatPrimeField f(g);
  GeneralizedFermatPrimeField r1(r2);
  a = b+b;

 // cout << (a.r - 9223372054034644992ULL==0) <<endl;
  if (a.number() - ((i + i)%p) == 0){
    cout << "Addition pass." << endl;
  }
  else{
    cout << "Addition error." << endl;
    exit(0);
  }
  

  //c=-d;
  //cout << c << endl;
  a = c-b;
  if( a == (k-i)%p)
  //mpz_class("26187125253387706293540091008842954506639669630097560675975923979393277866127545231347033330913504906656138429677493633359612409895490708092711084228610",10))
    cout << "subtraction pass." << endl;
  else{
    cout << "Subtraction error." << endl;
    exit(0);
  }
 // cout << "a = " << a << endl;
 // cout << "b = " << b << endl;
 // cout << "c = " << c << endl;
//  cout << "d = " << d << endl;
  a = c*b;
  // cout << (k*i)%p << endl;
  // cout << a << endl;
  if( a == (k*i)%p)
    //mpz_class("26187125253387706293540091008842954506639669630097560675975923979393277866127545231347033330913504906656138429677493633359612409895490708092711084228610",10))
    cout << "Multiplication pass." << endl;
  else{
    cout << "Multiplication error." << endl;
    exit(0);
  }

  for(int index = 0; index < 1; index ++){
    l = (l*r2)%p;
  }
  a = e.MulPowR(1);
  if(a == l)
      cout << "Multiplication by power of R pass." << endl;
  else{
    cout << "Multiplication by power of R error." << endl;
    exit(0);
  }
  //cout << "a = " << a << endl;
  a = i;
  a = a.inverse();
  mpz_t i_mpz;
  mpz_init (i_mpz);
  mpz_set(i_mpz,i.get_mpz_t());
  mpz_t p_mpz;
  mpz_init (p_mpz);
  mpz_set(p_mpz,p.get_mpz_t());
  mpz_invert (i_mpz, i_mpz, p_mpz);

  mpz_class i_inv(i_mpz);
  if(a == i_inv)
      cout << "Inversion pass." << endl;
  else{
    cout << "Inversion error." << endl;
    exit(0);
  }
 // cout << "a.inverse = " << a << endl;

//cout << a << endl;
//cout << c << endl;
  // a = i;
  // b = c / a;
  // cout << b << endl;
  // cout << (k*i_inv)%p << endl;
  // if(b == (k*i_inv)%p)
  //    cout << "Division pass." << endl;
  // else{
  //   cout << "Division error." << endl;
  //   exit(0);
  // }


  a = i;
  n = 1;
  for(int index = 0; index < 7; index++){
    n = (n*i)%p;
  }
  c = a^7;
  if(c == n)
     cout << "Exponentiation pass." << endl;
  else{
    cout << "Exponentiation error." << endl;
    exit(0);
  }
  
}

void GFPF(int k){

  mpz_t p_zz;
  mpz_init (p_zz);
  srgfn_prime * srgfn_ptr;


  switch (k)
    {
    // case 2:
    //   srgfn_ptr = p_list_1;
    //   break;

    case 4:
      srgfn_ptr = &p_list_2[0];
      break;

    case 8:
      srgfn_ptr = &p_list_3[0];
      break;

    case 16:
      srgfn_ptr = &p_list_4[0];
      break;

    case 32:
      srgfn_ptr = &p_list_5[0];
      break;

    case 64:
      srgfn_ptr = &p_list_6[0];
      break;

    case 128:
      srgfn_ptr = &p_list_7[0];
      break;

    default:
      {
  printf ("[ERROR: proper k is not specified!]\n");
  exit (EXIT_FAILURE);
      }
    }

 

//while(srgfn_ptr != NULL){
   srgfn_prime prime;

 //  //prime.k = k;
 prime.k = srgfn_ptr->k;
 prime.radix = srgfn_ptr->radix;
  gmp_compute_srgfn_u64 (p_zz, prime.radix, prime.k);

    usfixn64 conv_p1, conv_p2;
  get_proper_conv_p1p2 (&conv_p1, &conv_p2, prime);

  //conv_p1 = 4179340454199820289;
 // conv_p2 = 2485986994308513793;
  mpz_t p1p2_zz, r_zz;
  mpz_inits (p1p2_zz, r_zz, NULL);
  mpz_set_u64 (p1p2_zz, conv_p1);
  mpz_set_u64 (r_zz, conv_p2);
  mpz_mul (p1p2_zz, p1p2_zz, r_zz);
  mpz_set_u64 (r_zz, prime.radix);
  mpz_pow_ui (r_zz, r_zz, 2);
  mpz_t k_zz;
  mpz_init (k_zz);
  mpz_set_u64 (k_zz, prime.k);
  mpz_mul (r_zz, r_zz, k_zz);



//  if (((prime.k*prime.radix*prime.radix)< conv_p2) && ((prime.k*prime.radix*prime.radix)< conv_p1) &&
  if ((mpz_cmp (p1p2_zz, r_zz) <= 0))
    {
      printf (
    "[ERROR: conv primes must be twice the max value of the radix of SRGFN]\n");
      exit (EXIT_FAILURE);
    }



  init_fft_based_bigint_mult (conv_p1, conv_p2, k);
  init_gfpf_mult_data (&t_crt_data_global, &t_lhc_data_global, prime.radix, conv_p1, conv_p2);


  usfixn64* x = (usfixn64*) malloc(k * sizeof(usfixn64));

  usfixn64* y = (usfixn64*) malloc(k * sizeof(usfixn64));
  usfixn64* x_plain = (usfixn64*) malloc(k * sizeof(usfixn64));

  usfixn64* y_plain = (usfixn64*) malloc(k * sizeof(usfixn64));
  srand(time(NULL));
  
  for (int i = 0; i < k; i++){
      x[i] =  rand()%prime.radix;
      y[i] = rand()%prime.radix; 
      x_plain[i] = x[i];
      y_plain[i] = y[i];
  }
  //printf("radix %llu\n",prime.radix);
  mpz_t a_mpz;
  mpz_init (a_mpz);
  mpz_t b_mpz;
  mpz_init (b_mpz);
  usfixn64 radix = prime.radix;
  u64vector_to_radix_based_bigint_gmp (a_mpz, x, radix, k);
  u64vector_to_radix_based_bigint_gmp (b_mpz, y, radix, k);

  //cpu_timer gmp;
  //timer_record_start(&gmp);
  for(int i = 0; i < 1; ++i){//1000000
  //  BigPrimeField_Multiplication(a_mpz,b_mpz,p_zz,a_mpz);
    mpz_mul (a_mpz, a_mpz,b_mpz);
    mpz_mod (a_mpz, a_mpz, p_zz);
  }
  //timer_record_stop(&gmp);
  //timer_get_elapsed_time(&gmp, "gmp",1);
//  timer_get_elapsed_time(&gmp, "gmp_avg",1000000);

  //cpu_timer fft_based;
  //int n_chunks = 4;
  //usfixn64 * memory = (usfixn64*) malloc (n_chunks * sizeof(usfixn64) * k);
  //timer_record_start(&fft_based);
  for(int i = 0; i < 1; ++i){
   // memset(memory, 0x00, n_chunks*k*sizeof(usfixn64));
      GFPFMultiplication (x, y, k, &t_crt_data_global, &t_lhc_data_global );
  }
  
  //timer_record_stop(&fft_based);
  //timer_get_elapsed_time(&fft_based, "fft_based",1);

  mpz_t a_zz;
  mpz_init (a_zz);
  mpz_t b_zz;
  mpz_init (b_zz);

  // cpu_timer gmp_based;
  // timer_record_start(&gmp_based);


  for(int i = 0; i < 1; ++i){
  //  plain_mult_gmp_big_elements(x_plain, y_plain, k, R, p0);
    u64vector_to_radix_based_bigint_gmp (a_zz, x_plain, radix, k);
      //B=bigint(b)
    u64vector_to_radix_based_bigint_gmp (b_zz, y_plain, radix, k);

      mpz_mul (a_zz, a_zz, b_zz);
      mpz_mod (a_zz, a_zz, p_zz);

    bigint_to_u64vector_radixbased_gmp (x_plain, a_zz, radix, k);

  }
  // timer_record_stop(&gmp_based);
  // timer_get_elapsed_time(&gmp_based, "gmp_based",1);
//  timer_get_elapsed_time(&gmp_based, "gmp_based_avg",1000000);

 u64vector_to_radix_based_bigint_gmp (b_mpz, x, prime.radix, k);
 // radix_based_to_mpz_s64(b_mpz, x, prime.radix,k);
  if(mpz_cmp(a_mpz,b_mpz) == 0)
    printf("GFPFMultiplication: Pass\n");
  else{
    printf("GFPFMultiplication: Error\n");
  }
 u64vector_to_radix_based_bigint_gmp (b_mpz, x_plain, prime.radix, k);
 // radix_based_to_mpz_s64(b_mpz, x_plain, R,k);
  if(mpz_cmp(a_mpz,b_mpz) == 0)
    printf("GFPF_gmp: Pass\n");
  else{
    printf("GFPF_gmp: Error\n");
  }

  for (int i = 0; i < k; i++){
      x[i] = rand()%prime.radix;
      y[i] = rand()%prime.radix;
  }

  // radix_based_to_mpz_s64(a_mpz, x, R,k);
  // radix_based_to_mpz_s64(b_mpz, y, R,k);
  u64vector_to_radix_based_bigint_gmp (a_mpz, x, prime.radix, k);
  u64vector_to_radix_based_bigint_gmp (b_mpz, y, prime.radix, k);

  //   cpu_timer add_gmp;
  // timer_record_start(&add_gmp);
  for(int i = 0; i < 1; ++i){//1000000
  //  BigPrimeField_Multiplication(a_mpz,b_mpz,p_zz,a_mpz);
    mpz_add (a_mpz, a_mpz,b_mpz);
    mpz_mod (a_mpz, a_mpz, p_zz);
  }
  // timer_record_stop(&add_gmp);
  // timer_get_elapsed_time(&add_gmp, "add_gmp",1);

  //   cpu_timer add_gfpf;
  // timer_record_start(&add_gfpf);
  for(int i = 0; i < 1; ++i){//1000000
  //  BigPrimeField_Multiplication(a_mpz,b_mpz,p_zz,a_mpz);
    addition_big_elements (x, y, k, radix);
  }
  // timer_record_stop(&add_gfpf);
  // timer_get_elapsed_time(&add_gfpf, "add_gfpf",1);
  u64vector_to_radix_based_bigint_gmp (b_mpz, x, radix, k);
  if(mpz_cmp(a_mpz,b_mpz) == 0)
    printf("Addition: Pass\n");
  else{
    printf("Addition: Error\n");
  }

  for (int i = 0; i < k; i++){
      x[i] = rand()%prime.radix;
      y[i] = 0; 
  }
  y[1] = 1;

  // radix_based_to_mpz_s64(a_mpz, x, R,k);
  // radix_based_to_mpz_s64(b_mpz, y, R,k);
  u64vector_to_radix_based_bigint_gmp (a_mpz, x, prime.radix, k);
  u64vector_to_radix_based_bigint_gmp (b_mpz, y, prime.radix, k);
  mpz_t c_mpz;
  mpz_init (c_mpz);
  // cpu_timer gfpf_mul_pow_r;
  // timer_record_start(&gfpf_mul_pow_r);
  for(int i = 0; i < 1; ++i){
    mult_pow_R (x, i, k, prime.radix,0);
 //   mult_pow_R (x, i, k, R);
  }
  // timer_record_stop(&gfpf_mul_pow_r);
  // timer_get_elapsed_time(&gfpf_mul_pow_r, "gfpf_mul_pow_r",1);

  // cpu_timer gmp_pow;
  // timer_record_start(&gmp_pow);
  for(int i = 0; i < 1; ++i){
    mpz_pow_ui (c_mpz, b_mpz, i);
    mpz_mul (a_mpz, a_mpz,c_mpz);
    mpz_mod (a_mpz, a_mpz, p_zz);
   // BigPrimeField_Multiplication(a_mpz,c_mpz,p0,a_mpz);
  }
  // timer_record_stop(&gmp_pow);
  // timer_get_elapsed_time(&gmp_pow, "gmp_pow",1);

  // radix_based_to_mpz_s64(b_mpz, x, R,k);
   u64vector_to_radix_based_bigint_gmp (b_mpz, x, prime.radix, k);
  if(mpz_cmp(a_mpz,b_mpz) == 0)
    printf("mult_pow_R: Pass\n");
  else{
    printf("mult_pow_R: Error\n");
  }

}


