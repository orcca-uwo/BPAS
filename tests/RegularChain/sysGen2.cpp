#include <bpas.h>
#include <vector>
//#include "../MapleTestTool/MapleTestTool.hpp"
#include "sysGen.hpp"

void Sys3019Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-4*s^2*p+s^4-2*s^2+phi^2+1+4*p"),
		SMQP("s^3*m*phi^3-s^3*phi^3-5*phi^3*s^3*p-3*s*p*m*phi^3+3*s*p*phi^3+4*phi^3*s*p^2+phi^3*s^5-2*p^3*m*phi^3+2*p^3"),
		SMQP("-3*phi^3*s^2*p-2*s*p^3+phi^3*s^4-s^2*m*phi^3-s^2*phi^3+p*m*phi^3+p*phi^3")};
	std::vector<Symbol> vars = {Symbol("p"), Symbol("s"), Symbol("phi"), Symbol("m")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3019","4corps-1parameter-homog",showOutput,isLazard);
} // 4corps-1parameter-homog

// params := [x5,x4,x3,x2,x1];
void Sys3020Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x3*x12-x7*x8"),
		SMQP("x5*x11-x6*x10"),
		SMQP("x3*x9-x1*x9-x4*x8+x1*x8"),
		SMQP("x6*x12-x1*x12-x7*x11+x1*x11"),
		SMQP("x3*x10-x2*x10-x5*x8+x2*x8"),
		SMQP("x4*x11-x2*x11-x6*x9+x2*x9"),
		SMQP("x5*x12-x4*x12-x7*x10+x4*x10+x7*x9-x5*x9")};
	std::vector<Symbol> vars = {Symbol("x12"), Symbol("x11"), Symbol("x10"), Symbol("x9"), Symbol("x8"), Symbol("x7"), Symbol("x6"), Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3019","8-3-config-Li",showOutput,isLazard);
} // 8-3-config-Li

//params := [c_a,c_b,c_c];
void Sys3021Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a^2 - b^2 - c^2 + b*c*c_a"),
		SMQP("b^2 - c^2 - a^2 + c*a*c_b"),
		SMQP("c^2 - a^2 - b^2 + a*b*c_c"),
		SMQP("c_a^2 + s_a^2 -1"),
		SMQP("c_b^2 + s_b^2 -1"),
		SMQP("c_c^2 + s_c^2 -1")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("b"), Symbol("c"), Symbol("s_a"), Symbol("s_b"), Symbol("s_c"), Symbol("c_a"), Symbol("c_b"), Symbol("c_c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3021","AlKashiSinus",showOutput,isLazard);
} // AlKashiSinus

//params:=[r,t];
void Sys3023Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x*u-x*r-x*2 - u - 2*t^4 + 1"),
		SMQP("y*r - t^2*u - 1"),
		SMQP("z*t*u - 1"),
		SMQP("v - r*t")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("v"), Symbol("u"), Symbol("r"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3023","Alonso",showOutput,isLazard);
} // Alonso

//params:=[y,z,v];
void Sys3022Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x*u-x*r-x*2 - u - 2*t^4 + 1"),
		SMQP("y*r - t^2*u - 1"),
		SMQP("z*t*u - 1"),
		SMQP("v - r*t")};
	std::vector<Symbol> vars = {Symbol("r"), Symbol("u"), Symbol("t"), Symbol("x"), Symbol("y"), Symbol("z"), Symbol("v")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3024","Alonso-Li",showOutput,isLazard);
} // Alonso-Li

//params := [x,z];
void Sys3024Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*t*(t-1)^2 + (s-1)^3 +3*s - x"),
		SMQP("3*s*(s-1)^2 + t^3 + 3*t - y"),
		SMQP("-3*s*(s^2 -5*s+5*t)*t^3 -3*(s^3+6*s^2-9*s+1)*t^2 + t*(6*s^3 + 9*s^2 -18*s + 3) -3*s*(s-1) - z")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("s"), Symbol("y"), Symbol("x"), Symbol("z")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3024","Bezier",showOutput,isLazard);
} // Bezier

//params := [u,v];
void Sys3025Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x + y + z + t + u + v"),
		SMQP("x*y + y*z + z*t + t*u + x*v + u*v"),
		SMQP("x*y*z + y*z*t + z*t*u + x*y*v + x*u*v + t*u*v"),
		SMQP("x*y*z*t + y*z*t*u + x*y*z*v + x*y*u*v + x*t*u*v + z*t*u*v"),
		SMQP("x*y*z*t*u + x*y*z*t*v + x*y*z*u*v + x*y*t*u*v + x*z*t*u*v + y*z*t*u*v")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3025","Bjork60",showOutput,isLazard);
} // Bjork60

//params := [s];
void Sys3027Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("1/200*x*(s - x - z ) + y*(s - x - z ) - 35/2*x"), 
		SMQP("250*x*(s - x - z )*(z + 3/250) - 55/2*y"),
		SMQP("500*y + 25*x - 5*z")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("y"), Symbol("x"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3027","Blood-Coagulation",showOutput,isLazard);
} // Blood-Coagulation

//params := [s];
void Sys3026Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("1/200*x*s*(1 - 1/400*x) + y*s*(1 - 1/400*x) - 35/2*x"),
		SMQP("250*x*s*(1 - 1/600*y )*(z + 3/250) - 55/2*y"),
		SMQP("500*(y + 1/20*x)*(1 - 1/700*z) - 5*z")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("y"), Symbol("x"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3026","Blood-Coagulation-2",showOutput,isLazard);
} // Blood-Coagulation-2


//ineqs := [y]; pie := [1-x*y];
void Sys3157Test(bool showOutput, bool isLazard) {
//#original form:
//#(E x)(E y)[x^3-3*x*y^2+a*x+b, 3*x^2-y^2+a, y<>0, x*y-1<0];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^3-3*x*y^2+a*x+b"),
		SMQP("3*x^2-y^2+a")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("b"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3157","BM05-1",showOutput,isLazard);
} // BM05-1

//nie := [-w*x^2-v*x-u];
void Sys3158Test(bool showOutput, bool isLazard) {
//#original form:
//#(E x)[u*x^2+v*x+1=0, v*x^3+w*x+u=0, w*x^2+v*x+u<=0]
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("u*x^2+v*x+1"),
		SMQP("v*x^3+w*x+u")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("u"), Symbol("v"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3158","BM05-2",showOutput,isLazard);
} // BM05-2

void Sys3028Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*y^2*(y^2+x^2)+(b^2-3*a^2)*y^2-2*b*y^2*(y+x)+2*a^2*b*(y+x)-a^2*x^2+a^2*(a^2-b^2)"),
		SMQP("4*y^3+4*y*(y^2+x^2)-2*b*y^2-4*b*y*(y+x)+2*(b^2-3*a^2)*y+2*a^2*b"),
		SMQP("4*x*y^2-2*b*y^2-2*a^2*x+2*a^2*b")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3028","Bronstein-Wang",showOutput,isLazard);
} // Bronstein-Wang

void Sys3029Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("b1 + y + z - t - w"),
		SMQP("2*z*u + 2*y*v + 2*t*w - 2*w^2 - w - 1"),
		SMQP("3*z*u^2 + 3*y*v^2 - 3*t*w^2 + 3*w^3 + 3*w^2 - t + 4*w"),
		SMQP("6*x*z*v - 6*t*w^2 + 6*w^3 - 3*t*w + 6*w^2 - t + 4*w"),
		SMQP("4*z*u^3+ 4*y*v^3+ 4*t*w^3- 4*w^4 - 6*w^3+ 4*t*w- 10*w^2- w- 1"),
		SMQP("8*x*z*u*v +8*t*w^3 -8*w^4 +4*t*w^2 -12*w^3 +4*t*w -14*w^2 -3*w -1"),
		SMQP("12*x*z*v^2+12*t*w^3 -12*w^4 +12*t*w^2 -18*w^3 +8*t*w -14*w^2 -w -1"),
		SMQP("-24*t*w^3 + 24*w^4 - 24*t*w^2 + 36*w^3 - 8*t*w + 26*w^2 + 7*w + 1")};
	std::vector<Symbol> vars = {Symbol("b1"), Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("v"), Symbol("u"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3029","Butcher",showOutput,isLazard);
} // Butcher

//params := [v];
void Sys3030Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("30-30*x+v^4*(1-201*x)*y^4"),
		SMQP("1+x^4-(1+11*x^4)*y")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("v")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3030","CDC2-Cyclin",showOutput,isLazard);
} // CDC2-Cyclin

//params:=[c1, c2, c3, c5];
void Sys3031Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^3*y^2+c1*x^3*y+y^2+c2*x+c3"),
		SMQP("c4*x^4*y^2-x^2*y+y+c5"),
		SMQP("c4-1")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("c4"), Symbol("c1"), Symbol("c2"), Symbol("c3"), Symbol("c5")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3031","Cheaters-homotopy-easy",showOutput,isLazard);
} // Cheaters-homotopy-easy

//params:=[c1, c4, c2, c3, c5];
void Sys3032Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^3*y^2+c1*x^3*y+y^2+c2*x+c3"),
		SMQP("c4*x^4*y^2-x^2*y+y+c5")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("c1"), Symbol("c4"), Symbol("c2"), Symbol("c3"), Symbol("c5")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3032","Cheaters-homotopy-hard",showOutput,isLazard);
} // Cheaters-homotopy-hard

//params:=[t];
void Sys3033Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2 - 7 * x1 + x1 ^ 2 * x2 + t * x3 - t * x1"),
		SMQP("6 * x1 - x1 ^ 2 * x2 + 10 * t * x4 - t * x2"),
		SMQP("2 - 7 * x3 + x3 ^ 2 * x4 + t * x1 - t * x3"),
		SMQP("6 * x3 - x3 ^2 * x4  + 1 - t * x2 - t * x4")};
	std::vector<Symbol> vars = {Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3033","Chemical",showOutput,isLazard);
} // Chemical

void Sys3034Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("6*a33*a10*a20+10*a22*a10*a31+8*a32*a10*a21-162*a10^2*a21+16*a21*a30+14*a31*a20+48*a10*a30"),
		SMQP("15*a33*a10*a21-162*a10^2*a22-312*a10*a20+24*a10*a30+27*a31*a21+24*a32*a20+18*a22*a10*a32+30*a22*a30+84*a31*a10"),
		SMQP("-240*a10+420*a33-64*a22+112*a32"),
		SMQP("180*a33*a10-284*a22*a10-162*a10^2+60*a22*a32+50*a32*a10+70*a30+55*a33*a21+260*a31-112*a20"),
		SMQP("66*a33*a10+336*a32+90*a31+78*a22*a33-1056*a10-90*a21"),
		SMQP("136*a33-136"),
		SMQP("4*a22*a10*a30+2*a32*a10*a20+6*a20*a30-162*a10^2*a20+3*a31*a21*a10"),
		SMQP("28*a22*a10*a33+192*a30+128*a32*a10+36*a31*a20+36*a33*a20-300*a10*a21+40*a32*a21-648*a10^2+44*a22*a31")};
	std::vector<Symbol> vars = {Symbol("a10"), Symbol("a20"), Symbol("a21"), Symbol("a22"), Symbol("a30"), Symbol("a31"), Symbol("a32"), Symbol("a33")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3034","ChildDraw-1",showOutput,isLazard);
} // ChildDraw-1

void Sys3035Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("16*a20*a32+18*a21*a31+20*a22*a30"),
		SMQP("-80*a23+180*a34+855*a35"),
		SMQP("7*a20*a31+8*a21*a30"),
		SMQP("210*a35-210"),
		SMQP("40*a20*a34+44*a21*a33+48*a22*a32+52*a23*a31+280*a30"),
		SMQP("27*a20*a33+30*a21*a32+33*a22*a31+36*a23*a30"),
		SMQP("55*a20*a35+60*a21*a34+65*a22*a33+70*a23*a32+80*a30+375*a31"),
		SMQP("78*a21*a35+84*a22*a34+90*a23*a33-170*a20+102*a31+480*a32"),
		SMQP("136*a23*a35-114*a22+152*a33+720*a34"),
		SMQP("105*a22*a35+112*a23*a34-144*a21+126*a32+595*a33")};
	std::vector<Symbol> vars = {Symbol("a20"), Symbol("a32"), Symbol("a21"), Symbol("a31"), Symbol("a22"), Symbol("a30"), Symbol("a23"), Symbol("a35"), Symbol("a34"), Symbol("a33")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3035","ChildDraw-2",showOutput,isLazard);
} // ChildDraw-2

//params :=  [s];
void Sys3036Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s-x1*(1+x2+x3)"),
		SMQP("s-x2*(1+x1+x3)"),
		SMQP("s-x3*(1+x1+x2)")};
	std::vector<Symbol> vars = {Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3036","Cinquin-Demongeot-3-1",showOutput,isLazard);
} // Cinquin-Demongeot-3-1

//params :=  [s];
void Sys3037Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s-x1*(1+x2^2+x3^2)"),
		SMQP("s-x2*(1+x1^2+x3^2)"),
		SMQP("s-x3*(1+x1^2+x2^2)")};
	std::vector<Symbol> vars = {Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3037","Cinquin-Demongeot-3-2",showOutput,isLazard);
} // Cinquin-Demongeot-3-2

//params :=  [s];
void Sys3038Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s-x1*(1+x2^3+x3^3)"),
		SMQP("s-x2*(1+x1^3+x3^3)"),
		SMQP("s-x3*(1+x1^3+x2^3)")};
	std::vector<Symbol> vars = {Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3038","Cinquin-Demongeot-3-3",showOutput,isLazard);
} // Cinquin-Demongeot-3-3

//params :=  [s];
void Sys3039Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s-x1*(1+x2^4+x3^4)"),
		SMQP("s-x2*(1+x1^4+x3^4)"),
		SMQP("s-x3*(1+x1^4+x2^4)")};
	std::vector<Symbol> vars = {Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3039","Cinquin-Demongeot-3-4",showOutput,isLazard);
} // Cinquin-Demongeot-3-4

//params :=  [s];
void Sys3040Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s-x1*(1+x2+x3+x4+x5)"),
		SMQP("s-x2*(1+x1+x3+x4+x5)"),
		SMQP("s-x3*(1+x1+x2+x4+x5)"),
		SMQP("s-x4*(1+x1+x2+x3+x5)"),
		SMQP("s-x5*(1+x1+x2+x3+x4)")};
	std::vector<Symbol> vars = {Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3040","Cinquin-Demongeot-5-1",showOutput,isLazard);
} // Cinquin-Demongeot-5-1

//params := [r, a]; ineqs := [u+1,v-u,1-v,r-b,r-1];
void Sys3041Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("u^3+r*u^2+u^2-a*u+r*u+u-b-a+r+1"),
		SMQP("v^3+r*v^2-v^2-a*v-r*v+v-b+a+r-1"),
		SMQP("4*u^3+3*r*u^2-2*a*u-b"),
		SMQP("4*v^3+3*r*v^2-2*a*v-b")};
	std::vector<Symbol> vars = {Symbol("b"), Symbol("v"), Symbol("u"), Symbol("r"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3041","Collins-JSC02",showOutput,isLazard);
} // Collins-JSC02

//params := [u2, u1]; ineqs := [x2-2, u2];
void Sys3042Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1-u1)*(x1-1)+(x2-u2)*x2"),
		SMQP("(x1-1)^2+x2^2-1"),
		SMQP("u1*x4-u2*x3"),
		SMQP("x4*(x1-2)-x2*(x3-2)"),
		SMQP("(u1-x1)^2+(u2-x2)^2-(u1-x3)^2-(u2-x4)^2")};
	std::vector<Symbol> vars = {Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("u2"), Symbol("u1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3042","Cox-ISSAC07",showOutput,isLazard);
} // Cox-ISSAC07

void Sys3043Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a + b + c + d"),
		SMQP("a * b + b * c + c * d + d * a"),
		SMQP("a * b * c + b * c * d + c * d * a + d * a * b"),
		SMQP("a * b * c * d - 1")};
	std::vector<Symbol> vars = {Symbol("a"),Symbol("b"), Symbol("c"), Symbol("d")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3043","Cyclic-4",showOutput,isLazard);
} // Cyclic-4

//ineqs := [1+x1^3];
void Sys3179Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1-u1*(1+x1^3)"),
		SMQP("3*x1^2-u2*(1+x1^3)")};
	std::vector<Symbol> vars = {Symbol("x1"),Symbol("u1"), Symbol("u2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3179","DescartesFolium",showOutput,isLazard);
} // DescartesFolium

void Sys3156Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a*g"),
		SMQP("g*j + a*m + n*p + q"),
		SMQP("b*l"),
		SMQP("n*q"),
		SMQP("b*g + b*k + a*l + l*o + l*p + b + c"),
		SMQP("a*g + a*k + j*l + b*m + b*n + g*o + k*o + g*p + k*p + l*q + a + d + f + h + o + p"),
		SMQP("g*j + j*k + a*m + a*n + m*o + n*o + m*p + n*p + g*q + k*q + e + j + q + s - 1"),
		SMQP("j*m + j*n + m*q + n*q"),
		SMQP("j*n + m*q + 2*n*q"),
		SMQP("g*j + a*m + 2*a*n + n*o + n*p + 2*g*q + k*q + q + s"),
		SMQP("2*a*g + a*k + b*n + g*o + g*p + l*q + a + d"),
		SMQP("b*g + a*l"),
		SMQP("a*n + g*q"),
		SMQP("2*j*m + j*n + m*q"),
		SMQP("g*j + j*k + a*m + m*o + 2*m*p + n*p + e + 2*j + q"),
		SMQP("j*l + b*m + g*p + k*p + a + f + o + 2*p"),
		SMQP("l*p + b"),
		SMQP("j*n + m*q"),
		SMQP("g*p + a")};
	std::vector<Symbol> vars = {Symbol("s"), Symbol("q"), Symbol("p"), Symbol("o"), Symbol("n"), Symbol("m"), Symbol("l"), Symbol("k"), Symbol("j"), Symbol("h"), Symbol("g"), Symbol("f"), Symbol("e"), Symbol("d"), Symbol("c"), Symbol("b"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3156","DGP29",showOutput,isLazard);
} // DGP29

void Sys3155Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("s^15"),
		SMQP("t^15"),
		SMQP("u^15"),
		SMQP("u^5 - s^3*t*x + s^2*t^2*x + s^2*t^2*y - s*t^3*y")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("u"), Symbol("t"), Symbol("s")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3155","DGP6",showOutput,isLazard);
} // DGP6

//params :=[t];
void Sys3045Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^31 - x^6 - x - y"),
		SMQP("x^8  - z"),
		SMQP("x^10 - t")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3045","DonatiTraverso",showOutput,isLazard);
} // DonatiTraverso

//params:=[x];
void Sys3044Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("z^5 -t^4"),
		SMQP("t*z*y^2 +2*z^3*y -t^8 +2*t^5 +t^3 -t^2"),
		SMQP("t^4*x -t*x -t*y -z^2")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("z"), Symbol("y"), Symbol("x")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3044","DonatiTraverso-rev",showOutput,isLazard);
} // DonatiTraverso-rev

//nie := [x1, 2-x1, x2-2, 4-x2, y2+1, 1-y2, x+1, 9-x, y+6, 6-y];
void Sys3181Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x-x1*x2+y2"),
		SMQP("y-x1*y2-x2")};
	std::vector<Symbol> vars = {Symbol("y2"), Symbol("x2"), Symbol("x1"), Symbol("y"), Symbol("x")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3181","EdgeSquare",showOutput,isLazard);
} // EdgeSquare

//nie := [x^2 + y^2-1]; pie := [a, b];
void Sys3159Test(bool showOutput, bool isLazard) {
//#orginal form:
//# (E x)(E y)[b*(x - c)^2 + a*(y - d)^2 -a*b, x^2 + y^2-1>=0, a>0, b>0]
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("b*(x - c)^2 + a*(y - d)^2 -a*b")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("d"), Symbol("c"), Symbol("b"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3159","Ellipse",showOutput,isLazard);
} // Ellipse

//params := [u1, u2, u3];
void Sys3180Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1+3*x1*x2^2-x1^3-u1"),
		SMQP("3*x2+3*x1^2*x2-x2^3-u2"),
		SMQP("3*x1^2-3*x2^2-u3")};
	std::vector<Symbol> vars = {Symbol("x2"), Symbol("x1"), Symbol("u1"), Symbol("u2"), Symbol("u3")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3180","EnneperSurface",showOutput,isLazard);
} // EnneperSurface

void Sys3046Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*u6 + 2*u5 + 2*u4 + 2*u3 + 2*u2 + 1"),
		SMQP("8*U5*u6 + 8*U5*u5 + 8*U4*u6 +8*U4*u5 + 8*U4*u4 + 8*U3*u6 + 8*U3*u5 + 8*U3*u4 + 8*U3*u3 + 8*U2*u6 +8*U2*u5 + 8*U2*u4 + 8*U2*u3 + 8*U2*u2 -13"),
		SMQP("2*U6 + 2*U5 + 2*U4 + 2*U3 +2*U2 + 1"),
		SMQP("8*U6*u5 + 8*U6*u4 + 8*U6*u3 + 8*U6*u2 + 8*U5*u5 + 8*U5*u4 +8*U5*u3 + 8*U5*u2 + 8*U4*u4 + 8*U4*u3 + 8*U4*u2 + 8*U3*u3 + 8*U3*u2 +8*U2*u2 -13"),
		SMQP("U2*u2 -1"),
		SMQP("U3*u3 -1"),
		SMQP("U4*u4 -1"),
		SMQP("U5*u5 -1"),
		SMQP("U6*u6 -1")};
	std::vector<Symbol> vars = {Symbol("U6"),Symbol("U5"), Symbol("U4"), Symbol("U3"), Symbol("U2"), Symbol("u6"), Symbol("u5"), Symbol("u4"), Symbol("u3"), Symbol("u2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3046","F-633",showOutput,isLazard);
} // F-633

//params := [u2];
void Sys3047Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*u7 +2*u6 +2*u5 +2*u4 +2*u3 +2*u2 +1"),
		SMQP("u2*U2 -1"),
		SMQP("2*U7 +2*U6 +2*U5 +2*U4 +2*U3 +2*U2 +1"),
		SMQP("u3*U3 -1"),
		SMQP("u4*U4 -1"),
		SMQP("((16*u3 +16*u2 +8)*u5)*U4 +((16*u2 +8)*u5 +(16*u2 +8)*u4)*U3 +(8*u5 +8*u4 +8*u3)*U2 +18*u5 +18*u4 +18*u3 +18*u2 +11"),
		SMQP("u5*U5 -1"),
		SMQP("8*u6*U5 +(8*u6 +8*u5)*U4 +(8*u6 +8*u5 +8*u4)*U3 +(8*u6 +8*u5 +8*u4 +8*u3)*U2 +4*u6 +4*u5 +4*u4 +4*u3 +4*u2 +17"),
		SMQP("(16*u4*U3 +(16*u4 +16*u3)*U2 +8*u4 +8*u3 +8*u2 +18)*U5 +(16*u3*U2 +8*u3 +8*u2 +18)*U4 +(8*u2 +18)*U3 +18*U2 +11"),
		SMQP("(8*u5 +8*u4 +8*u3 +8*u2 +4)*U6 +(8*u4 +8*u3 +8*u2 +4)*U5 +(8*u3 +8*u2 +4)*U4 +(8*u2 +4)*U3 +4*U2 +17"),
		SMQP("u6*U6 -1"),
		SMQP("u7*U7 -1")};
	std::vector<Symbol> vars = {Symbol("U7"), Symbol("u7"), Symbol("U6"), Symbol("U5"), Symbol("U4"), Symbol("U3"), Symbol("U2"), Symbol("u6"), Symbol("u5"), Symbol("u4"), Symbol("u3"), Symbol("u2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3047","F-744",showOutput,isLazard);
} // F-744

void Sys3048Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a11 * x1 + a12 * x2 - b1"),
		SMQP("a21 * x1 + a22 * x2 - b2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("a11"), Symbol("a21"), Symbol("a12"), Symbol("a22"), Symbol("b1"), Symbol("b2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3048","GenLinSyst-2-2",showOutput,isLazard);
} // GenLinSyst-2-2

void Sys3049Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a*x+b*y-c"),
		SMQP("d*x+e*y-f"),
		SMQP("g*x+h*x-i")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d"), Symbol("e"), Symbol("f"), Symbol("g"), Symbol("h"), Symbol("i")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3049","GenLinSyst-3-2",showOutput,isLazard);
} // GenLinSyst-3-2

//params := [a11, a21, a31, a12, a22, a32, a13, a23, a33, b1, b2, b3];
void Sys3050Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a11 * x + a12 * y + a13 * z - b1"),
		SMQP("a21 * x + a22 * y + a23 * z - b2"),
		SMQP("a31 * x + a32 * y + a33 * z - b3")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("y"), Symbol("x"), Symbol("a11"), Symbol("a21"), Symbol("a31"), Symbol("a12"), Symbol("a22"), Symbol("a32"), Symbol("a13"), Symbol("a23"), Symbol("a33"), Symbol("b1"), Symbol("b2"), Symbol("b3")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3050","GenLinSyst-3-3",showOutput,isLazard);
} // GenLinSyst-3-3

//params:=[u1, u2, u3, u4];
void Sys3051Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x1*x8-x2*x8-x1*u2+x2*u1-x9*u1+x9*u2"),
		SMQP("-x5*x6+x4*x7+x5*x8-x7*x8-x4*x9+x6*x9"),
		SMQP("-x4*x6-x5*x7+x4*x8+x5*x9"),
		SMQP("-x3*x6+x3*u4+x7*u3-x7*u4"),
		SMQP("-x2*x4+x2*u4+x5*u2-x5*u4"),
		SMQP("x1*x4-x3*x4-x1*u3+x3*u1-x5*u1+x5*u3"),
		SMQP("x3^2+u3^2-u4^2"),
		SMQP("x2^2+u2^2-u4^2"),
		SMQP("x1^2+u1^2-u4^2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("x7"), Symbol("x8"), Symbol("x9"), Symbol("u1"), Symbol("u2"), Symbol("u3"), Symbol("u4")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3051","Geometry-Butterfly_3",showOutput,isLazard);
} // Geometry-Butterfly_3

void Sys3052Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*y*w-z*w+t*w"),
		SMQP("-2*u*w^2+10*v*w^2-20*w^3+7*t*u-35*t*v+10*t*w"),
		SMQP("2*y*w^2-2*z*w^2+6*t*w^2-7*y*t+7*z*t-21*t^2"),
		SMQP("-2*v^3+4*u*v*w+5*v^2*w-6*u*w^2-7*v*w^2+15*w^3+42*y*v-14*z*v-63*y*w+21*z*w-42*t*w+147*x"),
		SMQP("-9*u*w^3+45*v*w^3-135*w^4+14*z*v^2-14*t*v^2-28*z*u*w+70*t*u*w-14*z*v*w-28*7*t*v*w+4*7*z*w^2+86*7*t*w^2-42*7*y*z+14*7*z^2+42*7*y*t-14*7*z*t-21*7*x*u+105*7*x*v-315*7*x*w"),
		SMQP("6*y*w^3-9*z*w^3+36*t*w^3-2*7*x*v^2-4*7*y*t*w+6*7*z*t*w-24*7*t^2*w+4*7*x*u*w+2*7*x*v*w-4*7*x*w^2+ 56*7*x*y-35*7*x*z+84*7*x*t"),
		SMQP("2*u*v*w-6*v^2*w-u*w^2+13*v*w^2-5*w^3+14*y*w-28*t*w"),
		SMQP("u^2*w-3*u*v*w+5*u*w^2+14*y*w-28*t*w"),
		SMQP("-2*z*u*w-2*t*u*w+4*y*v*w+6*z*v*w-2*t*v*w-16*y*w^2-10*z*w^2+22*t*w^2+42*x*w"),
		SMQP("28*y*u*w+8*z*u*w-20*t*u*w-88*y*v*w-8*3*z*v*w+68*t*v*w+52*3*y*w^2+40*z*w^2-44*3*t*w^2-84*3*x*w"),
		SMQP("-4*y*z*w+10*y*t*w+8*z*t*w-20*t^2*w+12*x*u*w-30*x*v*w+15*x*w^2"),
		SMQP("-2*y^2*w+y*z*w+2*y*t*w-2*z*t*w+2*2*t^2*w-3*2*x*u*w+6*2*x*v*w-3*2*x*w^2"),
		SMQP("8*x*y*w-4*x*z*w+8*x*t*w")};
	std::vector<Symbol> vars = {'x','y','z','t','u','v','w'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3052","Gerdt",showOutput,isLazard);
} // Gerdt

void Sys3053Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a*g"),
		SMQP("d*f+b*g+a*h+w"),
		SMQP("c*i"),
		SMQP("w*f"),
		SMQP("c*g+a*i+d*i+e*i+c*j+t+c"),
		SMQP("c*f+a*g+d*g+e*g+c*h+b*i+a*j+d*j+e*j+y+u+v+i*w+a+d+e"),
		SMQP("a*f+d*f+e*f+w*g+b*g+a*h+d*h+e*h+w*j+b*j+x+z+w+b-1"),
		SMQP("w*f+b*f+w*h+b*h"),
		SMQP("2*w*f+b*f+w*h"),
		SMQP("2*a*f+d*f+e*f+2*w*g+b*g+a*h+w*j+x+w"),
		SMQP("c*f+2*a*g+d*g+e*g+w*i+a*j+y+a"),
		SMQP("c*g+a*i"),
		SMQP("a*f+w*g"),
		SMQP("b*f+w*h+2*b*h"),
		SMQP("d*f+b*g+a*h+2*d*h+e*h+b*j+z+w+2*b"),
		SMQP("d*g+c*h+b*i+d*j+u+a+2*d+e"),
		SMQP("d*i+c"),
		SMQP("b*f+w*h"),
		SMQP("d*g+a")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d"), Symbol("e"), Symbol("f"), Symbol("g"), Symbol("h"), Symbol("i"), Symbol("j")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3053","Gonnet",showOutput,isLazard);
} // Gonnet

//params := [a, b];
void Sys3190Test(bool showOutput, bool isLazard) {
//##- Family of *hard* problems (Haas parametric systems)
//#n := 5;
//#unknowns6 := [x,y];
//#sys6 := [ x^(2*n)+a*y^n-y, y^(2*n)+b*x^n-x ];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^10+a*y^5-y"),
		SMQP("y^10+b*x^5-x")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3190","Haas5",showOutput,isLazard);
} // Haas5

//params:=[C3,C2];
void Sys3054Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("B1+B2+B3+B4-1"),
		SMQP("2*B2*C2 + 2*B3*C3 + 2*B4*C4 - 1"),
		SMQP("3*B2*C2^2 +3*B3*C3^2 +3*B4*C4^2 -1"),
		SMQP("6*B3*A32*C2 +6*B4*A42*C2 +6*B4*A43*C3 -1"),
		SMQP("4*B2*C2^3 +4*B3*C3^3 +4*B4*C4^3 -1"),
		SMQP("8*B3*C3*A32*C2 +8*B4*C4*A42*C2 +8*B4*C4*A43*C3 -1"),
		SMQP("12*B3*A32*C2^2 +12*B4*A42*C2^2 +12*B4*A43*C3^2 -1"),
		SMQP("24*B4*A43*A32*C2 -1"),
		SMQP("-A21+C2"),
		SMQP("-A31-A32+C3"),
		SMQP("-A41-A42-A43+C4")};
	std::vector<Symbol> vars = {Symbol("A43"), Symbol("A42"), Symbol("A41"), Symbol("A32"), Symbol("A31"), Symbol("A21"), Symbol("B1"), Symbol("B2"), Symbol("B3"), Symbol("B4"), Symbol("C4"), Symbol("C3"), Symbol("C2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3054","Hairer-2-BGK",showOutput,isLazard);
} // Hairer-2-BGK

void Sys3055Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^9 - 3*x^6*z^3 + 3*x^3*z^6 - z^9 - 2*x^6*t^3 + 4*x^3*z^3*t^3 -2*z^6*t^3 + x^3*t^6 - z^3*t^6 + 3*x^7*u^2 - 6*x^4*z^3*u^2 +3*x*z^6*u^2 - 4*x^4*t^3*u^2 + 4*x*z^3*t^3*u^2 + x*t^6*u^2 + 3*x^5*u^4- 3*x^2*z^3*u^4 - 2*x^2*t^3*u^4 + x^3*u^6 + x^4 - 2*x^2*y^2 + y^4"),
		SMQP("84*x^6 - 60*x^3*z^3 + 3*z^6 - 40*x^3*t^3 + 4*z^3*t^3 + t^6 +105*x^4*u^2 - 24*x*z^3*u^2 - 16*x*t^3*u^2 + 30*x^2*u^4 + u^6 + 4*x"),
		SMQP("36*x^7 - 45*x^4*z^3 + 9*x*z^6 - 30*x^4*t^3 + 12*x*z^3*t^3 + 3*x*t^6 +63*x^5*u^2 - 36*x^2*z^3*u^2 - 24*x^2*t^3*u^2 + 30*x^3*u^4 - 3*z^3*u^4- 2*t^3*u^4 + 3*x*u^6 + 6*x^2 - 2*y^2"),
		SMQP("9*x^8 - 18*x^5*z^3 + 9*x^2*z^6 - 12*x^5*t^3 + 12*x^2*z^3*t^3 +3*x^2*t^6 + 21*x^6*u^2 - 24*x^3*z^3*u^2 + 3*z^6*u^2 - 16*x^3*t^3*u^2 +4*z^3*t^3*u^2 + t^6*u^2 + 15*x^4*u^4 - 6*x*z^3*u^4 - 4*x*t^3*u^4 +3*x^2*u^6 + 4*x^3 - 4*x*y^2")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3055","Hawes-1",showOutput,isLazard);
} // Hawes-1

void Sys3056Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a21*c1 + a11*c2"),
		SMQP("a22*c1 + a12*c2"),
		SMQP("a11*c1*(3*a12+2*a*c1^2)"),
		SMQP("a12*c1*(a12+4*a*c1^2)"),
		SMQP("a11*c1-3*a10*a11*c1+2*a*a11*c1^3+a21*c2"),
		SMQP("-3*a11^2*c1+2*12*c1-6*a10*a12*c1+16*a*a12*c1^3+2*a22*c2")};
	std::vector<Symbol> vars = {Symbol("a11"), Symbol("a21"), Symbol("a10"), Symbol("a22"), Symbol("a12"), Symbol("c2"), Symbol("c1"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3056","Hereman-2",showOutput,isLazard);
} // Hereman-2

//params := [b, a ];
void Sys3057Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a * a12^2 + 6*b*a12*c1^2 + 12*g*a12*c1^2+360*c1^4"),
		SMQP("a11*(a*a12^2+2*b*a12*c1^2+6*g*a12*c1^2+24*c1^4)"),
		SMQP("a11*(a*a10^2*c1-2*g*a10*c1^3+2*b*a12*c1^3+16*c1^5+c2)"),
		SMQP("a11*(a*a11^2+6*a*a10*a12+6*g*a10*c1^2-12*b*a12*c1^2-18*g*a12*c1^2-120*c1^4)"),
		SMQP("2*a*a11^2*a12 + 2*a*a10*a12^2 + b*a11^2*c1^2+3*g*a11^2*c1^2+12*g*a10*a12*c1^2-8*b*a12^2*c1^2-8*g*a12^2*c1^2-480*a12*c1^4"),
		SMQP("a*a10*a11^2*c1 + a*a10^2*a12*c1 - b*a11^2*c1^3-g*a11^2*c1^3 - 8*g*a10*a12*c1^3 + 2*b*a12^2*c1^3 + 136*a12*c1^5 + a12*c2")};
	std::vector<Symbol> vars = {Symbol("a11"), Symbol("a12"), Symbol("a10"), Symbol("c2"), Symbol("c1"), Symbol("g"), Symbol("b"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3057","Hereman-8",showOutput,isLazard);
} // Hereman-8

//params:=[d15,d14,d13,d12,d11,d10,d9,d8,d7,d6,d5,d4,d3,d2,d1,d0];
void Sys3058Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x2^2+x3^2-d2^2"),
		SMQP("x5^2+x6^2-d3^2"),
		SMQP("x8^2+x9^2-d4^2"),
		SMQP("x0^2-2*x0*x3+x1^2+x3^2-d6^2"),
		SMQP("x0^2-2*x0*x6+x4^2+x6^2-d7^2"),
		SMQP("x0^2-2*x0*x9+x7^2+x9^2-d8^2"),
		SMQP("x1^2-2*x1*x4-2*x2*x5-2*x3*x6+x4^2+d2^2+d3^2-d10^2"),
		SMQP("x1^2-2*x1*x7-2*x2*x8-2*x3*x9+x7^2+d2^2+d4^2-d11^2"),
		SMQP("x4^2-2*x4*x7-2*x5*x8-2*x6*x9+x7^2+d3^2+d4^2-d13^2"),
		SMQP("x1^2-2*x1*x10+x2^2-2*x2*x11+x10^2+x11^2-d12^2"),
		SMQP("x4^2-2*x4*x10+x5^2-2*x5*x11+x10^2+x11^2-d14^2"),
		SMQP("x7^2-2*x7*x10+x8^2-2*x8*x11+x10^2+x11^2-d15^2")};
	std::vector<Symbol> vars = {Symbol("x11"), Symbol("x10"), Symbol("x9"), Symbol("x8"), Symbol("x7"), Symbol("x6"), Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("x0"), Symbol("d15"), Symbol("d14"), Symbol("d13"), Symbol("d12"), Symbol("d11"), Symbol("d10"), Symbol("d9"), Symbol("d8"), Symbol("d7"), Symbol("d6"), Symbol("d5"), Symbol("d4"), Symbol("d3"), Symbol("d2"), Symbol("d1"), Symbol("d0")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3058","Hoffmann",showOutput,isLazard);
} // Hoffmann

//pie := [r1, -e1];
void Sys3160Test(bool showOutput, bool isLazard) {
//#original form
//#[IBVP]
//#(a1, a2, lambda1, lambda2, xi1, xi2, eta1, eta2)
//#2
//#(E lambda1)(E lambda2)(E xi1)(E xi2)(E eta1)(E eta2)
//#[
//#  lambda1^2 - lambda2^2 - eta1^2 + eta2^2 + xi1^2 - 2 xi1 xi2 + xi2^2 = 0
//#/\
//#  2 lambda1 lambda2 - 2 eta1 eta2 = 0
//#/\
//#  - lambda1 + eta1 + xi1 a2 - xi2 a2 = 0 
//#/\  
//#  - lambda2 + eta2 - xi1 a1 + xi2 a1 = 0
//#/\
//#  lambda1^2 + lambda2^2 + xi1^2 + xi2^2 + eta1^2 + eta2^2 - 1 = 0 
//#/\
//#  lambda1 > 0
//#/\ 
//#  - eta1 > 0
//#].
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("r1^2 - r2^2 - e1^2 + e2^2 + x1^2 - 2*x1*x2 + x2^2"),
		SMQP("2*r1*r2 - 2*e1*e2"),
		SMQP("-r1 + e1 + x1*a2 - x2*a2"),
		SMQP("-r2 + e2 - x1*a1 + x2*a1"),
		SMQP("r1^2 + r2^2 + x1^2 + x2^2 + e1^2 + e2^2 - 1")};
	std::vector<Symbol> vars = {Symbol("r1"), Symbol("r2"), Symbol("x1"), Symbol("x2"), Symbol("e1"), Symbol("e2"), Symbol("a1"), Symbol("a2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3160","IBVP",showOutput,isLazard);
} // IBVP

//nie := [u+1/2, 1/2-u];
void Sys3174Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-x1 + x2*u, -x2 + (1 + x1^2)*u + u^3")};
	std::vector<Symbol> vars = {Symbol("u"), Symbol("x1"), Symbol("x2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3174","Jirstrand22",showOutput,isLazard);
} // Jirstrand22

//nie := [1-u1^2, 1-u2^2, 1-u3^2];
void Sys3175Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-38*x2-170*x1*x2+148*x1^2*x2+4*x2^3+u1*(-52-2*x1+114*x1^2-79*x1^3+7*x2^2+14*x1*x2^2)+u3*(14-10*x1+37*x1^2-48*x1^3+8*x1^4-13*x2^2-13*x1*x2^2+20*x1^2*x2^2+11*x2^4)"),
		SMQP("-12-125*u2+u2^2+6*u2^3+95*x1-21*u2*x1+17*u2^2*x1-202*x1^2+81*u2*x1^2+139*x1^3"),
		SMQP("139*x2-112*x1*x2-388*x1^2*x2+215*x1^3*x2-38*x2^3+185*x1*x2^3+u1*(-11+35*x1-22*x1^2+5*x2^2+10*x1^3-17*x1*x2^2)+u3*(-44+3*x1-63*x1^2+34*x2^2+142*x1^3+63*x1*x2^2-54*x1^4-69*x1^2*x2^2-26*x2^4)")};
	std::vector<Symbol> vars = {Symbol("u3"), Symbol("u2"), Symbol("u1"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3175","Jirstrand23",showOutput,isLazard);
} // Jirstrand23

//pie := [2*x1+3*x1^2+9*x1^2*x2^2]; nie := [1-u^2];
void Sys3176Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-x1^3+x2"),
		SMQP("-x1^2-x2-x2^3+u")};
	std::vector<Symbol> vars = {Symbol("u"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3176","Jirstrand24",showOutput,isLazard);
} // Jirstrand24

//pie := [r]; nie := [t, 1-t, u+1, 1-u];
void Sys3177Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-t+2-r"),
		SMQP("-(3*t^2-2*t^3)-t^2+4*u-r*(6*t-6*t^2)")};
	std::vector<Symbol> vars = {Symbol("r"), Symbol("u"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3177","Jirstrand41",showOutput,isLazard);
} // Jirstrand41

//pie := [r]; nie := [t, 1-t, u+1, 1-u];
void Sys3178Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("t*x11+u-r*x11"),
		SMQP("(t*(x21-1)+1)^2-r*(x21-1)")};
	std::vector<Symbol> vars = {Symbol("r"), Symbol("u"), Symbol("t"), Symbol("x21"), Symbol("x11")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3178","Jirstrand42",showOutput,isLazard);
} // Jirstrand42

//params := [l, u, v, w];
void Sys3191Test(bool showOutput, bool isLazard) {
//##- A simple system with a big solution
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("l-a^2-b^2+2*a*b*w"),
		SMQP("l-a^2-c^2+2*a*c*v"),
		SMQP("l-b^2-c^2+2*b*c*u")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("b"), Symbol("c"), Symbol("l"), Symbol("u"), Symbol("v"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3191","John5",showOutput,isLazard);
} // John5

//params := [a, b, c, d];
void Sys3192Test(bool showOutput, bool isLazard) {
//##- A system from a geometric problem (Kahan ellipsoid)
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("b^2*(x-c)^2+a^2*(y-d)^2-a^2*b^2"),
		SMQP("x^2+y^2-1")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3192","Kahan-ellipsoid",showOutput,isLazard);
} // Kahan-ellipsoid

//params:=[b10, b9, b8, b7, b6, b5, b4, a10, a9, a8, a7, a6, a5, a4, a3];
void Sys3059Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x1*a4+x2*a5+x3*a6+x4*a8+x5*a9+x6*a10-1"),
		SMQP("y2*b4+y3*b5+y4*b7+y5*b8+y6*b9-1"),
		SMQP("b7"),
		SMQP("a3-b5"),
		SMQP("a7-a3"),
		SMQP("-7*b10*a5-14*a5*a9-14*a9^2+35*a5*a8+42*a6*a3-7*a6*a7-35*b9*a6-7*b10*a4+14*a8*a4-7*a8^2-42*a5*a4-14*a10*a3+21*a9*a8-14*a9*a4-42*b4*a6+21*b6*a4-28*b6*a8+7*b10*a8+21*b6*a9-7*b10*a9+14*b4*a10+21*b6*a5+28*b8*a6+14*b5*a10"),
		SMQP("-35*b5*a5-7*a9*a7-7*a4*a7+56*a5*a7-7*b9*a4-28*a8*a3+7*a8*a7+49*a4*a3-49*a5*a3-21*a9*a3-28*b5*a4+7*b7*a6+21*b5*a8+28*b5*a9-14*b8*a9-7*b4*a9-35*b8*a5+21*b9*a5+21*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("35*b5*a5-21*a9*a7-21*a4*a7+56*a5*a7-21*b9*a4+56*a8*a3+21*a8*a7+35*a4*a3-35*a5*a3-91*a9*a3+28*b5*a4+21*b7*a6-77*b5*a8+112*b5*a9-42*b8*a9-21*b4*a9+7*b8*a5-77*b9*a5+91*b4*a5+42*b8*a8-21*b9*a8+21*b9*a9"),
		SMQP("-7*b5*a5-7*a9*a7-7*a4*a7+56*a5*a7-7*b9*a4+7*a8*a7+49*a4*a3-49*a5*a3-49*a9*a3-28*b5*a4+7*b7*a6-7*b5*a8+56*b5*a9-14*b8*a9-7*b4*a9-35*b8*a5-7*b9*a5+49*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("42*b4*a8+28*b4*a4-35*b5*a5-7*a9*a7-7*a4*a7+28*a5*a7-21*b9*a4-14*a8*a3+7*a8*a7+7*a4*a3-7*a5*a3-7*a9*a3+14*b5*a4+49*b7*a6+14*b8*a4+21*b5*a8-14*b8*a9-91*b4*a9-7*b8*a5+21*b9*a5+35*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("14*b10*a5+28*a5*a9-14*a9^2-28*a5*a8-42*a6*a3+14*a6*a7+7*b9*a6-7*b10*a4-7*a8*a4-7*a8^2+21*a4^2+21*a5*a4+7*a10*a3+21*a9*a8-14*a9*a4+21*b4*a6+21*b5*a6-7*b6*a8+7*b10*a8+21*b6*a9-7*b10*a9-7*b4*a10-42*b6*a5-14*b8*a6-7*b5*a10"),
		SMQP("-7*b5*a5-7*a9*a7-7*a4*a7-7*b9*a4+7*a8*a7-7*a4*a3+7*a5*a3+7*a9*a3+28*b5*a4+7*b7*a6-7*b5*a8-14*b8*a9-7*b4*a9+21*b8*a5-7*b9*a5-7*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("7*b10*a5+14*a5*a9+14*a9^2-35*a5*a8-168*a6*a3+70*a6*a7+35*b9*a6+7*b10*a4-14*a8*a4+7*a8^2+63*a4^2+105*a5*a4+14*a10*a3-21*a9*a8-49*a9*a4+105*b4*a6+63*b5*a6-84*b6*a4+28*b6*a8-7*b10*a8+42*b6*a9+7*b10*a9-14*b4*a10-84*b6*a5-91*b8*a6-14*b5*a10"),
		SMQP("-14*a5^2+14*a5*a9-14*a5*a8+14*a6*a3-7*a6*a7+7*b9*a6+7*a8*a4-7*a4^2+7*a5*a4-14*b5*a6+7*b6*a4-7*b6*a8+7*b4*a10+7*b6*a5-7*b8*a6"),
		SMQP("252*a5^2-98*b10*a5-448*a5*a9-70*a9^2+616*a5*a8-42*a6*a3+154*a6*a7-49*b9*a6-35*b10*a4-371*a8*a4-35*a8^2+189*a4^2-399*a5*a4-133*a10*a3+105*a9*a8+182*a9*a4+105*b4*a6-63*b5*a6-84*b6*a4+301*b6*a8+35*b10*a8-147*b6*a9-35*b10*a9-371*b4*a10+42*b6*a5+266*b8*a6+133*b5*a10"),
		SMQP("-7*b6*a10+7*b10*a6+7*a10*a4-7*a6*a4+7*a6*a5-14*a6*a8+7*a6*a9"),
		SMQP("-7*b7*a5"),
		SMQP("7*b5*a5+7*a9*a7+7*a4*a7+7*b9*a4-7*a8*a7+7*a4*a3-7*a5*a3-7*a9*a3-28*b5*a4-7*b7*a6+7*b5*a8+14*b8*a9+7*b4*a9-21*b8*a5+7*b9*a5+7*b4*a5-14*b8*a8+7*b9*a8-7*b9*a9"),
		SMQP("63*b5*a5+7*a9*a7+7*a4*a7-56*a5*a7+7*b9*a4+56*a8*a3-7*a8*a7-49*a4*a3+49*a5*a3-7*a9*a3+28*b5*a4-7*b7*a6-49*b5*a8+14*b8*a9+7*b4*a9+35*b8*a5-49*b9*a5+7*b4*a5-14*b8*a8+7*b9*a8-7*b9*a9"),
		SMQP("-252*a5^2+70*b10*a5+392*a5*a9+14*a9^2-476*a5*a8+210*a6*a3-182*a6*a7+161*b9*a6+7*b10*a4+175*a8*a4+7*a8^2-189*a4^2+231*a5*a4+77*a10*a3-21*a9*a8+14*a9*a4-21*b4*a6-189*b5*a6+168*b6*a4-161*b6*a8-7*b10*a8-21*b6*a9+7*b10*a9+175*b4*a10+42*b6*a5-154*b8*a6-77*b5*a10"),
		SMQP("-7*b7*a5"),
		SMQP("-63*a5^2-7*b10*a5+49*a5*a9-14*a9^2-28*a5*a8+42*a6*a3-7*a6*a7+28*b9*a6-7*b10*a4+14*a8*a4-7*a8^2+21*a5*a4-14*a10*a3+21*a9*a8-14*a9*a4+21*b4*a6-63*b5*a6+21*b6*a4-28*b6*a8+7*b10*a8+21*b6*a9-7*b10*a9+14*b4*a10+21*b6*a5-35*b8*a6+14*b5*a10"),
		SMQP("112*b4*a8-175*b5*a5-7*a9*a7-7*a4*a7+168*a5*a7-119*b9*a4-56*a8*a3+7*a8*a7+49*a4*a3-49*a5*a3-105*a9*a3-28*b5*a4+119*b7*a6+112*b8*a4+161*b5*a8-14*b8*a9-231*b4*a9-147*b8*a5+161*b9*a5+105*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("56*b4*a8-133*b5*a5-21*a9*a7-21*a4*a7+112*a5*a7-77*b9*a4-56*a8*a3+21*a8*a7+35*a4*a3-35*a5*a3-35*a9*a3+28*b5*a4+77*b7*a6+56*b8*a4+91*b5*a8-42*b8*a9-133*b4*a9-49*b8*a5+91*b9*a5+35*b4*a5+42*b8*a8-21*b9*a8+21*b9*a9"),
		SMQP("7*b4*a8+7*b4*a4+7*b7*a6-14*b4*a9+7*b4*a5"),
		SMQP("7*b4*a8+7*b4*a4+7*b7*a6-14*b4*a9+7*b4*a5"),
		SMQP("28*b4*a8-35*b5*a5-7*a9*a7-7*a4*a7+28*a5*a7-35*b9*a4+7*a8*a7-7*a4*a3+7*a5*a3-21*a9*a3+28*b5*a4+35*b7*a6+28*b8*a4+21*b5*a8-14*b8*a9-63*b4*a9-7*b8*a5+21*b9*a5+21*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("14*b10*a5+28*a5*a9-14*a9^2-28*a5*a8-42*a6*a3+14*a6*a7+7*b9*a6-7*b10*a4-7*a8*a4-7*a8^2+21*a4^2+21*a5*a4+7*a10*a3+21*a9*a8-14*a9*a4+21*b4*a6+21*b5*a6-7*b6*a8+7*b10*a8+21*b6*a9-7*b10*a9-7*b4*a10-42*b6*a5-14*b8*a6-7*b5*a10"),
		SMQP("-14*b7*a5"),
		SMQP("-7*b10*a5-14*a5*a9-14*a9^2+35*a5*a8+14*a6*a7+7*b9*a6-7*b10*a4-28*a8*a4-7*a8^2+21*a4^2-21*a5*a4-14*a10*a3+21*a9*a8+7*a9*a4+21*b4*a6-21*b5*a6+14*b6*a8+7*b10*a8-7*b10*a9-28*b4*a10+7*b8*a6+14*b5*a10"),
		SMQP("-7*b6*a10+7*b10*a6+7*a10*a4-7*a6*a4+7*a6*a5-14*a6*a8+7*a6*a9"),
		SMQP("-7*b4*a8+7*b4*a4+14*b4*b6-7*b9*a4+14*b5*a4+7*b7*a6-7*b8*a4+14*b4*a9-7*b4*a5+7*b8*b6+7*b9*b6-14*b5*b6-14*b10*b4"),
		SMQP("-7*b7*a4+7*b7*a5-7*b4*a3+7*a3*b9-7*a7*b9+14*b8^2+7*b9^2+7*a7*b8-21*b5*a3+7*b5*a7+7*b7*b6+7*b5*b9+28*b4*b5-14*b4*b9-7*b8*b5-21*b8*b9"),
		SMQP("-7*b7*a4-7*b7*a5-14*b4*a3-7*a3*b8-7*a3*b9+7*a7*b9-14*b8^2-7*b9^2-7*b5^2-7*a7*b8+28*b5*a3-7*b5*a7+7*b7*b6-7*b5*b9-7*b4*b5-28*b4*b9+14*b8*b5+21*b8*b9+28*b4*b8+14*b4*a7"),
		SMQP("28*b4*a8-14*b4*a4-42*b5*a5+14*b4*b6+28*b9*a4-7*b5*a4-21*b7*a6-7*b8*a4-28*b5*a8+35*b5*a9+14*b8*a9-35*b4*a9+21*b8*a5-7*b9*a5-7*b4*a5-7*b8*a8+7*b9*a8-35*b8*b6+14*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9-14*b5*b6+14*b10*b4-14*b9*a9"),
		SMQP("14*b4^2-7*b7*a4-14*b7*a5-21*b4*a3-7*a3*b8-14*a3*b9+14*a7*b9-28*b8^2-14*b9^2-7*b5^2-14*a7*b8+49*b5*a3-14*b5*a7+7*b7*b6-14*b5*b9-35*b4*b5-21*b4*b9+21*b8*b5+42*b8*b9+28*b4*b8+21*b4*a7"),
		SMQP("28*b4^2-7*b7*a4-21*b7*a5-35*b4*a3-21*a3*b9+21*a7*b9-42*b8^2-21*b9^2-21*a7*b8+63*b5*a3-21*b5*a7+7*b7*b6-21*b5*b9-56*b4*b5-14*b4*b9+21*b8*b5+63*b8*b9+28*b4*b8+28*b4*a7"),
		SMQP("7*b4^2-7*b7*a5-7*b4*a3-7*a3*b9+7*a7*b9-14*b8^2-7*b9^2-7*a7*b8+21*b5*a3-7*b5*a7-7*b5*b9-21*b4*b5+7*b8*b5+21*b8*b9+7*b4*b8+7*b4*a7"),
		SMQP("-14*b4*a8+28*b4*a4+14*b5*a5+14*b4*b6-14*b9*a4-7*b5*a4+21*b7*a6-7*b8*a4+14*b5*a8-7*b5*a9+14*b8*a9+7*b4*a9-7*b8*a5+7*b9*a5+7*b4*a5-7*b8*a8+7*b9*a8+7*b8*b6+14*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9-14*b5*b6-28*b10*b4-14*b9*a9"),
		SMQP("7*b4*a3-7*a3*b8-7*b5^2+7*b5*a3-7*b4*b5+7*b8*b5"),
		SMQP("-7*b4*a8+7*b4*a4+7*b4*b6-7*b9*a4+7*b5*a4+7*b7*a6+7*b4*a9+7*b9*b6-7*b5*b6-7*b10*b4"),
		SMQP("-14*b4*a8+28*b4*a4+42*b5*a5-28*b4*b6-14*b9*a4-49*b5*a4+21*b7*a6-7*b8*a4+56*b5*a8-49*b5*a9+14*b8*a9+7*b4*a9+21*b8*a5-7*b9*a5-7*b4*a5-7*b8*a8+7*b9*a8+49*b8*b6-28*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9+28*b5*b6-28*b10*b4-14*b9*a9"),
		SMQP("70*b4*a8-56*b4*a4-42*b5*a5+14*b4*b6+70*b9*a4-7*b5*a4-63*b7*a6-91*b8*a4+14*b5*a8-7*b5*a9+14*b8*a9+7*b4*a9+105*b8*a5-49*b9*a5-133*b4*a5-7*b8*a8+7*b9*a8+91*b8*b6-70*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9-14*b5*b6-28*b10*b4-14*b9*a9"),
		SMQP("-7*b7*a4-7*b7*a5-28*b4*a3+7*a3*b8+21*a3*b9-7*a7*b9+14*b8^2+7*b9^2+21*b5^2+7*a7*b8-42*b5*a3+7*b5*a7+7*b7*b6-7*b5*b9+49*b4*b5-14*b4*b9-14*b8*b5-21*b8*b9"),
		SMQP("14*b4*a8-28*b4*a4-35*b4*b6+14*b9*a4-14*b5*a4-21*b7*a6+28*b8*a4-14*b5*a8+7*b5*a9-14*b8*a9-28*b4*a9-21*b8*a5+7*b9*a5+28*b4*a5+7*b8*a8-7*b9*a8-28*b8*b6-14*b9*b6-7*b8*b10-7*b10*b5+7*b10*b9+35*b5*b6+49*b10*b4+14*b9*a9"),
		SMQP("-7*b7*a4-7*b7*a5-14*b4*a3-7*a3*b8-7*a3*b9+7*a7*b9-14*b8^2-7*b9^2-7*b5^2-7*a7*b8+28*b5*a3-7*b5*a7+7*b7*b6-7*b5*b9-7*b4*b5-28*b4*b9+14*b8*b5+21*b8*b9+28*b4*b8+14*b4*a7"),
		SMQP("14*b8*a6-7*b4*a10+7*b4*a6-7*b5*a6-7*b9*a6"),
		SMQP("14*b4*a8-28*b4*a4+42*b5*a5-14*b4*b6+14*b9*a4+7*b5*a4-21*b7*a6+7*b8*a4-14*b5*a8+7*b5*a9-14*b8*a9-7*b4*a9-105*b8*a5+49*b9*a5+49*b4*a5+7*b8*a8-7*b9*a8-7*b8*b6-14*b9*b6-7*b8*b10-7*b10*b5+7*b10*b9+14*b5*b6+28*b10*b4+14*b9*a9"),
		SMQP("-7*a3*b7+7*b7*b5"),
		SMQP("-28*b4*a8+42*b4*a4+14*b5*a5+42*b4*b6-28*b9*a4+21*b5*a4+35*b7*a6-21*b8*a4+14*b5*a8-7*b5*a9+14*b8*a9+35*b4*a9-7*b8*a5+7*b9*a5-7*b4*a5-7*b8*a8+7*b9*a8+21*b8*b6+28*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9-42*b5*b6-56*b10*b4-14*b9*a9"),
		SMQP("-7*b4*a6+7*b4*a10+7*b5*a6+7*b9*a6-14*b8*a6"),
		SMQP("7*b4*a8-7*b4*a4-14*b4*b6+7*b9*a4-14*b5*a4-7*b7*a6+7*b8*a4-14*b4*a9+7*b4*a5-7*b8*b6-7*b9*b6+14*b5*b6+14*b10*b4"),
		SMQP("-7*a3*b7+7*b7*b5"),
		SMQP("-7*b4*a3+7*a3*b8+7*b5^2-7*b5*a3+7*b4*b5-7*b8*b5"),
		SMQP("-7*b7*a5-7*b4*a3+7*a3*b9+7*b5^2-7*b5*a3-7*b5*b9+7*b4*b5"),
		SMQP("-7*b4*a6+7*b4*a10+7*b5*a6+7*b9*a6-14*b8*a6"),
		SMQP("-7*a3*b7+7*b7*b5"),
		SMQP("-14*b4*a8+28*b4*a4+14*b4*b6-14*b9*a4-7*b5*a4+21*b7*a6-7*b8*a4+14*b5*a8-7*b5*a9+14*b8*a9+7*b4*a9+21*b8*a5-7*b9*a5-7*b4*a5-7*b8*a8+7*b9*a8+7*b8*b6+14*b9*b6+7*b8*b10+7*b10*b5-7*b10*b9-14*b5*b6-28*b10*b4-14*b9*a9"),
		SMQP("-7*b7*a5-7*b4*a3-7*a3*b9+7*a7*b9-14*b8^2-7*b9^2-7*a7*b8+21*b5*a3-7*b5*a7-7*b5*b9-14*b4*b5-7*b4*b9+7*b8*b5+21*b8*b9+14*b4*b8+7*b4*a7"),
		SMQP("14*b8*a6-7*b4*a10+7*b4*a6-7*b5*a6-7*b9*a6"),
		SMQP("7*b7*a4-35*b7*a5-7*b4*a3-14*a3*b8+21*a3*b9+7*a7*b9-14*b8^2-7*b9^2+14*b5^2-7*a7*b8+7*b5*a3-7*b5*a7-7*b7*b6-35*b5*b9-14*b4*b5+14*b4*b9+21*b8*b5+21*b8*b9"),
		SMQP("7*b7*a5"),
		SMQP("-84*a5^2-7*b10*a5+70*a5*a9-14*a9^2-49*a5*a8+84*a6*a3-28*a6*a7+49*b9*a6-7*b10*a4+14*a8*a4-7*a8^2-21*a4^2+21*a5*a4-14*a10*a3+21*a9*a8+7*a9*a4+21*b4*a6-105*b5*a6+42*b6*a4-28*b6*a8+7*b10*a8-7*b10*a9+14*b4*a10+42*b6*a5-35*b8*a6+14*b5*a10"),
		SMQP("-7*b6*a10+7*b10*a6+7*a10*a4-7*a6*a4+7*a6*a5-14*a6*a8+7*a6*a9"),
		SMQP("-7*b5*a5-7*a9*a7-7*a4*a7+56*a5*a7-7*b9*a4+7*a8*a7+49*a4*a3-49*a5*a3-49*a9*a3-28*b5*a4+7*b7*a6-7*b5*a8+56*b5*a9-14*b8*a9-7*b4*a9-35*b8*a5-7*b9*a5+49*b4*a5+14*b8*a8-7*b9*a8+7*b9*a9"),
		SMQP("-14*a5^2+14*a5*a9-14*a5*a8+14*a6*a3-7*a6*a7+7*b9*a6+7*a8*a4-7*a4^2+7*a5*a4-14*b5*a6+7*b6*a4-7*b6*a8+7*b4*a10+7*b6*a5-7*b8*a6"),
		SMQP("1701*a5*b5*a8-1071*b4*a5*a4-434*b8*b4*a10-966*b9*b6*a9-728*b10*b4*a5-504*b7*a10*a5-70*b10*b9*a9+840*b8*a6*a3-154*b8*a9*a4+504*b7*a10*a4-1323*b7*a6*a4+504*b7*a10*a9+2520*a5*a8*a7-315*a4*a3*a8+420*b8*b6*a5+2520*a8^2*a3-3108*b5*a4*a9+252*b7*a6*a5+1512*b8*a8*a9-224*b8*a8^2+1323*a6*b7*a8+126*a4*a7*a8-980*b8*a8*a4-2331*b5*a8^2+315*a8^2*a7-1834*b8*a10*a3+756*a7*a9^2-441*a7*a4^2+1092*b5*a4^2+546*b4*a4^2-2800*b4*a9^2+441*a8*b9*a9+1407*a8*b4*a9-1071*a8*a9*a7+6636*a8*b5*a9-4725*a8*a9*a3-2016*a7*a5*a4+756*b9*a10*a7-504*a8*b7*a10+3717*a9*a4*a3-1428*a9^2*b5+469*a9*b9*a4+315*a9*a4*a7+112*a9^2*b9+756*a9*a5*a3+770*b6*b9*a8+63*a4^2*a3-70*b4*a10*a3-70*b10*b9*a4-434*b4^2*a10-63*a3*a5*a4-1008*b4*b6*a4-504*b7*b6*a10+504*b7*b6*a6+966*b8*b6*a9-434*b8*b10*a9+630*b4*b6*a9-1561*a4*b4*a9+2590*b8*b5*a10+2268*b8*b9*a10+1869*b9*a4^2-756*b5^2*a6-1428*b4*a6*a3+3416*b8^2*a6+2296*b8*a6*a7-154*b4*b10*a8-84*b10*b5*a9-280*b9*a6*a7-574*b4*b9*a10-336*b9*a6*a3-2646*b9*b8*a6+252*b5^2*a10-434*b8*b10*a4+434*b4*b8*a6-2198*b4*b5*a10-5992*b8*a5*a9-868*b4*a5*a9-84*b10*b5*a4-686*b9*a10*a3-4564*b4*b9*a6-672*b8*b6*a4+742*b8*b6*a8-784*b9*b10*a5+742*b4*b6*a8+742*b4*a8*a4-728*b8*b10*a5+5670*b4*b5*a6-756*a9^2*a3-462*b4^2*a6-728*a5*b9*a9+1036*b4*a6*a7+504*a5^2*b4-252*b7*a6*a9+658*b9*b10*a8+588*b9*b6*a5-840*b9*b6*a4+5089*b4*a5*a8-399*b8*a5*a4-1512*b8^2*a10-910*b9*a8*a4-2919*b9*a5*a4-868*b8*a9^2+490*b4*b10*a4+504*b8*a5^2+714*b8*a4^2+420*b4*b6*a5+154*b4*a8^2-3024*b5*a6*a3+1260*b5*a6*a7+3402*b5*b9*a6-1869*b5*a8*a4+2016*b5*a10*a3-756*b5*a10*a7-826*b5*b9*a10+63*b5*a5*a4-1638*b5*b8*a6-756*a5*b5*a9-973*b9*a8^2+833*a5*b9*a8+4333*a5*b8*a8-756*b8*a10*a7+490*b10*b4*a9+2520*b9*a5^2+420*b5*b6*a4-266*b9^2*a6-756*b9^2*a10-154*b8*b10*a8+420*b5*b6*a9-1701*a5*a3*a8"),
		SMQP("74774*a5*b5*a8+9681*b4*a5*a4-5250*b8*b4*a10+15624*b9*b6*a9-7952*b10*b4*a5+1008*b7*a10*a5-7217*b10*b9*a9-1960*b8*a6*a3+315*b8*a9*a4+17836*a6*a3*a7+2520*b7*a10*a4-567*b7*a6*a4-7056*b7*a10*a9+9996*b6*a4*a7-15008*b6*a8*a7+12768*b6*a5*a3-11228*a5*a8*a7-5642*a4*a3*a8-19740*b6*a4*a3+22974*b8*b6*a5+25018*b6*a8*a3-14378*a8^2*a3+11242*b5*a4*a9+6552*b6*b5*a5-13944*b6*a5*a7-6258*b7*a6*a5-3087*b8*a8*a9+3927*b8*a8^2-3843*a6*b7*a8+2464*a4*a7*a8-2604*b8*a8*a4-4942*b5*a8^2+9100*a8^2*a7-2814*b8*a10*a3-2716*a7*a9^2-14112*a7*a4^2+19425*b5*a4^2+3318*b4*a4^2+15911*b4*a9^2-4305*a8*b9*a9-13377*a8*b4*a9-1162*a10*a3^2+7756*b4*a10*a7-8400*a8*a9*a7-8673*a8*b5*a9+24234*a8*a9*a3+1176*a7*a5*a4-8288*a7^2*a6-6552*b9*a10*a7+3836*a10*a3*a7+3024*a8*b7*a10-20188*a9*a4*a3+5929*a9^2*b5-1036*a9*b9*a4+10136*a9*a4*a7+8120*a9*a5*a7-2086*a9^2*b9-16660*a9*a5*a3-1008*a7^2*a10-5852*b6*b9*a8+14742*a4^2*a3-6972*a3^2*a6-8484*b4*a10*a3+2359*b10*b9*a4-3458*b4^2*a10-1218*a3*a5*a4-1890*b4*b6*a4-1008*b7*b6*a10-6552*b7*b6*a6+504*b7*b10*a10+8568*a5^2*a3-2919*b8*b6*a9+2793*b8*b10*a9-8617*a4*b4*a9-7266*b8*b5*a10-2520*b8*b9*a10+504*b7*b10*a6-4368*b9*a4^2-8722*b10*b5*a5+6804*b5^2*a6+8442*b4*a6*a3-2688*b8^2*a6-5236*b8*a6*a7-2240*b10*a9*a7+1778*b4*b10*a8+1288*b10*a4*a7+1232*b10*a8*a7+8351*b10*b5*a9+8400*b9*a6*a7+7042*b4*b9*a10-20902*b9*a6*a3+9002*b9*b8*a6+9968*b5^2*a10-231*b8*b10*a4+1778*b4*b8*a6-24738*b4*b5*a10+3633*b8*a5*a9-3437*b10*b5*a8+2261*b4*a5*a9-5257*b10*b5*a4+7070*b9*a10*a3-3724*b4*b9*a6-5817*b8*b6*a4-266*b10*a9*a3-1400*b10*a5*a3-5607*b8*b6*a8+1666*b9*b10*a5+2044*b4*b6*a8-4844*b4*a8*a4-3234*b6*a9*a3+6720*b6*a9*a7-742*b10*a8*a3-6090*b8*b10*a5-11130*b4*b5*a6+1736*a9^2*a3+1050*b4^2*a6+238*b10*a4*a3-7105*a5*b9*a9-9380*b4*a6*a7-2058*a5^2*b4+5733*b7*a6*a9+287*b9*b10*a8-10836*b9*b6*a5+378*b9*b6*a4+1015*b4*a5*a8-15603*b8*a5*a4+2016*b8^2*a10+4207*b9*a8*a4-1911*b9*a5*a4+42*b8*a9^2-266*b4*b10*a4+13482*b8*a5^2+6195*b8*a4^2+3556*b10*a5*a7-756*b4*b6*a5+364*b4*a8^2+16926*b5*a6*a3+12908*b5*a6*a7-20006*b5*b9*a6-18949*b5*a8*a4-19390*b5*a10*a3+1708*b5*a10*a7+10066*b5*b9*a10-39648*b5*a5*a4+41986*b5*b8*a6-51212*a5*b5*a9+4501*b9*a8^2+6097*a5*b9*a8-945*a5*b8*a8+6552*b8*a10*a7+1750*b10*b4*a9-3066*b9*a5^2-4662*b5*b6*a4-6706*b9^2*a6-1512*b9^2*a10+21854*b5*b6*a8+2625*b8*b10*a8+28476*a5^2*b5-17766*b5*b6*a9+12040*a5*a3*a8-2016*a5^2*a7"),
		SMQP("44114*a5*b5*a8+8715*b4*a5*a4+1106*b8*b4*a10+23772*b9*b6*a9-7056*b10*b4*a5+15120*b7*a10*a5-9667*b10*b9*a9+2296*b8*a6*a3+217*b8*a9*a4-6748*a6*a3*a7+7560*b7*a10*a4-11781*b7*a6*a4-15120*b7*a10*a9+6216*b6*a4*a7-20804*b6*a8*a7+15708*b6*a5*a3+5404*a5*a8*a7-34846*a4*a3*a8-5712*b6*a4*a3+16926*b8*b6*a5+12110*b6*a8*a3-10990*a8^2*a3+1750*b5*a4*a9+17052*b6*b5*a5-10920*b6*a5*a7-30030*b7*a6*a5-8925*b8*a8*a9+5621*b8*a8^2+13167*a6*b7*a8+16072*a4*a7*a8+4508*b8*a8*a4-6370*b5*a8^2+10612*a8^2*a7+4942*b8*a10*a3-9268*a7*a9^2-26712*a7*a4^2+6363*b5*a4^2+9954*b4*a4^2+38157*b4*a9^2-8547*a8*b9*a9-15435*a8*b4*a9-7742*a10*a3^2+18340*b4*a10*a7-17472*a8*a9*a7-10563*a8*b5*a9+46830*a8*a9*a3-31080*a7*a5*a4-224*a7^2*a6-15624*b9*a10*a7+14420*a10*a3*a7-1008*a8*b7*a10-38108*a9*a4*a3+16051*a9^2*b5+3220*a9*b9*a4+22736*a9*a4*a7+13160*a9*a5*a7-6482*a9^2*b9-32060*a9*a5*a3-5040*a7^2*a10-6664*b6*b9*a8+35658*a4^2*a3+9996*a3^2*a6-21196*b4*a10*a3+3437*b10*b9*a4-10374*b4^2*a10+40530*a3*a5*a4-5670*b4*b6*a4+9072*b7*b6*a10-15876*b7*b6*a6-2520*b7*b10*a10+19656*a5^2*a3+5187*b8*b6*a9-1813*b8*b10*a9-7308*b4*b6*a9-19299*a4*b4*a9-4942*b8*b5*a10-17640*b8*b9*a10+1512*b7*b10*a6-15120*b9*a4^2-12838*b10*b5*a5+21924*b5^2*a6+17094*b4*a6*a3-13664*b8^2*a6-13076*b8*a6*a7-5264*b10*a9*a7+3318*b4*b10*a8+3304*b10*a4*a7+224*b10*a8*a7+18893*b10*b5*a9+8680*b9*a6*a7+17822*b4*b9*a10-25802*b9*a6*a3+32942*b9*b8*a6-7672*b5^2*a10+707*b8*b10*a4+1134*b4*b8*a6-15022*b4*b5*a10+31003*b8*a5*a9-2135*b10*b5*a8-47817*b4*a5*a9-8323*b10*b5*a4+16450*b9*a10*a3-13356*b4*b9*a6-4011*b8*b6*a4+1106*b10*a9*a3-6328*b10*a5*a3-14749*b8*b6*a8+6566*b9*b10*a5+5124*b4*b6*a8-12516*b4*a8*a4-14154*b6*a9*a3+20580*b6*a9*a7+1918*b10*a8*a3-2926*b8*b10*a5-12390*b4*b5*a6+9016*a9^2*a3+3150*b4^2*a6-1414*b10*a4*a3-1547*a5*b9*a9-20412*b4*a6*a7+22722*a5^2*b4+12663*b7*a6*a9-3059*b9*b10*a8-11592*b9*b6*a5+4242*b9*b6*a4+15981*b4*a5*a8-21609*b8*a5*a4+6048*b8^2*a10+20909*b9*a8*a4+3003*b9*a5*a4+1918*b8*a9^2-798*b4*b10*a4-13986*b8*a5^2+945*b8*a4^2+8092*b10*a5*a7+168*b4*b6*a5-1932*b4*a8^2+66066*b5*a6*a3-13300*b5*a6*a7-51002*b5*b9*a6-8743*b5*a8*a4-34482*b5*a10*a3+14308*b5*a10*a7+22862*b5*b9*a10-23184*b5*a5*a4+22582*b5*b8*a6-33236*a5*b5*a9+1295*b9*a8^2+3563*a5*b9*a8-6139*a5*b8*a8+9576*b8*a10*a7+1218*b10*b4*a9-4494*b9*a5^2+6510*b5*b6*a4-3542*b9^2*a6+1512*b9^2*a10+22022*b5*b6*a8+5971*b8*b10*a8+12852*a5^2*b5-42378*b5*b6*a9-21784*a5*a3*a8-2016*a5^2*a7"),
		SMQP("1008*b9^2*a9+1008*b9^2*a4+1512*b9^2*a8+140*b7*b4*a10-3220*b7*b10*a9+8596*b7*a10*a3-924*b7*b4*a6+812*b7*b10*a4+18480*a9*a3*a7+8064*b4*a4*a3+10878*b8*a8*a3-672*b8*b9*a4-168*b8*a4*a3+2562*b8*a9*a7-13545*b7*a6*a3-7560*b4*a8*a3-504*b4*a8*a7+4032*b5^2*a5+5544*b4*b5*a4+7511*a7*b7*a6+30807*a7*a8*a3+9744*a7*b5*a9-1071*a7*b5*a5-3360*b8*b7*a10-7119*a7*b4*a5-6048*a7*b7*a10+15561*a7*b4*a9+6909*a7*b8*a5-6090*a7*b8*a8+546*b8*a9*a3-24843*b8*a5*a3+336*b8*b4*a9-3192*b8*b4*a5-1848*b8*b5*a4+2184*a7^2*a5-55062*a7*a4*a3+15057*a7^2*a4+1008*b4*b8*a8-1176*b4*b9*a8-5544*b4*a4*a7+8568*b4*b5*a8-5775*a7^2*a9+1344*a7*b8*a4-10017*a9*a3^2-4452*b7*a4^2+5684*b7*a4*a8+672*b8^2*a5+20097*a3*b4*a5+52857*a3^2*a4+7623*a3^2*a5-25095*a3*b4*a9-1232*b7*a5*a9-4536*b9^2*a5+13608*b7*a5*a8-7812*b5*a4*a3+5516*b7*b5*a10+2772*b5*a4*a7+504*b5^2*a4-336*b5^2*a8-8001*a5*a3*a7-2688*b5^2*a9-6321*a7^2*a8-7560*b4*b9*a4-41832*a3^2*a8-4536*b5*b9*a5-2184*b5*b8*a5-4536*b5*b4*a5-3528*b4^2*a9-1848*b5*b8*a8+5320*b7*b8*a6+504*b4^2*a5-2436*b7*a9*a8-7056*b7*b6*a5-3108*b7*b5*a6+8344*b7*a9*a4-2016*b7*b6*a4+5628*b7*b6*a9-812*b7*b10*a8-504*b5*b9*a8+3332*b7*b6*a8-3752*b7*a9^2-4200*b7*a5^2+2408*b7*b10*a5+5040*b5*b9*a9-1876*b7*a8^2-9660*a5*b7*a4-2016*b8^2*a9+1008*b8^2*a8-2016*b9*b5*a4-24528*b5*a9*a3+4032*b5*b8*a9-7728*b5*b4*a9+2016*b8^2*a4+3087*b9*a9*a3+9009*b9*a5*a3-4032*b9*b8*a9-1512*b9*b8*a8-1512*b9*b4*a5+7560*b9*b8*a5+8568*b9*b4*a9+17409*b5*a8*a3-6111*b5*a8*a7-63*b5*a5*a3+273*b9*a5*a7+7497*b9*a8*a3-10143*b9*a4*a3-3297*b9*a9*a7+3969*b9*a4*a7-3255*b9*a8*a7-672*b9*b7*a10+532*b9*b7*a6"),
		SMQP("-47271*a6*a3*a4+13335*a6*a3*a5-90132*a6*b5*b6-2331*a6*b5*a5-37632*a6*a3*b6-23436*a6*b4*b6-71757*a6*b4*a5-49350*a6*b5*a4+27216*a4*b6*a5-61320*a6*b4*a4-45528*a4^2*a5-6720*a4*b6^2-14280*a4*a5^2-30912*b6^2*a5-27888*b6*a5^2-23184*a4^3+4368*a4^2*b6-34272*a5^3+22400*b10*a5^2+55146*b10*b5*a6-16632*b8*a10*a5-31696*b10*b6*a9-28616*b10*b5*a10+9576*b10*b9*a10-13664*b10*a6*a7+17472*b10*a6*a3+5824*b6*a8*a4-2464*b6*b10*a4-17521*a6*a8*a7-31416*a9*b5*a10+11984*b6*b10*a8-37023*a6*a9*a7+50897*a6*a4*a7+5040*b10*a10*a7-6048*b6*a10*a7+3808*b6*a10*a3+58765*a6*b9*a4+8232*a6*a8*a3+4424*b10*a9*a8-3248*b10*a9*a4+22064*b6*a6*a7+5264*b10*a9^2-25704*a6*b7*a10-61173*a6*b9*a9+47803*a6*b9*a8-38808*b6*a9*a8-18536*b10*a8*a4-48874*b10*b9*a6+43568*b5*b6*a10+99435*a6*b4*a9-952*a8^3-2828*a6*b8*a8+7840*b6*a9^2+6552*b9*a10*a8-17528*b6*a8^2+5152*a8*a9^2+40712*a8*a4^2+39200*a8*b4*a10-23800*a8*b5*a10+42448*a9^2*a5+11088*a8*a10*a7+15064*a8^2*a4-12488*a8*a10*a3+25263*a6*a9*a3+8736*a9^3-24080*a9*a4^2-2856*a9*a8^2+14896*a9*a8*a4+1680*a9*a10*a3-46536*a9*b4*a10+9072*a10*a9*a7-28224*a10*a4*a7+40712*a10*b5*a4-44352*a10*b9*a4+43456*a10*a4*a3+4536*a10*b8*a8+6048*a10^2*b7+26208*a10*b9*a9+23744*b6*b10*a5-14042*b8*a6*a4+47488*a9*b6*a4-14196*a8*b4*a6-12488*a5*b10*a8+6608*b10*a4^2+5432*b10*a8^2-5600*b10^2*a5+42952*a5*b6*a9-27720*a4*b8*a10+65408*a4*b4*a10-13608*b9*a10*a5-39536*a9^2*a4+25144*a5*b10*a4+1456*b10^2*a4-17584*b6^2*a8+39648*b6^2*a9-3024*b8*a10*a9-13146*b8*a6*a9+1960*a5*b10*a9+9520*b10^2*a9-11088*a5*a10*a7-36358*b8*b10*a6-60592*a5^2*a8-4949*a5*b9*a6+63364*b9*b6*a6+87248*b4*a10*a5-35336*a5*a8^2+8358*a6*b5*a9+37149*a6*b5*a8+37485*a6^2*b7+12936*b4*b10*a6-4480*b10^2*a8-4144*b10*a10*a3+6104*a5*b5*a10+28168*a5*a10*a3+91280*a5*a9*a4+11242*b8*b6*a6-59192*a5*a9*a8+15624*b8*b10*a10-27797*b8*a6*a5+57680*b4*b6*a10-13608*b8*b6*a10-18032*b4*b10*a10+61600*a5*a8*a4-42056*a5*b6*a8+26656*a5^2*a9+3024*b9*b6*a10-4312*a5*a6*a7"),
		SMQP("-26082*a5*b5*a8-8715*b4*a5*a4-10514*b8*b4*a10-6972*b9*b6*a9-1456*b10*b4*a5-10080*b7*a10*a5+2331*b10*b9*a9-12936*b8*a6*a3-2065*b8*a9*a4-84*a6*a3*a7-8568*b7*a10*a4+8757*b7*a6*a4+1008*b7*a10*a9-1512*b6*a4*a7+10500*b6*a8*a7-5292*b6*a5*a3+9492*a5*a8*a7+41622*a4*a3*a8+1008*b6*a4*a3-12726*b8*b6*a5-9030*b6*a8*a3-8778*a8^2*a3+3906*b5*a4*a9+1092*b6*b5*a5+3528*b6*a5*a7+6174*b7*a6*a5+2037*b8*a8*a9+931*b8*a8^2-4599*a6*b7*a8-16464*a4*a7*a8-13076*b8*a8*a4+9114*b5*a8^2-84*a8^2*a7+2450*b8*a10*a3+3108*a7*a9^2+14112*a7*a4^2-819*b5*a4^2+10710*b4*a4^2-5285*b4*a9^2-357*a8*b9*a9-3045*a8*b4*a9+8022*a10*a3^2-5124*b4*a10*a7+3024*a8*a9*a7-3549*a8*b5*a9-8694*a8*a9*a3-1512*a7*a5*a4+672*a7^2*a6+6552*b9*a10*a7-9492*a10*a3*a7+7056*a8*b7*a10+1092*a9*a4*a3-5355*a9^2*b5+5796*a9*b9*a4-1176*a9*a4*a7-16296*a9*a5*a7+882*a9^2*b9+22764*a9*a5*a3+3024*a7^2*a10-1512*b6*b9*a8-25578*a4^2*a3-4284*a3^2*a6+7756*b4*a10*a3+483*b10*b9*a4-3682*b4^2*a10-3906*a3*a5*a4-7098*b4*b6*a4+2016*b7*b6*a10+6300*b7*b6*a6+504*b7*b10*a10-11592*a5^2*a3+21*b8*b6*a9-371*b8*b10*a9+3444*b4*b6*a9-1253*a4*b4*a9-3458*b8*b5*a10+8568*b8*b9*a10-5544*b7*b10*a6+6552*b9*a4^2+6846*b10*b5*a5-17388*b5^2*a6-31206*b4*a6*a3-5152*b8^2*a6+11452*b8*a6*a7+4704*b10*a9*a7-3710*b4*b10*a8-1848*b10*a4*a7-1680*b10*a8*a7-5229*b10*b5*a9-5040*b9*a6*a7-3486*b4*b9*a10+5754*b9*a6*a3-5446*b9*b8*a6+6384*b5^2*a10-35*b8*b10*a4+15106*b4*b8*a6+11774*b4*b5*a10-11683*b8*a5*a9-2121*b10*b5*a8-3143*b4*a5*a9+651*b10*b5*a4-5586*b9*a10*a3-2828*b4*b9*a6-2037*b8*b6*a4-2730*b10*a9*a3+6216*b10*a5*a3+8533*b8*b6*a8-5838*b9*b10*a5+10724*b4*b6*a8+6188*b4*a8*a4+6930*b6*a9*a3-9828*b6*a9*a7+1722*b10*a8*a3+4774*b8*b10*a5-4410*b4*b5*a6+840*a9^2*a3+11802*b4^2*a6+1806*b10*a4*a3-1197*a5*b9*a9+19292*b4*a6*a7+10374*a5^2*b4-6111*b7*a6*a9+3003*b9*b10*a8+5208*b9*b6*a5-1890*b9*b6*a4+10619*b4*a5*a8+21945*b8*a5*a4-6405*b9*a8*a4-14595*b9*a5*a4+1778*b8*a9^2-1498*b4*b10*a4+6426*b8*a5^2+2079*b8*a4^2-6636*b10*a5*a7-5544*b4*b6*a5-700*b4*a8^2+126*b5*a6*a3+588*b5*a6*a7+13146*b5*b9*a6-5313*b5*a8*a4+6258*b5*a10*a3-4116*b5*a10*a7-7518*b5*b9*a10+11256*b5*a5*a4+11634*b5*b8*a6+18732*a5*b5*a9+1785*b9*a8^2+5061*a5*b9*a8+91*a5*b8*a8-6552*b8*a10*a7+3878*b10*b4*a9+1302*b9*a5^2+2898*b5*b6*a4+462*b9^2*a6-2520*b9^2*a10-12222*b5*b6*a8-2443*b8*b10*a8-1764*a5^2*b5+13650*b5*b6*a9-5880*a5*a3*a8+6048*a5^2*a7"),
		SMQP("-16506*a5*b5*a8-15687*b4*a5*a4-5250*b8*b4*a10-8820*b9*b6*a9+10192*b10*b4*a5-11088*b7*a10*a5+1239*b10*b9*a9+22680*b8*a6*a3-1869*b8*a9*a4+29484*a6*a3*a7-17640*b7*a10*a4+26145*b7*a6*a4+5040*b7*a10*a9+504*b6*a4*a7+756*b6*a8*a7-5292*b6*a5*a3-18396*a5*a8*a7+38934*a4*a3*a8-3024*b6*a4*a3+10794*b8*b6*a5+20538*b6*a8*a3-19530*a8^2*a3+3066*b5*a4*a9-15036*b6*b5*a5+5544*b6*a5*a7+13230*b7*a6*a5+3465*b8*a8*a9+3759*b8*a8^2-11403*a6*b7*a8-28224*a4*a7*a8-12684*b8*a8*a4+5250*b5*a8^2+7308*a8^2*a7-14910*b8*a10*a3+5796*a7*a9^2+20160*a7*a4^2+1449*b5*a4^2-12474*b4*a4^2-15169*b4*a9^2+2583*a8*b9*a9-11865*a8*b4*a9-1386*a10*a3^2-12852*b4*a10*a7+5040*a8*a9*a7-6489*a8*b5*a9-4662*a8*a9*a3+18648*a7*a5*a4-10080*a7^2*a6+5544*b9*a10*a7-7812*a10*a3*a7+7056*a8*b7*a10+27972*a9*a4*a3-4431*a9^2*b5+2100*a9*b9*a4-7560*a9*a4*a7-7560*a9*a5*a7+2730*a9^2*b9-6804*a9*a5*a3+3024*a7^2*a10-3360*b6*b9*a8-37674*a4^2*a3-8316*a3^2*a6+8204*b4*a10*a3-5313*b10*b9*a4+12334*b4^2*a10-56322*a3*a5*a4-2226*b4*b6*a4-5040*b7*b6*a10+9324*b7*b6*a6+4536*b7*b10*a10-3528*a5^2*a3-10647*b8*b6*a9+4641*b8*b10*a9+14028*b4*b6*a9+22127*a4*b4*a9-3234*b8*b5*a10+5544*b8*b9*a10-7560*b7*b10*a6+2016*b9*a4^2+12894*b10*b5*a5-13356*b5^2*a6-23478*b4*a6*a3+7392*b8^2*a6-1092*b8*a6*a7-2590*b4*b10*a8-504*b10*a4*a7+3024*b10*a8*a7-8169*b10*b5*a9-11424*b9*a6*a7-13566*b4*b9*a10+21546*b9*a6*a3-15078*b9*b8*a6+6384*b5^2*a10-2415*b8*b10*a4-2926*b4*b8*a6+15694*b4*b5*a10-21735*b8*a5*a9-5061*b10*b5*a8+38717*b4*a5*a9+10479*b10*b5*a4-8610*b9*a10*a3+8036*b4*b9*a6+63*b8*b6*a4-1386*b10*a9*a3-504*b10*a5*a3+105*b8*b6*a8-3822*b9*b10*a5-6020*b4*b6*a8+7924*b4*a8*a4+2898*b6*a9*a3-10836*b6*a9*a7-3654*b10*a8*a3-6426*b8*b10*a5-3402*b4*b5*a6-12600*a9^2*a3-10038*b4^2*a6-882*b10*a4*a3-2793*a5*b9*a9+8764*b4*a6*a7-21882*a5^2*b4-8883*b7*a6*a9+1911*b9*b10*a8+7224*b9*b6*a5+10710*b9*b6*a4+1183*b4*a5*a8+12453*b8*a5*a4+2016*b8^2*a10-8841*b9*a8*a4-8799*b9*a5*a4-294*b8*a9^2+4102*b4*b10*a4+28602*b8*a5^2+8379*b8*a4^2-2268*b10*a5*a7-14280*b4*b6*a5+3724*b4*a8^2-47250*b5*a6*a3+19740*b5*a6*a7+14154*b5*b9*a6-861*b5*a8*a4+32802*b5*a10*a3-14868*b5*a10*a7-6510*b5*b9*a10+14784*b5*a5*a4+5586*b5*b8*a6+18228*a5*b5*a9+2877*b9*a8^2+5481*a5*b9*a8-6153*a5*b8*a8+6552*b8*a10*a7-1946*b10*b4*a9-4746*b9*a5^2-4662*b5*b6*a4+1470*b9^2*a6-1512*b9^2*a10-4326*b5*b6*a8+2793*b8*b10*a8-5796*a5^2*b5+21546*b5*b6*a9+68040*a5*a3*a8+6048*a5^2*a7"),
		SMQP("3766*a5*b5*a8+4977*b4*a5*a4-2450*b8*b4*a10+13944*b9*b6*a9-3248*b10*b4*a5+1008*b7*a10*a5-6657*b10*b9*a9+2856*b8*a6*a3+539*b8*a9*a4+11340*a6*a3*a7+2520*b7*a10*a4-1575*b7*a6*a4-7056*b7*a10*a9+7308*b6*a4*a7-12432*b6*a8*a7+10080*b6*a5*a3-3612*a5*a8*a7-2394*a4*a3*a8-14364*b6*a4*a3+23478*b8*b6*a5+21546*b6*a8*a3-14490*a8^2*a3+11690*b5*a4*a9+14280*b6*b5*a5-12600*b6*a5*a7-7770*b7*a6*a5-3087*b8*a8*a9+3871*b8*a8^2-3339*a6*b7*a8-336*a4*a7*a8+28*b8*a8*a4-3374*b5*a8^2+8988*a8^2*a7-1582*b8*a10*a3-2940*a7*a9^2-11088*a7*a4^2-5775*b5*a4^2+3318*b4*a4^2+15911*b4*a9^2-6321*a8*b9*a9-13209*a8*b4*a9-378*a10*a3^2+4956*b4*a10*a7-8064*a8*a9*a7-12537*a8*b5*a9+24570*a8*a9*a3-4536*a7*a5*a4-5376*a7^2*a6-6552*b9*a10*a7+2604*a10*a3*a7+3024*a8*b7*a10-18396*a9*a4*a3+8393*a9^2*b5+84*a9*b9*a4+9912*a9*a4*a7+1848*a9*a5*a7-966*a9^2*b9-10836*a9*a5*a3-1008*a7^2*a10-18396*b6*b9*a8+9702*a4^2*a3-4284*a3^2*a6-5236*b4*a10*a3+2919*b10*b9*a4-3458*b4^2*a10+5166*a3*a5*a4-1890*b4*b6*a4-1008*b7*b6*a10-5544*b7*b6*a6+504*b7*b10*a10+4536*a5^2*a3-3255*b8*b6*a9+2905*b8*b10*a9-8617*a4*b4*a9-8498*b8*b5*a10-2520*b8*b9*a10+504*b7*b10*a6-19488*b9*a4^2+1414*b10*b5*a5-15372*b5^2*a6+4746*b4*a6*a3-3136*b8^2*a6-4676*b8*a6*a7-2352*b10*a9*a7+1106*b4*b10*a8+1176*b10*a4*a7+1344*b10*a8*a7+9583*b10*b5*a9-9744*b9*a6*a7+21042*b4*b9*a10+378*b9*a6*a3-3766*b9*b8*a6-1568*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6+14*b4*b5*a10+9569*b8*a5*a9-4501*b10*b5*a8-10003*b4*a5*a9-4025*b10*b5*a4+13230*b9*a10*a3-7420*b4*b9*a6-3129*b8*b6*a4-378*b10*a9*a3-504*b10*a5*a3-8015*b8*b6*a8+7434*b9*b10*a5+2380*b4*b6*a8-4172*b4*a8*a4-4914*b6*a9*a3+7056*b6*a9*a7-630*b10*a8*a3-5138*b8*b10*a5-7770*b4*b5*a6+1512*a9^2*a3+1050*b4^2*a6+126*b10*a4*a3+24591*a5*b9*a9-9044*b4*a6*a7+3150*a5^2*b4+5733*b7*a6*a9-441*b9*b10*a8-11844*b9*b6*a5+13818*b9*b6*a4+8239*b4*a5*a8-11739*b8*a5*a4+2016*b8^2*a10+17871*b9*a8*a4+20937*b9*a5*a4+266*b8*a9^2-266*b4*b10*a4+6930*b8*a5^2+3171*b8*a4^2+2436*b10*a5*a7-1092*b4*b6*a5+28*b4*a8^2+36750*b5*a6*a3-9044*b5*a6*a7-11830*b5*b9*a6+5635*b5*a8*a4-8638*b5*a10*a3+2940*b5*a10*a7+3906*b5*b9*a10-6888*b5*a5*a4+11858*b5*b8*a6+5852*a5*b5*a9+5229*b9*a8^2-33327*a5*b9*a8-5873*a5*b8*a8+6552*b8*a10*a7+1750*b10*b4*a9-22050*b9*a5^2+16842*b5*b6*a4+7182*b9^2*a6-1512*b9^2*a10-770*b5*b6*a8+2681*b8*b10*a8-6804*a5^2*b5-19446*b5*b6*a9+5544*a5*a3*a8+2016*a5^2*a7"),
		SMQP("6846*a5*b5*a8+189*b4*a5*a4+3990*b8*b4*a10-6804*b9*b6*a9+8064*b10*b4*a5-1512*b7*a10*a5+1211*b10*b9*a9+8568*b8*a6*a3+1575*b8*a9*a4+8988*a6*a3*a7-9072*b7*a10*a4+11781*b7*a6*a4+5040*b7*a10*a9-168*b6*a4*a7-3724*b6*a8*a7-2268*b6*a5*a3-11788*a5*a8*a7+20286*a4*a3*a8+3024*b6*a4*a3+11298*b8*b6*a5+1890*b6*a8*a3-2898*a8^2*a3-1806*b5*a4*a9-9996*b6*b5*a5+840*b6*a5*a7+2646*b7*a6*a5+1701*b8*a8*a9+1995*b8*a8^2+5481*a6*b7*a8-6328*a4*a7*a8-420*b8*a8*a4-798*b5*a8^2+2156*a8^2*a7-6006*b8*a10*a3+3556*a7*a9^2+8568*a7*a4^2+5901*b5*a4^2-15666*b4*a4^2-12285*b4*a9^2+2667*a8*b9*a9-1197*a8*b4*a9+1134*a10*a3^2-4564*b4*a10*a7+2352*a8*a9*a7-1701*a8*b5*a9-6174*a8*a9*a3+4872*a7*a5*a4-5152*a7^2*a6+2016*b9*a10*a7-1988*a10*a3*a7-1008*a8*b7*a10+19404*a9*a4*a3-147*a9^2*b5-644*a9*b9*a4-9296*a9*a4*a7+1288*a9*a5*a7+2674*a9^2*b9+8316*a9*a5*a3+1008*a7^2*a10+3248*b6*b9*a8-26586*a4^2*a3+756*a3^2*a6+17556*b4*a10*a3-2317*b10*b9*a4+15078*b4^2*a10-8946*a3*a5*a4+5334*b4*b6*a4-1512*b7*b6*a10+4284*b7*b6*a6+2520*b7*b10*a10-10584*a5^2*a3-8043*b8*b6*a9+3885*b8*b10*a9+7980*b4*b6*a9+12747*a4*b4*a9-1554*b8*b5*a10+1008*b8*b9*a10-5544*b7*b10*a6+2856*b9*a4^2+3990*b10*b5*a5-756*b5^2*a6-3150*b4*a6*a3-336*b8^2*a6-8036*b8*a6*a7-1120*b10*a9*a7-1638*b4*b10*a8+392*b10*a4*a7+2128*b10*a8*a7-5901*b10*b5*a9-5376*b9*a6*a7-7126*b4*b9*a10+13650*b9*a6*a3-350*b9*b8*a6+504*b5^2*a10-147*b8*b10*a4-11046*b4*b8*a6-1554*b4*b5*a10-1995*b8*a5*a9+1743*b10*b5*a8+35385*b4*a5*a9+3675*b10*b5*a4-4466*b9*a10*a3+9828*b4*b9*a6+3171*b8*b6*a4+1134*b10*a9*a3+1512*b10*a5*a3-3843*b8*b6*a8+1274*b9*b10*a5-12180*b4*b6*a8+12180*b4*a8*a4+6426*b6*a9*a3-1428*b6*a9*a7-4158*b10*a8*a3-5250*b8*b10*a5-6426*b4*b5*a6-7560*a9^2*a3-8190*b4^2*a6-378*b10*a4*a3-2261*a5*b9*a9-4116*b4*a6*a7-23394*a5^2*b4-8631*b7*a6*a9-1085*b9*b10*a8+3528*b9*b6*a5+4158*b9*b6*a4-17829*b4*a5*a8+12201*b8*a5*a4+1008*b8^2*a10-5677*b9*a8*a4-3003*b9*a5*a4-1806*b8*a9^2+3150*b4*b10*a4+7938*b8*a5^2-2289*b8*a4^2+644*b10*a5*a7-6888*b4*b6*a5+1260*b4*a8^2-27594*b5*a6*a3+12348*b5*a6*a7+3402*b5*b9*a6-2961*b5*a8*a4+13482*b5*a10*a3-7084*b5*a10*a7-1582*b5*b9*a10+3696*b5*a5*a4-1638*b5*b8*a6-2604*a5*b5*a9-175*b9*a8^2-259*a5*b9*a8-11277*a5*b8*a8+4032*b8*a10*a7-882*b10*b4*a9-210*b9*a5^2-2814*b5*b6*a4-266*b9^2*a6+42*b5*b6*a8+525*b8*b10*a8-1764*a5^2*b5+8778*b5*b6*a9+12600*a5*a3*a8+2016*a5^2*a7"),
		SMQP("-315*b9^2*a9+63*b9^2*a4+231*b9^2*a8-896*b7*b4*a10+28*b7*b10*a9-364*b7*a10*a3+84*b7*b4*a6+196*b7*b10*a4-378*a9*a3*a7-315*b4*a4*a3-336*b8*a8*a3-273*b8*b9*a4-84*b8*a4*a3+21*b8*a9*a7-504*b7*a6*a3+1281*b4*a8*a3-252*b4*a8*a7-336*b5^2*a5+252*b4*b5*a4+826*a7*b7*a6-1764*a7*a8*a3-567*a7*b5*a9-336*a7*b5*a5+504*b8*b7*a10+546*a7*b4*a5+84*a7*b7*a10-630*a7*b4*a9-1176*a7*b8*a5+315*a7*b8*a8-336*b8*a9*a3+1596*b8*a5*a3+840*b8*b4*a5-231*b8*b5*a4+2394*a7*a4*a3-630*a7^2*a4-336*b4*b8*a8-42*b4*b9*a8+672*b4*a4*a7-84*b4*b5*a8+126*a7^2*a9+147*a7*b8*a4+252*a9*a3^2+777*b7*a4^2-1379*b7*a4*a8-1092*b8^2*a5-630*a3*b4*a5-2268*a3^2*a4-252*a3^2*a5-21*a3*b4*a9-1771*b7*a5*a9-210*b9^2*a5+1267*b7*a5*a8+231*b5*a4*a3+112*b7*b5*a10-189*b5*a4*a7+84*b5^2*a4-168*b5^2*a8+126*a5*a3*a7+336*b5^2*a9+378*a7^2*a8-126*b4*b9*a4+2016*a3^2*a8-210*b5*b9*a5+882*b5*b8*a5-672*b5*b4*a5+420*b4^2*a9+105*b5*b8*a8-112*b7*b8*a6+420*b4^2*a5+420*b7*a9*a8+42*b7*b6*a5+1428*b7*b5*a6-301*b7*a9*a4-861*b7*b6*a4-63*b7*b6*a9-28*b7*b10*a8+231*b5*b9*a8+679*b7*b6*a8-112*b7*a9^2+966*b7*a5^2-308*b7*b10*a5-315*b5*b9*a9+196*b7*a8^2-903*a5*b7*a4+42*b8^2*a9+126*b8^2*a8+63*b9*b5*a4+1113*b5*a9*a3-357*b5*b8*a9-84*b5*b4*a9+294*b8^2*a4-63*b9*a9*a3-966*b9*a5*a3+609*b9*b8*a9-525*b9*b8*a8-588*b9*b4*a5+966*b9*b8*a5+42*b9*b4*a9-1281*b5*a8*a3+609*b5*a8*a7+1218*b5*a5*a3+588*b9*a5*a7-21*b9*a8*a3+315*b9*a4*a3+189*b9*a9*a7-189*b9*a4*a7-105*b9*a8*a7-168*b9*b7*a10-700*b9*b7*a6"),
		SMQP("38787*a6*a3*a4+6237*a6*a3*a5+118524*a6*b5*b6+2079*a6*b5*a5+5376*a6*a3*b6+27300*a6*b4*b6+37653*a6*b4*a5+85218*a6*b5*a4-2016*a4*b6*a5+4704*a6*b4*a4+12432*a4^2*a5-9408*a4*b6^2+43848*a4*a5^2+14784*b6^2*a5+29232*b6*a5^2+23688*a4^3+12432*a4^2*b6+26208*a5^3-29568*b10*a5^2-107646*b10*b5*a6+12432*b8*a10*a5+20272*b10*b6*a9+38024*b10*b5*a10-16744*b10*b9*a10-1568*b10*a6*a7+15456*b10*a6*a3-2184*b6*a8*a4+5152*b6*b10*a4+10325*a6*a8*a7+18592*a9*b5*a10-6384*b6*b10*a8+15995*a6*a9*a7-32725*a6*a4*a7-5040*b10*a10*a7+15456*b6*a10*a7-17248*b6*a10*a3-20321*a6*b9*a4-38136*a6*a8*a3-4536*b10*a9*a8+10416*b10*a9*a4-8624*b6*a6*a7-3752*b10*a9^2+46032*a6*b7*a10+21049*a6*b9*a9+46837*a6*b9*a8+9464*b6*a9*a8+24136*b10*a8*a4+35966*b10*b9*a6-39088*b5*b6*a10-72891*a6*b4*a9+2576*a8^3+11032*a6*b8*a8+19880*b6*a9^2-17024*b9*a10*a8-18480*b6*a8^2-13160*a8*a9^2-53368*a8*a4^2+18984*a8*b4*a10+12040*a8*b5*a10+22960*a9^2*a5-6888*a8*a10*a7+27832*a8^2*a4+21616*a8*a10*a3-10731*a6*a9*a3+1456*a9^3-12152*a9*a4^2-952*a9*a8^2+17472*a9*a8*a4+12880*a9*a10*a3+42840*a9*b4*a10-13272*a10*a9*a7+15960*a10*a4*a7-56336*a10*b5*a4+8344*a10*b9*a4-15680*a10*a4*a3+728*a10*b8*a8-10080*a10^2*b7-20216*a10*b9*a9-6944*b6*b10*a5+15022*b8*a6*a4-32200*a9*b6*a4+9492*a8*b4*a6-1512*a5*b10*a8-14392*b10*a4^2-6832*b10*a8^2+1792*b10^2*a5-44632*a5*b6*a9-952*a4*b8*a10-29064*a4*b4*a10+13664*b9*a10*a5+7280*a9^2*a4-1960*a5*b10*a4+112*b10^2*a4+14896*b6^2*a8-23520*b6^2*a9+6272*b8*a10*a9-32102*b8*a6*a9+32200*a5*b10*a9-5264*b10^2*a9+15120*a5*a10*a7-2254*b8*b10*a6-30576*a5^2*a8-28371*a5*b9*a6-71932*b9*b6*a6-13048*b4*a10*a5-21784*a5*a8^2-29274*a6*b5*a9-4305*a6*b5*a8-54117*a6^2*b7+4872*b4*b10*a6+224*b10^2*a8+1904*b10*a10*a3-14112*a5*b5*a10-14784*a5*a10*a3-40712*a5*a9*a4-27790*b8*b6*a6+20664*a5*a9*a8-8456*b8*b10*a10-39291*b8*a6*a5-35280*b4*b6*a10+23464*b8*b6*a10+14448*b4*b10*a10+13496*a5*a8*a4+2296*a5*b6*a8-30912*a5^2*a9+40880*b9*b6*a10-21000*a5*a6*a7"),
		SMQP("-315*b9^2*a9-777*b9^2*a4-133*b9^2*a8-252*b7*b4*a10+280*b7*b10*a9-1344*b7*a10*a3-217*b7*b4*a6-1701*a9*a3*a7-714*b4*a4*a3-994*b8*a8*a3+959*b8*b9*a4-196*b8*a4*a3-35*b8*a9*a7+1379*b7*a6*a3+805*b4*a8*a3+315*b4*a8*a7-280*b5^2*a5-2520*b4*b5*a4-462*a7*b7*a6-3297*a7*a8*a3-903*a7*b5*a9+224*a7*b5*a5+168*b8*b7*a10+1106*a7*b4*a5+700*a7*b7*a10-2037*a7*b4*a9-1232*a7*b8*a5+567*a7*b8*a8-238*b8*a9*a3+2485*b8*a5*a3+42*b8*b4*a9-35*b8*b4*a5+889*b8*b5*a4+5103*a7*a4*a3-1470*a7^2*a4+14*b4*b8*a8-329*b4*b9*a8+161*b4*a4*a7-595*b4*b5*a8+294*a7^2*a9-77*a7*b8*a4+1071*a9*a3^2+1141*b7*a4^2-707*b7*a4*a8+476*b8^2*a5-2394*a3*b4*a5-5271*a3^2*a4-609*a3^2*a5+2639*a3*b4*a9-763*b7*a5*a9+518*b9^2*a5-385*b7*a5*a8+1715*b5*a4*a3+84*b7*b5*a10-1029*b5*a4*a7+196*b5^2*a4+112*b5^2*a8+798*a5*a3*a7+280*b5^2*a9+882*a7^2*a8+1827*b4*b9*a4+4200*a3^2*a8+518*b5*b9*a5+42*b5*b8*a5+1057*b5*b4*a5+245*b4^2*a9-91*b5*b8*a8+224*b7*b8*a6+77*b4^2*a5+336*b7*a9*a8+882*b7*b6*a5-280*b7*b5*a6-1029*b7*a9*a4-637*b7*b6*a4-707*b7*b6*a9+224*b7*b10*a8-133*b5*b9*a8-161*b7*b6*a8+168*b7*a9^2+910*b7*a5^2-420*b7*b10*a5-315*b5*b9*a9+168*b7*a8^2-147*a5*b7*a4-70*b8^2*a9-42*b8^2*a8+231*b9*b5*a4+2513*b5*a9*a3-245*b5*b8*a9+476*b5*b4*a9-826*b8^2*a4-252*b9*a9*a3-553*b9*a5*a3+665*b9*b8*a9+287*b9*b8*a8+77*b9*b4*a5-1274*b9*b8*a5-511*b9*b4*a9-1456*b5*a8*a3+749*b5*a8*a7-329*b5*a5*a3+364*b9*a5*a7-868*b9*a8*a3+1008*b9*a4*a3+21*b9*a9*a7-189*b9*a4*a7+427*b9*a8*a7+112*b9*b7*a10-392*b9*b7*a6"),
		SMQP("4389*a6*a3*a4+5523*a6*a3*a5+24024*a6*b5*b6+5229*a6*b5*a5-1848*a6*a3*b6+8064*a6*b4*b6+3339*a6*b4*a5+18690*a6*b5*a4-4788*a4*b6*a5+2688*a6*b4*a4+1092*a4^2*a5-4704*a4*b6^2+8652*a4*a5^2+2856*b6^2*a5+10920*b6*a5^2+5544*a4^3+5628*a4^2*b6+5040*a5^3-2632*b10*a5^2-23478*b10*b5*a6+7560*b8*a10*a5+3164*b10*b6*a9+8176*b10*b5*a10-4032*b10*b9*a10+448*b10*a6*a7+2352*b10*a6*a3-3836*b6*a8*a4-868*b6*b10*a4+6587*a6*a8*a7+3444*a9*b5*a10+308*b6*b10*a8+357*a6*a9*a7-5971*a6*a4*a7-1008*b10*a10*a7+3024*b6*a10*a7-5852*b6*a10*a3-3311*a6*b9*a4-6384*a6*a8*a3-1288*b10*a9*a8+1708*b10*a9*a4+1736*b6*a6*a7-1708*b10*a9^2+13608*a6*b7*a10+987*a6*b9*a9+5803*a6*b9*a8+2604*b6*a9*a8+3220*b10*a8*a4+6566*b10*b9*a6-1204*b5*b6*a10-13461*a6*b4*a9+476*a8^3+532*a6*b8*a8+4732*b6*a9^2-2520*b9*a10*a8-3416*b6*a8^2-1316*a8*a9^2-8092*a8*a4^2+812*a8*b4*a10+1316*a8*b5*a10+3640*a9^2*a5-2520*a8*a10*a7+2212*a8^2*a4+6244*a8*a10*a3-1701*a6*a9*a3-1176*a9^3-2492*a9*a4^2+3220*a9*a8*a4+84*a9*a10*a3+7980*a9*b4*a10-1512*a10*a9*a7+3528*a10*a4*a7-10360*a10*b5*a4-504*a10*b9*a4-5264*a10*a4*a3-3024*a10^2*b7-2520*a10*b9*a9-2632*b6*b10*a5+1750*b8*a6*a4-4256*a9*b6*a4+1176*a8*b4*a6+2380*a5*b10*a8-2632*b10*a4^2-28*b10*a8^2-560*b10^2*a5-16604*a5*b6*a9+1008*a4*b8*a10-5320*a4*b4*a10-1512*b9*a10*a5+1624*a9^2*a4-2324*a5*b10*a4+448*b10^2*a4+6188*b6^2*a8-4452*b6^2*a9+630*b8*a6*a9+3220*a5*b10*a9-560*b10^2*a9+1008*a5*a10*a7+938*b8*b10*a6-6496*a5^2*a8-7637*a5*b9*a6-17864*b9*b6*a6+476*b4*a10*a5-980*a5*a8^2-9030*a6*b5*a9+8421*a6*b5*a8-13167*a6^2*b7-2184*b4*b10*a6-448*b10^2*a8-112*b10*a10*a3-2044*a5*b5*a10+28*a5*a10*a3-15316*a5*a9*a4+322*b8*b6*a6+5740*a5*a9*a8-1008*b8*b10*a10-11585*b8*a6*a5-11284*b4*b6*a10+3024*b8*b6*a10+2128*b4*b10*a10+4648*a5*a8*a4+9604*a5*b6*a8-5264*a5^2*a9+9072*b9*b6*a10-5152*a5*a6*a7"),
		SMQP("90258*a5*b5*a8+18921*b4*a5*a4+23282*b8*b4*a10+6132*b9*b6*a9-7504*b10*b4*a5+14448*b7*a10*a5-5215*b10*b9*a9+30744*b8*a6*a3+7945*b8*a9*a4-18732*a6*a3*a7+7560*b7*a10*a4-13923*b7*a6*a4-3024*b7*a10*a9-2520*b6*a4*a7-24612*b6*a8*a7+15372*b6*a5*a3-3444*a5*a8*a7-62958*a4*a3*a8+7056*b6*a4*a3+39102*b8*b6*a5+5334*b6*a8*a3+24234*a8^2*a3-27930*b5*a4*a9-9492*b6*b5*a5-11592*b6*a5*a7-25830*b7*a6*a5+3423*b8*a8*a9+2849*b8*a8^2+46179*a6*b7*a8+29694*a4*a7*a8+26096*b8*a8*a4-25578*b5*a8^2+10668*a8^2*a7-3122*b8*a10*a3+1932*a7*a9^2-29862*a7*a4^2+27447*b5*a4^2-8358*b4*a4^2-3815*b4*a9^2-5271*a8*b9*a9+11361*a8*b4*a9-10374*a10*a3^2+17220*b4*a10*a7-26712*a8*a9*a7+27657*a8*b5*a9+9198*a8*a9*a3-19656*a7*a5*a4-672*a7^2*a6-7560*b9*a10*a7+23604*a10*a3*a7-17136*a8*b7*a10+3234*a9*a4*a3+13671*a9^2*b5+994*a9*b9*a4+10626*a9*a4*a7+30408*a9*a5*a7+2422*a9^2*b9-27132*a9*a5*a3-7056*a7^2*a10+7616*b6*b9*a8+30996*a4^2*a3+22428*a3^2*a6-4844*b4*a10*a3+1001*b10*b9*a4+8162*b4^2*a10+38808*a3*a5*a4+1890*b4*b6*a4+3024*b7*b6*a10-10332*b7*b6*a6-504*b7*b10*a10+13608*a5^2*a3-8673*b8*b6*a9+6503*b8*b10*a9-5964*b4*b6*a9-13181*a4*b4*a9+13202*b8*b5*a10-3528*b8*b9*a10-1512*b7*b10*a6-14406*b9*a4^2-5334*b10*b5*a5+40572*b5^2*a6+28518*b4*a6*a3-4592*b8^2*a6-30940*b8*a6*a7-7728*b10*a9*a7+1246*b4*b10*a8+8904*b10*a4*a7+672*b10*a8*a7+5985*b10*b5*a9+4592*b9*a6*a7+24542*b4*b9*a10-14490*b9*a6*a3+26726*b9*b8*a6-14112*b5^2*a10+3815*b8*b10*a4-20258*b4*b8*a6-57358*b4*b5*a10+24703*b8*a5*a9+7749*b10*b5*a8+4963*b4*a5*a9-6279*b10*b5*a4+7714*b9*a10*a3-29204*b4*b9*a6+9177*b8*b6*a4+5082*b10*a9*a3-4872*b10*a5*a3-19033*b8*b6*a8+6230*b9*b10*a5-11788*b4*b6*a8+2828*b4*a8*a4-9954*b6*a9*a3+23940*b6*a9*a7-1050*b10*a8*a3-10990*b8*b10*a5+20538*b4*b5*a6-7224*a9^2*a3-6090*b4^2*a6-5502*b10*a4*a3+49*a5*b9*a9-25564*b4*a6*a7-19614*a5^2*b4-693*b7*a6*a9-5495*b9*b10*a8-504*b9*b6*a5+2394*b9*b6*a4+3353*b4*a5*a8-7875*b8*a5*a4-6048*b8^2*a10+17507*b9*a8*a4+6321*b9*a5*a4-9674*b8*a9^2+1946*b4*b10*a4-25074*b8*a5^2-20307*b8*a4^2+10668*b10*a5*a7+6552*b4*b6*a5-2380*b4*a8^2-126*b5*a6*a3+12852*b5*a6*a7-34650*b5*b9*a6-11865*b5*a8*a4-1218*b5*a10*a3-924*b5*a10*a7+15470*b5*b9*a10-38766*b5*a5*a4-22050*b5*b8*a6-63084*a5*b5*a9-8869*b9*a8^2-13657*a5*b9*a8-20503*a5*b8*a8+7560*b8*a10*a7+6314*b10*b4*a9+10626*b9*a5^2-7098*b5*b6*a4+322*b9^2*a6-4536*b9^2*a10+24822*b5*b6*a8-2345*b8*b10*a8+12852*a5^2*b5-23562*b5*b6*a9-29568*a5*a3*a8-6048*a5^2*a7"),
		SMQP("1260*b9^2*a9+2016*b9^2*a4+2268*b9^2*a8+1120*b7*b4*a10-5096*b7*b10*a9+12488*b7*a10*a3-1029*b7*b4*a6+952*b7*b10*a4+25956*a9*a3*a7+9513*b4*a4*a3+18900*b8*a8*a3-2331*b8*b9*a4-2583*b8*a4*a3+2709*b8*a9*a7-18942*b7*a6*a3-8316*b4*a8*a3-1701*b4*a8*a7+6048*b5^2*a5+12600*b4*b5*a4+11662*a7*b7*a6+46746*a7*a8*a3+17388*a7*b5*a9-630*a7*b5*a5-5040*b8*b7*a10-12222*a7*b4*a5-9072*a7*b7*a10+23751*a7*b4*a9+8946*a7*b8*a5-8253*a7*b8*a8-441*b8*a9*a3-34839*b8*a5*a3+1071*b8*b4*a9-1890*b8*b4*a5-2520*b8*b5*a4+5040*a7^2*a5-80640*a7*a4*a3+22050*a7^2*a4+630*b4*b8*a8-819*b4*b9*a8-6363*b4*a4*a7+9765*b4*b5*a8-9198*a7^2*a9+1701*a7*b8*a4-14742*a9*a3^2-7056*b7*a4^2+9436*b7*a4*a8+2457*b8^2*a5+33705*a3*b4*a5+79758*a3^2*a4+10962*a3^2*a5-37611*a3*b4*a9+140*b7*a5*a9-6552*b9^2*a5+17164*b7*a5*a8-12852*b5*a4*a3+8680*b7*b5*a10+4032*b5*a4*a7+756*b5^2*a4-504*b5^2*a8-14490*a5*a3*a7-4032*b5^2*a9-8946*a7^2*a8-12411*b4*b9*a4-62496*a3^2*a8-6552*b5*b9*a5-1323*b5*b8*a5-2835*b5*b4*a5-4347*b4^2*a9-4599*b5*b8*a8+6083*b7*b8*a6+693*b4^2*a5-2856*b7*a9*a8-9912*b7*b6*a5-6552*b7*b5*a6+11732*b7*a9*a4-1848*b7*b6*a4+9492*b7*b6*a9-952*b7*b10*a8-756*b5*b9*a8+3556*b7*b6*a8-6160*b7*a9^2-8568*b7*a5^2+3472*b7*b10*a5+7308*b5*b9*a9-3080*b7*a8^2-14616*a5*b7*a4-3150*b8^2*a9+1134*b8^2*a8-2520*b9*b5*a4-35028*b5*a9*a3+6804*b5*b8*a9-13104*b5*b4*a9+4032*b8^2*a4+4662*b9*a9*a3+12222*b9*a5*a3-5481*b9*b8*a9-2079*b9*b8*a8-4851*b9*b4*a5+8757*b9*b8*a5+12915*b9*b4*a9+24570*b5*a8*a3-10962*b5*a8*a7+126*b5*a5*a3-1638*b9*a5*a7+10962*b9*a8*a3-13986*b9*a4*a3-4158*b9*a9*a7+5922*b9*a4*a7-5418*b9*a8*a7-1008*b9*b7*a10+1736*b9*b7*a6"),
		SMQP("315*b9^2*a9-63*b9^2*a4-147*b9^2*a8-756*b7*b4*a10-168*b7*b10*a9-252*b7*b4*a6+378*a9*a3*a7+567*b4*a4*a3+504*b8*a8*a3+105*b8*b9*a4-84*b8*a4*a3+147*b8*a9*a7-1008*b7*a6*a3-693*b4*a8*a3+168*b5^2*a5-252*b4*b5*a4+630*a7*b7*a6+1764*a7*a8*a3+567*a7*b5*a9-84*a7*b5*a5-168*b8*b7*a10-210*a7*b4*a5-252*a7*b7*a10+882*a7*b4*a9+588*a7*b8*a5-399*a7*b8*a8+168*b8*a9*a3-1428*b8*a5*a3-420*b8*b4*a5-105*b8*b5*a4-2394*a7*a4*a3+630*a7^2*a4+168*b4*b8*a8+126*b4*b9*a8-588*b4*a4*a7-126*a7^2*a9+21*a7*b8*a4-252*a9*a3^2+483*b7*a4^2-357*b7*a4*a8-84*b8^2*a5+630*a3*b4*a5+2268*a3^2*a4+252*a3^2*a5-903*a3*b4*a9-1701*b7*a5*a9-210*b9^2*a5+2793*b7*a5*a8-399*b5*a4*a3+588*b7*b5*a10+189*b5*a4*a7+84*b5^2*a4-126*a5*a3*a7-168*b5^2*a9-378*a7^2*a8-126*b4*b9*a4-2016*a3^2*a8-210*b5*b9*a5-126*b5*b8*a5-420*b5*b4*a5-336*b4^2*a9+147*b5*b8*a8+1344*b7*b8*a6+168*b4^2*a5-126*b7*b6*a5+504*b7*b5*a6+357*b7*a9*a4-651*b7*b6*a4+147*b7*b6*a9-147*b5*b9*a8+777*b7*b6*a8-168*b7*a9^2+966*b7*a5^2-252*b7*b10*a5+315*b5*b9*a9-168*b7*a8^2-1897*a5*b7*a4+294*b8^2*a9-294*b8^2*a8-63*b9*b5*a4-1281*b5*a9*a3+21*b5*b8*a9+420*b5*b4*a9+42*b8^2*a4+63*b9*a9*a3+546*b9*a5*a3-777*b9*b8*a9+441*b9*b8*a8+168*b9*b4*a5+462*b9*b8*a5+42*b9*b4*a9+1197*b5*a8*a3-525*b5*a8*a7+210*b5*a5*a3-168*b9*a5*a7+105*b9*a8*a3-315*b9*a4*a3-189*b9*a9*a7+189*b9*a4*a7+21*b9*a8*a7-840*b9*b7*a6"),
		SMQP("-3360*a6*a3*a4-3864*a6*a3*a5-21084*a6*b5*b6-8568*a6*b5*a5+7728*a6*a3*b6-1764*a6*b4*b6-5670*a6*b4*a5-462*a6*b5*a4+5040*a4*b6*a5-840*a6*b4*a4+4872*a4^2*a5+6888*a4*b6^2-11256*a4*a5^2+840*b6^2*a5-12768*b6*a5^2-3528*a4^3-7224*a4^2*b6+2016*a5^3-3556*b10*a5^2+1344*b10*b5*a6-5859*b8*a10*a5-1582*b10*b6*a9-4487*b10*b5*a10+1953*b10*b9*a10+196*b10*a6*a7-1428*b10*a6*a3+6160*b6*a8*a4+1568*b6*b10*a4-6412*a6*a8*a7-2121*a9*b5*a10-1246*b6*b10*a8+1596*a6*a9*a7+896*a6*a4*a7+504*b10*a10*a7-3024*b6*a10*a7+5824*b6*a10*a3+5572*a6*b9*a4+8484*a6*a8*a3+350*b10*a9*a8-518*b10*a9*a4-5656*b6*a6*a7+1106*b10*a9^2-5103*a6*b7*a10-882*a6*b9*a9+28*a6*b9*a8-2520*b6*a9*a8-644*b10*a8*a4-616*b10*b9*a6+6650*b5*b6*a10+3108*a6*b4*a9-238*a8^3-3290*a6*b8*a8-1106*b6*a9^2+1323*b9*a10*a8+1834*b6*a8^2+154*a8*a9^2+4634*a8*a4^2-532*a8*b4*a10-532*a8*b5*a10+112*a9^2*a5+1260*a8*a10*a7-1022*a8^2*a4-3122*a8*a10*a3-1764*a6*a9*a3+924*a9^3+1414*a9*a4^2+168*a9*a8^2-686*a9*a8*a4+294*a9*a10*a3-1239*a9*b4*a10+756*a10*a9*a7-252*a10*a4*a7-3367*a10*b5*a4-1386*a10*b9*a4+2044*a10*a4*a3-63*a10*b8*a8+1512*a10^2*b7+1134*a10*b9*a9+2828*b6*b10*a5-2534*b8*a6*a4+1624*a9*b6*a4+4326*a8*b4*a6-4844*a5*b10*a8+266*b10*a4^2+56*b10*a8^2+700*b10^2*a5+17752*a5*b6*a9-567*a4*b8*a10+3500*a4*b4*a10+1827*b9*a10*a5-1820*a9^2*a4+2674*a5*b10*a4-182*b10^2*a4-7504*b6^2*a8+1848*b6^2*a9+126*b8*a10*a9-2688*b8*a6*a9+3220*a5*b10*a9+322*b10^2*a9+1008*a5*a10*a7-4606*b8*b10*a6+13244*a5^2*a8-8918*a5*b9*a6+9940*b9*b6*a6-2653*b4*a10*a5-1736*a5*a8^2+168*a6*b5*a9+9870*a6*b5*a8+4410*a6^2*b7+210*b4*b10*a6+182*b10^2*a8+518*b10*a10*a3+1694*a5*b5*a10-4088*a5*a10*a3+8162*a5*a9*a4+1414*b8*b6*a6-4592*a5*a9*a8+567*b8*b10*a10+21406*b8*a6*a5+10430*b4*b6*a10-1449*b8*b6*a10-266*b4*b10*a10-2450*a5*a8*a4-14336*a5*b6*a8-5600*a5^2*a9-2898*b9*b6*a10+3752*a5*a6*a7"),
		SMQP("-76314*a5*b5*a8-36141*b4*a5*a4-20426*b8*b4*a10-19068*b9*b6*a9+560*b10*b4*a5-29232*b7*a10*a5+5719*b10*b9*a9-29512*b8*a6*a3-2737*b8*a9*a4+6076*a6*a3*a7+6552*b7*a10*a4-4977*b7*a6*a4+5040*b7*a10*a9-11592*b6*a4*a7+60228*b6*a8*a7-15708*b6*a5*a3+86436*a5*a8*a7+67102*a4*a3*a8-4368*b6*a4*a3-39942*b8*b6*a5-11102*b6*a8*a3-12194*a8^2*a3-2310*b5*a4*a9+7140*b6*b5*a5+15624*b6*a5*a7+46494*b7*a6*a5+7497*b8*a8*a9-4865*b8*a8^2-57267*a6*b7*a8-67914*a4*a7*a8-34328*b8*a8*a4+17682*b5*a8^2-8820*a8^2*a7-3766*b8*a10*a3+2772*a7*a9^2+44730*a7*a4^2-9975*b5*a4^2+19614*b4*a4^2-15449*b4*a9^2+3927*a8*b9*a9+567*a8*b4*a9+10766*a10*a3^2-50148*b4*a10*a7+20160*a8*a9*a7+2415*a8*b5*a9-33726*a8*a9*a3-39816*a7*a5*a4+16128*a7^2*a6+19656*b9*a10*a7-33012*a10*a3*a7+15120*a8*b7*a10+14042*a9*a4*a3-20055*a9^2*b5+9002*a9*b9*a4+14994*a9*a4*a7-90216*a9*a5*a7+602*a9^2*b9+42140*a9*a5*a3+7056*a7^2*a10-3080*b6*b9*a8-39816*a4^2*a3-20076*a3^2*a6+9548*b4*a10*a3+4207*b10*b9*a4-17290*b4^2*a10-26292*a3*a5*a4-7938*b4*b6*a4-17136*b7*b6*a10+16380*b7*b6*a6+504*b7*b10*a10-19656*a5^2*a3+6489*b8*b6*a9-2303*b8*b10*a9+5796*b4*b6*a9+6433*a4*b4*a9+1750*b8*b5*a10+31752*b8*b9*a10+8568*b7*b10*a6+11298*b9*a4^2+15582*b10*b5*a5-52668*b5^2*a6-52206*b4*a6*a3+23408*b8^2*a6+79324*b8*a6*a7+9072*b10*a9*a7-4886*b4*b10*a8-7560*b10*a4*a7-2016*b10*a8*a7-11697*b10*b5*a9-26768*b9*a6*a7-16646*b4*b9*a10+18242*b9*a6*a3-51086*b9*b8*a6+18144*b5^2*a10-6335*b8*b10*a4+36778*b4*b8*a6+22022*b4*b5*a10-56959*b8*a5*a9-7917*b10*b5*a8-20699*b4*a5*a9+1911*b10*b5*a4-15610*b9*a10*a3-22316*b4*b9*a6-9513*b8*b6*a4-4130*b10*a9*a3+8344*b10*a5*a3+29785*b8*b6*a8-16142*b9*b10*a5+23156*b4*b6*a8-6916*b4*a8*a4+13146*b6*a9*a3-44100*b6*a9*a7+5138*b10*a8*a3+10486*b8*b10*a5+40446*b4*b5*a6+3080*a9^2*a3+5250*b4^2*a6+2422*b10*a4*a3-4081*a5*b9*a9+39452*b4*a6*a7+40614*a5^2*b4+5733*b7*a6*a9+9863*b9*b10*a8+10920*b9*b6*a5-3066*b9*b6*a4+26327*b4*a5*a8+24087*b8*a5*a4-10080*b8^2*a10-20447*b9*a8*a4-21189*b9*a5*a4+3962*b8*a9^2-4354*b4*b10*a4+28602*b8*a5^2+21147*b8*a4^2-26460*b10*a5*a7-8568*b4*b6*a5+4004*b4*a8^2-26586*b5*a6*a3+4284*b5*a6*a7+57330*b5*b9*a6+10101*b5*a8*a4+27034*b5*a10*a3+252*b5*a10*a7-25718*b5*b9*a10+23226*b5*a5*a4+20538*b5*b8*a6+55860*a5*b5*a9+3997*b9*a8^2+12817*a5*b9*a8+38815*a5*b8*a8-19656*b8*a10*a7+3710*b10*b4*a9+7350*b9*a5^2+3738*b5*b6*a4+5894*b9^2*a6-7560*b9^2*a10-32550*b5*b6*a8-6727*b8*b10*a8-7812*a5^2*b5+34986*b5*b6*a9+13720*a5*a3*a8+42336*a5^2*a7"),
		SMQP("53991*a5*b5*a8+4347*b4*a5*a4-2450*b8*b4*a10+16359*b9*b6*a9-3563*b10*b4*a5+8988*b7*a10*a5-7462*b10*b9*a9-924*b8*a6*a3-721*b8*a9*a4-26460*a6*a3*a7+14364*b7*a10*a4-20097*b7*a6*a4-13104*b7*a10*a9-3843*b6*a4*a7-8337*b6*a8*a7+15183*b6*a5*a3+24864*a5*a8*a7-34461*a4*a3*a8+3213*b6*a4*a3+13083*b8*b6*a5-1134*b6*a8*a3+5166*a8^2*a3-20370*b5*a4*a9+14049*b6*b5*a5-12600*b6*a5*a7-19110*b7*a6*a5-567*b8*a8*a9+469*b8*a8^2+15624*a6*b7*a8+7035*a4*a7*a8+1288*b8*a8*a4-11025*b5*a8^2+6279*a8^2*a7-1834*b8*a10*a3-9555*a7*a9^2-17766*a7*a4^2+19929*b5*a4^2+3318*b4*a4^2+15344*b4*a9^2-8253*a8*b9*a9-3444*a8*b4*a9-7938*a10*a3^2+5964*b4*a10*a7-11844*a8*a9*a7+3864*a8*b5*a9+24381*a8*a9*a3-38052*a7*a5*a4+8736*a7^2*a6-7308*b9*a10*a7+15204*a10*a3*a7-5040*a8*b7*a10-20853*a9*a4*a3+18921*a9^2*b5+6979*a9*b9*a4+22827*a9*a4*a7+5124*a9*a5*a7-6671*a9^2*b9-15561*a9*a5*a3-6048*a7^2*a10-3787*b6*b9*a8+26208*a4^2*a3+18396*a3^2*a6-8260*b4*a10*a3+3500*b10*b9*a4-3458*b4^2*a10+51912*a3*a5*a4-1890*b4*b6*a4-756*b7*b6*a10-11025*b7*b6*a6-1512*b7*b10*a10+9324*a5^2*a3+1659*b8*b6*a9+259*b8*b10*a9-8127*b4*b6*a9-15295*a4*b4*a9+8134*b8*b5*a10+2268*b8*b9*a10+5355*b7*b10*a6-14070*b9*a4^2-8253*b10*b5*a5+26712*b5^2*a6+17598*b4*a6*a3+5936*b8^2*a6-2156*b8*a6*a7-5691*b10*a9*a7+1106*b4*b10*a8+4893*b10*a4*a7-357*b10*a8*a7+15855*b10*b5*a9+5768*b9*a6*a7+27818*b4*b9*a10-18102*b9*a6*a3+6258*b9*b8*a6-6804*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-43526*b4*b5*a10+5222*b8*a5*a9-3528*b10*b5*a8-22918*b4*a5*a9-5061*b10*b5*a4+9226*b9*a10*a3-28252*b4*b9*a6-3129*b8*b6*a4+945*b10*a9*a3-2205*b10*a5*a3-3857*b8*b6*a8-2569*b9*b10*a5+2380*b4*b6*a8-4172*b4*a8*a4-8883*b6*a9*a3+17073*b6*a9*a7+6930*b10*a8*a3-1673*b8*b10*a5+16506*b4*b5*a6+3087*a9^2*a3+1050*b4^2*a6-5733*b10*a4*a3-2660*a5*b9*a9-12068*b4*a6*a7+18018*a5^2*b4+13356*b7*a6*a9-1652*b9*b10*a8-987*b9*b6*a5-21*b9*b6*a4+6916*b4*a5*a8-14889*b8*a5*a4-6552*b8^2*a10+19166*b9*a8*a4+2247*b9*a5*a4+1148*b8*a9^2-266*b4*b10*a4-24570*b8*a5^2+3171*b8*a4^2+2436*b10*a5*a7-147*b4*b6*a5+28*b4*a8^2+25578*b5*a6*a3-756*b5*a6*a7-35910*b5*b9*a6-10815*b5*a8*a4-13230*b5*a10*a3+7224*b5*a10*a7+17486*b5*b9*a10-31374*b5*a5*a4+3150*b5*b8*a6-38115*a5*b5*a9-1246*b9*a8^2-2590*a5*b9*a8+15232*a5*b8*a8+4284*b8*a10*a7+4459*b10*b4*a9+8946*b9*a5^2-2226*b5*b6*a4+5110*b9^2*a6-7812*b9^2*a10+22743*b5*b6*a8+1295*b8*b10*a8+18144*a5^2*b5-35238*b5*b6*a9-44541*a5*a3*a8+4536*a5^2*a7"),
		SMQP("1260*b9^2*a9+2268*b9^2*a4+2016*b9^2*a8+3640*b7*b4*a10-4088*b7*b10*a9+15008*b7*a10*a3+168*b7*b4*a6+1456*b7*b10*a4+31878*a9*a3*a7+11340*b4*a4*a3+16002*b8*a8*a3-756*b8*b9*a4+504*b8*a4*a3+2646*b8*a9*a7-19383*b7*a6*a3-12348*b4*a8*a3-2016*b4*a8*a7+6804*b5^2*a5+9828*b4*b5*a4+8071*a7*b7*a6+43785*a7*a8*a3+14112*a7*b5*a9-2835*a7*b5*a5-4536*b8*b7*a10-13167*a7*b4*a5-9072*a7*b7*a10+25641*a7*b4*a9+13545*a7*b8*a5-8442*a7*b8*a8+1134*b8*a9*a3-39501*b8*a5*a3-252*b8*b4*a9-2772*b8*b4*a5-4284*b8*b5*a4+2520*a7^2*a5-81144*a7*a4*a3+21609*a7^2*a4+1008*b4*b8*a8-504*b4*b9*a8-7056*b4*a4*a7+12096*b4*b5*a8-9639*a7^2*a9+252*a7*b8*a4-19215*a9*a3^2-10080*b7*a4^2+10444*b7*a4*a8-2268*b8^2*a5+32823*a3*b4*a5+79695*a3^2*a4+11025*a3^2*a5-36477*a3*b4*a9+2912*b7*a5*a9-9072*b9^2*a5+13132*b7*a5*a8-12600*b5*a4*a3+6160*b7*b5*a10+5796*b5*a4*a7+756*b5^2*a4-1008*b5^2*a8-11529*a5*a3*a7-3780*b5^2*a9-8505*a7^2*a8-12096*b4*b9*a4-60480*a3^2*a8-6048*b5*b9*a5-2268*b5*b8*a5-11340*b5*b4*a5-4536*b4^2*a9-2016*b5*b8*a8+2492*b7*b8*a6-5880*b7*a9*a8-11676*b7*b6*a5-5544*b7*b5*a6+12992*b7*a9*a4-336*b7*b6*a4+8232*b7*b6*a9-1960*b7*b10*a8-1512*b5*b9*a8+3052*b7*b6*a8-5152*b7*a9^2-11844*b7*a5^2+4984*b7*b10*a5+8064*b5*b9*a9-1064*b7*a8^2-9828*a5*b7*a4-3024*b8^2*a9+2016*b8^2*a8-3528*b9*b5*a4-33012*b5*a9*a3+5292*b5*b8*a9-13356*b5*b4*a9+3024*b8^2*a4+2709*b9*a9*a3+13671*b9*a5*a3-4788*b9*b8*a9-3528*b9*b8*a8+252*b9*b4*a5+14868*b9*b8*a5+13356*b9*b4*a9+24507*b5*a8*a3-10395*b5*a8*a7+3087*b5*a5*a3-2331*b9*a5*a7+13419*b9*a8*a3-16821*b9*a4*a3-3717*b9*a9*a7+7749*b9*a4*a7-5607*b9*a8*a7-1512*b9*b7*a10+4508*b9*b7*a6"),
		SMQP("4347*b9^2*a9+11781*b9^2*a4+9765*b9^2*a8+7000*b7*b4*a10-13832*b7*b10*a9+54488*b7*a10*a3+8232*b7*b4*a6+3304*b7*b10*a4+128835*a9*a3*a7+34776*b4*a4*a3+47502*b8*a8*a3-7560*b8*b9*a4+9072*b8*a4*a3+15876*b8*a9*a7-77217*b7*a6*a3-37800*b4*a8*a3-9072*b4*a8*a7+30429*b5^2*a5+39312*b4*b5*a4+23842*a7*b7*a6+117999*a7*a8*a3+31437*a7*b5*a9-18522*a7*b5*a5-12096*b8*b7*a10-51282*a7*b4*a5-29232*a7*b7*a10+87822*a7*b4*a9+52038*a7*b8*a5-29988*a7*b8*a8+3906*b8*a9*a3-131355*b8*a5*a3-5040*b8*b4*a9-5040*b8*b4*a5-13608*b8*b5*a4-14112*a7^2*a5-260001*a7*a4*a3+66654*a7^2*a4+2016*b4*b8*a8-10080*b4*b9*a8-19152*b4*a4*a7+51408*b4*b5*a8-32130*a7^2*a9-5544*a7*b8*a4-73521*a9*a3^2-28224*b7*a4^2+15568*b7*a4*a8-20160*b8^2*a5+111825*a3*b4*a5+253449*a3^2*a4+18711*a3^2*a5-128583*a3*b4*a9+3416*b7*a5*a9-43659*b9^2*a5+40936*b7*a5*a8-60039*b5*a4*a3+9016*b7*b5*a10+38997*b5*a4*a7+2772*b5^2*a4-12915*b5^2*a8-630*a5*a3*a7-6048*b5^2*a9-22302*a7^2*a8-42336*b4*b9*a4-182952*a3^2*a8-27342*b5*b9*a5+6489*b5*b8*a5-65835*b5*b4*a5-13104*b4^2*a9-10962*b5*b8*a8+2912*b7*b8*a6-1008*b4^2*a5-23016*b7*a9*a8-54768*b7*b6*a5+6867*b7*b5*a6+48440*b7*a9*a4+168*b7*b6*a4+28896*b7*b6*a9-7336*b7*b10*a8-13230*b5*b9*a8+16240*b7*b6*a8-17584*b7*a9^2-44688*b7*a5^2+20944*b7*b10*a5+30051*b5*b9*a9+7336*b7*a8^2-23520*a5*b7*a4-9072*b8^2*a9+13104*b8^2*a8-7119*b9*b5*a4-114093*b5*a9*a3+11970*b5*b8*a9-30555*b5*b4*a9+11088*b8^2*a4+11466*b9*a9*a3+42084*b9*a5*a3-18270*b9*b8*a9-16002*b9*b8*a8+10269*b9*b4*a5+73521*b9*b8*a5+46557*b9*b4*a9+93177*b5*a8*a3-36351*b5*a8*a7+28476*b5*a5*a3-9450*b9*a5*a7+52353*b9*a8*a3-61362*b9*a4*a3-14553*b9*a9*a7+27027*b9*a4*a7-17703*b9*a8*a7-5040*b9*b7*a10+14315*b9*b7*a6"),
		SMQP("-35847*a6*a3*a4+15351*a6*a3*a5-20412*a6*b5*b6-10899*a6*b5*a5-23520*a6*a3*b6-46788*a6*b4*b6-29337*a6*b4*a5-77994*a6*b5*a4+56112*a4*b6*a5-61992*a6*b4*a4-69048*a4^2*a5+3360*a4*b6^2+26040*a4*a5^2-26880*b6^2*a5-27888*b6*a5^2+7056*a4^3-25200*a4^2*b6-34272*a5^3+22400*b10*a5^2+34398*b10*b5*a6-10584*b8*a10*a5-40768*b10*b6*a9-14952*b10*b5*a10+11928*b10*b9*a10-10080*b10*a6*a7+18144*b10*a6*a3+41216*b6*a8*a4+9632*b6*b10*a4+20335*a6*a8*a7-40376*a9*b5*a10+24640*b6*b10*a8-37919*a6*a9*a7+44625*a6*a4*a7+5040*b10*a10*a7-6048*b6*a10*a7+23296*b6*a10*a3+75957*a6*b9*a4-27384*a6*a8*a3+7224*b10*a9*a8+8848*b10*a9*a4+560*b6*a6*a7-3584*b10*a9^2+6552*a6*b7*a10-59437*a6*b9*a9+46991*a6*b9*a8-33992*b6*a9*a8-12824*b10*a8*a4-38430*b10*b9*a6+6608*b5*b6*a10+66255*a6*b4*a9-5432*a8^3-26656*a6*b8*a8+4816*b6*a9^2-4872*b9*a10*a8+17640*b6*a8^2-4144*a8*a9^2-1176*a8*a4^2+12544*a8*b4*a10+3640*a8*b5*a10+32592*a9^2*a5-12488*a8^2*a4-3304*a8*a10*a3+20559*a6*a9*a3+8960*a9^3+8736*a9*a4^2+10696*a9*a8^2+13440*a9*a8*a4-112*a9*a10*a3-49112*a9*b4*a10+8064*a10*a9*a7-11088*a10*a4*a7+51576*a10*b5*a4-40656*a10*b9*a4+5712*a10*a4*a3-3192*a10*b8*a8+31920*a10*b9*a9+31808*b6*b10*a5-6510*b8*a6*a4-224*a9*b6*a4-27804*a8*b4*a6-1624*a5*b10*a8+336*b10*a4^2-616*b10*a8^2-7056*b10^2*a5+65128*a5*b6*a9+2856*a4*b8*a10+41328*a4*b4*a10-16296*b9*a10*a5-27664*a9^2*a4-5992*a5*b10*a4-3024*b10^2*a4-37072*b6^2*a8+37632*b6^2*a9+4368*b8*a10*a9+5390*b8*a6*a9+168*a5*b10*a9+11088*b10^2*a9-11088*a5*a10*a7+14238*b8*b10*a6-60592*a5^2*a8+21343*a5*b9*a6+43372*b9*b6*a6+88592*b4*a10*a5+17304*a5*a8^2+28098*a6*b5*a9-6363*a6*b5*a8+5481*a6^2*b7+60984*b4*b10*a6-6048*b10^2*a8-10080*b10*a10*a3+13160*a5*b5*a10+29176*a5*a10*a3+34048*a5*a9*a4-40418*b8*b6*a6-79800*a5*a9*a8+1176*b8*b10*a10-30569*b8*a6*a5+25424*b4*b6*a10-7896*b8*b6*a10-18816*b4*b10*a10+84224*a5*a8*a4-96488*a5*b6*a8+38752*a5^2*a9-7728*b9*b6*a10-4312*a5*a6*a7"),
		SMQP("-23338*a5*b5*a8-26439*b4*a5*a4-7826*b8*b4*a10-7588*b9*b6*a9+336*b10*b4*a5-10752*b7*a10*a5+1967*b10*b9*a9-616*b8*a6*a3-693*b8*a9*a4+8428*a6*a3*a7-11256*b7*a10*a4+15393*b7*a6*a4+2352*b7*a10*a9-840*b6*a4*a7+7252*b6*a8*a7-6636*b6*a5*a3+196*a5*a8*a7+49462*a4*a3*a8-336*b6*a4*a3-5390*b8*b6*a5-2534*b6*a8*a3-13706*a8^2*a3+3682*b5*a4*a9-4172*b6*b5*a5+4200*b6*a5*a7+7350*b7*a6*a5+2513*b8*a8*a9+1855*b8*a8^2-4683*a6*b7*a8-20384*a4*a7*a8-15372*b8*a8*a4+7938*b5*a8^2+2380*a8^2*a7-5614*b8*a10*a3+4004*a7*a9^2+16128*a7*a4^2-1351*b5*a4^2+4438*b4*a4^2-9849*b4*a9^2+511*a8*b9*a9-4025*a8*b4*a9+4214*a10*a3^2-7700*b4*a10*a7+3696*a8*a9*a7-4585*a8*b5*a9-10038*a8*a9*a3+5208*a7*a5*a4-2912*a7^2*a6+6552*b9*a10*a7-7588*a10*a3*a7+7056*a8*b7*a10+5348*a9*a4*a3-5047*a9^2*b5+4676*a9*b9*a4-3304*a9*a4*a7-13384*a9*a5*a7+1498*a9^2*b9+16940*a9*a5*a3+3024*a7^2*a10-2688*b6*b9*a8-29610*a4^2*a3-6972*a3^2*a6+11340*b4*a10*a3+1295*b10*b9*a4-3458*b4^2*a10-17346*a3*a5*a4-658*b4*b6*a4+672*b7*b6*a10+6972*b7*b6*a6+1848*b7*b10*a10-11592*a5^2*a3-3647*b8*b6*a9+1337*b8*b10*a9+2156*b4*b6*a9+10591*a4*b4*a9-3458*b8*b5*a10+8568*b8*b9*a10-7560*b7*b10*a6+5824*b9*a4^2+8806*b10*b5*a5-16044*b5^2*a6-20454*b4*a6*a3-2464*b8^2*a6+5964*b8*a6*a7+3136*b10*a9*a7-2926*b4*b10*a8-1400*b10*a4*a7-112*b10*a8*a7-6209*b10*b5*a9-6832*b9*a6*a7-6174*b4*b9*a10+9338*b9*a6*a3-6790*b9*b8*a6+6384*b5^2*a10-3535*b8*b10*a4+17570*b4*b8*a6+14686*b4*b5*a10-13055*b8*a5*a9-3045*b10*b5*a8+2093*b4*a5*a9+2527*b10*b5*a4-5586*b9*a10*a3-1036*b4*b9*a6+4207*b8*b6*a4+406*b10*a9*a3+3304*b10*a5*a3+5593*b8*b6*a8-5110*b9*b10*a5+9100*b4*b6*a8+4676*b4*a8*a4+7602*b6*a9*a3-10164*b6*a9*a7-1414*b10*a8*a3+1358*b8*b10*a5-9786*b4*b5*a6-952*a9^2*a3+1050*b4^2*a6+910*b10*a4*a3-1617*a5*b9*a9+14812*b4*a6*a7+4718*a5^2*b4-8379*b7*a6*a9+2583*b9*b10*a8+6440*b9*b6*a5+182*b9*b6*a4+18095*b4*a5*a8+25277*b8*a5*a4-6713*b9*a8*a4-15519*b9*a5*a4+1162*b8*a9^2-1386*b4*b10*a4+11634*b8*a5^2+2555*b8*a4^2-5180*b10*a5*a7-2968*b4*b6*a5+700*b4*a8^2-13314*b5*a6*a3+6636*b5*a6*a7+13146*b5*b9*a6-2653*b5*a8*a4+14098*b5*a10*a3-8036*b5*a10*a7-7518*b5*b9*a10+10136*b5*a5*a4+8946*b5*b8*a6+18620*a5*b5*a9+2205*b9*a8^2+4753*a5*b9*a8-3633*a5*b8*a8-2520*b8*a10*a7+1302*b10*b4*a9-322*b9*a5^2-2198*b5*b6*a4+462*b9^2*a6-2520*b9^2*a10-9702*b5*b6*a8-679*b8*b10*a8-2772*a5^2*b5+16282*b5*b6*a9+12712*a5*a3*a8+6048*a5^2*a7"),
		SMQP("66192*a5*b5*a8+8337*b4*a5*a4-1526*b8*b4*a10+20244*b9*b6*a9-8400*b10*b4*a5+17248*b7*a10*a5-9051*b10*b9*a9-56*b8*a6*a3+1253*b8*a9*a4-36316*a6*a3*a7+21336*b7*a10*a4-32403*b7*a6*a4-15120*b7*a10*a9-4536*b6*a4*a7-11844*b6*a8*a7+17724*b6*a5*a3+25452*a5*a8*a7-72688*a4*a3*a8+4368*b6*a4*a3+7182*b8*b6*a5-4018*b6*a8*a3+19586*a8^2*a3-38178*b5*a4*a9+16716*b6*b5*a5-11592*b6*a5*a7-30534*b7*a6*a5+1099*b8*a8*a9-1855*b8*a8^2+26817*a6*b7*a8+25872*a4*a7*a8+5096*b8*a8*a4-20300*b5*a8^2+6342*a8^2*a7-1498*b8*a10*a3-9030*a7*a9^2-35154*a7*a4^2+23499*b5*a4^2+9954*b4*a4^2+19131*b4*a9^2-8827*a8*b9*a9+2723*a8*b4*a9-11774*a10*a3^2+17220*b4*a10*a7-22176*a8*a9*a7+27475*a8*b5*a9+19152*a8*a9*a3-46872*a7*a5*a4+10080*a7^2*a6-11256*b9*a10*a7+24276*a10*a3*a7-10416*a8*b7*a10-20132*a9*a4*a3+16611*a9^2*b5+5628*a9*b9*a4+30660*a9*a4*a7+20664*a9*a5*a7-9072*a9^2*b9-31514*a9*a5*a3-8400*a7^2*a10+112*b6*b9*a8+52416*a4^2*a3+22092*a3^2*a6-20524*b4*a10*a3+3381*b10*b9*a4-10374*b4^2*a10+68124*a3*a5*a4-5670*b4*b6*a4+3360*b7*b6*a10-14364*b7*b6*a6-3864*b7*b10*a10+23688*a5^2*a3+9723*b8*b6*a9-4109*b8*b10*a9-14028*b4*b6*a9-29085*a4*b4*a9+17290*b8*b5*a10+168*b8*b9*a10+7560*b7*b10*a6-13314*b9*a4^2-16758*b10*b5*a5+40572*b5^2*a6+21630*b4*a6*a3+12656*b8^2*a6-1148*b8*a6*a7-5712*b10*a9*a7+3766*b4*b10*a8+6888*b10*a4*a7-2688*b10*a8*a7+21189*b10*b5*a9+15120*b9*a6*a7+28854*b4*b9*a10-35378*b9*a6*a3+13454*b9*b8*a6-12432*b5^2*a10+427*b8*b10*a4+1974*b4*b8*a6-55398*b4*b5*a10-3871*b8*a5*a9-847*b10*b5*a8-46011*b4*a5*a9-8043*b10*b5*a4+13146*b9*a10*a3-46452*b4*b9*a6-6699*b8*b6*a4+1778*b10*a9*a3-7000*b10*a5*a3-3045*b8*b6*a8-1050*b9*b10*a5+5572*b4*b6*a8-12964*b4*a8*a4-16170*b6*a9*a3+21924*b6*a9*a7+7294*b10*a8*a3-2030*b8*b10*a5+39186*b4*b5*a6+7210*a9^2*a3+3150*b4^2*a6-4774*b10*a4*a3-3213*a5*b9*a9-17052*b4*a6*a7+28770*a5^2*b4+18249*b7*a6*a9-1211*b9*b10*a8-2520*b9*b6*a5-5838*b9*b6*a4+18907*b4*a5*a8-28707*b8*a5*a4-9408*b8^2*a10+23513*b9*a8*a4-231*b9*a5*a4-406*b8*a9^2-798*b4*b10*a4-32130*b8*a5^2+4473*b8*a4^2+8316*b10*a5*a7+4200*b4*b6*a5+308*b4*a8^2+39690*b5*a6*a3-4284*b5*a6*a7-40866*b5*b9*a6-11585*b5*a8*a4-22666*b5*a10*a3+10836*b5*a10*a7+22806*b5*b9*a10-37506*b5*a5*a4-4410*b5*b8*a6-55062*a5*b5*a9-5467*b9*a8^2+8281*a5*b9*a8+28987*a5*b8*a8+3192*b8*a10*a7+5250*b10*b4*a9+20370*b9*a5^2-2898*b5*b6*a4+1386*b9^2*a6-7560*b9^2*a10+31934*b5*b6*a8+1771*b8*b10*a8+18900*a5^2*b5-46578*b5*b6*a9-59542*a5*a3*a8-2016*a5^2*a7"),
		SMQP("63378*a5*b5*a8+4221*b4*a5*a4-2450*b8*b4*a10+16380*b9*b6*a9-3668*b10*b4*a5+11312*b7*a10*a5-8141*b10*b9*a9-168*b8*a6*a3-973*b8*a9*a4-39060*a6*a3*a7+18648*b7*a10*a4-28035*b7*a6*a4-15120*b7*a10*a9-7560*b6*a4*a7-6972*b6*a8*a7+16884*b6*a5*a3+33348*a5*a8*a7-44478*a4*a3*a8+9072*b6*a4*a3+9618*b8*b6*a5-8694*b6*a8*a3+11718*a8^2*a3-32214*b5*a4*a9+12852*b6*b5*a5-12600*b6*a5*a7-21546*b7*a6*a5+1449*b8*a8*a9-665*b8*a8^2+22617*a6*b7*a8+9408*a4*a7*a8+1540*b8*a8*a4-13986*b5*a8^2+5376*a8^2*a7-2590*b8*a10*a3-12348*a7*a9^2-19908*a7*a4^2+26481*b5*a4^2+3318*b4*a4^2+14567*b4*a9^2-8253*a8*b9*a9-189*a8*b4*a9-10458*a10*a3^2+6300*b4*a10*a7-12516*a8*a9*a7+16023*a8*b5*a9+18270*a8*a9*a3-48888*a7*a5*a4+13440*a7^2*a6-7224*b9*a10*a7+19404*a10*a3*a7-7728*a8*b7*a10-17640*a9*a4*a3+20265*a9^2*b5+9800*a9*b9*a4+26628*a9*a4*a7+11928*a9*a5*a7-8806*a9^2*b9-16884*a9*a5*a3-7728*a7^2*a10-896*b6*b9*a8+31122*a4^2*a3+25956*a3^2*a6-8260*b4*a10*a3+3283*b10*b9*a4-3458*b4^2*a10+66402*a3*a5*a4-1890*b4*b6*a4-1008*b7*b6*a10-11844*b7*b6*a6-2184*b7*b10*a10+12600*a5^2*a3+3297*b8*b6*a9-623*b8*b10*a9-10836*b4*b6*a9-17437*a4*b4*a9+14686*b8*b5*a10+4872*b8*b9*a10+6972*b7*b10*a6-14532*b9*a4^2-9534*b10*b5*a5+39060*b5^2*a6+18858*b4*a6*a3+10976*b8^2*a6-308*b8*a6*a7-6804*b10*a9*a7+1106*b4*b10*a8+6132*b10*a4*a7-924*b10*a8*a7+17535*b10*b5*a9+7168*b9*a6*a7+31906*b4*b9*a10-20790*b9*a6*a3+4186*b9*b8*a6-10752*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-58898*b4*b5*a10+161*b8*a5*a9-2793*b10*b5*a8-27475*b4*a5*a9-5817*b10*b5*a4+9086*b9*a10*a3-35980*b4*b9*a6-3129*b8*b6*a4+1386*b10*a9*a3-2772*b10*a5*a3-2471*b8*b6*a8-3962*b9*b10*a5+2380*b4*b6*a8-4172*b4*a8*a4-10206*b6*a9*a3+20412*b6*a9*a7+9450*b10*a8*a3-518*b8*b10*a5+31878*b4*b5*a6+5544*a9^2*a3+1050*b4^2*a6-7686*b10*a4*a3+2639*a5*b9*a9-13076*b4*a6*a7+19278*a5^2*b4+15813*b7*a6*a9-1645*b9*b10*a8+1512*b9*b6*a5-1050*b9*b6*a4+6475*b4*a5*a8-15519*b8*a5*a4-10080*b8^2*a10+20839*b9*a8*a4-315*b9*a5*a4+266*b8*a9^2-266*b4*b10*a4-33390*b8*a5^2+3171*b8*a4^2+2436*b10*a5*a7+168*b4*b6*a5+28*b4*a8^2+19278*b5*a6*a3+252*b5*a6*a7-42630*b5*b9*a6-15729*b5*a8*a4-11550*b5*a10*a3+8316*b5*a10*a7+20146*b5*b9*a10-33516*b5*a5*a4-5166*b5*b8*a6-56868*a5*b5*a9-3815*b9*a8^2-1043*a5*b9*a8+21595*a5*b8*a8+3192*b8*a10*a7+5362*b10*b4*a9+16926*b9*a5^2-4998*b5*b6*a4+9422*b9^2*a6-10248*b9^2*a10+28602*b5*b6*a8+833*b8*b10*a8+21420*a5^2*b5-41286*b5*b6*a9-63252*a5*a3*a8+2016*a5^2*a7"),
		SMQP("20055*a6*a3*a4+2457*a6*a3*a5+66948*a6*b5*b6+5859*a6*b5*a5+7392*a6*a3*b6+16380*a6*b4*b6+29757*a6*b4*a5+58422*a6*b5*a4-4704*a4*b6*a5-2688*a6*b4*a4+11256*a4^2*a5-5376*a4*b6^2+28224*a4*a5^2+12768*b6^2*a5+12096*b6*a5^2+12600*a4^3+7728*a4^2*b6+16128*a5^3-24864*b10*a5^2-70770*b10*b5*a6+8316*b8*a10*a5+15904*b10*b6*a9+22568*b10*b5*a10-11592*b10*b9*a10-4480*b10*a6*a7+14784*b10*a6*a3+2408*b6*a8*a4+5824*b6*b10*a4+1169*a6*a8*a7+16856*a9*b5*a10-11312*b6*b10*a8+13055*a6*a9*a7-22001*a6*a4*a7-3024*b10*a10*a7+9072*b6*a10*a7-11200*b6*a10*a3-11557*a6*b9*a4-12936*a6*a8*a3-4872*b10*a9*a8+4648*b10*a9*a4-5600*b6*a6*a7-1008*b10*a9^2+26964*a6*b7*a10+19957*a6*b9*a9+25165*a6*b9*a8+6104*b6*a9*a8+20160*b10*a8*a4+28210*b10*b9*a6-23072*b5*b6*a10-49203*a6*b4*a9+2744*a8^3+6244*a6*b8*a8+14336*b6*a9^2-9828*b9*a10*a8-15288*b6*a8^2-6272*a8*a9^2-34496*a8*a4^2+11648*a8*b4*a10+1820*a8*b5*a10+12544*a9^2*a5-3780*a8*a10*a7+18816*a8^2*a4+17080*a8*a10*a3-9471*a6*a9*a3-224*a9^3-4816*a9*a4^2-2296*a9*a8^2+9912*a9*a8*a4+4060*a9*a10*a3+27692*a9*b4*a10-8316*a10*a9*a7+9828*a10*a4*a7-36512*a10*b5*a4+7812*a10*b9*a4-8596*a10*a4*a3+2016*a10*b8*a8-6048*a10^2*b7-14868*a10*b9*a9-2240*b6*b10*a5+11858*b8*a6*a4-21784*a9*b6*a4+11844*a8*b4*a6-8960*a5*b10*a8-12488*b10*a4^2-3192*b10*a8^2+1568*b10^2*a5-24976*a5*b6*a9-2520*a4*b8*a10-16856*a4*b4*a10+5292*b9*a10*a5+4032*a9^2*a4+4816*a5*b10*a4+560*b10^2*a4+10864*b6^2*a8-16464*b6^2*a9+4536*b8*a10*a9-25046*b8*a6*a9+25536*a5*b10*a9-3472*b10^2*a9+8064*a5*a10*a7-9506*b8*b10*a6-18816*a5^2*a8-23499*a5*b9*a6-47572*b9*b6*a6-7476*b4*a10*a5-20944*a5*a8^2-15918*a6*b5*a9-10941*a6*b5*a8-34965*a6^2*b7+5208*b4*b10*a6+448*b10^2*a8+2128*b10*a10*a3-3948*a5*b5*a10-9156*a5*a10*a3-25984*a5*a9*a4-14098*b8*b6*a6+22176*a5*a9*a8-3528*b8*b10*a10-21315*b8*a6*a5-20048*b4*b6*a10+11592*b8*b6*a10+14000*b4*b10*a10+8680*a5*a8*a4+3472*a5*b6*a8-21504*a5^2*a9+26208*b9*b6*a10-13272*a5*a6*a7"),
		SMQP("6951*a6*a3*a4+3885*a6*a3*a5+33012*a6*b5*b6+4347*a6*b5*a5+2100*a6*a3*b6+10248*a6*b4*b6+4620*a6*b4*a5+28161*a6*b5*a4-3990*a4*b6*a5+2688*a6*b4*a4+798*a4^2*a5-2352*a4*b6^2+15330*a4*a5^2+4452*b6^2*a5+9660*b6*a5^2+8820*a4^3+3066*a4^2*b6+9576*a5^3-15316*b10*a5^2-35973*b10*b5*a6+4536*b8*a10*a5+3934*b10*b6*a9+9352*b10*b5*a10-6132*b10*b9*a10-6776*b10*a6*a7+12432*b10*a6*a3+2198*b6*a8*a4+5698*b6*b10*a4+2429*a6*a8*a7+8246*a9*b5*a10-6174*b6*b10*a8+7427*a6*a9*a7-7385*a6*a4*a7-1512*b10*a10*a7+4536*b6*a10*a7-5194*b6*a10*a3-8323*a6*b9*a4-10164*a6*a8*a3-3192*b10*a9*a8+3878*b10*a9*a4-2660*b6*a6*a7-1022*b10*a9^2+15120*a6*b7*a10+3745*a6*b9*a9+13720*a6*b9*a8+7658*b6*a9*a8+12782*b10*a8*a4+16709*b10*b9*a6-13118*b5*b6*a10-17304*a6*b4*a9+1358*a8^3-3353*a6*b8*a8+5222*b6*a9^2-5712*b9*a10*a8-7560*b6*a8^2-3122*a8*a9^2-17486*a8*a4^2+5054*a8*b4*a10+3122*a8*b5*a10+980*a9^2*a5-1764*a8*a10*a7+7098*a8^2*a4+7882*a8*a10*a3-2667*a6*a9*a3-140*a9^3-910*a9*a4^2-1120*a9*a8^2+1890*a9*a8*a4+490*a9*a10*a3+8162*a9*b4*a10-4284*a10*a9*a7+4788*a10*a4*a7-14672*a10*b5*a4+1596*a10*b9*a4-5824*a10*a4*a3+924*a10*b8*a8-3024*a10^2*b7-4956*a10*b9*a9+1540*b6*b10*a5+4025*b8*a6*a4-10864*a9*b6*a4+840*a8*b4*a6-14966*a5*b10*a8-11984*b10*a4^2-322*b10*a8^2+1792*b10^2*a5-15610*a5*b6*a9+1428*a4*b8*a10-9212*a4*b4*a10+2184*b9*a10*a5+3612*a9^2*a4+6510*a5*b10*a4+1036*b10^2*a4+4858*b6^2*a8-7518*b6^2*a9+672*b8*a10*a9+784*b8*a6*a9+19446*a5*b10*a9-476*b10^2*a9+4032*a5*a10*a7-5089*b8*b10*a6-8344*a5^2*a8-8792*a5*b9*a6-23380*b9*b6*a6-2338*b4*a10*a5-7742*a5*a8^2-14805*a6*b5*a9+315*a6*b5*a8-16506*a6^2*b7-2352*b4*b10*a6-1036*b10^2*a8+2828*b10*a10*a3-3010*a5*b5*a10-2534*a5*a10*a3-21882*a5*a9*a4-4165*b8*b6*a6+17010*a5*a9*a8-1428*b8*b10*a10-7280*b8*a6*a5-11774*b4*b6*a10+4620*b8*b6*a10+9940*b4*b10*a10+8736*a5*a8*a4+3654*a5*b6*a8-9464*a5^2*a9+13776*b9*b6*a10-8260*a5*a6*a7"),
		SMQP("-44247*a6*a3*a4+19383*a6*a3*a5-18060*a6*b5*b6-15939*a6*b5*a5-20496*a6*a3*b6-52164*a6*b4*b6-40257*a6*b4*a5-73290*a6*b5*a4+55440*a4*b6*a5-60816*a6*b4*a4-77280*a4^2*a5+6384*a4*b6^2+28056*a4*a5^2-26880*b6^2*a5-27888*b6*a5^2+17640*a4^3-31920*a4^2*b6-34272*a5^3+22400*b10*a5^2+31206*b10*b5*a6-15120*b8*a10*a5-41104*b10*b6*a9-15008*b10*b5*a10+14112*b10*b9*a10-13664*b10*a6*a7+25536*b10*a6*a3+54376*b6*a8*a4+15344*b6*b10*a4+15743*a6*a8*a7-53256*a9*b5*a10+21392*b6*b10*a8-48783*a6*a9*a7+48881*a6*a4*a7+5040*b10*a10*a7-6048*b6*a10*a7+22288*b6*a10*a3+89005*a6*b9*a4-25032*a6*a8*a3+10136*b10*a9*a8+10192*b10*a9*a4+560*b6*a6*a7-4984*b10*a9^2+12096*a6*b7*a10-67557*a6*b9*a9+76783*a6*b9*a8-30408*b6*a9*a8-10472*b10*a8*a4-49126*b10*b9*a6+7952*b5*b6*a10+60879*a6*b4*a9-5488*a8^3-39368*a6*b8*a8+7672*b6*a9^2-10080*b9*a10*a8-4256*b6*a8^2-9464*a8*a9^2-18760*a8*a4^2+37688*a8*b4*a10+7952*a8*b5*a10+49168*a9^2*a5-1008*a8*a10*a7+7000*a8^2*a4-896*a8*a10*a3+27279*a6*a9*a3+12432*a9^3+616*a9*a4^2+12600*a9*a8^2+20272*a9*a8*a4+5880*a9*a10*a3-44184*a9*b4*a10+9072*a10*a9*a7-10080*a10*a4*a7+48776*a10*b5*a4-48384*a10*b9*a4+5656*a10*a4*a3-3024*a10*b8*a8+35280*a10*b9*a9+32480*b6*b10*a5-4214*b8*a6*a4+4984*a9*b6*a4-22260*a8*b4*a6-392*a5*b10*a8-2968*b10*a4^2-2128*b10*a8^2-5600*b10^2*a5+65464*a5*b6*a9+2016*a4*b8*a10+37688*a4*b4*a10-12096*b9*a10*a5-36848*a9^2*a4-10136*a5*b10*a4-4592*b10^2*a4-37072*b6^2*a8+37632*b6^2*a9+6048*b8*a10*a9-13986*b8*a6*a9+8344*a5*b10*a9+9520*b10^2*a9-11088*a5*a10*a7+18326*b8*b10*a6-90832*a5^2*a8+11431*a5*b9*a6+44044*b9*b6*a6+97832*b4*a10*a5-20216*a5*a8^2+42882*a6*b5*a9-27363*a6*b5*a8-1071*a6^2*b7+60312*b4*b10*a6-4480*b10^2*a8-10192*b10*a10*a3+15176*a5*b5*a10+29176*a5*a10*a3+45416*a5*a9*a4-47642*b8*b6*a6-69944*a5*a9*a8-1008*b8*b10*a10-36113*b8*a6*a5+28112*b4*b6*a10-3024*b8*b6*a10-18032*b4*b10*a10+121576*a5*a8*a4-100856*a5*b6*a8+32704*a5^2*a9-3024*b9*b6*a10-4312*a5*a6*a7"),
		SMQP("-1575*b9^2*a9+1071*b9^2*a4+63*b9^2*a8+952*b7*b4*a10+1036*b7*b10*a9-3220*b7*a10*a3-1029*b7*b4*a6-476*b7*b10*a4-3402*a9*a3*a7-5418*b4*a4*a3-4536*b8*a8*a3-2457*b8*b9*a4+756*b8*a4*a3-63*b8*a9*a7+7896*b7*a6*a3+5985*b4*a8*a3+315*b4*a8*a7-1512*b5^2*a5+2016*b4*b5*a4-3374*a7*b7*a6-15876*a7*a8*a3-3843*a7*b5*a9+756*a7*b5*a5+1512*b8*b7*a10+1890*a7*b4*a5+2268*a7*b7*a10-7245*a7*b4*a9-5292*a7*b8*a5+2331*a7*b8*a8-1512*b8*a9*a3+12852*b8*a5*a3+378*b8*b4*a9+2205*b8*b4*a5+441*b8*b5*a4+21546*a7*a4*a3-5670*a7^2*a4-882*b4*b8*a8-945*b4*b9*a8+2961*b4*a4*a7-2835*b4*b5*a8+1134*a7^2*a9+315*a7*b8*a4+2268*a9*a3^2-63*b7*a4^2-2387*b7*a4*a8+756*b8^2*a5-7371*a3*b4*a5-20412*a3^2*a4-2268*a3^2*a5+9198*a3*b4*a9+2009*b7*a5*a9+1890*b9^2*a5-8141*b7*a5*a8+2079*b5*a4*a3-2072*b7*b5*a10-1197*b5*a4*a7-756*b5^2*a4+1134*a5*a3*a7+1512*b5^2*a9+3402*a7^2*a8+1827*b4*b9*a4+18144*a3^2*a8+1890*b5*b9*a5+1134*b5*b8*a5+4977*b5*b4*a5+1701*b4^2*a9-63*b5*b8*a8-2632*b7*b8*a6+189*b4^2*a5+1428*b7*a9*a8+3318*b7*b6*a5-756*b7*b5*a6-1393*b7*a9*a4+3003*b7*b6*a4-2667*b7*b6*a9+476*b7*b10*a8+63*b5*b9*a8-2345*b7*b6*a8+560*b7*a9^2+378*b7*a5^2-476*b7*b10*a5-1575*b5*b9*a9+1036*b7*a8^2+4893*a5*b7*a4-126*b8^2*a9+126*b8^2*a8+1071*b9*b5*a4+7749*b5*a9*a3-1449*b5*b8*a9+252*b5*b4*a9+630*b8^2*a4+693*b9*a9*a3-4914*b9*a5*a3+3213*b9*b8*a9-189*b9*b8*a8-1323*b9*b4*a5-4158*b9*b8*a5-1575*b9*b4*a9-6993*b5*a8*a3+3465*b5*a8*a7-1890*b5*a5*a3+1512*b9*a5*a7-2205*b9*a8*a3+3339*b9*a4*a3+441*b9*a9*a7-2205*b9*a4*a7+1071*b9*a8*a7+644*b9*b7*a6"),
		SMQP("-252*b9^2*a9+924*b9^2*a4+1092*b9^2*a8+2688*b7*b4*a10-1344*b7*b10*a9+5712*b7*a10*a3-1575*b7*b4*a6+672*b7*b10*a4+9828*a9*a3*a7+3507*b4*a4*a3+6972*b8*a8*a3-357*b8*b9*a4-1113*b8*a4*a3+1071*b8*a9*a7-5082*b7*a6*a3-2772*b4*a8*a3-231*b4*a8*a7+1008*b5^2*a5+3612*b4*b5*a4+2058*a7*b7*a6+15414*a7*a8*a3+3444*a7*b5*a9+126*a7*b5*a5-1344*b8*b7*a10-5250*a7*b4*a5-3024*a7*b7*a10+7581*a7*b4*a9+2478*a7*b8*a5-3255*a7*b8*a8-315*b8*a9*a3-10101*b8*a5*a3+861*b8*b4*a9-126*b8*b4*a5-1008*b8*b5*a4+336*a7^2*a5-28560*a7*a4*a3+7686*a7^2*a4+210*b4*b8*a8+63*b4*b9*a8-2121*b4*a4*a7+3255*b4*b5*a8-2730*a7^2*a9+1323*a7*b8*a4-4746*a9*a3^2-4956*b7*a4^2+5292*b7*a4*a8+651*b8^2*a5+10227*a3*b4*a5+26250*a3^2*a4+3990*a3^2*a5-12705*a3*b4*a9+3780*b7*a5*a9-1848*b9^2*a5-1372*b7*a5*a8-2940*b5*a4*a3+1344*b7*b5*a10+1764*b5*a4*a7-504*b5^2*a8-3318*a5*a3*a7-672*b5^2*a9-3318*a7^2*a8-3969*b4*b9*a4-20832*a3^2*a8-1848*b5*b9*a5+1071*b5*b8*a5+63*b5*b4*a5-1449*b4^2*a9-1701*b5*b8*a8-63*b7*b8*a6+231*b4^2*a5-672*b7*a9*a8-3528*b7*b6*a5-2016*b7*b5*a6+4284*b7*a9*a4+924*b7*b6*a4+2436*b7*b6*a9-672*b7*b10*a8+84*b5*b9*a8-252*b7*b6*a8-2016*b7*a9^2-3864*b7*a5^2+2352*b7*b10*a5+1764*b5*b9*a9-1344*b7*a8^2-2156*a5*b7*a4-1050*b8^2*a9+42*b8^2*a8-1092*b9*b5*a4-11340*b5*a9*a3+1932*b5*b8*a9-5376*b5*b4*a9+1176*b8^2*a4+1386*b9*a9*a3+4914*b9*a5*a3-315*b9*b8*a9-1533*b9*b8*a8-2625*b9*b4*a5+2079*b9*b8*a5+4305*b9*b4*a9+8358*b5*a8*a3-2646*b5*a8*a7-1470*b5*a5*a3-210*b9*a5*a7+3822*b9*a8*a3-5754*b9*a4*a3-1050*b9*a9*a7+2058*b9*a4*a7-1806*b9*a8*a7-672*b9*b7*a10+2016*b9*b7*a6"),
		SMQP("-693*b9^2*a9+1197*b9^2*a4-399*b9^2*a8-1456*b7*b4*a10+812*b7*b10*a9-4340*b7*a10*a3-504*b7*b4*a6-1036*b7*b10*a4-2646*a9*a3*a7-2709*b4*a4*a3-3864*b8*a8*a3-2667*b8*b9*a4+924*b8*a4*a3+147*b8*a9*a7+2856*b7*a6*a3+1911*b4*a8*a3+1008*b4*a8*a7-840*b5^2*a5+756*b4*b5*a4-658*a7*b7*a6-12348*a7*a8*a3-2457*a7*b5*a9+420*a7*b5*a5+504*b8*b7*a10+3822*a7*b4*a5+2100*a7*b7*a10-6174*a7*b4*a9-3948*a7*b8*a5+1701*a7*b8*a8-840*b8*a9*a3+9660*b8*a5*a3-420*b8*b4*a5+651*b8*b5*a4+16758*a7*a4*a3-4410*a7^2*a4+168*b4*b8*a8-1050*b4*b9*a8+1428*b4*a4*a7-1344*b4*b5*a8+882*a7^2*a9+273*a7*b8*a4+1764*a9*a3^2+2667*b7*a4^2-3829*b7*a4*a8+924*b8^2*a5-7938*a3*b4*a5-15876*a3^2*a4-1764*a3^2*a5+8169*a3*b4*a9-3605*b7*a5*a9+1302*b9^2*a5+749*b7*a5*a8+861*b5*a4*a3-448*b7*b5*a10-567*b5*a4*a7-924*b5^2*a4+336*b5^2*a8+882*a5*a3*a7+840*b5^2*a9+2646*a7^2*a8+1386*b4*b9*a4+14112*a3^2*a8+1302*b5*b9*a5+378*b5*b8*a5+1596*b5*b4*a5+672*b4^2*a9-273*b5*b8*a8+1288*b7*b8*a6+168*b4^2*a5+1092*b7*a9*a8+2226*b7*b6*a5+1932*b7*b5*a6-1379*b7*a9*a4+945*b7*b6*a4-1785*b7*b6*a9+700*b7*b10*a8-399*b5*b9*a8+161*b7*b6*a8+448*b7*a9^2+3486*b7*a5^2-1540*b7*b10*a5-693*b5*b9*a9+476*b7*a8^2+1407*a5*b7*a4+294*b8^2*a9-126*b8^2*a8+1197*b9*b5*a4+4767*b5*a9*a3-987*b5*b8*a9+2436*b5*b4*a9+546*b8^2*a4+1071*b9*a9*a3-3990*b9*a5*a3+1239*b9*b8*a9+861*b9*b8*a8+1176*b9*b4*a5-3066*b9*b8*a5-1974*b9*b4*a9-4431*b5*a8*a3+2247*b5*a8*a7-1302*b5*a5*a3+1344*b9*a5*a7-2163*b9*a8*a3+2961*b9*a4*a3-189*b9*a9*a7-2079*b9*a4*a7+1281*b9*a8*a7+336*b9*b7*a10-1820*b9*b7*a6"),
		SMQP("16877*a5*b5*a8+1428*b4*a5*a4+2016*b9*b6*a9-1617*b10*b4*a5+2184*b7*a10*a5-1099*b10*b9*a9-1960*b8*a6*a3-588*b8*a9*a4-13664*a6*a3*a7+5544*b7*a10*a4-8064*b7*a6*a4-4032*b7*a10*a9-3360*b6*a4*a7-728*b6*a8*a7+3696*b6*a5*a3+16282*a5*a8*a7-5201*a4*a3*a8+4704*b6*a4*a3+3444*b8*b6*a5-6104*b6*a8*a3+3136*a8^2*a3-8792*b5*a4*a9+1512*b6*b5*a5-4368*b6*a5*a7-1638*b7*a6*a5-84*b8*a8*a9-126*b8*a8^2+5985*a6*b7*a8-2093*a4*a7*a8+504*b8*a8*a4-2401*b5*a8^2+889*a8^2*a7-1008*b8*a10*a3-3073*a7*a9^2-756*a7*a4^2+6552*b5*a4^2+1155*b4*a9^2-2100*a8*b9*a9+903*a8*b4*a9-2296*a10*a3^2-1400*b4*a10*a7-1848*a8*a9*a7+2583*a8*b5*a9+4641*a8*a9*a3-16842*a7*a5*a4+5488*a7^2*a6+504*b9*a10*a7+4424*a10*a3*a7-2016*a8*b7*a10-5005*a9*a4*a3+5068*a9^2*b5+3283*a9*b9*a4+6755*a9*a4*a7-4522*a9*a5*a7-1127*a9^2*b9+2555*a9*a5*a3-2016*a7^2*a10+728*b6*b9*a8+3276*a4^2*a3+8400*a3^2*a6+2800*b4*a10*a3+1127*b10*b9*a4+20706*a3*a5*a4-2520*b7*b6*a10-1764*b7*b6*a6-2142*a5^2*a3+84*b8*b6*a9+210*b8*b10*a9-1848*b4*b6*a9-1932*a4*b4*a9+5544*b8*b5*a10+5544*b8*b9*a10+1953*b7*b10*a6-3780*b9*a4^2-2485*b10*b5*a5+8064*b5^2*a6+3444*b4*a6*a3+4536*b8^2*a6+1484*b8*a6*a7-3017*b10*a9*a7+1015*b10*a4*a7+1001*b10*a8*a7+4172*b10*b5*a9-980*b9*a6*a7+10472*b4*b9*a10-1036*b9*a6*a3-1988*b9*b8*a6-616*b5^2*a10-18760*b4*b5*a10+2667*b8*a5*a9-497*b10*b5*a8-1995*b4*a5*a9-1708*b10*b5*a4+112*b9*a10*a3-10752*b4*b9*a6+553*b10*a9*a3-455*b10*a5*a3-2317*b9*b10*a5+168*b6*a9*a3+6216*b6*a9*a7+1904*b10*a8*a3-273*b8*b10*a5+10416*b4*b5*a6-847*a9^2*a3-1841*b10*a4*a3-4907*a5*b9*a9-1848*b4*a6*a7+7266*a5^2*b4+4221*b7*a6*a9-1001*b9*b10*a8+5208*b9*b6*a5-672*b9*b6*a4-2793*b4*a5*a8-672*b8*a5*a4-5040*b8^2*a10+4235*b9*a8*a4-2142*b9*a5*a4+294*b8*a9^2-13230*b8*a5^2+952*b10*a5*a7-1176*b4*b6*a5-84*b5*a6*a3+2324*b5*a6*a7-8036*b5*b9*a6-2632*b5*a8*a4+1400*b5*a10*a3+112*b5*a10*a7+2912*b5*b9*a10-8400*b5*a5*a4-5012*b5*b8*a6-9317*a5*b5*a9-889*b9*a8^2-2695*a5*b9*a8+6699*a5*b8*a8-504*b8*a10*a7+735*b10*b4*a9+6258*b9*a5^2-672*b5*b6*a4+3556*b9^2*a6-4536*b9^2*a10+5768*b5*b6*a8-126*b8*b10*a8+7686*a5^2*b5-9408*b5*b6*a9-21287*a5*a3*a8+5040*a5^2*a7"),
		SMQP("-16842*a5*b5*a8-7959*b4*a5*a4-2450*b8*b4*a10-10500*b9*b6*a9+16464*b10*b4*a5-12096*b7*a10*a5+1071*b10*b9*a9+24024*b8*a6*a3+203*b8*a9*a4+25452*a6*a3*a7-16632*b7*a10*a4+25137*b7*a6*a4+5040*b7*a10*a9+504*b6*a4*a7+756*b6*a8*a7-9324*b6*a5*a3-18396*a5*a8*a7+65142*a4*a3*a8-3024*b6*a4*a3+9282*b8*b6*a5+10458*b6*a8*a3-23562*a8^2*a3+2898*b5*a4*a9-14700*b6*b5*a5+5544*b6*a5*a7+11718*b7*a6*a5+2961*b8*a8*a9+3871*b8*a8^2-3339*a6*b7*a8-28224*a4*a7*a8-17948*b8*a8*a4+5250*b5*a8^2+7308*a8^2*a7-21742*b8*a10*a3+5796*a7*a9^2+20160*a7*a4^2-2583*b5*a4^2-27594*b4*a4^2-13881*b4*a9^2+2751*a8*b9*a9-11865*a8*b4*a9-3402*a10*a3^2-12852*b4*a10*a7+5040*a8*a9*a7-6153*a8*b5*a9-12726*a8*a9*a3+18648*a7*a5*a4-10080*a7^2*a6+6552*b9*a10*a7-3780*a10*a3*a7+7056*a8*b7*a10+13860*a9*a4*a3-4599*a9^2*b5+3780*a9*b9*a4-7560*a9*a4*a7-7560*a9*a5*a7+2394*a9^2*b9+5292*a9*a5*a3+3024*a7^2*a10-3360*b6*b9*a8-37674*a4^2*a3-12348*a3^2*a6+33516*b4*a10*a3+2751*b10*b9*a4+30366*b4^2*a10-44226*a3*a5*a4+23310*b4*b6*a4-2016*b7*b6*a10+10332*b7*b6*a6+4536*b7*b10*a10-11592*a5^2*a3-8799*b8*b6*a9+4921*b8*b10*a9+12684*b4*b6*a9+17871*a4*b4*a9-3458*b8*b5*a10+8568*b8*b9*a10-11592*b7*b10*a6+6048*b9*a4^2+12726*b10*b5*a5-13356*b5^2*a6+378*b4*a6*a3+6944*b8^2*a6-2996*b8*a6*a7-3150*b4*b10*a8-504*b10*a4*a7+3024*b10*a8*a7-8001*b10*b5*a9-12432*b9*a6*a7-11550*b4*b9*a10+18522*b9*a6*a3-15526*b9*b8*a6+6384*b5^2*a10-10367*b8*b10*a4-19614*b4*b8*a6+5502*b4*b5*a10-13279*b8*a5*a9-5061*b10*b5*a8+74781*b4*a5*a9+6447*b10*b5*a4-5586*b9*a10*a3+22932*b4*b9*a6+18879*b8*b6*a4+6678*b10*a9*a3-2520*b10*a5*a3-2471*b8*b6*a8-3654*b9*b10*a5-22932*b4*b6*a8+38052*b4*a8*a4+8946*b6*a9*a3-10836*b6*a9*a7-7686*b10*a8*a3-5474*b8*b10*a5-22554*b4*b5*a6-4536*a9^2*a3-11718*b4^2*a6-882*b10*a4*a3-3633*a5*b9*a9-9828*b4*a6*a7-46914*a5^2*b4-14427*b7*a6*a9+1911*b9*b10*a8+8904*b9*b6*a5+2646*b9*b6*a4-46641*b4*a5*a8+34461*b8*a5*a4-8841*b9*a8*a4-18543*b9*a5*a4+266*b8*a9^2+630*b4*b10*a4+22050*b8*a5^2+1323*b8*a4^2-2268*b10*a5*a7-7224*b4*b6*a5+4284*b4*a8^2-46242*b5*a6*a3+20748*b5*a6*a7+15162*b5*b9*a6+3171*b5*a8*a4+29778*b5*a10*a3-15876*b5*a10*a7-7518*b5*b9*a10+6888*b5*a5*a4+1554*b5*b8*a6+17388*a5*b5*a9+2877*b9*a8^2+5313*a5*b9*a8-13601*a5*b8*a8+5544*b8*a10*a7-2058*b10*b4*a9-3570*b9*a5^2-12726*b5*b6*a4+2478*b9^2*a6-2520*b9^2*a10-4326*b5*b6*a8+2681*b8*b10*a8-4788*a5^2*b5+21210*b5*b6*a9+49896*a5*a3*a8+6048*a5^2*a7"),
		SMQP("7518*a5*b5*a8+1113*b4*a5*a4-994*b8*b4*a10-6636*b9*b6*a9+1904*b10*b4*a5-1008*b7*a10*a5+1155*b10*b9*a9+4200*b8*a6*a3-1505*b8*a9*a4+8988*a6*a3*a7-9576*b7*a10*a4+13545*b7*a6*a4+5040*b7*a10*a9-168*b6*a4*a7-3724*b6*a8*a7-2268*b6*a5*a3-11788*a5*a8*a7+20286*a4*a3*a8+3024*b6*a4*a3+12978*b8*b6*a5+1890*b6*a8*a3-2898*a8^2*a3+8274*b5*a4*a9-11676*b6*b5*a5+840*b6*a5*a7+3150*b7*a6*a5+2541*b8*a8*a9+1715*b8*a8^2+5481*a6*b7*a8-5068*a4*a7*a8-2884*b8*a8*a4-798*b5*a8^2+2156*a8^2*a7-9086*b8*a10*a3+3556*a7*a9^2+7308*a7*a4^2+861*b5*a4^2-10122*b4*a4^2-12845*b4*a9^2+2835*a8*b9*a9-357*a8*b4*a9+1134*a10*a3^2-4564*b4*a10*a7+2352*a8*a9*a7-1701*a8*b5*a9-6174*a8*a9*a3+14952*a7*a5*a4-5152*a7^2*a6+1512*b9*a10*a7-1988*a10*a3*a7-1008*a8*b7*a10+10584*a9*a4*a3-147*a9^2*b5+504*a9*b9*a4-10556*a9*a4*a7+1288*a9*a5*a7+2562*a9^2*b9+8316*a9*a5*a3+1008*a7^2*a10+4536*b6*b9*a8-17766*a4^2*a3+756*a3^2*a6+11956*b4*a10*a3-2373*b10*b9*a4+10094*b4^2*a10-17766*a3*a5*a4+630*b4*b6*a4-1008*b7*b6*a10+3780*b7*b6*a6+2520*b7*b10*a10-10584*a5^2*a3-7203*b8*b6*a9+3605*b8*b10*a9+8820*b4*b6*a9+10927*a4*b4*a9+1022*b8*b5*a10-504*b8*b9*a10-5544*b7*b10*a6+3108*b9*a4^2+4326*b10*b5*a5-756*b5^2*a6-2982*b4*a6*a3+2800*b8^2*a6-3276*b8*a6*a7-1120*b10*a9*a7-1358*b4*b10*a8+392*b10*a4*a7+2128*b10*a8*a7-5901*b10*b5*a9-3416*b9*a6*a7-9534*b4*b9*a10+11466*b9*a6*a3-2030*b9*b8*a6+504*b5^2*a10-427*b8*b10*a4-6062*b4*b8*a6+5558*b4*b5*a10-12971*b8*a5*a9+1743*b10*b5*a8+23065*b4*a5*a9+3675*b10*b5*a4-4578*b9*a10*a3+6580*b4*b9*a6-1533*b8*b6*a4+1134*b10*a9*a3+1512*b10*a5*a3+581*b8*b6*a8+378*b9*b10*a5-7756*b4*b6*a8+7196*b4*a8*a4+6426*b6*a9*a3-1428*b6*a9*a7-4158*b10*a8*a3-7714*b8*b10*a5-8946*b4*b5*a6-7560*a9^2*a3-7350*b4^2*a6-378*b10*a4*a3-6069*a5*b9*a9+1148*b4*a6*a7-18690*a5^2*b4-8631*b7*a6*a9-1029*b9*b10*a8+3864*b9*b6*a5+2814*b9*b6*a4+539*b4*a5*a8-2667*b8*a5*a4+2016*b8^2*a10-8337*b9*a8*a4-6783*b9*a5*a4-2366*b8*a9^2+2870*b4*b10*a4+14994*b8*a5^2+3255*b8*a4^2+644*b10*a5*a7-3864*b4*b6*a5+980*b4*a8^2-26082*b5*a6*a3+11844*b5*a6*a7+4410*b5*b9*a6-4221*b5*a8*a4+11970*b5*a10*a3-6580*b5*a10*a7-462*b5*b9*a10+3108*b5*a5*a4+2394*b5*b8*a6-2940*a5*b5*a9-231*b9*a8^2+3885*a5*b9*a8+3395*a5*b8*a8+4536*b8*a10*a7-1162*b10*b4*a9+1470*b9*a5^2-2814*b5*b6*a4-2058*b9^2*a6+504*b9^2*a10+42*b5*b6*a8+805*b8*b10*a8-1764*a5^2*b5+8778*b5*b6*a9+12600*a5*a3*a8+2016*a5^2*a7"),
		SMQP("66157*a5*b5*a8+7035*b4*a5*a4+2898*b8*b4*a10+14763*b9*b6*a9-6251*b10*b4*a5+9240*b7*a10*a5-6930*b10*b9*a9+7616*b8*a6*a3+1239*b8*a9*a4-21476*a6*a3*a7+14112*b7*a10*a4-20349*b7*a6*a4-11088*b7*a10*a9-2163*b6*a4*a7-12593*b6*a8*a7+16359*b6*a5*a3+15631*a5*a8*a7-42560*a4*a3*a8+861*b6*a4*a3+20139*b8*b6*a5+4354*b6*a8*a3+8806*a8^2*a3-17542*b5*a4*a9+6825*b6*b5*a5-13440*b6*a5*a7-17787*b7*a6*a5+21*b8*a8*a9+567*b8*a8^2+19971*a6*b7*a8+11872*a4*a7*a8+7140*b8*a8*a4-14336*b5*a8^2+7042*a8^2*a7-2394*b8*a10*a3-7714*a7*a9^2-20790*a7*a4^2+21693*b5*a4^2+3318*b4*a4^2+12257*b4*a9^2-8589*a8*b9*a9-231*a8*b4*a9-9226*a10*a3^2+9100*b4*a10*a7-14448*a8*a9*a7+8253*a8*b5*a9+21840*a8*a9*a3-31395*a7*a5*a4+6160*a7^2*a6-7056*b9*a10*a7+16604*a10*a3*a7-7056*a8*b7*a10-20188*a9*a4*a3+15953*a9^2*b5+5964*a9*b9*a4+20636*a9*a4*a7+11501*a9*a5*a7-4536*a9^2*b9-19747*a9*a5*a3-6048*a7^2*a10-4935*b6*b9*a8+31248*a4^2*a3+17220*a3^2*a6-11508*b4*a10*a3+4032*b10*b9*a4-3458*b4^2*a10+41937*a3*a5*a4-1890*b4*b6*a4-504*b7*b6*a10-11277*b7*b6*a6-1512*b7*b10*a10+15435*a5^2*a3-1281*b8*b6*a9+2247*b8*b10*a9-8127*b4*b6*a9-15295*a4*b4*a9+8946*b8*b5*a10+3024*b8*b9*a10+4851*b7*b10*a6-18354*b9*a4^2-5789*b10*b5*a5+29736*b5^2*a6+21294*b4*a6*a3+3024*b8^2*a6-8512*b8*a6*a7-6251*b10*a9*a7+1106*b4*b10*a8+5341*b10*a4*a7+203*b10*a8*a7+13363*b10*b5*a9+3500*b9*a6*a7+31542*b4*b9*a10-19474*b9*a6*a3+6762*b9*b8*a6-8372*b5^2*a10+861*b8*b10*a4-1498*b4*b8*a6-48118*b4*b5*a10+11613*b8*a5*a9-1540*b10*b5*a8-19852*b4*a5*a9-5537*b10*b5*a4+12306*b9*a10*a3-29344*b4*b9*a6-777*b8*b6*a4+2065*b10*a9*a3-2597*b10*a5*a3-5733*b8*b6*a8+1575*b9*b10*a5+1372*b4*b6*a8-4172*b4*a8*a4-10227*b6*a9*a3+18753*b6*a9*a7+5810*b10*a8*a3-3801*b8*b10*a5+16170*b4*b5*a6+3374*a9^2*a3+1050*b4^2*a6-6629*b10*a4*a3+7140*a5*b9*a9-13412*b4*a6*a7+7707*a5^2*b4+11403*b7*a6*a9-2688*b9*b10*a8-3507*b9*b6*a5+2667*b9*b6*a4+7777*b4*a5*a8-12957*b8*a5*a4-7056*b8^2*a10+21063*b9*a8*a4+8316*b9*a5*a4-1050*b8*a9^2-266*b4*b10*a4-23373*b8*a5^2-2121*b8*a4^2+3388*b10*a5*a7+2205*b4*b6*a5-980*b4*a8^2+20370*b5*a6*a3+784*b5*a6*a7-34762*b5*b9*a6-11375*b5*a8*a4-9618*b5*a10*a3+5572*b5*a10*a7+13902*b5*b9*a10-33390*b5*a5*a4+770*b5*b8*a6-49378*a5*b5*a9-1533*b9*a8^2-18816*a5*b9*a8+8967*a5*b8*a8+4032*b8*a10*a7+4459*b10*b4*a9+7203*b9*a5^2-2562*b5*b6*a4+7014*b9^2*a6-8064*b9^2*a10+22603*b5*b6*a8-189*b8*b10*a8+20601*a5^2*b5-32802*b5*b6*a9-34202*a5*a3*a8-1512*a5^2*a7"),
		SMQP("28035*a6*a3*a4+5229*a6*a3*a5+67452*a6*b5*b6+7119*a6*b5*a5+11088*a6*a3*b6+12180*a6*b4*b6+27993*a6*b4*a5+54390*a6*b5*a4+4536*a4*b6*a5-10416*a6*b4*a4+17640*a4^2*a5+20160*a4*a5^2+11088*b6^2*a5+10080*b6*a5^2+3528*a4^3+7560*a4^2*b6+16128*a5^3-36960*b10*a5^2-78834*b10*b5*a6+17640*b8*a10*a5+13944*b10*b6*a9+15624*b10*b5*a10-12600*b10*b9*a10-16128*b10*a6*a7+22176*b10*a6*a3-4368*b6*a8*a4+15960*b6*b10*a4+7301*a6*a8*a7+14616*a9*b5*a10-19992*b6*b10*a8+11403*a6*a9*a7-22869*a6*a4*a7-3024*b10*a10*a7+10080*b6*a10*a7-8904*b6*a10*a3-17913*a6*b9*a4-25032*a6*a8*a3-5992*b10*a9*a8+5544*b10*a9*a4-12768*b6*a6*a7-2016*b10*a9^2+24696*a6*b7*a10+26313*a6*b9*a9+13993*a6*b9*a8-1680*b6*a9*a8+26768*b10*a8*a4+38514*b10*b9*a6-28392*b5*b6*a10-46263*a6*b4*a9+1064*a8^3+15484*a6*b8*a8+20832*b6*a9^2-3528*b9*a10*a8-10528*b6*a8^2-3920*a8*a9^2-17136*a8*a4^2+2912*a8*b4*a10+2408*a8*b5*a10+24864*a9^2*a5-5040*a8*a10*a7+11984*a8^2*a4+13720*a8*a10*a3-2331*a6*a9*a3-4032*a9^3-8064*a9*a4^2+840*a9*a8^2+11704*a9*a8*a4+10080*a9*a10*a3+33768*a9*b4*a10-7056*a10*a9*a7+8064*a10*a4*a7-31248*a10*b5*a4+8064*a10*b9*a4-11592*a10*a4*a3-1512*a10*b8*a8-6048*a10^2*b7-18144*a10*b9*a9-2688*b6*b10*a5+12306*b8*a6*a4-9912*a9*b6*a4+10836*a8*b4*a6-44128*a5*b10*a8-22680*b10*a4^2-1064*b10*a8^2+9072*b10^2*a5-17472*a5*b6*a9-3528*a4*b8*a10-6552*a4*b4*a10+4536*b9*a10*a5+2016*a9^2*a4+19488*a5*b10*a4+1008*b10^2*a4+7896*b6^2*a8-20664*b6^2*a9+5040*b8*a10*a9-29022*b8*a6*a9+54768*a5*b10*a9-3024*b10^2*a9+5040*a5*a10*a7-18354*b8*b10*a6-10752*a5^2*a8-10143*a5*b9*a6-38556*b9*b6*a6-17304*b4*a10*a5-10528*a5*a8^2-34062*a6*b5*a9+5439*a6*b5*a8-30177*a6^2*b7+5880*b4*b10*a6+10080*b10*a10*a3-11760*a5*b5*a10-6384*a5*a10*a3-19488*a5*a9*a4-16674*b8*b6*a6+4144*a5*a9*a8-2520*b8*b10*a10-29127*b8*a6*a5-16296*b4*b6*a10+10584*b8*b6*a10+22176*b4*b10*a10-3864*a5*a8*a4-12096*a5*b6*a8-29568*a5^2*a9+27216*b9*b6*a10-15288*a5*a6*a7"),
		SMQP("1575*b9^2*a9+3297*b9^2*a4+3045*b9^2*a8+1512*b7*b4*a10-5376*b7*b10*a9+17304*b7*a10*a3+1512*b7*b4*a6+1008*b7*b10*a4+44058*a9*a3*a7+12726*b4*a4*a3+18018*b8*a8*a3-2079*b8*b9*a4+1575*b8*a4*a3+5019*b8*a9*a7-24969*b7*a6*a3-12978*b4*a8*a3-2688*b4*a8*a7+8568*b5^2*a5+12600*b4*b5*a4+10353*a7*b7*a6+46095*a7*a8*a3+14070*a7*b5*a9-4725*a7*b5*a5-5376*b8*b7*a10-16737*a7*b4*a5-10584*a7*b7*a10+30135*a7*b4*a9+14847*a7*b8*a5-9975*a7*b8*a8+987*b8*a9*a3-44730*b8*a5*a3-1113*b8*b4*a9-2625*b8*b4*a5-4158*b8*b5*a4+1008*a7^2*a5-95088*a7*a4*a3+25263*a7^2*a4+1008*b4*b8*a8-3276*b4*b9*a8-7560*b4*a4*a7+16884*b4*b5*a8-12201*a7^2*a9-903*a7*b8*a4-23121*a9*a3^2-9282*b7*a4^2+7602*b7*a4*a8-3381*b8^2*a5+37989*a3*b4*a5+92169*a3^2*a4+13671*a3^2*a5-45213*a3*b4*a9+1470*b7*a5*a9-13293*b9^2*a5+18326*b7*a5*a8-14994*b5*a4*a3+7392*b7*b5*a10+7518*b5*a4*a7+504*b5^2*a4-924*b5^2*a8-13335*a5*a3*a7-4368*b5^2*a9-8967*a7^2*a8-14364*b4*b9*a4-69552*a3^2*a8-7245*b5*b9*a5-1533*b5*b8*a5-15120*b5*b4*a5-5040*b4^2*a9-3675*b5*b8*a8+4473*b7*b8*a6-7056*b7*a9*a8-16548*b7*b6*a5-4704*b7*b5*a6+17178*b7*a9*a4-798*b7*b6*a4+9702*b7*b6*a9-1680*b7*b10*a8-3339*b5*b9*a8+5334*b7*b6*a8-7056*b7*a9^2-12684*b7*a5^2+6384*b7*b10*a5+10206*b5*b9*a9-12642*a5*b7*a4-3654*b8^2*a9+3822*b8^2*a8-2730*b9*b5*a4-42378*b5*a9*a3+5922*b5*b8*a9-14952*b5*b4*a9+3780*b8^2*a4+4998*b9*a9*a3+15414*b9*a5*a3-6867*b9*b8*a9-4977*b9*b8*a8-777*b9*b4*a5+21798*b9*b8*a5+18123*b9*b4*a9+28707*b5*a8*a3-10563*b5*a8*a7+5733*b5*a5*a3-1197*b9*a5*a7+18207*b9*a8*a3-21462*b9*a4*a3-5502*b9*a9*a7+9996*b9*a4*a7-7182*b9*a8*a7-1680*b9*b7*a10+4641*b9*b7*a6"),
		SMQP("64078*a5*b5*a8+1449*b4*a5*a4-2450*b8*b4*a10+22764*b9*b6*a9-1232*b10*b4*a5+11088*b7*a10*a5-10857*b10*b9*a9-8232*b8*a6*a3-2989*b8*a9*a4-43092*a6*a3*a7+24696*b7*a10*a4-27279*b7*a6*a4-23184*b7*a10*a9-5544*b6*a4*a7-16044*b6*a8*a7+24948*b6*a5*a3+68628*a5*a8*a7-29106*a4*a3*a8+11088*b6*a4*a3+28266*b8*b6*a5-7686*b6*a8*a3+630*a8^2*a3-23422*b5*a4*a9+8484*b6*b5*a5-22680*b6*a5*a7-20370*b7*a6*a5-4095*b8*a8*a9+2863*b8*a8^2+24381*a6*b7*a8-6720*a4*a7*a8+3052*b8*a8*a4-14630*b5*a8^2+12180*a8^2*a7-5614*b8*a10*a3-14952*a7*a9^2-15624*a7*a4^2+23457*b5*a4^2+3318*b4*a4^2+21203*b4*a9^2-14637*a8*b9*a9-9681*a8*b4*a9-9450*a10*a3^2+588*b4*a10*a7-15372*a8*a9*a7-1701*a8*b5*a9+41202*a8*a9*a3-71064*a7*a5*a4+15456*a7^2*a6-4536*b9*a10*a7+19068*a10*a3*a7-5040*a8*b7*a10-38304*a9*a4*a3+24689*a9^2*b5+14112*a9*b9*a4+36708*a9*a4*a7-28056*a9*a5*a7-6090*a9^2*b9+15624*a9*a5*a3-9072*a7^2*a10-14784*b6*b9*a8+23310*a4^2*a3+36036*a3^2*a6+1820*b4*a10*a3+7287*b10*b9*a4-3458*b4^2*a10+75222*a3*a5*a4-1890*b4*b6*a4-11088*b7*b6*a10-25956*b7*b6*a6+504*b7*b10*a10+1512*a5^2*a3-4263*b8*b6*a9+3409*b8*b10*a9-7812*b4*b6*a9-16177*a4*b4*a9+13678*b8*b5*a10+19656*b8*b9*a10+12600*b7*b10*a6-34104*b9*a4^2-4802*b10*b5*a5+22932*b5^2*a6+33978*b4*a6*a3+2912*b8^2*a6-2324*b8*a6*a7-14784*b10*a9*a7+1106*b4*b10*a8+4872*b10*a4*a7+5712*b10*a8*a7+26551*b10*b5*a9-5040*b9*a6*a7+61698*b4*b9*a10-27846*b9*a6*a3+12026*b9*b8*a6-7280*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-74242*b4*b5*a10+23933*b8*a5*a9-4837*b10*b5*a8-29407*b4*a5*a9-11753*b10*b5*a4+18942*b9*a10*a3-34636*b4*b9*a6-3129*b8*b6*a4+2646*b10*a9*a3-4536*b10*a5*a3-6503*b8*b6*a8+2226*b9*b10*a5+2380*b4*b6*a8-4172*b4*a8*a4-7182*b6*a9*a3+31500*b6*a9*a7+4410*b10*a8*a3-8666*b8*b10*a5-13818*b4*b5*a6+2772*a9^2*a3+1050*b4^2*a6-4914*b10*a4*a3-1701*a5*b9*a9-20132*b4*a6*a7+26838*a5^2*b4+14553*b7*a6*a9-4809*b9*b10*a8+4536*b9*b6*a5+9030*b9*b6*a4+3703*b4*a5*a8-15771*b8*a5*a4-18144*b8^2*a10+33327*b9*a8*a4+28497*b9*a5*a4+2786*b8*a9^2-266*b4*b10*a4-51030*b8*a5^2+3171*b8*a4^2+7140*b10*a5*a7-8904*b4*b6*a5+28*b4*a8^2+55902*b5*a6*a3-7700*b5*a6*a7-28294*b5*b9*a6-3605*b5*a8*a4+98*b5*a10*a3+4620*b5*a10*a7+10290*b5*b9*a10-54264*b5*a5*a4+2450*b5*b8*a6-28504*a5*b5*a9+2037*b9*a8^2-48111*a5*b9*a8+26383*a5*b8*a8+4536*b8*a10*a7+3766*b10*b4*a9+17766*b9*a5^2+5754*b5*b6*a4+5838*b9^2*a6-19656*b9^2*a10+31066*b5*b6*a8+2681*b8*b10*a8+37548*a5^2*b5-55734*b5*b6*a9-85680*a5*a3*a8+18144*a5^2*a7"),
		SMQP("-53634*a5*b5*a8-45843*b4*a5*a4-11466*b8*b4*a10-36372*b9*b6*a9+14448*b10*b4*a5-35280*b7*a10*a5+8323*b10*b9*a9+62552*b8*a6*a3-3969*b8*a9*a4+57820*a6*a3*a7-31752*b7*a10*a4+52101*b7*a6*a4+17136*b7*a10*a9-9576*b6*a4*a7+38388*b6*a8*a7-17724*b6*a5*a3+26292*a5*a8*a7+143038*a4*a3*a8-8400*b6*a4*a3+15330*b8*b6*a5+37618*b6*a8*a3-49154*a8^2*a3+27594*b5*a4*a9-36204*b6*b5*a5+17640*b6*a5*a7+56070*b7*a6*a5+10605*b8*a8*a9+9723*b8*a8^2-41895*a6*b7*a8-100968*a4*a7*a8-35700*b8*a8*a4+21882*b5*a8^2+8652*a8^2*a7-28854*b8*a10*a3+13524*a7*a9^2+78120*a7*a4^2-2331*b5*a4^2-28098*b4*a4^2-52773*b4*a9^2+7203*a8*b9*a9-16653*a8*b4*a9+350*a10*a3^2-61908*b4*a10*a7+24192*a8*a9*a7-16989*a8*b5*a9-35742*a8*a9*a3+13608*a7*a5*a4-6720*a7^2*a6+23688*b9*a10*a7-37380*a10*a3*a7+19152*a8*b7*a10+62972*a9*a4*a3-15435*a9^2*b5+15260*a9*b9*a4-16464*a9*a4*a7-73416*a9*a5*a7+9842*a9^2*b9+2492*a9*a5*a3+11088*a7^2*a10+112*b6*b9*a8-104202*a4^2*a3-26124*a3^2*a6+24556*b4*a10*a3-749*b10*b9*a4+34566*b4^2*a10-141330*a3*a5*a4+16758*b4*b6*a4-21168*b7*b6*a10+29484*b7*b6*a6+10584*b7*b10*a10-15624*a5^2*a3-27027*b8*b6*a9+11781*b8*b10*a9+26460*b4*b6*a9+50043*a4*b4*a9-13482*b8*b5*a10+27720*b8*b9*a10-13608*b7*b10*a6+20664*b9*a4^2+35238*b10*b5*a5-49644*b5^2*a6-36750*b4*a6*a3+4032*b8^2*a6+27132*b8*a6*a7+4368*b10*a9*a7-12054*b4*b10*a8-6216*b10*a4*a7+6720*b10*a8*a7-30429*b10*b5*a9-43232*b9*a6*a7-32774*b4*b9*a10+52178*b9*a6*a3-41230*b9*b8*a6+20160*b5^2*a10-10899*b8*b10*a4+1722*b4*b8*a6+39606*b4*b5*a10-33411*b8*a5*a9-5817*b10*b5*a8+113505*b4*a5*a9+14931*b10*b5*a4-25690*b9*a10*a3+32004*b4*b9*a6+22995*b8*b6*a4+1246*b10*a9*a3+616*b10*a5*a3+2877*b8*b6*a8-13622*b9*b10*a5-10836*b4*b6*a8+54516*b4*a8*a4+12138*b6*a9*a3-45108*b6*a9*a7-10318*b10*a8*a3-5586*b8*b10*a5-40194*b4*b5*a6-22456*a9^2*a3-15246*b4^2*a6-266*b10*a4*a3-11053*a5*b9*a9+27132*b4*a6*a7-70770*a5^2*b4-22239*b7*a6*a9+5747*b9*b10*a8+26040*b9*b6*a5+8862*b9*b6*a4-9429*b4*a5*a8+86793*b8*a5*a4+2016*b8^2*a10-35693*b9*a8*a4-54747*b9*a5*a4-126*b8*a9^2+2814*b4*b10*a4+47250*b8*a5^2-945*b8*a4^2-22092*b10*a5*a7-12936*b4*b6*a5+3612*b4*a8^2-119322*b5*a6*a3+47628*b5*a6*a7+54306*b5*b9*a6-16233*b5*a8*a4+83818*b5*a10*a3-27636*b5*a10*a7-25718*b5*b9*a10+36456*b5*a5*a4+7434*b5*b8*a6+49812*a5*b5*a9+6097*b9*a8^2+15757*a5*b9*a8-44877*a5*b8*a8+504*b8*a10*a7+2814*b10*b4*a9-3234*b9*a5^2-12222*b5*b6*a4+8918*b9^2*a6-7560*b9^2*a10-28686*b5*b6*a8-1155*b8*b10*a8-12852*a5^2*b5+66402*b5*b6*a9+145096*a5*a3*a8+42336*a5^2*a7"),
		SMQP("10563*a6*a3*a4+2877*a6*a3*a5+28140*a6*b5*b6-1953*a6*b5*a5+2100*a6*b4*b6+6909*a6*b4*a5+10794*a6*b5*a4+14784*a4*b6*a5+2016*a6*b4*a4-2688*a4^2*a5+15960*a4*a5^2+336*b6*a5^2+5544*a4^3-1680*a4^2*b6-2016*a5^3+2240*b10*a5^2-16086*b10*b5*a6+6888*b8*a10*a5+7112*b10*b5*a10-1400*b10*b9*a10+1344*b10*a6*a7+4760*b6*a8*a4+6741*a6*a8*a7-3248*a9*b5*a10+4704*b6*b10*a8-2373*a6*a9*a7-5845*a6*a4*a7-1008*b10*a10*a7+3360*b6*a10*a7+2191*a6*b9*a4-13944*a6*a8*a3+504*b10*a9*a8+3248*b10*a9*a4-5040*b6*a6*a7-1512*b10*a9^2+9912*a6*b7*a10+2793*a6*b9*a9+17325*a6*b9*a8-7560*b6*a9*a8+1624*b10*a8*a4-1386*b10*b9*a6-10864*b5*b6*a10-13251*a6*b4*a9-1008*a8^3-672*a6*b8*a8+9576*b6*a9^2-3304*b9*a10*a8-1008*b6*a8^2-3528*a8*a9^2-10360*a8*a4^2+5768*a8*b4*a10+3472*a8*b5*a10+19600*a9^2*a5-2352*a8*a10*a7+6776*a8^2*a4+5712*a8*a10*a3-1323*a6*a9*a3+1008*a9^3-4088*a9*a4^2+1512*a9*a8^2+9408*a9*a8*a4+6216*a9*a10*a3+14000*a9*b4*a10-1680*a10*a9*a7+3360*a10*a4*a7-6720*a10*b5*a4-3136*a10*b9*a4-6776*a10*a4*a3-3080*a10*b8*a8-2016*a10^2*b7-1120*a10*b9*a9-1610*b8*a6*a4-6216*a9*b6*a4+5124*a8*b4*a6+6664*a5*b10*a8-1288*b10*a4^2-2016*b10*a8^2+3864*a5*b6*a9+3640*a4*b8*a10-5768*a4*b4*a10-392*b9*a10*a5-4592*a9^2*a4-5096*a5*b10*a4-5040*b6^2*a8+1456*b8*a10*a9-14742*b8*a6*a9+2744*a5*b10*a9+336*a5*a10*a7+8778*b8*b10*a6-23632*a5^2*a8-203*a5*b9*a6-7980*b9*b6*a6+3024*b4*a10*a5-4648*a5*a8^2-5586*a6*b5*a9-3969*a6*b5*a8-10269*a6^2*b7+1848*b4*b10*a6-1008*b10^2*a8-3976*a5*b5*a10+2632*a5*a10*a3+3416*a5*a9*a4-18102*b8*b6*a6-13272*a5*a9*a8-3640*b8*b10*a10-17171*b8*a6*a5-1232*b4*b6*a10+6776*b8*b6*a10-224*b4*b10*a10+12040*a5*a8*a4-24472*a5*b6*a8-3584*a5^2*a9+7504*b9*b6*a10-5992*a5*a6*a7"),
		SMQP("-42189*a6*a3*a4+13293*a6*a3*a5-24276*a6*b5*b6+567*a6*b5*a5-20832*a6*a3*b6-45612*a6*b4*b6-27699*a6*b4*a5-53718*a6*b5*a4+17472*a4*b6*a5-67704*a6*b4*a4-37632*a4^2*a5-1344*a4*b6^2-4032*a4*a5^2-19488*b6^2*a5-14112*b6*a5^2-4032*a4^3-11424*a4^2*b6-16128*a5^3+10080*b10*a5^2+20202*b10*b5*a6-12348*b8*a10*a5-27608*b10*b6*a9-10444*b10*b5*a10+3444*b10*b9*a10-14896*b10*a6*a7+29904*b10*a6*a3+19712*b6*a8*a4+14224*b6*b10*a4+2709*a6*a8*a7-20692*a9*b5*a10+6664*b6*b10*a8-25669*a6*a9*a7+38395*a6*a4*a7+5040*b10*a10*a7-6048*b6*a10*a7+14672*b6*a10*a3+50519*a6*b9*a4-12600*a6*a8*a3+4144*b10*a9*a8+8008*b10*a9*a4+6832*b6*a6*a7-2856*b10*a9^2-1260*a6*b7*a10-33047*a6*b9*a9+15477*a6*b9*a8-17080*b6*a9*a8-8792*b10*a8*a4-17066*b10*b9*a6+12376*b5*b6*a10+51093*a6*b4*a9-4200*a6*b8*a8+1400*b6*a9^2+1596*b9*a10*a8-1232*b6*a8^2+2856*a8*a9^2+13552*a8*a4^2+10416*a8*b4*a10-10920*a8*b5*a10+18592*a9^2*a5+2016*a8*a10*a7-6272*a8^2*a4+23541*a6*a9*a3+3472*a9^3+1064*a9*a4^2+3752*a9*a8^2-280*a9*a8*a4-4088*a9*a10*a3-39676*a9*b4*a10+6048*a10*a9*a7-13104*a10*a4*a7+37156*a10*b5*a4-20328*a10*b9*a4+15680*a10*a4*a3+4452*a10*b8*a8+16968*a10*b9*a9+30352*b6*b10*a5+14462*b8*a6*a4+12656*a9*b6*a4-10668*a8*b4*a6-1120*a5*b10*a8-2744*b10*a4^2+1736*b10*a8^2-6832*b10^2*a5+30128*a5*b6*a9-24780*a4*b8*a10+50176*a4*b4*a10-11172*b9*a10*a5-16128*a9^2*a4-5432*a5*b10*a4-3304*b10^2*a4-15008*b6^2*a8+21840*b6^2*a9+5208*b8*a10*a9-5222*b8*a6*a9+2688*a5*b10*a9+6776*b10^2*a9-7056*a5*a10*a7-14462*b8*b10*a6-24192*a5^2*a8+16149*a5*b9*a6+23156*b9*b6*a6+61404*b4*a10*a5-14112*a5*a8^2+12054*a6*b5*a9+6447*a6*b5*a8+6363*a6^2*b7+46032*b4*b10*a6-1736*b10^2*a8-10136*b10*a10*a3+11592*a5*b5*a10+16128*a5*a10*a3+28616*a5*a9*a4-3598*b8*b6*a6-34384*a5*a9*a8+9660*b8*b10*a10-26523*b8*a6*a5+22792*b4*b6*a10-4452*b8*b6*a10-16408*b4*b10*a10+60592*a5*a8*a4-37184*a5*b6*a8+14112*a5^2*a9-5880*b9*b6*a10-2520*a5*a6*a7"),
		SMQP("-1575*b9^2*a9+1071*b9^2*a4+315*b9^2*a8-3248*b7*b4*a10+868*b7*b10*a9-5068*b7*a10*a3-525*b7*b4*a6-644*b7*b10*a4-3402*a9*a3*a7-5418*b4*a4*a3-4536*b8*a8*a3-2457*b8*b9*a4+756*b8*a4*a3-63*b8*a9*a7+2856*b7*a6*a3+5733*b4*a8*a3+315*b4*a8*a7-1512*b5^2*a5+2016*b4*b5*a4+994*a7*b7*a6-15876*a7*a8*a3-3843*a7*b5*a9-504*a7*b5*a5+1512*b8*b7*a10+2394*a7*b4*a5+2268*a7*b7*a10-7245*a7*b4*a9-6552*a7*b8*a5+2583*a7*b8*a8-1512*b8*a9*a3+12852*b8*a5*a3+378*b8*b4*a9+3213*b8*b4*a5+441*b8*b5*a4+21546*a7*a4*a3-5670*a7^2*a4-882*b4*b8*a8-1449*b4*b9*a8+2961*b4*a4*a7-1827*b4*b5*a8+1134*a7^2*a9+315*a7*b8*a4+2268*a9*a3^2+4473*b7*a4^2-6839*b7*a4*a8-1764*b8^2*a5-6615*a3*b4*a5-20412*a3^2*a4-2268*a3^2*a5+9198*a3*b4*a9-7399*b7*a5*a9+630*b9^2*a5+3535*b7*a5*a8+2079*b5*a4*a3-224*b7*b5*a10-1197*b5*a4*a7-756*b5^2*a4+1134*a5*a3*a7+1512*b5^2*a9+3402*a7^2*a8+1827*b4*b9*a4+18144*a3^2*a8+630*b5*b9*a5+2394*b5*b8*a5+945*b5*b4*a5+1701*b4^2*a9-315*b5*b8*a8+1064*b7*b8*a6+189*b4^2*a5+1932*b7*a9*a8+2058*b7*b6*a5+3780*b7*b5*a6-1729*b7*a9*a4-1029*b7*b6*a4-2163*b7*b6*a9+644*b7*b10*a8+315*b5*b9*a8+1771*b7*b6*a8+224*b7*a9^2+5166*b7*a5^2-2156*b7*b10*a5-1575*b5*b9*a9+868*b7*a8^2-399*a5*b7*a4-126*b8^2*a9+630*b8^2*a8+1071*b9*b5*a4+7749*b5*a9*a3-1449*b5*b8*a9+252*b5*b4*a9+630*b8^2*a4+693*b9*a9*a3-6174*b9*a5*a3+3213*b9*b8*a9-945*b9*b8*a8-819*b9*b4*a5-378*b9*b8*a5-1575*b9*b4*a9-7749*b5*a8*a3+3717*b5*a8*a7+1890*b5*a5*a3+2772*b9*a5*a7-1953*b9*a8*a3+3339*b9*a4*a3+441*b9*a9*a7-2205*b9*a4*a7+819*b9*a8*a7-3220*b9*b7*a6"),
		SMQP("69426*a5*b5*a8+6237*b4*a5*a4+5614*b8*b4*a10+15036*b9*b6*a9-5684*b10*b4*a5+11648*b7*a10*a5-7021*b10*b9*a9+12824*b8*a6*a3+1715*b8*a9*a4-34580*a6*a3*a7+18312*b7*a10*a4-28035*b7*a6*a4-12432*b7*a10*a9-6216*b6*a4*a7-11788*b6*a8*a7+17556*b6*a5*a3+27524*a5*a8*a7-51086*a4*a3*a8+7728*b6*a4*a3+19194*b8*b6*a5-2534*b6*a8*a3+14518*a8^2*a3-32214*b5*a4*a9+9828*b6*b5*a5-13272*b6*a5*a7-21042*b7*a6*a5+1561*b8*a8*a9-49*b8*a8^2+28161*a6*b7*a8+14672*a4*a7*a8+10220*b8*a8*a4-15890*b5*a8^2+6608*a8^2*a7-2926*b8*a10*a3-9884*a7*a9^2-22932*a7*a4^2+26481*b5*a4^2+3318*b4*a4^2+10535*b4*a9^2-9709*a8*b9*a9+3899*a8*b4*a9-11914*a10*a3^2+9548*b4*a10*a7-16212*a8*a9*a7+17983*a8*b5*a9+17262*a8*a9*a3-42504*a7*a5*a4+10976*a7^2*a6-6888*b9*a10*a7+20860*a10*a3*a7-10416*a8*b7*a10-14056*a9*a4*a3+22953*a9^2*b5+7336*a9*b9*a4+23716*a9*a4*a7+19768*a9*a5*a7-5222*a9^2*b9-27412*a9*a5*a3-7728*a7^2*a10-1680*b6*b9*a8+36162*a4^2*a3+25284*a3^2*a6-11508*b4*a10*a3+4403*b10*b9*a4-3458*b4^2*a10+58002*a3*a5*a4-1890*b4*b6*a4-672*b7*b6*a10-11844*b7*b6*a6-2184*b7*b10*a10+16632*a5^2*a3-735*b8*b6*a9+2065*b8*b10*a9-10836*b4*b6*a9-17437*a4*b4*a9+15358*b8*b5*a10+5880*b8*b9*a10+6300*b7*b10*a6-19236*b9*a4^2-8022*b10*b5*a5+39060*b5^2*a6+22554*b4*a6*a3+6944*b8^2*a6-8596*b8*a6*a7-7588*b10*a9*a7+882*b4*b10*a8+6692*b10*a4*a7-140*b10*a8*a7+14847*b10*b5*a9+4368*b9*a6*a7+37058*b4*b9*a10-18214*b9*a6*a3+5754*b9*b8*a6-11760*b5^2*a10+1225*b8*b10*a4-2590*b4*b8*a6-61250*b4*b5*a10+7217*b8*a5*a9-721*b10*b5*a8-21931*b4*a5*a9-5817*b10*b5*a4+10318*b9*a10*a3-39340*b4*b9*a6+903*b8*b6*a4+2842*b10*a9*a3-2996*b10*a5*a3-5775*b8*b6*a8-1666*b9*b10*a5+1148*b4*b6*a8-3948*b4*a8*a4-12558*b6*a9*a3+22764*b6*a9*a7+7994*b10*a8*a3-3038*b8*b10*a5+31878*b4*b5*a6+392*a9^2*a3+1050*b4^2*a6-8918*b10*a4*a3+5887*a5*b9*a9-14756*b4*a6*a7+12726*a5^2*b4+13797*b7*a6*a9-3493*b9*b10*a8+1848*b9*b6*a5+294*b9*b6*a4+11123*b4*a5*a8-12999*b8*a5*a4-10752*b8^2*a10+23191*b9*a8*a4+5397*b9*a5*a4-2422*b8*a9^2-266*b4*b10*a4-34902*b8*a5^2-4893*b8*a4^2+3332*b10*a5*a7+3192*b4*b6*a5-1428*b4*a8^2+16254*b5*a6*a3+1260*b5*a6*a7-41958*b5*b9*a6-15113*b5*a8*a4-8078*b5*a10*a3+6524*b5*a10*a7+18242*b5*b9*a10-35028*b5*a5*a4-13230*b5*b8*a6-58380*a5*b5*a9-3983*b9*a8^2-12243*a5*b9*a8+10171*a5*b8*a8+2856*b8*a10*a7+5362*b10*b4*a9+11718*b9*a5^2-4998*b5*b6*a4+9646*b9^2*a6-10584*b9^2*a10+25802*b5*b6*a8-1127*b8*b10*a8+20412*a5^2*b5-37254*b5*b6*a9-52052*a5*a3*a8-2016*a5^2*a7"),
		SMQP("7203*a6*a3*a4+2877*a6*a3*a5+27132*a6*b5*b6-4977*a6*b5*a5+3108*a6*b4*b6+5397*a6*b4*a5+12306*a6*b5*a4+14112*a4*b6*a5+3360*a6*b4*a4-7392*a4^2*a5+19992*a4*a5^2+336*b6*a5^2+9576*a4^3-4368*a4^2*b6-2016*a5^3+2240*b10*a5^2-15582*b10*b5*a6+1638*b8*a10*a5+6720*b10*b5*a10-336*b10*b9*a10+1344*b10*a6*a7+10024*b6*a8*a4+3493*a6*a8*a7-8120*a9*b5*a10+4368*b6*b10*a8-3605*a6*a9*a7-4613*a6*a4*a7-1008*b10*a10*a7+3024*b6*a10*a7+6671*a6*b9*a4-10920*a6*a8*a3+168*b10*a9*a8+3920*b10*a9*a4-4032*b6*a6*a7-392*b10*a9^2+12306*a6*b7*a10-5047*a6*b9*a9+28469*a6*b9*a8-4424*b6*a9*a8+3192*b10*a8*a4-5250*b10*b9*a6-9744*b5*b6*a10-16107*a6*b4*a9-224*a8^3-3640*a6*b8*a8+7224*b6*a9^2-6594*b9*a10*a8-5936*b6*a8^2-5320*a8*a9^2-19880*a8*a4^2+15064*a8*b4*a10+6958*a8*b5*a10+15792*a9^2*a5-2814*a8*a10*a7+12824*a8^2*a4+7280*a8*a10*a3-1995*a6*a9*a3+3248*a9^3-3976*a9*a4^2+280*a9*a8^2+9072*a9*a8*a4+5642*a9*a10*a3+9814*a9*b4*a10-1218*a10*a9*a7+4830*a10*a4*a7-8960*a10*b5*a4-6930*a10*b9*a4-7882*a10*a4*a3-2940*a10*b8*a8-2016*a10^2*b7+546*a10*b9*a9+1022*b8*a6*a4-9240*a9*b6*a4+2772*a8*b4*a6+6776*a5*b10*a8-1736*b10*a4^2-2800*b10*a8^2+4536*a5*b6*a9+5376*a4*b8*a10-10808*a4*b4*a10+3486*b9*a10*a5-4256*a9^2*a4-4536*a5*b10*a4-4032*b6^2*a8+2940*b8*a10*a9-14854*b8*a6*a9+2856*a5*b10*a9+336*a5*a10*a7+9282*b8*b10*a6-29008*a5^2*a8-6419*a5*b9*a6-7980*b9*b6*a6+7238*b4*a10*a5-11928*a5*a8^2+3990*a6*b5*a9-8673*a6*b5*a8-12789*a6^2*b7+1848*b4*b10*a6-1008*b10^2*a8-2170*a5*b5*a10+3850*a5*a10*a3-840*a5*a9*a4-20622*b8*b6*a6-6216*a5*a9*a8-4704*b8*b10*a10-17003*b8*a6*a5-1344*b4*b6*a10+8400*b8*b6*a10+672*b4*b10*a10+22008*a5*a8*a4-23464*a5*b6*a8+448*a5^2*a9+9072*b9*b6*a10-4648*a5*a6*a7"),
		SMQP("-2562*a5*b5*a8+861*b4*a5*a4-2450*b8*b4*a10+924*b9*b6*a9-896*b10*b4*a5+336*b7*a10*a5-637*b10*b9*a9-5992*b8*a6*a3-49*b8*a9*a4-3668*a6*a3*a7+3192*b7*a10*a4-4683*b7*a6*a4-1680*b7*a10*a9-168*b6*a4*a7+1820*b6*a8*a7+420*b6*a5*a3+5012*a5*a8*a7-5138*a4*a3*a8-336*b6*a4*a3-4578*b8*b6*a5-2366*b6*a8*a3+910*a8^2*a3-2730*b5*a4*a9+3948*b6*b5*a5-168*b6*a5*a7-798*b7*a6*a5-259*b8*a8*a9-777*b8*a8^2-2499*a6*b7*a8+1148*a4*a7*a8-2436*b8*a8*a4+266*b5*a8^2-700*a8^2*a7+1106*b8*a10*a3-1148*a7*a9^2-2268*a7*a4^2-231*b5*a4^2+3318*b4*a4^2+4571*b4*a9^2-245*a8*b9*a9+1015*a8*b4*a9+14*a10*a3^2+980*b4*a10*a7-840*a8*a9*a7+623*a8*b5*a9+2226*a8*a9*a3-5712*a7*a5*a4+1904*a7^2*a6-504*b9*a10*a7+532*a10*a3*a7+336*a8*b7*a10-4648*a9*a4*a3+21*a9^2*b5+1456*a9*b9*a4+3892*a9*a4*a7-1400*a9*a5*a7-1358*a9^2*b9+140*a9*a5*a3-336*a7^2*a10+1288*b6*b9*a8+5418*a4^2*a3+84*a3^2*a6-2268*b4*a10*a3+623*b10*b9*a4-3458*b4^2*a10+8106*a3*a5*a4-1890*b4*b6*a4+336*b7*b6*a10-1260*b7*b6*a6-840*b7*b10*a10+1512*a5^2*a3+3045*b8*b6*a9-1379*b8*b10*a9-1260*b4*b6*a9-5341*a4*b4*a9+1582*b8*b5*a10+168*b8*b9*a10+1848*b7*b10*a6+420*b9*a4^2-1470*b10*b5*a5+252*b5^2*a6-1638*b4*a6*a3+3920*b8^2*a6+5012*b8*a6*a7+392*b10*a9*a7-14*b4*b10*a8-112*b10*a4*a7-728*b10*a8*a7+1995*b10*b5*a9+2688*b9*a6*a7-742*b4*b9*a10-3766*b9*a6*a3-126*b9*b8*a6-168*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-1274*b4*b5*a10-4795*b8*a5*a9-581*b10*b5*a8-10339*b4*a5*a9-1281*b10*b5*a4+70*b9*a10*a3-5572*b4*b9*a6-3129*b8*b6*a4-434*b10*a9*a3-56*b10*a5*a3+2849*b8*b6*a8-1330*b9*b10*a5+2268*b4*b6*a8-3052*b4*a8*a4-966*b6*a9*a3+84*b6*a9*a7+1442*b10*a8*a3+826*b8*b10*a5+7182*b4*b5*a6+2408*a9^2*a3+1050*b4^2*a6+70*b10*a4*a3-4613*a5*b9*a9+868*b4*a6*a7+8778*a5^2*b4+2625*b7*a6*a9+567*b9*b10*a8-336*b9*b6*a5-1722*b9*b6*a4+1071*b4*a5*a8-6363*b8*a5*a4-672*b8^2*a10-525*b9*a8*a4-3003*b9*a5*a4+434*b8*a9^2-266*b4*b10*a4+630*b8*a5^2+3171*b8*a4^2-364*b10*a5*a7-532*b4*a8^2+8190*b5*a6*a3-3780*b5*a6*a7-126*b5*b9*a6+1295*b5*a8*a4-4382*b5*a10*a3+2324*b5*a10*a7+1610*b5*b9*a10-672*b5*a5*a4+3402*b5*b8*a6+840*a5*b5*a9-147*b9*a8^2+6853*a5*b9*a8+10507*a5*b8*a8-1512*b8*a10*a7+742*b10*b4*a9+3738*b9*a5^2+1050*b5*b6*a4-1778*b9^2*a6-168*b9^2*a10-14*b5*b6*a8-63*b8*b10*a8+756*a5^2*b5-2982*b5*b6*a9-7448*a5*a3*a8"),
		SMQP("1575*b9^2*a9-1071*b9^2*a4-63*b9^2*a8+56*b7*b4*a10-1036*b7*b10*a9+3220*b7*a10*a3+777*b7*b4*a6+476*b7*b10*a4+2772*a9*a3*a7+3654*b4*a4*a3+6300*b8*a8*a3+2520*b8*b9*a4-1197*b8*a4*a3+126*b8*a9*a7-5250*b7*a6*a3-3969*b4*a8*a3-567*b4*a8*a7+1512*b5^2*a5-1008*b4*b5*a4+2366*a7*b7*a6+16506*a7*a8*a3+3843*a7*b5*a9-756*a7*b5*a5-1512*b8*b7*a10-3402*a7*b4*a5-2268*a7*b7*a10+7497*a7*b4*a9+4788*a7*b8*a5-2394*a7*b8*a8+189*b8*a9*a3-15561*b8*a5*a3+189*b8*b4*a9-882*b8*b4*a5-189*b8*b5*a4-22176*a7*a4*a3+5670*a7^2*a4+378*b4*b8*a8+1197*b4*b9*a8-2709*b4*a4*a7+1071*b4*b5*a8-1134*a7^2*a9-252*a7*b8*a4-6678*a9*a3^2-945*b7*a4^2+3395*b7*a4*a8-441*b8^2*a5+13041*a3*b4*a5+24822*a3^2*a4-2142*a3^2*a5-10080*a3*b4*a9+7*b7*a5*a9-1890*b9^2*a5+6125*b7*a5*a8-4599*b5*a4*a3+2072*b7*b5*a10+1197*b5*a4*a7+756*b5^2*a4+3906*a5*a3*a7-1512*b5^2*a9-3402*a7^2*a8-1575*b4*b9*a4-18144*a3^2*a8-1890*b5*b9*a5-567*b5*b8*a5-5733*b5*b4*a5-1449*b4^2*a9-378*b5*b8*a8+1561*b7*b8*a6+63*b4^2*a5-1428*b7*a9*a8-2310*b7*b6*a5-1260*b7*b5*a6+1393*b7*a9*a4-1995*b7*b6*a4+2667*b7*b6*a9-476*b7*b10*a8-63*b5*b9*a8+1337*b7*b6*a8-560*b7*a9^2-2394*b7*a5^2+476*b7*b10*a5+1575*b5*b9*a9-1036*b7*a8^2-5565*a5*b7*a4+252*b8^2*a9-252*b8^2*a8-1071*b9*b5*a4-2709*b5*a9*a3+1449*b5*b8*a9-252*b5*b4*a9-630*b8^2*a4-63*b9*a9*a3+4284*b9*a5*a3-3276*b9*b8*a9+252*b9*b8*a8+567*b9*b4*a5+3717*b9*b8*a5+1323*b9*b4*a9+6363*b5*a8*a3-3465*b5*a8*a7+1260*b5*a5*a3-1512*b9*a5*a7+1575*b9*a8*a3-3969*b9*a4*a3-441*b9*a9*a7+2205*b9*a4*a7-1071*b9*a8*a7+364*b9*b7*a6"),
		SMQP("-45549*a6*a3*a4+13293*a6*a3*a5-22764*a6*b5*b6+63*a6*b5*a5-20832*a6*a3*b6-47124*a6*b4*b6-23919*a6*b4*a5-54474*a6*b5*a4+20832*a4*b6*a5-75432*a6*b4*a4-64512*a4^2*a5-1344*a4*b6^2+16128*a4*a5^2-19488*b6^2*a5-14112*b6*a5^2+6048*a4^3-18144*a4^2*b6-16128*a5^3+10080*b10*a5^2+19950*b10*b5*a6-13608*b8*a10*a5-27608*b10*b6*a9-10528*b10*b5*a10+3528*b10*b9*a10-14896*b10*a6*a7+29904*b10*a6*a3+37072*b6*a8*a4+14224*b6*b10*a4+8869*a6*a8*a7-20608*a9*b5*a10+6664*b6*b10*a8-25669*a6*a9*a7+47355*a6*a4*a7+5040*b10*a10*a7-6048*b6*a10*a7+14672*b6*a10*a3+44415*a6*b9*a4-14280*a6*a8*a3+2744*b10*a9*a8+6888*b10*a9*a4+6832*b6*a6*a7-2856*b10*a9^2-504*a6*b7*a10-30527*a6*b9*a9+9233*a6*b9*a8-22960*b6*a9*a8-9072*b10*a8*a4-14798*b10*b9*a6+12544*b5*b6*a10+58905*a6*b4*a9-1400*a8^3-9436*a6*b8*a8+1400*b6*a9^2+1512*b9*a10*a8+10808*b6*a8^2+56*a8*a9^2-1848*a8*a4^2-5264*a8*b4*a10-5768*a8*b5*a10+18592*a9^2*a5+2016*a8*a10*a7-22232*a8^2*a4-5320*a8*a10*a3+23541*a6*a9*a3+3472*a9^3+13944*a9*a4^2+7952*a9*a8^2+10360*a9*a8*a4-4088*a9*a10*a3-39760*a9*b4*a10+6048*a10*a9*a7-13104*a10*a4*a7+44520*a10*b5*a4-19152*a10*b9*a4+8400*a10*a4*a3+4536*a10*b8*a8+17136*a10*b9*a9+30352*b6*b10*a5+26250*b8*a6*a4+896*a9*b6*a4-23100*a8*b4*a6-5040*a5*b10*a8-3864*b10*a4^2+3136*b10*a8^2-6832*b10^2*a5+30128*a5*b6*a9-15624*a4*b8*a10+36960*a4*b4*a10-10584*b9*a10*a5-18368*a9^2*a4-11592*a5*b10*a4-3304*b10^2*a4-15008*b6^2*a8+21840*b6^2*a9+5040*b8*a10*a9+2338*b8*a6*a9+2688*a5*b10*a9+6776*b10^2*a9-7056*a5*a10*a7-6650*b8*b10*a6-14112*a5^2*a8+17913*a5*b9*a6+21644*b9*b6*a6+64008*b4*a10*a5+10528*a5*a8^2+12306*a6*b5*a9+3423*a6*b5*a8+4599*a6^2*b7+57120*b4*b10*a6-1736*b10^2*a8-10136*b10*a10*a3+12096*a5*b5*a10+16128*a5*a10*a3-3864*a5*a9*a4-3850*b8*b6*a6-52304*a5*a9*a8+9576*b8*b10*a10-18207*b8*a6*a5+13552*b4*b6*a10-13608*b8*b6*a10-16072*b4*b10*a10+85512*a5*a8*a4-35504*a5*b6*a8+14112*a5^2*a9-7056*b9*b6*a10-2520*a5*a6*a7"),
		SMQP("-1344*a5*b5*a8+483*b4*a5*a4-2450*b8*b4*a10+924*b9*b6*a9-896*b10*b4*a5-168*b7*a10*a5-637*b10*b9*a9-4984*b8*a6*a3+1211*b8*a9*a4-3668*a6*a3*a7+3696*b7*a10*a4-5481*b7*a6*a4-1008*b7*a10*a9-168*b6*a4*a7+1820*b6*a8*a7+420*b6*a5*a3+6692*a5*a8*a7-8288*a4*a3*a8-336*b6*a4*a3-4578*b8*b6*a5-2366*b6*a8*a3+3934*a8^2*a3-6762*b5*a4*a9+3948*b6*b5*a5-168*b6*a5*a7-630*b7*a6*a5+1281*b8*a8*a9-1757*b8*a8^2-1029*a6*b7*a8+1736*a4*a7*a8-1904*b8*a8*a4-756*b5*a8^2-910*a8^2*a7+98*b8*a10*a3-140*a7*a9^2-2646*a7*a4^2+273*b5*a4^2+3318*b4*a4^2+1547*b4*a9^2-567*a8*b9*a9+777*a8*b4*a9+14*a10*a3^2+980*b4*a10*a7-1638*a8*a9*a7+6447*a8*b5*a9-2688*a8*a9*a3-6720*a7*a5*a4+1904*a7^2*a6+336*b9*a10*a7+532*a10*a3*a7-336*a8*b7*a10-1246*a9*a4*a3-1995*a9^2*b5+154*a9*b9*a4+4522*a9*a4*a7-1400*a9*a5*a7-1022*a9^2*b9+1148*a9*a5*a3-336*a7^2*a10-280*b6*b9*a8+6048*a4^2*a3+84*a3^2*a6-1764*b4*a10*a3+623*b10*b9*a4-3458*b4^2*a10+7476*a3*a5*a4-1890*b4*b6*a4-168*b7*b6*a10-756*b7*b6*a6-840*b7*b10*a10+1512*a5^2*a3+3045*b8*b6*a9-1379*b8*b10*a9-1260*b4*b6*a9-5719*a4*b4*a9+3430*b8*b5*a10+2688*b8*b9*a10+1848*b7*b10*a6+42*b9*a4^2-1470*b10*b5*a5-756*b5^2*a6-2142*b4*a6*a3+6944*b8^2*a6+6524*b8*a6*a7+392*b10*a9*a7+434*b4*b10*a8-112*b10*a4*a7-728*b10*a8*a7+1995*b10*b5*a9+1176*b9*a6*a7+938*b4*b9*a10-4270*b9*a6*a3-4662*b9*b8*a6+168*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-4298*b4*b5*a10-9835*b8*a5*a9-693*b10*b5*a8-9331*b4*a5*a9-1281*b10*b5*a4+574*b9*a10*a3-8596*b4*b9*a6-3129*b8*b6*a4-434*b10*a9*a3-56*b10*a5*a3+3745*b8*b6*a8-1330*b9*b10*a5+2716*b4*b6*a8-3500*b4*a8*a4-966*b6*a9*a3+84*b6*a9*a7+1442*b10*a8*a3+826*b8*b10*a5+12222*b4*b5*a6+1400*a9^2*a3+1050*b4^2*a6+70*b10*a4*a3-245*a5*b9*a9+868*b4*a6*a7+8778*a5^2*b4+2625*b7*a6*a9+1351*b9*b10*a8-336*b9*b6*a5-1722*b9*b6*a4+3073*b4*a5*a8-6237*b8*a5*a4-2352*b8^2*a10+875*b9*a8*a4-357*b9*a5*a4+434*b8*a9^2-266*b4*b10*a4+630*b8*a5^2+3171*b8*a4^2-364*b10*a5*a7+1036*b4*a8^2+4662*b5*a6*a3-2268*b5*a6*a7+3402*b5*b9*a6+1197*b5*a8*a4-2198*b5*a10*a3+1484*b5*a10*a7-574*b5*b9*a10-2058*b5*a5*a4+882*b5*b8*a6-168*a5*b5*a9-721*b9*a8^2+3479*a5*b9*a8+11221*a5*b8*a8-2352*b8*a10*a7+742*b10*b4*a9+3738*b9*a5^2+1050*b5*b6*a4-266*b9^2*a6-1008*b9^2*a10+546*b5*b6*a8-847*b8*b10*a8+756*a5^2*b5-2982*b5*b6*a9-7322*a5*a3*a8"),
		SMQP("1701*b9^2*a9+3339*b9^2*a4+1995*b9^2*a8+1624*b7*b4*a10-3752*b7*b10*a9+16520*b7*a10*a3+2856*b7*b4*a6+952*b7*b10*a4+32004*a9*a3*a7+11340*b4*a4*a3+16296*b8*a8*a3-1029*b8*b9*a4+735*b8*a4*a3+3171*b8*a9*a7-23940*b7*a6*a3-12516*b4*a8*a3-3024*b4*a8*a7+7728*b5^2*a5+12600*b4*b5*a4+7504*a7*b7*a6+40572*a7*a8*a3+10836*a7*b5*a9-5880*a7*b5*a5-4032*b8*b7*a10-15456*a7*b4*a5-8736*a7*b7*a10+27216*a7*b4*a9+16464*a7*b8*a5-8883*a7*b8*a8+1617*b8*a9*a3-40803*b8*a5*a3-945*b8*b4*a9-1617*b8*b4*a5-4200*b8*b5*a4-4032*a7^2*a5-78876*a7*a4*a3+20160*a7^2*a4+672*b4*b8*a8-168*b4*b9*a8-6384*b4*a4*a7+12264*b4*b5*a8-8064*a7^2*a9-1365*a7*b8*a4-21924*a9*a3^2-8484*b7*a4^2+6748*b7*a4*a8-6069*b8^2*a5+37548*a3*b4*a5+76356*a3^2*a4+4284*a3^2*a5-38136*a3*b4*a9+644*b7*a5*a9-11361*b9^2*a5+12460*b7*a5*a8-17724*b5*a4*a3+1624*b7*b5*a10+13356*b5*a4*a7+840*b5^2*a4-4200*b5^2*a8+2016*a5*a3*a7-2184*b5^2*a9-8064*a7^2*a8-13608*b4*b9*a4-54432*a3^2*a8-9345*b5*b9*a5+3087*b5*b8*a5-22848*b5*b4*a5-4368*b4^2*a9-3549*b5*b8*a8+1841*b7*b8*a6-336*b4^2*a5-6888*b7*a9*a8-14952*b7*b6*a5+3192*b7*b5*a6+13244*b7*a9*a4-420*b7*b6*a4+8316*b7*b6*a9-2296*b7*b10*a8-2541*b5*b9*a8+4900*b7*b6*a8-4144*b7*a9^2-13944*b7*a5^2+5488*b7*b10*a5+9828*b5*b9*a9+952*b7*a8^2-8484*a5*b7*a4-1722*b8^2*a9+2394*b8^2*a8-4536*b9*b5*a4-30324*b5*a9*a3+4116*b5*b8*a9-6888*b5*b4*a9+3192*b8^2*a4+1953*b9*a9*a3+9429*b9*a5*a3-6573*b9*b8*a9-3171*b9*b8*a8+4767*b9*b4*a5+19866*b9*b8*a5+11319*b9*b4*a9+27384*b5*a8*a3-13188*b5*a8*a7+11676*b5*a5*a3-2688*b9*a5*a7+12264*b9*a8*a3-15561*b9*a4*a3-3717*b9*a9*a7+6867*b9*a4*a7-4011*b9*a8*a7-672*b9*b7*a10+3017*b9*b7*a6"),
		SMQP("35742*a5*b5*a8+8043*b4*a5*a4-994*b8*b4*a10+5796*b9*b6*a9-8288*b10*b4*a5+5040*b7*a10*a5-2989*b10*b9*a9+168*b8*a6*a3+259*b8*a9*a4+5628*a6*a3*a7+504*b7*a10*a4-945*b7*a6*a4-1008*b7*a10*a9+4200*b6*a4*a7-11228*b6*a8*a7+4788*b6*a5*a3-8204*a5*a8*a7-5922*a4*a3*a8-2016*b6*a4*a3+16170*b8*b6*a5+6930*b6*a8*a3-882*a8^2*a3+6594*b5*a4*a9-6972*b6*b5*a5-3864*b6*a5*a7-9114*b7*a6*a5-483*b8*a8*a9+1715*b8*a8^2+9513*a6*b7*a8+9898*a4*a7*a8+1400*b8*a8*a4-5166*b5*a8^2+3052*a8^2*a7-3038*b8*a10*a3-700*a7*a9^2-9450*a7*a4^2+10437*b5*a4^2+462*b4*a4^2+2387*b4*a9^2-861*a8*b9*a9-525*a8*b4*a9-1890*a10*a3^2+5740*b4*a10*a7-4368*a8*a9*a7-4725*a8*b5*a9+9954*a8*a9*a3+1176*a7*a5*a4-4256*a7^2*a6-4536*b9*a10*a7+3836*a10*a3*a7-3024*a8*b7*a10-7938*a9*a4*a3+6237*a9^2*b5-1106*a9*b9*a4-1330*a9*a4*a7+9128*a9*a5*a7+322*a9^2*b9-9828*a9*a5*a3-1008*a7^2*a10+3080*b6*b9*a8+6804*a4^2*a3+2772*a3^2*a6-4564*b4*a10*a3+1547*b10*b9*a4-602*b4^2*a10+2016*a3*a5*a4-2058*b4*b6*a4+1008*b7*b6*a10-3276*b7*b6*a6+504*b7*b10*a10+3528*a5^2*a3-4179*b8*b6*a9+2597*b8*b10*a9+1596*b4*b6*a9-4711*a4*b4*a9-994*b8*b5*a10-4536*b8*b9*a10-1512*b7*b10*a6+462*b9*a4^2-3906*b10*b5*a5+9324*b5^2*a6+5754*b4*a6*a3-1232*b8^2*a6-4844*b8*a6*a7-2240*b10*a9*a7+602*b4*b10*a8+1288*b10*a4*a7+1232*b10*a8*a7+2835*b10*b5*a9+5992*b9*a6*a7+770*b4*b9*a10-12054*b9*a6*a3+8498*b9*b8*a6-840*b5^2*a10+581*b8*b10*a4-1246*b4*b8*a6-6482*b4*b5*a10-1547*b8*a5*a9+2079*b10*b5*a8-847*b4*a5*a9-3717*b10*b5*a4+3262*b9*a10*a3-2884*b4*b9*a6-4557*b8*b6*a4+1134*b10*a9*a3-1512*b10*a5*a3-1435*b8*b6*a8+3122*b9*b10*a5-980*b4*b6*a8-1484*b4*a8*a4-1638*b6*a9*a3+6972*b6*a9*a7-2142*b10*a8*a3-7546*b8*b10*a5-4410*b4*b5*a6-1512*a9^2*a3-1470*b4^2*a6-378*b10*a4*a3-7637*a5*b9*a9-5012*b4*a6*a7-5754*a5^2*b4-1575*b7*a6*a9-1925*b9*b10*a8-3528*b9*b6*a5-2898*b9*b6*a4+8491*b4*a5*a8-15057*b8*a5*a4+2016*b8^2*a10-1687*b9*a8*a4-4389*b9*a5*a4-2366*b8*a9^2+910*b4*b10*a4+7434*b8*a5^2+2247*b8*a4^2+3556*b10*a5*a7+3864*b4*b6*a5-980*b4*a8^2+3150*b5*a6*a3+6468*b5*a6*a7-9366*b5*b9*a6-9219*b5*a8*a4-7854*b5*a10*a3+1708*b5*a10*a7+7826*b5*b9*a10-7770*b5*a5*a4+11802*b5*b8*a6-23940*a5*b5*a9+665*b9*a8^2+3605*a5*b9*a8+7259*a5*b8*a8+4536*b8*a10*a7+910*b10*b4*a9+2310*b9*a5^2-4830*b5*b6*a4-5642*b9^2*a6+504*b9^2*a10+9786*b5*b6*a8+805*b8*b10*a8+7308*a5^2*b5-5334*b5*b6*a9+1512*a5*a3*a8-2016*a5^2*a7"),
		SMQP("-5257*a5*b5*a8-231*b4*a5*a4-63*b9*b6*a9+168*b10*b4*a5-252*b7*a10*a5-14*b10*b9*a9-3640*b8*a6*a3-168*b8*a9*a4-3584*a6*a3*a7+252*b7*a10*a4+1323*b7*a6*a4-609*b6*a4*a7+581*b6*a8*a7+2625*b6*a5*a3+4109*a5*a8*a7+2212*a4*a3*a8+2247*b6*a4*a3+273*b8*b6*a5-2156*b6*a8*a3+28*a8^2*a3+238*b5*a4*a9-525*b6*b5*a5-672*b6*a5*a7-2541*b7*a6*a5-42*b8*a8*a9+1008*a6*b7*a8-637*a4*a7*a8+126*b8*a8*a4+140*b5*a8^2-28*a8^2*a7-504*b8*a10*a3-56*a7*a9^2+693*a7*a4^2-2016*b5*a4^2+42*b4*a9^2-42*a8*b9*a9-84*a8*b4*a9+308*a10*a3^2-700*b4*a10*a7+84*a8*a9*a7-336*a8*b5*a9-84*a8*a9*a3-1617*a7*a5*a4+728*a7^2*a6+252*b9*a10*a7-308*a10*a3*a7+119*a9*a4*a3+238*a9^2*b5+35*a9*b9*a4-119*a9*a4*a7-3773*a9*a5*a7-28*a9^2*b9+8813*a9*a5*a3-581*b6*b9*a8-2331*a4^2*a3+4872*a3^2*a6+2464*b4*a10*a3+28*b10*b9*a4+798*a3*a5*a4-252*b7*b6*a10-3339*b7*b6*a6-3339*a5^2*a3+168*b8*b6*a9+42*b8*b10*a9+147*b4*b6*a9+105*a4*b4*a9+756*b8*b5*a10+756*b8*b9*a10+1008*b7*b10*a6-819*b9*a4^2+1358*b10*b5*a5+252*b5^2*a6+2184*b4*a6*a3-2520*b8^2*a6+364*b8*a6*a7-28*b10*a9*a7-28*b10*a4*a7+28*b10*a8*a7+182*b10*b5*a9-112*b9*a6*a7+1204*b4*b9*a10+1232*b9*a6*a3+4172*b9*b8*a6-1036*b5^2*a10+980*b4*b5*a10+1176*b8*a5*a9-140*b10*b5*a8-525*b4*a5*a9+140*b10*b5*a4+1064*b9*a10*a3+3444*b4*b9*a6+28*b10*a9*a3+280*b10*a5*a3-126*b8*b6*a8+322*b9*b10*a5-147*b6*a9*a3+147*b6*a9*a7-28*b10*a8*a3-42*b8*b10*a5-9492*b4*b5*a6+56*a9^2*a3+28*b10*a4*a3-469*a5*b9*a9-924*b4*a6*a7-903*a5^2*b4-1890*b7*a6*a9-28*b9*b10*a8-357*b9*b6*a5+735*b9*b6*a4+84*b4*a5*a8-273*b8*a5*a4-504*b8^2*a10+637*b9*a8*a4+2268*b9*a5*a4+84*b8*a9^2-1575*b8*a5^2-280*b10*a5*a7-21*b4*b6*a5+2184*b5*a6*a3-3136*b5*a6*a7+2212*b5*b9*a6+1925*b5*a8*a4+1484*b5*a10*a3+56*b5*a10*a7-1568*b5*b9*a10+1323*b5*a5*a4+196*b5*b8*a6+3598*a5*b5*a9+28*b9*a8^2-2135*a5*b9*a8+420*a5*b8*a8-252*b8*a10*a7-168*b10*b4*a9+1617*b9*a5^2+1596*b5*b6*a4-2632*b9^2*a6-252*b9^2*a10-1645*b5*b6*a8+567*a5^2*b5-504*b5*b6*a9-7448*a5*a3*a8+504*a5^2*a7"),
		SMQP("-5943*a6*a3*a4+2247*a6*a3*a5-9492*a6*b5*b6+2961*a6*b5*a5-1596*a6*a3*b6-4032*a6*b4*b6-5796*a6*b4*a5-6195*a6*b5*a4-1638*a4*b6*a5-1974*a6*b4*a4+3108*a4^2*a5+336*a4*b6^2-14154*a4*a5^2-1932*b6^2*a5-924*b6*a5^2-3906*a4^3+462*a4^2*b6+3528*a5^3-448*b10*a5^2+10101*b10*b5*a6-3780*b8*a10*a5-1666*b10*b6*a9-4088*b10*b5*a10+2016*b10*b9*a10-224*b10*a6*a7+336*b10*a6*a3-140*b6*a8*a4+350*b6*b10*a4-4669*a6*a8*a7-1260*a9*b5*a10-70*b6*b10*a8+1197*a6*a9*a7+2177*a6*a4*a7+504*b10*a10*a7-1512*b6*a10*a7+490*b6*a10*a3+3157*a6*b9*a4+3696*a6*a8*a3+728*b10*a9*a8-854*b10*a9*a4+1316*b6*a6*a7+812*b10*a9^2-6804*a6*b7*a10-1491*a6*b9*a9-1778*a6*b9*a8-210*b6*a9*a8-1610*b10*a8*a4-3157*b10*b9*a6+3038*b5*b6*a10+7014*a6*b4*a9-196*a8^3-1631*a6*b8*a8-2408*b6*a9^2+1260*b9*a10*a8+658*b6*a8^2+868*a8*a9^2+3962*a8*a4^2+644*a8*b4*a10-1120*a8*b5*a10-4676*a9^2*a5+1260*a8*a10*a7-14*a8^2*a4-2660*a8*a10*a3+2331*a6*a9*a3+504*a9^3+2464*a9*a4^2-168*a9*a8^2-2702*a9*a8*a4-504*a9*a10*a3-5040*a9*b4*a10+756*a10*a9*a7-1764*a10*a4*a7+3206*a10*b5*a4+252*a10*b9*a4+4606*a10*a4*a3+1512*a10^2*b7+1260*a10*b9*a9+476*b6*b10*a5-1673*b8*a6*a4+826*a9*b6*a4-462*a8*b4*a6-518*a5*b10*a8+1358*b10*a4^2-28*b10*a8^2+280*b10^2*a5+1834*a5*b6*a9-504*a4*b8*a10+3710*a4*b4*a10+756*b9*a10*a5-812*a9^2*a4+1330*a5*b10*a4-224*b10^2*a4-1162*b6^2*a8+2478*b6^2*a9+924*b8*a6*a9-2282*a5*b10*a9+280*b10^2*a9-504*a5*a10*a7-595*b8*b10*a6+13580*a5^2*a8-56*a5*b9*a6+6748*b9*b6*a6-2758*b4*a10*a5-2618*a5*a8^2+2751*a6*b5*a9-4557*a6*b5*a8+5922*a6^2*b7+1596*b4*b10*a6+224*b10^2*a8+56*b10*a10*a3+3038*a5*b5*a10-2030*a5*a10*a3+10388*a5*a9*a4+1561*b8*b6*a6+3094*a5*a9*a8+504*b8*b10*a10+10108*b8*a6*a5+3542*b4*b6*a10-1512*b8*b6*a10-1064*b4*b10*a10-7826*a5*a8*a4+4690*a5*b6*a8-5432*a5^2*a9-4536*b9*b6*a10+812*a5*a6*a7"),
		SMQP("163800*a5*b5*a8+3171*b4*a5*a4-2450*b8*b4*a10+37212*b9*b6*a9-3584*b10*b4*a5+25536*b7*a10*a5-16457*b10*b9*a9+3192*b8*a6*a3-3577*b8*a9*a4-108948*a6*a3*a7+49896*b7*a10*a4-72513*b7*a6*a4-37296*b7*a10*a9-23688*b6*a4*a7-15372*b6*a8*a7+43092*b6*a5*a3+107352*a5*a8*a7-80724*a4*a3*a8+31248*b6*a4*a3+33642*b8*b6*a5-31206*b6*a8*a3+28518*a8^2*a3-78246*b5*a4*a9+27300*b6*b5*a5-36792*b6*a5*a7-35910*b7*a6*a5+1113*b8*a8*a9-1421*b8*a8^2+56511*a6*b7*a8+3780*a4*a7*a8+6328*b8*a8*a4-33432*b5*a8^2+12726*a8^2*a7-10654*b8*a10*a3-33894*a7*a9^2-31374*a7*a4^2+69321*b5*a4^2+3318*b4*a4^2+26453*b4*a9^2-23961*a8*b9*a9-147*a8*b4*a9-24234*a10*a3^2+3276*b4*a10*a7-23184*a8*a9*a7+27405*a8*b5*a9+48132*a8*a9*a3-135324*a7*a5*a4+40320*a7^2*a6-9576*b9*a10*a7+44604*a10*a3*a7-19152*a8*b7*a10-46536*a9*a4*a3+62097*a9^2*b5+26600*a9*b9*a4+64512*a9*a4*a7+9828*a9*a5*a7-19180*a9^2*b9-19446*a9*a5*a3-19152*a7^2*a10-8624*b6*b9*a8+52920*a4^2*a3+74340*a3^2*a6+476*b4*a10*a3+11095*b10*b9*a4-3458*b4^2*a10+177912*a3*a5*a4-1890*b4*b6*a4-12096*b7*b6*a10-29484*b7*b6*a6-3528*b7*b10*a10+12348*a5^2*a3-903*b8*b6*a9+2737*b8*b10*a9-21252*b4*b6*a9-28567*a4*b4*a9+39886*b8*b5*a10+34776*b8*b9*a10+21672*b7*b10*a6-48846*b9*a4^2-27426*b10*b5*a5+87444*b5^2*a6+41538*b4*a6*a3+32144*b8^2*a6+1540*b8*a6*a7-23184*b10*a9*a7+1106*b4*b10*a8+12600*b10*a4*a7+4032*b10*a8*a7+42231*b10*b5*a9+1624*b9*a6*a7+101458*b4*b9*a10-32718*b9*a6*a3-8750*b9*b8*a6-23856*b5^2*a10-119*b8*b10*a4+1442*b4*b8*a6-167426*b4*b5*a10+8351*b8*a5*a9-21*b10*b5*a8-42385*b4*a5*a9-19593*b10*b5*a4+18494*b9*a10*a3-96460*b4*b9*a6-3129*b8*b6*a4+6342*b10*a9*a3-10920*b10*a5*a3-2471*b8*b6*a8-6398*b9*b10*a5+2380*b4*b6*a8-4172*b4*a8*a4-17262*b6*a9*a3+55692*b6*a9*a7+13818*b10*a8*a3-10346*b8*b10*a5+80262*b4*b5*a6+2226*a9^2*a3+1050*b4^2*a6-13314*b10*a4*a3-511*a5*b9*a9-28196*b4*a6*a7+48930*a5^2*b4+38871*b7*a6*a9-8617*b9*b10*a8+22008*b9*b6*a5+1638*b9*b6*a4+3157*b4*a5*a8-23289*b8*a5*a4-36288*b8^2*a10+55447*b9*a8*a4+5943*b9*a5*a4+3206*b8*a9^2-266*b4*b10*a4-103194*b8*a5^2+3171*b8*a4^2+11844*b10*a5*a7-8232*b4*b6*a5+28*b4*a8^2+32886*b5*a6*a3+3108*b5*a6*a7-92190*b5*b9*a6-35427*b5*a8*a4-7854*b5*a10*a3+14364*b5*a10*a7+42994*b5*b9*a10-88158*b5*a5*a4-26166*b5*b8*a6-127302*a5*b5*a9-5117*b9*a8^2-33173*a5*b9*a8+57421*a5*b8*a8+3528*b8*a10*a7+7462*b10*b4*a9+45906*b9*a5^2-9030*b5*b6*a4+30422*b9^2*a6-36792*b9^2*a10+67242*b5*b6*a8+665*b8*b10*a8+58968*a5^2*b5-104454*b5*b6*a9-174846*a5*a3*a8+18144*a5^2*a7"),
		SMQP("-2772*b10^2*a6+2541*b6*a10*a5-371*a6^2*b8+959*a6*b10*a9+623*a6*b6*a8-6566*a6*b10*a8+18900*b6*b10*a6-9702*b6*b10*a10+3108*a6*b6*a5+5376*a6*b6*a9-14973*a6*b6*a4-1694*a6*a10*a3+455*a6*a10*a7-2429*a6*a5*a9+4298*a6*b10*a4+7616*a6*a5*a8-9492*a6^2*b4+3374*a6*b5*a10+9695*a6*b4*a10-7427*a6*b9*a10+5635*a6*a9*a4-1323*a6^2*b5-2373*a6*a9*a8+1295*a6*a8*a4+518*b10*a10*a9+1778*a6*b10*a5+2646*b6^2*a10+2226*a10*a9*a8-5145*a10*b6*a9-7525*a10^2*b4-252*a10^2*b8-217*a10^2*b5+2128*b10*a10*a8+1505*a10*b6*a8+8589*a10*b6*a4+91*a10*a9*a4+4529*a10*a5*a8-5390*a10*a5*a9-6391*a10*a8*a4-1050*a10*a5*a4-2128*a10*b10*a4-1673*a10^2*a3+3913*a10*b8*a6+2142*a10^2*b9-2317*a10*b10*a5+378*a10^2*a7-1358*a6*a9^2-5481*a6*a4^2+7189*a6^2*b9-525*a6^2*a3-364*a6*a8^2-6972*a6*a5*a4+2268*b10^2*a10-11340*b6^2*a6+3353*a6^2*a7+28*a10*a9^2+1638*a6*a5^2+3465*a10*a5^2-364*a10*a8^2+5355*a10*a4^2"),
		SMQP("-1008*b10^2*a6-189*b6*a10*a5-1771*a6^2*b8-182*a6*b10*a9+910*a6*b6*a8-7*a6*b10*a8+1008*b6*b10*a6-378*b6*b10*a10-210*a6*b6*a5+357*a6*b6*a9-399*a6*b6*a4+749*a6*a10*a3-427*a6*a10*a7+1757*a6*a5*a9+1519*a6*b10*a4-2240*a6*a5*a8-2247*a6^2*b4+91*a6*b5*a10+2590*a6*b4*a10-1421*a6*b9*a10-805*a6*a9*a4+336*a6^2*b5-1218*a6*a9*a8+2254*a6*a8*a4+14*b10*a10*a9+28*a6*b10*a5+1386*b6^2*a10+210*a10*a9*a8-357*a10*b6*a9-217*a10^2*b4-252*a10^2*b8-217*a10^2*b5+364*b10*a10*a8-259*a10*b6*a8+399*a10*b6*a4-161*a10*a9*a4-343*a10*a5*a8+238*a10*a5*a9-343*a10*a8*a4+378*a10*a5*a4-1120*a10*b10*a4-161*a10^2*a3+1057*a10*b8*a6+378*a10^2*b9-49*a10*b10*a5+126*a10^2*a7+140*a6*a9^2-1638*a6*a4^2+854*a6^2*b9+651*a6^2*a3+889*a6*a8^2+357*a6*a5*a4+252*b10^2*a10-1260*b6^2*a6-119*a6^2*a7+28*a10*a9^2-630*a6*a5^2-21*a10*a5^2-112*a10*a8^2+441*a10*a4^2"),
		SMQP("3213*b9^2*a9+1323*b9^2*a4+3003*b9^2*a8+560*b7*b4*a10-2716*b7*b10*a9+11284*b7*a10*a3+4032*b7*b4*a6+476*b7*b10*a4+11907*a9*a3*a7+5229*b4*a4*a3+8484*b8*a8*a3-1533*b8*b9*a4+1932*b8*a4*a3+2919*b8*a9*a7-16422*b7*a6*a3-6027*b4*a8*a3-1512*b4*a8*a7+6279*b5^2*a5+6804*b4*b5*a4+1169*a7*b7*a6+18774*a7*a8*a3-3024*a7*b5*a9-903*a7*b5*a5-504*b8*b7*a10-10605*a7*b4*a5-4452*a7*b7*a10+15687*a7*b4*a9+10605*a7*b8*a5-6363*a7*b8*a8+924*b8*a9*a3-24990*b8*a5*a3-1008*b8*b4*a9-420*b8*b4*a5-1995*b8*b5*a4-9576*a7^2*a5-35343*a7*a4*a3+9387*a7^2*a4+168*b4*b8*a8-2310*b4*b9*a8-2604*b4*a4*a7+9744*b4*b5*a8+63*a7^2*a9-1113*a7*b8*a4-7938*a9*a3^2-3507*b7*a4^2-679*b7*a4*a8-5124*b8^2*a5+17892*a3*b4*a5+38178*a3^2*a4-10458*a3^2*a5-23331*a3*b4*a9+1309*b7*a5*a9-8778*b9^2*a5+2303*b7*a5*a8-14574*b5*a4*a3-4480*b7*b5*a10+9954*b5*a4*a7-2688*b5^2*a4-4641*b5^2*a8+20853*a5*a3*a7+6384*b5^2*a9-5607*a7^2*a8-6930*b4*b9*a4-30240*a3^2*a8-7203*b5*b9*a5+1197*b5*b8*a5-11949*b5*b4*a5-1848*b4^2*a9-1281*b5*b8*a8-3752*b7*b8*a6-336*b4^2*a5-1428*b7*a9*a8-12390*b7*b6*a5+14469*b7*b5*a6+5803*b7*a9*a4+567*b7*b6*a4+6657*b7*b6*a9-1820*b7*b10*a8-2982*b5*b9*a8+3563*b7*b6*a8-1568*b7*a9^2-11298*b7*a5^2+5012*b7*b10*a5-882*b5*b9*a9+1988*b7*a8^2-819*a5*b7*a4+42*b8^2*a9+2646*b8^2*a8+1386*b9*b5*a4-17598*b5*a9*a3+777*b5*b8*a9+2499*b5*b4*a9+1806*b8^2*a4+1323*b9*a9*a3+10752*b9*a5*a3-5439*b9*b8*a9-3297*b9*b8*a8+2184*b9*b4*a5+16086*b9*b8*a5+3822*b9*b4*a9+18249*b5*a8*a3-2289*b5*a8*a7-735*b5*a5*a3-2499*b9*a5*a7+6909*b9*a8*a3-8379*b9*a4*a3+756*b9*a9*a7+2520*b9*a4*a7-2940*b9*a8*a7-672*b9*b7*a10+2716*b9*b7*a6"),
		SMQP("-21*b9*b8*a5+7*a7*b8*a5-14*b7*a6*a3+7*b7*b6*a5+28*b5*b4*a5-7*b5*b8*a5+7*b5*b9*a5-14*b9*b4*a5-7*a3*b4*a5-7*a5*b7*a4+14*b8^2*a5+14*b7*b5*a6+7*b7*a5^2+7*b9*a5*a3-7*b9*a5*a7+7*b9^2*a5+7*a7*b5*a5-21*b5*a5*a3"),
		SMQP("1575*b9^2*a9+4221*b9^2*a4+3885*b9^2*a8+1624*b7*b4*a10-4760*b7*b10*a9+18536*b7*a10*a3+3024*b7*b4*a6+952*b7*b10*a4+45360*a9*a3*a7+11592*b4*a4*a3+15792*b8*a8*a3-2856*b8*b9*a4+3696*b8*a4*a3+6636*b8*a9*a7-26544*b7*a6*a3-12264*b4*a8*a3-3024*b4*a8*a7+11949*b5^2*a5+13104*b4*b5*a4+7504*a7*b7*a6+38304*a7*a8*a3+8505*a7*b5*a9-7896*a7*b5*a5-4032*b8*b7*a10-17976*a7*b4*a5-9744*a7*b7*a10+29736*a7*b4*a9+16968*a7*b8*a5-10584*a7*b8*a8+1428*b8*a9*a3-45024*b8*a5*a3-2016*b8*b4*a9-1680*b8*b4*a5-4200*b8*b5*a4-6552*a7^2*a5-88200*a7*a4*a3+22680*a7^2*a4+672*b4*b8*a8-4704*b4*b9*a8-6384*b4*a4*a7+18816*b4*b5*a8-10584*a7^2*a9-2184*a7*b8*a4-25704*a9*a3^2-8736*b7*a4^2+3472*b7*a4*a8-7392*b8^2*a5+38304*a3*b4*a5+85680*a3^2*a4+5040*a3^2*a5-45192*a3*b4*a9-700*b7*a5*a9-17283*b9^2*a5+14392*b7*a5*a8-21819*b5*a4*a3+2632*b7*b5*a10+13293*b5*a4*a7+1092*b5^2*a4-5523*b5^2*a8+3024*a5*a3*a7-1428*b5^2*a9-7560*a7^2*a8-14112*b4*b9*a4-61488*a3^2*a8-9030*b5*b9*a5+3465*b5*b8*a5-23163*b5*b4*a5-4368*b4^2*a9-4242*b5*b8*a8+1232*b7*b8*a6-336*b4^2*a5-7896*b7*a9*a8-19488*b7*b6*a5+3675*b7*b5*a6+17276*b7*a9*a4-168*b7*b6*a4+9324*b7*b6*a9-2296*b7*b10*a8-6006*b5*b9*a8+6160*b7*b6*a8-6160*b7*a9^2-15344*b7*a5^2+7504*b7*b10*a5+10143*b5*b9*a9+2968*b7*a8^2-8064*a5*b7*a4-2856*b8^2*a9+5040*b8^2*a8-1071*b9*b5*a4-38829*b5*a9*a3+3738*b5*b8*a9-9471*b5*b4*a9+3696*b8^2*a4+4095*b9*a9*a3+15099*b9*a5*a3-6762*b9*b8*a9-6258*b9*b8*a8+3381*b9*b4*a5+27993*b9*b8*a5+16989*b9*b4*a9+32928*b5*a8*a3-11109*b5*a8*a7+10731*b5*a5*a3-3192*b9*a5*a7+19824*b9*a8*a3-22491*b9*a4*a3-5103*b9*a9*a7+9765*b9*a4*a7-6405*b9*a8*a7-1680*b9*b7*a10+4571*b9*b7*a6"),
		SMQP("777*b7*b4*a5+945*b7*b9*a5-1407*b7*b4*a9-1155*b7*b8*a5+63*b7^2*a6-126*b7*b8*a9-2016*b7*b5*a9+336*b7^2*a10+2919*a7*b7*a8+1008*a7*b7*a5-3255*a7*b7*a4-336*b7*b4*a8-11088*a7^2*a3-10080*a3^3+273*b7*b5*a8+1113*b7*a9*a7+2436*b7*b5*a4-6888*b7*a8*a3+4977*b7*a4*a3+119*b7*a5*a3-777*b7*a9*a3-903*b7*b9*a4+1134*b7*b8*a8+441*b7*b9*a8+231*b7*b9*a9+2016*a7^3+19152*a7*a3^2-672*b7*b8*a4+168*b7*b4*a4-63*b7*b5*a5"),
		SMQP("987*b7*b4*a5+2499*b7*b9*a5-1365*b7*b4*a9-3465*b7*b8*a5-315*b7^2*a6+966*b7*b8*a9-4704*b7*b5*a9+336*b7^2*a10+2709*a7*b7*a8+1848*a7*b7*a5-2709*a7*b7*a4-11088*a7^2*a3-10080*a3^3+3171*b7*b5*a8+1323*b7*a9*a7+588*b7*b5*a4-8904*b7*a8*a3+5187*b7*a4*a3-1603*b7*a5*a3+1029*b7*a9*a3-357*b7*b9*a4+714*b7*b8*a8-21*b7*b9*a8+21*b7*b9*a9+2016*a7^3+19152*a7*a3^2-672*b7*b8*a4+336*b7*b4*a4+1155*b7*b5*a5"),
		SMQP("-7*b7*a5*a3+7*b7*b5*a5"),
		SMQP("-257271*a6*a3*a4+66423*a6*a3*a5-84252*a6*b5*b6+16317*a6*b5*a5-107520*a6*a3*b6-201348*a6*b4*b6-143577*a6*b4*a5-170898*a6*b5*a4+76608*a4*b6*a5-313152*a6*b4*a4-329952*a4^2*a5-17472*a4*b6^2+129192*a4*a5^2-90048*b6^2*a5-56784*b6*a5^2+60984*a4^3-103152*a4^2*b6-62496*a5^3+38080*b10*a5^2+74382*b10*b5*a6-63504*b8*a10*a5-82208*b10*b6*a9-34048*b10*b5*a10-1008*b10*b9*a10-75712*b10*a6*a7+168000*b10*a6*a3+232232*b6*a8*a4+85120*b6*b10*a4-2177*a6*a8*a7-66360*a9*b5*a10-23744*b6*b10*a8-105231*a6*a9*a7+232897*a6*a4*a7+25200*b10*a10*a7-30240*b6*a10*a7+66752*b6*a10*a3+136493*a6*b9*a4-40488*a6*a8*a3-56*b10*a9*a8+34832*b10*a9*a4+40432*b6*a6*a7-12152*b10*a9^2-6048*a6*b7*a10-84357*a6*b9*a9+2135*a6*b9*a8-60984*b6*a9*a8-17752*b10*a8*a4-49742*b10*b9*a6+60256*b5*b6*a10+175623*a6*b4*a9+14896*a8^3-16912*a6*b8*a8+728*b6*a9^2+15120*b9*a10*a8-16912*b6*a8^2+9128*a8*a9^2-78344*a8*a4^2-24248*a8*b4*a10-68096*a8*b5*a10+62384*a9^2*a5+1008*a8*a10*a7-75544*a8^2*a4+20720*a8*a10*a3+121023*a6*a9*a3+8400*a9^3+63224*a9*a4^2+17976*a9*a8^2+27104*a9*a8*a4-34440*a9*a10*a3-119784*a9*b4*a10+27216*a10*a9*a7-50400*a10*a4*a7+197848*a10*b5*a4-53424*a10*b9*a4+22904*a10*a4*a3+46368*a10*b8*a8-6048*a10^2*b7+51408*a10*b9*a9+149632*b6*b10*a5+205394*b8*a6*a4-31864*a9*b6*a4+588*a8*b4*a6+12152*a5*b10*a8-25592*b10*a4^2+15232*b10*a8^2-35392*b10^2*a5+115304*a5*b6*a9-120960*a4*b8*a10+75880*a4*b4*a10-47376*b9*a10*a5-53872*a9^2*a4-88984*a5*b10*a4-21280*b10^2*a4-34832*b6^2*a8+75264*b6^2*a9+22176*b8*a10*a9-33474*b8*a6*a9+7784*a5*b10*a9+19040*b10^2*a9-31248*a5*a10*a7-54194*b8*b10*a6-78512*a5^2*a8+98063*a5*b9*a6+63644*b9*b6*a6+269080*b4*a10*a5-90328*a5*a8^2+54138*a6*b5*a9+1533*a6*b5*a8+4473*a6^2*b7+187656*b4*b10*a6+6160*b10^2*a8-56672*b10*a10*a3+63448*a5*b5*a10+68600*a5*a10*a3-151928*a5*a9*a4+1358*b8*b6*a6-99400*a5*a9*a8+54432*b8*b10*a10-121513*b8*a6*a5+84448*b4*b6*a10-49392*b8*b6*a10-56224*b4*b10*a10+546392*a5*a8*a4-120232*a5*b6*a8+51968*a5^2*a9-26208*b9*b6*a10-10808*a5*a6*a7"),
		SMQP("2940*a6*a3*a4+924*a6*a3*a5+10920*a6*b5*b6+1134*a6*b5*a5+672*a6*a3*b6+3528*a6*b4*b6+2268*a6*b4*a5+8736*a6*b5*a4-2016*a4*b6*a5+1050*a6*b4*a4+840*a4^2*a5-1176*a4*b6^2+3990*a4*a5^2+1848*b6^2*a5+3612*b6*a5^2+2898*a4^3+1596*a4^2*b6+3528*a5^3-4984*b10*a5^2-12012*b10*b5*a6+315*b8*a10*a5+1358*b10*b6*a9+3535*b10*b5*a10-2037*b10*b9*a10-1820*b10*a6*a7+4116*b10*a6*a3+238*b6*a8*a4+1484*b6*b10*a4+728*a6*a8*a7+3465*a9*b5*a10-1666*b6*b10*a8+2604*a6*a9*a7-2548*a6*a4*a7-504*b10*a10*a7+1512*b6*a10*a7-2156*b6*a10*a3-2702*a6*b9*a4-3360*a6*a8*a3-952*b10*a9*a8+1204*b10*a9*a4-448*b6*a6*a7-406*b10*a9^2+5103*a6*b7*a10+1134*a6*b9*a9+4606*a6*b9*a8+2646*b6*a9*a8+3892*b10*a8*a4+5180*b10*b9*a6-3934*b5*b6*a10-5712*a6*b4*a9+476*a8^3-1106*a6*b8*a8+1750*b6*a9^2-1995*b9*a10*a8-2576*b6*a8^2-938*a8*a9^2-5866*a8*a4^2+1778*a8*b4*a10+1358*a8*b5*a10+28*a9^2*a5-504*a8*a10*a7+2422*a8^2*a4+2212*a8*a10*a3-1260*a6*a9*a3-84*a9^3-182*a9*a4^2-462*a9*a8^2+616*a9*a8*a4-462*a9*a10*a3+2499*a9*b4*a10-1512*a10*a9*a7+1512*a10*a4*a7-5635*a10*b5*a4+462*a10*b9*a4-938*a10*a4*a3+483*a10*b8*a8-1008*a10^2*b7-1554*a10*b9*a9+476*b6*b10*a5+1162*b8*a6*a4-3710*a9*b6*a4+168*a8*b4*a6-3206*a5*b10*a8-3556*b10*a4^2-154*b10*a8^2+196*b10^2*a5-6482*a5*b6*a9+483*a4*b8*a10-3010*a4*b4*a10+1365*b9*a10*a5+1120*a9^2*a4+1204*a5*b10*a4+322*b10^2*a4+1988*b6^2*a8-2436*b6^2*a9+42*b8*a10*a9+336*b8*a6*a9+5194*a5*b10*a9-182*b10^2*a9+2520*a5*a10*a7-1330*b8*b10*a6-1540*a5^2*a8-3080*a5*b9*a6-8120*b9*b6*a6-637*b4*a10*a5-2954*a5*a8^2-4620*a6*b5*a9+168*a6*b5*a8-5418*a6^2*b7-798*b4*b10*a6-322*b10^2*a8+518*b10*a10*a3-1372*a5*b5*a10-2282*a5*a10*a3-6874*a5*a9*a4-1022*b8*b6*a6+6286*a5*a9*a8-483*b8*b10*a10-2114*b8*a6*a5-4354*b4*b6*a10+1533*b8*b6*a10+2926*b4*b10*a10+2254*a5*a8*a4+3178*a5*b6*a8-3920*a5^2*a9+4578*b9*b6*a10-2212*a5*a6*a7"),
		SMQP("910*a5*b5*a8-84*b4*a5*a4+84*b10*b4*a5-434*b8*a6*a3-364*a6*a3*a7+63*b7*a6*a4+84*b6*a5*a3+350*a4*a3*a8+336*b6*a4*a3-21*b8*b6*a5-322*b6*a8*a3+14*a8^2*a3-28*b5*a4*a9-42*b6*b5*a5-378*b7*a6*a5-14*b5*a8^2+378*b5*a4^2+154*a10*a3^2+42*a8*b5*a9-42*a8*a9*a3+28*a9*a4*a3-28*a9^2*b5+784*a9*a5*a3-378*a4^2*a3+420*a3^2*a6+350*b4*a10*a3+462*a3*a5*a4-63*b7*b6*a6-504*a5^2*a3-161*b10*b5*a5+504*b5^2*a6-105*b4*a6*a3-126*b8^2*a6-63*b8*a6*a7-14*b10*b5*a9+63*b9*a6*a7+511*b9*a6*a3+189*b9*b8*a6+154*b5^2*a10-350*b4*b5*a10-42*b8*a5*a9+14*b10*b5*a8-21*b4*a5*a9-14*b10*b5*a4+126*b4*b9*a6+14*b10*a9*a3+140*b10*a5*a3+21*b9*b10*a5-42*b6*a9*a3-14*b10*a8*a3-21*b8*b10*a5-84*b4*b5*a6+28*a9^2*a3+14*b10*a4*a3+42*a5*b9*a9+147*a5^2*b4-42*b9*b6*a5+42*b4*a5*a8+21*b8*a5*a4+42*b9*a5*a4-315*b8*a5^2-42*b4*b6*a5-735*b5*a6*a3+301*b5*a6*a7-637*b5*b9*a6-350*b5*a8*a4-308*b5*a10*a3-441*b5*a5*a4+497*b5*b8*a6-763*a5*b5*a9-21*a5*b9*a8+21*a5*b8*a8+147*b9*a5^2-336*b5*b6*a4-63*b9^2*a6+322*b5*b6*a8+630*a5^2*b5+42*b5*b6*a9-952*a5*a3*a8"),
		SMQP("308*b10^2*a6-301*b6*a10*a5+567*a6^2*b8-133*a6*b10*a9+399*a6*b6*a8+756*a6*b10*a8-2100*b6*b10*a6+1078*b6*b10*a10-252*a6*b6*a5-882*a6*b6*a9+1449*a6*b6*a4-56*a6*a10*a3-35*a6*a10*a7-679*a6*a5*a9-448*a6*b10*a4+168*a6*a5*a8+966*a6^2*b4-98*a6*b5*a10-1547*a6*b4*a10+777*a6*b9*a10-287*a6*a9*a4+357*a6^2*b5+315*a6*a9*a8-609*a6*a8*a4-28*b10*a10*a9-350*a6*b10*a5-294*b6^2*a10-336*a10*a9*a8+455*a10*b6*a9+735*a10^2*b4+28*a10^2*b8-21*a10^2*b5-266*b10*a10*a8-35*a10*b6*a8-1057*a10*b6*a4+77*a10*a9*a4-301*a10*a5*a8+322*a10*a5*a9+637*a10*a8*a4-84*a10*a5*a4+266*a10*b10*a4+231*a10^2*a3-301*a10*b8*a6-238*a10^2*b9+273*a10*b10*a5-42*a10^2*a7+182*a6*a9^2+735*a6*a4^2-1071*a6^2*b9-63*a6^2*a3-42*a6*a8^2+210*a6*a5*a4-252*b10^2*a10+1260*b6^2*a6-189*a6^2*a7+56*a10*a9^2+378*a6*a5^2-77*a10*a5^2+70*a10*a8^2-581*a10*a4^2"),
		SMQP("-7*a7*b8*a5-14*b8*a5*a3+21*b9*a5*a3-35*b5*b9*a5+7*b9*a5*a7-14*b8^2*a5+21*b9*b8*a5+7*b5*a5*a3+14*b9*b4*a5+21*b5*b8*a5-7*b7*b6*a5+7*a5*b7*a4-7*a7*b5*a5+14*b5^2*a5-35*b7*a5^2-14*b5*b4*a5-7*a3*b4*a5-7*b9^2*a5"),
		SMQP("2646*b9^2*a9+3654*b9^2*a4+4410*b9^2*a8+1372*b7*b4*a10-4592*b7*b10*a9+20552*b7*a10*a3+5460*b7*b4*a6+952*b7*b10*a4+35658*a9*a3*a7+11025*b4*a4*a3+15750*b8*a8*a3-2961*b8*b9*a4+3780*b8*a4*a3+6111*b8*a9*a7-29589*b7*a6*a3-12411*b4*a8*a3-3024*b4*a8*a7+14049*b5^2*a5+13356*b4*b5*a4+4921*a7*b7*a6+33327*a7*a8*a3+1890*a7*b5*a9-10395*a7*b5*a5-2520*b8*b7*a10-19341*a7*b4*a5-8820*a7*b7*a10+29295*a7*b4*a9+18585*a7*b8*a5-11403*a7*b8*a8+1386*b8*a9*a3-45675*b8*a5*a3-2016*b8*b4*a9-1260*b8*b4*a5-4095*b8*b5*a4-13608*a7^2*a5-77112*a7*a4*a3+19467*a7^2*a4+504*b4*b8*a8-4158*b4*b9*a8-5796*b4*a4*a7+18648*b4*b5*a8-5985*a7^2*a9-2205*a7*b8*a4-21609*a9*a3^2-7875*b7*a4^2+805*b7*a4*a8-9324*b8^2*a5+37863*a3*b4*a5+78561*a3^2*a4-10521*a3^2*a5-43596*a3*b4*a9+917*b7*a5*a9-16695*b9^2*a5+6223*b7*a5*a8-30996*b5*a4*a3-4676*b7*b5*a10+19656*b5*a4*a7+1008*b5^2*a4-10647*b5^2*a8+26397*a5*a3*a7+4032*b5^2*a9-7623*a7^2*a8-13986*b4*b9*a4-56952*a3^2*a8-13734*b5*b9*a5+7875*b5*b8*a5-28539*b5*b4*a5-4032*b4^2*a9-4725*b5*b8*a8-2800*b7*b8*a6-504*b4^2*a5-5880*b7*a9*a8-22386*b7*b6*a5+16695*b7*b5*a6+14315*b7*a9*a4+483*b7*b6*a4+10437*b7*b6*a9-2968*b7*b10*a8-6741*b5*b9*a8+6391*b7*b6*a8-4648*b7*a9^2-18102*b7*a5^2+8764*b7*b10*a5+7308*b5*b9*a9+4480*b7*a8^2-2919*a5*b7*a4-1638*b8^2*a9+4662*b8^2*a8+252*b9*b5*a4-34272*b5*a9*a3+2709*b5*b8*a9-1827*b5*b4*a9+3654*b8^2*a4+2205*b9*a9*a3+14868*b9*a5*a3-7497*b9*b8*a9-5103*b9*b8*a8+6615*b9*b4*a5+28413*b9*b8*a5+11781*b9*b4*a9+40824*b5*a8*a3-10395*b5*a8*a7+10962*b5*a5*a3-2583*b9*a5*a7+15372*b9*a8*a3-18081*b9*a4*a3-2709*b9*a9*a7+6363*b9*a4*a7-4347*b9*a8*a7-1008*b9*b7*a10+4193*b9*b7*a6"),
		SMQP("-5292*b10^2*a6+1785*b6*a10*a5-5383*a6^2*b8-623*a6*b10*a9+973*a6*b6*a8-2464*a6*b10*a8+16632*b6*b10*a6-7182*b6*b10*a10+924*a6*b6*a5+5334*a6*b6*a9-11613*a6*b6*a4+560*a6*a10*a3-1225*a6*a10*a7+5747*a6*a5*a9+7000*a6*b10*a4-2996*a6*a5*a8-9534*a6^2*b4+2128*a6*b5*a10+13741*a6*b4*a10-6755*a6*b9*a10+2219*a6*a9*a4-2457*a6^2*b5-5187*a6*a9*a8+7357*a6*a8*a4+854*b10*a10*a9+826*a6*b10*a5+4914*b6^2*a10+1218*a10*a9*a8-3885*a10*b6*a9-3409*a10^2*b4-252*a10^2*b8-1645*a10^2*b5+1792*b10*a10*a8-1435*a10*b6*a8+7581*a10*b6*a4-1505*a10*a9*a4-2443*a10*a5*a8-182*a10*a5*a9-2779*a10*a8*a4+3738*a10*a5*a4-4312*a10*b10*a4-245*a10^2*a3+805*a10*b8*a6+2142*a10^2*b9-1225*a10*b10*a5+378*a10^2*a7+14*a6*a9^2-7623*a6*a4^2+9107*a6^2*b9+651*a6^2*a3+1078*a6*a8^2-1470*a6*a5*a4+2268*b10^2*a10-11340*b6^2*a6+1393*a6^2*a7+700*a10*a9^2-2898*a6*a5^2+441*a10*a5^2-28*a10*a8^2+3087*a10*a4^2"),
		SMQP("203*b7*b4*a5+427*b7*b9*a5-413*b7*b4*a9-721*b7*b8*a5-147*b7^2*a6+406*b7*b8*a9-1456*b7*b5*a9+112*b7^2*a10+861*a7*b7*a8-56*a7*b7*a5-861*a7*b7*a4-3696*a7^2*a3-3360*a3^3+427*b7*b5*a8+483*b7*a9*a7+588*b7*b5*a4-2296*b7*a8*a3+1211*b7*a4*a3-91*b7*a5*a3+189*b7*a9*a3-77*b7*b9*a4+154*b7*b8*a8+35*b7*b9*a8-35*b7*b9*a9+672*a7^3+6384*a7*a3^2-224*b7*b8*a4+112*b7*b4*a4+875*b7*b5*a5"),
		SMQP("3542*a5*b5*a8+511*b4*a5*a4+651*b9*b6*a9-595*b10*b4*a5+644*b7*a10*a5-371*b10*b9*a9-1288*b8*a6*a3-252*b8*a9*a4-5264*a6*a3*a7+1932*b7*a10*a4-2331*b7*a6*a4-1344*b7*a10*a9-1323*b6*a4*a7-49*b6*a8*a7+1995*b6*a5*a3+6923*a5*a8*a7-1463*a4*a3*a8+1869*b6*a4*a3+1267*b8*b6*a5-2324*b6*a8*a3+1036*a8^2*a3-2814*b5*a4*a9+385*b6*b5*a5-1680*b6*a5*a7-763*b7*a6*a5-42*b8*a8*a9-42*b8*a8^2+2331*a6*b7*a8-910*a4*a7*a8+210*b8*a8*a4-735*b5*a8^2+287*a8^2*a7-504*b8*a10*a3-1043*a7*a9^2-21*a7*a4^2+1008*b5*a4^2+399*b4*a9^2-714*a8*b9*a9+273*a8*b4*a9-868*a10*a3^2-700*b4*a10*a7-588*a8*a9*a7+693*a8*b5*a9+1575*a8*a9*a3-6279*a7*a5*a4+2072*a7^2*a6+252*b9*a10*a7+1372*a10*a3*a7-672*a8*b7*a10-1666*a9*a4*a3+1806*a9^2*b5+1106*a9*b9*a4+2212*a9*a4*a7-2891*a9*a5*a7-385*a9^2*b9+2198*a9*a5*a3-672*a7^2*a10+49*b6*b9*a8+819*a4^2*a3+3864*a3^2*a6+1288*b4*a10*a3+385*b10*b9*a4+8106*a3*a5*a4-924*b7*b6*a10-1617*b7*b6*a6-2709*a5^2*a3+84*b8*b6*a9+84*b8*b10*a9-567*b4*b6*a9-609*a4*b4*a9+2100*b8*b5*a10+2100*b8*b9*a10+987*b7*b10*a6-1533*b9*a4^2-161*b10*b5*a5+2100*b5^2*a6+2016*b4*a6*a3+840*b8^2*a6+700*b8*a6*a7-1015*b10*a9*a7+329*b10*a4*a7+343*b10*a8*a7+1470*b10*b5*a9-448*b9*a6*a7+3892*b4*b9*a10-616*b9*a6*a3+476*b9*b8*a6-756*b5^2*a10-5460*b4*b5*a10+1085*b8*a5*a9-231*b10*b5*a8-938*b4*a5*a9-504*b10*b5*a4+392*b9*a10*a3-2604*b4*b9*a6+175*b10*a9*a3-245*b10*a5*a3-42*b8*b6*a8-693*b9*b10*a5+63*b6*a9*a3+2121*b6*a9*a7+644*b10*a8*a3-77*b8*b10*a5+420*b4*b5*a6-301*a9^2*a3-623*b10*a4*a3-1722*a5*b9*a9-924*b4*a6*a7+2471*a5^2*b4+777*b7*a6*a9-343*b9*b10*a8+1673*b9*b6*a5+21*b9*b6*a4-959*b4*a5*a8-343*b8*a5*a4-1848*b8^2*a10+1624*b9*a8*a4-140*b9*a5*a4+126*b8*a9^2-5817*b8*a5^2+224*b10*a5*a7-343*b4*b6*a5+1680*b5*a6*a3-672*b5*a6*a7-1092*b5*b9*a6+231*b5*a8*a4+1372*b5*a10*a3+56*b5*a10*a7+448*b5*b9*a10-2947*b5*a5*a4-2268*b5*b8*a6-217*a5*b5*a9-287*b9*a8^2-1708*a5*b9*a8+2597*a5*b8*a8-252*b8*a10*a7+189*b10*b4*a9+3311*b9*a5^2+756*b5*b6*a4+392*b9^2*a6-1596*b9^2*a10+945*b5*b6*a8-42*b8*b10*a8+777*a5^2*b5-3360*b5*b6*a9-9317*a5*a3*a8+3528*a5^2*a7"),
		SMQP("-3276*b10^2*a6+1659*b6*a10*a5-581*a6^2*b8-889*a6*b10*a9+4571*a6*b6*a8-1820*a6*b10*a8+13104*b6*b10*a6-6426*b6*b10*a10+588*a6*b6*a5+2478*a6*b6*a9-10815*a6*b6*a4-2324*a6*a10*a3-371*a6*a10*a7-4823*a6*a5*a9+4844*a6*b10*a4+9548*a6*a5*a8-8190*a6^2*b4+3332*a6*b5*a10+7427*a6*b4*a10-3913*a6*b9*a10+6097*a6*a9*a4+777*a6^2*b5-861*a6*a9*a8-2317*a6*a8*a4+826*b10*a10*a9+14*a6*b10*a5+3654*b6^2*a10+798*a10*a9*a8-2667*a10*b6*a9-2975*a10^2*b4+252*a10^2*b8-1211*a10^2*b5+1064*b10*a10*a8-1421*a10*b6*a8+4263*a10*b6*a4-1687*a10*a9*a4-1757*a10*a5*a8-658*a10*a5*a9-1589*a10*a8*a4+3486*a10*a5*a4-2072*a10*b10*a4+77*a10^2*a3-1309*a10*b8*a6+1386*a10^2*b9-1127*a10*b10*a5+126*a10^2*a7-1022*a6*a9^2-3465*a6*a4^2+5509*a6^2*b9-1407*a6^2*a3-1834*a6*a8^2-7350*a6*a5*a4+1764*b10^2*a10-8820*b6^2*a6+3395*a6^2*a7+644*a10*a9^2+3402*a6*a5^2+483*a10*a5^2+196*a10*a8^2+3213*a10*a4^2"),
		SMQP("-231*b7*b4*a5+441*b7*b9*a5-63*b7*b4*a9+189*b7*b8*a5+315*b7^2*a6-798*b7*b8*a9+1176*b7*b5*a9+231*a7*b7*a8+1708*a7*b7*a5-399*a7*b7*a4-252*b7*b4*a8-63*b7*b5*a8-231*b7*a9*a7+84*b7*b5*a4-504*b7*a8*a3+1113*b7*a4*a3+147*b7*a5*a3-609*b7*a9*a3-399*b7*b9*a4+462*b7*b8*a8+105*b7*b9*a8+231*b7*b9*a9-2079*b7*b5*a5"),
		SMQP("168*a6*a3*a4+294*a6*a3*a5+420*a6*b5*b6+189*a6*b5*a5-420*a6*a3*b6-210*a6*b4*a5-105*a6*b5*a4-294*a4*b6*a5-126*a6*b4*a4+273*a4^2*a5-336*a4*b6^2-105*a4*a5^2-84*b6^2*a5+462*b6*a5^2-378*a4^3+714*a4^2*b6-252*a5^3+98*b10*a5^2-273*b10*b5*a6-14*b10*b6*a9+252*b10*a6*a3-672*b6*a8*a4-14*b6*b10*a4+14*b6*b10*a8-364*a6*a4*a7-406*b6*a10*a3+364*a6*b9*a4-504*a6*a8*a3+14*b10*a9*a4+364*b6*a6*a7+42*a6*b9*a9-21*a6*b9*a8+42*b6*a9*a8-14*b10*a8*a4+21*b10*b9*a6+406*b5*b6*a10-21*a6*b4*a9+21*a6*b8*a8-28*b6*a9^2-14*b6*a8^2+350*a8*a4^2+70*a9^2*a5+14*a8^2*a4+252*a6*a9*a3+28*a9*a4^2-42*a9*a8*a4-406*a10*b5*a4+406*a10*a4*a3-140*b6*b10*a5-287*b8*a6*a4-70*a9*b6*a4+42*a8*b4*a6-35*a5*b10*a8+14*b10*a4^2-637*a5*b6*a9+350*a4*b4*a10+28*a9^2*a4+175*a5*b10*a4+322*b6^2*a8+42*b6^2*a9-42*b8*a6*a9+35*a5*b10*a9-21*b8*b10*a6-616*a5^2*a8+448*a5*b9*a6-364*b9*b6*a6+623*b4*a10*a5+35*a5*a8^2-231*a6*b5*a9+462*a6*b5*a8-63*a6^2*b7+84*b4*b10*a6-133*a5*b5*a10+133*a5*a10*a3+602*a5*a9*a4+287*b8*b6*a6-105*a5*a9*a8-1085*b8*a6*a5-350*b4*b6*a10-581*a5*a8*a4+651*a5*b6*a8+448*a5^2*a9-154*a5*a6*a7"),
		SMQP("588*b10^2*a6-441*b6*a10*a5+231*a6^2*b8+399*a6*b10*a9+567*a6*b6*a8-84*a6*b10*a8-2016*b6*b10*a6+798*b6*b10*a10-252*a6*b6*a5-1134*a6*b6*a9+1449*a6*b6*a4-84*a6*a10*a3+217*a6*a10*a7-735*a6*a5*a9-756*a6*b10*a4+714*a6^2*b4-714*a6*b5*a10-819*a6*b4*a10+665*a6*b9*a10-315*a6*a9*a4+609*a6^2*b5+147*a6*a9*a8-441*a6*a8*a4-224*b10*a10*a9+42*a6*b10*a5-378*b6^2*a10-84*a10*a9*a8+567*a10*b6*a9+259*a10^2*b4-84*a10^2*b8+175*a10^2*b5-70*b10*a10*a8+133*a10*b6*a8-945*a10*b6*a4+161*a10*a9*a4+7*a10*a5*a8+350*a10*a5*a9+133*a10*a8*a4-252*a10*a5*a4+406*a10*b10*a4+35*a10^2*a3+35*a10*b8*a6-126*a10^2*b9+217*a10*b10*a5-42*a10^2*a7+210*a6*a9^2+567*a6*a4^2-819*a6^2*b9-63*a6^2*a3+126*a6*a8^2+378*a6*a5*a4-252*b10^2*a10+1260*b6^2*a6-189*a6^2*a7-112*a10*a9^2+378*a6*a5^2-189*a10*a5^2-14*a10*a8^2-189*a10*a4^2"),
		SMQP("2709*b9^2*a9+3591*b9^2*a4+4599*b9^2*a8+1372*b7*b4*a10-4592*b7*b10*a9+20552*b7*a10*a3+5460*b7*b4*a6+952*b7*b10*a4+35532*a9*a3*a7+11025*b4*a4*a3+16380*b8*a8*a3-2961*b8*b9*a4+3780*b8*a4*a3+6237*b8*a9*a7-30030*b7*a6*a3-12159*b4*a8*a3-3024*b4*a8*a7+12033*b5^2*a5+13356*b4*b5*a4+4858*a7*b7*a6+33894*a7*a8*a3+1638*a7*b5*a9-7812*a7*b5*a5-2520*b8*b7*a10-19278*a7*b4*a5-8820*a7*b7*a10+29358*a7*b4*a9+18900*a7*b8*a5-11277*a7*b8*a8+1260*b8*a9*a3-46242*b8*a5*a3-2016*b8*b4*a9-1260*b8*b4*a5-4095*b8*b5*a4-14112*a7^2*a5-77616*a7*a4*a3+19530*a7^2*a4+504*b4*b8*a8-4662*b4*b9*a8-5796*b4*a4*a7+19152*b4*b5*a8-5922*a7^2*a9-2205*a7*b8*a4-21546*a9*a3^2-7875*b7*a4^2+553*b7*a4*a8-9324*b8^2*a5+37548*a3*b4*a5+79002*a3^2*a4-10962*a3^2*a5-44667*a3*b4*a9+245*b7*a5*a9-16254*b9^2*a5+7483*b7*a5*a8-29484*b5*a4*a3-4676*b7*b5*a10+19656*b5*a4*a7-9387*b5^2*a8+27342*a5*a3*a7+5040*b5^2*a9-7686*a7^2*a8-13986*b4*b9*a4-57456*a3^2*a8-13041*b5*b9*a5+6867*b5*b8*a5-28539*b5*b4*a5-4032*b4^2*a9-4977*b5*b8*a8-2800*b7*b8*a6-504*b4^2*a5-5880*b7*a9*a8-22890*b7*b6*a5+17451*b7*b5*a6+14315*b7*a9*a4+483*b7*b6*a4+10437*b7*b6*a9-2968*b7*b10*a8-5292*b5*b9*a8+6643*b7*b6*a8-4648*b7*a9^2-17262*b7*a5^2+8764*b7*b10*a5+6552*b5*b9*a9+4480*b7*a8^2-2415*a5*b7*a4-1638*b8^2*a9+5166*b8^2*a8-252*b9*b5*a4-35028*b5*a9*a3+2205*b5*b8*a9-1071*b5*b4*a9+3654*b8^2*a4+3339*b9*a9*a3+15372*b9*a5*a3-7623*b9*b8*a9-5733*b9*b8*a8+6552*b9*b4*a5+28098*b9*b8*a5+11718*b9*b4*a9+38997*b5*a8*a3-10332*b5*a8*a7+8379*b5*a5*a3-2520*b9*a5*a7+14049*b9*a8*a3-17703*b9*a4*a3-2835*b9*a9*a7+6363*b9*a4*a7-4473*b9*a8*a7-1008*b9*b7*a10+4256*b9*b7*a6"),
		SMQP("336*b9^2*a9+1764*b9^2*a4+1428*b9^2*a8+2408*b7*b4*a10-2380*b7*b10*a9+8932*b7*a10*a3-168*b7*b4*a6+812*b7*b10*a4+20496*a9*a3*a7+6363*b4*a4*a3+9828*b8*a8*a3-756*b8*b9*a4-189*b8*a4*a3+1974*b8*a9*a7-11046*b7*a6*a3-6405*b4*a8*a3-1176*b4*a8*a7+4032*b5^2*a5+6300*b4*b5*a4+4550*a7*b7*a6+23814*a7*a8*a3+7875*a7*b5*a9-1680*a7*b5*a5-2520*b8*b7*a10-7770*a7*b4*a5-5292*a7*b7*a10+14658*a7*b4*a9+6552*a7*b8*a5-4998*a7*b8*a8+357*b8*a9*a3-20853*b8*a5*a3-105*b8*b4*a9-1365*b8*b4*a5-1953*b8*b5*a4+1344*a7^2*a5-47124*a7*a4*a3+12558*a7^2*a4+504*b4*b8*a8-714*b4*b9*a8-3780*b4*a4*a7+7224*b4*b5*a8-6006*a7^2*a9+168*a7*b8*a4-11466*a9*a3^2-6237*b7*a4^2+5831*b7*a4*a8-777*b8^2*a5+18648*a3*b4*a5+45738*a3^2*a4+7182*a3^2*a5-21861*a3*b4*a9+2023*b7*a5*a9-5943*b9^2*a5+6265*b7*a5*a8-6867*b5*a4*a3+3416*b7*b5*a10+3465*b5*a4*a7+252*b5^2*a4-672*b5^2*a8-7686*a5*a3*a7-2016*b5^2*a9-4578*a7^2*a8-7182*b4*b9*a4-34272*a3^2*a8-2919*b5*b9*a5-147*b5*b8*a5-5796*b5*b4*a5-2520*b4^2*a9-2142*b5*b8*a8+1057*b7*b8*a6-3444*b7*a9*a8-7098*b7*b6*a5-3276*b7*b5*a6+8113*b7*a9*a4+273*b7*b6*a4+4683*b7*b6*a9-1148*b7*b10*a8-1428*b5*b9*a8+1589*b7*b6*a8-3248*b7*a9^2-7182*b7*a5^2+3164*b7*b10*a5+4599*b5*b9*a9-364*b7*a8^2-5229*a5*b7*a4-1764*b8^2*a9+1260*b8^2*a8-1575*b9*b5*a4-20349*b5*a9*a3+2961*b5*b8*a9-7812*b5*b4*a9+1890*b8^2*a4+2310*b9*a9*a3+8253*b9*a5*a3-2310*b9*b8*a9-2478*b9*b8*a8-609*b9*b4*a5+8652*b9*b8*a5+8337*b9*b4*a9+14091*b5*a8*a3-5733*b5*a8*a7+1764*b5*a5*a3-1764*b9*a5*a7+8799*b9*a8*a3-10962*b9*a4*a3-2394*b9*a9*a7+4704*b9*a4*a7-3402*b9*a8*a7-1008*b9*b7*a10+2821*b9*b7*a6"),
		SMQP("42*a5*b5*a8-1092*b4*a5*a4-364*b10*b4*a5-252*b7*a6*a4+273*b8*b6*a5-42*b6*b5*a5+63*b7*a6*a5+756*b4*a4^2-56*b4*a9^2+84*a8*b4*a9-308*b4*a10*a3-700*b4^2*a10-672*b4*b6*a4+252*b7*b6*a6+84*b4*b6*a9-56*a4*b4*a9+21*b10*b5*a5-1092*b4*a6*a3+504*b8^2*a6+252*b8*a6*a7+28*b4*b10*a8-252*b9*a6*a7+252*b9*a6*a3-756*b9*b8*a6+616*b4*b8*a6+308*b4*b5*a10+42*b8*a5*a9-1547*b4*a5*a9-1148*b4*b9*a6-21*b9*b10*a5+644*b4*b6*a8-700*b4*a8*a4+21*b8*b10*a5+1764*b4*b5*a6+84*b4^2*a6-42*a5*b9*a9+728*b4*a6*a7+609*a5^2*b4-210*b9*b6*a5+2114*b4*a5*a8-273*b8*a5*a4+210*b9*a5*a4-28*b4*b10*a4+315*b8*a5^2-126*b4*b6*a5-28*b4*a8^2-756*b5*a6*a3+252*b5*a6*a7+252*b5*b9*a6-21*b5*a5*a4-252*b5*b8*a6-21*a5*b5*a9+21*a5*b9*a8-21*a5*b8*a8-28*b10*b4*a9-147*b9*a5^2+252*b9^2*a6-126*a5^2*b5"),
		SMQP("-770*a6^2*b8+35*a6*b10*a9+203*a6*b6*a8-35*a6*b10*a8-252*b6*b10*a6-294*a6*b6*a5-105*a6*b6*a9+336*a6*b6*a4+133*a6*a10*a3+448*a6*a5*a9+287*a6*b10*a4-616*a6*a5*a8-357*a6^2*b4-133*a6*b5*a10+623*a6*b4*a10+70*a6*a9*a4+315*a6^2*b5-105*a6*a9*a8-133*a6*a8*a4+98*a6*b10*a5+252*b6^2*a10-504*a10*b6*a4+70*a6*a9^2-441*a6*a4^2+301*a6^2*b9+42*a6^2*a3+35*a6*a8^2+651*a6*a5*a4-154*a6^2*a7-252*a6*a5^2+252*a10*a4^2"),
		SMQP("-5565*a6*a3*a4+3717*a6*a3*a5-840*a6*b5*b6+1575*a6*b5*a5+1680*a6*a3*b6-840*a6*b4*b6-3255*a6*b4*a5+5250*a6*b5*a4+2436*a4*b6*a5-1260*a6*b4*a4-420*a4^2*a5+2352*a4*b6^2+2016*a4*a5^2-672*b6^2*a5-2520*b6*a5^2+252*a4^3-3444*a4^2*b6-4200*b10*a5^2-2730*b10*b5*a6+1260*b8*a10*a5-364*b10*b6*a9-1806*b10*b5*a10-42*b10*b9*a10-3360*b10*a6*a7+3024*b10*a6*a3+4312*b6*a8*a4+3164*b6*b10*a4-1323*a6*a8*a7-742*a9*b5*a10-2996*b6*b10*a8+1043*a6*a9*a7+707*a6*a4*a7+2044*b6*a10*a3-77*a6*b9*a4+1512*a6*a8*a3-476*b10*a9*a8+448*b10*a9*a4-2632*b6*a6*a7+308*b10*a9^2-539*a6*b9*a9+1113*a6*b9*a8+224*b6*a9*a8+3052*b10*a8*a4+3066*b10*b9*a6-2128*b5*b6*a10+189*a6*b4*a9-1428*a6*b8*a8+364*b6*a9^2+168*b9*a10*a8-868*b6*a8^2+84*a8*a9^2-1540*a8*a4^2+924*a8*b4*a10+210*a8*b5*a10-1904*a9^2*a5-126*a8*a10*a7+1148*a8^2*a4+2373*a6*a9*a3-56*a9^3+1036*a9*a4^2-28*a9*a8^2-1120*a9*a8*a4+574*a9*a10*a3-532*a9*b4*a10+126*a10*a9*a7+126*a10*a4*a7-574*a10*b5*a4+42*a10*b9*a4+154*a10*a4*a3-294*a10*b8*a8-210*a10*b9*a9+1904*b6*b10*a5+1162*b8*a6*a4-812*a9*b6*a4+420*a8*b4*a6-9240*a5*b10*a8-3388*b10*a4^2+168*b10*a8^2+1680*b10^2*a5+7168*a5*b6*a9-42*a4*b8*a10-364*a4*b4*a10-168*b9*a10*a5+224*a9^2*a4+5264*a5*b10*a4+168*b10^2*a4-2212*b6^2*a8-420*b6^2*a9+336*b8*a10*a9+826*b8*a6*a9+6944*a5*b10*a9+168*b10^2*a9-1008*a5*a10*a7-2898*b8*b10*a6-1176*a5^2*a8-567*a5*b9*a6+1792*b9*b6*a6-840*b4*a10*a5-2184*a5*a8^2-2562*a6*b5*a9-1617*a6*b5*a8-945*a6^2*b7+336*b4*b10*a6-168*b10^2*a8+1848*b10*a10*a3+210*a5*b5*a10+546*a5*a10*a3-2492*a5*a9*a4-2282*b8*b6*a6+4424*a5*a9*a8+42*b8*b10*a10+2121*b8*a6*a5+2576*b4*b6*a10+42*b8*b6*a10+3024*b4*b10*a10+1904*a5*a8*a4-9184*a5*b6*a8+672*a5^2*a9+84*b9*b6*a10-3696*a5*a6*a7"),
		SMQP("1617*b7*b4*a5+2457*b7*b9*a5-1743*b7*b4*a9-4011*b7*b8*a5-945*b7^2*a6+1890*b7*b8*a9-6048*b7*b5*a9+672*b7^2*a10+2751*a7*b7*a8+504*a7*b7*a5-2415*a7*b7*a4+336*b7*b4*a8-11088*a7^2*a3-10080*a3^3+4137*b7*b5*a8+1281*b7*a9*a7-252*b7*b5*a4-8904*b7*a8*a3+3801*b7*a4*a3-2905*b7*a5*a3+2079*b7*a9*a3-63*b7*b9*a4+126*b7*b8*a8-63*b7*b9*a8-609*b7*b9*a9+2016*a7^3+19152*a7*a3^2+3465*b7*b5*a5"),
		SMQP("-14*a6*b5*b6-42*a6*b5*a5+14*a6*b4*b6-301*a6*b4*a5-7*a6*b5*a4+112*a6*b4*a4+7*b10*b5*a6+70*a6*b9*a4-14*a6*b9*a9+7*a6*b9*a8-7*b10*b9*a6-161*a6*b4*a9-7*a6*b8*a8-91*b8*a6*a4+406*a8*b4*a6-168*a4*b4*a10+14*b8*a6*a9+7*b8*b10*a6-49*a5*b9*a6-70*b9*b6*a6-7*a6*b5*a9+14*a6*b5*a8-63*a6^2*b7-196*b4*b10*a6+91*b8*b6*a6+105*b8*a6*a5+168*b4*b6*a10"),
		SMQP("777*b7*b4*a5+945*b7*b9*a5-1911*b7*b4*a9-2163*b7*b8*a5-441*b7^2*a6+1890*b7*b8*a9-7056*b7*b5*a9+336*b7^2*a10+4935*a7*b7*a8-1512*a7*b7*a5-5271*a7*b7*a4-336*b7*b4*a8-22176*a7^2*a3-20160*a3^3-735*b7*b5*a8+3129*b7*a9*a7+6804*b7*b5*a4-11256*b7*a8*a3+5649*b7*a4*a3-1169*b7*a5*a3-441*b7*a9*a3-567*b7*b9*a4+1134*b7*b8*a8+441*b7*b9*a8+231*b7*b9*a9+4032*a7^3+38304*a7*a3^2-2016*b7*b8*a4+1008*b7*b4*a4+6993*b7*b5*a5"),
		SMQP("1260*b10^2*a6-126*b6*a10*a5-574*a6^2*b8+3073*a6*b10*a9+1267*a6*b6*a8-4459*a6*b10*a8+504*b6*b10*a6-1260*b6*b10*a10+798*a6*b6*a5-1911*a6*b6*a9-462*a6*b6*a4-217*a6*a10*a3+1806*a6*a10*a7-910*a6*a5*a9-1085*a6*b10*a4-686*a6*a5*a8-1785*a6^2*b4-1295*a6*b5*a10+1603*a6*b4*a10-2562*a6*b9*a10-154*a6*a9*a4+1449*a6^2*b5-2667*a6*a9*a8+1351*a6*a8*a4-840*b10*a10*a9+1624*a6*b10*a5-504*b6^2*a10+1008*a10*a9*a8-1260*a10*b6*a9-4116*a10^2*b4-504*a10^2*b8+1428*a10^2*b5+840*b10*a10*a8+3192*a10*b6*a8-126*a10*b6*a4+2100*a10*a9*a4+6720*a10*a5*a8-4704*a10*a5*a9-4620*a10*a8*a4-4788*a10*a5*a4+672*a10*b10*a4-1428*a10^2*a3+4872*a10*b8*a6+504*a10^2*b9-840*a10*b10*a5+1106*a6*a9^2-1575*a6*a4^2+623*a6^2*b9+336*a6^2*a3+1687*a6*a8^2+483*a6*a5*a4-140*a6^2*a7-672*a10*a9^2+756*a6*a5^2+3150*a10*a5^2-336*a10*a8^2+2394*a10*a4^2"),
		SMQP("399*b7*b4*a5+567*b7*b9*a5-609*b7*b4*a9-1029*b7*b8*a5-231*b7^2*a6+630*b7*b8*a9-2352*b7*b5*a9+168*b7^2*a10+1281*a7*b7*a8-168*a7*b7*a5-1281*a7*b7*a4-5544*a7^2*a3-5040*a3^3+735*b7*b5*a8+735*b7*a9*a7+924*b7*b5*a4-3528*b7*a8*a3+1743*b7*a4*a3-343*b7*a5*a3+441*b7*a9*a3-105*b7*b9*a4+210*b7*b8*a8+63*b7*b9*a8-63*b7*b9*a9+1008*a7^3+9576*a7*a3^2-336*b7*b8*a4+168*b7*b4*a4+1407*b7*b5*a5"),
		SMQP("-32949*a6*a3*a4-16779*a6*a3*a5-100884*a6*b5*b6+16695*a6*b5*a5+20832*a6*a3*b6-10668*a6*b4*b6-18291*a6*b4*a5-24990*a6*b5*a4-18480*a4*b6*a5-19488*a6*b4*a4+22176*a4^2*a5+17472*a4*b6^2-72744*a4*a5^2+3360*b6^2*a5-39984*b6*a5^2-20664*a4^3-12096*a4^2*b6+30240*a5^3-38080*b10*a5^2+29778*b10*b5*a6-20664*b8*a10*a5-784*b10*b6*a9-28728*b10*b5*a10+4200*b10*b9*a10-21504*b10*a6*a7+20160*b10*a6*a3+22568*b6*a8*a4+17360*b6*b10*a4-34419*a6*a8*a7+5152*a9*b5*a10-30800*b6*b10*a8+11683*a6*a9*a7+14259*a6*a4*a7+3024*b10*a10*a7-10080*b6*a10*a7+11536*b6*a10*a3+10983*a6*b9*a4+55944*a6*a8*a3-3640*b10*a9*a8-7616*b10*a9*a4-4816*b6*a6*a7+6664*b10*a9^2-29736*a6*b7*a10+5537*a6*b9*a9-32067*a6*b9*a8+15064*b6*a9*a8+12600*b10*a8*a4+15918*b10*b9*a6+21056*b5*b6*a10+21693*a6*b4*a9+3696*a8^3-13944*a6*b8*a8-16520*b6*a9^2+9912*b9*a10*a8-11200*b6*a8^2+9576*a8*a9^2+16968*a8*a4^2-4536*a8*b4*a10-17808*a8*b5*a10-41328*a9^2*a5+7056*a8*a10*a7-3528*a8^2*a4-9744*a8*a10*a3-5523*a6*a9*a3-1456*a9^3+19320*a9*a4^2-5768*a9*a8^2-23408*a9*a8*a4-14056*a9*a10*a3-30464*a9*b4*a10+5040*a10*a9*a7-10080*a10*a4*a7+20160*a10*b5*a4+5376*a10*b9*a4+20328*a10*a4*a3+9240*a10*b8*a8+6048*a10^2*b7+3360*a10*b9*a9+14336*b6*b10*a5+1134*b8*a6*a4+5656*a9*b6*a4+15540*a8*b4*a6-57848*a5*b10*a8-14280*b10*a4^2+6048*b10*a8^2+6720*b10^2*a5+33208*a5*b6*a9-6888*a4*b8*a10+17304*a4*b4*a10+1176*b9*a10*a5+3248*a9^2*a4+36344*a5*b10*a4+672*b10^2*a4-7168*b6^2*a8+2352*b6^2*a9-4368*b8*a10*a9+11690*b8*a6*a9+32088*a5*b10*a9+672*b10^2*a9-1008*a5*a10*a7-37422*b8*b10*a6+94640*a5^2*a8-28091*a5*b9*a6+39508*b9*b6*a6-24976*b4*a10*a5-32872*a5*a8^2+28182*a6*b5*a9-24297*a6*b5*a8+25011*a6^2*b7+5880*b4*b10*a6+2352*b10^2*a8+7392*b10*a10*a3+20216*a5*b5*a10-16184*a5*a10*a3+17416*a5*a9*a4+20818*b8*b6*a6+55720*a5*a9*a8+10920*b8*b10*a10+91021*b8*a6*a5+28448*b4*b6*a10-24360*b8*b6*a10+17472*b4*b10*a10-25816*a5*a8*a4+17416*a5*b6*a8-39872*a5^2*a9-18480*b9*b6*a10+20888*a5*a6*a7"),
		SMQP("-37149*a6*a3*a4-9891*a6*a3*a5-85428*a6*b5*b6-15057*a6*b5*a5-5292*a6*b4*b6-36099*a6*b4*a5-2646*a6*b5*a4-43680*a4*b6*a5-4368*a6*b4*a4-1344*a4^2*a5-24696*a4*a5^2-3024*b6*a5^2-11592*a4^3-336*a4^2*b6-6048*a5^3-2016*b10*a5^2+28602*b10*b5*a6-23184*b8*a10*a5-21504*b10*b5*a10+4368*b10*b9*a10-4032*b10*a6*a7+2632*b6*a8*a4-37387*a6*a8*a7-5208*a9*b5*a10-14112*b6*b10*a8-1701*a6*a9*a7+23947*a6*a4*a7+3024*b10*a10*a7-10080*b6*a10*a7+12383*a6*b9*a4+48552*a6*a8*a3-3304*b10*a9*a8-6608*b10*a9*a4+15120*b6*a6*a7+7560*b10*a9^2-30240*a6*b7*a10-9639*a6*b9*a9+6349*a6*b9*a8-4200*b6*a9*a8-3752*b10*a8*a4-16506*b10*b9*a6+53088*b5*b6*a10+48573*a6*b4*a9+4256*a8^3-20048*a6*b8*a8-7560*b6*a9^2+9744*b9*a10*a8-17248*b6*a8^2+3976*a8*a9^2+10696*a8*a4^2+5768*a8*b4*a10-20272*a8*b5*a10+5040*a9^2*a5+7056*a8*a10*a7+2520*a8^2*a4-7616*a8*a10*a3+693*a6*a9*a3+3024*a9^3-2632*a9*a4^2-5208*a9*a8^2+3136*a9*a8*a4-3528*a9*a10*a3-4872*a9*b4*a10+5040*a10*a9*a7-10080*a10*a4*a7+952*a10*b5*a4-10416*a10*b9*a4+19544*a10*a4*a3+9408*a10*b8*a8+6048*a10^2*b7+3696*a10*b9*a9-154*b8*a6*a4+24360*a9*b6*a4-12012*a8*b4*a6-13384*a5*b10*a8+3976*b10*a4^2+4816*b10*a8^2-10584*a5*b6*a9-10752*a4*b8*a10+13384*a4*b4*a10+2352*b9*a10*a5-10192*a9^2*a4+16072*a5*b10*a4+15120*b6^2*a8-4704*b8*a10*a9-4662*b8*a6*a9+5544*a5*b10*a9-1008*a5*a10*a7-25830*b8*b10*a6+13104*a5^2*a8-25515*a5*b9*a6+24948*b9*b6*a6+17976*b4*a10*a5-39928*a5*a8^2-882*a6*b5*a9+35343*a6*b5*a8+31059*a6^2*b7-7560*b4*b10*a6+3024*b10^2*a8+6552*a5*b5*a10-1512*a5*a10*a3+26600*a5*a9*a4+54810*b8*b6*a6+952*a5*a9*a8+10752*b8*b10*a10+48573*b8*a6*a5+3360*b4*b6*a10-20496*b8*b6*a10+1344*b4*b10*a10+28840*a5*a8*a4+50232*a5*b6*a8-4032*a5^2*a9-2688*b9*b6*a10+8568*a5*a6*a7"),
		SMQP("-14*a3*b7*b5+7*b7*b5^2+7*a3^2*b7"),
		SMQP("7*b7*a5*a3-7*b7*b5*a5"),
		SMQP("35*b7*b4*a5-42*b7*b5*a4+42*b7*a4*a3-14*b7*a5*a3-14*a7*a3*b8-28*a3*b8^2+42*b7*b5*b6-14*b5^2*b8-42*a3*b7*b6+42*b9*a3*b8-7*a7*b4*b9-14*a3*b5*a7+14*b8*b5*a7+42*a3^2*b5+14*b4*b8^2+14*a3*a7*b9-14*a3*b9^2-14*a3^2*b9-14*a7*b5*b9+7*a7*b4*b5+14*a3^2*b4+14*a7*b5^2-14*b4^2*b9-42*b8*b5*b9+7*b4*b9^2+42*b4*b5^2-42*a3*b5^2+7*b4^2*a3+14*b9^2*b5+14*a3*b8*b5+7*b5*b9*b4+7*a3*b4*b9+28*b8^2*b5+7*b4*a7*b8-77*a3*b4*b5+7*b7*b4*b6+14*b4^2*b5-7*b7*b4*a4+14*b5^2*b9-21*b4*b8*b5+14*b4*a3*b8-21*b4*b8*b9+14*b7*b5*a5"),
		SMQP("504*a7^2*b7-378*a7*b7*b4+168*b7^2*a9-1386*a3*a7*b7+1050*a3*b7*b4-1113*a3*b7*b8+525*b7*b4*b5-378*a3^2*b7+161*b7^2*a5-420*b7*b4*b9-399*a3*b7*b9+21*b7*b5^2-294*b7*b8*b5+189*a7*b7*b9+126*b7*b8^2+84*b7^2*b10+252*b7^2*a8-357*b7^2*a4-21*b7*b9^2-252*b7*b8*b4-231*b7*b5*b9+315*b7*b8*b9+42*b7*b4^2+3528*a3*b7*b5-1071*a7*b7*b5+399*a7*b7*b8-63*b7^2*b6"),
		SMQP("1351*a6*b5*b6+42*a6*b5*a5-8911*a6*b4*b6-1456*a6*b4*a5-721*a6*b5*a4-2114*a6*b4*a4-3780*b6^3+952*b10*b5*a6+210*b8*a10*a5-833*b10*b5*a10+875*b10*b9*a10-217*a9*b5*a10-140*a6*b9*a4-735*a6*b7*a10-476*a6*b9*a9-371*a6*b9*a8-1120*b10*b9*a6+1204*b5*b6*a10+1960*a6*b4*a9-133*a6*b8*a8+49*b9*a10*a8+616*a8*b4*a10+14*a8*b5*a10-812*a9*b4*a10+413*a10*b5*a4+616*a10*b9*a4-133*a10*b8*a8+126*a10^2*b7-14*a10*b9*a9+833*b8*a6*a4-56*a8*b4*a6-1645*a4*b8*a10+910*a4*b4*a10-259*b9*a10*a5+728*b8*a10*a9-385*b8*a6*a9-1547*b8*b10*a6-385*a5*b9*a6+5285*b9*b6*a6-364*b4*a10*a5+56*a6*b5*a9+770*a6*b5*a8+819*a6^2*b7+2681*b4*b10*a6-189*a5*b5*a10-2261*b8*b6*a6+637*b8*b10*a10+1323*b8*a6*a5+6503*b4*b6*a10+7*b8*b6*a10-2317*b4*b10*a10+756*b10^3-4158*b10^2*b6+7182*b6^2*b10-2506*b9*b6*a10"),
		SMQP("364*a6*b5*b6+378*a6*b5*a5-2044*a6*b4*b6+259*a6*b4*a5-245*a6*b5*a4-532*a6*b4*a4-840*b6^3-245*b10*b5*a6-42*b8*a10*a5-112*b10*b5*a10+196*b10*b9*a10+84*a9*b5*a10-658*a6*b9*a4-252*a6*b7*a10+98*a6*b9*a9-259*a6*b9*a8-707*b10*b9*a6+28*b5*b6*a10+287*a6*b4*a9+595*a6*b8*a8-28*b9*a10*a8+224*a8*b4*a10+56*a8*b5*a10-490*a9*b4*a10-112*a10*b5*a4+392*a10*b9*a4-140*a10*b8*a8+84*a10^2*b7-112*a10*b9*a9+1099*b8*a6*a4-322*a8*b4*a6-560*a4*b8*a10+392*a4*b4*a10+28*b9*a10*a5+364*b8*a10*a9-560*b8*a6*a9+581*b8*b10*a6+385*a5*b9*a6+784*b9*b6*a6+154*b4*a10*a5-7*a6*b5*a9-182*a6*b5*a8+483*a6^2*b7+1582*b4*b10*a6-98*a5*b5*a10-35*b8*b6*a6+140*b8*b10*a10-525*b8*a6*a5+1078*b4*b6*a10+84*b8*b6*a10-1106*b4*b10*a10+168*b10^3-924*b10^2*b6+1596*b6^2*b10-588*b9*b6*a10"),
		SMQP("-224*b9^2*a9+364*b9^2*a4-455*b9^2*a8-252*b4^2*b10-56*b7*b4*a10+28*b7*b10*a9-28*b7*a10*a3+525*b7*b4*a6+196*b7*b10*a4+532*b4*a4*a3+392*b8*a8*a3-315*b8*b9*a4-532*b8*a4*a3-658*b8*a9*a7+168*b7*a6*a3+70*b4*a8*a3-350*b4*a8*a7-189*b4*b5*a4+112*a7*b7*a6-490*a7*b5*a9-336*b8*b7*a10-616*a7*b4*a5+910*a7*b4*a9+588*a7*b8*a5-490*a7*b8*a8-28*b8*a9*a3-168*b8*a5*a3+847*b8*b4*a9+350*b8*b4*a5-7*b8*b5*a4-455*b4*b8*a8+1897*b4*b9*a8-98*b4*a4*a7-1764*b4*b5*a8+182*a7*b8*a4+126*b7*a4^2+490*b7*a4*a8+903*b8^2*a5+98*a3*b4*a5+280*a3*b4*a9-406*b7*a5*a9+203*b9^2*a5-686*b7*a5*a8+98*b5*a4*a3+28*b7*b5*a10+98*b5*a4*a7-168*b5^2*a4+84*b5^2*a8+336*b5^2*a9-392*b4*b9*a4-294*b5*b9*a5-42*b5*b8*a5-504*b5*b4*a5+147*b4^2*a9+812*b5*b8*a8+875*b7*b8*a6-357*b4^2*a5-84*b7*a9*a8+336*b7*b6*a5-252*b7*b5*a6+14*b7*a9*a4-378*b7*b6*a4-42*b7*b6*a9-28*b7*b10*a8-448*b5*b9*a8-434*b7*b6*a8+56*b7*a9^2-224*b7*b10*a5-469*b5*b9*a9+28*b7*a8^2+210*a5*b7*a4-154*b8^2*a9-1225*b8^2*a8-385*b9*b5*a4+1022*b5*a9*a3+329*b5*b8*a9-1785*b5*b4*a9-259*b9^2*b10+42*b4^2*b6+1638*b4*b5*b6+1309*b8^2*b6+364*a3*b8*b6-28*a3*b9*b6-224*a3*b8*b10+700*a3*b10*b5+56*a3*b10*b9-1148*a3*b5*b6+14*b9^2*b6+252*b4^2*a4+56*a3*b10*b4-42*b4^2*a8+168*b7*b6^2+798*b10*b9*b8+84*b5^2*b10-252*b5^2*b6-539*b8^2*b10-938*b5*b6*b8+1162*b5*b6*b9-581*b10*b9*b5-168*b7*b6*b10-280*a3*b4*b6-1155*b8*b9*b6+707*b4*b10*b9-987*b4*b10*b5-595*b4*b8*b10-952*b4*b9*b6-112*a7*b4*b6-28*a7*b9*b6+56*a7*b10*b9-224*a7*b10*b5-56*a7*b8*b10+196*a7*b8*b6-112*a7*b10*b4+280*a7*b5*b6-217*b8^2*a4+623*b4*b8*b6+343*b4*b8*a4+637*b5*b8*b10-602*b9*a9*a3+56*b9*a5*a3+546*b9*b8*a9+1512*b9*b8*a8+266*b9*b4*a5-994*b9*b8*a5-749*b9*b4*a9+1610*b5*a8*a3-490*b5*a8*a7-28*b9*a5*a7-854*b9*a8*a3+406*b9*a4*a3+490*b9*a9*a7-350*b9*a4*a7+658*b9*a8*a7+336*b9*b7*a10-847*b9*b7*a6"),
		SMQP("4354*a5*b5*a8+462*b4*a5*a4-2912*b8*b4*a10-3430*b9*b6*a9+3479*b10*b4*a5-672*b7*a10*a5+1876*b10*b9*a9+1848*b8*a6*a3-812*b8*a9*a4-336*b7*a10*a4+1176*b7*a6*a4+504*b7*a10*a9+8778*b8*b6*a5+700*b5*a4*a9-1176*b6*b5*a5+4620*b7*a6*a5-518*b8*a8*a9+609*b8*a8^2-3423*a6*b7*a8-1407*b8*a8*a4-1176*b5*a8^2-700*b8*a10*a3+882*b5*a4^2-2142*b4*a4^2-840*b4*a9^2+308*a8*b9*a9-595*a8*b4*a9-1008*b4*a10*a7+595*a8*b5*a9+672*b9*a10*a7+504*a8*b7*a10-392*a9^2*b5+364*a9*b9*a4-140*a9^2*b9-826*b6*b9*a8+2786*b4*a10*a3+560*b10*b9*a4+2590*b4^2*a10+1344*b4*b6*a4-3192*b7*b6*a10+9702*b7*b6*a6+1512*b7*b10*a10+1792*b8*b6*a9-1134*b8*b10*a9+7784*b4*b6*a9-84*a4*b4*a9-3332*b8*b5*a10+2016*b8*b9*a10-4515*b7*b10*a6+966*b9*a4^2-602*b10*b5*a5-3402*b5^2*a6-2940*b4*a6*a3+2996*b8^2*a6-56*b8*a6*a7+84*b4*b10*a8+903*b10*b5*a9-2100*b9*a6*a7-2394*b4*b9*a10+2604*b9*a6*a3-6748*b9*b8*a6+2086*b5^2*a10-343*b8*b10*a4-196*b4*b8*a6+728*b4*b5*a10+308*b8*a5*a9-2429*b10*b5*a8+4312*b4*a5*a9-469*b10*b5*a4-1386*b9*a10*a3+3752*b4*b9*a6+126*b8*b6*a4-2443*b8*b6*a8-497*b9*b10*a5-1344*b4*b6*a8+2198*b4*a8*a4-3059*b8*b10*a5+672*b4*b5*a6-1974*b4^2*a6+196*b10^2*b4-2380*a5*b9*a9+532*b4*a6*a7-3444*a5^2*b4+672*b7*a6*a9+658*b9*b10*a8+210*b9*b6*a5-2436*b9*b6*a4-1869*b4*a5*a8+2268*b8*a5*a4+672*b8^2*a10-224*b9*a8*a4-2478*b9*a5*a4+1120*b8*a9^2+462*b4*b10*a4+2016*b8*a5^2+504*b8*a4^2-10626*b4*b6*a5+728*b4*a8^2+4116*b9*b6^2+2303*b8*b10^2-924*b4*b6^2-3577*b5*b10^2+721*b10^2*b9-3164*b9*b6*b10-10087*b8*b10*b6-602*b4*b6*b10-14196*b6^2*b5+11634*b6^2*b8+13580*b5*b6*b10-7476*b5*a6*a3+4228*b5*a6*a7+2324*b5*b9*a6+231*b5*a8*a4+434*b5*a10*a3-504*b5*a10*a7-630*b5*b9*a10-2646*b5*a5*a4+5012*b5*b8*a6-1960*a5*b5*a9-371*b9*a8^2+2723*a5*b9*a8-3395*a5*b8*a8+336*b8*a10*a7-2863*b10*b4*a9+672*b9*a5^2+126*b5*b6*a4+2142*b9^2*a6-672*b9^2*a10+7028*b5*b6*a8+1120*b8*b10*a8+504*a5^2*b5-2660*b5*b6*a9"),
		SMQP("4606*a5*b5*a8-504*b4*a5*a4-2912*b8*b4*a10-4018*b9*b6*a9+3570*b10*b4*a5+84*b7*a10*a5+1967*b10*b9*a9+1848*b8*a6*a3-2177*b8*a9*a4-1092*b7*a10*a4+42*b7*a6*a4+504*b7*a10*a9+8085*b8*b6*a5+595*b5*a4*a9-294*b6*b5*a5+3045*b7*a6*a5-581*b8*a8*a9+588*b8*a8^2-3612*a6*b7*a8-1680*b8*a8*a4-1134*b5*a8^2-700*b8*a10*a3+882*b5*a4^2-3276*b4*a4^2-763*b4*a9^2+413*a8*b9*a9+518*a8*b4*a9-1008*b4*a10*a7+784*a8*b5*a9-84*b9*a10*a7+504*a8*b7*a10-497*a9^2*b5+1386*a9*b9*a4-378*a9^2*b9-924*b6*b9*a8+2128*b4*a10*a3+756*b10*b9*a4+2996*b4^2*a10+1680*b4*b6*a4-2436*b7*b6*a10+10206*b7*b6*a6+1512*b7*b10*a10+2737*b8*b6*a9-1239*b8*b10*a9+7070*b4*b6*a9-196*a4*b4*a9-4088*b8*b5*a10-252*b8*b9*a10-4830*b7*b10*a6+1344*b9*a4^2+217*b10*b5*a5-882*b5^2*a6-2982*b4*a6*a3+3248*b8^2*a6+70*b8*a6*a7+224*b4*b10*a8+1113*b10*b5*a9-1862*b9*a6*a7-1736*b4*b9*a10+2310*b9*a6*a3-11858*b9*b8*a6+2086*b5^2*a10-238*b8*b10*a4+1372*b4*b8*a6+6174*b4*b5*a10+2261*b8*a5*a9-2618*b10*b5*a8+2730*b4*a5*a9-364*b10*b5*a4-784*b9*a10*a3-812*b4*b9*a6+336*b8*b6*a4-1960*b8*b6*a8-91*b9*b10*a5-1316*b4*b6*a8+2436*b4*a8*a4-4445*b8*b10*a5-1470*b4*b5*a6-1932*b4^2*a6+616*b10^2*b4-4277*a5*b9*a9+392*b4*a6*a7-6531*a5^2*b4-273*b7*a6*a9+546*b9*b10*a8-294*b9*b6*a5-2352*b9*b6*a4-1190*b4*a5*a8-189*b8*a5*a4+2184*b8^2*a10-364*b9*a8*a4-1050*b9*a5*a4+1330*b8*a9^2+28*b4*b10*a4+4851*b8*a5^2+504*b8*a4^2-8358*b4*b6*a5+924*b4*a8^2+3696*b9*b6^2+2198*b8*b10^2-1344*b4*b6^2-3682*b5*b10^2+826*b10^2*b9-3164*b9*b6*b10-10402*b8*b10*b6+28*b4*b6*b10-13776*b6^2*b5+11424*b6^2*b8+13580*b5*b6*b10-7854*b5*a6*a3+4354*b5*a6*a7+7868*b5*b9*a6+210*b5*a8*a4-1834*b5*a10*a3+252*b5*a10*a7+280*b5*b9*a10-2835*b5*a5*a4-154*b5*b8*a6-2779*a5*b5*a9-364*b9*a8^2+3717*a5*b9*a8-3269*a5*b8*a8+1092*b8*a10*a7-3402*b10*b4*a9-147*b9*a5^2+336*b5*b6*a4+4466*b9^2*a6+84*b9^2*a10+6566*b5*b6*a8+1246*b8*b10*a8-630*a5^2*b5-2660*b5*b6*a9"),
		SMQP("-483*b9^2*a9-2121*b9^2*a4-3381*b9^2*a8+2786*b7*b4*a10-2590*b7*b10*a9+3010*b7*a10*a3-2373*b7*b4*a6+1946*b7*b10*a4+2121*b4*a4*a3-4578*b8*a8*a3+2793*b8*b9*a4+2352*b8*a4*a3-105*b8*a9*a7-210*b7*a6*a3-1659*b4*a8*a3+4200*b4*a8*a7-6678*b5^2*a5+1575*b4*b5*a4-4858*a7*b7*a6+1785*a7*b5*a9+4032*b8*b7*a10-546*a7*b4*a5-504*a7*b7*a10-2730*a7*b4*a9-378*a7*b8*a5+777*a7*b8*a8+840*b8*a9*a3-126*b8*a5*a3-4557*b8*b4*a9-4074*b8*b4*a5-4893*b8*b5*a4+1743*b4*b8*a8-3507*b4*b9*a8-3864*b4*a4*a7+1008*b4*b5*a8-2625*a7*b8*a4-3843*b7*a4^2+1715*b7*a4*a8-3213*b8^2*a5+2226*a3*b4*a5-1533*a3*b4*a9+5719*b7*a5*a9-840*b9^2*a5-8785*b7*a5*a8+3045*b5*a4*a3+14*b7*b5*a10-2877*b5*a4*a7-2835*b5^2*a4+7182*b5^2*a8-3717*b5^2*a9-672*b4*b9*a4-1953*b5*b9*a5+13041*b5*b8*a5-12285*b5*b4*a5-3465*b4^2*a9-3822*b5*b8*a8-8225*b7*b8*a6+2583*b4^2*a5-1806*b7*a9*a8+4998*b7*b6*a5+945*b7*b5*a6+1813*b7*a9*a4-2877*b7*b6*a4+3297*b7*b6*a9+322*b7*b10*a8-2730*b5*b9*a8-4375*b7*b6*a8-140*b7*a9^2-1764*b7*a5^2+3080*b7*b10*a5-1281*b5*b9*a9+2450*b7*a8^2+4431*a5*b7*a4-1176*b8^2*a9-735*b8^2*a8+2121*b9*b5*a4-1995*b5*a9*a3+4998*b5*b8*a9-882*b5*b4*a9-4242*b9^2*b10+630*b4^2*b6+15498*b4*b5*b6+1617*b8^2*b6+5208*a3*b8*b6-7098*a3*b9*b6+1554*a3*b8*b10-3234*a3*b10*b5+1218*a3*b10*b9-3864*a3*b5*b6+9618*b9^2*b6-420*a3*b10*b4+126*b4^2*a8+2394*b7*b6^2+252*b7*b10^2+105*b10*b9*b8+945*b5^2*b10-3276*b5^2*b6+2373*b8^2*b10+7287*b5*b6*b8-1680*b5*b6*b9+9597*b10*b9*b5-252*b7*b6*b10-1428*a3*b4*b6-10290*b8*b9*b6+2121*b4*b10*b9-4725*b4*b10*b5+1155*b4*b8*b10-11298*b4*b9*b6-924*a7*b4*b6-3570*a7*b9*b6-798*a7*b10*b9+42*a7*b10*b5+1050*a7*b8*b10+2184*a7*b8*b6+1092*a7*b10*b4+3318*a7*b5*b6+21*b8^2*a4-1743*b4*b8*b6-399*b4*b8*a4-8022*b5*b8*b10+2373*b9*a9*a3-672*b9*a5*a3-1491*b9*b8*a9+6762*b9*b8*a8+3759*b9*b4*a5+4809*b9*b8*a5+9912*b9*b4*a9-357*b5*a8*a3-861*b5*a8*a7-504*b5*a5*a3-294*b9*a5*a7+5649*b9*a8*a3-4809*b9*a4*a3-1533*b9*a9*a7+5649*b9*a4*a7-651*b9*a8*a7-4536*b9*b7*a10+7714*b9*b7*a6"),
		SMQP("2058*b7*b4*a5+1554*b7*b9*a5-7560*b7*b4*a9-11382*b7*b8*a5+3402*a7*b9^2+378*a7*b8^2+1512*b7^2*a6+2520*b7*b8*a9-31269*b7*b5*a9-252*b7^2*a10-2961*a7*b7*a8-6804*a7*b7*a5-5733*a7*b7*a4-5922*b7*b4*a8+5187*b7*b5*a8+11025*b7*a9*a7+25368*b7*b5*a4-9576*b7*a8*a3+6615*b7*a4*a3+1827*b7*a5*a3-693*b7*a9*a3-13041*b7*b9*a4+945*b7*b8*a8+6174*b7*b9*a8-4410*a7*a3*b8-1764*a3*b8^2+16737*b7*b5*b6+16254*b5^2*b8+9261*b7*b9*a9-5229*b5^3-210*a3*b7*b6-4536*b9*a3*b8-2646*a7*b4*b9+3276*a3^2*b8+10710*a3*b5*a7-7917*b8*b7*b6+3213*b7*b8*b10-4095*b8*b5*a7-9072*a3^2*b5-6426*a7*b7*b6-4284*b4*b8^2-10836*a3*a7*b9-4158*a3*b9^2+7686*a3^2*b9-2772*a7*b8*b9+6237*a7*b5*b9-378*a7*b4*b5+3528*a7*b7*b10-12474*a3^2*b4-9513*a7*b5^2+8316*b4^2*b9+11655*b8*b5*b9-18774*b4*b9^2-1008*a3*b7*b10-37233*b4*b5^2+18270*a7*b4*a3+18144*a3*b5^2+8694*b8^2*b9-1134*b4^2*a3-16443*b9^2*b5-2772*b8^3+1008*b4^2*a7-567*a3*b8*b5-11655*b7*b5*b10+38556*b5*b9*b4+13734*a3*b4*b9-3906*b9^2*b8-11340*b8^2*b5-126*b4*a7*b8-1134*a3*b4*b5+7497*a3*b5*b9+3570*b7*b4*b6-9828*b4^2*b5+1890*a7^2*b8-1890*a7^2*b5+1638*a7^2*b9-4536*a7^2*b4-567*b7*b8*a4+8442*b7*b4*a4-2772*b7*b4*b10-3213*b5^2*b9+5040*b9^3+1050*b9*b7*b6+1323*b9*b7*b10-126*b4*b8*b5-1386*b4*a3*b8+8190*b4*b8*b9+8050*b7*b5*a5"),
		SMQP("1911*b9^2*a9-2247*b9^2*a4-6153*b9^2*a8+5278*b7*b4*a10-5474*b7*b10*a9+5558*b7*a10*a3-9723*b7*b4*a6+5110*b7*b10*a4+7623*b4*a4*a3-630*b8*a8*a3+4893*b8*b9*a4-504*b8*a4*a3+2457*b8*a9*a7-294*b7*a6*a3-6237*b4*a8*a3+5292*b4*a8*a7-4158*b5^2*a5-3171*b4*b5*a4-182*a7*b7*a6+7623*a7*b5*a9-1764*a7*b5*a5+3024*b8*b7*a10-2142*a7*b4*a5-2520*a7*b7*a10-3402*a7*b4*a9-378*a7*b8*a5+567*a7*b8*a8+3024*b8*a9*a3-4914*b8*a5*a3-5733*b8*b4*a9-5544*b8*b4*a5-6783*b8*b5*a4+945*b4*b8*a8+903*b4*b9*a8-7308*b4*a4*a7+4704*b4*b5*a8+189*a7*b8*a4-5985*b7*a4^2+3325*b7*a4*a8-2457*b8^2*a5+10206*a3*b4*a5-8127*a3*b4*a9+10913*b7*a5*a9-3234*b9^2*a5-2303*b7*a5*a8+6615*b5*a4*a3+6538*b7*b5*a10-7623*b5*a4*a7-3759*b5^2*a4+3486*b5^2*a8-4389*b5^2*a9-2940*b4*b9*a4+2163*b5*b9*a5-3339*b5*b8*a5+3171*b5*b4*a5-9135*b4^2*a9-6720*b5*b8*a8-1981*b7*b8*a6+5985*b4^2*a5-4242*b7*a9*a8+1050*b7*b6*a5-26019*b7*b5*a6+2219*b7*a9*a4-10227*b7*b6*a4+12831*b7*b6*a9-1330*b7*b10*a8+3864*b5*b9*a8-2681*b7*b6*a8-1876*b7*a9^2-4032*b7*a5^2+9520*b7*b10*a5+7203*b5*b9*a9+3598*b7*a8^2-3927*a5*b7*a4+4788*b8^2*a9+441*b8^2*a8+5061*b9*b5*a4-16317*b5*a9*a3+3360*b5*b8*a9+168*b5*b4*a9-3864*b9^2*b10+1890*b4^2*b6+23394*b4*b5*b6-1827*b8^2*b6-16002*a3*b9*b6+4662*a3*b8*b10-9702*a3*b10*b5+3654*a3*b10*b9+5292*a3*b5*b6+12390*b9^2*b6-1260*a3*b10*b4+378*b4^2*a8-630*b7*b6^2-1260*b7*b10^2-11445*b10*b9*b8-9975*b5^2*b10+16296*b5^2*b6+5985*b8^2*b10-4557*b5*b6*b8-24780*b5*b6*b9+23667*b10*b9*b5+5292*b7*b6*b10+6048*a3*b4*b6+5628*b8*b9*b6+5523*b4*b10*b9-2247*b4*b10*b5+3969*b4*b8*b10-23898*b4*b9*b6-3780*a7*b4*b6-882*a7*b9*b6-3906*a7*b10*b9+2646*a7*b10*b5+2646*a7*b8*b10+252*a7*b8*b6+3276*a7*b10*b4+3150*a7*b5*b6+1197*b8^2*a4-7497*b4*b8*b6-1701*b4*b8*a4-4746*b5*b8*b10+8883*b9*a9*a3-3276*b9*a5*a3-13125*b9*b8*a9+4578*b9*b8*a8-10101*b9*b4*a5+13671*b9*b8*a5+13566*b9*b4*a9-5859*b5*a8*a3+189*b5*a8*a7+3276*b5*a5*a3+2394*b9*a5*a7+8127*b9*a8*a3-6867*b9*a4*a3-7371*b9*a9*a7+8379*b9*a4*a7-1197*b9*a8*a7-8568*b9*b7*a10+12068*b9*b7*a6"),
		SMQP("1764*b7*b4*a5+126*b7*b9*a5-1092*b7*b4*a9-1862*b7*b8*a5+672*a7*b9^2-588*a7*b8^2+4872*b7^2*a6+4368*b7*b8*a9-29736*b7*b5*a9-2268*b7^2*a10-5586*a7*b7*a8-7182*a7*b7*a5+2478*a7*b7*a4-4536*b7*b4*a8+4578*b7*b5*a8+8610*b7*a9*a7+14133*b7*b5*a4+336*b7*a8*a3+1386*b7*a4*a3-504*b7*a5*a3+3234*b7*a9*a3-9870*b7*b9*a4-2268*b7*b8*a8+3066*b7*b9*a8-3108*a7*a3*b8+588*a3*b8^2+11151*b7*b5*b6+14742*b5^2*b8+7686*b7*b9*a9-3213*b5^3-4116*a3*b7*b6+420*b9*a3*b8+1680*a7*b4*b9+1764*a3^2*b8+10458*a3*b5*a7-2940*b8*b7*b6+1176*b7*b8*b10-1155*b8*b5*a7-5418*a3^2*b5-6972*a7*b7*b6-7812*a3*a7*b9-2520*a3*b9^2+3780*a3^2*b9-2100*a7*b8*b9+9807*a7*b5*b9-4116*a7*b4*b5+2520*a7*b7*b10-7560*a3^2*b4-11781*a7*b5^2+20265*b8*b5*b9-12600*b4*b9^2+1680*a3*b7*b10-37989*b4*b5^2+16464*a7*b4*a3+24948*a3*b5^2+3276*b8^2*b9-16191*b9^2*b5+168*b8^3-10899*a3*b8*b5-7644*b7*b5*b10+36834*b5*b9*b4+2856*a3*b4*b9-3948*b9^2*b8-14742*b8^2*b5+1008*b4*a7*b8+1512*a3*b4*b5-861*a3*b5*b9+3990*b7*b4*b6+672*a7^2*b8-2436*a7^2*b5+1848*a7^2*b9-4536*a7^2*b4+504*b7*b8*a4+2856*b7*b4*a4-2184*b7*b4*b10-6825*b5^2*b9+5544*b9^3+4536*b9*b7*b6-672*b9*b7*b10+1092*b4*b8*b5-840*b4*a3*b8+840*b4*b8*b9+6181*b7*b5*a5"),
		SMQP("1274*a6*b5*b6-378*a6*b5*a5-1274*a6*b4*b6+511*a6*b4*a5-35*a6*b5*a4-1876*a6*b4*a4+203*b10*b5*a6+1470*b8*a10*a5+406*b10*b5*a10-238*b10*b9*a10-70*a9*b5*a10-406*a6*b9*a4-630*a6*b7*a10+770*a6*b9*a9-469*a6*b9*a8+1477*b10*b9*a6-308*b5*b6*a10-49*a6*b4*a9+1141*a6*b8*a8-98*b9*a10*a8+1372*a8*b4*a10+476*a8*b5*a10-854*a9*b4*a10-910*a10*b5*a4+1708*a10*b9*a4-238*a10*b8*a8-476*a10*b9*a9+1813*b8*a6*a4+266*a8*b4*a6-2086*a4*b8*a10-224*a4*b4*a10-490*b9*a10*a5+812*b8*a10*a9-1358*b8*a6*a9-2317*b8*b10*a6+427*a5*b9*a6-1274*b9*b6*a6-1582*b4*a10*a5-119*a6*b5*a9-98*a6*b5*a8-903*a6^2*b7+700*b4*b10*a6-504*a5*b5*a10+203*b8*b6*a6+238*b8*b10*a10-1323*b8*a6*a5+476*b4*b6*a10+1582*b8*b6*a10+56*b4*b10*a10-1204*b9*b6*a10"),
		SMQP("1764*a5*b5*a8-336*b4*a5*a4-3318*b8*b4*a10-5866*b9*b6*a9+2520*b10*b4*a5-672*b7*a10*a5+1617*b10*b9*a9+1008*b8*a6*a3-3493*b8*a9*a4-336*b7*a10*a4+1428*b7*a6*a4+504*b7*a10*a9+9786*b8*b6*a5-133*b5*a4*a9-1176*b6*b5*a5+2415*b7*a6*a5-497*b8*a8*a9+497*b8*a8^2-3801*a6*b7*a8-2695*b8*a8*a4-868*b5*a8^2-798*b8*a10*a3+378*b5*a4^2-2772*b4*a4^2-1183*b4*a9^2+329*a8*b9*a9+1463*a8*b4*a9-1008*b4*a10*a7+679*a8*b5*a9+672*b9*a10*a7+504*a8*b7*a10-469*a9^2*b5+2954*a9*b9*a4-574*a9^2*b9-2009*b6*b9*a8+2884*b4*a10*a3+490*b10*b9*a4+2996*b4^2*a10+1764*b4*b6*a4-3192*b7*b6*a10+8757*b7*b6*a6+1512*b7*b10*a10+4333*b8*b6*a9-889*b8*b10*a9+7385*b4*b6*a9-1876*a4*b4*a9-3234*b8*b5*a10+2016*b8*b9*a10-4452*b7*b10*a6+714*b9*a4^2-315*b10*b5*a5-630*b5^2*a6-2352*b4*a6*a3+924*b8^2*a6-672*b8*a6*a7+56*b4*b10*a8+1127*b10*b5*a9-1120*b9*a6*a7+182*b4*b9*a10+2520*b9*a6*a3-8330*b9*b8*a6+1722*b5^2*a10-350*b8*b10*a4+322*b4*b8*a6+3794*b4*b5*a10+2625*b8*a5*a9-2317*b10*b5*a8+1554*b4*a5*a9-434*b10*b5*a4-1442*b9*a10*a3+2506*b4*b9*a6-483*b8*b6*a4-1064*b8*b6*a8-777*b9*b10*a5-1904*b4*b6*a8+1848*b4*a8*a4-2751*b8*b10*a5-5922*b4*b5*a6-1932*b4^2*a6+112*b10^2*b4-3969*a5*b9*a9+392*b4*a6*a7-5943*a5^2*b4-2037*b7*a6*a9+581*b9*b10*a8-567*b9*b6*a5-1974*b9*b6*a4-2765*b4*a5*a8-147*b8*a5*a4+672*b8^2*a10+840*b9*a8*a4-1344*b9*a5*a4+1526*b8*a9^2+532*b4*b10*a4+4095*b8*a5^2+1134*b8*a4^2-9555*b4*b6*a5+1722*b4*a8^2+3990*b9*b6^2+2324*b8*b10^2-1050*b4*b6^2-3556*b5*b10^2+700*b10^2*b9-3059*b9*b6*b10-10129*b8*b10*b6-308*b4*b6*b10-14070*b6^2*b5+11571*b6^2*b8+13475*b5*b6*b10-4284*b5*a6*a3+2520*b5*a6*a7+6132*b5*b9*a6+245*b5*a8*a4+798*b5*a10*a3-504*b5*a10*a7-574*b5*b9*a10-1407*b5*a5*a4+1134*b5*b8*a6-1575*a5*b5*a9-273*b9*a8^2+2086*a5*b9*a8-798*a5*b8*a8+336*b8*a10*a7-3108*b10*b4*a9-63*b9*a5^2+525*b5*b6*a4+3262*b9^2*a6-672*b9^2*a10+6622*b5*b6*a8+1211*b8*b10*a8-1134*a5^2*b5-2387*b5*b6*a9"),
		SMQP("168*a7^2*b7-112*a7*b7*b4-476*a3*a7*b7+217*a3*b7*b4-308*a3*b7*b8+168*b7*b4*b5+56*a3^2*b7+91*b7^2*a5-98*b7*b4*b9-77*a3*b7*b9+224*b7*b5^2-91*b7*b8*b5+133*a7*b7*b9+14*b7*b8^2+28*b7^2*b10+140*b7^2*a8-91*b7^2*a4-21*b7*b9^2-56*b7*b8*b4-245*b7*b5*b9+147*b7*b8*b9+28*b7*b4^2+903*a3*b7*b5-385*a7*b7*b5+63*a7*b7*b8-49*b7^2*b6"),
		SMQP("168*a5*b5*a8+126*b4*a5*a4+504*b8*b4*a10+462*b9*b6*a9+595*b10*b4*a5+756*b7*a10*a5+196*b10*b9*a9+1008*b8*a6*a3-756*b7*a10*a4-1134*b7*a6*a4-2310*b8*b6*a5+840*b6*b5*a5+1386*b7*a6*a5+504*b8*a8*a4+126*b4*a4^2+140*b4*a9^2+42*a8*b9*a9-210*a8*b4*a9-756*b9*a10*a7-28*a9*b9*a4-28*a9^2*b9+112*b6*b9*a8-490*b4*a10*a3+196*b10*b9*a4-266*b4^2*a10-1176*b4*b6*a4+756*b7*b6*a10+504*b7*b6*a6-420*b8*b6*a9-210*b8*b10*a9-420*b4*b6*a9+140*a4*b4*a9-756*b8*b5*a10-2268*b8*b9*a10-315*b7*b10*a6+378*b9*a4^2+714*b10*b5*a5-2058*b4*a6*a3+1764*b8^2*a6+630*b8*a6*a7+140*b4*b10*a8+105*b10*b5*a9-770*b9*a6*a7-1862*b4*b9*a10+714*b9*a6*a3-2590*b9*b8*a6+105*b8*b10*a4+476*b4*b8*a6+2758*b4*b5*a10+1176*b8*a5*a9-210*b10*b5*a8-280*b4*a5*a9+105*b10*b5*a4+602*b9*a10*a3-2632*b4*b9*a6+714*b8*b6*a4-294*b8*b6*a8+511*b9*b10*a5+826*b4*b6*a8-266*b4*a8*a4-1491*b8*b10*a5+5670*b4*b5*a6-210*b4^2*a6+420*b10^2*b4-952*a5*b9*a9+196*b4*a6*a7+1176*a5^2*b4-91*b9*b10*a8+1554*b9*b6*a5+84*b9*b6*a4+112*b4*a5*a8+168*b8*a5*a4+1512*b8^2*a10-350*b9*a8*a4-630*b9*a5*a4-350*b4*b10*a4-756*b8*a5^2-504*b8*a4^2+546*b4*b6*a5+70*b4*a8^2-420*b9*b6^2-105*b8*b10^2-420*b4*b6^2-105*b5*b10^2+105*b10^2*b9-315*b8*b10*b6+630*b4*b6*b10+420*b6^2*b5-210*b6^2*b8-3402*b5*a6*a3+1134*b5*a6*a7+1512*b5*b9*a6-2268*b5*a10*a3+756*b5*a10*a7+910*b5*b9*a10+168*b5*a5*a4-2142*b5*b8*a6-84*a5*b5*a9-14*b9*a8^2+1036*a5*b9*a8-1092*a5*b8*a8+756*b8*a10*a7-35*b10*b4*a9+420*b9*a5^2+210*b5*b6*a4+812*b9^2*a6+756*b9^2*a10-420*b5*b6*a8+105*b8*b10*a8+210*b5*b6*a9"),
		SMQP("2114*a5*b5*a8-4935*b4*a5*a4-7882*b8*b4*a10-12264*b9*b6*a9+9086*b10*b4*a5-4788*b7*a10*a5+5754*b10*b9*a9+2436*b8*a6*a3-3122*b8*a9*a4-6300*b7*a10*a4+20223*b7*a6*a4+3024*b7*a10*a9+24234*b8*b6*a5-77*b5*a4*a9-2856*b6*b5*a5+7119*b7*a6*a5-126*b8*a8*a9-70*b8*a8^2-4788*a6*b7*a8-3199*b8*a8*a4+98*b5*a8^2-3206*b8*a10*a3-3087*b5*a4^2-2268*b4*a4^2-5264*b4*a9^2+966*a8*b9*a9-336*a8*b4*a9-1008*b4*a10*a7+378*a8*b5*a9+2772*b9*a10*a7+2016*a8*b7*a10-1148*a9^2*b5-630*a9*b9*a4+252*a9^2*b9-5355*b6*b9*a8+3136*b4*a10*a3-3087*b10*b9*a4+18284*b4^2*a10+5754*b4*b6*a4-8820*b7*b6*a10+33201*b7*b6*a6+5040*b7*b10*a10+2520*b8*b6*a9-2338*b8*b10*a9+29211*b4*b6*a9+3997*a4*b4*a9-4606*b8*b5*a10+8316*b8*b9*a10-15498*b7*b10*a6-4032*b9*a4^2+3689*b10*b5*a5+1386*b5^2*a6-5376*b4*a6*a3+1120*b8^2*a6-4060*b8*a6*a7-3416*b4*b10*a8-784*b10*b5*a9-3108*b9*a6*a7-23142*b4*b9*a10+5796*b9*a6*a3+7126*b9*b8*a6+5642*b5^2*a10-2191*b8*b10*a4-5558*b4*b8*a6+19390*b4*b5*a10-3038*b8*a5*a9-7742*b10*b5*a8+11053*b4*a5*a9+7553*b10*b5*a4-3318*b9*a10*a3+44506*b4*b9*a6+2730*b8*b6*a4-3157*b8*b6*a8-567*b9*b10*a5-10150*b4*b6*a8+3458*b4*a8*a4-11977*b8*b10*a5-33810*b4*b5*a6-19068*b4^2*a6-4872*b10^2*b4-2394*a5*b9*a9+4088*b4*a6*a7-10731*a5^2*b4-4788*b7*a6*a9+3780*b9*b10*a8+3129*b9*b6*a5-3780*b9*b6*a4+2534*b4*a5*a8+5796*b8*a5*a4-504*b8^2*a10-609*b9*a8*a4-1743*b9*a5*a4+3556*b8*a9^2-112*b4*b10*a4+6111*b8*a5^2+5229*b8*a4^2-25515*b4*b6*a5+728*b4*a8^2+29610*b9*b6^2+6258*b8*b10^2-8946*b4*b6^2-13902*b5*b10^2+5838*b10^2*b9-23331*b9*b6*b10-24423*b8*b10*b6+9240*b4*b6*b10-51534*b6^2*b5+26523*b6^2*b8+49539*b5*b6*b10-3612*b5*a6*a3+3668*b5*a6*a7+4228*b5*b9*a6+2828*b5*a8*a4+8722*b5*a10*a3-3780*b5*a10*a7-10794*b5*b9*a10+2835*b5*a5*a4+574*b5*b8*a6+7*a5*b5*a9-546*b9*a8^2+21*a5*b9*a8-4669*a5*b8*a8+252*b8*a10*a7-5194*b10*b4*a9+105*b9*a5^2+2415*b5*b6*a4-9870*b9^2*a6+1260*b9^2*a10+15344*b5*b6*a8+2884*b8*b10*a8-126*a5^2*b5-3717*b5*b6*a9"),
		SMQP("1498*b7*b4*a5-126*b7*b9*a5-735*b7*b4*a9-1932*b7*b8*a5+1197*a7*b9^2+651*a7*b8^2+2289*b7^2*a6+2058*b7*b8*a9-14532*b7*b5*a9-1008*b7^2*a10-2541*a7*b7*a8-3906*a7*b7*a5+231*a7*b7*a4-3066*b7*b4*a8+861*b7*b5*a8+4641*b7*a9*a7+9324*b7*b5*a4-84*b7*a8*a3+546*b7*a4*a3-126*b7*a5*a3+903*b7*a9*a3-4746*b7*b9*a4-42*b7*b8*a8+1449*b7*b9*a8-2793*a7*a3*b8+630*a3*b8^2+6552*b7*b5*b6+8358*b5^2*b8+3927*b7*b9*a9-1554*b5^3-1869*a3*b7*b6-714*b9*a3*b8-399*a7*b4*b9+1764*a3^2*b8+4137*a3*b5*a7-2247*b8*b7*b6+840*b7*b8*b10-1155*b8*b5*a7-3717*a3^2*b5-3402*a7*b7*b6+126*b4*b8^2-3045*a3*a7*b9-924*a3*b9^2+1911*a3^2*b9-2604*a7*b8*b9+5271*a7*b5*b9+483*a7*b4*b5+1428*a7*b7*b10-4767*a3^2*b4-6552*a7*b5^2+126*b4^2*b9+10626*b8*b5*b9-5859*b4*b9^2+756*a3*b7*b10-21462*b4*b5^2+9198*a7*b4*a3+13146*a3*b5^2+1869*b8^2*b9-315*b4^2*a3-9051*b9^2*b5+126*b8^3+252*b4^2*a7-5250*a3*b8*b5-4536*b7*b5*b10+20223*b5*b9*b4+1428*a3*b4*b9-1848*b9^2*b8-8463*b8^2*b5+315*b4*a7*b8+525*a3*b4*b5+252*a3*b5*b9+2835*b7*b4*b6+504*b4^2*b5+1302*a7^2*b8-630*a7^2*b5+126*a7^2*b9-2772*a7^2*b4-525*b7*b8*a4+1911*b7*b4*a4-1260*b7*b4*b10-3276*b5^2*b9+2709*b9^3+2079*b9*b7*b6-168*b9*b7*b10+777*b4*b8*b5-315*b4*a3*b8-567*b4*b8*b9+4501*b7*b5*a5"),
		SMQP("-231*b9^2*a9+3171*b9^2*a4-483*b9^2*a8-1512*b4^2*b10-5054*b7*b4*a10+1834*b7*b10*a9-3766*b7*a10*a3+8169*b7*b4*a6-3206*b7*b10*a4-5607*b4*a4*a3-4536*b8*a8*a3-3549*b8*b9*a4+4914*b8*a4*a3+567*b8*a9*a7-4200*b7*a6*a3+9513*b4*a8*a3-3024*b4*a8*a7+3402*b5^2*a5+4725*b4*b5*a4+1582*a7*b7*a6-1575*a7*b5*a9+2016*a7*b5*a5-378*a7*b4*a5+1512*a7*b7*a10+630*a7*b4*a9+1638*a7*b8*a5-63*a7*b8*a8-1260*b8*a9*a3+1260*b8*a5*a3-1323*b8*b4*a9+10962*b8*b4*a5+1197*b8*b5*a4-3087*b4*b8*a8+651*b4*b9*a8+2016*b4*a4*a7-8316*b4*b5*a8-1449*a7*b8*a4+5733*b7*a4^2-4109*b7*a4*a8+693*b8^2*a5-3024*a3*b4*a5+441*a3*b4*a9-8617*b7*a5*a9+1554*b9^2*a5+6391*b7*a5*a8-5985*b5*a4*a3-2282*b7*b5*a10+3843*b5*a4*a7+1701*b5^2*a4-882*b5^2*a8+819*b5^2*a9+1596*b4*b9*a4+3339*b5*b9*a5-4347*b5*b8*a5+3087*b5*b4*a5-2583*b4^2*a9+4158*b5*b8*a8+4697*b7*b8*a6-3087*b4^2*a5+1050*b7*a9*a8+1302*b7*b6*a5+12537*b7*b5*a6+1589*b7*a9*a4+2499*b7*b6*a4-4935*b7*b6*a9+1946*b7*b10*a8-2982*b5*b9*a8+217*b7*b6*a8+644*b7*a9^2+6300*b7*a5^2-6608*b7*b10*a5-1659*b5*b9*a9-1190*b7*a8^2-2793*a5*b7*a4-252*b8^2*a9-1197*b8^2*a8-2289*b9*b5*a4+3087*b5*a9*a3-2520*b5*b8*a9+4788*b5*b4*a9-84*b9^2*b10+4410*b4^2*b6-14490*b4*b5*b6-4977*b8^2*b6-378*a3*b8*b6+3402*a3*b9*b6-504*a3*b8*b10+2268*a3*b10*b5-1260*a3*b10*b9+3276*a3*b5*b6-3738*b9^2*b6+4032*b4^2*a4+2268*a3*b10*b4+3906*b4^2*a8+1386*b7*b6^2+252*b7*b10^2+5565*b10*b9*b8+6489*b5^2*b10-7308*b5^2*b6-1197*b8^2*b10+5859*b5*b6*b8+7644*b5*b6*b9-9177*b10*b9*b5-1764*b7*b6*b10-8568*a3*b4*b6+84*b8*b9*b6+1743*b4*b10*b9+189*b4*b10*b5-4851*b4*b8*b10+5838*b4*b9*b6+3780*a7*b4*b6+630*a7*b9*b6+2394*a7*b10*b9-1134*a7*b10*b5-2142*a7*b8*b10-1260*a7*b10*b4-2394*a7*b5*b6-441*b8^2*a4+9639*b4*b8*b6+63*b4*b8*a4-504*b5*b8*b10-567*b9*a9*a3+2646*b9*a5*a3+2373*b9*b8*a9+3318*b9*b8*a8-5061*b9*b4*a5-4851*b9*b8*a5-1554*b9*b4*a9+6111*b5*a8*a3-1701*b5*a8*a7-5292*b5*a5*a3-2646*b9*a5*a7-1953*b9*a8*a3+1827*b9*a4*a3+819*b9*a9*a7-2583*b9*a4*a7+1197*b9*a8*a7+1512*b9*b7*a10-5824*b9*b7*a6"),
		SMQP("231*b9^2*a9+1869*b9^2*a4-147*b9^2*a8+224*b4^2*b10-3178*b7*b4*a10+182*b7*b10*a9-266*b7*a10*a3+6069*b7*b4*a6-1666*b7*b10*a4-3129*b4*a4*a3-2562*b8*a8*a3-2961*b8*b9*a4+1092*b8*a4*a3-35*b8*a9*a7-546*b7*a6*a3+2751*b4*a8*a3-196*b4*a8*a7+462*b5^2*a5+1589*b4*b5*a4-1918*a7*b7*a6-413*a7*b5*a9+672*a7*b5*a5+1344*b8*b7*a10+1162*a7*b4*a5+504*a7*b7*a10+14*a7*b4*a9-126*a7*b8*a5+259*a7*b8*a8+336*b8*a9*a3+4074*b8*a5*a3-1771*b8*b4*a9+546*b8*b4*a5+357*b8*b5*a4+1085*b4*b8*a8-581*b4*b9*a8+980*b4*a4*a7-196*b4*b5*a8-707*a7*b8*a4+1239*b7*a4^2-2191*b7*a4*a8-651*b8^2*a5-5418*a3*b4*a5+273*a3*b4*a9-707*b7*a5*a9+924*b9^2*a5-2779*b7*a5*a8-4221*b5*a4*a3-2758*b7*b5*a10+2233*b5*a4*a7-441*b5^2*a4+1050*b5^2*a8-63*b5^2*a9-1652*b4*b9*a4+357*b5*b9*a5+1659*b5*b8*a5+21*b5*b4*a5+4081*b4^2*a9-42*b5*b8*a8-2135*b7*b8*a6-847*b4^2*a5+294*b7*a9*a8+714*b7*b6*a5+6867*b7*b5*a6-665*b7*a9*a4+3297*b7*b6*a4-1365*b7*b6*a9+1078*b7*b10*a8-2142*b5*b9*a8-1309*b7*b6*a8+700*b7*a9^2+1428*b7*a5^2-1288*b7*b10*a5-1071*b5*b9*a9+518*b7*a8^2+3213*a5*b7*a4-336*b8^2*a9-273*b8^2*a8+567*b9*b5*a4+651*b5*a9*a3-42*b5*b8*a9+1946*b5*b4*a9-42*b9^2*b10-1918*b4^2*b6-6818*b4*b5*b6-1785*b8^2*b6+3780*a3*b8*b6+1218*a3*b9*b6-1806*a3*b8*b10+1302*a3*b10*b5+42*a3*b10*b9-2856*a3*b5*b6-1386*b9^2*b6-224*b4^2*a4-252*a3*b10*b4-1862*b4^2*a8+966*b7*b6^2+756*b7*b10^2+1995*b10*b9*b8+3171*b5^2*b10-4452*b5^2*b6+819*b8^2*b10+5733*b5*b6*b8+4368*b5*b6*b9-3045*b10*b9*b5-2100*b7*b6*b10-1596*a3*b4*b6-3234*b8*b9*b6-2821*b4*b10*b9+4081*b4*b10*b5-791*b4*b8*b10+6482*b4*b9*b6+364*a7*b4*b6-686*a7*b9*b6+406*a7*b10*b9-658*a7*b10*b5+350*a7*b8*b10-448*a7*b8*b6-308*a7*b10*b4+602*a7*b5*b6+315*b8^2*a4+2443*b4*b8*b6+1043*b4*b8*a4-2478*b5*b8*b10-945*b9*a9*a3-420*b9*a5*a3+735*b9*b8*a9+2310*b9*b8*a8+2765*b9*b4*a5-2709*b9*b8*a5-2576*b9*b4*a9+3969*b5*a8*a3-1295*b5*a8*a7-2688*b5*a5*a3-770*b9*a5*a7-1113*b9*a8*a3+2205*b9*a4*a3+497*b9*a9*a7-1309*b9*a4*a7+791*b9*a8*a7+1176*b9*b7*a10-3374*b9*b7*a6"),
		SMQP("959*b9^2*a9-427*b9^2*a4-1099*b9^2*a8+462*b7*b4*a10-882*b7*b10*a9+798*b7*a10*a3-735*b7*b4*a6+630*b7*b10*a4+707*b4*a4*a3-1526*b8*a8*a3+1071*b8*b9*a4+784*b8*a4*a3+805*b8*a9*a7-630*b7*a6*a3-721*b4*a8*a3+1400*b4*a8*a7-546*b5^2*a5-595*b4*b5*a4-1134*a7*b7*a6+1435*a7*b5*a9+840*a7*b5*a5+1344*b8*b7*a10+154*a7*b4*a5-168*a7*b7*a10-910*a7*b4*a9+714*a7*b8*a5+427*a7*b8*a8+280*b8*a9*a3-42*b8*a5*a3-1519*b8*b4*a9-686*b8*b4*a5-1351*b8*b5*a4+581*b4*b8*a8-1225*b4*b9*a8-1288*b4*a4*a7+1568*b4*b5*a8-875*a7*b8*a4-777*b7*a4^2-63*b7*a4*a8+609*b8^2*a5-434*a3*b4*a5-1351*a3*b4*a9+1701*b7*a5*a9+1540*b9^2*a5-1491*b7*a5*a8+1015*b5*a4*a3+210*b7*b5*a10-959*b5*a4*a7-665*b5^2*a4+1834*b5^2*a8-959*b5^2*a9-784*b4*b9*a4+2989*b5*b9*a5-693*b5*b8*a5+1897*b5*b4*a5-1155*b4^2*a9-1162*b5*b8*a8-2331*b7*b8*a6+861*b4^2*a5-546*b7*a9*a8+3066*b7*b6*a5-21*b7*b5*a6-273*b7*a9*a4-1407*b7*b6*a4+1995*b7*b6*a9+126*b7*b10*a8-1302*b5*b9*a8-861*b7*b6*a8-84*b7*a9^2+924*b7*a5^2+840*b7*b10*a5+1113*b5*b9*a9+798*b7*a8^2-651*a5*b7*a4+1288*b8^2*a9+91*b8^2*a8+1407*b9*b5*a4-3185*b5*a9*a3+266*b5*b8*a9+2786*b5*b4*a9-1274*b9^2*b10+210*b4^2*b6+4606*b4*b5*b6+539*b8^2*b6+1736*a3*b8*b6-2366*a3*b9*b6+518*a3*b8*b10-1078*a3*b10*b5+406*a3*b10*b9-1288*a3*b5*b6+2926*b9^2*b6-140*a3*b10*b4+42*b4^2*a8+798*b7*b6^2+84*b7*b10^2-105*b10*b9*b8+35*b5^2*b10-532*b5^2*b6+791*b8^2*b10+2149*b5*b6*b8-840*b5*b6*b9+3339*b10*b9*b5-84*b7*b6*b10-476*a3*b4*b6-3570*b8*b9*b6+1267*b4*b10*b9-455*b4*b10*b5+385*b4*b8*b10-4046*b4*b9*b6-308*a7*b4*b6-1190*a7*b9*b6-266*a7*b10*b9+14*a7*b10*b5+350*a7*b8*b10+728*a7*b8*b6+364*a7*b10*b4+1106*a7*b5*b6+7*b8^2*a4-581*b4*b8*b6-133*b4*b8*a4-2954*b5*b8*b10+1631*b9*a9*a3+616*b9*a5*a3-3297*b9*b8*a9+1890*b9*b8*a8-791*b9*b4*a5-3017*b9*b8*a5+1484*b9*b4*a9-623*b5*a8*a3-119*b5*a8*a7-2688*b5*a5*a3-938*b9*a5*a7+2051*b9*a8*a3-1603*b9*a4*a3-1351*b9*a9*a7+1883*b9*a4*a7-385*b9*a8*a7-1512*b9*b7*a10+1722*b9*b7*a6"),
		SMQP("252*a7^2*b7-168*a7*b7*b4-924*a3*a7*b7+525*a3*b7*b4-756*a3*b7*b8-420*b7*b4*b5-126*a3^2*b7-245*b7^2*a5+168*b7*b4*b9-21*a3*b7*b9+294*b7*b5^2+315*b7*b8*b5+357*a7*b7*b9-294*b7*b8^2+42*b7^2*b10+210*b7^2*a8+21*b7^2*a4-189*b7*b9^2-84*b7*b8*b4-777*b7*b5*b9+693*b7*b8*b9+42*b7*b4^2+2079*a3*b7*b5-525*a7*b7*b5-63*a7*b7*b8-231*b7^2*b6"),
		SMQP("945*b9^2*a9+2331*b9^2*a4+1827*b9^2*a8+1344*b4^2*b10-5390*b7*b4*a10-14*b7*b10*a9+98*b7*a10*a3+10122*b7*b4*a6-1274*b7*b10*a4-3234*b4*a4*a3-546*b8*a8*a3-6237*b8*b9*a4+714*b8*a4*a3-483*b8*a9*a7+4494*b7*a6*a3-2415*b4*a8*a3+231*b4*a8*a7-882*b5^2*a5+2835*b4*b5*a4-3290*a7*b7*a6-483*a7*b5*a9-252*a7*b5*a5+1008*b8*b7*a10-42*a7*b4*a5+504*a7*b7*a10+231*a7*b4*a9+378*a7*b8*a5+651*a7*b8*a8-546*b8*a9*a3+1512*b8*a5*a3+903*b8*b4*a9-1827*b8*b4*a5+1785*b8*b5*a4+525*b4*b8*a8-2226*b4*b9*a8+1113*b4*a4*a7+6741*b4*b5*a8+1029*a7*b8*a4+945*b7*a4^2-2681*b7*a4*a8-3339*b8^2*a5-2499*a3*b4*a5+3822*a3*b4*a9+1043*b7*a5*a9+504*b9^2*a5-5285*b7*a5*a8-4263*b5*a4*a3-3122*b7*b5*a10+1911*b5*a4*a7-2121*b5^2*a4+2226*b5^2*a8+21*b5^2*a9-1911*b4*b9*a4-2247*b5*b9*a5+6489*b5*b8*a5-294*b5*b4*a5+7224*b4^2*a9-2310*b5*b8*a8+749*b7*b8*a6-168*b4^2*a5+546*b7*a9*a8+798*b7*b6*a5+3843*b7*b5*a6-3367*b7*a9*a4+3003*b7*b6*a4-147*b7*b6*a9+770*b7*b10*a8-2856*b5*b9*a8-875*b7*b6*a8+980*b7*a9^2+364*b7*b10*a5-1029*b5*b9*a9+994*b7*a8^2+6279*a5*b7*a4-1260*b8^2*a9+2205*b8^2*a8+2877*b9*b5*a4+777*b5*a9*a3+1470*b5*b8*a9-2268*b5*b4*a9+630*b9^2*b10-4074*b4^2*b6-5418*b4*b5*b6-63*b8^2*b6+4704*a3*b8*b6+3738*a3*b9*b6-2478*a3*b8*b10+3570*a3*b10*b5-1050*a3*b10*b9-6636*a3*b5*b6-2898*b9^2*b6-1344*b4^2*a4-1176*a3*b10*b4-3234*b4^2*a8+2142*b7*b6^2+1260*b7*b10^2+1953*b10*b9*b8+2247*b5^2*b10-5376*b5^2*b6+693*b8^2*b10+5145*b5*b6*b8+7392*b5*b6*b9-3129*b10*b9*b5-3528*b7*b6*b10-1176*a3*b4*b6-4914*b8*b9*b6-3759*b4*b10*b9+4851*b4*b10*b5-525*b4*b8*b10+7518*b4*b9*b6+84*a7*b4*b6-1302*a7*b9*b6+462*a7*b10*b9-1218*a7*b10*b5+798*a7*b8*b10-1092*a7*b8*b6+84*a7*b10*b4+2058*a7*b5*b6+1701*b8^2*a4+1617*b4*b8*b6+1281*b4*b8*a4-2604*b5*b8*b10-2415*b9*a9*a3+84*b9*a5*a3+567*b9*b8*a9-756*b9*b8*a8+7644*b9*b4*a5-2205*b9*b8*a5-6069*b9*b4*a9+903*b5*a8*a3-105*b5*a8*a7+756*b5*a5*a3-42*b9*a5*a7-1155*b9*a8*a3+2499*b9*a4*a3+1239*b9*a9*a7-2667*b9*a4*a7+105*b9*a8*a7+2520*b9*b7*a10-7882*b9*b7*a6"),
		SMQP("252*a7^2*b7-126*a7*b7*b4-84*b7^2*a9-504*a3*a7*b7+441*a3*b7*b4-546*a3*b7*b8-210*b7*b4*b5+378*a3^2*b7+329*b7^2*a5-84*b7*b4*b9-273*a3*b7*b9+336*b7*b5^2+21*b7*b8*b5+315*a7*b7*b9-126*b7*b8^2+84*b7^2*b10+252*b7^2*a8-105*b7^2*a4-147*b7*b9^2+42*b7*b8*b4-357*b7*b5*b9+441*b7*b8*b9+1449*a3*b7*b5-945*a7*b7*b5+21*a7*b7*b8-189*b7^2*b6"),
		SMQP("2037*b7*b4*a5-315*b7*b9*a5-273*b7*b4*a9-1911*b7*b8*a5+882*a7*b9^2-378*a7*b8^2+2457*b7^2*a6+1974*b7*b8*a9-14742*b7*b5*a9-1008*b7^2*a10-2205*a7*b7*a8-4536*a7*b7*a5-63*a7*b7*a4-4116*b7*b4*a8+609*b7*b5*a8+4725*b7*a9*a7+10332*b7*b5*a4-1008*b7*a8*a3+1827*b7*a4*a3-63*b7*a5*a3+819*b7*a9*a3-5313*b7*b9*a4+336*b7*b8*a8+1995*b7*b9*a8-3150*a7*a3*b8+252*a3*b8^2+7014*b7*b5*b6+9072*b5^2*b8+3759*b7*b9*a9-1890*b5^3-1974*a3*b7*b6-504*b9*a3*b8-126*a7*b4*b9+2016*a3^2*b8+5922*a3*b5*a7-2772*b8*b7*b6+1302*b7*b8*b10-1008*b8*b5*a7-3402*a3^2*b5-3780*a7*b7*b6+252*b4*b8^2-3654*a3*a7*b9-1260*a3*b9^2+1890*a3^2*b9-1008*a7*b8*b9+4536*a7*b5*b9-882*a7*b4*b5+1512*a7*b7*b10-4914*a3^2*b4-6678*a7*b5^2+252*b4^2*b9+11466*b8*b5*b9-5670*b4*b9^2+504*a3*b7*b10-22554*b4*b5^2+9324*a7*b4*a3+13860*a3*b5^2+1386*b8^2*b9-630*b4^2*a3-9576*b9^2*b5+252*b8^3+504*b4^2*a7-4914*a3*b8*b5-4746*b7*b5*b10+20034*b5*b9*b4+3024*a3*b4*b9-1008*b9^2*b8-9450*b8^2*b5+630*b4*a7*b8-1638*a3*b4*b5+630*a3*b5*b9+3150*b7*b4*b6+1008*b4^2*b5+756*a7^2*b8-1260*a7^2*b5+756*a7^2*b9-2520*a7^2*b4-924*b7*b8*a4+2562*b7*b4*a4-1680*b7*b4*b10-3150*b5^2*b9+2394*b9^3+2142*b9*b7*b6-294*b9*b7*b10+1386*b4*b8*b5-126*b4*a3*b8-1134*b4*b8*b9+4123*b7*b5*a5"),
		SMQP("-399*b4*a5*a4-42*b9*b6*a9-196*b10*b4*a5-42*b8*a9*a4+189*b7*a6*a4+315*b8*b6*a5+21*b5*a4*a9-126*b6*b5*a5+21*b8*a8*a4+21*b5*a4^2+546*b4*a4^2-140*b4*a9^2+210*a8*b4*a9+42*a9*b9*a4+21*b6*b9*a8-266*b4*a10*a3+21*b10*b9*a4-1246*b4^2*a10-378*b4*b6*a4-189*b7*b6*a6+42*b8*b6*a9-273*b4*b6*a9+343*a4*b4*a9-210*b9*a4^2-84*b4*a6*a3+70*b4*b10*a8-21*b8*b10*a4+1540*b4*b8*a6+266*b4*b5*a10-896*b4*a5*a9-21*b10*b5*a4-602*b4*b9*a6-546*b8*b6*a4-21*b8*b6*a8+812*b4*b6*a8-952*b4*a8*a4-630*b4*b5*a6+714*b4^2*a6+308*b4*a6*a7+504*a5^2*b4-147*b9*b6*a5+420*b9*b6*a4+1232*b4*a5*a8-315*b8*a5*a4-21*b9*a8*a4+147*b9*a5*a4+14*b4*b10*a4+273*b8*a4^2-315*b4*b6*a5-70*b4*a8^2-210*b9*b6^2+42*b4*b6^2-21*b9*b6*b10+21*b8*b10*b6-84*b4*b6*b10-42*b6^2*b5+273*b6^2*b8+21*b5*b6*b10-42*b5*a8*a4+126*b5*a5*a4-70*b10*b4*a9+21*b5*b6*a4+42*b5*b6*a8-21*b5*b6*a9"),
		SMQP("-84*b9^2*a4-56*b4^2*b10+42*b7*b4*a6+140*b4*a4*a3+7*b8*a8*a3+252*b8*b9*a4+91*b8*a4*a3+63*b7*a6*a3-70*b4*a8*a3-42*b5^2*a5-406*b4*b5*a4-14*b8*a9*a3-105*b8*a5*a3+28*b8*b4*a9+210*b8*b4*a5-7*b8*b5*a4-14*b4*b8*a8+14*b4*b9*a8+98*b4*b5*a8-84*a7*b8*a4+84*b7*a4^2+133*a3*b4*a5-7*a3*b4*a9+259*b5*a4*a3-84*b5*a4*a7-7*b5^2*a4+14*b5^2*a8-7*b5^2*a9+140*b4*b9*a4-49*b5*b9*a5+105*b5*b8*a5-217*b5*b4*a5+14*b4^2*a9-7*b5*b8*a8-98*b4^2*a5+84*b7*b6*a5-63*b7*b5*a6-168*b7*b6*a4+7*b5*b9*a8-14*b5*b9*a9-84*a5*b7*a4-14*b9*b5*a4+7*b5*a9*a3+14*b5*b8*a9-7*b5*b4*a9+28*b4^2*b6+322*b4*b5*b6+168*b8^2*b6-91*a3*b8*b6+154*a3*b9*b6-7*a3*b8*b10-7*a3*b10*b5+7*a3*b10*b9-238*a3*b5*b6+84*b9^2*b6+56*b4^2*a4+28*a3*b10*b4-28*b4^2*a8+84*b7*b6^2+7*b5^2*b10-14*b5^2*b6+7*b5*b6*b8+14*b5*b6*b9-7*b10*b9*b5-98*a3*b4*b6-252*b8*b9*b6-14*b4*b10*b9-14*b4*b10*b5+14*b4*b8*b10-140*b4*b9*b6-84*a7*b9*b6+84*a7*b8*b6+84*a7*b5*b6-168*b8^2*a4+14*b4*b8*b6-14*b4*b8*a4+7*b5*b8*b10+14*b9*a9*a3+49*b9*a5*a3-98*b9*b4*a5-28*b9*b4*a9-14*b5*a8*a3+42*b5*a5*a3-7*b9*a8*a3-154*b9*a4*a3+84*b9*a4*a7"),
		SMQP("56*a5*b5*a8-924*b4*a5*a4+560*b8*b4*a10+2898*b9*b6*a9-693*b10*b4*a5+42*b7*a10*a5-1582*b10*b9*a9+924*b8*a6*a3+1624*b8*a9*a4+966*b7*a10*a4-2856*b7*a6*a4-1176*b7*a10*a9-4431*b8*b6*a5-308*b5*a4*a9+42*b6*b5*a5-1932*b7*a6*a5-364*b8*a8*a9+448*b8*a8^2-84*a6*b7*a8+3136*b8*a8*a4-686*b5*a8^2+1540*b8*a10*a3+1386*b5*a4^2+2268*b4*a4^2+1512*b4*a9^2+280*a8*b9*a9-560*a8*b4*a9+1008*b4*a10*a7+182*a8*b5*a9+126*b9*a10*a7+168*a8*b7*a10+532*a9^2*b5+112*a9*b9*a4+112*a9^2*b9+1771*b6*b9*a8-1386*b4*a10*a3-742*b10*b9*a4-2520*b4^2*a10-1764*b4*b6*a4+2562*b7*b6*a10-8631*b7*b6*a6-1512*b7*b10*a10+1050*b8*b6*a9-406*b8*b10*a9-9597*b4*b6*a9-1428*a4*b4*a9+1274*b8*b5*a10+378*b8*b9*a10+4221*b7*b10*a6+840*b9*a4^2-1498*b10*b5*a5-1386*b5^2*a6-252*b4*a6*a3+3808*b8^2*a6+560*b8*a6*a7+854*b4*b10*a8+203*b10*b5*a9+896*b9*a6*a7+3164*b4*b9*a10-1512*b9*a6*a3+224*b9*b8*a6-1078*b5^2*a10+917*b8*b10*a4+2604*b4*b8*a6-10486*b4*b5*a10+3304*b8*a5*a9+896*b10*b5*a8-5628*b4*a5*a9+203*b10*b5*a4-182*b9*a10*a3-8316*b4*b9*a6+21*b8*b6*a4-3507*b8*b6*a8+707*b9*b10*a5+3962*b4*b6*a8-3080*b4*a8*a4+3101*b8*b10*a5+17346*b4*b5*a6+756*b4^2*a6+252*b10^2*b4-728*a5*b9*a9-1008*b4*a6*a7+5460*a5^2*b4+1344*b7*a6*a9-525*b9*b10*a8-3549*b9*b6*a5+2394*b9*b6*a4+2576*b4*a5*a8+1008*b8*a5*a4-924*b8^2*a10-672*b9*a8*a4+168*b9*a5*a4-1568*b8*a9^2-504*b4*b10*a4-4284*b8*a5^2-2772*b8*a4^2+3003*b4*b6*a5+28*b4*a8^2-9450*b9*b6^2-441*b8*b10^2+2646*b4*b6^2+4095*b5*b10^2-2583*b10^2*b9+9891*b9*b6*b10+2394*b8*b10*b6-630*b4*b6*b10+12474*b6^2*b5-4221*b6^2*b8-14427*b5*b6*b10-1344*b5*a6*a3+140*b5*a6*a7-4970*b5*b9*a6-1442*b5*a8*a4-56*b5*a10*a3+42*b5*a10*a7+2954*b5*b9*a10-882*b5*a5*a4-1064*b5*b8*a6-1484*a5*b5*a9+84*b9*a8^2+2128*a5*b9*a8-2296*a5*b8*a8-1134*b8*a10*a7+3843*b10*b4*a9+1260*b9*a5^2-1995*b5*b6*a4+784*b9^2*a6-1470*b9^2*a10-1512*b5*b6*a8-7*b8*b10*a8+1008*a5^2*b5+1239*b5*b6*a9"),
		SMQP("665*b7*b4*a5+105*b7*b9*a5-455*b7*b4*a9-1827*b7*b8*a5+504*a7*b9^2-840*a7*b8^2+2583*b7^2*a6+1106*b7*b8*a9-13384*b7*b5*a9-1008*b7^2*a10-2793*a7*b7*a8-1848*a7*b7*a5+1281*a7*b7*a4-1904*b7*b4*a8+4865*b7*b5*a8+3801*b7*a9*a7+5264*b7*b5*a4-1512*b7*a8*a3+945*b7*a4*a3+567*b7*a5*a3+1575*b7*a9*a3-5495*b7*b9*a4-1834*b7*b8*a8+2065*b7*b9*a8-168*a7*a3*b8+672*a3*b8^2+4396*b7*b5*b6+6216*b5^2*b8+3703*b7*b9*a9-1764*b5^3-1680*a3*b7*b6+168*b9*a3*b8+1344*a7*b4*b9+5208*a3*b5*a7-1568*b8*b7*b6+280*b7*b8*b10-420*b8*b5*a7-2520*a3^2*b5-2856*a7*b7*b6-4536*a3*a7*b9-1680*a3*b9^2+2016*a3^2*b9-672*a7*b8*b9+3780*a7*b5*b9-3192*a7*b4*b5+1008*a7*b7*b10-3024*a3^2*b4-4956*a7*b5^2+7308*b8*b5*b9-6720*b4*b9^2+672*a3*b7*b10-13860*b4*b5^2+6720*a7*b4*a3+10080*a3*b5^2+2184*b8^2*b9-6972*b9^2*b5-1008*b8^3-4452*a3*b8*b5-3080*b7*b5*b10+14616*b5*b9*b4+2016*a3*b4*b9-1512*b9^2*b8-4704*b8^2*b5+672*b4*a7*b8+1008*a3*b4*b5+924*a3*b5*b9+1064*b7*b4*b6-168*a7^2*b8-1176*a7^2*b5+1176*a7^2*b9-1680*a7^2*b4+2240*b7*b8*a4+784*b7*b4*a4-448*b7*b4*b10-2436*b5^2*b9+2352*b9^3+2240*b9*b7*b6-280*b9*b7*b10-2352*b4*b8*b5-1176*b4*a3*b8+3024*b4*b8*b9+413*b7*b5*a5"),
		SMQP("3885*b7*b4*a5+609*b7*b9*a5-7119*b7*b4*a9-4683*b7*b8*a5+1260*a7*b9^2+252*a7*b8^2-441*b7^2*a6+1890*b7*b8*a9-15456*b7*b5*a9+504*b7^2*a10-189*a7*b7*a8-4788*a7*b7*a5-4347*a7*b7*a4-2016*b7*b4*a8+21*b7*b5*a8+5229*b7*a9*a7+10878*b7*b5*a4-5544*b7*a8*a3+6237*b7*a4*a3+3591*b7*a5*a3+315*b7*a9*a3-5607*b7*b9*a4+1638*b7*b8*a8+3213*b7*b9*a8-2016*a7*a3*b8+3528*a3*b8^2+7854*b7*b5*b6+8316*b5^2*b8+3843*b7*b9*a9-2394*b5^3+3696*a3*b7*b6-10332*b9*a3*b8-1260*a7*b4*b9+3528*a3^2*b8+9324*a3*b5*a7-3444*b8*b7*b6+2016*b7*b8*b10-882*b8*b5*a7-11340*a3^2*b5-3528*a7*b7*b6-3528*b4*b8^2-7056*a3*a7*b9+1260*a3*b9^2+5040*a3^2*b9-504*a7*b8*b9+1890*a7*b5*b9-2268*a7*b4*b5+2016*a7*b7*b10-8568*a3^2*b4-4410*a7*b5^2+7560*b4^2*b9+6930*b8*b5*b9-10332*b4*b9^2-2016*a3*b7*b10-20034*b4*b5^2+9072*a7*b4*a3+7560*a3*b5^2+1764*b8^2*b9+252*b4^2*a3-8442*b9^2*b5+504*b8^3-3402*a3*b8*b5-5544*b7*b5*b10+23184*b5*b9*b4+504*a3*b4*b9-8064*b8^2*b5+252*b4*a7*b8+10836*a3*b4*b5+6930*a3*b5*b9+168*b7*b4*b6-11592*b4^2*b5+1008*a7^2*b8-1512*a7^2*b5+1008*a7^2*b9-2016*a7^2*b4-2772*b7*b8*a4+5796*b7*b4*a4-1008*b7*b4*b10-882*b5^2*b9+1764*b9^3-588*b9*b7*b6+1008*b9*b7*b10+3780*b4*b8*b5+756*b4*a3*b8+2772*b4*b8*b9+5131*b7*b5*a5"),
		SMQP("-448*a5*b5*a8-462*b4*a5*a4+707*b8*b4*a10+1701*b9*b6*a9-490*b10*b4*a5-903*b10*b9*a9+42*b8*a6*a3+238*b8*a9*a4+252*b7*a10*a4-819*b7*a6*a4-1974*b8*b6*a5-434*b5*a4*a9+42*b6*b5*a5-567*b7*a6*a5-105*b8*a8*a9+119*b8*a8^2+1008*a6*b7*a8+707*b8*a8*a4-217*b5*a8^2+805*b8*a10*a3+819*b5*a4^2+1386*b4*a4^2+196*b4*a9^2-63*a8*b9*a9+210*a8*b4*a9+504*b4*a10*a7+147*a8*b5*a9-504*b9*a10*a7-252*a8*b7*a10+574*a9^2*b5+210*a9*b9*a4-294*a9^2*b9+735*b6*b9*a8-1190*b4*a10*a3-651*b10*b9*a4-1330*b4^2*a10-420*b4*b6*a4+2016*b7*b6*a10-5355*b7*b6*a6-1008*b7*b10*a10+651*b8*b6*a9-133*b8*b10*a9-3948*b4*b6*a9-308*a4*b4*a9+455*b8*b5*a10-756*b8*b9*a10+2898*b7*b10*a6-315*b9*a4^2+56*b10*b5*a5+63*b5^2*a6-147*b4*a6*a3-224*b8^2*a6+245*b8*a6*a7-602*b4*b10*a8-91*b10*b5*a9+987*b9*a6*a7+861*b4*b9*a10-1323*b9*a6*a3+2254*b9*b8*a6-889*b5^2*a10+497*b8*b10*a4+511*b4*b8*a6-3479*b4*b5*a10+784*b8*a5*a9+973*b10*b5*a8-3122*b4*a5*a9+161*b10*b5*a4+651*b9*a10*a3-2219*b4*b9*a6-546*b8*b6*a4-1351*b8*b6*a8-42*b9*b10*a5+2282*b4*b6*a8-1330*b4*a8*a4+2030*b8*b10*a5+5817*b4*b5*a6+462*b4^2*a6+420*a5*b9*a9-280*b4*a6*a7+3318*a5^2*b4-504*b7*a6*a9-483*b9*b10*a8-1848*b9*b6*a5+1260*b9*b6*a4+980*b4*a5*a8+693*b8*a5*a4+357*b9*a8*a4+21*b9*a5*a4-266*b8*a9^2-154*b4*b10*a4-1638*b8*a5^2-693*b8*a4^2+630*b4*b6*a5-154*b4*a8^2-4284*b9*b6^2-378*b8*b10^2+1764*b4*b6^2+1890*b5*b10^2-1134*b10^2*b9+4410*b9*b6*b10+2016*b8*b10*b6-1260*b4*b6*b10+5796*b6^2*b5-2142*b6^2*b8-6678*b5*b6*b10+1365*b5*a6*a3-1351*b5*a6*a7-3227*b5*b9*a6-889*b5*a8*a4-1379*b5*a10*a3+756*b5*a10*a7+1113*b5*b9*a10-567*b5*a5*a4+1750*b5*b8*a6-1022*a5*b5*a9+105*b9*a8^2+378*a5*b9*a8+14*a5*b8*a8-252*b8*a10*a7+1988*b10*b4*a9+42*b9*a5^2-1050*b5*b6*a4-798*b9^2*a6-805*b5*b6*a8+259*b8*b10*a8+1008*a5^2*b5-231*b5*b6*a9"),
		SMQP("819*b9^2*a9+105*b9^2*a4-987*b9^2*a8+336*b4^2*b10+770*b7*b4*a10-2254*b7*b10*a9+2338*b7*a10*a3+378*b7*b4*a6+3290*b7*b10*a4+3850*b4*a4*a3+4046*b8*a8*a3-112*b8*b9*a4-2821*b8*a4*a3-2226*b8*a9*a7+3318*b7*a6*a3-6839*b4*a8*a3+1491*b4*a8*a7-2142*b5^2*a5+1295*b4*b5*a4-826*a7*b7*a6+861*a7*b5*a9-2520*a7*b5*a5-1008*b8*b7*a10-2058*a7*b4*a5-1176*a7*b7*a10+1995*a7*b4*a9+294*a7*b8*a5-294*a7*b8*a8+581*b8*a9*a3-3969*b8*a5*a3+4676*b8*b4*a9-2702*b8*b4*a5-973*b8*b5*a4-511*b4*b8*a8+2646*b4*b9*a8-3171*b4*a4*a7+3605*b4*b5*a8+2310*a7*b8*a4-2583*b7*a4^2+2639*b7*a4*a8-1386*b8^2*a5+5285*a3*b4*a5-1358*a3*b4*a9+4795*b7*a5*a9-1050*b9^2*a5-5677*b7*a5*a8+2681*b5*a4*a3+2702*b7*b5*a10-2457*b5*a4*a7-3479*b5^2*a4+2422*b5^2*a8-497*b5^2*a9-2667*b4*b9*a4-2975*b5*b9*a5+5712*b5*b8*a5-5012*b5*b4*a5-1176*b4^2*a9-2443*b5*b8*a8+4480*b7*b8*a6+1512*b4^2*a5-1470*b7*a9*a8-882*b7*b6*a5-12915*b7*b5*a6-1295*b7*a9*a4-4809*b7*b6*a4+5397*b7*b6*a9-1022*b7*b10*a8-364*b5*b9*a8-1267*b7*b6*a8-140*b7*a9^2-4284*b7*a5^2+4424*b7*b10*a5+1043*b5*b9*a9+1442*b7*a8^2+3339*a5*b7*a4-98*b8^2*a9-665*b8^2*a8+2093*b9*b5*a4-1183*b5*a9*a3+3374*b5*b8*a9-7714*b5*b4*a9-1386*b9^2*b10-882*b4^2*b6+12418*b4*b5*b6+4781*b8^2*b6+280*a3*b8*b6-5194*a3*b9*b6+490*a3*b8*b10+1078*a3*b10*b5+602*a3*b10*b9-896*a3*b5*b6+2226*b9^2*b6-336*b4^2*a4+140*a3*b10*b4-1050*b4^2*a8+42*b7*b6^2-84*b7*b10^2-875*b10*b9*b8-3619*b5^2*b10+4340*b5^2*b6-1183*b8^2*b10-5621*b5*b6*b8-1904*b5*b6*b9+4921*b10*b9*b5+420*b7*b6*b10+2660*a3*b4*b6-1442*b8*b9*b6+4557*b4*b10*b9-6041*b4*b10*b5+1183*b4*b8*b10-10122*b4*b9*b6-1932*a7*b4*b6+126*a7*b9*b6-798*a7*b10*b9-630*a7*b10*b5+714*a7*b8*b10-504*a7*b8*b6+1428*a7*b10*b4+1974*a7*b5*b6+1057*b8^2*a4-2891*b4*b8*b6-427*b4*b8*a4+1666*b5*b8*b10-245*b9*a9*a3-1918*b9*a5*a3-2926*b9*b8*a9+1337*b9*b8*a8+3822*b9*b4*a5+2086*b9*b8*a5-903*b9*b4*a9-2401*b5*a8*a3+231*b5*a8*a7+6048*b5*a5*a3+2562*b9*a5*a7+637*b9*a8*a3-119*b9*a4*a3+189*b9*a9*a7+651*b9*a4*a7+147*b9*a8*a7-168*b9*b7*a10-1988*b9*b7*a6"),
		SMQP("2247*b7*b4*a5+357*b7*b9*a5-2268*b7*b4*a9-4242*b7*b8*a5+1323*a7*b9^2-630*a7*b8^2+4536*b7^2*a6+4032*b7*b8*a9-30072*b7*b5*a9-1764*b7^2*a10-4914*a7*b7*a8-7056*a7*b7*a5+504*a7*b7*a4-5418*b7*b4*a8+5250*b7*b5*a8+9198*b7*a9*a7+17367*b7*b5*a4-2520*b7*a8*a3+2520*b7*a4*a3-126*b7*a5*a3+2142*b7*a9*a3-11025*b7*b9*a4-1512*b7*b8*a8+3906*b7*b9*a8-3654*a7*a3*b8-63*a3*b8^2+12201*b7*b5*b6+15183*b5^2*b8+7938*b7*b9*a9-3717*b5^3-3234*a3*b7*b6-630*b9*a3*b8+693*a7*b4*b9+2016*a3^2*b8+10836*a3*b5*a7-4200*b8*b7*b6+1764*b7*b8*b10-1701*b8*b5*a7-5922*a3^2*b5-6930*a7*b7*b6-378*b4*b8^2-8568*a3*a7*b9-3087*a3*b9^2+4662*a3^2*b9-1953*a7*b8*b9+8694*a7*b5*b9-3591*a7*b4*b5+2772*a7*b7*b10-8568*a3^2*b4-11277*a7*b5^2+1386*b4^2*b9+17892*b8*b5*b9-13671*b4*b9^2+1008*a3*b7*b10-37359*b4*b5^2+17136*a7*b4*a3+22932*a3*b5^2+4788*b8^2*b9-693*b4^2*a3-16254*b9^2*b5-756*b8^3+504*b4^2*a7-7938*a3*b8*b5-8568*b7*b5*b10+37485*b5*b9*b4+6300*a3*b4*b9-3843*b9^2*b8-13797*b8^2*b5+819*b4*a7*b8-1323*a3*b4*b5+1386*a3*b5*b9+4641*b7*b4*b6-630*b4^2*b5+882*a7^2*b8-2394*a7^2*b5+1890*a7^2*b9-4536*a7^2*b4+252*b7*b8*a4+4599*b7*b4*a4-2898*b7*b4*b10-5859*b5^2*b9+5355*b9^3+4011*b9*b7*b6-252*b9*b7*b10-378*b4*b8*b5-819*b4*a3*b8+1827*b4*b8*b9+6727*b7*b5*a5"),
		SMQP("1008*a7^2*b7-756*a7*b7*b4+168*b7^2*a9-3276*a3*a7*b7+819*a3*b7*b4-1428*a3*b7*b8+2184*b7*b4*b5+504*a3^2*b7+161*b7^2*a5-714*b7*b4*b9-399*a3*b7*b9+1848*b7*b5^2-1113*b7*b8*b5+567*a7*b7*b9+378*b7*b8^2+84*b7^2*b10+756*b7^2*a8-609*b7^2*a4+105*b7*b9^2-588*b7*b8*b4-1239*b7*b5*b9+441*b7*b8*b9+252*b7*b4^2+3969*a3*b7*b5-1575*a7*b7*b5+525*a7*b7*b8-63*b7^2*b6"),
		SMQP("497*b7*b4*a5-483*b7*b9*a5-735*b7*b4*a9-651*b7*b8*a5+1008*a7*b9^2+588*a7*b8^2+2247*b7^2*a6+2226*b7*b8*a9-14448*b7*b5*a9-1176*b7^2*a10-2877*a7*b7*a8-3780*a7*b7*a5+1029*a7*b7*a4-2352*b7*b4*a8+1029*b7*b5*a8+4557*b7*a9*a7+7938*b7*b5*a4+840*b7*a8*a3+357*b7*a4*a3-609*b7*a5*a3+987*b7*a9*a3-4179*b7*b9*a4-378*b7*b8*a8+1197*b7*b9*a8-2856*a7*a3*b8+336*a3*b8^2+6426*b7*b5*b6+8232*b5^2*b8+3843*b7*b9*a9-1386*b5^3-2520*a3*b7*b6+252*b9*a3*b8+336*a7*b4*b9+1512*a3^2*b8+3444*a3*b5*a7-1764*b8*b7*b6+672*b7*b8*b10-1302*b8*b5*a7-2772*a3^2*b5-3528*a7*b7*b6-2520*a3*a7*b9-1008*a3*b9^2+1512*a3^2*b9-2604*a7*b8*b9+5922*a7*b5*b9+168*a7*b4*b5+1344*a7*b7*b10-4032*a3^2*b4-6846*a7*b5^2+10962*b8*b5*b9-6048*b4*b9^2+1008*a3*b7*b10-21546*b4*b5^2+9408*a7*b4*a3+13860*a3*b5^2+2604*b8^2*b9-8946*b9^2*b5-168*b8^3-5082*a3*b8*b5-4200*b7*b5*b10+20412*b5*b9*b4+1008*a3*b4*b9-2772*b9^2*b8-8400*b8^2*b5+756*a3*b4*b5-1134*a3*b5*b9+2268*b7*b4*b6+1344*a7^2*b8-504*a7^2*b5-3024*a7^2*b4-420*b7*b8*a4+1344*b7*b4*a4-1008*b7*b4*b10-3654*b5^2*b9+3024*b9^3+2016*b9*b7*b6-336*b9*b7*b10+336*b4*b8*b5-588*b4*a3*b8+168*b4*b8*b9+4039*b7*b5*a5"),
		SMQP("1575*b7*b4*a5-3633*b7*b9*a5-945*b7*b4*a9-3297*b7*b8*a5+2268*a7*b9^2-567*b7^2*a6+126*b7*b8*a9+1680*b7*b5*a9+504*b7^2*a10+693*a7*b7*a8-1764*a7*b7*a5-1197*a7*b7*a4-861*b7*b5*a8+315*b7*a9*a7+1554*b7*b5*a4-1512*b7*a8*a3+819*b7*a4*a3-5103*b7*a5*a3-1827*b7*a9*a3+1575*b7*b9*a4+1386*b7*b8*a8+315*b7*b9*a8-3528*a7*a3*b8+1554*b7*b5*b6+1512*b5^2*b8-1323*b7*b9*a9-882*b5^3+168*a3*b7*b6-3024*b9*a3*b8-1008*a7*b4*b9+1512*a3^2*b8-1260*a3*b5*a7+336*b8*b7*b6+630*b8*b5*a7+252*a3^2*b5+504*a7*b7*b6+3528*a3*a7*b9+756*a3*b9^2-1512*a3^2*b9-1260*a7*b8*b9-4914*a7*b5*b9+1512*a7*b4*b5+1386*a7*b5^2-630*b8*b5*b9+3528*b4*b9^2+1134*b4*b5^2-3276*a3*b5^2-3528*b8^2*b9-1638*b9^2*b5+1890*a3*b8*b5-1512*b7*b5*b10-5796*b5*b9*b4+1260*a3*b4*b9+6804*b9^2*b8+252*b8^2*b5-2268*a3*b4*b5+8190*a3*b5*b9-84*b7*b4*b6+1008*a7^2*b8+504*a7^2*b5-1008*a7^2*b9-504*b7*b8*a4+882*b5^2*b9-3276*b9^3-3612*b9*b7*b6+1008*b9*b7*b10+9373*b7*b5*a5"),
		SMQP("13412*a5*b5*a8+2940*b4*a5*a4-1414*b8*b4*a10-1890*b9*b6*a9+3332*b10*b4*a5+378*b7*a10*a5+1827*b10*b9*a9+2184*b8*a6*a3+343*b8*a9*a4+630*b7*a10*a4-3276*b7*a6*a4+2058*b8*b6*a5+2275*b5*a4*a9+924*b6*b5*a5+8316*b7*a6*a5+273*b8*a8*a9-238*b8*a8^2-5040*a6*b7*a8-1414*b8*a8*a4-2296*b5*a8^2-1610*b8*a10*a3+4032*b5*a4^2-3528*b4*a4^2-287*b4*a9^2+189*a8*b9*a9-1302*a8*b4*a9-1008*b4*a10*a7+210*a8*b5*a9+630*b9*a10*a7+504*a8*b7*a10-497*a9^2*b5-1134*a9*b9*a4+630*a9^2*b9-882*b6*b9*a8+2926*b4*a10*a3+1260*b10*b9*a4+3248*b4^2*a10+1344*b4*b6*a4-3654*b7*b6*a10+11466*b7*b6*a6+2016*b7*b10*a10-1365*b8*b6*a9+203*b8*b10*a9+5628*b4*b6*a9+1288*a4*b4*a9-4312*b8*b5*a10+378*b8*b9*a10-4284*b7*b10*a6+1764*b9*a4^2-3976*b10*b5*a5-10332*b5^2*a6-4368*b4*a6*a3+7252*b8^2*a6+2912*b8*a6*a7+1120*b4*b10*a8+1925*b10*b5*a9-4284*b9*a6*a7-3528*b4*b9*a10+3276*b9*a6*a3-13790*b9*b8*a6+3836*b5^2*a10-994*b8*b10*a4-1358*b4*b8*a6-7028*b4*b5*a10-2513*b8*a5*a9-2240*b10*b5*a8+7609*b4*a5*a9-1540*b10*b5*a4-1386*b9*a10*a3-1820*b4*b9*a6-42*b8*b6*a4+56*b8*b6*a8-336*b9*b10*a5+308*b4*b6*a8+3248*b4*a8*a4-4060*b8*b10*a5+16212*b4*b5*a6-1176*b4^2*a6-2751*a5*b9*a9-112*b4*a6*a7-4368*a5^2*b4+3087*b7*a6*a9+1008*b9*b10*a8+3822*b9*b6*a5-1260*b9*b6*a4-5404*b4*a5*a8+882*b8*a5*a4+756*b8^2*a10-1764*b9*a8*a4-2184*b9*a5*a4+406*b8*a9^2+392*b4*b10*a4+3276*b8*a5^2+1386*b8*a4^2-7182*b4*b6*a5+392*b4*a8^2+4788*b9*b6^2+756*b8*b10^2-1260*b4*b6^2-3780*b5*b10^2+2268*b10^2*b9-5418*b9*b6*b10-7434*b8*b10*b6+2520*b4*b6*b10-13860*b6^2*b5+9954*b6^2*b8+11466*b5*b6*b10-15456*b5*a6*a3+9044*b5*a6*a7+1120*b5*b9*a6-700*b5*a8*a4-3458*b5*a10*a3+378*b5*a10*a7+126*b5*b9*a10-8820*b5*a5*a4+16870*b5*b8*a6-6818*a5*b5*a9-252*b9*a8^2+2100*a5*b9*a8-28*a5*b8*a8+882*b8*a10*a7-3640*b10*b4*a9+1428*b9*a5^2-1050*b5*b6*a4+4032*b9^2*a6+378*b9^2*a10+8456*b5*b6*a8-518*b8*b10*a8+2520*a5^2*b5-2184*b5*b6*a9"),
		SMQP("-42*a5*b5*a8+1092*b4*a5*a4+364*b10*b4*a5+252*b7*a6*a4-273*b8*b6*a5+42*b6*b5*a5-63*b7*a6*a5-756*b4*a4^2+56*b4*a9^2-84*a8*b4*a9+308*b4*a10*a3+700*b4^2*a10+672*b4*b6*a4-252*b7*b6*a6-84*b4*b6*a9+56*a4*b4*a9-21*b10*b5*a5+1092*b4*a6*a3-504*b8^2*a6-252*b8*a6*a7-28*b4*b10*a8+252*b9*a6*a7-252*b9*a6*a3+756*b9*b8*a6-616*b4*b8*a6-308*b4*b5*a10-42*b8*a5*a9+1547*b4*a5*a9+1148*b4*b9*a6+21*b9*b10*a5-644*b4*b6*a8+700*b4*a8*a4-21*b8*b10*a5-1764*b4*b5*a6-84*b4^2*a6+42*a5*b9*a9-728*b4*a6*a7-609*a5^2*b4+210*b9*b6*a5-2114*b4*a5*a8+273*b8*a5*a4-210*b9*a5*a4+28*b4*b10*a4-315*b8*a5^2+126*b4*b6*a5+28*b4*a8^2+756*b5*a6*a3-252*b5*a6*a7-252*b5*b9*a6+21*b5*a5*a4+252*b5*b8*a6+21*a5*b5*a9-21*a5*b9*a8+21*a5*b8*a8+28*b10*b4*a9+147*b9*a5^2-252*b9^2*a6+126*a5^2*b5"),
		SMQP("-2499*b9^2*a9+6447*b9^2*a4+525*b9^2*a8-1512*b4^2*b10-1694*b7*b4*a10+6202*b7*b10*a9-9646*b7*a10*a3-651*b7*b4*a6-6398*b7*b10*a4-13419*b4*a4*a3-9324*b8*a8*a3-5397*b8*b9*a4+5166*b8*a4*a3+1659*b8*a9*a7-14532*b7*a6*a3+27405*b4*a8*a3-3360*b4*a8*a7+8946*b5^2*a5+9681*b4*b5*a4-854*a7*b7*a6-6783*a7*b5*a9+3024*a7*b5*a5-1008*b8*b7*a10-9366*a7*b4*a5+3528*a7*b7*a10+2562*a7*b4*a9+12978*a7*b8*a5-1239*a7*b8*a8-3276*b8*a9*a3+9072*b8*a5*a3-9219*b8*b4*a9+17598*b8*b4*a5-147*b8*b5*a4-2163*b4*b8*a8+7455*b4*b9*a8+1176*b4*a4*a7-22764*b4*b5*a8-11193*a7*b8*a4+10521*b7*a4^2-5537*b7*a4*a8+9765*b8^2*a5-15372*a3*b4*a5+2457*a3*b4*a9-24325*b7*a5*a9+4830*b9^2*a5+18571*b7*a5*a8-13545*b5*a4*a3-5474*b7*b5*a10+8211*b5*a4*a7+8841*b5^2*a4-5082*b5^2*a8+3927*b5^2*a9-924*b4*b9*a4+8463*b5*b9*a5-13671*b5*b8*a5+7203*b5*b4*a5-1827*b4^2*a9+13902*b5*b8*a8-8239*b7*b8*a6-6867*b4^2*a5+4074*b7*a9*a8+8358*b7*b6*a5+45045*b7*b5*a6+4025*b7*a9*a4+7287*b7*b6*a4-15771*b7*b6*a9+3626*b7*b10*a8-5838*b5*b9*a8+1309*b7*b6*a8+1316*b7*a9^2+16884*b7*a5^2-19376*b7*b10*a5-6027*b5*b9*a9-4886*b7*a8^2-10101*a5*b7*a4-420*b8^2*a9-5901*b8^2*a8-9177*b9*b5*a4+12159*b5*a9*a3-7644*b5*b8*a9+8988*b5*b4*a9-840*b9^2*b10+3906*b4^2*b6-14154*b4*b5*b6+231*b8^2*b6+1890*a3*b8*b6+11970*a3*b9*b6-1764*a3*b8*b10+8568*a3*b10*b5-3528*a3*b10*b9-5796*a3*b5*b6-6762*b9^2*b6+1512*b4^2*a4+3780*a3*b10*b4+6930*b4^2*a8+2898*b7*b6^2+252*b7*b10^2+19257*b10*b9*b8+16485*b5^2*b10-20244*b5^2*b6-8085*b8^2*b10+8967*b5*b6*b8+25620*b5*b6*b9-26481*b10*b9*b5-4284*b7*b6*b10-16632*a3*b4*b6-13944*b8*b9*b6+3003*b4*b10*b9-4767*b4*b10*b5-6783*b4*b8*b10+11382*b4*b9*b6+5628*a7*b4*b6-8610*a7*b9*b6+6006*a7*b10*b9-2730*a7*b10*b5-5754*a7*b8*b10+12012*a7*b8*b6-4956*a7*b10*b4-2730*a7*b5*b6-105*b8^2*a4+15267*b4*b8*b6+1491*b4*b8*a4+2100*b5*b8*b10-3843*b9*a9*a3+5418*b9*a5*a3+9849*b9*b8*a9+6006*b9*b8*a8-8589*b9*b4*a5-17871*b9*b8*a5-2814*b9*b4*a9+15939*b5*a8*a3-2625*b5*a8*a7-18396*b5*a5*a3-11634*b9*a5*a7-5985*b9*a8*a3+7623*b9*a4*a3+2751*b9*a9*a7+357*b9*a4*a7+3381*b9*a8*a7+4536*b9*b7*a10-4396*b9*b7*a6"),
		SMQP("-147*a6*b5*b6-252*a6*b5*a5+2667*a6*b4*b6+441*a6*b4*a5+126*a6*b5*a4+126*a6*b4*a4+1260*b6^3-231*b10*b5*a6+105*b8*a10*a5+392*b10*b5*a10-350*b10*b9*a10+70*a9*b5*a10+210*a6*b9*a4+252*a6*b7*a10+210*a6*b9*a9+210*a6*b9*a8+735*b10*b9*a6-574*b5*b6*a10-651*a6*b4*a9+126*a6*b8*a8-70*b9*a10*a8+14*a8*b4*a10+112*a8*b5*a10+77*a9*b4*a10-308*a10*b5*a4+14*a10*b9*a4-14*a10*b8*a8-42*a10^2*b7-112*a10*b9*a9+126*a8*b4*a6+154*a4*b8*a10-154*a4*b4*a10+70*b9*a10*a5-14*b8*a10*a9-63*b8*a6*a9-126*b8*b10*a6-84*a5*b9*a6-2205*b9*b6*a6-119*b4*a10*a5-21*a6*b5*a9-252*a6*b5*a8-630*a6^2*b7-735*b4*b10*a6+21*a5*b5*a10+924*b8*b6*a6-154*b8*b10*a10-378*b8*a6*a5-1967*b4*b6*a10+224*b8*b6*a10+763*b4*b10*a10-252*b10^3+1386*b10^2*b6-2394*b6^2*b10+784*b9*b6*a10"),
		SMQP("-1967*b9^2*a9-1169*b9^2*a4-2345*b9^2*a8+1764*b7*b4*a10-4788*b7*b10*a9+588*b7*a10*a3-3108*b7*b4*a6+5292*b7*b10*a4+6916*b4*a4*a3+3500*b8*a8*a3+3213*b8*b9*a4-3409*b8*a4*a3-3031*b8*a9*a7+1764*b7*a6*a3-7154*b4*a8*a3+994*b4*a8*a7-1470*b5^2*a5+882*b4*b5*a4+2436*a7*b7*a6+1148*a7*b5*a9-4032*a7*b5*a5-336*b8*b7*a10-3388*a7*b4*a5-2016*a7*b7*a10+3346*a7*b4*a9-84*a7*b8*a5+287*a7*b8*a8+2177*b8*a9*a3-4515*b8*a5*a3-875*b8*b4*a9-2863*b8*b4*a5-3178*b8*b5*a4+1582*b4*b8*a8+364*b4*b9*a8-4802*b4*a4*a7+5166*b4*b5*a8+1841*a7*b8*a4-1554*b7*a4^2+546*b7*a4*a8-1659*b8^2*a5+8078*a3*b4*a5-4172*a3*b4*a9+210*b7*a5*a9-4123*b9^2*a5+2058*b7*a5*a8+4844*b5*a4*a3+8484*b7*b5*a10-4480*b5*a4*a7-3234*b5^2*a4+2898*b5^2*a8-1470*b5^2*a9-2450*b4*b9*a4-651*b5*b9*a5-1869*b5*b8*a5-6258*b5*b4*a5-4536*b4^2*a9-3913*b5*b8*a8+1995*b7*b8*a6+1512*b4^2*a5-756*b7*a9*a8-504*b7*b6*a5-19656*b7*b5*a6+1638*b7*a9*a4-10122*b7*b6*a4+7182*b7*b6*a9-756*b7*b10*a8+1799*b5*b9*a8+966*b7*b6*a8-1512*b7*a9^2-4032*b7*a5^2+3360*b7*b10*a5+1988*b5*b9*a9+1260*b7*a8^2-1806*a5*b7*a4-322*b8^2*a9-952*b8^2*a8+1400*b9*b5*a4-1456*b5*a9*a3+5264*b5*b8*a9-12222*b5*b4*a9-5068*b9^2*b10+420*b4^2*b6+20244*b4*b5*b6+5614*b8^2*b6+952*a3*b8*b6-10528*a3*b9*b6+1372*a3*b8*b10+532*a3*b10*b5+2828*a3*b10*b9-1820*a3*b5*b6+8792*b9^2*b6-1960*a3*b10*b4+84*b4^2*a8-168*b7*b6^2-504*b7*b10^2+294*b10*b9*b8-5922*b5^2*b10+11088*b5^2*b6-1778*b8^2*b10-9842*b5*b6*b8-9044*b5*b6*b9+10150*b10*b9*b5+1512*b7*b6*b10+5348*a3*b4*b6-3024*b8*b9*b6+8582*b4*b10*b9-10206*b4*b10*b5+770*b4*b8*b10-19684*b4*b9*b6-2632*a7*b4*b6+392*a7*b9*b6-1204*a7*b10*b9-980*a7*b10*b5+700*a7*b8*b10+700*a7*b8*b6+1400*a7*b10*b4+2464*a7*b5*b6+14*b8^2*a4-1162*b4*b8*b6-266*b4*b8*a4+2212*b5*b8*b10+2401*b9*a9*a3-3829*b9*a5*a3-1239*b9*b8*a9+1785*b9*b8*a8-259*b9*b4*a5+10598*b9*b8*a5+8743*b9*b4*a9-5614*b5*a8*a3+1316*b5*a8*a7+9786*b5*a5*a3+4172*b9*a5*a7+2254*b9*a8*a3-2345*b9*a4*a3-1841*b9*a9*a7+2443*b9*a4*a7+553*b9*a8*a7-2688*b9*b7*a10+3927*b9*b7*a6"),
		SMQP("623*a6*b5*b6-168*a6*b5*a5-6503*a6*b4*b6-1267*a6*b4*a5-1036*a6*b5*a4-1778*a6*b4*a4-2940*b6^3+497*b10*b5*a6-378*b8*a10*a5-763*b10*b5*a10+721*b10*b9*a10-259*a9*b5*a10+406*a6*b9*a4-105*a6*b7*a10-518*a6*b9*a9-140*a6*b9*a8-1225*b10*b9*a6+1092*b5*b6*a10+1477*a6*b4*a9-700*a6*b8*a8+35*b9*a10*a8-28*a8*b4*a10-126*a8*b5*a10-364*a9*b4*a10+735*a10*b5*a4-196*a10*b9*a4+49*a10*b8*a8+42*a10^2*b7+182*a10*b9*a9+350*b8*a6*a4-350*a8*b4*a6-539*a4*b8*a10+1106*a4*b4*a10+7*b9*a10*a5+280*b8*a10*a9+119*b8*a6*a9-476*b8*b10*a6-742*a5*b9*a6+3773*b9*b6*a6+364*b4*a10*a5+259*a6*b5*a9+560*a6*b5*a8+420*a6^2*b7+2723*b4*b10*a6+161*a5*b5*a10-1162*b8*b6*a6+455*b8*b10*a10+1764*b8*a6*a5+4249*b4*b6*a10-623*b8*b6*a10-1883*b4*b10*a10+588*b10^3-3234*b10^2*b6+5586*b6^2*b10-1498*b9*b6*a10"),
		SMQP("-21*b7*b4*a5+231*b7*b9*a5-5313*b7*b4*a9-8001*b7*b8*a5+2772*a7*b9^2-2268*a7*b8^2+3465*b7^2*a6+1470*b7*b8*a9-29484*b7*b5*a9-1008*b7^2*a10-4851*a7*b7*a8-3276*a7*b7*a5-945*a7*b7*a4-3360*b7*b4*a8+12075*b7*b5*a8+9387*b7*a9*a7+14616*b7*b5*a4-8064*b7*a8*a3+6363*b7*a4*a3+693*b7*a5*a3+1197*b7*a9*a3-12621*b7*b9*a4-3066*b7*b8*a8+5775*b7*b9*a8-756*a7*a3*b8-4032*a3*b8^2+11928*b7*b5*b6+11340*b5^2*b8+8925*b7*b9*a9-4788*b5^3-2856*a3*b7*b6+756*b9*a3*b8+1260*a7*b4*b9-504*a3^2*b8+9576*a3*b5*a7-6048*b8*b7*b6+1428*b7*b8*b10-3276*b8*b5*a7-3024*a3^2*b5-5796*a7*b7*b6-2520*b4*b8^2-10836*a3*a7*b9-7056*a3*b9^2+6048*a3^2*b9-1008*a7*b8*b9+5292*a7*b5*b9-5544*a7*b4*b5+2520*a7*b7*b10-7056*a3^2*b4-8316*a7*b5^2+4536*b4^2*b9+9828*b8*b5*b9-17388*b4*b9^2-27468*b4*b5^2+14616*a7*b4*a3+16380*a3*b5^2+9324*b8^2*b9-252*b4^2*a3-14112*b9^2*b5-5544*b8^3+2016*a3*b8*b5-8148*b7*b5*b10+27972*b5*b9*b4+14364*a3*b4*b9-3276*b9^2*b8-4284*b8^2*b5+756*b4*a7*b8-2520*a3*b4*b5+7056*a3*b5*b9+1260*b7*b4*b6-6552*b4^2*b5-252*a7^2*b8-2268*a7^2*b5+2772*a7^2*b9-3528*a7^2*b4+5376*b7*b8*a4+4452*b7*b4*a4-1176*b7*b4*b10-4032*b5^2*b9+4536*b9^3+2016*b9*b7*b6+588*b9*b7*b10-10332*b4*b8*b5-1764*b4*a3*b8+13356*b4*b8*b9+1099*b7*b5*a5"),
		SMQP("-105*b9^2*a9+357*b9^2*a4-315*b9^2*a8-364*b7*b4*a10-2044*b7*b10*a9+700*b7*a10*a3-483*b7*b4*a6+1988*b7*b10*a4+2079*b4*a4*a3+840*b8*a8*a3+357*b8*b9*a4-819*b8*a4*a3-735*b8*a9*a7+1848*b7*a6*a3-2184*b4*a8*a3-567*b4*a8*a7+1197*b5^2*a5+588*b4*b5*a4+2072*a7*b7*a6-231*a7*b5*a9-2100*a7*b5*a5-672*b8*b7*a10-1848*a7*b4*a5-672*a7*b7*a10+1911*a7*b4*a9-84*a7*b8*a5+315*a7*b8*a8+651*b8*a9*a3-2541*b8*a5*a3+147*b8*b4*a9-504*b8*b4*a5-336*b8*b5*a4+210*b4*b8*a8+567*b4*b9*a8-777*b4*a4*a7+1239*b4*b5*a8+357*a7*b8*a4+756*b7*a4^2-700*b7*a4*a8-2415*b8^2*a5+3549*a3*b4*a5-1491*a3*b4*a9-2660*b7*a5*a9-3087*b9^2*a5+3248*b7*a5*a8+1701*b5*a4*a3+2324*b7*b5*a10-2163*b5*a4*a7-1344*b5^2*a4+1197*b5^2*a8-84*b5^2*a9-777*b4*b9*a4+378*b5*b9*a5-1050*b5*b8*a5-5040*b5*b4*a5-777*b4^2*a9+399*b5*b8*a8+1939*b7*b8*a6-105*b4^2*a5+84*b7*a9*a8+924*b7*b6*a5-7245*b7*b5*a6-308*b7*a9*a4-4704*b7*b6*a4+2352*b7*b6*a9+28*b7*b10*a8-210*b5*b9*a8+1652*b7*b6*a8-56*b7*a9^2-1596*b7*a5^2+1064*b7*b10*a5+231*b5*b9*a9-28*b7*a8^2-336*a5*b7*a4+546*b8^2*a9-42*b8^2*a8-567*b9*b5*a4-609*b5*a9*a3+966*b5*b8*a9-2079*b5*b4*a9-2016*b9^2*b10+6720*b4*b5*b6+4032*b8^2*b6+1512*a3*b8*b6-5040*a3*b9*b6+3024*a3*b10*b5-1512*a3*b5*b6+2016*b9^2*b6+2016*b10*b9*b8-1932*b5^2*b10+2352*b5^2*b6-2016*b8^2*b10-6468*b5*b6*b8+3192*b5*b6*b9+1932*b10*b9*b5+2520*a3*b4*b6-3024*b8*b9*b6+4032*b4*b10*b9-6384*b4*b10*b5-7056*b4*b9*b6-1008*a7*b4*b6-1008*a7*b10*b5+1008*a7*b8*b6+1008*a7*b5*b6+2100*b5*b8*b10+903*b9*a9*a3-1449*b9*a5*a3-1071*b9*b8*a9-357*b9*b8*a8+2604*b9*b4*a5+6090*b9*b8*a5+630*b9*b4*a9-3192*b5*a8*a3+1827*b5*a8*a7+6447*b5*a5*a3+2436*b9*a5*a7+840*b9*a8*a3-819*b9*a4*a3-1239*b9*a9*a7+357*b9*a4*a7+315*b9*a8*a7-665*b9*b7*a6"),
		SMQP("-161*b7*b4*a5-469*b7*b9*a5-161*b7*b4*a9-497*b7*b8*a5+308*a7*b9^2-308*a7*b8^2+273*b7^2*a6-210*b7*b8*a9-168*b7*b5*a9-63*a7*b7*a8+336*a7*b7*a5-49*a7*b7*a4+791*b7*b5*a8+63*b7*a9*a7+28*b7*b5*a4-616*b7*a8*a3+231*b7*a4*a3-63*b7*a5*a3-63*b7*a9*a3-301*b7*b9*a4-70*b7*b8*a8+231*b7*b9*a8-168*a3*b8^2+56*b7*b5*b6-308*b5^2*b8+49*b7*b9*a9-140*b5^3-476*b9*a3*b8+224*a7*b4*b9+56*a3*b5*a7-140*b8*b7*b6-56*b7*b8*b10+28*b8*b5*a7-168*a3^2*b5-196*a3*b9^2-644*a7*b5*b9-336*a7*b4*b5+280*a7*b5^2-840*b8*b5*b9-56*b4*b9^2+1092*b4*b5^2-700*a3*b5^2+308*b8^2*b9-84*b9^2*b5-616*b8^3+896*a3*b8*b5-56*b7*b5*b10-1288*b5*b9*b4+700*a3*b4*b9+616*b9^2*b8+1148*b8^2*b5-84*a3*b4*b5+1484*a3*b5*b9-56*b7*b4*b6+588*b7*b8*a4+112*b5^2*b9-308*b9^3-196*b9*b7*b6+56*b9*b7*b10-1512*b4*b8*b5-84*b4*a3*b8+1064*b4*b8*b9-273*b7*b5*a5"),
		SMQP("504*b7*b4*a5+1554*b7*b9*a5-504*b7*b4*a9-1806*b7*b8*a5+504*a7*b9^2-252*b7^2*a6-1176*b7*b5*a9+252*b7^2*a10+378*a7*b7*a8-630*a7*b7*a5-630*a7*b7*a4-714*b7*b5*a8+126*b7*a9*a7-609*b7*b5*a4-504*b7*a8*a3+630*b7*a4*a3-252*b7*a5*a3+1134*b7*a9*a3+126*b7*b9*a4+756*b7*b8*a8+126*b7*b9*a8-1764*a7*a3*b8+2037*b7*b5*b6-3024*b5^2*b8-630*b7*b9*a9-441*b5^3+84*a3*b7*b6-252*b9*a3*b8-504*a7*b4*b9+756*a3^2*b8-630*a3*b5*a7+168*b8*b7*b6+1575*b8*b5*a7+126*a3^2*b5+252*a7*b7*b6+1764*a3*a7*b9-1512*a3*b9^2-756*a3^2*b9-3087*a7*b5*b9+756*a7*b4*b5+1953*a7*b5^2-5985*b8*b5*b9+504*b4*b9^2+5607*b4*b5^2-5418*a3*b5^2-504*b8^2*b9+3591*b9^2*b5+3465*a3*b8*b5-756*b7*b5*b10-4158*b5*b9*b4+1260*a3*b4*b9+1512*b9^2*b8+2646*b8^2*b5-2394*a3*b4*b5+2205*a3*b5*b9-42*b7*b4*b6+504*a7^2*b8+252*a7^2*b5-504*a7^2*b9-252*b7*b8*a4+2961*b5^2*b9-1008*b9^3-1176*b9*b7*b6+504*b9*b7*b10+5663*b7*b5*a5"),
		SMQP("-1582*a5*b5*a8-11823*b4*a5*a4+4340*b8*b4*a10-9954*b9*b6*a9+4312*b10*b4*a5-4284*b7*a10*a5+3059*b10*b9*a9-336*b8*a6*a3-3059*b8*a9*a4-4284*b7*a10*a4+13923*b7*a6*a4+3024*b7*a10*a9+13776*b8*b6*a5-602*b5*a4*a9-3864*b6*b5*a5+3528*b7*a6*a5+1113*b8*a8*a9-637*b8*a8^2-1449*a6*b7*a8-5194*b8*a8*a4+1694*b5*a8^2-1820*b8*a10*a3-2373*b5*a4^2+420*b4*a4^2-4207*b4*a9^2+315*a8*b9*a9+3171*a8*b4*a9+1008*b4*a10*a7-2415*a8*b5*a9+2268*b9*a10*a7+504*a8*b7*a10-497*a9^2*b5+196*a9*b9*a4+1414*a9^2*b9-5362*b6*b9*a8-2128*b4*a10*a3-1897*b10*b9*a4+18508*b4^2*a10-714*b4*b6*a4-6300*b7*b6*a10+18396*b7*b6*a6+3528*b7*b10*a10-1155*b8*b6*a9+833*b8*b10*a9+13062*b4*b6*a9+5537*a4*b4*a9+56*b8*b5*a10+6300*b8*b9*a10-8820*b7*b10*a6-3990*b9*a4^2+5600*b10*b5*a5+6300*b5^2*a6-3696*b4*a6*a3-13244*b8^2*a6-2296*b8*a6*a7-5320*b4*b10*a8-2443*b10*b5*a9-2632*b9*a6*a7-20692*b4*b9*a10+3696*b9*a6*a3+16688*b9*b8*a6+4340*b5^2*a10-3787*b8*b10*a4-21364*b4*b8*a6+22344*b4*b5*a10-1337*b8*a5*a9-3479*b10*b5*a8+2177*b4*a5*a9+6797*b10*b5*a4-1232*b9*a10*a3+43148*b4*b9*a6+7581*b8*b6*a4+3857*b8*b6*a8+1232*b9*b10*a5-350*b4*b6*a8+70*b4*a8*a4-10276*b8*b10*a5-30156*b4*b5*a6-20496*b4^2*a6-6384*b10^2*b4+1057*a5*b9*a9+5488*b4*a6*a7-5880*a5^2*b4-8253*b7*a6*a9+1855*b9*b10*a8+4368*b9*b6*a5+2478*b9*b6*a4+4417*b4*a5*a8+5103*b8*a5*a4-2520*b8^2*a10+371*b9*a8*a4+2709*b9*a5*a4+1330*b8*a9^2-1904*b4*b10*a4+4536*b8*a5^2+1407*b8*a4^2-11676*b4*b6*a5+238*b4*a8^2+23268*b9*b6^2+2352*b8*b10^2-13776*b4*b6^2-10752*b5*b10^2+6720*b10^2*b9-24948*b9*b6*b10-7812*b8*b10*b6+17640*b4*b6*b10-31584*b6^2*b5+8988*b6^2*b8+36288*b5*b6*b10+2688*b5*a6*a3-280*b5*a6*a7+7672*b5*b9*a6+5705*b5*a8*a4+8512*b5*a10*a3-3276*b5*a10*a7-13888*b5*b9*a10+3150*b5*a5*a4-3416*b5*b8*a6+2002*a5*b5*a9+497*b9*a8^2-8659*a5*b9*a8+1631*a5*b8*a8-756*b8*a10*a7-980*b10*b4*a9-1848*b9*a5^2-4578*b5*b6*a4-12908*b9^2*a6+3780*b9^2*a10+3290*b5*b6*a8+1309*b8*b10*a8-2016*a5^2*b5+3066*b5*b6*a9"),
		SMQP("231*b9^2*a9+15477*b9^2*a4-903*b9^2*a8-1680*b4^2*b10-22610*b7*b4*a10+3934*b7*b10*a9-9394*b7*a10*a3+27069*b7*b4*a6-11186*b7*b10*a4-23205*b4*a4*a3-16422*b8*a8*a3-25557*b8*b9*a4+17472*b8*a4*a3+861*b8*a9*a7-7350*b7*a6*a3+23415*b4*a8*a3-7728*b4*a8*a7+8694*b5^2*a5+23625*b4*b5*a4-602*a7*b7*a6-4305*a7*b5*a9+6552*a7*b5*a5+4032*b8*b7*a10+2310*a7*b4*a5+4536*a7*b7*a10+1470*a7*b4*a9+1638*a7*b8*a5+1239*a7*b8*a8-2184*b8*a9*a3+9702*b8*a5*a3-4767*b8*b4*a9+36918*b8*b4*a5+1029*b8*b5*a4-6783*b4*b8*a8-861*b4*b9*a8+7896*b4*a4*a7-13104*b4*b5*a8-399*a7*b8*a4+13167*b7*a4^2-18263*b7*a4*a8-567*b8^2*a5-14154*a3*b4*a5+3129*a3*b4*a9-22939*b7*a5*a9+5880*b9^2*a5+9877*b7*a5*a8-34881*b5*a4*a3-11774*b7*b5*a10+17997*b5*a4*a7+2163*b5^2*a4+1218*b5^2*a8+1533*b5^2*a9+1848*b4*b9*a4+8337*b5*b9*a5-4977*b5*b8*a5+2037*b5*b4*a5+1617*b4^2*a9+7770*b5*b8*a8+8141*b7*b8*a6-25935*b4^2*a5+3822*b7*a9*a8+4578*b7*b6*a5+43911*b7*b5*a6+2639*b7*a9*a4+18249*b7*b6*a4-15141*b7*b6*a9+7406*b7*b10*a8-11970*b5*b9*a8-581*b7*b6*a8+2828*b7*a9^2+18900*b7*a5^2-18368*b7*b10*a5-6363*b5*b9*a9-1106*b7*a8^2+4389*a5*b7*a4-1512*b8^2*a9-2205*b8^2*a8-189*b9*b5*a4+8463*b5*a9*a3-5334*b5*b8*a9+12978*b5*b4*a9-798*b9^2*b10+3738*b4^2*b6-52290*b4*b5*b6-15309*b8^2*b6+6216*a3*b8*b6+13314*a3*b9*b6-6762*a3*b8*b10+9618*a3*b10*b5-2562*a3*b10*b9+168*a3*b5*b6-12138*b9^2*b6+1680*b4^2*a4+5628*a3*b10*b4+16674*b4^2*a8+5166*b7*b6^2+2772*b7*b10^2+18627*b10*b9*b8+22407*b5^2*b10-27804*b5^2*b6-945*b8^2*b10+28329*b5*b6*b8+28728*b5*b6*b9-27909*b10*b9*b5-10332*b7*b6*b10-22092*a3*b4*b6-9030*b8*b9*b6-273*b4*b10*b9+9261*b4*b10*b5-15771*b4*b8*b10+19698*b4*b9*b6+8484*a7*b4*b6+546*a7*b9*b6+6594*a7*b10*b9-4830*a7*b10*b5-3822*a7*b8*b10-2940*a7*b8*b6-3108*a7*b10*b4-3318*a7*b5*b6+6111*b8^2*a4+40047*b4*b8*b6-5145*b4*b8*a4-8022*b5*b8*b10-3297*b9*a9*a3+6384*b9*a5*a3+6951*b9*b8*a9+12054*b9*b8*a8-7959*b9*b4*a5-17829*b9*b8*a5-11760*b9*b4*a9+21945*b5*a8*a3-6951*b5*a8*a7-16632*b5*a5*a3-7014*b9*a5*a7-7077*b9*a8*a3+10437*b9*a4*a3+3297*b9*a9*a7-13461*b9*a4*a7+4179*b9*a8*a7+6552*b9*b7*a10-24178*b9*b7*a6"),
		SMQP("-4851*b9^2*a9-1785*b9^2*a4-6069*b9^2*a8+6006*b7*b4*a10-7602*b7*b10*a9+4326*b7*a10*a3-9366*b7*b4*a6+10542*b7*b10*a4+16002*b4*a4*a3+11298*b8*a8*a3+6860*b8*b9*a4-9975*b8*a4*a3-8722*b8*a9*a7+3402*b7*a6*a3-13713*b4*a8*a3+1897*b4*a8*a7-3402*b5^2*a5+1617*b4*b5*a4+4746*a7*b7*a6+1211*a7*b5*a9-6552*a7*b5*a5-3024*b8*b7*a10-6454*a7*b4*a5-4200*a7*b7*a10+7609*a7*b4*a9+42*a7*b8*a5-1750*a7*b8*a8+4263*b8*a9*a3-9555*b8*a5*a3+2408*b8*b4*a9-7826*b8*b4*a5-7679*b8*b5*a4+2975*b4*b8*a8+9702*b4*b9*a8-11753*b4*a4*a7+1407*b4*b5*a8+7238*a7*b8*a4-5985*b7*a4^2+6657*b7*a4*a8+3486*b8^2*a5+16695*a3*b4*a5-8526*a3*b4*a9+3549*b7*a5*a9-5166*b9^2*a5-147*b7*a5*a8+9135*b5*a4*a3+14826*b7*b5*a10-8239*b5*a4*a7-5397*b5^2*a4+1218*b5^2*a8+21*b5^2*a9-6993*b4*b9*a4-4557*b5*b9*a5-6888*b5*b8*a5-5460*b5*b4*a5-13356*b4^2*a9-5621*b5*b8*a8+9492*b7*b8*a6+5796*b4^2*a5-2562*b7*a9*a8-2646*b7*b6*a5-46641*b7*b5*a6+4767*b7*a9*a4-19719*b7*b6*a4+13083*b7*b6*a9-3234*b7*b10*a8+5880*b5*b9*a8+483*b7*b6*a8-3444*b7*a9^2-8820*b7*a5^2+9912*b7*b10*a5+3549*b5*b9*a9+2142*b7*a8^2-3843*a5*b7*a4-1190*b8^2*a9-6671*b8^2*a8+2331*b9*b5*a4+567*b5*a9*a3+10738*b5*b8*a9-30618*b5*b4*a9-7518*b9^2*b10+1890*b4^2*b6+43806*b4*b5*b6+11795*b8^2*b6-2016*a3*b8*b6-21294*a3*b9*b6+3486*a3*b8*b10-3318*a3*b10*b5+6174*a3*b10*b9+1008*a3*b5*b6+15078*b9^2*b6-1596*a3*b10*b4+378*b4^2*a8-1890*b7*b6^2-1932*b7*b10^2-4781*b10*b9*b8-15225*b5^2*b10+26124*b5^2*b6-2401*b8^2*b10-25543*b5*b6*b8-22344*b5*b6*b9+22323*b10*b9*b5+6132*b7*b6*b10+10332*a3*b4*b6+154*b8*b9*b6+16107*b4*b10*b9-17031*b4*b10*b5+2905*b4*b8*b10-40950*b4*b9*b6-6580*a7*b4*b6+2282*a7*b9*b6-3850*a7*b10*b9-98*a7*b10*b5+1918*a7*b8*b10-560*a7*b8*b6+4508*a7*b10*b4+4690*a7*b5*b6+847*b8^2*a4-6797*b4*b8*b6-637*b4*b8*a4+7322*b5*b8*b10+2709*b9*a9*a3-7098*b9*a5*a3-238*b9*b8*a9+6755*b9*b8*a8-4410*b9*b4*a5+14518*b9*b8*a5+15099*b9*b4*a9-5607*b5*a8*a3-175*b5*a8*a7+15120*b5*a5*a3+8078*b9*a5*a7+315*b9*a8*a3-1785*b9*a4*a3-413*b9*a9*a7+3829*b9*a4*a7+2317*b9*a8*a7-5208*b9*b7*a10+6972*b9*b7*a6"),
		SMQP("-4627*b9^2*a9-217*b9^2*a4-3493*b9^2*a8+3066*b7*b4*a10-8694*b7*b10*a9+378*b7*a10*a3-5859*b7*b4*a6+10122*b7*b10*a4+12831*b4*a4*a3+8526*b8*a8*a3+5215*b8*b9*a4-7602*b8*a4*a3-6363*b8*a9*a7+4158*b7*a6*a3-13251*b4*a8*a3+546*b4*a8*a7-2142*b5^2*a5+1071*b4*b5*a4+6006*a7*b7*a6+1449*a7*b5*a9-7056*a7*b5*a5-2016*b8*b7*a10-7602*a7*b4*a5-3864*a7*b7*a10+7644*a7*b4*a9+1638*a7*b8*a5+147*a7*b8*a8+3402*b8*a9*a3-8316*b8*a5*a3-147*b8*b4*a9-4830*b8*b4*a5-5145*b8*b5*a4+2499*b4*b8*a8+3731*b4*b9*a8-8274*b4*a4*a7+10206*b4*b5*a8+4557*a7*b8*a4-2331*b7*a4^2+1155*b7*a4*a8+441*b8^2*a5+14700*a3*b4*a5-6867*a3*b4*a9-105*b7*a5*a9-5600*b9^2*a5+4935*b7*a5*a8+8085*b5*a4*a3+16758*b7*b5*a10-7917*b5*a4*a7-5607*b5^2*a4+1638*b5^2*a8-777*b5^2*a9-5530*b4*b9*a4-2877*b5*b9*a5-7749*b5*b8*a5-3843*b5*b4*a5-7875*b4^2*a9-6972*b5*b8*a8+6657*b7*b8*a6+2205*b4^2*a5-966*b7*a9*a8-1386*b7*b6*a5-40635*b7*b5*a6+3045*b7*a9*a4-19341*b7*b6*a4+12873*b7*b6*a9-1638*b7*b10*a8+3682*b5*b9*a8+2793*b7*b6*a8-2940*b7*a9^2-6804*b7*a5^2+5712*b7*b10*a5+4879*b5*b9*a9+1722*b7*a8^2-4473*a5*b7*a4-924*b8^2*a9-1995*b8^2*a8+2905*b9*b5*a4-651*b5*a9*a3+10878*b5*b8*a9-27678*b5*b4*a9-9170*b9^2*b10+630*b4^2*b6+39186*b4*b5*b6+11697*b8^2*b6+1176*a3*b8*b6-18186*a3*b9*b6+1890*a3*b8*b10+2982*a3*b10*b5+5082*a3*b10*b9-4872*a3*b5*b6+15442*b9^2*b6-3276*a3*b10*b4+126*b4^2*a8-630*b7*b6^2-1092*b7*b10^2+1337*b10*b9*b8-12075*b5^2*b10+23436*b5^2*b6-4683*b8^2*b10-25221*b5*b6*b8-16744*b5*b6*b9+16793*b10*b9*b5+2940*b7*b6*b10+9660*a3*b4*b6-5866*b8*b9*b6+17689*b4*b10*b9-19509*b4*b10*b5+1155*b4*b8*b10-36050*b4*b9*b6-4956*a7*b4*b6+1470*a7*b9*b6-1974*a7*b10*b9-2142*a7*b10*b5+882*a7*b8*b10+1176*a7*b8*b6+2436*a7*b10*b4+4326*a7*b5*b6+21*b8^2*a4-1743*b4*b8*b6-399*b4*b8*a4+8022*b5*b8*b10+3717*b9*a9*a3-6720*b9*a5*a3-497*b9*b8*a9+1624*b9*b8*a8-2723*b9*b4*a5+14847*b9*b8*a5+13034*b9*b4*a9-9933*b5*a8*a3+2667*b5*a8*a7+15624*b5*a5*a3+7098*b9*a5*a7+2121*b9*a8*a3-2793*b9*a4*a3-2877*b9*a9*a7+2961*b9*a4*a7+1533*b9*a8*a7-3864*b9*b7*a10+4242*b9*b7*a6"),
		SMQP("1540*a5*b5*a8+252*b4*a5*a4-938*b8*b4*a10+3472*b9*b6*a9-1106*b10*b4*a5-84*b7*a10*a5-1470*b10*b9*a9+1848*b8*a6*a3+1260*b8*a9*a4+1092*b7*a10*a4-3696*b7*a6*a4-1512*b7*a10*a9-3129*b8*b6*a5+336*b5*a4*a9-546*b6*b5*a5-588*b7*a6*a5-910*b8*a8*a9+938*b8*a8^2-1680*a6*b7*a8+3542*b8*a8*a4-1582*b5*a8^2+2030*b8*a10*a3+2142*b5*a4^2+1764*b4*a4^2+2604*b4*a9^2+742*a8*b9*a9-1484*a8*b4*a9+1008*b4*a10*a7+938*a8*b5*a9+336*b9*a10*a7+504*a8*b7*a10+420*a9^2*b5-56*a9*b9*a4-308*a9^2*b9+1883*b6*b9*a8-896*b4*a10*a3-1414*b10*b9*a4-2212*b4^2*a10-1596*b4*b6*a4+2436*b7*b6*a10-8253*b7*b6*a6-1512*b7*b10*a10+2828*b8*b6*a9-1162*b8*b10*a9-9359*b4*b6*a9-1176*a4*b4*a9-14*b8*b5*a10+1008*b8*b9*a10+4578*b7*b10*a6+1806*b9*a4^2-1316*b10*b5*a5-1890*b5^2*a6-1344*b4*a6*a3+5936*b8^2*a6+280*b8*a6*a7+672*b4*b10*a8+308*b10*b5*a9+1120*b9*a6*a7+5110*b4*b9*a10-1512*b9*a6*a3-4718*b9*b8*a6-602*b5^2*a10+1260*b8*b10*a4+4690*b4*b8*a6-12474*b4*b5*a10+5600*b8*a5*a9+714*b10*b5*a8-5992*b4*a5*a9+168*b10*b5*a4-826*b9*a10*a3-12026*b4*b9*a6+147*b8*b6*a4-5859*b8*b6*a8+126*b9*b10*a5+4746*b4*b6*a8-3332*b4*a8*a4+3682*b8*b10*a5+21462*b4*b5*a6+84*b4^2*a6+56*b10^2*b4-3024*a5*b9*a9-1288*b4*a6*a7+6132*a5^2*b4+2604*b7*a6*a9-560*b9*b10*a8-3843*b9*b6*a5+1974*b9*b6*a4+1512*b4*a5*a8+1134*b8*a5*a4-672*b8^2*a10-966*b9*a8*a4-1554*b9*a5*a4-1428*b8*a9^2-84*b4*b10*a4-4536*b8*a5^2-3654*b8*a4^2+1029*b4*b6*a5+532*b4*a8^2-12558*b9*b6^2+406*b8*b10^2+3066*b4*b6^2+4270*b5*b10^2-3430*b10^2*b9+13307*b9*b6*b10-1001*b8*b10*b6-1288*b4*b6*b10+12054*b6^2*b5+21*b6^2*b8-14903*b5*b6*b10-3192*b5*a6*a3+1204*b5*a6*a7-3388*b5*b9*a6-1834*b5*a8*a4-406*b5*a10*a3+4354*b5*b9*a10-1638*b5*a5*a4-1666*b5*b8*a6-3136*a5*b5*a9-42*b9*a8^2+4844*a5*b9*a8-4172*a5*b8*a8-1344*b8*a10*a7+3682*b10*b4*a9+2016*b9*a5^2-2373*b5*b6*a4+3290*b9^2*a6-2352*b9^2*a10+1344*b5*b6*a8+672*b8*b10*a8+2016*a5^2*b5-1099*b5*b6*a9"),
		SMQP("777*a6*b5*b6+420*a6*b5*a5-8337*a6*b4*b6-651*a6*b4*a5-714*a6*b5*a4-1638*a6*b4*a4-3780*b6^3+21*b10*b5*a6-420*b8*a10*a5-875*b10*b5*a10+917*b10*b9*a10-175*a9*b5*a10-882*a6*b9*a4-357*a6*b7*a10-126*a6*b9*a9-798*a6*b9*a8-1701*b10*b9*a6+1288*b5*b6*a10+1701*a6*b4*a9+294*a6*b8*a8+7*b9*a10*a8+196*a8*b4*a10-70*a8*b5*a10-854*a9*b4*a10+455*a10*b5*a4+364*a10*b9*a4-91*a10*b8*a8+126*a10^2*b7+70*a10*b9*a9+1932*b8*a6*a4-462*a8*b4*a6-1267*a4*b8*a10+1162*a4*b4*a10+35*b9*a10*a5+644*b8*a10*a9-735*b8*a6*a9-126*b8*b10*a6+252*a5*b9*a6+5019*b9*b6*a6+518*b4*a10*a5+147*a6*b5*a9+588*a6*b5*a8+1302*a6^2*b7+3885*b4*b10*a6+63*a5*b5*a10-1512*b8*b6*a6+595*b8*b10*a10+126*b8*a6*a5+5663*b4*b6*a10-371*b8*b6*a10-2989*b4*b10*a10+756*b10^3-4158*b10^2*b6+7182*b6^2*b10-2254*b9*b6*a10"),
		SMQP("280*a5*b5*a8-3318*b4*a5*a4+4123*b8*b4*a10+4809*b9*b6*a9+539*b10*b4*a5+1134*b7*a10*a5-1841*b10*b9*a9+1092*b8*a6*a3+1358*b8*a9*a4+378*b7*a10*a4-4095*b7*a6*a4-1512*b7*a10*a9-7140*b8*b6*a5-742*b5*a4*a9+798*b6*b5*a5-3087*b7*a6*a5-525*b8*a8*a9+679*b8*a8^2+1512*a6*b7*a8+4123*b8*a8*a4-371*b5*a8^2+2429*b8*a10*a3+2205*b5*a4^2+4410*b4*a4^2+1876*b4*a9^2+105*a8*b9*a9+210*a8*b4*a9+1512*b4*a10*a7+357*a8*b5*a9-1134*b9*a10*a7+770*a9^2*b5-70*a9*b9*a4-70*a9^2*b9+2737*b6*b9*a8-3920*b4*a10*a3-1841*b10*b9*a4-4774*b4^2*a10-3024*b4*b6*a4+4914*b7*b6*a10-12789*b7*b6*a6-2268*b7*b10*a10+2331*b8*b6*a9-1295*b8*b10*a9-14238*b4*b6*a9-1148*a4*b4*a9+973*b8*b5*a10-3402*b8*b9*a10+6489*b7*b10*a6+693*b9*a4^2-182*b10*b5*a5-819*b5^2*a6-1029*b4*a6*a3+1568*b8^2*a6+805*b8*a6*a7-1232*b4*b10*a8-518*b10*b5*a9+1099*b9*a6*a7+1645*b4*b9*a10-1995*b9*a6*a3+3094*b9*b8*a6-1463*b5^2*a10+2044*b8*b10*a4+1015*b4*b8*a6-9639*b4*b5*a10+5516*b8*a5*a9+2177*b10*b5*a8-10598*b4*a5*a9+994*b10*b5*a4+1253*b9*a10*a3-11711*b4*b9*a6-252*b8*b6*a4-4949*b8*b6*a8+2275*b9*b10*a5+8330*b4*b6*a8-4774*b4*a8*a4+2065*b8*b10*a5+25935*b4*b5*a6+1722*b4^2*a6+924*b10^2*b4-196*a5*b9*a9-700*b4*a6*a7+10206*a5^2*b4+1512*b7*a6*a9-1330*b9*b10*a8-5418*b9*b6*a5+4116*b9*b6*a4+5516*b4*a5*a8+2625*b8*a5*a4+756*b8^2*a10-623*b9*a8*a4+357*b9*a5*a4-1666*b8*a9^2-1498*b4*b10*a4-7434*b8*a5^2-3969*b8*a4^2+3444*b4*b6*a5-574*b4*a8^2-13776*b9*b6^2-1365*b8*b10^2+4368*b4*b6^2+5439*b5*b10^2-3171*b10^2*b9+13230*b9*b6*b10+5355*b8*b10*b6-2394*b4*b6*b10+18312*b6^2*b5-6888*b6^2*b8-20034*b5*b6*b10-3129*b5*a6*a3+637*b5*a6*a7-8407*b5*b9*a6-2219*b5*a8*a4-3451*b5*a10*a3+1134*b5*a10*a7+4795*b5*b9*a10-1029*b5*a5*a4+224*b5*b8*a6-3010*a5*b5*a9-35*b9*a8^2+2506*a5*b9*a8-3458*a5*b8*a8-378*b8*a10*a7+6377*b10*b4*a9+1134*b9*a5^2-3276*b5*b6*a4+518*b9^2*a6-378*b9^2*a10-3227*b5*b6*a8+686*b8*b10*a8+2016*a5^2*b5+1197*b5*b6*a9"),
		SMQP("504*a7^2*b7-336*a7*b7*b4-1848*a3*a7*b7+651*a3*b7*b4-1092*a3*b7*b8+126*b7*b4*b5-252*a3^2*b7-301*b7^2*a5-42*b7*b4*b9-21*a3*b7*b9+840*b7*b5^2+21*b7*b8*b5+525*a7*b7*b9-210*b7*b8^2+84*b7^2*b10+420*b7^2*a8-147*b7^2*a4-189*b7*b9^2-168*b7*b8*b4-1197*b7*b5*b9+819*b7*b8*b9+84*b7*b4^2+3339*a3*b7*b5-861*a7*b7*b5+63*a7*b7*b8-273*b7^2*b6"),
		SMQP("112*a5*b5*a8-987*b4*a5*a4+1414*b8*b4*a10+1575*b9*b6*a9-119*b10*b4*a5+63*b7*a10*a5-651*b10*b9*a9+462*b8*a6*a3+476*b8*a9*a4+441*b7*a10*a4-1260*b7*a6*a4-504*b7*a10*a9-1806*b8*b6*a5-238*b5*a4*a9-42*b6*b5*a5-1134*b7*a6*a5-210*b8*a8*a9+238*b8*a8^2+504*a6*b7*a8+1414*b8*a8*a4-119*b5*a8^2+728*b8*a10*a3+693*b5*a4^2+1449*b4*a4^2+602*b4*a9^2+63*a8*b9*a9+105*a8*b4*a9+504*b4*a10*a7+105*a8*b5*a9-63*b9*a10*a7+266*a9^2*b5-42*a9*b9*a4-42*a9^2*b9+987*b6*b9*a8-910*b4*a10*a3-651*b10*b9*a4-1547*b4^2*a10-840*b4*b6*a4+1323*b7*b6*a10-4284*b7*b6*a6-756*b7*b10*a10+798*b8*b6*a9-392*b8*b10*a9-4683*b4*b6*a9-406*a4*b4*a9+721*b8*b5*a10-189*b8*b9*a10+2205*b7*b10*a6+315*b9*a4^2-266*b10*b5*a5-315*b5^2*a6-168*b4*a6*a3+308*b8^2*a6+112*b8*a6*a7-427*b4*b10*a8-182*b10*b5*a9+546*b9*a6*a7+1113*b4*b9*a10-882*b9*a6*a3+1400*b9*b8*a6-581*b5^2*a10+679*b8*b10*a4+224*b4*b8*a6-4522*b4*b5*a10+2072*b8*a5*a9+749*b10*b5*a8-3514*b4*a5*a9+322*b10*b5*a4+84*b9*a10*a3-3598*b4*b9*a6-84*b8*b6*a4-1694*b8*b6*a8+567*b9*b10*a5+2653*b4*b6*a8-1547*b4*a8*a4+1099*b8*b10*a5+8106*b4*b5*a6+609*b4^2*a6+252*b10^2*b4-252*a5*b9*a9-266*b4*a6*a7+2982*a5^2*b4+504*b7*a6*a9-420*b9*b10*a8-2058*b9*b6*a5+1260*b9*b6*a4+1876*b4*a5*a8+756*b8*a5*a4-378*b8^2*a10-273*b9*a8*a4+231*b9*a5*a4-532*b8*a9^2-455*b4*b10*a4-2142*b8*a5^2-1386*b8*a4^2+1134*b4*b6*a5-203*b4*a8^2-4536*b9*b6^2-441*b8*b10^2+1512*b4*b6^2+1827*b5*b10^2-1071*b10^2*b9+4410*b9*b6*b10+1827*b8*b10*b6-882*b4*b6*b10+6048*b6^2*b5-2268*b6^2*b8-6678*b5*b6*b10-672*b5*a6*a3+70*b5*a6*a7-2800*b5*b9*a6-707*b5*a8*a4-112*b5*a10*a3+63*b5*a10*a7+1302*b5*b9*a10-441*b5*a5*a4+98*b5*b8*a6-910*a5*b5*a9-21*b9*a8^2+966*a5*b9*a8-1358*a5*b8*a8-441*b8*a10*a7+2128*b10*b4*a9+210*b9*a5^2-1092*b5*b6*a4+21*b9^2*a6-441*b9^2*a10-1043*b5*b6*a8+203*b8*b10*a8+504*a5^2*b5+357*b5*b6*a9"),
		SMQP("-7*b9*a5*a3+14*b7*a6*a3-28*b5*b4*a5+7*b5*b8*a5-7*b7*a5^2-14*b7*b5*a6-7*b5*b9*a5+21*b9*b8*a5+7*a5*b7*a4+21*b5*a5*a3-7*a7*b8*a5+7*a3*b4*a5-14*b8^2*a5-7*b7*b6*a5+14*b9*b4*a5-7*b9^2*a5+7*b9*a5*a7-7*a7*b5*a5"),
		SMQP("28*a5*b5*a8-7350*b4*a5*a4+3472*b8*b4*a10-4746*b9*b6*a9+2674*b10*b4*a5-2394*b7*a10*a5+1512*b10*b9*a9+210*b8*a6*a3-28*b8*a9*a4-2394*b7*a10*a4+8316*b7*a6*a4+1512*b7*a10*a9+6825*b8*b6*a5-364*b5*a4*a9-2562*b6*b5*a5+1764*b7*a6*a5+504*b8*a8*a9-308*b8*a8^2+1008*a6*b7*a8-1652*b8*a8*a4+1162*b5*a8^2-952*b8*a10*a3-1008*b5*a4^2+2772*b4*a4^2-2212*b4*a9^2-84*a8*b9*a9+2478*a8*b4*a9+1008*b4*a10*a7-1890*a8*b5*a9+1386*b9*a10*a7+252*a8*b7*a10-112*a9^2*b5-840*a9*b9*a4+1176*a9^2*b9-987*b6*b9*a8-1918*b4*a10*a3-546*b10*b9*a4+9226*b4^2*a10-840*b4*b6*a4-2898*b7*b6*a10+8379*b7*b6*a6+1764*b7*b10*a10-1638*b8*b6*a9+700*b8*b10*a9+5313*b4*b6*a9+3962*a4*b4*a9+826*b8*b5*a10+3402*b8*b9*a10-4536*b7*b10*a6-3024*b9*a4^2+3304*b10*b5*a5+3150*b5^2*a6-2436*b4*a6*a3-7168*b8^2*a6-854*b8*a6*a7-2926*b4*b10*a8-1484*b10*b5*a9-1848*b9*a6*a7-12054*b4*b9*a10+2016*b9*a6*a3+11186*b9*b8*a6+2422*b5^2*a10-2156*b8*b10*a4-11368*b4*b8*a6+11018*b4*b5*a10-616*b8*a5*a9-1582*b10*b5*a8+56*b4*a5*a9+3346*b10*b5*a4-546*b9*a10*a3+23492*b4*b9*a6+4893*b8*b6*a4+1015*b8*b6*a8+1470*b9*b10*a5+1498*b4*b6*a8-392*b4*a8*a4-6104*b8*b10*a5-13944*b4*b5*a6-10794*b4^2*a6-3192*b10^2*b4+1428*a5*b9*a9+3136*b4*a6*a7-1428*a5^2*b4-4284*b7*a6*a9+546*b9*b10*a8+1785*b9*b6*a5+2898*b9*b6*a4+5446*b4*a5*a8+3402*b8*a5*a4-1764*b8^2*a10-1260*b9*a8*a4+2184*b9*a5*a4+224*b8*a9^2-2870*b4*b10*a4+1260*b8*a5^2-252*b8*a4^2-6111*b4*b6*a5-1820*b4*a8^2+11718*b9*b6^2+1050*b8*b10^2-7182*b4*b6^2-5250*b5*b10^2+3486*b10^2*b9-13335*b9*b6*b10-3171*b8*b10*b6+9660*b4*b6*b10-15498*b6^2*b5+3591*b6^2*b8+18123*b5*b6*b10+1344*b5*a6*a3-140*b5*a6*a7+3584*b5*b9*a6+3976*b5*a8*a4+4760*b5*a10*a3-1890*b5*a10*a7-8022*b5*b9*a10+1638*b5*a5*a4-3346*b5*b8*a6+1316*a5*b5*a9+504*b9*a8^2-5040*a5*b9*a8-686*a5*b8*a8-630*b8*a10*a7-308*b10*b4*a9-924*b9*a5^2-3675*b5*b6*a4-8106*b9^2*a6+2142*b9^2*a10+406*b5*b6*a8+770*b8*b10*a8-1008*a5^2*b5+2709*b5*b6*a9"),
		SMQP("2527*b7*b4*a5+483*b7*b9*a5-1953*b7*b4*a9-4305*b7*b8*a5+1428*a7*b9^2+588*a7*b8^2+1281*b7^2*a6+1806*b7*b8*a9-14910*b7*b5*a9-504*b7^2*a10-1533*a7*b7*a8-4200*a7*b7*a5-1911*a7*b7*a4-3612*b7*b4*a8+609*b7*b5*a8+4893*b7*a9*a7+12348*b7*b5*a4-2856*b7*a8*a3+1113*b7*a4*a3+1323*b7*a5*a3+651*b7*a9*a3-5271*b7*b9*a4+840*b7*b8*a8+2331*b7*b9*a8-2604*a7*a3*b8+1512*a3*b8^2+7686*b7*b5*b6+8736*b5^2*b8+3927*b7*b9*a9-2058*b5^3+84*a3*b7*b6-4032*b9*a3*b8-1680*a7*b4*b9+2520*a3^2*b8+5964*a3*b5*a7-3150*b8*b7*b6+1470*b7*b8*b10-1218*b8*b5*a7-6552*a3^2*b5-3276*a7*b7*b6-336*b4*b8^2-4536*a3*a7*b9-252*a3*b9^2+3108*a3^2*b9-2016*a7*b8*b9+3822*a7*b5*b9+840*a7*b4*b5+1680*a7*b7*b10-6468*a3^2*b4-5670*a7*b5^2+1680*b4^2*b9+8778*b8*b5*b9-6216*b4*b9^2-20958*b4*b5^2+8988*a7*b4*a3+11004*a3*b5^2+1764*b8^2*b9-1092*b4^2*a3-8946*b9^2*b5+168*b8^3+672*b4^2*a7-4494*a3*b8*b5-5418*b7*b5*b10+20076*b5*b9*b4+2352*a3*b4*b9-756*b9^2*b8-8232*b8^2*b5+504*b4*a7*b8+1344*a3*b4*b5+3150*a3*b5*b9+3024*b7*b4*b6-420*b4^2*b5+1260*a7^2*b8-924*a7^2*b5+420*a7^2*b9-2352*a7^2*b4-1386*b7*b8*a4+3024*b7*b4*a4-1512*b7*b4*b10-2142*b5^2*b9+2184*b9^3+1092*b9*b7*b6+210*b9*b7*b10+1764*b4*b8*b5+336*b4*a3*b8-1008*b4*b8*b9+5131*b7*b5*a5"),
		SMQP("14385*b7*b4*a5+4305*b7*b9*a5-8043*b7*b4*a9-27363*b7*b8*a5+9828*a7*b9^2+2772*a7*b8^2+16443*b7^2*a6+14154*b7*b8*a9-119868*b7*b5*a9-7056*b7^2*a10-20097*a7*b7*a8-26208*a7*b7*a5-1071*a7*b7*a4-22344*b7*b4*a8+21105*b7*b5*a8+37233*b7*a9*a7+68628*b7*b5*a4-11592*b7*a8*a3+9261*b7*a4*a3+315*b7*a5*a3+8127*b7*a9*a3-43827*b7*b9*a4-6510*b7*b8*a8+17661*b7*b9*a8-13356*a7*a3*b8+4536*a3*b8^2+52668*b7*b5*b6+63504*b5^2*b8+32151*b7*b9*a9-17892*b5^3-12180*a3*b7*b6-4032*b9*a3*b8-4284*a7*b4*b9+7560*a3^2*b8+36540*a3*b5*a7-18480*b8*b7*b6+6636*b7*b8*b10-8568*b8*b5*a7-26964*a3^2*b5-24192*a7*b7*b6-1512*b4*b8^2-32508*a3*a7*b9-13608*a3*b9^2+19404*a3^2*b9-17640*a7*b8*b9+34272*a7*b5*b9-2268*a7*b4*b5+11088*a7*b7*b10-36036*a3^2*b4-44100*a7*b5^2+5544*b4^2*b9+70308*b8*b5*b9-53676*b4*b9^2+4032*a3*b7*b10-146412*b4*b5^2+64008*a7*b4*a3+92736*a3*b5^2+14364*b8^2*b9-2772*b4^2*a3-66024*b9^2*b5-2520*b8^3+2016*b4^2*a7-36540*a3*b8*b5-35700*b7*b5*b10+142884*b5*b9*b4+22680*a3*b4*b9-9072*b9^2*b8-56700*b8^2*b5+5292*b4*a7*b8-1260*a3*b4*b5+11340*a3*b5*b9+15708*b7*b4*b6-2520*b4^2*b5+6048*a7^2*b8-7056*a7^2*b5+5040*a7^2*b9-17136*a7^2*b4+4704*b7*b8*a4+14196*b7*b4*a4-7392*b7*b4*b10-20412*b5^2*b9+19404*b9^3+14196*b9*b7*b6-588*b9*b7*b10+1260*b4*b8*b5-5292*b4*a3*b8+8820*b4*b8*b9+33397*b7*b5*a5"),
		SMQP("-1197*a6*b5*b6-756*a6*b5*a5+8757*a6*b4*b6+819*a6*b4*a5+378*a6*b5*a4+882*a6*b4*a4+3780*b6^3-1617*b10*b5*a6-357*b8*a10*a5+728*b10*b5*a10-770*b10*b9*a10+322*a9*b5*a10+1218*a6*b9*a4+84*a6*b7*a10+546*a6*b9*a9+462*a6*b9*a8+777*b10*b9*a6-574*b5*b6*a10-2457*a6*b4*a9-126*a6*b8*a8-154*b9*a10*a8-70*a8*b4*a10-224*a8*b5*a10+371*a9*b4*a10-728*a10*b5*a4-154*a10*b9*a4+238*a10*b8*a8-126*a10^2*b7+224*a10*b9*a9-1008*b8*a6*a4+882*a8*b4*a6+1834*a4*b8*a10-1582*a4*b4*a10+238*b9*a10*a5-938*b8*a10*a9+483*b8*a6*a9+2310*b8*b10*a6-84*a5*b9*a6-6615*b9*b6*a6+343*b4*a10*a5-63*a6*b5*a9-756*a6*b5*a8-1890*a6^2*b7-1281*b4*b10*a6+231*a5*b5*a10+3276*b8*b6*a6-742*b8*b10*a10-630*b8*a6*a5-7553*b4*b6*a10-196*b8*b6*a10+2653*b4*b10*a10-756*b10^3+4158*b10^2*b6-7182*b6^2*b10+2044*b9*b6*a10"),
		SMQP("483*b7*b4*a5+567*b7*b9*a5-5061*b7*b4*a9-5313*b7*b8*a5+756*a7*b9^2-504*a7*b8^2-315*b7^2*a6+1974*b7*b8*a9-15498*b7*b5*a9+504*b7^2*a10-189*a7*b7*a8-5040*a7*b7*a5-4095*a7*b7*a4-3108*b7*b4*a8+1113*b7*b5*a8+5229*b7*a9*a7+15120*b7*b5*a4-6552*b7*a8*a3+3465*b7*a4*a3+2331*b7*a5*a3+315*b7*a9*a3-5943*b7*b9*a4+1596*b7*b8*a8+3255*b7*b9*a8-3276*a7*a3*b8+504*a3*b8^2+9534*b7*b5*b6+9072*b5^2*b8+3759*b7*b9*a9-2394*b5^3+420*a3*b7*b6-5796*b9*a3*b8-1008*a7*b4*b9+3528*a3^2*b8+8820*a3*b5*a7-3654*b8*b7*b6+2058*b7*b8*b10-1134*b8*b5*a7-7560*a3^2*b5-3780*a7*b7*b6-3024*b4*b8^2-6048*a3*a7*b9+252*a3*b9^2+3780*a3^2*b9+756*a7*b8*b9+2142*a7*b5*b9-2016*a7*b4*b5+2016*a7*b7*b10-7308*a3^2*b4-4662*a7*b5^2+6048*b4^2*b9+6174*b8*b5*b9-8568*b4*b9^2-1008*a3*b7*b10-21042*b4*b5^2+9324*a7*b4*a3+8316*a3*b5^2+3024*b8^2*b9-1008*b4^2*a3-8442*b9^2*b5+1008*b4^2*a7-1890*a3*b8*b5-6510*b7*b5*b10+19908*b5*b9*b4+6048*a3*b4*b9-1008*b9^2*b8-7812*b8^2*b5-504*b4*a7*b8+2268*a3*b4*b5+4662*a3*b5*b9+2520*b7*b4*b6-6552*b4^2*b5+756*a7^2*b8-1764*a7^2*b5+1260*a7^2*b9-2016*a7^2*b4-2562*b7*b8*a4+5712*b7*b4*a4-2184*b7*b4*b10-1134*b5^2*b9+2016*b9^3-252*b9*b7*b6+966*b9*b7*b10+5040*b4*b8*b5-504*b4*a3*b8+1512*b4*b8*b9+4627*b7*b5*a5"),
		SMQP("-567*b9^2*a9+1071*b9^2*a4-945*b9^2*a8-392*b7*b4*a10-6104*b7*b10*a9+1400*b7*a10*a3-1596*b7*b4*a6+5992*b7*b10*a4+5796*b4*a4*a3+2520*b8*a8*a3+1071*b8*b9*a4-2457*b8*a4*a3-2457*b8*a9*a7+3360*b7*a6*a3-6048*b4*a8*a3-1764*b4*a8*a7+2016*b5^2*a5+3024*b4*b5*a4+5488*a7*b7*a6+1512*a7*b5*a9-6552*a7*b5*a5-2016*b8*b7*a10-6048*a7*b4*a5-2016*a7*b7*a10+5796*a7*b4*a9+1008*a7*b8*a5+945*a7*b8*a8+1953*b8*a9*a3-4599*b8*a5*a3+567*b8*b4*a9-1197*b8*b4*a5-1260*b8*b5*a4+504*b4*b8*a8+1764*b4*b9*a8-2268*b4*a4*a7+2772*b4*b5*a8+1071*a7*b8*a4+1512*b7*a4^2-1400*b7*a4*a8-4725*b8^2*a5+8820*a3*b4*a5-4284*a3*b4*a9-6664*b7*a5*a9-8001*b9^2*a5+7840*b7*a5*a8+5040*b5*a4*a3+6664*b7*b5*a10-5040*b5*a4*a7-2520*b5^2*a4+4032*b5^2*a8-2520*b5^2*a9-2268*b4*b9*a4+567*b5*b9*a5-1449*b5*b8*a5-10836*b5*b4*a5-2268*b4^2*a9+63*b5*b8*a8+5201*b7*b8*a6-252*b4^2*a5+168*b7*a9*a8-840*b7*b6*a5-17136*b7*b5*a6-616*b7*a9*a4-13440*b7*b6*a4+8736*b7*b6*a9+56*b7*b10*a8+63*b5*b9*a8+3304*b7*b6*a8-112*b7*a9^2-2520*b7*a5^2+3472*b7*b10*a5+1512*b5*b9*a9-56*b7*a8^2-1344*a5*b7*a4+1134*b8^2*a9-126*b8^2*a8-2772*b9*b5*a4-2520*b5*a9*a3+3528*b5*b8*a9-6552*b5*b4*a9-6048*b9^2*b10+26208*b4*b5*b6+12096*b8^2*b6+3024*a3*b8*b6-10080*a3*b9*b6+9072*a3*b10*b5-9072*a3*b5*b6+8064*b9^2*b6+2016*b7*b6^2+6048*b10*b9*b8-5544*b5^2*b10+8064*b5^2*b6-6048*b8^2*b10-16632*b5*b6*b8+3024*b5*b6*b9+5544*b10*b9*b5-1008*b7*b6*b10+5040*a3*b4*b6-10080*b8*b9*b6+12096*b4*b10*b9-20160*b4*b10*b5-22176*b4*b9*b6-2016*a7*b4*b6-3024*a7*b10*b5+2016*a7*b8*b6+4032*a7*b5*b6+6552*b5*b8*b10+2457*b9*a9*a3-5103*b9*a5*a3-2457*b9*b8*a9-1071*b9*b8*a8+4851*b9*b4*a5+14490*b9*b8*a5+2331*b9*b4*a9-8064*b5*a8*a3+3024*b5*a8*a7+14616*b5*a5*a3+6048*b9*a5*a7+2520*b9*a8*a3-2457*b9*a4*a3-3465*b9*a9*a7+1071*b9*a4*a7+945*b9*a8*a7-1351*b9*b7*a6"),
		SMQP("-42*a7*b7*b4+168*b7^2*a9+714*a3*a7*b7+420*a3*b7*b4-21*a3*b7*b8+525*b7*b4*b5+126*a3^2*b7+553*b7^2*a5-462*b7*b4*b9-399*a3*b7*b9-903*b7*b5^2-378*b7*b8*b5-399*a7*b7*b9+462*b7*b8^2-168*b7^2*a8-273*b7^2*a4+231*b7*b9^2-126*b7*b8*b4+1113*b7*b5*b9-693*b7*b8*b9-168*a3*b7*b5-399*a7*b7*b5+399*a7*b7*b8+273*b7^2*b6"),
		SMQP("525*b9^2*a9+8967*b9^2*a4+147*b9^2*a8+336*b4^2*b10-12502*b7*b4*a10+266*b7*b10*a9-1862*b7*a10*a3+17535*b7*b4*a6-4774*b7*b10*a4-10311*b4*a4*a3-7266*b8*a8*a3-17367*b8*b9*a4+8736*b8*a4*a3-273*b8*a9*a7+1806*b7*a6*a3+3549*b4*a8*a3-1680*b4*a8*a7+1890*b5^2*a5+14931*b4*b5*a4-3766*a7*b7*a6-1155*a7*b5*a9+2520*a7*b5*a5+4032*b8*b7*a10+3066*a7*b4*a5+1512*a7*b7*a10+210*a7*b4*a9-1638*a7*b8*a5+1365*a7*b8*a8+168*b8*a9*a3+5922*b8*a5*a3-3717*b8*b4*a9+5250*b8*b4*a5-273*b8*b5*a4-693*b4*b8*a8-2079*b4*b9*a8+3864*b4*a4*a7+2016*b4*b5*a8+1491*a7*b8*a4+2709*b7*a4^2-10045*b7*a4*a8-3213*b8^2*a5-6510*a3*b4*a5+2163*a3*b4*a9-5705*b7*a5*a9+2184*b9^2*a5-2905*b7*a5*a8-19803*b5*a4*a3-7210*b7*b5*a10+9303*b5*a4*a7-1239*b5^2*a4+2982*b5^2*a8-105*b5^2*a9-7560*b4*b9*a4+1155*b5*b9*a5+4221*b5*b8*a5+399*b5*b4*a5+10059*b4^2*a9-714*b5*b8*a8-497*b7*b8*a6-9429*b4^2*a5+1722*b7*a9*a8+2982*b7*b6*a5+18837*b7*b5*a6-539*b7*a9*a4+11235*b7*b6*a4-5271*b7*b6*a9+3514*b7*b10*a8-5838*b5*b9*a8-1015*b7*b6*a8+1540*b7*a9^2+6300*b7*a5^2-5152*b7*b10*a5-3129*b5*b9*a9+1274*b7*a8^2+8967*a5*b7*a4-1176*b8^2*a9+273*b8^2*a8+3297*b9*b5*a4+2373*b5*a9*a3-210*b5*b8*a9+4158*b5*b4*a9-714*b9^2*b10-2562*b4^2*b6-21798*b4*b5*b6-4431*b8^2*b6+5880*a3*b8*b6+8358*a3*b9*b6-5838*a3*b8*b10+4998*a3*b10*b5+42*a3*b10*b9-9240*a3*b5*b6-4494*b9^2*b6-336*b4^2*a4+1428*a3*b10*b4+1302*b4^2*a8+3402*b7*b6^2+2268*b7*b10^2+7665*b10*b9*b8+9429*b5^2*b10-13188*b5^2*b6+1365*b8^2*b10+15771*b5*b6*b8+14280*b5*b6*b9-9471*b10*b9*b5-6804*b7*b6*b10-6132*a3*b4*b6-10290*b8*b9*b6-3339*b4*b10*b9+8127*b4*b10*b5-6489*b4*b8*b10+13734*b4*b9*b6+924*a7*b4*b6-1722*a7*b9*b6+1806*a7*b10*b9-2562*a7*b10*b5+462*a7*b8*b10-1932*a7*b8*b6-588*a7*b10*b4+2478*a7*b5*b6+6069*b8^2*a4+14805*b4*b8*b6+1197*b4*b8*a4-7098*b5*b8*b10-1995*b9*a9*a3+1680*b9*a5*a3+2541*b9*b8*a9+5250*b9*b8*a8+5859*b9*b4*a5-6279*b9*b8*a5-7056*b9*b4*a9+9555*b5*a8*a3-3549*b5*a8*a7-5544*b5*a5*a3-1722*b9*a5*a7-3255*b9*a8*a3+4935*b9*a4*a3+1659*b9*a9*a7-7287*b9*a4*a7+1785*b9*a8*a7+3528*b9*b7*a10-13286*b9*b7*a6"),
		SMQP("609*b9^2*a9-21*b9^2*a4+609*b9^2*a8+504*b4^2*b10-5390*b7*b4*a10-14*b7*b10*a9+98*b7*a10*a3+11067*b7*b4*a6-1526*b7*b10*a4-441*b4*a4*a3-378*b8*a8*a3+231*b8*b9*a4+504*b8*a4*a3-609*b8*a9*a7+3990*b7*a6*a3-315*b4*a8*a3+546*b4*a8*a7-882*b5^2*a5-6027*b4*b5*a4-3290*a7*b7*a6-735*a7*b5*a9-252*a7*b5*a5+1008*b8*b7*a10+1470*a7*b4*a5+504*a7*b7*a10-84*a7*b4*a9-630*a7*b8*a5-735*a7*b8*a8+1890*b8*a5*a3+1029*b8*b4*a9+420*b8*b4*a5+3675*b8*b5*a4+525*b4*b8*a8+105*b4*b9*a8+798*b4*a4*a7+966*b4*b5*a8-1113*a7*b8*a4+3213*b7*a4^2-1421*b7*a4*a8+441*b8^2*a5-4032*a3*b4*a5+1701*a3*b4*a9+791*b7*a5*a9+210*b9^2*a5-6545*b7*a5*a8+2205*b5*a4*a3-3122*b7*b5*a10-357*b5*a4*a7-1785*b5^2*a4+2562*b5^2*a8-147*b5^2*a9+2058*b4*b9*a4-2667*b5*b9*a5+5859*b5*b8*a5-2037*b5*b4*a5+7119*b4^2*a9-672*b5*b8*a8+1253*b7*b8*a6-945*b4^2*a5+546*b7*a9*a8+546*b7*b6*a5+4347*b7*b5*a6-3115*b7*a9*a4+987*b7*b6*a4-399*b7*b6*a9+770*b7*b10*a8-3864*b5*b9*a8-2135*b7*b6*a8+980*b7*a9^2+616*b7*b10*a5-1659*b5*b9*a9+994*b7*a8^2+4011*a5*b7*a4-1092*b8^2*a9-777*b8^2*a8+231*b9*b5*a4+1701*b5*a9*a3+840*b5*b8*a9-798*b5*b4*a9+840*b9^2*b10-3654*b4^2*b6-7518*b4*b5*b6-357*b8^2*b6+4032*a3*b8*b6+3150*a3*b9*b6-2646*a3*b8*b10+2646*a3*b10*b5-630*a3*b10*b9-5040*a3*b5*b6-3066*b9^2*b6-504*b4^2*a4-756*a3*b10*b4-3654*b4^2*a8+1890*b7*b6^2+1260*b7*b10^2+1029*b10*b9*b8+2415*b5^2*b10-6216*b5^2*b6+1407*b8^2*b10+5649*b5*b6*b8+7392*b5*b6*b9-3003*b10*b9*b5-3276*b7*b6*b10-252*a3*b4*b6-3696*b8*b9*b6-4641*b4*b10*b9+5397*b4*b10*b5-1155*b4*b8*b10+8526*b4*b9*b6+84*a7*b4*b6-1050*a7*b9*b6+210*a7*b10*b9-966*a7*b10*b5+1050*a7*b8*b10-1344*a7*b8*b6+84*a7*b10*b4+1806*a7*b5*b6-3045*b8^2*a4+2247*b4*b8*b6+1911*b4*b8*a4-2478*b5*b8*b10-2331*b9*a9*a3+252*b9*a5*a3+861*b9*b8*a9+3318*b9*b8*a8+4053*b9*b4*a5-2919*b9*b8*a5-5628*b9*b4*a9+4347*b5*a8*a3-1365*b5*a8*a7+756*b5*a5*a3-42*b9*a5*a7-2583*b9*a8*a3+567*b9*a4*a3+1491*b9*a9*a7-399*b9*a4*a7+1365*b9*a8*a7+2520*b9*b7*a10-7756*b9*b7*a6"),
		SMQP("-1533*b9^2*a9-1911*b9^2*a4-6069*b9^2*a8+7924*b7*b4*a10-13916*b7*b10*a9+3668*b7*a10*a3-11172*b7*b4*a6+16324*b7*b10*a4+20748*b4*a4*a3+11508*b8*a8*a3+7665*b8*b9*a4-10227*b8*a4*a3-7329*b8*a9*a7+6636*b7*a6*a3-24234*b4*a8*a3+3990*b4*a8*a7-4662*b5^2*a5+2814*b4*b5*a4+4732*a7*b7*a6+5208*a7*b5*a9-10836*a7*b5*a5-2016*b8*b7*a10-10164*a7*b4*a5-6048*a7*b7*a10+10038*a7*b4*a9+1008*a7*b8*a5+1617*a7*b8*a8+5523*b8*a9*a3-14049*b8*a5*a3+4809*b8*b4*a9-9219*b8*b4*a5-9450*b8*b5*a4+1974*b4*b8*a8+2184*b4*b9*a8-14406*b4*a4*a7+20454*b4*b5*a8+5523*a7*b8*a4-7686*b7*a4^2+5026*b7*a4*a8-2835*b8^2*a5+22470*a3*b4*a5-14280*a3*b4*a9+8330*b7*a5*a9-7035*b9^2*a5-1358*b7*a5*a8+14532*b5*a4*a3+23548*b7*b5*a10-13440*b5*a4*a7-11256*b5^2*a4+8778*b5^2*a8-2940*b5^2*a9-10542*b4*b9*a4-2247*b5*b9*a5+1575*b5*b8*a5-18564*b5*b4*a5-13608*b4^2*a9-11277*b5*b8*a8+7847*b7*b8*a6+4536*b4^2*a5-3612*b7*a9*a8-1092*b7*b6*a5-59850*b7*b5*a6+1022*b7*a9*a4-28686*b7*b6*a4+24990*b7*b6*a9-2716*b7*b10*a8+2079*b5*b9*a8+406*b7*b6*a8-3640*b7*a9^2-15372*b7*a5^2+11536*b7*b10*a5+8442*b5*b9*a9+4228*b7*a8^2-462*a5*b7*a4+2310*b8^2*a9-1722*b8^2*a8+7938*b9*b5*a4-9660*b5*a9*a3+13734*b5*b8*a9-31080*b5*b4*a9-14406*b9^2*b10+1260*b4^2*b6+63840*b4*b5*b6+21252*b8^2*b6+2856*a3*b8*b6-31584*a3*b9*b6+4116*a3*b8*b10+1596*a3*b10*b5+8484*a3*b10*b9-5460*a3*b5*b6+24780*b9^2*b6-5880*a3*b10*b4+252*b4^2*a8-504*b7*b6^2-1512*b7*b10^2+1722*b10*b9*b8-16212*b5^2*b10+30156*b5^2*b6-6972*b8^2*b10-30744*b5*b6*b8-22428*b5*b6*b9+28098*b10*b9*b5+4536*b7*b6*b10+16044*a3*b4*b6-13146*b8*b9*b6+28938*b4*b10*b9-30786*b4*b10*b5+2814*b4*b8*b10-60648*b4*b9*b6-7896*a7*b4*b6+1176*a7*b9*b6-3612*a7*b10*b9-2940*a7*b10*b5+2100*a7*b8*b10+2100*a7*b8*b6+4200*a7*b10*b4+7392*a7*b5*b6+1680*b8^2*a4-6762*b4*b8*b6-1302*b4*b8*a4+6552*b5*b8*b10+7959*b9*a9*a3-10731*b9*a5*a3-12369*b9*b8*a9+4263*b9*b8*a8+5313*b9*b4*a5+17934*b9*b8*a5+16863*b9*b4*a9-19110*b5*a8*a3+4704*b5*a8*a7+24570*b5*a5*a3+11760*b9*a5*a7+8526*b9*a8*a3-7035*b9*a4*a3-6279*b9*a9*a7+7329*b9*a4*a7-105*b9*a8*a7-7056*b9*b7*a10+8603*b9*b7*a6"),
		SMQP("896*a5*b5*a8+756*b4*a5*a4-1561*b8*b4*a10+2219*b9*b6*a9+210*b10*b4*a5+1344*b7*a10*a5-749*b10*b9*a9+168*b8*a6*a3+378*b8*a9*a4-588*b7*a10*a4-1911*b7*a6*a4-504*b7*a10*a9-2709*b8*b6*a5-378*b5*a4*a9+504*b6*b5*a5-987*b7*a6*a5-161*b8*a8*a9+259*b8*a8^2-336*a6*b7*a8+679*b8*a8*a4-581*b5*a8^2+721*b8*a10*a3+1071*b5*a4^2-252*b4*a4^2+980*b4*a9^2+119*a8*b9*a9-1120*a8*b4*a9+252*b4*a10*a7+763*a8*b5*a9-1092*b9*a10*a7+42*a9^2*b5+28*a9*b9*a4-602*a9^2*b9+434*b6*b9*a8-1106*b4*a10*a3-133*b10*b9*a4-70*b4^2*a10-1176*b4*b6*a4+2352*b7*b6*a10-4410*b7*b6*a6-756*b7*b10*a10+973*b8*b6*a9-665*b8*b10*a9-3913*b4*b6*a9-154*a4*b4*a9-973*b8*b5*a10-2268*b8*b9*a10+1596*b7*b10*a6+399*b9*a4^2-364*b10*b5*a5-315*b5^2*a6-273*b4*a6*a3+3808*b8^2*a6+497*b8*a6*a7+560*b4*b10*a8+931*b10*b5*a9+483*b9*a6*a7+1197*b4*b9*a10-483*b9*a6*a3-4424*b9*b8*a6-721*b5^2*a10+399*b8*b10*a4+2653*b4*b8*a6-1715*b4*b5*a10+1288*b8*a5*a9-63*b10*b5*a8-2142*b4*a5*a9+147*b10*b5*a4+819*b9*a10*a3-7175*b4*b9*a6-105*b8*b6*a4-1806*b8*b6*a8+1078*b9*b10*a5+1708*b4*b6*a8-630*b4*a8*a4-112*b8*b10*a5+10605*b4*b5*a6-672*b4^2*a6+448*b10^2*b4-532*a5*b9*a9-784*b4*a6*a7+2226*a5^2*b4+1428*b7*a6*a9-371*b9*b10*a8-1701*b9*b6*a5+1302*b9*b6*a4-98*b4*a5*a8+357*b8*a5*a4+1680*b8^2*a10+49*b9*a8*a4-777*b9*a5*a4-210*b8*a9^2-224*b4*b10*a4-1386*b8*a5^2-693*b8*a4^2+1617*b4*b6*a5+504*b4*a8^2-6174*b9*b6^2+224*b8*b10^2+882*b4*b6^2+2156*b5*b10^2-1736*b10^2*b9+6643*b9*b6*b10-931*b8*b10*b6+196*b4*b6*b10+6678*b6^2*b5-567*b6^2*b8-7819*b5*b6*b10-1785*b5*a6*a3+497*b5*a6*a7-2093*b5*b9*a6-1085*b5*a8*a4-3563*b5*a10*a3+1260*b5*a10*a7+3465*b5*b9*a10-1155*b5*a5*a4+28*b5*b8*a6-1358*a5*b5*a9+91*b9*a8^2+2072*a5*b9*a8-994*a5*b8*a8+588*b8*a10*a7+1008*b10*b4*a9+546*b9*a5^2-1113*b5*b6*a4+2772*b9^2*a6-420*b9^2*a10+105*b5*b6*a8+525*b8*b10*a8+1008*a5^2*b5-518*b5*b6*a9"),
		SMQP("581*b7*b4*a5+1897*b7*b9*a5-1715*b7*b4*a9-1071*b7*b8*a5-252*a7*b9^2-840*a7*b8^2+4179*b7^2*a6+3626*b7*b8*a9-26152*b7*b5*a9-2016*b7^2*a10-5229*a7*b7*a8-4200*a7*b7*a5+2709*a7*b7*a4-2744*b7*b4*a8+6293*b7*b5*a8+7245*b7*a9*a7+10388*b7*b5*a4-504*b7*a8*a3+1365*b7*a4*a3+2163*b7*a5*a3+4179*b7*a9*a3-9359*b7*b9*a4-3346*b7*b8*a8+2821*b7*b9*a8-168*a7*a3*b8+672*a3*b8^2+7336*b7*b5*b6+12264*b5^2*b8+6979*b7*b9*a9-3024*b5^3-3360*a3*b7*b6-84*b9*a3*b8+2352*a7*b4*b9+9744*a3*b5*a7-3584*b8*b7*b6+952*b7*b8*b10-1680*b8*b5*a7-5040*a3^2*b5-5880*a7*b7*b6-9072*a3*a7*b9-2436*a3*b9^2+4536*a3^2*b9-924*a7*b8*b9+10332*a7*b5*b9-5208*a7*b4*b5+2016*a7*b7*b10-6048*a3^2*b4-10752*a7*b5^2+17640*b8*b5*b9-13272*b4*b9^2+1344*a3*b7*b10-31248*b4*b5^2+13776*a7*b4*a3+22680*a3*b5^2+5712*b8^2*b9-13020*b9^2*b5-1008*b8^3-9240*a3*b8*b5-5432*b7*b5*b10+32508*b5*b9*b4+4032*a3*b4*b9-6804*b9^2*b8-11256*b8^2*b5+672*b4*a7*b8+2520*a3*b4*b5-2856*a3*b5*b9+3584*b7*b4*b6-168*a7^2*b8-2184*a7^2*b5+2184*a7^2*b9-3696*a7^2*b4+1904*b7*b8*a4+2128*b7*b4*a4-1624*b7*b4*b10-6972*b5^2*b9+6132*b9^3+5852*b9*b7*b6-952*b9*b7*b10-2352*b4*b8*b5-1176*b4*a3*b8+3024*b4*b8*b9-1715*b7*b5*a5"),
		SMQP("749*b9^2*a9+287*b9^2*a4-1729*b9^2*a8+406*b7*b4*a10-4970*b7*b10*a9-1162*b7*a10*a3-1659*b7*b4*a6+5446*b7*b10*a4+5159*b4*a4*a3+490*b8*a8*a3+1827*b8*b9*a4-1148*b8*a4*a3-623*b8*a9*a7+1050*b7*a6*a3-5425*b4*a8*a3+308*b4*a8*a7+1848*b5^2*a5+413*b4*b5*a4+2338*a7*b7*a6+973*a7*b5*a9-3360*a7*b5*a5-2870*a7*b4*a5-1512*a7*b7*a10+2870*a7*b4*a9+210*a7*b8*a5+1015*a7*b8*a8+1540*b8*a9*a3-4830*b8*a5*a3-1267*b8*b4*a9-1862*b8*b4*a5-1855*b8*b5*a4+1085*b4*b8*a8-133*b4*b9*a8-2884*b4*a4*a7+4340*b4*b5*a8-119*a7*b8*a4+63*b7*a4^2-791*b7*a4*a8-4011*b8^2*a5+6034*a3*b4*a5-4291*a3*b4*a9-2275*b7*a5*a9-4634*b9^2*a5+3661*b7*a5*a8+4417*b5*a4*a3+8218*b7*b5*a10-5285*b5*a4*a7-3353*b5^2*a4+4228*b5^2*a8-1127*b5^2*a9-2380*b4*b9*a4+3745*b5*b9*a5-2415*b5*b8*a5-10577*b5*b4*a5-2751*b4^2*a9-658*b5*b8*a8+833*b7*b8*a6+609*b4^2*a5-378*b7*a9*a8+546*b7*b6*a5-12495*b7*b5*a6-889*b7*a9*a4-8463*b7*b6*a4+6699*b7*b6*a9+182*b7*b10*a8-1722*b5*b9*a8+1771*b7*b6*a8-196*b7*a9^2-3612*b7*a5^2-1232*b7*b10*a5+1575*b5*b9*a9+742*b7*a8^2-651*a5*b7*a4+2464*b8^2*a9-77*b8^2*a8+273*b9*b5*a4-4403*b5*a9*a3+2198*b5*b8*a9-1372*b5*b4*a9-6146*b9^2*b10+210*b4^2*b6+11326*b4*b5*b6+5243*b8^2*b6+1400*a3*b8*b6-10766*a3*b9*b6-1162*a3*b8*b10+5810*a3*b10*b5+2926*a3*b10*b9+728*a3*b5*b6+5278*b9^2*b6-980*a3*b10*b4+42*b4^2*a8-882*b7*b6^2+84*b7*b10^2+6447*b10*b9*b8-2149*b5^2*b10+4172*b5^2*b6-4921*b8^2*b10-5747*b5*b6*b8+504*b5*b6*b9+3003*b10*b9*b5-924*b7*b6*b10+6244*a3*b4*b6-4578*b8*b9*b6+11011*b4*b10*b9-14903*b4*b10*b5+385*b4*b8*b10-14798*b4*b9*b6-2324*a7*b4*b6+490*a7*b9*b6+574*a7*b10*b9-2842*a7*b10*b5-490*a7*b8*b10+1064*a7*b8*b6+364*a7*b10*b4+1442*a7*b5*b6+7*b8^2*a4-581*b4*b8*b6-133*b4*b8*a4+3766*b5*b8*b10+3437*b9*a9*a3-2282*b9*a5*a3-5481*b9*b8*a9+1218*b9*b8*a8+5383*b9*b4*a5+8869*b9*b8*a5+2786*b9*b4*a9-7007*b5*a8*a3+3535*b5*a8*a7+10206*b5*a5*a3+3934*b9*a5*a7+3731*b9*a8*a3-3241*b9*a4*a3-3829*b9*a9*a7+2597*b9*a4*a7+245*b9*a8*a7-1512*b9*b7*a10+1064*b9*b7*a6"),
		SMQP("2121*b9^2*a9+16107*b9^2*a4+2751*b9^2*a8-5040*b4^2*b10-21686*b7*b4*a10+3850*b7*b10*a9-8806*b7*a10*a3+39795*b7*b4*a6-13790*b7*b10*a4-32403*b4*a4*a3-17934*b8*a8*a3-26523*b8*b9*a4+8652*b8*a4*a3+1155*b8*a9*a7-5838*b7*a6*a3+23793*b4*a8*a3-1848*b4*a8*a7+6426*b5^2*a5+14343*b4*b5*a4-12278*a7*b7*a6-5775*a7*b5*a9+3024*a7*b5*a5+6048*b8*b7*a10-2310*a7*b4*a5+5544*a7*b7*a10+2058*a7*b4*a9+7938*a7*b8*a5+1785*a7*b8*a8-672*b8*a9*a3+24822*b8*a5*a3-3549*b8*b4*a9+29190*b8*b4*a5+7287*b8*b5*a4-5817*b4*b8*a8-5691*b4*b9*a8+9744*b4*a4*a7+5208*b4*b5*a8-8673*a7*b8*a4+14553*b7*a4^2-22001*b7*a4*a8-7245*b8^2*a5-34818*a3*b4*a5+6531*a3*b4*a9-17941*b7*a5*a9+6888*b9^2*a5-3941*b7*a5*a8-32739*b5*a4*a3-18410*b7*b5*a10+17283*b5*a4*a7-1491*b5^2*a4+6510*b5^2*a8+1155*b5^2*a9-8400*b4*b9*a4+2919*b5*b9*a5+11025*b5*b8*a5-4893*b5*b4*a5+29799*b4^2*a9+2814*b5*b8*a8-11809*b7*b8*a6-28665*b4^2*a5+5082*b7*a9*a8+9870*b7*b6*a5+61425*b7*b5*a6-5215*b7*a9*a4+24927*b7*b6*a4-14259*b7*b6*a9+9002*b7*b10*a8-17262*b5*b9*a8-4067*b7*b6*a8+4676*b7*a9^2+16884*b7*a5^2-17192*b7*b10*a5-9261*b5*b9*a9+826*b7*a8^2+15099*a5*b7*a4-3696*b8^2*a9+1785*b8^2*a8+1701*b9*b5*a4+9597*b5*a9*a3-3486*b5*b8*a9+11886*b5*b4*a9+462*b9^2*b10-4410*b4^2*b6-54222*b4*b5*b6-8463*b8^2*b6+24108*a3*b8*b6+16590*a3*b9*b6-11298*a3*b8*b10+17178*a3*b10*b5-5082*a3*b10*b9-25032*a3*b5*b6-13902*b9^2*b6+5040*b4^2*a4+1596*a3*b10*b4+126*b4^2*a8+9450*b7*b6^2+5292*b7*b10^2+22617*b10*b9*b8+27321*b5^2*b10-40404*b5^2*b6+357*b8^2*b10+37191*b5*b6*b8+47376*b5*b6*b9-34587*b10*b9*b5-17388*b7*b6*b10-19572*a3*b4*b6-30114*b8*b9*b6-7791*b4*b10*b9+17283*b4*b10*b5-16485*b4*b8*b10+38766*b4*b9*b6+5628*a7*b4*b6-11130*a7*b9*b6+7014*a7*b10*b9-6762*a7*b10*b5-1722*a7*b8*b10+5460*a7*b8*b6-4956*a7*b10*b4+3822*a7*b5*b6+2541*b8^2*a4+41097*b4*b8*b6+2121*b4*b8*a4-12894*b5*b8*b10-8967*b9*a9*a3+1596*b9*a5*a3+8001*b9*b8*a9+10710*b9*b8*a8+9975*b9*b4*a5-21399*b9*b8*a5-23520*b9*b4*a9+24591*b5*a8*a3-6153*b5*a8*a7-19656*b5*a5*a3-10626*b9*a5*a7-8967*b9*a8*a3+19635*b9*a4*a3+4767*b9*a9*a7-9723*b9*a4*a7+4893*b9*a8*a7+11592*b9*b7*a10-29134*b9*b7*a6"),
		SMQP("14*a6*b5*b6+42*a6*b5*a5-14*a6*b4*b6+301*a6*b4*a5+7*a6*b5*a4-112*a6*b4*a4-7*b10*b5*a6-70*a6*b9*a4+14*a6*b9*a9-7*a6*b9*a8+7*b10*b9*a6+161*a6*b4*a9+7*a6*b8*a8+91*b8*a6*a4-406*a8*b4*a6+168*a4*b4*a10-14*b8*a6*a9-7*b8*b10*a6+49*a5*b9*a6+70*b9*b6*a6+7*a6*b5*a9-14*a6*b5*a8+63*a6^2*b7+196*b4*b10*a6-91*b8*b6*a6-105*b8*a6*a5-168*b4*b6*a10"),
		SMQP("-560*a5*b5*a8-9114*b4*a5*a4-3164*b8*b4*a10-3920*b9*b6*a9+1428*b10*b4*a5-2226*b7*a10*a5+1043*b10*b9*a9+84*b8*a6*a3-147*b8*a9*a4-4326*b7*a10*a4+12180*b7*a6*a4+2016*b7*a10*a9+7392*b8*b6*a5-567*b5*a4*a9-2436*b6*b5*a5+924*b7*a6*a5+581*b8*a8*a9+56*b8*a8^2+84*a6*b7*a8-2632*b8*a8*a4+182*b5*a8^2-532*b8*a10*a3-1386*b5*a4^2+1638*b4*a4^2-2009*b4*a9^2+721*a8*b9*a9+448*a8*b4*a9+1008*b4*a10*a7-1036*a8*b5*a9+714*b9*a10*a7+504*a8*b7*a10-231*a9^2*b5-1246*a9*b9*a4+518*a9^2*b9-2450*b6*b9*a8-2716*b4*a10*a3-182*b10*b9*a4+15862*b4^2*a10-2016*b4*b6*a4-2730*b7*b6*a10+9828*b7*b6*a6+2016*b7*b10*a10-1477*b8*b6*a9+35*b8*b10*a9+5992*b4*b6*a9+4984*a4*b4*a9-1106*b8*b5*a10+3654*b8*b9*a10-7140*b7*b10*a6-3486*b9*a4^2+3136*b10*b5*a5+4914*b5^2*a6-2688*b4*a6*a3-532*b8^2*a6-1568*b8*a6*a7-2030*b4*b10*a8-1015*b10*b5*a9-1260*b9*a6*a7-15078*b4*b9*a10+1932*b9*a6*a3+6020*b9*b8*a6+2590*b5^2*a10-4452*b8*b10*a4-8176*b4*b8*a6+16604*b4*b5*a10-1435*b8*a5*a9-3150*b10*b5*a8+147*b4*a5*a9+5166*b10*b5*a4-168*b9*a10*a3+27356*b4*b9*a6+7392*b8*b6*a4-924*b8*b6*a8+1064*b9*b10*a5+1190*b4*b6*a8+294*b4*a8*a4-5600*b8*b10*a5-13692*b4*b5*a6-20118*b4^2*a6-5320*b10^2*b4+931*a5*b9*a9+3220*b4*a6*a7-1932*a5^2*b4-4515*b7*a6*a9+1610*b9*b10*a8+1848*b9*b6*a5+4452*b9*b6*a4+7280*b4*a5*a8+5292*b8*a5*a4-420*b8^2*a10-1582*b9*a8*a4+798*b9*a5*a4+630*b8*a9^2-3178*b4*b10*a4+2268*b8*a5^2+1260*b8*a4^2-5880*b4*b6*a5-1218*b4*a8^2+11760*b9*b6^2+3136*b8*b10^2-12936*b4*b6^2-6104*b5*b10^2+2912*b10^2*b9-13552*b9*b6*b10-8288*b8*b10*b6+17528*b4*b6*b10-17304*b6^2*b5+6636*b6^2*b8+21532*b5*b6*b10+1344*b5*a6*a3-140*b5*a6*a7+4760*b5*b9*a6+4718*b5*a8*a4+3584*b5*a10*a3-1386*b5*a10*a7-7644*b5*b9*a10+1890*b5*a5*a4-7504*b5*b8*a6+1610*a5*b5*a9+518*b9*a8^2-3332*a5*b9*a8-2240*a5*b8*a8-210*b8*a10*a7-210*b10*b4*a9-756*b9*a5^2-6468*b5*b6*a4-7854*b9^2*a6+2310*b9^2*a10+2814*b5*b6*a8+1848*b8*b10*a8-1008*a5^2*b5+1652*b5*b6*a9"),
		SMQP("1855*b7*b4*a5-105*b7*b9*a5-1869*b7*b4*a9-2961*b7*b8*a5+924*a7*b9^2-336*a7*b8^2+1533*b7^2*a6+1974*b7*b8*a9-14994*b7*b5*a9-504*b7^2*a10-1533*a7*b7*a8-4704*a7*b7*a5-1407*a7*b7*a4-3780*b7*b4*a8+777*b7*b5*a8+4893*b7*a9*a7+12348*b7*b5*a4-2856*b7*a8*a3+1869*b7*a4*a3+903*b7*a5*a3+651*b7*a9*a3-5439*b7*b9*a4+756*b7*b8*a8+2415*b7*b9*a8-3024*a7*a3*b8+672*a3*b8^2+7434*b7*b5*b6+9156*b5^2*b8+3759*b7*b9*a9-2058*b5^3-672*a3*b7*b6-2688*b9*a3*b8-756*a7*b4*b9+2520*a3^2*b8+7056*a3*b5*a7-2982*b8*b7*b6+1554*b7*b8*b10-1050*b8*b5*a7-5292*a3^2*b5-3780*a7*b7*b6-168*b4*b8^2-4620*a3*a7*b9-672*a3*b9^2+2688*a3^2*b9-588*a7*b8*b9+3738*a7*b5*b9-924*a7*b4*b5+1680*a7*b7*b10-5880*a3^2*b4-6090*a7*b5^2+1512*b4^2*b9+10122*b8*b5*b9-6132*b4*b9^2-22050*b4*b5^2+9324*a7*b4*a3+12264*a3*b5^2+1512*b8^2*b9-1428*b4^2*a3-9366*b9^2*b5+336*b8^3+672*b4^2*a7-4326*a3*b8*b5-5334*b7*b5*b10+20160*b5*b9*b4+4116*a3*b4*b9-672*b9^2*b8-9156*b8^2*b5+588*b4*a7*b8-924*a3*b4*b5+2310*a3*b5*b9+3276*b7*b4*b6+168*b4^2*b5+756*a7^2*b8-1428*a7^2*b5+924*a7^2*b9-2352*a7^2*b4-1554*b7*b8*a4+3276*b7*b4*a4-1848*b7*b4*b10-2562*b5^2*b9+2184*b9^3+1260*b9*b7*b6+126*b9*b7*b10+2268*b4*b8*b5-1428*b4*b8*b9+4207*b7*b5*a5"),
		SMQP("-2611*b9^2*a9-889*b9^2*a4-3829*b9^2*a8+4354*b7*b4*a10-8414*b7*b10*a9+1442*b7*a10*a3-6699*b7*b4*a6+10402*b7*b10*a4+12831*b4*a4*a3+8190*b8*a8*a3+4039*b8*b9*a4-7602*b8*a4*a3-6363*b8*a9*a7+4494*b7*a6*a3-13923*b4*a8*a3+1218*b4*a8*a7-2646*b5^2*a5+2079*b4*b5*a4+4774*a7*b7*a6+1449*a7*b5*a9-7056*a7*b5*a5-2688*b8*b7*a10-7602*a7*b4*a5-3864*a7*b7*a10+7644*a7*b4*a9+1638*a7*b8*a5-189*a7*b8*a8+3402*b8*a9*a3-8316*b8*a5*a3+4809*b8*b4*a9-5250*b8*b4*a5-4809*b8*b5*a4+651*b4*b8*a8+4403*b4*b9*a8-8274*b4*a4*a7+8358*b4*b5*a8+4557*a7*b8*a4-3843*b7*a4^2+3787*b7*a4*a8+189*b8^2*a5+14700*a3*b4*a5-6867*a3*b4*a9+3479*b7*a5*a9-4928*b9^2*a5-329*b7*a5*a8+8085*b5*a4*a3+15694*b7*b5*a10-7917*b5*a4*a7-6363*b5^2*a4+2814*b5^2*a8-1533*b5^2*a9-6202*b4*b9*a4-2793*b5*b9*a5-1953*b5*b8*a5-7119*b5*b4*a5-7875*b4^2*a9-4536*b5*b8*a8+7973*b7*b8*a6+2205*b4^2*a5-1806*b7*a9*a8-1722*b7*b6*a5-39879*b7*b5*a6+1589*b7*a9*a4-18669*b7*b6*a4+14049*b7*b6*a9-1918*b7*b10*a8+1414*b5*b9*a8+721*b7*b6*a8-2380*b7*a9^2-8820*b7*a5^2+6496*b7*b10*a5+4711*b5*b9*a9+2002*b7*a8^2-1281*a5*b7*a4-1092*b8^2*a9-2919*b8^2*a8+4081*b9*b5*a4-651*b5*a9*a3+9450*b5*b8*a9-26922*b5*b4*a9-8498*b9^2*b10+630*b4^2*b6+40698*b4*b5*b6+14637*b8^2*b6+1176*a3*b8*b6-18186*a3*b9*b6+1890*a3*b8*b10+2982*a3*b10*b5+5082*a3*b10*b9-4872*a3*b5*b6+14098*b9^2*b6-3276*a3*b10*b4+126*b4^2*a8-630*b7*b6^2-1092*b7*b10^2+1757*b10*b9*b8-11319*b5^2*b10+21924*b5^2*b6-5775*b8^2*b10-24297*b5*b6*b8-13888*b5*b6*b9+15365*b10*b9*b5+2940*b7*b6*b10+9660*a3*b4*b6-6706*b8*b9*b6+18361*b4*b10*b9-20517*b4*b10*b5+1491*b4*b8*b10-37394*b4*b9*b6-4956*a7*b4*b6+1470*a7*b9*b6-1974*a7*b10*b9-2142*a7*b10*b5+882*a7*b8*b10+1176*a7*b8*b6+2436*a7*b10*b4+4326*a7*b5*b6+1113*b8^2*a4-3927*b4*b8*b6-735*b4*b8*a4+7686*b5*b8*b10+2373*b9*a9*a3-7728*b9*a5*a3-3017*b9*b8*a9+3556*b9*b8*a8+637*b9*b4*a5+11739*b9*b8*a5+11018*b9*b4*a9-8589*b5*a8*a3+2331*b5*a8*a7+15624*b5*a5*a3+7434*b9*a5*a7+3129*b9*a8*a3-2793*b9*a4*a3-2205*b9*a9*a7+2961*b9*a4*a7+1197*b9*a8*a7-3192*b9*b7*a10+3626*b9*b7*a6"),
		SMQP("168*a5*b5*a8-588*b4*a5*a4+924*b8*b4*a10-322*b9*b6*a9+287*b10*b4*a5-504*b7*a10*a5+154*b10*b9*a9+252*b8*a6*a3+224*b8*a9*a4+168*b7*a10*a4-42*b7*a6*a4-84*b7*a10*a9+189*b8*b6*a5+14*b5*a4*a9-378*b6*b5*a5+168*b7*a6*a5-210*b8*a8*a9+112*b8*a8^2-126*a6*b7*a8+616*b8*a8*a4-56*b5*a8^2+210*b8*a10*a3+504*b4*a4^2-42*b4*a9^2+84*a8*b9*a9+336*a8*b4*a9+336*b4*a10*a7-210*a8*b5*a9+462*b9*a10*a7+84*a8*b7*a10+56*a9^2*b5+112*a9*b9*a4+280*a9^2*b9+301*b6*b9*a8-406*b4*a10*a3+56*b10*b9*a4-518*b4^2*a10-336*b4*b6*a4-168*b7*b6*a10-189*b7*b6*a6-56*b8*b6*a9-28*b8*b10*a9-679*b4*b6*a9-210*a4*b4*a9+420*b8*b5*a10+882*b8*b9*a10+105*b7*b10*a6+42*b9*a4^2+210*b10*b5*a5-420*b4*a6*a3-126*b8^2*a6+266*b4*b10*a8-259*b10*b5*a9-126*b9*a6*a7-714*b4*b9*a10+42*b9*a6*a3+1554*b9*b8*a6+294*b5^2*a10+49*b8*b10*a4+224*b4*b8*a6-980*b4*b5*a10+588*b8*a5*a9-686*b4*a5*a9+7*b10*b5*a4-378*b9*a10*a3+770*b4*b9*a6+189*b8*b6*a4-511*b8*b6*a8+217*b9*b10*a5+308*b4*b6*a8-546*b4*a8*a4-273*b8*b10*a5+252*b4*b5*a6+168*b4^2*a6+112*b10^2*b4-112*a5*b9*a9+196*b4*a6*a7+504*a5^2*b4-210*b7*a6*a9-147*b9*b10*a8-231*b9*b6*a5+210*b9*b6*a4+1022*b4*a5*a8+126*b8*a5*a4-756*b8^2*a10-322*b9*a8*a4+378*b9*a5*a4-280*b8*a9^2-168*b4*b10*a4-504*b8*a5^2-378*b8*a4^2-273*b4*b6*a5-252*b4*a8^2-133*b8*b10^2+126*b4*b6^2-49*b5*b10^2+133*b10^2*b9-287*b9*b6*b10+224*b8*b10*b6-14*b4*b6*b10-126*b6^2*b5-189*b6^2*b8+203*b5*b6*b10-336*b5*b9*a6+112*b5*a8*a4+1092*b5*a10*a3-462*b5*a10*a7-714*b5*b9*a10-630*b5*b8*a6-84*a5*b5*a9+14*b9*a8^2-224*a5*b9*a8-588*a5*b8*a8-462*b8*a10*a7+329*b10*b4*a9+168*b9*a5^2-63*b5*b6*a4-924*b9^2*a6-126*b9^2*a10-28*b5*b6*a8+21*b8*b10*a8+595*b5*b6*a9"),
		SMQP("-399*b9^2*a9+399*b9^2*a4-525*b9^2*a8+1526*b7*b4*a10+14*b7*b10*a9-98*b7*a10*a3-2919*b7*b4*a6+518*b7*b10*a4+1239*b4*a4*a3+126*b8*a8*a3+833*b8*b9*a4-168*b8*a4*a3+203*b8*a9*a7-1386*b7*a6*a3+1197*b4*a8*a3-518*b4*a8*a7+378*b5^2*a5+399*b4*b5*a4+1190*a7*b7*a6+245*a7*b5*a9+84*a7*b5*a5-336*b8*b7*a10-490*a7*b4*a5-168*a7*b7*a10+28*a7*b4*a9+210*a7*b8*a5-175*a7*b8*a8-630*b8*a5*a3-1477*b8*b4*a9-1400*b8*b4*a5-1421*b8*b5*a4+203*b4*b8*a8+2331*b4*b9*a8-602*b4*a4*a7-2394*b4*b5*a8-49*a7*b8*a4-231*b7*a4^2+959*b7*a4*a8+1743*b8^2*a5+1344*a3*b4*a5-567*a3*b4*a9-749*b7*a5*a9-756*b9^2*a5+2303*b7*a5*a8+525*b5*a4*a3+1106*b7*b5*a10-301*b5*a4*a7+609*b5^2*a4-882*b5^2*a8+63*b5^2*a9-1218*b4*b9*a4+399*b5*b9*a5-2919*b5*b8*a5+1365*b5*b4*a5-2331*b4^2*a9+910*b5*b8*a8-413*b7*b8*a6+1701*b4^2*a5-210*b7*a9*a8-126*b7*b6*a5-1239*b7*b5*a6+1225*b7*a9*a4-1197*b7*b6*a4-63*b7*b6*a9-266*b7*b10*a8+1050*b5*b9*a8+245*b7*b6*a8-308*b7*a9^2+336*b7*a5^2-280*b7*b10*a5+483*b5*b9*a9-322*b7*a8^2-2121*a5*b7*a4+616*b8^2*a9-707*b8^2*a8-1071*b9*b5*a4-567*b5*a9*a3-434*b5*b8*a9-126*b5*b4*a9-378*b9^2*b10+630*b4^2*b6+3066*b4*b5*b6+245*b8^2*b6-1344*a3*b8*b6-1050*a3*b9*b6+882*a3*b8*b10-882*a3*b10*b5+210*a3*b10*b9+1680*a3*b5*b6+210*b9^2*b6+252*a3*b10*b4+126*b4^2*a8-630*b7*b6^2-420*b7*b10^2-371*b10*b9*b8-819*b5^2*b10+2100*b5^2*b6-343*b8^2*b10-1645*b5*b6*b8-2184*b5*b6*b9+1113*b10*b9*b5+1092*b7*b6*b10+84*a3*b4*b6+1582*b8*b9*b6+1449*b4*b10*b9-1701*b4*b10*b5+595*b4*b8*b10-2562*b4*b9*b6-28*a7*b4*b6+350*a7*b9*b6-70*a7*b10*b9+322*a7*b10*b5-350*a7*b8*b10+448*a7*b8*b6-28*a7*b10*b4-602*a7*b5*b6+49*b8^2*a4-1799*b4*b8*b6+161*b4*b8*a4+938*b5*b8*b10+777*b9*a9*a3-84*b9*a5*a3-343*b9*b8*a9+182*b9*b8*a8-2331*b9*b4*a5+1561*b9*b8*a5+1890*b9*b4*a9-189*b5*a8*a3+35*b5*a8*a7-252*b5*a5*a3+14*b9*a5*a7+441*b9*a8*a3-609*b9*a4*a3-497*b9*a9*a7+553*b9*a4*a7-35*b9*a8*a7-840*b9*b7*a10+1750*b9*b7*a6"),
		SMQP("2877*b7*b4*a5+3045*b7*b9*a5-9135*b7*b4*a9-3591*b7*b8*a5+2268*a7*b9^2+252*a7*b8^2-945*b7^2*a6-126*b7*b8*a9-19824*b7*b5*a9+1008*b7^2*a10-1449*a7*b7*a8-2016*a7*b7*a5-5607*a7*b7*a4-1008*b7*b4*a8+4641*b7*b5*a8+7497*b7*a9*a7+9156*b7*b5*a4-6552*b7*a8*a3+14049*b7*a4*a3+567*b7*a5*a3-1449*b7*a9*a3-10395*b7*b9*a4-882*b7*b8*a8+5481*b7*b9*a8-1008*a7*a3*b8-5544*a3*b8^2+10416*b7*b5*b6+4788*b5^2*b8+6615*b7*b9*a9-4284*b5^3-672*a3*b7*b6+756*b9*a3*b8-2268*a7*b4*b9-1008*a3^2*b8+3528*a3*b5*a7-5124*b8*b7*b6+2016*b7*b8*b10-3276*b8*b5*a7+1512*a3^2*b5-3024*a7*b7*b6-3528*b4*b8^2-8064*a3*a7*b9-8316*a3*b9^2+6048*a3^2*b9-504*a7*b8*b9+1260*a7*b5*b9-756*a7*b4*b5+2016*a7*b7*b10-5040*a3^2*b4-2016*a7*b5^2+7560*b4^2*b9+1512*b8*b5*b9-17388*b4*b9^2-2016*a3*b7*b10-15876*b4*b5^2+9072*a7*b4*a3+3276*a3*b5^2+8820*b8^2*b9+252*b4^2*a3-4788*b9^2*b5-3528*b8^3+10080*a3*b8*b5-7056*b7*b5*b10+25956*b5*b9*b4+16128*a3*b4*b9-5040*b9^2*b8-1764*b8^2*b5+252*b4*a7*b8-8568*a3*b4*b5+5292*a3*b5*b9+84*b7*b4*b6-11592*b4^2*b5-1008*a7^2*b5+2016*a7^2*b9-2016*a7^2*b4+1764*b7*b8*a4+5796*b7*b4*a4-1008*b7*b4*b10-504*b5^2*b9+3780*b9^3-924*b9*b7*b6+2016*b9*b7*b10-8316*b4*b8*b5+756*b4*a3*b8+10836*b4*b8*b9+833*b7*b5*a5"),
		SMQP("-1190*a5*b5*a8-12747*b4*a5*a4-1484*b8*b4*a10-3710*b9*b6*a9+1855*b10*b4*a5-2226*b7*a10*a5+1134*b10*b9*a9+84*b8*a6*a3+210*b8*a9*a4-4326*b7*a10*a4+10731*b7*a6*a4+2016*b7*a10*a9+7455*b8*b6*a5-861*b5*a4*a9-2562*b6*b5*a5+924*b7*a6*a5+1022*b8*a8*a9-133*b8*a8^2-861*a6*b7*a8-4606*b8*a8*a4+392*b5*a8^2-700*b8*a10*a3-1407*b5*a4^2+588*b4*a4^2-1680*b4*a9^2+322*a8*b9*a9+301*a8*b4*a9+1008*b4*a10*a7-1099*a8*b5*a9+714*b9*a10*a7+504*a8*b7*a10-252*a9^2*b5-1946*a9*b9*a4+616*a9^2*b9-4165*b6*b9*a8-2786*b4*a10*a3-595*b10*b9*a4+17696*b4^2*a10-2310*b4*b6*a4-2730*b7*b6*a10+9387*b7*b6*a6+2016*b7*b10*a10-1330*b8*b6*a9-70*b8*b10*a9+5005*b4*b6*a9+5901*a4*b4*a9-938*b8*b5*a10+3654*b8*b9*a10-6951*b7*b10*a6-2142*b9*a4^2+3262*b10*b5*a5+4914*b5^2*a6-2352*b4*a6*a3-4480*b8^2*a6-1400*b8*a6*a7-2814*b4*b10*a8-973*b10*b5*a9-1568*b9*a6*a7-13832*b4*b9*a10+2016*b9*a6*a3+6580*b9*b8*a6+2590*b5^2*a10-3990*b8*b10*a4-15092*b4*b8*a6+16674*b4*b5*a10-1288*b8*a5*a9-3087*b10*b5*a8-700*b4*a5*a9+5460*b10*b5*a4+98*b9*a10*a3+29680*b4*b9*a6+10248*b8*b6*a4+672*b8*b6*a8+1407*b9*b10*a5+2772*b4*b6*a8+490*b4*a8*a4-5999*b8*b10*a5-12306*b4*b5*a6-21840*b4^2*a6-5236*b10^2*b4+1680*a5*b9*a9+3500*b4*a6*a7-1428*a5^2*b4-4704*b7*a6*a9+1414*b9*b10*a8+1701*b9*b6*a5+3276*b9*b6*a4+5397*b4*a5*a8+8127*b8*a5*a4-420*b8^2*a10+483*b9*a8*a4+189*b9*a5*a4+504*b8*a9^2-4200*b4*b10*a4+2268*b8*a5^2-1281*b8*a4^2-6279*b4*b6*a5-14*b4*a8^2+11382*b9*b6^2+3115*b8*b10^2-13314*b4*b6^2-6125*b5*b10^2+2933*b10^2*b9-13279*b9*b6*b10-8624*b8*b10*b6+17990*b4*b6*b10-16926*b6^2*b5+6573*b6^2*b8+21511*b5*b6*b10+1344*b5*a6*a3-140*b5*a6*a7+5390*b5*b9*a6+5159*b5*a8*a4+3584*b5*a10*a3-1386*b5*a10*a7-7910*b5*b9*a10+756*b5*a5*a4-5740*b5*b8*a6+1484*a5*b5*a9+693*b9*a8^2-5299*a5*b9*a8-245*a5*b8*a8-210*b8*a10*a7-161*b10*b4*a9-1260*b9*a5^2-7455*b5*b6*a4-7252*b9^2*a6+2310*b9^2*a10+2478*b5*b6*a8+2058*b8*b10*a8-1008*a5^2*b5+1673*b5*b6*a9"),
		SMQP("-903*b9^2*a9+231*b9^2*a4-1197*b9^2*a8+1414*b7*b4*a10+70*b7*b10*a9-154*b7*a10*a3-2394*b7*b4*a6+826*b7*b10*a4+1554*b4*a4*a3+1302*b8*a8*a3+1491*b8*b9*a4-966*b8*a4*a3-1155*b8*a9*a7-1218*b7*a6*a3+357*b4*a8*a3-777*b4*a8*a7+378*b5^2*a5+651*b4*b5*a4+1414*a7*b7*a6-819*a7*b5*a9+84*a7*b5*a5-1008*b8*b7*a10-1218*a7*b4*a5-168*a7*b7*a10+1743*a7*b4*a9+1050*a7*b8*a5-945*a7*b8*a8+126*b8*a9*a3-840*b8*a5*a3-525*b8*b4*a9-903*b8*b4*a5-1365*b8*b5*a4+1365*b4*b8*a8+3150*b4*b9*a8-567*b4*a4*a7-4263*b4*b5*a8+273*a7*b8*a4-231*b7*a4^2+1351*b7*a4*a8+2289*b8^2*a5+1029*a3*b4*a5-714*a3*b4*a9-1645*b7*a5*a9+84*b9^2*a5+1183*b7*a5*a8+525*b5*a4*a3+1162*b7*b5*a10-189*b5*a4*a7+609*b5^2*a4-378*b5^2*a8+735*b5^2*a9-1407*b4*b9*a4+231*b5*b9*a5-2709*b5*b8*a5+672*b5*b4*a5-2520*b4^2*a9+1344*b5*b8*a8+329*b7*b8*a6+1512*b4^2*a5-378*b7*a9*a8+462*b7*b6*a5-1743*b7*b5*a6+1337*b7*a9*a4-1617*b7*b6*a4-231*b7*b6*a9-322*b7*b10*a8+378*b5*b9*a8-35*b7*b6*a8-196*b7*a9^2+336*b7*a5^2-644*b7*b10*a5-525*b5*b9*a9-266*b7*a8^2-1785*a5*b7*a4-420*b8^2*a9-2415*b8^2*a8-1407*b9*b5*a4+1785*b5*a9*a3+294*b5*b8*a9-2814*b5*b4*a9-798*b9^2*b10+630*b4^2*b6+5586*b4*b5*b6+2205*b8^2*b6-840*a3*b8*b6-1302*a3*b9*b6+378*a3*b8*b10+210*a3*b10*b5+462*a3*b10*b9-84*a3*b5*b6+966*b9^2*b6+504*a3*b10*b4+126*b4^2*a8-378*b7*b6^2-420*b7*b10^2+1281*b10*b9*b8-651*b5^2*b10+1596*b5^2*b6-1575*b8^2*b10-2625*b5*b6*b8-420*b5*b6*b9+21*b10*b9*b5+840*b7*b6*b10-168*a3*b4*b6-1470*b8*b9*b6+2289*b4*b10*b9-3213*b4*b10*b5+819*b4*b8*b10-4074*b4*b9*b6-252*a7*b4*b6+378*a7*b9*b6-42*a7*b10*b9-42*a7*b10*b5-378*a7*b8*b10+756*a7*b8*b6-252*a7*b10*b4-126*a7*b5*b6-63*b8^2*a4-1575*b4*b8*b6-63*b4*b8*a4+1806*b5*b8*b10-399*b9*a9*a3+84*b9*a5*a3+1575*b9*b8*a9+2184*b9*b8*a8-1176*b9*b4*a5-693*b9*b8*a5+1239*b9*b4*a9+1827*b5*a8*a3-693*b5*a8*a7-252*b5*a5*a3-42*b9*a5*a7-1071*b9*a8*a3+231*b9*a4*a3+567*b9*a9*a7-63*b9*a4*a7+1029*b9*a8*a7-168*b9*b7*a10+854*b9*b7*a6"),
		SMQP("-574*a5*b5*a8-31710*b4*a5*a4-12670*b8*b4*a10-5880*b9*b6*a9+2219*b10*b4*a5-4284*b7*a10*a5+1386*b10*b9*a9+924*b8*a6*a3+2128*b8*a9*a4-8316*b7*a10*a4+20286*b7*a6*a4+4536*b7*a10*a9+13146*b8*b6*a5-1274*b5*a4*a9-5124*b6*b5*a5+1008*b7*a6*a5+672*b8*a8*a9+2093*b8*a8^2-189*a6*b7*a8-10297*b8*a8*a4-994*b5*a8^2+322*b8*a10*a3-4914*b5*a4^2+1260*b4*a4^2-3668*b4*a9^2+1302*a8*b9*a9+2205*a8*b4*a9+2016*b4*a10*a7-1281*a8*b5*a9+756*b9*a10*a7+504*a8*b7*a10-644*a9^2*b5-3108*a9*b9*a4+672*a9^2*b9-5964*b6*b9*a8-6524*b4*a10*a3+2604*b10*b9*a4+28448*b4^2*a10-6972*b4*b6*a4-5292*b7*b6*a10+16758*b7*b6*a6+3528*b7*b10*a10-2226*b8*b6*a9-1372*b8*b10*a9+6930*b4*b6*a9+9310*a4*b4*a9-3598*b8*b5*a10+8316*b8*b9*a10-15183*b7*b10*a6-2268*b9*a4^2+7238*b10*b5*a5+11592*b5^2*a6-4368*b4*a6*a3+2128*b8^2*a6-3556*b8*a6*a7-3206*b4*b10*a8-1687*b10*b5*a9-2856*b9*a6*a7-22008*b4*b9*a10+4032*b9*a6*a3+322*b9*b8*a6+5348*b5^2*a10-10843*b8*b10*a4-8666*b4*b8*a6+36568*b4*b5*a10-644*b8*a5*a9-6545*b10*b5*a8-6020*b4*a5*a9+8435*b10*b5*a4+840*b9*a10*a3+48328*b4*b9*a6+25032*b8*b6*a4-4627*b8*b6*a8+3255*b9*b10*a5+6902*b4*b6*a8+2576*b4*a8*a4-10591*b8*b10*a5-28392*b4*b5*a6-42000*b4^2*a6-10500*b10^2*b4+3612*a5*b9*a9+6944*b4*a6*a7+168*a5^2*b4-9828*b7*a6*a9+1806*b9*b10*a8+3570*b9*b6*a5+3276*b9*b6*a4+20153*b4*a5*a8+22932*b8*a5*a4+504*b8^2*a10-252*b9*a8*a4-3318*b9*a5*a4+868*b8*a9^2-8008*b4*b10*a4+4032*b8*a5^2-3024*b8*a4^2-9702*b4*b6*a5-2926*b4*a8^2+17388*b9*b6^2+9303*b8*b10^2-29484*b4*b6^2-11361*b5*b10^2+2793*b10^2*b9-20832*b9*b6*b10-23835*b8*b10*b6+38262*b4*b6*b10-30996*b6^2*b5+15498*b6^2*b8+43764*b5*b6*b10+2688*b5*a6*a3-280*b5*a6*a7+15232*b5*b9*a6+11627*b5*a8*a4+5992*b5*a10*a3-2268*b5*a10*a7-14952*b5*b9*a10-2016*b5*a5*a4-21434*b5*b8*a6+3136*a5*b5*a9+1449*b9*a8^2-7707*a5*b9*a8-8449*a5*b8*a8-756*b8*a10*a7-721*b10*b4*a9-2856*b9*a5^2-15414*b5*b6*a4-10248*b9^2*a6+3276*b9^2*a10+6398*b5*b6*a8+6748*b8*b10*a8-2016*a5^2*b5+2562*b5*b6*a9"),
		SMQP("1071*b7*b9*a5-756*b7*b4*a9-315*b7*b8*a5+567*a7*b9^2+63*a7*b8^2-252*b7^2*a6-504*b7*b8*a9-1680*b7*b5*a9+252*b7^2*a10-126*a7*b7*a8+378*a7*b7*a5-630*a7*b7*a4+252*b7*b4*a8+798*b7*b5*a8+630*b7*a9*a7-420*b7*b5*a4-504*b7*a8*a3+1890*b7*a4*a3+1008*b7*a5*a3+126*b7*a9*a3-1071*b7*b9*a4-252*b7*b8*a8+630*b7*b9*a8-882*a3*b8^2+1344*b7*b5*b6-1323*b5^2*b8+378*b7*b9*a9-1071*b5^3-672*a3*b7*b6+441*b9*a3*b8-504*a7*b4*b9-1134*a3*b5*a7-273*b8*b7*b6-315*b8*b5*a7+1134*a3^2*b5+252*a7*b7*b6-2331*a3*b9^2-126*a7*b8*b9-1197*a7*b5*b9+756*a7*b4*b5+1008*a7*b5^2-2142*b8*b5*b9-1386*b4*b9^2+1197*b4*b5^2-1449*a3*b5^2+1197*b8^2*b9+1827*b9^2*b5-882*b8^3+3024*a3*b8*b5-756*b7*b5*b10-378*b5*b9*b4+3213*a3*b4*b9-252*b9^2*b8+1575*b8^2*b5-2331*a3*b4*b5+1575*a3*b5*b9-42*b7*b4*b6+252*a7^2*b5+945*b7*b8*a4+882*b5^2*b9-63*b9^3-735*b9*b7*b6+504*b9*b7*b10-2898*b4*b8*b5+63*b4*a3*b8+1890*b4*b8*b9+182*b7*b5*a5"),
		SMQP("105*b9^2*a9-2793*b9^2*a4-903*b9^2*a8+770*b7*b4*a10-1582*b7*b10*a9+2002*b7*a10*a3+903*b7*b4*a6+1946*b7*b10*a4+2625*b4*a4*a3+966*b8*a8*a3+3003*b8*b9*a4+336*b8*a4*a3-609*b8*a9*a7+4326*b7*a6*a3-3675*b4*a8*a3+2940*b4*a8*a7-3906*b5^2*a5-1785*b4*b5*a4-2842*a7*b7*a6+777*a7*b5*a9-1764*a7*b5*a5+1008*b8*b7*a10-546*a7*b4*a5-504*a7*b7*a10-1470*a7*b4*a9-1134*a7*b8*a5-735*a7*b8*a8+336*b8*a9*a3-1134*b8*a5*a3+1533*b8*b4*a9-2100*b8*b4*a5-693*b8*b5*a4-357*b4*b8*a8-651*b4*b9*a8-2604*b4*a4*a7+168*b4*b5*a8-357*a7*b8*a4-1575*b7*a4^2+2723*b7*a4*a8-567*b8^2*a5+3234*a3*b4*a5-273*a3*b4*a9+4711*b7*a5*a9-1806*b9^2*a5-6769*b7*a5*a8+5313*b5*a4*a3+1022*b7*b5*a10-4137*b5*a4*a7-861*b5^2*a4+4242*b5^2*a8-1743*b5^2*a9+1932*b4*b9*a4-3255*b5*b9*a5+7623*b5*b8*a5-8295*b5*b4*a5-693*b4^2*a9-2268*b5*b8*a8+469*b7*b8*a6+1323*b4^2*a5-798*b7*a9*a8-546*b7*b6*a5-4977*b7*b5*a6-1211*b7*a9*a4-3381*b7*b6*a4+3297*b7*b6*a9-686*b7*b10*a8-2268*b5*b9*a8-2359*b7*b6*a8-140*b7*a9^2-3528*b7*a5^2+4088*b7*b10*a5-315*b5*b9*a9+1442*b7*a8^2+3927*a5*b7*a4-1092*b8^2*a9-777*b8^2*a8+819*b9*b5*a4-483*b5*a9*a3+3024*b5*b8*a9-3360*b5*b4*a9-1176*b9^2*b10-378*b4^2*b6+6510*b4*b5*b6+1155*b8^2*b6+672*a3*b8*b6-3822*a3*b9*b6+1050*a3*b8*b10-714*a3*b10*b5-294*a3*b10*b9+1428*a3*b5*b6+2730*b9^2*b6-420*a3*b10*b4-882*b4^2*a8+630*b7*b6^2+252*b7*b10^2-1995*b10*b9*b8-2037*b5^2*b10-2352*b5^2*b6+1407*b8^2*b10+1953*b5*b6*b8-252*b5*b6*b9+6489*b10*b9*b5-252*b7*b6*b10+2352*a3*b4*b6-924*b8*b9*b6+777*b4*b10*b9-3381*b4*b10*b5+483*b4*b8*b10-4830*b4*b9*b6-924*a7*b4*b6-798*a7*b9*b6-798*a7*b10*b9+42*a7*b10*b5+1050*a7*b8*b10-588*a7*b8*b6+1092*a7*b10*b4+1554*a7*b5*b6-1533*b8^2*a4-2163*b4*b8*b6+273*b4*b8*a4-2898*b5*b8*b10-147*b9*a9*a3-1428*b9*a5*a3-1155*b9*b8*a9+3318*b9*b8*a8+2793*b9*b4*a5+4137*b9*b8*a5+3066*b9*b4*a9-1365*b5*a8*a3+147*b5*a8*a7+4788*b5*a5*a3+1470*b9*a5*a7+2121*b9*a8*a3-3549*b9*a4*a3-21*b9*a9*a7+3381*b9*a4*a7-147*b9*a8*a7-1512*b9*b7*a10+1036*b9*b7*a6"),
		SMQP("252*b9^2*a9-700*b7*b4*a10-28*b7*b10*a9+700*b7*a10*a3+147*b7*b4*a6-28*b7*b10*a4+441*b4*a4*a3+252*b8*a9*a7+3192*b7*a6*a3-504*b4*a8*a3+63*b4*a8*a7+2583*b5^2*a5-1260*b4*b5*a4+728*a7*b7*a6-2205*a7*b5*a9+2772*a7*b5*a5+504*a7*b4*a5-63*a7*b4*a9+1260*a7*b8*a5-126*b8*b4*a9-315*b8*b4*a5+252*b8*b5*a4+126*b4*b8*a8-63*b4*b9*a8-63*b4*a4*a7+945*b4*b5*a8+756*b7*a4^2-700*b7*a4*a8+2520*b8^2*a5-1701*a3*b4*a5-189*a3*b4*a9-1316*b7*a5*a9+1260*b9^2*a5+1904*b7*a5*a8+63*b5*a4*a3+308*b7*b5*a10-1449*b5*a4*a7-1512*b5^2*a4-441*b5^2*a8+2268*b5^2*a9-63*b4*b9*a4+5103*b5*b9*a5-7245*b5*b8*a5+6804*b5*b4*a5-63*b4^2*a9+1134*b5*b8*a8+616*b7*b8*a6-63*b4^2*a5+84*b7*a9*a8+6132*b7*b6*a5-5607*b7*b5*a6-308*b7*a9*a4-672*b7*b6*a4-1680*b7*b6*a9+28*b7*b10*a8-693*b5*b9*a8+1652*b7*b6*a8-56*b7*a9^2+2268*b7*a5^2-280*b7*b10*a5-819*b5*b9*a9-28*b7*a8^2-2184*a5*b7*a4+504*b8^2*a9+1071*b9*b5*a4+693*b5*a9*a3-630*b5*b8*a9+315*b5*b4*a9-6048*b4*b5*b6+1512*a3*b8*b6-5040*a3*b9*b6+4536*a3*b5*b6-2016*b9^2*b6-2016*b7*b6^2-252*b5^2*b10-1008*b5^2*b6-2772*b5*b6*b8+6552*b5*b6*b9+252*b10*b9*b5+1008*b7*b6*b10+2520*a3*b4*b6+1008*b8*b9*b6+1008*b4*b10*b5+1008*b4*b9*b6-1008*a7*b4*b6+1008*a7*b8*b6-1008*a7*b5*b6-252*b5*b8*b10+252*b9*a9*a3+1260*b9*a5*a3-756*b9*b8*a9-2079*b9*b4*a5-3780*b9*b8*a5-441*b9*b4*a9-1512*b5*a8*a3+2457*b5*a8*a7-3843*b5*a5*a3-1260*b9*a5*a7-252*b9*a9*a7-644*b9*b7*a6"),
		SMQP("-1330*a5*b5*a8+546*b4*a5*a4-364*b8*b4*a10-3276*b9*b6*a9+1652*b10*b4*a5+1764*b10*b9*a9+3444*b8*a6*a3-392*b8*a9*a4-504*b7*a10*a4+126*b7*a6*a4+6027*b8*b6*a5+952*b5*a4*a9+1050*b6*b5*a5+2079*b7*a6*a5+84*b8*a8*a9-196*b8*a8^2-2016*a6*b7*a8-364*b8*a8*a4+476*b5*a8^2-1148*b8*a10*a3-2772*b5*a4^2-3906*b4*a4^2-308*b4*a9^2+252*a8*b9*a9-546*a8*b4*a9-1008*b4*a10*a7-420*a8*b5*a9+1008*b9*a10*a7+504*a8*b7*a10-1064*a9^2*b5-504*a9*b9*a4+504*a9^2*b9-504*b6*b9*a8+2842*b4*a10*a3+1260*b10*b9*a4+3710*b4^2*a10+1848*b4*b6*a4-4032*b7*b6*a10+12222*b7*b6*a6+2016*b7*b10*a10-1428*b8*b6*a9+308*b8*b10*a9+7770*b4*b6*a9+700*a4*b4*a9-1372*b8*b5*a10+1512*b8*b9*a10-5796*b7*b10*a6+1764*b9*a4^2+623*b10*b5*a5-504*b5^2*a6-714*b4*a6*a3+2548*b8^2*a6-70*b8*a6*a7+1162*b4*b10*a8+224*b10*b5*a9-2394*b9*a6*a7-2772*b4*b9*a10+1386*b9*a6*a3-7154*b9*b8*a6+1316*b5^2*a10-952*b8*b10*a4-2072*b4*b8*a6+7546*b4*b5*a10+1414*b8*a5*a9-1988*b10*b5*a8+7399*b4*a5*a9-280*b10*b5*a4-1764*b9*a10*a3+2506*b4*b9*a6+2100*b8*b6*a4+1736*b8*b6*a8-651*b9*b10*a5-5530*b4*b6*a8+3710*b4*a8*a4-3325*b8*b10*a5-6090*b4*b5*a6-1050*b4^2*a6-3822*a5*b9*a9-532*b4*a6*a7-11865*a5^2*b4+1008*b7*a6*a9+1008*b9*b10*a8+1050*b9*b6*a5-3528*b9*b6*a4-2422*b4*a5*a8-1827*b8*a5*a4-1764*b9*a8*a4+966*b9*a5*a4+616*b8*a9^2+350*b4*b10*a4+6489*b8*a5^2+252*b8*a4^2-1890*b4*b6*a5+350*b4*a8^2+8568*b9*b6^2+756*b8*b10^2-3528*b4*b6^2-3780*b5*b10^2+2268*b10^2*b9-8820*b9*b6*b10-4032*b8*b10*b6+2520*b4*b6*b10-11592*b6^2*b5+4284*b6^2*b8+13356*b5*b6*b10-6762*b5*a6*a3+3122*b5*a6*a7+11578*b5*b9*a6+2828*b5*a8*a4+3220*b5*a10*a3-1512*b5*a10*a7-1764*b5*b9*a10+693*b5*a5*a4-9338*b5*b8*a6+4081*a5*b5*a9-252*b9*a8^2+2415*a5*b9*a8-3199*a5*b8*a8+504*b8*a10*a7-3934*b10*b4*a9-777*b9*a5^2+3108*b5*b6*a4+2142*b9^2*a6+644*b5*b6*a8-560*b8*b10*a8-5418*a5^2*b5+336*b5*b6*a9"),
		SMQP("-2989*b9^2*a9-343*b9^2*a4-2107*b9^2*a8+1526*b7*b4*a10-8218*b7*b10*a9-1946*b7*a10*a3-4074*b7*b4*a6+9926*b7*b10*a4+10094*b4*a4*a3+5362*b8*a8*a3+3444*b8*b9*a4-4991*b8*a4*a3-4634*b8*a9*a7+3066*b7*a6*a3-10969*b4*a8*a3-763*b4*a8*a7+336*b5^2*a5+889*b4*b5*a4+6482*a7*b7*a6+637*a7*b5*a9-7224*a7*b5*a5-2016*b8*b7*a10-7406*a7*b4*a5-3192*a7*b7*a10+7637*a7*b4*a9+1218*a7*b8*a5+994*a7*b8*a8+3031*b8*a9*a3-7707*b8*a5*a3+560*b8*b4*a9-3542*b8*b4*a5-3423*b8*b5*a4+1631*b4*b8*a8+770*b4*b9*a8-5509*b4*a4*a7+10339*b4*b5*a8+2590*a7*b8*a4+63*b7*a4^2-679*b7*a4*a8-2562*b8^2*a5+13027*a3*b4*a5-6202*a3*b4*a9-2723*b7*a5*a9-6860*b9^2*a5+6797*b7*a5*a8+6769*b5*a4*a3+16058*b7*b5*a10-7301*b5*a4*a7-5005*b5^2*a4+2492*b5^2*a8-931*b5^2*a9-3997*b4*b9*a4-1561*b5*b9*a5-4494*b5*b8*a5-8134*b5*b4*a5-5040*b4^2*a9-3465*b5*b8*a8+5404*b7*b8*a6+672*b4^2*a5-42*b7*a9*a8-1302*b7*b6*a5-31731*b7*b5*a6+1015*b7*a9*a4-17535*b7*b6*a4+11403*b7*b6*a9-602*b7*b10*a8+980*b5*b9*a8+3563*b7*b6*a8-1988*b7*a9^2-6468*b7*a5^2+224*b7*b10*a5+3899*b5*b9*a9+854*b7*a8^2-3171*a5*b7*a4-182*b8^2*a9-791*b8^2*a8+2429*b9*b5*a4-1211*b5*a9*a3+7350*b5*b8*a9-20384*b5*b4*a9-10262*b9^2*b10+210*b4^2*b6+32942*b4*b5*b6+13307*b8^2*b6+2576*a3*b8*b6-15806*a3*b9*b6-1498*a3*b8*b10+10178*a3*b10*b5+4606*a3*b10*b9-6496*a3*b5*b6+12838*b9^2*b6-1652*a3*b10*b4+42*b4^2*a8-210*b7*b6^2-588*b7*b10^2+9219*b10*b9*b8-9233*b5^2*b10+17836*b5^2*b6-8953*b8^2*b10-21735*b5*b6*b8-8120*b5*b6*b9+8491*b10*b9*b5+84*b7*b6*b10+8428*a3*b4*b6-9198*b8*b9*b6+18739*b4*b10*b9-24871*b4*b10*b5+385*b4*b8*b10-31094*b4*b9*b6-4004*a7*b4*b6+1162*a7*b9*b6+238*a7*b10*b9-4186*a7*b10*b5-826*a7*b8*b10+2072*a7*b8*b6+1036*a7*b10*b4+3794*a7*b5*b6+7*b8^2*a4-581*b4*b8*b6-133*b4*b8*a4+10458*b5*b8*b10+3227*b9*a9*a3-6356*b9*a5*a3-1890*b9*b8*a9+399*b9*b8*a8+1876*b9*b4*a5+15568*b9*b8*a5+9737*b9*b4*a9-10367*b5*a8*a3+4207*b5*a8*a7+17094*b5*a5*a3+7462*b9*a5*a7+3395*b9*a8*a3-3367*b9*a4*a3-3619*b9*a9*a7+2135*b9*a4*a7+1043*b9*a8*a7-2184*b9*b7*a10+2674*b9*b7*a6")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("b10"), Symbol("b9"), Symbol("b8"), Symbol("b7"), Symbol("b6"), Symbol("b5"), Symbol("b4"), Symbol("a10"), Symbol("a9"), Symbol("a8"), Symbol("a7"), Symbol("a6"), Symbol("a5"), Symbol("a4"), Symbol("a3")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3059","KdV",showOutput,isLazard);
} // KdV

void Sys3154Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("2*a*h*i + b*h^2+ 2*c*d*j - c*e*i - c*g*h -d*e*h"),
		SMQP("a*i^2 + 2*b*h*i + 2*c*f*j - c*g*i + d^2*j - d*e*i - d*g*h - e*f*h"),
		SMQP("b*i^2 + 2*d*r*j - d*g*i - e*l*i - f*g*h"),
		SMQP("f*(f*j - g*i)")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d"), Symbol("e"), Symbol("f"), Symbol("g"), Symbol("h"), Symbol("i"), Symbol("j"), Symbol("l"), Symbol("r")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3154","L",showOutput,isLazard);
} // L

//nie := [a, z-1];
void Sys3161Test(bool showOutput, bool isLazard) {
//#original form
//#(E a)(E z)[0<=a, z>=1, y1-2/3*a*(-z^4+z), y2*z^2-1/2*a*(z^4-1)]
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y1-2/3*a*(-z^4+z)"),
		SMQP("y2*z^2-1/2*a*(z^4-1)")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("z"), Symbol("y1"), Symbol("y2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3161","Lafferriere35",showOutput,isLazard);
} // Lafferriere35

//pie := [a];
void Sys3162Test(bool showOutput, bool isLazard) {
//#original form:
//#(E w)(E z)(E a)[a>0, 
//#w^2+z^2-1, 
//#w*((4*a^2-2)*z+2-a^2)+3*a, 
//#(a^2-2)*(w^2-z^2+z)-3*a
//#];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("w^2+z^2-1"),
		SMQP("w*((4*a^2-2)*z+2-a^2)+3*a"),
		SMQP("(a^2-2)*(w^2-z^2+z)-3*a")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("z"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3162","Lafferriere37",showOutput,isLazard);
} // Lafferriere37

void Sys3060Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a+b+c+d+e+f+g+h-1"),
		SMQP("-a^2*k-2*a*b*k-b^2*k-a*c*k-b*c*k-a*d*k-b*d*k-a*e*k-b*e*k-c*e*k-d*e*k-a*f*k-b*f*k-c*f*k-d*f*k+a+b"),
		SMQP("-a^2*l-a*b*l-a*c*l-a*d*l-a*e*l-b*e*l-c*e*l-d*e*l+a^2+2*a*b+b^2+a*e+b*e+a*f+b*f"),
		SMQP("a+c+e+g-m")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d"), Symbol("e"), Symbol("f"), Symbol("g"), Symbol("h"), Symbol("k"), Symbol("l"), Symbol("m")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3060","Lanconelli",showOutput,isLazard);
} // Lanconelli

//params := [k2,k1];
void Sys3062Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("20*k1+k1*y^4-20*k2*x-k2*x*y^4-2*x-4*x*y^4"),
		SMQP("2*x+4*x*y^4-1000*y-50*y^5")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("k2"), Symbol("k1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3062","Laurent-96",showOutput,isLazard);
} // Laurent-96

//params := [k4,a];
void Sys3061Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("600+50*y^2-156*x-13*x*y^2-1200*a*x-4000*a*x*y^2"),
		SMQP("12*a*x+40*a*x*y^2-12*k4*y-k4*y^3")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("k4"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3061","Laurent-96-2",showOutput,isLazard);
} // Laurent-96-2

//params := [m, c];
void Sys3063Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("250*m+250*x*y^2-189*x*m-189*x^2*y^2-25*c*x"),
		SMQP("6*x*m+6*x^2*y^2-100*c*x-75*y*m-75*x*y^3, -3*z")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("y"), Symbol("x"), Symbol("m"), Symbol("c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3063","Laurent-98",showOutput,isLazard);
} // Laurent-98

//params := [a,b,c];
void Sys3064Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(c^2-l^2)*y^2 +(b^2-l^2)*z^2 +(b^2+c^2-a^2-2*l^2)*y*z"),
		SMQP("(a^2-l^2)*z^2 +(c^2-l^2)*x^2 +(c^2+a^2-b^2-2*l^2)*z*x"),
		SMQP("(b^2-l^2)*x^2 +(a^2-l^2)*y^2 +(a^2+b^2-c^2-2*l^2)*x*y")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("l"), Symbol("a"), Symbol("b"), Symbol("c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"30","A",showOutput,isLazard);
} // Lazard-ASCM2001

//params:=[a4, a3, a2, a1];
void Sys3065Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x1*x1*x3+(-x1*x2*x2))*x4+(-x1*x2*x3*x3)+x2^3*x3"),
		SMQP("(-7*a4*x4)+(-3*a3*x3)+a2*x2+5*a1*x1"),
		SMQP("(-4*a4*x2*x4)+(-5*a4*x3*x3)+((-6*a3*x2)+(-5*a2*x1))*x3+2*a2*x2*x2"),
		SMQP("(-2*a4*x1*x4)+(5*a4*x2+2*a3*x1)*x3+a2*x1*x2"),
		SMQP("a4*x1*x3+2*a4*x2*x2+2*a3*x1*x2+a2*x1*x1"),
		SMQP("(7*a4*x3+10*a3*x2+5*a2*x1)*x4+(-2*a3*x3*x3)+4*a2*x2*x3+10*a1*x2*x1")};
	std::vector<Symbol> vars = {Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1"), Symbol("a4"), Symbol("a3"), Symbol("a2"), Symbol("a1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3065","Leykin-1",showOutput,isLazard);
} // Leykin-1

//params := [y];
void Sys3066Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x-110*t^2+495*t^3-1320*t^4+2772*t^5-5082*t^6+7590*t^7-8085*t^8+5555*t^9-2189*t^10+374*t^11"),
		SMQP("y-22*t+110*t^2-330*t^3+1848*t^5-3696*t^6+3300*t^7-1650*t^8+550*t^9-88*t^10-22*t^11")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("x"), Symbol("y")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3066","Lichtblau",showOutput,isLazard);
} // Lichtblau

void Sys3068Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y*z-y*t-x+a"),
		SMQP("z*t-z*x-y+a"),
		SMQP("t*x-t*y-z+a"),
		SMQP("x*y-x*z-t+a")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3068","Liu-Lorenz",showOutput,isLazard);
} // Liu-Lorenz

void Sys3067Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y*z-y*t-x+a"),
		SMQP("z*t-z*x-y+a"),
		SMQP("t*x-t*y-z+a"),
		SMQP("x*y-x*z-t+a")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("z"), Symbol("y"), Symbol("x"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3067","Liu-Lorenz-Li",showOutput,isLazard);
} // Liu-Lorenz-Li

//params := [w,a,b,c];
void Sys3069Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a*b - w*c"),
		SMQP("-x*y + y*u + x*v"),
		SMQP("x*t - t*b - x*c"),
		SMQP("y*z - y*w - z*a"),
		SMQP("-z*v + v*b + z*c - u*c"),
		SMQP("t*u - t*w + v*w - u*a")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("a"), Symbol("b"), Symbol("c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3069","MacLane",showOutput,isLazard);
} // MacLane

//params := [a, b, g]; nie := [g, 1-g]; pie := [a, b];
void Sys3183Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x1 -  ( (x1)^3 )/3 - y1 + g*(x1 - x2)"),
		SMQP("x2 -  ( (x2)^3 )/3 - y2 + g*(x2 - x1)"),
		SMQP("(x1 + a - b*y1)"),
		SMQP("(x2 + a - b*y2)")};
	std::vector<Symbol> vars = {Symbol("y1"), Symbol("y2"), Symbol("x1"), Symbol("x2"), Symbol("a"), Symbol("b"), Symbol("g")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3183","Mehta0",showOutput,isLazard);
} // Mehta0

//params := [a, b, g]; nie := [g, 1-g]; pie := [a, b];
void Sys3184Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1+6*g*x1-x1^3-3*g*x2-3*g*x3-3*y1"),
		SMQP("a+x1-b*y1"),
		SMQP("-3*g*x1+3*x2+6*g*x2-x2^3-3*g*x3-3*y2"),
		SMQP("a+x2-b*y2"),
		SMQP("-3*g*x1-3*g*x2+3*x3+6*g*x3-x3^3-3*y3"),
		SMQP("a+x3-b*y3")};
	std::vector<Symbol> vars = {Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("a"), Symbol("b"), Symbol("g")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3184","Mehta1",showOutput,isLazard);
} // Mehta1

//params := [a, b, g]; nie := [g, 1-g]; pie := [a, b];
void Sys3185Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1+9*g*x1-x1^3-3*g*x2-3*g*x3-3*g*x4-3*y1"),
		SMQP("a+x1-b*y1"),
		SMQP("-3*g*x1+3*x2+9*g*x2-x2^3-3*g*x3-3*g*x4-3*y2"),
		SMQP("a+x2-b*y2"),
		SMQP("-3*g*x1-3*g*x2+3*x3+9*g*x3-x3^3-3*g*x4-3*y3"),
		SMQP("a+x3-b*y3"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3+3*x4+9*g*x4-x4^3-3*y4"),
		SMQP("a+x4-b*y4")};
	std::vector<Symbol> vars = {Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("a"), Symbol("b"), Symbol("g")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3185","Mehta2",showOutput,isLazard);
} // Mehta2

//params := [a, b, g];  nie := [g, 1-g]; pie := [a, b];
void Sys3186Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1+12*g*x1-x1^3-3*g*x2-3*g*x3-3*g*x4-3*g*x5-3*y1"),
		SMQP("a+x1-b*y1"),
		SMQP("-3*g*x1+3*x2+12*g*x2-x2^3-3*g*x3-3*g*x4-3*g*x5-3*y2"),
		SMQP("a+x2-b*y2"),
		SMQP("-3*g*x1-3*g*x2+3*x3+12*g*x3-x3^3-3*g*x4-3*g*x5-3*y3"),
		SMQP("a+x3-b*y3"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3+3*x4+12*g*x4-x4^3-3*g*x5-3*y4"),
		SMQP("a+x4-b*y4"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3-3*g*x4+3*x5+12*g*x5-x5^3-3*y5"),
		SMQP("a+x5-b*y5")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("y1"), Symbol("x2"), Symbol("y2"), Symbol("x3"), Symbol("y3"), Symbol("x4"), Symbol("y4"), Symbol("x5"), Symbol("y5"), Symbol("a"), Symbol("b"), Symbol("g")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3186","Mehta3",showOutput,isLazard);
} // Mehta3

//params := [a, b, g]; nie := [g, 1-g]; pie := [a, b];
void Sys3187Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("3*x1+15*g*x1-x1^3-3*g*x2-3*g*x3-3*g*x4-3*g*x5-3*g*x6-3*y1"),
		SMQP("a+x1-b*y1"),
		SMQP("-3*g*x1+3*x2+15*g*x2-x2^3-3*g*x3-3*g*x4-3*g*x5-3*g*x6-3*y2"),
		SMQP("a+x2-b*y2"),
		SMQP("-3*g*x1-3*g*x2+3*x3+15*g*x3-x3^3-3*g*x4-3*g*x5-3*g*x6-3*y3"),
		SMQP("a+x3-b*y3"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3+3*x4+15*g*x4-x4^3-3*g*x5-3*g*x6-3*y4"),
		SMQP("a+x4-b*y4"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3-3*g*x4+3*x5+15*g*x5-x5^3-3*g*x6-3*y5"),
		SMQP("a+x5-b*y5"),
		SMQP("-3*g*x1-3*g*x2-3*g*x3-3*g*x4-3*g*x5+3*x6+15*g*x6-x6^3-3*y6"),
		SMQP("a+x6-b*y6")};
	std::vector<Symbol> vars = {Symbol("y1"), Symbol("y2"), Symbol("y3"), Symbol("y4"), Symbol("y5"), Symbol("y6"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("x5"), Symbol("x6"), Symbol("a"), Symbol("b"), Symbol("g")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3187","Mehta4",showOutput,isLazard);
} // Mehta4

//params := [a, g, x1, x2, x3, b]; nie := [g, 1-g]; pie := [a, b];
void Sys3188Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("27-81*b+81*b^2-27*b^3-162*b*g+324*b^2*g-162*b^3*g+243*b^2*g^2-243*b^3*g^2-108*t^2+189*b*t^2-90*b^2*t^2+9*b^3*t^2-324*g*t^2+702*b*g*t^2-324*b^2*g*t^2+18*b^3*g*t^2-243*g^2*t^2+729*b*g^2*t^2-243*b^2*g^2*t^2+36*t^4-27*b*t^4+3*b^2*t^4+108*g*t^4-54*b*g*t^4+81*g^2*t^4-t^6+27*b*x1^2-54*b^2*x1^2+27*b^3*x1^2-108*b^2*g*x1^2+108*b^3*g*x1^2+81*b^3*g^2*x1^2+54*t^2*x1^2-117*b*t^2*x1^2+54*b^2*t^2*x1^2-3*b^3*t^2*x1^2+108*g*t^2*x1^2-324*b*g*t^2*x1^2+108*b^2*g*t^2*x1^2-243*b*g^2*t^2*x1^2-18*t^4*x1^2+9*b*t^4*x1^2-36*g*t^4*x1^2+27*b*x2^2-54*b^2*x2^2+27*b^3*x2^2-108*b^2*g*x2^2+108*b^3*g*x2^2+81*b^3*g^2*x2^2+54*t^2*x2^2-117*b*t^2*x2^2+54*b^2*t^2*x2^2-3*b^3*t^2*x2^2+108*g*t^2*x2^2-324*b*g*t^2*x2^2+108*b^2*g*t^2*x2^2-243*b*g^2*t^2*x2^2-18*t^4*x2^2+9*b*t^4*x2^2-36*g*t^4*x2^2+27*b^2*x1^2*x2^2-27*b^3*x1^2*x2^2-54*b^3*g*x1^2*x2^2-27*t^2*x1^2*x2^2+81*b*t^2*x1^2*x2^2-27*b^2*t^2*x1^2*x2^2+162*b*g*t^2*x1^2*x2^2+9*t^4*x1^2*x2^2+27*b*x3^2-54*b^2*x3^2+27*b^3*x3^2-108*b^2*g*x3^2+108*b^3*g*x3^2+81*b^3*g^2*x3^2+54*t^2*x3^2-117*b*t^2*x3^2+54*b^2*t^2*x3^2-3*b^3*t^2*x3^2+108*g*t^2*x3^2-324*b*g*t^2*x3^2+108*b^2*g*t^2*x3^2-243*b*g^2*t^2*x3^2-18*t^4*x3^2+9*b*t^4*x3^2-36*g*t^4*x3^2+27*b^2*x1^2*x3^2-27*b^3*x1^2*x3^2-54*b^3*g*x1^2*x3^2-27*t^2*x1^2*x3^2+81*b*t^2*x1^2*x3^2-27*b^2*t^2*x1^2*x3^2+162*b*g*t^2*x1^2*x3^2+9*t^4*x1^2*x3^2+27*b^2*x2^2*x3^2-27*b^3*x2^2*x3^2-54*b^3*g*x2^2*x3^2-27*t^2*x2^2*x3^2+81*b*t^2*x2^2*x3^2-27*b^2*t^2*x2^2*x3^2+162*b*g*t^2*x2^2*x3^2+9*t^4*x2^2*x3^2+27*b^3*x1^2*x2^2*x3^2-81*b*t^2*x1^2*x2^2*x3^2"),
		SMQP("-81*t+189*b*t-135*b^2*t+27*b^3*t-162*g*t+648*b*g*t-594*b^2*g*t+108*b^3*g*t+486*b*g^2*t-729*b^2*g^2*t+81*b^3*g^2*t+81*t^3-99*b*t^3+27*b^2*t^3-b^3*t^3+270*g*t^3-324*b*g*t^3+54*b^2*g*t^3+243*g^2*t^3-243*b*g^2*t^3-9*t^5+3*b*t^5-18*g*t^5+27*t*x1^2-108*b*t*x1^2+99*b^2*t*x1^2-18*b^3*t*x1^2-216*b*g*t*x1^2+324*b^2*g*t*x1^2-36*b^3*g*t*x1^2+243*b^2*g^2*t*x1^2-45*t^3*x1^2+54*b*t^3*x1^2-9*b^2*t^3*x1^2-108*g*t^3*x1^2+108*b*g*t^3*x1^2-81*g^2*t^3*x1^2+3*t^5*x1^2+27*t*x2^2-108*b*t*x2^2+99*b^2*t*x2^2-18*b^3*t*x2^2-216*b*g*t*x2^2+324*b^2*g*t*x2^2-36*b^3*g*t*x2^2+243*b^2*g^2*t*x2^2-45*t^3*x2^2+54*b*t^3*x2^2-9*b^2*t^3*x2^2-108*g*t^3*x2^2+108*b*g*t^3*x2^2-81*g^2*t^3*x2^2+3*t^5*x2^2+54*b*t*x1^2*x2^2-81*b^2*t*x1^2*x2^2+9*b^3*t*x1^2*x2^2-162*b^2*g*t*x1^2*x2^2+27*t^3*x1^2*x2^2-27*b*t^3*x1^2*x2^2+54*g*t^3*x1^2*x2^2+27*t*x3^2-108*b*t*x3^2+99*b^2*t*x3^2-18*b^3*t*x3^2-216*b*g*t*x3^2+324*b^2*g*t*x3^2-36*b^3*g*t*x3^2+243*b^2*g^2*t*x3^2-45*t^3*x3^2+54*b*t^3*x3^2-9*b^2*t^3*x3^2-108*g*t^3*x3^2+108*b*g*t^3*x3^2-81*g^2*t^3*x3^2+3*t^5*x3^2+54*b*t*x1^2*x3^2-81*b^2*t*x1^2*x3^2+9*b^3*t*x1^2*x3^2-162*b^2*g*t*x1^2*x3^2+27*t^3*x1^2*x3^2-27*b*t^3*x1^2*x3^2+54*g*t^3*x1^2*x3^2+54*b*t*x2^2*x3^2-81*b^2*t*x2^2*x3^2+9*b^3*t*x2^2*x3^2-162*b^2*g*t*x2^2*x3^2+27*t^3*x2^2*x3^2-27*b*t^3*x2^2*x3^2+54*g*t^3*x2^2*x3^2+81*b^2*t*x1^2*x2^2*x3^2-27*t^3*x1^2*x2^2*x3^2")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("a"), Symbol("g"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3188","Mehta5",showOutput,isLazard);
} // Mehta5

//params := [a, x1, x2, x3, x4, g, b]; nie := [g, 1-g]; pie := [a, b];
void Sys3189Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("81-324*b+486*b^2-324*b^3+81*b^4-972*b*g+2916*b^2*g-2916*b^3*g+972*b^4*g+3888*b^2*g^2-7776*b^3*g^2+3888*b^4*g^2-5184*b^3*g^3+5184*b^4*g^3-594*t^2+1620*b*t^2-1512*b^2*t^2+540*b^3*t^2-54*b^4*t^2-2916*g*t^2+10692*b*g*t^2-11664*b^2*g*t^2+4212*b^3*g*t^2-324*b^4*g*t^2-3888*g^2*t^2+23328*b*g^2*t^2-31104*b^2*g^2*t^2+10368*b^3*g^2*t^2-432*b^4*g^2*t^2+15552*b*g^3*t^2-31104*b^2*g^3*t^2+6912*b^3*g^3*t^2+459*t^4-756*b*t^4+360*b^2*t^4-48*b^3*t^4+b^4*t^4+2916*g*t^4-4860*b*g*t^4+1944*b^2*g*t^4-144*b^3*g*t^4+6480*g^2*t^4-10368*b*g^2*t^4+2592*b^2*g^2*t^4+5184*g^3*t^4-6912*b*g^3*t^4-66*t^6+48*b*t^6-6*b^2*t^6-324*g*t^6+144*b*g*t^6-432*g^2*t^6+t^8+81*b*x1^2-243*b^2*x1^2+243*b^3*x1^2-81*b^4*x1^2-729*b^2*g*x1^2+1458*b^3*g*x1^2-729*b^4*g*x1^2+1944*b^3*g^2*x1^2-1944*b^4*g^2*x1^2-1296*b^4*g^3*x1^2+243*t^2*x1^2-891*b*t^2*x1^2+972*b^2*t^2*x1^2-351*b^3*t^2*x1^2+27*b^4*t^2*x1^2+729*g*t^2*x1^2-4374*b*g*t^2*x1^2+5832*b^2*g*t^2*x1^2-1944*b^3*g*t^2*x1^2+81*b^4*g*t^2*x1^2-5832*b*g^2*t^2*x1^2+11664*b^2*g^2*t^2*x1^2-2592*b^3*g^2*t^2*x1^2+7776*b^2*g^3*t^2*x1^2-243*t^4*x1^2+405*b*t^4*x1^2-162*b^2*t^4*x1^2+12*b^3*t^4*x1^2-1215*g*t^4*x1^2+1944*b*g*t^4*x1^2-486*b^2*g*t^4*x1^2-1944*g^2*t^4*x1^2+2592*b*g^2*t^4*x1^2-1296*g^3*t^4*x1^2+27*t^6*x1^2-12*b*t^6*x1^2+81*g*t^6*x1^2+81*b*x2^2-243*b^2*x2^2+243*b^3*x2^2-81*b^4*x2^2-729*b^2*g*x2^2+1458*b^3*g*x2^2-729*b^4*g*x2^2+1944*b^3*g^2*x2^2-1944*b^4*g^2*x2^2-1296*b^4*g^3*x2^2+243*t^2*x2^2-891*b*t^2*x2^2+972*b^2*t^2*x2^2-351*b^3*t^2*x2^2+27*b^4*t^2*x2^2+729*g*t^2*x2^2-4374*b*g*t^2*x2^2+5832*b^2*g*t^2*x2^2-1944*b^3*g*t^2*x2^2+81*b^4*g*t^2*x2^2-5832*b*g^2*t^2*x2^2+11664*b^2*g^2*t^2*x2^2-2592*b^3*g^2*t^2*x2^2+7776*b^2*g^3*t^2*x2^2-243*t^4*x2^2+405*b*t^4*x2^2-162*b^2*t^4*x2^2+12*b^3*t^4*x2^2-1215*g*t^4*x2^2+1944*b*g*t^4*x2^2-486*b^2*g*t^4*x2^2-1944*g^2*t^4*x2^2+2592*b*g^2*t^4*x2^2-1296*g^3*t^4*x2^2+27*t^6*x2^2-12*b*t^6*x2^2+81*g*t^6*x2^2+81*b^2*x1^2*x2^2-162*b^3*x1^2*x2^2+81*b^4*x1^2*x2^2-486*b^3*g*x1^2*x2^2+486*b^4*g*x1^2*x2^2+648*b^4*g^2*x1^2*x2^2-81*t^2*x1^2*x2^2+486*b*t^2*x1^2*x2^2-648*b^2*t^2*x1^2*x2^2+216*b^3*t^2*x1^2*x2^2-9*b^4*t^2*x1^2*x2^2+1458*b*g*t^2*x1^2*x2^2-2916*b^2*g*t^2*x1^2*x2^2+648*b^3*g*t^2*x1^2*x2^2-3888*b^2*g^2*t^2*x1^2*x2^2+135*t^4*x1^2*x2^2-216*b*t^4*x1^2*x2^2+54*b^2*t^4*x1^2*x2^2+486*g*t^4*x1^2*x2^2-648*b*g*t^4*x1^2*x2^2+648*g^2*t^4*x1^2*x2^2-9*t^6*x1^2*x2^2+81*b*x3^2-243*b^2*x3^2+243*b^3*x3^2-81*b^4*x3^2-729*b^2*g*x3^2+1458*b^3*g*x3^2-729*b^4*g*x3^2+1944*b^3*g^2*x3^2-1944*b^4*g^2*x3^2-1296*b^4*g^3*x3^2+243*t^2*x3^2-891*b*t^2*x3^2+972*b^2*t^2*x3^2-351*b^3*t^2*x3^2+27*b^4*t^2*x3^2+729*g*t^2*x3^2-4374*b*g*t^2*x3^2+5832*b^2*g*t^2*x3^2-1944*b^3*g*t^2*x3^2+81*b^4*g*t^2*x3^2-5832*b*g^2*t^2*x3^2+11664*b^2*g^2*t^2*x3^2-2592*b^3*g^2*t^2*x3^2+7776*b^2*g^3*t^2*x3^2-243*t^4*x3^2+405*b*t^4*x3^2-162*b^2*t^4*x3^2+12*b^3*t^4*x3^2-1215*g*t^4*x3^2+1944*b*g*t^4*x3^2-486*b^2*g*t^4*x3^2-1944*g^2*t^4*x3^2+2592*b*g^2*t^4*x3^2-1296*g^3*t^4*x3^2+27*t^6*x3^2-12*b*t^6*x3^2+81*g*t^6*x3^2+81*b^2*x1^2*x3^2-162*b^3*x1^2*x3^2+81*b^4*x1^2*x3^2-486*b^3*g*x1^2*x3^2+486*b^4*g*x1^2*x3^2+648*b^4*g^2*x1^2*x3^2-81*t^2*x1^2*x3^2+486*b*t^2*x1^2*x3^2-648*b^2*t^2*x1^2*x3^2+216*b^3*t^2*x1^2*x3^2-9*b^4*t^2*x1^2*x3^2+1458*b*g*t^2*x1^2*x3^2-2916*b^2*g*t^2*x1^2*x3^2+648*b^3*g*t^2*x1^2*x3^2-3888*b^2*g^2*t^2*x1^2*x3^2+135*t^4*x1^2*x3^2-216*b*t^4*x1^2*x3^2+54*b^2*t^4*x1^2*x3^2+486*g*t^4*x1^2*x3^2-648*b*g*t^4*x1^2*x3^2+648*g^2*t^4*x1^2*x3^2-9*t^6*x1^2*x3^2+81*b^2*x2^2*x3^2-162*b^3*x2^2*x3^2+81*b^4*x2^2*x3^2-486*b^3*g*x2^2*x3^2+486*b^4*g*x2^2*x3^2+648*b^4*g^2*x2^2*x3^2-81*t^2*x2^2*x3^2+486*b*t^2*x2^2*x3^2-648*b^2*t^2*x2^2*x3^2+216*b^3*t^2*x2^2*x3^2-9*b^4*t^2*x2^2*x3^2+1458*b*g*t^2*x2^2*x3^2-2916*b^2*g*t^2*x2^2*x3^2+648*b^3*g*t^2*x2^2*x3^2-3888*b^2*g^2*t^2*x2^2*x3^2+135*t^4*x2^2*x3^2-216*b*t^4*x2^2*x3^2+54*b^2*t^4*x2^2*x3^2+486*g*t^4*x2^2*x3^2-648*b*g*t^4*x2^2*x3^2+648*g^2*t^4*x2^2*x3^2-9*t^6*x2^2*x3^2+81*b^3*x1^2*x2^2*x3^2-81*b^4*x1^2*x2^2*x3^2-243*b^4*g*x1^2*x2^2*x3^2-243*b*t^2*x1^2*x2^2*x3^2+486*b^2*t^2*x1^2*x2^2*x3^2-108*b^3*t^2*x1^2*x2^2*x3^2+1458*b^2*g*t^2*x1^2*x2^2*x3^2-81*t^4*x1^2*x2^2*x3^2+108*b*t^4*x1^2*x2^2*x3^2-243*g*t^4*x1^2*x2^2*x3^2+81*b*x4^2-243*b^2*x4^2+243*b^3*x4^2-81*b^4*x4^2-729*b^2*g*x4^2+1458*b^3*g*x4^2-729*b^4*g*x4^2+1944*b^3*g^2*x4^2-1944*b^4*g^2*x4^2-1296*b^4*g^3*x4^2+243*t^2*x4^2-891*b*t^2*x4^2+972*b^2*t^2*x4^2-351*b^3*t^2*x4^2+27*b^4*t^2*x4^2+729*g*t^2*x4^2-4374*b*g*t^2*x4^2+5832*b^2*g*t^2*x4^2-1944*b^3*g*t^2*x4^2+81*b^4*g*t^2*x4^2-5832*b*g^2*t^2*x4^2+11664*b^2*g^2*t^2*x4^2-2592*b^3*g^2*t^2*x4^2+7776*b^2*g^3*t^2*x4^2-243*t^4*x4^2+405*b*t^4*x4^2-162*b^2*t^4*x4^2+12*b^3*t^4*x4^2-1215*g*t^4*x4^2+1944*b*g*t^4*x4^2-486*b^2*g*t^4*x4^2-1944*g^2*t^4*x4^2+2592*b*g^2*t^4*x4^2-1296*g^3*t^4*x4^2+27*t^6*x4^2-12*b*t^6*x4^2+81*g*t^6*x4^2+81*b^2*x1^2*x4^2-162*b^3*x1^2*x4^2+81*b^4*x1^2*x4^2-486*b^3*g*x1^2*x4^2+486*b^4*g*x1^2*x4^2+648*b^4*g^2*x1^2*x4^2-81*t^2*x1^2*x4^2+486*b*t^2*x1^2*x4^2-648*b^2*t^2*x1^2*x4^2+216*b^3*t^2*x1^2*x4^2-9*b^4*t^2*x1^2*x4^2+1458*b*g*t^2*x1^2*x4^2-2916*b^2*g*t^2*x1^2*x4^2+648*b^3*g*t^2*x1^2*x4^2-3888*b^2*g^2*t^2*x1^2*x4^2+135*t^4*x1^2*x4^2-216*b*t^4*x1^2*x4^2+54*b^2*t^4*x1^2*x4^2+486*g*t^4*x1^2*x4^2-648*b*g*t^4*x1^2*x4^2+648*g^2*t^4*x1^2*x4^2-9*t^6*x1^2*x4^2+81*b^2*x2^2*x4^2-162*b^3*x2^2*x4^2+81*b^4*x2^2*x4^2-486*b^3*g*x2^2*x4^2+486*b^4*g*x2^2*x4^2+648*b^4*g^2*x2^2*x4^2-81*t^2*x2^2*x4^2+486*b*t^2*x2^2*x4^2-648*b^2*t^2*x2^2*x4^2+216*b^3*t^2*x2^2*x4^2-9*b^4*t^2*x2^2*x4^2+1458*b*g*t^2*x2^2*x4^2-2916*b^2*g*t^2*x2^2*x4^2+648*b^3*g*t^2*x2^2*x4^2-3888*b^2*g^2*t^2*x2^2*x4^2+135*t^4*x2^2*x4^2-216*b*t^4*x2^2*x4^2+54*b^2*t^4*x2^2*x4^2+486*g*t^4*x2^2*x4^2-648*b*g*t^4*x2^2*x4^2+648*g^2*t^4*x2^2*x4^2-9*t^6*x2^2*x4^2+81*b^3*x1^2*x2^2*x4^2-81*b^4*x1^2*x2^2*x4^2-243*b^4*g*x1^2*x2^2*x4^2-243*b*t^2*x1^2*x2^2*x4^2+486*b^2*t^2*x1^2*x2^2*x4^2-108*b^3*t^2*x1^2*x2^2*x4^2+1458*b^2*g*t^2*x1^2*x2^2*x4^2-81*t^4*x1^2*x2^2*x4^2+108*b*t^4*x1^2*x2^2*x4^2-243*g*t^4*x1^2*x2^2*x4^2+81*b^2*x3^2*x4^2-162*b^3*x3^2*x4^2+81*b^4*x3^2*x4^2-486*b^3*g*x3^2*x4^2+486*b^4*g*x3^2*x4^2+648*b^4*g^2*x3^2*x4^2-81*t^2*x3^2*x4^2+486*b*t^2*x3^2*x4^2-648*b^2*t^2*x3^2*x4^2+216*b^3*t^2*x3^2*x4^2-9*b^4*t^2*x3^2*x4^2+1458*b*g*t^2*x3^2*x4^2-2916*b^2*g*t^2*x3^2*x4^2+648*b^3*g*t^2*x3^2*x4^2-3888*b^2*g^2*t^2*x3^2*x4^2+135*t^4*x3^2*x4^2-216*b*t^4*x3^2*x4^2+54*b^2*t^4*x3^2*x4^2+486*g*t^4*x3^2*x4^2-648*b*g*t^4*x3^2*x4^2+648*g^2*t^4*x3^2*x4^2-9*t^6*x3^2*x4^2+81*b^3*x1^2*x3^2*x4^2-81*b^4*x1^2*x3^2*x4^2-243*b^4*g*x1^2*x3^2*x4^2-243*b*t^2*x1^2*x3^2*x4^2+486*b^2*t^2*x1^2*x3^2*x4^2-108*b^3*t^2*x1^2*x3^2*x4^2+1458*b^2*g*t^2*x1^2*x3^2*x4^2-81*t^4*x1^2*x3^2*x4^2+108*b*t^4*x1^2*x3^2*x4^2-243*g*t^4*x1^2*x3^2*x4^2+81*b^3*x2^2*x3^2*x4^2-81*b^4*x2^2*x3^2*x4^2-243*b^4*g*x2^2*x3^2*x4^2-243*b*t^2*x2^2*x3^2*x4^2+486*b^2*t^2*x2^2*x3^2*x4^2-108*b^3*t^2*x2^2*x3^2*x4^2+1458*b^2*g*t^2*x2^2*x3^2*x4^2-81*t^4*x2^2*x3^2*x4^2+108*b*t^4*x2^2*x3^2*x4^2-243*g*t^4*x2^2*x3^2*x4^2+81*b^4*x1^2*x2^2*x3^2*x4^2-486*b^2*t^2*x1^2*x2^2*x3^2*x4^2+81*t^4*x1^2*x2^2*x3^2*x4^2"),
		SMQP("-324*t+1080*b*t-1296*b^2*t+648*b^3*t-108*b^4*t-972*g*t+5832*b*g*t-9720*b^2*g*t+5832*b^3*g*t-972*b^4*g*t+7776*b*g^2*t-23328*b^2*g^2*t+18144*b^3*g^2*t-2592*b^4*g^2*t-15552*b^2*g^3*t+20736*b^3*g^3*t-1728*b^4*g^3*t+648*t^3-1404*b*t^3+972*b^2*t^3-228*b^3*t^3+12*b^4*t^3+3888*g*t^3-9720*b*g*t^3+6804*b^2*g*t^3-1296*b^3*g*t^3+36*b^4*g*t^3+7776*g^2*t^3-23328*b*g^2*t^3+15552*b^2*g^2*t^3-1728*b^3*g^2*t^3+5184*g^3*t^3-20736*b*g^3*t^3+10368*b^2*g^3*t^3-216*t^5+252*b*t^5-72*b^2*t^5+4*b^3*t^5-1296*g*t^5+1296*b*g*t^5-216*b^2*g*t^5-2592*g^2*t^5+1728*b*g^2*t^5-1728*g^3*t^5+12*t^7-4*b*t^7+36*g*t^7+81*t*x1^2-486*b*t*x1^2+810*b^2*t*x1^2-486*b^3*t*x1^2+81*b^4*t*x1^2-1458*b*g*t*x1^2+4374*b^2*g*t*x1^2-3402*b^3*g*t*x1^2+486*b^4*g*t*x1^2+5832*b^2*g^2*t*x1^2-7776*b^3*g^2*t*x1^2+648*b^4*g^2*t*x1^2-5184*b^3*g^3*t*x1^2-324*t^3*x1^2+810*b*t^3*x1^2-567*b^2*t^3*x1^2+108*b^3*t^3*x1^2-3*b^4*t^3*x1^2-1458*g*t^3*x1^2+4374*b*g*t^3*x1^2-2916*b^2*g*t^3*x1^2+324*b^3*g*t^3*x1^2-1944*g^2*t^3*x1^2+7776*b*g^2*t^3*x1^2-3888*b^2*g^2*t^3*x1^2+5184*b*g^3*t^3*x1^2+108*t^5*x1^2-108*b*t^5*x1^2+18*b^2*t^5*x1^2+486*g*t^5*x1^2-324*b*g*t^5*x1^2+648*g^2*t^5*x1^2-3*t^7*x1^2+81*t*x2^2-486*b*t*x2^2+810*b^2*t*x2^2-486*b^3*t*x2^2+81*b^4*t*x2^2-1458*b*g*t*x2^2+4374*b^2*g*t*x2^2-3402*b^3*g*t*x2^2+486*b^4*g*t*x2^2+5832*b^2*g^2*t*x2^2-7776*b^3*g^2*t*x2^2+648*b^4*g^2*t*x2^2-5184*b^3*g^3*t*x2^2-324*t^3*x2^2+810*b*t^3*x2^2-567*b^2*t^3*x2^2+108*b^3*t^3*x2^2-3*b^4*t^3*x2^2-1458*g*t^3*x2^2+4374*b*g*t^3*x2^2-2916*b^2*g*t^3*x2^2+324*b^3*g*t^3*x2^2-1944*g^2*t^3*x2^2+7776*b*g^2*t^3*x2^2-3888*b^2*g^2*t^3*x2^2+5184*b*g^3*t^3*x2^2+108*t^5*x2^2-108*b*t^5*x2^2+18*b^2*t^5*x2^2+486*g*t^5*x2^2-324*b*g*t^5*x2^2+648*g^2*t^5*x2^2-3*t^7*x2^2+162*b*t*x1^2*x2^2-486*b^2*t*x1^2*x2^2+378*b^3*t*x1^2*x2^2-54*b^4*t*x1^2*x2^2-1458*b^2*g*t*x1^2*x2^2+1944*b^3*g*t*x1^2*x2^2-162*b^4*g*t*x1^2*x2^2+2592*b^3*g^2*t*x1^2*x2^2+162*t^3*x1^2*x2^2-486*b*t^3*x1^2*x2^2+324*b^2*t^3*x1^2*x2^2-36*b^3*t^3*x1^2*x2^2+486*g*t^3*x1^2*x2^2-1944*b*g*t^3*x1^2*x2^2+972*b^2*g*t^3*x1^2*x2^2-2592*b*g^2*t^3*x1^2*x2^2-54*t^5*x1^2*x2^2+36*b*t^5*x1^2*x2^2-162*g*t^5*x1^2*x2^2+81*t*x3^2-486*b*t*x3^2+810*b^2*t*x3^2-486*b^3*t*x3^2+81*b^4*t*x3^2-1458*b*g*t*x3^2+4374*b^2*g*t*x3^2-3402*b^3*g*t*x3^2+486*b^4*g*t*x3^2+5832*b^2*g^2*t*x3^2-7776*b^3*g^2*t*x3^2+648*b^4*g^2*t*x3^2-5184*b^3*g^3*t*x3^2-324*t^3*x3^2+810*b*t^3*x3^2-567*b^2*t^3*x3^2+108*b^3*t^3*x3^2-3*b^4*t^3*x3^2-1458*g*t^3*x3^2+4374*b*g*t^3*x3^2-2916*b^2*g*t^3*x3^2+324*b^3*g*t^3*x3^2-1944*g^2*t^3*x3^2+7776*b*g^2*t^3*x3^2-3888*b^2*g^2*t^3*x3^2+5184*b*g^3*t^3*x3^2+108*t^5*x3^2-108*b*t^5*x3^2+18*b^2*t^5*x3^2+486*g*t^5*x3^2-324*b*g*t^5*x3^2+648*g^2*t^5*x3^2-3*t^7*x3^2+162*b*t*x1^2*x3^2-486*b^2*t*x1^2*x3^2+378*b^3*t*x1^2*x3^2-54*b^4*t*x1^2*x3^2-1458*b^2*g*t*x1^2*x3^2+1944*b^3*g*t*x1^2*x3^2-162*b^4*g*t*x1^2*x3^2+2592*b^3*g^2*t*x1^2*x3^2+162*t^3*x1^2*x3^2-486*b*t^3*x1^2*x3^2+324*b^2*t^3*x1^2*x3^2-36*b^3*t^3*x1^2*x3^2+486*g*t^3*x1^2*x3^2-1944*b*g*t^3*x1^2*x3^2+972*b^2*g*t^3*x1^2*x3^2-2592*b*g^2*t^3*x1^2*x3^2-54*t^5*x1^2*x3^2+36*b*t^5*x1^2*x3^2-162*g*t^5*x1^2*x3^2+162*b*t*x2^2*x3^2-486*b^2*t*x2^2*x3^2+378*b^3*t*x2^2*x3^2-54*b^4*t*x2^2*x3^2-1458*b^2*g*t*x2^2*x3^2+1944*b^3*g*t*x2^2*x3^2-162*b^4*g*t*x2^2*x3^2+2592*b^3*g^2*t*x2^2*x3^2+162*t^3*x2^2*x3^2-486*b*t^3*x2^2*x3^2+324*b^2*t^3*x2^2*x3^2-36*b^3*t^3*x2^2*x3^2+486*g*t^3*x2^2*x3^2-1944*b*g*t^3*x2^2*x3^2+972*b^2*g*t^3*x2^2*x3^2-2592*b*g^2*t^3*x2^2*x3^2-54*t^5*x2^2*x3^2+36*b*t^5*x2^2*x3^2-162*g*t^5*x2^2*x3^2+243*b^2*t*x1^2*x2^2*x3^2-324*b^3*t*x1^2*x2^2*x3^2+27*b^4*t*x1^2*x2^2*x3^2-972*b^3*g*t*x1^2*x2^2*x3^2-81*t^3*x1^2*x2^2*x3^2+324*b*t^3*x1^2*x2^2*x3^2-162*b^2*t^3*x1^2*x2^2*x3^2+972*b*g*t^3*x1^2*x2^2*x3^2+27*t^5*x1^2*x2^2*x3^2+81*t*x4^2-486*b*t*x4^2+810*b^2*t*x4^2-486*b^3*t*x4^2+81*b^4*t*x4^2-1458*b*g*t*x4^2+4374*b^2*g*t*x4^2-3402*b^3*g*t*x4^2+486*b^4*g*t*x4^2+5832*b^2*g^2*t*x4^2-7776*b^3*g^2*t*x4^2+648*b^4*g^2*t*x4^2-5184*b^3*g^3*t*x4^2-324*t^3*x4^2+810*b*t^3*x4^2-567*b^2*t^3*x4^2+108*b^3*t^3*x4^2-3*b^4*t^3*x4^2-1458*g*t^3*x4^2+4374*b*g*t^3*x4^2-2916*b^2*g*t^3*x4^2+324*b^3*g*t^3*x4^2-1944*g^2*t^3*x4^2+7776*b*g^2*t^3*x4^2-3888*b^2*g^2*t^3*x4^2+5184*b*g^3*t^3*x4^2+108*t^5*x4^2-108*b*t^5*x4^2+18*b^2*t^5*x4^2+486*g*t^5*x4^2-324*b*g*t^5*x4^2+648*g^2*t^5*x4^2-3*t^7*x4^2+162*b*t*x1^2*x4^2-486*b^2*t*x1^2*x4^2+378*b^3*t*x1^2*x4^2-54*b^4*t*x1^2*x4^2-1458*b^2*g*t*x1^2*x4^2+1944*b^3*g*t*x1^2*x4^2-162*b^4*g*t*x1^2*x4^2+2592*b^3*g^2*t*x1^2*x4^2+162*t^3*x1^2*x4^2-486*b*t^3*x1^2*x4^2+324*b^2*t^3*x1^2*x4^2-36*b^3*t^3*x1^2*x4^2+486*g*t^3*x1^2*x4^2-1944*b*g*t^3*x1^2*x4^2+972*b^2*g*t^3*x1^2*x4^2-2592*b*g^2*t^3*x1^2*x4^2-54*t^5*x1^2*x4^2+36*b*t^5*x1^2*x4^2-162*g*t^5*x1^2*x4^2+162*b*t*x2^2*x4^2-486*b^2*t*x2^2*x4^2+378*b^3*t*x2^2*x4^2-54*b^4*t*x2^2*x4^2-1458*b^2*g*t*x2^2*x4^2+1944*b^3*g*t*x2^2*x4^2-162*b^4*g*t*x2^2*x4^2+2592*b^3*g^2*t*x2^2*x4^2+162*t^3*x2^2*x4^2-486*b*t^3*x2^2*x4^2+324*b^2*t^3*x2^2*x4^2-36*b^3*t^3*x2^2*x4^2+486*g*t^3*x2^2*x4^2-1944*b*g*t^3*x2^2*x4^2+972*b^2*g*t^3*x2^2*x4^2-2592*b*g^2*t^3*x2^2*x4^2-54*t^5*x2^2*x4^2+36*b*t^5*x2^2*x4^2-162*g*t^5*x2^2*x4^2+243*b^2*t*x1^2*x2^2*x4^2-324*b^3*t*x1^2*x2^2*x4^2+27*b^4*t*x1^2*x2^2*x4^2-972*b^3*g*t*x1^2*x2^2*x4^2-81*t^3*x1^2*x2^2*x4^2+324*b*t^3*x1^2*x2^2*x4^2-162*b^2*t^3*x1^2*x2^2*x4^2+972*b*g*t^3*x1^2*x2^2*x4^2+27*t^5*x1^2*x2^2*x4^2+162*b*t*x3^2*x4^2-486*b^2*t*x3^2*x4^2+378*b^3*t*x3^2*x4^2-54*b^4*t*x3^2*x4^2-1458*b^2*g*t*x3^2*x4^2+1944*b^3*g*t*x3^2*x4^2-162*b^4*g*t*x3^2*x4^2+2592*b^3*g^2*t*x3^2*x4^2+162*t^3*x3^2*x4^2-486*b*t^3*x3^2*x4^2+324*b^2*t^3*x3^2*x4^2-36*b^3*t^3*x3^2*x4^2+486*g*t^3*x3^2*x4^2-1944*b*g*t^3*x3^2*x4^2+972*b^2*g*t^3*x3^2*x4^2-2592*b*g^2*t^3*x3^2*x4^2-54*t^5*x3^2*x4^2+36*b*t^5*x3^2*x4^2-162*g*t^5*x3^2*x4^2+243*b^2*t*x1^2*x3^2*x4^2-324*b^3*t*x1^2*x3^2*x4^2+27*b^4*t*x1^2*x3^2*x4^2-972*b^3*g*t*x1^2*x3^2*x4^2-81*t^3*x1^2*x3^2*x4^2+324*b*t^3*x1^2*x3^2*x4^2-162*b^2*t^3*x1^2*x3^2*x4^2+972*b*g*t^3*x1^2*x3^2*x4^2+27*t^5*x1^2*x3^2*x4^2+243*b^2*t*x2^2*x3^2*x4^2-324*b^3*t*x2^2*x3^2*x4^2+27*b^4*t*x2^2*x3^2*x4^2-972*b^3*g*t*x2^2*x3^2*x4^2-81*t^3*x2^2*x3^2*x4^2+324*b*t^3*x2^2*x3^2*x4^2-162*b^2*t^3*x2^2*x3^2*x4^2+972*b*g*t^3*x2^2*x3^2*x4^2+27*t^5*x2^2*x3^2*x4^2+324*b^3*t*x1^2*x2^2*x3^2*x4^2-324*b*t^3*x1^2*x2^2*x3^2*x4^2")};
	std::vector<Symbol> vars = {Symbol("t"), Symbol("a"), Symbol("x1"), Symbol("x2"), Symbol("x3"), Symbol("x4"), Symbol("g"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3189","Mehta6",showOutput,isLazard);
} // Mehta6

//params:=[X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,X6,Y6,Z6,X1,Y1,Z1];
void Sys3070Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("c1^2+s1^2-1"),
		SMQP("c2^2+s2^2-1"),
		SMQP("c3^2+s3^2-1"),
		SMQP("c1+X2*s1+X3*c2+X4*s2+X5*c3+X6*s3-X1"),
		SMQP("c1+Y2*s1+Y3*c2+Y4*s2+Y5*c3+Y6*s3-Y1"),
		SMQP("c1+Z2*s1+Z3*c2+Z4*s2+Z5*c3+Z6*s3-Z1")};
	std::vector<Symbol> vars = {Symbol("c1"), Symbol("s2"), Symbol("c2"), Symbol("s3"), Symbol("c3"), Symbol("s1"), Symbol("X2"), Symbol("Y2"), Symbol("Z2"), Symbol("X3"), Symbol("Y3"), Symbol("Z3"), Symbol("X4"), Symbol("Y4"), Symbol("Z4"), Symbol("X5"), Symbol("Y5"), Symbol("Z5"), Symbol("X6"), Symbol("Y6"), Symbol("Z6"), Symbol("X1"), Symbol("Y1"), Symbol("Z1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3070","MESFET-3",showOutput,isLazard);
} // MESFET-3

//params:=[ a, b];
void Sys3071Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a - l3* c3 - l2*c1"),
		SMQP("b - l3*s3 - l2* s1"),
		SMQP("c1^2 + s1^2 - 1"),
		SMQP("c3^2 + s3^2 - 1")};
	std::vector<Symbol> vars = {Symbol("c1"), Symbol("c3"), Symbol("l2"), Symbol("l3"), Symbol("s1"), Symbol("s3"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3071","MontesS12",showOutput,isLazard);
} // MontesS12

//params:=[c];
void Sys3072Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("t^3 - c*u*t^2 - u*v^2 - u*w^2"),
		SMQP("t^3 - c*v*t^2 - v*u^2 - v*w^2"),
		SMQP("t^3 - c*w*t^2 - w*u^2 - w*v^2")};
	std::vector<Symbol> vars = {Symbol("u"), Symbol("v"), Symbol("w"), Symbol("t"), Symbol("c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3072","MontesS14",showOutput,isLazard);
} // MontesS14

//params :=[v,w,a,b];
void Sys3073Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-y*z+x*t"),
		SMQP("-y*u+x*v+y-v"),
		SMQP("z^2+t^2-w^2"),
		SMQP("u^2+v^2-a^2-2*u+1"),
		SMQP("z^2+t^2-2*z*u+u^2-2*t*v+v^2-b^2")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3073","Morgenstein",showOutput,isLazard);
} // Morgenstein

//pie := [a, b, c, r, h, a+b-c, b+c-a, c+a-b];
void Sys3163Test(bool showOutput, bool isLazard) {
//#original form:

//#(E s)(E b)(E c)[
//#a^2*h^2-4*s*(s-a)*(s-b)*(s-c), 
//#2*r*h-b*c,
//#2*s-a-b-c, 
//#a>0, b>0, c>0, r>0, h>0,
//#a+b-c>0, b+c-a>0, c+a-b>0
//#
//#]
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a^2*h^2-4*s*(s-a)*(s-b)*(s-c)"),
		SMQP("2*r*h-b*c"),
		SMQP("2*s-a-b-c")};
	std::vector<Symbol> vars = {Symbol("c"), Symbol("b"), Symbol("s"), Symbol("a"), Symbol("r"), Symbol("h")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3163","MPV89",showOutput,isLazard);
} // MPV89

void Sys3074Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2*z+y^2*z-z*a+1"),
		SMQP("x^2*y+y*z^2-y*a+1"),
		SMQP("x*y^2+x*z^2-x*a+1")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3074","Neural",showOutput,isLazard);
} // Neural

//params := [t];
void Sys3193Test(bool showOutput, bool isLazard) {
//##- Non-equidimensional example
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("t^2*x^2 - t*x"),
		SMQP("(t*x - 1)*y"),
		SMQP("t*x*z - y")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("z"), Symbol("x"), Symbol("t")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3193","Non-Equidim",showOutput,isLazard);
} // Non-Equidim

//params := [c];
void Sys3075Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("z*(x^2+y^2-c)+1"),
		SMQP("y*(x^2+z^2-c)+1"),
		SMQP("x*(y^2+z^2-c)+1")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("c")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3075","Noonburg",showOutput,isLazard);
} // Noonburg

//pie := [x, t, u+2, -u+2, v+2, -v+2, w+2, -w+2];
void Sys3165Test(bool showOutput, bool isLazard) {
//#original form:
//#A, B, C are unknowns
//#the others are parameters.
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("A^2+B^2-A*B*u-1"), 
		SMQP("B^2+C^2-B*C*v-t"),
		SMQP("A^2+C^2-A*C*w-x")};
	std::vector<Symbol> vars = {Symbol("A"), Symbol("B"), Symbol("C"), Symbol("x"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3165","P3P",showOutput,isLazard);
} // P3P

//pie := [ x, u+2, -u+2, v+2, -v+2, w+2, -w+2 ];
void Sys3164Test(bool showOutput, bool isLazard) {
//#original form:
//#A, B, C are unknowns, the others are parameters
//#sys := [A^2+B^2-A*B*u-1, 
//#B^2+C^2-B*C*v-1, 
//#A^2+C^2-A*C*w-x, 
//#x>0, u+2>0, -u+2>0,
//#v+2>0, -v+2>0,
//#w+2>0, -w+2>0
//#];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("A^2+B^2-A*B*u-1"),
		SMQP("B^2+C^2-B*C*v-1"),
		SMQP("A^2+C^2-A*C*w-x")};
	std::vector<Symbol> vars = {Symbol("A"), Symbol("B"), Symbol("C"), Symbol("x"), Symbol("u"), Symbol("v"), Symbol("w")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3164","P3P-Isosceles",showOutput,isLazard);
} // P3P-Isosceles

//params := [w,a,b,c,d,e];
void Sys3079Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x*u - u*w - x*a"),
		SMQP("y*t - t*w - y*a"),
		SMQP("x*v - v*b - x*c"),
		SMQP("z*t - t*b - z*c"),
		SMQP("y*v - v*d - y*e"),
		SMQP("z*u - u*d - z*e")};
	std::vector<Symbol> vars = {'x','y','z','t','u','v','w','a','b','c','d','e'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3079","Pappus",showOutput,isLazard);
} // Pappus

//params := [u,v,w,a];
void Sys3080Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x*y + x*z + x*t - u^2"),
		SMQP("x*y + y*z + y*t - v^2"),
		SMQP("x*z + y*z + z*t - w^2"),
		SMQP("x*t + y*t + z*t - a^2")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3080","Pavelle",showOutput,isLazard);
} // Pavelle

void Sys3081Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("-b0*a1 + a0*b1"),
		SMQP("(b1-b0)*a2 + (-a1-a0)*b2"),
		SMQP("17179869184*b1*l1 + (-17179869184*b2*l2) + 40392294970*a2 +      (-32678053073*a1)+172710882616*a0"),
		SMQP("(-17179869184*b0*l1) + (-17179869184*b2*l2) + 40392294970*a2 +     172710882616*a1 + (-32678053073*a0) - 43570737430"),
		SMQP("(8589934592*b1 + (-8589934592*b0))*l2 + 86355441308*a2 +     20196147485*a1 + 20196147485*a0 - 40392294970"),
		SMQP("(-17179869184*a1*l1) + (-17179869184*a2*l2) + 40392294970*b2 +     32678053073*b1 + 172710882616*b0"),
		SMQP("17179869184*a0*l1 + 17179869184*a2*l2 + (-40392294970*b2) +     172710882616*b1 + 32678053073*b0 + 43570737430"),
		SMQP("((-8589934592*a1) + (-8589934592*a0))*l2 + 86355441308*b2 +     (-20196147485*b1) + 20196147485*b0 - 40392294970")};
	std::vector<Symbol> vars = {Symbol("l1"),Symbol("l2"), Symbol("a0"), Symbol("a1"), Symbol("a2"), Symbol("b0"), Symbol("b1"), Symbol("b2")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3081","Pin1",showOutput,isLazard);
} // Pin1

//params := [d];
void Sys3082Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("35*y+35*y*x^2-106*x-6*x^3"),
		SMQP("5*x-y-y*x^2-5*y*z-5*y*z*x^2"),
		SMQP("z*(y-d)")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("y"), Symbol("x"), Symbol("d")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3082","Prey-Predator",showOutput,isLazard);
}  // Prey-Predator

void Sys3166Test(bool showOutput, bool isLazard) {
//#original form:
//#(E x1)(E y1)(E x2)(E y2) 
//# sys := [x1^2 + y1^2 - 1, (x2 - 10)^2 + y2^2 -9, 
//# x-(x1 + x2 )/2, y-(y1 + y2)/2];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x1^2 + y1^2 - 1"),
		SMQP("(x2 - 10)^2 + y2^2 -9"),
		SMQP("x-(x1 + x2 )/2, y-(y1 + y2)/2")};
	std::vector<Symbol> vars = {Symbol("x1"), Symbol("y1"), Symbol("x2"), Symbol("y2"), Symbol("x"), Symbol("y")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3166","Putnam",showOutput,isLazard);
} // Putnam

void Sys3083Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("t+v-a"),
		SMQP("x+y+z+t-u-w-a"),
		SMQP("x*z+y*z+x*t+z*t-u*w-u*a-w*a"),
		SMQP("x*z*t-u*w*a")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("t"), Symbol("u"), Symbol("v"), Symbol("w"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3083","Raksanyi",showOutput,isLazard);
} // Raksanyi

void Sys3084Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("d * x + c * y * z - x * (x^2 + a * y^2 + b * z^2)"),
		SMQP("d * y + c * z * x - x * (y^2 + a * z^2 + b * x^2)"),
		SMQP("d * z + c * x * y - x * (z^2 + a * x^2 + b * y^2)")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("z"), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3084","RationalInterp",showOutput,isLazard);
} // RationalInterp

// params := [a,b];
void Sys3086Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("c_1*c_2 - s_1*s_2 + c_1 - a"),
		SMQP("c_1*s_2 + c_2*s_1 + s_1 - b"),
		SMQP("c_1^2 + s_1^2 -1"),
		SMQP("c_2^2 + s_2^2 -1")};
	std::vector<Symbol> vars = {Symbol("s_1"), Symbol("c_1"), Symbol("s_2"), Symbol("c_2"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3086","RobotPlanoEasy",showOutput,isLazard);
} // RobotPlanoEasy

//params :=[l_2,l_3,a,b];
void Sys3087Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("l_3*c_1*c_2 - l_3*s_1*s_2 + l_2*c_1 - a"),
		SMQP("l_3*c_1*s_2 + l_3*c_2*s_1 + l_2*s_1 - b"),
		SMQP("c_1^2 + s_1^2 -1"),
		SMQP("c_2^2 + s_2^2 -1")};
	std::vector<Symbol> vars = {Symbol("s_1"), Symbol("c_1"), Symbol("s_2"), Symbol("c_2"), Symbol("l_2"), Symbol("l_3"), Symbol("a"), Symbol("b")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3087","RobotPlano",showOutput,isLazard);
} // RobotPlano

//pie := [ F, J, T, s, b1, b2, d, v, r1, r2, q, b1-b2 ];
void Sys3167Test(bool showOutput, bool isLazard) {
//#original form:
//#(E s)(E F)(E J)(E T)[
//#d -d*s - b1*J*s, b1*J*s + b2*J*T-(d + v + r1)*F + (1 - q)*r2*J,
//#v*F - (d + r2 )*J, -d*T + r1*F + q*r2*J - b2*T*J, 
//#F > 0, J > 0, T > 0, s > 0, b1 > 0, b2 > 0, d > 0,
//#v > 0, r1 > 0, r2 > 0, q > 0, b1-b2>0];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("d -d*s - b1*J*s"),
		SMQP("b1*J*s + b2*J*T-(d + v + r1)*F + (1 - q)*r2*J"),
		SMQP("v*F - (d + r2 )*J, -d*T + r1*F + q*r2*J - b2*T*J")};
	std::vector<Symbol> vars = {Symbol("F"), Symbol("J"), Symbol("T"), Symbol("s"), Symbol("b1"), Symbol("b2"), Symbol("d"), Symbol("q"), Symbol("r1"), Symbol("r2"), Symbol("v")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3167","SEIT",showOutput,isLazard);
} // SEIT

//pie := [u+1, v-u, 1-v, r-b, r-1, -r^2+24*r-16];
void Sys3168Test(bool showOutput, bool isLazard) {
//#original form:
//#(E b)(E u)(E v)[u^3+r*u^2+u^2-a*u+r*u+u-b-a+r+1,
//#v^3+r*v^2-v^2-a*v-r*v+v-b+a+r-1,
//#4*u^3+3*r*u^2-2*a*u-b,
//#4*v^3+3*r*v^2-2*a*v-b,
//#u+1>0,v-u>0,1-v>0,r-b>0,r-1>0, 
//#r^2-24*r+16<0
//#];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("u^3+r*u^2+u^2-a*u+r*u+u-b-a+r+1"),
		SMQP("v^3+r*v^2-v^2-a*v-r*v+v-b+a+r-1"),
		SMQP("4*u^3+3*r*u^2-2*a*u-b"),
		SMQP("4*v^3+3*r*v^2-2*a*v-b")};
	std::vector<Symbol> vars = {Symbol("v"), Symbol("u"), Symbol("b"), Symbol("a"), Symbol("r")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3168","Solotareff-4a",showOutput,isLazard);
} // Solotareff-4a

//pie := [u+1, v-u, 1-v, r-b, r-1, -r^2+24*r-16];
void Sys3169Test(bool showOutput, bool isLazard) {
//#original form:
//#(E a)(E u)(E v)[u^3+r*u^2+u^2-a*u+r*u+u-b-a+r+1,
//#v^3+r*v^2-v^2-a*v-r*v+v-b+a+r-1,
//#4*u^3+3*r*u^2-2*a*u-b,
//#4*v^3+3*r*v^2-2*a*v-b,
//#u+1>0,v-u>0,1-v>0,r-b>0,r-1>0, 
//#r^2-24*r+16<0
//#];
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("u^3+r*u^2+u^2-a*u+r*u+u-b-a+r+1"),
		SMQP("v^3+r*v^2-v^2-a*v-r*v+v-b+a+r-1"),
		SMQP("4*u^3+3*r*u^2-2*a*u-b"),
		SMQP("4*v^3+3*r*v^2-2*a*v-b")};
	std::vector<Symbol> vars = {Symbol("v"), Symbol("u"), Symbol("a"), Symbol("b"), Symbol("r")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3169","Solotareff-4b",showOutput,isLazard);
} // Solotareff-4b

void Sys3088Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("u*(u-1)"),
		SMQP("u * v - 1"),
		SMQP("(u-1) * w^2 + w + 1")};
	std::vector<Symbol> vars = {'w','v','u'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3088","Trivial-10",showOutput,isLazard);
} // trivial-10

void Sys3089Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2 + y^2 + z^2 -1"),
		SMQP("x^2 + y^2 -1")};
	std::vector<Symbol> vars = {'x','y','z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3089","Trivial-6",showOutput,isLazard);
} // trivial-6

void Sys3090Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x * (x-1) * (x-2)"),
		SMQP("x  * ((x-1) * (y-1) + (x-2) * y)"),
		SMQP("x * (x-1) * z")};
	std::vector<Symbol> vars = {'x','y','z'};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3090","Trivial-7",showOutput,isLazard);
} // trivial-7

// params := [x1];
void Sys3091Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x2^2 - x1"),
		SMQP("x3^2 -2 * x2 * x3 + x1"),
		SMQP("(x3 - x2) * x4")};
	std::vector<Symbol> vars = {Symbol("x4"),Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3091","Trivial-8",showOutput,isLazard);
} // trivial-8

//params := [x];
void Sys3092Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("v*v-2*y*v+u*u-2*x*u+y*y+x*x-1"),
		SMQP("v*v-u^3"),
		SMQP("-3*u*u*v-2*u*v+2*x*v+3*y*u*u"),
		SMQP("6*u*u*v*w*w-2*v*w-3*u*u*w+1")};
	std::vector<Symbol> vars = {Symbol("w"), Symbol("v"), Symbol("u"), Symbol("y"), Symbol("x")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3092","Vermeer",showOutput,isLazard);
} // Vermeer

//params:=[a11,a21];
void Sys3093Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a20*a11+a21+ a11*a02+3*a03"),
		SMQP("54*a20*a03+9*a20*a11*a02-9*a21*a02-9*a11*a12-18*a30*a11-2*a11^3"),
		SMQP("18*a30*a03-9*a20^2*a03+3*a30*a11*a02+3*a20*a02*a21+3*a20*a12*a11-3*a21*a12-3*a30*a21-2*a11^2*a21"),
		SMQP("3*a30*a21*a02+3*a30*a11*a12+3*a20*a21*a12-18*a20*a30*a03-2*a11*a21^2"),
		SMQP("9*a30*a21*a12- 27*a30^2*a03-2*a21^3")};
	std::vector<Symbol> vars = {Symbol("a12"),Symbol("a02"), Symbol("a03"), Symbol("a20"), Symbol("a30"), Symbol("a11"), Symbol("a21")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3093","Wang168",showOutput,isLazard);
} // Wang168

void Sys3094Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x1^2*x4 -1"),
		SMQP("x3^2*x5 -1"),
		SMQP("x2^2*x6 -1"),
		SMQP("x2^2*a1 +(2*x1^2*x2^2*x3^2*x5 +2*x1^2*x2^2*x3^2*x4)*x6 +x1^2*x2^2*x3^2*x5^2 +2*x1^2*x2^2*x3^2*x4*x5 +x1^2*x2^2*x3^2*x4^2"),
		SMQP("x2*x6*a1 +(x1^2*x2*x3^2*x5 +x1^2*x2*x3^2*x4)*x6^2 +(x1^2*x2*x3^2*x5^2 +2*x1^2*x2*x3^2*x4*x5 +x1^2*x2*x3^2*x4^2)*x6 +x1^2*x2*x3^2*x4*x5^2 +x1^2*x2*x3^2*x4^2*x5"),
		SMQP("x1^2*a2 +x1^2*x2^2*x3^2*x6^2 +(2*x1^2*x2^2*x3^2*x5 +2*x1^2*x2^2*x3^2*x4)*x6 +x1^2*x2^2*x3^2*x5^2 +2*x1^2*x2^2*x3^2*x4*x5"),
		SMQP("x1*x4*a2 +(x1*x2^2*x3^2*x5 +x1*x2^2*x3^2*x4)*x6^2 +(x1*x2^2*x3^2*x5^2 +2*x1*x2^2*x3^2*x4*x5 +x1*x2^2*x3^2*x4^2)*x6 +x1*x2^2*x3^2*x4*x5^2 +x1*x2^2*x3^2*x4^2*x5"),
		SMQP("x3^2*a3 +x1^2*x2^2*x3^2*x6^2 +(2*x1^2*x2^2*x3^2*x5 +2*x1^2*x2^2*x3^2*x4)*x6 +2*x1^2*x2^2*x3^2*x4*x5 +x1^2*x2^2*x3^2*x4^2"),
		SMQP("x3*x5*a3 +(x1^2*x2^2*x3*x5 +x1^2*x2^2*x3*x4)*x6^2 +(x1^2*x2^2*x3*x5^2 +2*x1^2*x2^2*x3*x4*x5 +x1^2*x2^2*x3*x4^2)*x6 +x1^2*x2^2*x3*x4*x5^2 +x1^2*x2^2*x3*x4^2*x5"),
		SMQP("y^2 +(-x1^2*x2^2*x3^2*x5 -x1^2*x2^2*x3^2*x4)*x6^2 +(-x1^2*x2^2*x3^2*x5^2 -2*x1^2*x2^2*x3^2*x4*x5 -x1^2*x2^2*x3^2*x4^2)*x6 -x1^2*x2^2*x3^2*x4*x5^2 -x1^2*x2^2*x3^2*x4^2*x5")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("a3"), Symbol("a2"), Symbol("a1"), Symbol("x6"), Symbol("x5"), Symbol("x4"), Symbol("x3"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3094","Wang-1991c",showOutput,isLazard);
} // Wang-1991c

//params:=[y, x];
void Sys3095Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("(x-u)^2+(y-v)^2-1"),
		SMQP("v^2-u^3"),
		SMQP("2*v*(x-u)+3*u^2*(y-v)"),
		SMQP("2*w*v-1")};
	std::vector<Symbol> vars = {Symbol("w"), Symbol("v"), Symbol("u"), Symbol("y"), Symbol("x")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3095","Wang93",showOutput,isLazard);
} // Wang93

void Sys3096Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x21-x12-x13"),
		SMQP("x22-x11-x13"),
		SMQP("x23-x11-x12"),
		SMQP("x30-x11^3-x12^3-x13^3"),
		SMQP("x21*x22*x23-x10*x30"),
		SMQP("-3*x11^2 * x104 - x103 - x102"),
		SMQP("-3*x12^2 * x104 - x103 - x101"),
		SMQP("-3*x13^2 * x104 - x103 - x101"),
		SMQP("x22*x23*x105 + x101"),
		SMQP("x22*x23*x105 + x101"),
		SMQP("x21*x23*x105 + x102"),
		SMQP("x21*x22*x105 + x103"),
		SMQP("x21*x22*x105 + x103"),
		SMQP("-x10*x105 + x104"),
		SMQP("-x30*x105 + 1")};
	std::vector<Symbol> vars = {Symbol("x105"), Symbol("x104"), Symbol("x103"), Symbol("x102"), Symbol("x101"), Symbol("x30"), Symbol("x23"), Symbol("x22"), Symbol("x21"), Symbol("x13"), Symbol("x12"), Symbol("x11"), Symbol("x10")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3096","Wu-Wang",showOutput,isLazard);
} // Wu-Wang

//params := [r, m, f]; ineqs := [f, m, r];
void Sys3097Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("m-x*y^2-x*(r+f)"),
		SMQP("r*x+x*y^2-f*y")};
	std::vector<Symbol> vars = {Symbol("x"), Symbol("y"), Symbol("r"), Symbol("m"), Symbol("f")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3097","WX-AB05-1",showOutput,isLazard);
} // WX-AB05-1

//params := [wb,hb,a]; pie := [ a+b-c, b+c-a, c+a-b, a, b, c, hb, wb ];
void Sys3170Test(bool showOutput, bool isLazard) {
//#original form
//#c, b, s are unknowns; wb, hb, a are parameters
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("wb^2*(a+c)^2-4*a*c*s*(s-b)"),
		SMQP("b^2*hb^2-4*s*(s-a)*(s-b)*(s-c)"),
		SMQP("2*s-a-b-c")};
	std::vector<Symbol> vars = {Symbol("c"), Symbol("b"), Symbol("s"), Symbol("wb"), Symbol("hb"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3170","Xia",showOutput,isLazard);
} // Xia

//params := [h, r, a]; ineqs := [a, b, c, a+b-c, b+c-a, c+a-b, h, r];
void Sys3098Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a^2*h^2-4*s*(s-a)*(s-b)*(s-c)"),
		SMQP("2*r*h-b*c"),
		SMQP("2*s-a-b-c")};
	std::vector<Symbol> vars = {Symbol("c"), Symbol("b"), Symbol("s"), Symbol("h"), Symbol("r"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3098","Xia-ISSAC07-1",showOutput,isLazard);
} // Xia-ISSAC07-1

//params := [z, y, x];
void Sys3099Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x-u*v"),
		SMQP("y-v"),
		SMQP("z-u^2")};
	std::vector<Symbol> vars = {Symbol("v"), Symbol("u"), Symbol("z"), Symbol("y"), Symbol("x")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3099","Xia-ISSAC07-2",showOutput,isLazard);
} // Xia-ISSAC07-2

//params := [y2, y1]; ineqs := [a, z-1];
void Sys3100Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("y1-2/3*a*(z-z^4)"),
		SMQP("y2*z^2-1/2*a*(z^4-1)")};
	std::vector<Symbol> vars = {Symbol("a"), Symbol("z"), Symbol("y2"), Symbol("y1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3100","Xia-Reachability",showOutput,isLazard);
} // Xia-Reachability

//params := [u2, u1]; ineqs := [x, x^2-y^2];
void Sys3101Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("16*x^2*u2^2-(u1^2+2*u1+1+u2^2)*(1-2*u1+u1^2+u2^2)"),
		SMQP("y^4*u2+(2-2*u2^2-2*u1^2)*y^3+u2*(u1^2-5+u2^2)*y^2")};
	std::vector<Symbol> vars = {Symbol("y"), Symbol("x"), Symbol("u2"), Symbol("u1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3101","XY-5-5-1",showOutput,isLazard);
} // XY-5-5-1

//params := [y2, y1, x2, x1]; ineqs := [y2];
void Sys3102Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("w^2 + z^2 - 1"),
		SMQP("x1*(z^2 - w^2) + 1/3*(3*x2 + 5)*z*w - 2/3*w - y1"),
		SMQP("1/3*(3*x2+5)*(z^2-w^2)-4*x1*z*w-5/3*z-y2")};
	std::vector<Symbol> vars = {Symbol("z"), Symbol("w"), Symbol("y2"), Symbol("y1"), Symbol("x2"), Symbol("x1")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3102","XY-5-7-2",showOutput,isLazard);
} // XY-5-7-2

//params := [ c, b, a];
void Sys3103Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("x^2+y^2-k^2*c^2"),
		SMQP("(1-x)^2+1-k^2*a^2"),
		SMQP("1+(1-y)^2-k^2*b^2")};
	std::vector<Symbol> vars = {Symbol("k"), Symbol("y"), Symbol("x"), Symbol("c"), Symbol("b"), Symbol("a")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3103","Yang011104",showOutput,isLazard);
} // Yang011104

void Sys3104Test(bool showOutput, bool isLazard) {
	std::vector<SparseMultivariateRationalPolynomial> polys = {
		SMQP("a*d*b"),
		SMQP("a*d*c"),
		SMQP("a*d*(a-d)"),
		SMQP("p^2*a - p*a^2 - a*b*c"),
		SMQP("q^2*a - q*a^2 - a*b*c"),
		SMQP("p^2*d - p*d^2 - d*b*c"),
		SMQP("q^2*d - q*d^2 - d*b*c")};
	std::vector<Symbol> vars = {Symbol("p"), Symbol("q"), Symbol("a"), Symbol("b"), Symbol("c"), Symbol("d")};
	for (size_t i = 0; i < polys.size(); ++i)
		polys[i].setRingVariables(vars);
	testSys(polys,vars,"3104","YangBaxterRosso",showOutput,isLazard);
} // YangBaxterRosso
