#ifndef SYSGEN_H
#define SYSGEN_H

#include <bpas.h>
using namespace std;


void testSys(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars, string sysNumber, string sysName, bool showOutput, bool isLazard, bool medianTiming = false);
void testSys(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& algvars, std::vector<Symbol>& transVars, string sysNumber, string sysName, bool showOutput, bool isLazard, bool medianTiming =false);

bool triangularizeValidate(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> vars, bool showOutput, bool isLazard);
bool triangularizeValidate(vector<SMQP> F, RegularChain<RN,SMQP> rc, vector<Symbol> algvars, vector<Symbol> transVars, bool showOutput, bool isLazard);

// Very Easy Tests //
void Sys126Test(bool showOutput, bool isLazard);
void Sys130Test(bool showOutput, bool isLazard);
void Sys132Test(bool showOutput, bool isLazard);
void Sys1000Test(bool showOutput, bool isLazard);
void Sys1001Test(bool showOutput, bool isLazard);
void Sys1002Test(bool showOutput, bool isLazard);
void Sys1003Test(bool showOutput, bool isLazard);
void Sys1004Test(bool showOutput, bool isLazard);
void Sys1005Test(bool showOutput, bool isLazard);
void Sys1006Test(bool showOutput, bool isLazard);
void Sys1007Test(bool showOutput, bool isLazard);
void Sys1008Test(bool showOutput, bool isLazard);
void Sys1009Test(bool showOutput, bool isLazard);
void Sys1010Test(bool showOutput, bool isLazard);
void Sys1011Test(bool showOutput, bool isLazard);
void Sys1012Test(bool showOutput, bool isLazard);
void Sys1013Test(bool showOutput, bool isLazard);
void Sys1014Test(bool showOutput, bool isLazard);
void Sys1015Test(bool showOutput, bool isLazard);
void Sys1016Test(bool showOutput, bool isLazard);
void Sys1221Test(bool showOutput, bool isLazard);
void Sys1255Test(bool showOutput, bool isLazard);
void Sys1289Test(bool showOutput, bool isLazard);
void Sys1302Test(bool showOutput, bool isLazard);
void Sys1303Test(bool showOutput, bool isLazard);
void Sys1304Test(bool showOutput, bool isLazard);
void Sys1305Test(bool showOutput, bool isLazard);
void Sys1364Test(bool showOutput, bool isLazard);
void Sys1366Test(bool showOutput, bool isLazard);
void Sys1373Test(bool showOutput, bool isLazard);
void Sys1397Test(bool showOutput, bool isLazard);
void Sys2852Test(bool showOutput, bool isLazard);
void Sys2985Test(bool showOutput, bool isLazard);

// Easy Tests //
void Sys2915Test(bool showOutput, bool isLazard);
void Sys2920Test(bool showOutput, bool isLazard);
void Sys2926Test(bool showOutput, bool isLazard);
void Sys2931Test(bool showOutput, bool isLazard);
void Sys2934Test(bool showOutput, bool isLazard);
void Sys3105Test(bool showOutput, bool isLazard);
void Sys3110Test(bool showOutput, bool isLazard);
void Sys3111Test(bool showOutput, bool isLazard);
void Sys3112Test(bool showOutput, bool isLazard);
void Sys3113Test(bool showOutput, bool isLazard);
void Sys3117Test(bool showOutput, bool isLazard);
void Sys3120Test(bool showOutput, bool isLazard);
void Sys3145Test(bool showOutput, bool isLazard);
void Sys3149Test(bool showOutput, bool isLazard);

// ISSAC2005 //
void Sys2887Test(bool showOutput, bool isLazard);
void Sys2888Test(bool showOutput, bool isLazard);
void Sys2889Test(bool showOutput, bool isLazard);
void Sys2890aTest(bool showOutput, bool isLazard);
void Sys2890bTest(bool showOutput, bool isLazard);
void Sys2890Test(bool showOutput, bool isLazard);
void Sys2891Test(bool showOutput, bool isLazard);
void Sys2892Test(bool showOutput, bool isLazard);
void Sys2893Test(bool showOutput, bool isLazard);
void Sys2894Test(bool showOutput, bool isLazard);
void Sys2895Test(bool showOutput, bool isLazard);
void Sys2896Test(bool showOutput, bool isLazard);
void Sys2897Test(bool showOutput, bool isLazard);
void Sys2898Test(bool showOutput, bool isLazard);
void Sys2899Test(bool showOutput, bool isLazard);
void Sys2900Test(bool showOutput, bool isLazard);

// ASCM09-zerodim //
void Sys3107Test(bool showOutput, bool isLazard);
void Sys3109Test(bool showOutput, bool isLazard);
void Sys3114Test(bool showOutput, bool isLazard);
void Sys3115Test(bool showOutput, bool isLazard);
void Sys3116Test(bool showOutput, bool isLazard);
void Sys3118Test(bool showOutput, bool isLazard);
void Sys3119Test(bool showOutput, bool isLazard);
void Sys3121Test(bool showOutput, bool isLazard);
void Sys3122Test(bool showOutput, bool isLazard);
void Sys3123Test(bool showOutput, bool isLazard);
void Sys3124Test(bool showOutput, bool isLazard);
void Sys3125Test(bool showOutput, bool isLazard);
void Sys3126Test(bool showOutput, bool isLazard);
void Sys3127Test(bool showOutput, bool isLazard);
void Sys3128Test(bool showOutput, bool isLazard);
void Sys3129Test(bool showOutput, bool isLazard);
void Sys3130Test(bool showOutput, bool isLazard);
void Sys3131Test(bool showOutput, bool isLazard);
void Sys3132Test(bool showOutput, bool isLazard);
void Sys3133Test(bool showOutput, bool isLazard);
void Sys3134Test(bool showOutput, bool isLazard);
void Sys3135Test(bool showOutput, bool isLazard);
void Sys3136Test(bool showOutput, bool isLazard);
void Sys3137Test(bool showOutput, bool isLazard);
void Sys3138Test(bool showOutput, bool isLazard);
void Sys3139Test(bool showOutput, bool isLazard);
void Sys3140Test(bool showOutput, bool isLazard);
void Sys3141Test(bool showOutput, bool isLazard);
void Sys3142Test(bool showOutput, bool isLazard);
void Sys3143Test(bool showOutput, bool isLazard);
void Sys3144Test(bool showOutput, bool isLazard);
void Sys3146Test(bool showOutput, bool isLazard);
void Sys3147Test(bool showOutput, bool isLazard);
void Sys3148Test(bool showOutput, bool isLazard);
void Sys3150Test(bool showOutput, bool isLazard);
void Sys3151Test(bool showOutput, bool isLazard);
void Sys3152Test(bool showOutput, bool isLazard);
void Sys3153Test(bool showOutput, bool isLazard);

// CM09-intersect-posdim //
void Sys3019Test(bool showOutput, bool isLazard); // 4corps-1parameter-homog
void Sys3020Test(bool showOutput, bool isLazard); // 8-3-config-Li
void Sys3021Test(bool showOutput, bool isLazard); // AlKashiSinus
void Sys3023Test(bool showOutput, bool isLazard); // Alonso
void Sys3022Test(bool showOutput, bool isLazard); // Alonso-Li
void Sys3024Test(bool showOutput, bool isLazard); // Bezier
void Sys3025Test(bool showOutput, bool isLazard); // Bjork60
void Sys3027Test(bool showOutput, bool isLazard); // Blood-coagulation
void Sys3026Test(bool showOutput, bool isLazard); // Blood-coagulation-2
void Sys3157Test(bool showOutput, bool isLazard); // BM05-1
void Sys3158Test(bool showOutput, bool isLazard); // BM05-2
void Sys3028Test(bool showOutput, bool isLazard); // Bronstein-Wang
void Sys3029Test(bool showOutput, bool isLazard); // Butcher
void Sys3030Test(bool showOutput, bool isLazard); // CDC2-cyclin
void Sys3031Test(bool showOutput, bool isLazard); // Cheaters-homotopy-easy
void Sys3032Test(bool showOutput, bool isLazard); // Cheaters-homotopy-hard
void Sys3033Test(bool showOutput, bool isLazard); // Chemical
void Sys3034Test(bool showOutput, bool isLazard); // ChildDraw-1
void Sys3035Test(bool showOutput, bool isLazard); // ChildDraw-2
void Sys3036Test(bool showOutput, bool isLazard); // Cinquin-Demongeot-3-1
void Sys3037Test(bool showOutput, bool isLazard); // Cinquin-Demongeot-3-2
void Sys3038Test(bool showOutput, bool isLazard); // Cinquin-Demongeot-3-3
void Sys3039Test(bool showOutput, bool isLazard); // Cinquin-Demongeot-3-4
void Sys3040Test(bool showOutput, bool isLazard); // Cinquin-Demongeot-5-1
void Sys3041Test(bool showOutput, bool isLazard); // Collins-jsc02
void Sys3042Test(bool showOutput, bool isLazard); // Cox-issac07
void Sys3043Test(bool showOutput, bool isLazard); // Cyclic-4
void Sys3179Test(bool showOutput, bool isLazard); // DescartesFolium
void Sys3156Test(bool showOutput, bool isLazard); // DGP29
void Sys3155Test(bool showOutput, bool isLazard); // DGP6
void Sys3045Test(bool showOutput, bool isLazard); // DonatiTraverso
void Sys3044Test(bool showOutput, bool isLazard); // DonatiTraverso-rev
void Sys3181Test(bool showOutput, bool isLazard); // EdgeSquare
void Sys3159Test(bool showOutput, bool isLazard); // Ellipse
void Sys3180Test(bool showOutput, bool isLazard); // EnneperSurface
void Sys3046Test(bool showOutput, bool isLazard); // F-633
void Sys3047Test(bool showOutput, bool isLazard); // F-744
void Sys3048Test(bool showOutput, bool isLazard); // GenLinSyst-2-2
void Sys3049Test(bool showOutput, bool isLazard); // GenLinSyst-3-2
void Sys3050Test(bool showOutput, bool isLazard); // GenLinSyst-3-3
void Sys3051Test(bool showOutput, bool isLazard); // Geometry-Butterfly_3
void Sys3052Test(bool showOutput, bool isLazard); // Gerdt
void Sys3053Test(bool showOutput, bool isLazard); // Gonnet
void Sys3190Test(bool showOutput, bool isLazard); // Haas5
void Sys3054Test(bool showOutput, bool isLazard); // Hairer-2-BGK
void Sys3055Test(bool showOutput, bool isLazard); // Hawes-1
void Sys3056Test(bool showOutput, bool isLazard); // Hereman-2
void Sys3057Test(bool showOutput, bool isLazard); // Hereman-8
void Sys3058Test(bool showOutput, bool isLazard); // Hoffmann
void Sys3160Test(bool showOutput, bool isLazard); // IBVP
void Sys3174Test(bool showOutput, bool isLazard); // Jirstrand22
void Sys3175Test(bool showOutput, bool isLazard); // Jirstrand23
void Sys3176Test(bool showOutput, bool isLazard); // Jirstrand24
void Sys3177Test(bool showOutput, bool isLazard); // Jirstrand41
void Sys3178Test(bool showOutput, bool isLazard); // Jirstrand42
void Sys3191Test(bool showOutput, bool isLazard); // John5
void Sys3192Test(bool showOutput, bool isLazard); // Kahan-ellipsoid
void Sys3059Test(bool showOutput, bool isLazard); // KdV
void Sys3154Test(bool showOutput, bool isLazard); // L
void Sys3161Test(bool showOutput, bool isLazard); // Lafferriere35
void Sys3162Test(bool showOutput, bool isLazard); // Lafferriere37
void Sys3060Test(bool showOutput, bool isLazard); // Lanconelli
void Sys3062Test(bool showOutput, bool isLazard); // laurent-96
void Sys3061Test(bool showOutput, bool isLazard); // laurent-96-2
void Sys3063Test(bool showOutput, bool isLazard); // laurent-98
void Sys3064Test(bool showOutput, bool isLazard); // Lazard-ascm2001
void Sys3065Test(bool showOutput, bool isLazard); // Leykin-1
void Sys3066Test(bool showOutput, bool isLazard); // Lichtblau
void Sys3068Test(bool showOutput, bool isLazard); // Liu-Lorenz
void Sys3067Test(bool showOutput, bool isLazard); // Liu-Lorenz-Li
void Sys3069Test(bool showOutput, bool isLazard); // MacLane
void Sys3183Test(bool showOutput, bool isLazard); // Mehta0
void Sys3184Test(bool showOutput, bool isLazard); // Mehta1
void Sys3185Test(bool showOutput, bool isLazard); // Mehta2
void Sys3186Test(bool showOutput, bool isLazard); // Mehta3
void Sys3187Test(bool showOutput, bool isLazard); // Mehta4
void Sys3188Test(bool showOutput, bool isLazard); // Mehta5
void Sys3189Test(bool showOutput, bool isLazard); // Mehta6
void Sys3070Test(bool showOutput, bool isLazard); // MESFET-3
void Sys3071Test(bool showOutput, bool isLazard); // MontesS12
void Sys3072Test(bool showOutput, bool isLazard); // MontesS14
void Sys3073Test(bool showOutput, bool isLazard); // Morgenstein
void Sys3163Test(bool showOutput, bool isLazard); // MPV89
void Sys3074Test(bool showOutput, bool isLazard); // Neural
void Sys3193Test(bool showOutput, bool isLazard); // Non-equidim
void Sys3075Test(bool showOutput, bool isLazard); // Noonburg
void Sys3165Test(bool showOutput, bool isLazard); // P3P
void Sys3164Test(bool showOutput, bool isLazard); // P3P-isosceles
void Sys3079Test(bool showOutput, bool isLazard); // Pappus
void Sys3080Test(bool showOutput, bool isLazard); // Pavelle
void Sys3081Test(bool showOutput, bool isLazard); // Pin1
void Sys3082Test(bool showOutput, bool isLazard); // Prey-predator
void Sys3166Test(bool showOutput, bool isLazard); // Putnam
void Sys3083Test(bool showOutput, bool isLazard); // Raksanyi
void Sys3084Test(bool showOutput, bool isLazard); // RationalInterp
void Sys3086Test(bool showOutput, bool isLazard); // RobotPlanoEasy
void Sys3087Test(bool showOutput, bool isLazard); // RobotPlano
void Sys3167Test(bool showOutput, bool isLazard); // SEIT
void Sys3168Test(bool showOutput, bool isLazard); // Solotareff-4a
void Sys3169Test(bool showOutput, bool isLazard); // Solotareff-4b
void Sys3088Test(bool showOutput, bool isLazard); // Trivial-10
void Sys3089Test(bool showOutput, bool isLazard); // Trivial-6
void Sys3090Test(bool showOutput, bool isLazard); // Trivial-7
void Sys3091Test(bool showOutput, bool isLazard); // Trivial-8
void Sys3092Test(bool showOutput, bool isLazard); // Vermeer
void Sys3093Test(bool showOutput, bool isLazard); // Wang168
void Sys3094Test(bool showOutput, bool isLazard); // Wang-1991c
void Sys3095Test(bool showOutput, bool isLazard); // Wang93
void Sys3096Test(bool showOutput, bool isLazard); // Wu-Wang
void Sys3097Test(bool showOutput, bool isLazard); // WX-ab05-1
void Sys3170Test(bool showOutput, bool isLazard); // Xia
void Sys3098Test(bool showOutput, bool isLazard); // Xia-ISSAC07-1
void Sys3099Test(bool showOutput, bool isLazard); // Xia-ISSAC07-2
void Sys3100Test(bool showOutput, bool isLazard); // Xia-reachability
void Sys3101Test(bool showOutput, bool isLazard); // XY-5-5-1
void Sys3102Test(bool showOutput, bool isLazard); // XY-5-7-2
void Sys3103Test(bool showOutput, bool isLazard); // Yang011104
void Sys3104Test(bool showOutput, bool isLazard); // YangBaxterRosso

// Uncategorized //
void Sys2250Test(bool showOutput, bool isLazard);


// Sample Get Functions //
void getSampleSys2915(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys2920(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys2926(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys2931(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys2934(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3105(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3110(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3111(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3112(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3113(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3117(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3120(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3145(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);
void getSampleSys3149(std::vector<SparseMultivariateRationalPolynomial>& polys, std::vector<Symbol>& vars);



#endif
