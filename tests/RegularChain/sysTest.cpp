#include <bpas.h>
#include <vector>
#include "sysGen.hpp"

using namespace std;

void triangularizeTests() {

    // Delete old timing data file
	remove("timing.txt");
	remove("rc-timing.txt");
	remove("copy-timing.txt");
	fstream fs;
	fs.open("timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	fs << "System\tprimitivePart\tsquareFreePart\tfactor\tpseudoDivide\tnormalForm\tsubresultantChain\tremoveRedundantComponents\tZDRC\tBPAS\tMaple\tVerification" << std::endl;
	fs.close();
	fs.open("rc-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	fs << "System\tintersect\tregularize\tregularGCD\tintersectFree\tintersectAlgebraic\textend\tcleanChain\tcleanSet\tGCDFreeFact\ttriangularize\tconstructChain\tconstructChains\ttotal\tBPAS\tMaple\tVerification" << std::endl;
	fs.close();
	fs.open("copy-timing.txt",std::fstream::in | std::fstream::out | std::fstream::app);
	fs << "System\ttsTime\trcTime\tzdrcTime\tinternalTotal\texternalTotal" << std::endl;
	fs.close();
	cout << "Starting triangularizeTests..." << endl;
	bool working(false);
	bool fixing(true);
	bool running(false);
	bool testing(true);
	bool laz(true);
	bool kalk(false);

//  	cerr << "\"VeryEasyTests\":" << endl;
	Sys126Test(working,laz);
	Sys130Test(working,laz);
	Sys132Test(working,laz);
	Sys1000Test(working,laz);
	Sys1001Test(working,laz);
	Sys1002Test(working,laz);
	Sys1003Test(working,laz);
	Sys1004Test(working,laz);
	Sys1005Test(working,laz);
	Sys1006Test(working,laz);
	Sys1007Test(working,laz);
	Sys1008Test(working,laz); // MD equivalent-non-equal reduction (ENER) of top polynomial
	Sys1009Test(working,laz);
	Sys1010Test(working,laz);
	Sys1011Test(working,laz);
	Sys1012Test(working,laz);
	Sys1013Test(working,laz);
	Sys1014Test(working,laz);
	Sys1015Test(working,laz);
	Sys1016Test(working,laz);
	Sys1221Test(working,laz); // MD ENER of top polynomial
	Sys1255Test(working,laz);

	Sys1289Test(working,laz); // redundant component in kalk
	Sys1302Test(working,laz);
	Sys1303Test(working,laz);
	Sys1304Test(working,laz);
	Sys1305Test(working,laz);
	Sys1364Test(working,laz);
	Sys1366Test(working,laz);
	Sys1373Test(working,laz);
	Sys1397Test(working,laz);
	Sys2852Test(working,laz);
	Sys2985Test(working,laz);

//  	cerr << "\"EasyTests\":" << endl;
 	Sys2915Test(working,laz); // Hairer-1: MD different but equivalent components, EXP SIZES NOT DEFINED FOR NVAR 0
// Sys2920Test(testing,laz); // Butcher: not clear if working or not because impeded by very slow computations (GCDFreeFactorization)
	Sys2926Test(working,laz); // Gerdt: MD ENER and factorization not done, repeated with differerent variable order in CM09-intersect-posdim (check that it's the same system!)
	// Sys2931Test(working,laz); // Pappus: pass, VERY slow to test (BPAS: 8s, Maple Test Total: 740s), repeated with differerent variable order in CM09-intersect-posdim, (~2s) kalk
	Sys2934Test(working,laz); // robotPlanoEasy: repeated with differerent variable order in CM09-intersect-posdim
	Sys3105Test(working,laz); // 4-body-homog: repeated in ASCM09-zerodim
	Sys3110Test(working,laz); // Arnborg-Lazard: repeated in ASCM09-zerodim
	Sys3111Test(working,laz); // Barry: repeated in ASCM09-zerodim
	Sys3112Test(working,laz); // Caprasse-Li: repeated in ASCM09-zerodim
	Sys3113Test(working,laz); // Caprasse: repeated in ASCM09-zerodim
	Sys3117Test(working,laz); // Czapor-Geddes-Wang: pass, but VERY slow to run (BPAS: 385s, Maple Test Total: 386s) [BPAS result is more complicated]
	Sys3120Test(working,laz); // GonzalezGonzalez: MD non-primitive polynomial and factorization not done?, repeated in ASCM09-zerodim
	Sys3145Test(working,laz); // Rose: slowish, repeated in ASCM09-zerodim
	Sys3149Test(working,laz); // Trinks-2: repeated in ASCM09-zerodim

// 	cerr << "\"ISSAC2005\":" << endl;
////	Sys2887Test(testing,laz); // Chemkin: very slow, gave up
	Sys2888Test(working,laz); // Cyclohexane: MD probably factorization
//	Sys2889Test(testing,laz); // Dessin-2: very slow, error, cannot add a constant to a TriangularSet (after 20m)
	Sys2890aTest(working,laz); // Eco4
	Sys2890bTest(working,laz); // Eco5
	Sys2890Test(working,laz); // Eco6
	// Sys2891Test(working,laz); // Eco7: slow (~16s) to run
	Sys2892Test(working,laz); // Fabfau: slow (~14s) to run
	Sys2893Test(working,laz); // Fee_1
//	Sys2894Test(testing,laz); // Gametwo5: very slow, massive coefficient swell
	Sys2895Test(working,laz); // Geneig
//	Sys2896Test(testing,laz); // ISSAC97: very slow, gave up after 4h
//	Sys2897Test(testing,laz); // Methan61: very slow, gave up after 4h
//	Sys2898Test(testing,laz); // Reimer-4: very slow, gave up after 4h
//	Sys2899Test(testing,laz); // Uteshev-Bikker: very slow, gave up after around 3h
	Sys2900Test(working,laz); // Weispfenning-94: slow (~10s), slower post-fac changes

//	cerr << "\"ASCM09-zerodim\":" << endl;
	Sys3105Test(working,laz); // 4-body-homog: slow to run laz (~50s) [BPAS result is more complicated]
	Sys3107Test(working,laz); // 5-body-homog: slow to run laz
	Sys3109Test(working,laz); // Arnborg-Lazard-rev: MD probably factorization
	Sys3110Test(working,laz); // Arnborg-Lazard
	Sys3111Test(working,laz); // Barry: Maple error lazBF: Error, invalid expression for eval
	Sys3112Test(working,laz); // Caprasse-Li
	Sys3113Test(working,laz); // Caprasse: segfault problem post-fac
	Sys3114Test(working,laz); // Chemical-reaction
	Sys3115Test(working,laz); // Circles
	Sys3116Test(working,laz); // Cyclic-5: kalk BPAS: error, pseudo-dividend is zero from SMQP.; Jan17/2020: No issues, very fast after cleanset and RittWu ordering
	Sys3117Test(working,laz); // Czapor-Geddes-Wang: MD, pass without GCDFreeFact, now slow (~6m) to run, exception post-fac
	Sys3118Test(working,laz); // Fabfaux: slow (~14s)
	Sys3119Test(working,laz); // Geometric-constraints
	Sys3120Test(working,laz); // GonzalezGonzalez
	Sys3121Test(working,laz); // Katsura-4
	Sys3122Test(working,laz); // LHLP1
	Sys3123Test(working,laz); // LHLP2
	Sys3124Test(working,laz); // LHLPv
	Sys3125Test(working,laz); // LHLP4
	Sys3126Test(working,laz); // LHLP5
	Sys3127Test(working,laz); // LHLP6
	Sys3128Test(working,laz); // Neural-network
	Sys3129Test(working,laz); // NLD-10-3
	Sys3130Test(working,laz); // NLD-3-4: 
	// Sys3131Test(working,laz); // NLD-3-5: slow (~1m) to run, segfault lazBF, (~46s) kalk
//	Sys3132Test(testing,laz); // NLD-4-5: very slow, gave up after 90m kalk FAILED TODO: verify this is still the case.
	Sys3133Test(working,laz); // NLD-7-3
	Sys3134Test(working,laz); // NLD-8-3: slow (~10s) to run, blad out of memory exception
	// Sys3135Test(working,laz); // NLD-9-3: was slow (~45s) to run, now (~6-7m) to run, (~17s) kalk
	// Sys3136Test(working,laz); // NQL-10-2: few seconds kalk; Jan15/2020: easy peasy
	// Sys3137Test(testing,laz); // NQL-10-4: slow, gave up kalk (expression swell); Jan15/2020: error in maple factorization
	// Sys3138Test(testing,laz); // NQL-15-2: very slow (~220s) to run, Maple dnf after 1h, (~3s) kalk but validation behaves strangely; Jan22/2020: validation probably fixed, but in lazard mode fails to factorize via maple.
	// Sys3139Test(working,laz); // NQL-5-4: blad out of memory exception, (<1s) kalk but slow (~1m) to validate
	Sys3140Test(working,laz); // P3P-special: segfault problem post-fac
	Sys3141Test(working,laz); // PlateForme2d-easy: slow to run, <10 s post-fac	
	Sys3142Test(working,laz); // R-5: MD factorization not done
	Sys3143Test(working,laz); // R-6
	// Sys3144Test(working,laz); // Reif: (156s) kalk
	Sys3145Test(working,laz); // Rose
	Sys3146Test(working,laz); // Simple-NQL-20-30: was much slower than Maple, segfault problem post-fac; now (~1s) kalk; Jan15/2020: easy-ish (4s laz)
	Sys3147Test(working,laz); // Takeuchi-Lu
//	Sys3148Test(testing,laz); // Themos-net: slow, gave up kalk at [SRC] calling a.subresultantChain(b,1):
	Sys3149Test(working,laz); // Trinks-2
	Sys3150Test(working,laz); // Trinks-difficult
//	Sys3151Test(testing,laz); // Uteshev-Bikker: slow, gave up (expression swell); Jan15/2020: error, cannot add a constant to a TriangularSet
	Sys3152Test(working,laz); // Wilkinson20	
	Sys3153Test(working,laz); // Wilkinsonxy

	// cerr << "\"CM09-intersect-posdim\":" << endl;
//	Sys3019Test(testing,laz); // 4corps-1parameter-homog:  slow kalk, gave up at [SRC] calling a.subresultantChain(b,1):
//	Sys3020Test(working,laz); // 8-3-config-Li: very slow (~30m) laz, (~24s) kalk (kalk PASSED with 90m validation)
	Sys3021Test(working,laz); // AlKashiSinus: slow laz val (~40s), slow laz (~140s), (~1s) kalk, (~30s) kalk val
	Sys3023Test(working,laz); // Alonso
////	Sys3022Test(fixing,laz); // Alonso-Li: (86m) kalk val TODO: fix can't add constant error [caused by improper sqf handling]
//	Sys3024Test(testing,laz); // Bezier:  slow kalk, gave up at [SRC] calling a.subresultantChain(b,1):
//	Sys3025Test(testing,laz); // Bjork60: TODO: kalk BPAS: error, a polynomial with the leading variable u already exists in the TriangularSet
	Sys3027Test(working,laz); // Blood-Coagulation
	Sys3026Test(working,laz); // Blood-Coagulation-2: segfault problem post-fac, (<~5s) kalk
	Sys3157Test(working,laz); // BM05-1
//	Sys3158Test(working,laz); // BM05-2: slow val (~95s) laz, very slow val lazB (>~5h)
	Sys3028Test(working,laz); // Bronstein-Wang: segfault lazBF
//	Sys3029Test(fixing,laz); // Butcher: Maple result cleaner, segfault lazBF; TODO: something is wrong in removeZero (prem to 0 not -8) :(
	Sys3030Test(working,laz); // CDC2-Cyclin
//	Sys3031Test(testing,laz); // Cheaters-homotopy-easy: slow kalk, gave up at [SRC] calling a.subresultantChain(b,1):
//	Sys3032Test(testing,laz); // Cheaters-homotopy-hard: (~s) kalk, TODO: complete kalk val
//	Sys3033Test(testing,laz); // Chemical: (~8s) kalk, TODO: complete kalk val
//	Sys3034Test(testing,laz); // ChildDraw-1: slow, gave up kalk (massive expression swell)
//	Sys3035Test(testing,laz); // ChildDraw-2: slow, gave up kalk (massive expression swell)
	Sys3036Test(working,laz); // Cinquin-Demongeot-3-1: exception problem post-fac, (<1s) kalk
	Sys3037Test(working,laz); // Cinquin-Demongeot-3-2: Maple result much simpler, (<1s) kalk
//	Sys3038Test(working,laz); // Cinquin-Demongeot-3-3: (1.5s) kalk, TODO: laz test
//	Sys3039Test(testing,laz); // Cinquin-Demongeot-3-4: (~5s) kalk, TODO: complete kalk val
//	Sys3040Test(working,laz); // Cinquin-Demongeot-5-1: passes lazBF, pre-fac segfault after 97m, (2s kalk), TODO test with primfac
//	Sys3041Test(working,laz); // Collins-JSC02: (<1s) kalk, TODO: laz test
	Sys3042Test(working,laz); // Cox-ISSAC07: exception problem lazBF, (<10s) kalk
	Sys3043Test(working,laz); // Cyclic-4
	Sys3179Test(working,laz); // DescartesFolium
//	Sys3156Test(testing,laz); // DGP29: *** Error in `./tests/RegularChain_test_sys_exe': corrupted size vs. prev_size: 0x000000000b6984d0 ***
	Sys3155Test(working,laz); // DGP6
	Sys3045Test(working,laz); // DonatiTraverso: segfault problem post-fac
	// Sys3044Test(working,kalk); // DonatiTraverso-rev: slow, gave up kalk at calling triangularSetNormalForm...
	Sys3181Test(working,laz); // EdgeSquare
//	Sys3159Test(working,laz); // Ellipse: note that this system is currently trivial without the constraints imposed
	Sys3180Test(working,laz); // EnneperSurface
//	Sys3046Test(working,laz); // F-633: slow (~120s) to run, exception lazBF
//	Sys3047Test(testing,laz); // F-744: slow, gave up kalk
	Sys3048Test(working,laz); // GenLinSyst-2-2: (<1s) kalk
	Sys3049Test(working,laz); // GenLinSyst-3-2
//	Sys3050Test(working,laz); // GenLinSyst-3-3: slow val lazB (~70s), super slow lazBF, (<1s) kalk
//	Sys3051Test(testing,laz); // Geometry-Butterfly_3: slow, eventually gave up kalk at [SRC] calling a.subresultantChain(b,1):
	Sys2926Test(working,laz); // Gerdt
	Sys3053Test(working,laz); // Gonnet: slow (~35s) to run
//	Sys3190Test(testing,laz); // Haas5: (~14s) kalk, TODO: kalk val
	Sys3054Test(working,laz); // Hairer-2-BGK: super slow lazBF
//	Sys3055Test(testing,laz); // Hawes-1: slow, quickly gave up kalk at [SRC] calling a.subresultantChain(b,1):
	Sys3056Test(working,laz); // Hereman-2
//	Sys3057Test(testing,laz); // Hereman-8: TODO: kalk BPAS: error, cannot add a constant to a TriangularSet
//	Sys3058Test(testing,laz); // Hoffmann: slow, quickly gave up kalk at [SRC] calling a.subresultantChain(b,1):
//	Sys3160Test(working,laz); // IBVP: (<1s) kalk, TODO: laz test
	Sys3174Test(working,laz); // Jirstrand22
	Sys3175Test(working,laz); // Jirstrand23: exception problem post-fac, (<1s) kalk, (3.4h) kalk val
	Sys3176Test(working,laz); // Jirstrand24
	Sys3177Test(working,laz); // Jirstrand41
	Sys3178Test(working,laz); // Jirstrand42: (<1s) kalk
//	Sys3191Test(working,laz); // John5: (<1s) kalk, (~900s) kalk val, TODO: laz test
	Sys3192Test(working,laz); // Kahan-ellipsoid
//	Sys3059Test(testing,laz); // KdV: *** Error in `./tests/RegularChain_test_sys_exe': corrupted size vs. prev_size: 0x0000000016cc4240 ***
	Sys3154Test(working,laz); // L: slow val (~40s), (<1s) kalk, kalk val UNKNOWN (killed after 26h)
	Sys3161Test(working,laz); // Lafferriere35
	Sys3162Test(working,laz); // Lafferriere37
	Sys3060Test(working,laz); // Lanconelli
	Sys3062Test(working,laz); // Laurent-96: (<1s) kalk
	Sys3061Test(working,laz); // Laurent-96-2: (<1s) kalk
	Sys3063Test(working,laz); // Laurent-98: Maple result slightly simpler, (<1s) kalk
//	Sys3064Test(working,laz); // Lazard-ASCM2001: (~2.5s) kalk, TODO: laz test, check that system is right, calls itself A, System 30
//	Sys3065Test(testing,laz); // Leykin-1: slow, quickly gave up after reduce call
//	Sys3066Test(working,laz); // Lichtblau: (<1s) kalk, TODO: laz test
//	Sys3068Test(testing,laz); // Liu-Lorenz: (~2s) kalk, TODO: complete kalk val
//	Sys3067Test(testing,laz); // Liu-Lorenz-Li: (~1.5s) kalk, TODO: complete kalk val
	// Sys3069Test(working,laz); // MacLane: slow to val (~60s), segfault lazBF, (<2s) kalk
	Sys3183Test(working,laz); // Mehta0
//	Sys3184Test(working,laz); // Mehta1: (~1s) kalk, TODO: laz test
//	Sys3185Test(testing,laz); // Mehta2: (~4s) kalk, TODO: complete very slow kalk val
//	Sys3186Test(testing,laz); // Mehta3: slow, gave up kalk at [SRC] calling a.subresultantChain(b,1): (equiv to later slow src comp in Mehta4?)
//	Sys3187Test(testing,laz); // Mehta4: slow, eventually gave up kalk at [SRC] calling a.subresultantChain(b,1):
//	Sys3188Test(testing,laz); // Mehta5: slow, quickly gave up kalk at [SRC] calling a.subresultantChain(b,1):
//	Sys3189Test(testing,laz); // Mehta6: slow, quickly gave up kalk at [SRC] calling a.subresultantChain(b,1):
//	Sys3070Test(testing,laz); // MESFET-3: BLAD error, factorization is wrong., TODO: ensure reliable factorization
	// Sys3071Test(working,laz); // MontesS12: was slow laz val (~30s), super slow lazBF (dnf), (<1s) kalk
	// Sys3072Test(working,laz); // MontesS14: slow to validate (~15s), (<1s) kalk, (>~1m) kalk val
//	Sys3073Test(working,laz); // Morgenstein: (~2s) kalk, TODO: laz test
	Sys3163Test(working,laz); // MPV89
	Sys3074Test(working,laz); // Neural: (<1s) kalk
	Sys3193Test(working,laz); // Non-Equidim: Maple result simpler
	Sys3075Test(working,laz); // Noonburg: Maple result simpler, (<~1s) kalk
//	Sys3165Test(testing,laz); // P3P: slow, gave up kalk at [SRC] calling a.subresultantChain(b,1):
//	Sys3164Test(testing,laz); // P3P-Isosceles: (~12s) kalk, TODO: kalk val, test laz
//	Sys3079Test(working,laz); // Pappus: was slow (~40s), now (~7s) lazB, very slow val (~2000s), (~4s) kalk
//	Sys3080Test(testing,laz); // Pavelle: slow, gave up, segault lazBF; (~20s) kalk, TODO: kalk val, test laz
//	Sys3081Test(testing,laz); // Pin1: kalk, segfault in SRC class
	Sys3082Test(working,laz); // Prey-Predator
	Sys3166Test(working,laz); // Putnam
	Sys3083Test(working,laz); // Raksanyi: (~2s) kalk
//	Sys3084Test(testing,laz); // RationalInterp: slow, gave up after 1h kalk
	Sys3087Test(working,laz); // RobotPlano: hanging problem lazBF, (<1s) kalk
	Sys3086Test(working,laz); // RobotPlanoEasy
	Sys3167Test(working,laz); // SEIT: was slow to validate (~50s)
//	Sys3168Test(working,laz); // Solotareff-4a: very slow (~1s) kalk, TODO: time laz
//	Sys3169Test(working,laz); // Solotareff-4b: (~1s) kalk, significant MD, TODO: test laz
	Sys3088Test(working,laz); // Trivial-10
	Sys3089Test(working,laz); // Trivial-6
	Sys3090Test(working,laz); // Trivial-7
	Sys3091Test(working,laz); // Trivial-8: MD factorization not done
//	Sys3092Test(working,laz); // Vermeer: (~2s) kalk, TODO: test laz
	Sys3094Test(working,laz); // Wang-1991c: super slow post-fac, (~1s) kalk
	Sys3095Test(working,laz); // Wang93: slow, HUGE MD (factorization should help), (~8s) kalk, (<1m) kalk val
//	Sys3096Test(fixing,laz); // Wu-Wang: error post-incremental change, corrupted size vs. prev_size: 0x0000000018e80680 or free(): invalid next size (fast): 0x00000000060b3fa0
	Sys3097Test(working,laz); // WX-AB05-1
//	Sys3170Test(testing,laz);  // Xia: (~10s) to run, slow val (~200s), (<1s) kalk
	Sys3098Test(working,laz); // Xia-ISSAC07-1
	Sys3099Test(working,laz); // Xia-ISSAC07-2
	Sys3100Test(working,laz); // Xia-reachability
	Sys3101Test(working,laz); // XY-5-5-1
//	Sys3102Test(fixing,laz); // XY-5-7-2: was slow (~33s), (~8s) lazB, much slower lazMF (coef-swell, dies after factor returns), (~300s) val, (~1s) kalk; now (~20s) kalk, TODO: test laz
	Sys3103Test(working,laz); // Yang011104: (<2s) kalk
	Sys3104Test(working,laz); // YangBaxterRosso

//	cerr << "\"Uncategorized\":" << endl;
//	Sys2250Test(working,laz);


/// 1303 PseudoDivide Issue ///
////SMQP p("c_1^2 + c_2*x_2 - 2*x_2^2");
//SMQP p("c_2^2*x_2 - c_2*x_1*x_2 - c_2*x_2^2");
//vector<Symbol> vars = {Symbol("c_1"), Symbol("c_2"), Symbol("x_1"), Symbol("x_2")};
//RegularChain<RN,SMQP> T(vars);
//SMQP q("x_1 - x_2");
//q.setRingVariables(vars);
//T += q;
//q = SMQP("c_1");
//q.setRingVariables(vars);
//T += q;
//q = T.pseudoDivide(p);
//cerr << "q = " << q << std::endl;
//q = T.reduce(p);
//cerr << "q = " << q << std::endl;
//q = T.normalForm(p);
//cerr << "q = " << q << std::endl;


}

int main() {
	// std::this_thread::sleep_for(std::chrono::milliseconds(100));
	triangularizeTests();
//	 triangularizeTestsParallel();

	// std::cerr << "Cumulative Maple worker time: " << g_MapleComputeTime << std::endl;
	// std::cerr << "Cumulative Maple wait time: " << g_mapleWaitTime << std::endl;
	return 0;
}
