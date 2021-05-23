

#include "MapleInterface/MapleInterface.h"
#include "Parser/bpas_parser.h"

#include "Utils/Unix_Timer.h"

void mapleKernelTextCB(void* data, int tag, const char* output) {
	fprintf(stdout, "Maple output: %s\n", output);
}

void mapleKernelErrorCB(void* data, long int tag, const char* output) {
	fprintf(stderr, "MAPLE ERROR:: %s\n", output);
	exit(1);
}

MKernelVector* mapleKVSingleton(int set, MKernelVector kv) {
	static MKernelVector singleton = NULL;
	if (set) {
		singleton = kv;
	}
	return &singleton;
}

MKernelVector* getMapleKVSingleton() {
	startMapleKernel(0, NULL);
	return mapleKVSingleton(0, NULL);
}

void startMapleKernel(int argc, char* argv[]) {
	MKernelVector* kv_p = mapleKVSingleton(0, NULL);
	if (*kv_p != NULL) {
		//already started
		return;
	}

	char err[2048];
	MCallBackVectorDesc cb = { mapleKernelTextCB,
								mapleKernelErrorCB,
								0,
								0,
								0,
								0,
								0,
								0
							};
	MKernelVector kv = NULL;
	if ( (kv=StartMaple(argc, argv, &cb, NULL, NULL, err)) == NULL ) {
		fprintf(stderr, "MapleInterface ERROR: Cannot start maple kernel: %s\n", err);
		exit(1);
	}
	mapleKVSingleton(1, kv);

//TODO!! for benchmarking always maple in serial.
// #if defined(SERIAL) && SERIAL
    char evalSerialStr[] = "kernelopts(numcpus=1):";
    EvalMapleStatement(*kv_p, evalSerialStr);
// #endif
}


void stopMapleKernel() {
	MKernelVector* kv_p = mapleKVSingleton(0, NULL);
	StopMaple(*kv_p);
}


void restartMapleKernel() {
	MKernelVector* kv_p = mapleKVSingleton(0, NULL);
	char err[2048];
	if (RestartMaple(*kv_p, err) == 0) {
		fprintf(stderr, "MapleInterface ERROR: Cannot restart the marple kernel %s\n", err);
		exit(1);
	}
//TODO!! for benchmarking always maple in serial.
// #if defined(SERIAL) && SERIAL
    char evalSerialStr[] = "kernelopts(numcpus=1):";
    EvalMapleStatement(*kv_p, evalSerialStr);
// #endif
}

/**
 * Helper function for ALGEB to string.
 */
char* algebToString_MplInt(ALGEB in) {
	MKernelVector* kv_p = mapleKVSingleton(0, NULL);

	char printStr[] = "sprintf:";
	ALGEB f = EvalMapleStatement(*kv_p, printStr);
	ALGEB resStrObj = EvalMapleProc(*kv_p, f, 2, ToMapleString(*kv_p, "%a"), in);
	char* resStr = MapleToString(*kv_p, resStrObj);
	size_t resSize = strlen(resStr);

	//+1 for null char.
	char* ret = (char*) malloc(sizeof(char)*(resSize+1));
	memcpy(ret, resStr, sizeof(char)*(resSize+1));
	return ret;
}

ALGEB convertToMaple_AAZ(const AltArrZ_t* aa, const char** syms) {
	char* buff;
	size_t size;
	FILE* stream;

	stream = open_memstream (&buff, &size);

	printPoly_AAZ(stream, aa, syms, aa->nvar);
	fprintf(stream, ":");


	fflush (stream);
	fclose (stream);

	MKernelVector* kv_p = getMapleKVSingleton();
	ALGEB f = EvalMapleStatement(*kv_p, buff);

	free(buff);

	return f;
}



char* gcd_MplInt_string(const char* ac, const char* bc) {
	// timer_id id = start_timer();
	MKernelVector* kv_p = getMapleKVSingleton();
    // timer_time elapsed1 = elapsed_time(&id);
    // double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Start kernel time: %f\n", time );


	// id = start_timer();
	ALGEB a = EvalMapleStatement(*kv_p, ac);
	ALGEB b = EvalMapleStatement(*kv_p, bc);

	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert to maple time: %f\n", time );


	// id = start_timer();
	char gcdStr[] = "gcd:";
	ALGEB gcdFunc = EvalMapleStatement(*kv_p, gcdStr);
	ALGEB g = EvalMapleProc(*kv_p, gcdFunc, 2, a, b);
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "gcd time: %f\n", time );


	// id = start_timer();
	char* g_str = algebToString_MplInt(g);

	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert from maple time: %f\n", time );

	return g_str;
}


AltArrZ_t* gcd_MplInt_AAZ(const AltArrZ_t* aa, const AltArrZ_t* ba, const char** syms) {


	// timer_id id = start_timer();
	MKernelVector* kv_p = getMapleKVSingleton();
    // timer_time elapsed1 = elapsed_time(&id);
    // double time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Start kernel time: %f\n", time );


	// id = start_timer();
	ALGEB a = convertToMaple_AAZ(aa, syms);
	ALGEB b = convertToMaple_AAZ(ba, syms);
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert to maple time: %f\n", time );


	// id = start_timer();
	char gcdStr[] = "gcd:";
	ALGEB gcdFunc = EvalMapleStatement(*kv_p, gcdStr);
	ALGEB g = EvalMapleProc(*kv_p, gcdFunc, 2, a, b);
	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "gcd time: %f\n", time );


	// id = start_timer();
	char* g_str = algebToString_MplInt(g);
	AltArr_t* retQ = generate_altarr_var_defined(g_str, syms, aa->nvar);
	AltArrZ_t* retZ = deepCopyPolynomial_AAZFromAA(retQ);
	freePolynomial_AA(retQ);
	free(g_str);

	// elapsed1 = elapsed_time(&id);
 //    time = (elapsed1.tv_sec + ((double)elapsed1.tv_usec / 1000000));
	// fprintf(stderr, "Convert from maple time: %f\n", time );

	return retZ;

}

void factorPolynomial_MplInt_string(const char* poly, int* numFacts, char*** retFactsPtr, char*** retExps, char** numericFact) {

	MKernelVector* kv_p = getMapleKVSingleton();
	restartMapleKernel();

	ALGEB a = EvalMapleStatement(*kv_p, poly);

	char factStr[] = "factors:";
	ALGEB factorFunc = EvalMapleStatement(*kv_p, factStr);
	ALGEB factList = EvalMapleProc(*kv_p, factorFunc, 1, a);

	// char verifyStr[] = "verify:";
	// ALGEB verifyFunc = EvalMapleStatement(*kv_p, verifyStr);
	// ALGEB didFactorAlg = EvalMapleProc(*kv_p, verifyFunc, 2, a, factList);
	// M_BOOL didFactor = MapleToM_BOOL(*kv_p, didFactorAlg);

	// if (didFactor) {
	// 	return 0;
	// }


	char nopsStr[] = "nops:";
	ALGEB nopsFunc = EvalMapleStatement(*kv_p, nopsStr);

	char opStr[] = "op:";
	ALGEB opFunc = EvalMapleStatement(*kv_p, opStr);

	ALGEB ringElemAlg = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 1), factList);
	if (numericFact != NULL) {
		*numericFact = algebToString_MplInt(ringElemAlg);
	}

	ALGEB facts = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 2), factList);
	ALGEB nopsAlg = EvalMapleProc(*kv_p, nopsFunc, 1, facts);
	M_INT nops = MapleToM_INT(*kv_p, nopsAlg);

	char** retFacts = (char**) malloc(sizeof(char*) * nops);
	char** exps = (char**) malloc(sizeof(char*) * nops);

	for (int i = 1; i <= nops; ++i) {
		ALGEB factAndExp = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, i), facts);

		ALGEB op1 = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 1), factAndExp);
		ALGEB op2 = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 2), factAndExp);

		retFacts[i-1] = algebToString_MplInt(op1);
		exps[i-1] = algebToString_MplInt(op2);
	}

	if (numFacts != NULL) {
		*numFacts = nops;
	}

	if (retFacts != NULL) {
		*retFactsPtr = retFacts;
	} else {
		for (int i = 0; i < nops; ++i) {
			free(retFacts[i]);
		}
		free(retFacts);
	}

	if (retExps != NULL) {
		*retExps = exps;
	} else {
		for (int i = 0; i < nops; ++i) {
			free(retExps[1]);
		}
		free(retExps);
	}
}

void factorPolynomial_MplInt_AAZ(const AltArrZ_t* aa, const char** syms, int* numFacts, AltArrZ_t*** retFactsPtr, int** retExps, mpz_t numericFact) {


	MKernelVector* kv_p = getMapleKVSingleton();
	restartMapleKernel();

	ALGEB a = convertToMaple_AAZ(aa, syms);

	char factStr[] = "factors:";
	ALGEB factorFunc = EvalMapleStatement(*kv_p, factStr);
	ALGEB factList = EvalMapleProc(*kv_p, factorFunc, 1, a);

	// char verifyStr[] = "verify:";
	// ALGEB verifyFunc = EvalMapleStatement(*kv_p, verifyStr);
	// ALGEB didFactorAlg = EvalMapleProc(*kv_p, verifyFunc, 2, a, factList);
	// M_BOOL didFactor = MapleToM_BOOL(*kv_p, didFactorAlg);

	// if (didFactor) {
	// 	return 0;
	// }



	char nopsStr[] = "nops:";
	ALGEB nopsFunc = EvalMapleStatement(*kv_p, nopsStr);

	char opStr[] = "op:";
	ALGEB opFunc = EvalMapleStatement(*kv_p, opStr);

	mpq_t ringElem;
	mpq_init(ringElem);
	ALGEB ringElemAlg = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 1), factList);
	char* reStr = algebToString_MplInt(ringElemAlg);
	mpz_set_str(mpq_numref(ringElem), reStr, 10);
	free(reStr);

	ALGEB facts = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 2), factList);
	ALGEB nopsAlg = EvalMapleProc(*kv_p, nopsFunc, 1, facts);
	M_INT nops = MapleToM_INT(*kv_p, nopsAlg);

	AltArrZ_t** retFacts = (AltArrZ_t**) malloc(sizeof(AltArrZ_t*) * nops);
	int* exps = (int*) malloc(sizeof(int) * nops);

	mpq_t cont;
	mpq_init(cont);
	for (int i = 1; i <= nops; ++i) {
		ALGEB factAndExp = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, i), facts);

		ALGEB op1 = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 1), factAndExp);
		ALGEB op2 = EvalMapleProc(*kv_p, opFunc, 2, ToMapleInteger(*kv_p, 2), factAndExp);

		char* op_str = algebToString_MplInt(op1);
		AltArr_t* opQ = generate_altarr_var_defined(op_str, syms, aa->nvar);

		AltArrZ_t* opZ = primitivePartAndContent_AAZFromAA(opQ, cont);

		// retFacts[i-1] = deepCopyPolynomial_AAZFromAA(opQ);
		retFacts[i-1] = opZ;
		exps[i-1] = MapleToM_INT(*kv_p, op2);
		for (int k = 0; k < exps[i-1]; ++k) {
			mpq_mul(ringElem, ringElem, cont);
		}

		freePolynomial_AA(opQ);
		free(op_str);
	}

	mpz_set(numericFact, mpq_numref(ringElem));
	mpq_clears(cont, ringElem, NULL);

	if (numFacts != NULL) {
		*numFacts = nops;
	}

	if (retFacts != NULL) {
		*retFactsPtr = retFacts;
	} else {
		for (int i = 0; i < nops; ++i) {
			freePolynomial_AAZ(retFacts[i]);
		}
		free(retFacts);
	}

	if (retExps != NULL) {
		*retExps = exps;
	} else {
		free(retExps);
	}

	return;
}


int triangularizeValidate_MplInt(const char** inputs, int nInputs, int isLazard) {

	MKernelVector* kv_p = getMapleKVSingleton();
	restartMapleKernel();

    char* evalSerialStr = "kernelopts(numcpus=1);";
    EvalMapleStatement(*kv_p, evalSerialStr);

	ALGEB algebList[nInputs];
	for (int i = 0; i < nInputs; ++i) {
        ALGEB res = EvalMapleStatement(*kv_p, inputs[i]);
        algebList[i] = (res);
    }

    char* procStrLazard = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R, 'output' = 'lazard','radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    char* procStrKalk = "TriangularizeSet := proc (F::list, in_rc::list, Rlist::list) local lrc, n, rc, R, i, results; R := RegularChains:-PolynomialRing(Rlist); rc := RegularChains:-ChainTools:-Chain(in_rc, RegularChains:-ChainTools:-Empty(R), R); lrc := RegularChains:-Triangularize(F, rc, R,'radical'='no'); n := nops(lrc); results := []; for i to n do results := [op(results), RegularChains:-Equations(op(i, lrc), R)] end do; results end proc:";
    ALGEB testProc;
    if(isLazard) {
    	testProc = EvalMapleStatement(*kv_p, procStrLazard);
    }
    else {
    	testProc = EvalMapleStatement(*kv_p, procStrKalk);
    }

    ALGEB result = EvalMapleProc(*kv_p, testProc, 3, algebList[0], algebList[1], algebList[2]);

    ALGEB algebList2[3];
    algebList2[0] = algebList[3];
    algebList2[1] = result;
    algebList2[2] = algebList[2];

    char* evalStr = "kernelopts(numcpus=24);";
    EvalMapleStatement(*kv_p, evalStr);

	char* easyValProc = "checkExactSets := proc(mrcs, brcs) local i,j,k,tlist,found,check,l; if (evalb(nops(mrcs) <> nops(brcs))) then  return false; end if; for i from 1 to nops(mrcs) do  found := false;  for j from 1 to nops(brcs) do  if (evalb(nops(op(j,brcs)) <> nops(op(i,mrcs)))) then   next;  end if;  tlist := [seq(evalb((expand(op(l,op(j,brcs))) - expand(op(l,op(i,mrcs)))=0)), l=1..nops(op(j,brcs)))];  check := true;  for k from 1 to nops(op(j,brcs)) do   check := check and tlist[k];   end do;  if (evalb(check)) then   found := true;   break;  end if;  end do;  if (evalb(found = false)) then  return false;  end if; end do; return true; end proc:";
	char* validateProcLazard = "EqualAsConstructibleSets := proc (dec1::list, dec2::list, Rlist::list) local n, lrc1, lrc2, lrs1, lrs2, cs1, cs2, R, i, rc; R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lrc1 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec1)), RegularChains:-ChainTools:-Empty(R), R); lrc1 := [op(lrc1), rc] end do; n := nops(dec2); lrc2 := []; for i to n do rc := RegularChains:-ChainTools:-Chain(ListTools:-Reverse(op(i, dec2)), RegularChains:-ChainTools:-Empty(R), R); lrc2 := [op(lrc2), rc] end do; if evalb(nops(lrc1) = 1) and evalb(nops(lrc2) = 1) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc1), R) and RegularChains:-ChainTools:-IsEmptyChain(op(1, lrc2), R) then return true end if; lrs1 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc1, [1], R); lrs2 := map(RegularChains:-ConstructibleSetTools:-RegularSystem, lrc2, [1], R); cs1 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs1, R); cs2 := RegularChains:-ConstructibleSetTools:-ConstructibleSet(lrs2, R); if RegularChains:-ConstructibleSetTools:-IsContained(cs1, cs2, R) and RegularChains:-ConstructibleSetTools:-IsContained(cs2, cs1, R) then return true else return false end if end proc:";
	char* validateProcKalk = "KalkbrenerTest := proc(dec1::list,dec2::list,Rlist::list) local g,G,rc,R,lsi1,lsi2,li,h,J,J1,J2,i,j,n,m,pass::boolean; if evalb(nops(dec1) = 0) and evalb(nops(dec2) = 0) then return true; end if;  R := RegularChains:-PolynomialRing(Rlist); n := nops(dec1); lsi1 := []; writeline(terminal,convert(ConstructSaturatedIdeals,string)); writeline(terminal,convert(SaturatedIdealsForFirstDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec1[i]; m := nops(rc); h := 1; for j to m-1 do  h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi1 := [op(lsi1), J]; end do; n := nops(dec2); lsi2 := []; writeline(terminal,convert(SaturatedIdealsForSecondDecomposition,string)); for i to n do writeline(terminal,convert(computingSaturatedIdeal,string)); rc := dec2[i]; m := nops(rc); h := 1; for j to m-1 do h := h * RegularChains:-Initial(rc[j],R); end do; J := PolynomialIdeals:-PolynomialIdeal(rc); J := PolynomialIdeals:-Saturate(J,h); lsi2 := [op(lsi2), J]; end do; writeline(terminal,convert(BeginMutualInclusionOfSaturatedIdealsTest,string)); n := nops(lsi1); m := nops(lsi2); for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(not pass) then writeline(terminal,convert(SaturatedIdealsFromFirstListNotIncludedInSecond,string)); end if; if evalb(pass) then  for i to n do for j to m do if PolynomialIdeals:-IdealContainment(lsi1[i],lsi2[j],lsi1[i]) then pass := true; break; else pass := false; end if; end do; if evalb(pass) then next; else break; end if; end do; if evalb(pass) then return true; end if; end if; writeline(terminal,convert(SaturatedIdealsFromSecondListNotIncludedInFirst,string)); writeline(terminal,convert(BeginIdentityOfIntersectionsOfSaturatedIdealsTest,string)); J1 := lsi1[1]; if n > 1 then for i from 2 to n do J1 := PolynomialIdeals:-Intersect(J1,lsi1[i]); end do; end if; J2 := lsi2[1]; if m > 1 then for i from 2 to m do J2 := PolynomialIdeals:-Intersect(J2,lsi2[i]); end do; end if; if PolynomialIdeals:-IdealContainment(J1,J2,J1) then return true; end if; for i from 1 to n do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi1[i])); for j from 1 to m do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi2[j]) then else pass := false; break; end if; end do; if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break;	end if; end do; if evalb(pass) = true then for j from 1 to m do G := PolynomialIdeals:-Generators(PolynomialIdeals:-Simplify(lsi2[j])); for i from 1 to n do pass := true; for g in G do if PolynomialIdeals:-RadicalMembership(g, lsi1[i]) then else pass := false; break;	end if; end do;	if evalb(pass) = true then break; end if; end do; if evalb(pass) = false then break; end if; end do; end if; if evalb(pass) = true then	return true; end if; J1 := PolynomialIdeals:-Radical(J1); J2 := PolynomialIdeals:-Radical(J2); if PolynomialIdeals:-IdealContainment(J1,J2,J1) then	return true; else return false;	end if; end proc:";

	ALGEB validateProc;

	validateProc = EvalMapleStatement(*kv_p, easyValProc);
	result = EvalMapleProc(*kv_p, validateProc, 2, algebList2[0], algebList2[1]);
	M_BOOL ret = MapleToM_BOOL(*kv_p, result);
	if (!ret) {

		if (isLazard) {
			validateProc = EvalMapleStatement(*kv_p, validateProcLazard);
	    }
	    else {
			validateProc = EvalMapleStatement(*kv_p, validateProcKalk);
	    }

	    result = EvalMapleProc(*kv_p, validateProc, 3, algebList2[0], algebList2[1], algebList2[2]);
	    ret = MapleToM_BOOL(*kv_p, result);
	    // ret = 1;
		// fprintf(stderr, "Did not pass exact set test!\n");
	} else {
		// fprintf(stderr, "Passed by exact set test!\n");
	}

    char* bpasstr = algebToString_MplInt(algebList2[0]);
    char* maplestr = algebToString_MplInt(algebList2[1]);
    fprintf(stderr, "BPAS result: \n");
    fprintf(stderr, "%s\n", bpasstr);
    fprintf(stderr, "Maple result: \n");
    fprintf(stderr, "%s\n", maplestr);
    fprintf(stderr, "BPAS size: %ld, Maple size: %ld\n", strlen(bpasstr), strlen(maplestr));

    return ret;
}