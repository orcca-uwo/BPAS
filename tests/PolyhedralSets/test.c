#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/PolyhedralSets/FME_Support_inequality.h"
#include "../../include/PolyhedralSets/FME_Support_unrolledll.h"
#include "../../include/PolyhedralSets/FME_Support_fme.h"


int run_via_maple(const char * maple_file_path)
{
//	char * cmd = (char*) malloc(128 + strlen(cmd_str_params));
//	sprintf(cmd, "%s %s\" %s", cmd_str, cmd_str_params, suffix);
//#if VERBOSE>=3
//	printf("executing [%s ]... \n", cmd);
//#endif

//	freopen("/dev/null", "w", stdout);

	const char pass_phrase[8] = "PASS";
	const char fail_phrase[8] = "FAIL";
	char cmd[256];
	sprintf(cmd, "which maple >/dev/null && echo '0'|| echo '1'");
//	printf("cmd=[%s]\n", cmd);
	FILE* pipe;
	pipe = popen(cmd, "r");
	if (pipe == NULL)
	{
		printf("ERROR: FAILED TO RUN [%s]\n", cmd);
		pclose(pipe);
		return (EXIT_FAILURE);
	}
	char buffer[1024];
	int value = 0;
	if (fgets(buffer, sizeof(buffer) - 1, pipe) != NULL)
	{
		value = atoi(buffer);
	}
	else
	{
		printf("ERROR: NOTHING READ!\n");
		value = -1;
	}

	pclose(pipe);

	if (value == 1)
	{
		printf("[WARNING: maple does NOT exist on this system!]\n");
		printf("[RESULT] ...... PASS\n");
		return (EXIT_FAILURE);
	}

	sprintf(cmd, "maple -q %s", maple_file_path);
	printf("cmd=[%s]\n", cmd);
	FILE* pipe_maple;
	pipe_maple = popen(cmd, "r");
	if (pipe_maple == NULL)
	{
		printf("[ERROR: FAILED TO RUN [%s]]\n", cmd);
		pclose(pipe_maple);
		return (EXIT_FAILURE);
	}

	while (fgets(buffer, sizeof(buffer) - 1, pipe) != NULL)
	{
		if (strstr(buffer, pass_phrase) != NULL)
			printf("[RESULT] ...... %s", buffer);
		if (strstr(buffer, fail_phrase) != NULL)
			printf("[RESULT] ...... %s", buffer);
	}
	pclose(pipe_maple);
	return EXIT_SUCCESS;
}

void makeMapleFile(inequality_t * input, unrolledLl_t * proj1, unrolledLl_t* proj2,
		unrolledLl_t* proj3, int varNum, int ineqNum, int elimNum)
{
	printf("Enter make maple function\n");
	FILE * file = fopen("mapleTest.mpl", "w");
	fprintf(file, "kernelopts(printbytes=false):\n" "with(PolyhedralSets):\n");
	fprintf(file, "elimNum := %d:", elimNum);
	fprintf(file, "l:=[");
	
//	printInequality(&(input[0]),'x');
	
	for (int i = 0; i < ineqNum; i++)
	{
		//printf("%d\n",i);
		for (int j = 0; j < varNum; j++)
		{
			fprintf(file, "+");
			fprintf(file, "( %s )", mpz_get_str(NULL, 10, input[i].coeff[j]));
			fprintf(file, "*x%d ", j);
		}
		//printf("finish j loop\n");	
		fprintf(file, "<= %s", mpz_get_str(NULL, 10, input[i].constant));
		//printf("const\n");
		if (i != ineqNum - 1)
			fprintf(file, ",");
	}
	printf("finish writing input polyhedron\n");
	
	fprintf(file, "]:\n");
	fprintf(file, "P := PolyhedralSet(l):\n");
	fprintf(file, "P1 := Project(P,[");
	for (int i = elimNum + 1; i < varNum; i++)
	{
		fprintf(file, "x%d", i);
		if (i != varNum - 1)
			fprintf(file, ",");
	}
	fprintf(file, "]):\n");
	fprintf(file, "nops(Relations(P1)):");

	fprintf(file, "l1:=[");

	
	
	inequalityNode_t *p1 = proj1->head->next;
	inequalityNode_t *p2 = proj2->head->next;
	inequalityNode_t *p3 = proj3->head->next;

	for (int i = 0; i < proj1->number; i++)
	{
		for (int j = 0; j < p1->fill; j++)
		{
			for (int k = 0; k < varNum; k++)
			{
				fprintf(file, "+");
				fprintf(file, "( %s )",
						mpz_get_str(NULL, 10, p1->data[j].coeff[k]));
				fprintf(file, "*x%d ", k);
			}
			fprintf(file, "<= %s", mpz_get_str(NULL, 10, p1->data[j].constant));
			fprintf(file, ",");
		}
		p1 = p1->next;
	}

	for (int i = 0; i < proj2->number; i++)
	{
		for (int j = 0; j < p2->fill; j++)
		{
			for (int k = 0; k < varNum; k++)
			{
				fprintf(file, "+");
				fprintf(file, "( %s )",
						mpz_get_str(NULL, 10, p2->data[j].coeff[k]));
				fprintf(file, "*x%d ", k);
			}
			fprintf(file, "<= %s", mpz_get_str(NULL, 10, p2->data[j].constant));
			fprintf(file, ",");
		}
		p2 = p2->next;
	}

	for (int i = 0; i < proj3->number; i++)
	{
		for (int j = 0; j < p3->fill; j++)
		{
			for (int k = 0; k < varNum; k++)
			{
				fprintf(file, "+");
				fprintf(file, "( %s )",
						mpz_get_str(NULL, 10, p3->data[j].coeff[k]));
				fprintf(file, "*x%d ", k);
			}
			fprintf(file, "<= %s", mpz_get_str(NULL, 10, p3->data[j].constant));
			fprintf(file, ",");
		}
		p3 = p3->next;
	}

	for (int i = 0; i <= elimNum; i++)
	{
		fprintf(file, "x%d = 0", i);
		if (i != elimNum)
			fprintf(file, ",");
	}

	const char * program_body = "]:\n"
			"Q := PolyhedralSet(l1):\n"
			"nops(Relations(Q)):\n"
			"Equal(P1 , Q):\n"
			"if nops(l1)=nops(Relations(P1)) then \n"
			"	if(Equal(P1,Q)=true) then \n"
			"		printf(\"LEVEL %d ..... PASS\\n\",elimNum): \n"
			"	else \n"
			"		printf(\"FAIL\\n\"):\n"
			"   fi: \n"
			"else \n"
			"	printf(\"FAIL\\n\"):\n"
			"fi: \n";
	fprintf(file, "%s", program_body);
	fclose(file);

}


void project_test(char * fileName, int varNum, int ineqNum)
{

	inequality_t * inputData = (inequality_t *) malloc(
			ineqNum * sizeof(inequality_t));
	getInputArray(fileName, inputData, varNum, ineqNum);

	matrix_t * initTestCone = (matrix_t *) malloc(sizeof(matrix_t));
	initialTestCone(inputData, ineqNum, initTestCone);

	FMEDS_t * inputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));
	FMEDS_t * outputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));

	allocFMEDS(varNum, inputFMEDS);

	allocFMEDS(varNum, outputFMEDS);


	for (int i = 0; i < ineqNum; i++)
	{
		addToFMEDS(inputFMEDS, &inputData[i], 0);
	}


//	for (int i = 0; i < ineqNum; i++)
//		freeInequality(&(inputData[i]));
//	free(inputData);

	matrix_t * theTestCone = (matrix_t *) malloc(sizeof(matrix_t));

	printf("********************FME data struct allocation Done\n");

	for (int i = 0; i < varNum - 1; i++)
	{
		printf("Step %d Statrs\n", i);

		allocMatrix(initTestCone->rowNum, initTestCone->colNum - i - 1,
				theTestCone);

		testCone(initTestCone, i + 1, theTestCone);

		oneStepFME(inputFMEDS, i, outputFMEDS, theTestCone);
                printf("fme done\n");
		
		makeMapleFile(inputData, outputFMEDS->posSubSys, outputFMEDS->negSubSys, outputFMEDS->zerSubSys, varNum, ineqNum, i);

		const char * cmd = "maple -q mapleTest.mpl";
		system(cmd);


		freeMatrix(theTestCone);

		freeFMEDS(inputFMEDS);

		inputFMEDS = outputFMEDS;

		outputFMEDS = (FMEDS_t *) malloc(sizeof(FMEDS_t));
		allocFMEDS(varNum , outputFMEDS);

	}

	freeMatrix(initTestCone);
	free(initTestCone);
	free(theTestCone);
	freeFMEDS(inputFMEDS);
	free(inputFMEDS);
	freeFMEDS(outputFMEDS);
	free(outputFMEDS);
}


int main(int argc, char ** argv)
{

	printf("\n ******** Testing %s ******** \n", argv[1]);

	int varNum = atoi(argv[2]);
	int ineqNum = atoi(argv[3]);

	project_test(argv[1], varNum, ineqNum);

	return 0;
}

