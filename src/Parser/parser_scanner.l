%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser_type.h"
#include "parser_grammar.tab.h"

%}

bpas_digit				[0-9]+
bpas_variables			[a-zA-Z_][a-zA-Z0-9_]*
bpas_whitespace			[ \t]*
bpas_newline			[\n]
bpas_period				.
bpas_ratnum				[0-9]+\/[0-9]+
bpas_minus				"-"
bpas_plus				"+"
bpas_divide				"/"
bpas_l_brace			"("
bpas_r_brace			")"
bpas_multiply			"*"
bpas_comma				","
bpas_l_square_brace		"["
bpas_r_square_brace		"]"
bpas_power				"^"
bpas_quit				"quit()"

%%
{bpas_variables}		{
							yylval.string_type = malloc(sizeof(char)*strlen(yytext)+1); 
							strncpy(yylval.string_type, yytext, strlen(yytext)+1);
							return VAR;
						}
{bpas_plus}				{return PLUS;}
{bpas_minus}			{return MINUS;}
{bpas_quit}				{
							exit(0);	
						}

{bpas_power}			{return POWER;}
{bpas_divide}			{return DIVIDE;}
{bpas_multiply}			{return MULTIPLY;}
{bpas_ratnum}			{
							yylval.string_type = malloc(sizeof(char)*strlen(yytext)+1); 
							strncpy(yylval.string_type, yytext, strlen(yytext)+1);
							return RATNUM;	
						}
{bpas_digit}			{
							yylval.string_type = malloc(sizeof(char)*strlen(yytext)+1); 
							strncpy(yylval.string_type, yytext, strlen(yytext)+1);
							return NUM;
						}

{bpas_l_brace}			{return L_BRACE;}
{bpas_r_brace}			{return R_BRACE;}
{bpas_comma}			{return COMMA;}
{bpas_l_square_brace}	{return LS_BRACE;}
{bpas_r_square_brace}	{return RS_BRACE;}
{bpas_whitespace}		{}
{bpas_period}			{fprintf(stderr, "BPAS ERROR: Unable to parse the character:  %s\n", yytext); exit(1);}
{bpas_newline}			{
							return *yytext;
						}

%%

int yywrap(){
	return 1;
}
