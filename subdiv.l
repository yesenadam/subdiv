/* subdiv scanner v1 */ 
%option noyywrap

%{ 
# include "subdiv.tab.h" 
#define YY_NO_UNPUT

#include <stdio.h>
#include <string.h>

int yyerror(char *s);
//int yylineno = 0;
int points=0;
#define MAXNUMPTS 50
#define MAXSHAPEDEFS 20
#define MAXPOINTDEFS 20
#define MAXSHAPEPOINTSLINES 20
struct pt {
float x;
float y;
} pts[50];
%} 

%% 

".."[0-9]+|"..."[0-9]+ {yylval.i = atoi(yytext+2);return DOTNUM;}
[\-]?[0-9]+"."[0-9]+|[\-]?"."[0-9]+|[\-]?[0-9]+ { yylval.f = atof(yytext);return NUM; }
",0"|","[1-9][0-9]*   {yylval.i = atoi(yytext+1);return COMMANUM;} /* ",0" allowed, means p0 */     
"levels"    {return LEVELS;}
"width"     {return WIDTH;}
"margin"    {return MARGIN;}
"noedges"|"no edges"    {return NOEDGES;}
"onlyedges"|"only edges"|"edgesonly"|"edges only"  {return ONLYEDGES;}
"pbox"      {return PBOX;}
"polygon"|"pgon"      {return PGON;}
"pcube"      {return PCUBE;}
"boxcolour"|"boxcolor"|"background" {return BOXCOLOUR;}
"def"       {return DEF;}
"par"       {return PAR;}
"left"|"l" {yylval.i=-1; return TURN;}
"right"|"r" {yylval.i=1; return TURN;}
"srand" {return SRAND;}
"crand"      {return CRAND;}
"!"|"not"       {return NOT;}
"rot2"       {return ROT2;}
"rot"       {return ROT;}
"="         {return EQUAL;}
">"         {return GT;}
"flip"      {return FLIP;}
"grid"      {return GRID;}
"crandom"      {return CRANDOM;}
"rand"      {return RAND;}
"wait" {return WAIT;}
"c"|"centre"|"center" {return CENTRE;}
"+"         {return PLUS;}
"- "   {return MINUS;} //trailing space!! I guess it will help... not to be confused with nums?
"/"         {return DIV;}
".."|"..."        {return DOTS;}
\" {return DITTO;}
[Pp][0-9]+  {yylval.i = atoi(yytext+1);return POINTNUM;} /*yytext+1 cuts initial P*/
[Ii][Pp][0-9]+  {yylval.i = atoi(yytext+2);return IPOINTNUM;} /*yytext+1 cuts initial P*/
^\n          {yylineno++;}/*ignore blank lines*/
\n         { yylineno++;return EOL; } 
":"[0-9c]+ {yylval.w=strdup(yytext+1);return PTSTR;} 
": "[0-9c]+ {yylval.w=strdup(yytext+2);return PTSTR;} 
[a-zA-Z][a-zA-Z0-9]*   {yylval.w=strdup(yytext);return WORD;} 
"//".*|"%".*|"#".*    /* ignore // or % or # comments */ 
[ \t]       { /* ignore whitespace */ } 
.           { yyerror("SCANNER ERROR");	}
%% 

