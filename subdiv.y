/* subdiv.y - bison parser for Ratsub graphics language - Adam Ponting, Nov 2019*/ 
%{ 
#define YY_NO_UNPUT // .. Unput = put chars back in input
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h>
#include <time.h>
#include <stdbool.h>
#define MAXINITPTS 100
#define MAXNUMPTS 100
#define MAXSHAPEDEFS 50
#define MAXPOINTDEFS 100 //per shape
#define MAXSHAPECALLS 50 //per shape
#define MAXCOLOURNAMES 200 //currently 146 colours, be sure to increase this if more colours added to struct
#define MAXSHAPESTOFLIP 20 //max length of shape list when "flip" command used.
#define PI 3.141592654
#define LEFT -1
#define RIGHT 1
#define YES 1
#define NO 0
#define EQUALS 1
#define GREATERTHAN 2
short ColNumber(char* cnamestr);
short GetShapeNum(char name[50]);
void BuildNewWaitShapes();
void CheckAndSetPoints(short cp1, short cp2, short cp3);
void CheckPoints(short cp1,short cp2,short cp3);
void CheckPtNumber();
void IncPointVars();
void IncShapePointVars();
void InitPtError(short ip);
void MakeLogFiles();
void PrintBranchedShapeCalls();
void PrintCentreProcedure();
void PrintColouringStuff();
void PrintDrawShape();
void PrintNewShapePoint(short j);
void PrintNewShapePoints();
void PrintShadingColProcedure();
void PrintShapeCall(short j);
void PrintShapeProcedure(short i);
void PrintTurnX();
void PrintTurnY();
void SetColNames();
void SetDecFract(int pfrom, int pto,float dfract);
void SetFromAndToPts(short fromP,short toP);
void SetGradientPtInfo(int gp, int ptnum, float c0, float c1, float c2);
void SetInitPt(float x, float y);
void SetPboxPoints(float x1, float y1, float x2, float y2);
void SetShadingCol(int cpn, float r,float g,float b);
void SetShapeCallCol(float r,float g,float b);
void SetShapeCol(float r,float g,float b);
int yyerror(char *s) __attribute__ ((noreturn));

struct {
    float x;
    float y;
    float col[3];
} initPt[MAXNUMPTS];

struct NewPtStruct{ //i.e. a "p3 p1 1/2 p2" line
    short   From;
    short   To;
    int     numer;
    int     denom;
    float   decFract;
    short   turnDirection;
    float   turnDist;
    float randHigh;
    float randLow;
    short Par[3]; //i.e. 0,1,2 in p3 par p0 p1 p2
    int   randToReuse;
    int   IsDecFract :1;
    int   IsRand :1;
    int   IsRandRange :1;
    int IsRandReuse :1;
    int   IsCentre :1;
    int   IsPar;
} *P;

struct ShapeCallStruct { //i.e. a "tri p0,5,6" / "tri :156" line
    char  name[50];
    short pts;
    short pt[MAXNUMPTS];
    float col[3]; //3 colours
    int   ColIsSet :1; //bool, 1 if cols have been set. Init to 0.
    int   ColIsRand :1; //bool, 1 if "colour rand" used. Init to 0.
    int   ColourAddIsUsed :1; //i.e. used in this particular shape call
    int   FlipIsUsed :1;
    int   NotIsUsed :1;
    int   RotIsUsed :1;
    int   Rot2IsUsed :1;
//    int   EqualIsUsed :1;
    int   IsDrawCall :1;
    short waitLevel; //used for later building wait shapes.
    char  waitShapeName[50];
};

struct colInfostruct {
    short ptnum; //number of InitPt
    float col[3]; //the colour this shape gives to the point
};

typedef struct { //stores everything in a shape definition - a "def quad 4" section
    char    name[50];
    short   N; //it's a shape with N vertices
    short   colPts; //for shading
    struct  colInfostruct cInfo[MAXINITPTS]; //for shading
    float   col[3]; //3 colours
    struct {
        float min;
        float max;
        } randRange[3];
    short   NewPts;
    int     ColIsSet :1; //bool, 1 if colours have been set. Init to 0.
    int     ColIsRand :1; //bool, 1 if "colour rand" used. Init to 0.
    int     RandRangeIsUsed :1;
    int     ShadingIsUsed :1;
    int     DrawNow :1; //bool, true if no shape calls - draw now
    int     IsWaitShape :1; 
    int     FlipIsUsed :1;
    int     GradientIsUsed :1; //hmm can use colourInfostruct to store? although its number of shape pt,
                                //not initpt, and always exactly 2 used.
    int     CentreIsUsed :1;
    short   centrePtNum;    //the point number to set to the average of the shape points; "c" means this point.
    struct  NewPtStruct NewPt[MAXPOINTDEFS];
    short   ShapeCalls;
    struct  ShapeCallStruct ShapeCall[MAXSHAPECALLS];
    short   waitLevel;
    short   BranchingIsUsed; //NB not boolean!! 0 means not set. 1 means >, 2 is = ?
    short   firstBranchedSCall;
    short   branchLevel;
} shape;
shape Shape[MAXSHAPEDEFS], *Si, *Ss; //Ss is Shape[s], the current Shape as data structure being built by parser

struct {
    char* name;
    short c[3]; //r,g,b
} col[MAXCOLOURNAMES];

struct {
        char  name[50];
        int num;
} shapesToFlip[MAXSHAPESTOFLIP];  

short numCols=0,currflipshape,numOfShapesToFlip=0;
int numer[100],n,imod; //for p4 p0 1,3,5,6,7/10 p2
short i,j,k,sp,maxLevel,initPts,s=-1,np,lnp,nss,nssp,shadepts,firstnp;
bool NoedgesIsSet, BoxColourIsSet, RandIsUsed, OnlyEdgesIsSet,GradientIsUsed, GridIsSet;
bool CentreIsUsed,ShadingIsUsed,ColourAddIsUsed,MarginArrayIsUsed,DrawIsUsed,FlipIsUsed;
bool OnList;
short ifrom,ito,pfrom,pto,waitnum,shapenum,turnDirection;
short MaxWaitLevel[MAXSHAPEDEFS]; //highest wait shape needed for each shape
float mod,gridTheta,fwidth=0.1,fmargin=1, boxcolour[3], marginArray[4],xsum,ysum,turnDist; //default width and margin
char *cnamestr,*cnamestr2, *ptlist;
char ch,errstr[100],shapename[100]; //no error strings longer than 100!

void Flippify(shape *T, shape *F);

void AddTurnInfo(){    
    lnp=firstnp;
    for (i=lnp;i<=lnp+pto-pfrom;i++){
        Ss->NewPt[i].turnDirection=turnDirection; //-1 for left, 1 for right
        Ss->NewPt[i].turnDist=turnDist;
fprintf(stderr,"AddTurnInfo: lnp %d i %d pto %d pfrom %d turnDir %d turnDist %g\n",lnp,i,pto,pfrom,Ss->NewPt[lnp].turnDirection,Ss->NewPt[lnp].turnDist);
    }
    turnDirection=0; //reset
}

%}

%union{
  int i;
  char*	w;
  float f;
} 

%token EOL DEF LEVELS WIDTH MARGIN DIV DOTS ONLYEDGES CENTRE PAR
%token WAIT GRID NOT EQUAL ROT ROT2 CRAND SRAND CRANDOM GT FLIP
%token NOEDGES BOXCOLOUR DITTO RAND PBOX PGON PCUBE DRAW PLUS MINUS 
%token <i> POINTNUM DOTNUM COMMANUM IPOINTNUM TURN
%token <w> WORD PTSTR
%token <f> NUM

%% 

sdvprogram: level settings optpointdefs defs 
;
level: LEVELS NUM EOL { maxLevel=$2; }
;
//======================================== SETTINGS ==================================
settings: | settings setting
; //settings are optional
setting: srandset | widthset | marginset | noedgesset | boxcolourset | onlyedgesset | gridset
;
widthset: WIDTH NUM EOL {fwidth=$2; }
;
marginset: MARGIN NUM EOL { fmargin=$2; }
| MARGIN NUM NUM EOL {// margin (left, right) (top,bottom)
     marginArray[0]=$2;
     marginArray[1]=$3;
     marginArray[2]=$2;
     marginArray[3]=$3;
     MarginArrayIsUsed=YES;
} | MARGIN NUM NUM NUM NUM EOL { //margin at left, top, right, bottom.
     marginArray[0]=$2;
     marginArray[1]=$3;
     marginArray[2]=$4;
     marginArray[3]=$5;
     MarginArrayIsUsed=YES;
};
srandset: SRAND NUM EOL {srandom($2);}
;
boxcolourset: BOXCOLOUR NUM NUM NUM EOL {
    BoxColourIsSet=1;
    boxcolour[0]=$2;
    boxcolour[1]=$3;
    boxcolour[2]=$4;
} | BOXCOLOUR colourname EOL {
    i=ColNumber(cnamestr); 
    BoxColourIsSet=1;
    boxcolour[0]=col[i].c[0]/255.0; 
    boxcolour[1]=col[i].c[1]/255.0;
    boxcolour[2]=col[i].c[2]/255.0;
};
noedgesset: NOEDGES EOL { NoedgesIsSet=1;}
;
onlyedgesset: ONLYEDGES EOL { OnlyEdgesIsSet=1;}
;
gridset: GRID NUM EOL { GridIsSet=1;
    gridTheta=$2;    }
;
//======================================= POINTDEFS ===================================
optpointdefs: {// if no init pts defined, defaults to "pbox", i.e. to pbox 0 0 100 100
    SetPboxPoints(0,0,100,100);
} | pointdefs
;
pointdefs: pointdef  | pointdefs pointdef //1 or more points?!
;
pointdef: cubedef | pgondef | pboxdef | POINTNUM NUM NUM EOL { 
    if ($1==1 && initPts==0)  //..if no P0 defined when P1 is defined, it defaults to P0 0 0
        SetInitPt(0,0);
    if ($1!=initPts) InitPtError($1);
    SetInitPt($2,$3);
} | POINTNUM CENTRE EOL { //e.g. p4 c
    if ($1!=initPts) InitPtError($1);
    xsum=0;
    ysum=0;
    for (i=0;i<initPts;i++) {
        xsum+=initPt[i].x;
        ysum+=initPt[i].y;
    }
    SetInitPt(xsum/initPts,ysum/initPts);
} | POINTNUM DITTO NUM EOL {
    if ($1!=initPts) InitPtError($1);
    if ($1==0) yyerror("Double quote (ditto) has nothing to repeat error.");
    SetInitPt(initPt[$1-1].x,$3); 
} | POINTNUM NUM DITTO EOL {
    if ($1!=initPts) InitPtError($1);
    if ($1==0) yyerror("Double quote (ditto) has nothing to repeat error.");
    SetInitPt($2,initPt[$1-1].y);
};
pboxdef: PBOX NUM NUM EOL {
    SetPboxPoints(0,0,$2,$3);
} | PBOX EOL {// defaults to pbox 0 0 100 100
    SetPboxPoints(0,0,100,100);
}  | PBOX NUM NUM NUM NUM EOL { //..pbox a,b c,d - make a rect of 4 points 
    SetPboxPoints($2,$3,$4,$5);
};
cubedef: PCUBE NUM NUM NUM NUM NUM NUM EOL {//cube 50 15 -5 0 -1000 -500
    struct pt {
        float x;
        float y;
        float z;
    } p[8], new[8], rot, xy[8];
    float DEGTORAD=2*PI/360.0;
    float rad=$2, eyez=$6, planez=$7;
    rot.x=$3*DEGTORAD;
    rot.y=$4*DEGTORAD;
    rot.z=$5*DEGTORAD;
    p[0].x=-rad; p[0].y=rad; p[0].z=-rad;
    p[1].x=-rad; p[1].y=rad; p[1].z=rad;
    p[2].x=rad; p[2].y=rad; p[2].z=rad;
    p[3].x=rad; p[3].y=rad; p[3].z=-rad;
    p[4].x=-rad; p[4].y=-rad; p[4].z=-rad;
    p[5].x=-rad; p[5].y=-rad; p[5].z=rad;
    p[6].x=rad; p[6].y=-rad; p[6].z=rad;
    p[7].x=rad; p[7].y=-rad; p[7].z=-rad;
    for (i=0;i<8;i++){
    //rotation about x axis
        new[i].x=p[i].x;
        new[i].y=cos(rot.x)*p[i].y-sin(rot.x)*p[i].z;
        new[i].z=sin(rot.x)*p[i].y+cos(rot.x)*p[i].z;
        p[i]=new[i];
    //rotation about y axis
        new[i].x=cos(rot.y)*p[i].x+sin(rot.y)*p[i].z;
        new[i].y=p[i].y;
        new[i].z=-sin(rot.y)*p[i].x+cos(rot.y)*p[i].z;
        p[i]=new[i];
    //rotation about z axis
        new[i].x=cos(rot.z)*p[i].x-sin(rot.z)*p[i].y;
        new[i].y=sin(rot.z)*p[i].x+cos(rot.z)*p[i].y;
        new[i].z=p[i].z;
    }
    //view from 0,0,eyez, viewplane is z=planez
    for (i=0;i<8;i++){
        if (new[i].z-eyez==0) yyerror("Divide by zero in pcube: a cube vertex is on the viewplane.");
        xy[i].x=(planez-eyez)*new[i].x/(new[i].z-eyez);
        xy[i].y=(planez-eyez)*new[i].y/(new[i].z-eyez);
        SetInitPt(xy[i].x,xy[i].y);
    }
};
pgondef: PGON NUM NUM NUM EOL {//pgon 5 50 0 : num of sides, radius, angular location of p0.
    for (i=0;i<$2;i++)
        SetInitPt($3*cos(2*PI*($4/360+i/$2)),$3*sin(2*PI*($4/360+i/$2)));
} | PGON NUM NUM EOL {//pgon 5 50 : num of sides, radius. p0 is at 0 degrees.
    for (i=0;i<$2;i++)
        SetInitPt($3*cos(2*PI*(i/$2)),$3*sin(2*PI*(i/$2)));
};

//=========================================== SHAPE DEFS =================================
defs: namedef | defs namedef 
;
namedef: defline optionals draworshapesets 
; 
defline: DEF WORD NUM EOL {
    s++;
    Ss=&Shape[s];
    nss=0; //reset num of subdiv shapes
    //Init shape def now, at its first use.
    Ss->col[0]=1; //white
    Ss->col[1]=1; 
    Ss->col[2]=1; 
    strcpy(Ss->name,$2);   //e.g. tri
    Ss->N=$3;      //e.g. 3
    //if it's the first shape, it's the main shape (first one subdivided)
    if (s==0)  
        if (initPts<Shape[0].N)
            yyerror("Not enough initial points for main shape.");
    np=Ss->N; // first new point num is same num as num of points in shape.
    P=&Ss->NewPt[np]; //i.e. P is the current new point in the current shape being parsed
    turnDirection=0;
};
//========================================== SHAPE OPTIONAL SETTINGS ========================
optionals: | optionals option
;
option: gradientline | shadingline | colourline | pointset// | turnpointset
//hmm but no - only 1 of some of these should be allowed (e.g. colour), many of others (e.g. ip0 red)
;
gradientline: POINTNUM NUM NUM NUM POINTNUM NUM NUM NUM EOL {
    SetGradientPtInfo(0,$1,$2,$3,$4);
    SetGradientPtInfo(1,$5,$6,$7,$8);
} | POINTNUM colourname POINTNUM colourname2 EOL {
    i=ColNumber(cnamestr);
    SetGradientPtInfo(0,$1,col[i].c[0]/255.0,col[i].c[1]/255.0,col[i].c[2]/255.0);
    i=ColNumber(cnamestr2);
    SetGradientPtInfo(1,$3,col[i].c[0]/255.0,col[i].c[1]/255.0,col[i].c[2]/255.0);
};
colourline: NUM NUM NUM EOL {
    SetShapeCol($1,$2,$3);
} | RAND EOL {
    Ss->ColIsRand=1;
    RandIsUsed=1;
    Ss->ColIsSet=1; //maybe not the best name for it now, as here colour is not set yet, but
    //there is a shape colour.
}| RAND NUM NUM NUM NUM NUM NUM EOL {// rand .2 .5 .3 .4 .6 .7
 //   Ss->ColIsRand=1;
    Ss->RandRangeIsUsed=1;
    Ss->ColIsSet=1;
    Ss->randRange[0].min=$2;
    Ss->randRange[0].max=$3;
    Ss->randRange[1].min=$4;
    Ss->randRange[1].max=$5;
    Ss->randRange[2].min=$6;
    Ss->randRange[2].max=$7;
} | colourname EOL {
    i=ColNumber(cnamestr); //cname sets cnamestr to colour name
    SetShapeCol(col[i].c[0]/255.0,col[i].c[1]/255.0,col[i].c[2]/255.0);
};
shadingline: shadings EOL
;
shadings: shading | shadings shading
;
shading: IPOINTNUM NUM NUM NUM { //ip0 1 1 0
    SetShadingCol($1,$2,$3,$4);
} | IPOINTNUM colourname { //ip0 light grey
    i=ColNumber(cnamestr); //cname sets cnamestr to colour name
    SetShadingCol($1,col[i].c[0]/255.0,col[i].c[1]/255.0,col[i].c[2]/255.0);
} | IPOINTNUM DOTNUM NUM NUM NUM { //ip0..3 1 1 0
    for (i=$1;i<=$2;i++)
        SetShadingCol(i,$3,$4,$5);
} | IPOINTNUM DOTNUM colourname { //ip0..3 red
    j=ColNumber(cnamestr); 
    for (i=$1;i<=$2;i++)
        SetShadingCol(i,col[j].c[0]/255.0,col[j].c[1]/255.0,col[j].c[2]/255.0);
};
colourname: WORD {cnamestr=$1;} 
    | colourname WORD {
        strcat(cnamestr," ");
        strcat(cnamestr,$2);
};
colourname2: WORD {cnamestr2=$1;} //for the second colourname on a line, used in gradient shading
    | colourname2 WORD {
        strcat(cnamestr2," ");
        strcat(cnamestr2,$2);
};
//================================= NEW SHAPE POINTS ==================================================
pointset:  singlepoint | multiplepoints
;
singlepoint: plainfract maybeturn EOL
| POINTNUM POINTNUM RAND POINTNUM EOL { //p3 p1 rand p2
    CheckAndSetPoints($1,$2,$4);
    P->IsRand=1; 
    RandIsUsed=1;
    IncPointVars();
} | POINTNUM POINTNUM CRAND POINTNUM EOL { //p3 p1 crand p2 - random number from 0 <=  <= 1
    CheckAndSetPoints($1,$2,$4);
    P->decFract=(random()%1001)/1000.0; 
    P->IsDecFract=1; 
    IncPointVars();
} | POINTNUM POINTNUM RAND NUM NUM POINTNUM EOL { //p3 p1 rand .1 .7 p2 ==> rand range
    if ($4>$5) yyerror("Low end of range is higher than high end");
    CheckAndSetPoints($1,$2,$6);
    P->IsRandRange=1;
    P->randLow=$4;
    P->randHigh=$5;
    RandIsUsed=1;
    IncPointVars();
} | POINTNUM POINTNUM CRAND NUM NUM POINTNUM EOL { //p3 p1 crand .1 .7 p2 ===> random number from .1 <=  <= .7
    if ($4>$5) yyerror("Low end of range is higher than high end");
    CheckAndSetPoints($1,$2,$6);
    mod=($5-$4)*1000.0;
    imod=mod;
    P->decFract=$4 + (( random()%imod ) + 1.0) / 1000.0; 
    fprintf(stderr,"$5 %g $4 %g mod %g imod %d decFract %g\n",$5,$4,mod,imod,P->decFract);
    P->IsDecFract=1; 
    IncPointVars();
} | plaindecfract maybeturn EOL
| POINTNUM POINTNUM DITTO POINTNUM EOL {//p3 p1 " p2 (with 0.6 above)
    if (Ss->NewPts==0) yyerror("Double quote (ditto) has nothing to repeat error.");
    CheckAndSetPoints($1,$2,$4);
    if (Ss->NewPt[np-1].IsDecFract){ 
        P->decFract=Ss->NewPt[np-1].decFract; //assuming np-1 is prev point
        P->IsDecFract=1; 
    } else if (Ss->NewPt[np-1].IsRandRange || Ss->NewPt[np-1].IsRand || Ss->NewPt[np-1].IsRandReuse){ 
    //ie prev line is p3 p1 rand .1 .7 p2 or p3 p1 rand p2 or p3 p1 " p2 with rand before the 1 or more dittos.
        P->IsRandReuse=1;
    } else {//p3 p1 " p2 (with 3/5 above)
        P->numer=Ss->NewPt[np-1].numer;
        P->denom=Ss->NewPt[np-1].denom;
    }
    IncPointVars();
} | POINTNUM CENTRE EOL { //p3 c
    if (np!=$1) yyerror("Bad new point number in new point centre definition"); 
    CentreIsUsed=YES;
    P->IsCentre=YES;
    IncPointVars();
} | POINTNUM PAR POINTNUM POINTNUM POINTNUM EOL { // p3 par p0 p1 p2
    P->IsPar=YES;
    P->Par[0]=$3;
    P->Par[1]=$4;
    P->Par[2]=$5;
    IncPointVars();
}; 
plaindecfract: POINTNUM POINTNUM NUM POINTNUM { //p3 p1 0.6 p2
    CheckAndSetPoints($1,$2,$4);
    P->decFract=$3; 
    P->IsDecFract=1; 
    IncPointVars();
} 
plainfract: POINTNUM POINTNUM NUM DIV NUM POINTNUM { // p3 p1 1/2 p2
    CheckAndSetPoints($1,$2,$6);
    P->numer=$3; // can't I work out the number now? than print error if not good. dont need numer&denom separate..
    P->denom=$5;
    IncPointVars();
};
//============================== ADDING MULTIPLE POINTS =========================
mayberangeendpoint: DOTS | DOTNUM { //.. or ..0 at end of  e.g. p4..7 p0.. .5 p1..0 <--
    Ss->NewPt[np-1].To=$1;
};
maybeturn: | TURN NUM { // left 0.5
    Ss->NewPt[np-1].turnDirection=$1; //-1 for left, 1 for right
    Ss->NewPt[np-1].turnDist=$2;
};
loopedmaybeturn: | TURN NUM { // left 0.5
    turnDirection=$1;
    turnDist=$2;
};
loopeddecfract: POINTNUM DOTNUM POINTNUM DOTS NUM POINTNUM {//p5..8 p0.. .95 p1..
    pfrom=$1;
    pto=$2;
///new error message, put elsewhere, make function etc. This checks if upper range is > lower range
    if ($2<=$1) {
        sprintf(errstr,"Bad new point range: p%d..%d ??",$1,$2);
        yyerror(errstr);
    }
    CheckPoints($1,$3,$6);
    firstnp=$1; //to use in loopedmaybeturn later
    for (i=0;i<=pto-pfrom;i++)
        SetDecFract($3+i,$6+i,$5);
};
loopedfract:POINTNUM DOTNUM POINTNUM DOTS NUM DIV NUM POINTNUM { // p3..5 p1.. 1/2 p2 (maybe .. or ..n next)
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$8);
    firstnp=$1; //to use in loopedmaybeturn later
    for (i=0;i<=pto-pfrom;i++){
        SetFromAndToPts($3+i,$8+i);
        P->numer=$5; // can't I work out the number now? than print error if not good. dont need numer&denom separate..
        P->denom=$7;
        IncPointVars();
    }
} 
loopedcrand:POINTNUM DOTNUM POINTNUM DOTS CRAND POINTNUM{//p5..8 p0.. crand p1 (maybe .. or ..n next)
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$6);
    firstnp=$1; //to use in loopedmaybeturn later
    for (i=0;i<=pto-pfrom;i++)
        SetDecFract($3+i,$6+i,(rand()%1001)/1000.0);
}
loopedrand:POINTNUM DOTNUM POINTNUM DOTS RAND POINTNUM{//p5..8 p0.. rand p1 (maybe .. or ..n next)
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$6);
    firstnp=$1; //to use in loopedmaybeturn later
    for (i=0;i<=pto-pfrom;i++){
        SetFromAndToPts($3+i,$6+i);
        P->IsRand=1; 
        IncPointVars();
    }
}
loopedtype: loopedfract | loopedcrand | loopedrand | loopeddecfract
;
multiplepoints: loopedtype mayberangeendpoint loopedmaybeturn EOL {
    if (turnDirection!=0)
        AddTurnInfo();
} | POINTNUM DOTNUM POINTNUM DOTS NUM DIV NUM POINTNUM loopedmaybeturn EOL { // p3..5 p1.. 1/2 p2 
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$8);
    firstnp=$1; //to use in loopedmaybeturn later
    for (i=0;i<=pto-pfrom;i++){
        SetFromAndToPts($3+i,$8);
        P->numer=$5; // can't I work out the number now? than print error if not good. dont need numer&denom separate..
        P->denom=$7;
        IncPointVars();
    }
    if (turnDirection!=0)
        AddTurnInfo();
} | POINTNUM DOTNUM POINTNUM NUM DOTNUM DIV NUM POINTNUM loopedmaybeturn EOL{//p3..6 p0 1..4/5 p1
    ifrom=$4;
    ito=$5;
    CheckPoints($1,$3,$8);
    pfrom=$1;
    pto=$2;
    firstnp=$1; //to use in loopedmaybeturn later
    if (pto-pfrom!=ito-ifrom) yyerror("Range mismatch between point numbers and fraction denominators");
    for (i=0;i<=ito-ifrom;i++){
        SetFromAndToPts($3,$8);
        P->numer=ifrom+i; // can't I work out the number now? than print error if not good. dont need numer&denom separate..
        P->denom=$7;
            fprintf(stderr,"i %d from %d to %d numer %d denom %d\n",i,P->From,P->To,P->numer,P->denom);

        IncPointVars();
    }
    if (turnDirection!=0)
        AddTurnInfo();
} | POINTNUM DOTNUM POINTNUM DOTS CRANDOM POINTNUM DOTS EOL {//p5..8 p0.. crandom p1..
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$6);
    srandomdev();
    for (i=0;i<=pto-pfrom;i++)
        SetDecFract($3+i,$6+i,(random()%1001)/1000.0);
} | POINTNUM DOTNUM POINTNUM DOTS NUM POINTNUM EOL {//p5..8 p1.. .95 p0
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$6);
    for (i=0;i<=pto-pfrom;i++)
        SetDecFract($3+i,$6,$5);
} | POINTNUM DOTNUM POINTNUM NUM POINTNUM DOTS EOL {//p5..8 p0 .95 p1..
    pfrom=$1;
    pto=$2;
    CheckPoints($1,$3,$5);
    for (i=0;i<=pto-pfrom;i++)
        SetDecFract($3,$5+i,$4);
}| POINTNUM DOTNUM POINTNUM commasepnumers DIV NUM POINTNUM EOL{//p3..6 p0 1,3,5,8/9 p1
//    ifrom=$4;
//    ito=$5;
    CheckPoints($1,$3,$7);
    pfrom=$1;
    pto=$2;
    fprintf(stderr,"pto %d pfrom %d n %d\n",pto,pfrom,n);
    if (pto-pfrom!=n-1) yyerror("Mismatch between point range and fraction numerators");
    for (i=0;i<n;i++){
        SetFromAndToPts($3,$7);
        P->numer=numer[i]; // can't I work out the number now? than print error if not good. dont need numer&denom separate..
        P->denom=$6;
        fprintf(stderr,"numer %d denom %d i %d\n",P->numer,P->denom,i);
        IncPointVars();
    }
};
commasepnumers: numer commanumers
;
numer: NUM {
    n=0;
    numer[n]=$1;
    n++;
};
commanumers: cnumer | commanumers cnumer
;
cnumer: COMMANUM {
    numer[n]=$1;
    n++;
};
//============================================ SHAPE CALLS ========================================
draworshapesets: | shapesets //| {Ss->DrawNow=1; }
//| shapesets //| draw
;
//draw: DRAW EOL { Ss->DrawNow=1; } // draw
//;
shapesets: shapeset | shapesets shapeset
;
shapeset: WORD verts colourset optflip optwait EOL { //sq p0,4,2,5 red
    if (waitnum==0) {
        strcpy(Ss->ShapeCall[nss].name,$1);
//        printf("%d %d %s %s %d\n",s,nss,$1,"draw",strcmp($1,"draw"));
        if (strcmp($1,"draw")==0) {
            DrawIsUsed=YES;
            Ss->ShapeCall[nss].IsDrawCall=YES;
        }
    }
    else {
        sprintf(shapename,"%s-wait-%d",$1,waitnum);
        strcpy(Ss->ShapeCall[nss].name,shapename);
        Ss->ShapeCall[nss].waitLevel=waitnum; //NB change this to only update if waitnum>current level.
        strcpy(Ss->ShapeCall[nss].waitShapeName,$1);
    }
    Ss->ShapeCalls++;
    nss++;
    nssp=0;//reset for next shape
    waitnum=0; //reset for next shape
} //"draw" EOL { printf("DRAW\n");
//} 
| GT NUM EOL { //branching, e.g. >3, >-3
    Ss->BranchingIsUsed=GREATERTHAN;
    Ss->firstBranchedSCall=nss;
    Ss->branchLevel=$2;
} | EQUAL NUM EOL { //branching, e.g. =3, =-3
    Ss->BranchingIsUsed=EQUALS;
    Ss->firstBranchedSCall=nss;
    Ss->branchLevel=$2;
};
optflip: | flip
;
flip: FLIP { FlipIsUsed=1;
    Ss->FlipIsUsed=1;
    Ss->ShapeCall[nss].FlipIsUsed=1;
};
optwait: | waitn | wait
;
wait: WAIT { waitnum=1;
};
waitn: WAIT NUM { waitnum=$2;
};
verts: dotverts | commasepvertices | pointstr
;
pointstr: PTSTR {// :230 as in sq :230 (or) sq : 230
    ptlist=$1;
    for (i=0;i<strlen(ptlist);i++){
        ch = ptlist[i];
        if (ch=='c') {//'c' for centre
            if (!Ss->CentreIsUsed) {//add new point for centre
        //so - make a NewPt for the shape, store the number of the centre Pt,
                Ss->CentreIsUsed=YES; //is this flag needed?
                CentreIsUsed=YES;
                Ss->centrePtNum=np;
                P->IsCentre=YES;
                IncPointVars();
            }
            Ss->ShapeCall[nss].pt[nssp]=Ss->centrePtNum;
        } else
            Ss->ShapeCall[nss].pt[nssp]=ch-'0';
        CheckPtNumber();
        IncShapePointVars();
    }
};
dotverts : POINTNUM DOTNUM { //p0..5
    for (i=$1;i<=$2;i++){
        Ss->ShapeCall[nss].pt[nssp]=i;
        IncShapePointVars();
    }
};
colourset: | NUM NUM NUM { // 1 0.5 .2
    Ss->ShapeCall[nss].ColIsSet=1;
    SetShapeCallCol($1,$2,$3);
} | RAND {
    Ss->ShapeCall[nss].ColIsRand=1;
    RandIsUsed=1;
} | colourname {
    j=ColNumber(cnamestr); //cname sets cnamestr to colour name
    Ss->ShapeCall[nss].ColIsSet=1;
    SetShapeCallCol(col[j].c[0]/255.0,col[j].c[1]/255.0,col[j].c[2]/255.0);
} | PLUS NUM NUM NUM { //Colour addition   + 0.1 0.1 -1
    ColourAddIsUsed=1;
    Ss->ShapeCall[nss].ColourAddIsUsed=1;
    SetShapeCallCol($2/10.0,$3/10.0,$4/10.0);
} | MINUS NUM NUM NUM { //Colour subtraction   - 0.1 0.1 -1 NB must be " - "
    ColourAddIsUsed=1;
    Ss->ShapeCall[nss].ColourAddIsUsed=1;
    SetShapeCallCol(-$2/10.0,-$3/10.0,-$4/10.0);
} | NOT { Ss->ShapeCall[nss].NotIsUsed=1;
} | ROT { Ss->ShapeCall[nss].RotIsUsed=1;
} | ROT2 { Ss->ShapeCall[nss].Rot2IsUsed=1;
//} | EQUAL{ Ss->ShapeCall[nss].EqualIsUsed=1;
}; //colourset is optional;
commasepvertices: vertex commavertices
;
vertex: POINTNUM {
    Ss->ShapeCall[nss].pt[nssp]=$1;
    CheckPtNumber();
    IncShapePointVars();
};
commavertices: cvertex | commavertices cvertex
;
cvertex: COMMANUM {
    Ss->ShapeCall[nss].pt[nssp]=$1;
    CheckPtNumber();
    IncShapePointVars();
};
%% 
//======================================== C ===================================
bool IsOnFlipList(char sh[50]){
    OnList=NO;
    for (k=0;k<numOfShapesToFlip;k++) {
        fprintf(stderr,"Comparing %s and %s\n",sh,shapesToFlip[k].name);
        if (strcmp(sh,shapesToFlip[k].name)==0) {
            fprintf(stderr,"SAME!\n");
            OnList=YES;
        }
    }
    return OnList;
}

void Flippify(shape *T, shape *F){
//makes one shape into 2, making a new flipped version, and changing the shapes called
//so that shape calls with FLIP and those without are opposite in the -flipped shape.
//also left & right in new pts are swapped.
    int i,j;
    struct ShapeCallStruct *Ti, *Fi;
    //no! put -flipped at end of name
    strcpy(T->name,F->name);
    strcat(T->name,"-flipped");
    T->N=F->N;
    T->colPts=F->colPts;
    for (i=0;i<F->colPts;i++){
        T->cInfo[i].ptnum=F->cInfo[i].ptnum;
        for (j=0;j<3;j++)
            T->cInfo[i].col[j]=F->cInfo[i].col[j];
    }
    for (i=0;i<3;i++)
        T->col[i]=F->col[i];
    T->NewPts=F->NewPts;
    T->ColIsSet=F->ColIsSet;
    T->ColIsRand=F->ColIsRand;
    T->ShadingIsUsed=F->ShadingIsUsed;
    T->DrawNow=F->DrawNow;
    T->IsWaitShape=F->IsWaitShape;
    T->FlipIsUsed=F->FlipIsUsed;
    T->GradientIsUsed=F->GradientIsUsed;
    for (i=F->N;i<F->N+F->NewPts;i++) {
        T->NewPt[i].From=F->NewPt[i].From;
        T->NewPt[i].To=F->NewPt[i].To;
        T->NewPt[i].numer=F->NewPt[i].numer;
        T->NewPt[i].denom=F->NewPt[i].denom;
        T->NewPt[i].decFract=F->NewPt[i].decFract;
        T->NewPt[i].turnDirection=-F->NewPt[i].turnDirection; //NB!! -1 x the old value.
        //i.e. swap left <--> right. "no turn" (turnDirection=0) stays the same.
        T->NewPt[i].turnDist=F->NewPt[i].turnDist;
        T->NewPt[i].IsDecFract=F->NewPt[i].IsDecFract;
        T->NewPt[i].IsRand=F->NewPt[i].IsRand;
        T->NewPt[i].IsCentre=F->NewPt[i].IsCentre;
    }
    T->ShapeCalls=F->ShapeCalls;
    for (i=0;i<F->ShapeCalls;i++) {
        Ti=&T->ShapeCall[i];
        Fi=&F->ShapeCall[i];
        Ti->FlipIsUsed=!Fi->FlipIsUsed; //NB!! the opposite
        //if flip is used in either, add "-flipped" to 1 of them
        strcpy(Ti->name,Fi->name);
        if (Fi->FlipIsUsed) //call flipped version in one of the 2 shape calls.
            strcat(Fi->name,"-flipped");
        else
            strcat(Ti->name,"-flipped");
        Ti->pts=Fi->pts;
        for (j=0;j<Fi->pts;j++)
            Ti->pt[j]=Fi->pt[j];
        for (j=0;j<3;j++) 
            Ti->col[j]=Fi->col[j];
        Ti->ColIsSet=Fi->ColIsSet;
        Ti->ColIsRand=Fi->ColIsRand;
        Ti->ColourAddIsUsed=Fi->ColourAddIsUsed;
        Ti->NotIsUsed=Fi->NotIsUsed;
        Ti->RotIsUsed=Fi->RotIsUsed;
//        Ti->EqualIsUsed=Fi->EqualIsUsed;
        Ti->IsDrawCall=Fi->IsDrawCall;
        Ti->waitLevel=Fi->waitLevel;
        strcpy(Ti->waitShapeName,Fi->waitShapeName);
    }
    T->waitLevel=F->waitLevel;
    T->BranchingIsUsed=F->BranchingIsUsed;
    T->firstBranchedSCall=F->firstBranchedSCall;
    T->branchLevel=F->branchLevel;
}

int main(int argc, char **argv) 
{ 
//** to debug, uncomment this and make with "bison -d subdiv.y --debug" instead of
//** the usual "bison -d subdiv.y" in the Makefile  (now can use "make debug"):
//  yydebug=1; 
    float minx,maxx,miny,maxy;
    SetColNames();
    yyparse(); 
    minx=initPt[0].x; //calc BoundingBox
    maxx=initPt[0].x;
    miny=initPt[0].y;
    maxy=initPt[0].y;
    for (i=1;i<initPts;i++) {
        if (initPt[i].x<minx) minx=initPt[i].x;
        if (initPt[i].x>maxx) maxx=initPt[i].x;
        if (initPt[i].y<miny) miny=initPt[i].y;
        if (initPt[i].y>maxy) maxy=initPt[i].y;
    }
    if (MarginArrayIsUsed){
        minx-=marginArray[0];
        maxy+=marginArray[1];
        maxx+=marginArray[2];
        miny-=marginArray[3];
    }
    else
        if (fmargin!=0) {
            minx-=fmargin;
            maxx+=fmargin;
            miny-=fmargin;
            maxy+=fmargin;
        }
    //in colour add mode, first shape is mid grey by default.
    if (ColourAddIsUsed && !Shape[0].ColIsSet){
        Shape[0].ColIsSet=YES;
        Shape[0].col[0]=0.5;
        Shape[0].col[1]=0.5;
        Shape[0].col[2]=0.5;
    }
    //Deal with flips.
    if (FlipIsUsed) {
        currflipshape=0;
        //could just loop this bit, until no new shapes added.
        for (i=0;i<=s;i++) {
            if (Shape[i].FlipIsUsed) 
                for (j=0;j<Shape[i].ShapeCalls;j++) 
                    if (Shape[i].ShapeCall[j].FlipIsUsed) {
                    //if not on shapesToFlip list, add to list
                        fprintf(stderr,"shape call flip: %d %d\n",i,j);
                        OnList=IsOnFlipList(Shape[i].ShapeCall[j].name);
                        if (!OnList) {
                            strcpy(shapesToFlip[numOfShapesToFlip].name,Shape[i].ShapeCall[j].name);
                            fprintf(stderr,"Added: %s\n",Shape[i].ShapeCall[j].name);
                            numOfShapesToFlip++;
                        }
                    }
            //now LINE is on list... go through its shape calls, add any to list
            //that arent already on it.
            
            //NB its this next bit that needs to be looped.
            int currshape;
            for (k=0;k<numOfShapesToFlip;k++) {
                currshape=GetShapeNum(shapesToFlip[k].name);
                for (i=0;i<Shape[currshape].ShapeCalls;i++) 
                    if (!IsOnFlipList(Shape[currshape].ShapeCall[i].name)) {
                        strcpy(shapesToFlip[numOfShapesToFlip].name,Shape[currshape].ShapeCall[i].name);
                        fprintf(stderr,"Added: %s\n",Shape[currshape].ShapeCall[i].name);
                        numOfShapesToFlip++;
                    }
            }
            fprintf(stderr,"Flip shapes added: %d\n",numOfShapesToFlip);
        }
        //that will handle maybe NOT ALL shapes that needs flips!
        //I can imagine programs it wont do properly.. shapes without flips, calling
        //shapes without flips, but THOSE call shapes with flips. Better to loop
        //and handle all possible programs with flips. (NOT DONE YET)
        
        //now just make flipped version of all shapes on the list, altering the originals
        //at the same time... :-)
        fprintf(stderr,"Num of shapes to flip: %d\n",numOfShapesToFlip);
        //now with list, loop through it and if any shape-flipped doesnt exist, make it.
        //and SomethingDone=YES. 
        for (i=0;i<numOfShapesToFlip;i++) {
        //so, need num of shape to copy..GetShapeNum(shapesToFlip[i].name)
        //and num of shape to copy to..
            s++; //now s is num of new shape.
            Flippify(&Shape[s],&Shape[GetShapeNum(shapesToFlip[i].name)]);
        }
        //if SomethingDone, do everything again. otherwise continue program.
    }
    BuildNewWaitShapes();
    // *** variable s ends up the number of the highest shape,
    // *** e.g. if 3 shapes, s=2, for Shape[0,...,2]
    for (i=0;i<=s;i++)
        if (Shape[i].ShapeCalls==0)
            Shape[i].DrawNow=1;
    MakeLogFiles();
    printf("%%!PS-Adobe-3.0 EPSF-3.0\n"); //just %! or %!PS-Adobe-3.0 dont work with GS.
    printf("%%%%BoundingBox: %g %g %g %g\n",minx,miny,maxx,maxy);
    printf("2 setlinejoin\n");  
    if (RandIsUsed)  //otherwise random numbers are the same every time
        printf("%u srand rand srand\n",(unsigned)time(NULL)%1000); 
    printf("/gr {grestore} bind def\n/gs {gsave} bind def\n"
           "/xy {aload pop} bind def\n" "/x {0 get} bind def\n/y {1 get} bind def\n"); 
    if (initPts<2) 
        yyerror("Not enough initial points defined"); //but maybe defining 1 would be useful...
    for (i=0;i<initPts;i++)
        printf("/P%d [%g %g] def\n",i,initPt[i].x,initPt[i].y);
    printf("0 0 0 setrgbcolor\n%g setlinewidth\n",fwidth);
    for (i=0;i<=s;i++) {
        Si=&Shape[i];
        if (Si->ColIsRand) {
            printf("/colour%d {",i);
            printf("rand 2147483647.0 div rand 2147483647.0 div rand 2147483647.0 div " 
                   "setrgbcolor} def\n");
        }
        else
            if (Si->ColIsSet) {
                printf("/colour%d {",i);
                printf("%g %g %g ",Si->col[0],Si->col[1],Si->col[2]); 
                printf("setrgbcolor} bind def\n");
            }
        if (ColourAddIsUsed && i==0) {
            printf("/colour0Array [" "%g %g %g ",Si->col[0],Si->col[1],Si->col[2]); 
            printf("] def\n");
        }
    }
    printf("/MAXLEVEL %d def\n",maxLevel);
    if (CentreIsUsed)
        PrintCentreProcedure();
    if (ShadingIsUsed)
        PrintShadingColProcedure();
    if (DrawIsUsed)
        PrintDrawShape();
    for (i=0;i<=s;i++)  //print the subdiv- procedures
        PrintShapeProcedure(i);
    printf("\n%%main\ngs\n");
    if (BoxColourIsSet) {
        printf("gs\nnewpath\n" "%g %g moveto\n",minx,miny);
        printf("%g %g lineto\n",maxx,miny);
        printf("%g %g lineto\n",maxx,maxy);
        printf("%g %g lineto\n",minx,maxy);
        printf("closepath\n" "%g %g %g setrgbcolor\n",boxcolour[0],boxcolour[1],boxcolour[2]);
        printf("fill\ngr\n");
    }
    if (Shape[0].ColIsSet) {
        if (Shape[0].ColIsRand)
            printf("rand 2147483647.0 div rand 2147483647.0 div rand 2147483647.0 div\nP0");
        else
            printf("%g %g %g P0",Shape[0].col[0],Shape[0].col[1],Shape[0].col[2]);
    } else
        printf("-101 0 0 P0");
    for (i=1;i<Shape[0].N;i++)
        printf(" P%d",i);
    printf(" -1\nsubdiv-%s\ngr\n",Shape[0].name);
    return 0;
} //main

int yyerror(char* s) {
    extern int yylineno;	// defined and maintained in lex.c
//    extern char *yytext;	// defined and maintained in lex.c
    fprintf (stderr,"ERROR on line %d: %s\n",yylineno,s);
    exit(1);
}

void PrintDrawShape() {
    printf("\n/draw { %% ( colR colG colB [pts] -- )\n");
    printf("    /pts exch def\n" "    /colB exch def\n");
    printf("    /colG exch def\n" "    /colR exch def\n");
    printf("    %%draw sq\n" "    newpath\n");
    printf("    pts 0 get xy moveto\n");
    printf("    1 1 pts length 1 sub {\n");
    printf("    %%or instead of these 2 lines: pts swap i get xy lineto\n");
    printf("            /i exch def \n" "            pts i get xy lineto\n");
    printf("    } for\n" "    closepath\n");
    if (!OnlyEdgesIsSet) {
        printf("    gs\n" "    colR -100 gt\n");
        printf("    {colR colG colB setrgbcolor} if\n" "    fill\n" "    gr\n");
    }
    if (!NoedgesIsSet)
        printf("    stroke\n");
    printf("} bind def\n");
}

void PrintNewShapePoint(short j){
    short fromMultiplier;
    P=&Si->NewPt[j]; // P is the jth new point defined in the ith Shape.
    if (P->IsPar) {
        printf("\t/P%d [ P%d x P%d x sub P%d x add\n",j,P->Par[1],P->Par[0],P->Par[2]);
        printf("\t\tP%d y P%d y sub P%d y add\n",P->Par[1],P->Par[0],P->Par[2]);
        printf("] def\n");
        return;
    }
    if (P->IsDecFract) { //e.g.   p5 p3 .6 p4
        printf("\t/P%d [ P%d x P%d x sub %g ",j,P->To,P->From,P->decFract);
        printf("mul P%d x add\n",P->From);
        if (P->turnDirection)  //turns left or right off line : do x coord
            PrintTurnX();                            
        printf("\t\tP%d y P%d y sub %g ",P->To,P->From,P->decFract);
        printf("mul P%d y add\n",P->From);
        if (P->turnDirection)  //turns left or right off line : do x coord
           PrintTurnY();                            
        printf("\t\t] def\n");
        return;
    }
    if (P->IsRandRange || P->IsRand || P->IsRandReuse) {//p5 p3 rand .3 .6 p4 OR //p5 p3 rand p4 OR p5 p3 " p4 with prev rand
        if (P->IsRandReuse)
            printf("        /r%d r%d def\n",j,j-1);
        else if (P->IsRandRange)
            printf("        /r%d rand 2147483647.0 div %g %g sub mul %g add def\n",j,P->randHigh,P->randLow,P->randLow);
        else
            printf("        /r%d rand 2147483647.0 div def\n",j);
        printf("        /P%d [ P%d x ",j,P->From);
        printf("1 r%d sub mul P%d x ",j,P->To);
        printf("r%d mul add\n",j);
        printf("            P%d y ",P->From);
        printf("1 r%d sub mul P%d y ",j,P->To);
        printf("r%d mul add ] def\n",j);
        return;
    }
    if (P->IsCentre) {// p4 c
        printf("       /P%d  [ [ ",j);
        for (k=0;k<Si->N;k++)
            printf("P%d ",k);
        printf(" ] find-centre\n" "            ]  def\n");
        return;
    }
//so then by elimination, must be a fraction. e.g. p5 p3 1/2 p4
    fromMultiplier=P->denom - P->numer;
//if new point 1/3 of the way from P0 to P2, it's (2xP0+P2)/3.
//so the 2 is the fromMultiplier.
    printf("        /P%d [ P%d x ",j,P->From);
     if (fromMultiplier!=1) //ie fraction isnt of form x/(x+1) like 1/2 or 2/3
        printf("%d mul ",fromMultiplier);
    printf("P%d x ",P->To);
    if (P->numer!=1) //not 1/x
       printf("%d mul ",P->numer);
    printf("add %d div\n",P->denom);
    if (P->turnDirection)  //turns left or right off line : do x coord
        PrintTurnX();                            
    printf("            P%d y ",P->From);
    if (fromMultiplier!=1)
        printf("%d mul ",fromMultiplier);
    printf("P%d y ",P->To);
    if (P->numer!=1)
        printf("%d mul ",P->numer);
    printf("add %d div\n",P->denom);
    if (P->turnDirection)  //turns left or right off line: do y coord
        PrintTurnY();                            
    printf("\t    ] def\n");
} //end PrintNewShapePoint()

void PrintNewShapePoints() {
    short start = Si->N;
    for (j=start;j<Si->NewPts+start;j++) 
        PrintNewShapePoint(j);
}

void PrintShapeCall(short j) {
//SC points to the jth shape call in Shape i.   
    struct ShapeCallStruct *SC = &Si->ShapeCall[j];
    printf("        %%stack for ");
    printf(SC->IsDrawCall ? "%s\n" : "subdiv-%s\n",SC->name); //i.e. "draw"
    if (Si->IsWaitShape)    printf("        colR colG colB \n");
    //if draw call, use SC colour, or shape colour, or param colour.
//    if (Si->RandRangeIsUsed) {
  //  }
    /*else if (SC->IsDrawCall) {*/
/*            if (SC->ColIsSet) printf("%g %g %g\n",SC->col[0],SC->col[1],SC->col[2]);*/
/*            else if (Si->ColIsRand) printf("randR randG randB\n");*/
/*            else if (Si->ColIsSet) printf("%g %g %g\n",Si->col[0],Si->col[1],Si->col[2]);*/
/*            else printf(" colR colG colB\n");*/
/*    } */
/*    */
    else if (SC->ColourAddIsUsed) printf(" colR %g add colG %g add colB %g add \n",SC->col[0],SC->col[1],SC->col[2]);
    else if (SC->ColIsRand) printf("rand 2147483647.0 div rand 2147483647.0 div rand 2147483647.0 div \n");
    else if (SC->ColIsSet) printf("%g %g %g \n",SC->col[0],SC->col[1],SC->col[2]);
    else if (SC->NotIsUsed) printf(" 1 colR sub 1 colG sub 1 colB sub\n");
    else if (SC->RotIsUsed) printf(" colB colR colG\n");
    else if (SC->Rot2IsUsed) printf("colG colB colR\n");
    else //if (SC->EqualIsUsed) {//if there is a shape colour, pass that. else pass param colour
        if (Si->ColIsRand || Si->RandRangeIsUsed)
        printf("randR randG randB\n");
       else if (Si->ColIsSet)
            printf("%g %g %g\n",Si->col[0],Si->col[1],Si->col[2]);
        else
            printf(" colR colG colB\n");
    //}
    //else printf(" colR colG colB\n");//printf("-101 0 0 \n");
    printf(" ");
    if (SC->IsDrawCall) //calling a draw shape, so pass array of pts
        printf("[ ");
    for (k=0;k<SC->pts;k++)
        printf("P%d ",SC->pt[k]);
    if (SC->IsDrawCall)
        printf(" ]\n");
    else //dont pass level if calling draw shape
        printf(" level\n");
}//end PrintShapeCall()

void PrintColouringStuff(){
    printf("        closepath\n");
    //So.. WHEN is gs and gr required?
    //- gradient colour seems to need it...
    // when stroking AND filling. either one alone doesnt need it
    //so.. dont print it if (no edges or only edges or !gradient)
    if (GradientIsUsed || !(OnlyEdgesIsSet || NoedgesIsSet))
        printf("        gs\n");
    if (Si->GradientIsUsed) {
        printf("		clip\n" "\t\t<</ShadingType 2 \n\t\t/ColorSpace [ /DeviceRGB ]\n");
        printf("\t\t/Coords [P%d x P%d y P%d x P%d y] \n",Si->cInfo[0].ptnum,Si->cInfo[0].ptnum,Si->cInfo[1].ptnum,Si->cInfo[1].ptnum);
        printf("\t\t/Extend [ true true ]\n\t\t/Function <<\n" 
               "\t\t/FunctionType 2\n\t\t/Domain [ 0 1 ]\n"); 
        printf("\t\t/C0 [ %g %g %g ]\n",Si->cInfo[0].col[0],Si->cInfo[0].col[1],Si->cInfo[0].col[2]); 
        printf("\t\t/C1 [ %g %g %g ]\n",Si->cInfo[1].col[0],Si->cInfo[1].col[1],Si->cInfo[1].col[2]); 
        printf("\t\t/N 1\n\t\t>>\n\t\t>> \n\t\tshfill\n"); 
    } else {
        if (Si->ShadingIsUsed) {
            printf("        /avex P0 x ");
            for (j=1;j<Si->N;j++)
                printf("P%d x add ",j);
            printf("%d div def\n",Si->N);
            printf("        /avey P0 y ");
            for (j=1;j<Si->N;j++)
                printf("P%d y add ",j);
            printf("%d div def\n",Si->N);
            printf("        %d [",Si->colPts); 
            for (j=0;j<Si->colPts;j++) 
                printf("[%g %g] ",initPt[Si->cInfo[j].ptnum].x,initPt[Si->cInfo[j].ptnum].y);
            printf("] [");
            for (j=0;j<Si->colPts;j++) 
                printf("[%g %g %g ] ",Si->cInfo[j].col[0],Si->cInfo[j].col[1],Si->cInfo[j].col[2]);
            printf("] [avex avey] \n" "        %%top of stack: N [ptlocs] [ptcolours] [x y]\n"
                   "        getshadingcolour setrgbcolor\n");
        }
        else 
            if (ColourAddIsUsed)
                if (i!=0 || maxLevel>0) //level is only ever 0 for shape 0
                    printf("        colR colG colB setrgbcolor\n");
                else { //shape is 0 so test if level==0..IF MAXLEVEL==0
                    printf("        level 0 eq\n" "        {colour0}\n"
                           "        {colR colG colB setrgbcolor} ifelse\n");
                }
            else 
                if (!OnlyEdgesIsSet) { //otherwise print this, but skip if "only edges"
                    printf("        colR -100 gt\n" "        {colR colG colB setrgbcolor}\n");
                    if (Si->ColIsSet || Si->ColIsRand)
                        printf("        {colour%d} ifelse\n",i);
                    else
                        printf("        if\n");
                }
        if (!OnlyEdgesIsSet)
            printf("        fill\n");
    }//!gradient
    if (GradientIsUsed || !(OnlyEdgesIsSet || NoedgesIsSet))
        printf("        gr\n");
}

void PrintShapeProcedure(short i){
    Si=&Shape[i]; // <<==== NB! Si is the ith Shape
    printf("\n/subdiv-%s { %% ( colR colG colB ",Si->name);
    for (j=0;j<Si->N;j++)
        printf("P%d ",j);
    printf("level -- )\n" "    /level exch 1 add def\n");
    for (j=Si->N-1;j>=0;j--) 
        printf("    /P%d exch def\n",j);
    printf("    /colB exch def\n    /colG exch def\n    /colR exch def\n");
    if (Si->RandRangeIsUsed) {
        printf("/randR rand 2147483647.0 div %g %g sub mul %g add def\n",Si->randRange[0].max,Si->randRange[0].min,Si->randRange[0].min);
        printf("/randG rand 2147483647.0 div %g %g sub mul %g add def\n",Si->randRange[1].max,Si->randRange[1].min,Si->randRange[1].min);
        printf("/randB rand 2147483647.0 div %g %g sub mul %g add def\n",Si->randRange[2].max,Si->randRange[2].min,Si->randRange[2].min);
//        printf("/randG rand 2147483647.0 div def\n");
  //      printf("/randB rand 2147483647.0 div def\n");
    } else if (Si->ColIsRand) {//need colour to pass, different each time through function.
        printf("/randR rand 2147483647.0 div def\n");
        printf("/randG rand 2147483647.0 div def\n");
        printf("/randB rand 2147483647.0 div def\n");
    }
    if (Si->DrawNow) 
        PrintNewShapePoints();
    else
       printf("    level //MAXLEVEL ge\n    {\n");
    if (Si->N>2 || !NoedgesIsSet) {//i.e. if a line and NO EDGES, dont bother making path.
        printf("        %%draw %s\n",Si->name);
        printf("        newpath\n" "        P0 xy moveto\n");
        for (j=1;j<Si->N;j++) 
            printf("        P%d xy lineto\n",j);
    }
    if (Si->N>2) //lines (N=2) dont need any of this colouring stuff, just stroke.
        PrintColouringStuff();
    if (!NoedgesIsSet)
        printf("        stroke\n");
    if (Si->DrawNow)
        printf("} bind def\n");
    else {
        printf("    }{\n"
               "        %%level<MAXLEVEL, so do %s subdivision\n",Si->name);
        PrintNewShapePoints();
        if (Si->BranchingIsUsed>0) 
            PrintBranchedShapeCalls();
        else {//no branches
            for (j=0;j<Si->ShapeCalls;j++) 
                PrintShapeCall(j);
            for (j=Si->ShapeCalls-1;j>=0;j--) 
                if (Si->ShapeCall[j].IsDrawCall)
                    printf("        draw\n");
                else
                    printf("        subdiv-%s\n",Si->ShapeCall[j].name);
        }
        printf("    } ifelse\n" "} bind def\n");
    }
}

void PrintBranchedShapeCalls() {
    printf("level ");
    if (Si->BranchingIsUsed==GREATERTHAN) 
        printf(Si->branchLevel>=0 ? "%d le {\n" : "MAXLEVEL %d add le {\n",Si->branchLevel);
    else //i.e. it uses EQUALS
        printf(Si->branchLevel>=0 ? "%d ne {\n" : "MAXLEVEL %d add ne {\n",Si->branchLevel);
    for (j=0;j<Si->firstBranchedSCall;j++) 
        PrintShapeCall(j);
    for (j=Si->firstBranchedSCall-1;j>=0;j--) 
        if (Si->ShapeCall[j].IsDrawCall)
            printf("        draw\n");
        else
            printf("        subdiv-%s\n",Si->ShapeCall[j].name);
    printf("}\n{\n");
    for (j=Si->firstBranchedSCall;j<Si->ShapeCalls;j++) 
        PrintShapeCall(j);
    for (j=Si->ShapeCalls-1;j>=Si->firstBranchedSCall;j--) 
        if (Si->ShapeCall[j].IsDrawCall)
            printf("        draw\n");
        else
            printf("        subdiv-%s\n",Si->ShapeCall[j].name);
    printf("} ifelse\n");
}

void MakeLogFiles(){
    fprintf(stderr,"WAIT LIST\n");
    for (i=0;i<=s;i++) 
        fprintf(stderr,"Shape %d max wait = %d\n",i,MaxWaitLevel[i]);
//FOR TESTING ONLY::
    fprintf(stderr,"NoedgesIsSet: %d\n",NoedgesIsSet);
    fprintf(stderr,"BoxColourIsSet: %d\n",BoxColourIsSet);
    fprintf(stderr,"RandIsUsed: %d\n",RandIsUsed);
    fprintf(stderr,"OnlyEdgesIsSet: %d\n",OnlyEdgesIsSet);
    fprintf(stderr,"GradientIsUsed: %d\n",GradientIsUsed);
    fprintf(stderr,"CentreIsUsed: %d\n",CentreIsUsed);
    fprintf(stderr,"ShadingIsUsed: %d\n",ShadingIsUsed);
    fprintf(stderr,"ColourAddIsUsed: %d\n",ColourAddIsUsed);    
    fprintf(stderr,"DrawIsUsed: %d\n",DrawIsUsed);
    fprintf(stderr,"MarginArrayIsUsed: %d\n\n",MarginArrayIsUsed);
    for (i=0;i<=s;i++) {
        Si=&Shape[i];
        fprintf(stderr,"======================================\nShape %d name: %s N: %d\n",i,Si->name,Si->N);
        fprintf(stderr,"ColourIsSet: %d ",Si->ColIsSet);
        fprintf(stderr,"ColourIsRand: %d\n",Si->ColIsRand);
        fprintf(stderr,"ShadingIsUsed: %d\n",Si->ShadingIsUsed);
        fprintf(stderr,"GradientIsUsed: %d\n",Si->GradientIsUsed);
            if (Si->GradientIsUsed){
                fprintf(stderr,"    Gradient info\n");
fprintf(stderr,"    Gradient point #0 is p%d\n",Si->cInfo[0].ptnum);
fprintf(stderr,"        Colour is %g %g %g\n",Si->cInfo[0].col[0],Si->cInfo[0].col[1],Si->cInfo[0].col[2]);
fprintf(stderr,"    Gradient point #1 is p%d\n",Si->cInfo[1].ptnum);
fprintf(stderr,"        Colour is %g %g %g\n",Si->cInfo[1].col[0],Si->cInfo[1].col[1],Si->cInfo[1].col[2]);
            }
       fprintf(stderr,"colourPts: %d ",Si->colPts);
        fprintf(stderr,"colour: %g %g %g\n",Si->col[0],Si->col[1],Si->col[2]);
        fprintf(stderr,"NewPts: %d\n",Si->NewPts);
        for (j=Si->N;j<Si->NewPts+Si->N;j++){
fprintf(stderr,"New Pt #%d : ",j);
fprintf(stderr,"\tFrom p%d To p%d\n",Si->NewPt[j].From,Si->NewPt[j].To);
fprintf(stderr,"\tnumer %d denom %d ",Si->NewPt[j].numer,Si->NewPt[j].denom);
fprintf(stderr,"\tdecFract %g\n",Si->NewPt[j].decFract);
fprintf(stderr,"\tturnDirection %d turnDist %g\n",Si->NewPt[j].turnDirection,Si->NewPt[j].turnDist);
fprintf(stderr,"\tIsDecFract %d IsRand %d IsCentre %d\n",Si->NewPt[j].IsDecFract,Si->NewPt[j].IsRand,Si->NewPt[j].IsCentre);
        }
        fprintf(stderr,"DrawNow: %d\n",Si->DrawNow);
        fprintf(stderr,"ShapeCalls: %d\n",Si->ShapeCalls);
        for (j=0;j<Si->ShapeCalls;j++){
            fprintf(stderr,"ShapeCall %d calls --> %s N: %d\n",j,Si->ShapeCall[j].name,Si->ShapeCall[j].pts);
            fprintf(stderr,"    ColourIsSet: %d\n",Si->ShapeCall[j].ColIsSet);
            fprintf(stderr,"    ColourIsRand: %d\n",Si->ShapeCall[j].ColIsRand);
            fprintf(stderr,"    IsDrawCall: %d\n",Si->ShapeCall[j].IsDrawCall);
            
fprintf(stderr,"    colour: %g %g %g\n",Si->ShapeCall[j].col[0],Si->ShapeCall[j].col[1],Si->ShapeCall[j].col[2]);
fprintf(stderr,"    ColourAddIsUsed: %d\n",Si->ShapeCall[j].ColourAddIsUsed);
fprintf(stderr,"    waitLevel: %d waitShapeName: %s\n",Si->ShapeCall[j].waitLevel,Si->ShapeCall[j].waitShapeName);
         }
    }
}

void BuildNewWaitShapes(){
    short newS;
    for (i=0;i<=s;i++) 
        for (j=0;j<Shape[i].ShapeCalls;j++)
            if (Shape[i].ShapeCall[j].waitLevel>MaxWaitLevel[i])
                MaxWaitLevel[GetShapeNum(Shape[i].ShapeCall[j].waitShapeName)]=Shape[i].ShapeCall[j].waitLevel;
    newS=s; //will be higher after new shapes added, but updating s would mess with i loop.
    for (i=0;i<=s;i++) {
        if (MaxWaitLevel[i]>0) {
            for (k=MaxWaitLevel[i];k>=1;k--){//add e.g. tri-wait-2,tri-wait-1
                newS++;
                sprintf(shapename,"%s-wait-%d",Shape[i].name,k);
                strcpy(Shape[newS].name,shapename);
                Shape[newS].N=Shape[i].N; //i.e. tri-wait-2 has N=3, same as tri
                Shape[newS].ShapeCalls=1;
                if (k>1) {
                    sprintf(shapename,"%s-wait-%d",Shape[i].name,k-1);
                    strcpy(Shape[newS].ShapeCall[0].name,shapename);
                } else //replace for -wait-1 the name of shape call with wait
                    strcpy(Shape[newS].ShapeCall[0].name,Shape[i].name);
                Shape[newS].ShapeCall[0].pts=Shape[newS].N;
                Shape[newS].IsWaitShape=1;
                for (sp=0;sp<Shape[newS].N;sp++)
                    Shape[newS].ShapeCall[0].pt[sp]=sp;
            }
        }
    }  
    s=newS; //the new highest shape number i.e. Shape[s]. There's no Shape[s+1].
}

void SetShapeCol(float r,float g,float b){
    Ss->col[0]=r;
    Ss->col[1]=g;
    Ss->col[2]=b; 
    Ss->ColIsSet=1;
}

void SetShadingCol(int cpn, float r,float g,float b){
    Ss->ShadingIsUsed=1;
    ShadingIsUsed=1;
    Ss->cInfo[Ss->colPts].ptnum=cpn;
    Ss->cInfo[Ss->colPts].col[0]=r;
    Ss->cInfo[Ss->colPts].col[1]=g;
    Ss->cInfo[Ss->colPts].col[2]=b;
    Ss->colPts++;
}

void SetShapeCallCol(float r,float g,float b){
    Ss->ShapeCall[nss].col[0]=r;
    Ss->ShapeCall[nss].col[1]=g;
    Ss->ShapeCall[nss].col[2]=b;
}

short GetShapeNum(char name[50]){
    short i;
    for (i=0;i<=s;i++)
        if (strcmp(name,Shape[i].name)==0)
            return i;
    yyerror("Shape not found!");
}

void IncPointVars(){
    Ss->NewPts++; //this is the first, second etc new point.
    np++;
    P=&Ss->NewPt[np];
}

void IncShapePointVars(){
    Ss->ShapeCall[nss].pts++;
    nssp++;
}

short ColNumber(char* cnamestr){
    short i;
    for (i=0;i<numCols;i++) 
        if (strcmp(cnamestr,col[i].name)==0) 
            return i; 
    yyerror("Colour name not found");
}

void PrintTurnX(){
    printf("\t\tP%d y P%d y sub ",P->To,P->From); 
    printf("%g mul ",P->turnDist);
    if (P->turnDirection==LEFT) 
        printf("sub\n");
    else  //RIGHT
        printf("add\n");
}

void PrintTurnY(){
    printf("\t\tP%d x P%d x sub ",P->To,P->From); 
    printf("%g mul ",P->turnDist);
    if (P->turnDirection==LEFT) 
        printf("add\n");
    else  //RIGHT
        printf("sub\n");
}

void SetFromAndToPts(short fromP,short toP){
    P->From=fromP;
    P->To=toP; 
}

void CheckPoints(short cp1,short cp2,short cp3){
    if (np!=cp1) {
        sprintf(errstr,"Bad definition of new point p%d should be p%d instead.",cp1,np);
        yyerror(errstr);
    }
    if (cp1==cp2 || cp1==cp3 || cp2==cp3) 
        yyerror("Bad new point definition. Duplicated point numbers.");
    if (cp3>np){
        sprintf(errstr,"Bad new point definition. There is no p%d.",cp3);
        yyerror(errstr);
    }
    if (cp2>np) {
        sprintf(errstr,"Bad new point definition. There is no p%d.",cp2);
        yyerror(errstr);
    }
}

void CheckAndSetPoints(short cp1, short cp2, short cp3) {
    CheckPoints(cp1,cp2,cp3);
    SetFromAndToPts(cp2,cp3);
}

void CheckPtNumber(){
    if (Ss->ShapeCall[nss].pt[nssp]>np-1){
        sprintf(errstr,"Bad point number - p%d not defined in %s.",Ss->ShapeCall[nss].pt[nssp],Ss->name);
        yyerror(errstr);
    }
}

void SetInitPt(float x, float y){
    if (GridIsSet) {
        initPt[initPts].x=x+y*cos(2*PI*gridTheta/360);
        initPt[initPts].y=y*sin(2*PI*gridTheta/360);
    }
    else {
        initPt[initPts].x=x;
        initPt[initPts].y=y;
    }
    initPts++;
}

void InitPtError(short ip)    {
        sprintf(errstr,"Incorrect initial point number p%d, should be p%d.",ip,initPts);
        yyerror(errstr);
}

void PrintCentreProcedure() {
    printf("\n/find-centre {1 dict begin %% [ P0 P1 .. Pn-1 ] -- x y\n");
    printf("        /pts exch def\n" "        /xsum 0 def\n");
    printf("        /ysum 0 def\n" "        /len pts length def\n");
    printf("        0 1 len 1 sub {\n");
    printf("                /i exch def\n");
    printf("                /xsum pts i get x xsum add def\n");
    printf("                /ysum pts i get y ysum add def\n");
    printf("        } for\n" "        xsum len div ysum len div\n");
    printf("end } bind def\n");
}

void PrintShadingColProcedure() {
    printf("\n/getshadingcolour {1 dict begin " 
        "%% ( numpts [Plocs] [Pcolours] [x y] -- r g b )\n");
    printf("    /P exch def\n" "    /Pcolours exch def\n");
    printf("    /Plocs exch def\n"  "    /N exch def\n");
    printf("    /ith {i get} def\n" "    /jth {j get} def\n");
    printf("    /distsum 0 def\n"   "    %%make array with dists to the points\n");
    printf("    /distAry [\n" "        0 1 N 1 sub {\n");
    printf("            /i exch def\n" "            /whichpt Plocs ith def\n");
    printf("            /px {whichpt x} def\n");
    printf("            /py {whichpt y} def\n");
    printf("            /distsq  P x px sub 2 exp  P y py sub 2 exp  add def\n");
    printf("            distsq 0 eq\n");
    printf("            {/invdistsq 10000 def } %%avoids divide-by-0\n");
    printf("            {/invdistsq 1 distsq div def } ifelse %%invert\n");
    printf("            /distsum distsum invdistsq add def\n");
    printf("            invdistsq  %%put in array\n");
    printf("        } for\n" "    ] def\n");
    printf("    /weightedAry [\n" "        0 1 N 1 sub {\n");
    printf("            /i exch def\n" "            distAry ith distsum div\n");
    printf("%%divide all by total to get numbers adding to 1.\n");
    printf("        } for\n" "    ] def\n");
    printf("    %%calculate colour\n" "    0 1 2 { %%loop over rgb: 0,1,2\n");
    printf("        /i exch def\n" "        /col 0 def\n");
    printf("        0 1 N 1 sub {\n" "            /j exch def\n");
    printf("            /product Pcolours jth ith weightedAry jth mul def\n");
    printf("            /col col product add def\n" "        } for\n");
    printf("        col %%put 1 colour number on stack\n");
    printf("    } for\n" "end } def\n\n");    
  }

void SetPboxPoints(float x1, float y1, float x2, float y2) {
    SetInitPt(x1,y1);
    SetInitPt(x1,y2);
    SetInitPt(x2,y2);
    SetInitPt(x2,y1);
}

void SetGradientPtInfo(int gp, int ptnum, float c0, float c1, float c2) {
    Ss->cInfo[gp].ptnum=ptnum;
    Ss->cInfo[gp].col[0]=c0;
    Ss->cInfo[gp].col[1]=c1;
    Ss->cInfo[gp].col[2]=c2;
    GradientIsUsed=YES;
    Ss->GradientIsUsed=YES;
}

void SetDecFract(int pfrom, int pto,float dfract){
    SetFromAndToPts(pfrom,pto);
    P->decFract=dfract; 
    P->IsDecFract=1; 
    IncPointVars();
}

void SetCol(char* name, int r, int g, int b){
    col[numCols].name=name;
    col[numCols].c[0]=r;
    col[numCols].c[1]=g;
    col[numCols].c[2]=b;
    numCols++;
}

void SetColNames() { // 146 colours <-- update this number as colours are added
    SetCol("alice blue",240,248,255);
    SetCol("antique white",250,235,215);
    SetCol("aqua",0,255,255);
    SetCol("aquamarine",127,255,212);
    SetCol("azure",240,255,255);
    SetCol("beige",245,245,220);
    SetCol("bisque",255,228,196);
    SetCol("black",0,0,0);
    SetCol("blanched almond",255,235,205);
    SetCol("blue violet",138,43,226);
    SetCol("blue",0,0,255);
    SetCol("brown",165,42,42);
    SetCol("burly wood",222,184,135);
    SetCol("cadet blue",95,158,160);
    SetCol("chartreuse",127,255,0);
    SetCol("chocolate",210,105,30);
    SetCol("coral",255,127,80);
    SetCol("corn silk",255,248,220);
    SetCol("cornflower blue",100,149,237);
    SetCol("crimson",220,20,60);
    SetCol("cyan",0,255,255);
    SetCol("dark blue",0,0,139);
    SetCol("dark cyan",0,139,139);
    SetCol("dark golden rod",184,134,11);
    SetCol("dark gray",169,169,169);
    SetCol("dark green",0,100,0);
    SetCol("dark grey",169,169,169);
    SetCol("dark khaki",189,183,107);
    SetCol("dark magenta",139,0,139);
    SetCol("dark olive green",85,107,47);
    SetCol("dark orange",255,140,0);
    SetCol("dark orchid",153,50,204);
    SetCol("dark red",139,0,0);
    SetCol("dark salmon",233,150,122);
    SetCol("dark sea green",143,188,143);
    SetCol("dark slate blue",72,61,139);
    SetCol("dark slate gray",47,79,79);
    SetCol("dark slate grey",47,79,79);
    SetCol("dark turquoise",0,206,209);
    SetCol("dark violet",148,0,211);
    SetCol("deep pink",255,20,147);
    SetCol("deep sky blue",0,191,255);
    SetCol("dim gray",105,105,105);
    SetCol("dim grey",105,105,105);
    SetCol("dodger blue",30,144,255);
    SetCol("firebrick",178,34,34);
    SetCol("floral white",255,250,240);
    SetCol("forest green",34,139,34);
    SetCol("gainsboro",220,220,220);
    SetCol("ghost white",248,248,255);
    SetCol("gold",255,215,0);
    SetCol("golden rod",218,165,32);
    SetCol("gray",128,128,128);
    SetCol("green yellow",173,255,47);
    SetCol("green",0,128,0);
    SetCol("grey",128,128,128);
    SetCol("honeydew",240,255,240);
    SetCol("hot pink",255,105,180);
    SetCol("indian red",205,92,92);
    SetCol("indigo",75,0,130);
    SetCol("ivory",255,255,240);
    SetCol("khaki",240,230,140);
    SetCol("lavender blush",255,240,245);
    SetCol("lavender",230,230,250);
    SetCol("lawn green",124,252,0);
    SetCol("lemon chiffon",255,250,205);
    SetCol("light blue",173,216,230);
    SetCol("light coral",240,128,128);
    SetCol("light cyan",224,255,255);
    SetCol("light golden rod yellow",250,250,210);
    SetCol("light gray",211,211,211);
    SetCol("light green",144,238,144);
    SetCol("light grey",211,211,211);
    SetCol("light pink",255,182,193);
    SetCol("light salmon",255,160,122);
    SetCol("light sea green",32,178,170);
    SetCol("light sky blue",135,206,250);
    SetCol("light slate gray",119,136,153);
    SetCol("light slate grey",119,136,153);
    SetCol("light steel blue",176,196,222);
    SetCol("light yellow",255,255,224);
    SetCol("lime green",50,205,50);
    SetCol("lime",0,255,0);
    SetCol("linen",250,240,230);
    SetCol("magenta",255,0,255);
    SetCol("maroon",128,0,0);
    SetCol("medium aquamarine",102,205,170);
    SetCol("medium blue",0,0,205);
    SetCol("medium orchid",186,85,211);
    SetCol("medium purple",147,112,219);
    SetCol("medium sea green",60,179,113);
    SetCol("medium slate blue",123,104,238);
    SetCol("medium spring green",0,250,154);
    SetCol("medium turquoise",72,209,204);
    SetCol("medium violet red",199,21,133);
    SetCol("midnight blue",25,25,112);
    SetCol("mint cream",245,255,250);
    SetCol("misty rose",255,228,225);
    SetCol("moccasin",255,228,181);
    SetCol("navajo white",255,222,173);
    SetCol("navy",0,0,128);
    SetCol("old lace",253,245,230);
    SetCol("olive drab",107,142,35);
    SetCol("olive",128,128,0);
    SetCol("orange red",255,69,0);
    SetCol("orange",255,165,0);
    SetCol("orchid",218,112,214);
    SetCol("pale golden rod",238,232,170);
    SetCol("pale green",152,251,152);
    SetCol("pale turquoise",175,238,238);
    SetCol("pale violet red",219,112,147);
    SetCol("papaya whip",255,239,213);
    SetCol("peach puff",255,218,185);
    SetCol("peru",205,133,63);
    SetCol("pink",255,192,203);
    SetCol("plum",221,160,221);
    SetCol("powder blue",176,224,230);
    SetCol("purple",128,0,128);
    SetCol("red",255,0,0);
    SetCol("rosy brown",188,143,143);
    SetCol("royal blue",65,105,225);
    SetCol("saddle brown",139,69,19);
    SetCol("salmon",250,128,114);
    SetCol("sandy brown",244,164,96);
    SetCol("sea green",46,139,87);
    SetCol("sea shell",255,245,238);
    SetCol("sienna",160,82,45);
    SetCol("silver",192,192,192);
    SetCol("sky blue",135,206,235);
    SetCol("slate blue",106,90,205);
    SetCol("slate gray",112,128,144);
    SetCol("slate grey",112,128,144);
    SetCol("snow",255,250,250);
    SetCol("spring green",0,255,127);
    SetCol("steel blue",70,130,180);
    SetCol("tan",210,180,140);
    SetCol("teal",0,128,128);
    SetCol("thistle",216,191,216);
    SetCol("tomato",255,99,71);
    SetCol("turquoise",64,224,208);
    SetCol("violet",238,130,238);
    SetCol("wheat",245,222,179);
    SetCol("white smoke",245,245,245);
    SetCol("white",255,255,255);
    SetCol("yellow green",154,205,50);
    SetCol("yellow",255,255,0);
}