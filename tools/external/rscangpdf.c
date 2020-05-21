/* performs r-scan statistics to find clusters and gaps in patloc output (.pll)
Copyright 2006 Jan Mrazek (mrazek@uga.edu)  */


/************************************************************************/
/* THIS PROGRAM IS DISTRIBUTED "AS IS" WITH NO WARRANTY,                */
/* OR SUPPORT, WITHOUT EVEN THE IMPLIED WARRANTY                        */
/* OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.              */
/* THE AUTHOR SHALL NOT BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       */
/* LIABILITY IN CONNECTION WITH THIS PROGRAM.                           */
/*                                                                      */
/* YOU CAN MODIFY AND REDISTRIBUTE THIS PROGRAM UNDER THE TERMS OF      */
/* THE GNU GENERAL PUBLIC LICENSE (http://www.gnu.org/licenses/gpl.txt) */
/*                                                                      */
/************************************************************************/

/********************************************************************************/
/* INSTALLATION INSTRUCTIONS:                                                   */
/* Compile with a command "gcc rscangpdf.c -lm -O3 -o ~/bin/rscangpdf" or       */
/* equivalent. cc compiler should work, too.                                    */
/* You may leave out the -O3 option (optimization) if it's causing problems     */
/* The program was tested on Redhat Linux.                                      */
/* Run as "rscangpdf <InputFile> <OutputFile>", where <InputFile> is the .pll   */
/* output from patloc, <OutputFile> is the output file name (extensions will be */
/* added by the program)                                                        */
/********************************************************************************/

/*
Known problems:
The program is stack memory-hungry. If you get a "segmentation fault" or other nondescript
runtime error, try some of the follwing:
1. increase the stack limit (ulimit -s unlimited)
2. lower the MAXMARKERS and MAXCLUSTERS values (a few lines below) and recompile the program
3. rewrite the program with dynamic array declarations
*/





#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXLINE 513
#define MAXMARKERS 700000
#define MAXCLUSTERS 20000
#define MINNTORRATIO 5.0



/********************************************************************/
/*                      Global declarations                         */
/********************************************************************/


long LSeq;





/********************************************************************/
/*                                                                  */
/*                      INCLUDED PROCEDURES                         */
/*                                                                  */
/********************************************************************/




void MakeLine (Start,End,Shift,Scale,OFi)
long Start;
long End;
long Shift;
double Scale;
FILE *OFi;
{
  long I,J,N;
  if (Start<=End)    {
    N=floor((double)(Start+1)/500000.0);
    I=Start-N*500000;
    J=End-N*500000;
    while (J>500000)    {
      fprintf (OFi,"%f %li M\n",I*Scale,Shift-N*2500);
      fprintf (OFi,"%f %li L S\n",500000*Scale,Shift-N*2500);
      ++N;
      I=0;
      J=End-N*500000;
    };
    fprintf (OFi,"%f %li M\n",I*Scale,Shift-N*2500);
    fprintf (OFi,"%f %li L S\n",J*Scale,Shift-N*2500);    }
  else    {
    MakeLine(1,End,Shift,Scale,OFi);
    MakeLine(Start,LSeq,Shift,Scale,OFi);
  };
}


double RSProbMax (X,N)
double X;
double N;
{
  double S=0.0;
  double P=1.0;
  double PP;
  int I,K;

  K=1;
  I=1;
  while (I<=N && I*X<1)    {
    P=P*(N+1-I)/I;
    K=-K;
    S-=P*K*pow(1.0-I*X,N-1);
    ++I;
  };
  return(S);
}



/***********************************************************************/
/***********************************************************************/
/***********************************************************************/




int main (argc, argv)
int argc;
char *argv[];
{
  FILE *IFi;
  FILE *OFi;
  long NMarkers;
  int I;
  long R,KK,LI,LP,L,LPP,MinDistMin1,MinDistMin5,MaxDistMax1,MaxDistMax5,J;
  long MinDistMax1[30],MinDistMax5[30],MaxDistMin1[30],MaxDistMin5[30],MinDist[30],MaxDist[30];
  static long MkrStart[MAXMARKERS],MkrEnd[MAXMARKERS];
  double F,P,X,PMax1,PMax5,PMin1,PMin5,A,B,C,PC;
  char Line[MAXLINE],OFiName[MAXLINE];
  long NClusters1,Start1[MAXCLUSTERS],End1[MAXCLUSTERS];
  long NClusters5,Start5[MAXCLUSTERS],End5[MAXCLUSTERS];
  long NGaps1,StartGaps1[MAXCLUSTERS],EndGaps1[MAXCLUSTERS];
  long NGaps5,StartGaps5[MAXCLUSTERS],EndGaps5[MAXCLUSTERS];
  char IQNone;

  if (argc!=4)    {
    printf ("\n\n\nUsage:  rscang InputFile OutputFile\n\n");
    exit(0);
  };


  NMarkers=-1;
  if ( (IFi=fopen(argv[1],"r")) == NULL )  {
    printf ("Unable to open input file %s.\n",argv[1]);
    exit (0);
  };
  strcpy(OFiName,argv[2]);
  //strcat(OFiName,".txt");
  if ( (OFi=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file %s.\n",OFiName);
    exit (0);
  };
  fprintf (OFi,"r-scan statistics on motifs found by Pattern Locator or Motif Locator.\n",argv[1]);
  if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
    printf ("Unable to read from file %s.\n",argv[1]);
    exit (0);
  };
  if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
    printf ("Unable to read from file %s.\n",argv[1]);
    exit (0);
  };
  I=0;
  while (Line[I]!='(')    ++I;
  LSeq=atol(Line+I+1);
  fprintf (OFi,"%s",Line);
  while (strncmp(Line,"     Start       End",19)!=0)    {
    if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
      printf ("Unable to read from file %s.\n",argv[1]);
      exit (0);
    };
    if (strncmp(Line,"     Start       End",19)!=0)    fprintf (OFi,"%s",Line);
  };
  while (feof(IFi)==0)    {
    if ( fgets(Line,MAXLINE-1,IFi)!=NULL )    {
      ++NMarkers;
      if (NMarkers>=MAXMARKERS)   {
        printf ("Motif(s) too general\n");
        exit (2);
      };
      MkrStart[NMarkers]=atol(Line);
      MkrEnd[NMarkers]=atol(Line+10);
      if (MkrStart[NMarkers]>MkrEnd[NMarkers])    {
        LI=MkrStart[NMarkers];
        MkrStart[NMarkers]=MkrEnd[NMarkers];
        MkrEnd[NMarkers]=LI;
      };
    };
  };
  fclose (IFi);


  /* sorting markers */
  for (L=0;L<=NMarkers;++L)    {
    LP=L-1;
    while (MkrStart[LP]>MkrStart[LP+1] && LP>=0)    {
      KK=MkrStart[LP];
      MkrStart[LP]=MkrStart[LP+1];
      MkrStart[LP+1]=KK;
      KK=MkrEnd[LP];
      MkrEnd[LP]=MkrEnd[LP+1];
      MkrEnd[LP+1]=KK;
      --LP;
    };
  };

  fprintf (OFi,"Number of patterns: %li\n",NMarkers+1);


  /* positions of markers are now known. NMarkers+1 markers were found and their starting and ending
     positions in the analyzed sequence are now in MkrStart[] and MkrEnd[] */

  /* combining overlaping markers */
  for (LP=1;LP<=NMarkers;++LP)    {
    while (MkrStart[LP]<=MkrEnd[LP-1] && LP<=NMarkers)    {
      if (MkrEnd[LP]>MkrEnd[LP-1])    MkrEnd[LP-1]=MkrEnd[LP];
      --NMarkers;
      for (LPP=LP;LPP<=NMarkers;++LPP)    {
        MkrStart[LPP]=MkrStart[LPP+1];
        MkrEnd[LPP]=MkrEnd[LPP+1];
      };
    };
  };
  fprintf (OFi,"Number of nonoverlapping patterns: %li\n\n",NMarkers+1);
  fprintf (OFi,"Locations of statistically significant clusters and overdispersions (gaps):\n");

  NClusters1=-1;
  NClusters5=-1;
  NGaps1=-1;
  NGaps5=-1;
  for (R=1;R<=20 && R<=NMarkers-1;++R)    {
    P=1.0;
    for (J=1;J<=R-1;++J)    {
      P=P*J;
    };
    L=LSeq;
    for (LP=0;LP<=NMarkers;++LP)    {
      L-=MkrEnd[LP]+1-MkrStart[LP];
    };
    F=NMarkers;
    X=log(F);
    X=X+(R-1)*log(X);
    PMax5=(X-log(-P*log(0.95)))/F;
    PMax1=(X-log(-P*log(0.99)))/F;
    PMin5=(X-log(-P*log(0.05)))/F;
    PMin1=(X-log(-P*log(0.01)))/F;
    /* exact formula */
    if (R==1 && NMarkers<=1000)    {
      B=1.0;
      A=1.0/F;
      while (fabs(A-B)>0.000001)    {
        C=(A+B)/2;
        PC=RSProbMax(C,F);
        if (PC>0.01)    {
          A=C;    }
        else    {
          B=C;
        };
      };
      PMax1=(A+B)/2;
      B=1.0;
      A=1.0/F;
      while (fabs(A-B)>0.000005)    {
        C=(A+B)/2;
        PC=RSProbMax(C,F);
        if (PC>0.05)    {
          A=C;    }
        else    {
          B=C;
        };
      };
      PMax5=(A+B)/2;
      B=1.0;
      A=1.0/F;
      while (fabs(A-B)>0.000005)    {
        C=(A+B)/2;
        PC=RSProbMax(C,F);
        if (PC>0.95)    {
          A=C;    }
        else    {
          B=C;
        };
      };
      PMin5=(A+B)/2;
      B=1.0;
      A=1.0/F;
      while (fabs(A-B)>0.000001)    {
        C=(A+B)/2;
        PC=RSProbMax(C,F);
        if (PC>0.99)    {
          A=C;    }
        else    {
          B=C;
        };
      };
      PMin1=(A+B)/2;
    };
    MaxDistMax5=floor(0.5+L*PMax5);
    MaxDistMax1=floor(0.5+L*PMax1);
    MaxDistMin5[R]=floor(0.5+L*PMin5);
    MaxDistMin1[R]=floor(0.5+L*PMin1);


    if (R==1)    {
      MinDistMin5=floor(0.5+L/F*(1.0-pow(0.95,1.0/(F-1.0))));
      MinDistMin1=floor(0.5+L/F*(1.0-pow(0.99,1.0/(F-1.0))));
      MinDistMax5[R]=floor(0.5+L/F*(1.0-pow(0.05,1.0/(F-1.0))));
      MinDistMax1[R]=floor(0.5+L/F*(1.0-pow(0.01,1.0/(F-1.0))));
    }
    else    {
      P=P*R;
      X=pow(F,1.0+1.0/(double)(R));
      MinDistMin5=floor(0.5+L*pow(-log(0.95)*P,1.0/(double)(R))/X);
      MinDistMin1=floor(0.5+L*pow(-log(0.99)*P,1.0/(double)(R))/X);
      MinDistMax5[R]=floor(0.5+L*pow(-log(0.05)*P,1.0/(double)(R))/X);
      MinDistMax1[R]=floor(0.5+L*pow(-log(0.01)*P,1.0/(double)(R))/X);
    };

    if ((double)(R)>log(F) && R>1)    {
      MinDistMin1=0;
      MinDistMin5=0;
      MinDistMax1[R]=LSeq;
      MinDistMax5[R]=LSeq;
    };

    X=log(P*F/MINNTORRATIO);
    if ((double)(R)>X/log(X) && R>1)    {
      MaxDistMin1[R]=0;
      MaxDistMin5[R]=0;
      MaxDistMax1=LSeq;
      MaxDistMax5=LSeq;
    };


    MinDist[R]=LSeq;
    MaxDist[R]=0;
    for (LP=0;LP<=NMarkers;++LP)    {
      L=0;
      for (LPP=LP-R;LPP<=LP-1;++LPP)    {
        if (LPP>=0)    {
          L+=MkrStart[LPP+1]-MkrEnd[LPP]-1;    }
        else if (LPP==-1)    {
          L+=MkrStart[0]-MkrEnd[NMarkers]+LSeq+1;    }
        else    {
          L+=MkrStart[LPP+2+NMarkers]-MkrEnd[LPP+NMarkers+1]-1;
        };
      };
      if (L<MinDist[R])    MinDist[R]=L;
      if (L>MaxDist[R])    MaxDist[R]=L;
      if (L<=MinDistMin5)    {
        LPP=LP-R;
        if (LPP<0)    LPP+=NMarkers+1;
        ++NClusters5;
        Start5[NClusters5]=MkrStart[LPP];
        End5[NClusters5]=MkrEnd[LP];
        if (L<=MinDistMin1)    {
          ++NClusters1;
          Start1[NClusters1]=MkrStart[LPP];
          End1[NClusters1]=MkrEnd[LP];
          fprintf (OFi,"%8li - %8li : Cluster, r=%li, probability <1%c\n",Start1[NClusters1],End1[NClusters1],R,'%');
        }
        else    {
          fprintf (OFi,"%8li - %8li : Cluster, r=%li, probability <5%c\n",Start5[NClusters5],End5[NClusters5],R,'%');
        };
      };
      if (L>=MaxDistMax5 /*&& R==1*/)    {
        LPP=LP-R;
        if (LPP<0)    LPP+=NMarkers+1;
        ++NGaps5;
        StartGaps5[NGaps5]=MkrEnd[LPP];
        EndGaps5[NGaps5]=MkrStart[LP];
        if (L>=MaxDistMax1)    {
          ++NGaps1;
          StartGaps1[NGaps1]=MkrEnd[LPP];
          EndGaps1[NGaps1]=MkrStart[LP];
          fprintf (OFi,"%8li - %8li : Gap or overdispersion, r=%li, probability <1%c\n",StartGaps1[NGaps1],EndGaps1[NGaps1],R,'%');
        }
        else    {
          fprintf (OFi,"%8li - %8li : Gap or overdispersion, r=%li, probability <5%c\n",StartGaps5[NGaps5],EndGaps5[NGaps5],R,'%');
        };
      };
    };
  };

  fprintf (OFi,"\nSignificantly even distribution of patterns:\n");
  IQNone=1;
  for (R=1;R<=20 && R<=NMarkers-1;++R)    {
    if (MinDist[R]>MinDistMax1[R])    {
      IQNone=0;
      fprintf (OFi,"Minimum observed r-scan length (%li) is higher than expected (%li): r=%li, probability <1%c\n",MinDist[R],MinDistMax1[R],R,'%');
    }
    else if (MinDist[R]>MinDistMax5[R])    {
      IQNone=0;
      fprintf (OFi,"Minimum observed r-scan length (%li) is higher than expected (%li): r=%li, probability <5%c\n",MinDist[R],MinDistMax5[R],R,'%');
    };
    if (MaxDist[R]<MaxDistMin1[R])    {
      IQNone=0;
      fprintf (OFi,"Maximum observed r-scan length (%li) is lower than expected (%li): r=%li, probability <1%c\n",MaxDist[R],MaxDistMin1[R],R,'%');
    }
    else if (MaxDist[R]<MaxDistMin5[R])    {
      IQNone=0;
      fprintf (OFi,"Maximum observed r-scan length (%li) is lower than expected (%li): r=%li, probability <5%c\n",MaxDist[R],MaxDistMin5[R],R,'%');
    }
  };
  if (IQNone==1)    {
    fprintf (OFi,"None detected.\n");
  };


  fclose (OFi);


/* Output circle */
{
  FILE *OFi;
  char Proc='%';
  double Scale;
  long I,J,KK,LI,LP,LII;

  strcpy(OFiName,argv[3]);
  //strcat(OFiName,".c.ps");
  if ( (OFi=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file.\n");
    exit(0);
  };

  Scale=360.0/(double)(LSeq);
  fprintf (OFi,"%c!PS-Adobe-2.0\n",Proc);
  fprintf (OFi,"%%BoundingBox: 60 70 540 730\n");
  fprintf (OFi,"gsave\n");
  fprintf (OFi,"1 dict begin\n");
  fprintf (OFi,"/showpage {} def\n");
  fprintf (OFi,"/lwidth 100 def\n");
  fprintf (OFi,"/vshift -23 def\n");
  fprintf (OFi,"/S { stroke } def\n");
  fprintf (OFi,"/Lshow { 10 10 scale 0 vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/Rshow { 10 10 scale dup stringwidth pop neg vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/Cshow { 10 10 scale dup stringwidth pop -2 div vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/WL { S lwidth 0.35 mul setlinewidth } def\n");
  fprintf (OFi,"/NL { S lwidth setlinewidth } def\n");
  fprintf (OFi,"/BL { S lwidth 1.6 mul setlinewidth } def\n");
  fprintf (OFi,"/BBL { S lwidth 2.4 mul setlinewidth } def\n");
  fprintf (OFi,"/TL0 { S lwidth 2.5 mul setlinewidth } def\n");
  fprintf (OFi,"/TL1 { S lwidth 4 mul setlinewidth } def\n");
  fprintf (OFi,"/TL2 { S lwidth 6 mul setlinewidth } def\n");
  fprintf (OFi,"/TL3 { S lwidth 8 mul setlinewidth } def\n");
  fprintf (OFi,"/M {moveto} def\n");
  fprintf (OFi,"/T {translate} def\n");
  fprintf (OFi,"/L {lineto} def\n");
  fprintf (OFi,"/RL { rlineto } def\n");
  fprintf (OFi,"0.01 0.01 scale\n");
  fprintf (OFi,"30000 48000 T\n");
  fprintf (OFi,"0.00 setgray\n");
/*  fprintf (OFi,"/Helvetica findfont 220 scalefont setfont\n");
  fprintf (OFi,"0 1000 M\n");
  fprintf (OFi,"(%s) Cshow\n",Genome);*/
  fprintf (OFi,"/Helvetica findfont 130 scalefont setfont\n");
  fprintf (OFi,"0 -1000 M\n");
  fprintf (OFi,"(%li bp) Cshow\n",LSeq);
  fprintf (OFi,"NL\n");
  fprintf (OFi,"0 0 25000 0 360 arc S\n");
  if (LSeq<=1000000)     {
    fprintf (OFi,"TL0\n");
    for (I=0;I<=LSeq;I+=50000)    {
      fprintf (OFi,"0 0 25000 %f %f arc S\n",90-Scale*I-0.1,90-Scale*I+0.1);
    };
  };
  fprintf (OFi,"TL1\n");
  for (I=0;I<=LSeq;I+=100000)    {
    fprintf (OFi,"0 0 25000 %f %f arc S\n",90-Scale*I-0.1,90-Scale*I+0.1);
  };
  fprintf (OFi,"TL2\n");
  for (I=0;I<=LSeq;I+=500000)    {
    fprintf (OFi,"0 0 25000 %f %f arc S\n",90-Scale*I-0.15,90-Scale*I+0.15);
  };
  fprintf (OFi,"TL3\n");
  for (I=0;I<=LSeq;I+=1000000)    {
    fprintf (OFi,"0 0 25000 %f %f arc S\n",90-Scale*I-0.15,90-Scale*I+0.15);
  };

  fprintf (OFi,"0.00 setgray\n");
  fprintf (OFi,"BL\n");
  KK=0;
  while (KK<=NMarkers)    {
    fprintf (OFi,"0 0 25800 %f %f arc S\n",I,90-Scale*MkrEnd[KK],90-Scale*MkrStart[KK]);
    ++KK;
  };


  fprintf (OFi,"1 0.75 0.25 setrgbcolor\n");
  fprintf (OFi,"BL\n");
  for (I=0;I<=NGaps5;++I)    {
    fprintf (OFi,"0 0 24200 %f %f arc S\n",90-Scale*EndGaps5[I],90-Scale*StartGaps5[I]);
  };
  fprintf (OFi,"TL2\n");
  for (I=0;I<=NGaps1;++I)    {
    fprintf (OFi,"0 0 24200 %f %f arc S\n",90-Scale*EndGaps1[I],90-Scale*StartGaps1[I]);
  };
  fprintf (OFi,"0 0.0 1 setrgbcolor\n");
  fprintf (OFi,"BL\n");
  for (I=0;I<=NClusters5;++I)    {
    fprintf (OFi,"0 0 24200 %f %f arc S\n",90-Scale*End5[I],90-Scale*Start5[I]);
  };
  fprintf (OFi,"TL2\n");
  for (I=0;I<=NClusters1;++I)    {
    fprintf (OFi,"0 0 24200 %f %f arc S\n",90-Scale*End1[I],90-Scale*Start1[I]);
  };

  fprintf (OFi,"stroke\n");
  fprintf (OFi,"grestore\n");
  fprintf (OFi,"end\n");
  fprintf (OFi,"showpage\n");
  fclose (OFi);

  strcpy (Line,"ps2pdf ");
  strcat(Line,argv[3]);
  //strcat(Line,".c.ps");
  system (Line);

};




/* Output linear */
{
  FILE *OFi;
  char Proc='%';
  double Scale;
  long I,J,KK,LI,LP,LII;

  strcpy(OFiName,argv[2]);
  strcat(OFiName,".l.ps");
  if ( (OFi=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file.\n");
    exit(0);
  };

  Scale=73000.0/500000.0;
  fprintf (OFi,"%c!PS-Adobe-2.0\n",Proc);
  fprintf (OFi,"%%BoundingBox: 60 70 540 730\n");
  fprintf (OFi,"gsave\n");
  fprintf (OFi,"1 dict begin\n");
  fprintf (OFi,"/showpage {} def\n");
  fprintf (OFi,"/lwidth 100 def\n");
  fprintf (OFi,"/vshift -23 def\n");
  fprintf (OFi,"/S { stroke } def\n");
  fprintf (OFi,"/Lshow { 10 10 scale 0 vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/Rshow { 10 10 scale dup stringwidth pop neg vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/Cshow { 10 10 scale dup stringwidth pop -2 div vshift rmoveto show 0.1 0.1 scale } def\n");
  fprintf (OFi,"/WL { S lwidth 0.8 mul setlinewidth } def\n");
  fprintf (OFi,"/NL { S lwidth setlinewidth } def\n");
  fprintf (OFi,"/BL { S lwidth 1.6 mul setlinewidth } def\n");
  fprintf (OFi,"/BBL { S lwidth 3 mul setlinewidth } def\n");
  fprintf (OFi,"/TL0 { S lwidth 2.5 mul setlinewidth } def\n");
  fprintf (OFi,"/TL1 { S lwidth 4 mul setlinewidth } def\n");
  fprintf (OFi,"/TL2 { S lwidth 6 mul setlinewidth } def\n");
  fprintf (OFi,"/TL3 { S lwidth 8 mul setlinewidth } def\n");
  fprintf (OFi,"/M {moveto} def\n");
  fprintf (OFi,"/T {translate} def\n");
  fprintf (OFi,"/L {lineto} def\n");
  fprintf (OFi,"/RL { rlineto } def\n");
  fprintf (OFi,"0.01 0.01 scale\n");
  fprintf (OFi,"90 rotate\n");
  fprintf (OFi,"3000 -7000 T\n");
  fprintf (OFi,"0.00 setgray\n");
/*  fprintf (OFi,"/Helvetica findfont 220 scalefont setfont\n");
  fprintf (OFi,"37000 3000 M\n");
  fprintf (OFi,"(%s) Cshow S\n",Genome);*/
  fprintf (OFi,"NL\n");
  MakeLine (0,LSeq,0,Scale,OFi);
  fprintf (OFi,"WL\n");
  fprintf (OFi,"0 -200 M\n");
  fprintf (OFi,"0 +200 L S\n");
  for (I=1;I*500000<=LSeq;++I)    {
    fprintf (OFi,"0 %li M\n",300-I*2500);
    fprintf (OFi,"0 %li L S\n",-300-I*2500);
    fprintf (OFi,"73000 %li M\n",300-(I-1)*2500);
    fprintf (OFi,"73000 %li L S\n",-300-(I-1)*2500);
  };
  fprintf (OFi,"TL0\n");
  for (I=10000;I<=LSeq;I+=10000)    {
    MakeLine(I-80,I+80,0,Scale,OFi);
  };
  fprintf (OFi,"TL1\n");
  for (I=50000;I<=LSeq;I+=50000)    {
    MakeLine(I-150,I+150,0,Scale,OFi);
  };
  fprintf (OFi,"TL2\n");
  for (I=100000;I<=LSeq;I+=100000)    {
    MakeLine(I-150,I+150,0,Scale,OFi);
  };

  fprintf (OFi,"0.00 setgray\n");
  fprintf (OFi,"BL\n");
  KK=0;
  while (KK<=NMarkers)    {
    LI=MkrStart[KK];
    LII=MkrEnd[KK];
    MakeLine(LI,LII,600,Scale,OFi);
    ++KK;
  };

  LP=-800;
  fprintf (OFi,"1 0.75 0.25 setrgbcolor\n");
  fprintf (OFi,"BL\n");
  for (I=0;I<=NGaps5;++I)    {
    MakeLine(StartGaps5[I],EndGaps5[I],LP,Scale,OFi);
  };
  fprintf (OFi,"BBL\n");
  for (I=0;I<=NGaps1;++I)    {
    MakeLine(StartGaps1[I],EndGaps1[I],LP,Scale,OFi);
  };
  fprintf (OFi,"0 0.0 1 setrgbcolor\n");
  fprintf (OFi,"BL\n");
  for (I=0;I<=NClusters5;++I)    {
    MakeLine(Start5[I],End5[I],LP,Scale,OFi);
  };
  fprintf (OFi,"BBL\n");
  for (I=0;I<=NClusters1;++I)    {
    MakeLine(Start1[I],End1[I],LP,Scale,OFi);
  };


  fprintf (OFi,"stroke\n");
  fprintf (OFi,"grestore\n");
  fprintf (OFi,"end\n");
  fprintf (OFi,"showpage\n");
  fclose (OFi);

  strcpy (Line,"ps2pdf ");
  strcat(Line,argv[2]);
  strcat(Line,".l.ps");
  system (Line);

};


}

