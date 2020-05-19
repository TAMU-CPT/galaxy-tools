/*********************************************
MOTIF LOCATOR
Copyright 2007 Jan Mrazek (mrazek@uga.edu) 
*********************************************/

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
/* Compile with a command "gcc motloc.c -lm -O3 -o ~/bin/motloc" or equivalent. */
/* cc compiler should work, too.                                                */
/* You may leave out the -O3 option (optimization) if it's causing problems     */
/* The program was tested on Redhat Linux.                                      */
/* Run as                                                                       */
/* "motloc <SequnceFile> <AlignmentFile> <OutputFile> <MinScore> <overlap> <strand>", */
/* where <SequenceFile> is an input sequence in GenBank or FASTA format,        */
/* <AlignmentFile> contains the aligned sequences (and nothing else), one       */
/* pattren per line, no spaces, <OutputFile> is the output file name            */
/* (extensions .pll and .plq  will be added by the program), <MinScore>         */
/* specifies minimum score with 0 defaulting to the 10th percentile among       */
/* scores for all sequences in the initial alignment but not lower than 0 and   */
/* not higher than log2(sequence length),                                       */
/* <overlap> 1 or 0 specifies treatment of overlapping matches (1 - leave,      */
/* 0 combine), and <strand> tells the program to search bothy DNA strands (1)   */
/* or only the direct strand (0).                                               */
/********************************************************************************/

/********************************************************************************
History and modifications of this program:

Created on 9/27/2007 by Jan Mrazek (mrazek@uga.edu) 
Modified on 1/9/2009 by Jan Mrazek (mrazek@uga.edu) 


*********************************************************************************/



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXLINE 1100
#define MAXSEQLEN 12000000      /* maximum length of the DNA sequence */
#define MAXMOTIFLENGTH 1000     /* maximum length of the aligned motifs */
#define MAXPATTERNS 200000      /* maximum number of patterns in the sequence. If more patterns are found,
                                   the program will stop and print a message that the pattern is too general. */
#define MAXALIGNED 50000        /* maximum number of sequences in the input alignment */


/*  Global declarations  */
static char Seq[MAXSEQLEN+2];
long LSeq;
char NTDef[16][5]=
        { {0,0,0,0,0},
          {1,1,0,0,0},    /* A (code 1) comprises 1 nucleotide A (code 1) */
          {1,2,0,0,0},    /* C (code 2) comprises 1 nucleotide C (code 2) */
          {1,3,0,0,0},    /* G (code 3) comprises 1 nucleotide G (code 3) */
          {1,4,0,0,0},    /* T (code 4) comprises 1 nucleotide T (code 4) */
          {2,2,3,0,0},    /* S (code 5) comprises 2 nucleotides C,G (codes 2,3) */
          {2,1,4,0,0},    /* W (code 6) comprises 2 nucleotides A,T (codes 1,4) */
          {2,1,3,0,0},    /* R (code 7) comprises 2 nucleotides A,G (codes 1,3) */
          {2,1,2,0,0},    /* M (code 8) comprises 2 nucleotides A,C (codes 1,2) */
          {2,3,4,0,0},    /* K (code 9) comprises 2 nucleotides G,T (codes 3,4) */
          {2,2,4,0,0},    /* Y (code 10) comprises 2 nucleotides C,T (codes 2,4) */
          {3,2,3,4,0},    /* B (code 11) comprises 3 nucleotides C,G,T (codes 2,3,4) */
          {3,1,3,4,0},    /* D (code 12) comprises 3 nucleotides A,G,T (codes 1,3,4) */
          {3,1,2,4,0},    /* H (code 13) comprises 3 nucleotides A,C,T (codes 1,2,4) */
          {3,1,2,3,0},    /* V (code 14) comprises 3 nucleotides A,C,G (codes 1,2,3) */
          {4,1,2,3,4} };  /* N (code 15) comprises 4 nucleotides A,C,G,T (codes 1,2,3,4) */



char NTtoI ( NT )
char NT;
{
  /* assigns integer numbers to nucleotides acording to the following code:
     A  ...  1
     C  ...  2
     G  ...  3
     T  ...  4
     S  ...  5
     W  ...  6
     R  ...  7
     M  ...  8
     K  ...  9
     Y  ... 10
     B  ... 11
     D  ... 12
     H  ... 13
     V  ... 14
     N  ... 15
     any other character is assigned 0. */

  if ( NT=='T' || NT=='t' ) return(4);
  if ( NT=='G' || NT=='g' ) return(3);
  if ( NT=='C' || NT=='c' ) return(2);
  if ( NT=='A' || NT=='a' ) return(1);
  if ( NT=='S' || NT=='s' ) return(5);
  if ( NT=='W' || NT=='w' ) return(6);
  if ( NT=='R' || NT=='r' ) return(7);
  if ( NT=='M' || NT=='m' ) return(8);
  if ( NT=='K' || NT=='k' ) return(9);
  if ( NT=='Y' || NT=='y' ) return(10);
  if ( NT=='B' || NT=='b' ) return(11);
  if ( NT=='D' || NT=='d' ) return(12);
  if ( NT=='H' || NT=='h' ) return(13);
  if ( NT=='V' || NT=='v' ) return(14);
  if ( NT=='N' || NT=='n' ) return(15);
  return(0);
}




char NTCompl ( NT )
char NT;
{
  /* Returns a nucleotide code complementary to the parameter */
  if (NT==0)    return(0);
  if (NT<=4)    return(5-NT);
  if (NT<=6)    return(NT);
  if (NT<=10)   return(17-NT);
  if (NT<=14)   return(25-NT);
  return(15);
}




int RSNTisNT (S,M)
char S;
int M;
{
  /* compares nucleotide S (A,C,G, or T) to nucleotide M (A,C,G,T,W,S,R,Y,M,K,B,D,H,V, or N) 
     Written by: Jan Mrazek, mrazek@uga.edu  */
  int I;
  if (M==15)    return(1);
  if (S==M)    return(1);
  for (I=1;I<=NTDef[M][0];++I)    {
    if (S==NTDef[M][I])    return(1);
  };
  return(0);
}



long IsMotif (L,S,LMot,MinScore,Score)
long L;
double S[MAXMOTIFLENGTH][5];
long LMot;
double MinScore;
double *Score;
{
  long LL;
  *Score=0;
  for (LL=0;LL<=LMot;++LL)    {
    if (Seq[L+LL]>0)    *Score+=S[LL][Seq[L+LL]];
  };
  if (*Score>=MinScore)    {
    return (L+LMot);
  }
  else    {
    return (0);
  };
}

long IsMotifC (L,S,LMot,MinScore,Score)
long L;
double S[MAXMOTIFLENGTH][5];
long LMot;
double MinScore;
double *Score;
{
  long LL;
  *Score=0;
  for (LL=0;LL<=LMot;++LL)    {
    if (Seq[L-LL]>0)    *Score+=S[LL][5-Seq[L-LL]];
  };
  if (*Score>=MinScore)    {
    return (L-LMot);
  }
  else    {
    return (0);
  };
}
      




int main (argc, argv)
int argc;
char *argv[];
{
  char Line[MAXLINE],PLine[MAXLINE],OFiName[MAXLINE],PatFiName[MAXLINE],SeqFiName[MAXLINE];
  FILE *IFi;      /* input file */
  FILE *OFiM;     /* list of pattern locations */
  FILE *OFiMS;    /* pattern sequences and flanks */
  long MkrStart[MAXPATTERNS],MkrEnd[MAXPATTERNS],PMkrStart[MAXPATTERNS],PMkrEnd[MAXPATTERNS];
  double F[MAXMOTIFLENGTH][5],S[MAXMOTIFLENGTH][5],P,PP,MinScore,PScore[MAXALIGNED+1],MinScorePercentile;
  int I,J,K,NMkr,IP,JP,JName,IMkr,IMkr0;
  long L,NPatterns,LP,LPP,LL,LI,LJ,LK,KK,NErr,NNOP;
  char Z;
  char IQStr,IQSeqPrint,IQOverlap,IQMinScorePercentile;
  char LineSeq[71],LineNotes[71],NT[17];
  long NNT[5],LMot,LastEnd;
  double FNT[5],MkrScore[MAXPATTERNS],PMkrScore[MAXPATTERNS];


  if (argc!=8)    {
    printf ("\n\n\nUsage:  motloc SequnceFile AlignmentFile OutputFile MinScore OverlapingPatterns Strand\n\n");
    printf ("SequenceFile: name of input file containing the analyzed sequence in GenBank or FASTA format.\n");
    printf ("AlignmentFile: name of input file containing the aligned motifs (1 motif per line and nothing else)\n");
    printf ("OutputFile: base for output file names (extensions will be added)\n");
    printf ("MinScore: minimum score of the matching motifs (0 will convert to default value)\n");
    printf ("OverlapingPatterns: 1=leave; 0=combine\n");
    printf ("Strand: 1=both strands; 0=only direct strand\n\n\n");
    exit(0);
  };

  strcpy (SeqFiName,argv[1]);
  strcpy (PatFiName,argv[2]);
  MinScore=atof(argv[5]);
  IQMinScorePercentile=0;
  MinScorePercentile=0.1;
  I=0;
  J=0;
  while (argv[5][I]!='\0')    {
    if (argv[5][I]=='%')    J=1;
printf ("%i %i\n",I,J);
    ++I;
  };
  if (J==1)    {
    IQMinScorePercentile=1;
    if (argv[5][0]=='%')    {
      MinScorePercentile=atof(argv[5]+1)/100.0;
    }
    else    {
      MinScorePercentile=atof(argv[5])/100.0;
    };
    if (MinScorePercentile<0.0)    MinScorePercentile=0.0;
    if (MinScorePercentile>1.0)    MinScorePercentile=1.0;
  };
printf ("%i %f\n",IQMinScorePercentile,MinScorePercentile);
  IQOverlap=atoi(argv[6]);
  IQStr=atoi(argv[7]);

  /* reading sequence */
  if ( (IFi=fopen(SeqFiName,"r")) == NULL )  {
    printf ("Unable to open input file %s.\n",SeqFiName);
    exit(0);
  };
  if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
    printf ("Unable to read from file %s.\n",SeqFiName);
  };
  IQSeqPrint=0;
  if (Line[0]=='>')    {
    printf ("Reading FASTA format\n");
    LSeq=0;
    Z=fgetc(IFi);
    while (Z!=EOF && IQSeqPrint==0)    {
      Z=NTtoI(Z);
      if (Z>0)    {
        if (LSeq<MAXSEQLEN)    {
          ++LSeq;
          Seq[LSeq]=Z;
          if (Z>4)    Seq[LSeq]=0;    }
        else    {
          printf ("Sequence is too long. Only the first %li bp will be searched.\n",LSeq);
          printf ("You can change MAXSEQLEN and recompile the program to allow searches in longer sequences.\n");
          IQSeqPrint=1;
        };
      };
      Z=fgetc(IFi);
    };
  }
  else if (Line[0]=='L' && Line[1]=='O' && Line[2]=='C' && Line[3]=='U' && Line[4]=='S' && Line[5]==' ' && Line[6]==' ')    {
    while (!feof(IFi) && (Line[0]!='O' || Line[1]!='R' || Line[2]!='I' || Line[3]!='G' || Line[4]!='I' || Line[5]!='N') )    {
      fgets(Line,MAXLINE-1,IFi);
    };
    if (feof(IFi))    {
      printf ("The input file %s does not appear to be in FASTA or GenBank format\n\n");
      printf ("FASTA: First line has to start with '>', sequence starts from second line\n");
      printf ("GenBank: First line has to start with 'LOCUS', sequence starts after a line starting with 'ORIGIN'\n");
      exit(0);
    }
    else    {
      printf ("Reading GenBank format\n");
      LSeq=0;
      Z=fgetc(IFi);
      while (Z!=EOF && IQSeqPrint==0)    {
        Z=NTtoI(Z);
        if (Z>0)    {
          if (LSeq<MAXSEQLEN)    {
            ++LSeq;
            Seq[LSeq]=Z;
            if (Z>4)    Seq[LSeq]=0;
          }
          else if (IQSeqPrint==0)    {
            printf ("Sequence is to long. Only the first %li bp will be searched.\n",LSeq);
            IQSeqPrint=1;
          };
        };
        Z=fgetc(IFi);
      };
    };
  }
  else    {
    printf ("The input file %s does not appear to be in FASTA or GenBank format\n\n");
    printf ("FASTA: First line has to start with '>', sequence starts from second line\n");
    printf ("GenBank: First line has to start with 'LOCUS', sequence starts after a line starting with 'ORIGIN'\n");
    exit(0);
  };
  fclose (IFi);


  /* Counting nucleotides (filling FNT[]) */
  for (L=1;L<=LSeq;++L)    {
    ++NNT[Seq[L]];
  };
  L=NNT[1]+NNT[2]+NNT[3]+NNT[4];
  FNT[1]=(double)NNT[1]/(double)(L);
  FNT[2]=(double)NNT[2]/(double)(L);
  FNT[3]=(double)NNT[3]/(double)(L);
  FNT[4]=(double)NNT[4]/(double)(L);



  /* Generating names and opening output files */
  strcpy (OFiName,argv[3]);
  //strcat (OFiName,".mll");
  if ( (OFiM=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file %s.\n",OFiName);
  };
  strcpy (OFiName,argv[4]);
  //strcat (OFiName,".mlq");
  if ( (OFiMS=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file %s.\n",OFiName);
  };



  /* reading alignment - filling weight matrix F[] and score matrix S[] */
  if ( (IFi=fopen(PatFiName,"r")) == NULL )  {
    printf ("Unable to open input file %s.\n",PatFiName);
    exit(0);
  };
  if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
    printf ("Unable to read from file %s.\n",PatFiName);
    exit(0);
  };


  LMot=-1;
  while (strchr("ACGTRYSWMKBDHVNacgtryswmkbdhvn",Line[LMot+1])!=NULL)    {
    ++LMot;
    for (I=0;I<=4;++I)    F[LMot][I]=FNT[I];   // pseudocounts equal to background frequencies
//    for (I=0;I<=4;++I)    F[LMot][I]=0.0;   // no pseudocounts
    Z=NTtoI(Line[LMot]);
    if (Z>4)    Z=0;
    F[LMot][Z]+=1.0;
  };

  printf ("Motif length: %li\n",LMot+1);
  if (LMot<3)    {
    printf ("!!! Error:  Motif too short  !!!\n");
    exit (0);
  };

  fgets(Line,MAXLINE-1,IFi);
  while (!feof(IFi))    {
    for (I=0;I<=LMot;++I)    {
      if (strchr("ACGTRYSWMKBDHVNacgtryswmkbdhvn",Line[I])==NULL)    {
        printf ("!!! Error:  Aligned sequences not the same lentgth  !!!\n");
        exit (0);
      };
      Z=NTtoI(Line[I]);
      if (Z>4)    Z=0;
      F[I][Z]+=1.0;
    };
    fgets(Line,MAXLINE-1,IFi);
  };
  for (I=0;I<=LMot;++I)    {
    P=F[I][1]+F[I][2]+F[I][3]+F[I][4];
    F[I][1]=F[I][1]/P;
    F[I][2]=F[I][2]/P;
    F[I][3]=F[I][3]/P;
    F[I][4]=F[I][4]/P;
    S[I][1]=log2(F[I][1]/FNT[1]);
    S[I][2]=log2(F[I][2]/FNT[2]);
    S[I][3]=log2(F[I][3]/FNT[3]);
    S[I][4]=log2(F[I][4]/FNT[4]);
  };

  fclose (IFi);

  /* Default minimum score or percentile minimum score */
  if (MinScore<=0.0 || IQMinScorePercentile==1)    {
    if ( (IFi=fopen(PatFiName,"r")) == NULL )  {
      printf ("Unable to open input file %s.\n",PatFiName);
      exit(0);
    };
    L=0;
    fgets(Line,MAXLINE-1,IFi);
    while (!feof(IFi) && L<MAXALIGNED)    {
      P=0;
      for (I=0;I<=LMot;++I)    {
        Z=NTtoI(Line[I]);
        if (Z>4)    Z=0;
        P+=S[I][Z];
      };
      ++L;
      PScore[L]=P;
      fgets(Line,MAXLINE-1,IFi);
    };
    fclose (IFi);
    for (LI=2;LI<=L;++LI)    {
      LJ=LI;
      while (LJ>=2 && PScore[LJ]<PScore[LJ-1])    {
        P=PScore[LJ];
        PScore[LJ]=PScore[LJ-1];
        PScore[LJ-1]=P;
        --LJ;
      };
    };
    LL=(double)(L)*MinScorePercentile;
    if (LL>L-1)    LL=L-1;
    MinScore=PScore[LL+1];
    if (IQMinScorePercentile==0)    {
      PP=0.0;
      for (I=0;I<=LMot;++I)    {
        P=S[I][1];
        if (P<S[I][2])    P=S[I][2];
        if (P<S[I][3])    P=S[I][3];
        if (P<S[I][4])    P=S[I][4];
        PP+=P;
      };
      PP=PP*0.0;
      if (PP>MinScore)    MinScore=PP;
      P=log2((double)(LSeq));
      if (P<MinScore)    MinScore=P;
    };
  };



  if (OFiM!=NULL)     {
    if (IQOverlap==1)    {
      fprintf (OFiM,"Locations of motifs in the analyzed sequence\nAnalyzed sequence: %s   (%li nucleotides)\n\n",SeqFiName,LSeq);
    }
    else    {
      fprintf (OFiM,"Locations of motifs in the analyzed sequence\nAnalyzed sequence: %s   (%li nucleotides)\n",SeqFiName,LSeq);
      fprintf (OFiM,"     Combining overlaping motifs ...\n");
    };
    fprintf (OFiM,"\n\nMotif scores and frequencies:\n\n");
    fprintf (OFiM,"Pos.  ------------Scores-------------   ------Frequencies------\n");
    fprintf (OFiM,"         A       C       G       T        A     C     G     T  \n");
    for (I=0;I<=LMot;++I)    {
      fprintf (OFiM,"%3i: %8.3f%8.3f%8.3f%8.3f  %6.3f%6.3f%6.3f%6.3f\n",I+1,S[I][1],S[I][2],S[I][3],S[I][4],F[I][1],F[I][2],F[I][3],F[I][4]);
    };
    fprintf (OFiM,"\nMinimum score: %f\n\n\n",MinScore);
    fprintf (OFiM,"     Start       End\n");
  };


  printf ("Searching for motifs ...\n");
  /* searching for motifs in the analyzed sequence - filling MkrStart[] and MkrEnd[] */
  NPatterns=-1;
  LastEnd=0;
    printf ("Direct strand ...\n");
    for (L=1;L<=LSeq-LMot;++L)    {
          if ( L%100000==0 || L==LSeq )    printf ("%li kb ,   %li patterns found\n",L/1000,NPatterns+1);
          LP=IsMotif(L,S,LMot,MinScore,&P);
          if (IQOverlap==0)    {
            if (LP>LastEnd)    {
              if (L<=LastEnd)    {
                MkrEnd[NPatterns]=LP;
              }
              else    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrStart[NPatterns]=L;
                MkrEnd[NPatterns]=LP;
              };
              LastEnd=LP;
            };
          }
          else    {
            if (LP>0)    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrStart[NPatterns]=L;
                MkrEnd[NPatterns]=LP;
                fprintf (OFiM,"%10i%10i\n",L,LP);
            };
          };  
    };

  if (IQStr==1)    {
    printf ("Complementary strand ...\n");
    LastEnd=LSeq;
    for (L=LSeq;L>=1+LMot;--L)    {
          if ( (LSeq-L)%100000==0 || L==1 )    printf ("%li kb ,   %li patterns found\n",(LSeq-L)/1000,NPatterns+1);
          LP=IsMotifC(L,S,LMot,MinScore,&P);
          if (IQOverlap==0)    {
            if (LP<LastEnd && LP>0)    {
              if (L>=LastEnd)    {
                MkrStart[NPatterns]=LP;
              }
              else    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrEnd[NPatterns]=L;
                MkrStart[NPatterns]=LP;
              };
              LastEnd=LP;
            };
          }
          else    {
            if (LP>0)    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrEnd[NPatterns]=L;
                MkrStart[NPatterns]=LP;
                fprintf (OFiM,"%10i%10i\n",LP,L);
            };
          };  
    };
  };
  if (IQOverlap==0)    {
    for (L=0;L<=NPatterns;++L)    {
      fprintf (OFiM,"%10i%10i\n",MkrStart[L],MkrEnd[L]);
    };
  };



/* sorting pattern locations */
  for (L=1;L<=NPatterns;++L)    {
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

  /* NPatterns+1 patterns were found and their starting and ending
     positions in are now in MkrStart[] and MkrEnd[] */

  /* counting non-overlaping patterns */
  for (LP=0;LP<=NPatterns;++LP)    {
    PMkrStart[LP]=MkrStart[LP];
    PMkrEnd[LP]=MkrEnd[LP];
  };
  NNOP=NPatterns;
  for (LP=1;LP<=NNOP;++LP)    {
    while (PMkrStart[LP]<=PMkrEnd[LP-1] && LP<=NNOP)    {
      PMkrEnd[LP-1]=PMkrEnd[LP];
      --NNOP;
      for (LPP=LP;LPP<=NNOP;++LPP)    {
        PMkrStart[LPP]=PMkrStart[LPP+1];
        PMkrEnd[LPP]=PMkrEnd[LPP+1];
      };
    };
  };


  /* Printing sequences and their flanks */
  LineSeq[70]='\0';
  LineNotes[70]='\0';
  strcpy(NT,"NACGTSWRMKYBDHVN");
  for (I=1;I<=50;++I)    Seq[LSeq+I]=0;
  fprintf (OFiMS,"Motif sequences and flanks\nAnalyzed sequence: %s   (%li nucleotides)\n\n",SeqFiName,LSeq);
  fprintf (OFiMS,"\n\nMotif scores and frequencies:\n\n");
  fprintf (OFiMS,"Pos.  ------------Scores-------------   ------Frequencies------\n");
  fprintf (OFiMS,"         A       C       G       T        A     C     G     T  \n");
  for (I=0;I<=LMot;++I)    {
    fprintf (OFiMS,"%3i: %8.3f%8.3f%8.3f%8.3f  %6.3f%6.3f%6.3f%6.3f\n",I+1,S[I][1],S[I][2],S[I][3],S[I][4],F[I][1],F[I][2],F[I][3],F[I][4]);
  };
  fprintf (OFiMS,"\nMinimum score: %f\n\n\n",MinScore);
  fprintf (OFiMS,"     Start       End\n");
  fprintf (OFiMS,"\nNumber of motifs in the analyzed sequence :  %li\n",NPatterns+1);
  fprintf (OFiMS,"Number of non-overlaping motifs :  %li\n\n\n",NNOP+1);
  LI=-1;
  while (LI<NPatterns)    {
    ++LI;
    LJ=LI;
    LK=MkrEnd[LI];
    while (LJ<NPatterns && MkrStart[LJ+1]<=LK)    {
      ++LJ;
      if (MkrEnd[LJ]>LK)    LK=MkrEnd[LJ];
    };
    fprintf (OFiMS,"\n");
    LP=10*(floor((double)(MkrStart[LI])/10.0)-1)+1;
    LPP=10*(floor((double)(LK+9-1)/10.0)+1);
    LL=(LPP-LP+1)%50;
    L=(MkrStart[LI]-LP)-(LPP-LK-1);
    if (LL==40)    {
      if (L>0)    {
        LPP+=10;    }
      else    {
        LP-=10;
      };    }
    else if (LL==30)    {
      LPP+=10;
      LP-=10;    }
    else if (LL==20)    {
      if (L>0)    {
        LPP+=10;    }
      else    {
        LP-=10;
      };
      LPP+=10;
      LP-=10;    }
    else if (LL==10)    {
      LPP+=20;
      LP-=20;
    };
    if (LP<1)    LP=1;
    if (LPP>LSeq)    LPP=LSeq;
    for (L=LP;L<=LPP;L+=50)    {
      sprintf(LineSeq,"%10li",L);
      for (I=10;I<=69;++I)    {
        LineSeq[I]=' ';
      };
      for (I=0;I<=69;++I)    {
        LineNotes[I]=' ';
      };
      LL=L-1;
      for (I=15;I<=60;I+=11)    {
        for (J=I;J<=I+9;++J)    {
          ++LL;
          LineSeq[J]=NT[Seq[LL]];
        };
      };
      for (LL=LI;LL<=LJ;++LL)    {
        if ( MkrStart[LL]>=L && MkrStart[LL]<=L+49 )    {
          I=MkrStart[LL]-L;
          I=15+I+(I/10);
          LineNotes[I]='>';
        };
      };
      for (LL=LI;LL<=LJ;++LL)    {
        if ( MkrEnd[LL]>=L && MkrEnd[LL]<=L+49 )    {
          I=MkrEnd[LL]-L;
          I=15+I+(I/10);
          if (LineNotes[I]=='>')    {
            LineNotes[I]='X';    }
          else    {
            LineNotes[I]='<';
          };
        };
      };
      fprintf (OFiMS,"%s\n%s\n",LineSeq,LineNotes);
    };
    LI=LJ;
  };
  fclose (OFiMS);

}

