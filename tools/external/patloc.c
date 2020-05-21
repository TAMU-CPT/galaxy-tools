/*********************************************
PATTERN LOCATOR
Copyright 2006 Jan Mrazek (mrazek@uga.edu) 
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
/* Please cite Mrazek, J. and Xie, S. (manuscript in preparation) if    */
/* you use this program in a published work.                            */
/************************************************************************/

/********************************************************************************/
/* INSTALLATION INSTRUCTIONS:                                                   */
/* Compile with a command "gcc patloc.c -lm -O3 -o ~/bin/patloc" or equivalent. */
/* cc compiler should work, too.                                                */
/* You may leave out the -O3 option (optimization) if it's causing problems     */
/* The program was tested on Redhat Linux.                                      */
/* Run as "patloc <SequnceFile> <PatternFile> <OutputFile> <overlap>", where    */
/* <SequenceFile> is an input sequence in GenBank or FASTA format,              */
/* <PatternFile> contains the pattern descriptions, one pattren per line,       */
/* no spaces, <OutputFile> is the output file name (extensions .pll and         */
/* .plq  will be added by the program), and <overlap> 1 or 0 specifies          */
/* treatment of overlapping patterns (1 - leave, 0 combine).                    */
/********************************************************************************/

/********************************************************************************
History and modifications of this program:

Modified on 7/14/2006 by Jan Mrazek (mrazek@uga.edu) 
   Minor mmodifications

Modified on 8/25/2006 by Jan Mrazek (mrazek@uga.edu) 
   Added the # symbol in the pattern definition + minor modifications

Modified on 3/31/2008 by Jan Mrazek (mrazek@uga.edu) 
   Changed declaration of large arrays to static to avoid stack oveflow
*********************************************************************************/



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXLINE 600
#define MAXSEQLEN 12000000      /* maximum length of the DNA sequence */
#define MAXPATTERNLENGTH 500    /* maximum length of the pattern description */
#define MAXPATTERNS 1000000     /* maximum number of patterns in the sequence. If more patterns are found,
                                   the program will stop and print a message that the pattern is too general. */
#define MAXNUMBERPATTERNS 200   /* maximum number of different patterns that can be searched in a single run of the program */
#define MAXGOODENDS 5000        /* This plays a role in the exhaustive search with variable-length repetitive segments 
                                   and allowing mismatches. 5000 should be more than sufficient even for complex patterns.
                                   (if it ever exceeds 5000 you probably don't want to wait for the result anyway) */


/*  Global declarations  */
char Seq[MAXSEQLEN+2];
static long LSeq;
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




long IsWord (L,MkrSeq,L0,MkrSeq0,MaxErr,NErr,NGoodEnds,GoodEnds,GoodEndsErr)
long L;
int *MkrSeq;
long L0;
int *MkrSeq0;
long MaxErr;
long *NErr;
long *NGoodEnds;
long *GoodEnds;
long *GoodEndsErr;

{
  /* if a pattern pointed to by MkrSeq (sequence of integers ending with 20000) starts at position L in Seq,
     then the ending position of the word in Seq is returned, otherwise 0. L0 is the original starting
     position or reference point as L can change in recursive calls. Allows maximum MaxErr errors (single base mismatches). 
     L,MkrSeq,L0,MkrSeq0,MaxErr are input parameters, NErr,NGoodEnds,GoodEnds,GoodEndsErr are output parameters.

     Written by: Jan Mrazek, mrazek@uga.edu  


Limitations:

- repeats and complements (+ and - characters in pattern definition) cannot refer to positions >= 1000 (e.g., +1000, -1002 etc. will not work).

- () can be nested with () and {} but {} cannot be nested with {}. For example, constructions like (()), ({}), {()}, or even ((({}){()}){})
  are ok but {{}} or {({})()} are not. 

*/

  int J,B,Len,Min,Max,I,IP,JP,I1,I2,I3,IG;
  long KK,LP,LI,LJ,KKK,PMaxErr,PNErr,PPNErr;
  int PMkrSeq[MAXPATTERNLENGTH];
  long LPP;
  long GoodEndsP[MAXGOODENDS],GoodEndsErrP[MAXGOODENDS],GoodEndsL0P[MAXGOODENDS],NGoodEndsP;
  long GoodEnds0[MAXGOODENDS],GoodEndsErr0[MAXGOODENDS],GoodEndsL00[MAXGOODENDS],NGoodEnds0;
  long GoodEnds1[MAXGOODENDS],GoodEndsErr1[MAXGOODENDS],GoodEndsL01[MAXGOODENDS],NGoodEnds1;
  long GoodEndsL0[MAXGOODENDS];
  int IQ,IQFound;


  *NGoodEnds=0;
  GoodEndsErr[0]=0;
  GoodEndsL0[0]=L0;
  GoodEnds[0]=L-1;
  J=0;
  while (1)    {
    B=MkrSeq[J];
    if ( B<0 )    {                        /* variable number of  repeats */
      Len=-B;
      Min=-MkrSeq[J+1]-1;
      Max=-MkrSeq[J+2]-1;
      JP=J+3;
      for (IP=0;IP<=Len-1;++IP)    {
        PMkrSeq[IP]=MkrSeq[JP+IP];
      };
      PMkrSeq[Len]=20000;    /*end the repetitive subpattern*/

      NGoodEnds0=*NGoodEnds;
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEnds0[IP]=GoodEnds[IP];
        GoodEndsErr0[IP]=GoodEndsErr[IP];
        GoodEndsL00[IP]=GoodEndsL0[IP];
      };


      if ( (PMkrSeq[0]==15 && PMkrSeq[1]==20000) || (PMkrSeq[0]==3002 && PMkrSeq[1]==1 && PMkrSeq[2]==15 && PMkrSeq[3]==20000) )    {     
                                                                                                      /* Special case: (N)[Min:Max] */
        *NGoodEnds=-1;
        for (IP=Min;IP<=Max;++IP)    {
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            ++*NGoodEnds;
            GoodEnds[*NGoodEnds]=GoodEnds0[IG]+IP;
            GoodEndsErr[*NGoodEnds]=GoodEndsErr0[IG];
            GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
          };
        };
        J+=3+Len;
      }
      else    {                                                     /* Everythig else; going from 1 to Min number of repeats*/
        for (I=1;I<=Min;++I)    {
          NGoodEnds1=-1;
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            PMaxErr=MaxErr-GoodEndsErr0[IG];
            KK=IsWord(GoodEnds0[IG]+1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);      /* */
            if (NGoodEndsP>=0)    {
              for (IP=0;IP<=NGoodEndsP;++IP)    {
                ++NGoodEnds1;
                GoodEnds1[NGoodEnds1]=GoodEndsP[IP];
                GoodEndsErr1[NGoodEnds1]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL01[NGoodEnds1]=GoodEndsL00[IG];
              };
              if (NGoodEnds1>=1)    {                /* eliminating duplicates */
                for (I1=1;I1<=NGoodEnds1;++I1)    {
                  I2=I1;
                  while (I2>=1 && (GoodEndsErr1[I2-1]>GoodEndsErr1[I2] || (GoodEndsErr1[I2-1]==GoodEndsErr1[I2] && GoodEnds1[I2-1]<GoodEnds1[I2]) ) )    {
                    LP=GoodEnds1[I2-1];
                    GoodEnds1[I2-1]=GoodEnds1[I2];
                    GoodEnds1[I2]=LP;
                    LP=GoodEndsErr1[I2-1];
                    GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                    GoodEndsErr1[I2]=LP;
                    LP=GoodEndsL01[I2-1];
                    GoodEndsL01[I2-1]=GoodEndsL01[I2];
                    GoodEndsL01[I2]=LP;
                    --I2;
                  };
                };
                for (I1=1;I1<=NGoodEnds1;++I1)    {
                  if (GoodEnds1[I1-1]==GoodEnds1[I1])    {
                    for (I2=I1;I2<=NGoodEnds1;++I2)    {
                      GoodEnds1[I2-1]=GoodEnds1[I2];
                      GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                      GoodEndsL01[I2-1]=GoodEndsL01[I2];
                    };
                    --NGoodEnds1;
                    --I1;
                  };
                };
              };
            };
          };
          NGoodEnds0=NGoodEnds1;
          for (IP=0;IP<=NGoodEnds1;++IP)     {
            GoodEnds0[IP]=GoodEnds1[IP];
            GoodEndsErr0[IP]=GoodEndsErr1[IP];
            GoodEndsL00[IP]=GoodEndsL01[IP];
          };
        };
        if (Min==0)    {
          NGoodEnds1=*NGoodEnds;
          for (IP=0;IP<=*NGoodEnds;++IP)     {
            GoodEnds1[IP]=GoodEnds[IP];
            GoodEndsErr1[IP]=GoodEndsErr[IP];
            GoodEndsL01[IP]=GoodEndsL0[IP];
          };
        };
        *NGoodEnds=NGoodEnds1;
        for (IP=0;IP<=NGoodEnds1;++IP)     {
          GoodEnds[IP]=GoodEnds1[IP];
          GoodEndsErr[IP]=GoodEndsErr1[IP];
          GoodEndsL0[IP]=GoodEndsL01[IP];
        };
  
  
        for (I=Min+1;I<=Max;++I)    {                                /* going from Min+1 to Max number of repeats */
          IQFound=0;
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            PMaxErr=MaxErr-GoodEndsErr0[IG];
            KK=IsWord(GoodEnds0[IG]+1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);      /* */
            if (NGoodEndsP>=0)    {
              IQFound=1;
              for (IP=0;IP<=NGoodEndsP;++IP)    {
                ++NGoodEnds1;
                GoodEnds1[NGoodEnds1]=GoodEndsP[IP];
                GoodEndsErr1[NGoodEnds1]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL01[NGoodEnds1]=GoodEndsL00[IG];
                ++*NGoodEnds;
                GoodEnds[*NGoodEnds]=GoodEndsP[IP];
                GoodEndsErr[*NGoodEnds]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
                if (*NGoodEnds>=1)    {
                  for (I1=1;I1<=*NGoodEnds;++I1)    {
                    I2=I1;
                    while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
                      LP=GoodEnds[I2-1];
                      GoodEnds[I2-1]=GoodEnds[I2];
                      GoodEnds[I2]=LP;
                      LP=GoodEndsErr[I2-1];
                      GoodEndsErr[I2-1]=GoodEndsErr[I2];
                      GoodEndsErr[I2]=LP;
                      LP=GoodEndsL0[I2-1];
                      GoodEndsL0[I2-1]=GoodEndsL0[I2];
                      GoodEndsL0[I2]=LP;
                      --I2;
                    };
                  };
                  for (I1=1;I1<=*NGoodEnds;++I1)    {
                    if (GoodEnds[I1-1]==GoodEnds[I1])    {
                      for (I2=I1;I2<=*NGoodEnds;++I2)    {
                        GoodEnds[I2-1]=GoodEnds[I2];
                        GoodEndsErr[I2-1]=GoodEndsErr[I2];
                        GoodEndsL0[I2-1]=GoodEndsL0[I2];
                      };
                      --*NGoodEnds;
                      --I1;
                    };
                  };
                };
                if (NGoodEnds1>=1)    {
                  for (I1=1;I1<=NGoodEnds1;++I1)    {
                    I2=I1;
                    while (I2>=1 && (GoodEndsErr1[I2-1]>GoodEndsErr1[I2] || (GoodEndsErr1[I2-1]==GoodEndsErr1[I2] && GoodEnds1[I2-1]<GoodEnds1[I2]) ) )    {
                      LP=GoodEnds1[I2-1];
                      GoodEnds1[I2-1]=GoodEnds1[I2];
                      GoodEnds1[I2]=LP;
                      LP=GoodEndsErr1[I2-1];
                      GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                      GoodEndsErr1[I2]=LP;
                      LP=GoodEndsL01[I2-1];
                      GoodEndsL01[I2-1]=GoodEndsL01[I2];
                      GoodEndsL01[I2]=LP;
                      --I2;
                    };
                  };
                  for (I1=1;I1<=NGoodEnds1;++I1)    {
                    if (GoodEnds1[I1-1]==GoodEnds1[I1])    {
                      for (I2=I1;I2<=NGoodEnds1;++I2)    {
                        GoodEnds1[I2-1]=GoodEnds1[I2];
                        GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                        GoodEndsL01[I2-1]=GoodEndsL01[I2];
                      };
                      --NGoodEnds1;
                      --I1;
                    };
                  };
                };
              };
            };
          };
          NGoodEnds0=NGoodEnds1;
          for (IP=0;IP<=NGoodEnds1;++IP)     {
            GoodEnds0[IP]=GoodEnds1[IP];
            GoodEndsErr0[IP]=GoodEndsErr1[IP];
            GoodEndsL00[IP]=GoodEndsL01[IP];
          };
          if (IQFound==0)    {
            I=Max+1;
          };
        };

        /* eliminating duplicates */
        if (*NGoodEnds>=1)    {
          for (I1=1;I1<=*NGoodEnds;++I1)    {
            I2=I1;
            while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
              LP=GoodEnds[I2-1];
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEnds[I2]=LP;
              LP=GoodEndsErr[I2-1];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsErr[I2]=LP;
              LP=GoodEndsL0[I2-1];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
              GoodEndsL0[I2]=LP;
              --I2;
            };
          };
          for (I1=1;I1<=*NGoodEnds;++I1)    {
            if (GoodEnds[I1-1]==GoodEnds[I1])    {
              for (I2=I1;I2<=*NGoodEnds;++I2)    {
                GoodEnds[I2-1]=GoodEnds[I2];
                GoodEndsErr[I2-1]=GoodEndsErr[I2];
                GoodEndsL0[I2-1]=GoodEndsL0[I2];
              };
              --*NGoodEnds;
              --I1;
            };
          };
        };
        if (*NGoodEnds>=0)    {
          J+=3+Len;
        }
        else    {
          return(0);
        };
      };
    }


    else if ( B<=15 )    {                        /* nucleotides (IUPAC code) */
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        ++GoodEnds[IG];
        if (RSNTisNT(Seq[GoodEnds[IG]],B)<=0)    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<=2000 )    {              /* direct repeats (+) */ 
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        ++GoodEnds[IG];
        LP=GoodEndsL0[IG]+B-1001;
        if ( Seq[GoodEnds[IG]]!=Seq[LP] || Seq[L]<=0 )    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<=3000 )    {              /* inverted repeats (-) */ 
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        ++GoodEnds[IG];
        LP=GoodEndsL0[IG]+B-2001;
        if ( Seq[GoodEnds[IG]]!=NTCompl(Seq[LP]) || Seq[L]<=0 )    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<4000 )    {         /* segments with errors */
      Len=B-3001;
      PMaxErr=MkrSeq[J+1]-1;
      JP=J+2;
      for (IP=0;IP<=Len-1;++IP)    {
        PMkrSeq[IP]=MkrSeq[JP+IP];
      };
      PMkrSeq[Len]=20000;

      NGoodEnds0=*NGoodEnds;
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEnds0[IP]=GoodEnds[IP];
        GoodEndsErr0[IP]=GoodEndsErr[IP];
        GoodEndsL00[IP]=GoodEndsL0[IP];
      };

      *NGoodEnds=-1;
      for (IG=0;IG<=NGoodEnds0;++IG)    {
        KK=IsWord(GoodEnds0[IG]+1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);    /* */
        if (NGoodEndsP>=0)    {
          for (IP=0;IP<=NGoodEndsP;++IP)    {
            ++*NGoodEnds;
            GoodEnds[*NGoodEnds]=GoodEndsP[IP];
/*            GoodEndsErr[*NGoodEnds]=GoodEndsErrP[IP];  */
            GoodEndsErr[*NGoodEnds]=0;
            GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
          };
        };
      };

      /* eliminating duplicates */
      if (*NGoodEnds>=1)    {
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          I2=I1;
          while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
            LP=GoodEnds[I2-1];
            GoodEnds[I2-1]=GoodEnds[I2];
            GoodEnds[I2]=LP;
            LP=GoodEndsErr[I2-1];
            GoodEndsErr[I2-1]=GoodEndsErr[I2];
            GoodEndsErr[I2]=LP;
            LP=GoodEndsL0[I2-1];
            GoodEndsL0[I2-1]=GoodEndsL0[I2];
            GoodEndsL0[I2]=LP;
            --I2;
          };
        };
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          if (GoodEnds[I1-1]==GoodEnds[I1])    {
            for (I2=I1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --I1;
          };
        };
      };
      if (*NGoodEnds>=0)    {
        J+=2+Len;
      }
      else    {
        return(0);
      };
    }
    else if ( B==4001 )    {         /* set new reference point (#) */
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEndsL0[IP]=GoodEnds[IP]+1;
      };
      ++J;
    }
    else    {                   /* end of the pattern */
      if (*NGoodEnds>=1)    {
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          I2=I1;
          while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
            LP=GoodEnds[I2-1];
            GoodEnds[I2-1]=GoodEnds[I2];
            GoodEnds[I2]=LP;
            LP=GoodEndsErr[I2-1];
            GoodEndsErr[I2-1]=GoodEndsErr[I2];
            GoodEndsErr[I2]=LP;
            LP=GoodEndsL0[I2-1];
            GoodEndsL0[I2-1]=GoodEndsL0[I2];
            GoodEndsL0[I2]=LP;
            --I2;
          };
        };
      };
      *NErr=GoodEndsErr[0];
      return(GoodEnds[0]);
    };
  };
}




long IsWordC (L,MkrSeq,L0,MkrSeq0,MaxErr,NErr,NGoodEnds,GoodEnds,GoodEndsErr)
long L;
int *MkrSeq;
long L0;
int *MkrSeq0;
long MaxErr;
long *NErr;
long *NGoodEnds;
long *GoodEnds;
long *GoodEndsErr;

{
  /* Same as IsWord but for complementary strand

     Written by: Jan Mrazek, mrazek@uga.edu  
  */

  int J,B,Len,Min,Max,I,IP,JP,I1,I2,I3,IG;
  long KK,LP,LI,LJ,KKK,PMaxErr,PNErr,PPNErr;
  int PMkrSeq[MAXPATTERNLENGTH];
  long LPP;
  long GoodEndsP[MAXGOODENDS],GoodEndsErrP[MAXGOODENDS],GoodEndsL0P[MAXGOODENDS],NGoodEndsP;
  long GoodEnds0[MAXGOODENDS],GoodEndsErr0[MAXGOODENDS],GoodEndsL00[MAXGOODENDS],NGoodEnds0;
  long GoodEnds1[MAXGOODENDS],GoodEndsErr1[MAXGOODENDS],GoodEndsL01[MAXGOODENDS],NGoodEnds1;
  long GoodEndsL0[MAXGOODENDS];
  int IQ,IQFound;


  *NGoodEnds=0;
  GoodEndsErr[0]=0;
  GoodEndsL0[0]=L0;
  GoodEnds[0]=L+1;
  J=0;
  while (1)    {
    B=MkrSeq[J];
    if ( B<0 )    {                        /* variable number of  repeats */
      Len=-B;
      Min=-MkrSeq[J+1]-1;
      Max=-MkrSeq[J+2]-1;
      JP=J+3;
      for (IP=0;IP<=Len-1;++IP)    {
        PMkrSeq[IP]=MkrSeq[JP+IP];
      };
      PMkrSeq[Len]=20000;    /*end the repetitive subpattern*/

      NGoodEnds0=*NGoodEnds;
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEnds0[IP]=GoodEnds[IP];
        GoodEndsErr0[IP]=GoodEndsErr[IP];
        GoodEndsL00[IP]=GoodEndsL0[IP];
      };


      if ( (PMkrSeq[0]==15 && PMkrSeq[1]==20000) || (PMkrSeq[0]==3002 && PMkrSeq[1]==1 && PMkrSeq[2]==15 && PMkrSeq[3]==20000) )    {     
                                                                                                      /* Special case: (N)[Min:Max] */
        *NGoodEnds=-1;
        for (IP=Min;IP<=Max;++IP)    {
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            ++*NGoodEnds;
            GoodEnds[*NGoodEnds]=GoodEnds0[IG]-IP;
            GoodEndsErr[*NGoodEnds]=GoodEndsErr0[IG];
            GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
          };
        };
        J+=3+Len;
      }
      else    {                                                     /* Everythig else */
        for (I=1;I<=Min;++I)    {
          NGoodEnds1=-1;
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            PMaxErr=MaxErr-GoodEndsErr0[IG];
            KK=IsWordC(GoodEnds0[IG]-1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);      /* */
            if (NGoodEndsP>=0)    {
              for (IP=0;IP<=NGoodEndsP;++IP)    {
                ++NGoodEnds1;
                GoodEnds1[NGoodEnds1]=GoodEndsP[IP];
                GoodEndsErr1[NGoodEnds1]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL01[NGoodEnds1]=GoodEndsL00[IG];
              };
              if (NGoodEnds1>=1)    {                /* eliminating duplicates */
                for (I1=1;I1<=NGoodEnds1;++I1)    {
                  I2=I1;
                  while (I2>=1 && (GoodEndsErr1[I2-1]>GoodEndsErr1[I2] || (GoodEndsErr1[I2-1]==GoodEndsErr1[I2] && GoodEnds1[I2-1]<GoodEnds1[I2]) ) )    {
                    LP=GoodEnds1[I2-1];
                    GoodEnds1[I2-1]=GoodEnds1[I2];
                    GoodEnds1[I2]=LP;
                    LP=GoodEndsErr1[I2-1];
                    GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                    GoodEndsErr1[I2]=LP;
                    LP=GoodEndsL01[I2-1];
                    GoodEndsL01[I2-1]=GoodEndsL01[I2];
                    GoodEndsL01[I2]=LP;
                    --I2;
                  };
                };
                for (I1=1;I1<=NGoodEnds1;++I1)    {
                  if (GoodEnds1[I1-1]==GoodEnds1[I1])    {
                    for (I2=I1;I2<=NGoodEnds1;++I2)    {
                      GoodEnds1[I2-1]=GoodEnds1[I2];
                      GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                      GoodEndsL01[I2-1]=GoodEndsL01[I2];
                    };
                    --NGoodEnds1;
                    --I1;
                  };
                };
              };
            };
          };
          NGoodEnds0=NGoodEnds1;
          for (IP=0;IP<=NGoodEnds1;++IP)     {
            GoodEnds0[IP]=GoodEnds1[IP];
            GoodEndsErr0[IP]=GoodEndsErr1[IP];
            GoodEndsL00[IP]=GoodEndsL01[IP];
          };
        };
        if (Min==0)    {
          NGoodEnds1=*NGoodEnds;
          for (IP=0;IP<=*NGoodEnds;++IP)     {
            GoodEnds1[IP]=GoodEnds[IP];
            GoodEndsErr1[IP]=GoodEndsErr[IP];
            GoodEndsL01[IP]=GoodEndsL0[IP];
          };
        };
        *NGoodEnds=NGoodEnds1;
        for (IP=0;IP<=NGoodEnds1;++IP)     {
          GoodEnds[IP]=GoodEnds1[IP];
          GoodEndsErr[IP]=GoodEndsErr1[IP];
          GoodEndsL0[IP]=GoodEndsL01[IP];
        };
  
  
        for (I=Min+1;I<=Max;++I)    {
          IQFound=0;
          for (IG=0;IG<=NGoodEnds0;++IG)    {
            PMaxErr=MaxErr-GoodEndsErr0[IG];
            KK=IsWordC(GoodEnds0[IG]-1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);      /* */
            if (NGoodEndsP>=0)    {
              IQFound=1;
              for (IP=0;IP<=NGoodEndsP;++IP)    {
                ++NGoodEnds1;
                GoodEnds1[NGoodEnds1]=GoodEndsP[IP];
                GoodEndsErr1[NGoodEnds1]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL01[NGoodEnds1]=GoodEndsL00[IG];
                ++*NGoodEnds;
                GoodEnds[*NGoodEnds]=GoodEndsP[IP];
                GoodEndsErr[*NGoodEnds]=GoodEndsErrP[IP]+GoodEndsErr0[IG];
                GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
                if (*NGoodEnds>=1)    {
                  for (I1=1;I1<=*NGoodEnds;++I1)    {
                    I2=I1;
                    while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
                      LP=GoodEnds[I2-1];
                      GoodEnds[I2-1]=GoodEnds[I2];
                      GoodEnds[I2]=LP;
                      LP=GoodEndsErr[I2-1];
                      GoodEndsErr[I2-1]=GoodEndsErr[I2];
                      GoodEndsErr[I2]=LP;
                      LP=GoodEndsL0[I2-1];
                      GoodEndsL0[I2-1]=GoodEndsL0[I2];
                      GoodEndsL0[I2]=LP;
                      --I2;
                    };
                  };
                  for (I1=1;I1<=*NGoodEnds;++I1)    {
                    if (GoodEnds[I1-1]==GoodEnds[I1])    {
                      for (I2=I1;I2<=*NGoodEnds;++I2)    {
                        GoodEnds[I2-1]=GoodEnds[I2];
                        GoodEndsErr[I2-1]=GoodEndsErr[I2];
                        GoodEndsL0[I2-1]=GoodEndsL0[I2];
                      };
                      --*NGoodEnds;
                      --I1;
                    };
                  };
                };
                if (NGoodEnds1>=1)    {
                  for (I1=1;I1<=NGoodEnds1;++I1)    {
                    I2=I1;
                    while (I2>=1 && (GoodEndsErr1[I2-1]>GoodEndsErr1[I2] || (GoodEndsErr1[I2-1]==GoodEndsErr1[I2] && GoodEnds1[I2-1]<GoodEnds1[I2]) ) )    {
                      LP=GoodEnds1[I2-1];
                      GoodEnds1[I2-1]=GoodEnds1[I2];
                      GoodEnds1[I2]=LP;
                      LP=GoodEndsErr1[I2-1];
                      GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                      GoodEndsErr1[I2]=LP;
                      LP=GoodEndsL01[I2-1];
                      GoodEndsL01[I2-1]=GoodEndsL01[I2];
                      GoodEndsL01[I2]=LP;
                      --I2;
                    };
                  };
                  for (I1=1;I1<=NGoodEnds1;++I1)    {
                    if (GoodEnds1[I1-1]==GoodEnds1[I1])    {
                      for (I2=I1;I2<=NGoodEnds1;++I2)    {
                        GoodEnds1[I2-1]=GoodEnds1[I2];
                        GoodEndsErr1[I2-1]=GoodEndsErr1[I2];
                        GoodEndsL01[I2-1]=GoodEndsL01[I2];
                      };
                      --NGoodEnds1;
                      --I1;
                    };
                  };
                };
              };
            };
          };
          NGoodEnds0=NGoodEnds1;
          for (IP=0;IP<=NGoodEnds1;++IP)     {
            GoodEnds0[IP]=GoodEnds1[IP];
            GoodEndsErr0[IP]=GoodEndsErr1[IP];
            GoodEndsL00[IP]=GoodEndsL01[IP];
          };
          if (IQFound==0)    {
            I=Max+1;
          };
        };

        /* eliminating duplicates */
        if (*NGoodEnds>=1)    {
          for (I1=1;I1<=*NGoodEnds;++I1)    {
            I2=I1;
            while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
              LP=GoodEnds[I2-1];
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEnds[I2]=LP;
              LP=GoodEndsErr[I2-1];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsErr[I2]=LP;
              LP=GoodEndsL0[I2-1];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
              GoodEndsL0[I2]=LP;
              --I2;
            };
          };
          for (I1=1;I1<=*NGoodEnds;++I1)    {
            if (GoodEnds[I1-1]==GoodEnds[I1])    {
              for (I2=I1;I2<=*NGoodEnds;++I2)    {
                GoodEnds[I2-1]=GoodEnds[I2];
                GoodEndsErr[I2-1]=GoodEndsErr[I2];
                GoodEndsL0[I2-1]=GoodEndsL0[I2];
              };
              --*NGoodEnds;
              --I1;
            };
          };
        };
        if (*NGoodEnds>=0)    {
          J+=3+Len;
        }
        else    {
          return(0);
        };
      };
    }


    else if ( B<=15 )    {                        /* nucleotides */
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        --GoodEnds[IG];
        if (RSNTisNT(NTCompl(Seq[GoodEnds[IG]]),B)<=0)    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<=2000 )    {              /* direct repeats (+) */ 
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        --GoodEnds[IG];
        LP=GoodEndsL0[IG]-B+1001;
        if ( Seq[GoodEnds[IG]]!=Seq[LP] || Seq[L]<=0 )    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<=3000 )    {              /* inverted repeats (-) */ 
      for (IG=0;IG<=*NGoodEnds;++IG)    {
        --GoodEnds[IG];
        LP=GoodEndsL0[IG]-B+2001;
        if ( Seq[GoodEnds[IG]]!=NTCompl(Seq[LP]) || Seq[L]<=0 )    {
          ++GoodEndsErr[IG];
          if (GoodEndsErr[IG]>MaxErr)     {
            for (I2=IG+1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --IG;
            if (*NGoodEnds<0)    {
              return(0);
            }
          };
        };
      };
      ++J;
    }


    else if ( B<4000 )    {         /* segments with errors */
      Len=B-3001;
      PMaxErr=MkrSeq[J+1]-1;
      JP=J+2;
      for (IP=0;IP<=Len-1;++IP)    {
        PMkrSeq[IP]=MkrSeq[JP+IP];
      };
      PMkrSeq[Len]=20000;

      NGoodEnds0=*NGoodEnds;
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEnds0[IP]=GoodEnds[IP];
        GoodEndsErr0[IP]=GoodEndsErr[IP];
        GoodEndsL00[IP]=GoodEndsL0[IP];
      };

      *NGoodEnds=-1;
      for (IG=0;IG<=NGoodEnds0;++IG)    {
        KK=IsWordC(GoodEnds0[IG]-1,PMkrSeq,GoodEndsL00[IG],MkrSeq0,PMaxErr,&PNErr,&NGoodEndsP,GoodEndsP,GoodEndsErrP);    /* */
        if (NGoodEndsP>=0)    {
          for (IP=0;IP<=NGoodEndsP;++IP)    {
            ++*NGoodEnds;
            GoodEnds[*NGoodEnds]=GoodEndsP[IP];
/*            GoodEndsErr[*NGoodEnds]=GoodEndsErrP[IP];  */
            GoodEndsErr[*NGoodEnds]=0;
            GoodEndsL0[*NGoodEnds]=GoodEndsL00[IG];
          };
        };
      };

      /* eliminating duplicates */
      if (*NGoodEnds>=1)    {
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          I2=I1;
          while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
            LP=GoodEnds[I2-1];
            GoodEnds[I2-1]=GoodEnds[I2];
            GoodEnds[I2]=LP;
            LP=GoodEndsErr[I2-1];
            GoodEndsErr[I2-1]=GoodEndsErr[I2];
            GoodEndsErr[I2]=LP;
            LP=GoodEndsL0[I2-1];
            GoodEndsL0[I2-1]=GoodEndsL0[I2];
            GoodEndsL0[I2]=LP;
            --I2;
          };
        };
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          if (GoodEnds[I1-1]==GoodEnds[I1])    {
            for (I2=I1;I2<=*NGoodEnds;++I2)    {
              GoodEnds[I2-1]=GoodEnds[I2];
              GoodEndsErr[I2-1]=GoodEndsErr[I2];
              GoodEndsL0[I2-1]=GoodEndsL0[I2];
            };
            --*NGoodEnds;
            --I1;
          };
        };
      };
      if (*NGoodEnds>=0)    {
        J+=2+Len;
      }
      else    {
        return(0);
      };
    }
    else if ( B==4001 )    {         /* set new reference point (#) */
      for (IP=0;IP<=*NGoodEnds;++IP)     {
        GoodEndsL0[IP]=GoodEnds[IP]-1;
      };
      ++J;
    }
    else    {                   /* end of the pattern */
      if (*NGoodEnds>=1)    {
        for (I1=1;I1<=*NGoodEnds;++I1)    {
          I2=I1;
          while (I2>=1 && (GoodEndsErr[I2-1]>GoodEndsErr[I2] || (GoodEndsErr[I2-1]==GoodEndsErr[I2] && GoodEnds[I2-1]<GoodEnds[I2]) ) )    {
            LP=GoodEnds[I2-1];
            GoodEnds[I2-1]=GoodEnds[I2];
            GoodEnds[I2]=LP;
            LP=GoodEndsErr[I2-1];
            GoodEndsErr[I2-1]=GoodEndsErr[I2];
            GoodEndsErr[I2]=LP;
            LP=GoodEndsL0[I2-1];
            GoodEndsL0[I2-1]=GoodEndsL0[I2];
            GoodEndsL0[I2]=LP;
            --I2;
          };
        };
      };
      *NErr=GoodEndsErr[0];
      return(GoodEnds[0]);
    };
  };
}







int main (argc, argv)
int argc;
char *argv[];
{
  struct Pattern { int Seq[MAXPATTERNLENGTH]; char Name[MAXLINE]; long LastEnd; char IQStrand; };
  char Line[MAXLINE],PLine[MAXLINE],OFiName[MAXLINE],PatFiName[MAXLINE],SeqFiName[MAXLINE];
  FILE *IFi;      /* input file */
  FILE *OFiM;     /* list of pattern locations */
  FILE *OFiMS;    /* pattern sequences and flanks */
  static long MkrStart[MAXPATTERNS],MkrEnd[MAXPATTERNS],PMkrStart[MAXPATTERNS],PMkrEnd[MAXPATTERNS];
  static char MkrType[MAXPATTERNS];
  struct Pattern Mkr[MAXNUMBERPATTERNS];
  int I,J,K,NMkr,IP,JP,JName,IMkr,IMkr0;
  long L,NPatterns,LP,LPP,LL,LI,LJ,LK,KK,NErr,NNOP;
  char Z;
  char IQStr,IQSeqPrint,IQOverlap;
  int IPar,Par[100],IBra,Bra[100];
  long GoodEnds[MAXGOODENDS],GoodEndsErr[MAXGOODENDS],NGoodEnds;
  char LineSeq[71],LineNotes[71],NT[17];


  if (argc!=6)    {
    printf ("\n\n\nUsage:  patloc SequnceFile PatternFile OutputFile OverlapingPatterns\n\n");
    printf ("SequenceFile: name of file containing the analyzed sequence in GenBank or FASTA format.\n");
    printf ("PatternFile: name of file containing the patterns to be located\n");
    printf ("OutputFile: base for output file names (extensions will be added)\n");
    printf ("OverlapingPatterns: 1=leave 0=combine\n\n\n");
    exit(0);
  };

  strcpy (SeqFiName,argv[1]);
  strcpy (PatFiName,argv[2]);
  IQOverlap=atoi(argv[5]);

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





  /* Generating names and opening output files */
  strcpy (OFiName,argv[3]);
  //strcat (OFiName,".pll");
  if ( (OFiM=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file %s.\n",OFiName);
  };
  strcpy (OFiName,argv[4]);
  //strcat (OFiName,".plq");
  if ( (OFiMS=fopen(OFiName,"w")) == NULL )  {
    printf ("Unable to open output file %s.\n",OFiName);
  };


  /* reading patterns - filling Mkr[].Seq */
  if ( (IFi=fopen(PatFiName,"r")) == NULL )  {
    printf ("Unable to open input file %s.\n",PatFiName);
    exit(0);
  };
  if ( fgets(Line,MAXLINE-1,IFi)==NULL )    {
    printf ("Unable to read from file %s.\n",PatFiName);
    exit(0);
  };

  NMkr=-1;
  IQStr=3;
  while (!feof(IFi) && NMkr<MAXNUMBERPATTERNS)    {
    ++NMkr;
    Mkr[NMkr].LastEnd=0;

    /* converting patterns to numerical codes */
    IPar=-1;
    IBra=-1;
    J=0;
    I=0;
    Mkr[NMkr].IQStrand=0;
    JName=0;
    if (Line[0]=='<')    {
      if (Line[1]=='>')    {
        Mkr[NMkr].IQStrand=2;
        I=2;    }
      else    {
        Mkr[NMkr].IQStrand=1;
        I=1;
      };    }
    else if (Line[0]=='>')    {
      if (Line[1]=='<')    {
        Mkr[NMkr].IQStrand=2;
        I=2;    }
      else    {
        Mkr[NMkr].IQStrand=0;
        I=1;
      };
    };
    if (Mkr[NMkr].IQStrand==2 || IQStr==2)    {
      IQStr=2;
    }
    else if (Mkr[NMkr].IQStrand==0)    {
      if (IQStr==1)    {
        IQStr=2;    }
      else    {
        IQStr=0;
      };
    }
    else if (Mkr[NMkr].IQStrand==1)    {
      if (IQStr==0)    {
        IQStr=2;    }
      else    {
        IQStr=1;
      };
    };

    while ( Line[I]!='\0' && Line[I]!='\n' && Line[I]!=' ' )    {
      Z=NTtoI(Line[I]);
      if (Z>0)    {
        Mkr[NMkr].Seq[J]=Z;
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++J;
        ++I;
      }
      else if (Line[I]=='#')    {
        Mkr[NMkr].Seq[J]=4001;
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++J;
        ++I;
      }
      else if (Line[I]=='(')    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++IPar;
        Par[IPar]=J;
        Mkr[NMkr].Seq[J]=19999;
        Mkr[NMkr].Seq[J+1]=19999;
        Mkr[NMkr].Seq[J+2]=19999;
        J+=3;
        ++I;    }
      else if (Line[I]==')' && IPar>=0)    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
        if (Line[I]!='[')    {
          printf ("Error in pattern no. %i, position %i\n",NMkr+1,I+1);
          ++I;
        }
        else    {
          Mkr[NMkr].Seq[Par[IPar]]=Par[IPar]+3-J;
          Mkr[NMkr].Name[JName]=Line[I];
          ++JName;
          ++I;
          K=1+atoi(Line+I);
          Mkr[NMkr].Seq[Par[IPar]+1]=-K;
          while ( Line[I]!=':' && Line[I]!='\n' )    {
            Mkr[NMkr].Name[JName]=Line[I];
            ++JName;
            ++I;
          };
          Mkr[NMkr].Name[JName]=Line[I];
          ++JName;
          ++I;
          if (Line[I]==']')    {
            K=10001;    }
          else    {
            K=1+atoi(Line+I);
          };
          Mkr[NMkr].Seq[Par[IPar]+2]=-K;
          while ( Line[I]!=']' && Line[I]!='\n' )    {
            Mkr[NMkr].Name[JName]=Line[I];
            ++JName;
            ++I;
          };
        };
        --IPar;
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
      }
      else if (Line[I]=='{')    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++IBra;
        Bra[IBra]=J;
        Mkr[NMkr].Seq[J]=19999;
        Mkr[NMkr].Seq[J+1]=19999;
        J+=2;
        ++I;
      }
      else if (Line[I]=='}' && IBra>=0)    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
        if (Line[I]!='[')    {
          printf ("Error in pattern no. %i, position %i\n",NMkr+1,I+1);
          ++I;
        }
        else    {
          Mkr[NMkr].Seq[Bra[IBra]]=3001+J-Bra[IBra]-2;
          Mkr[NMkr].Name[JName]=Line[I];
          ++JName;
          ++I;
          K=1+atoi(Line+I);
          Mkr[NMkr].Seq[Bra[IBra]+1]=K;
          while ( Line[I]!=']' && Line[I]!='\n' )    {
            Mkr[NMkr].Name[JName]=Line[I];
            ++JName;
            ++I;
          };
        };
        --IBra;
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
      }
      else if (Line[I]=='+')    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
        K=atoi(Line+I);
        if (K>0)    {
          Mkr[NMkr].Seq[J]=1000+K;
          ++J;
        };
        while (isdigit(Line[I]))    {
          Mkr[NMkr].Name[JName]=Line[I];
          ++JName;
          ++I;
        };    }
      else if (Line[I]=='-')    {
        Mkr[NMkr].Name[JName]=Line[I];
        ++JName;
        ++I;
        K=atoi(Line+I);
        if (K>0)    {
          Mkr[NMkr].Seq[J]=2000+K;
          ++J;
        };
        while (isdigit(Line[I]))    {
          Mkr[NMkr].Name[JName]=Line[I];
          ++JName;
          ++I;
        };
      }
      else    {
        printf ("Unexpected character %c in pattern no. %i\n",Line[I],NMkr+1);
        ++I;
      };
    };
    if (IPar!=-1)    printf ("Error in pattern no. %i, check parentheses\n",NMkr+1);
    if (IBra!=-1)    printf ("Error in pattern no. %i, check brackets\n",NMkr+1);
    if (Mkr[NMkr].Seq[0]==1001)    printf ("Error in pattern no. %i, starts with '+1'\n",NMkr+1);
    Mkr[NMkr].Seq[J]=20000;
    Mkr[NMkr].Name[JName]='\0';
    if (Mkr[NMkr].IQStrand==0)    strcat(Mkr[NMkr].Name,"   Direct strand");
    if (Mkr[NMkr].IQStrand==1)    strcat(Mkr[NMkr].Name,"   Complementary strand");
    if (Mkr[NMkr].IQStrand==2)    strcat(Mkr[NMkr].Name,"   Both strands");
    fgets(Line,MAXLINE-1,IFi);
  };

  /* end pattern conversion  */

/*for (I=0;I<50;++I)    {printf ("%i: %i\n",I,Mkr[0].Seq[I]);};*/

  for (I=0;I<=NMkr;++I)    {
    printf ("%s\n",Mkr[I].Name);
  };

  if (OFiM!=NULL)     {
    if (IQOverlap==1)    {
      fprintf (OFiM,"Locations of patterns in the analyzed sequence\nAnalyzed sequence: %s   (%li nucleotides)\n\n",SeqFiName,LSeq);
    }
    else    {
      fprintf (OFiM,"Locations of patterns in the analyzed sequence\nAnalyzed sequence: %s   (%li nucleotides)\n",SeqFiName,LSeq);
      fprintf (OFiM,"     Combining overlaping patterns of the same type ...\n");
    };
    fprintf (OFiM,"Patterns searched:\n");
    for (I=0;I<=NMkr;++I)    {
      fprintf (OFiM,"  %3i :  %s\n",I+1,Mkr[I].Name);
    };
    fprintf (OFiM,"\n     Start       End     Pattern\n");
  };

  printf ("Searching for patterns ...\n");
  /* searching for patterns in the analyzed sequence - filling MkrStart[] and MkrEnd[] */
  NPatterns=-1;
  if (IQStr==0 || IQStr==2)    {
    printf ("Direct strand ...\n");
    for (L=1;L<=LSeq;++L)    {
      if ( L%100000==0 || L==LSeq )    printf ("%li kb ,   %li patterns found\n",L/1000,NPatterns+1);
      for (IMkr=0;IMkr<=NMkr;++IMkr)    {
        if (Mkr[IMkr].IQStrand!=1)    {
          LP=IsWord(L,Mkr[IMkr].Seq,L,Mkr[IMkr].Seq,0,&NErr,&NGoodEnds,GoodEnds,GoodEndsErr);
          if (IQOverlap==0)    {
            if (LP>Mkr[IMkr].LastEnd)    {
              if (L<=Mkr[IMkr].LastEnd)    {
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
                MkrType[NPatterns]=IMkr+1;
              };
              Mkr[IMkr].LastEnd=LP;
            };
          }
          else    {
            if (LP>0)    {
              for (I=0;I<=NGoodEnds;++I)    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrStart[NPatterns]=L;
                MkrEnd[NPatterns]=GoodEnds[I];
                fprintf (OFiM,"%10i%10i      %3i\n",L,GoodEnds[I],IMkr+1);
              };
            };
          };  
        };
      };
    };
  };
  if (IQStr==1 || IQStr==2)    {
    printf ("Complementary strand ...\n");
    for (IMkr=0;IMkr<=NMkr;++IMkr)    {
      Mkr[IMkr].LastEnd=LSeq+1;
    };
    for (L=LSeq;L>=1;--L)    {
      if ( (LSeq-L)%100000==0 || L==1 )    printf ("%li kb ,   %li patterns found\n",(LSeq-L)/1000,NPatterns+1);
      for (IMkr=0;IMkr<=NMkr;++IMkr)    {
        if (Mkr[IMkr].IQStrand!=0)    {
          LP=IsWordC(L,Mkr[IMkr].Seq,L,Mkr[IMkr].Seq,0,&NErr,&NGoodEnds,GoodEnds,GoodEndsErr);
          if (IQOverlap==0)    {
            if (LP<Mkr[IMkr].LastEnd && LP>0)    {
              if (L>=Mkr[IMkr].LastEnd)    {
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
                MkrType[NPatterns]=IMkr+1;
              };
              Mkr[IMkr].LastEnd=LP;
            };
          }
          else    {
            if (LP>0)    {
              for (I=0;I<=NGoodEnds;++I)    {
                ++NPatterns;
                if (NPatterns>=MAXPATTERNS)   {
                  printf ("Pattern(s) too general\n");
                  exit(0);
                };
                MkrEnd[NPatterns]=L;
                MkrStart[NPatterns]=GoodEnds[I];
                fprintf (OFiM,"%10i%10i      %3i\n",GoodEnds[I],L,IMkr+1);
              };
            };
          };  
        };
      };
    };
  };
  if (IQOverlap==0)    {
    for (L=0;L<=NPatterns;++L)    {
      fprintf (OFiM,"%10i%10i      %3i\n",MkrStart[L],MkrEnd[L],MkrType[L]);
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
  fprintf (OFiMS,"Pattern sequences and flanks\nAnalyzed sequence: %s   (%li nucleotides)\n\n",SeqFiName,LSeq);
  fprintf (OFiMS,"Pattern(s) :\n");
  for (I=0;I<=NMkr;++I)    {
    fprintf (OFiMS,"    %3i :  %s\n",I+1,Mkr[I].Name);
  };
  fprintf (OFiMS,"Number of patterns in the analyzed sequence :  %li\n",NPatterns+1);
  fprintf (OFiMS,"Number of non-overlaping patterns :  %li\n\n\n",NNOP+1);
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

