/* rng.c
   This file contains the code for a high-quality random number
   generator written by Don Knuth.  The auxilliary routine
   ran_arr_cycle() has been modified slightly, and ran_init() is new.

   To use it:

      0.  #include "rng.h" (or "naututil.h" if you are using nauty)

      1.  Call ran_init(seed), where seed is any long integer.
          This step is optional, but if you don't use it you
          will always get the same sequence of random numbers.

      2.  For each random number, use the NEXTRAN macro.  It will
          give a random value in the range 0..2^30-1.  Alternatively,
          KRAN(k) will have a random value in the range 0..k-1.
          KRAN(k) actually gives you NEXTRAN mod k, so it is not
          totally uniform if k is very large.  In that case, you
          can use the slightly slower GETKRAN(k,var) to set the
          variable var to a better random number from 0..k-1.

    Brendan McKay, July 2002.  Fixed Nov 2002 on advice of DEK.

*/

/*    This program by D E Knuth is in the public domain and freely copyable
 *    AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to Volume 2 on pages 171 and following).              */

/*    N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
      included here; there's no backwards compatibility with the original. */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define MM (1L<<30)                 /* the modulus */
#define mod_diff(x,y) (((x)-(y))&(MM-1)) /* subtraction mod MM */

long ran_x[KK];                    /* the generator state */

#ifdef __STDC__
void ran_array(long aa[],int n)
#else
void ran_array(aa,n)    /* put n new random numbers in aa */
  long *aa;   /* destination */
  int n;      /* array length (must be at least KK) */
#endif
{
  int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_x[j];
  for (;j<n;j++) aa[j]=mod_diff(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_x[i]=mod_diff(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_x[i]=mod_diff(aa[j-KK],ran_x[i-LL]);
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
static long ran_arr_buf[QUALITY];
static long ran_arr_dummy=-1, ran_arr_started=-1;
long *ran_arr_ptr=&ran_arr_dummy; /* the next random number, or -1 */

#define TT  70   /* guaranteed separation between streams */
#define is_odd(x)  ((x)&1)          /* units bit of x */

#ifdef __STDC__
void ran_start(long seed)
#else
void ran_start(seed)    /* do this before using ran_array */
  long seed;            /* selector for different streams */
#endif
{
  int t,j;
  long x[KK+KK-1];              /* the preparation buffer */
  long ss=(seed+2)&(MM-2);
  for (j=0;j<KK;j++) {
    x[j]=ss;                      /* bootstrap the buffer */
    ss<<=1; if (ss>=MM) ss-=MM-2; /* cyclic shift 29 bits */
  }
  x[1]++;              /* make x[1] (and only x[1]) odd */
  for (ss=seed&(MM-1),t=TT-1; t; ) {       
    for (j=KK-1;j>0;j--) x[j+j]=x[j], x[j+j-1]=0; /* "square" */
    for (j=KK+KK-2;j>=KK;j--)
      x[j-(KK-LL)]=mod_diff(x[j-(KK-LL)],x[j]),
      x[j-KK]=mod_diff(x[j-KK],x[j]);
    if (is_odd(ss)) {              /* "multiply by z" */
      for (j=KK;j>0;j--)  x[j]=x[j-1];
      x[0]=x[KK];            /* shift the buffer cyclically */
      x[LL]=mod_diff(x[LL],x[KK]);
    }
    if (ss) ss>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_x[j+KK-LL]=x[j];
  for (;j<KK;j++) ran_x[j-LL]=x[j];
  for (j=0;j<10;j++) ran_array(x,KK+KK-1); /* warm things up */
  ran_arr_ptr=&ran_arr_started;
}

void
ran_init(long seed)    /* Added by BDM: use instead of ran_start. */
                       /*  But this is less important with this version */
{
    ran_start((unsigned long)seed % (MM-2));
}

#define ran_arr_next() (*ran_arr_ptr>=0? *ran_arr_ptr++: ran_arr_cycle())

long
ran_arr_cycle(void)
/* Modified by BDM to automatically initialise 
   if no explicit initialisation has been done */
{
   if (ran_arr_ptr==&ran_arr_dummy)
       ran_start(314159L); /* the user forgot to initialize */

  ran_array(ran_arr_buf,QUALITY);

  ran_arr_buf[KK]=-1;
  ran_arr_ptr=ran_arr_buf+1;
  return ran_arr_buf[0];
}
