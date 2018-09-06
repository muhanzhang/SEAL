/* sumlines.c - total the numbers appearing in various input lines. */
/* B. D. McKay.  Version of June 18, 2016. */

#ifndef GMP
#define GMP 1  /* Non-zero if gmp multi-precise integers are allowed.
                  In this case you need the GNU multi-precision library,
                  available with -lgmp if it is installed. */
#endif

#define USAGE \
"sumlines [-w|-W] [-v] [-d] [-n] [-f fmtfile]...  file file file ..."

#define HELPTEXT \
"   Sum lines matching specified formats.\n\
\n\
   Any number of input files can be given.  \"-\" means stdin.\n\
   If there are no files given, just stdin is assumed.\n\
   File names can contain wildcards, in which case all matching files\n\
      are used in numerically sorted order.\n\
\n\
   Formats are read from four sources in this order:\n\
   (1) Any files mentioned with -f on the command line (any number).\n\
   (2) The file named in the environment variable SUMLINES.FMT (if any)\n\
   (3) The file sumlines.fmt in the current directory (if it exists)\n\
   (4) The file sumlines.fmt in the home directory (if it exists)\n\
   All these are read if they exist and the results concatenated.\n\
   Formats exactly matching earlier formats (except perhaps for flags)\n\
        are not used.\n\
\n\
   Each format occupies exactly two lines.  The first line gives a\n\
   list of flags (DEFAULT FINAL ERROR UNIQUE COUNT CONTINUE NUMERIC\n\
   SILENT ENDFILE P=# separated by spaces, commas or |s).\n\
   The second line gives the format itself.\n\
\n\
   Example.  This totals the summary lines of autoson runs:\n\
     DEFAULT  # comment \n\
     cpu=%fu,%fs,%fx  pf=%d\n\
   There can also be blank lines and lines with only comments, but\n\
   not between the flags line and the format itself.\n\
\n\
   -d don't read sumlines.fmt or ~/sumlines.fmt or $SUMLINES.FMT \n\
   -w suppresses warning messages about no matching lines or no\n\
      matching final lines.\n\
   -W in addition, suppresses woarning about missing cases.\n\
   -n don't write the number of matching lines for each format.\n\
   -v produces a list of all the formats.\n"

#define DEFAULT   0  /* No special flags */
#define FINAL     1  /* At least one of these must be in each input file */
#define ERROR     2  /* Must be none of these */
#define UNIQUE    4  /* The %s and %c parts must be unique over all inputs */
#define COUNT     8  /* The output only states how many lines matched */
#define CONTINUE 16  /* Try to match later formats too */
#define NUMERIC  32  /* Use numerical comparison (see numstrcmp() below) */
#define SILENT   64  /* Don't report, just check */
#define ENDFILE 128  /* Usually appears at end of output */

/* The formats are tried against each input line one at a time, and the
   first one that matches is accepted.  The entire line must match.
   If the CONTINUE flag is present, the input line is also matched
   against further formats.

   Except in the case of formats with the COUNT flag, each format that
   matches any lines produces output giving the total value of each of the
   integers %d or real numbers %f in the lines which match.  If there are
   any %s or %c controls in the format, the output is given separately for
   each value of the matching strings which appear in the input lines.
   %m is like %d but allows arbitrarily large integers.
   %x is like %d but the maximum rather than the sum is accumulated.
   %p is like %d but the value is accumulated modulo some base.
   %h is like %d:%d:%f taken as h:m:s with a single floating value.

   In the case of the COUNT flag, the program only reports the number of
   input lines which matched the format.

   If a format has the UNIQUE flag, no two input lines may match with the
   same values of the %s and %c controls present.  Otherwise a warning
   message is written for each duplicate match.

   The sequence P=# where # is an integer value defines the base for the
   %p directive.  There can be no spaces in the sequence "P=#".  The
   default base is 2.

   %d  - matches an integer (small enough for 64 bits)
   %x  - same as %d but accumulates maximum rather than the sum
   %p  - same as %d but accumulates the value modulo a base
   %m  - matches a integer of unbounded size (if GMP!=0)
   %f  - matches a real number of the form ddddd.ddd or -ddddd.ddd
   %v  - same as %f but reports the average rather than the sum
   %h  - similar to %d:%d:%f taken as h:m:s with a single floating value
   %sx - matches a string, where 'x' is any character. 
         If 'x' is not a space, match zero or more characters from the
             current position up but not including the first 'x'.
         If 'x' is a space, match one or more characters from the current
             position up to and including the first non-space character
             which is followed by a space.
   %c  - matches a non-white character
   %%  - matches the character '%'
   %   - (with a space following the '%') matches zero or more spaces or
            tabs, as many as appear in the input.  In the output, this
            sequence appears as one space.
   %   - (appearing exactly at the end of the format) matches zero or
         more spaces at the end of the line.  In the output, nothing.
   %*d, %*m, %*x, %*p, %*f, %*sx, %*c - these are similar to the versions
         without the '*' except that the value is ignored (not used for
         summing, and not used to divide the output).  In the output,
         this field appears as a single character '*'.
   %#  - matches an unsigned integer.  For each format containing this
         control, a report is made of any breaks or duplicates in the
         sequence of matching numbers.  (So this is useful for checking a
         sequence of case numbers.)  At most one %# may appear in each format.
   %l  - matches a list of arbitrarily many (%d sized) integers

  At least one FINAL format must match in each file or a warning is given
  (unless -w is used, in which case no warning is given).

  A format marked ENDFILE will cause sumlines to act as if it started
  reading from a new input file.  This can have some effects on the
  order of output lines.
*/

#define HAS(i,flgs)  ((format[i].flags&(flgs)) != 0)

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <pwd.h>
#include <stdlib.h>
#include <glob.h>
#include <limits.h>
#include <unistd.h>

#if GMP
#include <gmp.h>
#endif

#if defined(__alpha)
typedef long integer;
#define DOUT "%ld"       /* Formats used to output %d/%x/%p,%f,%v quantities */
#define FOUT "%.2f"
#define VOUT "%.4f"
#define HMSOUT1 "%ld:%ld:%ld"
#define HMSOUT2 "%ld:%ld:%.2f"
#elif defined(__sun) || defined(__GNUC__) || (__STDC_VERSION__ > 199900L)
typedef long long integer;
#define DOUT "%lld"   
#define FOUT "%.2f"
#define VOUT "%.4f"
#define HMSOUT1 "%lld:%lld:%lld"
#define HMSOUT2 "%lld:%lld:%.2f"
#else
typedef long long integer;
#define DOUT "%Ld"   
#define FOUT "%.2f"
#define VOUT "%.4f"
#define HMSOUT1 "%Ld:%Ld:%Ld"
#define HMSOUT2 "%Ld:%Ld:%.2f"
#endif

static char *dout,*fout,*vout,*hmsout1,*hmsout2;
static integer maxint;   /* set by find_maxint() */


#define INCR(x,inc) \
     {if (((x) > 0 && maxint-(x) < (inc)) || ((x) < 0 && (maxint)+(x) < -(inc))) \
           {fprintf(stderr,">E overflow with %%d or %%p format\n"); exit(1);} \
      x += (inc);}    /*  x += inc   with safety check */

typedef int boolean;
#define FALSE 0
#define TRUE  1

typedef struct
{
    int nvals;
    integer *val;
} integerlist;

typedef union
{
    double f;
    integer d;
    integerlist *l;
#if GMP
    mpz_t *m;
#endif
} number;

#define D 0    /* Code for "integer" */
#define F 1    /* Code for "real" */
#define M 2    /* Code for "multiprecision integer" */
#define X 3    /* Code for "integer, take maximum" */
#define V 4    /* Code for "real, take average" */
#define P 5    /* Code for "integer, modulo some base" */
#define LD 6   /* Code for "list of integer" */
#define H 8    /* Code for "h:m:s" */

#define MAXLINELEN 100000   /* Maximum input line size
                              (longer lines are broken in bits) */
#define MAXVALUES 32  /* Maximum  total number of 
                             %d,%x,%p,%m,%v,%f,%h or %l items in a format */

#define MAXFORMATS 1000

static struct fmt_st
{
   integer pmod;
   int flags;
   char *fmt;
} format[MAXFORMATS];

typedef struct countrec
{
    struct countrec *left,*right,*parent;
    char *fmt;
    unsigned long count;
    number total[MAXVALUES];
} countnode;

static countnode *count_root[MAXFORMATS];
static unsigned long matching_lines[MAXFORMATS];
static integer total_position[MAXFORMATS];
static integer lastseq[MAXFORMATS];

#if GMP
static mpz_t mp_value[MAXVALUES];
#endif

static integerlist il[MAXVALUES];

#define A 0
#define L 1
#define R 2
#define LL 3
#define LR 4
#define RL 5
#define RR 6

#ifndef GLOB_BRACE       /* Allow {} processing -- Linux extension */
#define GLOB_BRACE 0
#endif

#ifndef GLOB_TILDE      /* Allow ~  processing -- Linux extension */
#define GLOB_TILDE 0
#endif

#ifndef GLOB_NOMATCH
#define GLOB_NOMATCH 0 /* Some versions don't have a special return for this */
#endif

#define GLOB_FLAGS (GLOB_ERR|GLOB_NOSORT|GLOB_BRACE|GLOB_TILDE)

#define HELP if (argc > 1 && (strcmp(argv[1],"-help")==0 \
                           || strcmp(argv[1],"--help")==0)) \
       { printf("\nUsage: %s\n\n%s",USAGE,HELPTEXT); return 0;}

/****************************************************************************/

static int
numstrcmp(char *s1, char *s2)
/* Same behaviour as strcmp(), except that when an unsigned integer is
   found in each string, the numerical values are compared instead of
   the ascii values.   Overflow is impossible.  Leading spaces before
   numbers are considered part of the numbers.  A number in one string
   is considered less than a non-number in the other string. */
{
    char *a1,*a2;

    while (1)
    {
        for (a1 = s1; *a1 == ' '; ++a1) {}
        if (isdigit(*a1))
        {
            for (s1 = a1+1; isdigit(*s1); ++s1) {}
        }
        else 
        {
            a1 = s1;
            ++s1;
        }

        for (a2 = s2; *a2 == ' '; ++a2) {}
        if (isdigit(*a2))
        {   
            for (s2 = a2+1; isdigit(*s2); ++s2) {}
        }  
        else
        {
            a2 = s2;
            ++s2;
        }

        if (!isdigit(*a1))
        {
            if (!isdigit(*a2))
            {
                if (*a1 < *a2) return -1;
                if (*a1 > *a2) return 1;
                if (*a1 == '\0') return 0;
            }
            else
                return 1;
        }
        else
        {
            if (!isdigit(*a2))
                return -1;
            else
            {
                for (; *a1 == '0'; ++a1) {}
                for (; *a2 == '0'; ++a2) {}

                if (s1-a1 < s2-a2) return -1;
                if (s1-a1 > s2-a2) return 1;
                for (; a1 < s1 && *a1 == *a2; ++a1, ++a2) {}
                if (a1 < s1)
                {
                    if (*a1 < *a2) return -1;
                    else           return 1;
                }
            }
        }
    }
}

/****************************************************************************/

static void
writeline(char *outf, number *val, unsigned long count)
/* Write an output line with the given format and values */
{
    int i,n;
    integer mins,nsecs;
    double secs,hms;
    boolean neg;

    n = 0;

    for (; *outf != '\0'; ++outf)
    {
        if (*outf == '%')
        {
            ++outf;
            if (*outf == '%' || *outf == '#')
                putchar(*outf);
            else if (*outf == 'd' || *outf == 'x' || *outf == 'p')
                printf(dout,val[n++].d);
            else if (*outf == 'f')
                printf(fout,val[n++].f);
            else if (*outf == 'v')
                printf(vout,val[n++].f/count);
	    else if (*outf == 'h')
	    {
		if (val[n].f < 0)
		{
		    neg = TRUE;
		    hms = -val[n].f;
		}
		else
		{
		    neg = FALSE;
		    hms = val[n].f;
		}
		mins = hms/60.0;
		secs = hms - 60*mins;
		nsecs = secs;
		++n;
		if (neg) printf("-");
		if (secs == nsecs)
		    printf(hmsout1,mins/60,mins%60,nsecs);
		else
		    printf(hmsout2,mins/60,mins%60,secs);
	    }
            else if (*outf == 'l')
            {
                for (i = 0; i < val[n].l->nvals; ++i)
                {
                    if (i > 0) printf(" ");
                    printf(dout,val[n].l->val[i]);
                }
                ++n;
            }
#if GMP
            else if (*outf == 'm')
                mpz_out_str(NULL,10,*(val[n++].m));
#endif
            else
            {
                fprintf(stderr,">E unknown output format %%%c\n",*outf);
                exit(1);
            }
        }
        else
            putchar(*outf);
    }
}

/*********************************************************************/

static void
print_counts(countnode *root, boolean printcounts)
/* Use a non-recursive inorder traversal to print the tree */
{
    int code;
    countnode *p;

    p = root;
    code = A;

    while (p)
    {
        switch (code)    /* deliberate flow-ons */
        {
         case A:
            if (p->left)
            {
                p = p->left;
                break;
            }
         case L:
            if (printcounts) printf("%5lu: ",p->count);
            writeline(p->fmt,p->total,p->count);
            if (p->right)
            {
                p = p->right;
                code = A;
                break;
            }
         case R:
            if (p->parent && p->parent->left == p) code = L;
            else                                   code = R;
            p = p->parent;
            break;
        }
    }
}

/*********************************************************************/

static void
print_common(countnode *root)
/* Print the common ends of the formats in the tree */
{
    int code;
    countnode *p;
    char *s0,*s1,*t0,*t1;
    int i,comm0,comm1,minlen,maxlen;

    if (root == NULL) return;

    p = root;
    code = A;

    s0 = s1 = p->fmt;
    while (*s1 != '\0') ++s1;
    comm0 = comm1 = minlen = maxlen = s1-s0;

    while (p)
    {   
        switch (code)    /* deliberate flow-ons */
        {
         case A:
            if (p->left)
            {   
                p = p->left;
                break;
            }
         case L:
            t0 = t1 = p->fmt;
            for (i = 0; i < comm0; ++i)
                if (s0[i] != t0[i]) break;
            comm0 = i;

            while (*t1 != '\0') ++t1;
            for (i = 1; i <= comm1; ++i)
                if (s1[-i] != t1[-i]) break;
            comm1 = i-1;
            if (t1-t0 < minlen) minlen = t1-t0;
            if (t1-t0 > maxlen) maxlen = t1-t0;

            if (p->right)
            {
                p = p->right;
                code = A;
                break;
            }
         case R:
            if (p->parent && p->parent->left == p) code = L;
            else                                   code = R;
            p = p->parent;
            break;
        }
    }

    if (comm0 + comm1 > minlen) comm1 = minlen - comm0;

    for (i = 0; i < comm0; ++i)
        printf("%c",s0[i]);
    if (comm0 + comm1 < maxlen) printf("*");
    for (i = comm1; i > 0; --i)
        printf("%c",s1[-i]);
}

/*********************************************************************/

static void
splay(countnode *p)
/* Splay the node p.  It becomes the new root. */
{
    countnode *q,*r,*s;
    countnode *a,*b,*c;
    int code;

#define LCHILD(x,y) {(x)->left = y; if (y) (y)->parent = x;}
#define RCHILD(x,y) {(x)->right = y; if (y) (y)->parent = x;}

    while (p->parent)
    {
        a = p->left;
        b = p->right;
        q = p->parent;
        if (q->left == p)
        {
            code = L;
            c = q->right;
        }
        else
        {
            code = R;
            c = q->left;
        }
        r = q->parent;
        if (r)
        {
            if (r->left == q) code = (code == L ? LL : LR);
            else              code = (code == L ? RL : RR);
            s = r->parent;
            p->parent = s;
            if (s)
            {
                if (s->left == r) s->left = p;
                else              s->right = p;
            }
        }
        else
        {
            p->parent = NULL;
        }
        
        switch (code)
        {
         case L:
            RCHILD(p,q); LCHILD(q,b); break;
         case R:
            LCHILD(p,q); RCHILD(q,a); break;
         case LL:
            RCHILD(p,q); RCHILD(q,r); LCHILD(q,b); LCHILD(r,c); break;
         case RR:
            LCHILD(p,q); LCHILD(q,r); RCHILD(r,c); RCHILD(q,a); break;
         case LR:
            LCHILD(p,q); RCHILD(p,r); RCHILD(q,a); LCHILD(r,b); break;
         case RL:
            LCHILD(p,r); RCHILD(p,q); RCHILD(r,a); LCHILD(q,b); break;
        }
    }
}

/*********************************************************************/

static void
add_one(countnode **to_root, char *fmt, integer pmod, int nval,
        number *val, int *valtype, int which, boolean numcompare)
/* Add one match to the node with the given format, creating it if it is new.
   The tree is then splayed to ensure good efficiency. */
{
    int i,j,cmp;
    countnode *p,*ppar,*new_node;
    integer w;

    p = *to_root;
    cmp = 0;

    while (p != NULL)
    {
        cmp = (numcompare ? numstrcmp(fmt,p->fmt) : strcmp(fmt,p->fmt));
        if (cmp == 0)
        {
            if (HAS(which,UNIQUE) && p->count == 1)
                printf("ERROR: Multiple matches for %s",fmt);
            for (i = 0; i < nval; ++i)
                if (valtype[i] == D)
                    {INCR(p->total[i].d,val[i].d);}
                else if (valtype[i] == X)
                    {if (val[i].d > p->total[i].d) p->total[i].d = val[i].d;}
                else if (valtype[i] == P)
                    {w = val[i].d % pmod; INCR(p->total[i].d,w);
                      p->total[i].d %= pmod;}
                else if (valtype[i] == LD)
                {
                    if (p->total[i].l->nvals < val[i].l->nvals)
                    {
                        if ((p->total[i].l->val
                           = (integer*)realloc(p->total[i].l->val,
                                          sizeof(integer)*val[i].l->nvals))
                                == NULL)
                        {
                            fprintf(stderr,"Malloc failed\n");
                            exit(1);
                        }
                    }
                    for (j = 0; j < p->total[i].l->nvals &&
                                j < val[i].l->nvals; ++j)
                        INCR(p->total[i].l->val[j],val[i].l->val[j]);
                    if (p->total[i].l->nvals < val[i].l->nvals)
                    {
                        for (j = p->total[i].l->nvals;
                                   j < val[i].l->nvals; ++j)
                            p->total[i].l->val[j] = val[i].l->val[j];
                        p->total[i].l->nvals = val[i].l->nvals;
                    }
                }
#if GMP
                else if (valtype[i] == M) 
                    mpz_add(*(p->total[i].m),*(p->total[i].m),*(val[i].m));
#endif
                else  
                    p->total[i].f += val[i].f;   /* F, V and H */
            ++p->count;
            splay(p);
            *to_root = p;
            return;
        }
        else if (cmp < 0)
        {
            ppar = p;
            p = p->left;
        }
        else
        {
            ppar = p;
            p = p->right;
        }
    }

    if ((new_node = (countnode*)malloc(sizeof(countnode))) == NULL)
    {
        fprintf(stderr,">E malloc failed in add_one()\n");
        exit(1);
    }

    if ((new_node->fmt = (char*)malloc(strlen(fmt)+1)) == NULL)
    {
        fprintf(stderr,">E malloc failed in add_one()\n");
        exit(1);
    }

    new_node->count = 1;
    strcpy(new_node->fmt,fmt);
    for (i = 0; i < nval; ++i)
    {
#if GMP
        if (valtype[i] == M)
        {
            if ((new_node->total[i].m
                                = (mpz_t*)malloc(sizeof(mpz_t))) == NULL)
            {
                fprintf(stderr,"Malloc failed\n");
                exit(1);
            }
            mpz_init_set(*(new_node->total[i].m),*(val[i].m));
        }
        else
#endif
        if (valtype[i] == LD)
        {
            if ((new_node->total[i].l
                        = (integerlist*)malloc(sizeof(integerlist))) == NULL)
            {  
                fprintf(stderr,"Malloc failed\n");
                exit(1);
            }
            if ((new_node->total[i].l->val
                 = (integer*)malloc(sizeof(integer)*val[i].l->nvals)) == NULL)
            {  
                fprintf(stderr,"Malloc failed\n");
                exit(1);
            }
            new_node->total[i].l->nvals = val[i].l->nvals;
            for (j = 0; j < val[i].l->nvals; ++j)
                new_node->total[i].l->val[j] = val[i].l->val[j];
        }
        else
            new_node->total[i] = val[i];
    }

    new_node->left = new_node->right = NULL;

    if (cmp == 0)
    {
        *to_root = new_node;
        new_node->parent = NULL;
    }
    else if (cmp < 0)
    {
        ppar->left = new_node;
        new_node->parent = ppar;
    }
    else
    {
        ppar->right = new_node;
        new_node->parent = ppar;
    }

    splay(new_node);
    *to_root = new_node;
}

/****************************************************************************/

static int
scanline(char *s, char *f, number *val, int *valtype,
         integer *seqno, char *outf)
/* Perform sscanf-like scan of line.
   The whole format must match.  outf is set to be an output format
   with unassigned values replaced by '*' and %s replaced by what
   it matches.  Assigned values except %s are put into val[] with
   their types in valtype[].  The number of values (not counting %#)
   is returned. 
   Integers matching %# are put into *seqno, with an error if there
   are more than one, and -1 if there are none.
   If the format doesn't match, -1 is returned.
   WARNING: the gmp and ilist values are pointers to static data,
   so they need to be copied if the values array is copied.
   See the comments at the start of the program for more information.
*/
{
    int n;                   /* Number of values assigned */
    int digit;
    boolean doass,neg,oflow,badgmp;
    integer ival;
    double dval,digval,comval;
    char ends,*saves;
    static boolean gmp_warning = FALSE;
    integer *ilist;
    size_t ilist_sz;
    int nilist;
#if GMP
    char mp_line[MAXLINELEN+1],*mp;
#endif
        
    n = 0;
    *seqno = -1;
    badgmp = oflow = FALSE;

    while (*f != '\0')
    {
        if (*f == '%')
        {
            ++f;
            if (*f == '*')
            {
                doass = FALSE;
                ++f;
            }
            else
                doass = TRUE;

            if (*f == '%')
            {
                if (!doass)
                {
                    fprintf(stderr,"Bad format item %%*\n");
                    exit(1);
                }
                if (*s++ != '%') return -1;
                ++f;
                *outf++ = '%';
                *outf++ = '%';
            }
            else if (*f == '\n')
            {
                if (!doass)
                {
                    fprintf(stderr,"Bad format item %%*\n");
                    exit(1);
                }
                while (*s != '\0')
                {
                    if (*s != ' ' && *s != '\n') return -1;
                    ++s;
                }
                --s;
            }
            else if (*f == 'c')
            {   
                if (*s == ' ' || *s == '\t' || *s == '\n') return -1;
                if (doass) *outf++ = *s;
                else       *outf++ = '*';
                ++f;
                ++s;
            }
            else if (*f == 's')
            {
                ends = *(f+1);
                if (ends == ' ')
                {
                    while (*s == ' ' || *s == '\t')
                    {
                        if (doass) *outf++ = *s;
                        ++s;
                    }
                }
                while (*s != '\n' && *s != ends)
                {
                    if (doass) *outf++ = *s;
                    ++s;
                }
                if (!doass) *outf++ = '*';
                ++f;
            }
#if GMP         
            else if (*f == 'd' || *f == 'x' || *f == 'p')
            {
#else
            else if (*f == 'd' || *f == 'x' || *f == 'p' || *f == 'm')
            {
                if (*f == 'm' && !gmp_warning)
                {
                    fprintf(stderr,
                     ">W not compiled with GMP, treating %%m like %%d\n");
                    gmp_warning = TRUE;
                }
#endif
                while (*s == ' ' || *s == '\t') ++s;
                if (!isdigit(*s) && *s != '-' && *s != '+') return -1;
                neg = (*s == '-');
                if (*s == '-' || *s == '+') ++s;
                ival = 0;
                while (isdigit(*s))
                {
                    digit =  *s++ - '0';
                    if (ival > (maxint-digit)/10)
                        oflow = TRUE;
                    else
                        ival = ival*10 + digit;
                }
                if (doass)
                {
                    *outf++ = '%';
                    if (*f == 'd' || *f == 'm')
                    {
                        *outf++ = 'd';
                        valtype[n] = D;
                    }
                    else if (*f == 'x')
                    {
                        *outf++ = 'x';
                        valtype[n] = X;
                    }
                    else
                    {
                        *outf++ = 'p';
                        valtype[n] = P;
                    }
                    val[n++].d = (neg ? -ival : ival);
                }
                else
                    *outf++ = '*';
                ++f;
            }
            else if (*f == 'l')
            {
                nilist = 0;
                if ((ilist = (integer*)malloc(200*sizeof(integer)))
                            == NULL)
                {
                    fprintf(stderr,"Malloc failed\n");
                    exit(1);
                }
                ilist_sz = 200;
                for (;;)
                {
                    saves = s;
                    while (*s == ' ' || *s == '\t') ++s;
                    if (!isdigit(*s) && *s != '-' && *s != '+')
                    {
                        s = saves;
                        break;
                    }
                    neg = (*s == '-');
                    if (*s == '-' || *s == '+') ++s;
                    ival = 0;
                    while (isdigit(*s))
                    {
                        digit =  *s++ - '0';
                        if (ival > (maxint-digit)/10)
                            oflow = TRUE;
                        else
                            ival = ival*10 + digit;
                    }
                    if (neg) ival = -ival;
                    if (nilist == ilist_sz)
                    {
                        if ((ilist
                            = (integer*)realloc((void*)ilist,
                                       (ilist_sz+500)*sizeof(integer)))
                                == NULL)
                        {
                            fprintf(stderr,"Malloc failed\n");
                            exit(1);
                        }
                        ilist_sz += 500;
                    }
                    ilist[nilist++] = ival;
                }
                if (doass)
                {
                    valtype[n] = LD;
                    val[n].l = &il[n];
                    val[n].l->nvals = nilist;
                    if (val[n].l->val) free(val[n].l->val);
                    val[n].l->val = ilist;
                    ++n;
                    *outf++ = '%';
                    *outf++ = 'l';
                }
                else
                {
                    free(ilist);
                    *outf++ = '*';
                }
                ++f;
            }
#if GMP
            else if (*f == 'm')
            {
                while (*s == ' ' || *s == '\t') ++s;        
                if (!isdigit(*s) && *s != '-' && *s != '+') return -1;
                mp = mp_line;
                if      (*s == '-') *mp++ = *s++;
                else if (*s == '+') s++;
                while (isdigit(*s)) *mp++ = *s++;
                *mp = '\0'; 
                if (doass)
                {
                    valtype[n] = M;
                    val[n].m = &mp_value[n];
                    if (mpz_set_str(mp_value[n],mp_line,10) < 0)
                        badgmp = TRUE;
                    ++n;
                    *outf++ = '%';
                    *outf++ = 'm';
                }
                else
                    *outf++ = '*';
                ++f;
            }
#endif
            else if (*f == '#')
            {
                while (*s == ' ' || *s == '\t') ++s;
                if (!isdigit(*s)) return -1;
                ival = 0;
                while (isdigit(*s))
                {
                    digit =  *s++ - '0';
                    if (ival > (maxint-digit)/10)
                        oflow = TRUE;
                    else
                        ival = ival*10 + digit;
                }
                if (*seqno >= 0)
                {
                    fprintf(stderr, 
                            ">E %%# can only be used once per format\n");
                    exit(1);
                }
                *seqno = ival;
                *outf++ = '#';
                ++f;
            }
            else if (*f == 'f' || *f == 'v')
            {
                while (*s == ' ' || *s == '\t') ++s;

                if (!isdigit(*s) && *s != '.' && *s != '-' && *s != '+')
                    return -1;
                neg = (*s == '-');
                if (*s == '-' || *s == '+') ++s;
                dval = 0.0;
                while (isdigit(*s)) dval = dval*10.0 + (*s++ - '0');
                if (*s == '.')
                {
                    digval = 1.0;
                    ++s;
                    while (isdigit(*s))
                    {
                        digval /= 10.0;
                        dval += (*s++ - '0') * digval;
                    }
                }
                if (doass)
                {
                    valtype[n] = (*f == 'f' ? F : V);
                    val[n++].f = (neg ? -dval : dval);
                    *outf++ = '%';
                    *outf++ = *f;
                }
                else
                    *outf++ = '*';
                ++f;
            }
            else if (*f == 'h')
            {
                while (*s == ' ' || *s == '\t') ++s;

                if (!isdigit(*s) && *s != '.' && *s != '-' && *s != '+' && *s != ':')
                    return -1;
                neg = (*s == '-');
                if (*s == '-' || *s == '+') ++s;
                dval = 0.0;
		comval = 0.0;
                while (isdigit(*s) || *s == ':')
		{
		    if (*s == ':')
		    {
			dval = dval*60.0 + comval;
			comval = 0.0;
			++s;
		    }
		    else
		        comval = comval*10.0 + (*s++ - '0');
		}
                if (*s == '.')
                {
                    digval = 1.0;
                    ++s;
                    while (isdigit(*s))
                    {
                        digval /= 10.0;
                        comval += (*s++ - '0') * digval;
                    }
                }
		dval = dval*60.0 + comval;
                if (doass)
                {
                    valtype[n] = H;
                    val[n++].f = (neg ? -dval : dval);
                    *outf++ = '%';
                    *outf++ = *f;
                }
                else
                    *outf++ = '*';
                ++f;
            }
            else if (*f == ' ')
            {
                while (*s == ' ' || *s == '\t') ++s;
                *outf++ = ' ';
                ++f;
            }
            else
            {
                fprintf(stderr,"Bad format item %%%c\n",*f);
                exit(1);
            }
        }
        else
        {
            if (*s != *f) return -1;
            *outf++ = *f;
            ++s;
            ++f;
        }
    }

    if (*s != '\0') return -1;

    *outf = '\0';

    if (oflow)
    {
        fprintf(stderr,"Integer too large\n");
        exit(1);
    }
    if (badgmp)
    {
        fprintf(stderr,"Illegal multiprecision integer\n");
        exit(1);
    }

    return n;
}

/****************************************************************************/

void
find_maxint(void)
{
/* Put the maximum possible integer value into maxint. */
/* New version with no integer overflow. */
    integer x,y;

    x = ((integer)1) << (8*sizeof(integer) - 2);
    y = x - 1;
    x += y;

    if (x <= 0)
    {
        fprintf(stderr,">E find_maxint() failed\n");
        exit(1);
    }

    maxint = x;
}

/****************************************************************************/

static void
sort_formats(int *order, int numformats)
/* Make order[0..numformats-1] a permutation of 0..numformats-1 being
   a good order to display the results. */
{
    double score[MAXFORMATS];
    int h,i,j,iw;

    for (i = 0; i < numformats; ++i)
    {
        if (matching_lines[i] == 0)
            score[i] = -1.0;
        else
            score[i] = i +
               ((100.0*total_position[i]) / matching_lines[i]) * numformats;
        order[i] = i;
    }

    j = numformats / 3;
    h = 1;
    do
        h = 3 * h + 1;
    while (h < j);

    do
    {
        for (i = h; i < numformats; ++i)
        {
            iw = order[i];
            for (j = i; score[order[j-h]] > score[iw]; )
            {
                order[j] = order[j-h];
                if ((j -= h) < h) break;
            }
            order[j] = iw;
        }
        h /= 3;
    }
    while (h > 0);
}

/****************************************************************************/

static void
read_formats(char *filename, int *numformatsp, boolean mustexist)
/* Read formats from the given file. */
{
    FILE *f;
    int i,c,flags,ignore;  
    char flagname[52];
    char line[MAXLINELEN+3];
    integer pmod;
    char *s;
    boolean oflow,badpmod;
    int digit;

    if (strcmp(filename,"-") == 0)
        f = stdin;
    else if ((f = fopen(filename,"r")) == NULL)
    {
        if (mustexist)
        {
            fprintf(stderr,">E Can't open %s for reading.\n",filename);
            exit(1);
        }
        return;
    }

    line[MAXLINELEN+2] = '\0';

    for (;;)
    {
        if ((c = getc(f)) == EOF) break;

        while (c == ' ' || c == '\t') c = getc(f);
        if (c == '\n') continue;
        if (c == EOF) break;
        
        if (c == '#')
        {
            while (c != '\n' && c != EOF) c = getc(f);
            continue;
        }

        ungetc(c,f);

        flags = 0;
        pmod = 2;
        for (;;)
        {
            while ((c = getc(f)) == ' '
                                || c == '|' || c == ',' || c == '\t') {}
            if (c == '#')
                while (c != '\n' && c != EOF) c = getc(f);
            if (c == '\n' || c == EOF) break;

            ungetc(c,f);

         /* There appear to be some issues with the [ flag in fscanf,
          * as to whether a null is appended.  We'll take no chances. */
            for (i = 0; i < 52; ++i) flagname[i] = '\0';
            ignore = fscanf(f,"%50[A-Za-z0-9=]",flagname);

            if      (strcmp(flagname,"DEFAULT") == 0)  {}
            else if (strcmp(flagname,"FINAL") == 0)    flags |= FINAL;
            else if (strcmp(flagname,"ERROR") == 0)    flags |= ERROR;
            else if (strcmp(flagname,"UNIQUE") == 0)   flags |= UNIQUE;
            else if (strcmp(flagname,"COUNT") == 0)    flags |= COUNT;
            else if (strcmp(flagname,"CONTINUE") == 0) flags |= CONTINUE;
            else if (strcmp(flagname,"NUMERIC") == 0)  flags |= NUMERIC;
            else if (strcmp(flagname,"SILENT") == 0)   flags |= SILENT;
            else if (strcmp(flagname,"ENDFILE") == 0)  flags |= ENDFILE;
            else if (flagname[0] == 'P' && flagname[1] == '=')
            {
                pmod = 0;
                oflow = FALSE;
                badpmod = (flagname[2] == '\0');
                for (s = flagname+2; *s != '\0'; ++s)
                {
                    if (isdigit(*s))
                    {
                        digit =  *s - '0';
                        if (pmod > (maxint-digit)/10)
                            oflow = TRUE;
                        else
                            pmod = pmod*10 + digit;
                    }
                    else
                        badpmod = TRUE;
                }
                if (badpmod)
                {
                    fprintf(stderr,">E Bad value for P= directive: %s\n",
                            flagname+2);
                    exit(1);
                }
                else if (oflow)
                {
                    fprintf(stderr,">E Value for P= is too large\n");
                    exit(1);
                }
    
            }
            else
            {
                fprintf(stderr,">E Unknown flag \"%s\" in %s\n",
                               flagname,filename);
                exit(1);
            }
        }

        if (fgets(line,MAXLINELEN,f) == NULL)
        {
            fprintf(stderr,">E Missing format in %s\n",filename);
            exit(1);
        }

        for (i = 0; i < *numformatsp; ++i)
            if (strcmp(line,format[i].fmt) == 0) break;
        if (i < *numformatsp) continue;

        if (*numformatsp == MAXFORMATS)
        {
            fprintf(stderr,">E Increase MAXFORMATS\n");
            exit(1);
        }

        format[*numformatsp].flags = flags;
        format[*numformatsp].pmod = pmod;
        if ((format[*numformatsp].fmt
                            = (char*)malloc(strlen(line)+1)) == NULL)
        {
            fprintf(stderr,">E malloc() failed in read_formats()\n");
            exit(1);
        }
        strcpy(format[*numformatsp].fmt,line);
        ++*numformatsp;
    }

    if (f != stdin) fclose(f);
}

/****************************************************************************/

static void
read_local_formats(int *numformatsp)
/* Read formats from sumlines.fmt in current directory */
{
        read_formats("sumlines.fmt",numformatsp,FALSE);
}

/****************************************************************************/

static void
read_global_formats(int *numformatsp)
/* Read formats from sumlines.fmt in home directory */
{
    struct passwd *pwd;
    char *homedir;
    char filename[4097];

    homedir = getenv("HOME");
    if (homedir == NULL && (pwd = getpwuid(getuid())) != NULL)
        homedir = pwd->pw_dir;

    if (homedir == NULL)
    {
        fprintf(stderr,">W Can't find home directory\n");
        return;
    }

    sprintf(filename,"%s/sumlines.fmt",homedir);
    read_formats(filename,numformatsp,FALSE);
}

/****************************************************************************/

static void
read_env_formats(int *numformatsp)
/* Read formats from $SUMLINES.FMT if it exists */
{
    char *filename;

    if ((filename = getenv("SUMLINES.FMT")) != 0)
        read_formats(filename,numformatsp,FALSE);
}

/****************************************************************************/

static boolean
readoneline(FILE *f, char *line, int size, int *nulls)
/* Get a line.  Read at most size-1 chars until EOF or \n.
   If \n is read, it is stored.  Then \0 is appended.
   *nulls is set to the number of NUL chars (which are also stored). */
{
    int i,c;

    *nulls = 0;
    for (i = 0; i < size-1; ++i)
    {
        c = getc(f);
        if (c == EOF) break;
        line[i] = c;
        if (c == '\0') ++*nulls;
        if (c == '\n') {++i; break;}
    }
    line[i] = '\0';

    return i > 0;
}

/****************************************************************************/

static int
pnumstrcmp(const void *a, const void *b)
/* numstrcmp on strings pointed at by a and b */
{
    return numstrcmp(*(char**)a,*(char**)b);
}

/****************************************************************************/

static void
doglob(char *patt, glob_t *globlk)
/* Find all files matching the given pattern, numeric sorting.
   Give a warning message if there are none. */
{
    int ret;

    ret = glob(patt,GLOB_FLAGS,NULL,globlk);

    if (ret != 0) globlk->gl_pathc = 0;

    if (ret == GLOB_NOSPACE)
    {
        fprintf(stderr,"ERROR: ran out of space during glob()\n");
        exit(1);
    }
    if (ret == GLOB_ERR)
    {
        fprintf(stderr,"ERROR: during glob(%s)\n",patt);
        exit(1);
    }
    if (ret != 0 && ret != GLOB_NOMATCH)
    {
        fprintf(stderr,"ERROR: value %d from glob(%s)\n",ret,patt);
        exit(1);
    }


    if (globlk->gl_pathc == 0) printf("WARNING: no files match %s\n",patt);
    
    if (globlk->gl_pathc >= 2)
        qsort(globlk->gl_pathv,globlk->gl_pathc,sizeof(char*),pnumstrcmp);
}

/****************************************************************************/

int
main(int argc, char *argv[])
{
    int i,j,nvals,argnum;
    number val[MAXVALUES];
    int valtype[MAXVALUES];
    char line[MAXLINELEN+2];
    char outf[MAXLINELEN+MAXVALUES+6];
    unsigned long matched,unmatched,finalmatched;
    unsigned long errorlines,totalerrorlines;
    unsigned long line_number,nullcount,numfiles,ifile;
    char *filename;
    FILE *infile;
    int numformats,firstarg,nulls;
    boolean havefinal,nowarn,noWarn,listformats,readfiles;
    integer seq;
    int order[MAXFORMATS];
    glob_t globlk,globlk_stdin,*pglob;
    char *glob_stdin_v[2];
    boolean printcounts;

    HELP;

    find_maxint();

    firstarg = 1;
    numformats = 0;
    nowarn = noWarn = FALSE;
    listformats = FALSE;
    readfiles = TRUE;
    printcounts = TRUE;

    globlk_stdin.gl_pathc = 1;
    globlk_stdin.gl_pathv = glob_stdin_v;
    glob_stdin_v[0] = "-";
    glob_stdin_v[1] = NULL;

    dout = DOUT;
    fout = FOUT;
    vout = VOUT;
    hmsout1 = HMSOUT1;
    hmsout2 = HMSOUT2;

    for (; firstarg < argc; ++firstarg)
    {
        if (argv[firstarg][0] == '-' && argv[firstarg][1] == 'f')
        {
            if (argv[firstarg][2] != '\0')
                read_formats(&argv[firstarg][2],&numformats,TRUE);
            else if (firstarg == argc - 1)
            {
                fprintf(stderr,">E No argument for -f\n");
                exit(1);
            }
            else
            {
                ++firstarg;
                read_formats(argv[firstarg],&numformats,TRUE);
            }
        }
        else if (strcmp(argv[firstarg],"-W") == 0)
            noWarn = TRUE;
        else if (strcmp(argv[firstarg],"-w") == 0)
            nowarn = TRUE;
        else if (strcmp(argv[firstarg],"-v") == 0)
            listformats = TRUE;
        else if (strcmp(argv[firstarg],"-d") == 0)
            readfiles = FALSE;
        else if (strcmp(argv[firstarg],"-n") == 0)
            printcounts = FALSE;
        else if (strcmp(argv[firstarg],"-V") == 0)
            vout = argv[++firstarg];
        else if (strcmp(argv[firstarg],"-F") == 0)
            fout = argv[++firstarg];
        else if (strcmp(argv[firstarg],"-D") == 0)
            dout = argv[++firstarg];
        else
            break;
    }

#if GMP
    for (i = 0; i < MAXVALUES; ++i) mpz_init(mp_value[i]);
#endif
    for (i = 0; i < MAXVALUES; ++i) 
    {
        il[i].nvals = 0;
        il[i].val = NULL;
    }

    if (noWarn) nowarn = TRUE;

    if (readfiles) read_local_formats(&numformats);
    if (readfiles) read_env_formats(&numformats);
    if (readfiles) read_global_formats(&numformats);

    if (listformats)
    {
        printf("%d formats:\n",numformats);
        for (i = 0; i < numformats; ++i)
            printf("%03x %s",format[i].flags,format[i].fmt);
    }

    if (numformats == 0)
    {
        fprintf(stderr,">E No formats\n");
        exit(1);
    }

    havefinal = FALSE;
    for (i = 0; i < numformats; ++i)
    {
        count_root[i] = NULL;
        matching_lines[i] = 0;
        total_position[i] = 0;
        if (HAS(i,FINAL)) havefinal = TRUE;
    }

    unmatched = totalerrorlines = 0;
    numfiles = 0;

    for (argnum = firstarg;
         argnum < (argc == firstarg ? argc+1 : argc); ++argnum)
    {
        if (argnum >= argc || strcmp(argv[argnum],"-") == 0)
            pglob = &globlk_stdin;
        else
        {
            pglob = &globlk;
            doglob(argv[argnum],pglob);
        }

        for (ifile = 0; ifile < pglob->gl_pathc; ++ifile)
        {
            matched = finalmatched = errorlines = 0;
            ++numfiles;

            if (strcmp(pglob->gl_pathv[ifile],"-") == 0)
            {
                filename = "stdin";
                infile = stdin;
            }
            else
            {
                filename = pglob->gl_pathv[ifile];
                if ((infile = fopen(filename,"r")) == NULL)
                {
                    fprintf(stderr,">E Can't open %s\n",filename);
                    exit(1);
                }
            }

            line_number = 0;
            nullcount = 0;
            while (readoneline(infile,line,MAXLINELEN,&nulls))
            {
                nullcount += nulls;
                line[MAXLINELEN] = '\n';
                line[MAXLINELEN+1] = '\0';
                if (line[0] == '\n') continue;
                ++line_number;
    
                for (i = 0; i < numformats; ++i)
                {
                    nvals
                      = scanline(line,format[i].fmt,val,valtype,&seq,outf);
                    if (nvals >= 0)
                    {
                        if (HAS(i,ENDFILE)) line_number = 0;
                        ++matched;
                        if (HAS(i,FINAL)) ++finalmatched;
                        if (HAS(i,ERROR)) ++errorlines;
                        ++matching_lines[i];
                        total_position[i] += line_number;
                        add_one(&count_root[i],outf,format[i].pmod,nvals,
                             val,valtype,i,HAS(i,NUMERIC));
                        if (!noWarn && matching_lines[i] > 1 && seq >= 0 
                                               && seq != lastseq[i]+1)
                        {
                            printf("WARNING: Sequence number");
                            if (seq == lastseq[i])
                            {   
                                printf(" ");
                                printf(dout,seq);
                                printf(" is repeated.\n");
                            }
                            else if (seq != lastseq[i]+2)
                            {
                                printf("s ");
                                printf(dout,lastseq[i]+1);
                                printf("-");
                                printf(dout,seq-1);
                                printf(" are missing.\n");
                            }
                            else
                            {
                                printf("  ");
                                printf(dout,seq-1);
                                printf(" is missing.\n");
                            }
                        }
                        lastseq[i] = seq;
                        if (!HAS(i,CONTINUE)) break;
                    }
                }
    
                if (i == numformats) ++unmatched;
            }
            if (errorlines != 0)
                printf("ERRORS: Error lines in file %s\n",filename);
            else if (matched == 0 && !nowarn)
                printf("WARNING: No matching lines in file %s\n",filename);
            else if (finalmatched == 0 && havefinal && !nowarn)
                printf("WARNING: No final lines in file %s\n",filename);
            if (nullcount > 0)
                printf("WARNING: %ld NULs found in file %s\n",
                                                        nullcount,filename);
            if (infile != stdin) fclose(infile);

            totalerrorlines += errorlines;
        }
        if (pglob == &globlk) globfree(pglob);
    }

    sort_formats(order,numformats);

    for (j = 0; j < numformats; ++j)
    {
        i = order[j];
        if (HAS(i,SILENT)) continue;

        if (HAS(i,COUNT))
        {
            if (matching_lines[i] > 0)
                printf("%5lu lines matched ",matching_lines[i]);
            print_common(count_root[i]);
        }
        else
            print_counts(count_root[i],printcounts);
    }

    if (unmatched > 0)
        printf("%5lu non-empty lines not matched\n",unmatched);
    if (argc > firstarg) printf("%5lu files read altogether\n",numfiles);
    if (totalerrorlines > 0) printf("%5lu errors found\n",totalerrorlines);

    exit(0);
}
