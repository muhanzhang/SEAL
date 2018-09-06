/* gtools.c : Common routines for gtools programs. */
/* Version 4.1, Jan 2016. */

/* Todo: size check if MAXN>0; option to free memory */

#include "gtools.h"

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

TLS_ATTR size_t ogf_linelen;
TLS_ATTR boolean is_pipe;

#if HAVE_FSEEKO
#define FSEEK_VER fseeko
#define FTELL_VER ftello
#define OFF_T_VER off_t
#else
#if !FTELL_DEC
extern long ftell(FILE*);
extern int fseek(FILE*,long,int);
#endif
#define FSEEK_VER fseek
#define FTELL_VER ftell
#define OFF_T_VER long
#endif

#if !POPEN_DEC
extern FILE *popen(const char*,const char*);
#endif

/*
  Version 1.1: Fixed sparse6 input for powers of 2.  May 9, 1998
  Version 1.2: Added "cmd: ..." option for opengraphfile(). 
               Fixed readg() bug (could not be invisible).  Oct 5, 1998
  Version 1.3: Added "is_pipe".  June 20, 2002
  Version 1.4: Stuff for autoconf.  August 30, 2002
  Version 1.5: Unlocked stdio for efficiency.  October 31, 2004
               Also fwrite() in place of fputs() for writeline().
  Version 1.6: i/o for sparsegraph; use of s6len; improve allocations
  Version 1.7: Add stringcounts()
  	       Add very long size code (see formats.txt)
  Version 1.8: Add gtools_check()
  Version 1.9: Add writepc_sg(), readpc_sg() and readpcle_sg()
               Add planar_code options to opengraphfile()
  Version 2.4: Add writeec_sg(), readec_sg()  (MISSING!)
               Add edge_code options to opengraphfile()
  Version 2.5: Remove sortints(), not used
  Version 2.6: Add sgtog6() and writeg6_sg()
  Version 2.7: Add lots of explicit casts.
               Fix planar code output for n > 255.
  Version 3.0: Procedures for incremental sparse6 format.
               Add checkgline()
  Version 4.0: Procedures for digraph6 format.
  Version 4.1: Made encodegraphsize() external.
*/

#define B(i) (1 << ((i)-1))
#define M(i) ((1 << (i))-1)

/*********************************************************************
opengraphfile(filename,codetype,assumefixed,position) 
          opens and positions a file for reading graphs.

  filename = the name of the file to open 
                (NULL means stdin, assumed already open)
             If filename starts with "cmd:", the remainder is taken
             to be a command to open a subshell for, using a pipe.
  codetype   = returns a code for the format.
                This is a combination of SPARSE6, GRAPH6, PLANARCODE,
                PLANARCODELE, PLANARCODEBE, EDGECODE, DIGRAPH6,
                UNKNOWN_TYPE and HAS_HEADER.
                If a header is present, that overrides the data.
                If there is no header, the first graph is examined.
  assumefixed = nonzero if files other than stdin or pipes should be
                assumed to be seekable and have equal record sizes.
                Ignored if there is a sparse6 header or the first
                graph has sparse6 format.
  position = the number of the record to position to
                (the first is number 1; 0 and -NOLIMIT also mean
                to position at start).  planar_code files can only
                be positioned at the start.

  If the file starts with ">", there must be a header.
  Otherwise opengraphfile() fails.

  The value returned is a file pointer or NULL.  
  If assumedfixed is not zero and position > 1, the global variable
  ogf_linelen is set to the length (including \n) of the length of the 
  first record.  UPDATE

  The global variable is_pipe is set to whether the input file is a pipe.

**********************************************************************/

FILE*
opengraphfile(char *filename, int *codetype, int assumefixed, long position)
{
    FILE *f;
    int c,bl,firstc;
    long i,l;
    OFF_T_VER pos,pos1,pos2;
    boolean bad_header;

    is_pipe = FALSE;

    if (filename == NULL)
    {
        f = stdin;
        assumefixed = FALSE;
    }
    else
    {
        if (filename[0] == 'c' && filename[1] == 'm'
            && filename[2] == 'd' && filename[3] == ':')
        {
#if !HAVE_POPEN
            gt_abort
               (">E The \"cmd:\" option is not available in this version.\n");
#else
            filename += 4;
            while (*filename == ' ') ++filename;
            f = popen(filename,"r");
#endif
            assumefixed = FALSE;
            is_pipe = TRUE;
        }
        else
            f = fopen(filename,"r");

        if (f == NULL)
        {
             fprintf(stderr,">E opengraphfile: can't open %s\n",filename);
            return NULL;
        }
    }

    FLOCKFILE(f);
    firstc = c = GETC(f);
    if (c == EOF)
    {
        *codetype = GRAPH6;
        FUNLOCKFILE(f);
        return f;
    }

    if (c != '>')
    {
        *codetype = firstc == ':' ? SPARSE6 : firstc == '&' ? DIGRAPH6 : GRAPH6;
        ungetc(c,f);
    }
    else
    {
        bad_header = FALSE;
        if ((c = GETC(f)) == EOF || c != '>')
            bad_header = TRUE;
        if (!bad_header && ((c = GETC(f)) == EOF || 
	         (c != 'g' && c != 's' && c != 'p')))
            bad_header = TRUE;        
        if (!bad_header && c == 'g')
        {
            if ((c = GETC(f)) == EOF || c != 'r' ||
                (c = GETC(f)) == EOF || c != 'a' ||
                (c = GETC(f)) == EOF || c != 'p' ||
                (c = GETC(f)) == EOF || c != 'h' ||
                (c = GETC(f)) == EOF || c != '6' ||
                (c = GETC(f)) == EOF || c != '<' ||
                (c = GETC(f)) == EOF || c != '<')
                    bad_header = TRUE;
            else
                *codetype = GRAPH6 | HAS_HEADER;
        }
        else if (!bad_header && c == 'd')
        {
            if ((c = GETC(f)) == EOF || c != 'i' ||
                (c = GETC(f)) == EOF || c != 'g' ||
                (c = GETC(f)) == EOF || c != 'r' ||
                (c = GETC(f)) == EOF || c != 'a' ||
                (c = GETC(f)) == EOF || c != 'p' ||
                (c = GETC(f)) == EOF || c != 'h' ||
                (c = GETC(f)) == EOF || c != '6' ||
                (c = GETC(f)) == EOF || c != '<' ||
                (c = GETC(f)) == EOF || c != '<')
                    bad_header = TRUE;
            else
                *codetype = DIGRAPH6 | HAS_HEADER;
        }
        else if (!bad_header && c == 'e')
        {
            if ((c = GETC(f)) == EOF || c != 'd' ||
                (c = GETC(f)) == EOF || c != 'g' ||
                (c = GETC(f)) == EOF || c != 'e' ||
                (c = GETC(f)) == EOF || c != '_' ||
                (c = GETC(f)) == EOF || c != 'c' ||
                (c = GETC(f)) == EOF || c != 'o' ||
                (c = GETC(f)) == EOF || c != 'd' ||
                (c = GETC(f)) == EOF || c != 'e' ||
                (c = GETC(f)) == EOF || c != '<' ||
                (c = GETC(f)) == EOF || c != '<')
                    bad_header = TRUE;
            else
                *codetype = EDGECODE | HAS_HEADER;
        }
        else if (!bad_header && c == 's')
        {
            if ((c = GETC(f)) == EOF || c != 'p' ||
                (c = GETC(f)) == EOF || c != 'a' ||
                (c = GETC(f)) == EOF || c != 'r' ||
                (c = GETC(f)) == EOF || c != 's' ||
                (c = GETC(f)) == EOF || c != 'e' ||
                (c = GETC(f)) == EOF || c != '6' ||
                (c = GETC(f)) == EOF || c != '<' ||
                (c = GETC(f)) == EOF || c != '<')
                    bad_header = TRUE;
            else
                *codetype = SPARSE6 | HAS_HEADER;
        }
        else if (!bad_header && c == 'p')
        {
            if ((c = GETC(f)) == EOF || c != 'l' ||
                (c = GETC(f)) == EOF || c != 'a' ||
                (c = GETC(f)) == EOF || c != 'n' ||
                (c = GETC(f)) == EOF || c != 'a' ||
                (c = GETC(f)) == EOF || c != 'r' ||
                (c = GETC(f)) == EOF || c != '_' ||
                (c = GETC(f)) == EOF || c != 'c' ||
                (c = GETC(f)) == EOF || c != 'o' ||
                (c = GETC(f)) == EOF || c != 'd' ||
                (c = GETC(f)) == EOF || c != 'e')
		    bad_header = TRUE;
	    else
	    {
		if ((c = GETC(f)) == EOF)
		    bad_header = TRUE;
		else if (c == ' ')
		{
		    if ((bl = GETC(f)) == EOF || (bl != 'l' && bl != 'b') ||
			(c = GETC(f)) == EOF || c != 'e' ||
			(c = GETC(f)) == EOF || c != '<' ||
			(c = GETC(f)) == EOF || c != '<')
			    bad_header = TRUE;
		    else if (bl == 'l')
                	*codetype = PLANARCODELE | HAS_HEADER;
		    else
                	*codetype = PLANARCODEBE | HAS_HEADER;
		}
		else if (c == '<')
		{
		    if ((c = GETC(f)) == EOF || c != '<')
			    bad_header = TRUE;
                    else
                        *codetype = PLANARCODE | HAS_HEADER;
		}
		else
                    bad_header = TRUE;
	    }
        }

        if (bad_header)
        {
            fprintf(stderr,">E opengraphfile: illegal header in %s\n",
                    filename == NULL ? "stdin" : filename);
            *codetype = UNKNOWN_TYPE | HAS_HEADER;
            FUNLOCKFILE(f);
            return NULL;
        }
    }

    if (position <= 1) return f;

    if (*codetype&PLANARCODEANY)
    {
	fprintf(stderr,
      ">E opengraphfile: planar_code files can only be opened at the start\n");
	*codetype = UNKNOWN_TYPE | HAS_HEADER;
            FUNLOCKFILE(f);
	    fclose(f);
            return NULL;
    }

    if (*codetype&EDGECODE)
    {
	fprintf(stderr,
      ">E opengraphfile: edge_code files can only be opened at the start\n");
	*codetype = UNKNOWN_TYPE | HAS_HEADER;
            FUNLOCKFILE(f);
	    fclose(f);
            return NULL;
    }

    if (!assumefixed || (*codetype&SPARSE6) || firstc == ':')
    {
        l = 1;
        while ((c = GETC(f)) != EOF)
        {
            if (c == '\n')
            {
                ++l;
                if (l == position) break;
            }
        }
        if (l == position) return f;

        fprintf(stderr,
           ">E opengraphfile: can't find line %ld in %s\n",position,
            filename == NULL ? "stdin" : filename);
        return NULL;
    }
    else
    {
        pos1 = FTELL_VER(f);
        if (pos1 < 0)
        {
            fprintf(stderr,">E opengraphfile: error on first ftell\n");
            return NULL;
        }

        for (i = 1; (c = GETC(f)) != EOF && c != '\n'; ++i) {}
        ogf_linelen = i;

        if (c == EOF)
        {
            fprintf(stderr,
                    ">E opengraphfile: required record no present\n");
            FUNLOCKFILE(f);
            return NULL;
        }
        
        pos2 = FTELL_VER(f);
        if (pos2 < 0)
        {
            fprintf(stderr,">E opengraphfile: error on second ftell\n");
            return NULL;
        }

        pos = pos1 + (position-1)*(pos2-pos1);
        if (FSEEK_VER(f,pos,SEEK_SET) < 0)
        {
            fprintf(stderr,">E opengraphfile: seek failed\n");
            return NULL;
        }
    }

    FUNLOCKFILE(f);
    return f;
}

/*********************************************************************/

void
writeline(FILE *f, char *s)
/* write a line with error checking */
/* \n is not appended automatically */
{
    size_t slen;

    slen = strlen(s);

    if (fwrite(s,1,slen,f) != slen || ferror(f))
        gt_abort(">E writeline : error on writing\n");
}

/*********************************************************************/
/* This function used to be called getline(), but this was changed due
   to too much confusion with the GNU function of that name.
*/

char*
gtools_getline(FILE *f)     /* read a line with error checking */
/* includes \n (if present) and \0.  Immediate EOF causes NULL return. */
{
    DYNALLSTAT(char,s,s_sz);
    int c;
    long i;

    DYNALLOC1(char,s,s_sz,5000,"gtools_getline");

    FLOCKFILE(f);
    i = 0;
    while ((c = GETC(f)) != EOF && c != '\n')
    {
        if (i == s_sz-3)
            DYNREALLOC(char,s,s_sz,3*(s_sz/2)+10000,"gtools_getline");
        s[i++] = (char)c;
    }
    FUNLOCKFILE(f);

    if (i == 0 && c == EOF) return NULL;

    if (c == '\n') s[i++] = '\n';
    s[i] = '\0';
    return s;
}

/****************************************************************************/

char*
getecline(FILE *f)     /* read an edge_code line */
/* No trailing \n or \0 is added.  Immediate EOF causes NULL return. */
{
    size_t headsize,bodysize;
    int sizesize,edgesize;
    int c1,c,i;
    DYNALLSTAT(unsigned char,s,s_sz);

    FLOCKFILE(f);
    if ((c1 = GETC(f)) == EOF) return NULL;

    if (c1 > 0)
    {
        bodysize = c1;
        edgesize = 1;
        headsize = 1;
    }
    else
    {
        if ((c = GETC(f)) == EOF)
            gt_abort(">E Incomplete edge_code line\n");
        else
        {
            sizesize = c >> 4;
            edgesize = c & 0xF;
            bodysize = 0;
            for (i = 0; i < sizesize; ++i)
            {
                if ((c = GETC(f)) == EOF)
                    gt_abort(">E Incomplete edge_code line\n");
                else
                    bodysize = (bodysize << 8) + c;
            }
            headsize = 2 + sizesize;
        }
    }

    DYNALLOC1(unsigned char,s,s_sz,headsize+bodysize,"getecline");

    s[0] = (unsigned char)c1;
    if (c1 == 0)
    {
        s[1] = (char)((sizesize << 4) + edgesize);
        for (i = 0; i < sizesize; ++i)
            s[headsize-1-i] = (bodysize >> 8*i) & 0xFF;
    }

    if (bodysize > 0 && fread(s+headsize,bodysize,1,f) != bodysize)
        gt_abort(">E Incomplete edge_code line\n");

    FUNLOCKFILE(f);
    return (char*)s;
}

int
graphsize(char *s)
/* Get size of graph out of graph6, digraph6 or sparse6 string. */
{
    char *p;
    int n;

    if (s[0] == ':' || s[0] == '&') p = s+1;
    else                            p = s;
    n = *p++ - BIAS6;

    if (n > SMALLN) 
    {
        n = *p++ - BIAS6;
        if (n > SMALLN)
        {
            n = *p++ - BIAS6;
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
        }
        else
        {
            n = (n << 6) | (*p++ - BIAS6);
            n = (n << 6) | (*p++ - BIAS6);
        }
    }
    return n;
}

/****************************************************************************/

void
encodegraphsize(int n, char **pp)
/* Encode the size n in a string starting at **p, and reset **p
   to point to the character after the size */
{
    char *p;

    p = *pp;
    if (n <= SMALLN) 
        *p++ = (char)(BIAS6 + n);
    else if (n <= SMALLISHN)
    {
        *p++ = MAXBYTE;
        *p++ = (char)(BIAS6 + (n >> 12));
        *p++ = (char)(BIAS6 + ((n >> 6) & C6MASK));
        *p++ = (char)(BIAS6 + (n & C6MASK));
    }
    else 
    {
        *p++ = MAXBYTE;
        *p++ = MAXBYTE;
        *p++ = (char)(BIAS6 + (n >> 30));
        *p++ = (char)(BIAS6 + ((n >> 24) & C6MASK));
        *p++ = (char)(BIAS6 + ((n >> 18) & C6MASK));
        *p++ = (char)(BIAS6 + ((n >> 12) & C6MASK));
        *p++ = (char)(BIAS6 + ((n >> 6) & C6MASK));
        *p++ = (char)(BIAS6 + (n & C6MASK));
    }

    *pp = p;
}

/****************************************************************************/

void
stringcounts(char *s, int *pn, size_t *pe)
/* Determine number of edges of graph6, digraph6 or sparse6 string */
{
    char *p;
    int i,j,k,x,nb,v,n,need;
    size_t count;
    boolean done;

    n = graphsize(s);
    *pn = n;

    p = s + (s[0] == ':' || s[0] == '&') + SIZELEN(n);
 
    if (s[0] == ':')  /* sparse6 */
    {
        count = 0;

        for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}
        k = 0;
        v = 0;
        done = FALSE;
        while (!done)
        {
            if (k == 0)
            {
                x = *(p++);
                if (x == '\n' || x == '\0')
                {
                    done = TRUE; continue;
                }
                else
                {
                    x -= BIAS6; k = 6;
                }
            }
            if ((x & B(k))) ++v;
            --k;

            need = nb;
            j = 0;
            while (need > 0 && !done)
            {
                if (k == 0)
                {
                    x = *(p++);
                    if (x == '\n' || x == '\0')
                    {
                        done = TRUE; continue;
                    }
                    else
                    {
                        x -= BIAS6; k = 6;
                    }
                }
                if (need >= k)
                {
                    j = (j << k) | (x & M(k));
                    need -= k; k = 0;
                }
                else
                {
                    k -= need;
                    j = (j << need) | ((x >> k) & M(need));
                    need = 0;
                }
            }
            if (done) continue;

            if (j > v)
                v = j;
            else if (v < n)
                ++count;
	}
    }
    else  /* graph6 or digraph6 */
    {
	count = 0;
	for (; *p != '\n' && *p != '\0'; ++p)
	    count += bytecount[*p - BIAS6];
    }

    *pe = count;
}

/****************************************************************************/

void
stringtograph(char *s, graph *g, int m)
/* Convert string (graph6, digraph6 or sparse6 format) to graph. */
/* Assumes g is big enough to hold it.   */
{
    char *p;
    int n,i,j,k,v,x,nb,need;
    size_t ii;
    set *gi,*gj;
    boolean done;

    n = graphsize(s);
    if (n == 0) return;

    p = s + (s[0] == ':' || s[0] == '&') + SIZELEN(n);

    if (TIMESWORDSIZE(m) < n)
        gt_abort(">E stringtograph: impossible m value\n");

    for (ii = m*(size_t)n; --ii > 0;) g[ii] = 0; g[0] = 0;

    if (s[0] != ':' && s[0] != '&')       /* graph6 format */
    {
        k = 1;
        for (j = 1; j < n; ++j)
        {
            gj = GRAPHROW(g,j,m);
    
            for (i = 0; i < j; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    gi = GRAPHROW(g,i,m);
                    ADDELEMENT(gi,j);
                    ADDELEMENT(gj,i);
                }
                x <<= 1;
            }
        }
    }
    else if (s[0] == '&')
    {
        k = 1;
        for (i = 0; i < n; ++i)
        {
            gi = GRAPHROW(g,i,m);
    
            for (j = 0; j < n; ++j)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    ADDELEMENT(gi,j);
                }
                x <<= 1;
            }
        } 
    }
    else    /* sparse6 format */
    {
        for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}

        k = 0;
        v = 0;
	done = FALSE;
        while (!done)
        {
	    if (k == 0)
	    {
		x = *(p++);
		if (x == '\n' || x == '\0')
		{
		    done = TRUE; continue;
		}
		else
		{
		    x -= BIAS6; k = 6;
		}
	    }
	    if ((x & B(k))) ++v;
	    --k;

	    need = nb;
	    j = 0;
	    while (need > 0 && !done)
	    {
		if (k == 0)
		{
		    x = *(p++);
		    if (x == '\n' || x == '\0')
		    {
		        done = TRUE; continue;
		    }
		    else
		    {
		       	x -= BIAS6; k = 6;
		    }
		}
		if (need >= k)
		{
		    j = (j << k) | (x & M(k));
		    need -= k; k = 0;
		}
		else
		{
		    k -= need;
		    j = (j << need) | ((x >> k) & M(need));
		    need = 0;
		}
	    }
	    if (done) continue;

	    if (j > v)
		v = j;
	    else if (v < n)
	    {
		ADDELEMENT(GRAPHROW(g,v,m),j);
		ADDELEMENT(GRAPHROW(g,j,m),v);
	    }
        }
    }
}

/****************************************************************************/

void
stringtograph_inc(char *s, graph *g, int m,
                  graph *prevg, int prevn)
/* Convert string (graph6, digraph6 or sparse6 format) to graph,
   allowing incremental sparse6 format with a prior graph assumed
   to have matching m,n values. 
   If prevg != NULL and type is is6, use prevg as prior graph.
   Assumes g is big enough to hold it. 
   *digraph is set according to the graph type.
*/
{
    char *p;
    int n,i,j,k,v,x,nb,need;
    size_t ii;
    set *gi,*gj;
    boolean done;

    if (s[0] == ';' && !prevg)
        gt_abort(">E stringtograph_inc missing prior graph\n");

    if (s[0] == ';')
    {
	n = prevn;
	p = s + 1;
        for (ii = m*(size_t)n; --ii > 0;) g[ii] = prevg[ii];
        g[0] = prevg[0];
    }
    else 
    { 
        n = graphsize(s);
        p = s + (s[0] == ':' || s[0] == '&') + SIZELEN(n);
        for (ii = m*(size_t)n; --ii > 0;) g[ii] = 0;
        g[0] = 0;
    }

    if (n == 0) return;

    if (TIMESWORDSIZE(m) < n)
        gt_abort(">E stringtograph_inc: impossible m value\n");

    if (s[0] != ':' && s[0] != ';' && s[0] != '&')   /* graph6 format */
    {
        k = 1;
        for (j = 1; j < n; ++j)
        {
            gj = GRAPHROW(g,j,m);
    
            for (i = 0; i < j; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    gi = GRAPHROW(g,i,m);
                    FLIPELEMENT(gi,j);
                    if (i != j) FLIPELEMENT(gj,i);
                }
                x <<= 1;
            }
        }
    }
    else if (s[0] == '&')       /* digraph6 format */
    {
        k = 1;
        for (j = 0; j < n; ++j)
        {
            gj = GRAPHROW(g,j,m);
    
            for (i = 0; i < n; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    FLIPELEMENT(gj,i);
                }
                x <<= 1;
            }
        }
    }
    else    /* sparse6 format */
    {
        for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}

        k = 0;
        v = 0;
	done = FALSE;
        while (!done)
        {
	    if (k == 0)
	    {
		x = *(p++);
		if (x == '\n' || x == '\0')
		{
		    done = TRUE; continue;
		}
		else
		{
		    x -= BIAS6; k = 6;
		}
	    }
	    if ((x & B(k))) ++v;
	    --k;

	    need = nb;
	    j = 0;
	    while (need > 0 && !done)
	    {
		if (k == 0)
		{
		    x = *(p++);
		    if (x == '\n' || x == '\0')
		    {
		        done = TRUE; continue;
		    }
		    else
		    {
		       	x -= BIAS6; k = 6;
		    }
		}
		if (need >= k)
		{
		    j = (j << k) | (x & M(k));
		    need -= k; k = 0;
		}
		else
		{
		    k -= need;
		    j = (j << need) | ((x >> k) & M(need));
		    need = 0;
		}
	    }
	    if (done) continue;

	    if (j > v)
		v = j;
	    else if (v < n)
	    {
		FLIPELEMENT(GRAPHROW(g,v,m),j);
		if (j != v) FLIPELEMENT(GRAPHROW(g,j,m),v);
	    }
        }
    }
}

/***********************************************************************/

graph*                 /* read graph into nauty format */
readgg(FILE *f, graph *g, int reqm, int *pm, int *pn, boolean *digraph) 
/* graph6, digraph6 and sparse6 formats are supported 
   f = an open file 
   g = place to put the answer (NULL for dynamic allocation) 
   reqm = the requested value of m (0 => compute from n) 
   *pm = the actual value of m 
   *pn = the value of n 
   *digraph = whether the input is a digraph
*/
{
    char *s,*p;
    int m,n;

    if ((readg_line = gtools_getline(f)) == NULL) return NULL;

    s = readg_line;
    if (s[0] == ':')
    {
        readg_code = SPARSE6;
        *digraph = FALSE;
        p = s + 1;
    }
    else if (s[0] == '&')
    {
	readg_code = DIGRAPH6;
        *digraph = TRUE;
	p = s + 1;
    }
    else
    {
        readg_code = GRAPH6;
        *digraph = FALSE;
        p = s;
    }

    while (*p >= BIAS6 && *p <= MAXBYTE) 
        ++p;
    if (*p == '\0')
        gt_abort(">E readgg: missing newline\n");
    else if (*p != '\n')
        gt_abort(">E readgg: illegal character\n");

    n = graphsize(s);
    if (readg_code == GRAPH6 && p - s != G6LEN(n))
        gt_abort(">E readgg: truncated graph6 line\n");
    if (readg_code == DIGRAPH6 && p - s != D6LEN(n))
        gt_abort(">E readgg: truncated digraph6 line\n");

    if (reqm > 0 && TIMESWORDSIZE(reqm) < n)
        gt_abort(">E readgg: reqm too small\n");
    else if (reqm > 0)
        m = reqm;
    else
        m = (n + WORDSIZE - 1) / WORDSIZE;

    if (g == NULL)
    {
        if ((g = (graph*)ALLOCS(n,m*sizeof(graph))) == NULL)
            gt_abort(">E readgg: malloc failed\n");
    }

    *pn = n;
    *pm = m;

    stringtograph(s,g,m);
    return g;
}

/***********************************************************************/

graph*                 /* read undirected graph into nauty format */
readg(FILE *f, graph *g, int reqm, int *pm, int *pn) 
/* graph6 and sparse6 formats are supported 
   f = an open file 
   g = place to put the answer (NULL for dynamic allocation) 
   reqm = the requested value of m (0 => compute from n) 
   *pm = the actual value of m 
   *pn = the value of n 

   Only allows undirected graphs.
*/
{
    boolean digraph;
    graph *gg;

    gg = readgg(f,g,reqm,pm,pn,&digraph);

    if (!gg) return NULL;
    if (digraph)
        gt_abort(">E readg() doesn't all digraphs; use readgg()\n");
    return gg;
}

/***********************************************************************/

int
checkgline(char *s)
/* Check if s[0..] appears to be a graph input line.  A complete check
   is not performed.  Note that graph input lines must end with \n.
   The value returned is 0 if no errors are found, otherwise:
     1 = missing newline
     2 = illegal character
     3 = graph6 or digraph6 line with wrong length
*/
{
    char *p;
    int n,t;

    if (s[0] == ':' || s[0] == ';')
    {
	t = SPARSE6;
	p = s + 1;
    }
    else if (s[0] == '&')
    {
	t = DIGRAPH6;
	p = s + 1;
    }
    else
    {
	t = GRAPH6;
	p = s;
    }

    while (*p >= BIAS6 && *p <= MAXBYTE)
        ++p;
    if (*p == '\0')
        return 1;
    else if (*p != '\n')
        return 2;

    if (t == GRAPH6)
    {
	n = graphsize(s);
	if (p - s != G6LEN(n)) return 3;
    }

    if (t == DIGRAPH6)
    {
	n = graphsize(s);
	if (p - s != D6LEN(n)) return 3;
    }

    return 0;
}

/***********************************************************************/

graph*                 /* read graph into nauty format */
readgg_inc(FILE *f, graph *g, int reqm, int *pm, int *pn,
 	  graph *prevg, int prevm, int prevn, boolean *digraph) 
/* graph6, digraph6 and sparse6 formats are supported 
   f = an open file 
   g = place to put the answer (NULL for dynamic allocation) 
   reqm = the requested value of m (0 => compute from n) 
          This is ignored for an incremental input.
   *pm = the actual value of m 
   *pn = the value of n 
   *digraph = whether the input is a digraph
   If prevg!=NULL, it is a prior graph for use in case the next
   input is a sparse6 increment.
*/
{
    char *s,*p;
    int m,n;

    if ((readg_line = gtools_getline(f)) == NULL) return NULL;

    s = readg_line;
    if (s[0] == ':')
    {
        readg_code = SPARSE6;
	*digraph = FALSE;
        p = s + 1;
    }
    else if (s[0] == ';')
    {
        readg_code = INCSPARSE6;
	*digraph = FALSE;
        p = s + 1;
    }
    else if (s[0] == '&')
    {
	readg_code = DIGRAPH6;
	*digraph = TRUE;
	p = s + 1;
    }
    else
    {
        readg_code = GRAPH6;
	*digraph = FALSE;
        p = s;
    }

    while (*p >= BIAS6 && *p <= MAXBYTE) 
        ++p;
    if (*p == '\0')
        gt_abort(">E readg_inc: missing newline\n");
    else if (*p != '\n')
        gt_abort(">E readg_inc: illegal character\n");

    if (readg_code == INCSPARSE6)
    {
	if (prevg == NULL) gt_abort(">E readg_inc: missing prior\n");
	n = prevn;
	m = prevm;
    }
    else
    {
        n = graphsize(s);
        if (readg_code == GRAPH6 && p - s != G6LEN(n))
            gt_abort(">E readg_inc: truncated graph6 line\n");
        if (readg_code == DIGRAPH6 && p - s != D6LEN(n))
            gt_abort(">E readg_inc: truncated digraph6 line\n");

        if (reqm > 0 && TIMESWORDSIZE(reqm) < n)
            gt_abort(">E readg_inc: reqm too small\n");
        else if (reqm > 0)
            m = reqm;
        else
            m = SETWORDSNEEDED(n);
    }

    if (g == NULL)
    {
        if ((g = (graph*)ALLOCS(n,m*sizeof(graph))) == NULL)
            gt_abort(">E readg_inc: malloc failed\n");
    }

    *pn = n;
    *pm = m;

    stringtograph_inc(s,g,m,prevg,prevn);

    return g;
}

/***********************************************************************/

graph*              /* read undirected graph into nauty format */
readg_inc(FILE *f, graph *g, int reqm, int *pm, int *pn,
 	  graph *prevg, int prevm, int prevn) 
/* graph6 and sparse6 formats are supported 
   f = an open file 
   g = place to put the answer (NULL for dynamic allocation) 
   reqm = the requested value of m (0 => compute from n) 
          This is ignored for an incremental input.
   *pm = the actual value of m 
   *pn = the value of n 
   *digraph = whether the input is a digraph
   If prevg!=NULL, it is a prior graph for use in case the next
   input is a sparse6 increment.
*/
{
    boolean digraph;
    graph *gg;

    gg = readgg_inc(f,g,reqm,pm,pn,prevg,prevm,prevn,&digraph);

    if (!gg) return NULL;
    if (digraph)
        gt_abort(">E readg_inc() doesn't all digraphs; use readgg_inc()\n");
    return gg;
}

/****************************************************************************/

void
stringtosparsegraph(char *s, sparsegraph *sg, int *nloops)
/* Convert string (graph6, digraph6 or sparse6 format)
 * to sparse graph.
 * Assumes sg exists and is initialised               
 * Also returns the number of loops  */
{
    char *p,*q;
    int n,nde,i,j,k,vv,x,nb,need;
    int *d,*e;
    size_t *v;
    int loops;
    boolean done;

    n = graphsize(s);

    q = s + (s[0] == ':' || s[0] == '&') + SIZELEN(n);

    sg->nv = n;

    DYNALLOC1(size_t,sg->v,sg->vlen,n,"stringtosparsegraph");
    DYNALLOC1(int,sg->d,sg->dlen,n,"stringtosparsegraph");

    v = sg->v;
    d = sg->d;
    for (i = 0; i < n; ++i) d[i] = 0;

    if (s[0] != ':' && s[0] != '&')       /* graph6 format */
    {
        p = q;
        k = 1;
        for (j = 1; j < n; ++j)
        {
            for (i = 0; i < j; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    d[i]++;
                    d[j]++;
                }
                x <<= 1;
            }
        }

        v[0] = 0;
        for (i = 1; i < n; ++i) v[i] = v[i-1]+d[i-1];
        nde = v[n-1]+d[n-1];
        for (i = 0; i < n; ++i) d[i] = 0;

        sg->nde = nde;
        DYNALLOC1(int,sg->e,sg->elen,nde,"stringtosparsegraph");
        e = sg->e;

        p = q;
        k = 1;

        for (j = 1; j < n; ++j)
        {
            for (i = 0; i < j; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    e[v[i]+d[i]++] = j;
                    e[v[j]+d[j]++] = i;
                }
                x <<= 1;
            }
        }

	*nloops = 0;
    }
    else if (s[0] == '&')    /* digraph6 */
    {
        p = q;
        k = 1;
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < n; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                    d[j]++;
                x <<= 1;
            }
        }

        v[0] = 0;
        for (i = 1; i < n; ++i) v[i] = v[i-1]+d[i-1];
        nde = v[n-1]+d[n-1];
        for (i = 0; i < n; ++i) d[i] = 0;

        sg->nde = nde;
        DYNALLOC1(int,sg->e,sg->elen,nde,"stringtosparsegraph");
        e = sg->e;

        p = q;
        k = 1;

	*nloops = 0;
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < n; ++i)
            {
                if (--k == 0)
                {
                    k = 6;
                    x = *(p++) - BIAS6;
                }
        
                if ((x & TOPBIT6))
                {
                    e[v[j]+d[j]++] = i;
		    if (i == j) ++*nloops;
                }
                x <<= 1;
            }
        }
    }
    else    /* sparse6 format */
    {
        for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}

        p = q;

        k = 0;
        vv = 0;
        done = FALSE;
	loops = 0;
        while (!done)
        {
            if (k == 0)
            {
                x = *(p++);
                if (x == '\n' || x == '\0')
                {
                    done = TRUE; continue;
                }
                else
                {
                    x -= BIAS6; k = 6;
                }
            }
            if ((x & B(k))) ++vv;
            --k;

            need = nb;
            j = 0;
            while (need > 0 && !done)
            {
                if (k == 0)
                {
                    x = *(p++);
                    if (x == '\n' || x == '\0')
                    {
                        done = TRUE; continue;
                    }
                    else
                    {
                        x -= BIAS6; k = 6;
                    }
                }
                if (need >= k)
                {
                    j = (j << k) | (x & M(k));
                    need -= k; k = 0;
                }
                else
                {
                    k -= need;
                    j = (j << need) | ((x >> k) & M(need));
                    need = 0;
                }
            }
            if (done) continue;

            if (j > vv)
                vv = j;
            else if (vv < n)
            {
                d[vv]++;
                if (vv != j) d[j]++;
		else         ++loops;
            }
        }

        v[0] = 0;
        for (i = 1; i < n; ++i) v[i] = v[i-1]+d[i-1];
        nde = v[n-1]+d[n-1];
        for (i = 0; i < n; ++i) d[i] = 0;

        sg->nde = nde;
        DYNALLOC1(int,sg->e,sg->elen,nde,"stringtosparsegraph");
        e = sg->e;

	p = q;

        k = 0;
        vv = 0;
        done = FALSE;
        while (!done)
        {
            if (k == 0)
            {
                x = *(p++);
                if (x == '\n' || x == '\0')
                {
                    done = TRUE; continue;
                }
                else
                {
                    x -= BIAS6; k = 6;
                }
            }
            if ((x & B(k))) ++vv;
            --k;

            need = nb;
            j = 0;
            while (need > 0 && !done)
            {
                if (k == 0)
                {
                    x = *(p++);
                    if (x == '\n' || x == '\0')
                    {
                        done = TRUE; continue;
                    }
                    else
                    {
                        x -= BIAS6; k = 6;
                    }
                }
                if (need >= k)
                {
                    j = (j << k) | (x & M(k));
                    need -= k; k = 0;
                }
                else
                {
                    k -= need;
                    j = (j << need) | ((x >> k) & M(need));
                    need = 0;
                }
            }
            if (done) continue;

            if (j > vv)
                vv = j;
            else if (vv < n)
            {
                e[v[vv]+d[vv]++] = j;
                if (vv != j) e[v[j]+d[j]++] = vv;
            }
        }
	*nloops = loops;
    }
}

/***********************************************************************/

sparsegraph*               /* read graph into sparsegraph format */
read_sgg_loops(FILE *f, sparsegraph *sg, int *nloops, boolean *digraph) 
/* graph6, digraph6 and sparse6 formats are supported
 * f = an open file
 * sg = place to put the answer (NULL for dynamic allocation) 
 *      - must be initialised if not NULL 
 * nloops := number of loops (each loop in a sparse6 string
 *        gives one loop in the sparse representation)
 */
{
    char *s,*p;
    int n,loops;

    if ((readg_line = gtools_getline(f)) == NULL) return NULL;

    s = readg_line;
    if (s[0] == ':')
    {
        readg_code = SPARSE6;
	*digraph = FALSE;
        p = s + 1;
    }
    else if (s[0] == '&')
    {
	readg_code = DIGRAPH6;
        *digraph = TRUE;
	p = s + 1;
    }
    else
    {
        readg_code = GRAPH6;
        *digraph = FALSE;
        p = s;
    }

    while (*p >= BIAS6 && *p <= MAXBYTE) 
        ++p;
    if (*p == '\0')
        gt_abort(">E read_sg: missing newline\n");
    else if (*p != '\n')
        gt_abort(">E read_sg: illegal character\n");

    n = graphsize(s);
    if (readg_code == GRAPH6 && p - s != G6LEN(n))
        gt_abort(">E read_sg: truncated graph6 line\n");
    if (readg_code == DIGRAPH6 && p - s != D6LEN(n))
        gt_abort(">E read_sg: truncated digraph6 line\n");

    if (sg == NULL)
    {
        if ((sg = (sparsegraph*)ALLOCS(1,sizeof(sparsegraph))) == NULL)
            gt_abort(">E read_sg: malloc failed\n");
        SG_INIT(*sg);
    }

    stringtosparsegraph(s,sg,&loops);
    *nloops = loops;

    return sg;
}

/***********************************************************************/

sparsegraph*          /* read undirected graph into sparsegraph format */
read_sg_loops(FILE *f, sparsegraph *sg, int *nloops) 
/* graph6 and sparse6 formats are supported
 * f = an open file
 * sg = place to put the answer (NULL for dynamic allocation) 
 *      - must be initialised if not NULL 
 * nloops := number of loops (each loop in a sparse6 string
 *        gives one loop in the sparse representation)
 * digraph = whether input line was a digraph
 */
{
    sparsegraph *sgg;
    boolean digraph;

    sgg = read_sgg_loops(f,sg,nloops,&digraph);
    if (!sgg) return NULL;
    if (digraph) gt_abort(">E read_sg_loops() can't handle digraphs,"
                  " use read_sgg_loops()\n");
    return sgg;
}

/***********************************************************************/

sparsegraph*                 /* read graph into sparsegraph format */
read_sg(FILE *f, sparsegraph *sg)
/* graph6 and sparse6 formats are supported
 *  *f = an open file
 *  *sg = place to put the answer (NULL for dynamic allocation) 
 *      - must be initialised if not NULL 
 */
{
    int loops;
    sparsegraph *sgg;
    boolean digraph;

    sgg = read_sgg_loops(f,sg,&loops,&digraph); 
    if (!sgg) return NULL;
    if (digraph) gt_abort(">E read_sg() can't handle digraphs,"
                  " use read_sgg_loops()\n");
    return sgg;
}

/****************************************************************************/

DYNALLSTAT(char,gcode,gcode_sz);  /* Used by ntog6, ntos6, ntod6 and sgtos6 */
TLS_ATTR size_t s6len;
TLS_ATTR int readg_code;
TLS_ATTR char *readg_line;

/****************************************************************************/

char*
ntod6(graph *g, int m, int n)
/* convert nauty graph to digraph6 string, including \n and \0 */
{
    int i,j,k;
    char *p,x;
    set *gj;
    size_t ii;

    ii = D6LEN(n)+3;

    DYNALLOC1(char,gcode,gcode_sz,ii,"ntod6");

    p = gcode;
    *p++ = '&';
    encodegraphsize(n,&p);

    k = 6;
    x = 0;

    for (j = 0; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        for (i = 0; i < n; ++i)
        {
            x <<= 1;
            if (ISELEMENT(gj,i)) x |= 1;
            if (--k == 0)
            {
                *p++ = (char)(BIAS6 + x);
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) *p++ = (char)(BIAS6 + (x << k));

    *p++ = '\n';
    *p = '\0';

    return gcode;
}

/****************************************************************************/

char*
ntog6(graph *g, int m, int n)
/* convert nauty graph to graph6 string, including \n and \0 */
{
    int i,j,k;
    char *p,x;
    set *gj;
    size_t ii;

    ii = G6LEN(n)+3;

    DYNALLOC1(char,gcode,gcode_sz,ii,"ntog6");

    p = gcode;
    encodegraphsize(n,&p);

    k = 6;
    x = 0;

    for (j = 1; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        for (i = 0; i < j; ++i)
        {
            x <<= 1;
            if (ISELEMENT(gj,i)) x |= 1;
            if (--k == 0)
            {
                *p++ = (char)(BIAS6 + x);
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) *p++ = (char)(BIAS6 + (x << k));

    *p++ = '\n';
    *p = '\0';

    return gcode;
}

/****************************************************************************/

char*
ntos6(graph *g, int m, int n)
/* convert nauty graph to sparse6 string, including \n and \0 */
{
    int i,j,k;
    char *p,x;
    set *gj;
    size_t ii;
    int r,rr,topbit,nb,lastj;
    char *plim;

    DYNALLOC1(char,gcode,gcode_sz,5000,"ntos6");

    plim = gcode + gcode_sz - 20;

    gcode[0] = ':';
    p = gcode+1;
    encodegraphsize(n,&p);

    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb)
    {}
    topbit = 1 << (nb-1);
    k = 6;
    x = 0;

    lastj = 0;
    for (j = 0; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        for (i = 0; i <= j; ++i)
        {
            if (ISELEMENT(gj,i))
            {
                if (p >= plim)
                {
                    ii = p - gcode;
                    DYNREALLOC(char,gcode,gcode_sz,
                                       3*(gcode_sz/2)+10000,"ntos6");
                    p = gcode + ii;
                    plim = gcode + gcode_sz - 20;
                }
                if (j == lastj)
                {
                    x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
                else
                {
                    x = (x << 1) | (char)1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                    if (j > lastj+1)
                    {
                        for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
                        {
                            if ((rr & topbit)) x = (x << 1) | (char)1;
                            else               x <<= 1;
                            if (--k == 0)
                            {
                                *p++ = (char)(BIAS6 + x);
                                k = 6;
                                x = 0;
                            }
                        }
                        x <<= 1;
                        if (--k == 0)
                        {
                            *p++ = (char)(BIAS6 + x);
                            k = 6;
                            x = 0;
                        }
                    }
                    lastj = j;
                }
                for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
                {
                    if ((rr & topbit)) x = (x << 1) | (char)1;
                    else               x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
            }
        }
    }

    if (k != 6)
    {
        if (k >= nb+1 && lastj == n-2 && n == (1<<nb))
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << (k-1)) - 1)));
        else
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << k) - 1)));
    }

    *p++ = '\n';
    *p = '\0';
    s6len = p - gcode;
    return gcode;
}

/****************************************************************************/

char*
ntois6(graph *g, graph *prevg, int m, int n)
/* convert nauty graph to incremental sparse6 string, including \n and \0.
   prevg == NULL implies there is no prior graph */
{
    int i,j,k;
    char *p,x;
    set *gj,*pgj;
    setword gdiff;
    size_t ii;
    int r,rr,topbit,nb,lastj,iw,nwords;
    char *plim;

    if (!prevg) return ntos6(g,m,n);

    DYNALLOC1(char,gcode,gcode_sz,5000,"ntois6");

    plim = gcode + gcode_sz - 20;

    gcode[0] = ';';
    p = gcode+1;

    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb)
    {}
    topbit = 1 << (nb-1);
    k = 6;
    x = 0;

    lastj = 0;
    for (j = 0; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        pgj = GRAPHROW(prevg,j,m); 
        nwords = SETWORDSNEEDED(j+1);
	for (iw = 0; iw < nwords; ++iw)
        {
	    gdiff = gj[iw] ^ pgj[iw];
            if (TIMESWORDSIZE(iw+1) > j+1) gdiff &= ALLMASK(SETBT(j+1));
	    while (gdiff)
	    {
		TAKEBIT(i,gdiff);
		i += TIMESWORDSIZE(iw);
   
                if (p >= plim)
                {
                    ii = p - gcode;
                    DYNREALLOC(char,gcode,gcode_sz,
                                       3*(gcode_sz/2)+10000,"ntois6");
                    p = gcode + ii;
                    plim = gcode + gcode_sz - 20;
                }
                if (j == lastj)
                {
                    x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
                else
                {
                    x = (x << 1) | (char)1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                    if (j > lastj+1)
                    {
                        for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
                        {
                            if ((rr & topbit)) x = (x << 1) | (char)1;
                            else               x <<= 1;
                            if (--k == 0)
                            {
                                *p++ = (char)(BIAS6 + x);
                                k = 6;
                                x = 0;
                            }
                        }
                        x <<= 1;
                        if (--k == 0)
                        {
                            *p++ = (char)(BIAS6 + x);
                            k = 6;
                            x = 0;
                        }
                    }
                    lastj = j;
                }
                for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
                {
                    if ((rr & topbit)) x = (x << 1) | (char)1;
                    else               x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
            }
        }
    }

    if (k != 6)
    {
        if (k >= nb+1 && lastj == n-2 && n == (1<<nb))
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << (k-1)) - 1)));
        else
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << k) - 1)));
    }

    *p++ = '\n';
    *p = '\0';
    s6len = p - gcode;
    return gcode;
}

/*************************************************************************/

char*
sgtos6(sparsegraph *sg)
/* Convert undirected sparse graph to sparse6 string including '\n'.
  It is null-terminated and its address (static memory) is returned.
  The length, not including the null, is put in s6len. */
{
    int *d,*e;
    int i,j,n;
    char *p,x,*plim;
    int nb,topbit;
    int dj,k,lastj;
    int r,rr;
    size_t ii,*v,vj,l;

    SG_VDE(sg,v,d,e);
    n = sg->nv;
    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}

    ii = (size_t)(nb+1)*(n/6+sg->nde/3);
    DYNALLOC1(char,gcode,gcode_sz,ii+1000,"sgtos6");
    plim = gcode + gcode_sz - 20;

    p = gcode;
    *p++ = ':';
    encodegraphsize(n,&p);

    topbit = 1 << (nb-1);
    k = 6;
    x = 0;

    lastj = 0;
    for (j = 0; j < n; ++j)
    {
        vj = v[j];
        dj = d[j];
        for (l = 0; l < dj; ++l)
        {
            i = e[vj+l];
            if (i <= j)
            {
                if (p >= plim)
                {
                    ii = p - gcode;
                    DYNREALLOC(char,
                               gcode,gcode_sz,5*(gcode_sz/4)+1000,"sgtos6");
                    p = gcode + ii;
                    plim = gcode + gcode_sz - 20;
                }
                if (j == lastj)
                {
                    x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
                else
                {
                    x = (x << 1) | (char)1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                    if (j > lastj+1)
                    {
                        for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
                        {
                            if ((rr & topbit)) x = (x << 1) | (char)1;
                            else               x <<= 1;
                            if (--k == 0)
                            {
                                *p++ = (char)(BIAS6 + x);
                                k = 6;
                                x = 0;
                            }
                        }
                            x <<= 1;
                        if (--k == 0)
                        {
                            *p++ = (char)(BIAS6 + x);
                            k = 6;
                            x = 0;
                        }
                    }
                    lastj = j;
                }
                for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
                {
                    if ((rr & topbit)) x = (x << 1) | (char)1;
                    else               x <<= 1;
                    if (--k == 0)
                    {
                        *p++ = (char)(BIAS6 + x);
                        k = 6;
                        x = 0;
                    }
                }
            }
        }
    }

    if (k != 6)
    {
        if (k >= nb+1 && lastj == n-2 && n == (1<<nb))
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << (k-1)) - 1)));
        else
	    *p++ = (char)(BIAS6 + ((x << k) | ((1 << k) - 1)));
    }

    *p++ = '\n';
    *p = '\0';
    s6len = p - gcode;
    return gcode;
}

/*************************************************************************/

char*
sgtog6(sparsegraph *sg)
/* Convert undirected sparse graph to graph6 string including '\n','\0'.
  It is null-terminated and its address (static memory) is returned. */
{
    int *d,*e,*ei;
    int i,j,n;
    char *p;
    size_t ii,*v,bodylen,org;
    static char g6bit[] = {32,16,8,4,2,1};

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    ii = G6LEN(n)+3;

    DYNALLOC1(char,gcode,gcode_sz,ii,"sgtog6");

    p = gcode;
    encodegraphsize(n,&p);

    bodylen = G6BODYLEN(n);
    for (ii = 0; ii < bodylen; ++ii) p[ii] = 0;
    p[bodylen] = '\n';
    p[bodylen+1] = '\0';

    for (i = 0, org = 0; i < n;  org += i, ++i)
    {
	ei = e + v[i];
	for (j = 0; j < d[i]; ++j)
	    if (ei[j] < i)
	    {
		ii = ei[j] + org;
		p[ii/6] |= g6bit[ii%6];
	    }
    }

    for (ii = 0; ii < bodylen; ++ii) p[ii] += BIAS6;

    return gcode;
}

/*************************************************************************/

char*
sgtod6(sparsegraph *sg)
/* Convert undirected sparse graph to digraph6 string including '\n','\0'.
  It is null-terminated and its address (static memory) is returned. */
{
    int *d,*e,*ei;
    int i,j,n;
    char *p;
    size_t ii,*v,bodylen,org;
    static char g6bit[] = {32,16,8,4,2,1};

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    ii = D6LEN(n)+3;

    DYNALLOC1(char,gcode,gcode_sz,ii,"sgtog6");

    p = gcode;
    *p++ = '&';
    encodegraphsize(n,&p);

    bodylen = D6BODYLEN(n);
    for (ii = 0; ii < bodylen; ++ii) p[ii] = 0;
    p[bodylen] = '\n';
    p[bodylen+1] = '\0';

    for (i = 0, org = 0; i < n;  org += n, ++i)
    {
	ei = e + v[i];
	for (j = 0; j < d[i]; ++j)
	{
	    ii = ei[j] + org;
	    p[ii/6] |= g6bit[ii%6];
	}
    }

    for (ii = 0; ii < bodylen; ++ii) p[ii] += BIAS6;

    return gcode;
}

/**************************************************************************/

void
writeg6(FILE *f, graph *g, int m, int n)
/* write graph to file in graph6 format */
{
    writeline(f,ntog6(g,m,n));
}

/**************************************************************************/

void
writed6(FILE *f, graph *g, int m, int n)
/* write graph to file in digraph6 format */
{
    writeline(f,ntod6(g,m,n));
}

/**************************************************************************/

void
writes6(FILE *f, graph *g, int m, int n)
/* write graph to file in sparse6 format */
{
    char *s;

    s = ntos6(g,m,n);

    if (fwrite(s,1,s6len,f) != s6len || ferror(f))
        gt_abort(">E writes6 : error on writing\n");
}

/**************************************************************************/

void
writeis6(FILE *f, graph *g, graph *prevg, int m, int n)
/* write graph to file in incremental sparse6 format
   prevg can be NULL if there is no previous graph */
{
    char *s;

    s = ntois6(g,prevg,m,n);

    if (fwrite(s,1,s6len,f) != s6len || ferror(f))
        gt_abort(">E writeis6 : error on writing\n");
}

/**************************************************************************/

void
writeg6_sg(FILE *f, sparsegraph *g)
/* write undirected sparse graph to file in sparse6 format */
{
    writeline(f,sgtog6(g));
}

/**************************************************************************/

void
writed6_sg(FILE *f, sparsegraph *g)
/* write undirected sparse graph to file in sparse6 format */
{
    writeline(f,sgtod6(g));
}

/**************************************************************************/

void
writes6_sg(FILE *f, sparsegraph *g)
/* write undirected sparse graph to file in sparse6 format */
{
    char *s;

    s = sgtos6(g);

    if (fwrite(s,1,s6len,f) != s6len || ferror(f))
        gt_abort(">E writes6 : error on writing\n");
}

/**************************************************************************/

DYNALLSTAT(unsigned char,buff,buff_sz);

void
writepc_sg(FILE *f, sparsegraph *sg)
/* write a sparse graph in planar_code format
    *f = an open file
    *sg = the graph to write
*/
{
    int bytes;
    size_t i,j,len,k,*v,vi;
    unsigned int w;
    int n,*d,*e,di;

#define BEPUT1(x) buff[j++]=(unsigned char)(x);
#define BEPUT2(x) w=(x); buff[j++]=(unsigned char)((w>>8)&0xFF); \
                buff[j++]=(unsigned char)(w&0xff);
#define BEPUT4(x) w=(x); buff[j++]=(unsigned char)((w>>24)&0xFF); \
          buff[j++]=(unsigned char)((w>>16)&0xff); \
          buff[j++]=(unsigned char)((w>>8)&0xFF); \
          buff[j++]=(unsigned char)(w&0xff);

    n = sg->nv;
    SG_VDE(sg,v,d,e);

    if (n <= 255)        bytes = 1;
    else if (n <= 65535) bytes = 2;
    else                 bytes = 4;

    len = bytes * (1 + n + (size_t)(sg->nde));
    if (bytes == 2)      len += 1;
    else if (bytes == 4) len += 3;

    DYNALLOC1(unsigned char,buff,buff_sz,len,"writepc_sg");

    if (bytes == 1)
    {
	j = 0;
	BEPUT1(n);
	for (i = 0; i < n; ++i)
	{
	    vi = v[i];
	    di = d[i];
	    for (k = 0; k < di; ++k) { BEPUT1(e[vi+k]+1); }
	    BEPUT1(0);
	}
    }
    else if (bytes == 2)
    {
	j = 0;
	BEPUT1(n);
	BEPUT2(n);
	for (i = 0; i < n; ++i)
	{
	    vi = v[i];
	    di = d[i];
	    for (k = 0; k < di; ++k) { BEPUT2(e[vi+k]+1); }
	    BEPUT2(0);
	}
    }
    else /* bytes==4 */
    {
	j = 0;
	BEPUT1(n);
	BEPUT2(n);
	BEPUT4(n);
	for (i = 0; i < n; ++i)
	{
	    vi = v[i];
	    di = d[i];
	    for (k = 0; k < di; ++k) { BEPUT4(e[vi+k]+1); }
	    BEPUT4(0);
	}
    }

    if (fwrite((void*)buff,1,j,f) != j)
        gt_abort(">E writepc_sg : error on writing\n");
}

/**************************************************************************/

sparsegraph*
readpc_sg(FILE *f,sparsegraph *sg)
/* read a planar_code graph into sparse graph format 
    *f = an open file
    *sg = place to put the answer (NULL for dynamic allocation) 
        - must be initialised if not NULL 
*/
{
#define BEGET1(x) { x = GETC(f); } 
#define BEGET2(x) { w1=GETC(f); w2=GETC(f); if (w2==EOF) x = EOF; else \
                 x = (w1<<8) | w2; }
#define BEGET4(x) { w1=GETC(f); w2=GETC(f); w3=GETC(f); w4=GETC(f); \
                 if (w4==EOF) x = EOF; \
                 else x = (w1<<24) | (w2<<16) | (w3<<8) | w4; }
    int w1,w2,w3,w4;
    int bytes,n;
    int i,j,*d,*e,di;
    size_t *v,vi;

    BEGET1(n);
    if (n == EOF || n < 0) return NULL;
    else if (n > 0)
	bytes = 1;
    else
    {
	BEGET2(n);
	if (n == EOF || n < 0)
	    gt_abort(">E readpc_sg : error 1 on reading\n");
	else if (n > 0)
	    bytes = 2;
	else
	{
	    BEGET4(n);
	    if (n == EOF || n < 0)
	        gt_abort(">E readpc_sg : error 2 on reading\n");
	    else if (n > 0)
		bytes = 4;
	    else
	        gt_abort(">E readpc_sg : error 3 on reading\n");
	}
    }
 
    if (sg == NULL)
    {
        if ((sg = (sparsegraph*)ALLOCS(1,sizeof(sparsegraph))) == NULL)
            gt_abort(">E readpc_sg: malloc failed\n");
        SG_INIT(*sg);
    }

    SG_ALLOC(*sg,n,2*(size_t)n,"readpc_sg");
    SG_VDE(sg,v,d,e);

    vi = 0;
    for (i = 0; i < n; ++i)
    {
	v[i] = vi;
	di = 0;
	do
	{
	    if      (bytes == 1) BEGET1(j)
            else if (bytes == 2) BEGET2(j)
            else                 BEGET4(j);
	    if (j == EOF) gt_abort(">E readpc_sg : error 4 on reading\n");

	    if (j > 0)
	    {
		if (vi == sg->elen)
		{
		    DYNREALLOC(int,sg->e,sg->elen,2*sg->elen,"readpc_sg");
		    e = sg->e;
		}
	        e[vi++] = j-1;
		++di;
	    }
	    else if (j == 0)
		d[i] = di;
	    else
		gt_abort(">E readpc_sg : error 5 on reading\n");
	} while (j != 0);
    }

    sg->nv = n;
    sg->nde = vi;
    return sg;
}

/**************************************************************************/

sparsegraph*
readpcle_sg(FILE *f,sparsegraph *sg)
/* read a planar_code graph into sparse graph format 
    *f = an open file
    *sg = place to put the answer (NULL for dynamic allocation) 
        - must be initialised if not NULL 
*/
{
#define LEGET1(x) { x = GETC(f); } 
#define LEGET2(x) { w2=GETC(f); w1=GETC(f); if (w1==EOF) x = EOF; else \
                 x = (w1<<8) | w2; }
#define LEGET4(x) { w4=GETC(f); w3=GETC(f); w2=GETC(f); w1=GETC(f); \
                 if (w1==EOF) x = EOF; \
                 else x = (w1<<24) | (w2<<16) | (w3<<8) | w4; }
    int w1,w2,w3,w4;
    int bytes,n;
    int i,j,*d,*e,di;
    size_t *v,vi;

    LEGET1(n);
    if (n == EOF || n < 0) return NULL;
    else if (n > 0)
	bytes = 1;
    else
    {
	LEGET2(n);
	if (n == EOF || n < 0)
	    gt_abort(">E readpcle_sg : error 1 on reading\n");
	else if (n > 0)
	    bytes = 2;
	else
	{
	    LEGET4(n);
	    if (n == EOF || n < 0)
	        gt_abort(">E readpcle_sg : error 2 on reading\n");
	    else if (n > 0)
		bytes = 4;
	    else
	        gt_abort(">E readpcle_sg : error 3 on reading\n");
	}
    }
 
    if (sg == NULL)
    {
        if ((sg = (sparsegraph*)ALLOCS(1,sizeof(sparsegraph))) == NULL)
            gt_abort(">E readpcle_sg: malloc failed\n");
        SG_INIT(*sg);
    }

    SG_ALLOC(*sg,n,2*(size_t)n,"readpcle_sg");
    SG_VDE(sg,v,d,e);

    vi = 0;
    for (i = 0; i < n; ++i)
    {
	v[i] = vi;
	di = 0;
	do
	{
	    if (bytes == 1)      LEGET1(j)
            else if (bytes == 2) LEGET2(j)
            else                 LEGET4(j);
	    if (j == EOF) gt_abort(">E readpcle_sg : error 4 on reading\n");

	    if (j > 0)
	    {
		if (vi == sg->elen)
		{
		    DYNREALLOC(int,sg->e,sg->elen,2*sg->elen,"readpcle_sg");
		    e = sg->e;
		}
	        e[vi++] = j-1;
		++di;
	    }
	    else if (j == 0)
		d[i] = di;
	    else
		gt_abort(">E readpcle_sg : error 5 on reading\n");
	} while (j != 0);
    }

    sg->nv = n;
    sg->nde = vi;
    return sg;
}

/**************************************************************************/

void
writelast(FILE *f)
/* write last graph read by readg() assuming no intervening gtools_getline() */
{
    writeline(f,readg_line);
}

/**************************************************************************/

int
longvalue(char **ps, long *l)
{
    boolean neg,pos;
    long sofar,last;
    char *s;

    s = *ps;
    pos = neg = FALSE;
    if (*s == '-')
    {
        neg = TRUE;
        ++s;
    }
    else if (*s == '+')
    {
        pos = TRUE;
        ++s;
    }

    if (*s < '0' || *s > '9') 
    {
        *ps = s;
        return (pos || neg) ? ARG_ILLEGAL : ARG_MISSING;
    }

    sofar = 0;

    for (; *s >= '0' && *s <= '9'; ++s)
    {
        last = sofar;
        sofar = sofar * 10 + (*s - '0');
        if (sofar < last || sofar > MAXARG)
        {
            *ps = s;
            return ARG_TOOBIG;
        }
    }
    *ps = s;
    *l = neg ? -sofar : sofar;
    return ARG_OK;
}
    
/**************************************************************************/

int
doublevalue(char **ps, double *l)
{
    boolean neg,pos;
    double sofar,weight;
    char *s;

    s = *ps;
    pos = neg = FALSE;
    if (*s == '-')
    {
        neg = TRUE;
        ++s;
    }
    else if (*s == '+')
    {
        pos = TRUE;
        ++s;
    }

    if ((*s < '0' || *s > '9') && *s != '.') 
    {
        *ps = s;
        return (pos || neg) ? ARG_ILLEGAL : ARG_MISSING;
    }

    sofar = 0.0;

    for (; *s >= '0' && *s <= '9'; ++s)
        sofar = sofar * 10 + (*s - '0');
 
    if (*s == '.')
    {
	weight = 1.0;
	for (++s; *s >= '0' && *s <= '9'; ++s)
	{
	    weight /= 10.0;
	    sofar += weight * (*s - '0');
	}
    }

    *ps = s;
    *l = neg ? -sofar : sofar;
    return ARG_OK;
}

/*************************************************************************/

void
arg_long(char **ps, long *val, char *id)
{
    int code;

    code = longvalue(ps,val);
    if (code == ARG_MISSING || code == ARG_ILLEGAL)
    {
        fprintf(stderr,">E %s: missing argument value\n",id);
        gt_abort(NULL);
    }
    else if (code == ARG_TOOBIG)
    {
        fprintf(stderr,">E %s: argument value too large\n",id);
        gt_abort(NULL);
    }
}

/*************************************************************************/

void
arg_int(char **ps, int *val, char *id)
{
    int code;
    long longval;

    code = longvalue(ps,&longval);
    *val = longval;
    if (code == ARG_MISSING || code == ARG_ILLEGAL)
    {
        fprintf(stderr,">E %s: missing argument value\n",id);
        gt_abort(NULL);
    }
    else if (code == ARG_TOOBIG || *val != longval)
    {
        fprintf(stderr,">E %s: argument value too large\n",id);
        gt_abort(NULL);
    }
}

/*************************************************************************/

void
arg_double(char **ps, double *val, char *id)
{
    int code;

    code = doublevalue(ps,val);
    if (code == ARG_MISSING || code == ARG_ILLEGAL)
    {
        fprintf(stderr,">E %s: missing argument value\n",id);
        gt_abort(NULL);
    }
}

/************************************************************************/

boolean
strhaschar(char *s, int c)
/* Check if s contains c.  Saves the bother of figuring out whether
  strchr() is available, or index() or whatever.  */
{
    int i;

    for (i = 0; s[i] != '\0'; ++i)
        if (s[i] == c) return TRUE;

    return FALSE;
}

/************************************************************************/

void
arg_range(char **ps, char *sep, long *val1, long *val2, char *id)
{
    int code;
    char *s;

    s = *ps;
    code = longvalue(&s,val1);
    if (code != ARG_MISSING)
    {
        if (code == ARG_ILLEGAL)
        {
            fprintf(stderr,">E %s: bad range\n",id);
            gt_abort(NULL);
        }
        else if (code == ARG_TOOBIG)
        {
            fprintf(stderr,">E %s: value too big\n",id);
            gt_abort(NULL);
        }
    }
    else if (*s == '\0' || !strhaschar(sep,*s))
    {
        fprintf(stderr,">E %s: missing value\n",id);
        gt_abort(NULL);
    }
    else
        *val1 = -NOLIMIT;

    if (*s != '\0' && strhaschar(sep,*s))
    {
        ++s;
        code = longvalue(&s,val2);
        if (code == ARG_MISSING)
            *val2 = NOLIMIT;
        else if (code == ARG_TOOBIG)
        {
            fprintf(stderr,">E %s: value too big\n",id);
            gt_abort(NULL);
        }
        else if (code == ARG_ILLEGAL)
        {
            fprintf(stderr,">E %s: illegal range\n",id);
            gt_abort(NULL);
        }
    }
    else
        *val2 = *val1;

    *ps = s;
}

/************************************************************************/

void
arg_sequence(char **ps, char *sep,
             long *val, int maxvals, int *numvals, char *id)
{
    int code,ival;
    char *s;

    s = *ps;

    for (ival = 0; ival < maxvals; ++ival)
    {
        code = longvalue(&s,&val[ival]);
        if (code == ARG_ILLEGAL)
        {
            fprintf(stderr,">E %s: illegal value\n",id);
            gt_abort(NULL);
        }
        else if (code == ARG_TOOBIG)
        {
            fprintf(stderr,">E %s: value too big\n",id);
            gt_abort(NULL);
        }
	else if (code == ARG_MISSING)
        {
            fprintf(stderr,">E %s: value missing\n",id);
            gt_abort(NULL);
        }

	if (*s == '\0' || !strhaschar(sep,*s))
	{
	    *numvals = ival+1;
	    *ps = s;
	    return;
	}
	++s;
    }
    fprintf(stderr,">E %s: too many values\n",id);
    gt_abort(NULL);
}

/************************************************************************/

void
arg_doublerange(char **ps, char *sep, double *val1, double *val2, char *id)
{
    int code;
    char *s;

    s = *ps;
    code = doublevalue(&s,val1);
    if (code != ARG_MISSING)
    {
        if (code == ARG_ILLEGAL)
        {
            fprintf(stderr,">E %s: bad range\n",id);
            gt_abort(NULL);
        }
    }
    else if (*s == '\0' || !strhaschar(sep,*s))
    {
        fprintf(stderr,">E %s: missing value\n",id);
        gt_abort(NULL);
    }
    else
        *val1 = -NOLIMIT;

    if (*s != '\0' && strhaschar(sep,*s))
    {
        ++s;
        code = doublevalue(&s,val2);
        if (code == ARG_MISSING)
            *val2 = NOLIMIT;
        else if (code == ARG_ILLEGAL)
        {
            fprintf(stderr,">E %s: illegal range\n",id);
            gt_abort(NULL);
        }
    }
    else
        *val2 = *val1;

    *ps = s;
}

/***********************************************************************/

void
writerange(FILE *f, int c, long lo, long hi)    /* Write a range. */
{
    if (c != '\0') fprintf(f,"%c",c);
    if (lo != -NOLIMIT) fprintf(f,"%ld",lo);
    if (lo != hi)
    {
        fprintf(stderr,":");
        if (hi != NOLIMIT) fprintf(f,"%ld",hi);
    }
}

/************************************************************************/

void
gt_abort(const char *msg)     /* Write message and halt. */
{
    if (msg) fprintf(stderr,"%s",msg);
    ABORT(">E gtools");
}

/************************************************************************/

char*
stringcopy(char *s)   /* duplicate string */
{
    char *scopy;
    size_t i,len;

    for (len = 0; s[len] != '\0'; ++len)
    {}

    if ((scopy = (char*)ALLOCS(len+1,1)) == NULL)
        gt_abort(">E stringcopy: malloc failed\n");

    for (i = 0; i <= len; ++i)
        scopy[i] = s[i];

    return scopy;
}

/*****************************************************************************
*                                                                            *
*  gtools_check() checks that this file is compiled compatibly with the      *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
gtools_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE,"Error: WORDSIZE mismatch in gtools.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE,"Error: MAXM inadequate in gtools.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE,"Error: MAXN inadequate in gtools.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE,"Error: gtools.c version mismatch\n");
        exit(1);
    }
}
