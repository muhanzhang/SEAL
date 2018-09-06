/* naugstrings.c : Write graph6 or sparse6 strings into array. */
/* Version 1.1, Jun 2015. */

#include "naugstrings.h"

/****************************************************************************/

static void
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
gtog6string(graph *g, char **pp, int m, int n)
/* convert dense graph to graph6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int i,j,k;
    char *p,x;
    set *gj;
    size_t ii;

    ii = G6LEN(n)+3;

    p = *pp;
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
    *pp = p;
}

/****************************************************************************/

void
gtod6string(graph *g, char **pp, int m, int n)
/* convert dense graph to digraph6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int i,j,k;
    char *p,x;
    set *gj;
    size_t ii;

    ii = D6LEN(n)+3;

    p = *pp;
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
    *pp = p;
}

/****************************************************************************/

void
gtois6string(graph *g, graph *prevg, char **pp, int m, int n)
/* convert dense graph to incremental spargse6 string, including \n
   but not \0.
   Start at *pp and set *pp to point one char after the \n */
{
    int i,j,k;
    char *p,x;
    set *gj,*pgj;
    setword gdiff;
    int r,rr,topbit,nb,lastj,iw,nwords;

    if (!prevg) 
    {
	gtos6string(g,pp,m,n);
	return;
    }

    **pp = ';';
    p = *pp+1;

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
                            if (rr & topbit) x = (x << 1) | (char)1;
                            else             x <<= 1;
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
                    if (rr & topbit) x = (x << 1) | (char)1;
                    else             x <<= 1;
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
    *pp = p;
}

/****************************************************************************/

void
gtos6string(graph *g, char **pp, int m, int n)
/* convert dense graph to sparse6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int i,j,k;
    char *p,x;
    set *gj;
    int r,rr,topbit,nb,lastj;

    **pp = ':';
    p = *pp+1;
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
                            if (rr & topbit) x = (x << 1) | (char)1;
                            else             x <<= 1;
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
                    if (rr & topbit) x = (x << 1) | (char)1;
                    else             x <<= 1;
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
    *pp = p;
}

/*************************************************************************/

void
sgtos6string(sparsegraph *sg, char **pp)
/* convert sparse graph to graph6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int *d,*e;
    int i,j,n;
    char *p,x;
    int nb,topbit;
    int dj,k,lastj;
    int r,rr;
    size_t ii,*v,vj,l;

    SG_VDE(sg,v,d,e);
    n = sg->nv;
    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}

    ii = (size_t)(nb+1)*(n/6+sg->nde/3);

    p = *pp;
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
                            if (rr & topbit) x = (x << 1) | (char)1;
                            else             x <<= 1;
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
                    if (rr & topbit) x = (x << 1) | (char)1;
                    else             x <<= 1;
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
    *pp = p;
}

/*************************************************************************/

void
sgtog6string(sparsegraph *sg, char **pp)
/* convert dense graph to graph6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int *d,*e,*ei;
    int i,j,n;
    char *p;
    size_t ii,*v,bodylen,org;
    static char g6bit[] = {32,16,8,4,2,1};

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    ii = G6LEN(n)+3;

    p = *pp;
    encodegraphsize(n,&p);

    bodylen = G6BODYLEN(n);
    for (ii = 0; ii < bodylen; ++ii) p[ii] = 0;
    p[bodylen] = '\n';

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

    *pp = p + bodylen + 1;
}

/*************************************************************************/

void
sgtod6string(sparsegraph *sg, char **pp)
/* convert dense graph to digraph6 string, including \n but not \0
   Start at *pp and set *pp to point one char after the \n */
{
    int *d,*e,*ei;
    int i,j,n;
    char *p;
    size_t ii,*v,bodylen,org;
    static char g6bit[] = {32,16,8,4,2,1};

    SG_VDE(sg,v,d,e);
    n = sg->nv;

    ii = G6LEN(n)+3;

    p = *pp;
    encodegraphsize(n,&p);

    bodylen = D6BODYLEN(n);
    for (ii = 0; ii < bodylen; ++ii) p[ii] = 0;
    p[bodylen] = '\n';

    for (i = 0, org = 0; i < n;  org += i, ++i)
    {
	ei = e + v[i];
	for (j = 0; j < d[i]; ++j)
	{
	    ii = ei[j] + org;
	    p[ii/6] |= g6bit[ii%6];
	}
    }

    for (ii = 0; ii < bodylen; ++ii) p[ii] += BIAS6;

    *pp = p + bodylen + 1;
}
