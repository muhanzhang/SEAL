/* naucompare - quick and dirty text file comparitor to compare
 * two program outputs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXLINELEN 100000
#define MAXDIFFS 5

static char line1[MAXLINELEN+1];
static char line2[MAXLINELEN+1];
static char line1mod[MAXLINELEN+100];
static char line2mod[MAXLINELEN+100];

/* special[] contains formats for lines that might not compare exactly.
 * See scanline() for the meaning of the format characters.
 * These strings should have \n at the end.
 */

static char* special[] =   /* Make sure to use \n */
{
    "Dreadnaut version %*f (%*d bits).\n",
    "cpu time = %*f seconds\n",
    "canupdates=%d; cpu time = %*f seconds\n",
    "Mode=%s m=%*d n=%d labelorg=%d edges=%d options=(%s))\n",
    "group time: %*f, %*f, %*f, group order time: %*f;%\n",
    ">Z %d graphs generated in %*f sec\n",
    ">Z %d graphs labelled from %s to %s in %*f sec.\n",
    "%d graphs altogether; cpu=%*f sec\n",
    "group time: %*f,%*f,%*f, order:%*f total:%*f "
       "(Schreier fails: %d); exp_paths time:%*f; aut_check time:%*f\n",
    "%d cell; code = %d; cpu time = %*f seconds\n",
    "%d cells; code = %d; cpu time = %*f seconds\n",
    ">Z %d graphs read from stdin; %d coloured graphs written to stdout; %*f sec\n",
    ">Z %d graphs read from stdin; %d written to stdout; %*f sec\n",
    "%d graphs altogether from %d read; cpu=%*f sec\n",
    ">A %*s %s\n"
};
#define NUMSPECIALS (sizeof(special)/sizeof(*special))

/***********************************************************************/

static int
scanline(char *in, char *fmt, char *out)
/* Compare input string against format and return 1 if it matches or
 * 0 if it doesn't match.  In the case of matching, construct an output
 * string with volatile fields replaced by "*" and some standardisation
 * of the non-volatile fields.  Format specifiers are:
 
   %d  - matches an integer (maybe with sign)
   %f  - matches a real number of the form ddd.ddd (maybe with sign)
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
   %*d, %*f, %*sx, %*c - these are similar to the versions
         without the '*' except that the fields are volatile.
 * The entire line must match.  */
{
    int doass,dots;
    char ends,*outf,*s,*f;

    outf = out;
    s = in;
    f = fmt;

    while (*f != '\0')
    {
        if (*f == '%')
        {
            ++f;
            if (*f == '*')
            {
                doass = 0;
                ++f;
            }
            else
                doass = 1;

            if (*f == '%')
            {
                if (*s++ != '%') return 0;
                ++f;
                *outf++ = '%';
            }
            else if (*f == '\n')
            {
                while (*s != '\0')
                {
                    if (*s != ' ' && *s != '\n') return 0;
                    ++s;
                }
                --s;
            }
            else if (*f == 'c')
            {
                if (*s == ' ' || *s == '\t' || *s == '\n') return 0;
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
		    dots = 0;
                    while (*s == ' ' || *s == '\t')
		    {
			if (doass && dots == 0) *outf++ = ' ';
                        ++s;
			++dots;
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
            else if (*f == 'd')
            {
                while (*s == ' ' || *s == '\t') ++s;
                if (!isdigit(*s) && *s != '-' && *s != '+') return 0;
                if (*s == '-' || *s == '+')
		{
		    if (doass) *outf++ = *s;
		    ++s;
		}
                while (isdigit(*s))
                {
                    if (doass) *outf++ = *s;
		    ++s;
                }
                if (!doass) *outf++ = '*';
                ++f;
            }
            else if (*f == 'f')
            {
                while (*s == ' ' || *s == '\t') ++s;
                if (!isdigit(*s) && *s != '.' && *s != '-' && *s != '+')
                    return 0;
                if (*s == '-' || *s == '+')
		{
		    if (doass) *outf++ = *s;
		    ++s;
		}
		dots = 0;
                while (isdigit(*s) || (dots == 0 && *s == '.'))
                {
		    if (*s == '.') ++dots;
                    if (doass) *outf++ = *s;
		    ++s;
                }
                if (!doass) *outf++ = '*';
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
            if (*s != *f) return 0;
            *outf++ = *f;
            ++s;
            ++f;
        }
    }

    if (*s != '\0') return 0;

    *outf = '\0';

    return 1;
}

/***********************************************************************/

int
main(int argc, char *argv[])
{
    int i,diffs,lineno;
    char *l1,*l2;
    FILE *f1,*f2;

    if (argc != 3)
    {
	fprintf(stderr,"Usage: naucompare file1 file2\n");
	exit(0);
    }

    if (strcmp(argv[1],"-") == 0)
	f1 = stdin;
    else if ((f1 = fopen(argv[1],"r")) == NULL)
    {
	fprintf(stderr,">E naucompare can't open file 1\n");
	exit(1);
    }

    if (strcmp(argv[2],"-") == 0)
	f2 = stdin;
    else if ((f2 = fopen(argv[2],"r")) == NULL)
    {
	fprintf(stderr,">E naucompare can't open file 2\n");
	exit(1);
    }

    diffs = lineno = 0;

    for (;;)
    {
	l1 = fgets(line1,MAXLINELEN,f1);
	l2 = fgets(line2,MAXLINELEN,f2);

	if (l1 == NULL || l2 == NULL) break;

	++lineno;
	if (strcmp(l1,l2) != 0)
	{
	    for (i = 0; i < NUMSPECIALS; ++i)
		if (scanline(l1,special[i],line1mod)
		              && scanline(l2,special[i],line2mod)
                              && strcmp(line1mod,line2mod) == 0)
		    break;
	    if (i == NUMSPECIALS)
	    {
	        if (diffs >= MAXDIFFS)
		{
		    printf("----And more differences not reported.\n");
		    break;
		}
		printf("----Difference on line %d:\n",lineno);
	        printf("%s%s",l1,l2);
	        ++diffs;
	    }
	}
    }

    if (l1 != NULL || l2 != NULL) ++diffs;
    if (l1 != NULL && l2 == NULL) printf("File 1 is longer\n");
    if (l1 == NULL && l2 != NULL) printf("File 2 is longer\n");

    if (diffs == 0) printf("OK\n");

    exit(diffs != 0);
}
