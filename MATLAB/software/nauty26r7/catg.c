/* catg.c version 1.3; B D McKay, Jan 2013 */

#define USAGE "catg [-xv] [infile]..."

#define HELPTEXT \
"  Copy files to stdout with all but the first header removed.\n\
\n\
     -x  Don't write a header.\n\
         In the absence of -x, a header is written if\n\
         there is one in the first input file.\n\
\n\
     -v  Summarize to stderr.\n"

#include "gtools.h"

/***********************************************************************/

static FILE*
openfile_head(char *filename, char **header)
/* Open the file for reading and return a FILE*.  If open fails, or
   the header is malformed, write a message and return NULL.  Set
   *header to point to the header string (statically allocated).
   A filename of NULL or "-" means stdin. */
{
    int c,i;
    char *actname;
    FILE *f;
    DYNALLSTAT(char,head,head_sz);

    if (filename == NULL || strcmp(filename,"-") == 0)
    {
        f = stdin;
        actname = "stdin";
    }
    else 
    {
        f = fopen(filename,"r");
        actname = filename;
    }

    if (f == NULL)
    {
        fprintf(stderr,">E catg: can't open file %s\n",actname);
        return NULL;
    }

    DYNALLOC1(char,head,head_sz,100,"catg");

    c = getc(f);
    if (c == '>')
    {
        i = 0;
        head[i++] = (char)c;

        c = getc(f);
        if (c != '>')
        {
            fprintf(stderr,">E catg: bad header in %s\n",actname);
            fclose(f);
            return NULL;
        }
        head[i++] = (char)c;

        do
        {
            c = getc(f);
            if (c == EOF)
            {
                fprintf(stderr,">E catg: bad header in %s\n",actname);
                fclose(f);
                return NULL;
            }

            if (i >= head_sz-1)
                DYNREALLOC(char,head,head_sz,head_sz+100,"catg");
            head[i++] = (char)c;
        }
        while (c != '<' || head[i-2] != '<');

        head[i] = '\0';
    }
    else
    {
	ungetc(c,f);
        head[0] = '\0';
    }

    *header = head;

    return f;
}

/***********************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*head;
    FILE *infile;
    int nfiles,i,j;
    char *arg,sw;
    boolean vswitch,xswitch;
    boolean badargs;
    char buff[1024];
    size_t nr;
    DYNALLSTAT(char*,filename,filename_sz);

    HELP; PUTVERSION;

    DYNALLOC1(char*,filename,filename_sz,200,"catg");

    xswitch = vswitch = FALSE;

    nfiles = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('v',vswitch)
                else SWBOOLEAN('x',xswitch)
                else badargs = TRUE;
            }
        }
        else
        {
            if (nfiles >= filename_sz)
                DYNREALLOC(char*,filename,filename_sz,
                                      filename_sz+200,"catg");
            filename[nfiles++] = arg;
        }
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    for (i = 0; i < (nfiles == 0 ? 1 : nfiles); ++i)
    {
        infilename = (nfiles == 0 ? NULL : filename[i]);
        infile = openfile_head(infilename,&head);
        if (infile == NULL) ABORT("catg");

        if (i == 0 && !xswitch) fprintf(stdout,"%s",head);

        while ((nr = fread(buff,1,1024,infile)) > 0)
        {
            fwrite(buff,1,nr,stdout);
            if (ferror(stdout))
            {
                fprintf(stderr,">E catg: error in writing to stdout\n");
                ABORT("catg");
            }
        }
        if (ferror(infile))
        {
            fprintf(stderr,">E catg: error in reading from %s\n",
                    infilename == NULL ? "stdin" : infilename);
            ABORT("catg");
        }
        fclose(infile);
    }

    if (vswitch)  { }

    exit(0);
}
