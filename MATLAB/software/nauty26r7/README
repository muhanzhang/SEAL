README file for nauty 2.6 

Brendan McKay, bdm@cs.anu.edu.au
Adolfo Piperno, piperno@di.uniroma1.it

------------------------------------------------------------

The most recent distribution of nauty and Traces can be found at
http://cs.anu.edu.au/~bdm/nauty and http://pallini.di.uniroma1.it .

The manual nug26.pdf is available at that site and is also included
in the distribution package.

Note that nauty and Traces are copyright but free to use for most
purposes. The details are in the file COPYRIGHT.

------------------------------------------------------------

INSTALLATION.

See the manual for more information.

If you have a working shell, and "make", you can run
   ./configure
followed by
   make 
to compile nauty and Traces for your system.

If that succeeds without problem, you will have have the
program dreadnaut ready to run.

There are some options that can be specified at the ./configure
step; see the manual.

If you don't have a shell or make, manually edit the files nauty.h,
naututil.h and gtools.h as distributed.  The parts between the lines
======= near the start are the main things to look at.  After this
manual editing, you can use makefile.basic as a guide to compilation.

Programs which use an older version of nauty need to be
recompiled (** not just relinked **). Make sure they define
the options structure using one of
DEFAULTOPTIONS_GRAPH
DEFAULTOPTIONS_SPARSEGRAPH
DEFAULTOPTIONS_DIGRAPH
DEFAULTOPTIONS_SPARSEDIGRAPH
DEFAULTOPTIONS_TRACES

------------------------------------------------------------

TESTING.

After compiling nauty successfully, it is recommended that you run
the included test programs.  The simplest way is
    make checks

------------------------------------------------------------

MAILING LIST.

There is a mailing list for announcements and discussion about
nauty and related topics.  You can subscribe at
http://mailman.anu.edu.au/mailman/listinfo/nauty

------------------------------------------------------------

OTHER FILES IN THE PACKAGE.

Also in the package (documentation at the start of each source file).

sumlines.c  -  This is a program designed to digest the outputs from
  multiple runs of a program (such as a computation split into multiple
  parts).  Lines matching given patterns can be counted and checked,
  and numbers appearing in them can be accumulated.  Instructions appear
  in the source file.  See the option GMP near the head of the program
  before trying to compile.

sorttemplates.c  - Some carefully tuned generic quicksort procedures.

bliss2dre.c - A program which reads one file in Bliss format and writes
  it in dreadnaut format.

blisstog.c - A program which reads one or more files in Bliss format
  and writes all the graphs in sparse6 format.

poptest.c - A program for testing the POPCOUNT macro.

dretodot.c - A program that reads files in dreadaut format and writes
  dot files suitable for drawing with graphviz.

------------------------------------------------------------

Windows.

For running nauty in Windows, Cygwin is recommended.

If configure gives an error message similar to this:
   can not guess host type: you must specify one
then try
   ./configure --build=unknown

------------------------------------------------------------

Making 32-bit executables on 64-bit Linux systems.

(In bash or sh:)
CFLAGS=-m32 CXXFLAGS=-m32 LDFLAGS=-m32 ./configure
make clean; make

This requires 32-bit libraries to be available.  On Ubuntu
they are called ia32-libs and libc6-dev-i386.

------------------------------------------------------------

RECENT CHANGES.

See the file changes24-26.txt for a longer list.
See the file README_24 for a list of older changes.
