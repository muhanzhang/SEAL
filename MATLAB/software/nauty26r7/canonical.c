#include "mex.h"
#include "nausparse.h"
#include <stdint.h>
#include <assert.h>
/*#include <iostream>
//using namespace std
// For invoking nauty to get the canonical labeling of a graph
// MATLAB usage: canon_labels = canonical(subgraph, num_edges, degrees, colors)
// compile using mex: 
       mex canonical.c nauty.c nautil.c naugraph.c schreier.c naurng.c nausparse.c
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int n = mxGetM(prhs[0]);
  int ne = (int)*mxGetPr(prhs[1]);


  plhs[0] = mxCreateNumericMatrix(n, 1, mxINT32_CLASS,  mxREAL);
  int *output = mxGetData(plhs[0]);
  double* gr = mxGetPr(prhs[0]);
  double* deg = mxGetPr(prhs[2]);
  double* color = mxGetPr(prhs[3]);

  
  DYNALLSTAT(int, orbits, orbits_sz);
  /*DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, lab, lab_sz);
  */
  int i, j;
  int lab[n];
  for (i = 0; i < n; ++i){
    /*lab[i] = (int)(lab0[i] + 0.5);*/
    lab[i] = i;
  }
  int ptn[n];
  for (i = 0; i < n; ++i){
    ptn[i] = (int)(color[i] + 0.5);
  }
  
  static DEFAULTOPTIONS_SPARSEGRAPH(options);
  statsblk stats;
  sparsegraph sg;
  sparsegraph canong;

 
  options.writeautoms = FALSE;
  options.getcanon = TRUE;
  options.defaultptn = FALSE;
  options.schreier = TRUE;

  SG_INIT(sg);
  SG_INIT(canong);
  int m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
  /*
  int *lab = mxGetData(plhs[0]);
  */

  /*DYNALLOC1(int,lab,lab_sz,n,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,n,"malloc");*/
  DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
  
  SG_ALLOC(sg,n,ne,"malloc");
  sg.nv = n;
  sg.nde = ne;

  int degcount = 0;
  int k = 0;
  for (i = 0; i < n; ++i){
    sg.v[i] = degcount;
    sg.d[i] = (int)(deg[i] + 0.5);
    k = 0;
    for (j = 0; j < n; ++j){
      if ((int)(gr[j * n + i] + 0.5) == 0) continue;
      else {
        sg.e[degcount + k] = j;
        ++k;
      }
    }
    assert((int)(deg[i] + 0.5) == k);
    degcount += k; 
  }

  /*
  for (i = 0; i < degcount; ++i){
    printf("%d", sg.e[i]);
  }
  printf("\n");
  for (i = 0; i < n; ++i){
    printf("%d", sg.d[i]);
  }
  printf("\n");
  */


  sparsenauty(&sg, lab, ptn, orbits, &options, &stats, &canong);
  

  for (i = 0; i < n; ++i){
    output[i] = lab[i];
  }

  /*
  for (i = 0; i < degcount; ++i){
    printf("%d", canong.e[i]);
  }
  printf("\n");
  for (i = 0; i < n; ++i){
    printf("%d", canong.d[i]);
  }
  printf("\n");
  */

}
