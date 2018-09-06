/* rng.h : definitions for using Don Knuth's random number generator.

   To use it:
     1.  Call ran_init(seed) with any long seed.  (Optional,
	   but you will always get the same sequence otherwise.)
     2.  Use NEXTRAN to get the next number (0..2^30-1).
         Alternatively, use KRAN(k) to get a random number 0..k-1.
	 For large k, KRAN(k) is not quite uniform.  In that case
         use GETKRAN(k,var) to set the variable var to a better
         random number 0..k-1.

   Some of these definitions are also in naututil.h.
*/

#ifdef __cplusplus
extern "C" {
#endif

extern long *ran_arr_ptr;
long ran_arr_cycle(void);
void ran_init(long seed);
void ran_array(long *aa, int n);

#ifdef __cplusplus
}
#endif

#define MAXRAN (0x3fffffffL)    /* Values are 0..MAXRAN */
#define NEXTRAN (*ran_arr_ptr>=0 ? *ran_arr_ptr++ : ran_arr_cycle())
#define KRAN(k) (NEXTRAN%(k))
#define RANREAL ((NEXTRAN+0.5)/(MAXRAN+1.0))  /* Uniform (0,1) */

#define MAXSAFE(k) (((MAXRAN+1)/(k))*(k))
#define GETKRAN(k,var) {long __getkran; \
    do {__getkran = NEXTRAN;} while (__getkran >= MAXSAFE(k)); \
    var = __getkran % (k);}
#define INITRANBYTIME ran_init((long)time(NULL))
