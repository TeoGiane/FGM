#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ggm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_exp_mc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void omp_set_num_cores(void *);
extern void rgwish_c(void *, void *, void *, void *, void *, void *);
extern void rwish_c(void *, void *, void *, void *);
extern void scale_free(void *, void *);
extern void transfer_data(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ggm_bdmcmc_map",                         (DL_FUNC) &ggm_bdmcmc_map,                         18},
    {"log_exp_mc",                             (DL_FUNC) &log_exp_mc,                              8},
    {"omp_set_num_cores",                      (DL_FUNC) &omp_set_num_cores,                       1},
    {"rgwish_c",                               (DL_FUNC) &rgwish_c,                                6},
    {"rwish_c",                                (DL_FUNC) &rwish_c,                                 4},
    {"scale_free",                             (DL_FUNC) &scale_free,                              2},
    {"transfer_data",                          (DL_FUNC) &transfer_data,                           5},
    {NULL, NULL, 0}
};

void R_init_FGM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
