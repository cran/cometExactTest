#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

// RegisteringDynamic Symbols
void R_init_cometexacttest(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
