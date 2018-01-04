#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

// Constants used by exact test
static double const OVER_THRESH  = -1;
static double const UNDER_THRESH = 0;

///////////////////////////////////////////////////////////////////////////////
// Utilities

// Store precomputed factorials
double *lnfacs;

double lnfac(int a){
    double z;
    int y;
    if ((a == 1) || a == 0) return 0.0;
    else{
        z = 0;
        for (y = 2; y<=a; y++ ) z = log(y)+z;
        return z;
    }
}

void precompute_factorials(int N){
    int i;
    lnfacs = malloc(sizeof(double) * (N+1));
    for (i=0; i < N+1; i++) {
		lnfacs[i] = lnfac(i);
	}
}

// Numeric functions
int min(int x, int y){ return x < y ? x : y; }
int max(int x, int y){ return x > y ? x : y; }

/* Free the pointers in the array, then free the array. */
void free_ptr_array(void **array, int len) {
    int i;
    for (i=0; i < len; i++) {
        free(array[i]);
    }
    free(array);
}

int sum_cells( int *tbl, int *cells, int num_entries ){
  int i, sum = 0;
  for (i = 0; i < num_entries; i++){
    sum += tbl[cells[i]];
  }
  return sum;
}

///////////////////////////////////////////////////////////////////////////////
// Functions for computing the exact test
void fixed_cells( int k, int pos, int val, int *cells ){
  // Increment pos so counting starts at one
  pos += 1;

  // Define variables
  int *before, *after;
  int i, j, cell_num;

  int num_before  = pow(2, pos-1);
  int num_after   = pow(2, k-pos);

  // Initialize arrays
  before = malloc(sizeof(int) * num_before);
  after = malloc(sizeof(int) * num_after);

  // Compute all integers of k bits before and after the fixed position
  for (i=0; i < num_before; i++) before[i] = i;
  for (i=0; i < num_after; i++) after[i]  = i << pos;

  // Combine the binary strings and set the bit in the given position to 1
  cell_num = 0;
  for (i = 0; i < num_before; i++){
    for (j=0; j < num_after; j++){
      cells[cell_num] = (before[i] + after[j]) | (val << (pos-1));
      cell_num += 1;
    }
  }

  // Clean up memory
  free(before);
  free(after);

}


double denom(int k, int num_entries, int N, int *tbl){
  int i;
  double total = (k-1)*lnfacs[N];
  for (i=0; i < num_entries; i++){ total += lnfacs[tbl[i]]; }
  return total;
}

// Create a list of all k-bit binary strings with one 1
int *get_ex_cells(int k){
  int *ex_cells;
  int i;
  ex_cells = malloc(sizeof(int) * k);
  for (i=0; i < k; i++) ex_cells[i] = 1 << i;
  return ex_cells;
}

// Counts the number of ones in the binary representation of an integer
int num_ones(int n){
  int count = 0;
  while (n > 0){
    count += 1;
    n &= n - 1;
  }
  return count;
}

// Create a list of all k-bit binary strings with >1 1s
int *get_co_cells(int k){
  int *co_cells;
  int i, cell_num, num_co_cells = pow(2, k)-k-1;
  co_cells = malloc(sizeof(int) * num_co_cells);
  cell_num = 0;
  for (i=pow(2, k); i > 0; i--){
    if( num_ones(i) > 1){
      co_cells[cell_num] = i;
      cell_num++;
    }
  }
  return co_cells;
}

void derive_remaining_cells(int k, int N, int *margins, int *ex_cells, int *tbl, int *mar_rems){
  int i;
  for (i = 0; i < k; i++){
    tbl[ex_cells[i]] = mar_rems[i];
  }
  tbl[0] = N;
  for (i = 1; i < pow(2, k); i++) tbl[0] -= tbl[i];
}

int min_affected_margin(int k, int cell, int *mar_rems){
  int i, min_val;
  min_val = INT_MAX;
  for (i = 0; i < k; i++){
    // Check if the ith variable is one for the given cell
    // and replace if it is less than the current min
    if ( cell & (1 << i) && mar_rems[i] <= min_val ) min_val = mar_rems[i];
  }
  return min_val;
}

// Determines if any integer in an array is negative
int contains_negative(int *arr, int len){
  int i;
  for (i = 0; i < len; i++){
    if (arr[i] < 0) return 1;
  }
  return 0;
}

int exact_test_helper(double *pval, int *num_tbls, int k, double pvalthresh, int num_entries,
                      int N, double numerator, int *margins, int *ex_cells, int *co_cells,
                      int num_co_cells, int *tbl, int **mar_stack, int co_in, int T_rem, int T_obs){
  int res = UNDER_THRESH;
  if (co_in >= num_co_cells){
    derive_remaining_cells( k, N, margins, ex_cells, tbl, mar_stack[co_in] );
    if (contains_negative(tbl, num_entries) == 0){
      double addp = exp( numerator - denom( k, num_entries, N,  tbl) );
      pval[0] += addp;
      if (T_obs < sum_cells(tbl, ex_cells, k)){ // T > T_x
        pval[1] += addp;
      }
      num_tbls[0] += 1;
    }
  if ((pval[0]+pval[1])/2 > pvalthresh) {
      res = OVER_THRESH;
    }
  }
  else {
    // Define required variables
    int i, cell, val, MarRem;
    double coef;
    int *mar_rems;

    cell     = co_cells[co_in];
    coef     = num_ones( cell );
    mar_rems = mar_stack[co_in];

    // Determine which variables are in the margin
    MarRem = min_affected_margin( k, cell, mar_rems );

    // Iterate over the possible values the current cell can take
    for (val = 0; val < min(MarRem, (int) floor(T_rem/coef)) + 1; val++){
      // Update margins
      for (i=0; i < k; i++){
        if (cell & (1 << i)) mar_stack[co_in+1][i] = mar_rems[i] - val;
        else mar_stack[co_in+1][i] = mar_rems[i];
      }

      // Create new table using the current value
      tbl[cell] = val;
      res = exact_test_helper( pval, num_tbls, k, pvalthresh, num_entries, N, numerator,
                               margins, ex_cells, co_cells, num_co_cells, tbl,
                               mar_stack, co_in + 1, T_rem-coef*val, T_obs);
      if (res < 0) {
        break;
      }
    }
  }
  return res;
}

double comet_exact_test(int k, int N, int *ctbl, int *final_num_tbls, double pvalthresh){
  // Computation variables
  int kbar, Tobs, T, margin_size, num_co_cells;
  double numerator;
  double *pval;
  int *ex_cells, *co_cells, *blank_tbl, *margins, *num_tbls, *cells;
  int **mar_stack;

  int num_entries, i, res; // Helper variables
  double final_p_value;

  num_entries = 1 << k;

  // Compute binary representation of each cell
  ex_cells     = get_ex_cells(k);
  co_cells     = get_co_cells(k);
  num_co_cells = pow(2, k)-k-1;

  // Allocate memory stack for marginal remainders
  mar_stack = malloc(sizeof(int *) * (num_co_cells + 1));
  for (i=0; i < num_co_cells + 1; i++) mar_stack[i] = malloc(sizeof(int) * k);
  cells = malloc(sizeof(int) * 1 << (k-1));

  /* Compute additional fixed values */
  // Margins
  margins = malloc(sizeof(int) * 2 * k);
  margin_size = pow(2, k-1);
  for (i = 0; i < k; i++){
    // Negative margin
    fixed_cells(k, i, 0, cells);
    margins[i]   = sum_cells( ctbl, cells, margin_size );

    // Positive margin
    fixed_cells(k, i, 1, cells);
    margins[i+k] = sum_cells( ctbl, cells, margin_size );

    // Initialize margin stack
    mar_stack[0][i] = margins[i+k];
  }
  // Numerator
  numerator = 0.0;
  for (i = 0; i < 2*k; i++) numerator += lnfacs[margins[i]];

  // Observed
  Tobs = sum_cells(ctbl, ex_cells, k);
  kbar = 0;
  for (i = 0; i < k; i++) kbar += margins[i+k];

  /* Set up recursion */
  // Set remaining co-occurrences allowed
  T = kbar - Tobs;
  pval = malloc(sizeof(double) * 2);
  num_tbls = malloc(sizeof(int));
  pval[0] = 0.0;
  pval[1] = 0.0; // mid-pvalue
  num_tbls[0] = 0;

  // Construct a blank table to start with
  blank_tbl = malloc(sizeof(int) * num_entries);

  // Run the exact test recursion which will update the results array
  // with the number of extreme tables and the pval
  res = exact_test_helper( pval, num_tbls, k, pvalthresh, num_entries, N, numerator,
                           margins, ex_cells, co_cells, num_co_cells, blank_tbl,
                           mar_stack, 0, T, Tobs);
  *final_num_tbls = num_tbls[0];
  final_p_value = (pval[0]+pval[1])/2;

  // Free memory
  free(cells);
  free(margins);
  free(pval);
  free(num_tbls);
  free(co_cells);
  free(ex_cells);
  free(blank_tbl);
  free_ptr_array((void **) mar_stack, num_co_cells + 1);

  return (res == OVER_THRESH) ? res : final_p_value;

}

///////////////////////////////////////////////////////////////////////////////
// Wrapper function for calling the exact test from R
SEXP cometexacttest(SEXP k, SEXP N, SEXP tbl_, SEXP pvalthresh){
	SEXP result = PROTECT(allocVector(REALSXP, 1));
	int *num_tbls, *tbl;
	int i, num_entries = length(tbl_);
	tbl = malloc(sizeof(int) * num_entries);
  num_tbls = malloc(sizeof(int));
	for (i = 0; i < num_entries; i++){ tbl[i] = INTEGER(tbl_)[i]; }

  precompute_factorials(asInteger(N));
	REAL(result)[0] = comet_exact_test(asInteger(k), asInteger(N), tbl, num_tbls, asReal(pvalthresh));
	UNPROTECT(1);
  free(tbl);
  free(num_tbls);
	return result;
}
