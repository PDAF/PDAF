#ifdef SNGLPREC
#define gemmTYPE sgemm
#define gemvTYPE sgemv
#define gesvTYPE sgesv
#define syevTYPE ssyev
#define gesvdTYPE sgesvd
#define syevxTYPE ssyevx
#define potrfTYPE spotrf
#define trtrsTYPE strtrs
#define larnvTYPE slarnv
#define MPI_REALTYPE MPI_REAL
#define WORDLENGTH_REAL 1
#else
#define gemmTYPE dgemm
#define gemvTYPE dgemv
#define gesvTYPE dgesv
#define syevTYPE dsyev
#define gesvdTYPE dgesvd
#define syevxTYPE dsyevx
#define potrfTYPE dpotrf
#define trtrsTYPE dtrtrs
#define larnvTYPE dlarnv
#define MPI_REALTYPE MPI_DOUBLE_PRECISION
#define WORDLENGTH_REAL 2
#endif

