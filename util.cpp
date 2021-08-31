#include "util.h"

// *** Utils, CitySim        *** //
// *** jerome.kaempf@epfl.ch *** //

// extern use of LAPACK in Fortran
extern "C" {
    // LU decomposition and backsubstitution of a general matrix
    void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
    // LU decomposition of a general matrix
    void dgetrf_(int* m, int *n, double* a, int* lda, int* ipiv, int* info);
    // generates inverse of a matrix given its LU decomposition
    void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
    // computes all eigenvalues and, optionally, eigenvectors of an n-by-n real symmetric matrix A
    void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
}

// prints a matrix
void print_matrix(string desc, int m, int n, double* a, int lda) {
    int i, j;
    cout << "\n " << desc << endl;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) cout << a[i+j*lda] << " ";
        cout << endl;
    }
    return;
}

// prints a vector of int
void print_int_vector(string desc, int n, int* a) {
    int j;
    cout << "\n " << desc << endl;
    for (j = 0; j < n; j++) cout << a[j] << " ";
    cout << endl;
    return;
}

// inverse a square matrix
void inverse_square_matrix(double* a, int n)
{
    int ipiv[n];
    int lwork = n*n;
    double work[lwork];
    int info;
    //print_matrix("Matrix A", n, n, a, n);
    dgetrf_(&n,&n,a,&n,ipiv,&info);
    if (info<0) throw(string("LAPACK dgetrf: argument" + toString(-info) + " had an illegal value"));
    else if (info>0) throw(string("LAPACK dgetrf: U(" + toString(info) + "," + toString(info) + ") is exactly zero. The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.had an illegal value"));

    dgetri_(&n,a,&n,ipiv,work,&lwork,&info);
    if (info<0) throw(string("LAPACK dgetri: argument" + toString(-info) + " had an illegal value"));
    else if (info>0) throw(string("LAPACK dgetri: U(" + toString(info) + "," + toString(info) + ") is exactly zero; the matrix is singular and its inverse could not be computed."));

    //print_matrix("Matrix A^(-1)", n, n, a, n);
    return;
}

// solve for the problem Ax=b using dgesv
void solve_Ax_equal_b(double* a, double *b, int n)
{
    // run of DGESV - A * x = b (notation: row x column) NPxNP NPx1 = NPx1
    int one = 1;
    int ipiv[n]; // pivot indices
    int info;
    //print_matrix("Matrix A", n, n, a, n);
    //print_matrix("Vector b", n, 1, b, 1);

    dgesv_(&n, &one, a, &n, ipiv, b, &n, &info);
    if (info<0) throw(string("LAPACK dgesv: argument" + toString(-info) + " had an illegal value"));
    else if (info>0) throw(string("LAPACK dgesv: U(" + toString(info) + "," + toString(info) + ") is exactly zero.  The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed."));

    // Print LU factorization
    //print_matrix( "LU factorization", n, n, a, n );
    // Print pivot indices
    //print_int_vector( "Pivot indices", n, ipiv );
    // Print solution
    //print_matrix("Vector x", n, 1, b, 1);
    return;
}

// solve for the eigenvalues of the real symmetric square matrix A
void eigenvalues_A(double *a, double *w, int n)
{
    // run of DSYEV
    char jobz = 'V', uplo = 'U';
    int info;
    int lwork = n*n;
    double work[lwork];
    //print_matrix("Matrix A", n, n, a, n);
    // compute the eigenvalues
    dsyev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, &info);
    // print the eigenvalues
    //print_matrix("Vector w", 1, n, w, 1);
    return;
}

// associate the logstreams
void associate(ostream* pLogStream, ostream& logStream) {
    // the read buffers are associated
    if(pLogStream!=NULL) logStream.rdbuf(pLogStream->rdbuf());
    if (!logStream.good()) throw(string("Unable to define correctly the logStream."));
}

string tabs(unsigned int number) {
    string myTabs = "";
    for (unsigned int i=0; i<number; ++i) myTabs += "\t";
    return myTabs;
}

void save(const string filename, ostringstream &oss, bool overwrite) {

  fstream output;

  if (overwrite)
    output.open(filename.c_str(), ios::out | ios::binary);
  else
    output.open(filename.c_str(), ios::out | ios::binary | ios::app);

  output << oss.str();
  output.close();

  return;

}

//calculates absolute value
uint32_t Abs(int32_t a) {
 if(a < 0)
  return -1*a;
 else
  return a;
}

float nfix(void)
{
const float r = 3.442620f;     /* The start of the right tail */
static float x, y;
 for(;;)
  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }
                         /* iz>0, handle the wedges of other strips */
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

     /* initiate, try to exit for(;;) for loop*/
      hz=SHR3;
      iz=hz&127;
      if( Abs(hz) < kn[iz] ) return (hz*wn[iz]);
  }

}

void zigset(uint32_t jsrseed)
{  const double m1 = 2147483648.0;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;

   jsr^=jsrseed;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(uint32_t) ((dn/q)*m1);
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(unsigned short int i=126;i>=1;i--) // put the index i in the loop declaration, JK - 22/02/09
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(uint32_t) ((dn/tn)*m1);
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }
}

double normallyDistributedSPRNG_Ziggurat() {

    return RNOR;

}

double randomUniform(double minValue,double maxValue)
{
	return minValue + UNI * (maxValue - minValue);
}

float cosAngleBetween(const GENPoint& a, const GENPoint& b) {

    return GEN::dot_product(a,b)/(a.Radius()*b.Radius());

}
