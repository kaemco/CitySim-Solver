#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdint.h> // for usage of uint32_t

using namespace std;

// *** Utils, CitySim        *** //
// *** jerome.kaempf@epfl.ch *** //

// Print a matrix
void print_matrix( string desc, int m, int n, double* a, int lda );

// Print a vector of integers
void print_int_vector( string desc, int n, int* a );

// inverse a square matrix
void inverse_square_matrix(double* a, int n);

// solve for the problem Ax=b using dgesv
void solve_Ax_equal_b(double* a, double *b, int n);

// solve for the eigenvalues of the real symmetric square matrix A
void eigenvalues_A(double *a, double *w, int n);

// associate the logstreams
void associate(ostream* pLogStream, ostream& logStream);

string tabs(unsigned int number);

template <typename T> inline T to(const string &s) {

    T value;
    std::stringstream ss(s);
    ss >> value;

    return value;

}

template <typename T> inline string toString(const T &s) {

    std::stringstream ss;
    ss << fixed << s;

    return ss.str();

}

template <class T> ostream& operator<<(ostream &os, const pair<T,T>& p)
{
    os << p.first << "\t" << p.second;
    return os;
}

template <class T> ostream& operator<<(ostream &os, const vector<T>& v)
{
    for (unsigned int i=0; i<v.size()-1; i++) {
        os << v[i] << "\t";
    }
    os << v[v.size()-1];
    return os;
}

template <class T> void save(const string filename, const vector<T> &vector) {

  fstream output (filename.c_str(), ios::out | ios::binary);
  output.setf(ios_base::fixed);
  output.precision(12);

  for (unsigned int i=0;i<vector.size();i++) output << vector[i] << endl;

  output.close();

  return;

}

template <class T> void save(const string filename, const T &value) {

  fstream output (filename.c_str(), ios::out | ios::binary);
  output.setf(ios_base::fixed);
  output.precision(12);

  output << value << endl;
  output.close();

  return;

}

template <class T> void load(const string filename, vector<T> &vector) {

    fstream input (filename.c_str(), ios::in | ios::binary);
    if (!input.is_open()) throw(string("Error loading file: " + filename));

    string buffer;
    do {
        input >> buffer;
        if (input.eof()) break;
        vector.push_back(to<T>(buffer));
    } while (!buffer.empty());

    input.close();

    return;

}

void save(const string filename, ostringstream &oss, bool overwrite=true);

// ******************************** Ziggurat Method for Randn **********************************
// changed for all platforms -> uint32_t, JK - 22/02/09
// note: on a 32-bits machine, the unsigned long int and the unsigned int are of the same size
// (maximum of the 32-bits available = 0..2^32-1)
// HOWEVER, on a 64-bits machine, the unsigned int is still in the range 0..2^32-1, but the
// unsigned long int is now between 0..2^64-1
//
// THEREFORE, we need to use the following definition for the unsigned int:
// typedef unsigned long int uint32_t
// as this method was designed using unsigned int in 32 bits
// *********************************************************************************************

#pragma GCC diagnostic ignored "-Wunused-variable"

static uint32_t jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static int32_t hz;
static uint32_t iz, kn[128];
static float wn[128],fn[128];

#define RNOR (hz=SHR3, iz=hz&127, (Abs(hz) < kn[iz])? hz*wn[iz] : nfix())

#pragma GCC diagnostic warning "-Wunused-variable"

/* nfix() generates variates from the residue when rejection in RNOR occurs. */

float nfix(void);

/*--------This procedure sets the seed and creates the tables------*/

void zigset(uint32_t jsrseed);

// ******************************** End of Ziggurat Method for Randn *************************** //

double normallyDistributedSPRNG_Ziggurat();

double randomUniform(double minValue,double maxValue);

#endif
