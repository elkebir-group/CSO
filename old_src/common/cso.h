/*
 * cso.h
 *
 *  Created on: 17-feb-2009
 *      Author: M. El-Kebir
 */

#ifndef CSO_H_
#define CSO_H_

//#define TIXML_USE_STL
#define GET_BIT(n, x, i) (((x) >> ((n) - (i) - 1)) & 1)
#define GENERATE_AND_MASK(i, n) (((1 << (n)) - 1) << (i))

#include <assert.h>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <stdlib.h>
#include <math.h>

#ifdef _MSC_VER
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include <limits.h>
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <stack>

// forward declarations
class Genotype;
class GenotypeGamete;
class Data;

typedef std::vector<float> FloatVector;
typedef std::vector<double> DoubleVector;

typedef std::vector<FloatVector> FloatMatrix;
typedef std::vector<DoubleVector> DoubleMatrix;

typedef std::vector<Genotype> GenotypeVector;
typedef std::vector<const Genotype*> GenotypePointerVector;

typedef std::set<Genotype> GenotypeSet;
typedef std::set<const Genotype*> GenotypePointerSet;

typedef std::list<Genotype> GenotypeList;
typedef std::list<const Genotype*> GenotypePointerList;

typedef std::tr1::unordered_map<int, std::tr1::unordered_map<int, GenotypeGamete> > GenotypeGameteMatrix;

typedef struct
{
  int _c;
  double _prob;
} Gamete;

extern const Gamete g_InvalidGamete;
extern const Data* g_pData;
#define INVALID_GAMETE -1

typedef std::list<Gamete> GameteList;
typedef std::vector<Gamete> GameteVector;

typedef enum {Linked, WeaklyLinked, Unlinked} LinkageType;

// function prototypes
void toBitstring(int val, int n, char* buf);
int fromBitstring(int n, const char* buf);
int numberOfDifferences(int n, int val1, int val2);
int chromosomeCompare(int n, int val1, int val2, std::vector<int>& homyzygous, std::vector<int>& heterozygous);
int randInt(int n);
double randDouble();
int randChromosome(int nLoci);
unsigned long probToPop(double p, double gamma);
int largestSubGroupSize(int nLoci, int val, int target);
bool gamete_eq(const Gamete& gamete1, const Gamete& gamete2);
bool gamete_lt(const Gamete& gamete1, const Gamete& gamete2);


#endif /* CSO_H_ */
