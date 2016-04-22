#ifndef TaxonSetHDR
#define TaxonSetHDR

#include "boost/dynamic_bitset.hpp"
using namespace boost;

#include<iostream>
#include<cmath>
using namespace std;

// General function to return log( prod_{i=0}^{n-1} (i+x) )
//   which is used in the new model repeatedly.
//   Assume that x > 0.
double logA(double x,int n);

// General function to return log( number of unrooted binary trees with n taxa)
double logUB(unsigned int n);

class TaxonSet {
public:
  static TaxonSet getTaxaWithInitialValue(int bitToSet, double weight);

  static TaxonSet getTaxa(int taxacnt);

  static TaxonSet get() {
    return t;
  }

  TaxonSet() {
    taxa.resize(numTaxa, true);
  }

  TaxonSet flip();
  bool leaf() {
    return taxa.count() == 1;
  }

  bool all() {
    return taxa.count() == taxa.size();
  }

  void setAll() {
    taxa.set();
  }

  double getPriorProbability() {
    return priorProbability;
  }

  void updatePriorProbability();
  TaxonSet excludeTaxon(int num);
  unsigned int minus(int num);
  // Two sets are compatible if one is the subset of the other or a subset of the complement of the other or vice versa
  bool isCompatible(TaxonSet z);
  bool isSubsetOf(const TaxonSet& p) {
    return taxa.is_subset_of(p.taxa);
  }

  bool isSupersetOf(const TaxonSet& p) {
    return (taxa | p.taxa) == taxa;
  }

  void setWeight(double x) {
    weight = x;
  }

  double getWeight() {
    return weight;
  }

  int getNumTaxa() {
    return numTaxa;
  }
  void resize(int num, bool set);
  void print(ostream& f) const;

  bool operator<(const TaxonSet& p) {
    return taxa < p.taxa;
  }

  bool operator>(const TaxonSet& p) {
    return taxa > p.taxa;
  }

  bool operator==(const TaxonSet& p) {
    return taxa == p.taxa;
  }

  bool operator!=(const TaxonSet& p) {
    return taxa != p.taxa;
  }

  TaxonSet operator-(const TaxonSet& other) {
    return TaxonSet(dynamic_bitset<>(taxa - other.taxa), weight);
  }

  TaxonSet operator&(TaxonSet x) {
    return TaxonSet(taxa & x.taxa, weight);
  }

  TaxonSet operator|(TaxonSet x) {
    return TaxonSet(taxa | x.taxa, weight);
  }

  bool operator[] (const int nIndex) {
    return taxa[nIndex];
  }

  friend ostream& operator<< (ostream &out, TaxonSet &t) {
    out << t.taxa;
    return out;
  }

  friend bool operator<(const TaxonSet& foo,const TaxonSet& foo1) {     return foo.taxa < foo1.taxa;  }

  friend bool operator>(const TaxonSet& foo,const TaxonSet& foo1) {     return foo.taxa > foo1.taxa;  }

private:

  TaxonSet(int value) {
    taxa.resize(numTaxa, false);
    taxa[value] = 1;
  }

  TaxonSet(dynamic_bitset<> t, double wt) {
    taxa = t;
    weight = wt;
  }

  dynamic_bitset<> taxa;
  double weight;
  double priorProbability;
  static int numTaxa;
  static TaxonSet t;
};

#endif
