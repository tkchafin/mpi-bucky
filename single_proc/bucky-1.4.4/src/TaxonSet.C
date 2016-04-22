#include"TaxonSet.h"
double logA(double x,int n)
{
  if(n<0) {
    cerr << "Error: logA called with negative n = " << n << endl;
    exit(1);
  }
  if(n==0)
    return 0.0;
  double answer = log(x);
  for(int i=1;i<n;i++)
    answer += log(i + x);
  return answer;
}

// General function to return log( number of unrooted binary trees with n taxa)
double logUB(unsigned int n){
  if(n<1){
    cerr << "Error: logUB called with non-positive n = " << n << endl;
    exit(1);
  }
  if (n<4) return 0.0;
  double answer = 0.0;
  for (unsigned int i=4;i<n+1;i++)
    answer += log(static_cast<double>(2*i - 5));
  return answer;
}

int TaxonSet::numTaxa = 0;
TaxonSet TaxonSet::t;

TaxonSet TaxonSet::getTaxaWithInitialValue(int bitToSet, double weight) {
  TaxonSet t(bitToSet);
  t.weight = weight;
  return t;
}

TaxonSet TaxonSet::getTaxa(int taxacnt) {
  numTaxa = taxacnt;
  t.resize(numTaxa, true);
  return t;
}

TaxonSet TaxonSet::flip() {
  dynamic_bitset<> temp = taxa;
  temp.flip();
  return TaxonSet(temp, weight);
}

void TaxonSet::updatePriorProbability() {
  vector<int> left, right;
  for(int i=1;i<=numTaxa;i++) {
    if(taxa[i - 1]) left.push_back(i);
    else right.push_back(i);
  }
  priorProbability = exp(logUB(left.size()+1) + logUB(right.size()+1)-logUB(numTaxa));
}

TaxonSet TaxonSet::excludeTaxon(int num) {
  assert (num >= 0 && num < taxa.size());
  if (!taxa[num]) {
    return TaxonSet(taxa, weight);
  }

  dynamic_bitset<> result(taxa);
  result = result.flip(num);
  return TaxonSet(result, weight);
}

unsigned int TaxonSet::minus(int num) {
  dynamic_bitset<> tosub(taxa.size(), num);
  dynamic_bitset<> diff = taxa - tosub;
  return diff.to_ulong();
}

// Two sets are compatible if one is the subset of the other or a subset of the complement of the other or vice versa
bool TaxonSet::isCompatible(TaxonSet z) {
  dynamic_bitset<> yz = taxa & z.taxa;
  if (yz == taxa || yz == z.taxa) {
    return true;
  }

  TaxonSet yflip = flip();
  yz = yflip.taxa & z.taxa;
  if (yz == yflip.taxa || yz == z.taxa) {
    return true;
  }

  return false;
}

void TaxonSet::resize(int num, bool set) {
  numTaxa = num;
  taxa.resize(numTaxa, set);
}

void TaxonSet::print(ostream& f) const {
  vector<int> left;
  vector<int> right;
  for(int i=1;i<=numTaxa;i++) {
    if(taxa[i-1])
      left.push_back(i);
    else
      right.push_back(i);
  }
  f << "{";
  if(left.size()>0)
    f << left[0];
  for(int i=1;i<left.size();i++)
    f << "," << left[i];
  f << "|";
  if(right.size()>0)
    f << right[0];
  for(int i=1;i<right.size();i++)
    f << "," << right[i];
  f << "}";
}
