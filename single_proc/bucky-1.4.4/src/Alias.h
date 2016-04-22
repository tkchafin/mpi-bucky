#ifndef ALIAS_H_
#define ALIAS_H_
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

typedef struct cell {
    int alias;
    int index;
    double cutoff;
} Cell;

class Alias {
public:
  Alias(vector<double>& ps, vector<int>& is);

  virtual ~Alias() {
  }

  int pick(Rand& rg);
  void printTwoPointDistributions();
private:
  vector<Cell> cells;
  void computeTwoPointDistributions();
};

#endif /* ALIAS_H_ */
