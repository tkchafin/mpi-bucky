#include<string>
#include<vector>
#include<iostream>
#include <iomanip>
#include<cstdlib>
#include<sstream>
#include<cassert>
#include<algorithm>
#include<functional>
#include<set>
using namespace std;
#include "boost/unordered_map.hpp"
using namespace boost;
#include "TGM.h"
bool cmpTreeWeights(TreeWeight x,TreeWeight y) {  double diff = x.getWeight() - y.getWeight(); if (diff < 0) diff = -diff; if (diff < .00000001) return x.getIndex() < y.getIndex(); return x.getWeight() > y.getWeight(); }
