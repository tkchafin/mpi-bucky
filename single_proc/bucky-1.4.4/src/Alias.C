#include"Rand.h"
#include "Alias.h"

//// ps - non-zero probabilities for topologies
//// is - corresponding indices
Alias::Alias(vector<double>& ps, vector<int>& is) {
  int index = 0;
  int nCells = ps.size();
  cells.resize(nCells);
  for(vector<double>::iterator pit = ps.begin(); pit != ps.end(); pit++, index++) {
    double cutoff = *pit * nCells;
    cells[index].cutoff = cutoff;
    cells[index].alias = is[index];
    cells[index].index = is[index];
  }

  computeTwoPointDistributions();
}

int Alias::pick(Rand& rg) {
  double rand = rg.runif(); // rand is in [0, 1)
  double kfr = cells.size() * rand;
  int k = (int) floor(kfr); // k is in [0, n)
  double r = kfr - (double) k;
  double diff = r - cells[k].cutoff;
  if (diff > 0.000000001) {
    return cells[k].alias;
  }

  return cells[k].index;
}

void Alias::computeTwoPointDistributions() {
    vector<Cell>::iterator litr = cells.end() - 1, ritr = cells.begin();

    // sort the vector of multiplied probabilities such that all cells
    // having probability less than 1 are at the end, cells with probabilities
    // greater than 1 are at the beginning
    // e.g., if initial cells are 0.8, 0.9, 1.1, 1.2
    // after this while loop, cells become 1.1, 1.2, 0.8, 0.9
    while (true) {
        while (ritr != cells.end() && ritr->cutoff >= 1.00000) {
            ritr++;
        }
        if (ritr == cells.end()) {
            break;
        }

        while (litr != cells.begin() && litr->cutoff <= 1.00000) {
            litr--;
        }

        if (litr == cells.begin() && litr->cutoff <= 1.00000) {
            break;
        }

        if (litr < ritr)
            break;

        Cell temp = *ritr;
        *ritr = *litr;
        *litr = temp;
        litr--; ritr++;
    }

    vector<Cell>::iterator li = cells.end() - 1, ri = cells.end() - 1;
    while (true) {
        //find a cell with probability initial cutoff greater than 1.0
        while (ri != cells.begin() && ri->cutoff <= 1.00000) {
            ri--;
        }
        if (ri == cells.begin() && ri->cutoff <= 1.00000) {
            break;
        }

        //find a cell with probability initial cutoff less than 1.0
        while (li != cells.begin() && li->cutoff >= 1.00000) {
            li--;
        }

        if (li == cells.begin() && li->cutoff >= 1.00000) {
            break;
        }

        li->alias = ri->index;
        ri->cutoff -= (1 - li->cutoff);
        li--;
    }
}

void Alias::printTwoPointDistributions() {
  cerr << "Two point distribution\n";
  cerr << "Cutoff:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].cutoff << "\t";
  }
  cerr << endl;
  cerr << "Index:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].index << "\t";
  }
  cerr << endl;
  cerr << "Alias:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].alias << "\t";
  }
  cerr << endl;
}

/*// simple test
int main(int argc, char *argv[]) {
  vector<double> ps;



//     ps.push_back(.0231);
//     ps.push_back(.94205);
//     ps.push_back(.03485);


      ps.push_back(0.5);
      ps.push_back(0.5);

//  ps.push_back(1.0095/15);
//  ps.push_back(.9795/15);
//  ps.push_back(.98775/15);
//  ps.push_back(1.023/15);
//  ps.push_back(1.0365/15);
//  ps.push_back(1.005/15);
//  ps.push_back(.96375/15);
//  ps.push_back(1.00575/15);
//  ps.push_back(1.02975/15);
//  ps.push_back(.96825/15);
//  ps.push_back(1.014/15);
//  ps.push_back(.999/15);
//  ps.push_back(.97575/15);
//  ps.push_back(1.01775/15);
//  ps.push_back(.98925/15);

     vector<int> is;
   for (int i = 0; i < ps.size(); i++) {
      is.push_back(i);
  }

  Alias a(ps, is);
  a.printTwoPointDistributions();

  int sum = 0;
  Rand r;
  for (int i = 0; i < 1000000000; i++) {
      sum += a.pick(r);
  }

  cerr << sum << endl;

  ps.clear();
  is.clear();
    ps.push_back(.17);
    ps.push_back(.02);
    ps.push_back(.15);
    ps.push_back(.01);
    ps.push_back(.04);
    ps.push_back(.25);
    ps.push_back(.05);
    ps.push_back(.03);
    ps.push_back(.2);
    ps.push_back(.08);
    for (int i = 0; i < ps.size(); i++) {
       is.push_back(i);
   }

    a = Alias(ps, is);
    a.printTwoPointDistributions();
    sum = 0;
    Rand r1;
    for (int i = 0; i < 100000000; i++) {
        sum += a.pick(r1);
    }
    cerr << sum << endl;
}*/
