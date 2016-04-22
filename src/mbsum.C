// BUCKy - Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)
// Copyright (C) 2006-2102 by Bret Larget

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License, version 2, as
// published by the Free Software Foundation.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License (the file gpl.txt included with this
// distribution or http://www.gnu.org/licenses/gpl.txt) for more
// details.

// version 1.2b: multiple files from the same gene can be combined.
// version 1.3.0: Translate table, if present, are copied verbatim from
//              the first file, assumed the same in other files.
// version 1.4.0: requires translate tables, outputs unrooted formatted trees.
// version 1.4.2: removes new extra characters between `=' and tree in MrBayes 3.2 output

// mbsum.C

// Input:    one or more MrBayes .t files
//           Lines with the format `tree <name> = <parenthetic tree representation>'
//             are read as trees.  Other lines are ignored.
//
// Options:  [-n number-of-skipped-trees] will skip over the first number-of-skipped-trees trees in each input file
//           [--skip number-of-skipped-trees] will skip over the first number-of-skipped-trees trees in each input file
//           [-o outfile] name of file for writing output
//           [-h] prints the brief usage message
//           [--help] prints the brief usage message
//           [--version] prints the version number
//
// Output:   a single file with the name <filename - .extension> + .in (by default)
//             (Output file name is the original file name, stripped of the last extension if there is one, plus .in .)
//           The output file contains one line per sampled tree topology with two fields,
//             the topology (parenthetic representation) and the count.
//           Tree topologies are sorted by count, from highest to lowest.
//           The output file is in the appropriate input format for BUCKy.
//
// Usage:    mbsum [--help || -h] [<--skip || -n> number-of-skipped-trees] [<--out || -o> output-file] [--version] [input filename(s)]

#define VERSION "1.4.4"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include "mbsumtree.h"
using namespace std;
using namespace mbsumtree;

void usageError(char *name)
{
  cerr << "Usage: " << name << " mbsum [--help || -h] [<--skip || -n> number-of-skipped-trees] [<--out || -o> output-file] [--version] [input filename(s)]" << endl;
  exit(1);
}

class TopCountNode {
public:
  TopCountNode() {}
  TopCountNode(string t) : count(1), top(t), left(NULL), right(NULL), parent(NULL) {} // create the root node
  TopCountNode(string t,TopCountNode* p) : count(1), top(t), left(NULL), right(NULL), parent(p) {} // create a child node
  ~TopCountNode() {
    if(left != NULL) {
      left->~TopCountNode();
      delete left;
    }
    if(right != NULL) {
      right->~TopCountNode();
      delete right;
    }
  }
  void createLeft(string t) {
    left = new TopCountNode(t);
  }
  void createRight(string t) {
    right = new TopCountNode(t);
  }
  void pass(string t,TopCountNode* d) {
    d->check(t);
  }
  void check(string t) {
    if(t == top)
      count++;
    else if(t < top) {
      if(left==NULL)
	createLeft(t);
      else
	pass(t,left);
    }
    else { // t > top
      if(right==NULL)
	createRight(t);
      else
	pass(t,right);
    }
  }
  void printSubtree(ostream& f) {
    if(left!=NULL)
      left->printSubtree(f);
    if(right!=NULL)
      right->printSubtree(f);
    print(f);
  }
  void print(ostream& f) {
    f << top << " " << count << endl;
  }
  void collect(vector<TopCountNode*>& cnodes) {
    if(left != NULL)
      left->collect(cnodes);
    if(right != NULL)
      right->collect(cnodes);
    cnodes.push_back(this);
  }
  int getCount() { return count; }
  string getTop() { return top; }
private:
  TopCountNode* left;
  TopCountNode* right;
  TopCountNode* parent;
  int count;
  string top;
};

bool cmpTCNodes(TopCountNode* x,TopCountNode* y) {
  if(x->getCount() != y->getCount())
    return( x->getCount() > y->getCount() );
  else
    return( x->getTop() < y-> getTop() );
}

int main(int argc, char *argv[])
{
  if(argc==1)
    usageError(argv[0]);

  int numSkip = 0;
  string outFile;
  vector<string> fileNames;

  for(int i=1;i<argc;) {
    string a = (string)(argv[i]);
    if(a == "-h" || a == "--help")
      usageError(argv[0]);
    if(a == "--version") {
      cout << "mbsum version " << VERSION << endl;
      exit(1);
    }
    if(a == "-n" || a == "--skip") {
      if(i == argc-1)
	usageError(argv[0]);
      istringstream f((string)(argv[++i]));
      if(f.fail())
	usageError(argv[0]);
      f >> numSkip;
      if(f.fail())
	usageError(argv[0]);
      if(numSkip < 0) {
	cerr << "Error: Cannot skip a negative number of trees." << endl;
	usageError(argv[0]);
      }
    }
    else if(a == "-o" || a == "--out") {
      if(i == argc-1)
	usageError(argv[0]);
      istringstream f((string)(argv[++i]));
      if(f.fail())
	usageError(argv[0]);
      outFile.clear();
      f >> outFile;
      if(f.fail())
	usageError(argv[0]);
    }
    else
      fileNames.push_back(a);
    i++;
  }

  if(outFile.size() == 0) {
    if(fileNames.size() > 0) {
      string fileRoot = (string)(fileNames[0]);
      int ext = fileRoot.rfind('.');
      if(ext!=string::npos) // '.' in filename, erase from last '.'
	fileRoot.erase(ext,fileRoot.length() - ext);
      // erase directory
      ext = fileRoot.rfind('/');
      if(ext!=string::npos)
	fileRoot.erase(0,ext+1);
      outFile = fileRoot + ".in";
    }
    else
      outFile = "bucky.in";
  }

  ofstream sumOut(outFile.c_str());
  if(sumOut.fail()){
    cerr << "Error: Could not open file " << outFile << " ." << endl;
    exit(1);
  }

  int numFiles = fileNames.size();
  vector<int> numTrees(numFiles,0);
  vector<int> numTaxa(numFiles,0);
  TopCountNode* root;
  Pruner* mbsumPruner = new Pruner();
  for(int i=0;i<numFiles;i++) {
    ifstream f(fileNames[i].c_str());
    if(f.fail())
      cerr << "Error: Could not open file " << fileNames[i] << " ." << endl;
    else
      cout << "Reading file " << fileNames[i] << ": " << flush;



    string line;
    int lineNumber=0;
    int numActuallySkipped=0;
    int skipInFile = numSkip;
    bool readTtable = false;
    // determine EOL character
    char foo;
    do {
      foo = f.get();
    } while( (foo != '\n') && (foo != '\r') && !f.fail() );
    if(foo == '\r') {
      foo = f.peek();
      if(foo != '\n')
	foo = '\r';
    }
    // Now use foo as the EOL character

    cerr << endl; // added for debugging!

    while(getline(f,line,foo)) {
      lineNumber++;
      istringstream s(line);
      string keyTree,name,equalSign;

      if (readTtable){
	char ch;
	s >> noskipws >> ch;
	while (ch != ';' && ch != EOF && s.good()) {
	  sumOut << ch;
	  s >> noskipws >> ch;
	}
	if (ch == ';'){
	  sumOut << ";\n";
	  readTtable = false;
	} else {
	  sumOut << "\n";
	  continue;
	}
      }

      s >> keyTree;
      if (i==0){
	// read translate table only in the first file, assuming the same in other files.
	if (keyTree=="translate" || keyTree=="TRANSLATE" || keyTree=="Translate"){
	  readTtable = true;
	  sumOut << "translate";
	  char ch;
	  s >> noskipws >> ch;
	  while (ch != ';' && ch != EOF && s.good()) {
	    sumOut << ch;
	    s >> noskipws >> ch;
	  }
	  if (ch == ';'){
	    sumOut << ";\n";
	    readTtable = false;
	  } else {
	    sumOut << "\n";
	    continue;
	  }
	}
      }

      // skip if the remaining line is not in format "  tree name = treeRep"
      if(keyTree != "tree")
	continue;

      s >> name >> equalSign;
      if(equalSign != "=")
	continue;
      // The rest should be a parenthetic representation of a tree.  Assume no spaces!
      // Skip the first numSkip trees
      if(skipInFile-- > 0) {
	numActuallySkipped++;
	continue;
      }

      // Skip extra output in MrBayes3.2 [&U] before tree string

      while ( !s.fail() ) {
	char c;
	c = s.get();
	if ( c == '(' ) {
	  s.putback(c);
	  break;
	}
      }
      
      string treeString;
      s >> treeString;
      Tree tree(treeString,lineNumber, mbsumPruner);
      numTrees[i]++;

      ostringstream g;
      tree.printTop(g); // changed from mb2badger to not include an endline
      string top = g.str();

      if((i==0) && numTrees[0] == 1) // first tree
	root = new TopCountNode(top);
      else
	root->check(top);
    }
    cout << "Skipped " << numActuallySkipped << " trees, read " << numTrees[i] << " trees." << endl;
  }

  vector<TopCountNode*> cnodes;
  root->collect(cnodes);
  sort(cnodes.begin(),cnodes.end(),cmpTCNodes);

  for(vector<TopCountNode*>::iterator n=cnodes.begin();n!=cnodes.end();n++)
    (*n)->print(sumOut);

  int total = 0;
  for(int i=0;i<numFiles;i++)
    total += numTrees[i];
  cout << "Read " << total << " total trees of which " << cnodes.size() << " are distinct." << endl;
  cout << "Output written to file " << outFile << endl;
  return(0);
}
