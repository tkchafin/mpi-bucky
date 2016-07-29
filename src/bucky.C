// MPI-BUCKY
// Tyler K. Chafin - tkchafin@uark.edu
//
// Modification of: 
//   BUCKy 1.4.4 Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)
//   Copyright (C) 2006-2014 by Bret Larget and Cecile Ane and Riley Larget

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License, version 2, as
// published by the Free Software Foundation.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License (the file gpl.txt included with this
// distribution or http://www.gnu.org/licenses/gpl.txt) for more
// details.

// Version History

// First alpha version 7 March, 2005
// Version 1.1 released 30 October, 2006
// Version 1.2a 21 August, 2007
// Version 1.2b 17 January, 2008
// Version 1.3.0  27 May, 2009
// Version 1.3.1  29 October, 2009
// Version 1.4.3   9 July, 2014
// Version 1.4.4  22 June, 2015
//------MPI------
// Version 2.0.0 2016 Tyler K. Chafin

// File:     bucky.C

// Input:    [options] <summary files>
//
// Options:
//            [-a alpha-value]                --- Dirichlet process hyper-parameter
//            [-n number-of-MCMC-updates]     --- number of MCMC updates per chain
//            [-c number-of-chains]           --- number of separate chains for MCMCMC
//            [-r MCMCMC-rate]                --- rate in updates to propose MCMCMC swap
//            [-m alpha-multiplier]           --- multiple used to determine alpha for heated chains
//            [-s subsample-rate]             --- subsample rate for output file (summaries do not subsample)
//            [-o output-file-root]           --- root name of all output files
//            [-s1 seed1]                     --- 0 < seed1 < 4294967295
//            [-s2 seed2]                     --- 0 < seed1 < 4294967295
//            [--use-independence-prior]
//            [--calculate-pairs]
//            [--create-sample-file]
//            [--create-single-file]
//            [--create-joint-file]
//            [-h] prints the brief usage message
//            [--help] prints the brief usage message
//            [--version] prints brief intro and version number
//
// Output:
//
//  fileRoot.cluster       --- summary of cluster size (number of groups)
//  fileRoot.concordance   --- split-by-split summary of concordance factor
//  fileRoot.gene          --- gene by gene summary of joint posterior information
//  fileRoot.joint         --- table of raw second stage MCMC counts for joint posterior
//  fileRoot.input         --- names of input files
//  fileRoot.out           --- output file
//  fileRoot.pairs         --- output if pairwise analysis is desired (very slow!)
//  fileRoot.sample        --- output of raw sampled GTMs (very big!)
//  fileRoot.single        --- summary of table of single gene posterior probabilities

// Changes in version 1.1
// --- added probabilities to concordance file output
// --- fixed bug with split index that assumed no more than 8 taxa
// --- implemented new and better updateGroups; use it as the default
// --- added runtime argument ---no-use-update-groups

// Changes in version 1.2
// --- Added more information to the --help output including default values
// --- Used Huffman coding idea to create a binary tree for proposing trees that minimizes
//     the average number of comparisons necessary to map the uniform number to the tree
//     which should increase the speed of the program, especially for data sets with many poorly resolved trees.
// --- Added probabilities to the .cluster file.
// --- Changed function getProb which speeds up retrieving probabilities for genes (by a lot!) by using a lookup table.
// --- Fixed error in writing clades in primary concordance tree
// --- Concordance summary includes primary concordance tree
// --- Eliminated some output files, and do not print some output by default. Add controls to print more output.

// Changes in version 1.2b
// --- Added the distribution of genome-wide concordance factors for splits in the primary concordance tree
// --- Fixed error in updatePairCounts, count at the bottom-right
// --- Fixed uninitialized default values for createXxxFile, fixed error in setFileNames

// Changes in version 1.3.0
// --- Use git repository: use git commands to get more details on previous versions.
// --- Added independent runs are diagnostic summaries to check convergence.
// --- Fixed bug in the group update.
// --- Translate tables are read and used to check all taxa are the same.

// Changes in version 1.4.3
// --- Updated included boost library to BOOST 1.55.0
// --- Fixed incorrect return type in set.all()
// --- Added virtual destructor for Table class
// --- Added option to turn off the population tree calculation

// Changes in version 1.4.3
// --- Fixed error with genome-wide CFs, that occurred when grid<ngenes and for low-CF splits.
//

// Changes in version 2.0.0 
// MPI parallel implementation 
// ---Major revision to main() 
// ---multiple-core enabled MCMCMC
// ---Changed way track alphas (query through State::getAlpha()) 
// ---MPI serialization and gather of relevant tables
// ---Serialized many multi-dimensional vectors (states etc)
// Memory Optimization 
// ---Converted TGMTable to ints (32 bits) rather than doubles (64 bits)
// ---Removed --opt-space option (deprecated Table class)
 


//Developer notes
//Known issues:
//--- Threads not closed properly on calling help



#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <limits>
#include <map>
#include <set>
#include "bucky.h"
#include "mbsumtree.h"
#include "mpi.h"

using namespace std;
string PRE="MPI-BUCKy (BUCKy 2.0.0) using OpenMPI Message Passing Interface";
string PRE2="Written by Tyler K. Chafin, last modified June 2016";
string PRE3="Based on :";
string VERSION = "1.4.4";
string DATE = "22 June 2015";
string COPYRIGHT = "Copyright (C) 2006-2015 by Bret Larget, Cecile Ane and Riley Larget\nPlease cite original authors";

int countTaxa(string top) {
  int countTaxa=0;
  istringstream s(top);
  char c;
  int z;
  do {
    c = s.peek();
    if(!isdigit(c))
      s >> c;
    else {
      s >> z;
      countTaxa++;
    }
  } while(c != ';' && c != EOF && s.good());

  return countTaxa;
}

bool getTaxa(string& file, set<string>& taxaInIFile, map<string, int>& translateMap) {
  ifstream f(file.c_str());
  if(f.fail()) {
    cerr <<"Error: Cannot open file " << file << "." << endl;
    MPI_Finalize();
    exit(1);
  }

  // search for the gene's translate table
  bool hasTtable = false;

  string keyword;
  f >> keyword;
  if(keyword=="translate" || keyword=="TRANSLATE" || keyword=="Translate") {
    hasTtable =true;
    // read the gene's translate table
    string taxonname;
    int taxonnumber;
    bool done=false;  // will be done when semicolon is reached
    while(f.good() && !f.eof() && !done){
      f >> taxonnumber >> taxonname;
      if(f.fail()){
        cerr << "Error in reading the translate table for file: "<< file << ". This is probably because of incorrect format of translate table. Please make sure it ends with a ';'." <<endl;
        exit(1);
      }
      size_t pos=taxonname.rfind(",");
      if (pos!=string::npos)
        taxonname.replace(pos,1,"");
      else {
        pos = taxonname.rfind(";");
        if (pos!=string::npos){
          taxonname.replace(pos,1,"");
          done=true;
        }
        else {
          f >> keyword; // take out the comma
          if (keyword==";") done=true;
        }
      }

      taxaInIFile.insert(taxonname);
      if (!translateMap[taxonname]) {
        translateMap[taxonname] = translateMap.size();
      }

    }
  }
  f.close();
  return hasTtable;

}

void getTaxaSubset(int rank, vector<string>& inputFiles, map<string, int>& prunedTranslateMap, vector<string>& translateTable, string pruneFile, bool shouldPruneGene, bool& changed) {
  vector<int> genesToIgnore;
  map<string, int> translateMap;
  set<string> taxaNames;
  size_t fileNum = 1;
  unsigned int maxTaxa = 0;
  string firstFile = inputFiles[0];
  if (!pruneFile.empty()) {
    bool hasTable = getTaxa(pruneFile, taxaNames, translateMap);
    if (!hasTable) {
	  if (rank==0){
        cerr << "\nBucky cannot find the translate table in " << pruneFile
	     << ". The prune file should start with keyword 'translate'."
	     << "Please cross check translate table format. "<< endl;
      }
      exit(0);
    }
    fileNum = 0;
    maxTaxa = taxaNames.size();
    firstFile = pruneFile;
  }
  else {
    bool hasTable = getTaxa(inputFiles[0], taxaNames, translateMap);
    if (!hasTable) {
      if (rank==0){
        cerr << "\nBucky cannot find the translate table in " << inputFiles[0]
	     << ".\nPlease double check the presence or format\n"
	     << "of the translate table in this file.\n"
	     << "All input files should start with keyword 'translate'."<<endl;
	   }
      exit(0);
    }
  }

  for (;fileNum<inputFiles.size();fileNum++) {
    set<string> taxaInIFile;
    bool hasTable = getTaxa(inputFiles[fileNum], taxaInIFile, translateMap);
    if (!hasTable) {
	  if (rank==0){
        cerr << "\nBucky cannot find translate table in " << inputFiles[fileNum]
	     << ".\nPlease double check the presence or format\n"
	     << "of the translate table in this file.\n"
	     << "All input files should start with keyword 'translate'."<< endl;
	   }
      exit(0);
    }

    set<string> output;
    int count = taxaInIFile.size();
    set_intersection(taxaNames.begin(), taxaNames.end(), taxaInIFile.begin(), taxaInIFile.end(), inserter(output, output.begin()));
    if (shouldPruneGene && (taxaNames.size() != output.size() || taxaNames.size() > taxaInIFile.size())) {
        if (rank==0)
        cerr << "Skipping gene " << inputFiles[fileNum] << " as translate table does not match with that of reference: " << firstFile << endl;
        genesToIgnore.push_back(fileNum);
    }
    else {
      if (count != taxaInIFile.size() || count != taxaNames.size()) {
        changed = true;
      }
      if (output.size() < maxTaxa) {
		if(rank==0){
          cerr << "\nCannot prune to taxa subset specified in " << firstFile << endl;
          cerr << inputFiles[fileNum] << " does not have all taxa specified in the prune file" << endl;
	    }
        for (set<string>::iterator itr = output.begin(); itr != output.end(); itr++) {
          taxaNames.erase(*itr);
        }
        if (rank==0){
          cerr << "Missing taxa for this locus:\n";
          for (set<string>::iterator itr = taxaNames.begin(); itr != taxaNames.end(); itr++) {
            cerr << *itr << "\n";
          }
	    }

        exit(0);
      }

      taxaNames = output;
    }
  }

  // Give ids to all taxa left after intersection.
  // Make sure these numbers fall in continuous range.
  // if a taxa ti is given number 33 before and taxa 32 gets pruned ti should get 32 now.
  vector<string> unmappedTaxa;
  int taxaCnt = taxaNames.size();
  translateTable.resize(taxaCnt);
  for (set<string>::iterator itr = taxaNames.begin(); itr != taxaNames.end();itr++) {
    int num = translateMap[*itr];
    translateMap.erase(*itr);
    if (num <= taxaCnt) {
      prunedTranslateMap[*itr] = num;
      translateTable[num - 1] = *itr;
    }
    else {
      unmappedTaxa.push_back(*itr);
    }
  }

  int i = 0;
  for (map<string, int>::iterator itr = translateMap.begin(); itr != translateMap.end(); itr++) {
    if (itr->second <= taxaCnt) {
      prunedTranslateMap[unmappedTaxa[i]] = itr->second;
      translateTable[itr->second - 1] = unmappedTaxa[i];
      i++;
    }
  }

  // if pg option is specified, loci not having all taxa should be ignored, erase them from inputFiles
  for (int i = genesToIgnore.size() - 1; i >= 0; i--) {
      inputFiles.erase(inputFiles.begin() + genesToIgnore[i]);
  }
}

string pruneTop(string top, int lineNum, mbsumtree::Pruner* p) {
  mbsumtree::Tree t(top, lineNum, p);
  ostringstream newtop;
  t.printTop(newtop);
  return newtop.str();
}

void normalize(vector<double>& table) {
  double sum = 0;
  for(int j=0;j<table.size();j++)
    sum += table[j];

  if(sum > 0)
    for(int j=0;j<table.size();j++)
      table[j] /= sum;
}

//
// readFile
//
// input file is expected to have each line contain two fields separated by white space
//   the first field is a topology string (parenthetic representation)
//   the second field is a positive weight (usually a count from an MCMC sample)
// add read topologies to row i of the table
// when new topologies are discovered, add a new column to the table

void readFile(string filename, int i, TGM*& tgm, int &max,
              vector<int>& taxid, map<string, int>& translateMap, bool changed)
{
  vector<double> table;
  vector<string> topologies;
  ifstream f(filename.c_str());
  if(f.fail()) {
    cerr <<"Error: Cannot open file " << filename << "." << endl;
    MPI_Finalize();
    exit(1);
  }

  map<int, int> translateID;
  // in case the numbers in the translate table are not from 1 to Ntax
  // translateID[ taxon number given in the input file] = ID given in
  // bucky's translate table used, common for all genes.

  string keyword;
  int numTaxa=translateMap.size();
  mbsumtree::Pruner* thePruner = new BuckyPruner(numTaxa);
  int lineNum = 1;
  // read the gene's translate table
    string taxonname;
    int taxonnumber;
    f >> keyword; // take out 'translate'
    bool done=false;  // will be done when semicolon is reached
    while(f.good() && !f.eof() && !done){
      f >> taxonnumber >> taxonname;
      if(f.fail()){
	cerr << "Error in reading the translate table for gene "<<i<<endl;
	exit(1);
      }
      lineNum++;
      size_t pos=taxonname.rfind(",");
      if (pos!=string::npos)
	taxonname.replace(pos,1,"");
      else {
	pos = taxonname.rfind(";");
	if (pos!=string::npos){
	  taxonname.replace(pos,1,"");
	  done=true;
	}
	else {
	  f >> keyword; // take out the comma
	  if (keyword==";") done=true;
	}
      }
      int taxonid = translateMap[taxonname];
      if (taxonid != 0){
        translateID[taxonnumber] = taxonid;
        taxid.push_back(taxonid);
      }
    }
    numTaxa=taxid.size();

  // set zeros for ith row of table
  while(f.good() && !f.eof()) {
    string top;
    double weight;
    f >> top >> weight;
    if(f.fail()) {
      if (!top.empty() && top.find("(") == string::npos) {
        cerr << "Invalid topology string " << top << " in file:" << filename
            <<". Please make sure translate table does not contain multiple ';'?"
            << endl;

        exit(0);
      }
      break;
    }

    // replace taxon numbers by taxid in the topology
      istringstream itop(top);
      ostringstream otop;
      char c;
      int z;
      do {
	c = itop.peek();
	if(!isdigit(c)){
	  itop >> c;
	  otop << c;
	}
	else {
	  itop >> z;
	  int taxonid = translateID[z];
	  // if translate table does not have mapping for z we replace it by zero
	  // and zeroes will be pruned later.
	  otop << taxonid;
	}
      } while(c != ';' && c != EOF && itop.good());
      top = otop.str();

    top = pruneTop(top, lineNum, thePruner);

    // determine longest length of topology string
    if(top.length()>max)
      max = top.length();

    topologies.push_back(top);
    table.push_back(weight);
    lineNum++;
  }
  f.close();
  delete thePruner;

  normalize(table);
  vector<double>::iterator citr;
  vector<string>::iterator titr;
  for (titr = topologies.begin(), citr = table.begin(); titr != topologies.end(); titr++, citr++) {
    tgm->addGeneCount(*titr, i, *citr);
  }
}

double State::calculateLogPriorProb()
{
  logPriorProb=0;
  for(int i=0;i<numGroups;i++)
    logPriorProb += logA(alphaOverTop,counts[indices[i]]);
  logPriorProb -= logA(alpha,genes.size());
  return logPriorProb;

}

int State::updateOne(int i,Rand& rand) {
  int oldTop = tops[i];
  int newTop;
  if(i < genes.size())
    //    newTop = genes[i]->pickTree(rand);
    newTop = genes[i]->pickTreeFast(rand);
  else {
    cerr << "Internal Error: i = " << i << " too large in updateOne." << endl;
    exit(1);
  }

  if(newTop==oldTop)
    return 1;
  int n1 = counts[oldTop];
  int n2 = counts[newTop];
  double logHR=0;
  if(!useIndependencePrior) { // logHR=0 if independence prior
    logHR += log(n2 + alphaOverTop);
    logHR -= log(n1-1 + alphaOverTop);
  }
  if(log(rand.runif()) > logHR) // reject
      return 0;
  else { // accept
    tops[i] = newTop;
    counts[oldTop]--;
    counts[newTop]++;
    if(n2==0) {
      numGroups++;
      indices.push_back(newTop);
      sort(indices.begin(),indices.end());
    }
    if(n1==1) {
      numGroups--;
      int j=0;
      while(indices[j]!=oldTop)
	j++;
      indices.erase(indices.begin()+j,indices.begin()+j+1);
    }
    logPriorProb += logHR;
    if(i < genes.size())
      logPosteriorProbProduct += (genes[i]->getLogProb(newTop) - genes[i]->getLogProb(oldTop));
    return 1;
  }
}

/********************************************************************************

This is the old very slow algorithm.

void State::updateOneGroup(int group,Rand& rand) {
  // group is index in indices of the group to update
  // find list of topologies that:
  //  (1) are not topologies for a different group; and
  //  (2) have positive single gene prior for all genes in the group
  // current topology is in this set
  // pick one with prob. proportional to prod_{set} prob(top)

  int oldTop = indices[group];

  // find set of genes currently with oldTop
  vector<int> set(counts[oldTop]); // indices of genes with oldTop
  int j=0;
  for(int i=0;i<genes.size();i++)
    if(tops[i] == oldTop)
      set[j++] = i;

  // loop through topologies to find possible tops and weights
  vector<int> possibleTops(0);
  vector<double> logWeights(0);
  for(int top=0;top<numTreesSampled;top++) {
    if(top == oldTop) {
      possibleTops.push_back(top);
      double w = 0.0;
      for(int i=0;i<set.size();i++) {
	double p = genes[set[i]]->getProb(top);
	if(p>0)
	  w += log(p);
	else {
	  cerr << "Error: Current top has zero probability!" << endl;
	  exit(1);
	}
      }
      logWeights.push_back(w);
    }
    else {
      if(!(counts[top]>0)) { // top not already used by another group
	double w = 0.0;
	bool noZeroProbs = true;
	for(int i=0;i<set.size() && noZeroProbs;i++) {
	  double p = genes[set[i]]->getProb(top);
	  if(p>0)
	    w += log(p);
	  else
	    noZeroProbs = false;
	}
	if(noZeroProbs) {
	  possibleTops.push_back(top);
	  logWeights.push_back(w);
	}
      }
    }
  }

  if(possibleTops.size()==1) // oldTop is the only possible choice
    return;

  // subtract max from each log weight for numerical accuracy
  double m = logWeights[0];
  for(int i=1;i<logWeights.size();i++)
    if(logWeights[i] > m)
      m = logWeights[i];
  double totalWeight=0.0;
  for(int i=0;i<logWeights.size();i++) {
    logWeights[i] = exp(logWeights[i] - m); // set logWeights[i] to weight
    totalWeight += logWeights[i];
  }
  double x = rand.runif() * totalWeight;
  int newTopIndex = 0;
  x -= logWeights[newTopIndex];
  while(x > 0) {
    x -= logWeights[++newTopIndex];
  }
  int newTop = possibleTops[newTopIndex];
  if(newTop == oldTop) // picked oldTop
    return;
  // make changes to newTop
  for(int i=0;i<set.size();i++)
    tops[set[i]] = newTop;
  counts[newTop] = counts[oldTop];
  counts[oldTop] = 0;
  j=0;
  while(indices[j]!=oldTop)
    j++;
  indices.erase(indices.begin()+j,indices.begin()+j+1);
  indices.push_back(newTop);
  sort(indices.begin(),indices.end());
}

*******************************************************************************/

// New, faster update group algorithm
int State::updateOneGroup(int gene,Rand& rand) {
  // pick a random topology for the gene
  // check that it:
  //  (1) is not a topology for another group
  //  (2) has positive single gene prior for all genes in the group
  // return 1 if change accepted, 0 otherwise.

  int oldTop = tops[gene];
  int newTop = genes[gene]->pickTreeFast(rand);
  if(newTop==oldTop)
    return 0;

  // check if newTop is already in use
  if(counts[newTop]>0)
    return 0;

  // possible new topology for group
  vector<int> set(counts[oldTop]); // indices of genes with oldTop
  int j=0;
  double oldLogProb = 0.0;
  double newLogProb = 0.0;
  for(int i=0;i<genes.size();i++) {
    if(tops[i] == oldTop) {
      oldLogProb += genes[i]->getLogProb(oldTop);
      set[j++] = i;
      if (genes[i]->hasProb(newTop))
        newLogProb += genes[i]->getLogProb(newTop);
      else // acceptance probability is 0
        return 0;
    }
  }

  double geneNewLogProb = 0.0;
  if (genes[gene]->hasProb(newTop)) {
    geneNewLogProb = genes[gene]->getLogProb(newTop);
  }

  if(log(rand.runif()) < newLogProb - oldLogProb + genes[gene]->getLogProb(oldTop) - geneNewLogProb) {
    // make changes to newTop
    for(int i=0;i<set.size();i++)
      tops[set[i]] = newTop;
    counts[newTop] = counts[oldTop];
    counts[oldTop] = 0;
    j=0;
    while(indices[j]!=oldTop)
      j++;
    indices.erase(indices.begin()+j,indices.begin()+j+1);
    indices.push_back(newTop);
    sort(indices.begin(),indices.end());
    return 1;
  }
  else
    return 0;
}

int State::update(Rand& rand) {
  int sum=0;
  for(int i=0;i<genes.size();i++)
    sum += updateOne(i,rand);
  return sum;
}

void State::print(ostream& f) {
  f << "State:" << endl;
  f << "  alpha = " << alpha << endl;
  f << "  logAlpha = " << logAlpha << endl;
  f << "  numGroups = " << numGroups << endl;
  f << "  length of indices = " << indices.size() << endl; // numGroups
  f << "  logPriorProb = " << logPriorProb << endl;
  f << "  logPosteriorProbProduct = " << logPosteriorProbProduct << endl;
  f << "  logH = " << getLogH() << endl;
  f << "  counts:" << endl;
  f << "    Top Count" << endl;
  for(int i=0;i<numGroups;i++)
    f << setw(7) << indices[i] << setw(6) << counts[indices[i]] << endl;
  f << "  tops:" << endl;
  f << "   gene   top" << endl;
  for(int i=0;i<genes.size();i++)
    f << setw(7) << i << setw(6) << tops[i] << endl;

}

void State::updateSplits(vector<int>& totals,vector<vector<int> >& topSplitsIndex, int numSplits) {

  int numSplitsPerTree = topSplitsIndex[0].size();
  vector<int> counts(numSplits);
  for(int i=0;i<numSplits;i++)
    counts[i] = 0;
  for(int i=0;i<genes.size();i++)
    for(int j=0;j<numSplitsPerTree;j++) {
      counts[topSplitsIndex[tops[i]][j]]++;
    }
  for(int i=0;i<numSplits;i++)
    totals[(i*numSplits)+counts[i]]++;
}

void State::updatePairCounts(vector<int>& counts)
{
  // slow way to do this.  maybe we should have State maintain a vector of groups....
  // only fill in upper triangle (i <= j)
  for(int i=0;i<genes.size();i++)
    for(int j=i;j<genes.size();j++)
      if( tops[i] == tops[j] )
	counts[(i*genes.size())+j]++;
}

Node* Node::getNeighbor(int i) const { return edges[i]->getOtherNode(this); }

TaxonSet Node::setAllTaxa(Node* parent,TaxonSet all) {
  int ep;
  TaxonSet allMx = all;
  Edge* e;
  if(leaf) {
    allMx = all.excludeTaxon(number);
    taxa[0] = allMx;
    e = edges[0];
  }
  else {
    for(int i=0;i<3;i++) {
      Node* n = getNeighbor(i);
      if(n != parent) {
	TaxonSet y = n->setAllTaxa(this,all);
	setTaxa(i, all - y);
	allMx = allMx & y;
      }
      else {
	e = edges[i];
	ep = i;
      }
    }
  }
  TaxonSet s = (allMx[0] == 0 ? all - allMx : allMx);
  if(!leaf)
    setTaxa(ep,allMx);
  e->setSplit(s);
  return allMx;
}

void Tree::setAllTaxa() {
  for(int i=0;i<3;i++)
    root->setTaxa(i,all - root->getNeighbor(i)->setAllTaxa(root,all));
}

void Node::print(ostream& f) const
{
  f << setw(3) << number << ":";
  TaxonSet t = taxa[0];
  f << " taxa = (" << t;
  if(!leaf)
    for(int i=1;i<3;i++) {
        t = taxa[i];
      f << "," << t;
    }

  f << "), ";
  for(int i=0;i<3;i++) {
    if(i>0 && !leaf)
      f << ",";
    if(i==0 || !leaf)
      f << " e[" << edges[i]->getNumber() << "] -> " << "n[" << getNeighbor(i)->getNumber() << "]";
  }
  f << endl;
}

void Edge::print(ostream& f,int numTaxa) const
{
  f << setw(3) << number << ":";
  f << " n[" << nodes[0]->getNumber() << "] <-> " << "n[" << nodes[1]->getNumber() << "], split = ";
  split.print(f);
  f << endl;
}

void Tree::print(ostream& f) const
{
  f << top << endl;
  f << "numTaxa = " << numTaxa << ", numNodes = " << numNodes << ", numEdges = " << numEdges << endl;
  f << "Nodes:" << endl;
  for(int i=0;i<numNodes;i++)
    nodes[i]->print(f);
  f << "Edges:" << endl;
  for(int i=0;i<numEdges;i++)
    edges[i]->print(f,numTaxa);
}

void Tree::construct(string& top) {
    // top is unrooted tree topology of form (X, Y, Z) X, Y, Z can be treated as rooted subtree topologies.
    istringstream topStr(top);
    int currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
    char ch;
    topStr >> ch;
    Node *n1 = connectInt(topStr, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    Node *n2 = connectInt(topStr, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    Node *n3 = connectInt(topStr, currentInternalNode, currentLeafEdge, currentInternalEdge);
    root = nodes[currentInternalNode++];

    Edge *e1, *e2, *e3;
    if (n1->isLeaf())
        e1 = edges[currentLeafEdge++];
    else
        e1 = edges[currentInternalEdge++];

    if (n2->isLeaf())
        e2 = edges[currentLeafEdge++];
    else
        e2 = edges[currentInternalEdge++];

    if (n3->isLeaf())
        e3 = edges[currentLeafEdge++];
    else
        e3 = edges[currentInternalEdge++];

    root->setEdge(0, e1);
    root->setEdge(1, e2);
    root->setEdge(2, e3);

    n1->setEdge(0, e1);
    n2->setEdge(0, e2);
    n3->setEdge(0, e3);

    e1->setNode(0, root);
    e1->setNode(1, n1);
    e2->setNode(0, root);
    e2->setNode(1, n2);
    e3->setNode(0, root);
    e3->setNode(1, n3);
}

Node* Tree::connectInt(istream& topStr,int& currentInternalNode,int& currentLeafEdge,int& currentInternalEdge)
{
  int n;
  char ch = topStr.peek();
  if(ch!='(') {
    topStr >> n;
    return nodes[n-1];
  }
  else {
    topStr >> ch;
    Node *node1,*node2;
    node1 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=',') {
      cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl;
      exit(1);
    }
    node2 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=')') {
      cerr << "Error: Cannot parse topology string, expected ')', found " << ch << endl;
      exit(1);
    }
    return connectThreeNodes(node1,node2,nodes[currentInternalNode++],currentLeafEdge,currentInternalEdge);
  }
}

Node *Tree::connectThreeNodes(Node* node1,Node* node2,Node* node3,int& currentLeafEdge,int& currentInternalEdge)
{
  // node3 must be an internal node.
  Edge *edge1,*edge2;
  if(node1->isLeaf())
    edge1 = edges[currentLeafEdge++];
  else
    edge1 = edges[currentInternalEdge++];
  if(node2->isLeaf())
    edge2 = edges[currentLeafEdge++];
  else
    edge2 = edges[currentInternalEdge++];
  node1->setEdge(0,edge1);
  node2->setEdge(0,edge2);
  node3->setEdge(1,edge1);
  node3->setEdge(2,edge2);

  edge1->setNode(0,node1);
  edge1->setNode(1,node3);
  edge2->setNode(0,node2);
  edge2->setNode(1,node3);
  return node3;
}

void Tree::getSplits(SplitSet& s) {
    s.resize(numSplits);
    s.setAll(all);
    s.setNumTaxa(numTaxa);
    for(int i=0,j=numTaxa;i<numSplits;i++)
      s.setSplit(i,edges[j++]->getSplit());
  }






// assumes that numChains > 1 and useIndependencePrior == false
void MPI_mcmcmc(vector<State*>& states, vector<double>& alphas, 
        Rand* rand, vector<int>& accepts, vector<int>& proposals, 
        int rank, const vector<int>& runs, const vector<int>& chains , 
        int numChains , int runNum, vector<int>& rank_idx, MPI_Comm comm)
{
    
  //rand is swap_rand- all ranks in sync
  //Because all ranks received same swap seeds
  //!!!ALL RANKS MUST CALL MCMCMC ONCE FOR EACH RUN THEY HAVE CHAINS FOR!!!
  
  int a = (int) ( (numChains-1)*rand->runif() );
  //cout << rank << " has chosen chain "<<a<<" for exchange with chain "<<a+1<<"(run "<<runNum<<")\n";
  int b = a+1;
  int local_a = a+(numChains*runNum);
  int local_b = a+(numChains*runNum)+1;
  
  int send_rank, recv_rank; 

  //Query "global" tables to capture ranks involved in send/recv
  send_rank = rank_idx[local_a]; //Is this correct?? 
  recv_rank = rank_idx[local_b];  

  double aLogPriorProb, aLogPostProb, bLogPriorProb, bLogPostProb;
  //All ranks sample an acceptance probability, advance RNG in sync
  double logAccept = log(rand->runif()); 
  bool accept; 
  
  if ((rank == send_rank) || (rank == recv_rank)){
    proposals[local_a]++;
    if (send_rank != recv_rank){
      if (rank == send_rank){
		//Exchange ALPHA between send & recv ranks
		double myAlpha = states[local_a]->getAlpha();
		double otherAlpha;
		//cout << rank << " sending "<<myAlpha<<" to rank "<< recv_rank<<endl;
		MPI_Sendrecv(&myAlpha, 1, MPI_DOUBLE, recv_rank, 3, &otherAlpha, 1, MPI_DOUBLE, recv_rank, 3, comm, MPI_STATUS_IGNORE);
        //cout <<flush; 
        //cout << "Rank "<<send_rank<<" will exchange chain "
        //     << a << " (index "<<local_a<<") with chain "<<b<<
        //        " (index "<<local_b<<") at rank "<<recv_rank<<endl; 
        //Capture probabilities to send
        aLogPriorProb = states[local_a]->getLogPriorProb(); 
        //cout << aLogPriorProb<<"\n";
        states[local_a]->setAlpha(otherAlpha);
        aLogPostProb = states[local_a]->getLogPriorProb(); 
        //cout << aLogPostProb<<"\n";
        //Send logPriorProb and logPostProb to recv_rank
        MPI_Send(&aLogPriorProb, 1, MPI_DOUBLE, recv_rank, 0, comm);
        MPI_Send(&aLogPostProb, 1, MPI_DOUBLE, recv_rank, 1, comm); 
        //recv_rank will calculate and return if swap is accepted
        MPI_Recv(&accept, 1, MPI::BOOL, recv_rank, 2, comm, MPI_STATUS_IGNORE);
        
        if (accept == true){ 
            accepts[local_a]++; 
        }else{
            states[local_a]->setAlpha(myAlpha);
        }
      }else if (rank == recv_rank){ 
		//Exchange ALPHA between send & recv ranks
		double myAlpha = states[local_b]->getAlpha();
		double otherAlpha;
		//cout << rank << " sending "<<myAlpha<<" to rank "<< send_rank<<endl;
		MPI_Sendrecv(&myAlpha, 1, MPI_DOUBLE, send_rank, 3, &otherAlpha, 1, MPI_DOUBLE, send_rank, 3, comm, MPI_STATUS_IGNORE);
        //cout << "Rank "<<rank<<" is receiving from rank "<<send_rank<<endl;
        //Calculate log probabilities
        bLogPriorProb = states[local_b]->getLogPriorProb(); 
        //cout << bLogPriorProb<<"\n";
        states[local_b]->setAlpha(otherAlpha);
        bLogPostProb = states[local_b]->getLogPriorProb(); 
        //cout << bLogPostProb<<"\n";
        //Send logPriorProb and logPostProb from send_rank
        MPI_Recv(&aLogPriorProb, 1, MPI_DOUBLE, send_rank, 0, comm, MPI_STATUS_IGNORE); 
        MPI_Recv(&aLogPostProb, 1, MPI_DOUBLE, send_rank, 1, comm, MPI_STATUS_IGNORE); 
        //Receiving process calculates if proposed change is accepted
        double logProbDiff = ((-aLogPriorProb - bLogPriorProb)+(aLogPostProb + bLogPostProb));
        if (logAccept < logProbDiff){
            accept = true;
            accepts[local_a]++; 
        }else{
            accept = false; 
            states[local_b]->setAlpha(myAlpha);
        }
        //Return bool to send_rank
        MPI_Send(&accept, 1, MPI::BOOL, send_rank, 2, comm);
      }
    }else if ((send_rank == recv_rank)&&(rank == send_rank)){ 
	  double aAlpha = states[local_a]->getAlpha();
	  double bAlpha = states[local_b]->getAlpha();
      //cout <<flush; 
      //cout << "Rank "<<send_rank<<" will exchange chain "
      //     << a << " (index "<<local_a<<") with chain "<<b
      //     << " (index "<<local_b<<") at rank "<<recv_rank<<endl; 
      aLogPriorProb = states[local_a]->getLogPriorProb(); 
      bLogPriorProb = states[local_b]->getLogPriorProb(); 
      states[local_a]->setAlpha(bAlpha);
      states[local_b]->setAlpha(aAlpha);
      bLogPostProb = states[local_b]->getLogPriorProb();
      aLogPostProb = states[local_a]->getLogPriorProb(); 
      //Receiving process calculates if proposed change is accepted
      double logProbDiff = ((-aLogPriorProb - bLogPriorProb)+(aLogPostProb + bLogPostProb));
      if (logAccept < logProbDiff){ //accept
          accept = true; 
          accepts[local_a]++; 
      }else{ //reject
          accept = false; 
          states[local_a]->setAlpha(aAlpha);
          states[local_b]->setAlpha(bAlpha);
      }
    }
  }
}





void usage(int rank, Defaults defaults)
{
  if (rank == 0){
  cerr << "Usage: bucky [options] <input files>" << endl << endl;
  cerr << "  Options:" << endl;
  cerr << "  Parameter                      | Usage                      | Default Value" << endl;
  cerr << "  -------------------------------------------------------------------" << endl;
  cerr << "  alpha                          | -a number                  | " << defaults.getAlpha() << endl;
  cerr << "  # of runs                      | -k integer                 | " << defaults.getNumRuns() << endl;
  cerr << "  # of MCMC updates              | -n integer                 | " << defaults.getNumUpdates() << endl;
  cerr << "  # of chains                    | -c integer                 | " << defaults.getNumChains() << endl;
  cerr << "  MCMCMC Rate                    | -r integer                 | " << defaults.getMCMCMCRate() << endl;
  cerr << "  alpha multiplier               | -m number                  | " << defaults.getAlphaMultiplier() << endl;
  cerr << "  subsample rate                 | -s integer                 | " << defaults.getSubsampleRate() << endl;
  cerr << "  output root file name          | -o name                    | " << defaults.getRootFileName() << endl;
  cerr << "  input file list file           | -i filename                | " << defaults.getInputListFileName() << endl;
  cerr << "  random seed 1                  | -s1 integer                | " << defaults.getSeed1() << endl;
  cerr << "  random seed 2                  | -s2 integer                | " << defaults.getSeed2() << endl;
  cerr << "  CF cutoff for display          | -cf number                 | " << defaults.getSwCFcutoff() << endl;
  cerr << "  create sample file             | --create-sample-file       | " << (defaults.getCreateSampleFile() == true ? "true" : "false") << endl;
  cerr << "  create joint file              | --create-joint-file        | " << (defaults.getCreateJointFile() == true ? "true" : "false") << endl;
  cerr << "  create single file             | --create-single-file       | " << (defaults.getCreateSingleFile() == true ? "true" : "false") << endl;
  cerr << "  use independence prior         | --use-independence-prior   | " << (defaults.getUseIndependencePrior() == true ? "true" : "false") << endl;
  cerr << "  calculate pairs                | --calculate-pairs          | " << (defaults.getCalculatePairs() == true ? "true" : "false") << endl;
  cerr << "  taxon set                      | -p prune-file              | common taxa" << endl;
  cerr << "  skip genes with fewer taxa     | -sg                        | false" << endl;
  cerr << "  use update groups              | --use-update-groups        | " << (defaults.getUseUpdateGroups() == true ? "true" : "false") << endl;
  cerr << "  use update groups              | --do-not-use-update-groups | " << endl;
  cerr << "  do not build population tree?  | --no-population-tree       | false"  << endl;
  cerr << "  grid-size for genomewide CF    | --genomewide-grid-size     | " << defaults.getNumGenomewideGrid() << endl;
  cerr << "  help                           | -h OR --help               |" << endl;
  cerr << "  version                        | --version                  |" << endl;
  cerr << "  -------------------------------------------------------------------" << endl << endl;
  }
  
  MPI_Finalize();
  exit(0);
}

void showParameters(ostream& f,FileNames& fn,Defaults defaults,ModelParameters& mp,RunParameters& rp)
{
  f << "  Parameter              | Usage                    | Default Value | Value Used" << endl;
  f << "  ------------------------------------------------------------------------------" << endl;
  f << "  alpha                  | -a number                | " << left << setw(14) << defaults.getAlpha()                                             << "| " << mp.getAlpha() << endl;
  f << "  # of runs              | -k integer               | " << left << setw(14) << defaults.getNumRuns()                                           << "| " << rp.getNumRuns() << endl;
  f << "  # of MCMC updates      | -n integer               | " << left << setw(14) << defaults.getNumUpdates()                                        << "| " << rp.getNumUpdates() << endl;
  f << "  # of chains            | -c integer               | " << left << setw(14) << defaults.getNumChains()                                         << "| " << rp.getNumChains() << endl;
  f << "  MCMCMC Rate            | -r integer               | " << left << setw(14) << defaults.getMCMCMCRate()                                        << "| " << rp.getMCMCMCRate() << endl;
  f << "  alpha multiplier       | -m number                | " << left << setw(14) << defaults.getAlphaMultiplier()                                   << "| " << rp.getAlphaMultiplier() << endl;
  f << "  subsample rate         | -s integer               | " << left << setw(14) << defaults.getSubsampleRate()                                     << "| " << rp.getSubsampleRate() << endl;
  f << "  output root file name  | -o name                  | " << left << setw(14) << defaults.getRootFileName()                                      << "| " << fn.getRootFileName() << endl;
  f << "  input file list file   | -i filename              | " << left << setw(14) << defaults.getInputListFileName()                                 << "| " << fn.getInputListFileName() << endl;
  f << "  random seed 1          | -s1 integer              | " << left << setw(14) << defaults.getSeed1()                                             << "| " << rp.getSeed1() << endl;
  f << "  random seed 2          | -s2 integer              | " << left << setw(14) << defaults.getSeed2()                                             << "| " << rp.getSeed2() << endl;
  f << "  CF cutoff for display  | -cf number               | " << left << setw(14) << defaults.getSwCFcutoff()                                        << "| " << rp.getSwCFcutoff() << endl;
  f << "  create sample file     | --create-sample-file     | " << left << setw(14) << (defaults.getCreateSampleFile() == true ? "true" : "false")     << "| " << (rp.getCreateSampleFile() ? "true" : "false") << endl;
  f << "  create joint file      | --create-joint-file      | " << left << setw(14) << (defaults.getCreateJointFile() == true ? "true" : "false")      << "| " << (rp.getCreateJointFile() ? "true" : "false") << endl;
  f << "  create single file     | --create-single-file     | " << left << setw(14) << (defaults.getCreateSingleFile() == true ? "true" : "false")     << "| " << (rp.getCreateSingleFile() ? "true" : "false") << endl;
  f << "  use independence prior | --use-independence-prior | " << left << setw(14) << (defaults.getUseIndependencePrior() == true ? "true" : "false") << "| " << (mp.getUseIndependencePrior() ? "true" : "false") << endl;
  f << "  calculate pairs        | --calculate-pairs        | " << left << setw(14) << (defaults.getCalculatePairs() == true ? "true" : "false")       << "| " << (rp.getCalculatePairs() ? "true" : "false" ) << endl;
  f << "  use update groups      | --use-update-groups      | " << left << setw(14) << (defaults.getUseUpdateGroups() == true ? "true" : "false")      << "| " << (rp.getUseUpdateGroups() ? "true" : "false")<< endl;
  f << "  File with prune list   | -p pruneFile             | " << left << setw(14) << ""                                                              << "| " << rp.getPruneFile() << endl;
  f << "  skip genes             | -sg                      | " << left << setw(14) <<"false" << "| " << (rp.getPruneGene() ? "true" : "false")<< endl;
  f << "  no population tree?    | --no-population-tree     | " << left << setw(14) <<"false" << "| " << (rp.shouldBuildPopulationTree() ? "false" : "true") << endl;
  f << "  genomewide CF grid size| --genomewide-grid-size   | " << left << setw(14) << defaults.getNumGenomewideGrid()                                 << "| " << rp.getNumGenomewideGrid() << endl;
  f << "  ------------------------------------------------------------------------------" << endl;
}

void intro(ostream& f) {
  f << endl << endl; 
  f << "Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)" << endl;
  f << endl<< PRE << endl << PRE2 << endl << endl<< PRE3 << endl ;
  f << "BUCKy version " << VERSION << ", " << DATE << endl;
  f << COPYRIGHT << endl << endl;
  f << "This is free software; see the source for copying conditions.  There is NO" << endl;
  f << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl << endl;
}

int readArguments(int rank, int argc, char *argv[],FileNames& fn,ModelParameters& mp,RunParameters& rp,Defaults& defaults)
{
  int k=1;
  bool done = (argc>1 ? false : true);
  while(!done && k<argc) {
    string flag=argv[k];
    if(flag=="--version")
      exit(0);
    if(flag=="-h" || flag=="--help")
      usage(rank, defaults);
    else if(flag=="-a") {
      double alpha;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> alpha) )
	usage(rank, defaults);
      if(alpha <= 0.0)
	cerr << "Warning: parameter alpha must be positive.  Ignoring argument -a " << alpha << "." << endl;
      else
	mp.setAlpha(alpha);
      k++;
    }
    else if(flag=="-n") {
      unsigned int numUpdates;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numUpdates) )
	usage(rank, defaults);
      rp.setNumUpdates(numUpdates);
      k++;
    }
    else if(flag=="-k") {
      unsigned int numRuns;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numRuns) )
	usage(rank,defaults);
      if (numRuns==0){
	cerr << "Warning: parameter number-of-runs must be at least one. Setting to default: -k 2." << endl;
	numRuns=2;
      }
      rp.setNumRuns(numRuns);
      k++;
    }
    else if(flag=="-c") {
      unsigned int numChains;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numChains) )
	usage(rank,defaults);
      if(numChains==0) {
	cerr << "Warning: parameter number-of-chains must be at least one. Setting to default -c 1." << endl;
	numChains=1;
      }
      if(mp.getUseIndependencePrior()){
	cerr << "Warning: parameter number-of-chains can only be used with non-independence prior. Setting to default -c 1." << endl;
	numChains=1;
      }
      rp.setNumChains(numChains);
      k++;
    }
    else if(flag=="-m") {
      double alphaMultiplier;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> alphaMultiplier) )
	usage(rank,defaults);
      if(alphaMultiplier<=0)
	cerr << "Warning: parameter alpha-multiplier must be positive. Ignoring argument -m " << alphaMultiplier << "." << endl;
      else
	rp.setAlphaMultiplier(alphaMultiplier);
      k++;
    }
    else if(flag=="-r") {
      unsigned int mcmcmcRate;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> mcmcmcRate) )
	usage(rank,defaults);
      if(mcmcmcRate==0) {
	cerr << "Warning: parameter MCMCMC-rate must be at least one. Setting to default -r 1." << endl;
	mcmcmcRate=1;
      }
      rp.setMCMCMCRate(mcmcmcRate);
      k++;
    }
    else if(flag=="-s") {
      unsigned int subsampleRate;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> subsampleRate) )
	usage(rank,defaults);
      if(subsampleRate==0) {
	cerr << "Warning: parameter subsample-rate must be at least one. Setting to default -s 1." << endl;
	subsampleRate=1;
      }
      rp.setSubsampleRate(subsampleRate);
      k++;
    }
    else if(flag=="-o") {
      string root = argv[++k];
      fn.setFileNames(root);
      k++;
    }
    else if(flag=="-i") {
      string filename = argv[++k];
      fn.setInputListFileName(filename);
      k++;
    }
    else if(flag=="-s1") {
      unsigned int seed1;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> seed1) )
	usage(rank,defaults);
      if(seed1==0)
	cerr << "Warning: parameter seed1 must be at least one. Ignorning command -s1 " << seed1 << "." << endl;
      else
	rp.setSeed1(seed1);
      k++;
    }
    else if(flag=="-s2") {
      unsigned int seed2;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> seed2) )
	usage(rank,defaults);
      if(seed2==0)
	cerr << "Warning: parameter seed2 must be at least one. Ignorning command -s2 " << seed2 << "." << endl;
      else
	rp.setSeed2(seed2);
      k++;
    }
    else if(flag=="-cf") {
      double cutoff;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> cutoff) )
	usage(rank,defaults);
      rp.setSwCFcutoff(cutoff);
      k++;
    }
    else if(flag=="--create-sample-file") {
      rp.setCreateSampleFile(true);
      k++;
    }
    else if(flag=="--create-joint-file") {
      rp.setCreateJointFile(true);
      k++;
    }
    else if(flag=="--create-single-file") {
      rp.setCreateSingleFile(true);
      k++;
    }
    else if(flag=="--use-independence-prior") {
      mp.setUseIndependencePrior(true);
      // to make sure alpha=Inf and is not used...
      mp.setAlpha(numeric_limits<double>::infinity());
      cerr << "Using independence prior: setting alpha=" << mp.getAlpha() << endl;
      if(rp.getNumChains()>1){
	cerr << "Warning: Independence prior can only be used with 1 chain. Setting number of chains to default -c 1." << endl;
	rp.setNumChains(1);
      }
      k++;
    }
    else if(flag=="--calculate-pairs") {
      rp.setCalculatePairs(true);
      k++;
    }
    else if(flag=="--use-update-groups") {
      rp.setUseUpdateGroups(true);
      k++;
    }
    else if(flag=="--do-not-use-update-groups") {
      rp.setUseUpdateGroups(false);
      k++;
    }
    else if (flag =="-p") {
      rp.setPruneFile(argv[++k]);
      k++;
    }
    else if (flag =="-sg") {
      rp.setPruneGene(true);
      k++;
    }
    else if (flag =="--genomewide-grid-size") {
      unsigned int ngrid;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> ngrid) )
	usage(rank,defaults);
      if(ngrid==0)
	cerr << "Warning: parameter grid-size must be at least one. Ignoring option genomewide-grid-size " << ngrid << "." << endl;
      else
	rp.setNumGenomewideGrid(ngrid);
      k++;
    }
    else if (flag =="--opt-space") {
	  if (rank==0)
	    cerr << "Warning: --opt-space not available yet in MPI-BUCKy, ignoring flag!" << endl;
      rp.setOptSpace(false);
      k++;
    }
    else if (flag =="--no-population-tree") {
      rp.setBuildPopulationTree(false);
      k++;
    }
    else
      done = true;
  }

  if(fn.getInputListFileName().empty() && argc-k<1)
    usage(rank,defaults);
  return k;
}

static const string WHITESPACE_CHARS(" \t\n\r\f\v");

// Returns the string defined by S, with all characters in CHARS
// stripped from its righthand side
string stripRight(const string& s,
                  const string& chars = WHITESPACE_CHARS) {
  string::size_type pos = s.find_last_not_of(chars);
  return (pos == string::npos ? "" : s.substr(0, pos + 1));
}

void readInputFileList(string inputListFileName, vector<string>& inputFiles) {
  ifstream f(inputListFileName.c_str());
  string filename;
  while (getline(f, filename)) {
    inputFiles.push_back(stripRight(filename));
  }
}

void MPI_readInputFiles(int rank, vector<string>& inputFiles,TGM* &tgm, vector<string>& translateTable,
                    int& max, ostream& fout, vector<vector<int> >& taxid, string prunefile, bool shouldPruneGene)
{
  bool changed = false;
  map<string, int> translateMap;
  getTaxaSubset(rank, inputFiles, translateMap, translateTable, prunefile, shouldPruneGene, changed);
  if (inputFiles.size() == 0) {
      if (rank ==0)
        cerr << "All genes skipped by the pruning step, Exiting..." << endl;
      MPI_Finalize(); 
      exit(0);
  }
  if (translateTable.size() < 4) {
    cerr << "Too few taxa (" << translateTable.size() << ") left after pruning process, Taxa common to all genes:" << endl;
    for (map<string,int>::iterator itr = translateMap.begin(); itr != translateMap.end(); itr++) {
      cerr << itr->first << "  ";
    }
    cerr << endl;
    exit(0);
  }


  int part = inputFiles.size() / 50;

  taxid.resize(inputFiles.size());
  for (size_t i=0;i<inputFiles.size();i++){
    readFile(inputFiles[i],i,tgm,max,taxid[i],translateMap, changed);
  }

}

void Defaults::print(ostream& f) {
  f << "BUCKy default values:" << endl;
  f << "-a alpha-value = " << alpha << endl;
  f << "-k number-of-chains = " << numRuns << endl;
  f << "-n number-of-MCMC-updates = " << numUpdates << endl;
  f << "-c number-of-chains = " << numChains << endl;
  f << "-r MCMCMC-rate = " << mcmcmcRate << endl;
  f << "-m alpha-multiplier = " << alphaMultiplier << endl;
  f << "-s subsample-rate = " << subsampleRate << endl;
  f << "-f output-file-root = " << rootFileName << endl;
  f << "-s1 seed1 = " << seed1 << endl;
  f << "-s2 seed2 = " << seed2 << endl;
  f << "--create-sample-file = " << (createSampleFile ? "true" : "false") << endl;
  f << "--create-joint-file = " << (createJointFile ? "true" : "false") << endl;
  f << "--create-single-file = " << (createSingleFile ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (useIndependencePrior ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (calculatePairs ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (useUpdateGroups ? "true" : "false") << endl;
}

void printParameters(ostream& f,FileNames& fn,ModelParameters& mp,RunParameters& rp) {
//this function seems to be obsolete and not used.
  f << "Parameters used:" << endl;
  f << "-a alpha-value =" << mp.getAlpha() << endl;
  f << "-k number-of-chains = " << rp.getNumRuns() << endl;
  f << "-n number-of-MCMC-updates = " << rp.getNumUpdates() << endl;
  f << "-c number-of-chains = " << rp.getNumChains() << endl;
  f << "-r MCMCMC-rate = " << rp.getMCMCMCRate() << endl;
  f << "-m alpha-multiplier = " << rp.getAlphaMultiplier() << endl;
  f << "-s subsample-rate = " << rp.getSubsampleRate() << endl;
  f << "-f output-file-root = " << fn.getRootFileName() << endl;
  f << "-s1 seed1 = " << rp.getSeed1() << endl;
  f << "-s2 seed2 = " << rp.getSeed2() << endl;
  f << "--create-sample-file = " << (rp.getCreateSampleFile() ? "true" : "false") << endl;
  f << "--create-joint-file = " << (rp.getCreateJointFile() ? "true" : "false") << endl;
  f << "--create-single-file = " << (rp.getCreateSingleFile() ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (mp.getUseIndependencePrior() ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (rp.getCalculatePairs() ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (rp.getUseUpdateGroups() ? "true" : "false") << endl;
}

// Fix error in finding concordance tree
// Two sets are compatible if one is the subset of the other or a subset of the complement of the other or vice versa
bool isCompatible(unsigned int x,unsigned int y,unsigned int all) {
  if( x == (x & y) )
    return true;
  if( y == (x & y) )
    return true;
  unsigned int xc = x^all;
  if( xc == (xc & y) )
    return true;
  if( y == (xc & y) )
    return true;
  return false;
}
void GenomewideDistribution::updateSamplewide(vector<int> &splitGeneCount, unsigned int numUpdates) {
  vector<double> splitGenePP(splitGeneCount.size());
  for (int j=0;j<splitGeneCount.size();j++)
    splitGenePP[j] = (double)(splitGeneCount[j])/numUpdates;
  this->updateSamplewide(splitGenePP);
}

void GenomewideDistribution::updateSamplewide(vector<double> &splitGenePP){
  if (splitGenePP.size() != samplewide.size()){
      cerr << "InternalError: wrong number of genes in GenomewideDistribution::updateSamplewide;" << endl;
      exit(1);
  }
  int ngenes=samplewide.size()-1;
  vector<double> csum(samplewide.size());
  for(int j=0;j<splitGenePP.size();j++) {
    samplewide[j] = splitGenePP[j];
    if (j>0) csum[j] =  csum[j-1]+samplewide[j];
    else csum[j] = samplewide[j];
  }
  int j=0;
  while(csum[j] < 0.005) j++;
  samplewideCredibilityInterval[0]= ((double) j)/ngenes;
  while(csum[j] < 0.025) j++;
  samplewideCredibilityInterval[1]= ((double) j)/ngenes;
  while(csum[j] < 0.05) j++;
  samplewideCredibilityInterval[2]= ((double) j)/ngenes;
  while(csum[j] < 0.95) j++;
  samplewideCredibilityInterval[3]= ((double) j)/ngenes;
  while(csum[j] < 0.975) j++;
  samplewideCredibilityInterval[4]= ((double) j)/ngenes;
  while(csum[j] < 0.995 && j<ngenes) j++;
  samplewideCredibilityInterval[5]= ((double) j)/ngenes;
}

void GenomewideDistribution::updateConvolutionWeight(double alpha){
  int ngenes = samplewide.size()-1;
  int ngrid  = genomewide.size()-1;
  double priorQ = 1 - priorProbability;
  vector<vector<double> > logweight;
  logweight.resize(ngrid+1);
  for (int i=0; i<ngrid+1; i++) logweight[i].resize(ngenes+1);
  vector<double> logalphap(ngenes+1);
  for (int j=0; j<ngenes+1; j++) logalphap[j] = log(alpha*priorProbability + j);
  vector<double> logalphaq(ngenes+1);
  for (int j=0; j<ngenes+1; j++) logalphaq[j] = log(alpha*priorQ + j);
  vector<double> xinside(ngrid-1);
  vector<double> u(ngrid-1);
  for (int i=1; i<ngrid; i++){
    xinside[i-1] = ((double) (i))/ngrid;
    u[i-1] = log(xinside[i-1])-log(1-xinside[i-1]);
  }
  for (int i=1; i<ngrid; i++){
    logweight[i][1]=alpha*priorProbability*log(xinside[i-1])+(alpha*priorQ+ngenes-2)*log(1-xinside[i-1]);
    logweight[i][0]   = logweight[i][1] - u[i-1] + logalphap[0] - logalphaq[ngenes-1];
  }
  for (int j=2; j<ngenes+1; j++){
    for (int i=1; i<ngrid; i++)
      logweight[i][j] = logweight[i][j-1]+u[i-1] - logalphap[j-1]+logalphaq[ngenes-j];
  }
  // Next: calculate the normalizing constant Z=Beta(alpha*priorP+1, alpha*priorQ+ngenes-1).
  // logZ ~ log of sum(exp(logweight[,j]))/ngrid, for any j, avoiding the need for
  // a library for the Beta or the Gamma function. Checked with R: very accurate
  // approximation as soon as ngrid > 100 and j well chosen (so that weight[,j] not too skewed).
  int bestj = (int) ( ((double)ngenes)/2 + alpha*(0.5-priorProbability));
  if (bestj>ngenes-1) bestj=ngenes-1;
  else if (bestj<1) bestj=1;
  double logZ = 0;
  for (int i=1; i<ngrid; i++) logZ += exp(logweight[i][bestj]);
  logZ = log( logZ/ngrid );
  for (int i=1; i<ngrid; i++)
    for (int j=0; j<ngenes+1; j++)
      convolutionWeight[i][j] = exp(logweight[i][j] -logZ);
  for (int j=0; j<ngenes+1; j++){
    convolutionWeight[0][j] = 0.0;
    convolutionWeight[ngrid][j] = 0.0;
  }
  convolutionWeight[0][0]=(double) ngrid;
  for (int i=1; i<ngrid+1; i++)
    convolutionWeight[0][0] -= convolutionWeight[i][0];
  convolutionWeight[ngrid][ngenes]=(double) ngrid;
  for (int i=0; i<ngrid; i++)
    convolutionWeight[ngrid][ngenes] -= convolutionWeight[i][ngenes];
}

void GenomewideDistribution::updateGenomewide(double alpha) {
  int ngenes = samplewide.size()-1;
  int ngrid  = genomewide.size()-1;
  genomewidePosteriorMean = alpha/(alpha+ngenes)*priorProbability + ngenes/(alpha+ngenes)*samplewidePosteriorMean;

  vector<double> csum(genomewide.size());
  for(int i=0;i<ngrid+1;i++) {
    // matrix multiplication: genomewide = convolutionWeights * samplewide (density)
    double sumprod=0;
    for (int j=0; j<ngenes+1; j++)
      sumprod += convolutionWeight[i][j]*samplewide[j]/ngrid;
    genomewide[i]=sumprod;
    if (i>0) csum[i] += csum[i-1]+genomewide[i];
    else csum[i] = genomewide[i];
  }

  int i=0;
  while(csum[i] < 0.005 && i<ngrid) i++;
  /* the condition i<ngrid is not needed in theory, but in practice it is because
     cumsum[ngrid] may be < 1.0, due to the approximation in Z (see above) */
  genomewideCredibilityInterval[0]= ((double) i)/ngrid;
  while(csum[i] < 0.025 && i<ngrid) i++;
  genomewideCredibilityInterval[1]= ((double) i)/ngrid;
  while(csum[i] < 0.05 && i<ngrid) i++;
  genomewideCredibilityInterval[2]= ((double) i)/ngrid;
  while(csum[i] < 0.95 && i<ngrid) i++;
  genomewideCredibilityInterval[3]= ((double) i)/ngrid;
  while(csum[i] < 0.975 && i<ngrid) i++;
  genomewideCredibilityInterval[4]= ((double) i)/ngrid;
  while(csum[i] < 0.995 && i<ngrid) i++;
  genomewideCredibilityInterval[5]= ((double) i)/ngrid;
}


// Currently one function to write all output
// Should have separate functions for each type of output
void writeOutput(ostream& fout,FileNames& fileNames,int max,int numTrees,int numTaxa,vector<string> topologies,
		 int numGenes,RunParameters& rp,ModelParameters& mp,TGMTable *newTable,
		 vector<vector<int> >& clusterCount, vector<TaxonSet>& splits,  vector<vector<vector<int> > >& splitsGeneMatrix,
		 vector<vector<int> >& pairCounts,   vector<Gene*>& genes, vector<double>& alphas,
		 vector<vector<int> >& mcmcmcAccepts,vector<vector<int> >& mcmcmcProposals, vector<string>& translateTable)
{
  // .joint
  if(rp.getCreateJointFile()) {
    cout << "Writing joint posterior table to file " << fileNames.getJointFile() << "...." << flush;
    ofstream jointStr(fileNames.getJointFile().c_str());
    if(jointStr.fail()) {
      cerr <<"Error: Cannot open file " << fileNames.getJointFile() << "." << endl;
    }
    jointStr.setf(ios::fixed, ios::floatfield);
    jointStr.setf(ios::showpoint);
    for(int i=0;i<numTrees;i++) {
      jointStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
      for(int j=0;j<numGenes;j++) {
          jointStr << " " << setw(8) << setprecision(5) << (double) newTable->getCounts(i, j) / ((double)rp.getNumRuns()*rp.getNumUpdates());
      }
      jointStr << endl;
    }
    cout << "done." << endl;
  }

  // .cluster
  cout << "Writing cluster summary to file " << fileNames.getClusterFile() << "...." << flush;
  ofstream clusterStr(fileNames.getClusterFile().c_str());
  if(clusterStr.fail()) {
    cerr <<"Error: Cannot open file " << fileNames.getClusterFile() << "." << endl;
  }
  clusterStr.setf(ios::fixed, ios::floatfield);
  clusterStr.setf(ios::showpoint);

  vector<double> clusterPP(numGenes+1);
  vector<double> wsum(rp.getNumRuns());
  double wsumAvg=0.0, wsumSD=0.0;
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    wsum[irun]=0.0;

  for (int i=0; i<numGenes+1; i++){
    clusterPP[i]=0.0;
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      clusterPP[i] += (double)clusterCount[irun][i];
    wsumAvg += i * clusterPP[i];
    clusterPP[i] /= ((double)rp.getNumUpdates() * rp.getNumRuns());
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      wsum[irun] += i* (double)clusterCount[irun][i];
  }
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    wsum[irun] /= (double)rp.getNumUpdates();
  wsumAvg /= ((double)rp.getNumUpdates() * rp.getNumRuns());

  clusterStr << "mean #groups = " << setprecision(3) << wsumAvg << endl ;
  if (rp.getNumRuns()>1){
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      wsumSD += (wsum[irun]-wsumAvg)*(wsum[irun]-wsumAvg);
    wsumSD = sqrt(wsumSD / (rp.getNumRuns()-1));
    clusterStr << "SD across runs = " << wsumSD <<endl;
  }
  clusterStr << endl ;

  clusterStr << "credible regions for # of groups" << endl;
  clusterStr << "probability region" << endl;
  clusterStr << "------------------" << endl;

  int a=numGenes,b=0;
  for(int i=0;i<numGenes+1;i++)
    if(clusterPP[i]>0) {
      if(a>i)	a = i;
      if(b<i)	b = i;
    }

  int lo=a,hi;
  double sum = clusterPP[lo];
  while(sum < .005)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .995)    sum += clusterPP[++hi];
  clusterStr << "  0.99      (" << lo << "," << hi << ")" << endl;

  lo=a;  sum = clusterPP[lo];
  while(sum < .025)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .975)    sum += clusterPP[++hi];
  clusterStr << "  0.95      (" << lo << "," << hi << ")" << endl;

  lo=a;  sum = clusterPP[lo];
  while(sum < .050)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .950)    sum += clusterPP[++hi];
  clusterStr << "  0.90      (" << lo << "," << hi << ")" << endl;
  clusterStr << "------------------" << endl << endl;

  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    clusterStr << "Distribution of cluster number in run " << irun+1 << ":" << endl;
    clusterStr << " # of    raw    posterior" << endl;
    clusterStr << "groups  counts probability" << endl;
    clusterStr << "--------------------------" << endl;
    for(int i=a;i<b+1;i++)
      clusterStr << setw(3) << i << " " << setw(10) << clusterCount[irun][i]
		 << setw(12) << setprecision(8) << clusterCount[irun][i] / (double)rp.getNumUpdates() << endl;
    clusterStr << "--------------------------" << endl << endl;
  }
  clusterStr.close();
  cout << "done." << endl;


  // .concordance
  cout << "Writing concordance factors to " << fileNames.getConcordanceFile() << "...." << flush;
  ofstream concordanceStr(fileNames.getConcordanceFile().c_str());
  if(concordanceStr.fail()) {
    cerr <<"Error: Cannot open file " << fileNames.getConcordanceFile() << "." << endl;
  } else{
  concordanceStr.setf(ios::fixed, ios::floatfield);
  concordanceStr.setf(ios::showpoint);

  vector<vector<double> > splitsGeneMatrixPP(splits.size());
  for(int i=0;i<splits.size();i++) {
    splitsGeneMatrixPP[i].resize(numGenes+1);
    for(int j=0;j<numGenes+1;j++){
      splitsGeneMatrixPP[i][j]=0.0;
      for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
	    splitsGeneMatrixPP[i][j] += (double)splitsGeneMatrix[irun][i][j];
      splitsGeneMatrixPP[i][j] /= ((double)rp.getNumUpdates() * rp.getNumRuns());
    }
  }

  // use TreeWeight class for sorting splits by mean number of genes
  vector<TreeWeight> sw(splits.size());
  for(int i=0;i<splits.size();i++) {
    sw[i].setWeight(0);
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrixPP[i][j]>0)
	sw[i].addWeight(j* splitsGeneMatrixPP[i][j]);
    sw[i].setIndex(i);
  }
  sort(sw.begin(),sw.end(),cmpTreeWeights);

  // ctree is the vector of splits in the primary concordance tree
  // gwDistr is the vector of GenomewideDistribution for splits in the primary concordance tree
  vector<TaxonSet> ctree;
  vector<GenomewideDistribution*> gwDistr;
  // otherClade is the vector of splits with CF>.10 but not in the concordance tree
  // otherGwDistr is the vector of GenomewideDistribution for those splits
  vector<TaxonSet> otherClade;
  vector<GenomewideDistribution*> otherGwDistr;

  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    TaxonSet y = splits[i];
    y.setWeight(sw[w].getWeight());
    bool keep=true;
    if (sw[w].getWeight() <= rp.getSwCFcutoff() * numGenes)
      keep=false;
    bool add=true;
    for(int j=0;j<ctree.size();j++) {
      TaxonSet z = ctree[j];
      if(!y.isCompatible(z)) {
	add = false;
	break;
      }
    }
    if(add || keep){
      y.updatePriorProbability();
      if (add)
	    ctree.push_back(y);
      else
	    otherClade.push_back(y);
      GenomewideDistribution* g;
      //TKC: I don't g is ever freed from memory?
      //Need to go through gwDistr and otherGwDistr and delete
      g = new GenomewideDistribution(numGenes, rp.getNumGenomewideGrid());
      g->updatePriorProbability(y);
      g->updateSamplewide(splitsGeneMatrixPP[i]);
      if (!mp.getUseIndependencePrior()){
	// under independence: genomewide CF is concentrated on priorProbability.
	g->updateConvolutionWeight(alphas[0]);
	g->updateGenomewide(alphas[0]);
      }
      double meanpostCFsd=0.0; // SD of estimated sample-wide CF, across runs
      if (rp.getNumRuns()>1){
	    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
	      double meanpostCF =0.0;
	      for(int j=0;j<numGenes+1;j++)
	        meanpostCF += j * (double)splitsGeneMatrix[irun][i][j];
	      meanpostCF /= (double)rp.getNumUpdates();
	      meanpostCFsd += (meanpostCF - y.getWeight()) * (meanpostCF - y.getWeight());
	    }
	    meanpostCFsd = sqrt(meanpostCFsd / (rp.getNumRuns()-1));  // number of genes
	    g->setSamplewidePosteriorMeanSD(meanpostCFsd / numGenes); // proportion of genes
      }
      if (add)
	    gwDistr.push_back(g);
      else
	    otherGwDistr.push_back(g);
    }
  }

  // Now, use ctree to actually build the tree
  // Begin by creating a TaxonSet for each split and then partially ordering them so that subsets precede supersets

  double AvgOfSDacrossSplits = 0.0;

  vector<TaxonSet> tset;
  for(vector<TaxonSet>::iterator p=ctree.begin();p!=ctree.end();p++) {
    TaxonSet a = *p;
    // Currently, splits are represented by the part containing the first taxon.
    // We will instead represent a taxon set for a split as the taxa excluding an outgroup,
    // presumed to be the last taxon.

    if (!a.all() && a[numTaxa - 1]) {
       a = a.flip();
    }

    tset.push_back(a);
  }

  sort(tset.begin(),tset.end());
  tset.push_back(TaxonSet::get());

  concordanceStr << "translate" << endl;
  int taxaNum = 1;
  int taxaSize = translateTable.size();
  for (vector<string>::iterator itr = translateTable.begin(); itr != translateTable.end(); itr++) {
    concordanceStr << " " << taxaNum << " " << *itr;
    if (taxaNum < taxaSize)
      concordanceStr << ",";
    else
      concordanceStr << ";\n";

    concordanceStr << "\n";
    taxaNum++;
  }

  quartet::TreeBuilder tb;
  string quartetTree, quartetTreeWithWts;
  if (rp.shouldBuildPopulationTree()){
    tb.getTree(newTable, numTaxa, quartetTree, quartetTreeWithWts);
    concordanceStr << "Population Tree:" << endl;
    concordanceStr << quartetTree << endl << endl;
  }
  ConcordanceTree z(tset,numGenes);
  concordanceStr << "Primary Concordance Tree Topology:" << endl;
  z.printTopology(concordanceStr);

  if(rp.shouldBuildPopulationTree()){
    concordanceStr<< "Population Tree, With Branch Lengths In Estimated Coalescent Units:" << endl;
    concordanceStr << quartetTreeWithWts << endl << endl;
  }
  concordanceStr << "Primary Concordance Tree with Sample Concordance Factors:" << endl;
  z.print(concordanceStr, numGenes);

  if(rp.shouldBuildPopulationTree()){
    tb.printTies(concordanceStr);
    concordanceStr << endl;
  }
  concordanceStr << "Splits in the Primary Concordance Tree: sample-wide ";
  if (!mp.getUseIndependencePrior())
    concordanceStr << "and genome-wide ";
  concordanceStr << "mean CF (95% credibility)";
  if (rp.getNumRuns()>1)
    concordanceStr << ", SD of mean sample-wide CF across runs";
  concordanceStr << endl;

  for(int w=0;w<ctree.size();w++) {
    ctree[w].print(concordanceStr);
    gwDistr[w]->printSampleCF(concordanceStr);
    if (!mp.getUseIndependencePrior())
      gwDistr[w]->printGenomeCF(concordanceStr);
    if (rp.getNumRuns()>1){
      concordanceStr << "\t" << gwDistr[w]->getSamplewidePosteriorMeanSD();
      AvgOfSDacrossSplits += gwDistr[w]->getSamplewidePosteriorMeanSD();
    }
    concordanceStr << endl;
  }
  concordanceStr << endl;

  concordanceStr << "Splits NOT in the Primary Concordance Tree but with estimated CF > "
		 << rp.getSwCFcutoff() <<":"<<endl;
  for(int w=0;w<otherClade.size();w++) {
    otherClade[w].print(concordanceStr);
    otherGwDistr[w]->printSampleCF(concordanceStr);
    if (!mp.getUseIndependencePrior())
      otherGwDistr[w]->printGenomeCF(concordanceStr);
    if (rp.getNumRuns()>1){
      concordanceStr << "\t" << otherGwDistr[w]->getSamplewidePosteriorMeanSD();
      AvgOfSDacrossSplits += otherGwDistr[w]->getSamplewidePosteriorMeanSD();
    }
    concordanceStr << endl;
  }
  concordanceStr << endl;

  AvgOfSDacrossSplits /= (gwDistr.size()+otherGwDistr.size());
  concordanceStr<<"Average SD of mean sample-wide CF: " << AvgOfSDacrossSplits<<endl<<endl;

  concordanceStr << "All Splits:" << endl;

  //TKC: Free memory from GenomeWideDistribution objects
  //Fixed memory leak
  for (int i = 0; i <gwDistr.size(); i++)
    delete gwDistr[i];
  for (int i = 0; i <otherGwDistr.size(); i++)
    delete otherGwDistr[i];

  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    splits[i].print(concordanceStr);
    concordanceStr << endl;

    int a=numGenes,b=0;
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrixPP[i][j]>0) {
	if(a>j)	a = j;
	if(b<j)	b = j;
      }

    concordanceStr << "#Genes      count in run(s) 1";
    if (rp.getNumRuns()>1)
      concordanceStr << " through "<< rp.getNumRuns();
    concordanceStr << ", Overall probability, Overall cumulative probability" << endl;
    double csum=0;
    for(int j=a;j<b+1;j++) {
      concordanceStr << setw(6) << j ;
      for (int irun=0; irun<rp.getNumRuns(); irun++)
	    concordanceStr << " " << setw(10) << splitsGeneMatrix[irun][i][j];
      csum += splitsGeneMatrixPP[i][j];
      concordanceStr << " " << setw(12) << setprecision(6) << splitsGeneMatrixPP[i][j]
		     << " " << setw(11) << setprecision(6) << csum << endl;
    }
    concordanceStr << endl;
    concordanceStr << "mean CF = " << setw(7) << setprecision(3) << sw[w].getWeight() / numGenes << " (proportion of loci)" << endl;
    concordanceStr << "        = " << setw(7) << setprecision(3) << sw[w].getWeight()            << " (number of loci)" << endl;

    int lo=a,hi;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .005)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .995)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "99% CI for CF = (" << lo << "," << hi << ")" << endl;

    lo=a;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .025)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .975)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "95% CI for CF = (" << lo << "," << hi << ")" << endl;

    lo=a;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .050)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .950)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "90% CI for CF = (" << lo << "," << hi << ")" << endl << endl;
  }

  cout << "done." << endl;

  cout << "Average SD of mean sample-wide CF: " << AvgOfSDacrossSplits<<endl;
  }
  concordanceStr.close();

  // .pairs
  if(rp.getCalculatePairs()) {
    cout << "Writing tree pair data to " << fileNames.getPairTreeFile() << "...." << flush;
    ofstream pairsStr(fileNames.getPairTreeFile().c_str());
    if(pairsStr.fail()) {
      cerr <<"Error: Cannot open file " << fileNames.getPairTreeFile() << "." << endl;
    } else {
      pairsStr.setf(ios::fixed, ios::floatfield);
      pairsStr.setf(ios::showpoint);
      double total = pairCounts[0][0];
      for(int i=0;i<numGenes;i++) {
	pairsStr << setw(4) << i << " ";
	for(int j=0;j<numGenes;j++)
	  pairsStr << setw(7) << setprecision(4) << (double) (i < j ? pairCounts[i][j] : pairCounts[j][i]) / total;
	pairsStr << endl;
      }
      pairsStr.close();
      cout << "done." << endl;
    }
  }

// .topologies --- eliminated
//  cout << "Writing topology single and joint posterior distribution to " << fileNames.getTreePosteriorFile() << "...." << flush;
//  ofstream treePosteriorStr(fileNames.getTreePosteriorFile().c_str());
//  treePosteriorStr.setf(ios::fixed, ios::floatfield);
//  treePosteriorStr.setf(ios::showpoint);
//  treePosteriorStr << "index "
//		   << setw(max) << left << "topology" << right
//		   << setw(10) << "single"
//		   << setw(10) << "joint" << endl;
//  for(int i=0;i<numTrees;i++) {
//    double single=0.0;
//    int joint=0;
//    for(int j=0;j<numGenes;j++) {
//      single += genes[j]->getProb(i);
//      joint += newTable[i][j];
//    }
//    treePosteriorStr << setw(5) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right
//		     << setw(10) << setprecision(6) << (double) single / (double) numGenes
//		     << setw(10) << setprecision(6) << (double) joint / (double) rp.getNumUpdates() / (double) numGenes
//		     << endl;
//  }
//  treePosteriorStr.close();
//  cout << "done." << endl;

// .gene

  cout << "Writing single and joint gene posteriors to " << fileNames.getGenePosteriorFile() << "...." << flush;
  ofstream genePostStr(fileNames.getGenePosteriorFile().c_str());
  if(genePostStr.fail()) {
    cerr <<"Error: Cannot open file " << fileNames.getGenePosteriorFile() << "." << endl;
  } else {
    for(int i=0;i<numGenes;i++)
      genes[i]->print(genePostStr,newTable,rp.getNumUpdates()*rp.getNumRuns(),topologies,max);
    genePostStr.close();
    cout << "done." << endl;
  }

  if(rp.getNumChains()>1) {
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    for (int irun=0; irun<rp.getNumRuns(); irun++){
      cout << "\nMCMCMC acceptance statistics in run " << irun+1 <<":" << endl;
      cout << "         alpha1 <-->         alpha2 accepted proposed proportion" << endl;
      for(int i=0;i<rp.getNumChains()-1;i++) {
	cout << setw(15) << alphas[i] << " <-->" << setw(15) << alphas[i+1];
	cout << setw(9) << mcmcmcAccepts[irun][i] << setw(9) << mcmcmcProposals[irun][i];
	if(mcmcmcProposals[irun][i] > 0) {
	  cout << setw(11) << setprecision(6) << (double) (mcmcmcAccepts[irun][i]) / (double) (mcmcmcProposals[irun][i]);
	}
	cout << endl;
      }
    }
  }

}


void MPI_seed(int s1, int s2, int procs, std::vector<int> &samples){
    srand(s1); 
    for (int i=0; i <procs*2; i+=2){
        samples[i]=rand(); 
        //cout << "Random num " << i << " is " << samples[i] << endl;
    }
    srand(s2); 
    for (int i=1; i <= procs*2; i+=2){
        samples[i]=rand(); 
        //cout << "Random num " << i << " is " << samples[i] << endl;
    }
}








/* main */







//TKC:
//Note: Valgrind shows a ton of memory leaks. 

/*Known issues
 * --mpi process aborts when proc has more than 1 run
 * ----Double free, may have fixed with most recent edits (24 Apr)
 * ----Check further, run tests
 * --Error on calling --help
 * --Dividing runs/chains not correct. 5 runs on 2 processes, each proc gets 3
 * --create-sample-file needs to be called per run (e.g. for.. i<my_runs)
 * ----Data structure needs to match
 * --*******Seg fault when a run is split across two threads********
 * ----Occurs AFTER burn-in
 * ----Not in MPI_MCMCMC
 * 
 * --A lot of this would probably be easier if we made it more object-oriented
 * ----e.g. "Chain" object which has its own members for id, alpha, state
 * ----Rather than having all of these independent collections (more error prone)
 * ----Consider re-write, once bugs with current version are worked out.
 * */



int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  int my_rank, p;
  my_rank = MPI::COMM_WORLD.Get_rank(); //Get rank
  MPI_Comm_size(MPI::COMM_WORLD, &p); //Get number of processes
  
  //cout << flush; 
  cout << "\nHello from rank " << my_rank <<"\n";
  MPI_Barrier(MPI_COMM_WORLD); 
   
  time_t beginTime;
  time(&beginTime);
  if (my_rank == 0)
    intro(cout);
  // Set Default Parameters
  int debug=0;
  FileNames fileNames("run1");
  ModelParameters mp;
  RunParameters rp;
  Defaults defaults(fileNames,mp,rp);
  
  
  // Read arguments from the command line
  // return value k is the position of the first input file
  //   i.e. argv[k] is the first input file
  int k=readArguments(my_rank,argc,argv,fileNames,mp,rp,defaults);

  // Populate input files vector from remaining command line arguments
  vector<string> inputFiles(argv + k, argv + argc);

  // If input list filename is given as argument, read input filenames
  if (!(fileNames.getInputListFileName().empty())) {
    readInputFileList(fileNames.getInputListFileName(), inputFiles);
  }

  //TKC: Calc number of runs each process rank is responsible for

  int total = rp.getNumChains()*rp.getNumRuns(); 
  int start, end; 

  //Get group of processes in world communicator
  MPI_Group world_group; 
  MPI_Comm_group(MPI_COMM_WORLD, &world_group); 
  MPI_Comm MPI_New_World; 
  
  if (total < p){
      if (my_rank <= total){
        start = my_rank; 
        end = my_rank; 
      }
      MPI_Group new_group; 
      vector<int> subset(total);
      for (int i=0; i<total; i++){
          subset[i]=i; 
      }
      //Populate new world communicator with the subset that is needed
      MPI_Group_incl(world_group, total, subset.data() , &new_group);  
      MPI_Comm_create(MPI_COMM_WORLD, new_group, &MPI_New_World); 
      
  }else{
    start = my_rank*(total/p); 
    end = ((total/p)+my_rank*(total/p))-1; 
    if (my_rank == p-1)
      end += total%p; 
    //Populate new world communicator
    MPI_Comm_create(MPI_COMM_WORLD, world_group, &MPI_New_World); 
  }

  //Terminate unused processes 
  if (MPI_New_World == MPI_COMM_NULL){
      //Bye cruel world 
      cout << "Warning: Process rank "<<my_rank<<" is terminating early. "
              "Please check number of processes requested."<<endl; 
      MPI_Finalize(); 
      exit(0); 
  }
  //Recapture number of processes and rank
  //cout <<"Old world had: "<<p<<" processes\n";
  MPI_Comm_size(MPI_New_World, &p);
  my_rank = MPI::COMM_WORLD.Get_rank();
  //cout <<"New world has: "<<p<<" processes\n";

  MPI_Barrier(MPI_New_World); 
  //cout << flush;   
  //cout << "Rank " << my_rank <<" has chains: ";
  //for (int i=start; i <= end; i++){
   //   cout<< i << " ";
  //}
  //cout << endl;  
  
  
  int numGenes;
  int max=0;
  int numTaxa;
  vector<string> translateTable;
  vector<vector<int> > taxid;
  ofstream fout(fileNames.getOutFile().c_str());
  // Open outfile to save all window output
  if (my_rank==0){
    cout << "Screen output written to file " << fileNames.getOutFile() << endl;
    ofstream fout(fileNames.getOutFile().c_str());
    if(fout.fail()) {
        cerr <<"Error: Cannot open file " << fileNames.getOutFile() << "." << endl;
        MPI_Finalize();
        exit(1);
    }
    intro(fout);

    cout << "Program initiated at " << ctime(&beginTime) << endl << endl;

    for(int i=0;i<argc;i++)
        fout << " " << argv[i];
    fout << endl << endl;
    showParameters(fout,fileNames,defaults,mp,rp);
    fout << endl;

    cout << "Reading in summary files...";
  }

  TGM *topToGeneMap = new TGM();
  MPI_readInputFiles(my_rank, inputFiles,topToGeneMap,translateTable,max,fout,taxid, rp.getPruneFile(), rp.getPruneGene());
  numGenes = inputFiles.size();
  numTaxa = translateTable.size();
  if (my_rank==0){
    cout << "done." << endl;
  }
  
  // check that genome-wide grid size > number of genes: otherwise the genome-wide estimates are innacurate.
  if (rp.getNumGenomewideGrid() < 3*numGenes){
    if (my_rank ==0){
      cout << "The grid size to get genome-wide CFs is too small compared to the number of sampled genes...";
      cout << " changing this grid size to "<< 3*numGenes << " (3 * # genes)"<<endl;
    }
    rp.setNumGenomewideGrid(3*numGenes);
  }

  // Check for missing taxa.
  for(int i=0; i<inputFiles.size(); i++) {
    bool problem = false;
    if (taxid[i].size() != translateTable.size())
      problem=true;
    if ((problem) && (my_rank==0)) {
      cerr <<"\nError: Expecting " << translateTable.size() << " taxa."<<endl;
      cerr <<  "       File " << inputFiles[i] << " contains trees with the following " << taxid[i].size() << " taxa:" <<endl;
      for (int j=0; j<taxid[i].size(); j++){
		cerr << setw(4) << j+1 << " " << setw(4) << taxid[i][j];
		if (taxid[i][j] <= translateTable.size()){
		  cerr << " " << translateTable[taxid[i][j]-1] << endl;
	    }else{
		  cerr << " " << "(outside range of translate table)" <<endl;
        }
      exit(1);
	  }
	}
  }

  // Find taxa sequenced in all genes
  TaxonList allTaxList(translateTable);
  allTaxList.setNumberGenes(taxid);
  vector<int> taxlist = allTaxList.getNonMissing();
  allTaxList.setInclude(taxlist);
  if (my_rank == 0){
    allTaxList.print(cout);
  }
  // fixit:
  // prune the missing taxa from all topologies, which is a vector<string>
  // find the unique topologies
  // collapse the table to only unique topologies, with updated counts
  // ??update the translate table and the taxid vectors??

  int numTrees = topToGeneMap->getNumTrees();
  
  //Rank 0 only
  if (my_rank ==0){
    cout << "Read " << numGenes << " genes with a total of " << numTrees << " different sampled tree topologies" << endl;
    cout << "Writing input file names to file " << fileNames.getInputFile() << "...." << flush;
    
    ofstream inputStr(fileNames.getInputFile().c_str());
    inputStr << "Gene Filename" << endl;
    inputStr << "============================================" << endl;
    for(int i=0;i<inputFiles.size();i++)
        inputStr << setw(4) << i << " " << inputFiles[i] << endl;
    inputStr << "============================================" << endl;
    cout << "done." << endl;

    cout << "Sorting trees by average posterior probability..." << flush;
  }
  
  // -------------------------------!!!------------------------------
  //Need to check that each proc sorts topToGeneMap the same. Could alternativey  
  //create in Rank 0 and bcast it                                    
  //--------------------------------!!!------------------------------
  
  topToGeneMap->reorder(); // sort according to topology counts
  vector<string> topologies = topToGeneMap->getTopologies();
  if (my_rank ==0){
    cout << "done." << endl;
  }
  // Print table of individual gene posteriors

  // .single
  // Rank 0 only
  if((rp.getCreateSingleFile()) && (my_rank==0)) {
    cout << "Writing single gene posterior distribution table to file " << fileNames.getSingleFile() << "...." << flush;
    ofstream singleStr(fileNames.getSingleFile().c_str());
    singleStr.setf(ios::fixed, ios::floatfield);
    singleStr.setf(ios::showpoint);

    for(int i=0;i<numTrees;i++) {
      singleStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
      for(int j=0;j<numGenes;j++)
        singleStr << " " << setw(8) << topToGeneMap->getCounts(topologies[i], j);
        singleStr << endl;
    }
    singleStr.close();
    cout << "done." << endl;
  }

  // Initialize random number generator.

  
  /*TKC: Master takes user seeds and samples np*2 random ints 
   Then MPI_Scatter random ints to processes, which use these as local seeds*/
  
  if (my_rank ==0){
    cout << "Initializing random number generator...";
  }
  
  vector<int> global_seeds; 
  unsigned int swap_seed; 
  if (my_rank == 0){
    global_seeds.resize(p*2);
    MPI_seed(rp.getSeed1(), rp.getSeed2(), p, global_seeds);
    //cout << "Starting with seeds "<<rp.getSeed1()<<" and "<<rp.getSeed2()<<endl;
    swap_seed = rp.getSeed1()+rp.getSeed2();

  }
  int local_seeds[2]; 
  //Send seed for swap accepts to all processes
  MPI_Bcast(&swap_seed, 1, MPI_INT, 0, MPI_New_World); 
  srand(swap_seed); 
  vector <Rand*> swap_rand; 
  swap_rand.resize(rp.getNumRuns());
  for (int irun=0; irun<rp.getNumRuns(); irun++){
      swap_rand[irun] = new Rand(rand(), rand()); 
      unsigned int s1, s2; 
      swap_rand[irun]->getSeed(s1, s2); 
      //cout << my_rank <<" has seed run "<<irun<<" of "<<s1<<" - "<<s2<<endl; 
  }

  //cout << my_rank<<" has swap seed of "<<swap_seed<<endl; 

  //for (int i=0; i<10; i++){
  //    cout << my_rank<<" (swap) sample "<<i<<" is "<<int(swap_rand.runif()*rp.getNumChains())<<endl;
  //}
  
  Rand mcmc_rand(rp.getSeed1(),rp.getSeed2());
  if (my_rank ==0){
    cout << "done." << endl;
  }
  MPI_Scatter(global_seeds.data(), 2, MPI_INT, local_seeds, 2, MPI_INT, 0, MPI_New_World);
  rp.setSeed1(local_seeds[0]);
  rp.setSeed2(local_seeds[1]);
  //MPI_Barrier(MPI_New_World);
  //cout << flush; 
  //cout << "Rank "<<my_rank << " has seed " << rp.getSeed1() << " - " << rp.getSeed2() << endl;


  // Create genes
  if (my_rank ==0){
    cout << "Initializing gene information..." << endl << flush;
  }
  
  vector<Gene*> genes(numGenes);
  vector<double> counts(numTrees);

  // create gene objects
  for(int i=0;i<numGenes;i++) {
    genes[i] = new Gene(i);
  }

  // update counts for gene objects
  int topIndex = 0;
  for (vector<string>::iterator titr = topologies.begin(); titr != topologies.end(); titr++, topIndex++) {
    const vector<GeneCounts>& geneCounts = topToGeneMap->getCounts(*titr);
    for (vector<GeneCounts>::const_iterator gitr = geneCounts.begin(); gitr != geneCounts.end(); gitr++) {
      genes[gitr->getGene()]->addCount(topIndex, gitr->getCount());
    }
  }

  // update probabilities and create cells using Alias method for faster sampling
  for(int i=0;i<numGenes;i++) {
    genes[i]->updateState(topologies.size());
  }
  //cout << "done." << endl << flush;
  // Finished with table and taxid, so delete them.
  delete topToGeneMap;
  for(int i=0;i<taxid.size();i++){
    taxid[i].clear();
  }
  taxid.clear();

  vector<vector<TaxonSet> > topologySplitsMatrix(topologies.size());
  if(numTaxa>3)
    for(int i=0;i<topologies.size();i++)
      topologySplitsMatrix[i].resize(numTaxa-3);
  vector<TaxonSet> splits; // sorted list of realized splits
  
  for(int i=0;i<topologies.size();i++) {
    TaxonSet all = TaxonSet::getTaxa(numTaxa);
    Tree t(all,topologies[i]);
    SplitSet s;
    t.getSplits(s);

    for(int j=0;j<numTaxa-3;j++) {
      TaxonSet x = s.getSplit(j);
      topologySplitsMatrix[i][j] = x;
      int n = find(splits.begin(),splits.end(),x) - splits.begin();
      if(n==splits.size())
	    splits.push_back(x);
    }
  }
//  topologyStr.close();
//  cout << "done." << endl;
  if (my_rank == 0){
    cout << "Initializing splits counts (found " << splits.size() << " distinct splits)..." << flush;
  }
  sort(splits.begin(),splits.end());

  // topologySplitsIndexMatrix[i][j] = index into splits vector of the jth split of topology i
  vector<vector<int> > topologySplitsIndexMatrix(topologies.size());
  if(numTaxa>3) {
    for(int i=0;i<topologies.size();i++) {
      topologySplitsIndexMatrix[i].resize(numTaxa-3);
      for(int j=0;j<numTaxa-3;j++)
	topologySplitsIndexMatrix[i][j] = find(splits.begin(),splits.end(),topologySplitsMatrix[i][j]) - splits.begin();
    }
  }

  // splitsGeneMatrix[irun][i][j] = the number of times that split i is in exactly j genes, in run 'irun'.
  //TKC: Change to 1 serialized vector per run, easier for MPI later
  vector<vector<int>> splitsGeneMatrix(rp.getNumRuns());
  for (int irun=0; irun<rp.getNumRuns(); irun++){
    splitsGeneMatrix[irun].resize((splits.size()*numGenes+1));
    for(int i=0;i<splits.size();i++) {
      for(int j=0;j<numGenes+1;j++)
	    splitsGeneMatrix[irun][(i*numGenes+1)+j] = 0;
    }
  }
  if (my_rank == 0){
    cout << "done." << endl;
  }
  
  if(rp.getNumUpdates()==0){
	  MPI_Finalize();
    return 0;
  }
  if (my_rank ==0){
    cout << "Setting initial MCMC state..." << flush;
  }

  //Allocate space for alphas
  vector<double> global_alphas(rp.getNumChains());
  global_alphas.resize(rp.getNumChains()*rp.getNumRuns());
    
  //Allocate space for collecting states later (have each process make a local version)
  if (my_rank==0)
    vector<vector<State*> > global_states(rp.getNumRuns());

  //Make serial list of indices
  vector<int> global_index(total);
  vector<int> global_runs(total);
  vector<int> global_ranks(total);
    
  //Resize index vectors
  global_index.resize(total);
  global_runs.resize(total);
  int position=0;
  int numPer; 
  if (total < p){
    numPer=1; 
  }else{
    numPer=total/p;
  }
  int rankIdx=0; 
  int count=0; 
  for (int irun=0; irun<rp.getNumRuns(); irun++){
    //Allocate space for later collecting final states
    //global_states[irun].resize(rp.getNumChains()); 
    for(int i=0;i<rp.getNumChains();i++){  
      if (i==0){
        global_alphas[irun*rp.getNumChains()] = mp.getAlpha();
      }else{
        global_alphas[i+(irun*rp.getNumChains())]=rp.getAlphaMultiplier()*global_alphas[i-1];
      }
      //cout << "Global alpha "<<i<<" is "<<global_alphas[i]<<endl;
      //cout << "Run " << irun+1 <<" chain "<<i+1<<" has alpha "<<global_alphas[i]<<endl;
            
      //Place default objects in global_states -> need to make a copy constructor for later collecting (?)
      //global_states[irun][i] = new State(global_alphas[i],numTaxa,numTrees,genes,mp.getUseIndependencePrior(),rand);
      //Store chain indices
            
      if (rankIdx == p-1){
        numPer += total%p; 
      }else if (count == numPer){
        count = 0; 
        rankIdx++; 
      }
      global_ranks[position] = rankIdx; 
      count++; 
                 
      global_index[position] = i;
      global_runs[position] = irun; 
      position++; 
    }
  }
    
  int my_runs = (global_runs[end]+1)-(global_runs[start]); 
  int my_chains = end+1 - start;
    
  vector<State*>local_states(total); 
    
  for (int i=start; i<=end; i++){
    //Each local_states will only fill chains under each rank
    //Place default objects in global_states -> need to make a copy constructor for later collecting (?)
    local_states[i] = new State(global_alphas[i],numTaxa,numTrees,genes,mp.getUseIndependencePrior(),mcmc_rand);
  }
  if (my_rank ==0){
    cout << "done." << endl <<flush;
  }

  cout << flush;      
 
  if (my_rank == 0){
    cout << "Global indices: "; 
    for (int i=0; i < global_index.size(); i++)
      cout << global_index[i] << " ";
    cout << endl << "Global run indices: ";
    for (int k=0; k < global_runs.size(); k++)
      cout << global_runs[k] << " ";
    cout << endl;
    cout << endl << "Global rank indices: ";
    for (int j=0; j < global_ranks.size(); j++)
      cout << global_ranks[j] << " ";
    cout << endl;
      
  }
  
    //cout << "Rank "<<my_rank<<" starts at "<<start<<" and ends at "<<end << endl <<flush;
    //cout << my_rank << " has " << my_chains << " chains\n";
    //cout << my_rank << " has " << my_runs << " runs\n";
    //cout << my_rank << " has chains: ";
    //for (int i=start; i<=end; i++)
    //  cout << global_index[i] << " ";
    //cout << endl;
  
  
  if (my_rank==0){
    cout << "Initializing MCMCMC acceptance counters and pairwise counters..." << flush;
  }  

    
  //TKC: Again, for MPI version only need 1D vectors
  vector<int> local_mcmcmcAccepts(total);
  vector<int> local_mcmcmcProposals(total);
  
  for (int i=0; i<total; i++){
    local_mcmcmcAccepts[i] = local_mcmcmcProposals[i] = 0;
  }

  vector<int> local_pairCounts(numGenes*numGenes); // not adapted to several runs. Runs will be pooled.
  for(int i=0;i<numGenes;i++) {   
    for(int j=0;j<numGenes;j++)
      local_pairCounts[(i*numGenes)+j] = 0;
  }
  
  if (my_rank==0){
    cout << "done." << endl;
  }
  
  MPI_Barrier(MPI_New_World); 
  time_t beginMCMCtime;
  if (my_rank==0){
    cout << "MCMC initiated at " << ctime(&beginMCMCtime) << endl;
  }
  int numBurn = rp.getNumUpdates()/10;
  if (my_rank==0){
    cout << "Beginning burn-in with " << numBurn << " updates (10% extra of desired updates)...";
  }

  

  int part = numBurn / 50;
  int thisRun; 
  
  //BEGIN BURN IN CYCLE

  for(int cycle=0;cycle<numBurn;cycle++) {
    //TKC: Need to take a look at updateOneGroup()   
    //Update each chain
    for(int i=start;i<=end;i++) {
      local_states[i]->update(mcmc_rand);
	  if(rp.getUseUpdateGroups()) {
		int gene = (int)(mcmc_rand.runif()*genes.size());
	      local_states[i]->updateOneGroup(gene,mcmc_rand);
	  }
	}	
    //TKC: If 2 or more chains, and mcmcmc interval is met, 
    // then for each run, swap chain states
    if(cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1){ 
	  int begin = global_runs[start];
      for(int irun=begin; irun < (begin + my_runs); irun++){      

        //MCMCMC
        //thisRun = global_runs[i];
        MPI_mcmcmc(local_states, global_alphas, swap_rand[irun], 
          local_mcmcmcAccepts, local_mcmcmcProposals, my_rank, 
          global_runs, global_index, rp.getNumChains(), 
          irun, global_ranks, MPI_New_World);   
      }
	} 
  }

  MPI_Barrier(MPI_New_World); 

  if (my_rank==0){
    //cout << "*" << endl;
    cout << "done." << endl << flush;
  }
  
  if (my_rank==0){
    cout << "Initializing summary tables..." << flush;
  }
    
  vector<int> local_clusterCount(rp.getNumRuns()*(numGenes+1));
  for (int irun=0; irun<rp.getNumRuns(); irun++){
    for(int i=0;i<numGenes+1;i++)
      local_clusterCount[(irun*numGenes+1)+i] = 0;
  }
  
  MPI_Barrier(MPI_New_World); 
  if (my_rank==0){
    cout << "done." << endl << flush;
  }

  
  if (my_rank==0){
    cout << "Beginning " << rp.getNumUpdates() << " MCMC updates...";
  }
  //int part = rp.getNumUpdates() / 50;
  

  vector<ofstream*> sampleFileStr(my_runs); //Declare vector of ofstream
  //Some local processes will not use
  
  for (int i=0; i<sampleFileStr.size(); i++){
	  //cout << "Rank " << my_rank << " has: " << i;
  }
  //cout << flush; 

//Not currently working!!!
  /*if(rp.getCreateSampleFile()) {
      //for (int i=start; i <= end; i++){
        //if (local_ind
		
      for (unsigned int irun=0; irun<my_runs; irun++){
		thisRun = global_runs[((irun*rp.getNumChains())+start)];
		cout << my_rank << " has file for run " << thisRun << endl;

        sampleFileStr[irun] = new ofstream(fileNames.getSampleFile(thisRun+1).c_str());
        cout << fileNames.getSampleFile(thisRun+1) << endl;
        if (sampleFileStr[irun]->fail()){
	      cerr << "Error: could not open file " << fileNames.getSampleFile(thisRun+1) << endl;
	      exit(1);
        }
        sampleFileStr[thisRun]->setf(ios::fixed, ios::floatfield);
        sampleFileStr[thisRun]->setf(ios::showpoint);
      }
    } 
  }*/
  
  
  //TKC: Serialized 2D vector
  vector<int> local_accept(total);
  for (int i=0; i<total; i++){
    local_mcmcmcAccepts[i] = local_mcmcmcProposals[i] = local_accept[i] = 0;
  }
  
  TGMTable* localTable;
  localTable = new TGMTable(numGenes, topologies);

  //FULL MCMC
  for(int cycle=0;cycle<rp.getNumUpdates();cycle++) {  
    //Update each chain
    for(int i=start;i<=end;i++) {
	  //Update accepts
      local_accept[i] += local_states[i]->update(mcmc_rand);
	  if(rp.getUseUpdateGroups()) {
	    int gene = (int)(mcmc_rand.runif()*genes.size());
	    local_states[i]->updateOneGroup(gene,mcmc_rand);
	  }
	}	
	
    //TKC: If 2 or more chains, and mcmcmc interval is met, 
    // then for each run, swap chain states
    if (cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1){
	  int begin = global_runs[start];
      for(int irun=begin; irun < (begin + my_runs); irun++){     
        //MCMCMC
        MPI_mcmcmc(local_states, global_alphas, swap_rand[irun], 
          local_mcmcmcAccepts, local_mcmcmcProposals, my_rank, 
          global_runs, global_index, rp.getNumChains(), 
          irun, global_ranks, MPI_New_World);
	  }
    }
    
    //Execute the following if process controls a rank 0 chain
    //NEED TO LOOK INTO CONTENTS OF CLUSTERCOUNT AND SPLITSGENEMATRIX
    //ALSO WHY ARE SPLITS, COUNTS, AND TABLE UPDATED EVERY SINGLE ITERATION? 
    //SHOULD THEY BE SUBSAMPLED??
	for (int i=start; i<= end; i++){
	//If chain is cold
	  if (local_states[i]->getAlpha() == mp.getAlpha()){
	    local_states[i]->updateTable(localTable);
		local_states[i]->updateSplits(splitsGeneMatrix[global_runs[i]],topologySplitsIndexMatrix, splits.size());
		local_clusterCount[(global_runs[i]*numGenes)+(local_states[i]->getNumGroups())]++;
		if (rp.getCalculatePairs() && cycle % rp.getSubsampleRate() == 0){
		  local_states[i]->updatePairCounts(local_pairCounts); 
		/*if (rp.getCreateSampleFile() && cycle % rp.getSubsampleRate() == 0){
		  *(sampleFileStr[global_runs[i]]) << setw(8) << local_accept[i];
		  local_accept[i] = 0;
		  local_states[i]->sample(*sampleFileStr[global_runs[i]]);*/
		}
	  }
	} 
  }  

  MPI_Barrier(MPI_New_World); 
  
  if (my_rank == 0)
    cout << "done." << endl << flush;  


  //Free up memory used by the ofstreams
  /*
  if(rp.getCreateSampleFile() && hasCold==true){
    for (int i=0; i<sampleFileStr.size(); i++){
      if (sampleFileStr[i] != NULL){
		  //cout << "-----Rank " << my_rank << " deleting index: " << i << endl;
	    sampleFileStr[i]->close();
	    //delete sampleFileStr[i];
	  }
    }
  }*/
  
  
  //Check final states
  string name = to_string(my_rank) + ".state";
  ofstream f(name);
  for (int i=start; i<end; i++)
    local_states[i]->print(f);
  
      
  //Check final table
  string name2 = to_string(my_rank) + ".table";
  ofstream g(name2);
  localTable->print(g);
  

  //Doughter processes send, master collects all final data structures

  vector<int> serialTable;
  serialTable = localTable->getSerialTable();
  int sz = serialTable.size();
  //If master, collect 
  if (my_rank == 0){
	//Foreach daughter process, receive
	for (int i=1; i<p; i++){
	  vector<int> tempTable;	  
	  tempTable.resize(sz);
	  int tgmTag = 100 + i;
	  
	  //Receive serialized TGM table from daughter
	  MPI_Recv(&tempTable[0], sz, MPI_INT, i, tgmTag, MPI_New_World, MPI_STATUS_IGNORE);
	  localTable->readSerialTable(tempTable);
	  tempTable.clear();
	  
	  //Receive serialized clusterCount from daughter
	  vector<int> tempClust;
	  tempClust.resize(local_clusterCount.size());
	  int clustTag = 200 + i;
	  MPI_Recv(&tempClust[0], tempClust.size(), MPI_INT, i, clustTag, MPI_New_World, MPI_STATUS_IGNORE);
      for (int k=0; k < tempClust.size(); k++){
        local_clusterCount[k] += tempClust[k];
        tempClust.clear(); 
	  }
      
      //Receive serialized splitsGeneMatrix vectors for each run
      for (int irun=0; irun < rp.getNumRuns(); irun++){
		vector<int> tempSplits((splits.size()*numGenes+1));
		int tag = ((30 + i)*10) + irun;
        MPI_Recv(&tempSplits[0], tempSplits.size(), MPI_INT, i, tag, MPI_New_World, MPI_STATUS_IGNORE);
	    for (int j=0; j<tempSplits.size(); j++){
		  splitsGeneMatrix[irun][j] += tempSplits[j];  
	    }
	    tempSplits.clear();
	  }
	  
	  //Receive pairCounts
	  int pairTag = 400 + i;
	  vector<int> tempCounts(local_pairCounts.size());
	  MPI_Recv(&tempCounts[0], tempCounts.size(), MPI_INT, i, pairTag, MPI_New_World, MPI_STATUS_IGNORE);
	  for (int k=0; k<tempCounts.size(); k++)
	    local_pairCounts[k] += tempCounts[k];
	  tempCounts.clear();
        
      //Receive alphas
      int alphaTag = 500 + i;
      vector<double> tempAlphas(total);
      MPI_Recv(&tempAlphas[0], total, MPI_DOUBLE, i, alphaTag, MPI_New_World, MPI_STATUS_IGNORE);
      for (int l=0; l<total; l++){
	    if (tempAlphas[l] != 0){
	      global_alphas[l] = tempAlphas[l];
	    }
	  }
	  //Update final alpha values for master
	  for (int m=start; m<=end; m++)
	    global_alphas[m] = local_states[m]->getAlpha();    
	}	
	//Else, send to master	
  }else{
	  
	//Send TGMTable
	int tgmTag = 100 + my_rank;
	MPI_Send(&serialTable[0], sz, MPI_INT, 0, tgmTag, MPI_New_World);
	serialTable.clear();
	
	//Send local_clusterCount
	int clustTag = 200 + my_rank;
	MPI_Send(&local_clusterCount[0], local_clusterCount.size(), MPI_INT, 0, clustTag, MPI_New_World);
    local_clusterCount.clear();
    
    //For each run, send serialized splitsGeneMatrix
    int splitTag = (30 + my_rank)*10;
    for (int irun=0; irun < rp.getNumRuns(); irun++){
	  int tempTag = splitTag + irun;
	  vector<int> tempSplits(splits.size());  
	  MPI_Send(&splitsGeneMatrix[irun][0], splitsGeneMatrix[irun].size(), MPI_INT, 0, tempTag, MPI_New_World);
	  splitsGeneMatrix[irun].clear();
	}
	//Send local_pairCounts
	int pairTag = 400 + my_rank;
	MPI_Send(&local_pairCounts[0], local_pairCounts.size(), MPI_INT, 0, pairTag, MPI_New_World);
	local_pairCounts.clear();
	
	//Send alphas
	vector<double> local_alphas(total);
	int alphaTag = 500 + my_rank;
	for(int i=start;i<=end;i++) 
	  local_alphas[i] = local_states[i]->getAlpha();
	 MPI_Send(&local_alphas[0], total, MPI_DOUBLE, 0, alphaTag, MPI_New_World);
	 local_alphas.clear();
	 global_alphas.clear();
  }
  
  //Exit, below not tested yet
  //MPI_Finalize();
  //exit(0);
  
  MPI_Barrier(MPI_New_World); 
   
   if (my_rank == 0){
     string name2 = "final.table";
     ofstream g(name2);
     localTable->print(g);
     cout << "Final alphas: ";
     for (int i=0; i<total; i++)
       cout<< global_alphas[i] << " ";
     cout << endl;
   }
   
   //NOTE:  
   //topologies -- same
   //local_table -- serialized, send done, incorporated by master
   //topSplitsIndexMatrix -- same
   //clusterCOunt -- Serialized, send done, incorporated by master
   //splitsGeneMatrix -- Serialized, send done, incorporated
   //topologySplitsIndexMatrix-- Same
   //local_pairCounts -- Serialized, sent, incorporated
   //alphas -- All collected and incorporated
   //accepts
   //proposals
   //translateTable
   //splits
   
   //DO ALL MPI_GATHERV'S TO RANK 0 IN HERE, ONLY RANK 0 CREATES FINAL OUTS
    
   /* Need to gather:
    * alphas
    * accepts
    * proposals
    * */
  
    
  //writeOutput(fout,fileNames,max,numTrees,numTaxa,topologies,numGenes,rp,mp,
	//      newTable,clusterCount,splits,splitsGeneMatrix,
	//      pairCounts,genes,alphas,mcmcmcAccepts,mcmcmcProposals, translateTable);

        
  for(int i=0;i<numGenes;i++)
    delete genes[i];

  for(int i=start;i<=end;i++)
    delete local_states[i];
  
  delete localTable;
      
//^^Make sure to do the same for global_states (rank 0)

  time_t endTime;
  time(&endTime);

  if (my_rank==0){
    cout << "Program ended at " << ctime(&endTime) << endl << endl;
    int diff=endTime-beginTime,days=diff/(24*60*60),hours=diff%(24*60*60)/(60*60);
    int minutes=diff%(60*60)/60,seconds=diff%60;
    fout.close();
  }
  
 //} //END RANK 0 WORK
  MPI_Finalize();
  return 0;
} 
