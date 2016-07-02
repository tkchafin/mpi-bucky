#ifndef BUCKYHDR
#define BUCKYHDR

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include "boost/dynamic_bitset.hpp"
#include "boost/unordered_map.hpp"
#include "TaxonSet.h"
#include "mbsumtree.h"
#include "Rand.h"
#include "Alias.h"
#include "TGM.h"
#include "Quartets.h"
using namespace std;

static double LOG_ZERO = log(double(0));

class Gene {
public:
  Gene(int n) {
    number = n;
    total = 0.0;
    numTrees = 0;
  }

  void addCount(int topIndex, double count) {
    if (count > 0) {
      numTrees++;
      total += count;
      indices.push_back(topIndex);
      counts.push_back(count);
    }
  }

  void updateState(int size) {
    vector<double> probs(indices.size(), 0.0);
    logProbs.resize(size,LOG_ZERO);
    hasNonZeroProb.resize(size, false);
    for(int i=0;i<indices.size();i++) {
      double prob = counts[i]/total;
      probs[i] = prob;
      logProbs[indices[i]] = log(prob);
      hasNonZeroProb[indices[i]] = true;
    }
    //added for faster probability retreival
    alias = new Alias(probs, indices);
  }
  ~Gene() {
    indices.clear();
    counts.clear();
    logProbs.clear();
    delete alias;
  }
  int getNumber() const { return number; }
  int getNumTrees() const { return numTrees; }
  double getTotal() const { return total; }
  int getIndex(int i) const { return indices[i]; }
  double getCount(int i) const { return counts[i]; }
  int pickTreeFast(Rand& r) const {
    return alias->pick(r);
  }

  // Can return log(0) if the gene prior does not have top.
  // Use this in combination with hasProb to not get log(0)
  double getLogProb(int top) const { return logProbs[top]; }

  bool hasProb(int top) const { return hasNonZeroProb[top]; }

//  double getProb(int top) const {
//    int i=0;
//    while(indices[i]!=top && i<indices.size())
//      i++;
//    if(i==indices.size())
//      return 0.0;
//    else
//      return (double)(counts[i])/(double)(total);
//  }
  void print(ostream& f,string file) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);

    f << "Gene " << number << " (from file " << file << "):" << endl;
    f << "  numTrees = " << numTrees << endl;
    f << "    i indices probability" << endl;
    for(int i=0;i<numTrees;i++)
      f << setw(5) << i << ":" << setw(7) << indices[i] << setw(12) << setprecision(8) << counts[i] << endl;
    f << endl;
  }
  void print(ostream& f, TGMTable* table, int numCycles,const vector<string> &topologies,int max) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);
    f << "Gene " << number << ":" << endl;
    f << "  numTrees = " << numTrees << endl;
    f << "  index "
      << setw(max) << left << "topology" << right
      << setw(10) << "single"
      << setw(10) << "joint" << endl;
    for(int i=0;i<numTrees;i++) {
        f << setw(7) << indices[i] << " "
        << setw(max) << left << topologies[indices[i]].substr(0,max) << right
        << setw(10) << setprecision(6) << counts[i] / total
        << setw(10) << setprecision(6) << (double) table->getCounts(indices[i], number) / (double) numCycles << endl;
    }
    f << endl;
  }
private:
  int number;
  int numTrees;
  double total;
  vector<int> indices;      // length is the number of topologies with positive probability
  vector<double> counts;    // length is the number of topologies with positive probability
  vector<double> logProbs;     // length is the number of trees in the whole data set
                            // added to make the retrieval of probability information very fast (at the cost of space)
  vector<bool> hasNonZeroProb;  // added to check if probs are valid. If counts[i] = 0, probs[i]=log(0) and isCountZero[i] = true
  Alias *alias;
};

class State {
public:
  State(double a,int numTaxa,int numTrees,vector<Gene*> g,bool useIP,Rand& rand) {
    alpha = a;
    logAlpha = log(alpha);
    numTreesSampled = numTrees;
    logNumTopologies = 0;
    for(int i=4;i<=numTaxa;i++)
      logNumTopologies += log(double(2.0*i-5));
    genes = g;
    alphaOverTop = exp(logAlpha - logNumTopologies);
    useIndependencePrior = useIP;
    numGroups = 0;
    counts.resize(numTrees);
    for(int i=0;i<numTrees;i++)
      counts[i] = 0;
    for(int i=0;i<genes.size();i++) {
      int top;
      top = genes[i]->pickTreeFast(rand);
      if(counts[top]==0) {
	numGroups++;
	indices.push_back(top);
      }
      counts[top]++;
      tops.push_back(top);
    }
    int numDataGenes = genes.size();

    sort(indices.begin(),indices.end());
    calculateLogPriorProb();
    calculateLogPosteriorProbProduct();
  }
  double getAlpha() { return alpha; }
  int getNumGroups() { return numGroups; }
  double getLogNumTopologies() { return logNumTopologies; }
  double getAlphaOverTop() { return alphaOverTop; }
  int getCounts(int i) { return counts[i]; }
  int getTops(int i) { return tops[i]; }
  void setAlpha(double a) {
    alpha = a;
    logAlpha = log(a);
    alphaOverTop = exp(logAlpha - logNumTopologies);
    calculateLogPriorProb();
  }
  void setTop(int i,int newTop) {
    int oldTop = tops[i];
    int n1 = counts[oldTop];
    int n2 = counts[newTop];
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
    double logHR=0;
    if(!useIndependencePrior) { // logHR=0 if independence prior
      logHR += log(double(n2 + alphaOverTop));
      logHR -= log(double(n1-1 + alphaOverTop));
    }
    logPriorProb += logHR;
    if(i < genes.size())
      logPosteriorProbProduct += (genes[i]->getLogProb(newTop) - genes[i]->getLogProb(oldTop));

  }
  double getLogPriorProb() { return logPriorProb; }
  double calculateLogPriorProb();
  double getLogPosteriorProbProduct() { return logPosteriorProbProduct; }
  double calculateLogPosteriorProbProduct() {
    logPosteriorProbProduct = 0;
    for(int i=0;i<genes.size();i++)
      logPosteriorProbProduct += genes[i]->getLogProb(tops[i]);
    return logPosteriorProbProduct;
  }
  double getLogH() { return logPriorProb + logPosteriorProbProduct; }
  int updateOne(int,Rand&);
  int update(Rand&);
  void sample(ostream& f) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);
    f << setw(6) << numGroups << setw(12) << setprecision(6) << getLogH();
    for(int i=0;i<genes.size();i++)
      f << setw(4) << tops[i];
    f << endl;
  }

  
  void updateTable(TGMTable* &table) {
    for(int j=0;j<genes.size();j++)
      table->addGeneCount(tops[j], j);
  }
  void print(ostream&);
  void testPrintState(){
      cout << "alpha= " << alpha << "; "
           << "logAlpha = " << logAlpha << "; "
           << "numGroups = " << numGroups << "; "
           << "logPriorProb = " << logPriorProb << "; "
           << "posteriorProb = " << logPosteriorProbProduct << "\n";
  }
  void updateSplits(vector<vector<int> >&,vector<vector<int> >&);
  int updateOneGroup(int,Rand&);
  void updatePairCounts(vector<vector<int> >&);
 private:
  double alpha,logAlpha,alphaOverTop;
  bool useIndependencePrior;
  vector<Gene*> genes; // vector of pointers to genes
  int numGroups,numTreesSampled;
  double logNumTopologies;
  vector<int> indices; // indices[i] = index of ith positive entry in counts
  vector<int> counts; // counts[i] = number of genes with topology i
  vector<int> tops; // tops[i] = topology for gene i
  double logPriorProb;
  double logPosteriorProbProduct;
};

class Edge;

class Node {
public:
  Node(int n,int lf) : number(n) {
    if(lf)
      leaf = true;
    else
      leaf = false;
  }
  ~Node() {}
  int getNumber() const { return number; }
  Edge* getEdge(int n) const { return edges[n]; }
  bool isLeaf() const { return leaf; }
  TaxonSet getTaxa(int i) const { return taxa[i]; }
  void setEdge(int n,Edge *e) { edges[n]=e; }
  void setTaxa(int i,TaxonSet x) { taxa[i] = x; }
  TaxonSet setAllTaxa(Node*,TaxonSet);
  void print(ostream&) const;
  Node* getNeighbor(int) const;
private:
  int number;
  TaxonSet taxa[3]; // integer representing bit vector of taxa in each of three subtrees for internal nodes
  bool leaf;
  Edge* edges[3];
};

class Edge {
public:
  Edge(int n) : number(n) {}
  int getNumber() const { return number; }
  Node* getNode(int i) { return nodes[i]; }
  Node* getOtherNode(const Node *n) const { return (n==nodes[0] ? nodes[1] : nodes[0]); }
  TaxonSet getSplit() { return split; }
  void setSplit(TaxonSet s) {split = s; }
  void setNode(int i,Node *n) { nodes[i] = n; }
  void print(ostream&,int) const;
  ~Edge() {}
private:
  int number;
  TaxonSet split;
  Node* nodes[2];
};

class SplitSet;

class Tree {
public:
  Tree(TaxonSet allTaxa,string top) : all(allTaxa) {
    numTaxa = all.getNumTaxa();
    numNodes = 2*numTaxa - 2;
    numEdges = 2 * numTaxa - 3;
    numSplits = numTaxa - 3;
    nodes.resize(numNodes);
    edges.resize(numEdges);
    for(int i=0;i<numTaxa;i++)
      nodes[i] = new Node(i,1);
    for(int i=numTaxa;i<numNodes;i++)
      nodes[i] = new Node(i,0);
    for(int i=0;i<numEdges;i++)
      edges[i] = new Edge(i);
    construct(top);
    setAllTaxa();
  }
  ~Tree() {
    for(int i=0;i<numNodes;i++)
      delete nodes[i];
    for(int i=0;i<numEdges;i++)
      delete edges[i];
  }
  void construct(string&);
  Node* connectInt(istream&,int&,int&,int&);
  Node* connectThreeNodes(Node*,Node*,Node*,int&,int&);
  void print(ostream&) const;
  void setAllTaxa();
  void getSplits(SplitSet&);
  int getNumTaxa() const { return numTaxa; }
  int getNumNodes() const { return numNodes; }
  int getNumEdges() const { return numEdges; }
  int getNumSplits() const { return numSplits; }
  TaxonSet getAll() const { return all; }
  Node* getNode(int i) const { return nodes[i]; }
  Edge* getEdge(int i) const { return edges[i]; }
  string getTop() const { return top; }
private:
  int numTaxa,numNodes,numEdges,numSplits;
  TaxonSet all;
  vector<Node *> nodes;
  vector<Edge *> edges;
  string top;
  Node *root;
};

class SplitSet {
 public:
  SplitSet() {}
  void print(ostream& f) {
    for(int i=0;i<splits.size();i++) {
      splits[i].print(f);
      f << endl;
    }
  }
  void printShort(ostream& f) {
    for(int i=0;i<splits.size();i++) {
      TaxonSet t = splits[i];
      f << setw(5) << t << " ";
    }
  }
  int getNumTaxa() { return numTaxa; }
  TaxonSet getSplit(int i) { return splits[i]; }
  void setNumTaxa(int x) { numTaxa = x; }
  void setAll(TaxonSet a) { all = a; }
  void setSplit(int i,TaxonSet c) {
    splits[i] = c;
  }
  void copySplits(Tree& t) { t.getSplits(*this); }
  void resize(int n) { splits.resize(n); }
 private:
  TaxonSet all;
  int numTaxa;
  vector<TaxonSet> splits;
};

class ConcordanceEdge;

class ConcordanceNode {
 public:
  ConcordanceNode(int numTaxa,int i,int numGenes) {
    taxonSet = TaxonSet::getTaxaWithInitialValue(i, numGenes);
    leaf = true;
    number = minTaxa = i;
  }
  ConcordanceNode(TaxonSet ts,int num) {
    taxonSet = ts;
    leaf = ts.leaf();
    number = num;
  }
  ~ConcordanceNode() { edges.clear(); }
  TaxonSet getTaxonSet() { return taxonSet; }
  int getNumber() { return number; }
  bool all() { return taxonSet.all(); }
  void addEdge(ConcordanceEdge* e) { edges.push_back(e); }
  int getDegree() { return edges.size(); }
  ConcordanceEdge* getEdge(int i) { return edges[i]; }
  ConcordanceNode* getNeighbor(ConcordanceEdge*);
  void printRoot(ostream&, int);
  void print(ostream&,ConcordanceEdge*, int);
  void printTopologyRoot(ostream&);
  void printTopology(ostream&,ConcordanceEdge*);
  int setMinTaxa(ConcordanceEdge*);
  void setMinTaxaRoot();
  int getMinTaxa() { return minTaxa; }
 private:
  int number;
  int minTaxa;
  TaxonSet taxonSet;
  bool leaf;
  vector<ConcordanceEdge*> edges;
};

class ConcordanceEdge {
 public:
  ConcordanceEdge() {}
  ConcordanceEdge(ConcordanceNode* child,ConcordanceNode* parent,int num) {
    nodes[0] = child;
    nodes[1] = parent;
    number = num;
  }
  ~ConcordanceEdge() {}
  void addNodes(ConcordanceNode* child,ConcordanceNode* parent) {
    nodes[0] = child;
    nodes[1] = parent;
  }
  ConcordanceNode* getChild() { return nodes[0]; }
  ConcordanceNode* getParent() { return nodes[1]; }
  ConcordanceNode* getOtherNode(ConcordanceNode* n) {
    if(n == nodes[0])
      return nodes[1];
    if(n == nodes[1])
      return nodes[0];
    cerr << "Internal Error: getOtherNode called on edge " << number << " and node " << n->getNumber() << endl;
    exit(1);
  }
  int getNumber() { return number; }
  void setNumber(int x) { number = x; }
 private:
  int number;
  ConcordanceNode* nodes[2];
};

class ConcordanceTree {
 public:
  ConcordanceTree() {}
  ConcordanceTree(vector<TaxonSet>,int);
  ~ConcordanceTree() {
    nodes.clear();
    edges.clear();
  }
  void setMinTaxa();
  void print(ostream&, int);
  void printFull(ostream&);
  void printTopology(ostream&);
 private:
  int numTaxa,numNodes,numEdges;
  vector<ConcordanceNode*> nodes;
  vector<ConcordanceEdge*> edges;
};

ConcordanceNode* ConcordanceNode::getNeighbor(ConcordanceEdge* e) {
  for(vector<ConcordanceEdge*>::iterator p = edges.begin();p != edges.end(); p++)
    if((*p) == e)
      return e->getOtherNode(this);
  cerr << "Internal Error: called getNeighbor on node " << number << " and edge " << e->getNumber() << endl;
  exit(1);
}

class ConcordanceEdgeCompare : public binary_function<ConcordanceEdge*,ConcordanceEdge*,bool> {
public:
  ConcordanceEdgeCompare(ConcordanceNode* n) {
    node = n;
  }

  bool operator()(ConcordanceEdge* x,ConcordanceEdge* y) {
    return (x->getOtherNode(node)->getMinTaxa() < y->getOtherNode(node)->getMinTaxa());
  }

private:
  ConcordanceNode* node;
};

int ConcordanceNode::setMinTaxa(ConcordanceEdge* parent) {
  minTaxa = number;
  if(!leaf) {
    for(vector<ConcordanceEdge*>::iterator e=edges.begin();e!=edges.end();e++)
      if(*e != parent) {
	int m = getNeighbor(*e)->setMinTaxa(*e);
	if(m < minTaxa)
	  minTaxa = m;
      }
  }
  sort(edges.begin(),edges.end(),ConcordanceEdgeCompare(this));
  return minTaxa;
}

void ConcordanceNode::setMinTaxaRoot() {
  minTaxa = number;
  for(vector<ConcordanceEdge*>::iterator e=edges.begin();e!=edges.end();e++) {
    int m = getNeighbor(*e)->setMinTaxa(*e);
    if(m < minTaxa)
      minTaxa = m;
  }
  sort(edges.begin(),edges.end(),ConcordanceEdgeCompare(this));
}

void ConcordanceTree::setMinTaxa() {
  nodes.back()->setMinTaxaRoot();
}

// Create a concordance tree from vector<TaxonSet>
//   assumes that taxon sets have the same number of taxa
//   will add taxon sets corresponding to leaves
//   assumes that the taxon set is partially sorted so subsets precede supersets
//   assumes that taxon sets for splits do not include the outgroup
//     except for the last taxon set in the list which is all of the taxa (and the root of the concordance tree)

ConcordanceTree::ConcordanceTree(vector<TaxonSet> ts,int numGenes) {
  numTaxa = ts.front().getNumTaxa();
  // create ConcordanceNodes for the leaves
  for(int i=0;i<numTaxa;i++)
    nodes.push_back(new ConcordanceNode(numTaxa,i,numGenes));
  // create ConcordanceNodes for the internal nodes
  int nextNode = numTaxa;
  for(vector<TaxonSet>::iterator p=ts.begin();p!=ts.end();p++)
    nodes.push_back(new ConcordanceNode(*p,numTaxa++));
  int nextEdge = 0;
  for(vector<ConcordanceNode*>::iterator c=nodes.begin();c!=nodes.end() && !(*c)->all();c++) {
    vector<ConcordanceNode*>::iterator p = c+1;
    TaxonSet tset = (*c)->getTaxonSet();
    while(p!= nodes.end() && !(tset.isSubsetOf( (*p)->getTaxonSet() )))
      p++;
    if(p==nodes.end()) {
      cerr << "InternalError: taxon has no superset in ConcordanceTree::ConcordanceTree(vector<TaxonSet>);" << endl;
      exit(1);
    }
    edges.push_back(new ConcordanceEdge(*c,*p,nextEdge++));
    (*c)->addEdge(edges.back());
    (*p)->addEdge(edges.back());
  }
  setMinTaxa();
}

void ConcordanceTree::printFull(ostream& f) {
  f << "Nodes:" << endl << endl;
  f << "  Number Degree Subset (Edge,Neighbor)" << endl;
  for(vector<ConcordanceNode*>::iterator n=nodes.begin();n!=nodes.end();n++) {
    f << setw(8) << (*n)->getNumber() << setw(7) << (*n)->getDegree() << " ";
    TaxonSet t = (*n)->getTaxonSet();
    f << t;
    for(int i=0;i<(*n)->getDegree();i++) {
      ConcordanceEdge* e = (*n)->getEdge(i);
      f << " (" << e->getNumber() << "," << e->getOtherNode(*n)->getNumber() << ")";
    }
    f << endl;
  }
  f << endl << endl << "Edges:" << endl << endl;
  f << "  Number Child Parent" << endl;
  for(vector<ConcordanceEdge*>::iterator e=edges.begin();e!=edges.end();e++)
    f << setw(8) << (*e)->getNumber() << setw(6) << (*e)->getChild()->getNumber() << setw(7) << (*e)->getParent()->getNumber() << endl;
  f << endl << endl;
}

void ConcordanceNode::printTopology(ostream& f,ConcordanceEdge* parentEdge) {
  if(leaf)
    f << number + 1;
  else {
    f << "(";
    int numChildren = edges.size() - 1;
    for(vector<ConcordanceEdge*>::iterator e=edges.begin(); e!=edges.end();e++)
      if(*e != parentEdge) {
	getNeighbor(*e)->printTopology(f,*e);
	if(--numChildren > 0)
	  f << ",";
      }
    f << ")";
  }
}

void ConcordanceNode::print(ostream& f,ConcordanceEdge* parentEdge, int numGenes) {
  if(leaf)
    f << number + 1;
  else {
    f << "(";
    int numChildren = edges.size() - 1;
    for(vector<ConcordanceEdge*>::iterator e=edges.begin(); e!=edges.end();e++)
      if(*e != parentEdge) {
	getNeighbor(*e)->print(f,*e, numGenes);
	if(--numChildren > 0)
	  f << ",";
      }
    f << ")";
  }
  f << ":" << setprecision(3) << taxonSet.getWeight() / numGenes;
}

void ConcordanceNode::printTopologyRoot(ostream& f) {
  f << "(";
  for(vector<ConcordanceEdge*>::iterator e=edges.begin(); e!=edges.end();e++) {
    getNeighbor(*e)->printTopology(f,*e);
    if(*e != edges.back())
      f << ",";
  }
  f << ");" << endl;
}

void ConcordanceNode::printRoot(ostream& f, int numGenes) {
  f << "(";
  for(vector<ConcordanceEdge*>::iterator e=edges.begin(); e!=edges.end();e++) {
    getNeighbor(*e)->print(f,*e, numGenes);
    if(*e != edges.back())
      f << ",";
  }
  f << ");" << endl;
}

void ConcordanceTree::printTopology(ostream& f) {
  f.setf(ios::fixed, ios::floatfield);
  f.setf(ios::showpoint);
  nodes.back()->printTopologyRoot(f);
  f << endl;
}

void ConcordanceTree::print(ostream& f, int numGenes) {
  f.setf(ios::fixed, ios::floatfield);
  f.setf(ios::showpoint);
  nodes.back()->printRoot(f, numGenes);
  f << endl;
}

class FileNames {
public:
  FileNames(string fileRoot) :
    rootFileName(fileRoot),
    clusterFile(fileRoot + ".cluster"),
    concordanceFile(fileRoot + ".concordance"),
    //    geneFile(fileRoot + ".gene"),
    //    genePosteriorFile(fileRoot + ".genepost"),
    genePosteriorFile(fileRoot + ".gene"),
    jointFile(fileRoot + ".joint"),
    inputFile(fileRoot + ".input"),
    outFile(fileRoot + ".out"),
    pairTreeFile(fileRoot + ".pairs"),
    sampleFile(fileRoot + ".sample"),
    singleFile(fileRoot + ".single")
    //    splitsFile(fileRoot + ".splits"),
    //    topologyFile(fileRoot + ".top"),
    //    treePosteriorFile(fileRoot + ".topologies")
  {}
  void setFileNames(string fileRoot) {
    rootFileName = fileRoot;
    clusterFile = fileRoot + ".cluster";
    concordanceFile = fileRoot + ".concordance";
    //    geneFile = fileRoot + ".gene";
    //    genePosteriorFile = fileRoot + ".genepost";
    genePosteriorFile = fileRoot + ".gene";
    jointFile = fileRoot + ".joint";
    inputFile = fileRoot + ".input";
    outFile = fileRoot + ".out";
    pairTreeFile = fileRoot + ".pairs";
    sampleFile = fileRoot + ".sample";
    singleFile = fileRoot + ".single";
    //    splitsFile = fileRoot + ".splits";
    //    topologyFile = fileRoot + ".top";
    //    treePosteriorFile = fileRoot + ".topologies";
  }
  string getRootFileName() { return rootFileName; }
  string getClusterFile() { return clusterFile; }
  string getConcordanceFile() { return concordanceFile; }
  //  string getGeneFile() { return geneFile; }
  string getGenePosteriorFile() { return genePosteriorFile; }
  string getJointFile() { return jointFile; }
  string getInputFile() { return inputFile; }
  string getOutFile() { return outFile; }
  string getPairTreeFile() { return pairTreeFile; }
  string getSampleFile() { return sampleFile; }
  string getSampleFile(unsigned int k) {
    ostringstream run_k;
    run_k << k;
    return(sampleFile + ".run" + run_k.str());
  }
  string getSingleFile() { return singleFile; }
  //  string getSplitsFile() { return splitsFile; }
  //  string getTopologyFile() { return topologyFile; }
  //  string getTreePosteriorFile() { return treePosteriorFile; }
  string getInputListFileName() { return inputListFileName; }
  void setInputListFileName(string filename) { inputListFileName = filename; }
private:
  string rootFileName;
  string clusterFile;
  string concordanceFile;
  //  string geneFile;
  string genePosteriorFile;
  string jointFile;
  string inputFile;
  string outFile;
  string pairTreeFile;
  string sampleFile;
  string singleFile;
  //  string splitsFile;
  //  string topologyFile;
  //  string treePosteriorFile;
  string inputListFileName;
};

class ModelParameters {
public:
  ModelParameters() : alpha(1.0), useIndependencePrior(false) {}
  double getAlpha() { return alpha; }
  void setAlpha(double a) { alpha = a; }
  bool getUseIndependencePrior() { return useIndependencePrior; }
  void setUseIndependencePrior(bool iprior) { useIndependencePrior = iprior; }
private:
  double alpha;
  bool useIndependencePrior;
};

class RunParameters {
public:
  RunParameters() :
    alphaMultiplier(10.0),
    seed1(1234),
    seed2(5678),
    numUpdates(100000),
    subsampleRate(1),
    numRuns(2),
    numChains(1),
    mcmcmcRate(100),
    calculatePairs(false),
    useUpdateGroups(true),
    createSampleFile(false),
    createJointFile(false),
    createSingleFile(false),
    numGenomewideGrid(1000),
    swCFcutoff(.05),
    optSpace(false),
    shouldPruneGene(false),
    buildPopulationTree(true)
  {}
  double getAlphaMultiplier() { return alphaMultiplier; }
  void setAlphaMultiplier(double x) { alphaMultiplier = x; }
  unsigned int getSeed1() { return seed1; }
  void setSeed1(unsigned int x) { seed1 = x; }
  unsigned int getSeed2() { return seed2; }
  void setSeed2(unsigned int x) { seed2 = x; }
  unsigned int getNumUpdates() { return numUpdates; }
  void setNumUpdates(unsigned int x) { numUpdates = x; }
  unsigned int getSubsampleRate() { return subsampleRate; }
  void setSubsampleRate(unsigned int x) { subsampleRate = x; }
  unsigned int getNumRuns() { return numRuns; }
  void setNumRuns(unsigned int x) { numRuns = x; }
  unsigned int getNumChains() { return numChains; }
  void setNumChains(unsigned int x) { numChains = x; }
  unsigned int getMCMCMCRate() { return mcmcmcRate; }
  void setMCMCMCRate(unsigned int x) { mcmcmcRate = x; }
  unsigned int getNumGenomewideGrid(){ return numGenomewideGrid;}
  void setNumGenomewideGrid(unsigned int x){ numGenomewideGrid = x;}
  bool getCalculatePairs() { return calculatePairs; }
  void setCalculatePairs(bool x) { calculatePairs = x; }
  bool getUseUpdateGroups() { return useUpdateGroups; }
  void setUseUpdateGroups(bool x) { useUpdateGroups = x; }
  bool getCreateSampleFile() { return createSampleFile; }
  void setCreateSampleFile(bool x) { createSampleFile = x; }
  bool getCreateJointFile() { return createJointFile; }
  void setCreateJointFile(bool x) { createJointFile = x; }
  bool getCreateSingleFile() { return createSingleFile; }
  void setCreateSingleFile(bool x) { createSingleFile = x; }
  double getSwCFcutoff() { return  swCFcutoff; }
  void   setSwCFcutoff(double x) { swCFcutoff = x; }
  string getPruneFile() { return pruneFile; }
  void setPruneFile(string pFile) { pruneFile = pFile; }
  bool getPruneGene() { return shouldPruneGene; }
  void setPruneGene(bool pg) { shouldPruneGene = pg; }
  bool shouldOptSpace() { return optSpace;}
  void setOptSpace(bool o) { optSpace = o; }
  bool shouldBuildPopulationTree() { return buildPopulationTree;}
  void setBuildPopulationTree(bool x) { buildPopulationTree = x; }
  /*void MPI_calcMyChains(int procs, int rank){
    int total = numChains*numRuns; 
    int myChains; 
    if (procs >= total){ //If more procs than total chains
        if (total)
        if (rank > procs-total-1){
            myChains = 0; 
        }else if (rank == total-(procs-1)){
            myChains = total%procs;
        }else{
            myChains = total/procs;
        }
    }else{ //If more chains than procs
        
    }
  }*/
     
private:
  double alphaMultiplier,swCFcutoff;  // cutoff on sample-wide CF to display splits
  unsigned int seed1,seed2,numUpdates,subsampleRate,numRuns,numChains,mcmcmcRate,numGenomewideGrid;
  bool calculatePairs;
  bool useUpdateGroups;
  bool createSampleFile;
  bool createJointFile;
  bool createSingleFile;
  string pruneFile;
  bool shouldPruneGene;
  bool optSpace; // set this if gene topology counts forms a sparse matrix. We do not store zeros in the 'table' in this case.
                 // default is false because this is not time-efficient.
  bool buildPopulationTree;
};

class Defaults {
public:
  Defaults(FileNames &fn,ModelParameters &mp,RunParameters &rp) :
    alpha(mp.getAlpha()),
    numUpdates(rp.getNumUpdates()),
    numRuns(rp.getNumRuns()),
    numChains(rp.getNumChains()),
    mcmcmcRate(rp.getMCMCMCRate()),
    alphaMultiplier(rp.getAlphaMultiplier()),
    subsampleRate(rp.getSubsampleRate()),
    rootFileName(fn.getRootFileName()),
    inputListFileName(fn.getInputListFileName()),
    seed1(rp.getSeed1()),
    seed2(rp.getSeed2()),
    useIndependencePrior(mp.getUseIndependencePrior()),
    calculatePairs(rp.getCalculatePairs()),
    useUpdateGroups(rp.getUseUpdateGroups()),
    createSampleFile(rp.getCreateSampleFile()),
    createJointFile(rp.getCreateJointFile()),
    createSingleFile(rp.getCreateSingleFile()),
    numGenomewideGrid(rp.getNumGenomewideGrid()),
    swCFcutoff(rp.getSwCFcutoff()),
    optSpace(rp.shouldOptSpace())
  {}
  void print(ostream&);
  double getAlpha() { return alpha; }
  unsigned int getNumUpdates() { return numUpdates; }
  unsigned int getNumRuns() { return numRuns; }
  unsigned int getNumChains() { return numChains; }
  unsigned int getMCMCMCRate() { return mcmcmcRate; }
  double getAlphaMultiplier() { return alphaMultiplier; }
  unsigned int getSubsampleRate() { return subsampleRate; }
  string getRootFileName() { return rootFileName; }
  string getInputListFileName() { return inputListFileName; }
  unsigned int getSeed1() { return seed1; }
  unsigned int getSeed2() { return seed2; }
  unsigned int getNumGenomewideGrid() { return numGenomewideGrid; }
  bool getUseIndependencePrior() { return useIndependencePrior; }
  bool getCalculatePairs() { return calculatePairs; }
  bool getUseUpdateGroups() { return useUpdateGroups; }
  bool getCreateSampleFile() { return createSampleFile; }
  bool getCreateJointFile() { return createJointFile; }
  bool getCreateSingleFile() { return createSingleFile; }
  double getSwCFcutoff() { return swCFcutoff; }
  bool shouldOptSpace() { return optSpace;}
  void setOptSpace(bool o) { optSpace = o; }
private:
  double alpha;
  unsigned int numUpdates;
  unsigned int numChains, numRuns;
  unsigned int mcmcmcRate;
  double alphaMultiplier, swCFcutoff;
  unsigned int subsampleRate;
  string rootFileName;
  string inputListFileName;
  unsigned int seed1;
  unsigned int seed2;
  unsigned int numGenomewideGrid;
  bool useIndependencePrior;
  bool calculatePairs;
  bool useUpdateGroups;
  bool createSampleFile;
  bool createJointFile;
  bool createSingleFile;
  bool optSpace;
};

class GenomewideDistribution {
 public:
  GenomewideDistribution(int ngene, unsigned int ngrid){
    samplewide.resize(ngene+1);
    genomewide.resize(ngrid+1);
    convolutionWeight.resize(ngrid+1);
    for (int i=0; i<convolutionWeight.size(); i++)
      convolutionWeight[i].resize(ngene+1);
    samplewideCredibilityInterval.resize(6);
    genomewideCredibilityInterval.resize(6);
  }
  ~GenomewideDistribution() {
    samplewide.clear();
    genomewide.clear();
    convolutionWeight.clear();
    samplewideCredibilityInterval.clear();
    genomewideCredibilityInterval.clear();
  }
  double getPriorProbability(){return priorProbability;}
  void updatePriorProbability(TaxonSet s){
    s.updatePriorProbability();
    priorProbability = s.getPriorProbability();
    samplewidePosteriorMean = s.getWeight()/(samplewide.size()-1.0);
  };
  void updateConvolutionWeight(double);
  void updateSamplewide(vector<double> &);
  void updateSamplewide(vector<int> &, unsigned int);
  void setSamplewidePosteriorMeanSD(double x){ samplewidePosteriorMeanSD = x;}
  double getSamplewidePosteriorMeanSD(){ return samplewidePosteriorMeanSD;}
  void updateGenomewide(double);
  void printSampleCF(ostream&);
  void printGenomeCF(ostream&);
  void printfull();
 private:
  double priorProbability;
  vector<vector<double> > convolutionWeight;
  vector<double> samplewide;
  vector<double> genomewide;
  double genomewidePosteriorMean,samplewidePosteriorMean, samplewidePosteriorMeanSD;
  vector<double> samplewideCredibilityInterval, genomewideCredibilityInterval;
  // vectors of size 6: lo(99%), lo(95%), lo(90%), hi(90%), hi(95%), hi(99%).
};


void GenomewideDistribution::printfull(){
  cerr << "Prior Probability:  " << priorProbability << endl;
  cerr << "Sample post. mean:  " << samplewidePosteriorMean << endl;
  cerr << "99% interval:       ("<< samplewideCredibilityInterval[0]<<","<<samplewideCredibilityInterval[5]<<")"<<endl;
  cerr << "95% interval:       ("<< samplewideCredibilityInterval[1]<<","<<samplewideCredibilityInterval[4]<<")"<<endl;
  cerr << "90% interval:       ("<< samplewideCredibilityInterval[2]<<","<<samplewideCredibilityInterval[3]<<")"<<endl;
  cerr << "Genome post. mean:  " << genomewidePosteriorMean << endl;
  cerr << "99% interval:       ("<< genomewideCredibilityInterval[0]<<","<<genomewideCredibilityInterval[5]<<")"<<endl;
  cerr << "95% interval:       ("<< genomewideCredibilityInterval[1]<<","<<genomewideCredibilityInterval[4]<<")"<<endl;
  cerr << "90% interval:       ("<< genomewideCredibilityInterval[2]<<","<<genomewideCredibilityInterval[3]<<")"<<endl;
  cerr << "Sample wide distr.:  ";
  for (int j=0;j<samplewide.size();j++) cerr << setw(8) << setprecision(2) << j;
  cerr << endl << "                     ";
  for (int j=0;j<samplewide.size();j++) cerr << setw(8) << setprecision(2) << samplewide[j];
  cerr << endl;

  cerr << "Genome wide distr.:  ";
  for (int i=0;i<genomewide.size();i++) cerr << setw(8) << setprecision(2) << i;
  cerr << endl << "                     ";
  for (int i=0;i<genomewide.size();i++) cerr << setw(8) << setprecision(2) << genomewide[i];
  cerr << endl;

  cerr << "Convolution Weights (first and last 10 lines):\n";
  for (int i=0;i<convolutionWeight.size();i++){
    if (i<10 || i>convolutionWeight.size()-11){
      for (int j=0;j<samplewide.size();j++)
	cerr << setw(8) << setprecision(2) << convolutionWeight[i][j] ;
      cerr << endl;
    }
    if (i==10) cerr<<"..."<<endl;
  }
  cerr << endl;
}
void GenomewideDistribution::printSampleCF(ostream& f){ // 95% credibility only
  f << " " << setprecision(3) << samplewidePosteriorMean
    <<"("<<samplewideCredibilityInterval[1]<<","<<samplewideCredibilityInterval[4]<<")";
}
void GenomewideDistribution::printGenomeCF(ostream& f){
  f << " " << setprecision(3) << genomewidePosteriorMean
    <<"("<<genomewideCredibilityInterval[1]<<","<<genomewideCredibilityInterval[4]<<")";
}


class TaxonList {
 public:
  TaxonList(vector<string>& tTable){
    Nall = tTable.size();
    for (int i=0; i<Nall; i++)
      name.push_back(tTable[i]);
    alleleOf.resize(Nall); // for now only one allele allowed.
    Ntax = Nall;
    isMissingOne.resize(Ntax);
    numberGenes.resize(Ntax);
  }
  ~TaxonList(){
    name.clear();
    isMissingOne.clear();
    include.clear();
    numberGenes.clear();
    alleleOf.clear();
  }
  void setNumberGenes(vector<vector<int> >& taxid){
    for (int j=0; j<Ntax; j++)
      numberGenes[j]=0;
    for (int i=0; i<taxid.size(); i++){
      for (int j=0; j<taxid[i].size(); j++){
	if (taxid[i][j]-1<numberGenes.size())
	  numberGenes[taxid[i][j]-1]++;
	else
	  cerr << "Warning: taxon ID ("<<taxid[i][j]
	       <<") exceeded the expected max Taxon number."<<endl;
      }
    }
    for (int j=0; j<Ntax; j++)
      if (numberGenes[j]<taxid.size())
	isMissingOne[j]=true;
      else
	isMissingOne[j]=false;
  }
  vector<vector<bool> > getAvailabilityMatrix(vector<vector<int> >& taxid); //fixit: for later may be
  vector<int> getNonMissing(){
    vector<int> mytax;
    for (int j=0; j<Ntax; j++)
      if (!isMissingOne[j])
	mytax.push_back(j+1);
    return mytax;
  }
  void setInclude(vector<int>& taxid){
    for (int i=0; i<Ntax; i++)
      include.push_back(false);
    for (int i=0; i<taxid.size(); i++){
      if (taxid[i]-1<include.size())
	include[taxid[i]-1]=true;
      else
	cerr << "Warning: taxon ID ("<<taxid[i]
	     <<") exceeded the expected max Taxon number ("
	     <<include.size()<<")"<<endl;
    }
  }
  void print(ostream& f){
    //f << "All taxa:"<< endl;
    //f << "\ntranslate" << endl;
    //for (int i=0; i<Ntax; i++)
    //  f << setw(4) << i+1 << " "<<name[i]<< (i<(Ntax-1)?",":";") <<endl;
    //f << "Taxa with missing data:"<< endl;
    //for (int i=0; i<Ntax; i++)
    //  if (isMissingOne[i])
    //	f << setw(4) << i+1 <<endl;
    f << "Number of genes sequenced for each taxon:"<<endl;
    for (int i=0; i<Ntax; i++)
      f << setw(4) << i+1 << " "<<name[i]<<" "<< setw(6)<<numberGenes[i]<<endl;
    //f << "Taxa included in the analysis:\n";
    //for (int i=0; i<Ntax; i++)
    //  if (include[i])
    //	f << setw(4) << i+1 <<endl;
  }
 private:
  int Ntax;             // total number of individuals
  int Nall;             // total number of sequences/alleles
  vector<string> name;  // taxon or allele names
  vector<bool> isMissingOne; // is the taxon missing for at least one gene?
  vector<bool> include;      // should the taxon be included in the analysis?
  vector<int> numberGenes;   // how many genes have a sequence for this taxon?
  vector<int> alleleOf;      // if not an allele: 0. Otherwise: ID of the individual
};

// Bucky Pruner
class BuckyPruner : public mbsumtree::Pruner {
public:
  BuckyPruner(int maxTaxa) {
    mTaxa = maxTaxa;
  }

  virtual bool prune(int num) {
    return num <= 0 || num > mTaxa;
  }
private:
  int mTaxa;
};
#endif
