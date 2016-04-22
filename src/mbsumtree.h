#ifndef MBSUMTREE_H_
#define MBSUMTREE_H_

using namespace std;
namespace mbsumtree {
class Node;

class Edge {
public:
  Edge() {}
  ~Edge() {
    for(int i=0;i<2;i++)
      nodes[i] = NULL;
  }
  Edge(int n,double x) :
    number(n),
    length(x)
    {}
  int getNumber() { return number; }
  void setLength(double x) { length=x; }
  double getLength() { return length; }
  Node* getOtherNode(Node *n) { return (nodes[0] == n ? nodes[1] : nodes[0]); }
  void setNode(Node *n,int i) { nodes[i]=n; }
  void setNodes(Node *m,Node *n) { nodes[0]=m; nodes[1]=n; }
  Node* getNode(int n) { return nodes[n]; }
  void swapNodes() {
    Node *temp=nodes[0];
    nodes[0] = nodes[1];
    nodes[1] = temp;
  }
private:
  int number;
  double length;
  Node* nodes[2];
};

class Node {
public:
  Node() { number = -1; }
  ~Node() { edges.clear(); }
  Node(int n,bool l) : number(n), leaf(l) {}
  int getNumber() const { return number; }
  void setNumber(int x) { number = x; }
  bool isLeaf() { return leaf; }
  void setLeaf(bool x) { leaf = x; }
  int getNumEdges() { return edges.size(); }
  Edge* getEdge(int n) { return edges[n]; }
  void setEdge(Edge* e) { edges.push_back(e); }
  void setEdge(Edge *e,int n) { edges[n] = e; }
  void removeEdge(Edge* e) {
    for (vector<Edge*>::iterator itr = edges.begin(); itr != edges.end();)
      if (*itr == e) {
        edges.erase(itr);
        break;
      }
      else {
        itr++;
      }
    }
  void clearEdges() { edges.clear(); }
  void reorder(Edge*);
  void rootReorder();
  int getMinTaxa() { return minTaxa; }
  int setMinTaxa(int x) { minTaxa = x; return x; }
  int setMinTaxa(Edge*);
  void setNumbers(int&,Edge*);
  void printTop(ostream&,Edge*);
private:
  int number;
  bool leaf;
  vector<Edge*> edges;
  int minTaxa;
};

class EdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  EdgeCompare(Edge* e,Node* n) : parent(e),node(n) {}

  bool operator()(Edge* x,Edge* y) {
    if(x==parent)
      return false;
    if(y==parent)
      return true;
    return (x->getOtherNode(node)->getMinTaxa() < y->getOtherNode(node)->getMinTaxa());
  }
private:
  const Edge* parent;
  Node* node;
};

class RootEdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  RootEdgeCompare(Node* r) : root(r) {}

  bool operator()(Edge* x,Edge* y) {
    return (x->getOtherNode(root)->getMinTaxa() < y->getOtherNode(root)->getMinTaxa());
  }
private:
  Node* root;
};

class Pruner {
public:
  virtual bool prune(int num) {
    return false;
  }
};

class Tree {
public:
  Tree(string,int,Pruner*);
  ~Tree() {
    for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
      delete *e;
    for(vector<Node*>::iterator n=nodes.begin();n!=nodes.end();n++)
      delete *n;
    nodes.clear();
    edges.clear();
  }
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return numNodes; }
  int getNumEdges() { return numEdges; }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  void printTop(ostream&);
  void reorder(Node*);
  bool readSubtree(istringstream&,Node*,int&,int);
  void setNumbers(Node*);
  void setMinTaxa(Node*);
private:
  int numTaxa;
  int numNodes;
  int numEdges;
  vector<Node*> nodes;
  vector<Edge*> edges;
  int maxTaxNumber;
  Pruner *thePruner;
};
}
#endif
