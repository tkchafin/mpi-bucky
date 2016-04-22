#ifndef QUARTETS_H_
#define QUARTETS_H_

namespace quartet {

class Quartets {
public:
    Quartets(string top);
    virtual ~Quartets();
private:
    vector<int> qs;
};

class SuperNode;

class TaxonSet {
public:
    TaxonSet() {
    }

    void add(int taxon) {
        tset.push_back(taxon);
    }

    void merge(TaxonSet*& other) {
        //merge two sorted ranges, result is also sorted.
        int initialSize = tset.size();
        tset.resize(initialSize + other->tset.size());
        vector<int> otherset = other->getTSet();
        copy(otherset.begin(), otherset.begin() + otherset.size(), tset.begin() + initialSize);
        inplace_merge(tset.begin(), tset.begin() + initialSize, tset.end());
    }

    vector<int>& getTSet() {
        return tset;
    }

    void remove(TaxonSet*& other) {
        for (vector<int>::iterator itr = other->tset.begin(); itr != other->tset.end(); itr++) {
            vector<int>::iterator pos = find(tset.begin(), tset.end(), *itr);
            if (pos != tset.end())
                tset.erase(pos);
        }
    }

    int intersect(TaxonSet*& other) {
        vector<int>& t2 = other->tset;
        vector<int> v(tset.size() + t2.size());
        vector<int>::iterator it = set_intersection(tset.begin(), tset.end(), t2.begin(), t2.end(), v.begin());
        if (it -v.begin() > 0 && tset.size() == t2.size() && t2.size() == v.size()) {
            return numeric_limits<int>::max();
        }
        return int(it - v.begin());
    }

    void printWithTaxonSet(ostream &str, TaxonSet*& other, int numTaxa);

    friend ostream& operator<<(ostream& f, TaxonSet& t) {
        for (vector<int>::iterator itr = t.tset.begin(); itr != t.tset.end(); itr++) {
            f << *itr << ",";
        }

        long pos = f.tellp();
        f.seekp(pos - 1);//remove last ','
        return f;
    }

private:
    vector<int> tset;
};

class Edge;

class Node {
public:
  Node(int n,int lf) : number(n) {
    for (int i = 0; i < 3; i++)
      t[i] = new TaxonSet();

    if(lf)
      leaf = true;
    else
      leaf = false;
  }
  ~Node() { /*delete[] t;*/}
  int getNumber() const { return number; }
  Edge* getEdge(int n) const { return edges[n]; }
  bool isLeaf() const { return leaf; }
  void setEdge(int n,Edge *e) { edges[n]=e; }
  void print(ostream&, ostream&, const Node *, const Node *, int) const;
  Node* getNeighbor(int) const;
  void mergeTaxa(int index, Node *n, int nIndex) {
      t[index]->merge(n->t[nIndex]);
  }
  void addTaxa(int index, int taxon) {
      t[index]->add(taxon);
  }

  vector<int>& getTaxa(int index) {
      return t[index]->getTSet();
  }

  TaxonSet*& getTset(int index) {
      return t[index];
  }
  void print(ostream& f) const;

private:
  int number;
  bool leaf;
  Edge* edges[3];
  TaxonSet* t[3];
};

class Edge {
public:
  Edge(int n) : number(n) {weight  = 0.0;}
  int getNumber() const { return number; }
  Node* getNode(int i) { return nodes[i]; }
  Node* getOtherNode(const Node *n) const { return (n==nodes[0] ? nodes[1] : nodes[0]); }
  void setNode(int i,Node *n) { nodes[i] = n; }
  void print(ostream&,int) const;
  void addWeight(double wt) { weight += wt; }
  double getWeight() { return weight; }
  ~Edge() {}
private:
  int number;
  double weight;
  Node* nodes[2];
};

class Tree {
public:
  Tree(int nTaxa,string t) {
    top = t;
    numTaxa = nTaxa;
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
  }
  ~Tree() {
    for(int i=0;i<numNodes;i++)
      delete nodes[i];
    for(int i=0;i<numEdges;i++)
      delete edges[i];
  }
  void print(ostream& f) const
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
  void construct(string&);
  Node* connectInt(istream&,int&,int&,int&,int&);
  Node* connectThreeNodes(Node*,Node*,Node*,int&,int&);
  void modifyOutgroup(int,string&,string&) const;
  void setAllTaxa();
  void getQuartets(vector<int>& rInd, vector<int>& cInd);
  int getQuartets(vector<int>& rInd, vector<int>& cInd, Node *n1, Node *n2, string& quartet);
  int getNumTaxa() const { return numTaxa; }
  int getNumNodes() const { return numNodes; }
  int getNumEdges() const { return numEdges; }
  int getNumSplits() const { return numSplits; }
  Node* getNode(int i) const { return nodes[i]; }
  Edge* getEdge(int i) const { return edges[i]; }
  string getTop() const { return top; }
  int intersect(vector<int>& t1, vector<int>& t2);
  void getTsets(vector<TaxonSet*>& tsets, Node *n1, Node *n2);
private:
  int numTaxa,numNodes,numEdges,numSplits;
  vector<Node *> nodes;
  vector<Edge *> edges;
  string top;
  Node *root;
};

class SuperNode {
public:
    SuperNode(int n, int lf) {
        number = n;
        tset = new TaxonSet();
        if (lf) {
            numLeafs = 1;
            stringstream str;
            str << n;
            topology = str.str();
            tset->add(n);
        }
        else {
            numLeafs = 0;
        }
    }

    void add(SuperNode *lc, SuperNode *rc) {
        numLeafs = lc->numLeafs + rc->numLeafs;
        if (lc->getLowestTaxon() < rc->getLowestTaxon()) {
            left = lc;
            right = rc;
        }
        else {
            left = rc;
            right = lc;
        }

        tset->merge(lc->tset);
        tset->merge(rc->tset);
        topology = "(" + left->topology + "," + right->topology + ")";
    }

    void setLeft(SuperNode *l) {
        left = l;
    }

    void setRight(SuperNode *r) {
        right = r;
    }

    SuperNode *getLeft() {
        return left;
    }

    SuperNode *getRight() {
        return right;
    }

    int getNumNodes() {
        return numLeafs;
    }

    int getLowestTaxon() {
        return tset->getTSet()[0];
    }

    string getTopology() {
        return topology;
    }

    TaxonSet* getTaxa() {
        return tset;
    }

    void print(ostream& f) {
        f << topology;
    }


    // if this represents (2,4) and n2 represents (1,3)
    // output will be ((1,3),(2,4))
    string printWithSuperNode(SuperNode *n2) {
        stringstream mStream;
        if (getLowestTaxon() < n2->getLowestTaxon()) {
            mStream << "(" << topology << "," << n2->topology << ")";
        }
        else {
            mStream << "(" << n2->topology << "," << topology << ")";
        }

        return mStream.str();
    }

    // if n1 represents (2,5) and n2 represents (1,3) and numTaxa is 6
    // output will be 1,2,3,5 | 4,6
    // required for printing tie information at the end of .concordance file.
    string printUngroupedWithSuperNode(SuperNode *n2, int numTaxa) {
        stringstream str;
        tset->printWithTaxonSet(str, n2->tset, numTaxa);
        return str.str();
    }

private:
    int number;
    int numLeafs;
    // each supernode can refer to two other supernodes at max except for root
    SuperNode *left;
    SuperNode *right;
    TaxonSet *tset;
    string topology;//topology represented by this supernode. For example if node has 2,1 as childs, topology would be (1,2)

};

class TieInfo {
public:
    TieInfo(){
    }

    void print(ostream& f) {
        f << tiedQuartet << "\t" << support << ", " << cUnits << ",  ";
        for (set<string>::iterator itr1 = ties.begin(); itr1 != ties.end();) {
            f << *itr1;
            if (++itr1 != ties.end()) f << ", " ;
        }
        f << endl;
    }

    void addTieInfo(int node1, int node2, vector<SuperNode *>& superNodes, int numTaxa) {
        ties.insert(superNodes[node1-1]->printUngroupedWithSuperNode(superNodes[node2-1], numTaxa));
    }

    void setSupport(double s) {
        support = s;
    }

    void setCUnits(double cu) {
        cUnits = cu;
    }

    void setTiedQuartet(string q) {
        tiedQuartet = q;
    }

    double getSupport() {
        return support;
    }

    double getCUnits() {
        return cUnits;
    }

    string getTiedQuartet() {
        return tiedQuartet;
    }

private:
    set<string> ties;
    double support; //tied support value
    double cUnits; // support value in coalescent units
    string tiedQuartet;
};

class TreeBuilder {
public:
    void getTree(Table* newTable, int numTaxa, string& top, string& topWithWts);
    void printTies(ostream& f);
private:
    string getTreeFromQuartetCounts(vector<vector<double> >& counts, int numTaxa, map<string, TieInfo*>& ties);
    double computeConfidence(int m, int n, vector<int>& activeNodes, vector<vector<double> >& counts);
    double computeNewConfidence(int i, int j, int b, vector<int>& activeNodes, vector<vector<double> >& counts, vector<vector<double> >& confidence);
    double computeNewCardinality(int maxI, int maxJ, int node1, int numTaxa, vector<int>& activeNodes, vector<vector<double> >& size);
    vector<SuperNode*> superNodes;
    vector<TieInfo*> qlist; // list of quartets, their supports n tie info
};
}
#endif /* QUARTETS_H_ */
