#ifndef __TGM__h
#define __TGM__h

class TreeWeight {
public:
    TreeWeight() {}
    void setWeight(double x) { weight=x; }
    void addWeight(double x) { weight += x; }
    double getWeight() { return weight; }
    void setIndex(int x) { index=x; }
    int getIndex() { return index; }
private:
    double weight;
    int index;
};

bool cmpTreeWeights(TreeWeight x,TreeWeight y);

class GeneCounts {
public:
    GeneCounts(int index, double c) {
        gene = index;
        count = c;
    }

    int getGene() const {
        return gene;
    }

    double getCount() const {
        return count;
    }

    void addCount(double c) {
        count += c;
    }

private:
    int gene;
    double count;
};

class Table {
public:
  virtual ~Table() {}
  
  virtual void addGeneCount(string top, int gene, double count) = 0;

  virtual void addGeneCount(int topIndex, int gene, double count) = 0;
  
  virtual double getCounts(int top, int gene) = 0;
  
  virtual double getCounts(string top, int gene) = 0;
  
  virtual const vector<GeneCounts> getCounts(string top) = 0;
  
  virtual double getTotalCounts(int top) = 0;
  
  virtual void reorder() = 0;
  
  virtual int getNumTrees() = 0;
  
  virtual vector<string> getTopologies() = 0;
  
  virtual void print(ostream &f) = 0;
  
  virtual vector<double> getSerialTable() = 0;
  
  virtual void readSerialTable(vector<double>& nt) = 0;
  
};

class TGMTable : public Table
{
public:
    TGMTable() {nGenes = 0;}

    ~TGMTable(){
		table.clear();
		tops.clear();
	}

    TGMTable(int numGenes, vector<string>& topologies) {
        int numTrees = topologies.size();
        tops = topologies;
        nGenes = numGenes;
        table.resize(numTrees);
        for (int i = 0; i < numTrees; i++) {
            table[i].resize(nGenes, 0);
        }
    }
    
    virtual vector<double> getSerialTable(){
	    vector<double> retTable;
	    int sz = this->table.size();
		for(int t=0; t < this->tops.size(); t++){
			for(int g=0; g < sz; g++){
				retTable.push_back(this->table[t][g]);
			}
		}
		return retTable;
	}

	virtual void readSerialTable(vector<double>& nt) {
		int sz = this->table.size();
		for(int t=0; t < this->tops.size(); t++){
			for(int g=0; g < sz; g++){
				double temp = nt[t*sz+g];
				if (temp > 0 && temp < 1)
				  temp = 0;
				this->table[t][g] += temp;
			}
		}
	}

    virtual void addGeneCount(string top, int gene, double count) {
        vector<string>::iterator it = find(tops.begin(), tops.end(), top);
        int index;
        if (it == tops.end()) {
            index = tops.size();
            tops.push_back(top);
            table.push_back(vector<double>(nGenes, 0));
        }
        else {
            index = it - tops.begin();
        }

        if (table[index].size() <= gene) {
            table[index].resize(gene + 1, 0.0);
        }
        table[index][gene] += count;
    }

    virtual const vector<GeneCounts> getCounts(string tp) {
        int index = find(tops.begin(), tops.end(), tp) - tops.begin();
        vector<GeneCounts> counts;
        for(int i = 0 ; i < table[index].size(); i++) {
            counts.push_back(GeneCounts(i, table[index][i]));
        }

        return counts;
    }

    virtual void reorder() {
        vector<TreeWeight> topWts(tops.size());
        for (int i = 0;  i < tops.size(); i++) {
            topWts[i].setIndex(i);
            topWts[i].setWeight(0);
        }

        int topIndex = 0;
        for (vector<vector<double> >::iterator itr = table.begin(); itr != table.end(); itr++, topIndex++) {
            for (vector<double>::iterator git = itr->begin(); git != itr->end(); git++) {
                topWts[topIndex].addWeight(*git);
            }
        }

        sort(topWts.begin(), topWts.end(), cmpTreeWeights);
        vector<string> topCopy = tops;
        int j = 0;
        for (vector<TreeWeight>::iterator it = topWts.begin(); it != topWts.end(); it++, j++) {
            tops[j] = topCopy[it->getIndex()];
        }

        vector<vector<double> > copy;
        copy.resize(table.size());
        for (int i = 0; i < table.size();i++) {
            copy[i].resize(table[i].size());
            assert(table[0].size() == table[i].size());
            for (int j = 0; j < table[i].size(); j++) {
                copy[i][j] = table[i][j];
            }
        }

        for(int i=0;i<table.size();i++) {
            for (int j = 0; j < table[i].size(); j++) {
                table[i][j] = copy[topWts[i].getIndex()][j];
            }
        }

        for (int i = 0; i < table.size(); i++) {
            for (int j = 0; j < table[i].size(); j++) {
                assert(table[i][j] == copy[topWts[i].getIndex()][j]);
            }
        }
    }

    virtual vector<string> getTopologies() {
        return tops;
    }

    virtual void addGeneCount(int topIndex, int gene, double count) {
        table[topIndex][gene] += count;
    }

    virtual void getCounts(int topIndex, unordered_map<int, double>& geneCounts) {
        for (int i = 0; i < table[topIndex].size(); i++) {
            geneCounts[i] = table[topIndex][i];
        }
    }

    virtual double getCounts(string top, int gene) {
        int n = find(tops.begin(),tops.end(),top) - tops.begin();
        if(n==tops.size()) {
            cerr << "Internal error in finding topology" << endl;
            exit(0);
        }
        return table[n][gene];
    }

    virtual double getTotalCounts(int top) {
        double total = 0;
        for (int i = 0; i < table[top].size(); i++) {
            total += table[top][i];
        }

        return total;
    }

    virtual double getCounts(int topIndex, int gene) {
        return table[topIndex][gene];
    }

    virtual int getNumTrees() {
        return table.size();
    }
    
    virtual void print(ostream &f)
    {
        for (int i = 0; i < tops.size(); i++) {
            f << "Topology: " << tops[i] << endl;
            for (int j = 0; j < table.size(); j++) {
                f << "Gene " << j << ": Count "<< table[i][j] << endl;
            }
            f << endl;
        }
    }
    


private:
//Why is this a double? Would take less mem for ints
    vector<vector<double> > table;
    vector<string> tops;
    int nGenes;
};


class TGM : public Table
{
public:
    TGM() {}

    ~TGM() {}

    TGM(vector<string>& topologies) {
        tops = topologies;
        int topIndex = 1;
        for (vector<string>::iterator itr = tops.begin(); itr != tops.end(); itr++, topIndex++) {
            topIndexMap[*itr] = topIndex;
        }
        topGeneCounts.resize(tops.size());
    }

    void addGeneCount(string top, int gene, double count) {
        int index = topIndexMap[top];
        if (index == 0) {
            index = tops.size() + 1;
            topIndexMap[top] = index;
            tops.push_back(top);
            vector<GeneCounts> geneCounts;
            geneCounts.push_back(GeneCounts(gene, count));
            topGeneCounts.push_back(geneCounts);
        }
        else {
            vector<GeneCounts>& geneCounts = topGeneCounts[index - 1];
            for (vector<GeneCounts>::iterator itr = geneCounts.begin(); itr != geneCounts.end(); itr++) {
                if (itr->getGene() == gene) {
                    itr->addCount(count);
                    return;
                }
            }

            geneCounts.push_back(GeneCounts(gene, count));
        }
    }

    virtual void addGeneCount(int topIndex, int gene, double count) {
        vector<GeneCounts>& geneCounts = topGeneCounts[topIndex];
        for (vector<GeneCounts>::iterator itr = geneCounts.begin(); itr != geneCounts.end(); itr++) {
            if (itr->getGene() == gene) {
                itr->addCount(count);
                return;
            }
        }

        geneCounts.push_back(GeneCounts(gene, count));
    }

    virtual const vector<GeneCounts> getCounts(string top) {
        return topGeneCounts[topIndexMap[top] - 1];
    }

    const vector<GeneCounts>& getCounts(int topIndex) {
        return topGeneCounts[topIndex];
    }

    virtual double getCounts(string top, int gene) {
        return getCounts(topIndexMap[top] - 1, gene);
    }

    virtual double getCounts(int top, int gene) {
        vector<GeneCounts> geneCounts = topGeneCounts[top];
        for (vector<GeneCounts>::iterator itr = geneCounts.begin(); itr != geneCounts.end(); itr++) {
            if (itr->getGene() == gene) {
                return itr->getCount();
            }
        }

        return 0.0;
    }

    virtual double getTotalCounts(int top) {
        double total = 0;
        for (int i = 0; i < topGeneCounts[top].size(); i++) {
            total += topGeneCounts[top][i].getCount();
        }

        return total;
    }

    virtual int getNumTrees() {
        return tops.size();
    }

    virtual void printAvgTopsGenes() {
        int total = 0;
        for (vector<vector<GeneCounts> >::iterator titr = topGeneCounts.begin(); titr != topGeneCounts.end(); titr++) {
            total += titr->size();
        }

        cerr << "Total: " << total << " " << topGeneCounts.size() << endl;
        cerr << "Average genes for topology: " << ((double) total/topGeneCounts.size()) << endl;
    }

    virtual void reorder()
    {
        vector<TreeWeight> topWts(tops.size());
        for (int i = 0;  i < tops.size(); i++) {
            topWts[i].setIndex(i);
            topWts[i].setWeight(0);
        }

        int topIndex = 0;
        for (vector<vector<GeneCounts> >::iterator itr = topGeneCounts.begin(); itr != topGeneCounts.end(); itr++, topIndex++) {
            for (vector<GeneCounts>::iterator git = itr->begin(); git != itr->end(); git++) {
                topWts[topIndex].addWeight(git->getCount());
            }
        }

        sort(topWts.begin(), topWts.end(), cmpTreeWeights);
        vector<string> topCopy = tops;
        vector<vector<GeneCounts> > topGeneCountsCopy = topGeneCounts;
        int j = 0;
        for (vector<TreeWeight>::iterator it = topWts.begin(); it != topWts.end(); it++, j++) {
            tops[j] = topCopy[it->getIndex()];
            topGeneCounts[j] = topGeneCountsCopy[it->getIndex()];
            topIndexMap[tops[j]] = j + 1;
        }
    }

    virtual vector<string> getTopologies()
    {
        return tops;
    }

    virtual void print(ostream &f)
    {
        for (vector<string>::iterator titr = tops.begin();titr != tops.end(); titr++) {
            f << "Topology : " << *titr << endl;
            vector<GeneCounts>& geneCounts = topGeneCounts[topIndexMap[*titr] - 1];
            for (vector<GeneCounts>::iterator git = geneCounts.begin(); git != geneCounts.end(); git++) {
                f << "Gene " << git->getGene() << ": count " << git->getCount() << endl;
            }

            f << endl;
        }
    }
    
    virtual vector<double> getSerialTable(){}
	virtual void readSerialTable(vector<double>& nt){}

private:
    vector<string> tops;
    vector<vector<GeneCounts> > topGeneCounts;
    unordered_map<string, int> topIndexMap;
};


#endif
