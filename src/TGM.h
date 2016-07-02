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


class TGMTable 
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
            table[i].resize(nGenes);
            for (int k=0; k <nGenes; k++)
              table[i][k] = 0;
        }
    }
    
    vector<int> getSerialTable(){
	    vector<int> retTable;
	    int sz = this->table.size();
		for(int t=0; t < sz; t++){
			for(int g=0; g < this->table[t].size(); g++){
				retTable.push_back(this->table[t][g]);
			}
		}
		return retTable;
	}

	void readSerialTable(vector<int>& nt) {
		int sz = this->table.size();
		int ts = this->tops.size();
		cout << "Num tops: "<<ts<<" Num genes = "<<this->nGenes<<endl;
		int t, g;
		for (int i=0; i< nt.size(); i++){
			this->table[t][g] += nt[i];
			g++;
			if (g >= this->nGenes){
				t++;
				g=0;
			}
		}
		cout << "Last values: t="<<t<<" & g="<<g<<endl;
	}


    const vector<GeneCounts> getCounts(string tp) {
        int index = find(tops.begin(), tops.end(), tp) - tops.begin();
        vector<GeneCounts> counts;
        for(int i = 0 ; i < table[index].size(); i++) {
            counts.push_back(GeneCounts(i, table[index][i]));
        }

        return counts;
    }

    void reorder() {
        vector<TreeWeight> topWts(tops.size());
        for (int i = 0;  i < tops.size(); i++) {
            topWts[i].setIndex(i);
            topWts[i].setWeight(0);
        }

        int topIndex = 0;
        for (vector<vector<int> >::iterator itr = table.begin(); itr != table.end(); itr++, topIndex++) {
            for (vector<int>::iterator git = itr->begin(); git != itr->end(); git++) {
                topWts[topIndex].addWeight(*git);
            }
        }

        sort(topWts.begin(), topWts.end(), cmpTreeWeights);
        vector<string> topCopy = tops;
        int j = 0;
        for (vector<TreeWeight>::iterator it = topWts.begin(); it != topWts.end(); it++, j++) {
            tops[j] = topCopy[it->getIndex()];
        }

        vector<vector<int> > copy;
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

    vector<string> getTopologies() {
        return tops;
    }

    
    void addGeneCount(int topIndex, int gene) {
        table[topIndex][gene]++;
    }
    
    void addGeneCount(string top, int gene) {
        vector<string>::iterator it = find(tops.begin(), tops.end(), top);
        int index;
        if (it == tops.end()) {
            index = tops.size();
            tops.push_back(top);
            table.push_back(vector<int>(nGenes, 0));
        }
        else {
            index = it - tops.begin();
        }

        if (table[index].size() <= gene) {
            table[index].resize(gene + 1, 0.0);
        }
        table[index][gene]++;
    }


    void getCounts(int topIndex, unordered_map<int, int>& geneCounts) {
        for (int i = 0; i < table[topIndex].size(); i++) {
            geneCounts[i] = table[topIndex][i];
        }
    }

    int getCounts(string top, int gene) {
        int n = find(tops.begin(),tops.end(),top) - tops.begin();
        if(n==tops.size()) {
            cerr << "Internal error in finding topology" << endl;
            exit(0);
        }
        return table[n][gene];
    }

   int getTotalCounts(int top) {
        int total = 0;
        for (int i = 0; i < table[top].size(); i++) {
            total += table[top][i];
        }

        return total;
    }

    int getCounts(int topIndex, int gene) {
        return table[topIndex][gene];
    }

    int getNumTrees() {
        return table.size();
    }
    
    void print(ostream &f)
    {
        for (int i = 0; i < tops.size(); i++) {
            f << "Topology: " << tops[i] << endl;
            for (int j = 0; j < nGenes; j++) {
                f << "Gene " << j << ": Count "<< table[i][j] << endl;
            }
            f << endl;
        }
    }
    


private:
//Convert to ints cuts memory usage of table in half
    vector<vector<int> > table;
    vector<string> tops;
    int nGenes;
};


class TGM 
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

    void addGeneCount(int topIndex, int gene, double count) {
        vector<GeneCounts>& geneCounts = topGeneCounts[topIndex];
        for (vector<GeneCounts>::iterator itr = geneCounts.begin(); itr != geneCounts.end(); itr++) {
            if (itr->getGene() == gene) {
                itr->addCount(count);
                return;
            }
        }

        geneCounts.push_back(GeneCounts(gene, count));
    }

    const vector<GeneCounts> getCounts(string top) {
        return topGeneCounts[topIndexMap[top] - 1];
    }

    const vector<GeneCounts>& getCounts(int topIndex) {
        return topGeneCounts[topIndex];
    }

    double getCounts(string top, int gene) {
        return getCounts(topIndexMap[top] - 1, gene);
    }

    double getCounts(int top, int gene) {
        vector<GeneCounts> geneCounts = topGeneCounts[top];
        for (vector<GeneCounts>::iterator itr = geneCounts.begin(); itr != geneCounts.end(); itr++) {
            if (itr->getGene() == gene) {
                return itr->getCount();
            }
        }

        return 0.0;
    }

    double getTotalCounts(int top) {
        double total = 0;
        for (int i = 0; i < topGeneCounts[top].size(); i++) {
            total += topGeneCounts[top][i].getCount();
        }

        return total;
    }

    int getNumTrees() {
        return tops.size();
    }

    void printAvgTopsGenes() {
        double total = 0;
        for (vector<vector<GeneCounts> >::iterator titr = topGeneCounts.begin(); titr != topGeneCounts.end(); titr++) {
            total += titr->size();
        }

        cerr << "Total: " << total << " " << topGeneCounts.size() << endl;
        cerr << "Average genes for topology: " << ((double) total/topGeneCounts.size()) << endl;
    }

    void reorder()
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

    vector<string> getTopologies()
    {
        return tops;
    }

    void print(ostream &f)
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
    

private:
    vector<string> tops;
    vector<vector<GeneCounts> > topGeneCounts;
    unordered_map<string, int> topIndexMap;
};


#endif
