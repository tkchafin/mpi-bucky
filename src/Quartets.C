#define NDEBUG
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
#include<map>
using namespace std;
#include "boost/unordered_map.hpp"
using namespace boost;
#include "TGM.h"
#include "Quartets.h"
using namespace quartet;

string DUMMY = "NONEED";
vector<vector<int> > nCr;
void precomputeNcR(int numTaxa) {
    // compute nC1, nC2, nC3, nC4 for n in [0, numTaxa - 1], this is used for ranking quartets
    int numSuperNodes = 2 * numTaxa - 3;
    nCr.resize(numSuperNodes, vector<int>(4, 0));
    for (int i = 0; i < numSuperNodes; i++) {
        nCr[i][0] = i;
    }

    for (int i = 1; i < numSuperNodes; i++) {
        for (int j = 1; j < 4; j++) {
            nCr[i][j] = nCr[i-1][j-1] + nCr[i-1][j];
        }
    }
}

int computeRIndex(int t1, int t2, int t3, int t4) {
    //compute subset rank. assume sorted order is t1, t2, t3, t4 -
    // rank for a quartet can be obtained by computing t1 Choose 1 + t2 Choose 2 + t3 Choose 3 + t4 Choose 4
    return nCr[t1 - 1][0] + nCr[t2 - 1][1] + nCr[t3 - 1][2] + nCr[t4 - 1][3];
}

void getQuartetRowColumnIndex(int t1, int t2, int t3, int t4, int& rIndex, int& cIndex) {
    // make t3 smallest of t3, t4
    if (t3 > t4) {
        int temp = t4;
        t4 = t3;
        t3 = temp;
    }

    // make t1 smallest of t1, t2
    if (t1 > t2) {
        int temp = t1;
        t1 = t2;
        t2 = temp;
    }

    if (t1 > t3) {
        // swap t1 with t3, t2 with t4 - makes t1 smallest
        int temp = t1;
        t1 = t3;
        t3 =temp
                ;
        temp = t2;
        t2 = t4;
        t4 = temp;
    }

    if (t2 < t3) {
        //t2 2nd smallest - order t1, t2, t3, t4
        assert (t1 < t2);
        assert (t2< t3);
        assert (t3< t4);
        cIndex = 0;
        rIndex = computeRIndex(t1, t2, t3, t4);
    }
    else {
        //t3 2nd smallest
        if (t2 < t4) {
            //  order t1, t3, t2, t4
            assert (t1 < t3);
            assert (t3 < t2);
            assert (t2< t4);
            cIndex = 1;
            rIndex = computeRIndex(t1, t3, t2, t4);
        }
        else {
            // order t1, t3, t4, t2
            assert (t1 < t3);
            assert (t3 < t4);
            assert (t4 < t2);
            cIndex = 2;
            rIndex = computeRIndex(t1, t3, t4, t2);
        }
    }
}

void printQuartet(string& quartet, TaxonSet*& t1, TaxonSet*& t2, TaxonSet*& t3, TaxonSet*& t4) {
    stringstream str1, str2;
    int min1 = t1->getTSet()[0];
    int min2 = t2->getTSet()[0];
    int min3 = t3->getTSet()[0];
    int min4 = t4->getTSet()[0];

    int min5, min6;
    if (min1 < min2) {
        str1 << *t1;
        str1 << "; ";
        str1 << *t2;
        min5 = min1;
    }
    else {
        str1 << *t2;
        str1 << "; ";
        str1 << *t1;
        min5 = min2;
    }

    if (min3 < min4) {
        str2 << *t3;
        str2 << "; ";
        str2 << *t4;
        min6 = min3;
    }
    else {
        str2 << *t4;
        str2 << "; ";
        str2 << *t3;
        min6 = min4;
    }

    if (min5 < min6) {
        str1 << "|";
        str2 << "}";
        quartet = "{" + str1.str() + str2.str();
    }
    else {
        str1 << "}";
        str2 << "|";
        quartet = "{" + str2.str() + str1.str();
    }
}

void printSet(ostream &str, vector<int>& result) {
    for (vector<int>::iterator it = result.begin(); it != result.end(); it++) {
         str << *it;
         if (it + 1 != result.end()) str << ",";
     }
}

void printSetComplement(ostream &str, vector<int>& result, int numTaxa) {
    int lowestTaxon = 1;
    vector<int>::iterator it = result.begin();
    for (; it != result.end() && lowestTaxon <= numTaxa;) {
        if (*it != lowestTaxon) {
            str << lowestTaxon << ",";
        }
        else {
            it++;
        }

        lowestTaxon++;
    }

    if (lowestTaxon > numTaxa) {
        for (;it != result.end(); it++) {
            str << *it << ",";
        }
    }
    else if (it == result.end()) {
        for (;lowestTaxon <= numTaxa; lowestTaxon++) {
            str << lowestTaxon << ",";
        }
    }
    long pos = str.tellp();
    str.seekp(pos - 1); //remove last ','
}

TaxonSet* getSetComplement(TaxonSet *t, int numTaxa) {
    vector<int>& result = t->getTSet();
    TaxonSet *complement = new TaxonSet();
    int lowestTaxon = 1;
    vector<int>::iterator it = result.begin();
    for (; it != result.end() && lowestTaxon <= numTaxa;) {
        if (*it != lowestTaxon) {
            complement->add(lowestTaxon);
        }
        else {
            it++;
        }

        lowestTaxon++;
    }

    if (lowestTaxon > numTaxa) {
        for (;it != result.end(); it++) {
            complement->add(*it);
        }
    }
    else if (it == result.end()) {
        for (;lowestTaxon <= numTaxa; lowestTaxon++) {
            complement->add(lowestTaxon);
        }
    }

    return complement;
}

void TreeBuilder::getTree(TGMTable* newTable, int numTaxa, string& top, string& topWithWts) {
    vector<string> topologies = newTable->getTopologies();
    precomputeNcR(numTaxa);
    int numNodes = numTaxa + numTaxa - 3;
    int numQuartets = numNodes * (numNodes - 1) * (numNodes - 2) * (numNodes -3) / 24;
    vector<vector<double> > counts;
    counts.resize(numQuartets);//(numTaxa -3) new super nodes are created
    for (int i = 0; i < counts.size(); i++) {
        counts[i].resize(3, 0.0);
    }

    for (int i = 0; i < topologies.size(); i++) {
        double totalCount = (double) newTable->getTotalCounts(i);
        Tree t(numTaxa, topologies[i]);
        vector<int> rind, cind;
        t.getQuartets(rind, cind);
        for (int i = 0; i < rind.size(); i++) {
            counts[rind[i]][cind[i]] += totalCount;
        }
    }

    for (int i = 0; i < counts.size(); i++) {
        double total = counts[i][0] + counts[i][1] + counts[i][2];
        if (total != 0.0) {
            counts[i][0] /= total;
            counts[i][1] /= total;
            counts[i][2] /= total;
        }
    }
    map<string, TieInfo*> ties;
    top = getTreeFromQuartetCounts(counts, numTaxa, ties);
    Tree t(numTaxa, top);
    for (int i = 0; i < numTaxa; i++) {
        t.getEdge(i)->addWeight(10.0);//set prob 1.0, branch length in #coalescent units as 10 for all leaf edges
    }

    //compute weights for internal edges using formula W(AB|CD)/|A||B||C||D|
    for (int i = numTaxa; i < t.getNumEdges(); i++) {
        Edge *e = t.getEdge(i);
        Node *n1 = e->getNode(0);
        Node *n2 = e->getNode(1);
        vector<int> rInd, cInd;
        string quartet;
        int abcd = t.getQuartets(rInd, cInd, n1, n2, quartet); // returns |A||B||C||D|
        double wt = 0.0;
        for (int i = 0; i < rInd.size(); i++) {
            wt += counts[rInd[i]][cInd[i]];
        }

        wt /= (double) abcd;
        double s = wt;
        //convert prob to # of coalescent units.
        if (wt > 0.99997)
            wt = 10.0;
        else if (wt <= .33334)
            wt = 0.0;
        else
            wt = -log(3.0/2.0 * (1-wt));

        e->addWeight(wt);

        // update quartet weight, #coalescent units
        string top;
        TieInfo* info;
        if (n1->getNeighbor(0) == n2) {
            stringstream str;
            n1->print(str);
            top = str.str();
            map<string, TieInfo*>::iterator it = ties.find(top);
            if (it != ties.end()) {
                info = it->second;
            }
            else {
                info = new TieInfo();
            }
        }
        else {
            stringstream str;
            n2->print(str);
            top = str.str();
            map<string, TieInfo*>::iterator it = ties.find(top);
            if (it != ties.end()) {
                info = it->second;
            }
            else {
                info = new TieInfo();
            }
        }

        info->setSupport(s);
        info->setCUnits(wt);
        info->setTiedQuartet(quartet);
        qlist.push_back(info);
    }

    t.modifyOutgroup(numTaxa, top, topWithWts);//choose last taxon as outgroup and print topologies to top, topWithWts.
}

// returns smallest number in range [1, numTaxa] that is not present in t1, t2
int getComplement(vector<int>& t1, vector<int>& t2, int numTaxa)
{
    vector<int>::iterator it1 = t1.begin();
    vector<int>::iterator it2 = t2.begin();
    for (int i = 1; i <= numTaxa; i++) {
        if (i == *it1) {
            it1++;
            continue;
        }

        if (i == *it2) {
            it2++;
            continue;
        }

        return i;
    }
    cerr << "Warning: function 'getComplement' called on an inappropriate set of 2 clades.\n";
    cout << "Warning: function 'getComplement' called on an inappropriate set of 2 clades.\n";
    return numTaxa+1; // this case should never occur
}

string TreeBuilder::getTreeFromQuartetCounts(vector<vector<double> >& counts, int numTaxa, map<string, TieInfo*>& ties) {
    vector<vector<double> > normCounts(counts.size());// normalized counts
    for (int i = 0; i < counts.size(); i++) {
        normCounts[i] = counts[i];
    }

    // normalize counts.
    int numOfCurrentQuartets = nCr[numTaxa - 1][2] + nCr[numTaxa - 1][3]; //C(numTaxa,4)  = C(numTaxa - 1, 3) + C(numTaxa - 1, 4)
    for (int i = 0; i < numOfCurrentQuartets; i++) {
        if (normCounts[i][0] == normCounts[i][1]) {
            if (normCounts[i][1] == normCounts[i][2]) {
                normCounts[i][0] = .333334;
                normCounts[i][1] = .333333;
                normCounts[i][2] = .333333;
            }
            else if (normCounts[i][2] < normCounts[i][1]){
                normCounts[i][0] = 0.5;
                normCounts[i][1] = 0.5;
                normCounts[i][2] = 0.0;
            }
            else {
                normCounts[i][0] = 0.0;
                normCounts[i][1] = 0.0;
                normCounts[i][2] = 1.0;
            }
        }
        else if (normCounts[i][0] < normCounts[i][1]) {
            if (normCounts[i][1] == normCounts[i][2]) {
                normCounts[i][0] = 0.0;
                normCounts[i][1] = 0.5;
                normCounts[i][2] = 0.5;
            }
            else if (normCounts[i][1] < normCounts[i][2]) {
                normCounts[i][0] = 0.0;
                normCounts[i][1] = 0.0;
                normCounts[i][2] = 1.0;
            }
            else {
                normCounts[i][0] = 0.0;
                normCounts[i][1] = 1.0;
                normCounts[i][2] = 0.0;
            }
        }
        else {
            if (normCounts[i][0] == normCounts[i][2]) {
                normCounts[i][0] = 0.5;
                normCounts[i][1] = 0.0;
                normCounts[i][2] = 0.5;
            }
            else if (normCounts[i][0] > normCounts[i][2]) {
                normCounts[i][0] = 1.0;
                normCounts[i][1] = 0.0;
                normCounts[i][2] = 0.0;
            }
            else {
                normCounts[i][0] = 0.0;
                normCounts[i][1] = 0.0;
                normCounts[i][2] = 1.0;
            }
        }
    }

    int numNodes = numTaxa + numTaxa - 3;
    for (int i = 0; i < numNodes; i++) {
        superNodes.push_back(new SuperNode(i + 1, i < numTaxa));
    }

    vector<int> activeNodes; // keeps track of list of super nodes active in current iteration
    for (int i = 0; i < numTaxa; i++) {
        activeNodes.push_back(i + 1);
    }

    vector<vector<double> > confidence(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > size(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > support(numNodes, vector<double>(numNodes, 0));
    for (int j = 1; j <= numTaxa; j++) {
        for (int i = 1; i < j; i++) {
            confidence[i - 1][j -1] = computeConfidence(i, j, activeNodes, normCounts);
            size[i - 1][j - 1] = nCr[numTaxa - 2][1]; //initial cardinality is always numTaxa-2 choose 2, no need to do summation
            support[i -1][j - 1] = confidence[i - 1][j -1]/size[i - 1][j - 1];
        }
    }

    int currentNode = numTaxa;
    while (true) {
        // find nodes that have maximum support, call them maxI, maxJ
        int maxI = -1, maxJ = -1;
        double maxSupport = -1;
        vector<int> maxIs, maxJs; // In case of tie, these keep track of nodes that have maximum support
        for (int j = 0; j < activeNodes.size(); j++) {
            int node1 = activeNodes[j];
            for (int i = j + 1; i < activeNodes.size(); i++) {
                int node2 = activeNodes[i];
                if (support[node1 - 1][node2 - 1] == maxSupport) {
                    if (maxI == node1 || maxI == node2 || maxJ == node1 || maxJ == node2) {
                        // add only if there is a conflict between maxI, maxJ and node1, node2
                        maxIs.push_back(node1);
                        maxJs.push_back(node2);
                    }
                }
                else if (support[node1 - 1][node2 - 1] > maxSupport) {
                    maxSupport = support[node1 - 1][node2 - 1];
                    maxI = node1;
                    maxJ = node2;
                    maxIs.clear();
                    maxJs.clear();
                }
            }
        }

        // if there are multiple maxima, break tie by computing supports
        // with original counts
        if (maxIs.size() > 0) {
            double conf = computeConfidence(maxI, maxJ, activeNodes, counts);
            double sz = size[maxI - 1][maxJ - 1];
            maxSupport = conf / sz;
            string mStr = superNodes[maxI - 1]->printWithSuperNode(superNodes[maxJ - 1]);

            for (int j = 0; j < maxIs.size(); j++) {
                int node1 = maxIs[j];
                int node2 = maxJs[j];
                double conf = computeConfidence(node1, node2, activeNodes, counts);
                double sz = size[node1 - 1][node2 - 1];
                double support = conf / sz;

                if (support > maxSupport) {
                    //maxI, maxJ is no longer the chosen one. Delete all its tieinfo.
                    map<string, TieInfo*>::iterator it = ties.find(mStr);
                    if (it != ties.end()) {
                        TieInfo*& info = it->second;
                        ties.erase(mStr);
                        delete info;
                    }
                    maxI = node1;
                    maxJ = node2;
                    maxSupport = support;
                    mStr = superNodes[maxI - 1]->printWithSuperNode(superNodes[maxJ - 1]);
                }
                else if (support == maxSupport) {
                    if (maxI == node1 || maxI == node2 || maxJ == node1 || maxJ == node2) {
                        map<string, TieInfo*>::iterator it = ties.find(mStr);
                        TieInfo* info;
                        if (it == ties.end()) {
                            info = new TieInfo();
                            ties[mStr] = info;
                        }
                        else {
                            info = it->second;
                        }

                        info->addTieInfo(node1, node2, superNodes, numTaxa);
                    }
                }
            }
        }

        currentNode++;
        superNodes[currentNode - 1]->add(superNodes[maxI - 1], superNodes[maxJ - 1]);
        vector<int>::iterator itr = find(activeNodes.begin(), activeNodes.end(), maxI);
        assert(itr != activeNodes.end());
        activeNodes.erase(itr);
        itr = find(activeNodes.begin(), activeNodes.end(), maxJ);
        assert(itr != activeNodes.end());
        activeNodes.erase(itr);

        if (activeNodes.size() == 2) {
            activeNodes.push_back(currentNode); //quit when there are 3 active nodes(super nodes)
            break;
        }

        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            for (int j = i + 1; j < activeNodes.size(); j++) {
                int node2 = activeNodes[j];
                for (int k = j + 1; k < activeNodes.size(); k++) {
                    int node3 = activeNodes[k];
                    //find weights of 3 different resolutions for node1, node2, node3 and newly created node(currentNode)
                    int rInd1, cInd1, rInd2, cInd2;
                    int rIndResult, cIndResult;

                    getQuartetRowColumnIndex(node1, node2, node3, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node1, node2, node3, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node1, node2, node3, currentNode, rIndResult, cIndResult);
                    normCounts[rIndResult][cIndResult] = normCounts[rInd1][cInd1] + normCounts[rInd2][cInd2];
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];

                    getQuartetRowColumnIndex(node1, node3, node2, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node1, node3, node2, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node1, node3, node2, currentNode, rIndResult, cIndResult);
                    normCounts[rIndResult][cIndResult] = normCounts[rInd1][cInd1] + normCounts[rInd2][cInd2];
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];

                    getQuartetRowColumnIndex(node3, node2, node1, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node3, node2, node1, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node3, node2, node1, currentNode, rIndResult, cIndResult);
                    normCounts[rIndResult][cIndResult] = normCounts[rInd1][cInd1] + normCounts[rInd2][cInd2];
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];
                }
            }

        }

        //compute confidence, support for currentNode and every other active node
        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            confidence[node1 - 1][currentNode - 1] = computeNewConfidence(maxI, maxJ, node1,activeNodes,normCounts,confidence);
            size[node1 - 1][currentNode - 1] = computeNewCardinality(maxI, maxJ, node1, numTaxa, activeNodes, size);
            support[node1 - 1][currentNode - 1] = confidence[node1 - 1][currentNode - 1] / size[node1 - 1][currentNode - 1];
        }

        // modify confidence, support of active node pairs from previous iteration
        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            for (int j = i + 1; j < activeNodes.size(); j++) {
                int node2 = activeNodes[j];
                int rind, cind;
                getQuartetRowColumnIndex(maxI, maxJ, node1, node2, rind, cind);
                confidence[node1 - 1][node2 - 1] -=  normCounts[rind][cind];
                size[node1 - 1][node2 - 1] -= superNodes[maxI - 1]->getNumNodes()
                        * superNodes[maxJ - 1]->getNumNodes()
                        * superNodes[node1 - 1]->getNumNodes()
                        * superNodes[node2-1]->getNumNodes();

                support[node1 -1][node2 - 1] = confidence[node1 - 1][node2 - 1] / size [node1 - 1][node2 - 1];
            }
        }
        activeNodes.push_back(currentNode);
    }

    SuperNode *node1 = superNodes[activeNodes[0] - 1];
    SuperNode *node2 = superNodes[activeNodes[1] - 1];
    SuperNode *node3 = superNodes[activeNodes[2] - 1];;

    stringstream top;
    top << "(";
    node1->print(top);
    top << ",";
    node2->print(top);
    top << ",";
    node3->print(top);
    top << ");";
    return top.str();
}

//compute initial confidence, can also be used for computing confidence in later steps but inefficient,  use computeNewConfidence instead
double TreeBuilder::computeConfidence(int m, int n, vector<int>& activeNodes, vector<vector<double> >& counts) {
    double confidence = 0.0;
    for (int i = 0; i < activeNodes.size(); i++) {
        int node1 = activeNodes[i];
        if (node1 == m || node1 == n)
            continue;
        for (int j = i + 1; j < activeNodes.size(); j++) {
            int node2 = activeNodes[j];
            if (node2 == m || node2 == n)
                continue;

            int rIndex, cIndex;
            getQuartetRowColumnIndex(m, n, node1, node2, rIndex, cIndex);
            confidence += counts[rIndex][cIndex];
        }
    }

    return confidence;
}

// compute confidence of pairs that contain new super node, formula(D!=K,B): C(K, B) = C(I,B) + C(J,B) - Sum(counts[I][B][J][D] + counts[J][B][I][D])
double TreeBuilder::computeNewConfidence(int i, int j, int b, vector<int>& activeNodes, vector<vector<double> >& counts, vector<vector<double> >& confidence) {
    double conf = 0.0;
    if (i < b) {
        conf += confidence[i - 1][b - 1];
    }
    else {
        conf += confidence[b - 1][i - 1];
    }

    if (j < b) {
        conf += confidence[j - 1][b - 1];
    }
    else {
        conf += confidence[b - 1][j - 1];
    }

    for (int ii = 0; ii < activeNodes.size(); ii++) {
        int d = activeNodes[ii];
        if (d == b)
            continue;

        int rind, cind;
        getQuartetRowColumnIndex(i, b, j, d, rind, cind);
        conf -= counts[rind][cind];

        getQuartetRowColumnIndex(j, b, i, d, rind, cind);
        conf -= counts[rind][cind];
    }

    return conf;
}

// compute cardinality of pairs that contain new super node, formula T(K, B) = T(I,B) + T(J,B) - 2 * size(I) * size(J) * size(B) * (ntaxa - size(B) - size(K))
double TreeBuilder::computeNewCardinality(int maxI, int maxJ, int node1, int numTaxa, vector<int>& activeNodes, vector<vector<double> >& size) {
    double sz = 0.0;
    if (maxI < node1) {
        sz += size[maxI - 1][node1 - 1];
    }
    else {
        sz += size[node1 - 1][maxI - 1];
    }

    if (maxJ < node1) {
        sz += size[maxJ - 1][node1 - 1];
    }
    else {
        sz += size[node1 - 1][maxJ - 1];
    }

    int sizeI = superNodes[maxI - 1]->getNumNodes();
    int sizeJ = superNodes[maxJ - 1]->getNumNodes();
    int sizeB = superNodes[node1 - 1]->getNumNodes();
    int sizeK = sizeI + sizeJ;
    sz -= 2 * sizeI * sizeJ * sizeB * (numTaxa - sizeB - sizeK);
    return sz;
}

void TreeBuilder::printTies(ostream& f) {
    if (qlist.empty())
        return;

    vector<TieInfo*>::iterator tieitr = qlist.begin();
    f << "Four-way partitions in the Population Tree:";
    f << " sample-wide CF, coalescent units and Ties(if present)" << endl;
    for(;tieitr != qlist.end(); tieitr++) {
        TieInfo *info = *tieitr;
        info->print(f);
    }
}

Node* Node::getNeighbor(int i) const { return edges[i]->getOtherNode(this); }

void Node::print(ostream& f) const
{
    if (leaf) {
        f << t[0]->getTSet()[0];
        return;
    }

    f << "(";
    getNeighbor(1)->print(f);
    f << ",";
    getNeighbor(2)->print(f);
    f << ")";
}

void Node::print(ostream& f, ostream& topStr, const Node *caller, const Node *root, int numTaxa) const
{
    if (leaf) {
        f << t[0]->getTSet()[0];
        topStr << t[0]->getTSet()[0];
        return;
    }

    int rank[2];
    Node *n[2];
    double weight[2];

    // print is called from one of the three neighbors, find other two neighbors
    for (int i = 2, pos = 1; pos >= 0; i--) {
        Node *temp = getNeighbor(i);
        if (temp != caller) {
            n[pos] = temp;
            weight[pos] =edges[i]->getWeight();
            if (i != 0 || root == this) {
                rank[pos--] = t[i]->getTSet()[0];
            }
            else {
                rank[pos--] = getComplement(t[1]->getTSet(), t[2]->getTSet(), numTaxa);
            }
        }
    }

    // check which neighbor has lower taxa and print it first
    f << "(";
    topStr << "(";
    if (rank[0] < rank[1]) {
        n[0]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[0] << ",";
        topStr << ",";
        n[1]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[1];
    }
    else {
        n[1]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[1] << ",";
        topStr << ",";
        n[0]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[0];

    }
    f << ")";
    topStr << ")";

}

void Edge::print(ostream& f,int numTaxa) const
{
    f << setw(3) << number << "(" << weight << ")" << ":";
    f << " n[" << nodes[0]->getNumber() << "] <-> " << "n[" << nodes[1]->getNumber() << "], split = ";
    f << endl;
}

// reroots the tree to parent node of node representing taxa
// Note: does not change the tree structure, just returns a 'new' topology string
void Tree::modifyOutgroup(int taxon, string& top, string& topWithWts) const
{
    stringstream f;// stream to print topology with weights
    f.precision(3);
    stringstream topStr;// stream to print topology without weights
    int index = 0;
    //identify index of leaf node representing this taxa.
    for (;index < numTaxa; index++) {
        if (nodes[index]->getTset(0)->getTSet()[0] == taxon) {
            break;
        }
    }

    Node *newRoot = nodes[index]->getNeighbor(0);
    Node *nbrs[3];
    int tset[3];
    for (int i = 0; i < 3; i++) {
        nbrs[i] = newRoot->getNeighbor(i);
        tset[i] = newRoot->getTset(i)->getTSet()[0];
    }

    if (newRoot != root) {
        // since new root is not actually a root, taxonsets on 3 edges will be {a}U{b}, {a}, {b}
        // change the first set to: compliment of {a}U{b}
        tset[0] = getComplement(newRoot->getTset(1)->getTSet(), newRoot->getTset(2)->getTSet(), numTaxa);
    }

    f << "(";
    topStr << "(";
    if (tset[0] < tset[1]) {
        if (tset[0] < tset[2]) {
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            if (tset[1] < tset[2]) {
                nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
                topStr << ",";
                nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(2)->getWeight();
            }
            else {
                nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
                topStr << ",";
                nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(1)->getWeight();
            }
        }
        else {
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
            topStr << ",";
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(1)->getWeight();
        }
    }
    else if (tset[1] < tset[2]) {
        nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
        topStr << ",";
        if (tset[0] < tset[2]) {
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight();
        }
        else {
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
            topStr << ",";
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight();
        }
    }
    else {
        nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
        topStr << ",";
        nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
        topStr << ",";
        nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(0)->getWeight();
    }

    f << ");";
    topStr << ");";

    topWithWts = f.str();
    top = topStr.str();
}

void Tree::construct(string& top) {
    int currentLeafNode = 0, currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
    Node *n1,*n2,*n3;
    char ch;
    istringstream topStr(top);
    topStr >> ch;
    n1 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    n2 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    n3 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
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

    root->mergeTaxa(0, n1, 0);
    root->mergeTaxa(1, n2, 0);
    root->mergeTaxa(2, n3, 0);
}

Node* Tree::connectInt(istream& topStr,int& currentLeafNode,int& currentInternalNode,int& currentLeafEdge,int& currentInternalEdge)
{
    int n;
    char ch = topStr.peek();
    if(ch!='(') {
        topStr >> n;
        nodes[currentLeafNode]->addTaxa(0, n);
        return nodes[currentLeafNode++];
    }
    else {
        topStr >> ch;
        Node *node1,*node2;
        node1 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
        topStr >> ch;
        if(ch!=',') {
            cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl;
            exit(1);
        }
        node2 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
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

    node3->mergeTaxa(1, node1, 0);
    node3->mergeTaxa(2, node2, 0);
    node3->mergeTaxa(0, node1, 0);
    node3->mergeTaxa(0, node2, 0);
    return node3;
}

// returns tsets of size 4. Taxonsets on path from n1 to n2 are skipped
void Tree::getTsets(vector<TaxonSet*>& tsets, Node *n1, Node *n2) {
    int i, j;
    int matchI = 0, matchJ = 0;
    int greatestMatch = 0;
    for (i = 0; i < 3; i++) {
        TaxonSet*& t1 = n1->getTset(i);
        for (j = 0; j < 3; j++) {
            TaxonSet*& t2 = n2->getTset(j);
            int match = t1->intersect(t2);
            if (match != 0 && match >= greatestMatch) {
                greatestMatch = match;
                matchI = i;
                matchJ= j;
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (i != matchI) {
            TaxonSet *tset = n1->getTset(i);
            if (n1 != root && i == 0) {
                tsets.push_back(getSetComplement(tset, numTaxa));
            }
            else {
                tsets.push_back(tset);
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (i != matchJ) {
            TaxonSet *tset = n2->getTset(i);
            if (n2 != root && i == 0) {
                tsets.push_back(getSetComplement(tset, numTaxa));
            }
            else {
                tsets.push_back(tset);
            }
        }
    }
}

void Tree::getQuartets(vector<int>& rInd, vector<int>& cInd) {
    for (int i = numTaxa; i < numNodes; i++) {
        Node *n1 = nodes[i];
        for (int j = i + 1; j < numNodes; j++) {
            Node *n2 = nodes[j];
            getQuartets(rInd, cInd, n1 , n2, DUMMY);
        }
    }
}

int Tree::getQuartets(vector<int>& rInd, vector<int>& cInd, Node *n1, Node *n2, string& quartet) {
//    for (int i = numTaxa; i < numNodes; i++) {

//        for (int j = i + 1; j < numNodes; j++) {

            vector<TaxonSet*> txnsets;
            getTsets(txnsets, n1, n2);
            vector<vector<int> > tsets(4);
            for (int i = 0; i < txnsets.size(); i++) {
                tsets[i] = txnsets[i]->getTSet();
            }

            assert(tsets.size() == 4);
            for (int j1 = 0; j1 < tsets[0].size(); j1++) {
                for (int j2 = 0; j2 < tsets[1].size(); j2++) {
                    for (int j3 = 0; j3 < tsets[2].size(); j3++) {
                        for (int j4 = 0; j4 < tsets[3].size(); j4++) {
                            int rIndex, cIndex;
                            getQuartetRowColumnIndex(tsets[0][j1], tsets[1][j2], tsets[2][j3], tsets[3][j4], rIndex, cIndex);
                            rInd.push_back(rIndex);
                            cInd.push_back(cIndex);
                        }
                    }
                }
            }

            if (quartet != DUMMY) {
                printQuartet(quartet, txnsets[0], txnsets[1], txnsets[2], txnsets[3]);
            }
//        }
//    }

            return tsets[0].size() * tsets[1].size() * tsets[2].size() * tsets[3].size();
}

void TaxonSet::printWithTaxonSet(ostream &str, TaxonSet*& other, int numTaxa) {
    str << "{";
    vector<int> result(tset.size() + (other->tset).size());
    set_union(tset.begin(), tset.begin() + tset.size(), (other->tset).begin(), (other->tset).begin() + (other->tset).size(), result.begin());
    if (result[0] == 1) {
        printSet(str, result);
        str << "|";
        printSetComplement(str, result, numTaxa);
    }
    else {
        printSetComplement(str, result, numTaxa);
        str << "|";
        printSet(str, result);

    }
    str << "}";
}

/*int main(int argc, char *argv[]) {
    string top = "(1,(2,(3,(4,(5,6)))));";
//      string top = "(((1,4),2),(3,(5,6)));";
//    string top = "((1,2),((3,4),(5,6)));";
//    string top = "(((((1,2),3),4),5),6);";
//        string top = "(1,((2,3),4));";


      vector<string> topologies;
      topologies.push_back(top);
      Table *newTable = new TGM(topologies);
      newTable->addGeneCount(0,0,10);
      newTable->addGeneCount(0,1,10);
      newTable->addGeneCount(0,2,10);
      newTable->addGeneCount(0,3,10);

      TreeBuilder t;
      top = t.getTree(newTable, 6);
      cerr << "Final topology: " << top<< endl;
}*/
