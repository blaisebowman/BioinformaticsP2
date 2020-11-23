#include<bits/stdc++.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <random>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment #2, Part 3.
//Gibbs Sampling Algorithm, using Laplace's Rule of Succession

int randomInd(int s, int e) {
    static random_device rng;
    static std::mt19937 eng(std::time(nullptr));
    std::uniform_int_distribution<int> val(s, e);
    return val(eng);
}

pair<vector<string>, double> gibbs_sampling_algorithm(vector<string> DNA, int k, int d, int t, double num) {
    vector<string> motifs, removedMotifs;
    vector<int> startingPositions;
    static random_device randomGenerator;
    static mt19937 generator(randomGenerator());
    for (int i = 0; i < DNA.size(); i++) {
        transform(DNA[i].begin(), DNA[i].end(), DNA[i].begin(), ::tolower);
    }
    for (int i = 0; i < DNA.size(); i++) {
        int start = randomInd(0, DNA[0].size() - k - 1);
        startingPositions.push_back(start); //picks random starting index for motifs for each sequence
    }
    for (int i = 0; i < DNA.size(); i++) {
        string m = DNA[i].substr(startingPositions[i], k);
        string bef = DNA[i].substr(0, startingPositions[i]);
        string after = DNA[i].substr(startingPositions[i] + k, DNA[i].size() - startingPositions[i] + k + 1);
        transform(m.begin(), m.end(), m.begin(), ::toupper);
        DNA[i] = bef + m + after; //rejoins string with new motif in uppercase
        motifs.push_back(m); //assigns the motifs
    }
    static uniform_int_distribution<int> str(0, DNA.size() - 1);
    int randomKmer = str(generator);
    vector<string> removedSequences;
    removedSequences.push_back(DNA[randomKmer]);
    removedMotifs.push_back(motifs[randomKmer]);
    motifs.erase(motifs.begin() + randomKmer); //erases random k-mer from motif vector
    vector<int> aNum, cNum, gNum, tNum;
    for (int i = 0; i < motifs[0].size(); i++) {
        int a = 0;
        int cV = 0;
        int g = 0;
        int tt = 0;
        for (int j = 0; j < motifs.size(); j++) {
            if (motifs[j][i] == 'a' || motifs[j][i] == 'A') {
                a++;
            }
            if (motifs[j][i] == 'c' || motifs[j][i] == 'C') {
                cV++;
            }
            if (motifs[j][i] == 'g' || motifs[j][i] == 'G') {
                g++;
            }
            if (motifs[j][i] == 't' || motifs[j][i] == 'T') {
                tt++;
            }
        }
        aNum.push_back(a);
        cNum.push_back(cV);
        gNum.push_back(g);
        tNum.push_back(tt);
    }
    /*cout << "Score Matrix: " << endl;
    cout << "A: ";*/
    for (int i = 0; i < motifs[0].size(); i++) {
        aNum[i]++;
        //cout << aNum[i] << "  ";
    }
    /*cout << endl;
    cout << "C: ";*/
    for (int i = 0; i < motifs[0].size(); i++) {
        cNum[i]++;
        //cout << cNum[i] << "  ";
    }
    /*cout << endl;
    cout << "G: ";*/
    for (int i = 0; i < motifs[0].size(); i++) {
        gNum[i]++;
        //cout << gNum[i] << "  ";
    }
    /*cout << endl;
    cout << "T: ";*/
    for (int i = 0; i < motifs[0].size(); i++) {
        tNum[i]++;
        //cout << tNum[i] << "  ";
    }
    /*cout << endl;
    cout << "Profile Matrix: " << endl;
    cout << "A: ";*/
    /*for (int i = 0; i < motifs[0].size(); i++) {
        if (aNum[i] != 0) {
            //cout << ((double)(aNum[i]) / 4) << " ";
            cout << ((double) (aNum[i]) / 8) << " ";
        } else {
            cout << 0 << "  ";
        }
    }
    cout << endl;
    cout << "C: ";
    for (int i = 0; i < motifs[0].size(); i++) {
        if (cNum[i] != 0) {
            //cout << ((double)(cNum[i]) / 4) << " ";
            cout << ((double) (cNum[i]) / 8) << " ";
        } else {
            cout << 0 << "  ";
        }
    }
    cout << endl;
    cout << "G: ";
    for (int i = 0; i < motifs[0].size(); i++) {
        if (gNum[i] != 0) {
            //cout << ((double)(gNum[i]) / 4) << " ";
            cout << ((double) (gNum[i]) / 8) << " ";
        } else {
            cout << 0 << "  ";
        }
    }
    cout << endl;
    cout << "T: ";
    for (int i = 0; i < motifs[0].size(); i++) {
        if (tNum[i] != 0) {
            //cout << ((double)(tNum[i]) / 4) << " ";
            cout << ((double) (tNum[i]) / 8) << " ";
        } else {
            cout << 0 << "  ";
        }
    }
    cout << endl;*/

    vector<string> splitTheMotifs;
    for (int i = 0; i < removedMotifs.size(); i++) {
        for (int r = 0; r <= DNA[0].size() - k; r++) {
            splitTheMotifs.push_back(removedSequences[i].substr(r, k)); //splits the motifs of
        }
    }
    vector<double> probabilities;
    map<string, double> map;
    for (int i = 0; i < splitTheMotifs.size(); i++) {
        double prob = 1.0; //use a double to keep track of the score 
        for (int j = 0; j < k; j++) {
            if (splitTheMotifs[i][j] == 'a' || splitTheMotifs[i][j] == 'A') {
                if (aNum[j] != 0 && prob != 0.0) {
                    //prob *= ((double)(aNum[j]) / 4);
                    prob *= ((double) (aNum[j]) / 8);
                } else {
                    prob = 0;
                    break;
                }
            }
            if (splitTheMotifs[i][j] == 'c' || splitTheMotifs[i][j] == 'C') {
                if (cNum[j] != 0 && prob != 0.0) {
                    //prob *= ((double)(cNum[j]) / 4);
                    prob *= ((double) (cNum[j]) / 8);
                } else {
                    prob = 0;
                    break;
                }
            }
            if (splitTheMotifs[i][j] == 'g' || splitTheMotifs[i][j] == 'G') {
                if (gNum[j] != 0 && prob != 0.0) {
                    //prob *= ((double)(gNum[j]) / 4);
                    prob *= ((double) (gNum[j]) / 8);

                } else {
                    prob = 0;
                    break;
                }
            }
            if (splitTheMotifs[i][j] == 't' || splitTheMotifs[i][j] == 'T') {
                if (tNum[j] != 0 && prob != 0.0) {
                    //prob *= ((double)(tNum[j]) / 4);
                    prob *= ((double) (tNum[j]) / 8);

                } else {
                    prob = 0;
                    break;
                }
            }
        }
        probabilities.push_back(prob);
        map.insert(pair<string, double>(splitTheMotifs[i], prob));
    }

    /*for (int i = 0; i < probabilities.size(); i++) {
        cout << "(" << probabilities[i] << ") " << splitTheMotifs[i] << "  ";
    }
    cout << endl;
    cout << "Before insertion, the motifs are: " << endl;
    for (int i = 0; i < motifs.size(); i++) {
        cout << motifs[i] << endl;
    }*/

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> dis(probabilities.begin(), probabilities.end());
    std::map<int, int> mp;
    int val;
    for (int i = 0; i < 1; i++) {
        ++mp[dis(gen)];
    }
    string newS;
    vector<int> addedIn;
    string updated = "";
    transform(removedSequences[0].begin(), removedSequences[0].end(), removedSequences[0].begin(), ::tolower);
    //cout << removedMotifs[0] << endl;
    for (int i = 0; i < splitTheMotifs.size(); i++) {
        transform(splitTheMotifs[i].begin(), splitTheMotifs[i].end(), splitTheMotifs[i].begin(), ::tolower);
    }
    double score;
    for (auto p: mp) {
       /* cout << p.first << " generated " << p.second << " times\n";
        cout << "The motif associated with index " << p.first << " is " << splitTheMotifs[p.first];
        cout << "\n" << "The probability/score of this motif is: " << probabilities[p.first] << " " << endl;*/
        score = probabilities[p.first];
        transform(splitTheMotifs[p.first].begin(), splitTheMotifs[p.first].end(), splitTheMotifs[p.first].begin(),
                  ::toupper);
        motifs.insert(motifs.begin() + randomKmer, (splitTheMotifs[p.first]));
    }

    if (score != num) { //if score does not remain the same, repeat.
        gibbs_sampling_algorithm(DNA, k, d, t, score);
    }
    for (int i = 0; i < splitTheMotifs.size(); i++) {
        if (i == 0) {
            newS.append(splitTheMotifs[i].substr(0, splitTheMotifs[i].size()));
        } else if (i == splitTheMotifs.size()) {
            newS.append(splitTheMotifs[i].substr(0, splitTheMotifs[i].size()));
        } else {
            string tmp = string{splitTheMotifs[i].back()};
            newS.append(tmp);
        }
    }
    string m2 = motifs[randomKmer];
    transform(m2.begin(), m2.end(), m2.begin(), ::toupper);
    string bef2 = removedSequences[0].substr(0, startingPositions[randomKmer]);
    string after2 = removedSequences[0].substr(startingPositions[randomKmer] + k,
                                               removedSequences[0].size() - startingPositions[randomKmer] + k + 1);
    DNA[randomKmer] = newS;
    transform(DNA[randomKmer].begin(), DNA[randomKmer].end(), DNA[randomKmer].begin(), ::toupper);
    return make_pair(motifs, score);
}

int main(int argc, char *argv[]) {
    //TAKE AS COMMAND LINE ARGUMENTS!
    //ifstream file ("assignment2Updated.fasta");
    ifstream file(argv[1]); //--->> for CMD/LINUX
    string line, line2;
    map<string, string> lines;
    vector<string> sequence;
    vector<string> data;
    if (file.is_open()) {
        string val;
        while (getline(file, line)) {
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.find('>') != string::npos) { //if the line begins with '>'
                sequence.push_back(line);
            } else {
                data.push_back(line);
            }
        }
        for (int i = 0; i < sequence.size(); i++) {
            lines.insert(pair<string, string>(sequence[i], data[i]));
        }
        file.close();
    }
    vector<string> sequences;
    for (const auto &s  : lines) {
        sequences.push_back(s.second);
    }
    long kL = strtol(argv[2], NULL, 10); //---->>> number of sequences
    int k = kL;
    long dL = strtol(argv[3], NULL, 10); //---->>> number of characters in sequence
    int d = dL;
    double num;
    pair<vector<string>, double> result;
    result = gibbs_sampling_algorithm(sequences, k, d, sequences.size(), num);
    ofstream out;
    ofstream ofile("assignment2.o3");
    if (ofile.is_open()) {
        for (int i = 0; i < result.first.size(); i++) {
            ofile << result.first[i] << endl;
        }
        ofile << result.second << endl;
    }
    ofile.close();
    ofile.clear();
    return 0;
}