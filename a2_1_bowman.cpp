#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <random>
#include <cstdlib>
#include <regex>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment #2

vector<string> generatorSeq(const int sequenceLength, const int numSequences) {
    //given a sequence length, generate a pair of DNA sequences of length sequenceLength
    //generate numSequences DNA sequences
    vector<string> generated;
    const string dnaChars = "actg";
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    uniform_int_distribution<> pos(0, dnaChars.size() - 1);
    string sequence1;
    for (int a = 0; a < numSequences; a++) {
        sequence1 = "";
        for (int i = 0; i < sequenceLength; i++) {
            sequence1 += dnaChars[pos(randomGenerator)];
        }
        generated.push_back(sequence1);
    }
    return generated;
}

string motifGenerator(int k) {
    const string dnaChars = "ATCG";
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    uniform_int_distribution<> pos(0, dnaChars.size() - 1);
    string motifString;
    for (int i = 0; i < k; i++) {
        motifString += dnaChars[pos(generator)];
    }
    return motifString;
}

string motifMismatch(string motif, int mismatch) {
    const string dnaChars = "ATCG";
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    uniform_int_distribution<> pos(0, dnaChars.size() - 1);
    string motifString;
    int track;
    vector<int> avoid;
    while (track < mismatch) {
        uniform_int_distribution<int> dist(0, motif.size() - 1); //generate random index between zero and
        int randIndex = dist(generator);
        char c = dnaChars[pos(generator)];
        if ((count(avoid.begin(), avoid.end(), randIndex))) {
            if (motif.at(randIndex) != c) {
                motif.at(randIndex) = c;
                track++;
                avoid.push_back(randIndex);
            }
        } else {
            motif.at(randIndex) = c;
            track++;
            avoid.push_back(randIndex);
        }
    }
    return motif;
}

vector<string> insertMotif(vector<string> DNA, int motifLength, int mismatch) {
    vector<string> updatedSequence; // the updated DNA sequence with the inserted MOTIF
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    string motifString;
    string usedMotif = motifGenerator(motifLength);
    for (int i = 0; i < DNA.size(); i++) {
        uniform_int_distribution<int> dist(0, DNA.at(i).size() - motifLength);
        int randIndex = dist(generator);
        string thisMotif = motifMismatch(usedMotif, mismatch);
        DNA.at(i).replace(randIndex, motifLength, thisMotif);
    }
    return DNA;
}

int main(int argc, char *argv[]) {
    int t = strtol(argv[1], nullptr, 0); //---->>> number of sequences
    int n = strtol(argv[2], nullptr, 0); //---->>> number of characters in sequences
    int k = strtol(argv[3], nullptr, 0); //---->>> length of motif
    int d = strtol(argv[4], nullptr, 0); //---->>> number of mismatches
    cout << t << " " << n << " " << k << " " << d <<endl;
    ofstream out;
    ofstream ofile("assignment2.fasta");
    vector<string> sequences = generatorSeq(n, t);
    if (ofile.is_open()) {
        for (int i = 0; i < sequences.size(); i++) {
            ofile << ">seq" + to_string(i + 1) << endl;
            ofile << sequences[i] << endl;
        }
    }
    ofile.close();
    ofile.clear();
    return 0;
}
