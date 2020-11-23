#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <set>
#include <random>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment 2 Sequence Generator
//Modified from my version of Assignment 1.

vector <string> algorithmAnalysisv2 (int numSequences, int sequenceLength) {
    vector <string> sequences;
    const string dnaChars = "ATCG";
    random_device randomGenerator;
    mt19937 generator(randomGenerator());
    uniform_int_distribution<> pos(0, dnaChars.size() - 1);
    for (int seqNum = 0; seqNum <= numSequences; seqNum++) {
        string dnaSeq;
        for (int i = 0; i < sequenceLength; i++) {
            dnaSeq += dnaChars[pos(generator)];
        }
        sequences.push_back(dnaSeq);
    }
return sequences;
}

int main(int argc, char *argv []) {
    // Companion to a2_2_bowman.cpp and a2_3_bowman.cpp, used for Question 4 part A and part B to generate random string sequence pairs

    //START OF PART A
    ofstream ofile ("10seq25char.o1");
    vector <string> result1 = algorithmAnalysisv2(10, 25);
    if(ofile.is_open()) {
        for (int i = 1; i < result1.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result1[i] << endl;
        }
    }
    ofile.close();
    ofile.clear();

    ofile.open("10seq100char.o2");
    vector <string> result2 = algorithmAnalysisv2(10, 100);
    if(ofile.is_open()) {
        for (int i = 1; i < result2.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result2[i] << endl;

        }
    }
    ofile.close();
    ofile.clear();
    ofile.open("10seq250char.o3");
    vector <string> result3 = algorithmAnalysisv2(10, 250);
    if(ofile.is_open()) {
        for (int i = 1; i < result3.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result3[i] << endl;

        }
    }
    ofile.close();
    ofile.clear();

    //Start of Part B:

    ofile.open("5seq25char.o4");
    vector <string> result4 = algorithmAnalysisv2(5, 25);
    if(ofile.is_open()) {
        for (int i = 1; i < result4.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result4[i] << endl;

        }
    }
    ofile.close();
    ofile.clear();

    ofile.open("25seq25char.o5");
    vector <string> result5 = algorithmAnalysisv2(25, 25);
    if(ofile.is_open()) {
        for (int i = 1; i < result5.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result5[i] << endl;

        }
    }
    ofile.close();
    ofile.clear();

    ofile.open("50seq25char.o6");
    vector <string> result6 = algorithmAnalysisv2(50, 25);
    if(ofile.is_open()) {
        for (int i = 1; i < result6.size(); i++) {
            ofile << ">seq" << i << endl;
            ofile << result6[i] << endl;
        }
    }
    ofile.close();
    ofile.clear();
    return 0;
}

