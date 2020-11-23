#include<bits/stdc++.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>

using namespace std;

//Blaise Bowman, CIS 4930 Bioinformatics Assignment #2, Part 2.
//Brute Force Motif Search

vector <int> next_leaf(vector<int> digits,int t,int k){
    for(int i=t; i>=0; i--){
        if(digits[i] < k){
            digits[i]++; //increments least significant digit
            return digits;
        }
        digits[i] = 1;
    }
    return digits;
}

int characterCounter(char c,int j,vector<int> s,vector <string> DNA,int t){
    int cnt = 0;
    for(int i=0; i<t; i++){
        if((DNA[i][s[i]+j] == c) || ((DNA[i][s[i]+j] == c-32))){ //c-32 is equal to A, C, G, or T
            cnt++; //increments the count of the character passed in as an argument.
        }
    }
    return cnt;
}

pair<string,int> score(vector<int> s,vector<string> DNA,int t,int k){
    int sc = 0;
    string consensus_motif;
    pair <string,int> p;
    for(int col = 0; col < k; col++){
        int max_count = 0;
        char max_char;
        int cnt = characterCounter('a',col,s,DNA,t);
        if(cnt > max_count){
            max_count = cnt;
            max_char = 'a';
        }
        cnt = characterCounter('c',col,s,DNA,t);
        if(cnt > max_count){
            max_count = cnt;
            max_char = 'c';
        }
        cnt = characterCounter('g',col,s,DNA,t);
        if(cnt > max_count){
            max_count = cnt;
            max_char = 'g';
        }
        cnt = characterCounter('t',col,s,DNA,t);
        if(cnt > max_count){
            max_count = cnt;
            max_char = 't';
        }
        sc += max_count;
        consensus_motif.push_back(max_char);
    }
    p = make_pair(consensus_motif,sc);
    return p;
}

pair<string,vector<int>> brute_force_motif_search(vector <string> DNA,int t,int n,int k){
    int best_score;
    string consensus_motif;
    vector <int> s(t);
    vector <int> z(t);
    vector <int> all_motifs(t);
    pair <string,int> p;
    pair <string,vector<int>> final_result;
    for(int i=0; i< t; i++) {
        s[i] = z[i] = 1;
    }
    p = score(s,DNA,t,k);
    consensus_motif = p.first;
    best_score = p.second;
    all_motifs = s;
    while(true){
        if(s == z) {
            break;
        }
        vector <int> tmp = next_leaf(s,t,n-k+1);
        s.insert( s.end(), tmp.begin(), tmp.end());
        //s += next_leaf(s,t,n-k+1);
        pair<string,int> q;
        q = score(s,DNA,t,k);
        if(q.second > best_score){
            consensus_motif = q.first;
            best_score = q.second;
            all_motifs = s;
        }
        if(s == z) {
            break;
        }
    }
    final_result = make_pair(consensus_motif,all_motifs);
    return final_result;
}

vector <string> perString (pair <string, vector <int>> result, int t,  vector <string> DNA, int k){
    vector <string> motifs;
    for (int i = 0; i < t; i++){
        motifs.push_back(DNA[i].substr(result.second[i],k));
    }
    return motifs;
}

int main(int argc, char *argv[]){
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
    pair<string,vector<int>> result;
    long kL = strtol(argv[2], NULL, 10); //---->>> number of sequences
    int k = kL;
    long dL = strtol(argv[3], NULL, 10); //---->>> number of characters in sequences
    int d = dL;
    result = brute_force_motif_search(sequences,sequences.size(),sequences[0].size(),k);
    vector <string> motifs = perString(result,sequences.size(),sequences,k);
    ofstream out;
    ofstream ofile("assignment2.o2");
    if(ofile.is_open()){
        for (const auto &x: motifs){
            ofile << x << endl;
        }
    }
    ofile.close();
    ofile.clear();
    return 0;
}