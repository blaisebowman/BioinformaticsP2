Blaise Bowman, CIS 4930: Bioinformatics Project #2

To run this program, I transferred a2_1_bowman.cpp, a2_2_bowman.cpp, a2_3_bowman.cpp, and sequenceGenerator.cpp to Linux via WinSCP. 

Within Linux, I used the following commands, in order:

g++ sequenceGenerator.cpp
./a.out
	//The above generates 10 sequences of length 25, 100, and 250, which are written to 10seq25char.o1, 10seq100char.o2, and 10seq250char.o3 respectively. 
	//The above generates 5 sequences of length 25, which is written to 5seq25char.o4
	//The above generates 25 sequences of length 25, which is written to 25seq25char.o5
	//The above generates 50 sequences of length 25, which is written to 50seq25char.o6

g++ a2_1_bowman.cpp 
./a.out
	//The above generates a pair of t DNA string of length n, with a (k,d)-motif at a random location in each DNA string. 
	//Writes output to assignment2.fasta

g++ a2_2_bowman.cpp 
	//Compiles the program containing the Brute Force Motif Search Algorithm. 

./a.out assignment2.fasta 3 1
	//replace "assignment2.fasta" with a .fasta file to cross-check.
	//replace 3 (which is k) with a specified k-value.
	//replace 1 (which is d) with a specified number of mismatches.
	//writes output to a2.o2, which contains the generated motifs and the associated score.

g++ a2_3_bowman.cpp 
	//Compiles the program containing the Gibbs Sampling Algorithm. 

./a.out assignment2.fasta 3 1
	//replace "assignment2.fasta" with a .fasta file to cross-check.
	//replace 3 (which is k) with a specified k-value.
	//replace 1 (which is d) with a specified number of mismatches.
	//writes output to assignment2.o3, which contains the generated motifs and the associated score. 

g++ a2_2_bowman.cpp
	//recompile a2_2_bowman.cpp in order to run algorithm analysis for Brute Force Motif Search

time ./a.out 10seq25char.o1 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 25 characters in length is passed 

time ./a.out 10seq100char.o2 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 100 characters in length is passed 

time ./a.out 10seq250char.o3 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 250 characters in length is passed 

time ./a.out 5seq25char.o4 15 4
	//The above command determines the execution time of the algorithm with the 5 sequences of 25 characters in length is passed 

time ./a.out 25seq25char.o5 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 25 characters in length is passed 

time ./a.out 50seq25char.o6 15 14
	//The above command determines the execution time of the algorithm with the 50 sequences of 25 characters in length is passed

g++ a2_3_bowman.cpp
	//recompile a2_3_bowman.cpp in order to run algorithm analysis for Gibbs Sampling

time ./a.out 10seq25char.o1 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 25 characters in length is passed 

time ./a.out 10seq100char.o2 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 100 characters in length is passed 

time ./a.out 10seq250char.o3 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 250 characters in length is passed 

time ./a.out 5seq25char.o4 15 4
	//The above command determines the execution time of the algorithm with the 5 sequences of 25 characters in length is passed 

time ./a.out 25seq25char.o5 15 4
	//The above command determines the execution time of the algorithm with the 10 sequences of 25 characters in length is passed 

time ./a.out 50seq25char.o6 15 4
	//The above command determines the execution time of the algorithm with the 50 sequences of 25 characters in length is passed  
 
