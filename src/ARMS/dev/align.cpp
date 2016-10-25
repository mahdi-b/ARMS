#include <iostream>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <cstdlib>
#include <string>
#include <unordered_map>
using namespace seqan;

int main(int argc, char const ** argv) {

    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;     // align type

    if (argc != 3) {
	std::cout<<"Mising params";
        return 1;  // Invalid number of arguments.
    }
    std::unordered_map<std::string,TSequence> seqs;
    std::string input_file = argv[1];
    std::string seq_name = argv[2];
    TSequence seq = "";

    std::cout<< "Reading " + input_file + ".\n";
    // Open file and create RecordReader.
    std::fstream in(input_file, std::ios::binary | std::ios::in);
    RecordReader<std::fstream, SinglePass<> > reader(in);

    // Read file record-wise.
    std::string id;
    while (!atEnd(reader)) {
        if (!readRecord(id, seq, reader, Fasta())) {
            std::cout << "Error reading fasta";
    }
        seqs[id] = seq;
    }
    for ( auto it = seqs.begin(); it != seqs.end(); ++it ) {
        std::cout << it->first << "\n" << it->second;
        std::cout << std::endl;
    } 


    return 0;
}
