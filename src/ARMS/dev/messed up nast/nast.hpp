#ifndef NAST_HPP
#define NAST_HPP
#include <iostream>
#include <fstream>
#include <string>
#include<cstdlib>
#include<vector>
using namespace std;
/*
 *  nast.hpp
 *
 *
 *  Created by Pat Schloss on 12/17/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is my implementation of the NAST (nearest alignment space termination) algorithm as described in:
 *
 *	DeSantis TZ, Hugenholtz P, Keller K, Brodie EL, Larsen N, Piceno YM, Phan R, & Anderson GL.  2006.  NAST: a multiple
 *		sequence alignment server for comparative analysis of 16S rRNA genes.  Nucleic Acids Research.  34:W394-9.
 *
 *	To construct an object one needs to provide a method of getting a pairwise alignment (alignment) and the template
 *	and candidate sequence that are to be aligned to each other.
 *
 */



/**************************************************************************************************/

class Nast {

public:
	Nast(string, string);
	float getSimilarityScore();
	int getMaxInsertLength();

	string alignment;
	string candidateSeq;
	string templateSeq;
private:
	void pairwiseAlignSeqs();
	void regapSequences();
	void removeExtraGaps(string&, string, string);


	int maxInsertLength;
};

/**************************************************************************************************/

#endif