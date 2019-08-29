#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <thread>
#include <cassert>


struct score{
	std::size_t val;         // score for the current position
	std::vector<std::pair<std::size_t, std::size_t>> preds;  // the predecessors of the current position
};


class nussinov_dp{
	
	std::map<char, std::map<char, std::size_t>> score_func;     // score function
	score **scores;        // scores mat
	std::string seq, seq_header, structure;   // the sequence, it's header and the predicted structure
	std::size_t hl_size;    // the minimum allowed hairpin-loops size
	std::size_t n_basepairs;  // number of base pairs



public:
	nussinov_dp(int, char**);
	void _predict(void);
	void _traceback(std::size_t, std::size_t);
	void _write(std::string);
	void _report(void);
	virtual ~nussinov_dp();


};
