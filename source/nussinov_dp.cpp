#include "nussinov_dp.hpp"
nussinov_dp::nussinov_dp(int argv, char** argc) : hl_size{ 0 }, n_basepairs{ 0 } {
	if (argv != 3) {
		std::cerr << "argument error!\n";
		exit(EXIT_FAILURE);}
	std::ifstream seq_if_stream(argc[1], std::ios::in);
	std::ifstream func_if_stream(argc[2], std::ios::in);
	if (!seq_if_stream || !func_if_stream) {
		std::cerr << "file could not be opened!\n";
		exit(EXIT_FAILURE);}
	func_if_stream >> hl_size;
	char c1, c2;
	std::size_t val;
	while (func_if_stream >> c1 >> c2 >> val) 
		score_func[c1][c2] = val;
	std::string record;
	while(std::getline(seq_if_stream, record))
		if (record[0] == '>')
			seq_header = record;
		else
			seq += record;
	seq_if_stream.close();
	func_if_stream.close();
	scores = new score*[seq.size()];
	for (std::size_t i = 0; i < seq.size()-hl_size; ++i) {
		scores[i] = new score[seq.size()];
		scores[i][i+hl_size].val = 0;
		if(i!=0)
			scores[i][i+hl_size-1].val = 0; }
	assert(hl_size >= 0 && hl_size <= seq.size() - 2);
	structure = seq;
}
nussinov_dp::~nussinov_dp() {
	for (std::size_t i = 0; i < seq.size() - hl_size; ++i)
		delete[] scores[i];
	delete[] scores;
}
void nussinov_dp::_predict(void){
	std::size_t ind1, ind2, _max;
	for (std::size_t n = 1 + hl_size; n < seq.size(); ++n) { 
		for (std::size_t j = n; j < seq.size(); ++j) {
			std::size_t i = j - n ;  
			ind1 = i; ind2 = j - 1; _max = scores[i][j - 1].val;
			_max = scores[i + 1][j - 1].val + score_func[seq[i]][seq[j]] > _max ? ind1 = i + 1, ind2 = j - 1, (scores[i + 1][j - 1].val + score_func[seq[i]][seq[j]]) : _max;
			_max = scores[i + 1][j].val > _max ? ind1 = i+1, ind2 = j, scores[i+1][j].val : _max;
			scores[i][j].val = _max;
			scores[i][j].preds.push_back(std::pair<std::size_t, std::size_t>(ind1, ind2));
			for(std::size_t k = i + 1 + hl_size; k < j - hl_size; ++k)  
				if (scores[i][k].val + scores[k + 1][j].val > scores[i][j].val) { 
					scores[i][j].val = scores[i][k].val + scores[k + 1][j].val;
					scores[i][j].preds[0] = std::pair<std::size_t, std::size_t>(i, k);
					if (scores[i][j].preds.size() > 1)
						scores[i][j].preds[1] = std::pair<std::size_t, std::size_t>(k + 1, j);
					else scores[i][j].preds.push_back(std::pair<std::size_t, std::size_t>(k + 1, j));
				}
		}
	}
	_traceback(0, seq.size() - 1);   //   calling traceback from the upper right corner
	for (std::size_t i = 0; i < seq.size(); ++i)
		if (structure[i] != '(' && structure[i] != ')')
			structure[i] = '.';
}
void nussinov_dp::_traceback(std::size_t i, std::size_t j){
	if (i + hl_size < j) { 
		if (scores[i][j].preds.size() > 1) {
			_traceback(scores[i][j].preds[0].first, scores[i][j].preds[0].second);
			_traceback(scores[i][j].preds[1].first, scores[i][j].preds[1].second);  }
		else if (scores[i][j].preds[0].first == i && scores[i][j].preds[0].second == j - 1) 
			_traceback(i, j - 1);
		else if (scores[i][j].preds[0].first == i + 1 && scores[i][j].preds[0].second == j) 
			_traceback(i + 1, j);
		else {
			structure[i] = '(';
			structure[j] = ')';
			++n_basepairs;
			_traceback(i + 1, j - 1);  }
	}	
}
void nussinov_dp::_report() {
	std::cout << "\n\n\na sequence of size " << seq.size() << " has been read\n\n";
	std::cout << "minimum size of hairpin-loops: " << hl_size << '\n';
	std::cout << "\nscoring-function: \n";
	for (std::map<char, std::map<char, std::size_t>>::iterator it1 = score_func.begin(); it1 != score_func.end(); ++it1)
		for (std::map<char, std::size_t>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			std::cout << "   " <<it1->first << " - " << it2->first << " => " << it2->second << '\n';
	std::cout << "\n\nmax score: " << scores[0][seq.size() - 1].val << "\n";
	std::cout << "\nnumber of base pairs: " << n_basepairs << "\n";
}
void nussinov_dp::_write(std::string f_out){
	std::ofstream of_stream(f_out, std::ios::out);
	if (!of_stream) {
		std::cerr << "file could not be created: " << f_out << "\n";
		exit(EXIT_FAILURE); }
	of_stream << seq_header << "\n" << seq << "\n" << structure;
	of_stream.close();
	_report();
	std::cout << "\nthe structure as a dot-bracket format has been written to '" << f_out << "'\n\n\n";
}
