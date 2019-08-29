#include "nussinov_dp.hpp"

void b_task(nussinov_dp& _seq) {
	_seq._predict();
}

int main(int arcv, char* argc[]){
	std::ios_base::sync_with_stdio(false);
	nussinov_dp _seq(arcv, argc);
	std::thread sub_task(b_task, std::ref(_seq));
	std::string f_out;
	std::cout << "\nEnter a name for output file: ";
	std::cin >> f_out;
	sub_task.join();

	_seq._write(f_out);
	return EXIT_SUCCESS;
}
