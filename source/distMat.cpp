
/*
* main.cpp
*
*  Created on: Feb 25, 2019
*      Author: mogsem07
*/

/*
	** this program calculate a distance matrix between a sample in the first file and other samples in the second file.
	** the program considers only the SNPs on the chromosome 1 'chr1'.  but it can be easily modified to consider all chromosomes together, or an other arbitrary chromosome

	* the firs argument is a vcf file contains the SNPs of only one sample, if the file contains SNPs of more than one sample, the SNPs of the first sample will be read, and the rest will be ignored 
	* the second argument is a vcf file contains the SNPs of other samples, the number of samples in this file is unlimited 
	* the third argument is a penalty-function. we give this as an option, so we can apply weighted distance penalty, e.g there is distance between 1|1 and 0|1 , and also a distance between 1|1 and 0|0 , but the two distances are not equal. this can be defined as desired from the last column in the penalty-function file.  an example for how a penalty-function should look like is also submitted
	* the fourth argument is the output file
	* the fifth argument is optional, and it defines the number of threads to be used. More than 221 will be ignored. Default is 21


*/


#include "commons.hpp"
#include "_thread.hpp"

extern void execution_time(std::chrono::time_point<std::chrono::steady_clock> &, const std::string&);
extern std::vector<std::tuple<std::size_t, std::size_t>> thread_intervals(unsigned ,unsigned, unsigned);
extern void* _task(void*);
extern unsigned tasks_num(const unsigned&);

int main(int argc, char **argv)
{
	std::ios_base::sync_with_stdio(false);
	std::size_t thrN{DEF_THR_NUM};
	if (argc < 5 || argc > 6){                        // check validation of arguments
		std::cerr << "argument error!" << std::endl;
		exit(EXIT_FAILURE);
	}
	else if (argc == 6){
		std::string thr = argv[5];
		for (char c : thr)
			if (!isdigit(c)){
				std::cout << "Invalid argument: " << thr << std::endl;
				exit(EXIT_FAILURE);
			}
		thrN = atoi( argv[5] );
		thrN = thrN > 0 ? ( thrN < MAX_THR_NUM ? thrN : MAX_THR_NUM ) : DEF_THR_NUM ;
	}

	std::vector<std::ifstream> inputs(3);
	for (auto i = 1; i < 4; ++i)                // check if the files are readable
	{
		inputs[i - 1].open(argv[i], std::ios::in);
		if (!inputs[i - 1])
		{
			std::cerr << "file could not be opened: " << argv[i] << "\n";
			exit(EXIT_FAILURE);
		}
	}

	auto start_t = std::chrono::steady_clock::now();

	std::string record;
	std::size_t size_{0}, i, j;
	std::vector<std::string> headers;
	headers.push_back("Reference");
	while (std::getline(inputs[1], record))
	{
		if (record.substr(0, 6) == "#CHROM")     // count the number of samples 
		{
			std::stringstream st_stream(record);
			std::string temp;
			for (i = 1; st_stream.good(); ++i)
			{
				st_stream >> temp;
				if (i >= 10)
				{
					headers.push_back(temp);
					++size_;
				}
			}
			break;
		}
	}
	

	std::map<std::string, std::string> ref;
	std::string chrom, pos, id, rf, alt, snp;
	std::stringstream st_stream;

	while (std::getline(inputs[0], record))            // read the SNPs in first file
	{
		st_stream.str(record);
		st_stream >> chrom;
		if (chrom == "chr1")
		{
			st_stream >> pos >> id >> rf >> alt;
			if(rf.size()==1 && alt.size() == 1 && rf != "." && alt != ".") // ignore INDELs
			{
				for (i = 6; st_stream.good(); ++i)
					{
						st_stream >> snp;
						if (i == 10)         //  the program expect only one sample in the first file, if it contains more than one sample, the SNPs of the first sample will be read, and the rest will be ignored
							{
								ref.insert(std::pair<std::string, std::string>(chrom + pos, snp.substr(0, 3).erase(1, 1)));
								break;
							}
					}
			}
		}
		st_stream.clear();
	}
	inputs[0].close();
	size_+=1;
	std::vector<std::string> data(size_, "00");
	std::map<std::string, std::vector<std::string>> samples; 
	std::map<std::string,std::string>::iterator query1;
	bool flag{false};
	while (std::getline(inputs[1], record)){       // read the SNPs of other samples
		st_stream.str(record);		
		st_stream >> chrom ;
		if (chrom == "chr1"){
			st_stream >> pos >> id >> rf >> alt;
			if(rf.size()==1 && alt.size() == 1 && rf != "." && alt != "."){
				for (i = 6; st_stream.good(); ++i){
					st_stream >> snp;
					if (i >= 10){                // the number of samples in the second file can be any, they all will be considered 
						data[i - 9] = snp.substr(0, 3).erase(1, 1);
						if(data[i-9] != "00")    // check if all samples don't have the current SNP
							flag = true;
					}
				}
				query1 = ref.find(chrom+pos);
				if(!flag && query1 == ref.end())  // check at the reference as well, if it does not have it either, then skip it!
					flag = false;
				else{
					samples.insert(std::pair<std::string, std::vector<std::string>>(chrom + pos, data));
					flag = false;
				}
			}
		}
		st_stream.clear();
	}
	inputs[1].close();


	std::map<std::string, std::vector<std::string>>::iterator query2;
	std::vector<std::string> temp(size_,"00");

	for(std::pair<const std::string, std::string> it:ref){         //  merge SNPs of first and second files
		query2 = samples.find(it.first);
		if(query2 == samples.end()){
			temp[0] = it.second;
			samples.insert(std::pair<std::string, std::vector<std::string>>(it.first, temp));
		}else{
			query2->second[0] = it.second;
		}
	}


	std::map<std::string, std::map<std::string, std::size_t>> penalty;  //  load the penalty function
	std::string sam1, sam2;
	std::size_t pen;
	while (inputs[2] >> sam1 >> sam2 >> pen)
		penalty[sam1][sam2] = pen;
	inputs[2].close();


	// ==========================================

	if(thrN > std::floor((2.5*(size_))/4))
        thrN = std::floor((2.5*(size_))/4);
	
	std::size_t tskN = tasks_num(size_);
	std::vector<std::tuple<std::size_t, std::size_t>> intervals = thread_intervals(tskN/thrN, size_, thrN);  // create an interval of jobs for each thread to be assigned to.
	thrN = intervals.size();
    
	// ==========================================

	uint64_t **distMat = new uint64_t*[size_];  // create a distance matrix
	for(i = 0; i < size_; i++){
		distMat[i] = new uint64_t[size_];
		memset(distMat[i], 0, sizeof(uint64_t)*(size_));
	}

	// =================================
	std::mutex synchro;
	_thread *ths = new _thread[thrN];    // create N threads
	for(i =0; i < thrN; ++i)
		ths[i] = _thread(_task, 6, &samples, distMat, &size_, &penalty, &intervals[i], &synchro);
	
	for(i =0; i < thrN; ++i)
		ths[i].join(NULL);

	// ======================================================

	std::ofstream of_stream(argv[4], std::ios::out);      // writ the calculated distance matrix into the output file
	of_stream << "            \t";
	for (i = 0; i <  size_; ++i)
		of_stream << std::setw(12) << headers[i] << "\t";
	of_stream << std::endl;
	for (i = 0; i <  size_; ++i) {
		of_stream << std::setw(12) <<  headers[i] <<"\t";
		for (j = 0; j <  size_; ++j) {
			of_stream <<std::setw(12) <<  distMat[i][j] << "\t";
		}
		of_stream << std::endl;
	}
	of_stream.close();
 
	// ====================================

	execution_time(start_t, "Calculation-time");

	//  Clean up!
	for(i =0; i < thrN; ++i)
	{
		ths[i].cleanUP();
	}
	delete[] ths;

	for(i = 0; i < size_; i++)
		delete[] distMat[i];
	delete[] distMat;

    // ======================

	return EXIT_SUCCESS;
}
