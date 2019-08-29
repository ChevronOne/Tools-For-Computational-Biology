
#include "commons.hpp"

std::vector<std::tuple<std::size_t, std::size_t>> thread_intervals(unsigned tks_per_thr,unsigned samS, unsigned thrN){
    std::vector<std::tuple<std::size_t, std::size_t>> intervals;
    std::size_t sum{0};
    std::size_t t{0};

    std::tuple<std::size_t, std::size_t> tup;
    std::get<0>(tup) = t;
    intervals.push_back(tup);

    for(unsigned i = 0, j = samS - 1; i < samS -1; ++i, j-=1){
        sum += j;
        if(i == samS-2){
            std::get<1>(intervals[t]) = i + 1;
            break;
        }
        if (sum >= tks_per_thr){
            std::get<1>(intervals[t]) = i;
            t += 1;
            std::get<0>(tup) = i+1;
            intervals.push_back(tup);
            sum = 0;
        }
    }
    return intervals;
}



void* _task(void* argv)
{

    std::map<std::string, std::vector<std::string>>* samples =  (std::map<std::string, std::vector<std::string>>*) (((void**)argv)[0]);
    uint64_t** distMat =  (uint64_t**) (((void**)argv)[1]);
    uint64_t* size_ =  (uint64_t*)  (((void**)argv)[2]);
    std::map<std::string, std::map<std::string, std::size_t>>* penalty =  (std::map<std::string, std::map<std::string, std::size_t>>*) (((void**)argv)[3]);
    std::tuple<std::size_t, std::size_t>* range =  (std::tuple<std::size_t, std::size_t>*) (((void**)argv)[4]);
    std::mutex* synchro =  (std::mutex*) (((void**)argv)[5]);

	for (auto i = std::get<0>(*range); i <= std::get<1>(*range); ++i)
	{
		for(auto j = i+1; j < *size_; ++j)
		{
			for(std::map<std::string, std::vector<std::string>>::iterator it = samples->begin(); it != samples->end(); ++it){
                
				distMat[i][j]+=(*penalty)[it->second[i]][it->second[j]];
			}
            distMat[j][i] = distMat[i][j];
			(*synchro).lock();
			std::cout << "dist < "<< std::setw(2) << i << " , " << std::setw(2) << j << " >: " << distMat[i][j] << "\n";
			(*synchro).unlock();
		}
	}

    pthread_exit(0);
}


void execution_time(std::chrono::time_point<std::chrono::steady_clock> &start_t, const std::string &st)
{
	std::cout << "\n"<< st <<":\n";
	std::cout << std::setw(7)
			  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_t).count()
			  << " Milli-Seconds" << std::endl;
	std::cout << std::setw(7)
			  << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_t).count()
			  << " Seconds" << std::endl;
	std::cout << std::setw(7)
			  << std::chrono::duration_cast<std::chrono::minutes>(std::chrono::steady_clock::now() - start_t).count()
			  << " Minutes" << std::endl;
}

unsigned tasks_num(const unsigned &thr){
    unsigned tsk{0};
    for(std::size_t i = 1; i < thr; ++i)
        tsk+=i;
    return tsk;
}


