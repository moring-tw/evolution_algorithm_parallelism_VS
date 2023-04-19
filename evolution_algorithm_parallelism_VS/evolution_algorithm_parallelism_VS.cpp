// evolution_algorithm_parallelism_VS.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//
#include <omp.h>
#include <stdlib.h>

#include <map>

#include "headers/header.h"
#include "headers/cell.h"
#include "headers/publicfunction.h"



int main()
{
	size_t num_of_run = 1;

	evolutionalgorithm::Variables v;
	evolutionalgorithm::ReadConfig rc;
	std::ifstream config_file_path;
	string config_file_path_string;
	config_file_path.open("option.txt", std::ios::in);
	std::getline(config_file_path, config_file_path_string);
	std::cout << config_file_path_string << std::endl;
	system("pause");
	rc.NextRun(v, config_file_path_string);
	system("pause");
	double start, end;
	int problem_counter = 0;
	double* allBest;
	double meanBest = 0;
	string filepath = CreateFolder();
	_CrtCheckMemory();

	for (problem_counter = 0; problem_counter < v.GetProblemAmount(); ++problem_counter) {
		{
			if (v.IsProblemExecuting(problem_counter)) {
				evolutionalgorithm::ExecutiveVariables ev(v); // clone the variables v
				num_of_run = ev.GetRun();
				allBest = (double*)malloc(sizeof(double) * num_of_run);
				meanBest = 0;
				int thread_id = omp_get_thread_num();
				std::cout << "P" << problem_counter - 28 << " @ T" << thread_id << " " << std::endl;
				Benchmarks* fp;
				if (problem_counter > 28 && problem_counter < 49) {
					fp = FunctionChooser(problem_counter, ev);
				}
				else {
					fp = FunctionChooser(problem_counter, ev);
				}
				double now_best = DBL_MAX;
				start = clock();
				evolutionalgorithm::differentialevolution::ProcessController pc;
				std::vector<size_t> success_run_counter(5, 0);
				std::vector<size_t> total_found_optima(5, 0);
				const std::vector<size_t> num_of_optima{ 2, 5, 1, 4, 2, 18, 36, 81, 216, 12, 6, 8, 6, 6, 8, 6, 8, 6, 8, 8 };
				int run_counter = 0;
				for (run_counter = 0; run_counter < num_of_run; ++run_counter) {
					std::cout << "Problem :" << problem_counter - 28 << " @ Run " << run_counter + 1 << ": ";
					if (run_counter == num_of_run) break;
					now_best = DBL_MAX;
					pc.SetExecutiveVariables(ev);
					pc.SetOutputFileStream(filepath, std::to_string(problem_counter) + "R" + std::to_string(run_counter + 1), "csv");
					fp->nextRun();
					pc.NextRun();
					pc.SetBenchmark(fp);
					pc.Execute();
					double best = pc.GetBest();
					if (allBest)
						allBest[run_counter] = best;
					else
						throw invalid_argument("out of memory.");
					meanBest += best;
					if (best < now_best) {
						now_best = best;
					}
					std::vector<size_t>& res = pc.GetHowManyOptima();
					std::cout << std::scientific << now_best << std::defaultfloat << " ";
					for (size_t result_counter = 0; result_counter < res.size(); ++result_counter) {
						double accuracy = 0.1 / pow(10, result_counter);
						cout << accuracy << ":"
							<< res[result_counter]
							<< " ";
						if (num_of_optima[problem_counter - 29] == res[result_counter])
							++success_run_counter[result_counter];
						total_found_optima[result_counter] += res[result_counter];
					}
					std::cout << std::endl;
				}
				end = clock();
				std::cout << "Problem :" << problem_counter - 28 << std::endl;
				meanBest /= num_of_run;
				for (size_t result_counter = 0; result_counter < success_run_counter.size(); ++result_counter) {
					std::cout << "PR: " << (double)total_found_optima[result_counter] / ((double)num_of_run * (double)num_of_optima[problem_counter - 29]) << " SR: " << (double)success_run_counter[result_counter] / (double)num_of_run << std::endl;
				}
				std::cout << std::scientific;
				std::cout << "Mean: " << meanBest << " STD: " << evolutionalgorithm::PublicFunction::CalculateStandardDeviation(allBest, meanBest, num_of_run) << std::endl << std::endl;
				std::cout << std::defaultfloat;
				//std::cout << "Time(ms):" << (end - start) << std::endl << std::endl;
				delete fp;
				free(allBest);
			}
		}
	}
	std::system("PAUSE");
}


// 執行程式: Ctrl + F5 或 [偵錯] > [啟動但不偵錯] 功能表
// 偵錯程式: F5 或 [偵錯] > [啟動偵錯] 功能表

// 開始使用的提示: 
//   1. 使用 [方案總管] 視窗，新增/管理檔案
//   2. 使用 [Team Explorer] 視窗，連線到原始檔控制
//   3. 使用 [輸出] 視窗，參閱組建輸出與其他訊息
//   4. 使用 [錯誤清單] 視窗，檢視錯誤
//   5. 前往 [專案] > [新增項目]，建立新的程式碼檔案，或是前往 [專案] > [新增現有項目]，將現有程式碼檔案新增至專案
//   6. 之後要再次開啟此專案時，請前往 [檔案] > [開啟] > [專案]，然後選取 .sln 檔案

Benchmarks* FunctionChooser(unsigned serial_number, evolutionalgorithm::ExecutiveVariables& ev) {
	//Benchmarks* fp;
	if (serial_number > 0 && serial_number < 29) {
		Benchmarks* fp = new CEC2013Benchmarks(ev.GetDE_Dim(), serial_number);
		return fp;
	}
	else if (serial_number >= 29 && serial_number <= 48) {
		CEC2013* fp;
		fp = new CEC2013(serial_number - 28);
		ev.UpdateDE_Dim(fp->get_dimension());
		ev.SetDE_MaxFE(fp->get_maxfes());
		for (size_t dim_counter = 0; dim_counter < ev.GetDE_Dim(); ++dim_counter) {
			ev.SetBoundMin(dim_counter, fp->get_lbound(dim_counter));
			ev.SetBoundMax(dim_counter, fp->get_ubound(dim_counter));
		}
		return fp;
	}
	Benchmarks* fp;
	switch (serial_number) {
	case 0:
		fp = new F1();
		return fp;
		break;
	default:
		fp = nullptr;
		std::cout << "Unknown function number." << std::endl;
		return fp;
	}
}

string CreateFolder() {
	char dir_buf[std::size("Www Mmm dd hh:mm:ss yyyy\n")];
	auto t = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(t);
	std::tm* now_time = new tm;
	ctime_s(&dir_buf[0], 26, &end_time);

	localtime_s(now_time, &end_time);
	std::strftime(&dir_buf[0], 26, "./%Y%m%dT%H%M%S", now_time);
	wchar_t wtext[23];
	size_t out_size;
	mbstowcs_s(&out_size, wtext, dir_buf, 23);
	if (!CreateDirectory(wtext, NULL) && ERROR_ALREADY_EXISTS != GetLastError()) {
		std::cout << GetLastError() << " " << out_size << " ERROR occured when trying to create a folder. " << wtext << std::endl;
		exit(56789);
	}
	delete now_time;
	string res = string(dir_buf);
	std::cout << res << std::endl;
	return res;
}