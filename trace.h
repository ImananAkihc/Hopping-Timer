#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
using namespace std;

#define READ_LEN 8

void loadTrace(vector<pair<uint32_t, double>>&vec, const char*filename, int totalflow)
{
	set<uint32_t> idset;
	double time;
	cout << "open " << filename << endl;
	FILE* pf = fopen(filename, "rb");
	if (!pf) {
		cout << filename << "not found." << endl;
		system("pause");
		exit(-1);
	}
	vec.clear();
	char read_key[16];
	uint32_t ret = 0;
	while (fread(read_key, 1, 16, pf))
	{
		time = *(double*)read_key;
		uint32_t key = *(uint32_t*)(read_key+8);
		idset.insert(key);
		vec.push_back(pair<uint32_t, double>(key, time));
		ret++;
		if (ret == totalflow)
			break;
	}
	fclose(pf);
	cout << "load " << ret << " items, " << idset.size() << " distinct items." << endl;
}

void loadTimeStamp(vector<pair<uint32_t, double>>&vec, const char*filename, int totalflow)
{
	set<uint32_t> idset;
	double time;
	cout << "open " << filename << endl;
	FILE* pf = fopen(filename, "rb");
	if (!pf) {
		cout << filename << "not found." << endl;
		system("pause");
		exit(-1);
	}
	vec.clear();
	char read_key[21];
	uint32_t ret = 0;
	double prevtime = 0;
	while (fread(read_key, 1, 21, pf))
	{
		time = *(double*)(read_key + 13);
		//cout << fixed << setprecision(8) << time << endl;
		if (prevtime == 0)
			prevtime = time;
		else{
			if (time < 0){
				cout << "CAIDA timestamp error, time < 0" << endl;
				system("pause");
			}
			if (time < prevtime){
				cout << "CAIDA timestamp error, time < prevtime" << endl;
				system("pause");
			}
			prevtime = time;
		}

		uint32_t key = *(uint32_t*)read_key;
		idset.insert(key);
		vec.push_back(pair<uint32_t, double>(key, time));
		ret++;
		if (ret == totalflow)
			break;
	}
	fclose(pf);
	cout << "load " << ret << " items, " << idset.size() << " distinct items." << endl;
}

void loadWebpage(vector<pair<uint32_t, double>>&vec, const char*filename, int totalflow)
{
	set<uint32_t> idset;
	double time;
	cout << "open " << filename << endl;
	FILE* pf = fopen(filename, "rb");
	if (!pf) {
		cout << filename << "not found." << endl;
		system("pause");
		exit(-1);
	}
	vec.clear();
	char read_key[13];
	
	uint32_t ret = 0;
	while (fread(read_key, 1, 13, pf))
	{
		uint32_t key = *(uint32_t*)read_key;
		time = *(double*)(read_key + 4);
		idset.insert(key);
		vec.push_back(pair<uint32_t, double>(key, time));
		ret++;
		if (ret == totalflow)
			break;
	}
	fclose(pf);
	cout << "load " << ret << " items, " << idset.size() << " distinct items." << endl;
}

void getUnique(vector<pair<uint32_t, double>>&unique, vector<pair<uint32_t, double>>vec){
	set<uint32_t> tmp;
	for (int i = 0; i < vec.size(); i++){
		if (tmp.count(vec[i].first) == 0){
			tmp.insert(vec[i].first);
			unique.push_back(vec[i]);
		}
	}
}