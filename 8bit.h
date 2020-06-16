#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
using namespace std;

#define GROUP 0
#define HLLBUCKETS 128
#define speed 0

int countFirst1(int n){
	if (n == 0)
		return 10;
	if (n % 2 == 1)
		return 1;
	return countFirst1(n / 2) + 1;
}

class T_8bit
{
public:
	int* timestamp;
	int* counters;
	BOBHash32* hash;
	int hashnum;
	double t;
	int size;
	int NCC;
	int bits;
	int cells;
	int HALF;
	T_8bit(int BITS, int _hashnum, double _t, int _size, int _NCC, BOBHash32* _hash)
	{
		bits = BITS;
		cells = int(pow(2, BITS)) - 1;
		HALF = ((cells + 1) / 2);
		hashnum = _hashnum;
		t = _t / double(HALF);
		size = _size;
		NCC = _NCC;
		hash = _hash;
		timestamp = new int[size];
		counters = new int[cells];

		for (int i = 0; i < size; i++)
			timestamp[i] = -1;
		memset(counters, 0, cells * sizeof(int));
		cout << "create " << bits << "bits, tableSize " << size << ", hashnum " << hashnum << ", memory " << size << ", NCC " << NCC << endl;
	}

	~T_8bit(){
		delete timestamp;
		delete counters;
	}

	inline bool CheckSameWindow(int prev, int curr){
		if ((curr - prev) >= HALF ||
			(prev - curr < HALF && prev > curr))
			return false;
		if (prev == -1 || curr == -1){
			return false;
		}
		return true;
	}

	

	bool CheckAndClear(int index, double time, bool ifClear){
		int label = int(time / t) % cells;

		if (!ifClear)
			return CheckSameWindow(timestamp[index], label);
		int head = (index / NCC) * NCC;
		for (int x = 0; x < NCC; x++){
			uint32_t y = (head + x) % size;
			if (timestamp[y] != -1 && !CheckSameWindow(timestamp[y], label)){
				timestamp[y] = -1;
			}
		}

		return false;
	}
	//check the whether the timestamp[index] is avaliable and clear the bits nearby

	void insert(uint32_t id, double time)
	{
		int label = unsigned int(time / t) % cells;
		int min = cells;
		for (int i = 0; i<hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			if (!speed){
				if (min == cells)
					min = timestamp[s];
				else if (!CheckSameWindow(timestamp[s], label))
					min = -1;
				else{
					int oldgap = (label + cells - min) % cells;
					int newgap = (label + cells - timestamp[s]) % cells;
					if (newgap > oldgap && min != -1)
						min = timestamp[s];
				}
			}
			timestamp[s] = label;
			CheckAndClear(s, time, 1);
		}
		//cout << hashnum << ": min " << min << " label " << label << endl;
		if (!speed){
			int begin, end;
			begin = label;
			end = (label + HALF) % cells;
			if (CheckSameWindow(min, label)){
				begin = (min + HALF) % cells;
			}
			int cnt = 0;
			//cout << "begin: " << begin << " end: " << end << endl;
			while (begin != end){
				counters[begin]++;
				begin = (begin + 1) % cells;
				cnt++;
			}

			//cout << "+" << cnt << endl;
			for (int i = 1; i < HALF; i++){
				counters[((label + cells) - i) % cells] = 0;
			}
		}
	}

	int countflow(double time)
	{
		int empty = 0;
		for (int i = 0; i < size; i++){
			if (!CheckAndClear(i, time, 0))
				empty++;
		}
		return int(-log(double(empty) / size) * size);
	}

	bool query(uint32_t id, double time)
	{
		for (int i = 0; i<hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			if (!CheckAndClear(s,time, 0))
				return 0;
		}
		return 1;
	}

	int getDistinctItems(double time){
		int label = int(time / t) % cells;
		//cout << "label: " << label << endl;
		return counters[label];
	}
	

	int getDistinctItems_HLL(double time){
		double sum = 0;
		double KMax[HLLBUCKETS] = { 0 };
		for (int i = 0; i < size; i++){
			if (!CheckAndClear(i, time, 0))
				continue;
			if (KMax[i%HLLBUCKETS] < countFirst1(i / HLLBUCKETS))
				KMax[i%HLLBUCKETS] = countFirst1(i / HLLBUCKETS);
		}
		for (int i = 0; i < HLLBUCKETS; i++){
			if (KMax[i] == 0)
				KMax[i] = 0.00001;
			sum += 1.0 / KMax[i];
		}
		sum = HLLBUCKETS / sum;
		double constant = 1;
		/*switch ((int)log(HLLBUCKETS)) {
		case 4:
		constant = 0.673 * HLLBUCKETS * HLLBUCKETS;
		case 5:
		constant = 0.697 * HLLBUCKETS * HLLBUCKETS;
		case 6:
		constant = 0.709 * HLLBUCKETS * HLLBUCKETS;
		default:
		constant = (0.7213 / (1 + 1.079 / HLLBUCKETS)) * HLLBUCKETS * HLLBUCKETS;
		}*/
		return (int)(pow(2, sum) * HLLBUCKETS * constant / hashnum);
	}
	// hyperloglog, not used at present

	double setBitRate(){
		int n = 0;
		for (int i = 0; i < size; i++)
			if (timestamp[i] != -1)
				n++;
		return n * 1.0 / size;
	}
};