#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
using namespace std;

class TSV
{
public:
	int* group;
	BOBHash32* hash;
	int t;
	int size;
	TSV(double _t, int _size, BOBHash32* _hash)
	{
		t = int(_t * 10);
		size = _size / 8;
		hash = _hash;
		group = new int[size];
		for (int i = 0; i < size; i++){
			group[i] = -1;
		}
		cout << "create TSV, tableSize " << size << ", memory " << _size << endl;
	}

	~TSV(){
		delete group;
	}

	void insert(uint32_t id, int time)
	{
		uint32_t s = hash[0].run((char*)&id, sizeof(uint32_t)) % size;
		group[s] = time;
	}

	bool query(uint32_t id, int time)
	{
		//cout << "query time " << time << endl;
		uint32_t s = hash[0].run((char*)&id, sizeof(uint32_t)) % size;
		if (time - group[s] < t && group[s] != -1)
			return 1;
		return 0;
	}

	int countflow(int time)
	{
		int empty = 0;
		for (int i = 0; i<size; i++)
		{
			if (!query(i, time))
				empty++;
		}
		return int(-log(double(empty) / double(size))*double(size));
	}
};