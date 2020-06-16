#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
#include <ctime>
using namespace std;

class CVS
{
public:
	int* group;
	BOBHash32* hash;
	int k;
	int C;
	int z;
	int size;
	CVS(double _t, int _size, BOBHash32* _hash)
	{
		C = 15;
		size = _size * 2;
		k = C * size / 10000;
		z = size;
		hash = _hash;
		group = new int[size];
		for (int i = 0; i < size; i++){
			group[i] = 0;
		}
		srand((time(0)));
		cout << "create CVS, tableSize " << size << ", memory " << _size << ", k = " << k << " C = " << C << endl;
	}

	~CVS(){
		delete group;
	}

	void insert(uint32_t id)
	{
		uint32_t s = hash[0].run((char*)&id, sizeof(uint32_t)) % size;
		if (group[s] == 0)
			z--;
		group[s] = C;
		for (int i = 0; i < k; i++){
			uint32_t r1 = rand();
			uint32_t r2 = rand();
			uint32_t r3 = rand();
			uint32_t t = (r1 * r2 + r3) % size;
			if (group[t] == 1)
				z++;
			if (group[t] >= 1)
				group[t]--;
		}
	}

	bool query(uint32_t id, int time)
	{
		uint32_t s = hash[0].run((char*)&id, sizeof(uint32_t)) % size;
		if (group[s] != 0)
			return 1;
		return 0;
	}

	int countflow()
	{
		//cout << "z: " << z << endl;
		return int(-log(double(z) / double(size))*double(size));
	}
};