#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
using namespace std;

class TBF
{
public:
	int* group;
	BOBHash32* hash;
	int hashnum;
	double t;
	int size;
	TBF(int _hashnum, double _t, int _size, BOBHash32* _hash)
	{
		hashnum = _hashnum;
		t = _t * 10 + 1;
		size = _size / 2;
		hash = _hash;
		group = new int[size];
		for (int i = 0; i < size; i++){
			group[i] = -1;
		}
		cout << "create TBF, tableSize " << size << ", hashnum " << hashnum << ", memory " << _size << endl;
	}

	~TBF(){
		delete group;
	}

	void insert(uint32_t id, uint64_t time)
	{
		int label = int(time * 10) % int(t);
		for (int i = 0; i<hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			group[s] = label;
		}
		for (int i = 0; i < size; i++){
			if (group[i] == (int(label + 1) % int(t)))
				group[i] = -1;
		}
	}

	bool query(uint32_t id, uint64_t time)
	{
		int label = int(time * 10) % int(t);
		for (int i = 0; i < hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			if (group[s] == -1 || group[s] == (int(label + 1) % int(t)))
				return 0;
		}
		return 1;
	}

};