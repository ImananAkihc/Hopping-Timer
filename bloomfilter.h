#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
using namespace std;


class BloomFilter
{
public:
	int* table;
	BOBHash32* hash;
	int hashnum;
	int size;
	BloomFilter(int _hashnum, int _size, BOBHash32* _hash)
	{
		hashnum = _hashnum;
		size = _size * 8;
		hash = _hash;
		table = new int[size];
		memset(table, 0, size * sizeof(int));
		cout << "create standard bloom filter, tableSize " << size << ", hashnum " << hashnum << ", memory " << size/8 << endl;
	}

	~BloomFilter(){
		delete table;
	}

	void insert(uint32_t id)
	{
		for (int i = 0; i < hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			table[s] = 1;
		}
	}

	bool query(uint32_t id)
	{
		for (int i = 0; i < hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			if (!table[s])
				return 0;
		}
		return 1;
	}

	void clear(){
		memset(table, 0, size * sizeof(int));
		return;
	}
};