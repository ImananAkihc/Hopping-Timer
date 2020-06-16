#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include "BOBHash32.h"
using namespace std;

class Swamp
{
public:
	BOBHash32* hash;
	int hashnum;
	int interval;
	int tableSize;
	int* FreqTable;
	uint32_t** Window;
	int curr;
	int s;
	int distinctItems;
	Swamp(int _hashnum, int _interval, int _size, BOBHash32* _hash)
	{
		
		hashnum = _hashnum;
		interval = 10 * _interval; //每秒有10个元素插入
		hash = _hash;

		Window = new uint32_t*[hashnum];
		for (int i = 0; i < hashnum; i++)
			Window[i] = new uint32_t[interval];

		tableSize = _size - int(interval * hashnum * 1.5); //设counter都是8位的，再设window里面每个指纹是12位，window占去的内存要从table大小中减去
		if (tableSize < 0){
			cout << "size error" << endl;
			system("pause");
		}

		distinctItems = 0;
		s = 0;
		FreqTable = new int[tableSize];

		for (int i = 0; i < hashnum; i++)
			memset(Window[i], 0, interval * sizeof(uint32_t));
		
		memset(FreqTable, 0, tableSize * sizeof(uint32_t));
		cout << "create swamp, tableSize " << tableSize << ", hashnum " << hashnum << ", memory " << _size << endl;
	}

	~Swamp(){
		for (int i = 0; i < hashnum; i++)
			delete Window[i];
		delete Window;
		delete FreqTable;
	}

	void updateD(int prev, int cur){
		if (prev == 1)
			distinctItems--;
		if (cur == 1)
			distinctItems++;
	}

	void insert(uint32_t id){
		if (s < interval)
		{
			curr = s;
			int currmin = interval;
			for (int i = 0; i < hashnum; i++){
				uint32_t h = (hash[i].run((char*)&id, sizeof(uint32_t))) % tableSize;
				Window[i][s] = h;
				FreqTable[Window[i][curr]]++;
				if (FreqTable[Window[i][curr]] < currmin)
					currmin = FreqTable[Window[i][curr]];
			}
			updateD(0, currmin);
			s++;
		}
		else{
			int prevmin = interval;
			int currmin = interval;
			for (int i = 0; i < hashnum; i++){
				if (FreqTable[Window[i][curr]] < prevmin)
					prevmin = FreqTable[Window[i][curr]];
				FreqTable[Window[i][curr]]--;
				uint32_t h = (hash[i].run((char*)&id, sizeof(uint32_t))) % tableSize;
				Window[i][curr] = h;
				FreqTable[Window[i][curr]]++;
				if (FreqTable[Window[i][curr]] < currmin)
					currmin = FreqTable[Window[i][curr]];
			}
			updateD(prevmin, currmin);
			curr++;
			if (curr == interval)
				curr = 0;
		}
	}

	bool query(uint32_t id)
	{
		for (int i = 0; i < hashnum; i++)
		{
			uint32_t h = (hash[i].run((char*)&id, sizeof(uint32_t))) % tableSize;
			if (FreqTable[h] == 0)
				return 0;
		}
		return 1;

	}

	int getDistinctItems(){
		return distinctItems;
	}
};