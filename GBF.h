#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include "BOBHash32.h"
using namespace std;

class GBF
{
public:
	int** group;
	BOBHash32* hash;
	int hashnum;
	double t;
	int size;
	int curr;
	int bits;
	int cells;
	GBF(int BITS, int _hashnum, double _t, int _size, BOBHash32* _hash)
	{
		bits = BITS;
		cells = int(pow(2, bits));
		hashnum = _hashnum;
		t = _t / double(cells /2);
		size = _size * 16 / cells;
		hash = _hash;
		curr = 0;
		group = new int * [size];
		for (int i = 0; i < size; i++){
			group[i] = new int[cells / 32];
			memset(group[i], 0, (cells / 32) * sizeof(int));
		}
		cout << "create GBF, tableSize " << size << ", hashnum " << hashnum << ", memory " << _size << ", sub-window " << cells << endl;
	}

	~GBF(){
		for (int i = 0; i < size; i++){
			delete group[i];
		}
		delete group;
	}
	void set(int row, int label){
		//cout << "set " << row << " " << label << endl;
		group[row][label / 32] |= (1 << (label % 32));
		return;
	}

	void clear(int row, int label){
		//cout << "clear " << row << " " << label << endl;
		group[row][label / 32] &= ~(1 << (label % 32));
		return;
	}

	void insert(uint32_t id, uint64_t time)
	{
		int label = int(time / t) % cells;
		if (label != curr){
			//cout << "new label " << label << endl;
			curr = label;
			int TBC = (curr + 1) % cells;
			for (int i = 0; i < size; i++)
				clear(i, (TBC + cells / 2 - 1) % cells);
		}
		for (int i = 0; i<hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			set(s, curr);
		}
	}

	bool query(uint32_t id, int show)
	{
		int* ans = new int[cells / 32];
		for (int i = 0; i < cells / 32; i++)
			ans[i] = -1;
		for (int i = 0; i < hashnum; i++)
		{
			uint32_t s = hash[i].run((char*)&id, sizeof(uint32_t)) % size;
			for (int j = 0; j < cells / 32; j++)
				ans[j] = ans[j] & group[s][j];
			
		}
		for (int i = 0; i < cells / 32; i++)
			if (ans[i])
				return 1;
		return 0;
	}

};