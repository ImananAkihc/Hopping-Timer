#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdint.h>
#include <vector>
#include <ctime>
#include <set>
#include "math.h"
#include "time.h"
#include <algorithm>
#include <fstream>
#include "BOBHash32.h"
#include "trace.h"
#include "8bit.h"
#include "swamp.h"
#include "bloomfilter.h"
#include "GBF.h"
#include "TBF.h"
#include "TSV.h"
#include "CVS.h"

using namespace std;

#define QueryNum 100000
#define HashNum 8
#define BitLength 8
#define Size 131072 //Size of bf = 128KB
#define Interval 6553.6 // propose that every second, 10 items arrived. 10000 elements in the window
#define TimeInterval 0.1 // Used for time based
#define DataSet 0 // 0=CAIDA 1=DISTINCT 2=TIME 3=WEBPAGE 4=CAMPUS


BOBHash32 hashes[1000];
vector<uint32_t> QuerySet;

ofstream fp;

void CreateQuerySet(int querynum){
	for (int i = 0; i < querynum; i++){
		uint32_t a1 = rand();
		uint32_t a2 = rand() % 2;
		uint32_t b1 = rand();
		uint32_t b2 = rand() % 2;
		uint32_t com = ((a1 | (a2 << 15)) << 16) | (b1 | (b2 << 15));
		QuerySet.push_back(com);
	}
}

void InitBOBHash(int hashnum){
	set<int>seeds;
	for (int i = 0; i < hashnum; i++)
	{
		uint32_t seed = rand() % MAX_PRIME32;
		while (seeds.find(seed) != seeds.end())
			seed = rand() % MAX_PRIME32;
		seeds.insert(seed);
		hashes[i].initialize(seed);
	}
}

void InitDataSet(vector<pair<uint32_t, double>>&vec, int totalflow){
	if (DataSet == 1)
	{
		vector<pair<uint32_t, double>> flow;
		loadTrace(flow, "formatted00.dat", totalflow * 10);
		cout << "deleting ";
		getUnique(vec, flow);
		cout << "remain " << vec.size() << endl;
	}
	else if (DataSet == 2)
		loadTimeStamp(vec, "130000.dat", totalflow);
	else if (DataSet == 3)
		loadWebpage(vec, "webdocs_form00.dat", totalflow);
	else if (DataSet == 4)
		loadWebpage(vec, "campus.dat", totalflow);
	else
		loadTrace(vec, "formatted00.dat", totalflow);
}


void CompareFPRwithTime(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_32K(8, HashNum, Interval, Size / 4, 64, hashes);
	T_8bit t_48K(8, HashNum, Interval, Size * 3 / 8, 64, hashes);
	T_8bit t_64K(8, HashNum, Interval, Size / 2, 64, hashes);

	set<uint32_t> unique;
	double fal, neg, pos;
	double t_32KFPR = 0, t_64KFPR = 0, t_48KFPR = 0;

	int gap = 10 * Interval;
	if (DataSet == 1)
		gap /= 2;
	int r = 0;
	cout << "test every " << gap << " insertion" << endl;

	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		now = double(i) / 10;
		t_32K.insert(flow[i].first, now);
		t_64K.insert(flow[i].first, now);
		t_48K.insert(flow[i].first, now);
		if ((i + 1) % gap == 0){
			fal = neg = pos = 0;
			unique.clear();
			for (int j = i; j > max(0, i - int(Interval * 10)); j--)
				unique.insert(flow[j].first);

			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_32K.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_32KFPR += pos * 1.0 / fal;

			fal = neg = pos = 0;

			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_64K.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_64KFPR += pos * 1.0 / fal;


			fal = neg = pos = 0;
			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_48K.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_48KFPR += pos * 1.0 / fal;


			if (r % 4 == 0){
				fp << t_32KFPR / 4 << " ";
				cout << "the FPR of 32K: " << t_32KFPR / 4 << " ";

				fp << t_48KFPR / 4 << " ";
				cout << "the FPR of 48K: " << t_48KFPR / 4 << " ";
				
				fp << t_64KFPR / 4 << endl;
				cout << "the FPR of 64K: " << t_64KFPR / 4 << endl;

				t_32KFPR = t_64KFPR = t_48KFPR = 0;
			}
			r++;
		}

	}
}

void CompareDistinctwithTime(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_32K(8, HashNum, Interval, Size / 4, 64, hashes);
	T_8bit t_64K(8, HashNum, Interval, Size / 2, 64, hashes);
	T_8bit t_48K(8, HashNum, Interval, Size * 3/ 8, 64, hashes);
	set<uint32_t> unique;
	double t32KError = 0, t64KError = 0, t48KError = 0;

	int gap = 10 * Interval;
	if (DataSet == 1)
		gap /= 2;
	int r = 0;
	cout << "test every " << gap << " insertion" << endl;
	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		now = double(i) / 10;
		t_32K.insert(flow[i].first, now);
		t_64K.insert(flow[i].first, now);
		t_48K.insert(flow[i].first, now);
		if ((i + 1) % gap == 0){
			unique.clear();
			for (int j = i; j > max(0, i - int(Interval * 10)); j--)
				unique.insert(flow[j].first);

			cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;

			t32KError += abs(double(unique.size()) - double(t_32K.getDistinctItems(now))) / (unique.size() * 1.0);
			t64KError += abs(double(unique.size()) - double(t_64K.getDistinctItems(now))) / (unique.size() * 1.0);
			t48KError += abs(double(unique.size()) - double(t_48K.getDistinctItems(now))) / (unique.size() * 1.0);

			if (r % 4 == 0){
				fp << t32KError / 4 << " ";
				cout << "the RE of 32K: " << t32KError / 4 << " ";

				fp << t48KError / 4 << " ";
				cout << "the RE of 48K: " << t48KError / 4 << " ";

				fp << t64KError / 4 << endl;
				cout << "the RE of 64K: " << t64KError / 4 << endl;

				t32KError = t48KError = t64KError = 0;
			}
			r++;
		}

	}
}

void CompareBMDistinctwithTime(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_32K(8, 1, Interval, Size / 4, 64, hashes);
	T_8bit t_64K(8, 1, Interval, Size / 2, 64, hashes);
	T_8bit t_48K(8, 1, Interval, Size * 3 / 8, 64, hashes);
	set<uint32_t> unique;
	double t32KError = 0, t64KError = 0, t48KError = 0;

	int gap = 10 * Interval;
	if (DataSet == 1)
		gap /= 2;
	int r = 0;
	cout << "test every " << gap << " insertion" << endl;
	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		now = double(i) / 10;
		t_32K.insert(flow[i].first, now);
		t_64K.insert(flow[i].first, now);
		t_48K.insert(flow[i].first, now);
		if ((i + 1) % gap == 0){
			unique.clear();
			for (int j = i; j > max(0, i - int(Interval * 10)); j--)
				unique.insert(flow[j].first);

			cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;

			t32KError += abs(double(unique.size()) - double(t_32K.countflow(now))) / (unique.size() * 1.0);
			t64KError += abs(double(unique.size()) - double(t_64K.countflow(now))) / (unique.size() * 1.0);
			t48KError += abs(double(unique.size()) - double(t_48K.countflow(now))) / (unique.size() * 1.0);

			if (r % 4 == 0){
				fp << t32KError / 4 << " ";
				cout << "the RE of 32K: " << t32KError / 4 << " ";

				fp << t48KError / 4 << " ";
				cout << "the RE of 48K: " << t48KError / 4 << " ";

				fp << t64KError / 4 << endl;
				cout << "the RE of 64K: " << t64KError / 4 << endl;

				t32KError = t48KError = t64KError = 0;
			}
			r++;
		}

	}
}

void TimeBasedCompareFPRwithTime(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_64K(8, HashNum, TimeInterval, Size / 2, 4, hashes);
	T_8bit t_128K(8, HashNum, TimeInterval, Size, 4, hashes);
	T_8bit t_64KO(8, HashNum, TimeInterval, Size / 2, 4, hashes);
	T_8bit t_128KO(8, HashNum, TimeInterval, Size, 4, hashes);

	set<uint32_t> unique;
	double fal, neg, pos;
	double t_64KFPR = 0, t_128KFPR = 0, t_64KOFPR = 0, t_128KOFPR = 0;

	double gap = TimeInterval / 2.0;
	int r = 0;
	int flag = 0;
	cout << "test every " << gap << " second" << endl;
	double start = flow[0].second;
	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		if (flag == 1){
			for (int j = 0; j < t_64KO.size; j += t_64KO.NCC)
				t_64KO.CheckAndClear(j, now, 1);
			for (int j = 0; j < t_128KO.size; j += t_128KO.NCC)
				t_128KO.CheckAndClear(j, now, 1);
			flag = 0;
		}
		now = flow[i].second - start;
		t_64K.insert(flow[i].first, now);
		t_128K.insert(flow[i].first, now);
		t_64KO.insert(flow[i].first, now);
		t_128KO.insert(flow[i].first, now);
		if (i == flow.size() - 1 || int(now / gap) != int((flow[i + 1].second - start) / gap)){
			cout << "test: " << r << " " << now << endl;
			flag = 1;
			unique.clear();
			for (int j = i; j > 0; j--){
				if (now - (flow[j].second - start) > TimeInterval)
					break;
				unique.insert(flow[j].first);
			}
			cout << unique.size() << " items in 0.05s" << endl;

			fal = neg = pos = 0;
			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_64K.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_64KFPR += pos * 1.0 / fal;
			cout << "64: " << pos * 1.0 / fal << endl;

			fal = neg = pos = 0;
			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_128K.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_128KFPR += pos * 1.0 / fal;
			cout << "128: " << pos * 1.0 / fal << endl;

			fal = neg = pos = 0;
			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_64KO.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_64KOFPR += pos * 1.0 / fal;
			cout << "64O: " << pos * 1.0 / fal << endl;

			fal = neg = pos = 0;
			for (int j = 0; j < querySize; j++)
				if (unique.count(QuerySet[j]) == 0)
				{
					fal++;
					if (t_128KO.query(QuerySet[j], now) == 1)
						pos++;
				}

			t_128KOFPR += pos * 1.0 / fal;
			cout << "128O: " << pos * 1.0 / fal << endl;


			fp << t_64KFPR << " ";
			cout << "the FPR of 64K: " << t_64KFPR << " ";

			fp << t_128KFPR << " ";
			cout << "the FPR of 128K: " << t_128KFPR << " ";

			fp << t_64KOFPR << " ";
			cout << "the FPR of 64KO: " << t_64KOFPR << " ";

			fp << t_128KOFPR << endl;
			cout << "the FPR of 128KO: " << t_128KOFPR << endl;

			t_64KFPR = t_128KFPR = t_64KOFPR = t_128KOFPR = 0;
		}

	}
}

void TimeBasedCompareDistinctwithTime(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_32K(8, HashNum, TimeInterval, Size / 4, 4, hashes);
	T_8bit t_64K(8, HashNum, TimeInterval, Size / 2, 4, hashes);
	T_8bit t_32KO(8, HashNum, TimeInterval, Size / 4, 4, hashes);
	T_8bit t_64KO(8, HashNum, TimeInterval, Size / 2, 4, hashes);
	set<uint32_t> unique;
	double t32KError = 0, t64KError = 0, t32KOError = 0, t64KOError = 0;

	double gap = TimeInterval / 2;
	int r = 0;
	int flag = 0;
	cout << "test every " << gap << " second" << endl;
	double start = flow[0].second;
	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		if (flag == 1){
			for (int j = 0; j < t_64KO.size; j += t_64KO.NCC)
				t_64KO.CheckAndClear(j, now, 1);
			for (int j = 0; j < t_32KO.size; j += t_32KO.NCC)
				t_32KO.CheckAndClear(j, now, 1);
			flag = 0;
		}
		now = flow[i].second - start;
		t_32K.insert(flow[i].first, now);
		t_64K.insert(flow[i].first, now);
		t_32KO.insert(flow[i].first, now);
		t_64KO.insert(flow[i].first, now);
		if (i == flow.size() - 1 || int(now / gap) != int((flow[i + 1].second - start) / gap)){
			flag = 1;
			unique.clear();
			for (int j = i; j > 0; j--){
				if (now - (flow[j].second - start) > TimeInterval)
					break;
				unique.insert(flow[j].first);
			}
			cout << unique.size() << " distinct items in recent " << 0.05 << "s" << endl;

			t32KError += abs(double(unique.size()) - double(t_32K.getDistinctItems(now))) / (unique.size() * 1.0);
			t64KError += abs(double(unique.size()) - double(t_64K.getDistinctItems(now))) / (unique.size() * 1.0);
			t32KOError += abs(double(unique.size()) - double(t_32KO.getDistinctItems(now))) / (unique.size() * 1.0);
			t64KOError += abs(double(unique.size()) - double(t_64KO.getDistinctItems(now))) / (unique.size() * 1.0);

			fp << t32KError << " ";
			cout << "the RE of 32K: " << t32KError << " ";

			fp << t64KError << " ";
			cout << "the RE of 64K: " << t64KError << endl;

			fp << t32KOError << " ";
			cout << "the RE of 32KO: " << t32KOError << " ";

			fp << t64KOError << endl;
			cout << "the RE of 64KO: " << t64KOError << endl;

			t32KError = t64KError = t32KOError = t64KOError = 0;
		}

	}
}

void TimeBasedCompareBMDistinctwithTime(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	T_8bit t_4K(8, 1, TimeInterval, Size / 32, 4, hashes);
	T_8bit t_8K(8, 1, TimeInterval, Size / 16, 4, hashes);
	T_8bit t_4KO(8, 1, TimeInterval, Size / 32, 4, hashes);
	T_8bit t_8KO(8, 1, TimeInterval, Size / 16, 4, hashes);
	set<uint32_t> unique;
	double t4KError = 0, t8KError = 0, t4KOError = 0, t8KOError = 0;

	double gap = TimeInterval / 2;
	int r = 0;
	int flag = 0;
	cout << "test every " << gap << " second" << endl;
	double start = flow[0].second;
	for (int i = 0; i < flow.size(); i++)
	{
		if (i % 10000 == 0)
			cout << i << endl;
		if (flag == 1){
			for (int j = 0; j < t_4KO.size; j += t_4KO.NCC)
				t_4KO.CheckAndClear(j, now, 1);
			for (int j = 0; j < t_8KO.size; j += t_8KO.NCC)
				t_8KO.CheckAndClear(j, now, 1);
			flag = 0;
		}
		now = flow[i].second - start;
		t_4K.insert(flow[i].first, now);
		t_8K.insert(flow[i].first, now);
		t_4KO.insert(flow[i].first, now);
		t_8KO.insert(flow[i].first, now);
		if (i == flow.size() - 1 || int(now / gap) != int((flow[i + 1].second - start) / gap)){
			unique.clear();
			for (int j = i; j > 0; j--){
				if (now - (flow[j].second - start) > TimeInterval)
					break;
				unique.insert(flow[j].first);
			}
			flag = 1;
			cout << unique.size() << " distinct items in recent " << 0.05 << "s" << endl;

			t4KError += abs(double(unique.size()) - double(t_4K.countflow(now))) / (unique.size() * 1.0);
			t8KError += abs(double(unique.size()) - double(t_8K.countflow(now))) / (unique.size() * 1.0);
			t4KOError += abs(double(unique.size()) - double(t_4KO.countflow(now))) / (unique.size() * 1.0);
			t8KOError += abs(double(unique.size()) - double(t_8KO.countflow(now))) / (unique.size() * 1.0);

			fp << t4KError << " ";
			cout << "the RE of 4K: " << t4KError << " ";

			fp << t8KError << " ";
			cout << "the RE of 8K: " << t8KError << endl;

			fp << t4KOError << " ";
			cout << "the RE of 4KO: " << t4KOError << " ";

			fp << t8KOError << endl;
			cout << "the RE of 8KO: " << t8KOError << endl;
			t4KError = t8KError = t4KOError = t8KOError = 0;
		}

	}
}

void CompareFPRwithDifferentMemory(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int interv = Interval; interv <= (2 * Interval); interv += (Interval / 4)){
		T_8bit t_16K(8, HashNum, interv, Size / 8, 64, hashes);
		T_8bit t_32K(8, HashNum, interv, Size / 4, 64, hashes);
		T_8bit t_64K(8, HashNum, interv, Size / 2, 64, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double t_32KFPR = 0, t_64KFPR = 0, t_16KFPR = 0;
		int r = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_32K.insert(flow[i].first, now);
			t_64K.insert(flow[i].first, now);
			t_16K.insert(flow[i].first, now);

			if ((i + 1) % gap == 0){
				r++;
				fal = neg = pos = 0;
				unique.clear();
				for (int j = i; j > max(0, i - interv * 10); j--)
					unique.insert(flow[j].first);


				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_32K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_32KFPR += pos * 1.0 / fal;

				fal = neg = pos = 0;

				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_64K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_64KFPR += pos * 1.0 / fal;

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_16K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_16KFPR += pos * 1.0 / fal;

			}
		}
		cout << "average of " << r << " times" << endl;
		fp << t_16KFPR / r << " ";
		cout << "the FPR of 16K: " << t_16KFPR / r << " " << endl;

		fp << t_32KFPR / r << " ";
		cout << "the FPR of 32K: " << t_32KFPR / r << " " << endl;

		fp << t_64KFPR / r << endl;
		cout << "the FPR of 64K: " << t_64KFPR / r << " " << endl;

		
	}
}

void CompareDistinctElements(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int interv = Interval; interv <= (2 * Interval); interv += (Interval / 4)){
		cout << "interv: " << interv << endl;
		T_8bit t_32K(8, HashNum, interv, Size / 4, 64, hashes);
		T_8bit t_64K(8, HashNum, interv, Size / 2, 64, hashes);
		T_8bit t_48K(8, HashNum, interv, Size * 3 / 8, 64, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t32KError = 0, t64KError = 0, t48KError = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_32K.insert(flow[i].first, now);
			t_64K.insert(flow[i].first, now);
			t_48K.insert(flow[i].first, now);

			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - interv * 10); j--){
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " distinct items in recent " << interv * 10 << " items" << endl;

				t32KError += abs(double(unique.size()) - double(t_32K.getDistinctItems(now))) / (unique.size() * 1.0);
				t64KError += abs(double(unique.size()) - double(t_64K.getDistinctItems(now))) / (unique.size() * 1.0);
				t48KError += abs(double(unique.size()) - double(t_48K.getDistinctItems(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t32KError / r << " ";
		cout << "the RE of 32K: " << t32KError / r << " " << endl;

		fp << t48KError / r << " ";
		cout << "the RE of 48K: " << t48KError / r << " " << endl;

		fp << t64KError / r << endl;
		cout << "the RE of 64K: " << t64KError / r << " " << endl;
	}
}

void CompareBMDistinctElements(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int interv = Interval; interv <= (2 * Interval); interv += (Interval / 4)){
		cout << "interv: " << interv << endl;
		T_8bit t_32K(8, 1, interv, Size / 32, 64, hashes);
		T_8bit t_64K(8, 1, interv, Size / 8, 64, hashes);
		T_8bit t_48K(8, 1, interv, Size / 16, 64, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t32KError = 0, t64KError = 0, t48KError = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_32K.insert(flow[i].first, now);
			t_64K.insert(flow[i].first, now);
			t_48K.insert(flow[i].first, now);

			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - interv * 10); j--){
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " distinct items in recent " << interv * 10 << " items" << endl;
				cout << "t32K: " << t_32K.countflow(now) << endl;
				t32KError += abs(double(unique.size()) - double(t_32K.countflow(now))) / (unique.size() * 1.0);
				cout << "t64K: " << t_64K.countflow(now) << endl;
				t64KError += abs(double(unique.size()) - double(t_64K.countflow(now))) / (unique.size() * 1.0);
				cout << "t48K: " << t_48K.countflow(now) << endl;
				t48KError += abs(double(unique.size()) - double(t_48K.countflow(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t32KError / r << " ";
		cout << "the RE of 32K: " << t32KError / r << " " << endl;

		fp << t48KError / r << " ";
		cout << "the RE of 48K: " << t48KError / r << " " << endl;

		fp << t64KError / r << endl;
		cout << "the RE of 64K: " << t64KError / r << " " << endl;
	}
}

void TimeBasedCompareFPRwithDifferentMemory(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);
	
	for (double interv = TimeInterval / 4; interv <= (4 * TimeInterval); interv *= 2){
		cout << "interv " << interv << endl;
		T_8bit t_64K(8, HashNum, interv, Size / 2, 4, hashes);
		T_8bit t_128K(8, HashNum, interv, Size, 4, hashes);
		T_8bit t_64KO(8, HashNum, interv, Size / 2, 4, hashes);
		T_8bit t_128KO(8, HashNum, interv, Size, 4, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double t_64KFPR = 0, t_128KFPR = 0, t_64KOFPR = 0, t_128KOFPR = 0;
		double gap = TimeInterval * 2.0;
		int r = 0;
		int flag = 0;
		cout << "test every " << gap << " second" << endl;
		double start = flow[0].second;

		for (int i = 0; i < flow.size(); i++)
		{
			now = flow[i].second - start;
			if (flag == 1){
				//cout << "clear at " << now << endl;
				//cout << int(now / interv) << " " << int((flow[i - 1].second - start) / interv) << endl;
				for (int j = 0; j < t_64KO.size; j += t_64KO.NCC)
					t_64KO.CheckAndClear(j, now, 1);
				for (int j = 0; j < t_128KO.size; j += t_128KO.NCC)
					t_128KO.CheckAndClear(j, now, 1);
				flag = 0;
			}
			if (i % 10000 == 0){
				cout << i << endl;
				//cout << now << endl;
			}
			t_64K.insert(flow[i].first, now);
			t_128K.insert(flow[i].first, now);
			t_64KO.insert(flow[i].first, now);
			t_128KO.insert(flow[i].first, now);

			if (i == flow.size() - 1 || int(now / gap) != int((flow[i + 1].second - start) / gap)){
				r++;
				flag = 1;
				unique.clear();
				for (int j = i; j > 0; j--){
					if (now - (flow[j].second - start) > interv)
					{	
						cout << i - j << " different items" << endl;
						break;
					}
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " items in the window" << endl;


				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_64K.query(QuerySet[j], now) == 1)
							pos++;
					}
				t_64KFPR += pos * 1.0 / fal;
				cout << "64: " << pos * 1.0 / fal << " ";

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_128K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_128KFPR += pos * 1.0 / fal;
				cout << "128: " << pos * 1.0 / fal << endl;

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_64KO.query(QuerySet[j], now) == 1)
							pos++;
					}
				t_64KOFPR += pos * 1.0 / fal;
				cout << "64: " << pos * 1.0 / fal << " ";

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_128KO.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_128KOFPR += pos * 1.0 / fal;
				cout << "128: " << pos * 1.0 / fal << endl;
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t_64KFPR / r << " ";
		cout << "the FPR of 64K: " << t_64KFPR / r << " ";

		fp << t_128KFPR / r << " ";
		cout << "the FPR of 128K: " << t_128KFPR / r << " " << endl;

		fp << t_64KOFPR / r << " ";
		cout << "the FPR of 64KO: " << t_64KOFPR / r << " ";

		fp << t_128KOFPR / r << endl;
		cout << "the FPR of 128KO: " << t_128KOFPR / r << " " << endl;


	}
}

void TimeBasedCompareDistinctElements(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (double interv = TimeInterval / 4; interv <= (4 * TimeInterval); interv *= 2){
		cout << "interv: " << interv << endl;
		T_8bit t_32K(8, HashNum, interv, Size / 4, 4, hashes);
		T_8bit t_64K(8, HashNum, interv, Size / 2, 4, hashes);
		T_8bit t_32KO(8, HashNum, interv, Size / 4, 4, hashes);
		T_8bit t_64KO(8, HashNum, interv, Size / 2, 4, hashes);
		set<uint32_t> unique;
		int r = 0;
		int flag = 0;
		double t32KError = 0, t64KError = 0, t32KOError = 0, t64KOError = 0;
		double gap = interv / 2.0;
		double start = flow[0].second;

		for (int i = 0; i < flow.size(); i++)
		{
			if (flag == 1){
				for (int j = 0; j < t_32KO.size; j += t_32KO.NCC)
					t_32KO.CheckAndClear(j, now, 1);
				for (int j = 0; j < t_64KO.size; j += t_64KO.NCC)
					t_64KO.CheckAndClear(j, now, 1);
				flag = 0;
			}
			now = flow[i].second - start;
			t_32K.insert(flow[i].first, now);
			t_64K.insert(flow[i].first, now);
			t_32KO.insert(flow[i].first, now);
			t_64KO.insert(flow[i].first, now);

			if (i == flow.size() - 1 || int(now / gap) != int((flow[i + 1].second - start) / gap)){
				r++;
				flag = 1;
				unique.clear();
				for (int j = i; j > 0; j--){
					if (now - (flow[j].second - start) > interv)
						break;
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " distinct items in the window" << endl;

				t32KError += abs(double(unique.size()) - double(t_32K.getDistinctItems(now))) / (unique.size() * 1.0);
				t64KError += abs(double(unique.size()) - double(t_64K.getDistinctItems(now))) / (unique.size() * 1.0);
				t32KOError += abs(double(unique.size()) - double(t_32KO.getDistinctItems(now))) / (unique.size() * 1.0);
				t64KOError += abs(double(unique.size()) - double(t_64KO.getDistinctItems(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t32KError / r << " ";
		cout << "the RE of 32K: " << t32KError / r << " " << endl;

		fp << t64KError / r << " ";
		cout << "the RE of 64K: " << t64KError / r << " " << endl;

		fp << t32KOError / r << " ";
		cout << "the RE of 32KO: " << t32KOError / r << " " << endl;

		fp << t64KOError / r << endl;
		cout << "the RE of 64KO: " << t64KOError / r << " " << endl;
	}
}

void TimeBasedCompareBMDistinctElements(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (double interv = TimeInterval / 4; interv <= (4 * TimeInterval); interv *= 2){
		cout << "interv: " << interv << endl;
		T_8bit t_4K(8, 1, interv, Size / 32, 4, hashes);
		T_8bit t_8K(8, 1, interv, Size / 16, 4, hashes);
		T_8bit t_4KO(8, 1, interv, Size / 32, 4, hashes);
		T_8bit t_8KO(8, 1, interv, Size / 16, 4, hashes);
		set<uint32_t> unique;
		int r = 0;
		int flag = 0;
		double t4KError = 0, t8KError = 0, t4KOError = 0, t8KOError = 0;
		int gap = flow.size() / 8;
		double start = flow[0].second;

		for (int i = 0; i < flow.size(); i++)
		{
			now = flow[i].second - start;
			if (flag == 1){
				for (int j = 0; j < t_4KO.size; j += t_4KO.NCC)
					t_4KO.CheckAndClear(j, now, 1);
				for (int j = 0; j < t_8KO.size; j += t_8KO.NCC)
					t_8KO.CheckAndClear(j, now, 1);
				flag = 0;
			}
			t_4K.insert(flow[i].first, now);
			t_8K.insert(flow[i].first, now);
			t_4KO.insert(flow[i].first, now);
			t_8KO.insert(flow[i].first, now);

			if ((i + 1) % gap == 0){
				r++;
				flag = 1;
				unique.clear();
				for (int j = i; j > 0; j--){
					if (now - (flow[j].second - start) > interv)
						break;
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " distinct items in the window" << endl;

				t4KError += abs(double(unique.size()) - double(t_4K.countflow(now))) / (unique.size() * 1.0);
				t8KError += abs(double(unique.size()) - double(t_8K.countflow(now))) / (unique.size() * 1.0);
				t4KOError += abs(double(unique.size()) - double(t_4KO.countflow(now))) / (unique.size() * 1.0);
				t8KOError += abs(double(unique.size()) - double(t_8KO.countflow(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t4KError / r << " ";
		cout << "the RE of 4K: " << t4KError / r << " " << endl;

		fp << t8KError / r << " ";
		cout << "the RE of 8K: " << t8KError / r << " " << endl;

		fp << t4KOError / r << " ";
		cout << "the RE of 4KO: " << t4KOError / r << " " << endl;

		fp << t8KOError / r << endl;
		cout << "the RE of 8KO: " << t8KOError / r << " " << endl;

	}
}

/*void CompareDistinctNoError(char* filename, int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	loadTrace(flow, filename, flowSize);

	for (int interv = Interval / 2; interv <= (3 * Interval); interv += (Interval / 2)){
		cout << "interv: " << interv << endl;
		T_8bit t_32K(8, HashNum, interv, Size/2, 100, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t32K = 0, real = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_32K.insert(flow[i].first, now);

			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0,i - interv * 10); j--){
					unique.insert(flow[j].first);
				}
				cout << unique.size() << " distinct items in recent " << interv * 10 << " items ";
				cout << "8bit answer " << t_32K.getDistinctItems(now) << endl;

				t32K += t_32K.getDistinctItems(now);
				real += unique.size();
			}
		}
		cout << "average of " << r << " times" << endl;

		t_8bitDistinct_NoError << t32K / r << " ";
		cout << "the count of 32K: " << t32K / r << " " ;

		t_8bitDistinct_NoError << real / r << endl;
		cout << "the count of real: " << real / r << endl;
	}
}*/

void CompareFPRWithOther(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int m = 0; m <= 30; m += 2){
		T_8bit t_8bit(8, HashNum, Interval, (Size * (2 + m)) / 16, 4, hashes);
		Swamp swamp(HashNum, Interval, (Size * (18 + m)) / 16, hashes);
		TBF tbf(HashNum, Interval, (Size * (2 + m)) / 16, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double EbitFPR = 0, swampFPR = 0, TBFFPR = 0;
		int r = 0;
		int gap = 10 * Interval;
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			if (i % 10000 == 0)
				cout << i << endl;
			now = double(i) / 10;
			t_8bit.insert(flow[i].first, now);
			swamp.insert(flow[i].first);
			tbf.insert(flow[i].first, now);
			if ((i+1) % gap == 0){
				r++;
				fal = neg = pos = 0;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--){
					unique.insert(flow[j].first);
				}
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_8bit.query(QuerySet[j], now) == 1)
							pos++;
					}

				EbitFPR += pos * 1.0 / fal;

				fal = neg = pos = 0;

				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (swamp.query(QuerySet[j]) == 1)
							pos++;
					}
				swampFPR += pos * 1.0 / fal;

				fal = neg = pos = 0;

				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (tbf.query(QuerySet[j],0) == 1)
							pos++;
					}

				TBFFPR += pos * 1.0 / fal;

			}
		}
		cout << "average of " << r << " times" << endl;

		fp << EbitFPR / r << " ";
		cout << "the FPR of 8bit: " << EbitFPR / r << " ";

		fp << swampFPR / r << " ";
		cout << "the FPR of swamp: " << swampFPR / r << " ";

		fp << TBFFPR / r << endl;
		cout << "the FPR of tbf: " << TBFFPR / r << endl << endl;
		
	}
}

void CompareDistinctElementsWithOther(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int m = 0; m <= 30; m += 2){
		T_8bit t_8bit(8, HashNum, Interval, (Size * (2 + m)) / 16 - 512, 4, hashes);
		T_8bit t_8bitBM(8, 1, Interval, (Size * (2 + m)) / 16, 64, hashes);
		Swamp swamp(HashNum, Interval, (Size * (18 + m)) / 16, hashes);
		TSV tsv(Interval, (Size * (2 + m)) / 16, hashes);
		CVS cvs(Interval, (Size * (2 + m)) / 16, hashes);
		int r = 0;
		double t_8bitError = 0, swamp_Error = 0, TSV_Error = 0, CVS_Error = 0, t_8bitBMError = 0;
		set<uint32_t> unique;
		int gap = flow.size() / 8;
		gap = (gap / int(Interval * 10 / pow(2, BitLength - 1))) * int(Interval * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_8bit.insert(flow[i].first, now);
			t_8bitBM.insert(flow[i].first, now);
			swamp.insert(flow[i].first);
			tsv.insert(flow[i].first, i);
			cvs.insert(flow[i].first);
			
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0,i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;

				t_8bitError += abs(double(unique.size()) - double(t_8bit.getDistinctItems(now))) / (unique.size() * 1.0);
				//cout << t_8bit.getDistinctItems(now) << endl;
				t_8bitBMError += abs(double(unique.size()) - double(t_8bitBM.countflow(now))) / (unique.size() * 1.0);
				//cout << t_8bitBM.countflow(now) << endl;
				swamp_Error += abs(double(unique.size()) - double(swamp.getDistinctItems())) / (unique.size() * 1.0);
				//cout << swamp.getDistinctItems() << endl;
				TSV_Error += abs(double(unique.size()) - double(tsv.countflow(i))) / (unique.size() * 1.0);
				//cout << tsv.countflow(i) << endl;
				CVS_Error += abs(double(unique.size()) - double(cvs.countflow())) / (unique.size() * 1.0);
				//cout << cvs.countflow() << endl;
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t_8bitError / r << " ";
		cout << "the RE of 8bit " << t_8bitError / r << " " << endl;

		fp << t_8bitBMError / r << " ";
		cout << "the RE of 8bitBM " << t_8bitBMError / r << " " << endl;

		fp << swamp_Error / r << " ";
		cout << "the RE of swamp: " << swamp_Error / r << " " << endl;

		fp << TSV_Error / r << " ";
		cout << "the RE of TSV: " << TSV_Error / r << " " << endl;

		fp << CVS_Error / r << endl;
		cout << "the RE of CVS: " << CVS_Error / r << " " << endl;

	}
}

void CompareBits(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);
	for (int size = Size / 8; size <= Size; size += Size / 8){
		T_8bit t_4bit(4, HashNum, Interval, size * 2, 4, hashes);
		T_8bit t_6bit(6, HashNum, Interval, size * 3 / 2, 4, hashes);
		T_8bit t_8bit(8, HashNum, Interval, size, 4, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double t_4FPR = 0, t_6FPR = 0, t_8FPR = 0;
		int r = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(Interval * 10 / pow(2, BitLength - 1))) * int(Interval * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_4bit.insert(flow[i].first, now);
			t_6bit.insert(flow[i].first, now);
			t_8bit.insert(flow[i].first, now);
			
			if ((i + 1) % gap == 0){
				r++;
				fal = neg = pos = 0;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);


				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_4bit.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_4FPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;

				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_6bit.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_6FPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_8bit.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_8FPR += (pos * 1.0 / fal);

			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t_4FPR / r << " ";
		cout << "the FPR of 4bit: " << t_4FPR / r << " " << endl;

		fp << t_6FPR / r << " ";
		cout << "the FPR of 6bit: " << t_6FPR / r << " " << endl;

		fp << t_8FPR / r << endl;
		cout << "the FPR of 8bit: " << t_8FPR / r << " " << endl;
	}
}

void CompareBitsDistinct(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int size = Size / 8; size <= Size; size += Size / 8){
		T_8bit t_4(4, HashNum, Interval, size * 2, 4, hashes);
		T_8bit t_6(6, HashNum, Interval, size * 3 / 2, 4, hashes);
		T_8bit t_8(8, HashNum, Interval, size, 4, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t4Error = 0, t6Error = 0, t8Error = 0;
		int gap = 12345;
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_4.insert(flow[i].first, now);
			t_6.insert(flow[i].first, now);
			t_8.insert(flow[i].first, now);
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;

				t4Error += abs(double(unique.size()) - double(t_4.getDistinctItems(now))) / (unique.size() * 1.0);
				t6Error += abs(double(unique.size()) - double(t_6.getDistinctItems(now))) / (unique.size() * 1.0);
				t8Error += abs(double(unique.size()) - double(t_8.getDistinctItems(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t4Error / r << " ";
		cout << "the RE of 4: " << t4Error / r << " " << endl;

		fp << t6Error / r << " ";
		cout << "the RE of 6: " << t6Error / r << " " << endl;

		fp << t8Error / r << endl;
		cout << "the RE of 8: " << t8Error / r << " " << endl;
	}
}

void CompareBitsDistinctBM(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int size = Size / 64; size <= Size / 8; size += Size / 64){
		T_8bit t_4(4, 1, Interval, size * 2, 4, hashes);
		T_8bit t_6(6, 1, Interval, size * 3 / 2, 4, hashes);
		T_8bit t_8(8, 1, Interval, size, 4, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t4Error = 0, t6Error = 0, t8Error = 0;
		int gap = 12345;
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_4.insert(flow[i].first, now);
			t_6.insert(flow[i].first, now);
			t_8.insert(flow[i].first, now);
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;

				t4Error += abs(double(unique.size()) - double(t_4.countflow(now))) / (unique.size() * 1.0);
				t6Error += abs(double(unique.size()) - double(t_6.countflow(now))) / (unique.size() * 1.0);
				t8Error += abs(double(unique.size()) - double(t_8.countflow(now))) / (unique.size() * 1.0);
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t4Error / r << " ";
		cout << "the RE of 4: " << t4Error / r << " " << endl;

		fp << t6Error / r << " ";
		cout << "the RE of 6: " << t6Error / r << " " << endl;

		fp << t8Error / r << endl;
		cout << "the RE of 8: " << t8Error / r << " " << endl;
	}
}

void CompareNCC(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int size = Size; size <= Size * 2; size += Size / 8){
		T_8bit t_5NCC(8, 8, Interval, size, 4, hashes);
		T_8bit t_10NCC(8, 8, Interval, size, 6, hashes);
		T_8bit t_20NCC(8, 8, Interval, size, 8, hashes);
		T_8bit t_ONCC(8, 8, Interval, size, 32, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double t_5NCCFPR = 0, t_10NCCFPR = 0, t_20NCCFPR = 0, t_ONCCFPR = 0;
		int r = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(Interval * 10 / pow(2, BitLength - 1))) * int(Interval * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_5NCC.insert(flow[i].first, now);
			t_10NCC.insert(flow[i].first, now);
			t_20NCC.insert(flow[i].first, now);
			t_ONCC.insert(flow[i].first, now);
			
			if ((i + 1) % gap == 0){
				r++;
				for (int j = 0; j < t_ONCC.size; j += t_ONCC.NCC)
					t_ONCC.CheckAndClear(j, now, 1);
				fal = neg = pos = 0;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);


				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_5NCC.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_5NCCFPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;

				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_10NCC.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_10NCCFPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_20NCC.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_20NCCFPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_ONCC.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_ONCCFPR += (pos * 1.0 / fal);

			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t_5NCCFPR / r << " ";
		cout << "the FPR of 5: " << t_5NCCFPR / r << " " << endl;

		fp << t_10NCCFPR / r << " ";
		cout << "the FPR of 10: " << t_10NCCFPR / r << " " << endl;

		fp << t_20NCCFPR / r << " ";
		cout << "the FPR of 20: " << t_20NCCFPR / r << " " << endl;

		fp << t_ONCCFPR / r << endl;
		cout << "the FPR of O: " << t_ONCCFPR / r << " " << endl;
	}
}

void CompareNCCDistinct(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int interv = Interval; interv <= (2 * Interval); interv += (Interval / 4)){
		cout << "interv: " << interv << endl;
		T_8bit t_5NCC(8, 2, interv, Size/2, 16, hashes);
		T_8bit t_10NCC(8, 2, interv, Size/2, 32, hashes);
		T_8bit t_20NCC(8, 2, interv, Size/2, 64, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t5NCCError = 0, t10NCCError = 0, t20NCCError = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_5NCC.insert(flow[i].first, now);
			t_10NCC.insert(flow[i].first, now);
			t_20NCC.insert(flow[i].first, now);
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0,i - interv * 10); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << interv * 10 << " items" << endl;

				t5NCCError += abs(double(unique.size()) - double(t_5NCC.getDistinctItems(now))) / (unique.size() * 1.0);
				cout << t_5NCC.getDistinctItems(now) << endl;
				t10NCCError += abs(double(unique.size()) - double(t_10NCC.getDistinctItems(now))) / (unique.size() * 1.0);
				cout << t_10NCC.getDistinctItems(now) << endl;
				t20NCCError += abs(double(unique.size()) - double(t_20NCC.getDistinctItems(now))) / (unique.size() * 1.0);
				cout << t_20NCC.getDistinctItems(now) << endl;
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t5NCCError / r << " ";
		cout << "the RE of 5: " << t5NCCError / r << " " << endl;

		fp << t10NCCError / r << " ";
		cout << "the RE of 10: " << t10NCCError / r << " " << endl;

		fp << t20NCCError / r << endl;
		cout << "the RE of 20: " << t20NCCError / r << " " << endl;
	}
}

void CompareNCCDistinctBM(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int size = Size / 2; size <= Size * 2; size += Size / 2){
		T_8bit t_5NCC(8, 1, Interval, size, 8, hashes);
		T_8bit t_10NCC(8, 1, Interval, size, 16, hashes);
		T_8bit t_20NCC(8, 1, Interval, size, 32, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t5NCCError = 0, t10NCCError = 0, t20NCCError = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(Interval * 10 / pow(2, BitLength - 1))) * int(Interval * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_5NCC.insert(flow[i].first, now);
			t_10NCC.insert(flow[i].first, now);
			t_20NCC.insert(flow[i].first, now);
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << Interval * 10 << " items" << endl;
				//cout << "r: " << r << endl;
				t5NCCError += abs(double(unique.size()) - double(t_5NCC.countflow(now))) / (unique.size() * 1.0);
				//cout << "5: " << abs(double(unique.size()) - double(t_5NCC.countflow(now))) << endl;
				t10NCCError += abs(double(unique.size()) - double(t_10NCC.countflow(now))) / (unique.size() * 1.0);
				//cout << "10: " << abs(double(unique.size()) - double(t_10NCC.countflow(now))) << endl;
				t20NCCError += abs(double(unique.size()) - double(t_20NCC.countflow(now))) / (unique.size() * 1.0);
				//cout << "20: " << abs(double(unique.size()) - double(t_20NCC.countflow(now))) << endl;
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t5NCCError / r << " ";
		cout << "the RE of 5: " << t5NCCError / r << " " << endl;

		fp << t10NCCError / r << " ";
		cout << "the RE of 10: " << t10NCCError / r << " " << endl;

		fp << t20NCCError / r << endl;
		cout << "the RE of 20: " << t20NCCError / r << " " << endl;
	}
}


void CompareHash(int flowSize, int querySize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int h = 2; h <= 16; h+= 2){
		T_8bit t_64K(8, h, Interval, Size / 2, 4, hashes);
		T_8bit t_128K(8, h, Interval, Size, 4, hashes);
		T_8bit t_96K(8, h, Interval, Size * 3 / 4, 4, hashes);
		set<uint32_t> unique;
		double fal, neg, pos;
		double t_64KFPR = 0, t_128KFPR = 0, t_96KFPR = 0;
		int r = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(Interval * 10 / pow(2, BitLength - 1))) * int(Interval * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_64K.insert(flow[i].first, now);
			t_128K.insert(flow[i].first, now);
			t_96K.insert(flow[i].first, now);
			
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0, i - int(Interval * 10)); j--)
					unique.insert(flow[j].first);

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_64K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_64KFPR += (pos * 1.0 / fal);


				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_128K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_128KFPR += (pos * 1.0 / fal);

				fal = neg = pos = 0;
				for (int j = 0; j < querySize; j++)
					if (unique.count(QuerySet[j]) == 0)
					{
						fal++;
						if (t_96K.query(QuerySet[j], now) == 1)
							pos++;
					}

				t_96KFPR += (pos * 1.0 / fal);

			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t_64KFPR / r << " ";
		cout << "the FPR of 64: " << t_64KFPR / r << " " << endl;

		fp << t_96KFPR / r << " ";
		cout << "the FPR of 96: " << t_96KFPR / r << " " << endl;

		fp << t_128KFPR / r << endl;
		cout << "the FPR of 128: " << t_128KFPR / r << " " << endl;
	}
}

void CompareHashDistinct(int flowSize){
	double now = 0;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);

	for (int interv = Interval; interv <= (2 * Interval); interv += (Interval / 4)){
		cout << "interv: " << interv << endl;
		T_8bit t_4h(8, HashNum/2, interv, Size/2, 64, hashes);
		T_8bit t_6h(8, HashNum*3/4, interv, Size/2, 64, hashes);
		T_8bit t_8h(8, HashNum, interv, Size/2, 64, hashes);
		set<uint32_t> unique;
		int r = 0;
		double t4hError = 0, t6hError = 0, t8hError = 0;
		int gap = flow.size() / 8;
		gap = (gap / int(interv * 10 / pow(2, BitLength - 1))) * int(interv * 10 / pow(2, BitLength - 1));
		cout << "test every " << gap << " insertion" << endl;
		for (int i = 0; i < flow.size(); i++)
		{
			now = double(i) / 10;
			t_4h.insert(flow[i].first, now);
			t_6h.insert(flow[i].first, now);
			t_8h.insert(flow[i].first, now);
			
			if ((i + 1) % gap == 0){
				r++;
				unique.clear();
				for (int j = i; j > max(0,i - interv * 10); j--)
					unique.insert(flow[j].first);
				cout << unique.size() << " distinct items in recent " << interv * 10 << " items" << endl;

				t4hError += (abs(double(unique.size()) - double(t_4h.getDistinctItems(now))) / (unique.size() * 1.0));
				t6hError += (abs(double(unique.size()) - double(t_6h.getDistinctItems(now))) / (unique.size() * 1.0));
				t8hError += (abs(double(unique.size()) - double(t_8h.getDistinctItems(now))) / (unique.size() * 1.0));
			}
		}
		cout << "average of " << r << " times" << endl;

		fp << t4hError / r << " ";
		cout << "the RE of 4: " << t4hError / r << " " << endl;

		fp << t6hError / r << " ";
		cout << "the RE of 6: " << t6hError / r << " " << endl;

		fp << t8hError / r << endl;
		cout << "the RE of 8: " << t8hError / r << " " << endl;
	}
}

void CompareSpeedWithOther(int flowSize, int querySize){
	double now = 0;
	clock_t start, finish;
	double duration;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);
	T_8bit t_8bit(8, HashNum, Interval, Size, 64, hashes);
	Swamp swamp(HashNum, Interval, Size, hashes);
	TBF tbf(HashNum, Interval, Size, hashes);
	
	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		t_8bit.insert(flow[i].first, now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		swamp.insert(flow[i].first);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		tbf.insert(flow[i].first, now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << endl;



	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		t_8bit.query(QuerySet[i], now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		swamp.query(QuerySet[i]);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		tbf.query(QuerySet[i], now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl; 
	fp << 0.1 / duration << endl;
}

void CompareDistinctSpeedWithOther(int flowSize, int querySize){
	double now = 0;
	clock_t start, finish;
	double duration;
	vector<pair<uint32_t, double>> flow;
	InitDataSet(flow, flowSize);
	T_8bit t_8bit(8, HashNum, Interval, Size, 64, hashes);
	T_8bit t_8bitBM(8, 1, Interval, Size, 64, hashes);
	Swamp swamp(HashNum, Interval, Size, hashes);
	TSV tsv(Interval, Size, hashes);
	CVS cvs(Interval, Size, hashes);

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		t_8bit.insert(flow[i].first, now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		t_8bitBM.insert(flow[i].first, now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		swamp.insert(flow[i].first);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		tsv.insert(flow[i].first, i);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << " ";

	start = clock();
	for (int i = 0; i < flow.size(); i++)
	{
		now = double(i) / 10;
		cvs.insert(flow[i].first);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.1 / duration << endl;





	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		now = double(i) / 10;
		t_8bit.getDistinctItems(now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 1 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		now = double(i) / 10;
		t_8bitBM.countflow(now);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 1 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		now = double(i) / 10;
		swamp.getDistinctItems();
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 1 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize / 100; i++)
	{
		now = double(i) / 10;
		tsv.countflow(i);
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 0.01 / duration << " ";

	start = clock();
	for (int i = 0; i < querySize; i++)
	{
		now = double(i) / 10;
		cvs.countflow();
	}
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "duration: " << duration << endl;
	fp << 1 / duration << endl;
}


int main(){

	printf("Main running, init random reed\n");
	srand((time(0)));

	InitBOBHash(100);
	printf("Hash initialized\n");

	CreateQuerySet(QueryNum);
	printf("Query set created\n");

	fp.open("result.txt");
	printf("Init ofstream\n");
	
	//CompareFPRwithTime(500000, QueryNum);
	//CompareDistinctwithTime(500000);
	//CompareBMDistinctwithTime(500000);
	//TimeBasedCompareFPRwithTime(500000, QueryNum);
	//TimeBasedCompareDistinctwithTime(500000);
	//TimeBasedCompareBMDistinctwithTime(500000);
	//CompareFPRwithDifferentMemory(120000, QueryNum);
	//CompareDistinctElements(120000);
	//CompareBMDistinctElements(120000);
	//TimeBasedCompareFPRwithDifferentMemory(500000, QueryNum);
	//TimeBasedCompareDistinctElements(500000);
	//TimeBasedCompareBMDistinctElements(500000);
	//CompareFPRWithOther(120000, QueryNum);
	//CompareDistinctElementsWithOther(120000);
	//CompareBits(120000, QueryNum);
	//CompareBitsDistinct(120000);
	//CompareBitsDistinctBM(120000);
	CompareNCC(240000, QueryNum);
	//CompareNCCDistinct(1200000);
	//CompareNCCDistinctBM(1200000);
	//CompareHash(120000, QueryNum);
	//CompareHashDistinct(120000);
	//CompareSpeedWithOther(100000, QueryNum);
	//CompareDistinctSpeedWithOther(100000, QueryNum * 10);
	system("pause");
}