#pragma once
struct mysignal
{
	double loss;
	double distance;
	mysignal* next = nullptr;
};
class mypoint
{
public:
	mypoint *parents = nullptr;
	int depth;
	int n_wall;//该点所对应的墙的序号从0取到19
	double loss_db;
	double x;
	double y;
	static int pointcount;
	mypoint *children[20];//20是墙个数上限，修改墙的上限时需要修改这个
	int n_children = 0;
	bool can_reach_target = 0;
	mypoint();
	mypoint( int depth, double loss_db, double x, double y);
	~mypoint();
	double get_lossdb_by_distant(mypoint* p);
	static void release_the_tree(mypoint* p);
	static void ergodic_the_tree(mypoint* p, mysignal** s, mypoint target);
	static void release_the_signal(mysignal* s);
	static double* generate_signal(int* sequence, mysignal* signal_struct);
	static double* generate_signal_without_noize(int* sequence, mysignal* signal_struct);
};

void test();


