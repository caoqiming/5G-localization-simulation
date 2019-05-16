#pragma once
class mywall
{
public:
	double start_x;
	double start_y;
	double end_x;
	double end_y;
	double loss_db;
	mywall(double x1, double y1, double x2, double y2, double loss);
	mywall();
	~mywall();
	void copy(mywall copy_from);
	void changedata(double x1, double y1, double x2, double y2, double loss);
};

