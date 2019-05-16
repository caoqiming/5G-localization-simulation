#include "mywall.h"



mywall::mywall(double x1, double y1, double x2, double y2, double loss)
{
	this->loss_db = loss;
	this->start_x = x1;
	this->start_y = y1;
	this->end_y = y2;
	this->end_x = x2;
}

mywall::mywall()
{

}

mywall::~mywall()
{
}

void mywall::copy(mywall copy_from)
{
	this->end_x = copy_from.end_x;
	this->end_y = copy_from.end_y;
	this->start_x = copy_from.start_x;
	this->start_y = copy_from.start_y;
	this->loss_db = copy_from.loss_db;
}
void mywall::changedata(double x1, double y1,double x2, double y2, double loss)
{
	this->end_x = x2;
	this->end_y = y2;
	this->start_x = x1;
	this->start_y = y1;
	this->loss_db = loss;

}