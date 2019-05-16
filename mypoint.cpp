#include<math.h>
#include "mypoint.h"
#include"localization.h"
#include<iostream>
#include<time.h>

using namespace std;
mypoint::mypoint(int depth, double loss_db, double x, double y)
{
	this->pointcount++;
	this->depth = depth;
	this->loss_db = loss_db;
	this->x = x;
	this->y = y;
}
mypoint::mypoint()
{
	this->n_wall = 0;
	this->pointcount++;
	this->depth = 0;
	this->loss_db = 5;
	this->x = -1;
	this->y = -1;
}
mypoint::~mypoint()
{
}

double mypoint::get_lossdb_by_distant(mypoint* p)
{
	double distance = pow(pow(p->x - this->x,2)+pow(p->y - this->y,2),0.5);
	if (distance == 0)
		return 1000000000000;
	return 20 * log10(distance/1000);
}

void mypoint::release_the_tree(mypoint* p)
{
	for(int k=0;k<20;k++)
	{
		if (p->children[k])
		{
			release_the_tree(p->children[k]);
		}
	}
	for (int j = 0;j < 20;j++)
	{
		p->children[j] = nullptr;
	}
	p->n_children = 0;
	if(p->parents)
		delete p;
}

void mypoint::ergodic_the_tree(mypoint * p, mysignal** s,mypoint target)
{
	mysignal* temp=*s;
	if (p->can_reach_target)
	{
		mysignal* newsignal= new mysignal;
		temp->next = newsignal;
		newsignal->distance = pow(pow(target.x - p->x, 2) + pow(target.y - p->y, 2), 0.5);
		newsignal->loss = p->loss_db;
		*s = newsignal;
	}
	for (int i = 0;i < p->n_children;i++)
	{
		ergodic_the_tree(p->children[i], s, target);
	}
}

void mypoint::release_the_signal(mysignal * s)
{

	if (s->next)
	{
		release_the_signal(s->next);
	}

	delete s;

}

double* mypoint::generate_signal(int* sequence, mysignal* signal_struct)//signal是单链表，不是数组
{
	if (!signal_struct->next)//若为空
	{
		double* empty = new double[RECEIVE_LENTH];
		for (int i = 0;i < RECEIVE_LENTH;i++)
		{
			empty[i] = 0;
		}
		return empty;
	}
	mysignal* signal_struct_now = signal_struct->next;//第一个节点不包含信息也不计入个数
	if (!signal_struct_now)
	{
		double* signal = new double[RECEIVE_LENTH];
		for (int i = 0;i < RECEIVE_LENTH;i++)
		{
			signal[i] = gaussrand()*NOIZE_BZC;
		}
		
		return signal;
	}
	double* signal = new double[RECEIVE_LENTH];//若场景变大或者最低db有变这个数组需要扩大
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal[i] = 0;
	}
	double*  sequence_double = new double[MSEQUENCE_LENTH];
	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		sequence_double[i] = double(sequence[i]);
	}
	while (true)
	{
		
		int delay = (int)floor(signal_struct_now->distance / 2.44);
		double loss_db_all = signal_struct_now->loss + 20 * log10(signal_struct_now->distance / 1000) + 32.4478 + 20 * log10(FREQUENCE);
		double loss_D = pow(10, 0.1*loss_db_all);

		for (int j = 0;j < MSEQUENCE_LENTH;j++)
		{
			signal[j + delay] += sequence_double[j] / loss_D;
		}

		if (signal_struct_now->next)//后面还有
		{
			signal_struct_now = signal_struct_now->next;
		}
		else
		{
			break;
		}
		

	}

	delete[] sequence_double;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal[i] += gaussrand()*NOIZE_BZC;
	}
	return signal;
}
double* mypoint::generate_signal_without_noize(int* sequence, mysignal* signal_struct)//signal是单链表，不是数组
{
	if (!signal_struct->next)//若为空
	{
		double* empty = new double[RECEIVE_LENTH];
		for (int i = 0;i < RECEIVE_LENTH;i++)
		{
			empty[i] = 0;
		}
		return empty;
	}

	mysignal* signal_struct_now = signal_struct->next;//第一个节点不包含信息也不计入个数
	double* signal = new double[RECEIVE_LENTH];//若场景变大或者最低db有变这个数组需要扩大
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal[i] = 0;
	}
	double*  sequence_double = new double[MSEQUENCE_LENTH];
	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		sequence_double[i] = double(sequence[i]);
	}
	while (true)
	{

		int delay = (int)floor(signal_struct_now->distance / 2.44);
		double loss_db_all = signal_struct_now->loss + 20 * log10(signal_struct_now->distance / 1000) + 32.4478 + 20 * log10(FREQUENCE);
		double loss_D = pow(10, 0.1*loss_db_all);

		for (int j = 0;j < MSEQUENCE_LENTH;j++)
		{
			signal[j + delay] += sequence_double[j] / loss_D;
		}

		if (signal_struct_now->next)//后面还有
		{
			signal_struct_now = signal_struct_now->next;
		}
		else
		{
			break;
		}


	}

	delete[] sequence_double;
	return signal;
}


void test()
{
	WALLNUMBER = 0;
}