#include <math.h>
#include<iostream>
#include <algorithm>
#include<Dense>
#include "localization.h"

using namespace std;
int WALLNUMBER = 4;
int WALLNUMBER_SAVE= WALLNUMBER;
int WALLNUMBER_SPECULATE = WALLNUMBER;
double NOIZE_BZC = 1.0e-10;
int* m_sequence()
{
	const int xishu1[10] = { 0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 };
	int seg[11] = { 1 };
	int backq = 0;
	int* sequence = new int[MSEQUENCE_LENTH];//2^10-1
	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		sequence[i] = seg[10];
		backq = 0;
		for (int j = 0;j < 10;j++)
		{
			backq += seg[j] * xishu1[j];
		}
		backq += seg[10];
		for (int k = 10;k > 0;k--)
			seg[k] = seg[k - 1];

		if (backq % 2 == 0)
			seg[0] = 0;
		else
			seg[0] = 1;
	}

	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		if (sequence[i] == 0)
		{
			sequence[i] = -1;
		}
	}

	return sequence;
}
mypoint* get_symmetric_point(mypoint* p, mywall ab)
{
	double vector_ab_x = ab.end_x - ab.start_x;
	double vector_ab_y = ab.end_y - ab.start_y;

	double a = 0;//ab的方程ax+by+c=0
	double b = 0;

	if (vector_ab_y == 0)
	{
		a = 0;
		b = 1;
	}
	else
	{
		a = 1;
		b = vector_ab_x / vector_ab_y;
		double l = pow(a * a + b * b, 0.5);
		a = a / l;
		b = b / l;
	}
	double c = -(a*ab.start_x + b * ab.start_y);
	double distance = abs(a*p->x + b * p->y + c);//点到线的距离
	if (distance == 0)
	{
		mypoint* p_new = new mypoint(p->depth + 1, p->loss_db + ab.loss_db, p->x, p->y);
		p_new->parents = p;
		return p_new;//点在直线上就直接返回p
	}

	if (a*(ab.start_x - p->x) + b * (ab.start_y - p->y) < 0)//法向量需要反向
	{
		a = -a;
		b = -b;
	}

	double c_x;
	double c_y;
	c_x = p->x + 2 * distance*a;
	c_y = p->y + 2 * distance*b;
	mypoint* p_new = new mypoint(p->depth + 1, p->loss_db + ab.loss_db, c_x, c_y);
	p_new->parents = p;
	return p_new;
}
bool checkblock1(mypoint p, mypoint q, mywall ab)
{//由p到q障碍物有ab时返回真 19.5.10更新 若点在墙上不算被挡住
	double a, b;


	double ax = p.x;
	double ay = p.y;
	double bx = q.x;
	double by = q.y;
	double cx = ab.start_x;
	double cy = ab.start_y;
	double dx = ab.end_x;
	double dy = ab.end_y;
	if (cx == dx)
	{
		if (ax == cx or bx == cx)
			return false;
	}
	else
	{
		a = (cy - dy) / (cx - dx);
		b = cy - a * cx;
		if (abs(ay - a * ax - b )<1e-8 or abs(by - a * bx - b)<1e-8)
			return false;
	}





	if (!(min(ax, bx) <= max(cx, dx) && min(cy, dy) <= max(ay, by) && min(cx, dx) <= max(ax, bx) && min(ay, by) <= max(cy, dy)))
		return false;

	double u, v, w, z;//分别记录两个向量
	u = (cx - ax)*(by - ay) - (bx - ax)*(cy - ay);
	v = (dx - ax)*(by - ay) - (bx - ax)*(dy - ay);
	w = (ax - cx)*(dy - cy) - (dx - cx)*(ay - cy);
	z = (bx - cx)*(dy - cy) - (dx - cx)*(by - cy);
	return (u*v <= 0.00000001 && w*z <= 0.00000001);

}
bool checkblock2(mypoint p, mypoint q, mywall wall[])
{//由p到q障碍物只有1个（前面已经判断有ab）时返回真
	int n = 0;
	for (int i = 0;i < WALLNUMBER;i++)
	{
		if (checkblock1(p, q, wall[i]))
			n++;
	}
	return n == 1;
}
bool checkblock3(mypoint p, mypoint q, mywall wall[])
{//由p到q没有障碍物时返回真
	int n = 0;
	for (int i = 0;i < WALLNUMBER;i++)
	{
		if (checkblock1(p, q, wall[i]))
			n++;
	}
	return n == 0;
}
just_point* get_intersection(just_point p1, just_point p2, just_point p3, just_point p4)
{
	double a1,a2,b1,b2;
	if (p1.x == p2.x)
	{
		a2 = (p3.y - p4.y) / (p3.x - p4.x);
		b2 = p3.y - a2 * p3.x;
		just_point* ans=new just_point;
		ans->x = p1.x;
		ans->y = a2 * ans->x + b2;
		return ans;
	}
	else 
	{
		a1 = (p1.y - p2.y) / (p1.x - p2.x);
		b1 = p1.y - a1 * p1.x;
		if (p3.x == p4.x)
		{
			just_point* ans = new just_point;
			ans->x = p3.x;
			ans->y = a1 * ans->x + b1;
			return ans;
		}
	}

	a2 = (p3.y - p4.y) / (p3.x - p4.x);
	b2 = p3.y - a2 * p3.x;
	just_point* ans = new just_point;
	ans->x = (b2 - b1) / (a1 - a2);
	ans->y = a1* ans->x+b1;
	return ans;
}
bool super_checkblock(mypoint* p, mypoint* target,mywall* wall)
{
	//因为是迭代算法，传进来的anchor并不是真的anchor 所以这里要自己找真的anchor
	mypoint* r_anchor = p;
	while(r_anchor->parents)
	{
		r_anchor = r_anchor->parents;
	}



	bool result = true;
	mypoint* pp = p;
	just_point *ppp=new just_point;
	ppp->front = nullptr;//第一个指针不记录数据
	ppp->x = -1;
	ppp->y = -1;
	int miao = 0;
	while (1)
	{
		just_point p1, p2, p3, p4;
		p1.x = pp->x;
		p1.y = pp->y;
		if (miao == 0)
		{
			p2.x = target->x;
			p2.y = target->y;
		}
		else
		{
			p2.x = ppp->x;
			p2.y = ppp->y;
		}


		p3.x = wall[pp->n_wall].start_x;
		p3.y = wall[pp->n_wall].start_y;
		p4.x = wall[pp->n_wall].end_x;
		p4.y = wall[pp->n_wall].end_y;
		ppp->next=get_intersection(p1,p2,p3,p4);
		ppp->next->front = ppp;
		ppp = ppp->next;
		pp = pp->parents;
		if (pp == r_anchor)
			break;

		miao++;
	}
	miao = 0;
	while(1)
	{
		mypoint temp;
		temp.x = ppp->x;
		temp.y = ppp->y;

		if (miao == 0)//第一次运行是考虑基站与点ppp
		{
			if (!checkblock3(temp, *r_anchor, wall))//有障碍物
			{
				result = false;
				break;
			}
		}

		if (ppp->front and ppp->front->x==-1 and ppp->front->y == -1)
		{
			//判断最后一个是ppp与目标
			if (!checkblock3(temp, *target, wall))//有障碍物
			{
				result = false;
				break;
			}
			break;
		}
		else
		{//中间的判断部分
			mypoint temp2;
			temp2.x = ppp->front->x;
			temp2.y = ppp->front->y;
			if (!checkblock3(temp, temp2, wall))//有障碍物
			{
				result = false;
				//cout << "miao?";
				break;
			}

		}
		
		ppp = ppp->front;
		miao++;
	} 
	ppp = ppp->front;
	while (ppp->next)
	{
		ppp=ppp->next;
		delete ppp->front;
	}
	delete ppp;
	return result;
}

bool check_loss(double loss_fanshe, double loss_juli)
{
	if (loss_fanshe + loss_juli + 32.4478 + 20 * log10(FREQUENCE) < HIGHESTDB)
	{

		//cout <<"当前loss"<< loss_fanshe + loss_juli + 32.4478 + 20 * log10(FREQUENCE)<<endl;
		return 1;
	}
	else
		return 0;
}
void build_the_tree(mywall wall[], mypoint* target, mypoint* anchor)//anchor是根节点
{
	for (int i = 0;i < WALLNUMBER;i++)
	{
		mypoint* p_child = get_symmetric_point(anchor, wall[i]);
		p_child->n_wall = i;
		//cout<<"当前的点在： ("<< p_child->x <<","<< p_child ->y <<")当前对称的墙是："<<i<<endl;
		if (check_loss(p_child->loss_db, target->get_lossdb_by_distant(p_child)) && checkblock1(*p_child, *target, wall[i]))
		{
			//cout << "能量合格，check1有效";
			anchor->children[anchor->n_children] = p_child;
			anchor->n_children++;


			if (checkblock2(*p_child, *target, wall) and super_checkblock(p_child, target, wall))
			{
				p_child->can_reach_target = 1;
				//cout << " 无其他阻挡";
			}
			else
			{
				p_child->can_reach_target = 0;
				//cout << " 有其他阻挡";
			}
			//cout << "继续迭代" << endl << endl;
			build_the_tree(wall, target, p_child);
		}
		else
		{
			//cout << "无效"<<endl<<endl;
			delete p_child;
			mypoint::pointcount--;
		}
	}

}



TOA_MPC* get_one_mpc(double simulation_signal[RECEIVE_LENTH],int* sequence)//获取最大的mpc
{
	const int delay_rang = RECEIVE_LENTH - MSEQUENCE_LENTH;
	double inner_product[delay_rang] = { 0 };
	double max_inner_product =-1;
	int delay=0;
	double amplitude=0;
	for (int delay_int = 0;delay_int < delay_rang;delay_int++)
	{
		for (int i = 0;i < MSEQUENCE_LENTH;i++)
		{
			inner_product[delay_int] += simulation_signal[i+ delay_int] * sequence[i];
		}
		if (inner_product[delay_int] > max_inner_product)
		{
			max_inner_product = inner_product[delay_int];
			delay = delay_int;
		}
	}
	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		amplitude += simulation_signal[i + delay] * sequence[i];
	}
	amplitude = amplitude / MSEQUENCE_LENTH;
	TOA_MPC* mpc = new TOA_MPC;
	mpc->amplitude = amplitude;
	mpc->delay = double(delay);
	return mpc;
}
TOA_MPC* seperate_mpc(double* simulation_signal, int* sequence)
{
	double signal[RECEIVE_LENTH];
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal[i] = simulation_signal[i];
	}

	TOA_MPC* mpc_root = new TOA_MPC;
	TOA_MPC** mpc_temp = &mpc_root;
	TOA_MPC* mpc_return = mpc_root;
	TOA_MPC* current_mpc = 0;
	while(1)
	{	
		current_mpc = get_one_mpc(signal, sequence);
		if (current_mpc->amplitude > NOIZE_BZC and current_mpc->amplitude>1e-10)//该mpc分量在识别范围内
		{
			TOA_MPC* temp = *mpc_temp;
			temp->next = current_mpc;
			*mpc_temp = current_mpc;
			for (int i = 0;i < MSEQUENCE_LENTH;i++)
			{
				signal[i + int(current_mpc->delay)] -= sequence[i]* current_mpc->amplitude;
			}

		}
		else
		{
			delete current_mpc;
			break;
		} 
	}

	return mpc_return;
}
mysignal* mpc_to_mysignal(TOA_MPC* p)//mpc格式转mysignal格式 只在分析mpc性能里用了
{
	mysignal* s = new mysignal;
	mysignal* ss = s;


	if (! p->next)
	{
		//cout<<"mpc_to_mysignal输入mpc为空"<<endl;
		return s;
	}
	

	while (true)
	{
		if (p->next)
		{
			p = p->next;
			s->next = new mysignal;
			s = s->next;
			s->distance = p->delay * 2.4414;
			s->loss = -10*log10(p->amplitude)-103.33-20*log10(s->distance /1000);
		}
		else
		{
			break;
		}
	}
	return ss;
}
void print_mysignal(mysignal* s)//输出多径信息
{
	cout.setf(std::ios::left);
	cout.width(15);
	cout << s->distance;
	cout<<",";
	cout.width(5);
	cout << s->loss+20 * log10(s->distance / 1000) + 32.4478 + 20 * log10(FREQUENCE) <<endl;
	if (s->next)
	{
		print_mysignal(s->next);
	}
}
void release_toa_mpc(TOA_MPC* p)
{
	if (p->next)
	{
		release_toa_mpc(p->next);
	}

	delete p;

}
void release_toa_mpc(TOA_MPC* p,int depth)
{

	if (p->next != nullptr)
	{
		release_toa_mpc(p->next, depth+1);
	}
	if (depth != 0)
	{
		delete p;
		TOA_MPC::number--;
	}

}
TRAN_LOS_DATA GET_LOS_DATA(double* simulation_signal)
{
	//一下两个参数的单位不是秒而分别是8.138*10^-9秒和(8.138*10^-9)^2秒^2
	/*  Mean excess delay : It is an indicator of how quickly
		the signal decays after the MPC corresponding to the
		LoS path.A small value is indicative of LoS conditions  */
	double mean_excess_delay;
	double Numerator = 0;
	double Denominator = 0;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		double square= pow(simulation_signal[i], 2);
		Numerator += i * square;
		Denominator += square;
	}
	mean_excess_delay = Numerator / Denominator;


	/*
		Root mean square (RMS) delay spread: This metric
	    is the standard deviation of the excess delay, relative
	    to the LoS component and like τm, a small value is a
	    strong indicator of LoS
	*/
	double RMS_delay_spread;
	Numerator = 0;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		Numerator += pow(i- mean_excess_delay,2) * pow(simulation_signal[i],2);
	}
	RMS_delay_spread = Numerator / Denominator;

	/*
	    Amplitude kurtosis,This is a measure of the peakiness 
		of hi(t),a large value being indicative ofLoS conditions
	*/

	double Amplitude_kurtosis;
	double mean=0;
	double standard_deviation=0;
	double signal_abs[RECEIVE_LENTH];
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal_abs[i] = abs(simulation_signal[i]);
	}

	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		mean += signal_abs[i];
	}
	mean = mean / RECEIVE_LENTH;//均值
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		standard_deviation += pow(signal_abs[i]- mean,2);
	}
	standard_deviation = pow(standard_deviation / (RECEIVE_LENTH-1),0.5);//标准差
	Numerator = 0;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		Numerator += pow(signal_abs[i] - mean, 4);
	}
	Amplitude_kurtosis=Numerator / (RECEIVE_LENTH*standard_deviation);

	/*Total received power:*/
	double power= Denominator/ RECEIVE_LENTH;

	/*Maximum signal amplitude*/
	double amplitude_max = 0;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		if (signal_abs[i] > amplitude_max)
			amplitude_max = signal_abs[i];
	}
	/*ToA of the first MPC*/
	//mpc->next->delay
	TRAN_LOS_DATA data;
	data.Amplitude_kurtosis = Amplitude_kurtosis;
	data.amplitude_max = amplitude_max;
	data.mean_excess_delay = mean_excess_delay;
	data.power = power;
	data.RMS_delay_spread = RMS_delay_spread;
	return data;
}
TOA_MPC** get_fingerprints(mypoint* anchor,mywall* wall)
{
	TOA_MPC** fingerprints=new TOA_MPC*[max_x];
	for (int i = 0;i < max_x;i++)
	{
		fingerprints[i] = new TOA_MPC[max_y];
	}

	mypoint** v_target=new mypoint*[max_x];
	for (int i = 0;i < max_x;i++)
	{
		v_target[i] = new mypoint[max_y];
	}

	int* sequence = m_sequence();
	for (int i = 0;i < max_x;i++)//开始获取指纹信息
	{
		for (int j = 0;j < max_y;j++)
		{
			cout <<"指纹信息已获取:"<< double(i* max_y +j)/(max_x*max_y)*100+0.01<<"%\r";
			v_target[i][j].x = i;
			v_target[i][j].y = j;
			if (checkblock3(v_target[i][j], *anchor, wall))
				anchor->can_reach_target = 1;//有line of sight
			build_the_tree(wall, &v_target[i][j], anchor);
			mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
			mysignal* signal = signal_temp;
			mypoint::ergodic_the_tree(anchor, &signal_temp, v_target[i][j]);
			double* simulation_signal = mypoint::generate_signal_without_noize(sequence, signal);//指纹信息本来是要多取几个平均的，这里直接不要噪声效果一样，复杂度更低
			TOA_MPC* mpc = seperate_mpc(simulation_signal, sequence);
			fingerprints[i][j] = *mpc;
			delete[] simulation_signal;
			mypoint::release_the_tree(anchor);
			mypoint::release_the_signal(signal);
		}
	}
	delete[] sequence;
	//释放new出来的格点
	for (int i = 0; i < max_x; ++i)
	{
		delete[] v_target[i];
	}
	delete[] v_target;

	mypoint::pointcount -= 10000;
	return fingerprints;
}
TOA_MPC* simulate_once(mypoint* anchor, mypoint* target, mywall* wall)
{//返回仿真一次之后分离出的多径信息。 用了之后需要release_toa_mpc
	if (checkblock3(*target, *anchor, wall))
		anchor->can_reach_target = 1;//有line of sight
	build_the_tree(wall, target, anchor);
	mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
	mysignal* signal = signal_temp;
	mypoint::ergodic_the_tree(anchor, &signal_temp, *target);
	int* sequence = m_sequence();
	double* simulation_signal = mypoint::generate_signal(sequence, signal);
	TOA_MPC* mpc = seperate_mpc(simulation_signal, sequence);
	delete[] simulation_signal;
	delete[] sequence;
	mypoint::release_the_tree(anchor);
	mypoint::release_the_signal(signal);
	anchor->can_reach_target = false;
	return mpc;
}
double* simulate_once_receive(mypoint* anchor, mypoint* target, mywall* wall)//模拟收到信号
{
	if (checkblock3(*target, *anchor, wall))
		anchor->can_reach_target = 1;//有line of sight
	build_the_tree(wall, target, anchor);
	mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
	mysignal* signal = signal_temp;
	mypoint::ergodic_the_tree(anchor, &signal_temp, *target);
	int* sequence = m_sequence();
	double* simulation_signal = mypoint::generate_signal(sequence, signal);
	delete[] sequence;
	mypoint::release_the_tree(anchor);
	mypoint::release_the_signal(signal);
	return simulation_signal;
}
double* generate_signal_by_mpc(TOA_MPC mpc, int* sequence)//与树生成信号不同，这里利用mpc生成信号。
{
	TOA_MPC* p= &mpc;
	double* signal = new double[RECEIVE_LENTH];
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		signal[i] = 0;
	}
	double*  sequence_double = new double[MSEQUENCE_LENTH];
	for (int i = 0;i < MSEQUENCE_LENTH;i++)
	{
		sequence_double[i] = double(sequence[i]);
	}

	while(p->next)
	{
		p = p->next;
		for (int j = 0;j < MSEQUENCE_LENTH;j++)
		{
			signal[j+int(p->delay)] = sequence_double[j] * p->amplitude;
		}
	}

	delete[]sequence_double;
	return signal;
}
just_point localization_by_fingerprints(TOA_MPC** fingerprints, TOA_MPC* target_mpc)
{
	int best_i = 0;
	int best_j = 0;
	double best_match = 0;
	int* sequence=m_sequence();
	double* target_signal = generate_signal_by_mpc(*target_mpc ,sequence);

	for (int i = 0;i < max_x;i++)
	{
		for (int j = 0;j < max_y;j++) 
		{
			double* current_signal = generate_signal_by_mpc(fingerprints[i][j], sequence);
			double current_match = 0;
			double current_power = 0;
			for (int k = 0;k <RECEIVE_LENTH ;k++)
			{
				current_power += pow(current_signal[k], 2);
				current_match += current_signal[k] * target_signal[k];
			}

			current_match = current_match / pow(current_power,0.5);
			if (current_match >= best_match)
			{
				best_match = current_match;
				best_i = i;
				best_j = j;
			}

			delete[] current_signal;
		}
	}

	just_point ans;
	ans.x = best_i;
	ans.y = best_j;
	delete[] target_signal;
	delete[] sequence;
	return ans;
}
mywall* add_tentative_wall(double x1,double y1,double x2,double y2,mywall* current_wall,double loss)
{
	mywall* new_wall = new mywall[WALLNUMBER_SPECULATE+1];
	for (int i = 0;i < WALLNUMBER_SPECULATE;i++)
	{
		new_wall[i].copy(current_wall[i]);
	}

	new_wall[WALLNUMBER_SPECULATE].changedata(x1,y1,x2,y2,loss);
	WALLNUMBER_SPECULATE++;
	delete[] current_wall;
	return new_wall;

}
just_point localize_by_last_square_method(double r1, double r2, double r3, double x1, double y1, double x2, double y2, double x3, double y3)
{
	Eigen::Matrix<double, 3, 1> R(pow(r1, 2) - pow(x1, 2) - pow(y1, 2), pow(r2, 2) - pow(x2, 2) - pow(y2, 2), pow(r3, 2) - pow(x3, 2) - pow(y3, 2));
	Eigen::Matrix<double, 3, 3> A;
	A <<-2 * x1, -2 * y1, 1,
		-2 * x2, -2 * y2, 1,
		-2 * x3, -2 * y3, 1;
	
	Eigen::Matrix<double, 3, 1> ANS;
	ANS = (A.transpose()*A).inverse()*A.transpose()*R;
	just_point p;
	p.x = ANS(0);
	p.y = ANS(1);

	return p;
}
//tracing
movement::movement()
{
	this->run_time = 0;
	this->a_x=0;
	this->a_y=0;
	this->v_x=0;
	this->v_y=0;
	this->x=0;//预测x
	this->y=0;//预测y
	this->x1=0;//上一个定位 定位与预测不一样
	this->y1=0;//
	this->x2=0;//上上个定位
	this->y2=0;
}
void movement::update(double x,double y)
{
	this->run_time++;

	if (this->run_time == 1)
	{
		this->x2 = x;
		this->y2 = y;
		return;
	}
	else if (this->run_time == 2)
	{
		this->x1 = x;
		this->y1 = y;
		this->v_x = x1 - x2;
		this->v_y = y1 - y2;
		this->a_x = 0;
		this->a_y = 0;
		return;
	}

	this->v_x = x - x1;
	this->v_y = y - y1;
	this->a_x = x - 2 * this->x1 + this->x2;
	this->a_y = y - 2 * this->y1 + this->y2;
	this->x2 = this->x1;
	this->y2 = this->y1;
	this->x1 = x;
	this->y1 = y;
}
void movement::predict() 
{
	if (run_time == 1)
	{
		cout<<"初始化低于2次不能进行预测！"<<endl;
	}
	//之前update时已经将x、x1、x2移位了
	this->x = this->x1 + this->v_x + this->a_x ;
	this->y = this->y1 + this->v_y + this->a_y;
}
double** movement::output_predict()
{
	double deta = 100;//方差

	double** pdf = new double*[max_x];
	for (int i = 0;i < max_x;i++)
	{
		pdf[i] = new double[max_y];
	}

	for (double i = 0;i < max_x;i++)
	{
		for (double j = 0;j < max_y;j++)
		{
			pdf[(int)i][(int)j] = exp(-((pow(i-this->x, 2)+ pow(j - this->y, 2))/deta));
		}
	}
	return pdf;
}
just_point trace_by_fingerprints(TOA_MPC** fingerprints, TOA_MPC* target_mpc,double** pdf)
{
	int best_i = 0;
	int best_j = 0;
	double best_match = 0;
	int* sequence = m_sequence();
	double* target_signal = generate_signal_by_mpc(*target_mpc, sequence);

	for (int i = 0;i < max_x;i++)
	{
		for (int j = 0;j < max_y;j++)
		{
			double* current_signal = generate_signal_by_mpc(fingerprints[i][j], sequence);
			double current_match = 0;
			double current_power = 0;
			for (int k = 0;k < RECEIVE_LENTH;k++)
			{
				current_power += pow(current_signal[k], 2);
				current_match += current_signal[k] * target_signal[k];
			}

			current_match = current_match / pow(current_power, 0.5)*pdf[i][j];
			if (current_match > best_match)
			{
				best_match = current_match;
				best_i = i;
				best_j = j;
			}

			delete[] current_signal;
		}
	}

	just_point ans;
	ans.x = best_i;
	ans.y = best_j;
	delete[] target_signal;
	delete[] sequence;
	return ans;
}
double compare_signal(double* s1, double* s2)
{
	double neiji = 0;
	double p1 = 0;
	double p2 = 0;
	for (int i = 0;i < RECEIVE_LENTH;i++)
	{
		neiji += s1[i] * s2[i];
		p1 += s1[i] * s1[i];
		p2 += s2[i] * s2[i];
	}
	/*
	cout << "*****" << endl;
	cout<< neiji <<endl;
	cout << p1 << endl;
	cout << p2 << endl;
	cout << "......." << endl;
   */
	if (p1 == 0 or p2 == 0)
	{
		return 0;
	}
	return abs(neiji) / pow(p1*p2,0.5);
}
double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}