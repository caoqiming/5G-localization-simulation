#pragma once
#include"mypoint.h"
#include"mywall.h"
const int RECEIVE_LENTH = 1500;
const int MSEQUENCE_LENTH = 1023;
const double FREQUENCE = 3500;//3.5GHz
const double HIGHESTDB = 100;//能分辨的最高损耗
extern int WALLNUMBER;//修改的时候一定要把mypoint.h里的*children也改了
extern int WALLNUMBER_SAVE;
extern int WALLNUMBER_SPECULATE;
extern double NOIZE_BZC;//噪声标准差
const int max_x = 100;//地图边界大小
const int max_y = 100;

struct just_point 
{
	double x;
	double y;
	just_point* next = nullptr;
	just_point* front = nullptr;
};
struct TOA_MPC 
{
	static int number;
	double amplitude;
	double delay;//取整数以表示采样点延迟
	TOA_MPC* next = nullptr;
};
struct TRAN_LOS_DATA {
	double mean_excess_delay;
	double RMS_delay_spread;
	double Amplitude_kurtosis;
	double power;
	double amplitude_max;
};
class movement 
{//运动模型预测 
public:
	int run_time;//2次更新之后开始预测，3次更新之后才有加速度
	double a_x;
	double a_y;
	double v_x;
	double v_y;
	double x;//预测x
	double y;//预测y
	double x1;//上一个定位 定位与预测不一样
	double y1;//
	double x2;//上上个定位
	double y2;

	movement();
	void update(double x, double y);//这里的x y 是上一次最终预测结果,更新速度与加速度信息
	void predict();//预测信息储存在x y 中
	double **output_predict();
};

TOA_MPC* seperate_mpc(double* simulation_signal,int* sequence);
TOA_MPC* get_one_mpc(double simulation_signal[RECEIVE_LENTH], int* sequence);
mysignal* mpc_to_mysignal(TOA_MPC* p);
void print_mysignal(mysignal* s);
void release_toa_mpc(TOA_MPC* p);
void release_toa_mpc(TOA_MPC* p, int depth);
TRAN_LOS_DATA GET_LOS_DATA(double* simulation_signal);
int* m_sequence();
mypoint* get_symmetric_point(mypoint* p, mywall ab);
bool checkblock1(mypoint p, mypoint q, mywall ab);
bool checkblock2(mypoint p, mypoint q, mywall wall[]);
bool checkblock3(mypoint p, mypoint q, mywall wall[]);
just_point* get_intersection(just_point p1, just_point p2, just_point p3, just_point p4);
bool super_checkblock(mypoint* p, mypoint* target, mywall* wall);
bool check_loss(double loss_fanshe, double loss_juli);
void build_the_tree(mywall wall[], mypoint* target, mypoint* anchor);
TOA_MPC** get_fingerprints(mypoint* anchor, mywall* wall);
double* simulate_once_receive(mypoint* anchor, mypoint* target, mywall* wall);
TOA_MPC* simulate_once(mypoint* anchor, mypoint* target, mywall* wall);
just_point localization_by_fingerprints(TOA_MPC** fingerprints, TOA_MPC* target_mpc);
double* generate_signal_by_mpc(TOA_MPC mpc, int* sequence);
mywall* add_tentative_wall(double x1, double y1, double x2, double y2, mywall* current_wall, double loss = 5);
just_point localize_by_last_square_method(double r1, double r2, double r3, double x1, double y1, double x2, double y2, double x3, double y3);
just_point trace_by_fingerprints(TOA_MPC** fingerprints, TOA_MPC* target_mpc, double** pdf);
double compare_signal(double* s1, double* s2);
double gaussrand();

