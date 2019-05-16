#include<iostream>
#include"mypoint.h"
#include"mywall.h"
#include"localization.h"
#include<Python.h>
#include <fstream>

using namespace std;
int mypoint::pointcount = 0;
int TOA_MPC::number = 0;
//m序列生成

//mian8里面用到的
double mpc_test(mywall* wall,mypoint anchor,mypoint target,int* sequence)
{
	double accuracy=0;
	target.x = ((double)rand() / RAND_MAX)*(100 - 0) + 0;
	target.y = ((double)rand() / RAND_MAX)*(100 - 0) + 0;

	build_the_tree(wall, &target, &anchor);
	mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
	mysignal* signal = signal_temp;
	mypoint::ergodic_the_tree(&anchor, &signal_temp, target);
	double* simulation_signal = mypoint::generate_signal(sequence, signal);
	TOA_MPC* mpc = seperate_mpc(simulation_signal, sequence);
	delete[] simulation_signal;
	//分析signal 与 mpc 的区别
	mysignal* mpc_ans = mpc_to_mysignal(mpc);//将分解的多径向量转格式
	double* signal_real = mypoint::generate_signal_without_noize(sequence, signal);
	double* signal_mpc = mypoint::generate_signal_without_noize(sequence, mpc_ans);
	accuracy = compare_signal(signal_mpc, signal_real);

	
	cout<<"分解得到的"<<endl;
	print_mysignal(mpc_ans->next);
	cout << "事实上的" << endl;
	print_mysignal(signal->next);
	/*
	for (int j = 0;j < RECEIVE_LENTH;j++)
	{
		cout<< signal_mpc[j]<<endl;
	}
	*/
	delete[] signal_real;
	delete[] signal_mpc;
	release_toa_mpc(mpc);
	mypoint::release_the_tree(&anchor);
	mypoint::release_the_signal(signal);
	mypoint::release_the_signal(mpc_ans);
	return accuracy;
}

int main1()//用指纹信息定位性能测试
{
	WALLNUMBER++;
	//定义墙
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	wall[4] = mywall(60, 0, 60, 30, 5);
	//定义target和anchor
	mypoint target(0, 0, 5, 4);
	mypoint anchor(0, 0, 5, 95);

	just_point ans;
	int testtime = 500;
	double accuracy = 0;
	TOA_MPC** fingerprints = get_fingerprints(&anchor, wall);
	for (int fc_change = 0;fc_change < 1;fc_change++)//改变噪声的循环
	{
		NOIZE_BZC = 8e-10;//+1e-10;
		accuracy = 0;

		for (int runtime = 0;runtime < testtime;runtime++)//求均值的循环
		{
			target.x = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
			target.y = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
			TOA_MPC* simulate_seperate_mpc = simulate_once(&anchor, &target, wall);
			ans=localization_by_fingerprints(fingerprints, simulate_seperate_mpc);
			//释放第定位所用mpc
			release_toa_mpc(simulate_seperate_mpc);
			accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
		}


		cout << "噪声标准差为 " << NOIZE_BZC << "时：" << endl;
		cout << accuracy / testtime << "," << endl;
	}
		
	//释放指纹信息
	for (int i = 0;i < max_x;i++)
	{
		for (int j = 0;j < max_y;j++)
		{
			release_toa_mpc(&fingerprints[i][j], 0);
		}
	}
	for (int i = 0; i < max_x; ++i)
	{
		delete[] fingerprints[i];
	}
	delete[] fingerprints;

	delete[] wall;
	return 0;
}
int main2()//用于svm收集数据
{
	ofstream fout("data_new_y.json");
	const int data_number=1000;//数据量
	fout << "{\"number\":"<< data_number<<",\"data\":[";
	//动态定义墙
	WALLNUMBER++;
	mywall* wall=new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	wall[4] = mywall(60, 0, 60, 30, 5);
	//定义target和anchor
	mypoint target(0, 0, 50, 40);
	mypoint anchor(0, 0, 5, 95);
	srand((int)time(NULL));
	for (int i = 0;i < data_number;i++)
	{
		target.x = ((double)rand() / RAND_MAX)*(100 - 0) + 0;
		target.y = ((double)rand() / RAND_MAX)*(100 - 0) + 0;


		
		if (checkblock3(target, anchor, wall))
			anchor.can_reach_target = 1;//有line of sight
		else
			anchor.can_reach_target = 0;//没有los



		build_the_tree(wall, &target, &anchor);
		mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
		mysignal* signal = signal_temp;
		mypoint::ergodic_the_tree(&anchor, &signal_temp, target);
		int* sequence = m_sequence();
		double* simulation_signal = mypoint::generate_signal(sequence, signal);


		TRAN_LOS_DATA train_data= GET_LOS_DATA(simulation_signal);



		fout << "[" << anchor.can_reach_target<<",";
		fout << train_data.Amplitude_kurtosis*1e28 << ",";
		fout << train_data.amplitude_max*1e9 << ",";
		fout << train_data.mean_excess_delay*1e-8 << ",";
		fout << train_data.power*1e19 << ",";
		fout << train_data.RMS_delay_spread*1e-8 << "]";
		if (i != data_number - 1)
			fout << ",";

		delete[] simulation_signal;
		delete[] sequence;
		mypoint::release_the_tree(&anchor);
		mypoint::release_the_signal(signal);
	}
	fout << "]}";
	fout.close();
	delete[] wall;
	return 0;
}
int main3()//调用python函数
{
	Py_Initialize();
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('./')");
	PyObject* pmod = PyImport_ImportModule("demo");
	PyObject* pfunc = PyObject_GetAttrString(pmod, "test");
	PyObject *pArgs = PyTuple_New(2);               //函数调用的参数传递均是以元组的形式打包的,2表示参数个数  
	PyTuple_SetItem(pArgs, 0, Py_BuildValue("d", 5.6));//0---序号  i表示创建int型变量  
	PyTuple_SetItem(pArgs, 1, Py_BuildValue("d", 7.2));//1---序号 
	//PyTuple_SetItem(pArgs, 2, Py_BuildValue("d", 7.2));
	//PyTuple_SetItem(pArgs, 3, Py_BuildValue("d", 7.2)); 
	//PyTuple_SetItem(pArgs, 4, Py_BuildValue("d", 7.2)); 
	PyObject* presult = PyObject_CallObject(pfunc, pArgs);
	int ans = _PyLong_AsInt(presult);
	cout << ans << endl;
	Py_Finalize();
	system("pause");
	return 0;
}
int main5()
{//仿真三个基站的toa定位算法
	//WALLNUMBER++;
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	//wall[4] = mywall(50, 0, 50, 30, 5);
	//定义target和anchor
	mypoint target(0, 0, 34.5, 57.9);
	mypoint anchor1(0, 0, 5, 5);
	mypoint anchor2(0, 0, 95, 5);
	mypoint anchor3(0, 0, 5, 95);
	
	int testtime = 1000;
	TOA_MPC* receive_1,* receive_2,* receive_3;
	double r1, r2, r3;
	double accuracy = 0;
	for (int fc_change = 0;fc_change < 20;fc_change++)//改变噪声的循环
	{
		NOIZE_BZC = NOIZE_BZC+1e-10;//+1e-10;
		accuracy = 0;
		
		for (int runtime = 0;runtime < testtime;runtime++)//求均值的循环
		{
			target.x = ((double)rand() / RAND_MAX)*(100 - 0) + 0;
			target.y = ((double)rand() / RAND_MAX)*(100 - 0) + 0;
			receive_1 = simulate_once(&anchor1, &target, wall);
			receive_2 = simulate_once(&anchor2, &target, wall);
			receive_3 = simulate_once(&anchor3, &target, wall);
			if (!receive_1->next or !receive_2->next or !receive_3->next)
			{
				just_point ans;
				ans.x= ((double)rand() / RAND_MAX)*(100 - 0) + 0;
				ans.y= ((double)rand() / RAND_MAX)*(100 - 0) + 0;
				accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
			}
			else
			{
				r1 = receive_1->next->delay*2.4414;
				r2 = receive_2->next->delay*2.4414;
				r3 = receive_3->next->delay*2.4414;
				just_point ans;
				ans = localize_by_last_square_method(r1, r2, r3, 5, 5, 95, 5, 5, 95);
				accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
			}

			release_toa_mpc(receive_1);
			release_toa_mpc(receive_2);
			release_toa_mpc(receive_3);
			receive_1->next = nullptr;
			receive_2->next = nullptr;
			receive_3->next = nullptr;
		}
		//cout << "噪声标准差为 " << NOIZE_BZC << "时：" << endl;
		cout << accuracy / testtime << ","<<endl;

	}






	delete[] wall;
	return 0;
}
int main4()//自动更新障碍信息
{
	//动态定义墙
	WALLNUMBER++;
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	wall[4] = mywall(50, 0, 50, 50, 5);
	//定义target和anchor
	mypoint target(0, 0, 50, 40);
	mypoint anchor(0, 0, 5, 5);

	mywall* tentative_wall = new mywall[WALLNUMBER_SPECULATE];
	for (int i = 0;i < 5;i++)//猜测墙的位置
	{

	}




	delete[] wall;
	return 0;
}
int main6() 
{//移动追踪算法
		//定义墙
	WALLNUMBER++;
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	wall[4] = mywall(60, 80, 80, 60, 5);
	//定义target和anchor
	mypoint anchor(0, 0, 5, 5);
	TOA_MPC** fingerprints = get_fingerprints(&anchor, wall);



	just_point ans;//用于向运动预测模型输入参数
	movement trace;
	mypoint target(0, 0, 87.4, 90.7);

	for (int runtime = 0;runtime < 10;runtime++)//里面更改目标的位置进行连续追踪
	{
		switch (runtime)
		{
		case 0:
			target.x = 87.4;
			break;
		case 1:
			target.x = 78.8;
			break;
		case 2:
			target.x = 70.2;
			break;
		case 3:
			target.x = 61.6;
			break;
		case 4:
			target.x = 53;
			break;
		case 5:
			target.x = 44.4;
			break;
		case 6:
			target.x = 35.8;
			break;
		case 7:
			target.x = 27.2;
			break;
		case 8:
			target.x = 18.6;
			break;
		case 9:
			target.x = 10;
			break;
		}

		TOA_MPC* simulate_seperate_mpc = simulate_once(&anchor, &target, wall);

		if (runtime < 2)
		{
			ans = localization_by_fingerprints(fingerprints, simulate_seperate_mpc);
			if (runtime < 1)
			{
				trace.update(ans.x, ans.y);
			}
		}
		else
		{
			trace.update(ans.x, ans.y);
			trace.predict();
			double** pdf = trace.output_predict();
			ans = trace_by_fingerprints(fingerprints, simulate_seperate_mpc, pdf); //trace_by_fingerprints(fingerprints, simulate_seperate_mpc, pdf);
			for (int i = 0; i < max_x; ++i)
			{
				delete[] pdf[i];
			}
			delete[] pdf;
		}
		cout << "runtime:" << runtime << endl;
		cout << ans.x << "," << ans.y << endl;
		//释放第定位所用mpc
		release_toa_mpc(simulate_seperate_mpc);
	}

	
	//释放指纹信息
	for (int i = 0;i < max_x;i++)
	{
		for (int j = 0;j < max_y;j++)
		{
			release_toa_mpc(&fingerprints[i][j], 0);
		}
	}
	for (int i = 0; i < max_x; ++i)
	{
		delete[] fingerprints[i];
	}
	delete[] fingerprints;
	TOA_MPC::number = TOA_MPC::number - 10000;
	delete[] wall;
	return 0;
}
int main7()
{
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	//定义target和anchor
	mypoint anchor(0, 0, 5, 5);
	mypoint target1(0, 0, 5, 5);
	mypoint target2(0, 0, 5, 5);
	int* sequence = m_sequence();
	TOA_MPC*ans1=simulate_once( &anchor, &target1,  wall);
	TOA_MPC*ans2 = simulate_once(&anchor, &target2, wall);

	return 0;
}
int main8()//多径分解算法性能分析 
{
	//定义墙
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	//定义target和anchor
	mypoint target(0, 0, 5, 4);
	mypoint anchor(0, 0, 5, 5);
	int* sequence = m_sequence();
	srand((int)time(NULL));
	double accuracy = 0;
	for (int fc_change = 0;fc_change < 1;fc_change++)//改变噪声的循环
	{
		NOIZE_BZC = NOIZE_BZC;//+1e-10;
		accuracy = 0;
		int testtime = 1;
		for (int runtime = 0;runtime < testtime;runtime++)//求均值的循环
		{
			accuracy+=mpc_test(wall, anchor, target, sequence);
		}
		cout << "噪声标准差为 " << NOIZE_BZC << "时：" << endl;
		cout << accuracy/ testtime << endl;

	}
	


	delete[] sequence;
	return 0;
}
int main9()//基于LOS检测的性能分析
{
	WALLNUMBER++;
	mywall* wall = new mywall[WALLNUMBER];
	wall[0] = mywall(0, 0, 100, 0, 5);
	wall[1] = mywall(100, 0, 100, 100, 5);
	wall[2] = mywall(0, 0, 0, 100, 5);
	wall[3] = mywall(0, 100, 100, 100, 5);
	wall[4] = mywall(60, 0, 60, 30, 5);

	//定义target和anchor
	mypoint target(0, 0, 5, 4);
	mypoint anchor1(0, 0, 5, 5);
	mypoint anchor2(0, 0, 95, 5);
	mypoint anchor3(0, 0, 5, 95);
	TOA_MPC* receive_1, *receive_2, *receive_3;
	double r1, r2, r3;

	just_point ans;
	int testtime = 100;
	double accuracy = 0;
	TOA_MPC** fingerprints = get_fingerprints(&anchor3, wall);
	for (int fc_change = 0;fc_change < 10;fc_change++)//改变噪声的循环
	{
		NOIZE_BZC = 1e-10;
		accuracy = 0;

		for (int runtime = 0;runtime < testtime;runtime++)//求均值的循环
		{
			target.x = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
			target.y = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
			if (checkblock3(target, anchor1, wall) and checkblock3(target, anchor2, wall) and checkblock3(target, anchor3, wall))
			{//有视距，TOA定位
				receive_1 = simulate_once(&anchor1, &target, wall);
				receive_2 = simulate_once(&anchor2, &target, wall);
				receive_3 = simulate_once(&anchor3, &target, wall);
				if (!receive_1->next or !receive_2->next or !receive_3->next)
				{
					just_point ans;
					ans.x = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
					ans.y = ((double)rand() / RAND_MAX)*(70 - 30) + 30;
					cout<<"不该出现这个的";
					accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
				}
				else
				{
					r1 = receive_1->next->delay*2.4414;
					r2 = receive_2->next->delay*2.4414;
					r3 = receive_3->next->delay*2.4414;
					just_point ans;
					ans = localize_by_last_square_method(r1, r2, r3, 5, 5, 95, 5, 5, 95);
					accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
				}

				release_toa_mpc(receive_1);
				release_toa_mpc(receive_2);
				release_toa_mpc(receive_3);
				receive_1->next = nullptr;
				receive_2->next = nullptr;
				receive_3->next = nullptr;
			}
			else
			{//无视距，指纹定位
				TOA_MPC* simulate_seperate_mpc = simulate_once(&anchor3, &target, wall);
				ans = localization_by_fingerprints(fingerprints, simulate_seperate_mpc);
				//释放第定位所用mpc
				release_toa_mpc(simulate_seperate_mpc);
				accuracy += pow(pow(ans.x - target.x, 2) + pow(ans.y - target.y, 2), 0.5);
			}


		}


		cout << "噪声标准差为 " << NOIZE_BZC << "时：" << endl;
		cout << accuracy / testtime << "," << endl;
	}

	//释放指纹信息
	for (int i = 0;i < max_x;i++)
	{
		for (int j = 0;j < max_y;j++)
		{
			release_toa_mpc(&fingerprints[i][j], 0);
		}
	}
	for (int i = 0; i < max_x; ++i)
	{
		delete[] fingerprints[i];
	}
	delete[] fingerprints;

	delete[] wall;
	return 0;



}
int main()
{
	srand((int)time(NULL));
	main9();
	system("pause");
	return 0;
}

/*
	if (checkblock3(target,anchor, wall))
		anchor.can_reach_target = 1;//有line of sight
	build_the_tree(wall, &target, &anchor);
	mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
	mysignal* signal = signal_temp;
	mypoint::ergodic_the_tree(&anchor, &signal_temp, target);


	mypoint::release_the_tree(&anchor);
	mypoint::release_the_signal(signal);
	system("pause");
*/