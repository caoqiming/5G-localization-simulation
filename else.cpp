/*显示mpc结果
for (int i = 0;i < mpc->number;i++)
{
	mpc_temp = mpc_temp->next;
	cout << "分离出的第" << i << "个mpc衰减系数为" << mpc_temp->amplitude << endl;
	cout << "分离出的第" << i << "个mpc延迟为" << mpc_temp->delay << endl;
}
*/

/*部分源

	if (checkblock3(target, anchor, wall))
		anchor.can_reach_target = 1;//有line of sight
	build_the_tree(wall, &target, &anchor);
	mysignal* signal_temp = new mysignal;//第一个signal是空的，不算进去
	mysignal::signal_number++;
	mysignal* signal = signal_temp;
	mypoint::ergodic_the_tree(&anchor, &signal_temp, target);
	int* sequence=m_sequence();
	double* simulation_signal= mypoint::generate_signal(sequence, signal);




TOA_MPC* mpc = seperate_mpc(simulation_signal, sequence);


//CHECK_LOS(simulation_signal);
//cout<<mysignal::signal_number<<endl;

release_toa_mpc(mpc);


delete[] simulation_signal;
delete[] sequence;
mypoint::release_the_tree(&anchor);
mypoint::release_the_signal(signal);
*/

/*  指纹定位
//定义墙
mywall wall[4];
wall[0] = mywall(0, 0, 20, 0, 5);
wall[1] = mywall(20, 0, 20, 10, 5);
wall[2] = mywall(0, 0, 0, 10, 5);
wall[3] = mywall(0, 10, 20, 10, 5);
//定义target和anchor
mypoint target(0, 0, 5, 4);
mypoint anchor(0, 0, 1, 1);

TOA_MPC* simulate_seperate_mpc = simulate_once(&anchor, &target, wall);

TOA_MPC** fingerprints = get_fingerprints(&anchor, wall);
localization_by_fingerprints(fingerprints, simulate_seperate_mpc);



//释放第定位所用mpc
release_toa_mpc(simulate_seperate_mpc);
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


system("pause");

*/