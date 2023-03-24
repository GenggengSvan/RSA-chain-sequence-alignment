#include<iostream>
#include<algorithm>
#include<string>
#include <iomanip> 
using namespace std;

//定义两组RNA最大条数为300，这个值可以更改
const int max_num = 300;

//三个数寻找最大值
int Max_three(int a, int b, int c) {
	int max = a;
	if (b > max) max = b;
	if (c > max) max = c;
	return max;
}

//大小写忽略函数
int ignore(char a) {
	int result = 0;
	if ((int)a > 95) result = (int)a - 32;
	else result = (int)a;
	return a;
}

//碱基比较赋分函数
int judge(char a, char b) {
	if (a == b) return 5;
	else if (a != b) return -4;
}

//两条RNA进行序列比对函数
int Sequence_Alignment(string s1, string s2) {
	int k = -3;
	int num_s1 = s1.length() + 1;
	int num_s2 = s2.length() + 1;
	int A[2][100000];
	int line = 1;
	for (int i = 0;i < num_s2;i++) {
		A[0][i] = i * k;
	}
	A[1][0] = (line++) * k;
	for (int i = 1;i < num_s1;i++) {
		for (int j = 1;j < num_s2;j++) {
			A[1][j] = Max_three((A[0][j - 1] + judge(s1[i - 1], s2[j - 1])), (A[0][j] + k), (A[1][j - 1] + k));
		}
		for (int j = 0;j < num_s2;j++) {
			A[0][j] = A[1][j];
			A[1][j] = 0;
		}
		A[1][0] = (line++) * k;
	}
	return A[0][num_s2 - 1];
}

//寻找增广路径函数
int findpath(int x, bool* visx, bool* visy, int* la, int* lb, int** w, int* match, int n, int* slack,int *fa)
{
	int tempDelta;
	visx[x] = true;
	for (int y = 0;y < n;y++) {
		if (visy[y])continue;
		tempDelta = la[x] + lb[y] - w[x][y];
		if (tempDelta == 0) {
			visy[y] = true;
			fa[y + n] = x;
			if (match[y] == -1) {
				return y + n;
			}
			//记录交替树的父节点信息（为了区别X，Y集合，Y的点都映射成n+y）
			fa[match[y]] =n+y;
			int res = findpath(match[y],visx,visy,la,lb,w,match,n,slack,fa);
			if (res > 0)return res;//返回增广路径的末端叶子节点
		}
		else if (slack[x] > tempDelta)//统计以x为准的slack值。
			slack[x] = tempDelta;
	}
	return -1;
}
//最大流寻找最佳匹配
int max_match(int** w,int n)
{
	int la[max_num], lb[max_num];//标杆数组
	for (int i = 0; i < n; i++) {
		la[i] = -1000;
		lb[i] = 0;
		for (int j = 0; j < n; j++)
			la[i] = max(la[i], w[i][j]);
	}
	bool visx[max_num], visy[max_num];//左右访问标记
	int delta;
	int fa[max_num], slack[max_num];//fa统计match的答案
	int match[max_num];
	memset(match, -1, sizeof(match));
	for (int x = 0; x < n; ++x) {
		for (int i = 0; i < n; ++i) slack[i] = 1<<30;
		for (int i = 0;i < n + n;i++)fa[i] = -1;
		memset(visx, false, sizeof(visx));
		memset(visy, false, sizeof(visy));
		int fir = 1;int leaf = -1;
		while (true) {
			if (fir == 1) {
				leaf = findpath(x, visx, visy, la, lb, w, match, n, slack, fa);
				fir = 0;
			}
			else {
				for (int i = 0;i <n;i++) {
					if (slack[i] == 0) {
						slack[i] = 1<<30;//slack设为无穷小
						leaf = findpath(i,visx, visy, la, lb, w, match, n, slack, fa);
						if (leaf > 0)break;
					}
				}
			}
			if (leaf > 0) {
				int p = leaf;
				while (p > 0) {
					match[p - n] = fa[p];
					p = fa[fa[p]];//回溯
				}
				break;
			}
			else {
				int delta = 1<<30;
				for (int i = 0; i < n; ++i)
					if (visx[i] && delta > slack[i])
						delta = slack[i];
				for (int i = 0; i < n; ++i)
					if (visx[i]) { la[i] -= delta;slack[i] -= delta; }//X点的slack要响应改变，slack变0说明有新边加入
				for (int j = 0; j <n; ++j) {
					if (visy[j])
						lb[j] += delta;
				}
			}
		}
	}

	//得到答案
	int ans = 0;
	for (int i = 0; i < n; i++)
		if(w[match[i]][i]!=-1000) //如果是-1000,相当于这里其实并没有RNA匹配
			ans += w[match[i]][i];
	return ans;
}

//贪心算法寻找最佳匹配
int match_tanxin(int num_1,int num_2,string* p_1,string* p_2){
	int max_all = -500;
	int* p_2_judge = new int[num_2];
	for (int i = 0;i < num_2;i++) {
		p_2_judge[i] = 0;
	}
	for (int t = 0;t < num_2;t++) {
		int max_all_temp = 0;
		for (int i = 0;i < num_2;i++) {
			p_2_judge[i] = 0;
		}
		for (int i = 0;i < num_1;i++) {
			int max_temp = -500;
			int index = 0;
			for (int j = t;j < num_2;j++) {
				if (p_2_judge[j]== 0) {
					int temp = Sequence_Alignment(p_1[i],p_2[j]);
					if (temp > max_temp) {
						max_temp = temp;
						index = j;
					}
				}
			}
			for (int j = 0;j < t;j++) {
				if (p_2_judge[j] == 0) {
					int temp = Sequence_Alignment(p_1[i],p_2[j]);
					if (temp > max_temp) {
						max_temp = temp;
						index = j;
					}
				}
			}
			p_2_judge[index]  = 1;
			max_all_temp += max_temp;
		}
		if (max_all_temp > max_all) {
			max_all = max_all_temp;

		}
	}
	return max_all;
}

int main() {
	int num_1, num_2;
	cin >> num_1 >> num_2;
	string* p_1 = new string[num_1];
	for (int i = 0;i < num_1;i++) {
		cin >> p_1[i];
	}
	string* p_2 = new string[num_2];
	for (int i = 0;i < num_2;i++) {
		cin>> p_2[i];
	}
	int n = max(num_1, num_2);
	int** allscores = new int*[n];
	for (int i = 0;i < n;i++) {
		allscores[i] = new int[n];
		for (int j = 0;j < n;j++) {
			allscores[i][j] = -1000;
		}
	}
	for (int i = 0;i < num_1;i++) {
		for (int j = 0;j < num_2;j++) {
			allscores[i][j] = Sequence_Alignment(p_1[i], p_2[j]);
		}
	}
	//最大流计算最佳匹配
	int result = max_match(allscores,n);
	//如果使用贪心算法，注释上一行，使用本行
	//int result = match_tanxin(num_1, num_2, p_1, p_2);
	cout << result<<endl;


	for (int i = 0;i < num_1;i++) {
		delete[]allscores[i];
	}
	delete[]allscores;
	delete[]p_1;
	delete[]p_2;
	return 0;
}