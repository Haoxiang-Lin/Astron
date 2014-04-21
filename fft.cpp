#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

const long double pi = 3.14159265358979323846264338327950;

class Fourier 
{
    public: 
        //构造函数, 输入时域数组头尾
        Fourier(double *begin, double *end);
        //对位序颠倒的时域数组进行Fourier变换. +1 正变换 ; -1 逆变换
        void Transform(const int isign);
        //位序颠倒操作
        void SwapData();
        //打印数据
        void PrintData();
    private:
        vector<double> data;
        //n是数组长度, nn是点数, nn是n的一半
        int n, nn;
};

Fourier::Fourier(double *begin, double *end)
{
    while (begin != end) {
       data.push_back(*begin++);
    }
    nn = data.size()/2;
    n = nn << 1;
}

void Fourier::Transform(const int isign)
{
    //循环相关
    int i, j, m, mmax, istep;
    //角度相关
    double theta, wr, wi, wpr, wpi;
    //临时变量
    double wtemp, tempr, tempi;

    mmax = 2;
    //长度循环: 合并长度的循环, 表示合并操作的长度为mmax
	while (mmax < n) {
		istep = mmax << 1;
		theta = - isign * (2*pi/mmax);
        //用于递增角度
        wtemp = sin(0.5*theta);
		wpr   = -2.0 * wtemp * wtemp;
		wpi   = sin(theta);
        //初始角度为0
		wr    = 1.0;
		wi    = 0.0;
        //角度循环: 循环角度
		for (m = 1; m < mmax; m += 2) {
            //取样循环: 固定角度, 取样时域, 进行奇偶合并操作
			for (i = m; i <= n; i += istep) {
                //生成i和j的取样数列
				j = i + mmax;

                //奇取样点的操作
				tempr     =  wr * data[j-1] - wi * data[j];
				tempi     =  wr * data[j]   + wi * data[j-1];

                //共轭合并 (巧妙地利用了单位根的周期性, 
                 //结合角度循环直接算出所有系数, 无需重头开始计算)
				data[j-1] =  data[i-1] - tempr;
				data[j]   =  data[i]   - tempi;

                //合并
				data[i-1] =  data[i-1] + tempr;
				data[i]   =  data[i] + tempi;

                //归一化
                if (isign == -1 && mmax != n/2) {
                    data[j-1] /= 2.0;
                    data[j]   /= 2.0;
                    data[i-1] /= 2.0;
                    data[i]   /= 2.0;
                }
			}
            //递增角度
            wtemp = wr;
			wr = wtemp * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
        //长度翻倍
		mmax <<= 1;
    }
}


void Fourier::SwapData()
{
    int i, j, m;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (i < j) {
			swap(data[j-1], data[i-1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
}

void Fourier::PrintData()
{
    string file_name = "log.txt";
    ofstream outdata;
    outdata.open(file_name.c_str(), ios::out);
    outdata << fixed << showpoint << right << setprecision(6);

    for (vector<double>::iterator iter = data.begin(); iter != data.end(); iter += 2) {
        outdata << x << setw(30) << *iter << endl;
    }
}

