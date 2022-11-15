#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

int main()
{
    const char * fn = "F:\\C++Code\\RiemannSolver\\ExactRS\\ExactRSExactRS-density.plt";
    ifstream fin;
    ofstream fout;
    fout.open("F:\\C++Code\\LagrangeDG\\sod_shock_tube_polar\\output\\rho1dpolar.plt", ios::trunc);
    fin.open(fn, ios::in);
    if (!fin.is_open())
	{
		cout << "无法找到这个文件！" << endl;
	}
    double X[80], rho[80];
    for (int i=0; i<80; i++)
    {
        X[i] = 0;
        rho[i] = 0;
    }
    double temp;
    char title[20];
    fin.getline(title,20);
    fout<<title<<endl;
    cout<<title<<endl;
    int t = 0;
double x;
    while(fin>>temp)
    {
        if (t%2 == 0)
        {
            x = temp;
            fout<<temp<<"\t";
        }
        else{
            if (x<0.241801)
            {
                fout<<temp<<endl;
            }
            else if (x<0.46612)
            {
                fout<<x*temp<<endl;
            }
            else if (x < 0.668239)
            {
                fout<<x*temp<<endl;
            }
            else if (x < 0.868895)
            {
                fout<<temp*pow(x,3.0/5.0)<<endl;
            }
            else{
                fout<<temp;
            }
        }
        t++;
    }
    fin.close();
    fout.close();
    system("pause");
}