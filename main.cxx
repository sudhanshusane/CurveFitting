#include <iostream>
#include "spline.h"
#include "cubichermite.h"
#include <vector>

using namespace std;

int main()
{
vector<double> time;
vector<double> x,y,z;

time.push_back(0);
time.push_back(5);
time.push_back(10);

x.push_back(2);
x.push_back(1.8);
x.push_back(0);

y.push_back(0);
y.push_back(1.8);
y.push_back(2);

z.push_back(0);
z.push_back(0);
z.push_back(0);

/*
tk::spline s_x;
s_x.set_points(time, x);
tk::spline s_y;
s_y.set_points(time, y);

for(int j = 0; j <= 20; j++)
{
  double t = j*1.0;
  cout << s_x(t) << "," << s_y(t) << ",0" <<endl;
}
*/

curve::CubicHermite spline;
spline.SetPoints(time,x,y,z);

for(int j = 0; j <= 10; j++)
{
  double t = j*1.0;
  double pt[3];
  spline.GetLocationAtTime(t, pt);
  std::cout << pt[0] << "," << pt[1] << "," << pt[2] << std::endl;  
}


}
