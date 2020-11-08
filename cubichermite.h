#ifndef CURVE_CUBICHERMITE_H
#define CURVE_CUBICHERMITE_H

#include <iostream>
#include <vector>

namespace
{

namespace curve
{

class CubicHermite
{
public:
CubicHermite()
{
  tension = -1.0;
  bias = 0.0;
}

void SetTension(double t)
{
  tension = t;
}

void SetBias(double b)
{
  bias = b;
}

void SetPoints(std::vector<double> in_time, std::vector<double> x1, std::vector<double> y1, std::vector<double> z1)
{
  if(x1.size() == in_time.size() && x1.size() == y1.size() 
    && x1.size() == z1.size() && x1.size() >= 2)
  {
    time.assign(in_time.begin(), in_time.end());
    x.assign(x1.begin(), x1.end());
    y.assign(y1.begin(), y1.end());
    z.assign(z1.begin(), z1.end());
  }
  else
  {
    std::cout << "Error: Vector sizes don't align in SetPoints()" << std::endl;
  }
}

void GetLocationAtTime(double t, double *pt)
{
  if(time.size() >= 2)
  {
    if(t >= time.front() && t <= time.back())
    { 
        int index[4] = {0};
        // index[0], index[3] can be set afterwards. Calculate index 1,2 based on where t lies.
        for(int i = 0; i < time.size()-1; i++)
        {
          if(t >= time[i] && t <= time[i+1])
          {
            index[1] = i;
            index[2] = i+1;
            break;
          }
        }
        if(index[1] == 0)
          index[0] = index[1];
        else
          index[0] = index[1] - 1;
        
        if(index[2] == time.size()-1)
          index[3] = index[2];
        else
          index[3] = index[2]+1;  

        double mu = (t - time[index[1]])/(time[index[2]] - time[index[1]]);  
        pt[0] = HermiteInterpolate(x[index[0]], x[index[1]], x[index[2]], x[index[3]], mu);
        pt[1] = HermiteInterpolate(y[index[0]], y[index[1]], y[index[2]], y[index[3]], mu);
        pt[2] = HermiteInterpolate(z[index[0]], z[index[1]], z[index[2]], z[index[3]], mu);
      }
    else
    {
      std::cout << "Requested time outside range. Extrapolation not supported." << std::endl;
    }
  }
}

private:

double HermiteInterpolate(double y0, double y1, double y2, double y3, double mu)
{
  double m0,m1,mu2,mu3;
  double a0,a1,a2,a3;

  mu2 = mu * mu;
  mu3 = mu2 * mu;
  m0  = (y1-y0)*(1+bias)*(1-tension)/2;
  m0 += (y2-y1)*(1-bias)*(1-tension)/2;
  m1  = (y2-y1)*(1+bias)*(1-tension)/2;
  m1 += (y3-y2)*(1-bias)*(1-tension)/2;
  a0 =  2*mu3 - 3*mu2 + 1;
  a1 =    mu3 - 2*mu2 + mu;
  a2 =    mu3 -   mu2;
  a3 = -2*mu3 + 3*mu2;

  return(a0*y1+a1*m0+a2*m1+a3*y2);
}

double tension, bias;
std::vector<double> time, x, y, z;
};

} // namespace spline

} // namespace

#endif /* SPLINE_CUBICHERMITE_H */
