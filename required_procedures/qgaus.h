//qgaus.h -- routines for gaussian quadrature

template <class T>
double qgaus(T &func, const double a, const double b)
{
  static const double x[] = {0.1488743389816312,0.4333953941292472,
			     0.6794095682990244,0.8650633666889845,
			     0.9739065285171717};
  static const double w[] = {0.2955242247147529,0.2692667193099963,
			     0.2190863625159821,0.1494513491505806,
			     0.0666713443086881};
  double xm = 0.5*(b+a);
  double xr = 0.5*(b-a);
  double s = 0;
  for (int j = 0; j < 5; j++)
    {
      double dx = xr * x[j];
      s += w[j]*(func(xm - dx) + func(xm + dx));
    }
  return s *= xr;
};

