#ifndef _STG_INTERFACE_H
#define _STG_INTERFACE_H

#include<string>
#include<vector> 
#include <bits/stdc++.h>

template <typename T>
class STGInterface{
   int npts;
   std::string turb_model;
   std::vector<T> x, y, z;
   std::vector<T> deltag;
   std::vector<T> tke_int, eps_int;
   std::vector<T> sigma;
public:
   STGInterface(const std::string);
   void info();
   T get_sigma_max();
   T get_sigma_inv3();
   T get_maxdel(const int);
   int get_npts(){return npts;}
   void get_cords(const int, T &, T &, T &);
   void get_turb_prop(const int, T &, T &);
};

#endif // _STG_INTERFACE_H

