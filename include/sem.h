#ifndef _SEM_H
#define _SEM_H
#include "stg_interface.h"

template <typename T>
struct STGgpuData{
   T *xsem, *ysem, *zsem;
   T *eps_sem_x, *eps_sem_y, *eps_sem_z;
   T *uvec;
   T sigma;
   T *xx, *a_ij;
};

template <typename T>
class SEM{
    int nsem;
    STGgpuData<T> _gpu;
    T delta0, ub0, sigma_max; 
    T lx_sem, ly_sem, lz_sem;
    std::vector<T>  xsem, ysem, zsem;
    std::vector<T>  eps_sem_x, eps_sem_y, eps_sem_z;
public:
    SEM(STGInterface<T> &);
    ~SEM();
    void info();
    void generate_sem(STGInterface<T> &);
    void write_sem();
    void read_sem();
    void update_sem(const T);
    void get_sem_fluc (const std::vector<T> &,const  T, const T, std::vector<T> &);
};

#endif // _SEM_H
