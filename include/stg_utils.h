#ifndef _STG_UTILS_H
#define _STG_UTILS_H
#include<vector>
#include<cmath>


namespace STGUtils 
{

template <typename T>
void get_re_stress(const T & tke, std::vector<T> & Rij)
{
Rij[0] = 4.0/9.0 * 2.0 * tke;
Rij[1] = 2.0/9.0 * 2.0 * tke; 
Rij[2] = 1.0/3.0 * 2.0 * tke; 
Rij[3] = - 0.3 * tke; 
Rij[4] = 0.0; 
Rij[5] = 0.0; 
}

template <typename T>
void re_cholesky_decomp(const std::vector<T> & Rij, std::vector<T> & aij)
{
    aij[0] = std::sqrt(std::max(Rij[0], 1.0e-10));
    aij[1] = 0.0;
    aij[2] = 0.0;
    aij[3] = Rij[3]/aij[0];
    aij[4] = std::sqrt(std::max(Rij[1]-std::pow(aij[3],2), 0.0));
    aij[5] = 0.0; 
    aij[6] = 0.0;
    aij[7] = 0.0;
    aij[8] = std::sqrt(std::max(Rij[2]- std::pow(aij[6],2)
                       - std::pow(aij[7],2), 0.0));
}

template <typename T>
T get_sigma(const T & tke, const T & eps, const T & delta0, const T& delm) 
{
    return std::max(std::min(std::pow(tke,1.5)/eps, 0.41*delta0), delm);
}

template <typename T>
T f_1d(const T  x)
{
    T _f1d; 
    if(std::abs(x)< 1.0)
    _f1d = std::sqrt(1.5) * (1.0 - std::abs(x));
    else
    _f1d = 0.0;
    return _f1d;
}

template <typename T>
struct Array2D{
    int Nx_, Ny_;
    T** arr;
    Array2D();
    Array2D(const int Nx, const int Ny)
    {
        arr = new T *[Nx];
        for(auto i=0; i<Nx; ++i)
        {
            arr[i] = new T [Ny]();
        }
        Nx_ = Nx;
        Ny_ = Ny;

        for(auto i=0; i<Nx; ++i){
            for(auto j=0; j<Ny; ++j)
            {
                arr[i][j] = 0.0;
            }
        }

    }

    ~Array2D()
    {
        for(auto i=0; i<Nx_; ++i){
            delete[] arr[i];
        }
        delete[] arr;
    }

    T* operator [](const int & i)
    {
        return arr[i];
    }

    // copy constructor 
    Array2D(const Array2D & rhs)
    {
        Nx_ = rhs.Nx_;
        Ny_ = rhs.Ny_;
        arr = new T *[Nx_];
        for(auto i =0; i<Nx_; ++i)
        {
            arr[i] = new T [Ny_]();
            for(auto j=0; j<Ny_; ++j)
            {
                arr[i][j] = rhs.arr[i][j];
            }
        }
    }

    // move constructor
    Array2D(const Array2D &&rhs) noexcept
    {
        Nx_ = rhs.Nx_;
        Ny_ = rhs.Ny_;
        arr = rhs.arr;
        rhs.arr = nullptr;
    }
};
}
#endif // _STG_UTILS_H
