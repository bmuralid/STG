#ifndef _STG_H
#define _STG_H

#include "stg_interface.h"
#include "sem.h"
#include "stg_utils.h"
#include <memory>


template <typename T>
class STG{
    int nt, npts;
    STGInterface<T> stg_int;
    SEM<T> sem;
    std::unique_ptr<STGUtils::Array2D<T>> uprime, vprime, wprime;
public:
STG();
~STG();
void update_stg();
void get_stg_fluc(const int, const T);
void finalize_stg();
};

#endif

