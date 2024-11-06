#include "stg.h"
#include "stg_inputs.h"
#include <cstdio>
#include <fstream>
#include <iostream>

template class STG<double>;

template <typename T>
STG<T>::STG() : stg_int("stg_interface_01.dat"), sem(stg_int) {
  std::cout << "In STG constructor" << std::endl;
  stg_int.info();
  npts = stg_int.get_npts();
  sem.generate_sem(stg_int);
  // sem.write_sem();
  // sem.info();
  uprime = nullptr;
  vprime = nullptr;
  wprime = nullptr;
  std::cout << "End STG constructor" << std::endl;
}

template <typename T> void STG<T>::get_stg_fluc(const int nt_, const T dt) {
  npts = stg_int.get_npts();
  nt = nt_;
  T x_, y_, z_, delm_;
  T tke_, eps_;

  std::vector<T> xx, uvec;
  xx.assign(4, 0.0);
  uvec.assign(3, 0.0);

  uprime = std::make_unique<STGUtils::Array2D<T>>(nt, npts);
  vprime = std::make_unique<STGUtils::Array2D<T>>(nt, npts);
  wprime = std::make_unique<STGUtils::Array2D<T>>(nt, npts);

  for (auto i = 0; i < nt; ++i) {
    if ((i - 1) % 10 == 0)
      printf("\nUpdating time step %d, time = %e", i, i * dt);
    for (auto j = 0; j < npts; ++j) {
      stg_int.get_cords(j, x_, y_, z_);
      delm_ = stg_int.get_maxdel(j);
      stg_int.get_turb_prop(j, tke_, eps_);

      xx[0] = x_;
      xx[1] = y_;
      xx[2] = z_;
      xx[3] = delm_;
      sem.get_sem_fluc(xx, tke_, eps_, uvec);
      uprime->arr[i][j] = uvec[0];
      vprime->arr[i][j] = uvec[1];
      wprime->arr[i][j] = uvec[2];
    }
    sem.update_sem(dt);
  }
}

template <typename T> STG<T>::~STG() {
  std::fstream f;
  f.open("uprime.dat", std::ios::out);
  if (f.is_open()) {
    for (auto i = 0; i < npts; ++i) {
      for (auto j = 0; j < nt; ++j) {
        f << std::scientific << uprime->arr[j][i] << "\t";
      }
      f << std::endl;
    }
    f.close();
  }

  f.open("vprime.dat", std::ios::out);
  if (f.is_open()) {
    for (auto i = 0; i < npts; ++i) {
      for (auto j = 0; j < nt; ++j) {
        f << std::scientific << vprime->arr[j][i] << "\t";
      }
      f << std::endl;
    }
    f.close();
  }
  f.open("wprime.dat", std::ios::out);
  if (f.is_open()) {
    for (auto i = 0; i < npts; ++i) {
      for (auto j = 0; j < nt; ++j) {
        f << std::scientific << wprime->arr[j][i] << "\t";
      }
      f << std::endl;
    }
    f.close();
  }
}

int main(int argc, char *argv[]) {
  int nt;

  if (argc > 1) {
    nt = atoi(argv[1]);
  } else {
    nt = 10;
  }

  double dt = 2.0e-6;

  if (argc > 2) {
    nt = std::atoi(argv[1]);
    dt = std::atof(argv[2]);
  }

  STGInputs::read_input();

  STG<double> stg;

  stg.get_stg_fluc(nt, dt);

  return 0;
}
