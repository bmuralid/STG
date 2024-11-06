#include "sem.h"
#include "stg_inputs.h"
#include "stg_utils.h"
#include <fstream>
#include <omp.h>
#include <vector>

#define THREAD_NUM 8

template class SEM<double>;

namespace Random {
std::mt19937_64 mt{std::random_device{}()};
template <typename T> void getrand(std::vector<T> &arr) {
  std::uniform_real_distribution<T> unif{0, 1};
  for (std::size_t i = 0; i < arr.size(); ++i)
    arr[i] = unif(mt);
}
} // namespace Random

template <typename T> SEM<T>::SEM(STGInterface<T> &stg_int) {
  T sigma_max = stg_int.get_sigma_max();
  lx_sem = STGInputs::lx + 2.0 * sigma_max;
  ly_sem = STGInputs::ly + 2.0 * sigma_max;
  lz_sem = STGInputs::lz + 2.0 * sigma_max;

  auto vol_sem = lx_sem * ly_sem * lz_sem;
  nsem = int(vol_sem * stg_int.get_sigma_inv3());

  // Assign the vectors
  xsem.assign(nsem, 0.0);
  ysem.assign(nsem, 0.0);
  zsem.assign(nsem, 0.0);

  eps_sem_x.assign(nsem, 0.0);
  eps_sem_y.assign(nsem, 0.0);
  eps_sem_z.assign(nsem, 0.0);

  if (STGInputs::lGenerateSEM)
    generate_sem(stg_int);
  else
    read_sem();
}

template <typename T> void SEM<T>::generate_sem(STGInterface<T> &stg_int) {
  sigma_max = stg_int.get_sigma_max();

  Random::getrand(xsem);
  Random::getrand(ysem);
  Random::getrand(zsem);

  Random::getrand(eps_sem_x);
  Random::getrand(eps_sem_y);
  Random::getrand(eps_sem_z);

  for (auto i = 0; i < nsem; ++i) {
    xsem[i] = -sigma_max + xsem[i] * lx_sem;
    ysem[i] = -sigma_max + ysem[i] * ly_sem;
    zsem[i] = -sigma_max + zsem[i] * lz_sem;

    eps_sem_x[i] = -1.0 + eps_sem_x[i] * 2.0;
    eps_sem_y[i] = -1.0 + eps_sem_y[i] * 2.0;
    eps_sem_z[i] = -1.0 + eps_sem_z[i] * 2.0;
  }
  write_sem();
}

template <typename T> void SEM<T>::write_sem() {
  std::string fname = "sem.dat";
  std::fstream f;
  f.open(fname, std::ios::out);
  if (f.is_open()) {
    f << nsem << "\n";
    f.precision(6);
    for (auto i = 0; i < nsem; ++i) {
      f << std::scientific << xsem[i] << "\t" << ysem[i] << "\t" << zsem[i]
        << "\t" << eps_sem_x[i] << "\t" << eps_sem_y[i] << "\t" << eps_sem_z[i]
        << "\n";
    }

    f.close();
  }
}

template <typename T> void SEM<T>::read_sem() {
  std::string fname = "sem.dat";
  std::fstream f;
  f.open(fname, std::ios::in);
  if (f.is_open()) {
    f >> nsem;
    for (auto i = 0; i < nsem; ++i)
      f >> std::scientific >> xsem[i] >> ysem[i] >> zsem[i] >> eps_sem_x[i] >>
          eps_sem_y[i] >> eps_sem_z[i];
    f.close();
  }
}

template <typename T> void SEM<T>::update_sem(const T dt) {

  std::vector<T> rnd;
  rnd.assign(5, 0.0);

  for (auto i = 0; i < nsem; ++i) {
    xsem[i] += dt * STGInputs::ub0;
    if (xsem[i] > sigma_max) {
      Random::getrand(rnd);
      xsem[i] = -sigma_max;
      ysem[i] = -sigma_max + rnd[0] * ly_sem;
      zsem[i] = -sigma_max + rnd[1] * lz_sem;

      if (rnd[2] < 0.5)
        eps_sem_x[i] = -1.0;
      else
        eps_sem_x[i] = 1.0;

      if (rnd[3] < 0.5)
        eps_sem_y[i] = -1.0;
      else
        eps_sem_y[i] = 1.0;

      if (rnd[4] < 0.5)
        eps_sem_z[i] = -1.0;
      else
        eps_sem_z[i] = 1.0;
    }
  }
}

template <typename T>
void SEM<T>::get_sem_fluc(const std::vector<T> &xx, const T tke, const T eps,
                          std::vector<T> &uvecSEM) {
  std::vector<T> Rij, aij;

  Rij.assign(6, 0.0);
  aij.assign(9, 0.0);

  // stg_int.get_turb_prop(xx[0], xx[1], xx[2], tke, eps);

  STGUtils::get_re_stress(tke, Rij);

  STGUtils::re_cholesky_decomp(Rij, aij);

  T sigma = STGUtils::get_sigma(tke, eps, STGInputs::delta0, xx[3]);
  T uSEM = 0.0;
  T vSEM = 0.0;
  T wSEM = 0.0;

// SEM calculation loop
#pragma omp parallel for reduction(+ : uSEM) reduction(+ : vSEM)               \
    reduction(+ : wSEM)
  for (auto k = 0; k < nsem; ++k) {
    auto f1 = STGUtils::f_1d((xx[0] - xsem[k]) / sigma);
    auto f2 = STGUtils::f_1d((xx[1] - ysem[k]) / sigma);
    auto f3 = STGUtils::f_1d((xx[2] - zsem[k]) / sigma);

    auto s1 =
        (aij[0] * eps_sem_x[k] + aij[1] * eps_sem_y[k] + aij[2] * eps_sem_z[k]);
    auto s2 =
        (aij[3] * eps_sem_x[k] + aij[4] * eps_sem_y[k] + aij[5] * eps_sem_z[k]);
    auto s3 =
        (aij[6] * eps_sem_x[k] + aij[7] * eps_sem_y[k] + aij[8] * eps_sem_z[k]);

    uSEM = uSEM + f1 * f2 * f3 * s1;
    vSEM = vSEM + f1 * f2 * f3 * s2;
    wSEM = wSEM + f1 * f2 * f3 * s3;
  }
  uvecSEM[0] = uSEM *
               std::sqrt((lx_sem * ly_sem * lz_sem) / std::pow(sigma, 3)) /
               std::sqrt(nsem * 1.0);
  uvecSEM[1] = vSEM *
               std::sqrt((lx_sem * ly_sem * lz_sem) / std::pow(sigma, 3)) /
               std::sqrt(nsem * 1.0);
  uvecSEM[2] = wSEM *
               std::sqrt((lx_sem * ly_sem * lz_sem) / std::pow(sigma, 3)) /
               std::sqrt(nsem * 1.0);
}

template <typename T> void SEM<T>::info() {
  std::cout << "\n"
            << "SEM information"
            << "\n";
  std::cout << "-----------------------------";
  std::cout << "No. of SEM pts.:"
            << "\t" << nsem << "\n";
}

template <typename T> SEM<T>::~SEM() { write_sem(); }
