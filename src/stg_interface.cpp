#include "stg_interface.h"
#include "stg_inputs.h"
#include <fstream>
#include <iostream>

template class STGInterface<double>;

template <typename T> STGInterface<T>::STGInterface(const std::string fname) {
  std::fstream f;
  std::string str;
  f.open(fname, std::ios::in);
  if (f.is_open()) {
    f >> npts >> str >> turb_model;

    x.assign(npts, 0.0);
    y.assign(npts, 0.0);
    z.assign(npts, 0.0);
    deltag.assign(npts, 0.0);
    tke_int.assign(npts, 0.0);
    eps_int.assign(npts, 0.0);
    sigma.assign(npts, 0.0);

    for (auto i = 0; i < npts; ++i) {
      f >> std::scientific >> x[i] >> y[i] >> z[i] >> deltag[i] >> tke_int[i] >>
          eps_int[i];
    }

    f.close();
  } else {
    printf("Unable to open file");
  }

  for (auto i = 0; i < npts; ++i) {
    sigma[i] = std::max(std::min(std::pow(tke_int[i], 1.5) / eps_int[i],
                                 0.41 * STGInputs::delta0),
                        deltag[i]);
  }
}

template <typename T>
void STGInterface<T>::get_cords(const int k, T &x_, T &y_, T &z_) {
  x_ = x[k];
  y_ = y[k];
  z_ = z[k];
}

template <typename T> T STGInterface<T>::get_maxdel(const int k) {
  return deltag[k];
}

template <typename T>
void STGInterface<T>::get_turb_prop(const int k, T &tke_, T &eps_) {
  tke_ = tke_int[k];
  eps_ = eps_int[k];
}

template <typename T> T get_max(const std::vector<T> &arr) {
  T maxval = *std::max_element(arr.begin(), arr.end());
  return maxval;
}

template <typename T> T get_min(const std::vector<T> &arr) {
  T minval = *std::min_element(arr.begin(), arr.end());
  return minval;
}

template <typename T> void STGInterface<T>::info() {
  auto xmin = get_min<T>(x);
  auto xmax = get_max<T>(x);

  auto ymin = get_min<T>(y);
  auto ymax = get_max<T>(y);

  auto zmin = get_min<T>(z);
  auto zmax = get_max<T>(z);

  auto sigma_min = get_min<T>(sigma);
  auto sigma_max = get_max<T>(sigma);

  std::cout << "\n-------------------------------------------------------------"
               "-----------\n";
  std::cout << npts << "\t" << turb_model << "\n";
  std::cout << std::scientific << sigma_min << "\t" << sigma_max << "\n";
  std::cout << std::scientific << xmin << '\t' << xmax << '\t' << ymin << '\t'
            << ymax << '\t' << zmin << '\t' << zmax << '\n';
}

template <typename T> T STGInterface<T>::get_sigma_max() {
  return get_max<T>(sigma);
}

template <typename T> T STGInterface<T>::get_sigma_inv3() {
  T sigma_inv = 0.0;
  for (auto i = 0; i < npts; ++i)
    sigma_inv = std::max(sigma_inv, std::pow(1.0 / sigma[i], 3));
  return sigma_inv;
}
