#include "curand_kernel.h"
#include "sem.h"
#include "stg_inputs.h"
#include "stg_utils.h"
#include <fstream>
#include <vector>

template class SEM<double>;

template <typename T> __device__ T f_1d_gpu(const T x) {
  T _f1d;
  if (std::abs(x) < 1.0)
    _f1d = std::sqrt(1.5) * (1.0 - std::abs(x));
  else
    _f1d = 0.0;
  return _f1d;
}

__global__ void setup_kernel(curandState *state) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  /* Each thread gets same seed , a different sequence
  number , no offset */
  curand_init(1234, idx, 0, &state[idx]);
}

template <typename T>
__global__ void update_sem_gpu(T *xsem, T *ysem, T *zsem, T *eps_sem_x,
                               T *eps_sem_y, T *eps_sem_z, const T ub0,
                               const T dt, const T sigma_max, const T ly_sem,
                               const T lz_sem, curandState *state,
                               const int nsem) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = gridDim.x * blockDim.x;
  /* Copy state to local memory for efficiency */
  curandState localState = state[idx];
  T rnd[5];

  for (int i = idx; i < nsem; i += stride) {
    xsem[i] += dt * ub0;
    if (xsem[i] > sigma_max) {

      for (auto k = 0; k < 5; ++k)
        rnd[k] = curand_uniform(&localState);
      xsem[i] = -sigma_max;
      ysem[i] = -sigma_max + rnd[0] * ly_sem;
      zsem[i] = -sigma_max + rnd[1] * lz_sem;

      eps_sem_x[i] = rnd[2] < 0.5 ? -1.0 : 1.0;
      eps_sem_y[i] = rnd[3] < 0.5 ? -1.0 : 1.0;
      eps_sem_z[i] = rnd[4] < 0.5 ? -1.0 : 1.0;
    }
  }
}

// Using atomic reduction
template <typename T>
__global__ void get_sem_fluc_gpu(const T *xsem, const T *ysem, const T *zsem,
                                 const T *eps_sem_x, const T *eps_sem_y,
                                 const T *eps_sem_z, T *uvec, const T sigma,
                                 const T *xx, const T *a_ij, const int nsem) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = gridDim.x * blockDim.x;
  T usem = 0.0, vsem = 0.0, wsem = 0.0;

  for (int i = idx; i < nsem; i += stride) {
    T f1 = f_1d_gpu((xx[0] - xsem[i]) / sigma);
    T f2 = f_1d_gpu((xx[1] - ysem[i]) / sigma);
    T f3 = f_1d_gpu((xx[2] - zsem[i]) / sigma);

    T s1 = (a_ij[0] * eps_sem_x[i] + a_ij[1] * eps_sem_y[i] +
            a_ij[2] * eps_sem_z[i]);
    T s2 = (a_ij[3] * eps_sem_x[i] + a_ij[4] * eps_sem_y[i] +
            a_ij[5] * eps_sem_z[i]);
    T s3 = (a_ij[6] * eps_sem_x[i] + a_ij[7] * eps_sem_y[i] +
            a_ij[8] * eps_sem_z[i]);

    usem = usem + f1 * f2 * f3 * s1;
    vsem = vsem + f1 * f2 * f3 * s2;
    wsem = wsem + f1 * f2 * f3 * s3;
  }

  atomicAdd(&uvec[0], usem);
  atomicAdd(&uvec[1], vsem);
  atomicAdd(&uvec[2], wsem);
}

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

  std::cout << "Sigma max: " << sigma_max << std::endl;
  if (STGInputs::lGenerateSEM) {
      std::cout << "Generating SEM" << std::endl;
    generate_sem(stg_int);
   } else
    read_sem();
  // CUDA initialization
  size_t size = nsem * sizeof(T);

  cudaMalloc((void **)&_gpu.xsem, size);
  cudaMalloc((void **)&_gpu.ysem, size);
  cudaMalloc((void **)&_gpu.zsem, size);

  cudaMalloc((void **)&_gpu.uvec, 3 * sizeof(T));

  cudaMalloc((void **)&_gpu.eps_sem_x, size);
  cudaMalloc((void **)&_gpu.eps_sem_y, size);
  cudaMalloc((void **)&_gpu.eps_sem_z, size);

  cudaMalloc((void **)&_gpu.xx, 3 * sizeof(T));
  cudaMalloc((void **)&_gpu.a_ij, 9 * sizeof(T));
  cudaMemcpy(_gpu.xsem, xsem.data(), size, cudaMemcpyHostToDevice);
  cudaMemcpy(_gpu.ysem, ysem.data(), size, cudaMemcpyHostToDevice);
  cudaMemcpy(_gpu.zsem, zsem.data(), size, cudaMemcpyHostToDevice);

  cudaMemcpy(_gpu.eps_sem_x, eps_sem_x.data(), size, cudaMemcpyHostToDevice);
  cudaMemcpy(_gpu.eps_sem_y, eps_sem_y.data(), size, cudaMemcpyHostToDevice);
  cudaMemcpy(_gpu.eps_sem_z, eps_sem_z.data(), size, cudaMemcpyHostToDevice);

  _gpu.sigma = 0.0;
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

    eps_sem_x[i] = eps_sem_x[i] < 0.5 ? -1.0 : 1.0;
    eps_sem_y[i] = eps_sem_y[i] < 0.5 ? -1.0 : 1.0;
    eps_sem_z[i] = eps_sem_z[i] < 0.5 ? -1.0 : 1.0;
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
  size_t threads_per_block = 1024;
  size_t number_of_blocks = (nsem + threads_per_block -1) /threads_per_block;

  cudaError_t err;
  curandState *devStates;
  cudaMalloc((void **)&devStates,
             threads_per_block * number_of_blocks * sizeof(curandState));

  setup_kernel<<<number_of_blocks, threads_per_block>>>(devStates);

  update_sem_gpu<<<number_of_blocks, threads_per_block>>>(
      _gpu.xsem, _gpu.ysem, _gpu.zsem, _gpu.eps_sem_x, _gpu.eps_sem_y,
      _gpu.eps_sem_z, STGInputs::ub0, dt, sigma_max, ly_sem, lz_sem, devStates,
      nsem);
  err = cudaGetLastError();
  if (err != cudaSuccess)
    printf("Error : %s \n", cudaGetErrorString(err));

  cudaDeviceSynchronize();
  cudaFree(devStates);
}

template <typename T>
void SEM<T>::get_sem_fluc(const std::vector<T> &xx, const T tke, const T eps,
                          std::vector<T> &uvecSEM) {
  std::vector<T> Rij, aij;

  Rij.assign(6, 0.0);
  aij.assign(9, 0.0);


  STGUtils::get_re_stress(tke, Rij);

  STGUtils::re_cholesky_decomp(Rij, aij);

  T sigma = STGUtils::get_sigma(tke, eps, STGInputs::delta0, xx[3]);
  T uSEM = 0.0;
  T vSEM = 0.0;
  T wSEM = 0.0;
  cudaError_t err;

  cudaMemcpy(_gpu.a_ij, aij.data(), 9 * sizeof(T), cudaMemcpyHostToDevice);
  cudaMemcpy(_gpu.xx, xx.data(), 3 * sizeof(T), cudaMemcpyHostToDevice);

  cudaMemset(_gpu.uvec, 0, 3 * sizeof(T));

  size_t threads_per_block = 256;  //1024;
  size_t number_of_blocks = 32; //(nsem +threads_per_block -1)/threads_per_block;

  get_sem_fluc_gpu<<<number_of_blocks, threads_per_block>>>(
      _gpu.xsem, _gpu.ysem, _gpu.zsem, _gpu.eps_sem_x, _gpu.eps_sem_y,
      _gpu.eps_sem_z, _gpu.uvec, sigma, _gpu.xx, _gpu.a_ij, nsem);

  err = cudaGetLastError();
  if (err != cudaSuccess)
    printf("Error : %s \n", cudaGetErrorString(err));

  cudaDeviceSynchronize();

  cudaMemcpy(uvecSEM.data(), _gpu.uvec, 3 * sizeof(T), cudaMemcpyDeviceToHost);

  uvecSEM[0] = uvecSEM[0] *
               std::sqrt((lx_sem * ly_sem * lz_sem) / std::pow(sigma, 3)) /
               std::sqrt(nsem * 1.0);
  uvecSEM[1] = uvecSEM[1] *
               std::sqrt((lx_sem * ly_sem * lz_sem) / std::pow(sigma, 3)) /
               std::sqrt(nsem * 1.0);
  uvecSEM[2] = uvecSEM[2] *
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

template <typename T> SEM<T>::~SEM() {
  size_t size = nsem * sizeof(T);

  cudaMemcpy(xsem.data(), _gpu.xsem, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(ysem.data(), _gpu.ysem, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(zsem.data(), _gpu.zsem, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(eps_sem_x.data(), _gpu.eps_sem_x, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(eps_sem_y.data(), _gpu.eps_sem_y, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(eps_sem_z.data(), _gpu.eps_sem_z, size, cudaMemcpyDeviceToHost);


  write_sem();
  // Deallocate  GPU field variables
  std::cout << "Deallocation GPU" << std::endl;
  cudaFree(_gpu.xsem);
  cudaFree(_gpu.ysem);
  cudaFree(_gpu.zsem);
  cudaFree(_gpu.uvec);
  cudaFree(_gpu.eps_sem_x);
  cudaFree(_gpu.eps_sem_y);
  cudaFree(_gpu.eps_sem_z);
  cudaFree(_gpu.a_ij);
  cudaFree(_gpu.xx);
}
