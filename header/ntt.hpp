#ifndef NTT_H
#define NTT_H

#include "finite_field.hpp"

template<int k, int p>
class ntt
{
    private:
    finite_field<p> w = 0, w_inv = 0;


    public:
    constexpr ntt();
    void operator()(finite_field<p>* ary, 
    int k_idx, bool inverse) const;
};

template<int k, int p>
class convolution
{
    private:
    static constexpr ntt<k, 998244353> ntt_1 = ntt<k, 998244353>();
    static constexpr ntt<k, 897581057> ntt_2 = ntt<k, 897581057>();
    static constexpr ntt<k, 880803841> ntt_3 = ntt<k, 880803841>();

    void convolution_kernel_1(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r, 
    int size_a, int size_b, int k_idx);
    void convolution_kernel_2(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
    int size_a, int size_b, int k_idx);
    void convolution_kernel_3(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
    int size_a, int size_b, int k_idx);

    public:
    constexpr convolution() {};
    void operator()(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
    int size_a, int size_b, int k_idx);
};

template<int k>
class arb_ntt
{
    private:
    static constexpr int p1 = 998244353, p2 = 897581057, p3 = 880803841; 
    static constexpr __int128_t pqr = (__int128_t)p1*p2*p3;
    static constexpr __int128_t factor_p1 = ((__int128_t)p2*p3*
    (finite_field<p1>(1)/finite_field<p1>((1ll*p2*p3) % p1)).get_val()) \
    % pqr;
    static constexpr __int128_t factor_p2 = ((__int128_t)p1*p3*
    (finite_field<p2>(1)/finite_field<p2>((1ll*p1*p3) % p2)).get_val()) \
    % pqr;
    static constexpr __int128_t factor_p3 = ((__int128_t)p1*p2*
    (finite_field<p3>(1)/finite_field<p3>((1ll*p1*p2) % p3)).get_val()) \
    % pqr;
    static constexpr ntt<k, p1> ntt_p1;
    static constexpr ntt<k, p2> ntt_p2;
    static constexpr ntt<k, p3> ntt_p3;
    int allocated_idx = -1;
    int *space_i;
    finite_field<p1> *space_p1;
    finite_field<p2> *space_p2;
    finite_field<p3> *space_p3;
    int *b, *b_inv;
    finite_field<p1> *a1, *b1;
    finite_field<p2> *a2, *b2;
    finite_field<p3> *a3, *b3;
    
    int pow(long long int base, int idx, int p);
    int wisdom(int size, int p);
    void fft(int* a, int* b, int* b_inv, int size, int k_idx, int p);

    public:
    arb_ntt() {};
    void allocate(int k_idx);
    void allocate(int* out_space_i, 
    finite_field<p1>* out_space_p1, finite_field<p2>* out_space_p2,
    finite_field<p3>* out_space_p3, int k_idx);
    void operator()(int* ary, int size, int p, bool inverse);
    ~arb_ntt();
};

#include "ntt.tpp"

#endif