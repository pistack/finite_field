template<int k, int p>
constexpr ntt<k, p>::ntt()
{
    int idx_2 = 0;
    int tmp = p-1;
    finite_field<p> candidate = quad_non_res<p>();

    while(tmp % 2 == 0)
    {idx_2++; tmp >>= 1;}
    candidate = pow(candidate, tmp<<(idx_2-k));
    w = candidate; w_inv = finite_field<p>(1)/w;
    return;
}

template<int k, int p>
void ntt<k, p>::operator()(finite_field<p>* ary, int k_idx, bool inverse) const
{
    finite_field<p> w_idx = (inverse ? pow(w_inv, 1<<(k-k_idx)) : 
    pow(w, 1<<(k-k_idx)));
    finite_field<p> inv_size, tmp_1, tmp_2;
    finite_field<p> factor_curr, factor_accum;
    int size = 1<<k_idx;
    int dist, mask, half_curr_size;
    int j = 0;

    for(int i=1; i<size; i++)
    {
        mask = (size >> 1);
        while(j >= mask)
        {j -= mask; mask >>= 1;}
        j += mask;
        if(i < j)
        {tmp_1 = ary[i]; ary[i] = ary[j]; ary[j] = tmp_1;}
    }
    dist = size;
    for(int curr_size = 2; curr_size <= size; curr_size <<= 1)
    {
        dist >>= 1;
        factor_curr = pow(w_idx, dist);
        for(int from=0; from<size; from += curr_size)
        {
            factor_accum = 1; half_curr_size = curr_size >>1;
            for(int i=0; i<half_curr_size; i++)
            {
                tmp_1 = ary[from+i];
                tmp_2 = factor_accum*ary[from+half_curr_size+i];
                ary[from+i] = tmp_1 + tmp_2;
                ary[from+half_curr_size+i] = tmp_1 - tmp_2;
                factor_accum *= factor_curr;
            }
        }
    }

    if(inverse)
    {
        inv_size = finite_field<p>(1)/finite_field<p>(size);
        for(int i=0; i<size; i++)
        ary[i] *= inv_size;
    }
    return;
}

template<int k, int p>
constexpr ntt<k, 998244353> convolution<k, p>::ntt_1;

template<int k, int p>
constexpr ntt<k, 897581057> convolution<k, p>::ntt_2;

template<int k, int p>
constexpr ntt<k, 880803841> convolution<k, p>::ntt_3;

template<int k, int p>
void convolution<k, p>::convolution_kernel_1(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
int size_a, int size_b, int k_idx)
{
    finite_field<998244353> *ext_a, *ext_b;
    ext_a = new finite_field<998244353>[1<<k_idx];
    ext_b = new finite_field<998244353>[1<<k_idx];
    for(int i=0; i<size_a; i++)
    ext_a[i] = a[i].get_val() % 998244353;
    for(int i=0; i<size_b; i++)
    ext_b[i] = b[i].get_val() % 998244353;
    for(int i=size_a; i<(1<<k_idx); i++)
    ext_a[i] = 0;
    for(int i=size_b; i<(1<<k_idx); i++)
    ext_b[i] = 0;
    ntt_1(ext_a, k_idx, false); ntt_1(ext_b, k_idx, false);
    for(int i=0; i<(1<<k_idx); i++)
    ext_a[i] *= ext_b[i];
    ntt_1(ext_a, k_idx, true);
    for(int i=0; i<size_a+size_b-1; i++)
    r[i] = ext_a[i].get_val() % p;
    delete ext_a; delete ext_b;
    return;
}

template<int k, int p>
void convolution<k, p>::convolution_kernel_2(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
int size_a, int size_b, int k_idx)
{
    constexpr long long int n = 1ll*998244353*897581057;
    constexpr long long int f1 = (1ll*897581057*415935157);
    constexpr long long int f2 = (1ll*998244353*523588941);
    finite_field<p> *result = new finite_field<p>[size_a+size_b-1];
    finite_field<998244353> *ext_a1, *ext_b1;
    finite_field<897581057> *ext_a2, *ext_b2;

    ext_a1 = new finite_field<998244353>[1<<k_idx];
    ext_b1 = new finite_field<998244353>[1<<k_idx];
    ext_a2 = new finite_field<897581057>[1<<k_idx];
    ext_b2 = new finite_field<897581057>[1<<k_idx];

    for(int i=0; i<size_a; i++)
    {ext_a1[i] = a[i].get_val() % 998244353; ext_a2[i] = a[i].get_val() % 897581057;}
    for(int i=0; i<size_b; i++)
    {ext_b1[i] = b[i].get_val() % 998244353; ext_b2[i] = b[i].get_val() % 897581057;}
    for(int i=size_a; i<(1<<k_idx); i++)
    {ext_a1[i] = 0; ext_a2[i] = 0;}
    for(int i=size_b; i<(1<<k_idx); i++)
    {ext_b1[i] = 0; ext_b2[i] = 0;}
    ntt_1(ext_a1, k_idx, false); ntt_1(ext_b1, k_idx, false);
    ntt_2(ext_a2, k_idx, false); ntt_2(ext_b2, k_idx, false);
    for(int i=0; i<(1<<k_idx); i++)
    {ext_a1[i] *= ext_b1[i]; ext_a2[i] *= ext_b2[i];}
    ntt_1(ext_a1, k_idx, true); ntt_2(ext_a2, k_idx, true);
    for(int i=0; i<size_a+size_b-1; i++)
    {
        r[i] = (((__int128_t)f1*ext_a1[i].get_val()+
        (__int128_t)f2*ext_a2[i].get_val()) % n) % p;
    }
    delete ext_a1; delete ext_b1;
    delete ext_a2; delete ext_b2;
    return;
}

template<int k, int p>
void convolution<k, p>::convolution_kernel_3(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
int size_a, int size_b, int k_idx)
{
    constexpr __int128_t n = (__int128_t)998244353*897581057*880803841;
    constexpr __int128_t f1 = (__int128_t)41593599*897581057*880803841;
    constexpr __int128_t f2 = (__int128_t)998244353*635786105*880803841;
    constexpr __int128_t f3 = (__int128_t)998244353*897581057*220201354;

    finite_field<p> *result = new finite_field<p>[size_a+size_b-1];
    finite_field<998244353> *ext_a1, *ext_b1;
    finite_field<897581057> *ext_a2, *ext_b2;
    finite_field<880803841> *ext_a3, *ext_b3;

    ext_a1 = new finite_field<998244353>[1<<k_idx];
    ext_b1 = new finite_field<998244353>[1<<k_idx];
    ext_a2 = new finite_field<897581057>[1<<k_idx];
    ext_b2 = new finite_field<897581057>[1<<k_idx];
    ext_a3 = new finite_field<880803841>[1<<k_idx];
    ext_b3 = new finite_field<880803841>[1<<k_idx];

    for(int i=0; i<size_a; i++)
    {
        ext_a1[i] = a[i].get_val() % 998244353; 
        ext_a2[i] = a[i].get_val() % 897581057;
        ext_a3[i] = a[i].get_val() % 880803841;
    }
    for(int i=0; i<size_b; i++)
    {
        ext_b1[i] = b[i].get_val() % 998244353; 
        ext_b2[i] = b[i].get_val() % 897581057;
        ext_b3[i] = b[i].get_val() % 880803841;
    }

    for(int i=size_a; i<(1<<k_idx); i++)
    {ext_a1[i] = 0; ext_a2[i] = 0; ext_a3[i] = 0;}
    for(int i=size_b; i<(1<<k_idx); i++)
    {ext_b1[i] = 0; ext_b2[i] = 0; ext_b3[i] = 0;}

    ntt_1(ext_a1, k_idx, false); ntt_1(ext_b1, k_idx, false);
    ntt_2(ext_a2, k_idx, false); ntt_2(ext_b2, k_idx, false);
    ntt_3(ext_a3, k_idx, false); ntt_3(ext_b3, k_idx, false);

    for(int i=0; i<(1<<k_idx); i++)
    {ext_a1[i] *= ext_b1[i]; ext_a2[i] *= ext_b2[i]; ext_a3[i] *= ext_b3[i];}
    ntt_1(ext_a1, k_idx, true); ntt_2(ext_a2, k_idx, true); ntt_3(ext_a3, k_idx, true);
    for(int i=0; i<size_a+size_b-1; i++)
    {
        r[i] = (((__int128_t)f1*ext_a1[i].get_val()+
        (__int128_t)f2*ext_a2[i].get_val()+
        (__int128_t)f3*ext_a3[i].get_val()) % n) % p;
    }
    delete ext_a1; delete ext_b1;
    delete ext_a2; delete ext_b2;
    delete ext_a3; delete ext_b3;
    return;
}

template<int k, int p>
void convolution<k, p>::operator()(finite_field<p>* a, finite_field<p>* b, finite_field<p>* r,
int size_a, int size_b, int k_idx)
{
    if(1ll*(p-1)*(p-1) < 998244353>>k)
    return convolution_kernel_1(a, b, r, size_a, size_b, k_idx);
    if(1ll*(p-1)*(p-1)<(1ll*998244353*897581057)>>k)
    return convolution_kernel_2(a, b, r, size_a, size_b, k_idx);
    return convolution_kernel_3(a, b, r, size_a, size_b, k_idx);
}

template<int k>
constexpr int arb_ntt<k>::p1;

template<int k>
constexpr int arb_ntt<k>::p2;

template<int k>
constexpr int arb_ntt<k>::p3;

template<int k>
constexpr __int128_t arb_ntt<k>::pqr;

template<int k>
constexpr __int128_t arb_ntt<k>::factor_p1;

template<int k>
constexpr __int128_t arb_ntt<k>::factor_p2;

template<int k>
constexpr __int128_t arb_ntt<k>::factor_p3;

template<int k>
constexpr ntt<k, 998244353> arb_ntt<k>::ntt_p1;

template<int k>
constexpr ntt<k, 897581057> arb_ntt<k>::ntt_p2;

template<int k>
constexpr ntt<k, 880803841> arb_ntt<k>::ntt_p3;

template<int k>
arb_ntt<k>::~arb_ntt()
{
    if(allocated_idx>-1)
    {
        delete space_i; delete space_p1;
        delete space_p2; delete space_p3;
    }
}

template<int k>
void arb_ntt<k>::allocate(int k_idx)
{
    int size = 1<<k_idx;
    if(allocated_idx>=k_idx)
    return;
    if(allocated_idx>-1)
    {
        delete space_i; delete space_p1;
        delete space_p2; delete space_p3;
    }
    space_i = new int[size];
    space_p1 = new finite_field<p1>[2*size];
    space_p2 = new finite_field<p2>[2*size];
    space_p3 = new finite_field<p3>[2*size];
    b = space_i; b_inv = space_i+size/2;
    a1 = space_p1; b1 = space_p1+size;
    a2 = space_p2; b2 = space_p2+size;
    a3 = space_p3; b3 = space_p3+size;
    return;
}

template<int k>
void arb_ntt<k>::allocate(int* out_space_i,
finite_field<p1>* out_space_p1, finite_field<p2>* out_space_p2,
finite_field<p3>* out_space_p3, int k_idx)
{
    int size = 1<<k_idx;
    b = out_space_i; b_inv = out_space_i+size/2;
    a1 = out_space_p1; b1 = out_space_p1+size;
    a2 = out_space_p2; b2 = out_space_p2+size;
    a3 = out_space_p3; b3 = out_space_p3+size;
    return;
}

template<int k>
int arb_ntt<k>::pow(long long int base, int idx, int p)
{
    long long int result = 1;
    while(idx > 0)
    {
        if(idx % 2 == 1)
        {result *= base; result %= p;}
        base *= base; base %= p; idx >>= 1;
    }
    return result;
}

template<int k>
int arb_ntt<k>::wisdom(int size, int p)
{
    bool is_primitive;
    int tmp, num_fact, factor[20], g;
    if(p == 2)
    return 1;
    tmp = p-1; num_fact = 0;
    for(int i=2; i*i<=tmp; i++)
    {
        if(tmp % i == 0)
        {
            factor[num_fact] = i;
            num_fact++;
            while(tmp % i == 0)
            tmp /= i;
        }
    }
    for(g = 2; g <= p-1; g++)
    {
        is_primitive = true;
        for(int i=0; i<num_fact; i++)
        {
            if(pow(g, (p-1)/factor[i], p) == 1)
            {
                is_primitive = false;
                break;
            }
        }
        if(is_primitive)
        break;
    }
    g = pow(g, (p-1)/(2*size), p);
    return g;
}

template<int k>
void arb_ntt<k>::fft(int* a, int* b, int* b_inv, int size, int k_idx, int p)
{
    int fft_size = 1<<k_idx;
    int tmp;
    __int128_t tmp_128;
    for(int i=0; i<size; i++)
    {
        b1[i] = b[i] % p1; b2[i] = b[i] % p2; b3[i] = b[i] % p3;
        tmp = ((1ll)*a[i]*b_inv[i]) % p;
        a1[i] = tmp % p1; a2[i] = tmp % p2; a3[i] = tmp % p3;
    }
    for(int i=size; i<fft_size; i++)
    {b1[i] = 0; b2[i] = 0; b3[i] = 0; a1[i] = 0; a2[i] = 0; a3[i] = 0;}
    for(int i=1; i<size; i++)
    {b1[fft_size-i] = b1[i]; b2[fft_size-i] = b2[i]; b3[fft_size-i] = b3[i];}
    ntt_p1(a1, k_idx, false); ntt_p2(a2, k_idx, false); ntt_p3(a3, k_idx, false);
    ntt_p1(b1, k_idx, false); ntt_p2(b2, k_idx, false); ntt_p3(b3, k_idx, false);
    for(int i=0; i<fft_size; i++)
    {a1[i] *= b1[i]; a2[i] *= b2[i]; a3[i] *= b3[i];}
    ntt_p1(a1, k_idx, true); ntt_p2(a2, k_idx, true); ntt_p3(a3, k_idx, true);
    for(int i=0; i<size; i++)
    {
        tmp_128 = (((((a1[i].get_val()*factor_p1) % pqr)+
        a2[i].get_val()*factor_p2) % pqr)+a3[i].get_val()*factor_p3) % pqr;
        tmp_128 = (tmp_128 < 0 ? pqr+tmp_128 : tmp_128);
        a[i] = tmp_128 % p; a[i] = (1ll*a[i]*b_inv[i]) % p;  
    }
    return;
}

template<int k>
void arb_ntt<k>::operator()(int* ary, int size, int p, bool inverse)
{
    int fft_size;
    int k_idx;
    int w, w_inv, size_inv, factor;
    int tmp1, tmp2;

    if((p-1)/size % 2 == 0)
    {
        w = wisdom(size, p);
        w_inv = pow(w, p-2, p);
        fft_size = 1; k_idx = 0;
        while(fft_size < 2*size)
        {fft_size <<= 1; k_idx++;}
        for(int i=0; i<size; i++)
        {
            b[i] = pow(w_inv, (1ll*i*i) % (2*size), p);
            b_inv[i] = pow(w, (1ll*i*i) % (2*size), p);
        }
        if(inverse)
        {
            fft(ary, b_inv, b, size, k_idx, p);
            size_inv = pow(size, p-2, p);
            for(int i=0; i<size; i++)
            ary[i] = (1ll*ary[i]*size_inv) % p;
        }
        else
        fft(ary, b, b_inv, size, k_idx, p);
    }
    else
    {
        w = wisdom(size/2, p);
        w_inv = pow(w, p-2, p);
        fft_size = 1; k_idx = 0;
        while(fft_size < size)
        {fft_size <<= 1; k_idx++;}
    
        for(int i=0; i<size/2; i++)
        b[i] = ary[2*i+1];
        for(int i=0; i<size/2; i++)
        ary[i] = ary[2*i];
        for(int i=0; i<size/2; i++)
        ary[size/2+i] = b[i];
    
        for(int i=0; i<size/2; i++)
        {
            b[i] = pow(w_inv, (1ll*i*i) % size, p);
            b_inv[i] = pow(w, (1ll*i*i) % size, p);
        }
        factor = 1;
        if(inverse)
        {
            size_inv = pow(size, p-2, p);
            fft(ary, b_inv, b, size/2, k_idx, p);
            fft(ary+size/2, b_inv, b, size/2, k_idx, p);
            for(int i=0; i<size/2; i++)
            {
                tmp1 = ary[i]; tmp2 = (1ll*ary[i+size/2]*factor) % p;
                ary[i] = (tmp1+tmp2) % p;
                ary[i+size/2] = (tmp1+p-tmp2) % p;
                factor = (1ll*factor*w_inv) % p; 
            }
            for(int i=0; i<size; i++)
            ary[i] = (1ll*ary[i]*size_inv) % p;
        }
        else
        {
            fft(ary, b, b_inv, size/2, k_idx, p);
            fft(ary+size/2, b, b_inv, size/2, k_idx, p);
            for(int i=0; i<size/2; i++)
            {
                tmp1 = ary[i]; tmp2 = (1ll*ary[i+size/2]*factor) % p;
                ary[i] = (tmp1+tmp2) % p;
                ary[i+size/2] = (tmp1+p-tmp2) % p;
                factor = (1ll*factor*w) % p; 
            }
        }
    }
    return;
}

