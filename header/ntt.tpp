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

template<int k>
constexpr ntt<k, 998244353> arb_ntt<k>::ntt_p1;

template<int k>
constexpr ntt<k, 897581057> arb_ntt<k>::ntt_p2;

template<int k>
constexpr ntt<k, 880803841> arb_ntt<k>::ntt_p3;

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
void arb_ntt<k>::fft(int* a, finite_field<998244353>* a1, finite_field<897581057>* a2, finite_field<880803841>* a3,
int* b, finite_field<998244353>* b1, finite_field<897581057>* b2, finite_field<880803841>* b3,
int* b_inv, int size, int k_idx, int p)
{
    int fft_size = 1<<k_idx;
    int tmp;
    long long int c1, c2, c3;
    for(int i=0; i<size; i++)
    {
        b1[i] = b[i] % 998244353; b2[i] = b[i] % 897581057; b3[i] = b[i] % 880803841;
        tmp = (1ll*a[i]*b_inv[i]) % p;
        a1[i] = tmp % 998244353; a2[i] = tmp % 897581057; a3[i] = tmp % 880803841;
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
        c1 = a1[i].get_val();
        c2 = (897581057 - ((523588941*(c1 - a2[i].get_val())) % 897581057)) % 897581057;
        c3 = (880803841 - ((220201354*((c1+998244353*c2 - a3[i].get_val()) % 880803841)) % 880803841)) % 880803841;
        a[i] = (b_inv[i]*(c1+998244353*((c2+897581057*c3) % p) % p)) % p;
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
    int *b, *b_inv;
    finite_field<998244353> *a1, *b1;
    finite_field<897581057> *a2, *b2;
    finite_field<880803841> *a3, *b3;

    if((p-1)/size % 2 == 0)
    {
        w = wisdom(size, p);
        w_inv = pow(w, p-2, p);
        fft_size = 1; k_idx = 0;
        while(fft_size < 2*size)
        {fft_size <<= 1; k_idx++;}

        b = new int[size]; b_inv = new int[size];
        a1 = new finite_field<998244353>[fft_size]; 
        b1 = new finite_field<998244353>[fft_size];
        a2 = new finite_field<897581057>[fft_size]; 
        b2 = new finite_field<897581057>[fft_size];
        a3 = new finite_field<880803841>[fft_size]; 
        b3 = new finite_field<880803841>[fft_size];

        for(int i=0; i<size; i++)
        {
            b[i] = pow(w_inv, (1ll*i*i) % (2*size), p);
            b_inv[i] = pow(w, (1ll*i*i) % (2*size), p);
        }
        if(inverse)
        {
            fft(ary, a1, a2, a3, b_inv, b1, b2, b3, b, size, k_idx, p);
            size_inv = pow(size, p-2, p);
            for(int i=0; i<size; i++)
            ary[i] = (1ll*ary[i]*size_inv) % p;
        }
        else
        fft(ary, a1, a2, a3, b, b1, b2, b3, b_inv, size, k_idx, p);
    }
    else
    {
        w = wisdom(size/2, p);
        w_inv = pow(w, p-2, p);
        fft_size = 1; k_idx = 0;
        while(fft_size < size)
        {fft_size <<= 1; k_idx++;}

        b = new int[size/2]; b_inv = new int[size/2];
        a1 = new finite_field<998244353>[fft_size]; 
        b1 = new finite_field<998244353>[fft_size];
        a2 = new finite_field<897581057>[fft_size]; 
        b2 = new finite_field<897581057>[fft_size];
        a3 = new finite_field<880803841>[fft_size]; 
        b3 = new finite_field<880803841>[fft_size];

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
            fft(ary, a1, a2, a3, b_inv, b1, b2, b3, b, size/2, k_idx, p);
            fft(ary+size/2, a1, a2, a3, b_inv, b1, b2, b3, b, size/2, k_idx, p);
            for(int i=0; i<size/2; i++)
            {
                tmp1 = ary[i]; tmp2 = (1ll*ary[i+size/2]*factor) % p;
                ary[i] = (tmp1+tmp2) % p;
                ary[i+size/2] = (tmp1+(p-tmp2)) % p;
                factor = (1ll*factor*w_inv) % p; 
            }
            for(int i=0; i<size; i++)
            ary[i] = (1ll*ary[i]*size_inv) % p;
        }
        else
        {
            fft(ary, a1, a2, a3, b, b1, b2, b3, b_inv, size/2, k_idx, p);
            fft(ary+size/2, a1, a2, a3, b, b1, b2, b3, b_inv, size/2, k_idx, p);
            for(int i=0; i<size/2; i++)
            {
                tmp1 = ary[i]; tmp2 = (1ll*ary[i+size/2]*factor) % p;
                ary[i] = (tmp1+tmp2) % p;
                ary[i+size/2] = (tmp1+(p-tmp2)) % p;
                factor = (1ll*factor*w) % p; 
            }
        }
    }
    delete b; delete b_inv;
    delete a1; delete a2; delete a3;
    delete b1; delete b2; delete b3;
    return;
}

