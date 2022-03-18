template<int k, int p>
constexpr ntt<k+2, p> multipoint_evaluation<k, p>::ntt_p;

template<int k, int p>
multipoint_evaluation<k, p>::~multipoint_evaluation()
{
    if(allocated_idx>-1)
    delete space;
}

template<int k, int p>
void multipoint_evaluation<k, p>::allocate(int k_idx)
{
    if(allocated_idx>=k_idx)
    return;
    if(allocated_idx>-1)
    delete space;
    space = new finite_field<p>[(k_idx+21)<<k_idx];
    poly_tree = space;
    collect = space+(k_idx<<k_idx);
    ext_f = space+((k_idx+1)<<k_idx);
    ext_g = space+((k_idx+5)<<k_idx);
    tmp = space+((k_idx+9)<<k_idx);
    tmp2 = space+((k_idx+13)<<k_idx);
    tmp3 = space+((k_idx+17)<<k_idx);
    return;
}

template<int k, int p>
void multipoint_evaluation<k, p>::allocate(finite_field<p>* out_space, int k_idx)
{
    poly_tree = out_space;
    collect = out_space+(k_idx<<k_idx);
    ext_f = out_space+((k_idx+1)<<k_idx);
    ext_g = out_space+((k_idx+5)<<k_idx);
    tmp = out_space+((k_idx+9)<<k_idx);
    tmp2 = out_space+((k_idx+13)<<k_idx);
    tmp3 = out_space+((k_idx+17)<<k_idx);
    return;    
}

template<int k, int p>
int multipoint_evaluation<k, p>::construct(finite_field<p>* x, int size)
{
    int k_idx = 0;
    int pos1, pos2;
    int half_size, curr_size, tot_size = 1;
    int fft_curr_size, fft_tot_size;

    while(tot_size < size)
    {k_idx++; tot_size <<= 1;}
    fft_tot_size = tot_size << 2;
    for(int i=0; i<size; i++)
    poly_tree[i] = finite_field<p>(-1)*x[i];
    for(int i=size; i<tot_size; i++)
    poly_tree[i] = 0;
    curr_size = 1; pos1 = 0;
    for(int k_curr=1; k_curr < k_idx; k_curr++)
    {
        half_size = curr_size;
        curr_size <<= 1; 
        fft_curr_size = curr_size << 1;
        pos2 = 0;
        for(int i=0; i<fft_tot_size; i += fft_curr_size)
        {
            for(int j=0; j<half_size; j++)
            tmp[i+j] = poly_tree[pos1+pos2+j];
            tmp[i+half_size] = 1;
            for(int j=half_size+1; j<fft_curr_size; j++)
            tmp[i+j] = 0;
            pos2 += half_size;
        }
        pos1 += tot_size; pos2 = 0;
        for(int i=0; i<fft_tot_size; i += (fft_curr_size<<1))
        {
            ntt_p(tmp+i, k_curr+1, false);
            ntt_p(tmp+fft_curr_size+i, k_curr+1, false);
            for(int j=0; j<fft_curr_size; j++)
            tmp[i+j] *= tmp[i+fft_curr_size+j];
            ntt_p(tmp+i, k_curr+1, true);
            for(int j=0; j<curr_size; j++)
            poly_tree[pos1+pos2+j] = tmp[i+j];
            pos2 += curr_size;
        }
    }
    return k_idx;
}

template<int k, int p>
void multipoint_evaluation<k, p>::inverse(finite_field<p>* p1,
finite_field<p>* p1_inv, finite_field<p>* p1_rev,
finite_field<p>* p1_inv_pre, int deg_p1, int deg_inv)
{
    int deg_next, deg_curr, deg_rev;
    int k_idx, half_size_inv, size_inv;
    deg_curr = 0; k_idx = 1; size_inv = 2;
    p1_inv[0] = 1; p1_inv[1] = 0;
    while(deg_curr < deg_inv)
    {
        deg_next = (deg_curr<<1)+1;
        deg_next = (deg_next < deg_inv ? deg_next : deg_inv);
        half_size_inv = size_inv; size_inv <<= 1; k_idx++;
        for(int i=0; i<=deg_curr; i++)
        p1_inv_pre[i] = p1_inv[i];
        for(int i=half_size_inv; i<size_inv; i++)
        p1_inv[i] = 0;
        deg_rev = (deg_next < deg_p1 ? deg_next : deg_p1);
        for(int i=0; i<=deg_rev; i++)
        p1_rev[i] = p1[deg_p1-i];
        for(int i=deg_rev+1; i<size_inv; i++)
        p1_rev[i] = 0;
        ntt_p(p1_rev, k_idx, false);
        ntt_p(p1_inv, k_idx, false);
        for(int i=0; i<size_inv; i++)
        {
            p1_rev[i] *= p1_inv[i];
            p1_rev[i] = finite_field<p>(1) - p1_rev[i];
            p1_inv[i] *= p1_rev[i];
        }
        ntt_p(p1_inv, k_idx, true);
        for(int i=0; i<=deg_curr; i++)
        p1_inv[i] = p1_inv_pre[i];
        for(int i=deg_next+1; i<size_inv; i++)
        p1_inv[i] = 0;
        deg_curr = deg_next;
    }
}

template<int k, int p>
void multipoint_evaluation<k, p>::remainder(finite_field<p>* f, 
finite_field<p>* g, finite_field<p>* r, 
finite_field<p>* supp1, finite_field<p>* supp2,
int deg_f, int deg_g, int k_idx)
{
    int size = 1<<k_idx;
    int deg_inv = deg_f - deg_g;
    finite_field<p> swp;

    if(deg_g > deg_f)
    {
        for(int i=0; i<size; i++)
        r[i] = f[i];
        return;
    }
    inverse(g, r, supp1, supp2, deg_g, deg_inv);
    for(int i=0; i<=deg_inv; i++)
    supp1[i] = f[deg_f-i];
    for(int i=deg_inv+1; i<size; i++)
    {supp1[i] = 0; r[i] = 0;}
    ntt_p(supp1, k_idx, false);
    ntt_p(g, k_idx, false);
    ntt_p(r, k_idx, false);
    for(int i=0; i<size; i++)
    r[i] *= supp1[i];
    ntt_p(r, k_idx, true);
    for(int i=deg_inv+1; i<size; i++)
    r[i] = 0;
    for(int i=0; i<=deg_inv/2; i++)
    {swp = r[i]; r[i] = r[deg_inv-i]; r[deg_inv-i] = swp;}
    ntt_p(r, k_idx, false);
    for(int i=0; i<size; i++)
    g[i] *= r[i];
    ntt_p(g, k_idx, true);
    for(int i=0; i<size; i++)
    r[i] = f[i]-g[i];
    return;
}

template<int k, int p>
finite_field<p>* multipoint_evaluation<k, p>::operator()(finite_field<p>* f, int deg, 
finite_field<p>* x, int num_pts)
{
    int k_idx, fft_tot_size, fft_f_size, tot_size, curr_size;
    int fft_size, fft_k_idx, start, pos, deg_f;
    finite_field<p>* ev = new finite_field<p>[num_pts];

    if(deg == 0)
    {
        for(int i=0; i<num_pts; i++)
        ev[i] = f[0];
        return ev;
    }
    if(deg == 1)
    {
        for(int i=0; i<num_pts; i++)
        ev[i] = f[0]+f[1]*x[i];
        return ev;
    }
    if(num_pts == 1)
    {
        ev[0] = (f[deg]*x[0]+f[deg-1]);
        for(int i=deg-2; i>=0; i--)
        {ev[0] *= x[0]; ev[0] += f[i];}
        return ev;
    }

    // construct polynomial sagment tree
    k_idx = construct(x, num_pts);

    // first division step
    tot_size = 1<<k_idx; curr_size = tot_size>>1;
    fft_size = tot_size<<1; fft_k_idx = k_idx+1;
    fft_f_size = 2*deg+1;

    while(fft_size < fft_f_size)
    {fft_size <<= 1; fft_k_idx++;}

    start = (k_idx-1)*tot_size;

    // extend polynomial f 
    for(int i=0; i<=deg; i++)
    ext_f[i] = f[i];
    for(int i=deg+1; i<fft_size; i++)
    ext_f[i] = 0;

    // extend n/2 degree polynomial g1
    for(int i=0; i<curr_size; i++)
    ext_g[i] = poly_tree[start+i];
    ext_g[curr_size] = 1;
    for(int i=curr_size+1;i<fft_size; i++)
    ext_g[i] = 0;

    // f % g1
    remainder(ext_f, ext_g, 
    tmp, tmp2, tmp3, 
    deg, curr_size, fft_k_idx);

    // collect f % g1
    for(int i=0; i<curr_size; i++)
    collect[i] = tmp[i];

    // extend n/2 degree polynomial g2
    for(int i=0; i<curr_size; i++)
    ext_g[i] = poly_tree[start+curr_size+i];
    ext_g[curr_size] = 1;
    for(int i=curr_size+1; i<fft_size; i++)
    ext_g[i] = 0;

    // f % g2
    remainder(ext_f, ext_g,
    tmp, tmp2, tmp3,
    deg, curr_size, fft_k_idx);

    // collect f % g2
    for(int i=0; i<curr_size; i++)
    collect[i+curr_size] = tmp[i];
    fft_tot_size = tot_size << 2; fft_f_size = tot_size << 1;
    curr_size = tot_size >> 1; fft_size = tot_size;
    fft_k_idx = k_idx;

    // visits polynomial sagment tree
    for(int k_i = k_idx-1; k_i > 0; k_i--)
    {
        pos = 0; deg_f = curr_size-1;
        for(int i=0; i<fft_f_size; i += fft_size)
        {
            for(int j=0; j<curr_size; j++)
            ext_f[i+j] = collect[pos+j];
            for(int j=curr_size; j<fft_size; j++)
            ext_f[i+j] = 0;
            pos += curr_size;
        }
        start -= tot_size; curr_size >>= 1; pos = start;
        for(int i=0; i<fft_tot_size; i += fft_size)
        {
            for(int j=0; j<curr_size; j++)
            ext_g[i+j] = poly_tree[pos+j];
            ext_g[i+curr_size] = 1;
            for(int j=curr_size+1; j<fft_size; j++)
            ext_g[i+j] = 0;
            pos += curr_size;
        }
        pos = 0;
        for(int i=0; i<fft_f_size; i += fft_size)
        {
            remainder(ext_f+i, ext_g+pos,
            tmp+pos, tmp2+pos, tmp3+pos, 
            deg_f, curr_size, fft_k_idx);
            remainder(ext_f+i, ext_g+pos+fft_size,
            tmp+pos+fft_size, tmp2+pos+fft_size, tmp3+pos+fft_size,
            deg_f, curr_size, fft_k_idx);
            pos += fft_size<<1;
        }
        pos = 0;
        for(int i=0; i<tot_size; i += curr_size)
        {
            for(int j=0; j<curr_size; j++)
            collect[i+j] = tmp[pos+j];
            pos += fft_size;
        }
        fft_size >>= 1; fft_k_idx--;
    }
    for(int i=0; i<num_pts; i++)
    ev[i] = collect[i];
    return ev;
}