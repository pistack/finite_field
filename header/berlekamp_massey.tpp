template<typename field>
berlekamp_massey<field>::~berlekamp_massey()
{
    if(allocated_size>0)
    delete space;
}

template<typename field>
void berlekamp_massey<field>::allocate(int size)
{
    if(allocated_size>=size)
    return;
    if(allocated_size>0)
    delete space;
    space = new field[size];
    c = space; t1 = space+size/3; t2 = space+2*(size/3);
    return;
}

template<typename field>
void berlekamp_massey<field>::allocate(field* out_space, int size)
{
    c = out_space; 
    t1 = out_space+size/3; 
    t2 = out_space+2*(size/3);
}

template<typename field>
field* berlekamp_massey<field>::operator()(field* s, int* deg, int size)
{
    field b(1), d, b_inv;
    field *r;
    t1[0] = 1; c[0] = 1;
    int L = 0, deg_t1 = 0, m = 1;
    for(int n = 0; n<size; n++)
    {
        d = s[n];
        for(int i=1; i<=L; i++)
        d += c[i]*s[n-i];
        if(d == 0)
        m++;
        else if((2*L) <= n)
        {
            for(int i=0; i<=L; i++)
            t2[i] = c[i];
            for(int i=L+1; i<=n+1-L; i++)
            c[i] = 0;
            b_inv = field(1)/b;
            for(int i=0; i<=deg_t1; i++)
            c[i+m] -= (d*b_inv)*t1[i]; 
            for(int i=0; i<=L; i++)
            t1[i] = t2[i];
            deg_t1 = L;
            b = d; m = 1; L = n+1-L;
        }
        else
        {
            b_inv = field(1)/b;
            for(int i=0; i<= deg_t1; i++)
            c[i+m] -= (d*b_inv)*t1[i];
            m++;
        }
    }
    *deg = L;
    r = new field[L+1];
    for(int i=0; i<=L; i++)
    r[i] = c[i];
    return r;
}