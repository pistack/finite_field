template<typename field>
field* berlekamp_massey(field* s, int* deg, int size)
{
    field b(1), d, b_inv;
    field *r, *c, *t1, *t2;
    c = new field[size]; t1 = new field[size]; t2 = new field[size];
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
    delete c; delete t1; delete t2;
    return r;
}