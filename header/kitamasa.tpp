template<typename ring>
kitamasa<ring>::~kitamasa()
{
    if(allocated_space>0)
    delete space;
}

template<typename ring>
void kitamasa<ring>::allocate(int size)
{
    if(allocated_space>=size)
    return;
    if(allocated_space>0)
    delete space;
    allocated_space = size;
    space = new ring[allocated_space];
    a = space; b = space+(size/4);
    r = space+2*(size/4);
    return;
}

template<typename ring>
void kitamasa<ring>::allocate(ring* out_space, int size)
{
    a = out_space; b = out_space+(size/4);
    r = out_space+2*(size/4);
    return;
}

template<typename ring>
int kitamasa<ring>::mult(ring* x, ring* y, ring* mod,
int deg_x, int deg_y, int deg_mod)
{
    ring tmp;
    int deg_new = deg_x+deg_y;
    for(int i=0; i<=deg_new; i++)
    r[i] = 0;
    for(int i=0; i<=deg_x; i++)
    {
        for(int j=0; j<=deg_y; j++)
        r[i+j] += x[i]*y[j];
    }
    for(int i=deg_new; i>=deg_mod; i--)
    {
        tmp = r[i]; r[i] = 0;
        for(int j=1; j<=deg_mod; j++)
        r[i-j] -= tmp*mod[j];
    }
    while(r[deg_new] == 0)
    deg_new--;
    if(deg_new < 0) deg_new = 0;
    for(int i=0; i<=deg_new; i++)
    x[i] = r[i];
    return deg_new;
}

template<typename ring>
ring kitamasa<ring>::operator()(unsigned long long int n, 
ring* init_term, ring* char_poly, int deg)
{
    int deg_a = 0, deg_b = 1;
    ring result(0);

    a[0] = 1;
    b[0] = 0; b[1] = 1;
    while(n>0)
    {
        if(n % 2 == 1)
        deg_a = mult(a, b, char_poly, deg_a, deg_b, deg);
        deg_b = mult(b, b, char_poly, deg_b, deg_b, deg);
        n >>= 1;
    }
    for(int i=0; i<=deg_a; i++)
    result += init_term[i]*a[i];
    return result;
}