template<int p>
constexpr finite_field<p> & finite_field<p>::operator/=(const finite_field<p> & div)
{
    int a = (div.value<0 ? p+div.value : div.value); 
    int b = p, x0 = 0, x1 = 1, q = 0, t = 0;
    while(a > 0)
    {
        q = b/a; t = x0; x0 = x1; x1 = t - q*x1;
        t = a; a = b % a; b = t;
    }
    x0 = (x0 < 0 ? p+x0 : x0);
    this -> value = (1ll*(this->value)*x0) % p;
    return *this;
}

template<typename ED>
ED gcd(ED a, ED b)
{
    ED tmp;
    while(b != 0)
    {tmp = b; b = a % b; a = tmp;}
    return a;
}

template<typename ED, typename... EDs>
ED gcd(ED a, ED b, EDs... cs)
{
    return gcd(a, gcd(b, cs...));
}

template<typename mon>
constexpr mon pow(mon base, unsigned long long int idx)
{
    mon result = 1;
    while(idx > 0)
    {
        if(idx % 2 == 1)
        result *= base;
        base *= base; idx >>= 1;
    }
    return result;
}

template<int p>
constexpr finite_field<p> quad_non_res()
{
    finite_field<p> val(0);
    for(int i=2; i<p; i++)
    {
        val = i;
        if(pow(val, (p-1)>>1) != 1)
        break;
    }
    return val;
}

template<int p>
constexpr finite_field<p> tonelli_shanks<p>::operator()(finite_field<p> x)
{
    int q = p-1, m = 0, i = 0;
    finite_field<p> b(0), c(0), t(0), r(0);
    if(!(pow(x, (p-1)>>1)==1))
    return finite_field<p>(0);
    while(m % 2 == 0)
    {q >>= 1; m++;}
    c = pow(non_residual, q);
    t = pow(x, q);
    r = pow(x, (q+1)>>1);
    while(!(t == 1))
    {
        b = c;
        for(int j=1; j<m; j++)
        {
            t *= t;
            if(t == 1)
            {i = j; break;}
        }
        for(int j=1; j<m-i; j++)
        b *= b;
        m = i; c = b*b; t *= c;
        r *= b;
    }
    return r;
}