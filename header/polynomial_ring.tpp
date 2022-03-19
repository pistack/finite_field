template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(const int & c)
{
    coeff = new finite_field<p>[max_deg+1]; deg = 0;
    coeff[0] = c;
}

template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(const finite_field<p> & c)
{
    coeff = new finite_field<p>[max_deg+1]; deg = 0;
    coeff[0] = c;
}

template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(int d, finite_field<p>* c)
{
    deg = d; coeff = new finite_field<p>[max_deg+1];
    for(int i=0; i<=deg; i++)
    coeff[i] = c[i];
}

template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(const polynomial_ring<max_deg, p> & poly)
{
    deg = poly.deg; coeff = new finite_field<p>[max_deg+1];
    for(int i=0; i<=deg; i++)
    coeff[i] = poly.coeff[i];
}

template<int max_deg, int p>
polynomial_ring<max_deg, p>::~polynomial_ring()
{
    if(this->deg > -1)
    delete coeff;
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator<(const polynomial_ring<max_deg, p> & poly) const
{
    if(deg != poly.deg)
    return (deg < poly.deg);
    for(int i=deg; i>=0; i--)
    {
        if(this -> coeff[i] != poly.coeff[i])
        return (this -> coeff[i] < poly.coeff[i]);
    }
    return false;
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator>(const polynomial_ring<max_deg, p> & poly) const
{
    if(deg != poly.deg)
    return (deg > poly.deg);
    for(int i=deg; i>=0; i--)
    {
        if(this -> coeff[i] != poly.coeff[i])
        return (this -> coeff[i] > poly.coeff[i]);
    }
    return false;
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator==(const polynomial_ring<max_deg, p> & poly) const
{
    if(deg != poly.deg)
    return false;
    for(int i=deg; i>=0; i--)
    {
        if(coeff[i] != poly.coeff[i])
        return false;
    }
    return true;
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator!=(const polynomial_ring<max_deg, p> & poly) const
{
    if(deg == poly.deg)
    return false;
    for(int i=deg; i>=0; i--)
    {
        if(coeff[i] == poly.coeff[i])
        return false;
    }
    return true;
}

template<int max_deg, int p>
polynomial_ring<max_deg , p> & polynomial_ring<max_deg, p>::operator=(const polynomial_ring<max_deg, p> & assgin)
{
    if(this == &assgin)
    return *this;
    if(this -> deg == -1)
    this -> coeff = new finite_field<p>[max_deg+1];

    this->deg = assgin.deg;
    for(int i=0; i<=(this->deg); i++)
    this->coeff[i] = assgin.coeff[i];
    return *this;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator+(const polynomial_ring<max_deg, p> & add) const
{
    polynomial_ring<max_deg, p> poly;
    int common_deg;

    poly.coeff = new finite_field<p>[max_deg+1];
    if(this->deg > add.deg)
    {
        poly.deg = this->deg;
        for(int i=poly.deg; i>add.deg; i--)
        poly.coeff[i] = this -> coeff[i];
        common_deg = add.deg;
    }
    else
    {
        poly.deg = add.deg;
        for(int i=poly.deg; i>this->deg; i--)
        poly.coeff[i] = add.coeff[i];
        common_deg = this->deg;
    }
    for(int i=0; i<=common_deg; i++)
    poly.coeff[i] = this->coeff[i] + poly.coeff[i];
    while(poly.coeff[poly.deg] == 0)
    poly.deg--;
    if(poly.deg < 0)
    poly.deg = 0;
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator-(const polynomial_ring<max_deg, p> & sub) const
{
    polynomial_ring<max_deg, p> poly;
    int common_deg;
    poly.coeff = new finite_field<p>[max_deg+1];
    if(this->deg > sub.deg)
    {
        poly.deg = this->deg;
        for(int i=poly.deg; i>sub.deg; i--)
        poly.coeff[i] = this->coeff[i];
        common_deg = sub.deg;
    }
    else
    {
        poly.deg = sub.deg;
        for(int i=poly.deg; i>this->deg; i--)
        poly.coeff[i] = sub.coeff[i];
        common_deg = this->deg;
    }
    for(int i=0; i<=common_deg; i++)
    poly.coeff[i] = this->coeff[i] - sub.coeff[i];
    while(poly.coeff[poly.deg] == 0)
    poly.deg--;
    if(poly.deg < 0)
    poly.deg = 0;
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator*(const polynomial_ring<max_deg, p> & mult) const
{
    int tmp = 1, max_deg_idx = 0;
    polynomial_ring<max_deg, p> poly;
    poly.coeff = new finite_field<p>[max_deg+1];
    if(this -> deg == 0 && this -> coeff[0] == 0 || mult.deg == 0 && mult.coeff[0] == 0)
    {poly.coeff[0] = 0; poly.deg = 0; return poly;}
    poly.deg = this->deg + mult.deg;
    for(int i=0; i<=poly.deg; i++)
    poly.coeff[i] = 0;
    for(int i=0; i<=this->deg; i++)
    {
        for(int j=0; j<=mult.deg; j++)
        poly.coeff[i+j] += this -> coeff[i]*mult.coeff[j];
    }
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator%(const polynomial_ring<max_deg, p> & mod) const
{
    finite_field<p> inv_b;
    finite_field<p> *tmp;
    polynomial_ring<max_deg, p> poly;
    poly.coeff = new finite_field<p>[max_deg+1];
    if(this -> deg < mod.deg)
    {
        poly.deg = this -> deg;
        for(int i=0; i<=poly.deg; i++)
        poly.coeff[i] = this -> coeff[i];
        return poly;
    }
    poly.deg = mod.deg-1;
    inv_b = finite_field<p>(1)/mod.coeff[mod.deg];
    tmp = new finite_field<p>[mod.deg+1];
    for(int i=0; i<=(this->deg); i++)
    poly.coeff[i] = this -> coeff[i];
    for(int i=0; i<=mod.deg; i++)
    tmp[i] = inv_b*mod.coeff[mod.deg-i];
    for(int i=this->deg; i>=mod.deg; i--)
    {
        for(int j=mod.deg; j>=0; j--)
        poly.coeff[i-j] -= tmp[j]*poly.coeff[i];
    }
    while(poly.coeff[poly.deg] == 0)
    poly.deg--;
    if(poly.deg < 0)
    poly.deg = 0;
    delete tmp;
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator/(const polynomial_ring<max_deg, p> & div) const
{
    finite_field<p> inv_b;
    finite_field<p> *tmp1, *tmp2;
    polynomial_ring<max_deg, p> poly;
    poly.coeff = new finite_field<p>[max_deg+1];
    if(this -> deg < div.deg)
    {
        poly.deg = 0; poly.coeff[0] = 0;
        return poly;
    }
    poly.deg = this -> deg - div.deg;
    inv_b = finite_field<p>(1)/div.coeff[div.deg];
    tmp1 = new finite_field<p>[this->deg+1];
    tmp2 = new finite_field<p>[div.deg+1];
    for(int i=0; i<=this->deg; i++)
    tmp1[i] = this -> coeff[i];
    for(int i=0; i<=div.deg; i++)
    tmp2[i] = inv_b*div.coeff[div.deg-i];
    for(int i=this->deg; i>=div.deg; i--)
    {
        poly.coeff[i-div.deg] = inv_b*tmp1[i];
        for(int j=div.deg; j>=0; j--)
        tmp1[i-j] -= tmp2[j]*tmp1[i];
    }
    delete tmp1; delete tmp2;
    return poly;
}

template<int max_deg, int p>
finite_field<p>* polynomial_ring<max_deg, p>::to_ary(int* d) const
{
    finite_field<p>* ary = new finite_field<p>[this->deg+1];
    *d = this->deg;
    for(int i=0; i<=this->deg; i++)
    ary[i] = this->coeff[i];
    return ary;
}

template<int max_deg, int p>
polynomial_ring<max_deg, p> mod_pow(polynomial_ring<max_deg, p> base, unsigned long long int idx, polynomial_ring<max_deg, p> & mod)
{
    finite_field<p> c[1] = {finite_field<p>(1)};
    polynomial_ring<max_deg, p> r(0, c);
    while(idx > 0)
    {
        if(idx % 2 == 1)
        {r = base*r; r = r % mod;}
        base = base*base; base = base % mod; idx >>= 1;
    }
    return r;
}

template<int max_deg, int p>
bool is_prime(polynomial_ring<max_deg, p> poly)
{
    finite_field<p> c[2] = {finite_field<p>(0), finite_field<p>(1)};
    polynomial_ring<max_deg, p> x(1, c), pow_x = x;
    int primes[20], num_prime = 0;
    int n = poly.get_deg(), tmp = n;
    if(n == 0)
    return false;
    if(n == 1)
    return true;
    for(int i=2; i*i<=tmp; i++)
    {
        if(tmp % i == 0)
        {
            primes[num_prime] = i; num_prime++;
            while(tmp % i == 0)
            tmp /= i;
        }
    }
    if(tmp > 1)
    {primes[num_prime] = tmp; num_prime++;}
    for(int i=0; i<n/primes[num_prime-1]; i++)
    pow_x = mod_pow(pow_x, p, poly);
    if(!(gcd(pow_x-x, poly).is_unit()))
    return false;
    for(int i=num_prime-1; i>0; i--)
    {
        for(int j=n/primes[i]; j<n/primes[i-1]; j++)
        pow_x = mod_pow(pow_x, p, poly);
        if(!(gcd(pow_x-x, poly).is_unit()))
        return false;
    }
    for(int i=n/primes[0]; i<n; i++)
    pow_x = mod_pow(pow_x, p, poly);
    return (pow_x == x);
}

template<int max_deg>
bool is_prime(polynomial_ring<max_deg, 2> poly)
{
    finite_field<2> c[2] = {finite_field<2>(0), finite_field<2>(1)};
    polynomial_ring<max_deg, 2> x(1, c), pow_x = x;
    int primes[20], num_prime = 0;
    int n = poly.get_deg(), tmp = n;
    if(n == 0)
    return false;
    if(n == 1)
    return true;
    for(int i=2; i*i<=tmp; i++)
    {
        if(tmp % i == 0)
        {
            primes[num_prime] = i; num_prime++;
            while(tmp % i == 0)
            tmp /= i;
        }
    }
    if(tmp > 1)
    {primes[num_prime] = tmp; num_prime++;}
    for(int i=0; i<n/primes[num_prime-1]; i++)
    {pow_x = pow_x*pow_x; pow_x = pow_x % poly;}
    if(gcd(pow_x-x, poly) != 1)
    return false;
    for(int i=num_prime-1; i>0; i--)
    {
        for(int j=n/primes[i]; j<n/primes[i-1]; j++)
        {pow_x = pow_x*pow_x; pow_x = pow_x % poly;}
        if(gcd(pow_x-x, poly) != 1)
        return false;
    }
    for(int i=n/primes[0]; i<n; i++)
    {pow_x = pow_x*pow_x; pow_x = pow_x % poly;}
    return (pow_x == x);
}