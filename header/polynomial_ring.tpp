template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(int d, finite_field<p>* c)
{
    deg = d;
    for(int i=0; i<=deg; i++)
    coeff[i] = c[i];
}

template<int max_deg, int p>
polynomial_ring<max_deg, p>::polynomial_ring(const polynomial_ring<max_deg, p> & poly)
{
    deg = poly.deg;
    for(int i=0; i<=deg; i++)
    coeff[i] = poly.coeff[i];
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator<(const polynomial_ring<max_deg, p> & poly)
{
    if(deg < poly.deg)
    return true;
    if(poly.deg < deg)
    return false;
    for(int i=deg; i>=0; i--)
    {
        if(coeff[i]<poly.coeff[i])
        return true;
        if(poly.coeff[i]<coeff[i])
        return false;
    }
    return false;
}

template<int max_deg, int p>
bool polynomial_ring<max_deg, p>::operator==(const polynomial_ring<max_deg, p> & poly)
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
polynomial_ring<max_deg , p> & polynomial_ring<max_deg, p>::operator=(const polynomial_ring<max_deg, p> & assgin)
{
    if(this == &assgin)
    return;
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
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator-(const polynomial_ring<max_deg, p> & sub) const
{
    polynomial_ring<max_deg, p> poly;
    int common_deg;
    if(this->deg > sub.deg)
    {
        poly.deg = this->deg;
        for(int i=poly.deg; i>add.deg; i--)
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
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator*(const polynomial_ring<max_deg, p> & mult) const
{
    polynomial_ring<max_deg, p> poly;
    poly.deg = mult.deg + this->deg;
    for(int i=0; i<=poly.deg; i++)
    poly.coeff[i] = 0;
    for(int i=0; i<=this->deg; i++)
    {
        for(int j=0; j<=mult.deg; j++)
        poly.coeff[i+j] += this->coeff[i]+mult.coeff[j];
    }
    return poly;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator%(const polynomial_ring<max_deg, p> & mod) const
{
    finite_field<p> inv_coeff;
    finite_field<p> reduce[max_deg+1];
    polynomial_ring<max_deg, p> rem;
    inv_coeff = 1/mod.coeff[mod.deg];
    for(int i=0; i<=mod.deg; i--)
    reduce[i] = inv_coeff*mod.coeff[mod.deg-i];
    for(int i=0; i<=this->deg; i++)
    rem[i] = this->coeff[i];
    for(int i=this->deg; i>=mod.coeff; i--)
    {
        for(int j=mod.coeff; j>=0; j--)
        rem[i-j] -= reduce[j]*rem[i];
    }
    rem.deg = (this->deg > mod.deg ? mod.deg : this->deg);
    while(rem.coeff[rem.deg] == 0)
    rem.deg--;
    return rem;
}

template<int max_deg, int p>
const polynomial_ring<max_deg, p> polynomial_ring<max_deg, p>::operator/(const polynomial_ring<max_deg, p> & div) const
{
    int pos;
    finite_field<p> inv_coeff;
    finite_field<p> reduce[max_deg+1], rem[max_deg+1];
    polynomial_ring<max_deg, p> q;
    if(this->deg < div->deg)
    {q.deg = 0; q.coeff[0] = 0; return q;}
    inv_coeff = 1/div.coeff[div.deg];
    for(int i=0; i<=div.deg; i--)
    reduce[i] = inv_coeff*div.coeff[div.deg-i];
    for(int i=0; i<=this->deg; i++)
    rem[i] = this->coeff[i];
    q.deg = this -> deg - div.deg;
    pos = q.deg;
    for(int i=this->deg; i>=mod.coeff; i--)
    {
        q[pos] = inv_coeff*rem[i];
        for(int j=mod.coeff; j>=0; j--)
        rem[i-j] -= reduce[j]*rem[i];
        pos--;
    }
    return q;
}