#include <iostream>
#include "../header/ntt.hpp"

typedef long long int ll;

ll pow(ll base, ll idx, ll mod)
{
    ll base_tmp = base % mod;
    ll idx_tmp = idx;
    ll result=1;
    
    while(idx_tmp > 0)
    {
        if(idx_tmp % 2 == 1)
        {
            result *= base_tmp;
            result %= mod;
        }
        base_tmp *= base_tmp;
        base_tmp %= mod;
        idx_tmp >>= 1;
    }
    return result;
}

int a[300000];

int main(void)
{
    int n, p, pow_n;
    ll m, tmp;
    arb_ntt<20> arb;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(); std::cout.tie();
    std::cin >> n >> p >> m;
    for(int i=0; i<n; i++)
    std::cin >> a[i];
    pow_n = pow(n, m/2, p);
    if(m % 4 == 1)
    arb(a, n, p, false);
    else if(m % 4 > 1)
    {
        for(int i=1; i<=n/2; i++)
        {tmp = a[i]; a[i] = a[n-i]; a[n-i] = tmp;}
        if(m % 4 == 3)
        arb(a, n, p, false);   
    }
    for(int i=0; i<n; i++)
    a[i] = (1ll*pow_n*a[i]) % p;
    for(int i=0; i<n; i++)
    std::cout << a[i] << ' ';
    std::cout << '\n';
    return 0;
}

