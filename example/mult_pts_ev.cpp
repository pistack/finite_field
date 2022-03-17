#include <iostream>
#include "../header/multipoint_evaluation.hpp"

const int p = 998244353;
const int k_max = 17;
finite_field<p> b[100000], f[200001];
finite_field<p> space[(k_max+21)<<k_max];
int main(void)
{
    int n, q;
    finite_field<p>* result;
    multipoint_evaluation<k_max, p> m_ev;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(); std::cout.tie();
    m_ev.allocate(space, k_max);
    std::cin >> n >> q;
    for(int i=0; i<=n; i++)
    std::cin >> f[n-i];
    for(int i=0; i<q; i++)
    std::cin >> b[i];
    result = m_ev(f, n, b, q);
    for(int i=0; i<q; i++)
    std::cout << result[i] << '\n';
    return 0;
}