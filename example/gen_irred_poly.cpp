#include <random>
#include <iostream>
#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"

const int max_deg = 2047;
finite_field<2> coeff[max_deg+1];
int main(void)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> dist(0, 2);
    int deg; polynomial_ring<max_deg, 2> poly;
    int count = 0;
    std::cin >> deg;
    do
    {
        count++;
        coeff[0] = 1;
        for(int i=1; i<deg; i++)
        coeff[i] = dist(gen);
        coeff[deg] = 1;
        poly = polynomial_ring<max_deg, 2>(deg, coeff);
    }
    while(!is_prime(poly));
    std::cout << poly << '\n';
}

