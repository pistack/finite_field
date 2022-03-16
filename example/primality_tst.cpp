#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"
#include <iostream>

const int length = 62;

int poly_str_to_poly(char* poly_str, finite_field<2>* buffer)
{
    int i=0, tmp = 0, deg = 0;
    while(poly_str[i] != '\0')
    {
        if(poly_str[i] != 'x' && poly_str[i] != '^')
        {
            if(poly_str[i+1] == '+' || poly_str[i+1] == '\0')
            {
                tmp = poly_str[i] - '0';
                buffer[tmp] = 1;
                if(deg < tmp)
                deg = tmp;
                if(poly_str[i+1] == '\0')
                break;
                i += 2;
            }
            else
            {
                tmp = 10*(poly_str[i]-'0')+(poly_str[i+1]-'0');
                buffer[tmp] = 1;
                if(deg < tmp)
                deg = tmp;
                if(poly_str[i+2] == '\0')
                break;
                i += 3;
            }
        }
        i++;
    }
    return deg;
}

finite_field<2> buffer[33];
char poly_str[161];


int main(void)
{
    int deg;
    polynomial_ring<length, 2> P;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie();
    for(int i=0; i<=32; i++)
    buffer[i] = 0;
    std::cin >> poly_str;
    deg = poly_str_to_poly(poly_str, buffer);
    P = polynomial_ring<length, 2>(deg, buffer);
    std::cout << is_prime(P) << '\n';
    return 0;
}