#include <stdio.h>
#include<gmp.h>
#include <stdbool.h>
#include <time.h>

void xorshift (int nbits,mpz_t rng) {


    /*limite de valores representáveis em nbits (Teto)*/
    mpz_t ceiling;
    mpz_init(ceiling);
    mpz_ui_pow_ui(ceiling,2,nbits);

    /*gerando uma seed a partir do tempo atual*/
    mpz_t seed;
    mpz_init(seed);
    gmp_randstate_t state;
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomm (seed, state, ceiling);



    mpz_set(rng,seed); /*rng (valor grande aleatório) inicialmente recebe valor da seed*/
    mpz_t temp1; /*temp1 guardara numeros grandes temporariamente */
    mpz_init(temp1);


    /*seed ^=seed<<13;seed ^=seed>>17;seed ^=seed<<5;*/

    mpz_mul_2exp (temp1,rng,13); /*temp1=rng<<13*/
    mpz_powm(rng,rng,temp1,ceiling); /*rng=rng^(rng<<13)*/


    mpz_tdiv_q_2exp (temp1,rng,17);/*temp1=rng>>17*/
    mpz_powm(rng,rng,temp1,ceiling);/*rng=rng^(rng>>17)*/

    mpz_mul_2exp (temp1,rng,5);/*temp1=rng<<5*/
    mpz_powm(rng,rng,temp1,ceiling);/*rng=rng^(rng<<5)*/
    return;
}


int main () {
    mpz_t rng; /*estrutura de tamanho arbitrário a ser usada para
 representar numero aleatorio gerado*/
    mpz_init(rng);
    mpz_set_ui(rng,1);
    while (mpz_cmp_ui(rng,1)==0) {
    /*CONFIGURAR NÚMERO DE BITS NA LINHA ABAIXO(primeiro param)*/
        xorshift(40,rng);
    }
    gmp_printf("%Zd \n",rng);
    return 0;
}