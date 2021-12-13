#include <stdio.h>
#include<gmp.h>
#include <stdbool.h>
#include <time.h>

void lgc (int nbits,int i,mpz_t x) {

    mpz_t m; /*Módulo M determinado pelo número de bits*/
    mpz_init(m);
    mpz_ui_pow_ui(m,2,nbits); /*m=2^nbits*/

    /*gerando uma seed a partir do tempo atual*/
    mpz_t seed;
    mpz_init(seed);
    gmp_randstate_t state;
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomm (seed, state, m);

    mpz_t a,c,temp;
    mpz_init (a);
    mpz_init(c);
    mpz_init(temp); /*valor temporário usado para segurar números grandes*/
    mpz_set_ui(a,6364136223846793005); /*Ṕara esta implementação,Multiplicador= 2^64 */
    mpz_set_ui(c,1); /*Para esta implementação, incremento =1 */

    for (int j=0;j < i;j++) {
        /* a partir de agora,seed será usado como xn e x como xn+1*/
        mpz_mul(temp,seed,a); /* a*x */
        mpz_add(temp,temp,c); /* a*x + c */
        mpz_mod(temp,temp,m); /* a*x + c mod m */
        mpz_set(seed,x); /*x = x n+1 */
        mpz_set(x,temp); /* x n+1 = a*x + c mod m  */

    }
}
int main() {
    mpz_t rng; /*note que como a biblioteca usa estruturas para segurar os valores,não será necessário um valor de retorno*/
    mpz_init(rng); /*rng será usado para apresentar o valor final gerado*/
    /*CONFIGURAR NÚMERO DE BITS NA LINHA ABAIXO(primeiro parametro)*/
    lgc(40,1,rng);
    /*SEGUNDO PARAM= I,número de iterações,pode ser determinado
    como qualquer número ímpar */
    gmp_printf("%Zd \n",rng);
}
