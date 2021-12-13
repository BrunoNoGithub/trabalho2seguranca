#include <stdio.h>
#include<gmp.h>
#include <stdbool.h>
#include <time.h>
/*Método que retorna,por recursividade, se o MDC de 2 números é 1 */
bool gdcisone(mpz_t x,mpz_t y) {
    if ( (mpz_cmp_ui(x,0)==0) || ( mpz_cmp_ui(y,0)==0) ){
        return ( (mpz_cmp_ui(x,1)==0) || ( mpz_cmp_ui(y,0)==1) );
    }
    if ( (mpz_cmp(x,y)==0 )) {
        if (mpz_cmp_ui(x,1)==0) {
            return true;
        } else {
            return false;
        }
    };
    if (mpz_cmp(x,y)>0) {
        mpz_sub(x,x,y);
        return gdcisone(x,y);
    }
    mpz_sub(y,y,x);
    return gdcisone(x,y);
    }




bool fermat (mpz_t n,int precision) {
    mpz_t resto; 
    mpz_init(resto);
    mpz_mod_ui(resto,n,2);

    if (mpz_cmp_ui(resto,0)==0)    { /*testar se é impar*/
        return false;
    }

    mpz_t nmenos1,seed,a,rnglimit,a_elevado,limitebit,x,y;
    mpz_init(nmenos1);
    mpz_sub_ui(nmenos1,n,1); /*variavel que guarda n-1*/

    mpz_init(seed); /*seed aleatória*/
    gmp_randstate_t state;
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));


    mpz_init(a); /*guarda a aleatório gerado*/
    mpz_init(a_elevado); /*guardará a^(n-1)*/

    mpz_init(x);
    mpz_init(y);
    /*x e y serão iguais a a e n, mas serão passadas e alterados para o método de MDC*/
    mpz_init(rnglimit);
    mpz_sub_ui(rnglimit,nmenos1,2); /*limite de aleatoridade para a (n-3)*/

    mpz_init(limitebit);
    mpz_ui_pow_ui(limitebit,2,4096);

    for (long j = 0; j < (precision);j=j+1) {
        mpz_urandomm (a, state, rnglimit); /*gera numero aleatório até n-3*/
        mpz_add_ui(a,a,2); /*incrementa 2 ao numb aleatorio gerado,formando range [2,n-1]*/

        mpz_set(x,a);
        mpz_set(y,n);
        if (gdcisone(x,y)){ /*se o mdc de a e n for 1*/
            mpz_powm(a_elevado,a,nmenos1,limitebit); /*a^(n-1)*/
            mpz_mod(resto,a_elevado,n); /*resto=a^(n-1) mod n*/
            if (mpz_cmp_ui(resto,1)!=0){ /*se resto for diferente de 1,é composto*/
                return false;
            }
        } else { /*se mdc não for 1,é composto*/
            return false;
        }
    } /*se chegou aqui, pode ser primo*/
    return true;
}

int main() {
    mpz_t number;
    mpz_init(number);
    mpz_set_ui(number,2333);
    if (fermat(number,100)){
        printf("Primo(?)");
    } else {
        printf("Composto");
    };
    return 0;
}