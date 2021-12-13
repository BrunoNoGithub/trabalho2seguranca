#include <stdio.h>
#include<gmp.h>
#include <stdbool.h>
#include <time.h>



bool miller (mpz_t var) {
    /*var equivale a n*/

    mpz_t resto; /*variável resto para guardar valores */
    mpz_init(resto);
    mpz_mod_ui(resto,var,2); /*resto=var mod 2*/
    /*Conferindo se var é ímpar*/
    if (mpz_cmp_ui(resto,0)==0)    {
        return false;
    }
    /*n-1 será utilizado bastante, então estará numa variável*/
    mpz_t nmenos1;
    mpz_init(nmenos1);
    mpz_sub_ui(nmenos1,var,1);
    /*valor k (expoente)*/
    long k = 0;
    /*valor q inicalizado*/
    mpz_t q;
    mpz_init(q);
    /*flag que será usada para concrolar laço*/
    bool flag = false;
    /*temp1 e temp2 para guardar valores temporarios*/
    mpz_t temp1;
    mpz_init(temp1);
    mpz_t temp2;
    mpz_init(temp2);
    /*teto de valores= 2^4096*/
    mpz_t ceiling;
    mpz_init(ceiling);
    mpz_ui_pow_ui(ceiling,2,4096);

    /*laço para achar valores de k e q*/
    while (!flag) {
        k=k+1; /*para cada loop,incrementa k
        primeiro loop->k=1*/
        mpz_ui_pow_ui(temp1,2,k); /* temp1=2^k */
        mpz_mod(temp2,nmenos1,temp1); /* temp2=(n-1) mod (2^k) */
        if (mpz_cmp_ui(temp2,0)!=0) { /*se o resultado da divisão for não inteiro */
            /*se o resultado da divisão por 2^k for não inteiro,então o k correto será
            o valor de k anterior (decrementa o k) e q será o resultado da divisão*/
            flag=true; /*determina fim do loop*/
            k=k-1; /*decrementa o k*/
            mpz_ui_pow_ui(temp1,2,k); /*temp1=2^k*/
            mpz_div(q,nmenos1,temp1); /*q=(n-1)/2^k*/
        }
        }
    /*a=valor aleatório*/
    mpz_t a;
    mpz_init(a);
    mpz_sub_ui(temp2,nmenos1,3); /*temp2=n-4*/

    gmp_randstate_t state;
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));

    mpz_urandomm (a, state, temp2); /*a em alcance 0 a n-4*/


    mpz_add_ui(a,a,2); /*desloca alcance para 2 a n-2*/


    mpz_powm(temp1,a,q,ceiling); /*temp1=a^q */
    mpz_mod(temp2,temp1,var); /*temp2= a^q mod n */
    if ((mpz_cmp_ui(temp2,1))==0 || (mpz_cmp_ui(temp2,-1))==0) { /*se a^q mod n =1...*/
        return true; /*pode ser primo*/
    }
    for (long j = 0; j < (k);j=j+1) {   /*para j =0 até k-1*/
        mpz_pow_ui(temp1,a,2); 
        mpz_pow_ui(temp1,temp1,j);
        mpz_powm(temp1,temp1,q,ceiling); /*temp1=a^q2j*/
        mpz_mod(temp2,temp1,var); /*temp2 = a^q2j mod n*/
        if (mpz_cmp(temp2,nmenos1)==0) { /* se a^2qj mod n = n-1 */
            return true;
        }
    }
    return false;
    }




int main () {
    mpz_t n;
    mpz_init(n);
    /*Inserir valor a ser testado abaixo*/
    mpz_set_ui(n, 2333);
    if (miller(n)) {
        printf("Talvez Primo");
    }else {
        printf("Composto");
    }
}