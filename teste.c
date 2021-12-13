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




















int main(){
    bool mil= false;
    mpz_t numb;
    mpz_init(numb);
    bool prime=false;
    int nbits=2048;
    if (mil) {
        while (!prime){
            lgc(nbits,1,numb);
            prime=miller(numb);
            }
    } else {
        while (!prime){
            lgc(nbits,1,numb);
            prime=fermat(numb,2);
        }
    }
    gmp_printf("%i : %Zd",nbits,numb);
}