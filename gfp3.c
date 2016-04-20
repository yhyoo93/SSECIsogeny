/*
    Copyright (c) 2011 Luca De Feo.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
  Modifed by Dieter Fishbein
*/

#include <pthread.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <openssl/rand.h>
#include <sys/time.h>
#include <string.h>
#define DEBUG
#define GF_TMP_REGS 9
#define MAX(a,b) (((a)>(b))? (a):(b))

#include "keccak.c"
#include <stdarg.h>
#include <limits.h>


/****************** TYPES *****************/

struct GF_params;
typedef struct GF_params GF_params;
int count=1;
int dieter=0;
int version=2;
int count1=0;
int addCount=0;
int count2=0;
int count3=0;
int count2iso=0;
int count2mul=0;
int count3iso=0;
int count3mul=0;
int countMul = 0;

// Elements of GF(p^2)
typedef struct {
  GF_params* parent;
  mpz_t a, b;
} GF;


// The field GF(p^2)
// basically its characteristic and some work registers
struct GF_params {
  mpz_t p, tmp1, tmp2, tmp3;
  GF GFtmp[GF_TMP_REGS];
  gmp_randstate_t state;
  int initialized;
};

void add_GF(GF *res, const GF x, const GF y);

//threading stuff
typedef struct {
    GF *a, *x, *y;
    long int *u;
} arith;

void *add_t(void *ref){
    arith *arithInfo;
    arithInfo = (arith *) ref;
    add_GF(arithInfo->a, *(arithInfo->x), *(arithInfo->y));
}

/****some global variables for threading*****/

typedef struct queue_point queue_point;

typedef struct mt_info {
    
    queue_point *qp2, *qp3;
    
}mt_info;

mt_info *t;
pthread_t threads[2];



/******** IMPLEMENTATION OF GF(p^2) ********/

// Memory management
void init_GF(GF* x, GF_params* parent) {
  mpz_init(x->a);
  mpz_init(x->b);
  x->parent = parent;
}

void clear_GF(GF *x) {
  mpz_clear(x->a);
  mpz_clear(x->b);
}

// Initialization of GF(p,2)
int setup_GF(GF_params* field, const char* characteristic) {
  if (!characteristic) {
    // If p is NULL, use 2^387 * 3^242 - 1 
    // as default value
    mpz_init_set_ui(field->p, 3);
    mpz_pow_ui(field->p, field->p, 242);
    mpz_mul_2exp(field->p, field->p, 387);
    mpz_sub_ui(field->p, field->p, 1);
  } else {
  mpz_init( field->p );
    mpz_set_str(field->p, characteristic, 0);
  }
  // Check that the Legendre symbol of -1 is -1 (p = 3 mod 4)
  if (mpz_fdiv_ui(field->p, 4) != 3) {
    mpz_clear(field->p);
    return 0;
  }

  gmp_randinit_default(field->state);
  mpz_init(field->tmp1); mpz_init(field->tmp2); mpz_init(field->tmp3);
  int i;
  for (i = 0 ; i < GF_TMP_REGS ; i++)
    init_GF(&field->GFtmp[i], field);
  field->initialized = 1;
  
  return 1;
}

void free_GF(GF_params* field) {
  if (field->initialized) {
    int i;
    for (i = 0 ; i < GF_TMP_REGS ; i++)
      clear_GF(&field->GFtmp[i]);
    mpz_clear(field->p); mpz_clear(field->tmp1); mpz_clear(field->tmp2); mpz_clear(field->tmp3);
    field->initialized = 0;
  }
}

// IO of elements
// outputs are strings in base 16
// (to save some bytes)
void set_GF(GF* x, const char* a, const char* b) {
  mpz_set_str(x->a, a, 0);
  mpz_set_str(x->b, b, 0);
}

void set_GFc(GF* x, const char* p) {
  mpz_set_str(x->parent->p, p, 0);  
}


void get_GF(char *a, char *b, const GF x) {
  gmp_sprintf(a, "%#Zx", x.a);
  gmp_sprintf(b, "%#Zx", x.b);
}


// Arithmetic modulo X^2 + 1

/*
  There seems to be a bug in GMP 4.2.1 that makes mpz_mod give
  unpredictable results when the mpz_t holding the result is the same
  as one of the operands.
*/

void copy_GF(GF* res, const GF x){
  mpz_set(res->a, x.a);
  mpz_set(res->b, x.b);
  res->parent = x.parent;
}

void add_GF(GF *res, const GF x, const GF y) {


 /*   printf("\n***** %d ********\n", addCount);
    gmp_printf("\nx.a: %Zd\n", x.a);
    gmp_printf("\ny.a: %Zd\n", y.a);
    printf("\nx.a size: %d\n", (x.a)->_mp_size);
    printf("\ny.a size: %d\n", (y.a)->_mp_size);
*/
    mpz_add(x.parent->tmp1, x.a, y.a);
    mpz_mod(res->a, x.parent->tmp1, x.parent->p);
/*
    gmp_printf("\nres: %Zd\n", res->a);
    gmp_printf("\nres: %Zx\n", res->a);
    printf("\nres size: %d\n", (res->a)->_mp_size);
    printf("\n******************\n", addCount);
    
      /*  gmp_printf("\nres: %Zd\n", res->a);
    printf("\nres size: %d\n", (res->a)->_mp_size);
    printf("\n******************\n", addCount);    
    */
    mpz_add(x.parent->tmp1, x.b, y.b);
    mpz_mod(res->b, x.parent->tmp1, x.parent->p);
    res->parent = x.parent;
    
/*    printf("\n\nInput/output after execution count: %d\n", countAdd);
    print_GF(x,"\nx");
    print_GF(y,"\ny");
    print_GF(*res,"\nres");
    gmp_printf("\nres.a: %Zx\nres.b: %Zx\n\n", (res->a), res->b);
    printf("\n*****************************\n");
    
    
    */
    
    addCount++;
}

void add_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
/*    mpz_t tmp1;
    mpz_init(tmp1); */ 
    
    mpz_add(tmp[0], x.a, y.a);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_add(tmp[0], x.b, y.b);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    
  //  mpz_clear(tmp1);
}

void add_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_add_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}

void sub_GF(GF *res, const GF x, const GF y) {
  mpz_sub(x.parent->tmp1, x.a, y.a);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_sub(x.parent->tmp1, x.b, y.b);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

//specify mpz tmp variables in tmp
void sub_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
 //   mpz_t tmp1;
 //   mpz_init(tmp1);    
    
    mpz_sub(tmp[0], x.a, y.a);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_sub(tmp[0], x.b, y.b);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    
  //  mpz_clear(tmp1);
    
}

void sub_GF_ui(GF *res, const GF x, unsigned long int u) {
  mpz_sub_ui(x.parent->tmp1, x.b, u);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  mpz_set(res->a, x.a);
  res->parent = x.parent;
}
//specify tmps
void sub_GF_ui_t(GF *res, const GF x, unsigned long int u, mpz_t *tmp) {
    
   // mpz_t tmp1;
   // mpz_init(tmp1);
    
    mpz_sub_ui(tmp[0], x.b, u);
    mpz_mod(res->b, tmp[0], x.parent->p);
    mpz_set(res->a, x.a);
    res->parent = x.parent;
    
 //   mpz_clear(tmp1);
}

void neg_GF(GF *res, const GF x) {
  if (mpz_sgn(x.a) == 0)
    mpz_set(res->a, x.a);
  else
     mpz_sub(res->a, x.parent->p, x.a); 
                                
  if (mpz_sgn(x.b) == 0)
    mpz_set(res->b, x.b);
  else
    mpz_sub(res->b, x.parent->p, x.b);
  res->parent = x.parent;
}

void scalar_GF(GF *res, const GF x, mpz_t s) {
  mpz_mul(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void scalar_GF_si(GF *res, const GF x, long int s) {
  mpz_mul_si(x.parent->tmp1, x.a, s);
  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mul_si(x.parent->tmp1, x.b, s);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);
  res->parent = x.parent;
}

void scalar_GF_si_t(GF *res, const GF x, long int s, mpz_t *tmp) {

    mpz_mul_si(tmp[0], x.a, s);
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mul_si(tmp[0], x.b, s);
    mpz_mod(res->b, tmp[0], x.parent->p);
    res->parent = x.parent;
    

}


void mul_GF(GF *res, const GF x, const GF y) {
  
 /*   printf("\n****** %d ********", countMul) ;
    printf("\nInputs to mul_GF ") ; 
    gmp_printf("\n x.a: %Zd", x.a);
    gmp_printf("\n x.b: %Zd", x.b);
    gmp_printf("\n y.a: %Zd", y.a);  
    gmp_printf("\n y.b: %Zd\n", y.b);  
    printf("\nInputs to mpz_mul ", countMul) ;   
*/
      
  mpz_add(y.parent->tmp1, x.a, x.b);
  mpz_sub(y.parent->tmp2, y.b, y.a);
    
   // printf("\n*************** %d **********\n", countMul);
  /*  gmp_printf("\n tmp1: %Zd\n", y.parent->tmp1);
    gmp_printf("\n tmp2: %Zd\n", y.parent->tmp2);
    printf("\n tmp1 size: %d\n", (y.parent->tmp1)->_mp_size);
    printf("\n tmp2 size: %d\n", (y.parent->tmp2)->_mp_size);
   */
    
    
    
  mpz_mul(y.parent->tmp3, y.parent->tmp1, y.parent->tmp2);
    
 /*   gmp_printf("\n res: %Zd\n", y.parent->tmp3);
    printf("\n res (hex): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", ((y.parent->tmp3)->_mp_d)[25], ((y.parent->tmp3)->_mp_d)[24], ((y.parent->tmp3)->_mp_d)[23], ((y.parent->tmp3)->_mp_d)[22], ((y.parent->tmp3)->_mp_d)[21], ((y.parent->tmp3)->_mp_d)[20], ((y.parent->tmp3)->_mp_d)[19], ((y.parent->tmp3)->_mp_d)[18], ((y.parent->tmp3)->_mp_d)[17], ((y.parent->tmp3)->_mp_d)[16], ((y.parent->tmp3)->_mp_d)[15], ((y.parent->tmp3)->_mp_d)[14], ((y.parent->tmp3)->_mp_d)[13], ((y.parent->tmp3)->_mp_d)[12], ((y.parent->tmp3)->_mp_d)[11], ((y.parent->tmp3)->_mp_d)[10], ((y.parent->tmp3)->_mp_d)[9], ((y.parent->tmp3)->_mp_d)[8], ((y.parent->tmp3)->_mp_d)[7], ((y.parent->tmp3)->_mp_d)[6], ((y.parent->tmp3)->_mp_d)[5], ((y.parent->tmp3)->_mp_d)[4], ((y.parent->tmp3)->_mp_d)[3], ((y.parent->tmp3)->_mp_d)[2], ((y.parent->tmp3)->_mp_d)[1], ((y.parent->tmp3)->_mp_d)[0]);

    printf("\n res size: %d\n", (y.parent->tmp3)->_mp_size);
   // printf("\n*****************************\n");
   */ 
  
  //  printf("\n*************** %d **********\n", countMul);
 /*   gmp_printf("\n x.a: %Zd\n", x.a);
    gmp_printf("\n y.b: %Zd\n", y.b);
    
    gmp_printf("\n x.a (h): %Zx\n", x.a);
    gmp_printf("\n y.b (h): %Zx\n", y.b);
    
    printf("\n x.a size: %d\n", (x.a)->_mp_size);
    printf("\n y.b size: %d\n", (y.b)->_mp_size);
   */ 
    
    
 /*   printf("\n*************** %d **********\n", countMul);
    gmp_printf("\n x.a: %Zx\n", x.a);
    gmp_printf("\n y.b: %Zx\n", y.b);
    
    gmp_printf("\n y.a (h): %Zx\n", y.a);
    gmp_printf("\n x.b (h): %Zx\n", x.b);
    
    printf("\n x.a size: %d\n", (x.a)->_mp_size);
    printf("\n y.b size: %d\n", (y.b)->_mp_size);
    printf("\n y.a size: %d\n", (y.a)->_mp_size);
    printf("\n x.b size: %d\n", (x.b)->_mp_size);
*/
    
  mpz_mul(y.parent->tmp1, x.a, y.b);
    
  //  gmp_printf("\n res: %Zd\n", y.parent->tmp1);
    // printf("\n res (hex): %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx %lx\n\n", ((y.parent->tmp3)->_mp_d)[25], ((y.parent->tmp3)->_mp_d)[24], ((y.parent->tmp3)->_mp_d)[23], ((y.parent->tmp3)->_mp_d)[22], ((y.parent->tmp3)->_mp_d)[21], ((y.parent->tmp3)->_mp_d)[20], ((y.parent->tmp3)->_mp_d)[19], ((y.parent->tmp3)->_mp_d)[18], ((y.parent->tmp3)->_mp_d)[17], ((y.parent->tmp3)->_mp_d)[16], ((y.parent->tmp3)->_mp_d)[15], ((y.parent->tmp3)->_mp_d)[14], ((y.parent->tmp3)->_mp_d)[13], ((y.parent->tmp3)->_mp_d)[12], ((y.parent->tmp3)->_mp_d)[11], ((y.parent->tmp3)->_mp_d)[10], ((y.parent->tmp3)->_mp_d)[9], ((y.parent->tmp3)->_mp_d)[8], ((y.parent->tmp3)->_mp_d)[7], ((y.parent->tmp3)->_mp_d)[6], ((y.parent->tmp3)->_mp_d)[5], ((y.parent->tmp3)->_mp_d)[4], ((y.parent->tmp3)->_mp_d)[3], ((y.parent->tmp3)->_mp_d)[2], ((y.parent->tmp3)->_mp_d)[1], ((y.parent->tmp3)->_mp_d)[0]);
 //   printf("\n res size: %d\n", (y.parent->tmp1)->_mp_size);
    //printf("\n*****************************\n");
    
        
  mpz_mul(y.parent->tmp2, y.a, x.b);
    
/*    gmp_printf("\n\n tmp1: %Zx\n", y.parent->tmp1);
    gmp_printf("\n tmp2: %Zx\n", y.parent->tmp2);
    
    printf("\n tmp1 size: %d\n", (y.parent->tmp1)->_mp_size);
    printf("\n tmp2 size: %d\n", (y.parent->tmp2)->_mp_size);
    printf("\n******************************\n", countMul);
 
*/     

  mpz_sub(y.parent->tmp3, y.parent->tmp3, y.parent->tmp1);
    


    
       
  mpz_add(y.parent->tmp1, y.parent->tmp1, y.parent->tmp2); 
    
        mpz_mod(res->a, y.parent->tmp1, y.parent->p);
    
      mpz_add(y.parent->tmp3, y.parent->tmp3, y.parent->tmp2);
    
        mpz_mod(res->b, y.parent->tmp3, y.parent->p);
    
      
  res->parent = x.parent;
    
  /*  printf("\n\n Outputs from mul_GF ") ; 
    gmp_printf("\n res.a: %Zd", res->a);
    gmp_printf("\n res.b: %Zd", res->b);
    printf("\n*****************", countMul) ; 
  */
    countMul++;  
}

//mul version where you can specify the temp variables. It will take the tmp veriables in the parent of the fourth object
void mul_GF_t(GF *res, const GF x, const GF y, mpz_t *tmp) {
    
   // mpz_t tmp1, tmp2, tmp3;
   // mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
    
    mpz_add(tmp[0], x.a, x.b);  
    
    mpz_sub(tmp[1], y.b, y.a);   
  
    mpz_mul(tmp[2], tmp[0], tmp[1]);
    
    mpz_mul(tmp[0], x.a, y.b);
    mpz_mul(tmp[1], y.a, x.b);
    
    mpz_sub(tmp[2], tmp[2], tmp[0]);
    mpz_add(tmp[2], tmp[2], tmp[1]);
    mpz_add(tmp[0], tmp[0], tmp[1]);
    
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mod(res->b, tmp[2], x.parent->p);
    res->parent = x.parent;
    
 //   mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); 
    
}

void sqr_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.b);
    
 //   if(dieter==1)  
 //       gmp_printf("\n\ntmp1_sqr: %Zd\n\n", x.parent->tmp1);    
    
  mpz_add(x.parent->tmp1, x.parent->tmp1, x.parent->tmp1);
   
  

  mpz_add(x.parent->tmp2, x.b, x.a);
  mpz_sub(x.parent->tmp3, x.b, x.a);
  mpz_mul(x.parent->tmp2, x.parent->tmp2, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp1, x.parent->p);
  mpz_mod(res->b, x.parent->tmp2, x.parent->p);
  res->parent = x.parent;
}

//specify the GF that gives you your mpz temp variables
void sqr_GF_t(GF *res, const GF x, mpz_t *tmp) {
    
  //  mpz_t tmp1, tmp2, tmp3;
   // mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
    
    mpz_mul(tmp[0], x.a, x.b);
    mpz_add(tmp[0], tmp[0], tmp[0]);
    
    mpz_add(tmp[1], x.b, x.a);
    mpz_sub(tmp[2], x.b, x.a);
    mpz_mul(tmp[1], tmp[1], tmp[2]);
    
    mpz_mod(res->a, tmp[0], x.parent->p);
    mpz_mod(res->b, tmp[1], x.parent->p);
    res->parent = x.parent;
    
   // mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); 
}


int inv_GF(GF *res, const GF x) {
  mpz_mul(x.parent->tmp1, x.a, x.a);
  mpz_addmul(x.parent->tmp1, x.b, x.b);
  if (!mpz_invert(x.parent->tmp3, x.parent->tmp1, x.parent->p))
    return 0;

  mpz_mul(x.parent->tmp1, x.b, x.parent->tmp3);
  mpz_neg(x.parent->tmp3, x.parent->tmp3);
  mpz_mul(x.parent->tmp2, x.a, x.parent->tmp3);

  mpz_mod(res->a, x.parent->tmp2, x.parent->p);
  mpz_mod(res->b, x.parent->tmp1, x.parent->p);

  res->parent = x.parent;
  return 1;
}

int div_GF(GF *res, const GF x, const GF y) {
 if (!inv_GF(&x.parent->GFtmp[0], y)) return 0;
  mul_GF(res, x, x.parent->GFtmp[0]);
  return 1;
}

// Miscellaneaous
int cmp_GF(const GF x, const GF y) {
  int c = mpz_cmp(x.a, y.a);
  if (c == 0) c = mpz_cmp(x.b, y.b);
  return c;
}

int is_one_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_cmp_ui
                                 (x.b, 1) == 0);
}

int is_zero_GF(const GF x) {
  return (mpz_sgn(x.a) == 0) && (mpz_sgn(x.b) == 0);
}

void random_GF(GF *res) {
  mpz_urandomm(res->a, res->parent->state, res->parent->p);
  mpz_urandomm(res->b, res->parent->state, res->parent->p);
}

void print_GF(const GF x, char * a) {
  gmp_printf("%s: %Zd*x + %Zd\n",a, x.a, x.b);
}

//returns the characteristic of the field containing the given GF element
char * getChar( GF *res) {
  char * sa1= malloc( mpz_sizeinbase(res->parent->p, 10)+10 );
    mpz_get_str(sa1, 10, res->parent->p);
  return sa1;
}


char * getA(GF *x) {
    char * num;
    num  = malloc(sizeof(char)*1000);
    mpz_get_str(num,10,(*x).a);
    return num;
}

char * getB(  GF *x ) {
    char * num;
    num  = malloc(sizeof(char)*1000);
    mpz_get_str(num,10,(*x).b);
    return num;
}

//Returns 1 if a=b and another int otherwise
int equals(GF *a, GF *b){
    GF tmp1;
    
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &tmp1, parent );
    sub_GF(&tmp1, *a, *b);
    return is_zero_GF(tmp1);
    
}


/***** ALGORITHMS RELATED TO ECC **********/

// One step of Montgomery ladder
void mont_ladder(GF *res1x, GF *res1z,
     GF *res2x, GF *res2z,
     const GF x1, const GF z1,
     const GF x2, const GF z2,
     const GF dx, const GF dz,
     const GF A24) {
  GF* tmp = x1.parent->GFtmp;

  add_GF(&tmp[4], x1, z1);         // a = (self.x + self.z)
  sub_GF(&tmp[5], x1, z1);         // b = (self.x - self.z)
  sub_GF(&tmp[6], x2, z2);
  add_GF(&tmp[7], x2, z2);
  sqr_GF(&tmp[1], tmp[4]);         // aa = a.square()
  sqr_GF(&tmp[3], tmp[5]);         // bb = b.square()
  sub_GF(&tmp[0], tmp[1], tmp[3]); // e = aa - bb
  mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
  mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b

  add_GF(&tmp[2], tmp[6], tmp[7]);
  sqr_GF(&tmp[2], tmp[2]);
  mul_GF(&tmp[2], tmp[2], dz);     // x2 = diff.z*(da + cb).square()

  sub_GF(&tmp[8], tmp[6], tmp[7]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(&tmp[8], tmp[8], dx);     // z2 = diff.x*(da - cb).square()

  mul_GF(res1z, A24, tmp[0]);
  add_GF(res1z, *res1z, tmp[3]);
  mul_GF(res1z, *res1z, tmp[0]);   // z1 = e*(bb + self.curve.A24*e))
  mul_GF(res1x, tmp[1], tmp[3]);  // x1 = aa*bb

  copy_GF(res2x, tmp[2]);
  copy_GF(res2z, tmp[8]);
}

/* Montgomery point doubling */
void mont_double(GF *resx, GF *resz,
     const GF x, const GF z,
     const GF A24) {

  GF* tmp = x.parent->GFtmp;
  add_GF(&tmp[0], x, z);           // a = (x + z)
  sqr_GF(&tmp[1], tmp[0]);         // aa = a^2
  sub_GF(&tmp[2], x, z);           // b = (x - z)
  sqr_GF(&tmp[3], tmp[2]);         // bb = b^2
  sub_GF(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  mul_GF(resz, A24, tmp[4]);   
  add_GF(resz, *resz, tmp[3]);      
  mul_GF(resz, *resz, tmp[4]);   // Z = c (bb + A24 c)) 
  mul_GF(resx, tmp[1], tmp[3]);  // X = aa bb  
}

/* Montgomery point tripling */
void mont_triple(GF *resx, GF *resz,
     const GF x, const GF z,
     const GF A24) {
  GF* tmp = x.parent->GFtmp;

  // Very dirty function, assuming that mont_double uses
  // registers 0 and 2 to store resp. x+z and x-z
  mont_double(&tmp[5], &tmp[6], x, z, A24);

  sub_GF(&tmp[7], tmp[5], tmp[6]);
  add_GF(&tmp[8], tmp[5], tmp[6]);
  mul_GF(&tmp[5], tmp[7], tmp[0]); // da = (x2 - z2)*a
  mul_GF(&tmp[6], tmp[8], tmp[2]); // cb = (x2 + z2)*b

  add_GF(&tmp[7], tmp[5], tmp[6]);
  sqr_GF(&tmp[7], tmp[7]);
  mul_GF(&tmp[7], tmp[7], z);     // X = z*(da + cb)^2

  sub_GF(&tmp[8], tmp[5], tmp[6]);
  sqr_GF(&tmp[8], tmp[8]);
  mul_GF(resz, tmp[8], x);     // Z = x*(da - cb)^2

  copy_GF(resx, tmp[7]);
}

/*
  Converts a Montgomery point to an Edwards point.
  Guarantee: avoids the temporary registers 0-4.
*/
void mont_to_ed(GF* Rx, GF* Ry,
    const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  /*
    X = x(x+z) / y(x+z)
    Y = y(x-z) / y(x+z)
  */
  add_GF(&tmp[5], x, z);
  sub_GF(&tmp[6], x, z);
  mul_GF(&tmp[7], y, tmp[5]);
  inv_GF(&tmp[8], tmp[7]);
  mul_GF(&tmp[7], tmp[5], tmp[8]);
  mul_GF(Rx, x, tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[8]);
  mul_GF(Ry, y, tmp[7]);
}

typedef struct {
    GF tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, dz, dx, x, z;
} thread_p1;

/*void *p3(void *thr){

    thread_p1 *t = (thread_p1 *) thr;
    
   // print_GF(t->tmp2,"tmp2");
    
    /* P3 */
  //  mul_GF(&(t->tmp2), t->tmp0, t->tmp4); // da = (P.x - P.z)*a
/*    mul_GF(&(t->tmp3), t->tmp1, t->tmp5); // cb = (P.x + P.z)*b
    
    add_GF(&(t->tmp8), t->tmp2, t->tmp3);
    sqr_GF(&(t->tmp0), t->tmp8);
    mul_GF(&(t->x), t->tmp0, t->dz);     // x2 = diff.z*(da + cb).square()
   
  /*  sub_GF(&(t->tmp8), t->tmp2, t->tmp3);
    sqr_GF(&(t->tmp0), t->tmp8);
    mul_GF(&(t->z), t->tmp0, t->dx);     // z2 = diff.x*(da - cb).square()
  */  
    
 //  pthread_exit(0);
//}

// Three-point ladder addition step:
//   P1 = 2 P1
//   P2 = dadd(P1, P2, D2)
//   P3 = dadd(P1, P3, D3)
void mont_tradd(GF *x1, GF *z1,
    GF *x2, GF *z2,
    GF *x3, GF *z3,
    const GF dx2, const GF dz2,
    const GF dx3, const GF dz3,
    const GF A24) {
        dieter++;
    
    /*   GF_params *parent1, *parent2, *parent3;
     char * prime = malloc(sizeof(char)*5000);
     mpz_get_str (prime, 10, (x1->parent)->p);
     parent1 = malloc(sizeof(GF_params));
     parent2 = malloc(sizeof(GF_params));
     parent3 = malloc(sizeof(GF_params));
     setup_GF(parent1, prime);
     setup_GF(parent2, prime);
     setup_GF(parent3, prime);
     
     GF* tmpB = parent1->GFtmp;
     GF* tmpC = parent2->GFtmp;*/
    GF* tmp = x1->parent->GFtmp;
    
   /* if(dieter==version){ 
        print_GF(*x3,"\n\nx3\n\n"); 
        print_GF(*z3,"\n\nz3\n\n"); 
        print_GF(*x1,"\n\nx1\n\n"); 
        print_GF(*z1,"\n\nz1\n\n");
        print_GF(*x2,"\n\nx2\n\n"); 
        print_GF(*z2,"\n\nz2\n\n");
    }
    */
    
    // The use of temporary registers in this function is tailored
    // so that it is compatible with
    //       {x1, x2} == {&tmp[0], &tmp[2]}
    //       {z1, z2} == {&tmp[1], &tmp[3]}
    // so that mont_3ladder may safely call it.
    // Be careful when you change indices!
    
    add_GF(&tmp[4], *x1, *z1);         // a = (self.x + self.z)
    sub_GF(&tmp[5], *x1, *z1);         // b = (self.x - self.z)
    sub_GF(&tmp[6], *x2, *z2);
    add_GF(&tmp[7], *x2, *z2);
    sub_GF(&tmp[0], *x3, *z3);
    add_GF(&tmp[1], *x3, *z3); 
    
   /*if(count1==1600){  
        
        print_GF(tmp[6],"\n\ntmp6");
        print_GF(tmp[4],"\n\ntmp4");
        print_GF(tmp[7],"\n\ntmp7");
        print_GF(tmp[5],"\n\ntmp5");
    }*/
    
    // print_GF(tmp[1], "\n\ntmp1");
    
  //  if(dieter==version)
  //      print_GF((tmp[4]),"\n\ntmp4_bt\n\n"); 
    
    /*  copy_GF(&tmpB[4],tmpA[4]);
     copy_GF(&tmpB[5],tmpA[5]);
     
     copy_GF(&tmpC[4],tmpA[4]);
     copy_GF(&tmpC[5],tmpA[5]);
     
     pthread_t thread1;
     thread_p1 thread_params = {tmpA[0], tmpA[1], tmpA[2], tmpA[3], tmpA[4], tmpA[5], tmpA[6], tmpA[7], tmpA[8], dz3, dx3, *x3, *z3}; 
     int rc = pthread_create(&thread1, NULL, p3, (void *) &thread_params); */ 
    
    /* P3 */
    
   /* if(dieter==version){   
        print_GF((tmp[0]),"\n\ntmp0\n\n");
        print_GF((tmp[4]),"\n\ntmp4\n\n");
        print_GF((tmp[1]),"\n\ntmp1\n\n");
        print_GF((tmp[5]),"\n\ntmp5\n\n");
    }*/
    
    mul_GF(&tmp[2], tmp[0], tmp[4]); // da = (P.x - P.z)*a
    mul_GF(&tmp[3], tmp[1], tmp[5]); // cb = (P.x + P.z)*b
    
/*    if(dieter==version){   
        print_GF(tmp[2],"\n\ntmp2\n\n");
        print_GF(tmp[3],"\n\ntmp3\n\n");
    }  */
    
    add_GF(&tmp[8], tmp[2], tmp[3]);
    sqr_GF(&tmp[0], tmp[8]);
    mul_GF(x3, tmp[0], dz3);     // x2 = diff.z*(da + cb).square()
    
    sub_GF(&tmp[8], tmp[2], tmp[3]);
    sqr_GF(&tmp[0], tmp[8]);
    mul_GF(z3, tmp[0], dx3);     // z2 = diff.x*(da - cb).square()
    
    /* P2 */
    mul_GF(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
    mul_GF(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b
    
    add_GF(&tmp[2], tmp[6], tmp[7]);
    sqr_GF(&tmp[3], tmp[2]);
    mul_GF(x2, tmp[3], dz2);     // x2 = diff.z*(da + cb).square()
    
    sub_GF(&tmp[3], tmp[6], tmp[7]);
    sqr_GF(&tmp[8], tmp[3]);
    mul_GF(z2, tmp[8], dx2);     // z2 = diff.x*(da - cb).square()
    
    
    /* P1 */
    sqr_GF(&tmp[6], tmp[4]);         // aa = a.square()
    sqr_GF(&tmp[7], tmp[5]);         // bb = b.square()
    sub_GF(&tmp[8], tmp[6], tmp[7]); // e = aa - bb
    mul_GF(&tmp[4], A24, tmp[8]);
    add_GF(&tmp[5], tmp[4], tmp[7]);
    mul_GF(z1, tmp[5], tmp[8]);      // z1 = e*(bb + self.curve.A24*e))
    mul_GF(x1, tmp[6], tmp[7]);      // x1 = aa*bb
    
     count1++;
    
 /*   if(count1==1000){ 
         printf("\n********VALUES TO CHECK***********%d\n",count1);
        print_GF(*(x1), "\n\nt.x1");
        print_GF(*(z1), "\n\nt.z1");
        print_GF(*(x2), "\n\nt.x2");
        print_GF(*(z2), "\n\nt.z2");
        print_GF(*x3,"\n\nt.x3\n\n");
        print_GF(*z3,"\n\nt.z3\n\n"); 
        printf("\n*****************************\n");

        
    }*/
   
   
}


// 3-point ladder to compute P + [t]Q
// Inputs: t, P, Q, Q - P
void mont_3ladder(GF* Rx, GF* Rz,
                  const mpz_t t,
                  const GF Px, const GF Pz,
                  const GF Qx, const GF Qz,
                  const GF QPx, const GF QPz,
                  const GF A24) {
    GF* tmp = Px.parent->GFtmp;
    set_GF(&tmp[0], "0", "1"); 
    set_GF(&tmp[1], "0", "0");
    copy_GF(&tmp[2], Qx); copy_GF(&tmp[3], Qz);
    copy_GF(Rx, Px); copy_GF(Rz, Pz); 
    int bit = mpz_sizeinbase(t, 2) - 1;
    for ( ; bit >=0 ; bit--) {
        if (mpz_tstbit(t, bit) == 0) {
            mont_tradd(&tmp[0], &tmp[1], &tmp[2], &tmp[3], Rx, Rz,
                       Qx, Qz, Px, Pz, A24);
        } else {
            mont_tradd(&tmp[2], &tmp[3], &tmp[0], &tmp[1], Rx, Rz,
                       Qx, Qz, QPx, QPz, A24);
        }
    }
}




/*void *EdwardsCompute(void *thread){
    
    thread_params *thread1;
    thread1 = (thread_params *) thread; 
    
    init_GF(thread1->aPx, (thread1->field)); init_GF(thread1->aPy, (thread1->field));
    mont_to_ed(thread1->aPx, thread1->aPy, *(thread1->Px), *(thread1->Py), *(thread1->Pz));
    init_GF(thread1->aQx, (thread1->field)); init_GF(thread1->aQy, (thread1->field));
    mont_to_ed(thread1->aQx, thread1->aQy, *(thread1->Qx), *(thread1->Qy), *(thread1->Qz));
    
    mul_GF((thread1->tmp4), *(thread1->aPx), *(thread1->aQx)); // tmp4 = C = aPx * aQx
    mul_GF((thread1->tmp5), *(thread1->aPy), *(thread1->aQy)); // tmp5 = D = aPy * aQy
    add_GF((thread1->tmp0), *(thread1->aPx), *(thread1->aPy)); // tmp0 = A = aPx + aPy
    add_GF((thread1->tmp2), *(thread1->aQx), *(thread1->aQy)); // tmp2 = B = aQx + aQy
    mul_GF((thread1->tmp7), *(thread1->tmp4), *(thread1->tmp5));
    mul_GF((thread1->tmp6), *(thread1->d), *(thread1->tmp7)); // tmp6 = E = d * aPx * aQx * aPy * aQy
    sqr_GF((thread1->tmp8), *(thread1->tmp6));
    neg_GF((thread1->tmp7), *(thread1->tmp8));
  /*  add_GF_ui(tmp[7], tmp[7], 1);
    inv_GF(tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
    add_GF_ui(tmp[6], tmp[6], 1);
    mul_GF(tmp[7], a, tmp[4]);
    sub_GF(tmp[1], tmp[5], tmp[7]);
    mul_GF(tmp[7], tmp[6], tmp[1]);
    init_GF(PQy, field);
    mul_GF(PQy, tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
    neg_GF(tmp[6], tmp[6]);
    add_GF_ui(tmp[6], tmp[6], 2);
    mul_GF(tmp[1], tmp[0], tmp[2]);
    sub_GF(tmp[3], tmp[1], tmp[4]);
    sub_GF(tmp[1], tmp[3], tmp[5]);
    mul_GF(tmp[7], tmp[6], tmp[1]);
    init_GF(PQx, field);
    mul_GF(PQx, tmp[7], tmp[8]);*/ // PQx = (1-E)(A B - C - D) / (1-E^2)

//}

/*
  Computes [m]P + [n]Q, with P and Q points on the Montgomery curve
  with parameters A,B.  Uses Edwards' coordinates for
  calculations.  */
void shamir(GF* Rx, GF* Ry, GF* Rz,
      const GF A, const GF B,
      const GF Px, const GF Py, const GF Pz,
      const GF Qx, const GF Qy, const GF Qz,
      const mpz_t m, const mpz_t n) {
  
  
   
  // some temporary registers
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;
  // some other dynamically allocated registers 
  GF a, d, aPx, aPy, aQx, aQy, PQx, PQy;
  
 //  thread_params thread1 = {&aPx, &aPy, &aQx, &aQy, &Px, &Py, &Pz, &Qx, &Qy, &Qz, &field, &d, &a, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7], &tmp[8]};
    
 // pthread_t threads;  
 // int rc = pthread_create(&threads, NULL, EdwardsCompute, (void *) &thread1);   
  
  /* 
     Parameters of the Edwards curve equivalent to this one:
       a = (A+2)/B
       d = (A-2)/B
  */
  inv_GF(&tmp[0], B);
  copy_GF(&tmp[1], A);
  add_GF_ui(&tmp[1], tmp[1], 2);
  init_GF(&a, field);
  mul_GF(&a, tmp[1], tmp[0]);
  sub_GF_ui(&tmp[1], tmp[1], 4);
  init_GF(&d, field);
  mul_GF(&d, tmp[1], tmp[0]);

  /*
    Computing the Ewdards coordinates of P and Q:
      aPx, aPy = Edwards(P)
      aQx, aQy = Edwards(Q)
  */
  init_GF(&aPx, field); init_GF(&aPy, field);
  mont_to_ed(&aPx, &aPy, Px, Py, Pz);
  init_GF(&aQx, field); init_GF(&aQy, field);
  mont_to_ed(&aQx, &aQy, Qx, Qy, Qz);

  /*
    Computing P+Q using affine Edwards.
  */
  mul_GF(&tmp[4], aPx, aQx); // tmp4 = C = aPx * aQx
  mul_GF(&tmp[5], aPy, aQy); // tmp5 = D = aPy * aQy
  add_GF(&tmp[0], aPx, aPy); // tmp0 = A = aPx + aPy
  add_GF(&tmp[2], aQx, aQy); // tmp2 = B = aQx + aQy
  mul_GF(&tmp[7], tmp[4], tmp[5]);
  mul_GF(&tmp[6], d, tmp[7]); // tmp6 = E = d * aPx * aQx * aPy * aQy
  sqr_GF(&tmp[8], tmp[6]);
  neg_GF(&tmp[7], tmp[8]);
  add_GF_ui(&tmp[7], tmp[7], 1);
  inv_GF(&tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
  add_GF_ui(&tmp[6], tmp[6], 1);
  mul_GF(&tmp[7], a, tmp[4]);
  sub_GF(&tmp[1], tmp[5], tmp[7]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQy, field);
  mul_GF(&PQy, tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
  neg_GF(&tmp[6], tmp[6]);
  add_GF_ui(&tmp[6], tmp[6], 2);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  sub_GF(&tmp[3], tmp[1], tmp[4]);
  sub_GF(&tmp[1], tmp[3], tmp[5]);
  mul_GF(&tmp[7], tmp[6], tmp[1]);
  init_GF(&PQx, field);
  mul_GF(&PQx, tmp[7], tmp[8]); // PQx = (1-E)(A B - C - D) / (1-E^2)*/
  
  int bit = MAX(mpz_sizeinbase(m, 2), mpz_sizeinbase(n, 2)) - 1;
  mpz_set_ui(Rx->a, 0); mpz_set_ui(Ry->a, 0); mpz_set_ui(Rz->a, 0);   
  mpz_set_ui(Rx->b, 0); mpz_set_ui(Ry->b, 1); mpz_set_ui(Rz->b, 1);
  Rx->parent = Ry->parent = Rz->parent = Px.parent;

  for ( ; bit >=0 ; bit--){
    /* Double, using projective Edwards */
    add_GF(&tmp[1], *Rx, *Ry);
    sqr_GF(&tmp[0], tmp[1]); // tmp0 = B = (Rx + Ry)^2
    sqr_GF(&tmp[1], *Rx); // tmp1 = C = Rx^2
    sqr_GF(&tmp[2], *Ry); // tmp2 = D = Ry^2
    mul_GF(&tmp[3], a, tmp[1]); // tmp3 = E = a C
    add_GF(&tmp[4], tmp[3], tmp[2]); // tmp4 = F = E + D
  sqr_GF(&tmp[5], *Rz); // tmp5 = H = Rz^2
    scalar_GF_si(&tmp[7], tmp[5], 2);
    sub_GF(&tmp[6], tmp[4], tmp[7]); // tmp6 = J = F - 2H
    sub_GF(&tmp[7], tmp[0], tmp[1]);
    sub_GF(&tmp[8], tmp[7], tmp[2]);
    mul_GF(Rx, tmp[8], tmp[6]); // Rx = (B-C-D) J
    sub_GF(&tmp[7], tmp[3], tmp[2]);
    mul_GF(Ry, tmp[7], tmp[4]); // Ry = (E-D) F
    mul_GF(Rz, tmp[4], tmp[6]); // Rz = F J

    /* Double and Add, using projective Edwards */
    int r = mpz_tstbit(m, bit) | (mpz_tstbit(n, bit) << 1);
    if (r) {
      if (r == 1) {
  mul_GF(&tmp[0], *Rx, aPx); // tmp0 = C = Rx aPx
  mul_GF(&tmp[1], *Ry, aPy); // tmp1 = D = Ry aPy
  add_GF(&tmp[2], aPx, aPy);  // tmp2 = H = aPx + aPy
      } else if (r == 2) {
  mul_GF(&tmp[0], *Rx, aQx); // tmp0 = C = Rx aQx
  mul_GF(&tmp[1], *Ry, aQy); // tmp1 = D = Ry aQy
  add_GF(&tmp[2], aQx, aQy);  // tmp2 = H = aQx + aQy
      } else {
  mul_GF(&tmp[0], *Rx, PQx); // tmp0 = C = Rx PQx
  mul_GF(&tmp[1], *Ry, PQy); // tmp1 = D = Ry PQy
  add_GF(&tmp[2], PQx, PQy);  // tmp2 = H = PQx + PQy
      }
      sqr_GF(&tmp[3], *Rz); // tmp3 = B = Rz^2
      mul_GF(&tmp[5], tmp[0], tmp[1]);
      mul_GF(&tmp[4], d, tmp[5]); // tmp4 = E = d C D
      sub_GF(&tmp[5], tmp[3], tmp[4]); // tmp5 = F = B - E
      add_GF(&tmp[6], tmp[3], tmp[4]); // tmp6 = G = B + E
      add_GF(&tmp[7], *Rx, *Ry);
      mul_GF(&tmp[8], tmp[7], tmp[2]);
      sub_GF(&tmp[7], tmp[8], tmp[0]);
      sub_GF(&tmp[8], tmp[7], tmp[1]);
      mul_GF(&tmp[7], tmp[8], tmp[5]);
      mul_GF(Rx, *Rz, tmp[7]); // Rx = Rz F ((Rx+Ry)H - C - D)
      mul_GF(&tmp[7], a, tmp[0]);
      sub_GF(&tmp[8], tmp[1], tmp[7]);
      mul_GF(&tmp[7], tmp[6], tmp[8]);
      mul_GF(Ry, *Rz, tmp[7]); // Ry = Rz G (D - a C)
      mul_GF(Rz, tmp[5], tmp[6]); // Rz = F G
    }
  }

  /* Convert to Montgomery */
  add_GF(&tmp[0], *Rz, *Ry);
  sub_GF(&tmp[1], *Rz, *Ry);
  mul_GF(Ry, tmp[0], *Rz); // Ry = (Rz+Ry)Rz
  mul_GF(Rz, tmp[1], *Rx); // Rz = (Rz-Ry)Rx
  mul_GF(Rx, tmp[0], *Rx); // Rx = (Rz+Ry)Rx

  clear_GF(&a); clear_GF(&d);
  clear_GF(&aPx); clear_GF(&aPy);
  clear_GF(&aQx); clear_GF(&aQy);
  clear_GF(&PQx); clear_GF(&PQy);
}




/************* ISOGENIES ******************/
typedef struct {
  GF u, r;
} iso;

typedef GF iso2;

typedef struct {
  GF p, p2;
} iso3;

void init_iso3( iso3 *iso, GF_params *parent){
 
    init_GF(&(iso->p), parent);
    init_GF(&(iso->p2), parent);
    
}

void copy_iso3( iso3 *res, iso3 iso){
    
    copy_GF(&(res->p), iso.p);
    copy_GF(&(res->p2), iso.p2);
    
}

typedef struct {
  GF Ap2;
} iso4;

/* Utility routine to compute (A+2)/4 */
void a24(GF* A24, const GF A) {
  GF_params* field = A.parent;
  GF* tmp = field->GFtmp;

  add_GF_ui(&tmp[0], A, 2);
  mpz_set_ui(field->tmp1, 4);
  mpz_invert(field->tmp2, field->tmp1, field->p);
  scalar_GF(A24, tmp[0], field->tmp2);
}

/*
  Compute an isomorphism of the montgomery curve
  sending (x,z) to (0,0).
*/
void isom_comp(iso* iso, GF* iA, GF* iB, GF* iA24,
         const GF A, const GF B, const GF A24,
         const GF x, const GF z) {
  GF* tmp = A.parent->GFtmp;
  
  mont_double(&tmp[1], &tmp[2], x, z, A24);
  div_GF(&tmp[0], tmp[1], tmp[2]);  // P2x = x([2]P) / z([2]P)
  neg_GF(&iso->r, tmp[0]);  // r = -P2x
  scalar_GF_si(&tmp[1], tmp[0], 3);
  add_GF(&tmp[1], tmp[1], A);  // a2 = 3 P2x + A
  mul_GF(&tmp[2], iso->r, z);
  add_GF(&tmp[3], tmp[2], x);
  div_GF(&iso->u, z, tmp[3]);  // u = z / (x - z P2x)
  mul_GF(iA, tmp[1], iso->u);  // iA = a2 u
  mul_GF(iB, B, iso->u);  // iB = B u
  a24(iA24, *iA);
}

/* Apply an isomorphism of Montgomery curves */
void isom_apply(GF* X, GF* Y, GF* unused,
    const iso iso,
    const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], iso.r, z);
  add_GF(&tmp[1], x, tmp[0]);
  if (Y)
    mul_GF(Y, y, iso.u); // Y = y u
  mul_GF(X, tmp[1], iso.u); // X = (x + r z) u
}


/*
  Compute a 2-isogeny of the montgomery curve
  sending (x,z) to (1,...).
*/
void iso2_comp(iso2* iso, GF* iA, GF* iB, GF* iA24,
         const GF A, const GF B,
         const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  sub_GF(&tmp[0], x, z);
  sqr_GF(&tmp[1], tmp[0]);
  inv_GF(&tmp[0], tmp[1]);
  mul_GF(&tmp[1], tmp[0], z);
  mul_GF(iso, tmp[1], x); // iA2 = x z / (x-z)^2
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(iB, B, *iso); // iB = B iA2
  mul_GF(iA, tmp[0], *iso); // iA = (A+6) iA2
  a24(iA24, *iA);
}

/* Apply a 2-isogeny of Montgomery curves */
void iso2_apply(GF* X, GF* Y, GF* Z,
    const iso2 iso,
    const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;
  
  sub_GF(&tmp[3], x, z);
  sqr_GF(&tmp[4], tmp[3]);
  if (Y) {
    mul_GF(&tmp[4], x, tmp[4]); // ... X = x iA2 (x - z)^2
    sqr_GF(&tmp[0], x); // Px2 = x^2
    sqr_GF(&tmp[1], z);
    sub_GF(&tmp[2], tmp[0], tmp[1]);
    mul_GF(&tmp[1], y, tmp[2]);
    mul_GF(Y, iso, tmp[1]); // Y = iA2 y (x^2 - z^2)
    mul_GF(Z, z, tmp[0]); // Z = z x^2
  } else {
    mul_GF(Z, z, x); // Z = z x
  }
  mul_GF(X, iso, tmp[4]); // X = iA2 (x - z)^2
    
    count2iso++;
}

/*
  Compute a 3-isogeny of the montgomery curve
*/
void iso3_comp(iso3* iso, GF* iA, GF* iB, GF* iA24,
         const GF A, const GF B,
       const GF x, const GF z) {
  GF* tmp = x.parent->GFtmp;

  div_GF(&iso->p, x, z);             // p
  sqr_GF(&iso->p2, iso->p);          // p^2

  scalar_GF_si(&tmp[3], iso->p, -6); 
  add_GF(&tmp[4], tmp[3], A);
  mul_GF(&tmp[3], tmp[4], iso->p);
  add_GF_ui(&tmp[4], tmp[3], 6);     // (-6p + A)p + 6
  
  mul_GF(iB, B, iso->p2);      // iB = B p^2
  mul_GF(iA, tmp[4], iso->p);  // iA = ((-6p + A)p + 6)p

  a24(iA24, *iA);
}

/* Apply a 3-isogeny of Montgomery curves */
void iso3_apply(GF* X, GF* Y, GF* Z,
    const iso3 iso,
    const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], z, iso.p);
  sub_GF(&tmp[1], x, tmp[0]); // h = x - p z
                              // if zero, P is in the kernel
  mul_GF(&tmp[2], x, iso.p);
  sub_GF(&tmp[0], tmp[2], z); // rh = x p - z
  sqr_GF(&tmp[3], tmp[0]);
  mul_GF(&tmp[2], x, tmp[3]); // X0 = x (x p - z)^2

  if (Y) {
    mul_GF(&tmp[2], tmp[2], tmp[1]); // X0 *= h
    mul_GF(&tmp[3], x, z);
    sub_GF_ui(&tmp[4], iso.p2, 1);
    mul_GF(&tmp[5], tmp[3], tmp[4]);
    scalar_GF_si(&tmp[3], tmp[5], -2);
    mul_GF(&tmp[4], tmp[0], tmp[1]);
    add_GF(&tmp[5], tmp[3], tmp[4]);
    mul_GF(&tmp[3], tmp[5], tmp[0]); // (rh (rh h + 2xz(1-p^2)))
    sqr_GF(&tmp[7], tmp[1]);
    mul_GF(&tmp[8], tmp[7], z);
    mul_GF(Y, y, tmp[3]); // Y = y (rh (rh h + 2xz(1-p^2)))
    mul_GF(Z, tmp[8], tmp[1]); // Z = h^2 h z
  } else {
    sqr_GF(&tmp[3], tmp[1]);
    mul_GF(Z, tmp[3], z); // Z = h^2 z
  }
  copy_GF(X, tmp[2]);
    
  count3iso++;
}



/*
  Compute a 4-isogeny of the Montgomery curve
  sending (1,...) to infinity.
*/
void iso4_comp(iso4* iso, GF* iA, GF* iB, GF* iA24,
         const GF A, const GF B) {
  GF* tmp = A.parent->GFtmp;

  add_GF_ui(&iso->Ap2, A, 2); // Ap2 = A + 2
  sub_GF_ui(&tmp[0], A, 2);
  neg_GF(&tmp[0], tmp[0]);
  inv_GF(&tmp[2], tmp[0]); // iAm2 = 1 / (2-A)
  add_GF_ui(&tmp[0], A, 6);
  mul_GF(&tmp[1], tmp[0], tmp[2]);
  mul_GF(iB, B, tmp[2]); // iB = B iAm2
  scalar_GF_si(iA, tmp[1], -2); // iA = -2 (A+6) iAm2
  a24(iA24, *iA);
}

/* Apply a 4-isogeny of Montgomery curves */
void iso4_apply(GF* X, GF* Y, GF* Z,
    const iso4 iso,
    const GF x, const GF y, const GF z) {
  GF* tmp = x.parent->GFtmp;

  mul_GF(&tmp[0], x, z); // z1 = x z
  sub_GF(&tmp[2], x, z);
  sqr_GF(&tmp[1], tmp[2]); // x1 = (x - z)^2
  mul_GF(&tmp[2], tmp[0], iso.Ap2); // zA2 = z1 Ap2
  scalar_GF_si(&tmp[3], tmp[0], 4); // fourz = 4 z1
  add_GF(&tmp[6], tmp[1], tmp[2]);
  add_GF(&tmp[4], tmp[1], tmp[3]);
  mul_GF(&tmp[5], tmp[6], tmp[4]); // x0 = (x1+zA2)(x1+fourz)

  if (Y) {
    mul_GF(&tmp[4], x, tmp[1]); // B = x x1
    mul_GF(&tmp[5], tmp[5], tmp[4]); // x0 *= B
    sqr_GF(&tmp[6], z);
    sub_GF(&tmp[7], tmp[0], tmp[6]);
    scalar_GF_si(&tmp[6], tmp[7], 2);
    add_GF(&tmp[7], tmp[6], tmp[1]); // C = x1 + 2(z1 - z^2)
    sqr_GF(&tmp[6], tmp[1]);
    mul_GF(&tmp[8], tmp[2], tmp[3]);
    sub_GF(&tmp[6], tmp[6], tmp[8]); // D = x1^2 - zA2 fourz
    mul_GF(&tmp[8], tmp[7], tmp[6]);
    mul_GF(Y, y, tmp[8]); // Y = y C D
    sqr_GF(&tmp[6], tmp[4]);
    sub_GF_ui(&tmp[7], iso.Ap2, 4);
      
    //The code below implements neg_GF(&tmp[7], tmp[7]). For some reason just calling the function here has unpredictable behavious due to a bug in mpz_sub
      if (mpz_sgn(tmp[7].a) == 0)
          mpz_set(tmp[7].a, tmp[7].a);
      else
        mpz_sub(tmp[7].a, tmp[7].parent->p, tmp[7].a); 
                                    
      if (mpz_sgn(tmp[7].b) == 0)
          mpz_set(tmp[7].b, x.b);
      else
          mpz_sub(tmp[7].b, tmp[7].parent->p, tmp[7].b);
      
    mul_GF(&tmp[8], tmp[6], tmp[7]);
    mul_GF(Z, z, tmp[8]); // Z = z B^2 (4 - Ap2)
  } else {
    sub_GF(&tmp[4], tmp[3], tmp[2]);
    mul_GF(Z, tmp[1], tmp[4]); // Z = x1 (fourz - zA2)
  }
  copy_GF(X, tmp[5]);
}

/* Apply a 4-isogeny of Montgomery curves (hopefully thread safe since we can specify temp variables */
 void iso4_apply_t(GF* X, GF* Y, GF* Z,
                 const iso4 iso,
                 const GF x, const GF y, const GF z, GF * tmp) {
     
     mpz_t tmpA[3];
     mpz_init(tmpA[0]);
     mpz_init(tmpA[1]);
     mpz_init(tmpA[2]);
     
    // dieter=1; 
    mul_GF_t(&tmp[0], x, z, tmpA); // z1 = x z
    sub_GF_t(&tmp[2], x, z, tmpA);
    sqr_GF_t(&tmp[1], tmp[2], tmpA); // x1 = (x - z)^2  //probably the buggy line
     
  /*   print_GF(tmp[1], "tmp1");
     print_GF(*Y, "Y");
     print_GF(*Z, "Z");
     print_GF(x, "x");
     print_GF(y, "y");
     print_GF(z, "z");  */
     
    mul_GF_t(&tmp[2], tmp[0], iso.Ap2, tmpA); // zA2 = z1 Ap2
    scalar_GF_si_t(&tmp[3], tmp[0], 4, tmpA); // fourz = 4 z1
    add_GF_t(&tmp[6], tmp[1], tmp[2], tmpA);
    add_GF_t(&tmp[4], tmp[1], tmp[3], tmpA);
    mul_GF_t(&tmp[5], tmp[6], tmp[4], tmpA); // x0 = (x1+zA2)(x1+fourz)
 
    if (Y) {
        mul_GF_t(&tmp[4], x, tmp[1], tmpA); // B = x x1
        mul_GF_t(&tmp[5], tmp[5], tmp[4], tmpA); // x0 *= B
        sqr_GF_t(&tmp[6], z, tmpA);
        sub_GF_t(&tmp[7], tmp[0], tmp[6], tmpA);
        scalar_GF_si_t(&tmp[6], tmp[7], 2, tmpA);
        add_GF_t(&tmp[7], tmp[6], tmp[1], tmpA); // C = x1 + 2(z1 - z^2)
        sqr_GF_t(&tmp[6], tmp[1], tmpA);
        
      //  printf("\ndieter\n");
        
        mul_GF_t(&tmp[8], tmp[2], tmp[3], tmpA);
        sub_GF_t(&tmp[6], tmp[6], tmp[8], tmpA); // D = x1^2 - zA2 fourz
        mul_GF_t(&tmp[8], tmp[7], tmp[6], tmpA);
        mul_GF_t(Y, y, tmp[8], tmpA); // Y = y C D
        sqr_GF_t(&tmp[6], tmp[4], tmpA);
        sub_GF_ui_t(&tmp[7], iso.Ap2, 4, tmpA);
 
        //The code below implements neg_GF(&tmp[7], tmp[7]). For some reason just calling the function here has unpredictable behavious due to a bug in mpz_sub
        if (mpz_sgn(tmp[7].a) == 0)
            mpz_set(tmp[7].a, tmp[7].a);
        else
            mpz_sub(tmp[7].a, tmp[7].parent->p, tmp[7].a); 
 
        if (mpz_sgn(tmp[7].b) == 0)
            mpz_set(tmp[7].b, x.b);
        else
            mpz_sub(tmp[7].b, tmp[7].parent->p, tmp[7].b);
 
            mul_GF_t(&tmp[8], tmp[6], tmp[7], tmpA);
            mul_GF_t(Z, z, tmp[8], tmpA); // Z = z B^2 (4 - Ap2)
        } else {
            sub_GF_t(&tmp[4], tmp[3], tmp[2], tmpA);
            mul_GF_t(Z, tmp[1], tmp[4], tmpA); // Z = x1 (fourz - zA2)
          }
        copy_GF(X, tmp[5]);
      //  dieter=0;
     
    // mpz_clear(tmpA[0]); mpz_clear(tmpA[1]); mpz_clear(tmpA[2]); 
    }
 
 
/******* COMPOSITE ISOGENIES **************/

/* Implementation of a queue */

typedef struct queue_point {
    GF x, z;
    int h;
    struct queue_point *next, *prev;
};


void copy_QP(queue_point *res, queue_point qp){
   // printf("dieter00\n");
    copy_GF(&(res->x), qp.x); //printf("dieter--\n");
    copy_GF(&(res->z), qp.z); //printf("dieter++\n");
    res->h = qp.h; // printf("dieter3\n");
    res->next = qp.next; //printf("dieter4\n");
    res->prev = qp.prev; //printf("dieter5\n");
    
}

void print_QP(queue_point qp, char * str){
    
    printf("\n%s\n",str);
    print_GF(qp.x, "x");
    print_GF(qp.z, "z");
    printf("\nh: %d", qp.h);
    
}

void init_QP( queue_point *qp, GF_params *parent ){
    //printf("dieter++\n");
    init_GF(&(qp->x), parent);
    init_GF(&(qp->z), parent);
    qp->h=0; //printf("dieter2\n");
    qp->next = malloc(sizeof(queue_point));
    qp->prev = malloc(sizeof(queue_point)); //printf("dieter3\n");
}

typedef struct {
    GF *X, *Y, *Z, x, y, z;
    iso4  d;
} thread_p2;

void * apply_t(void * thr){
    
    thread_p2 *t = (thread_p2 *) thr;
    
    int i=0;
    GF tmp[GF_TMP_REGS];// = malloc(sizeof(GF)*GF_TMP_REGS);
    for (i = 0 ; i < GF_TMP_REGS ; i++)
        init_GF(&tmp[i], t->X->parent);
   
    iso4_apply_t(t->X, t->Y, t->Z, t->d, t->x, t->y, t->z, tmp);
}

#define Q_INIT(q,field) do {       \
    q = malloc(sizeof(queue_point)); \
    if (q) {           \
      q->next = q->prev = NULL;      \
      init_GF(&q->x, field);       \
      init_GF(&q->z, field);       \
      q->h = 0;          \
    }            \
  } while(0)
#define Q_CLEAR(q) do {   \
      clear_GF(&q->x);    \
      clear_GF(&q->z);    \
      free(q);      \
    } while(0)
#define Q_PUSH(tail,q) do {     \
    tail->next = q;       \
    q->prev = tail;       \
    tail = q;         \
  } while(0)
#define Q_POP(tail,q) do {    \
    q = tail;           \
    tail = tail->prev;    \
    if (tail) {           \
      tail->next = NULL;    \
    }             \
  } while(0)
#define Q_NEXT(q) (q->next)
#define Q_PREV(q) (q->prev)
#define Q_ISHEAD(q) (q->prev==NULL)
#define Q_ISTAIL(q) (q->next==NULL)

// These bits of code are almost identical for 1, 2, 3, 4
// isogenies, thus we "template" them.
#define APPLY_ISOG(apply,obj,lower) do {    \
    for ( tmp = tail ; tmp ; tmp = Q_PREV(tmp)) { \
      apply(&tmp->x, NULL, &tmp->z, obj,    \
      tmp->x, tmp->x, tmp->z);  \
      tmp->h = tmp->h - lower;        \
    }       \
    if (Px && Py && Pz){  \
      apply(Px, Py, Pz, obj, *Px, *Py, *Pz);\
     }                           \
    if (Qx && Qy && Qz)         \
      apply(Qx, Qy, Qz, obj, *Qx, *Qy, *Qz);    \
  } while (0)              
#define COMP_ISOG(comp,obj) do {    \
    Q_POP(tail, tmp);         \
    comp(&obj, A, B, A24, *A, *B, tmp->x, tmp->z);  \
    Q_CLEAR(tmp);         \
  } while (0)


union isogenies {
  struct {
    iso d1;
    iso2 d2;
    iso4 d4;
  };
  iso3 d3;
};




void *p2(){
    
}

void *p3(){
    
}
    

/* Push (Px, Py, Pz) through the isogeny of kernel
 generated by (Rx, Rz) using the given strategy. */
void push_through_iso(GF *A, GF *B, GF *A24,
                      const GF Rx, const GF Rz,
                      const int ell, int *strategy, int h,
                      GF *Px, GF *Py, GF *Pz,
                      GF *Qx, GF *Qy, GF *Qz, int e) {
  
    GF_params* field = A->parent;
    int split, i, first = 1, first2=1;
    union isogenies phi;
    queue_point *tail, *tmp;
  /*  count2iso = 0;
    count2mul = 0;
    count3iso = 0;
    count3mul = 0;
    */
    if (ell == 2) {
        init_GF(&phi.d1.u, field);
        init_GF(&phi.d1.r, field);
        init_GF(&phi.d2, field);
        init_GF(&phi.d4.Ap2, field);
    } else {
        init_GF(&phi.d3.p, field);
        init_GF(&phi.d3.p2, field);
    }
    
    Q_INIT(tail, field);
    copy_GF(&tail->x, Rx);
    copy_GF(&tail->z, Rz);
    tail->h = h;
    while (tail) {
        h = tail->h;
        split = strategy[h];
        // Descend to the floor
        while (h > 1) {
            Q_INIT(tmp, field);
            copy_GF(&tmp->x, tail->x);
            copy_GF(&tmp->z, tail->z);
            
            for ( i=0 ; i < h - split ; i++) {    
                if (ell == 2){
                    mont_double(&tmp->x, &tmp->z,
                                tmp->x, tmp->z, *A24); 
                    
                   // count2mul++;
                }
                else{ 
                    mont_triple(&tmp->x, &tmp->z,
                                tmp->x, tmp->z, *A24);
                   // count3mul++;
                }
                
            }
            tmp->h = split;
            
            Q_PUSH(tail, tmp);
            h = split;
            split = strategy[h];
        }
        // For ell=2, at the first iteration, bring the
        // 2-torsion point to (0,0)
        if (ell == 2 && first) {
            first = 0;
            Q_INIT(tmp, field); // slight abuse
            mont_double(&tmp->x, &tmp->z, tail->x, tail->z, *A24);
            isom_comp(&phi.d1, A, B, A24,
                      *A, *B, *A24, tmp->x, tmp->z);
            Q_CLEAR(tmp);
            APPLY_ISOG(isom_apply, phi.d1, 0);
        }
        
        
        // Compute and apply the isogeny
        if (ell == 2) {

            COMP_ISOG(iso2_comp, phi.d2);
            APPLY_ISOG(iso2_apply, phi.d2, 1);
                
        } else {
            COMP_ISOG(iso3_comp, phi.d3);
            
           // print_GF(phi.d3.p, "phi.d3.p"); 
           // print_GF(phi.d3.p2, "phi.d3.p2"); 
            APPLY_ISOG(iso3_apply, phi.d3, 1);
        }
       // count3++;   
        
    }
    // For ell=2 there is still a 4-isogeny to apply
    if (ell == 2) {
        iso4_comp(&phi.d4, A, B, A24, *A, *B);
        // This works because the queue is empty
        APPLY_ISOG(iso4_apply, phi.d4, 2);

    }
    
    if (ell == 2) {
        clear_GF(&phi.d1.u);
        clear_GF(&phi.d1.r);
        clear_GF(&phi.d2);
        clear_GF(&phi.d4.Ap2);
    } else {
        clear_GF(&phi.d3.p);
        clear_GF(&phi.d3.p2);
    } 
    
}





////////////////////dieter////////////////////////////////////////////////

//threading stuff
    

/*    
void *add_ui_t(void *arithInfo){
    add_GF_ui(arithInfo->a, arithInfo->(*x), arithInfo->(*u));
}

void *sub_t(void *arithInfo){
     sub_GF(arithInfo->a, arithInfo->(*x), arithInfo->(*y));
}
    
void *sub_ui(void *arithInfo){
    sub_GF_ui(arithInfo->a, arithInfo->(*x), arithInfo->(*u));
}

void *neg(void *arithInfo){
    neg_GF(arithInfo->a,  arithInfo->(*x));
        
}
    
void *scalar(void *arithInfo){
    ;
        
}

void *scalar_si(void *arithInfo){
    ;   
        
}
    
void *mul(void *arithInfo){
    ;   
        
}
    
void *sqr(void *arithInfo){
    ;   
        
}

void *inv(void *arithInfo){
    ;   
        
}
    
void *div(void *arithInfo){
    ;   
        
}
 
 */
    
//Montgomery Curve
typedef struct {
    GF A, B, A24;
} MC;

void copy_MC(MC *res, MC curve){
    copy_GF( &(*res).A, curve.A );
    copy_GF( &(*res).B, curve.B );
    copy_GF( &(*res).A24, curve.A24 );
}

void init_MC(MC* x) {
    GF_params* parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &(*x).A, parent);
    init_GF( &(*x).B, parent);
    init_GF( &(*x).A24, parent);
}


void clear_MC(MC *x) {
    clear_GF(&(*x).A);
    clear_GF(&(*x).B);
    clear_GF(&(*x).A24);
}

void print_Curve(MC *x){
    print_GF( (*x).A,"A" );
    print_GF( (*x).B,"B" );
    print_GF( (*x).A24,"A24" ); 
}

void set_Curve( MC *curve, GF A, GF B, GF A24 ){
    (*curve).A = A;
    (*curve).B = B;
    (*curve).A24 = A24;
}


void j_invariant(GF *final, MC *curve){
    GF tmp1, tmp2, tmp3, denom, num;
    
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    init_GF( &tmp1, parent );
    init_GF( &tmp2, parent );
    init_GF( &tmp3, parent );
    init_GF( &denom, parent );
    init_GF( &num, parent );
    sqr_GF(&tmp1, (*curve).A);
    sub_GF_ui(&denom, tmp1, 4);
    sub_GF_ui(&tmp2, tmp1, 3 );
    sqr_GF( &tmp3, tmp2);
    mul_GF( &tmp1, tmp3, tmp2);
    scalar_GF_si(&num, tmp1, 256);
    div_GF(final, num, denom);
}


//Montgomery Point
typedef struct {
    GF x,y,z;
    MC curve;
} MP;

void init_MP(MP* a) {
    GF_params *parent;
    MC *curve;
    curve = malloc(sizeof(MC));
    parent = malloc(sizeof(GF_params));
    setup_GF(parent,""); 
    
    init_GF( &(*a).x, parent );
    init_GF( &(*a).y, parent );
    init_GF( &(*a).z, parent );
    init_MC(&a->curve);
}
    
void copy_MP( MP *res, MP P ){
    copy_GF(&(*res).x,P.x);
    copy_GF(&(*res).y,P.y);
    copy_GF(&(*res).z,P.z);  
    copy_MC(&(*res).curve,P.curve);
}
    
void set_Curve_MC( MP *res, MC curve ){
    (*res).curve.A = curve.A;
    (*res).curve.B = curve.B;
    (*res).curve.A24 = curve.A24;
}

void set_MP( MP *res, GF x, GF y, GF z, MC curve ){
    (*res).x = x;
    (*res).y = y;
    (*res).z = z;
    
    (*res).curve.A = curve.A;
    (*res).curve.B = curve.B;
    (*res).curve.A24 = curve.A24;
}


void clear_MP(MP *a) {
    clear_GF(&(*a).x);
    clear_GF(&(*a).y);
    clear_GF(&(*a).z);
}

void print_MP(MP *a, char * string){
    printf("%s: \n",string);
    print_GF( (*a).x,"x" );
    print_GF( (*a).y,"y" );
    print_GF( (*a).z,"z" );
    print_Curve( &((*a).curve) );
}


//P+Q make sure both curves are the same!!!!
void  add(MP *res, MP Q, MP P){

    GF *x, *y, *z;
    
    GF* pxtmp = P.x.parent->GFtmp;
    GF* pytmp = P.y.parent->GFtmp;
    GF* pztmp = P.y.parent->GFtmp;

    x = malloc(sizeof(GF));
    z = malloc(sizeof(GF));
    y = malloc(sizeof(GF));
    
    init_GF(x, P.x.parent);
    init_GF(z, P.x.parent);
    init_GF(y, P.x.parent);
    
    mul_GF( &pxtmp[1], Q.y, P.z );
    mul_GF( &pxtmp[2] , Q.x, P.z );
    
  //  printf("dietergeorge-1\n");
   // print_GF(pxtmp[2],"\npx2");
    
    mul_GF( &pxtmp[3] , Q.z, P.z );

    mul_GF( &pxtmp[4], P.y ,Q.z );
    
  //  printf("dietergeorge-3\n");
  //  print_GF(pxtmp[4],"\npx4");
    
    sub_GF( &pxtmp[5], pxtmp[4], pxtmp[1]);
    sqr_GF( &pxtmp[6], pxtmp[5] );

    mul_GF( &pxtmp[4], P.x, Q.z );
    sub_GF( &pxtmp[7], pxtmp[4], pxtmp[2] );
    
    sqr_GF( &pxtmp[8], pxtmp[7] );
    
   /* printf("dietergeorge\n");
    print_GF(pxtmp[7],"\nx");
    print_GF(pxtmp[8],"\ny");
   */ 
    mul_GF( &pxtmp[0], pxtmp[7], pxtmp[8] );
    mul_GF( &pytmp[2], pxtmp[8], pxtmp[2] );

    mul_GF( &pxtmp[4], Q.curve.B, pxtmp[6] );
    mul_GF( &pytmp[1], Q.curve.A, pxtmp[8] );
    sub_GF( &pytmp[8], pxtmp[4], pytmp[1] );
    mul_GF( &pztmp[1], pytmp[8], pxtmp[3] );
    sub_GF( &pxtmp[4], pytmp[1], pxtmp[0] );
    scalar_GF_si( &pytmp[1], pytmp[2], 2 );
    sub_GF( &pytmp[4], pxtmp[4], pytmp[1] ); 

    mul_GF(x, pxtmp[7], pytmp[4] );
    mul_GF(z, pxtmp[0], pxtmp[3] );

    sub_GF( &pxtmp[4], pytmp[2], pytmp[4] );
    mul_GF( &pytmp[1], pxtmp[5], pxtmp[4] );   
    mul_GF( &pytmp[6], pxtmp[0], pxtmp[1] );
    sub_GF( y, pytmp[1], pytmp[6] );
    
    (*res).x = *x;
    (*res).y = *y;
    (*res).z = *z;
    (*res).curve = P.curve;
    
}

//Returns the negative of P
void neg(MP *res, MP P){
    GF *yN;
    
    yN = malloc(sizeof(GF));
    init_GF( yN, (P).x.parent);
    neg_GF( yN , (P).y );

    (*res).x = (P).x;
    (*res).y = *yN;
    (*res).z = (P).z;
    (*res).curve = P.curve;
    free(yN);
}



//P-Q make sure both curves are the same!!
void  subtract(MP *res, MP P, MP Q){
    MP *point;
    point = malloc(sizeof(MP));
    init_MP(point);
    
    neg(point,Q);
    add(res,P,*point);
    free(point);

}


//returns a random integer from 0 to m (inclusive)
void rand_range(mpz_t *num, mpz_t m){
    
    unsigned long int bytes;
    mpz_t *tmp1, *tmp2;
    tmp1 = malloc(sizeof(mpz_t));
    tmp2 = malloc(sizeof(mpz_t));
    mpz_init(tmp1);
    mpz_init(tmp2);
    bytes = ceil( mpz_sizeinbase(m, 2)/8.0 ) + 20;
    unsigned char buf[bytes];
    RAND_bytes(buf, bytes);
    
    unsigned char buf1[bytes];
    char* buf2 = buf1;

    int k=0;
    for (k; k < bytes; k++)
        buf2 += sprintf(buf2, "%02X", buf[k]);
    
    sprintf(buf2,"\n");
    *(buf2 + 1) = '\0';
    mpz_set_str(tmp1, buf1, 16);
    mpz_add_ui(tmp2, m, 1);
    
    mpz_mod(num, tmp1, tmp2); 
    free(tmp1);free(tmp2);
}



void rand_subgroup(mpz_t *m, mpz_t *n, char * l, char * e){
  
    mpz_t *num, *_l, *_e, *tmp1, *le, *le1, *l1;
    num = malloc(sizeof(mpz_t));
    _l = malloc(sizeof(mpz_t));
    _e = malloc(sizeof(mpz_t));
    tmp1 = malloc(sizeof(mpz_t));
    le = malloc(sizeof(mpz_t));
    le1 = malloc(sizeof(mpz_t));
    l1 = malloc(sizeof(mpz_t));
     
    long int e_int = 0;
    e_int = atoi(e);
    mpz_init(le);
    mpz_init(le1);
    mpz_init(l1);
    mpz_init(tmp1);
    mpz_init(num);
    
    mpz_init_set_str(_l, l, 10);
    mpz_init_set_str(_e, e, 10);
    mpz_add_ui(l1, _l, 1);
    rand_range(num, *l1);
        
    if ( mpz_cmp_ui(*num, 1)==0 ){
        mpz_set_ui(m, 1);
        mpz_pow_ui(le, _l, e_int);
        rand_range(n, *le);
    }else{
        mpz_pow_ui(le1, _l, e_int-1);
        rand_range(tmp1, *le1);
        mpz_mul(m, tmp1, _l);
        mpz_set_ui(n, 1);
    }
    
    free(num);free(_l);free(_e);free(tmp1);free(le);free(le1);free(l1);
    
}

double run_ZKP(double *time, char * eA_str, char * eB_str, char * lA_str, char * lB_str, int *strA, int lenA, int *strB, int lenB, 
                    MP *PA, MP *QA, MP *PB, MP *QB, int rounds, MP **R_array, MP **phiR_array, MP **psiS_array, MC **E_R_array, MC **E_RS_array, MC **E_SR_array){
  int good=0;

  int eA = atoi(eA_str);
  int eB = atoi(eB_str);

  int lA = atoi(lA_str);
  int lB = atoi(lB_str);


  // E
  MC *E;
  E = malloc(sizeof(MC));
  init_MC(E);
  copy_MC(E, PA->curve);

  printf("******************** E ********************\n");
  print_Curve(E);


  


  printf("---------------------Computing Peggy's Secret S\n");
  
  MP *S;
  S = malloc(sizeof(MP));
  init_MP(S);

  mpz_t *mA, *nA;
  mA = malloc(sizeof(mpz_t));
  nA = malloc(sizeof(mpz_t));
  mpz_init(*mA);
  mpz_init(*nA);

  rand_subgroup(mA,nA,lA_str,eA_str);

  shamir(&S->x, &S->y, &S->z, E->A, E->B, PA->x, PA->y, PA->z, QA->x, QA->y, QA->z, *mA, *nA);

  copy_MC(&S->curve, *E);


  printf("------------------Computing phi: E -> E/<S>,  phi(P_B), phi(Q_B)\n");

  MC *E_S;
  E_S = malloc(sizeof(MC));
  init_MC(E_S);
  copy_MC(E_S, PA->curve);

  MP *phiPB, *phiQB;
  phiPB = malloc(sizeof(MP));
  phiQB = malloc(sizeof(MP));
  init_MP(phiPB);
  init_MP(phiQB);

  copy_GF(&phiPB->x, PB->x);
  copy_GF(&phiPB->y, PB->y);
  copy_GF(&phiPB->z, PB->z);
  copy_GF(&phiQB->x, QB->x);
  copy_GF(&phiQB->y, QB->y);
  copy_GF(&phiQB->z, QB->z);

  push_through_iso(&E_S->A, &E_S->B, &E_S->A24, S->x, S->z, lA, strA, lenA-1, &phiPB->x, &phiPB->y, &phiPB->z, &phiQB->x, &phiQB->y, &phiQB->z, eA);

  /* not necessary
  copy_MC(&phiPB->curve, *E_S);
  copy_MC(&phiQB->curve, *E_S);
  */

  print_Curve(E_S);


  /////////////////////////////////////////////////////////////////////////////////
  // run the ZKP rounds

  for(int r=0; r < rounds; r++) {
    printf("round %d\n",r);

    printf("----------Choosing random R and computing phi(R)\n");
    
    mpz_t *mB, *nB;
    mB = malloc(sizeof(mpz_t));
    nB = malloc(sizeof(mpz_t));
    mpz_init(*mB);
    mpz_init(*nB);

    rand_subgroup(mB,nB,lB_str,eB_str);

    R_array[r] = malloc(sizeof(MP));
    init_MP(R_array[r]);

    shamir(&R_array[r]->x, &R_array[r]->y, &R_array[r]->z, E->A, E->B, PB->x, PB->y, PB->z, QB->x, QB->y, QB->z, *mB, *nB);

    copy_MC(&R_array[r]->curve, *E);

    // Compute phi(R) by mB*phi(PB)+nB*phi(QB) instead of pushing through isogeny
    phiR_array[r] = malloc(sizeof(MP));
    init_MP(phiR_array[r]);

    shamir(&phiR_array[r]->x, &phiR_array[r]->y, &phiR_array[r]->z, E_S->A, E_S->B, phiPB->x, phiPB->y, phiPB->z, phiQB->x, phiQB->y, phiQB->z, *mB, *nB);
    
    copy_MC(&phiR_array[r]->curve, *E_S);


    /////////////////////////////////////////////////////////
    printf("----------Computing E/<R> and psi(S)\n");

    E_R_array[r] = malloc(sizeof(MC));
    init_MC(E_R_array[r]);
    copy_MC(E_R_array[r], *E);

    psiS_array[r] = malloc(sizeof(MP));
    init_MP(psiS_array[r]);
    copy_MP(psiS_array[r], *S);
    
    push_through_iso(&E_R_array[r]->A, &E_R_array[r]->B, &E_R_array[r]->A24, R_array[r]->x, R_array[r]->z, lB, strB, lenB-1, &psiS_array[r]->x, &psiS_array[r]->y, &psiS_array[r]->z, NULL, NULL, NULL, eB);

    copy_MC(&psiS_array[r]->curve, *E_R_array[r]);


    //////////////////////////////////////////////////////
    printf("------------Computing E/<R,S>\n");

    E_RS_array[r] = malloc(sizeof(MC));
    init_MC(E_RS_array[r]);
    copy_MC(E_RS_array[r], *E_R_array[r]);

    push_through_iso(&E_RS_array[r]->A, &E_RS_array[r]->B, &E_RS_array[r]->A24, psiS_array[r]->x, psiS_array[r]->z, lA, strA, lenA-1, NULL, NULL, NULL, NULL, NULL, NULL, eA);


    
    free(mB); free(nB);



    //print_Curve(E_RS_array[r]);


    // check with E/<S,R>
    E_SR_array[r] = malloc(sizeof(MC));
    init_MC(E_SR_array[r]);
    copy_MC(E_SR_array[r], *E_S);

    push_through_iso(&E_SR_array[r]->A, &E_SR_array[r]->B, &E_SR_array[r]->A24, phiR_array[r]->x, phiR_array[r]->z, lB, strB, lenB-1, NULL, NULL, NULL, NULL, NULL, NULL, eB);



  }

  free(mA); free(nA);
  free(phiPB); free(phiQB);

  return 0;
}




//reads public parameters for use with ss_isogeny_exchange_dfc in from file. Note that this file must match the format and naming conventions
//of the one generated by ss_isogeny_gen_file().
void params_from_file(char * p, char *eA, char *eB, char *lA, char *lB, int *strA, int *lenA, int *strB, int *lenB, MP *PA, MP *QA, MP *PB, MP *QB, char * file ){

    FILE *fr; 
    fr = fopen (file, "rt");
    MC curve;
    init_MC(&curve);
    int lineMax = 10000;
    char line[lineMax];

    //taking in p,lA, lB, eA, eB
    fgets(line, lineMax, fr);
    sscanf (line, "%s", p);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", lA);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", lB);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", eA);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", eB);
    
    //taking in the curve
    GF A, B, A24;
    GF_params *parent;
    parent = malloc(sizeof(GF_params));
    setup_GF(parent, p);
   
    char *a;
    char *b;
    a = malloc(sizeof(char)*lineMax);
    b = malloc(sizeof(char)*lineMax);

    init_GF(&A, parent);
    init_GF(&B, parent);
    init_GF(&A24, parent);
 
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&A, a, b);
   
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&B, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
   
    set_GF(&A24, a, b);

    set_Curve( &curve, A, B, A24);
  
    //taking in PA
    GF x1,y1,z1;
    
    init_GF(&x1, parent);
    init_GF(&y1, parent);
    init_GF(&z1, parent);
   
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x1, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y1, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z1, a, b);
    
    set_MP( PA, x1, y1, z1, curve );
    
    //taking in QA
    GF x2,y2,z2;
    
    init_GF(&x2, parent);
    init_GF(&y2, parent);
    init_GF(&z2, parent);
    
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x2, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y2, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z2, a, b);
    
    set_MP( QA, x2, y2, z2, curve );
    
    //taking in PB
    GF x3,y3,z3;
    
    init_GF(&x3, parent);
    init_GF(&y3, parent);
    init_GF(&z3, parent);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x3, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y3, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z3, a, b);
    
    set_MP( PB, x3, y3, z3, curve );
    
    //taking in QB
    GF x4,y4,z4;
    
    init_GF(&x4, parent);
    init_GF(&y4, parent);
    init_GF(&z4, parent);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&x4, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&y4, a, b);
    
    fgets(line, lineMax, fr);
    sscanf (line, "%s", a);
    fgets(line, lineMax, fr);
    sscanf (line, "%s", b);
    
    set_GF(&z4, a, b);
    
    set_MP( QB, x4, y4, z4, curve );
    
    //taking in strA
    fgets(line, lineMax, fr);
    sscanf (line, "%d", lenA);
    
    int k=0;
    for(; k<*lenA; k++){
        fgets(line, lineMax, fr);
        sscanf (line, "%d", &strA[k]);
    }
    
    //taking in strB
    fgets(line, lineMax, fr);
    sscanf (line, "%d", lenB);
    
    int r=0;
    for(; r<*lenB; r++){
        fgets(line, lineMax, fr);
        sscanf (line, "%d", &strB[r]);
    }
    
    free(a);
    free(b);
}


char *concat(const char *s1, ...)
{
  va_list args;
  const char *s;
  char *p, *result;
  unsigned long l, m, n;
 
  m = n = strlen(s1);
  va_start(args, s1);
  while ((s = va_arg(args, char *))) {
    l = strlen(s);
    if ((m += l) < l) break;
  }
  va_end(args);
  if (s || m >= INT_MAX) return NULL;
 
  result = (char *)malloc(m + 1);
  if (!result) return NULL;
 
  memcpy(p = result, s1, n);
  p += n;
  va_start(args, s1);
  while ((s = va_arg(args, char *))) {
    l = strlen(s);
    if ((n += l) < l || n > m) break;
    memcpy(p, s, l);
    p += l;
  }
  va_end(args);
  if (s || m != n || p != result + n) {
    free(result);
    return NULL;
  }
 
  *p = 0;
  return result;
}


char* join_strings(char* strings[], char* seperator, int count) {
    char* str = NULL;             /* Pointer to the joined strings  */
    size_t total_length = 0;      /* Total length of joined strings */
    int i = 0;                    /* Loop counter                   */

    /* Find total length of joined strings */
    for (i = 0; i < count; i++) total_length += strlen(strings[i]);
    total_length++;     /* For joined string terminator */
    total_length += strlen(seperator) * (count - 1); // for seperators

    str = (char*) malloc(total_length);  /* Allocate memory for joined strings */
    str[0] = '\0';                      /* Empty string we can append to      */

    /* Append all the strings */
    for (i = 0; i < count; i++) {
        strcat(str, strings[i]);
        if (i < (count - 1)) strcat(str, seperator);
    }

    return str;
}


char *MC_str(MC curve, int base) {
  return concat(mpz_get_str(NULL, base, curve.A.a), mpz_get_str(NULL, base, curve.A.b),
                mpz_get_str(NULL, base, curve.B.a), mpz_get_str(NULL, base, curve.B.b),
                mpz_get_str(NULL, base, curve.A24.a), mpz_get_str(NULL, base, curve.A24.b), NULL);
}

char *MP_str(MP *point, int base) {
  print_MP(point,"point");
  char *st = concat(mpz_get_str(NULL, base, point->x.a), mpz_get_str(NULL, base, point->x.b),
                    mpz_get_str(NULL, base, point->y.a), mpz_get_str(NULL, base, point->y.b),
                    mpz_get_str(NULL, base, point->z.a), mpz_get_str(NULL, base, point->z.b),
                    mpz_get_str(NULL, base, point->curve.A.a), mpz_get_str(NULL, base, point->curve.A.b),
                    mpz_get_str(NULL, base, point->curve.B.a), mpz_get_str(NULL, base, point->curve.B.b),
                    mpz_get_str(NULL, base, point->curve.A24.a), mpz_get_str(NULL, base, point->curve.A24.b), NULL);
  return st;
}


void string_data(char** data, int rounds, MC **com1, MC **com2, int *chal, MP **resp1, MP **resp2, char **hresp, int hrlen) {
  char **rdata;
  rdata = malloc(rounds * sizeof(char*));

  int base = 10;
  for(int r=0; r<rounds; r++) {
      char *E1 = concat(mpz_get_str(NULL, base, com1[r]->A.a), mpz_get_str(NULL, base, com1[r]->A.b),
                        mpz_get_str(NULL, base, com1[r]->B.a), mpz_get_str(NULL, base, com1[r]->B.b),
                        mpz_get_str(NULL, base, com1[r]->A24.a), mpz_get_str(NULL, base, com1[r]->A24.b), NULL);
      //printf("E1: %s\n", E1);  

      char *E2 = concat(mpz_get_str(NULL, base, com2[r]->A.a), mpz_get_str(NULL, base, com2[r]->A.b),
                        mpz_get_str(NULL, base, com2[r]->B.a), mpz_get_str(NULL, base, com2[r]->B.b),
                        mpz_get_str(NULL, base, com2[r]->A24.a), mpz_get_str(NULL, base, com2[r]->A24.b), NULL);
      //printf("E2: %s\n", E2);  

      char *ch;
      if (chal[r] == 0) ch = "01";
      else ch = "10";
      //printf("ch: %s\n", ch);

      char *ch1R1 = concat(mpz_get_str(NULL, base, resp1[2*r]->x.a), mpz_get_str(NULL, base, resp1[2*r]->x.b),
                        mpz_get_str(NULL, base, resp1[2*r]->y.a), mpz_get_str(NULL, base, resp1[2*r]->y.b),
                        mpz_get_str(NULL, base, resp1[2*r]->z.a), mpz_get_str(NULL, base, resp1[2*r]->z.b),
                        mpz_get_str(NULL, base, resp1[2*r]->curve.A.a), mpz_get_str(NULL, base, resp1[2*r]->curve.A.b),
                        mpz_get_str(NULL, base, resp1[2*r]->curve.B.a), mpz_get_str(NULL, base, resp1[2*r]->curve.B.b),
                        mpz_get_str(NULL, base, resp1[2*r]->curve.A24.a), mpz_get_str(NULL, base, resp1[2*r]->curve.A24.b), NULL);
      //printf("ch1R1: %s\n", ch1R1);  

      char *ch1R2 = concat(mpz_get_str(NULL, base, resp2[2*r]->x.a), mpz_get_str(NULL, base, resp2[2*r]->x.b),
                        mpz_get_str(NULL, base, resp2[2*r]->y.a), mpz_get_str(NULL, base, resp2[2*r]->y.b),
                        mpz_get_str(NULL, base, resp2[2*r]->z.a), mpz_get_str(NULL, base, resp2[2*r]->z.b),
                        mpz_get_str(NULL, base, resp2[2*r]->curve.A.a), mpz_get_str(NULL, base, resp2[2*r]->curve.A.b),
                        mpz_get_str(NULL, base, resp2[2*r]->curve.B.a), mpz_get_str(NULL, base, resp2[2*r]->curve.B.b),
                        mpz_get_str(NULL, base, resp2[2*r]->curve.A24.a), mpz_get_str(NULL, base, resp2[2*r]->curve.A24.b), NULL);
      //printf("ch1R2: %s\n", ch1R2);  

      char *ch2R1 = concat(mpz_get_str(NULL, base, resp1[2*r+1]->x.a), mpz_get_str(NULL, base, resp1[2*r+1]->x.b),
                        mpz_get_str(NULL, base, resp1[2*r+1]->y.a), mpz_get_str(NULL, base, resp1[2*r+1]->y.b),
                        mpz_get_str(NULL, base, resp1[2*r+1]->z.a), mpz_get_str(NULL, base, resp1[2*r+1]->z.b),
                        mpz_get_str(NULL, base, resp1[2*r+1]->curve.A.a), mpz_get_str(NULL, base, resp1[2*r+1]->curve.A.b),
                        mpz_get_str(NULL, base, resp1[2*r+1]->curve.B.a), mpz_get_str(NULL, base, resp1[2*r+1]->curve.B.b),
                        mpz_get_str(NULL, base, resp1[2*r+1]->curve.A24.a), mpz_get_str(NULL, base, resp1[2*r+1]->curve.A24.b), NULL);
      //printf("ch2R1: %s\n", ch2R1);  

      char *ch2R2 = concat(mpz_get_str(NULL, base, resp2[2*r+1]->x.a), mpz_get_str(NULL, base, resp2[2*r+1]->x.b),
                        mpz_get_str(NULL, base, resp2[2*r+1]->y.a), mpz_get_str(NULL, base, resp2[2*r+1]->y.b),
                        mpz_get_str(NULL, base, resp2[2*r+1]->z.a), mpz_get_str(NULL, base, resp2[2*r+1]->z.b),
                        mpz_get_str(NULL, base, resp2[2*r+1]->curve.A.a), mpz_get_str(NULL, base, resp2[2*r+1]->curve.A.b),
                        mpz_get_str(NULL, base, resp2[2*r+1]->curve.B.a), mpz_get_str(NULL, base, resp2[2*r+1]->curve.B.b),
                        mpz_get_str(NULL, base, resp2[2*r+1]->curve.A24.a), mpz_get_str(NULL, base, resp2[2*r+1]->curve.A24.b), NULL);
      //printf("ch2R2: %s\n", ch2R2);  

      rdata[r] = concat(E1, E2, ch, ch1R1, ch1R2, ch2R1, ch2R2, hresp[2*r], hresp[2*r+1], NULL);
      //printf("\nround %d data: %s\n", r, rdata[r]);

    }
    *data = join_strings(rdata, "", rounds);
    //printf("\n\nall data:\n%s\n", *data);
}

//1st argument specifies file with parameters, second # of times to run the key exchange. If second argument is null, only 1 iteration is performed
int main(int argc, char *argv[]) {
    int iterations;
    srand(time);
    
    if( argc < 2 ){
        printf("ERROR: Must specify file name in 1st command line agrument.\n");
        
    }else{
        printf("File where Public parameters are: '%s'\n", argv[1]);
    }
    
    if( !argv[2] )
        iterations = 1;
    else{
        iterations = atoi(argv[2]);
        printf("Number of Iterations: %s\n", argv[2]);
    }    
    
    int MAX_LENGTH = 10000;
    int *strA, *strA_t, *strB, *strB_t; 
    int lenA=0;
    int lenB=0;
    MP *PA;
    MP *QA;
    MP *PB;
    MP *QB;
    char *p, *eA, *eB, *lA, *lB;
    
    p = malloc(sizeof(char)*MAX_LENGTH);
    eA = malloc(sizeof(char)*MAX_LENGTH);
    eB = malloc(sizeof(char)*MAX_LENGTH);
    lA = malloc(sizeof(char)*MAX_LENGTH);
    lB = malloc(sizeof(char)*MAX_LENGTH);
    PA = malloc(sizeof(MP));
    QA = malloc(sizeof(MP));
    PB = malloc(sizeof(MP));
    QB = malloc(sizeof(MP));
    strA_t = malloc(MAX_LENGTH*sizeof(int));
    strB_t = malloc(MAX_LENGTH*sizeof(int));
    
    params_from_file( p, eA, eB, lA, lB, strA_t, &lenA, strB_t, &lenB, PA, QA, PB, QB, argv[1] );

    strA = malloc(lenA*sizeof(int));
    strB = malloc(lenB*sizeof(int));
    
    //Putting strategies into arrays of the exact required length
    int r=0;
    for(r; r<lenA; r++)
        strA[r] = strA_t[r];
        
    int k=0;
    for(k; k<lenB; k++)
        strB[k] = strB_t[k]; 

    

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    int rounds=32; //also equal to the bit length of hash output (must be a multiple of 8)

    MP *R_array[rounds];
    MP *phiR_array[rounds];
    MP *psiS_array[rounds];

    MC *E_R_array[rounds];
    MC *E_RS_array[rounds];
    MC *E_SR_array[rounds];
    
    run_ZKP(time, eA, eB, lA, lB, strA, lenA, strB, lenB, PA, QA, PB, QB, rounds, R_array, phiR_array, psiS_array, E_R_array, E_RS_array, E_SR_array);

    for (int r=0; r<rounds; r++) {
      printf("\n\nround %d \n", r+1);
      print_MP(R_array[r], "R");
      print_MP(phiR_array[r], "phi(R)");
      print_MP(psiS_array[r], "psi(S)");
      printf("********* E/<R> *********\n");
      print_Curve(E_R_array[r]);
      printf("********* E/<R,S> *********\n");
      print_Curve(E_RS_array[r]);
      printf("********* E/<S,R> *********\n");
      print_Curve(E_SR_array[r]);
      
    }



/*
    //compute commitment/challenge/responses
    

    MC **com1, **com2;
    com1 = malloc(rounds * sizeof(MC*));
    com2 = malloc(rounds * sizeof(MC*));
    printf("com1: %p, com2: %p\n", com1, com2);
    
    int *chal;
    chal = malloc(rounds * sizeof(int)); //first challenge bits
    
    MP **resp1, **resp2;
    resp1 = malloc(2 * rounds * sizeof(MP*));
    resp2 = malloc(2 * rounds * sizeof(MP*));

    for(int r=0; r<rounds; r++) {
      printf("\nbegin round %d\n",r);
      ZKP_identity(time, eA, eB, lA, lB, strA, lenA, strB, lenB, PA, QA, PB, QB, r, com1, com2, chal, resp1, resp2);
      printf("\nfinished round %d\n\n",r);
    }

    for(int r=0; r<rounds; r++) {
      printf("Round %d\n\n", r);
      printf("E1 = E/<R>------------------------\n");
      print_Curve(com1[r]);
      printf("E2 = E/<S,R>------------------------\n");
      print_Curve(com2[r]);
      printf("\nChallenge bit order: %d, %d\n\n", chal[r], 1-chal[r]);
      printf("Challenge 1  ----------------\n");
      print_MP(resp1[2*r], "resp1");
      print_MP(resp2[2*r], "resp2");
      printf("Challenge 2  ----------------\n");
      print_MP(resp1[2*r +1], "resp1");
      print_MP(resp2[2*r +1], "resp2");
      printf("\n\n");
    }

    
    // put all the data into a string to compute hash
    char *data = NULL;
    int base = 10;

    uint8_t **hresp;
    hresp = malloc(2 * rounds * sizeof(char*));
    int hrlen = 64; // in bytes

    for(int r=0; r<rounds; r++) {
      uint8_t hr1[hrlen];
      uint8_t hr2[hrlen];

      char *ch1resp = concat(MP_str(resp1[2*r], base), MP_str(resp2[2*r], base), NULL);
      char *ch2resp = concat(MP_str(resp1[2*r+1], base), MP_str(resp2[2*r+1], base), NULL);
      

      keccak((uint8_t*) ch1resp, strlen(ch1resp), hr1, hrlen);
      keccak((uint8_t*) ch2resp, strlen(ch2resp), hr2, hrlen);
      
      /*
      printf("hash 1: ");
      for(int i=0; i<hrlen; i++) {
        printf("%02X", hr1[i]);
      }
      printf("\nhash 2: ");
      for(int i=0; i<hrlen; i++) {
        printf("%02X", hr2[i]);
      }
      printf("\n");
      */ /*

      //converting uint8_t array into char array
      char h1[hrlen+1];
      char h2[hrlen+1];
      h1[hrlen] = '\0';
      h2[hrlen] = '\0';
      memcpy(h1, hr1, hrlen);
      memcpy(h2, hr2, hrlen);

/*
      for(int i=0; i<hrlen; i++) {
        printf("%d ", (uint8_t)h1[i]);
      }
      printf("\n");
      for(int i=0; i<hrlen; i++) {
        printf("%d ", (uint8_t)h2[i]);
      }
      printf("\nstring1: %s\nstring2: %s\n", h1, h2);
*/ /*

      hresp[2*r] = h1;
      hresp[2*r+1] = h2;

    }

    string_data(&data, rounds, com1, com2, chal, resp1, resp2, hresp, hrlen);
    //printf("data: %s\n",data);
    printf("\ndata length %lu\n", strlen(data));


    int hashlength = rounds/8;
    uint8_t hash[hashlength];
    printf("hashlength: %d bytes (%d bits)\n", hashlength, rounds);

    keccak((uint8_t*) data, strlen(data), hash, hashlength);
    
    printf("hash (Hex): ");
    for(int i=0; i<hashlength; i++) {
      printf("%02X", hash[i]);
    }
    printf("\n");

    


    /*
    for(i; i<iterations; i++){
        
      //  pthread_create(&threads[0], NULL, p2, NULL);
      //  pthread_create(&threads[1], NULL, p3, NULL);
        
        double *time;
        time = malloc(sizeof(double));
        good = ss_isogeny_exchange_dfc(time, eA, eB, lA, lB, strA, lenA, strB, lenB, PA, QA, PB, QB);
        totalTime += *time;
        
        if (!good)
            errors +=1;
    }
    */
    /*
    avgTime = totalTime/(iterations*1.0);
    printf ("\nAverage Time for %d iteration(s) (sec) : %f \n", iterations, avgTime);
    printf ("\n Number of key exchanges NOT completed successfully : %d \n", errors);
    printf("\nWITH_FIX\n");
    printf("\nParameters: %s\n",  argv[1]);
    
 /*   
    arith *info = malloc(sizeof(arith));
    info->a = &((*PA).x); 
    info->x = &((*PA).y);
    info->y = &((*PA).z);
    
    add_GF(&((*QA).x), (*PA).y, (*PA).z);
    print_GF((*QA).x, "res");
    
    pthread_t threads[10];  
    int rc = pthread_create(&threads[0], NULL, add_t, (void *) info); 
    
    print_GF((*PA).x, "3tmp4");  
    printf("\n\n*************\n\n");   
    
    */
    
    //pthread_exit(NULL);
    
    free(p); free(eA); free(eB); free(lA); free(lB); free(PA); free(QA); free(PB); free(QB); free(strA_t); free(strB_t); free(strA); free(strB);

}  

  
 
    
      
