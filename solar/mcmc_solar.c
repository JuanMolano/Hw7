//
//  mcmc_solar.c
//
//
//  Created by Juan Pablo Molano on 11/11/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#define USAGE "./mcmc_solar.x n_burn n_steps"

//------------------------------------------------------------------------------------------------------------

float *load_matrix(char *filename, int *n, int *m);
void get_data( float *array , float *x, int n_row);
float *reserva(int n);
void print_data(float *array, float *array2, float *array3,  float *array4 , float *array5 , float *array6 , float *array7 , float *array8,  int n_puntos);
void get_patch( float *array,  float *x , int n_row);
void my_model( float *x_obs, float *array, float a, float b, float c , float d,  int n_row);
float likelihood(float *y_obs, float *y_model, int n_row);
int min_likelihood( float *array , int n_row);
void print_array(float *array, int n_puntos);

//------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv){
    float *matrix , *tim, *man, *count, *a_walk, *b_walk, *c_walk, *d_walk, *l_walk, *y_init, *y_prime, *best_y;
    srand48(time(NULL));
    float a_prime, b_prime, c_prime, d_prime, l_prime, l_init, gamma, alpha, beta , best_a , best_b, best_c, best_d;
    int n_row, n_cols , i, rlin = 0;
    
    
    a_walk = reserva(atof(argv[1])+atof(argv[2]));
    b_walk = reserva(atof(argv[1])+atof(argv[2]));
    c_walk = reserva(atof(argv[1])+atof(argv[2]));
    d_walk = reserva(atof(argv[1])+atof(argv[2]));
    l_walk = reserva(atof(argv[1])+atof(argv[2]));
    y_init = reserva(n_row);
    y_prime = reserva(n_row);

    
    a_walk[0] = drand48();
    b_walk[0] = drand48();
    c_walk[0] = drand48();
    d_walk[0] = drand48();
    
    

    matrix = load_matrix( "sol.dat" , &n_row, &n_cols);
    
    man = reserva( n_row);
    tim = reserva(n_row);
    man = reserva(n_row);
    best_y = reserva( n_row);
    get_data( matrix , tim , n_row);
    get_patch( matrix ,  man  ,  n_row);
    
    for (i=0;i<n_row;i++) {
        if (man[i] == -99) {
            man[i] = 0.0;
        }
    }
    
    my_model( tim , y_init , a_walk[0] , b_walk[0] , c_walk[0] , d_walk[0] , n_row );
    
    l_walk[0] = likelihood( man , y_init, n_row);
    
    

    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
  
    for (i=0; i < atof(argv[1]) + atof(argv[2]); i++) {
        
        a_prime = gsl_ran_gaussian(r, 0.1) + a_walk[i];
        b_prime = gsl_ran_gaussian(r, 0.1) + b_walk[i];
        c_prime = gsl_ran_gaussian(r, 0.1) + c_walk[i];
        d_prime = gsl_ran_gaussian(r, 0.1) + d_walk[i];
        
        my_model( tim, y_init, a_walk[i], b_walk[i], c_walk[i], d_walk[i], n_row);
        my_model(tim, y_prime ,a_prime, b_prime, c_prime, d_prime, n_row);
        
        l_prime = likelihood( man, y_prime , n_row);
        l_init = likelihood( man, y_init ,  n_row);
        
        gamma = l_init - l_prime;
        if (gamma>=0.0) {
            a_walk[i+1] = a_prime;
            b_walk[i+1] = b_prime;
            c_walk[i+1] = c_prime;
            d_walk[i+1] = d_prime;
            l_walk[i+1] = l_prime;
        }
        else{
            
            beta = drand48();
            alpha = exp(gamma);
            if (beta <= alpha) {
                
                a_walk[i+1]  = a_prime;
                b_walk[i+1]  = b_prime;
                c_walk[i+1]  = c_prime;
                d_walk[i+1]  = d_prime;
                l_walk[i+1]  = l_prime;
            }
            else{
                
                a_walk[i+1] = a_walk[i];
                b_walk[i+1] = b_walk[i];
                c_walk[i+1] = c_walk[i];
                d_walk[i+1] = d_walk[i];
                l_walk[i+1] = l_init;
                
            }
            
        }
    }
    
    i = min_likelihood( l_walk , n_row);
    best_a = a_walk[i];
    best_b = b_walk[i];
    best_c = c_walk[i];
    best_d = d_walk[i];
    
    my_model( tim , best_y , best_a, best_b , best_c , best_d , n_row);
    
    print_data( tim , man,  l_walk , a_walk , b_walk , c_walk , d_walk , best_y, n_row);
    
    gsl_rng_free (r);

    
    return 0;
}

//------------------------------------------------------------------------------------------------------------

float *load_matrix(char *filename, int *n, int *m){
    float *matrix;
    FILE *in;
    int n_row = 4632, n_cols = 5;
    int i;
    int j;
    
    if(!(in=fopen(filename, "r"))){
        printf("Problem opening file %s\n", filename);
        exit(1);
    }

    
    matrix = malloc(n_row * n_cols * sizeof(float));
    
    for(i=0;i<n_row;i++){
        for(j=0;j<n_cols;j++){
            fscanf(in, "%f", &matrix[i*n_cols + j]);
        }
    }    
    *n = n_row;
    *m = n_cols;
    return matrix;
}

//------------------------------------------------------------------------------------------------------------

void get_data( float *array , float *x , int n_row){
    
    int i,j;
    
    for(i=0;i<n_row;i++){
        
            x[i] = array[i*5] + array[i*5 + 1]/12 + array[i*5 + 2]/365;
       
    }

}

//------------------------------------------------------------------------------------------------------------

void get_patch( float *array,  float *x , int n_row){
    
    int i,j;
    
    for(i=0;i<n_row;i++){
        
        x[i] = array[i*5 + 3];
    }
    

}

//------------------------------------------------------------------------------------------------------------

float *reserva(int n){
    float *array;
    int i;
    if(!(array = malloc(n * sizeof(float)))){
        printf("Problema en reserva\n");
        exit(1);
    }
    for(i=0;i<n ;i++){
        array[i] = 0.0;
    }
    return array;
}

//------------------------------------------------------------------------------------------------------------

void print_data(float *array, float *array2 , float *array3 , float *array4 , float *array5 , float *array6, float *array7, float *array8 , int n_puntos){
    int i;
    for(i=0;i<n_puntos;i++){
        printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", array[i] , array2[i] , array3[i] , array4[i] , array5[i] , array6[i], array7[i], array8[i]);
    }
}

//------------------------------------------------------------------------------------------------------------

float likelihood(float *y_obs, float *y_model, int n_row){
    
    int i;
    float c=0;
    
    for (i=0; i<n_row; i++) {

        c = c + ( y_obs[i] - y_model[i]);
    }
    c = (1.0/2.0)*pow(c,2);
    

    return c;
    
}

//------------------------------------------------------------------------------------------------------------

void my_model( float *x_obs, float *array , float a, float b, float c , float d,  int n_row){
    

    int i;
    for (i=0; i<n_row; i++) {
        array[i] = a*cos((2*3.141592653589793/d)*x_obs[i] + b) + c;
    }
}

//------------------------------------------------------------------------------------------------------------

int min_likelihood( float *array , int n_row){
    
    int i;
    int ind = 0;
    
    for(i = 1; i < n_row; i ++){
        if ( array[i] < ind) {
            ind = i;
        }
        
    }

    return ind;
}

void print_array(float *array, int n_puntos){
    int i;
    for(i=0;i<n_puntos;i++){
        printf("%f \n", array[i]);
    }
}