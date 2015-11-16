//
//  mcmc_lotkavolterra.c
//  
//
//  Created by Juan Pablo Molano on 15/11/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#define USAGE "./mcmc_lotkavolterra.x n_burn n_steps"

//------------------------------------------------------------------------------------------------------------

void load_matrix(char *filename, float *time, float *pre, float *dep , int n);
void print_data(float *array, float *array2 , float *array3 , float *array4 , float *array5 , float *array6, float *array7, float *array8 , float *array9 , float *array10,  int n_puntos);
float *reserva(int n);
float likelihood(float *y_obs, float *y_model, int n_row);
void ec_dep(float *array, float x, float y, float gamma, float delta);
void ec_pre(float *array, float x, float y, float alpha, float beta);
int min_likelihood( float *array , int n_row);


//------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv){
    
    float *tim, *pre, *dep;
    float  m_x , m_y , k1_x , k1_y , k2_x , k2_y ,  k3_x , k3_y ,  k4_x , k4_y;
    float x1 , x2 , x3, x4, y1, y2, y3, y4 , h;
    float alpha_prime, beta_prime, gamma_prime, delta_prime, l_prime, l_init, gam, alph, bet , best_alpha , best_beta, best_gamma, best_delta;
    float *fit_pre, *fit_dep , *dep_init, *dep_prime , *pre_init , *pre_prime;
    int i, j , n_row = 96;
    float *alpha_a, *beta_a , *gamma_a , *delta_a , *walk_a, *best_dep, *best_pre;
    srand48(time(NULL));

    alpha_a = reserva( atof(argv[1]) + atof(argv[2]));
    beta_a = reserva( atof(argv[1]) + atof(argv[2]));
    gamma_a= reserva( atof(argv[1]) + atof(argv[2]));
    delta_a = reserva( atof(argv[1]) + atof(argv[2]));
    walk_a = reserva( atof(argv[1]) + atof(argv[2]));
    
    alpha_a[0] = drand48();
    beta_a[0] = drand48();
    gamma_a[0] = drand48();
    delta_a[0] = drand48();

    tim = reserva(96);
    pre = reserva(96);
    dep = reserva(96);
    fit_dep = reserva(96);
    fit_pre = reserva(96);
    dep_init = reserva(96);
    dep_prime = reserva(96);
    pre_init = reserva(96);
    pre_prime = reserva(96);
    best_dep = reserva(96);
    best_pre = reserva(96);
    
    load_matrix( "new_data.txt" , tim , pre , dep, 96);

    h = tim[1] - tim[0];
    
    fit_pre[0] = pre[0];
    fit_dep[0] = dep[0];
  
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    for(i=1;i<96;i++){
        
        ec_pre( &k1_x , fit_pre[i-1] , fit_dep[i-1] , alpha_a[0] , beta_a[0]);
        ec_dep( &k1_y , fit_pre[i-1] , fit_dep[i-1] , gamma_a[0] , delta_a[0]);
    
        x1 = fit_pre[i-1] + (h/2.0) * k1_x;
        y1 = fit_dep[i-1] + (h/2.0) * k1_y;
    

        ec_pre( &k2_x , fit_pre[i-1] , fit_dep[i-1] , alpha_a[0] , beta_a[0]);
        ec_dep( &k2_y , fit_pre[i-1] , fit_dep[i-1] , gamma_a[0] , delta_a[0]);
    

        x2 = fit_pre[i-1] + (h/2.0) * k2_x;
        y2 = fit_dep[i-1] + (h/2.0) * k2_y;
    

        ec_pre( &k3_x , fit_pre[i-1] , fit_dep[i-1] , alpha_a[0] , beta_a[0]);
        ec_dep( &k3_y , fit_pre[i-1] , fit_dep[i-1] , gamma_a[0] , delta_a[0]);
    

        x3 = fit_pre[i-1] + h * k3_x;
        y3 = fit_dep[i-1] + h * k3_y;
    
        ec_pre( &k4_x , fit_pre[i-1] , fit_dep[i-1] , alpha_a[0] , beta_a[0]);
        ec_dep( &k4_y , fit_pre[i-1] , fit_dep[i-1] , gamma_a[0] , delta_a[0]);
    
    
        m_x = (1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
        m_y = (1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
    
        fit_pre[i] = fit_pre[i-1] + h * m_x;
        fit_dep[i] = fit_dep[i-1] + h * m_y;
    }

    for (i=0; i < atof(argv[1]) + atof(argv[2]); i++) {
        
        alpha_prime = gsl_ran_gaussian(r, 0.1) + alpha_a[i];
        beta_prime = gsl_ran_gaussian(r, 0.1) + beta_a[i];
        gamma_prime = gsl_ran_gaussian(r, 0.1) + gamma_a[i];
        delta_prime = gsl_ran_gaussian(r, 0.1) + delta_a[i];
        
        
        for(j=1;j<96;j++){
            
            dep_init[0] = dep[0];
            pre_init[0] = pre[0];
            
            
            ec_pre( &k1_x , pre_init[j-1] , dep_init[j-1] , alpha_a[i] , beta_a[i]);
            ec_dep( &k1_y , pre_init[j-1] , dep_init[j-1] , gamma_a[i] , delta_a[i]);
            
            x1 = pre_init[j-1] + (h/2.0) * k1_x;
            y1 = dep_init[j-1] + (h/2.0) * k1_y;
            
            
            ec_pre( &k2_x , pre_init[j-1] , dep_init[j-1] , alpha_a[i] , beta_a[i]);
            ec_dep( &k2_y , pre_init[j-1] , dep_init[j-1] , gamma_a[i] , delta_a[i]);
            
            
            x2 = pre_init[j-1] + (h/2.0) * k2_x;
            y2 = dep_init[j-1] + (h/2.0) * k2_y;
            
            
            ec_pre( &k3_x , pre_init[j-1] , dep_init[j-1] , alpha_a[i] , beta_a[i]);
            ec_dep( &k3_y , pre_init[j-1] , dep_init[j-1] , gamma_a[i] , delta_a[i]);
            
            
            x3 = pre_init[j-1] + h * k3_x;
            y3 = dep_init[j-1] + h * k3_y;
            
            ec_pre( &k4_x , pre_init[j-1] , dep_init[j-1] , alpha_a[i] , beta_a[i]);
            ec_dep( &k4_y , pre_init[j-1] , dep_init[j-1] , gamma_a[i] , delta_a[i]);
            
            
            m_x = (1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
            m_y = (1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
            
            pre_init[j] = pre_init[i-1] + h * m_x;
            dep_init[j] = dep_init[i-1] + h * m_y;
        }
        
        for(j=1;j<96;j++){
            
            dep_prime[0] = dep[0];
            pre_prime[0] = pre[0];
            
            
            ec_pre( &k1_x , pre_prime[j-1] , dep_prime[j-1] , alpha_prime , beta_prime);
            ec_dep( &k1_y , pre_prime[j-1] , dep_prime[j-1] , gamma_prime , delta_prime);
            
            x1 = pre_prime[j-1] + (h/2.0) * k1_x;
            y1 = dep_prime[j-1] + (h/2.0) * k1_y;
            
            
            ec_pre( &k2_x , pre_prime[j-1] , dep_prime[j-1] ,  alpha_prime , beta_prime);
            ec_dep( &k2_y , pre_prime[j-1] , dep_prime[j-1] , gamma_prime , delta_prime);
            
            
            x2 = pre_prime[j-1] + (h/2.0) * k2_x;
            y2 = dep_prime[j-1] + (h/2.0) * k2_y;
            
            
            ec_pre( &k3_x , pre_prime[j-1] , dep_prime[j-1] ,  alpha_prime , beta_prime);
            ec_dep( &k3_y , pre_prime[j-1] , dep_prime[j-1] , gamma_prime , delta_prime);
            
            
            x3 = pre_prime[j-1] + h * k3_x;
            y3 = dep_prime[j-1] + h * k3_y;
            
            ec_pre( &k4_x , pre_prime[j-1] , dep_prime[j-1] ,  alpha_prime , beta_prime);
            ec_dep( &k4_y , pre_prime[j-1] , dep_prime[j-1] , gamma_prime , delta_prime);
            
            
            m_x = (1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
            m_y = (1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
            
            pre_prime[j] = pre_prime[j-1] + h * m_x;
            dep_prime[j] = dep_prime[j-1] + h * m_y;
        }

        
        l_prime = likelihood( dep, dep_prime , n_row) + likelihood( pre , pre_prime , n_row);
        l_init = likelihood( dep , dep_init ,  n_row) + likelihood( pre , pre_init , n_row);
        
        gam = l_init - l_prime;
        
        if (gam>=0.0) {
            alpha_a[i+1] = alpha_prime;
            beta_a[i+1] = beta_prime;
            gamma_a[i+1] = gamma_prime;
            delta_a[i+1] = delta_prime;
            walk_a[i+1] = l_prime;
        }
        else{
            
            bet = drand48();
            alph = exp(gam);
            
            if (bet <= alph) {
                
                alpha_a[i+1] = alpha_prime;
                beta_a[i+1] = beta_prime;
                gamma_a[i+1] = gamma_prime;
                delta_a[i+1] = delta_prime;
                walk_a[i+1] = l_prime;

            }
            else{
                
                alpha_a[i+1] = alpha_a[i];
                beta_a[i+1] = beta_a[i];
                gamma_a[i+1] = gamma_a[i];
                delta_a[i+1] = delta_a[i];
                walk_a[i+1] = l_prime;
                
            }
            
        }
    }
    
    j = min_likelihood( walk_a , n_row);
    best_alpha = alpha_a[j];
    best_beta = beta_a[j];
    best_gamma = gamma_a[j];
    best_delta = delta_a[j];
    
    for(i=1;i<96;i++){
        
        best_dep[0] = dep[0];
        best_pre[0] = pre[0];
        
        ec_pre( &k1_x , best_pre[i-1] , best_dep[i-1] , best_alpha , best_beta);
        ec_dep( &k1_y ,  best_pre[i-1] , best_dep[i-1], best_gamma , best_delta);
        
        x1 = fit_pre[i-1] + (h/2.0) * k1_x;
        y1 = fit_dep[i-1] + (h/2.0) * k1_y;
        
        
        ec_pre( &k2_x ,  best_pre[i-1] , best_dep[i-1] ,  best_alpha , best_beta);
        ec_dep( &k2_y ,  best_pre[i-1] , best_dep[i-1] , best_gamma , best_delta);
        
        
        x2 = fit_pre[i-1] + (h/2.0) * k2_x;
        y2 = fit_dep[i-1] + (h/2.0) * k2_y;
        
        
        ec_pre( &k3_x ,  best_pre[i-1] , best_dep[i-1],  best_alpha , best_beta);
        ec_dep( &k3_y ,  best_pre[i-1] , best_dep[i-1], best_gamma , best_delta);
        
        
        x3 = fit_pre[i-1] + h * k3_x;
        y3 = fit_dep[i-1] + h * k3_y;
        
        ec_pre( &k4_x ,  best_pre[i-1] , best_dep[i-1],  best_alpha , best_beta);
        ec_dep( &k4_y ,  best_pre[i-1] , best_dep[i-1] , best_gamma , best_delta);
        
        
        m_x = (1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
        m_y = (1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
        
        best_pre[i] = best_pre[i-1] + h * m_x;
        best_dep[i-1] = best_dep[i-1] + h * m_y;
    }
    
    
    print_data( tim , pre, dep, best_pre , best_pre, walk_a , alpha_a , beta_a , gamma_a , delta_a , n_row);
    
    gsl_rng_free (r);

    
    return 0;

}





//------------------------------------------------------------------------------------------------------------

void load_matrix(char *filename, float *time, float *pre, float *dep , int n){
    
    FILE *in;
    int i;
    in=fopen(filename,"r");
    if(!in){
        printf("problems opening the file %s\n", filename);
        exit(1);
    }
    for(i=0;i<(n);i++){
        
        fscanf(in, "%f %f %f\n", &time[i], &pre[i], &dep[i]);
        
    }
    fclose(in);
}

//------------------------------------------------------------------------------------------------------------

void print_data(float *array, float *array2 , float *array3 , float *array4 , float *array5 , float *array6, float *array7, float *array8 , float *array9 ,  float *array10 ,   int n_puntos){
    int i;
    for(i=0;i<n_puntos;i++){
        printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", array[i] , array2[i] , array3[i] , array4[i] , array5[i] , array6[i], array7[i] , array8[i] , array9[i] , array10[i]);
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


//------------------------------------------------------------------------------------------------------------

void ec_pre(float *array, float x, float y, float alpha, float beta){
    
    *array = x*(alpha - beta*y);
    
}

//------------------------------------------------------------------------------------------------------------

void ec_dep(float *array, float x, float y, float gamma, float delta){
    
    *array = -y*(gamma - delta*x);
    
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

