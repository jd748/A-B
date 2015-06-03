#include <random>
#include "Eigen/Eigen"
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "DP.hh"

int main()
{
    // Patients (n), attributes (p), iterations in solving DP (N1), iterations in estimation of ratio (N2);
    int p = 2;
    int n = 6;
    int N1 = 100;
    int N2 = 100;

    std::mt19937 generator1 (61245);
    std::mt19937 generator2 (16746);
    std::mt19937 generator3 (27351);
    std::mt19937 generator4 (66459);
    std::mt19937 generator5 (12612);
    std::normal_distribution <double> nd(0.0, 1.0);
    std::chi_squared_distribution <double> xs(p - 2);

    //s is Cholesky factorization of Covar matrix
    //mu is vector of patient attribute means
    Eigen::MatrixXd s(p-1,p-1);
    Eigen::MatrixXd si(p-1,p-1);
    Eigen::MatrixXd zee(n,p-1);
    Eigen::MatrixXd Z(n,p-1);
    Eigen::MatrixXd Z2(n,p);
    //Update mu as necessary
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(p-1);
    
    //Generates matrix of standard normals, then uses cholesky factorization to get correct covar matrix
    for (int i = 0; i < n; i++){
        for (int j = 0; j < p-1; j++){
            zee(i,j) = nd(generator1);
        }
    }

    for (int i = 0; i < p-1; i++){
        for (int j = 0; j < p-1; j++){
            if (i==j)
            {
                s(i,j) = 1.0;
            }
            else
            {
                s(i,j) = 0.1;
            }
        }
    }

    si = s.inverse();

    double temp = 0;
    
    //Z2 is matrix of patient attributes (n x p)
    Z = zee*s;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){
            if(j > 0){
                Z2(i,j) = Z(i,j-1) + mu(j-1);
            } else{
                Z2(i,j) = 1;
            }
        }
    }
            
    //Eta + Xi computation
    
    double eta [N1];
    double xi [N1];

    for (int i = 0; i < N1; i++){
        eta[i] = nd(generator3);
        xi[i] = xs(generator4);
    }

    //Solving 2-D DP using mesh
    double M = 2000.0;
    double delta = 0.2;
    int steps = ceil(M/delta)+1;

    DP table(n, delta, M);

    //syntax helpers
    double lambda = 0;
    int m = 0;

    double lambda_up = 0;
    double lambda_down = 0;
    double roundup_plus = 0;
    double rounddown_plus = 0;
    double roundup_minus = 0;
    double rounddown_minus = 0;

    double distdown_minus = 0;
    double distup_minus = 0;
    double distdown_plus = 0;
    double distup_plus = 0;

    double plus = 0;
    double minus = 0;
    int screwups = 0;


    //Boundary condition
    for (int i = 0; i < 2*n + 1; i++){
        for (int j = 0; j < steps;j++){
            m = i-n;
            lambda = delta*j;
            table.set(0, i, j, pow(m, 2) + lambda);
        }
    }
    
    //Solving mesh
    for (int l = 1;l < n; l++){
        std::cout << l << "\n";
        for (int i = l-1; i < 2*n+1-l; i++){
            m = i-n;
            for (int j = 0; j < steps; j++){
                lambda = delta*j;
                temp = 0;
                for (int iter = 0; iter < N1; iter++){ 
                    lambda_up = pow(sqrt(lambda)+eta[iter],2)+xi[iter];
                    lambda_down = pow(sqrt(lambda)-eta[iter],2)+xi[iter];
                    
                    if (lambda_down > M){
                        lambda_down = M;
                    }
                    if (lambda_up > M){
                        lambda_up = M;
                    }

                    roundup_plus = ceil(lambda_up/delta);
                    rounddown_plus = floor(lambda_up/delta); 
                    distdown_plus  = lambda_up - rounddown_plus*delta;
                    distup_plus = roundup_plus*delta - lambda_up;

                    roundup_minus = ceil(lambda_down/delta);
                    rounddown_minus = floor(lambda_down/delta);
                    distdown_minus = lambda_down - rounddown_minus*delta;
                    distup_minus = roundup_minus*delta - lambda_down;
                    
                    try{ 
                        plus = (1/delta)*(distdown_plus*table.at(l-1, i+1, rounddown_plus) + distup_plus*table.at(l-1, i+1, roundup_plus));
                        minus = (1/delta)*(distdown_minus*table.at(l-1,i-1,rounddown_minus) + distup_minus*table.at(l-1,i-1,roundup_minus));
                   }
                    catch (const std::out_of_range& e){
                        std::cout << "roundup/down_plus " << roundup_plus << ", " << rounddown_plus << "\n";
                        std::cout << "roundup/down_minus " << roundup_minus << ", " << rounddown_minus << "\n";   
                        std::cout << "Line 151 \n";
                        return 0;
                    }
                    temp = temp*iter/(iter+1) + std::min(plus,minus)/(iter+1);

                    /*if(temp < 0){
                        std::cout << lambda_down << " " << lambda_up << "\n";
                        std::cout << distdown_plus << " " << distup_plus << "\n";
                        std::cout << distdown_minus << " " << distup_minus << "\n";    
                        return 0;
                    }*/
                }
                
                try{
                    table.set(l,i,j,temp);          
                } 
                catch (const std::out_of_range& e){
                    std::cout << "(" << l << ", " << i << ", " << j <<") \n";
                }

            }
        }
    }

    //Vs Naive random allocation
    //Eff = x^T P_Z x where x is allocations, P_Z = I - Z(Z^T Z)^(-1) Z^T

    Eigen::MatrixXd PZ;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Eigen::MatrixXd Z2I;
    
    //Generates matrix of standard normals, then uses cholesky factorization to get correct covar matrix
    //Vector of n/2 1's and -1's
    int rand_x[n];
    for (int i = 0; i < n; i++){
        rand_x[i] = -1 + 2*floor(2*i/n);
    }

    std::vector <int> myvector (rand_x, rand_x+n);
    Eigen::VectorXd random_x(n);
    Eigen::VectorXd dp_x(n);

    double eff_r = 0;
    double eff_dp = 0;
    double eff = 0;

    //Tracks variable Delta as in paper
    Eigen::VectorXd var_Delta(p-1);
    var_Delta.setZero();
    Eigen::VectorXd var_Delta_up(p-1);
    Eigen::VectorXd var_Delta_down(p-1);
    Eigen::MatrixXd condition(n,n); //Used to prevent floating point problems
    int var_delta;
    double mahalanobis_plus = 0;
    double mahalanobis_minus = 0;
    //Eigen::EigenSolver<Eigen::MatrixXd> es;
    
    //Large loop to test policies
    for (int asdf = 0; asdf < N2; asdf++){
        //std::cout << "asdf = " << asdf << "\n";
        var_delta = n+1;
        
        for (int i = 0; i < n; i++){
            for (int j = 0; j < p-1; j++){
                zee(i,j) = nd(generator1);
            }
        }    
        
        //Z is matrix of patient attributes (n x p)  
        Z = zee*s;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < p; j++){
                if(j > 0){
                    Z2(i,j) = Z(i,j-1) + mu(j-1);
                } else{
                    Z2(i,j) = 1;
                }
            }
        }

        Z2I = Z2.transpose()*Z2;
        PZ = I - Z2*(Z2I.inverse())*Z2.transpose(); 
        //es.compute(PZ, false);
        //temp = 0;
        //for (int i = 0; i < n; i++){
        //    if (temp > (double)(es.eigenvalues()(i,0).real())){
        //        temp = (double)(es.eigenvalues()(i,0).real());
        //    }
        //}
        //PZ = PZ + Eigen::MatrixXd::Identity(n,n)*temp;
        //This is where randomized sampling is done 
        for (int i = 0; i < N1; i++){
            std::random_shuffle ( myvector.begin(), myvector.end());
            for (int j = 0; j < n; j++){
                random_x(j) = myvector[j];
            }
            temp = random_x.transpose() * PZ * random_x;
            eff_r = eff_r*i/(i+1) + temp/(i+1);
        }

        //This is where we do "optimal" sampling
        for (int i = 0; i < n; i++){
            var_Delta_up = var_Delta + Z2.block(i,1,1,p-1).transpose();
            var_Delta_down = var_Delta - Z2.block(i,1,1,p-1).transpose();
           
            mahalanobis_plus = var_Delta_up.transpose()*si*var_Delta_up;
            roundup_plus = ceil(mahalanobis_plus/delta);
            rounddown_plus = floor(mahalanobis_plus/delta);
            distup_plus = roundup_plus*delta - mahalanobis_plus;
            distdown_plus = mahalanobis_plus - rounddown_plus*delta;
            //std::cout << var_delta << "\n";
            //std::cout << rounddown_plus << " " << roundup_plus << "\n";
            //std::cout << "(" << n-i-1 << "," << var_delta+1 << ", " << rounddown_plus << ") \n"; 
            //std::cout << "plus done \n";

            mahalanobis_minus = var_Delta_down.transpose()*si*var_Delta_down;
            roundup_minus = ceil(mahalanobis_minus/delta);
            rounddown_minus = floor(mahalanobis_minus/delta);
            distup_minus = roundup_minus*delta - mahalanobis_minus;
            distdown_minus = mahalanobis_minus - rounddown_minus*delta;
            //std::cout << "(" << n-i-1 << "," << var_delta-1 << ", " << rounddown_plus << ") \n";
           
            try{
            plus = (1/delta)*(distdown_plus*table.at(n-i-1,var_delta+1, rounddown_plus) + distup_plus*table.at(n-i-1,var_delta+1,roundup_plus));
           minus = (1/delta)*(distdown_minus*table.at(n-i-1,var_delta-1, rounddown_minus) + distup_minus*table.at(n-i-1,var_delta-1,roundup_minus));
            }
            catch (const std::out_of_range& e){
                std::cout << "mahalanobis_plus = " << mahalanobis_plus << "\n";
                std::cout << "mahalanobis_minus = " << mahalanobis_minus << "\n";
                std::cout << "distdown/up plus =" << distdown_plus << ", " << distup_plus << "\n";
                std::cout << "distdown/up minus =" << distdown_minus << ", " << distup_minus << "\n";
                std::cout << "line 273 \n";
                std::cout << asdf << "\n";
                screwups++;
            }
            if (minus >= plus){
                var_delta++;
                var_Delta = var_Delta_up;
                dp_x(i) = 1;
            } else{
                var_delta--;
                var_Delta = var_Delta_down;
                dp_x(i) = -1;
            }
        }
    
        eff_dp = dp_x.transpose()*PZ*dp_x;
        //std::cout << eff_r << ", " << eff_dp << "\n";
        temp += eff_r/eff_dp;
        eff = eff*asdf/(asdf+1) + (eff_r/eff_dp)/(asdf+1);
    }

    for (int l = 0; l < n; l++){
        for (int i=0; i < 2*n+1;i++){
            std::cout << table.at(l,i,0) << " ";
        }
        std::cout << "\n";
    } 

    std::cout << temp/N1 << "\n";
    std::cout << eff << "\n";
    std::cout << "\n" << PZ << "\n";
    std::cout << "\n" << eff_r << "\n \n" << eff_dp << "\n";
    std::cout << "\n" << dp_x << "\n \n" << dp_x.transpose()*PZ*dp_x;
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(PZ); 
    std::cout << "\n" << es.eigenvalues();
}

