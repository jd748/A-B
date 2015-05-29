#include "DP.hh"
#include <cmath>

DP::DP(int n, double delta, double M)
    : n(n), 
    delta(delta),
    M(M){

    this->num_i = 2*n+1;
    this->num_j = ceil(M/delta)+1;

    int array_size = this->num_i*this->num_j*n;

    this->arr.resize(array_size, 0);

}

int DP::index(int l, int i, int j){
    return l*this->num_i*this->num_j+i*this->num_j+j;
}

double DP::at(int l, int i, int j){
    return this->arr.at(this->index(l,i,j));
} 

void DP::set(int l, int i, int j, double x){
    this->arr.at(this->index(l,i,j)) = x;
}

