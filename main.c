#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
/*
This program takes in a signal l of 2**n length and returns the FFT(l). My main resource for learning about 
FFT has been the book "Lectures on the Fourier Transform and its Application" by Brad G. Osgood. This resource
has been very helpful in understanding the DFT, CFT, and FFT. 

*/

/*
Break down of the FFT.

Divide and conquer nature. O(n*log(n))

Take in a signal l of radix-2 (length(l) = 2**n for some natural number n).
Sort the indices using bit reversal


*/

//struct to hold complex numbers.
typedef struct { 
    double real;
    double img;
} complex;


//helper functions for complex numbers.
complex ADD_C(complex a, complex b){
    complex ADD_COMP;
    ADD_COMP.real = a.real + b.real;
    ADD_COMP.img = a.img + b.img;
    return ADD_COMP;
}

complex SUB_C(complex a, complex b){
    complex SUB_COMP;
    SUB_COMP.real = a.real - b.real;
    SUB_COMP.img = a.img - b.img;
    return SUB_COMP;
}

complex MUL_C(complex a, complex b){
    //(a.real + a.img i)(b.real + b.img i)
    //(a.real * b.real) + (a.real*b.img i) + (a.img * b.real i) + (a.img * b.img * i**2)
    //(a.real * b.real- a.img * b.img) + (a.real*b.img + a.img * b.real) i 
    complex MUL_COMP;
    MUL_COMP.real = (a.real*b.real-a.img*b.img);
    MUL_COMP.img = (a.real*b.img + a.img*b.real);
    return MUL_COMP;
}

double MAG_C(complex a){
    //get the magnitude of complex number a.
    double mag_a = sqrt(pow(a.real,2)+pow(a.img,2));
    return mag_a;
}

double PHASE_C(complex a){
    //get the phase of a complex number a. TOA: tangent = opposite / adjacent. 
    double toa = a.img / a.real;
    double phase_a = atan(toa);
    return phase_a;
}

int reverse_bits(int num, int n_powers) {
    int reverse_num = 0;
    int i;

    for (i = 0; i < n_powers; i++) {
        if ((num >> i) & 1) { // Check the i-th bit from the right
            reverse_num |= (1 << ((n_powers - 1) - i)); // Set the mirrored bit in reverse_num
        }
    }
    return reverse_num;
}



int main(){
    int bits = 1;
    int powers = 3; //2**3 = 8

    printf("Reverse bit of %u is %u\n", bits, reverse_bits(bits, powers));

    return 0;
}