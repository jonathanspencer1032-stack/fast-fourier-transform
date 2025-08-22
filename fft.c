#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
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
    /*
    Inputs: 2 complex number structs
    Outputs: The sum of the two input complex numbers
    */
    complex ADD_COMP;
    ADD_COMP.real = a.real + b.real;
    ADD_COMP.img = a.img + b.img;
    return ADD_COMP;
}

complex SUB_C(complex a, complex b){
    /*
    Inputs: 2 complex number structs
    Outputs: The difference of the two input complex numbers
    */
    complex SUB_COMP;
    SUB_COMP.real = a.real - b.real;
    SUB_COMP.img = a.img - b.img;
    return SUB_COMP;
}

complex MUL_C(complex a, complex b){
    /*
    Inputs: 2 complex number structs
    Outputs: The product of the two input complex numbers
    */

    //(a.real + a.img i)(b.real + b.img i)
    //(a.real * b.real) + (a.real*b.img i) + (a.img * b.real i) + (a.img * b.img * i**2)
    //(a.real * b.real- a.img * b.img) + (a.real*b.img + a.img * b.real) i 
    complex MUL_COMP;
    MUL_COMP.real = (a.real*b.real-a.img*b.img);
    MUL_COMP.img = (a.real*b.img + a.img*b.real);
    return MUL_COMP;
}

double MAG_C(complex a){
    /*
    Inputs: 1 complex number struct
    Outputs: The magnitude of that complex number
    */
    double mag_a = sqrt(pow(a.real,2)+pow(a.img,2));
    return mag_a;
}

double PHASE_C(complex a){
    /*
    Inputs: 1 complex number
    Outputs: The phase of the complex number
    */
    
    double phase_a = atan2(a.img, a.real);
    return phase_a;
}

int reverse_bits(int num, int lgn) {
    /*
    Inputs: -integer to be bit flipped.
            -log(n) where n is the size of the dft.
    Outputs: integer representation of the flipped bits for the input int and size.
    */
    int reverse_num = 0;
    int i;

    for (i = 0; i < lgn; i++) {
        if ((num >> i) & 1) { // Check the i-th bit from the right
            reverse_num |= (1 << ((lgn - 1) - i)); // Set the mirrored bit in reverse_num
        }
    }
    return reverse_num;
}

complex* compute_twiddle_factors(int N){ // compute twiddle factors for roots of unity N.
    /*
    Inputs: size of the DFT.
    Outputs: pointer to the beginning of an array of complex numbers on the heap.
    Note: you must free this at the end.
    */

    complex* W = (complex *)malloc(sizeof(complex)*N); //get a pointer to the front of an N-length array of complex structs.
    if (W == NULL) {
        perror("malloc failed");
        exit(1);
    }

    for(int k = 0; k < N; k++){
        W[k].real = cos(2*M_PI*k/N); 
        W[k].img = -sin(2*M_PI*k/N); //negative for the fft, possitive for ifft.
    }

    return W;
}

complex* FFT(complex* x, int lgn, complex* W){ // take in a vector of length 2**lgn and return its DFT.
    /*
    Inputs: - x: the complex samples.
            - lgn: log2(n) where n is the size of x.
            - W: pointer to an array of twiddle factors.
    Outputs: Complex Array with the fft output.
    Note on memory management: you must free x_out at the end of program.
    Note on scalling the output. This fft scales the output by 1/sqrt(n). The ifft is also scalled by 1/sqrt(n).
    */  
    int n = 1 << lgn; // 2^n get the number of samples.
    int one_sqrt_n = 1/sqrt(n); // for scaling the output.

    complex* x_out = malloc(n * sizeof(complex)); //allocate heap memory for the output.
    
    int bit_reverse_index; 
    for(int i = 0; i < n; i++){ // get the bit reversal lookup table.
        bit_reverse_index = reverse_bits(i, lgn);
        x_out[bit_reverse_index] = x[i]; // map bit reversed array to the input array.
    }

    //Think of these for loops as matrix multiplications. The inmost for loop is the butterfly.
    for(int i = 0; i< lgn; i++){ // O(lgn) loop
        int m = 1 << (i + 1); // 2^(i+1)

        // Get th root of unity twiddle factor

        for(int k = 0; k<n; k = k+m){ // iterate by m.

            for(int j = 0; j < m/2; j++){ //Butterfly opperation.
                complex t = MUL_C(x_out[k+j+m/2], W[j * n / m]); 
                complex u = x_out[k+j];
                
                //how should I scale?
                x_out[k + j] = ADD_C(u, t);
                x_out[k + j + m/2] = SUB_C(u, t);
            }
        }
    }

    return x_out;
}

int main() {
    int N = 8;
    complex* twiddle = compute_twiddle_factors(N);

    complex x[8] = { {0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0} };
    complex* x_fft = FFT(x, 3, twiddle);

    //Check results
    printf("FFT output:\n");
    for(int i = 0; i < N; i++){
        printf("x_fft[%d]: %.5f + %.5fi, |%.5f|\n", 
               i, x_fft[i].real, x_fft[i].img, MAG_C(x_fft[i]));
    }
    free(twiddle);
    free(x_fft);
    return 0;
}