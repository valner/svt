#include <fstream>
#include <iostream>
#include <string.h>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

// scalar multiply of two matrix rows
double 
mult_row_by_row(double *row1, double *row2, int size) {
    double sum = 0.;
    
    for (int i = 0; i < size; i++)
        sum += row1[i] * row2[i];

    return sum;
}

// generate square positive definite matrix dim (A) = size and store it on disk 
int 
generate_matrix(const int size, char* outfile, const int seed) {
    int i, j;
    double **matrix, **output;
    double sum;

    if (0 == size) return 0; // empty matrix -- do nothing

    matrix = new double*[size]; // allocate space for initial matrix
    for (i = 0; i < size; i++)
        matrix[i] = new double [size];
    // pseudo-random initialization
    srand48(seed);
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
	    matrix[i][j] = drand48();

    // compute C = A A^T
    output = new double*[size];
    for (i = 0; i < size; i++)
        output[i] = new double [size];
    
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++) {
            output[i][j] = mult_row_by_row(matrix[i], matrix[j], size); 
        }
    for (i = 0; i < size; i++)
        delete matrix[i];

    delete[] matrix;
    // store matrix as binary file
    ofstream file(outfile);
    if (file.fail()) {
        throw std::logic_error("Cannot open file for write or smth.");
    }
    // format: size; {double values for rows}
    file.write((const char*) &size, sizeof(size));
    for (i = 0; i < size; i++)
        file.write ((const char*) output[i], sizeof(output[i][0])* size);
    file.close();
    for (i = 0; i < size; i++)
        delete output[i];
    delete[] output; 
    return 0;
}
// function to read matrix from file and return pointer
/* 
double** 
read_matrix(const char *infile, int &outsize) {
    int i,j;
    int size;
    FILE *fd;
    double **outmatrix;

    fd = fopen(infile, "r");
    fread(&size, sizeof(size), 1, fd);
    outsize = size;
    outmatrix = (double**) malloc (size * sizeof(double*));
    
    for (i = 0; i < size; i++) {
        outmatrix[i] = (double*) malloc (size * sizeof(double));
        fread(outmatrix, sizeof(double), size, fd);
    }

    return outmatrix; 

}
*/
 // ===========================================================================
int main(int argc, char *argv[]) {
    if (argc != 3) {
        throw std::logic_error("Wrong command line. Usage: ./generate_spd size outfile.");
    }
    generate_matrix(atoi(argv[1]), argv[2], 100);    
    return 0;
}
