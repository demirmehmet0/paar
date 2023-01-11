// The source code of Paar's algorithm(Paar1-Paar2) in our framework is provided by the repository given in https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string>
#include "random"

#define SIZE 32
#define TIME_LIMIT 90000
using namespace std;


mt19937 rand_generator;
default_random_engine generator;
string file;
ifstream TheMatrix;
int DIM;

int NumInputs;
int NumMatrices;
char *flag;

void binprint(long long int x); // output the last NumInputs bits of x

int hamming_weight(uint64_t input)
{
    return __builtin_popcount(input & 0xffffffff) + __builtin_popcount((input >> 32) & 0xffffffff);
}

int paar_algorithm1(std::fstream *f, uint64_t *input_matrix)
{

    int xor_count = 0;
    uint64_t number_of_columns = DIM;
    int hw_max;
    int i_max = 0, j_max = 0, k_max = 0;
    uint64_t tmp;
    int hw;
    uint64_t new_column;
    vector<pair<int, int> > program;

    // compute naive xor count
    for (uint64_t i = 0; i < DIM; i++)
    {
        xor_count += hamming_weight(input_matrix[i]);
    }
    xor_count -= DIM;
    *f << "Naive XOR count: " << xor_count << endl;
    /* cout << "Naive XOR count: " << xor_count << endl; */
    *f << "SLP:" << endl
        << endl;
    /* cout << "SLP:" << endl
         << endl; */
    int say = 0;
    int ctr = DIM;
    do
    {
        std::random_device rd;
        int rand_number = rd() % 5;
        // refresh distribution
        if(rand_number < -1){
            hw_max = 0;
            for (uint64_t i = 0; i < number_of_columns; i++)
            {
                for (uint64_t j = i + 1; j < number_of_columns; j++)
                {
                    for (uint64_t k = j + 1; k < number_of_columns; k++)
                    {
                        tmp = input_matrix[i] & input_matrix[j] & input_matrix[k];
                        hw = hamming_weight(tmp);
                        
                        if (hw > hw_max )
                        {
                            hw_max = hw;
                            i_max = i;
                            j_max = j;
                            k_max = k;
                        }
                    }
                }
            }
            
            if (hw_max > 1)
            {
                new_column = input_matrix[i_max] & input_matrix[j_max] & input_matrix[k_max];
                input_matrix[number_of_columns] = new_column;
                input_matrix[i_max] = (new_column ^ ((1llu << DIM) - 1)) & input_matrix[i_max];
                input_matrix[j_max] = (new_column ^ ((1llu << DIM) - 1)) & input_matrix[j_max];
                input_matrix[k_max] = (new_column ^ ((1llu << DIM) - 1)) & input_matrix[k_max];
                cout << "x" << ctr << " = x" << i_max << " + x" << j_max << " + x" << k_max << endl;
                xor_count -= (hw_max - 1);
                number_of_columns++;
                ctr++;
            }
        }else{
            hw_max = 0;
            for (uint64_t i=0; i<number_of_columns; i++) {
                for (uint64_t j=i+1; j<number_of_columns; j++) {
                    tmp = input_matrix[i] & input_matrix[j];
                    hw = hamming_weight(tmp);
                    if (hw > hw_max) {
                        hw_max = hw;
                        i_max = i;
                        j_max = j;
                    }
                }
            }
            
            if (hw_max > 1) {
                new_column = input_matrix[i_max] & input_matrix[j_max];
                input_matrix[number_of_columns] = new_column;
                input_matrix[i_max] = ( new_column^((1llu << DIM)-1) ) & input_matrix[i_max];
                input_matrix[j_max] = ( new_column^((1llu << DIM)-1) ) & input_matrix[j_max];
                cout << "x" << ctr << " = x" << i_max << " + x" << j_max << endl;
                xor_count -= (hw_max-1);
                number_of_columns++;
                ctr++;
            }
        }
        say++;
    } while (hw_max > 1);
    
    
    for (uint64_t i = 0; i < DIM; i++)
    {
        bool plus_flag = 0;
        *f << endl
            << "y" << i;
        cout << endl
             << "y" << i;
        for (uint64_t j = 0; j < number_of_columns; j++)
        {
            if ((input_matrix[j] & (1ll << (DIM - 1 - i))) != 0)
            {
                if (plus_flag == 0)
                {
                    *f << " = x" << j;
                    cout << " = x" << j;
                    plus_flag = 1;
                }
                else
                {
                    *f << " + x" << j;
                    cout << " + x" << j;
                }
            }
        }
    }

    *f << endl
        << endl;
    /* cout << endl
         << endl; */
    return xor_count;
}

void ReadMatrix(uint64_t *input_matrix)
{
    int NumTargets, NumInputs;
    TheMatrix >> NumTargets;
    TheMatrix >> NumInputs;
    DIM = NumTargets;
    // check that NumInputs is < wordsize
    if (NumInputs >= 8 * sizeof(long long int))
    {
        cout << "too many inputs" << endl;
        exit(0);
    }
    int bit;
    long long int PowerOfTwo = pow(2.0, NumTargets - 1);
    for (int i = 0; i < NumTargets; i++)
    {
        for (int j = 0; j < NumInputs; j++)
        {
            TheMatrix >> bit;
            if (i == 0)
            {
                input_matrix[j] = 0;
            }
            if (bit)
            {
                input_matrix[j] += PowerOfTwo;
            }
        }
        PowerOfTwo = PowerOfTwo / 2;
        TheMatrix.get();
    }
}

int Threshold;

int proccessMatrix(string mode, string file){
    TheMatrix.open(file);
    NumMatrices = 1;

    uint64_t *input_matrix;
    input_matrix = (uint64_t *)malloc((64 + 200) * sizeof(uint64_t));
    string filePath = file;

    string s;
    s.append("matrices/result/");
    s.append("result.");
    s.append(mode);
    s.append(".");
    // split filepath to get filename
    string filename = filePath.substr(filePath.find_last_of("/") + 1);
    s.append(filename);

    //cout << "File: " << s << endl;
    std::fstream f(s, std::fstream::out);
    int result;
    for (int i = 0; i < NumMatrices; i++)
    {
        ReadMatrix(input_matrix);
        
        if (mode.compare("paar1") == 0)
        {
            result = paar_algorithm1(&f, input_matrix);
            f << "Paar1 XOR count:" << result << endl;
            cout << "Paar1 XOR count:" << result << endl;
        }
        f << endl
            << endl
            << endl;
        /* cout << endl
            << endl
             << endl; */
    }
    f.close();
    free(input_matrix);
    return result;
}

int main(int argc, char *argv[])
{

    int i;
    string mode;
    clock_t t1 = clock();
    if (argc < 3)
    {
        cout << "Yanlis kullanim!" << endl;
        return 0;
    }
    mode = argv[1];
    file = argv[2];
    cout << "File: " << file << endl;
    cout << "Mode: " << mode << endl;
    int min = 9999;
    for (int i = 0; i < 5000; i++)
    {
        srand(time(NULL));
        //clean variables
        mode = "paar1";
        file = "./matrix.txt";
        TheMatrix = ifstream();
        DIM = NULL;
        NumInputs;
        NumMatrices;
        //reset random seed
        srand(time(NULL));

        int result = proccessMatrix(mode, file);
        if (result < min)
        {
            min = result;
            if (min == 88)
            {
                exit(0);
            }
        }
    }
    cout << "Min: " << min << endl;
    


    return 0;
}
