#pragma once
#include "examples.h"

using namespace std;
using namespace seal;

void test_matrix_initialization();
void test_matrix_init_partial();
void test_matconv();
void test_print_iter();
void test_innerprod();
void test_conviter();
void test_innerprod_vector();
void test_init_matrix();
void test_init_matrix_uint();
void test_init_matrix_uint_by_kernel();
void test_matrix_conversion_with_coeffiter();
void test_secret_product();
void test_matrix_conversion_with_rnsiter();
void test_matrix_dot_product();
void test_matrix_dot_product_uint(uint64_t coeff_degree);
void test_matrix_dot_product_uint_t(uint64_t coeff_degree);
void test_conv_nega();
