#include <stdlib.h>
#include <stdint.h>
#include <string.h> // memcmp
#include "bn.h"
#include "bn_lcl.h"
#include "ec_lcl.h"
#include "dut.h"
#include "random.h"

static EC_GROUP *group;
static BIGNUM *order;

const size_t chunk_size = 20;
const size_t number_measurements = 10000; // per test

uint8_t do_one_computation(uint8_t *data) {
    EC_POINT *r = EC_POINT_new(group);
    BIGNUM *scalar = BN_bin2bn(data, chunk_size, NULL);
    BN_nnmod(scalar, scalar, order, NULL);
    ec_GF2m_simple_mul(group, r, scalar, 0, NULL, NULL, NULL);
    return 0;
}

void init_dut(void) {
    int function_status = -1;
    group = EC_GROUP_new_by_curve_name(NID_sect163r1);
    if (NULL == group)
    {
        printf("Failed to create new EC Group\n");
        function_status = -1;
    }
    
    order = BN_new();
    EC_GROUP_get_order(group, order, NULL);
}

void prepare_inputs(uint8_t *input_data, uint8_t *classes) {
  randombytes(input_data, number_measurements * chunk_size);
  for (size_t i = 0; i < number_measurements; i++) {
    classes[i] = randombit();
  }
}

