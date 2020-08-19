#include <stdlib.h>
#include <stdint.h>
#include <string.h>				// memcmp
#include "bn.h"
#include "bn_lcl.h"
#include "ec_lcl.h"
#include "dut.h"
#include "random.h"

static EC_GROUP *group;
static BIGNUM *order;

const size_t chunk_size = 21;
const size_t number_measurements = 1000;	// per test

uint8_t do_one_computation(uint8_t *data) {
	EC_POINT *r = EC_POINT_new(group);
	BIGNUM *scalar = BN_bin2bn(data, chunk_size, NULL);
	ec_GF2m_simple_mul(group, r, scalar, 0, NULL, NULL, NULL);
	return 0;
}

void init_dut(void) {
}

void prepare_inputs(uint8_t *input_data, uint8_t *classes) {
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *r = BN_new();
	BIGNUM *s = BN_new();
	group = EC_GROUP_new_by_curve_name(NID_sect163r1);
	order = BN_new();
	EC_GROUP_get_order(group, order, NULL);
	printf("\n");

	randombytes(input_data, number_measurements * chunk_size);
	for (size_t i = 0; i < number_measurements; i++) {
		BIGNUM *scalar =
				BN_bin2bn(input_data + (size_t)i * chunk_size, chunk_size,
				NULL);
		BN_nnmod(scalar, scalar, order, ctx);
		BN_add(r, scalar, order);
		BN_add(s, r, order);
		BN_copy(scalar, BN_num_bits(r) > BN_num_bits(order) ? r : s);

		classes[i] = randombit();
		if (classes[i] == 0) {
			BN_set_bit(scalar, BN_num_bits(order) - 1);
		} else {
			BN_clear_bit(scalar, BN_num_bits(order) - 1);
		}

		BN_bn2bin(scalar, input_data + (size_t)i * chunk_size);
		BN_free(scalar);
	}
	BN_free(r);
	BN_free(s);
}
