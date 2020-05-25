
void ep_mul_monty(ep_t r, const ep_t p, const bn_t k) {
	int i, j, bits;
	fp_t tmp, rand;
	ep_t t[2];
	bn_t n, l;

	bn_null(n);
	bn_null(l);
	fp_null(tmp);
	fp_null(rand);
	ep_null(t[0]);
	ep_null(t[1]);

	if (bn_is_zero(k) || ep_is_infty(p)) {
		ep_set_infty(r);
		return;
	}

	TRY {
		bn_new(n);
		bn_new(l);
		fp_new(tmp);
		fp_new(rand);
		ep_new(t[0]);
		ep_new(t[1]);

		ep_curve_get_ord(n);
		bits = bn_bits(n);
		bn_abs(l, k);
		bn_add(l, l, n);
		bn_add(n, l, n);
		dv_swap_cond(l->dp, n->dp, RLC_MAX(l->used, n->used),
			bn_get_bit(l, bits) == 0);
		l->used = RLC_SEL(l->used, n->used, bn_get_bit(l, bits) == 0);

		ep_norm(t[0], p);
		ep_dbl(t[1], t[0]);

#if EP_ADD == PROJC
		/* Blind both points independently. */
		for (i = 0; i < 2; i++) {
			fp_rand(rand);
			fp_mul(t[i]->z, t[i]->z, rand);
			fp_sqr(tmp, rand);
			fp_mul(t[i]->x, t[i]->x, tmp);
			fp_mul(rand, rand, tmp);
			fp_mul(t[i]->y, t[i]->y, rand);
		}
#endif

		for (i = bits - 1; i >= 0; i--) {
			j = bn_get_bit(l, i);
			dv_swap_cond(t[0]->x, t[1]->x, RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y, t[1]->y, RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z, t[1]->z, RLC_FP_DIGS, j ^ 1);
			ep_add(t[0], t[0], t[1]);
			ep_dbl(t[1], t[1]);
			dv_swap_cond(t[0]->x, t[1]->x, RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->y, t[1]->y, RLC_FP_DIGS, j ^ 1);
			dv_swap_cond(t[0]->z, t[1]->z, RLC_FP_DIGS, j ^ 1);
		}

		ep_norm(r, t[0]);
		if (bn_sign(k) == RLC_NEG) {
			ep_neg(r, r);
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(n);
		bn_free(l);
		fp_free(tmp);
		fp_free(rand);
		ep_free(t[1]);
		ep_free(t[0]);
	}
}
