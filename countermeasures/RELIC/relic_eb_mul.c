
void eb_mul_lodah(eb_t r, const eb_t p, const bn_t k) {
	int bits, i, j;
	dv_t x1, z1, x2, z2, r1, r2, r3, r4, r5;
	const dig_t *b;
	bn_t t, n;

	if (bn_is_zero(k)) {
		eb_set_infty(r);
		return;
	}

	bn_null(n);
	bn_null(t);
	dv_null(x1);
	dv_null(z1);
	dv_null(x2);
	dv_null(z2);
	dv_null(r1);
	dv_null(r2);
	dv_null(r3);
	dv_null(r4);
	dv_null(r5);

	TRY {
		bn_new(n);
		bn_new(t);
		dv_new(x1);
		dv_new(z1);
		dv_new(x2);
		dv_new(z2);
		dv_new(r1);
		dv_new(r2);
		dv_new(r3);
		dv_new(r4);
		dv_new(r5);

		fb_sqr(z2, p->x);
		fb_sqr(x2, z2);
		dv_zero(r5, 2 * RLC_FB_DIGS);

		b = eb_curve_get_b();
		eb_curve_get_ord(n);
		bits = bn_bits(n);

		bn_abs(t, k);
		bn_add(t, t, n);
		bn_add(n, t, n);
		dv_swap_cond(t->dp, n->dp, RLC_MAX(t->used, n->used),
			bn_get_bit(t, bits) == 0);
		t->used = RLC_SEL(t->used, n->used, bn_get_bit(t, bits) == 0);

		switch (eb_curve_opt_b()) {
			case RLC_ZERO:
				break;
			case RLC_ONE:
				fb_add_dig(x2, x2, (dig_t)1);
				break;
			case RLC_TINY:
				fb_add_dig(x2, x2, b[0]);
				break;
			default:
				fb_addn_low(x2, x2, b);
				break;
		}

		/* Blind both points independently. */
		fb_rand(z1);
		fb_mul(x1, z1, p->x);
		fb_rand(r1);
		fb_mul(z2, z2, r1);
		fb_mul(x2, x2, r1);

		for (i = bits - 1; i >= 0; i--) {
			j = bn_get_bit(t, i);
			fb_mul(r1, x1, z2);
			fb_mul(r2, x2, z1);
			fb_add(r3, r1, r2);
			fb_muln_low(r4, r1, r2);
			dv_swap_cond(x1, x2, RLC_FB_DIGS, j ^ 1);
			dv_swap_cond(z1, z2, RLC_FB_DIGS, j ^ 1);
			fb_sqr(z1, r3);
			fb_muln_low(r1, z1, p->x);
			fb_addd_low(x1, r1, r4, 2 * RLC_FB_DIGS);
			fb_rdcn_low(x1, x1);
			fb_sqr(r1, z2);
			fb_sqr(r2, x2);
			fb_mul(z2, r1, r2);
			switch (eb_curve_opt_b()) {
				case RLC_ZERO:
					fb_sqr(x2, r2);
					break;
				case RLC_ONE:
					fb_add(r1, r1, r2);
					fb_sqr(x2, r1);
					break;
				case RLC_TINY:
					fb_sqr(r1, r1);
					fb_sqrl_low(x2, r2);
					fb_mul1_low(r5, r1, b[0]);
					fb_addd_low(x2, x2, r5, RLC_FB_DIGS + 1);
					fb_rdcn_low(x2, x2);
					break;
				default:
					fb_sqr(r1, r1);
					fb_sqrl_low(x2, r2);
					fb_muln_low(r5, r1, b);
					fb_addd_low(x2, x2, r5, 2 * RLC_FB_DIGS);
					fb_rdcn_low(x2, x2);
					break;
			}
			dv_swap_cond(x1, x2, RLC_FB_DIGS, j ^ 1);
			dv_swap_cond(z1, z2, RLC_FB_DIGS, j ^ 1);
		}

		if (fb_is_zero(z1)) {
			/* The point q is at infinity. */
			eb_set_infty(r);
		} else {
			if (fb_is_zero(z2)) {
				fb_copy(r->x, p->x);
				fb_add(r->y, p->x, p->y);
				fb_set_dig(r->z, 1);
			} else {
				/* r3 = z1 * z2. */
				fb_mul(r3, z1, z2);
				/* z1 = (x1 + x * z1). */
				fb_mul(z1, z1, p->x);
				fb_add(z1, z1, x1);
				/* z2 = x * z2. */
				fb_mul(z2, z2, p->x);
				/* x1 = x1 * z2. */
				fb_mul(x1, x1, z2);
				/* z2 = (x2 + x * z2)(x1 + x * z1). */
				fb_add(z2, z2, x2);
				fb_mul(z2, z2, z1);

				/* r4 = (x^2 + y) * z1 * z2 + (x2 + x * z2)(x1 + x * z1). */
				fb_sqr(r4, p->x);
				fb_add(r4, r4, p->y);
				fb_mul(r4, r4, r3);
				fb_add(r4, r4, z2);

				/* r3 = (z1 * z2 * x)^{-1}. */
				fb_mul(r3, r3, p->x);
				fb_inv(r3, r3);
				/* r4 = (x^2 + y) * z1 * z2 + (x2 + x * z2)(x1 + x * z1) * r3. */
				fb_mul(r4, r4, r3);
				/* x2 = x1 * x * z2 * (z1 * z2 * x)^{-1} = x1/z1. */
				fb_mul(x2, x1, r3);
				/* z2 = x + x1/z1. */
				fb_add(z2, x2, p->x);

				/* z2 = z2 * r4 + y. */
				fb_mul(z2, z2, r4);
				fb_add(z2, z2, p->y);

				fb_copy(r->x, x2);
				fb_copy(r->y, z2);
				fb_set_dig(r->z, 1);
			}
		}

		r->norm = 1;
		if (bn_sign(k) == RLC_NEG) {
			eb_neg(r, r);
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(n);
		bn_free(t);
		dv_free(x1);
		dv_free(z1);
		dv_free(x2);
		dv_free(z2);
		dv_free(r1);
		dv_free(r2);
		dv_free(r3);
		dv_free(r4);
		dv_free(r5);
	}
}
