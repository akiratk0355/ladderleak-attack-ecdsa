
/*-
 * This functions computes (in constant time) a point multiplication over the
 * EC group.
 *
 * At a high level, it is Montgomery ladder with conditional swaps.
 *
 * It performs either a fixed scalar point multiplication
 *          (scalar * generator)
 * when point is NULL, or a generic scalar point multiplication
 *          (scalar * point)
 * when point is not NULL.
 *
 * scalar should be in the range [0,n) otherwise all constant time bets are off.
 *
 * NB: This says nothing about EC_POINT_add and EC_POINT_dbl,
 * which of course are not constant time themselves.
 *
 * The product is stored in r.
 *
 * Returns 1 on success, 0 otherwise.
 */
static int ec_mul_consttime(const EC_GROUP *group, EC_POINT *r,
                            const BIGNUM *scalar, const EC_POINT *point,
                            BN_CTX *ctx)
{
    int i, cardinality_bits, group_top, kbit, pbit, Z_is_one;
    EC_POINT *s = NULL;
    BIGNUM *k = NULL;
    BIGNUM *lambda = NULL;
    BIGNUM *cardinality = NULL;
    BN_CTX *new_ctx = NULL;
    int ret = 0;

    if (ctx == NULL && (ctx = new_ctx = BN_CTX_new()) == NULL)
        return 0;

    BN_CTX_start(ctx);

    s = EC_POINT_new(group);
    if (s == NULL)
        goto err;

    if (point == NULL) {
        if (!EC_POINT_copy(s, group->generator))
            goto err;
    } else {
        if (!EC_POINT_copy(s, point))
            goto err;
    }

    EC_POINT_BN_set_flags(s, BN_FLG_CONSTTIME);

    cardinality = BN_CTX_get(ctx);
    lambda = BN_CTX_get(ctx);
    k = BN_CTX_get(ctx);
    if (k == NULL || !BN_mul(cardinality, &group->order, &group->cofactor, ctx))
        goto err;

    /*
     * Group cardinalities are often on a word boundary.
     * So when we pad the scalar, some timing diff might
     * pop if it needs to be expanded due to carries.
     * So expand ahead of time.
     */
    cardinality_bits = BN_num_bits(cardinality);
    group_top = cardinality->top;
    if ((bn_wexpand(k, group_top + 2) == NULL)
        || (bn_wexpand(lambda, group_top + 2) == NULL))
        goto err;

    if (!BN_copy(k, scalar))
        goto err;

    BN_set_flags(k, BN_FLG_CONSTTIME);

    if ((BN_num_bits(k) > cardinality_bits) || (BN_is_negative(k))) {
        /*-
         * this is an unusual input, and we don't guarantee
         * constant-timeness
         */
        if (!BN_nnmod(k, k, cardinality, ctx))
            goto err;
    }

    if (!BN_add(lambda, k, cardinality))
        goto err;
    BN_set_flags(lambda, BN_FLG_CONSTTIME);
    if (!BN_add(k, lambda, cardinality))
        goto err;
    /*
     * lambda := scalar + cardinality
     * k := scalar + 2*cardinality
     */
    kbit = BN_is_bit_set(lambda, cardinality_bits);
    BN_consttime_swap(kbit, k, lambda, group_top + 2);

    group_top = group->field.top;
    if ((bn_wexpand(&s->X, group_top) == NULL)
        || (bn_wexpand(&s->Y, group_top) == NULL)
        || (bn_wexpand(&s->Z, group_top) == NULL)
        || (bn_wexpand(&r->X, group_top) == NULL)
        || (bn_wexpand(&r->Y, group_top) == NULL)
        || (bn_wexpand(&r->Z, group_top) == NULL))
        goto err;

    if (group->meth->field_type == NID_X9_62_prime_field) {
        /*
         * Apply blinding independently to r and s: use point r as temporary
         * variable.
         *
         * Blinding requires using a projective representation, and later we
         * implement the ladder step using EC_POINT_add() and EC_POINT_dbl():
         * this requires EC_METHODs that support projective coordinates in their
         * point addition and point doubling implementations.
         *
         * All the built-in EC_METHODs for prime curves at the moment satisfy
         * this condition, while the EC_METHOD for binary curves does not.
         *
         * This is why we check for a prime field to apply blinding.
         */

        /* first randomize r->Z to blind s. */
        do {
            if (!BN_rand_range(&r->Z, &group->field)) {
                ECerr(EC_F_EC_MUL_CONSTTIME, ERR_R_BN_LIB);
                goto err;
            }
        } while (BN_is_zero(&r->Z));

        /* convert r->Z to the correct field representation. */
        if (group->meth->field_encode != NULL
            && !group->meth->field_encode(group, &r->Z, &r->Z, ctx)) {
            goto err;
        }

        /* scale s->X and s->Y by r->Z^2 and r->Z^3, respectively. */
        if (!group->meth->field_sqr(group, &r->X, &r->Z, ctx)
            || !group->meth->field_mul(group, &r->Y, &r->X, &r->Z, ctx)) {
            goto err;
        }
        if (!group->meth->field_mul(group, &s->X, &s->X, &r->X, ctx)
            || !group->meth->field_mul(group, &s->Y, &s->Y, &r->Y, ctx)
            || !group->meth->field_mul(group, &s->Z, &s->Z, &r->Z, ctx)) {
            goto err;
        }
        /* mark the flag in s to full projective coordinates. */
        s->Z_is_one = 0;

        /* blinding: now rerandomize r->Z to make r a blinded copy of s. */
        do {
            if (!BN_rand_range(&r->Z, &group->field)) {
                ECerr(EC_F_EC_MUL_CONSTTIME, ERR_R_BN_LIB);
                goto err;
            }
        } while (BN_is_zero(&r->Z));

        /* convert r->Z to the correct field representation. */
        if (group->meth->field_encode != NULL
            && !group->meth->field_encode(group, &r->Z, &r->Z, ctx)) {
            goto err;
        }

        /* scale r->X and r->Y by r->Z^2 and r->Z^3, respectively. */
        if (!group->meth->field_sqr(group, &r->X, &r->Z, ctx)
            || !group->meth->field_mul(group, &r->Y, &r->X, &r->Z, ctx)) {
            goto err;
        }
        if (!group->meth->field_mul(group, &r->X, &s->X, &r->X, ctx)
            || !group->meth->field_mul(group, &r->Y, &s->Y, &r->Y, ctx)
            || !group->meth->field_mul(group, &r->Z, &s->Z, &r->Z, ctx)) {
            goto err;
        }
        /* mark the flag in s to full projective coordinates. */
        r->Z_is_one = 0;
    } else {
        /* top bit is a 1, in a fixed pos; binary curves are not blinded here */
        if (!EC_POINT_copy(r, s))
            goto err;
    }

    EC_POINT_BN_set_flags(r, BN_FLG_CONSTTIME);

    if (!EC_POINT_dbl(group, s, s, ctx))
        goto err;

    pbit = 0;

#define EC_POINT_CSWAP(c, a, b, w, t) do {         \
        BN_consttime_swap(c, &(a)->X, &(b)->X, w); \
        BN_consttime_swap(c, &(a)->Y, &(b)->Y, w); \
        BN_consttime_swap(c, &(a)->Z, &(b)->Z, w); \
        t = ((a)->Z_is_one ^ (b)->Z_is_one) & (c); \
        (a)->Z_is_one ^= (t);                      \
        (b)->Z_is_one ^= (t);                      \
} while(0)

    /*-
     * The ladder step, with branches, is
     *
     * k[i] == 0: S = add(R, S), R = dbl(R)
     * k[i] == 1: R = add(S, R), S = dbl(S)
     *
     * Swapping R, S conditionally on k[i] leaves you with state
     *
     * k[i] == 0: T, U = R, S
     * k[i] == 1: T, U = S, R
     *
     * Then perform the ECC ops.
     *
     * U = add(T, U)
     * T = dbl(T)
     *
     * Which leaves you with state
     *
     * k[i] == 0: U = add(R, S), T = dbl(R)
     * k[i] == 1: U = add(S, R), T = dbl(S)
     *
     * Swapping T, U conditionally on k[i] leaves you with state
     *
     * k[i] == 0: R, S = T, U
     * k[i] == 1: R, S = U, T
     *
     * Which leaves you with state
     *
     * k[i] == 0: S = add(R, S), R = dbl(R)
     * k[i] == 1: R = add(S, R), S = dbl(S)
     *
     * So we get the same logic, but instead of a branch it's a
     * conditional swap, followed by ECC ops, then another conditional swap.
     *
     * Optimization: The end of iteration i and start of i-1 looks like
     *
     * ...
     * CSWAP(k[i], R, S)
     * ECC
     * CSWAP(k[i], R, S)
     * (next iteration)
     * CSWAP(k[i-1], R, S)
     * ECC
     * CSWAP(k[i-1], R, S)
     * ...
     *
     * So instead of two contiguous swaps, you can merge the condition
     * bits and do a single swap.
     *
     * k[i]   k[i-1]    Outcome
     * 0      0         No Swap
     * 0      1         Swap
     * 1      0         Swap
     * 1      1         No Swap
     *
     * This is XOR. pbit tracks the previous bit of k.
     */

    for (i = cardinality_bits - 1; i >= 0; i--) {
        kbit = BN_is_bit_set(k, i) ^ pbit;
        EC_POINT_CSWAP(kbit, r, s, group_top, Z_is_one);
        if (!EC_POINT_add(group, s, r, s, ctx))
            goto err;
        if (!EC_POINT_dbl(group, r, r, ctx))
            goto err;
        /*
         * pbit logic merges this cswap with that of the
         * next iteration
         */
        pbit ^= kbit;
    }
    /* one final cswap to move the right value into r */
    EC_POINT_CSWAP(pbit, r, s, group_top, Z_is_one);
#undef EC_POINT_CSWAP

    ret = 1;

 err:
    EC_POINT_clear_free(s);
    BN_CTX_end(ctx);
    BN_CTX_free(new_ctx);

    return ret;
}
