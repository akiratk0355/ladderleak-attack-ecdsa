
/*-
 * Computes scalar*point and stores the result in r.
 * point can not equal r.
 * Uses a modified algorithm 2P of
 *     Lopez, J. and Dahab, R.  "Fast multiplication on elliptic curves over
 *     GF(2^m) without precomputation" (CHES '99, LNCS 1717).
 *
 * To protect against side-channel attack the function uses constant time swap,
 * avoiding conditional branches.
 */
static int ec_GF2m_montgomery_point_multiply(const EC_GROUP *group,
                                             EC_POINT *r,
                                             const BIGNUM *scalar,
                                             const EC_POINT *point,
                                             BN_CTX *ctx)
{
    BIGNUM *x1, *x2, *z1, *z2;
    int ret = 0, i, group_top;
    BN_ULONG mask, word;

    if (r == point) {
        ECerr(EC_F_EC_GF2M_MONTGOMERY_POINT_MULTIPLY, EC_R_INVALID_ARGUMENT);
        return 0;
    }

    /* if result should be point at infinity */
    if ((scalar == NULL) || BN_is_zero(scalar) || (point == NULL) ||
        EC_POINT_is_at_infinity(group, point)) {
        return EC_POINT_set_to_infinity(group, r);
    }

    /* only support affine coordinates */
    if (!point->Z_is_one)
        return 0;

    /*
     * Since point_multiply is static we can guarantee that ctx != NULL.
     */
    BN_CTX_start(ctx);
    x1 = BN_CTX_get(ctx);
    z1 = BN_CTX_get(ctx);
    if (z1 == NULL)
        goto err;

    x2 = &r->X;
    z2 = &r->Y;

    group_top = group->field.top;
    if (bn_wexpand(x1, group_top) == NULL
        || bn_wexpand(z1, group_top) == NULL
        || bn_wexpand(x2, group_top) == NULL
        || bn_wexpand(z2, group_top) == NULL)
        goto err;

    if (!BN_GF2m_mod_arr(x1, &point->X, group->poly))
        goto err;               /* x1 = x */
    if (!BN_one(z1))
        goto err;               /* z1 = 1 */
    if (!group->meth->field_sqr(group, z2, x1, ctx))
        goto err;               /* z2 = x1^2 = x^2 */
    if (!group->meth->field_sqr(group, x2, z2, ctx))
        goto err;
    if (!BN_GF2m_add(x2, x2, &group->b))
        goto err;               /* x2 = x^4 + b */

    /* blinding: make sure z1 and z2 are independently blinded. */
    do {
        if (!BN_rand(z1, BN_num_bits(&group->field) - 1, -1, 0)) {
            ECerr(EC_F_EC_GF2M_MONTGOMERY_POINT_MULTIPLY, ERR_R_BN_LIB);
            return 0;
        }
    } while (BN_is_zero(z1));

    /* first blind (x2,z2) using z1 as the random field element. */
    if ((group->meth->field_encode != NULL
         && !group->meth->field_encode(group, z1, z1, ctx))
        || !group->meth->field_mul(group, x2, x2, z1, ctx)
        || !group->meth->field_mul(group, z2, z2, z1, ctx))
        return 0;

    /* now generate another random field element to blind (x1,z1) */
    do {
        if (!BN_rand(z1, BN_num_bits(&group->field) - 1, -1, 0)) {
            ECerr(EC_F_EC_GF2M_MONTGOMERY_POINT_MULTIPLY, ERR_R_BN_LIB);
            return 0;
        }
    } while (BN_is_zero(z1));

    if ((group->meth->field_encode != NULL
         && !group->meth->field_encode(group, z1, z1, ctx))
        || !group->meth->field_mul(group, x1, x1, z1, ctx))
        return 0;

    /* find top most bit and go one past it */
    i = scalar->top - 1;
    mask = BN_TBIT;
    word = scalar->d[i];
    while (!(word & mask))
        mask >>= 1;
    mask >>= 1;
    /* if top most bit was at word break, go to next word */
    if (!mask) {
        i--;
        mask = BN_TBIT;
    }

    for (; i >= 0; i--) {
        word = scalar->d[i];
        while (mask) {
            BN_consttime_swap(word & mask, x1, x2, group_top);
            BN_consttime_swap(word & mask, z1, z2, group_top);
            if (!gf2m_Madd(group, &point->X, x2, z2, x1, z1, ctx))
                goto err;
            if (!gf2m_Mdouble(group, x1, z1, ctx))
                goto err;
            BN_consttime_swap(word & mask, x1, x2, group_top);
            BN_consttime_swap(word & mask, z1, z2, group_top);
            mask >>= 1;
        }
        mask = BN_TBIT;
    }

    /* convert out of "projective" coordinates */
    i = gf2m_Mxy(group, &point->X, &point->Y, x1, z1, x2, z2, ctx);
    if (i == 0)
        goto err;
    else if (i == 1) {
        if (!EC_POINT_set_to_infinity(group, r))
            goto err;
    } else {
        if (!BN_one(&r->Z))
            goto err;
        r->Z_is_one = 1;
    }

    /* GF(2^m) field elements should always have BIGNUM::neg = 0 */
    BN_set_negative(&r->X, 0);
    BN_set_negative(&r->Y, 0);

    ret = 1;

 err:
    BN_CTX_end(ctx);
    return ret;
}
