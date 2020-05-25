#include <inttypes.h>
#include <stdio.h>
#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "mocksig.h"

Domain mock::setup(uint32_t n_bit, mpz_class n) {
	Domain pp = { n_bit, n, };

	return pp;
}

mpz_class mock::keygen(Domain pp, gmp_randclass& rand) {
	mpz_class d = rand.get_z_range(pp.n);
	return d;
}

SignatureLeak mock::sign(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand) {
	mpz_class r = rand.get_z_range(pp.n);
	mpz_class rr;
	mpz_mod_ui(rr.get_mpz_t(), r.get_mpz_t(), 1 << leak);

	mpz_class h = rand.get_z_range(pp.n);
	mpz_class s;
	mpz_class tmp = r - h * d;
	mpz_mod(s.get_mpz_t(), tmp.get_mpz_t(), pp.n.get_mpz_t());

	SignatureLeak sigma(h, s, rr);

	return sigma;
}

SignatureLeak mock::sign_msb_leak(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand) {
	mpz_class r = rand.get_z_range(pp.n);
	mpz_class rr;
	if (r < (pp.n >> leak))
		rr = mpz_class(0);
	else
		rr = mpz_class(1);

	mpz_class h = rand.get_z_range(pp.n);
	mpz_class s;
	mpz_class tmp = r - h * d;
	mpz_mod(s.get_mpz_t(), tmp.get_mpz_t(), pp.n.get_mpz_t());

	SignatureLeak sigma(h, s, rr);

	return sigma;
}

SignatureSimple mock::sign_filter(Domain pp, mpz_class d, int M, uint32_t leak, uint32_t filter, gmp_randclass& rand) {
	mpz_class nonce_lim = pp.n >> leak;
	mpz_class hlim = pp.n >> filter;
	mpz_class r, h, s, stmp;
	r = rand.get_z_range(nonce_lim);
	h = rand.get_z_range(hlim);
	stmp = r - h * d;
	mpz_mod(s.get_mpz_t(), stmp.get_mpz_t(), pp.n.get_mpz_t());

	SignatureSimple sigma(h, s);

	return sigma;
}

/* ---------------- utilities ---------------- */

/* Counting cycles */
unsigned long long rdtsc(void) {
        unsigned hi, lo;
        __asm__ volatile ("rdtsc" : "=a"(lo), "=d"(hi));
        return ((unsigned long long)lo)|(((unsigned long long)hi)<<32);
}

void sigprint(SignatureSimple sig) {
	std::cout << "(h=" << sig.h << ", s=" << sig.s << ")" << std::endl;
}

void sigprint(SignatureSC25519 sig) {
	mpz_class h, s;
	gs_to_mpz(sig.h, h);
	gs_to_mpz(sig.s, s);
	gmp_printf("(h=%Zd, s=%Zd)\n", h.get_mpz_t(), s.get_mpz_t());
}

void sigvprint(std::vector<SignatureSimple>& sigs, uint32_t idx_start, uint32_t idx_end) {
	unsigned int pos = 0;
	for (const SignatureSimple& s : sigs) {
		if (idx_start <= pos && pos <= idx_end)
			sigprint(s);
		pos += 1;
	}
}

void sigsave(std::vector<SignatureSimple>& sigs, std::string filename, bool str) {
	FILE *fp = fopen(filename.c_str(), "wb");

	for (auto& sig : sigs) {
		if (!str) {
			mpz_out_raw(fp, sig.h.get_mpz_t());
			mpz_out_raw(fp, sig.s.get_mpz_t());
		} else {
			mpz_out_str(fp, 10, sig.h.get_mpz_t());
			fprintf(fp, ",");
			mpz_out_str(fp, 10, sig.s.get_mpz_t());
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void sigsave(std::vector<SignatureSC25519>& sigs, std::string filename, bool str) {
	FILE *fp = fopen(filename.c_str(), "wb");

	mpz_class h, s;
	for (auto& sig : sigs) {
		gs_to_mpz(sig.h, h);
		gs_to_mpz(sig.s, s);
		if (!str) {
			mpz_out_raw(fp, h.get_mpz_t());
			mpz_out_raw(fp, s.get_mpz_t());
		} else {
			mpz_out_str(fp, 10, h.get_mpz_t());
			fprintf(fp, ",");
			mpz_out_str(fp, 10, s.get_mpz_t());
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void sigload(std::vector<SignatureSimple>& sigs, std::string filename, uint32_t lim, bool str, mpz_class bound, uint32_t resume_offset) {
	FILE *fp = fopen(filename.c_str(), "rb");
	mpz_class h, s;
	uint32_t size = sigs.size();
	uint32_t ctr = 0;

	if (!str) {
		while (size < lim) {
			if (mpz_inp_raw(h.get_mpz_t(), fp) == MPZ_LOAD_ERROR)
				break;
			if (mpz_inp_raw(s.get_mpz_t(), fp) == MPZ_LOAD_ERROR)
				break;
			if (ctr < resume_offset) {
				ctr++;
				continue;
			}
			sigs.emplace_back(h, s);
			if (bound != 0 && h >= bound)
				gmp_printf("WARNING: found h larger than ub=%Zd at %lu\n", bound.get_mpz_t(), ctr);
			size++;
			ctr++;
			//sigprint(SignatureSimple(h, s));
		}
	} else {
		while (size < lim) {
			if (mpz_inp_str(h.get_mpz_t(), fp, 10) == MPZ_LOAD_ERROR)
				break;
			if (mpz_inp_str(s.get_mpz_t(), fp, 10) == MPZ_LOAD_ERROR)
				break;
			if (ctr < resume_offset) {
				ctr++;
				continue;
			}
			sigs.emplace_back(h, s);
			if (bound != 0 && h >= bound)
				gmp_printf("WARNING: found h larger than ub=%Zd at %lu\n", bound.get_mpz_t(), ctr);
			size++;
			ctr++;
			//sigprint(SignatureSimple(h, s));
		}
	}
	fclose(fp);
}
;

void sigload(std::vector<SignatureSC25519>& sigs, std::string filename) {
	FILE *fp = fopen(filename.c_str(), "rb");
	mpz_class h, s;

	while (1) {
		if (mpz_inp_raw(h.get_mpz_t(), fp) == MPZ_LOAD_ERROR)
			break;
		if (mpz_inp_raw(s.get_mpz_t(), fp) == MPZ_LOAD_ERROR)
			break;
		sc25519 temph, temps;
		mpz_to_gs(temph, h);
		mpz_to_gs(temps, s);

		sigs.emplace_back(temph, temps);
		//sigprint(SignatureSimple(h, s));
	}
	fclose(fp);
}
;

void countbias(std::vector<SignatureSimple>& sigs, Domain pp, mpz_class sk, mpz_class expected_lb, mpz_class expected_ub) {
	size_t bad = 0;
	for (auto& sig : sigs) {
		mpz_class k_tmp = sk * sig.h + sig.s;
		mpz_class k;
		mpz_mod(k.get_mpz_t(), k_tmp.get_mpz_t(), pp.n.get_mpz_t());
		if (k < expected_lb || expected_ub < k)
			bad += 1;
	}
	double error_rate = 100 * ((double) bad) / ((double) sigs.size());
	gmp_printf("%.2f percent of k are in the range (%Zd,%Zd)\n", 100 - error_rate, expected_lb.get_mpz_t(), expected_ub.get_mpz_t());
}

std::vector<uint8_t> mpz_to_vector(const mpz_class x) {
	size_t size = MPZ_BYTE_SIZE;
	std::vector<uint8_t> v(size);
	mpz_export(&v[0], &size, -1, 1, 0, 0, x.get_mpz_t());
	if (sgn(x) < 0)
		v.back() = 1; // last byte indicates the sign

	return v;
}

mpz_class vector_to_mpz(std::vector<uint8_t>& v) {
	mpz_class x;
	bool is_negative = false;
	if (v.back() == 1) {
		is_negative = true;
		v.back() = 0;
	}
	mpz_import(x.get_mpz_t(), v.size(), -1, 1, 0, 0, &v[0]);

	if (is_negative)
		x = -x;
	return x;
}
/*
 * Offset should start from at least 1. Otherwise we assume that the bit bound of original
 * integers + margin is at most 127, and the integers should not be truncated at all.
 * Instead we calculate how many left shift is needed to adjust the most significant bit
 * to 128-th bit.
 * */
void compute_ofst_uint128(int& ofst, size_t& pad, const int& bit_bound, int margin) {
	ofst = bit_bound + margin - 127;
	if (ofst < 1) {
		pad = (size_t)(-ofst + 1);
	} else {
		pad = 0;
	}
}

/* Extract most significant 128 bits */
void mpz_to_uint128_opt(uint128_t& r, mpz_class& x, const int& ofst) {
	// ofst is the bit position that uint128_t should start at
	if (ofst < 1) {
		mpz_to_uint128(r, x, 0);
		return;
	}
	mpz_class truncated = x >> (ofst - 1);
	mpz_to_uint128(r, truncated, 0);
#if 0
	printf("%lu\n", mpz_sizeinbase(x.get_mpz_t(),2));
	printf("%lu\n", mpz_sizeinbase(truncated.get_mpz_t(),2));
#endif
}

void mpz_to_gs(gs& sc, mpz_class& x) {
	for (auto& vv : sc.v)
		vv = 0;
	mpz_export(sc.v, NULL, -1, 8, 0, 0, x.get_mpz_t());
}

void mpz_to_gs(sc25519& sc, mpz_class& x) {
	sc = sc25519( { .v = { 0, 0, 0, 0 } });
	mpz_export(sc.v, NULL, -1, 8, 0, 0, x.get_mpz_t());
}

void mpz_to_gss(gss& sc, mpz_class& x, const int& offset) {
	gs temp;
	mpz_to_gs(temp, x);
	sc.v[0] = temp.v[offset];
	sc.v[1] = temp.v[offset + 1];
}

void mpz_to_uint64(uint64_t& rhigh, uint64_t& rlow, mpz_class& x, const int& offset) {
	gs temp;
	mpz_to_gs(temp, x);
	uint64_t low = temp.v[offset];
	uint64_t high = temp.v[offset + 1];
	rlow = low;
	rhigh = high;
}

void mpz_to_uint128(uint128_t& r, mpz_class& x, const int& offset) {
	gs temp;
	mpz_to_gs(temp, x);
	uint128_t low = temp.v[offset];
	uint128_t high = temp.v[offset + 1];
	r = low + (high << 64);
}
void mpz_to_uint96(uint96_t& r, mpz_class& x, const int& offset) {
	gs temp;
	mpz_to_gs(temp, x);
	uint32_t low = (uint32_t)(temp.v[offset] >> 32);
	uint64_t high = temp.v[offset + 1];
	r.l = low;
	r.h = high;
}
void uint128_to_uint96(uint96_t& r, const uint128_t& x) {
	uint64_t low = (uint64_t) x;
	r.l = (uint32_t)(low >> 32);
	r.h = (uint64_t)(x >> 64);
}
void uint128v_to_uint96v(std::vector<uint96_t>& rlist, std::vector<uint128_t>& xlist) {
	rlist.reserve(xlist.size());
	uint96_t tmp;
	for (uint128_t &x : xlist) {
		uint128_to_uint96(tmp, x);
		rlist.emplace_back(tmp);
	}
}

void gs_to_mpz(gs& sc, mpz_class& x) {
	mpz_import(x.get_mpz_t(), 4, -1, 8, 0, 0, sc.v);
}

void gs_to_mpz(sc25519& sc, mpz_class& x) {
	mpz_import(x.get_mpz_t(), 4, -1, 8, 0, 0, sc.v);
}

void uint128_to_mpz(uint128_t& r, mpz_class& x) {
	mpz_import(x.get_mpz_t(), 1, 1, sizeof(r), 0, 0, &r);
}

void print_uint128(uint128_t& r) {
	mpz_class tmp;
	uint128_to_mpz(r, tmp);
	gmp_printf("%Zd\n", tmp.get_mpz_t());
}

void gssv_to_uintv(std::vector<uint64_t>& h0list, std::vector<uint64_t>& h1list, std::vector<gss>& hlist) {
	h0list.reserve(hlist.size());
	h1list.reserve(hlist.size());
	for (auto &h : hlist) {
		h0list.emplace_back(h.v[0]);
		h1list.emplace_back(h.v[1]);
	}
}

void uintv_to_gssv(std::vector<uint64_t>& h0list, std::vector<uint64_t>& h1list, std::vector<gss>& hlist) {
	if (h0list.size() != h1list.size()) {
		printf("WARNING: uintv_to_gssv detected unequal vector sizes\n");
	}

	hlist.reserve(h0list.size());
	for (size_t i = 0; i < h0list.size(); i++) {
		gss temp_sc;
		temp_sc.v[0] = h0list[i];
		temp_sc.v[1] = h1list[i];
		hlist.emplace_back(temp_sc);
	}
}

void uint64v_to_uint128v(std::vector<uint64_t>& highlist, std::vector<uint64_t>& lowlist, std::vector<uint128_t>& list) {
	if (highlist.size() != lowlist.size()) {
		printf("WARNING: uint64v_to_uint128v detected unequal vector sizes\n");
	}

	list.reserve(highlist.size());
	for (size_t i = 0; i < highlist.size(); i++) {
		uint128_t high = highlist[i];
		uint128_t low = lowlist[i];
		list.emplace_back((high << 64) + low);
	}
}

void pack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr) {
	hlist.reserve(sigsptr->size());
	for (auto& sig : *sigsptr) {
		gs hscalar;
		mpz_to_gs(hscalar, sig.h);
		hlist.emplace_back(hscalar);
	}
}

void packhalf(std::vector<gss>& hlist, std::vector<SignatureSimple>* sigsptr, const int& offset) {
	hlist.reserve(sigsptr->size());
	for (auto& sig : *sigsptr) {
		gss half_h;
		mpz_to_gss(half_h, sig.h, offset);
		hlist.emplace_back(half_h);
	}
}

void unpack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr) {
	mpz_class x;
	for (auto& h : hlist) {
		gs_to_mpz(h, x);
		sigsptr->emplace_back(SignatureSimple(x, 0));
	}
}

void packedprint(std::vector<gs>& hlist) {
	for (auto& h : hlist) {
		gsprint(h);
	}
}

void gsprint(gs& sc) {
	mpz_class x;
	gs_to_mpz(sc, x);
	gmp_printf("%Zd\n", x.get_mpz_t());
}

void gsprint(sc25519& sc) {
	mpz_class x;
	gs_to_mpz(sc, x);
	gmp_printf("%Zd\n", x.get_mpz_t());
}

/* arithmetic for short group scalar (=gss) */
void gss_add(gss *r, const gss *x, const gss *y) {
	r->v[1] = x->v[1] + y->v[1];
	if ((UINT64_MAX - x->v[0]) < y->v[0]) {
		r->v[1] = r->v[1] + 1;
		r->v[0] = y->v[0] - (UINT64_MAX - x->v[0]) - 1;
	} else {
		r->v[0] = x->v[0] + y->v[0];
	}
}

// assume x > y
void gss_sub(gss *r, const gss *x, const gss *y) {
	r->v[1] = x->v[1] - y->v[1];
	if (x->v[0] < y->v[0]) {
		r->v[1] = r->v[1] - 1;
		r->v[0] = UINT64_MAX - (y->v[0] - x->v[0] - 1);
	} else {
		r->v[0] = x->v[0] - y->v[0];
	}
}

int gss_lt(const gss *a, const gss *b) {
	if (a->v[1] < b->v[1])
		return 1;
	else if (a->v[1] > b->v[1])
		return 0;
	else if (a->v[0] < b->v[0])
		return 1;
	else if (a->v[0] > b->v[0])
		return 0;
	return 0;
}

int gss_lteq(const gss *a, const gss *b) {
	if (a->v[0] == b->v[0] && a->v[1] == b->v[1])
		return 1;
	return gss_lt(a, b);
}
