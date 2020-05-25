#ifndef MOCKSIG_H
#define MOCKSIG_H

#define MPZ_LOAD_ERROR 0
#define MPZ_BYTE_SIZE 34

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

extern "C" {

#include "qDSA/Curve25519-asm/sc25519.h"
#include "qDSA/Curve25519-asm/scalar.h"

}
#include <vector>

struct Domain {
	uint32_t n_bit;
	mpz_class n;
};

class SignatureLeak {
public:
	mpz_class h;
	mpz_class s;
	mpz_class rr; // leaked nonce

	SignatureLeak(mpz_class h, mpz_class s, mpz_class rr) :
			h(h), s(s), rr(rr) {
	}

	bool operator <(const SignatureLeak& sig) const {
		return (h < sig.h);
	}

	bool operator >(const SignatureLeak& sig) const {
		return (h > sig.h);
	}
};

class SignatureSimple {
public:
	mpz_class h;
	mpz_class s;

	SignatureSimple() :
			h(0), s(0) {
	}
	SignatureSimple(mpz_class h, mpz_class s) :
			h(h), s(s) {
	}

	bool operator <(const SignatureSimple& sig) const {
		return (h < sig.h);
	}

	bool operator >(const SignatureSimple& sig) const {
		return (h > sig.h);
	}
};

class SignatureSC25519 {
public:
	sc25519 h;sc25519 s;

	SignatureSC25519() :
			h( { .v = { 0, 0, 0, 0 } }), s( { .v = { 0, 0, 0, 0 } }) {
	}
	SignatureSC25519(sc25519 h, sc25519 s) :
			h(h), s(s) {
	}

	bool operator <(const SignatureSC25519& sig) const {
		return (1 == sc25519_lt(&h, &sig.h));
	}

	bool operator >(const SignatureSC25519& sig) const {
		return (0 == sc25519_lt(&h, &sig.h));
	}
};

typedef struct {
	uint64_t v[4];
} gs;

typedef struct {
	uint64_t v[2];
} gss;

struct uint96_t {
	uint64_t h;
	uint32_t l;

	// default + parameterized constructor
	uint96_t(uint64_t h = 0, uint32_t l = 0) :
			h(h), l(l) {
	}

	// assignment operator modifies object, therefore non-const
	uint96_t& operator=(const uint96_t& a) {
		h = a.h;
		l = a.l;
		return *this;
	}
	// doesn't modify object. therefore const.
	uint96_t operator+(const uint96_t& a) const {
		return uint96_t(h + a.h + ((l + a.l) < l), l + a.l);
	}
	uint96_t operator-(const uint96_t& a) const {
		return uint96_t(h - a.h - ((l - a.l) > l), l - a.l);
	}
	inline bool operator==(const uint96_t& a) const {
			return ((h == a.h) & (h == a.l));
	}
	inline bool operator<(const uint96_t& a) const {
		return ((h < a.h) | ((h == a.h) & (l < a.l)));
	}
	inline bool operator<=(const uint96_t& a) const {
		return ((*this < a) | (*this == a));
	}
	inline bool operator>(const uint96_t& a) const {
		return ((h > a.h) | ((h == a.h) & (l > a.l)));
	}
	inline bool operator>=(const uint96_t& a) const {
		return ((*this > a) | (*this == a));
	}
	uint96_t operator >>(const unsigned offset) {
		if (offset >= 96) {
			return uint96_t(0, 0);
		} else if (offset == 64) {
			return uint96_t(0, (uint32_t) (h >> 32));
		} else if (offset == 32) {
			return uint96_t(h >> 32, (uint32_t) h);
		} else if (offset == 0) {
			return *this;
		} else if (offset < 32) {
			return uint96_t(h >> offset, (uint32_t) ((h << (64 - offset)) >> (64 - offset)) + (l >> offset));
		} else if (offset < 64) {
			return uint96_t(h >> offset, (uint32_t) (h << (64 - offset) >> 32));
		} else if (offset < 96) {
			return uint96_t(0, (uint32_t) (uint32_t) (h >> (offset - 32)));
		} else {
			return uint96_t(0, 0);
		}
	}
};

namespace mock {
Domain setup(uint32_t n_bit, mpz_class n);
mpz_class keygen(Domain pp, gmp_randclass& rand);
SignatureLeak sign(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand);
SignatureLeak sign_msb_leak(Domain pp, mpz_class d, int M, uint32_t leak, gmp_randclass& rand);
SignatureSimple sign_filter(Domain pp, mpz_class d, int M, uint32_t leak, uint32_t filter, gmp_randclass& rand);
}

unsigned long long rdtsc(void);
void sigprint(SignatureSimple sig);
void sigprint(SignatureSC25519 sig);
void sigvprint(std::vector<SignatureSimple>& sigs, uint32_t idx_start, uint32_t idx_end);
void sigsave(std::vector<SignatureSimple>& sigs, std::string filename, bool str = false);
void sigsave(std::vector<SignatureSC25519>& sigs, std::string filename, bool str = false);
template<class Iterator>
void sigsave_it(Iterator start, Iterator end, std::string filename) {
	FILE *fp = fopen(filename.c_str(), "wb");

	for (auto it = start; it != end; ++it) {
		mpz_out_raw(fp, it->h.get_mpz_t());
		mpz_out_raw(fp, it->s.get_mpz_t());
	}
	fclose(fp);
}
void sigload(std::vector<SignatureSimple>& sigs, std::string filename, uint32_t lim, bool str = false, mpz_class bound = 0, uint32_t resume_offset = 0);
void sigload(std::vector<SignatureSC25519>& sigs, std::string filename);
std::vector<uint8_t> mpz_to_vector(const mpz_class x);
mpz_class vector_to_mpz(std::vector<uint8_t>& v);
void countbias(std::vector<SignatureSimple>& sigs, Domain pp, mpz_class sk, mpz_class expected_lb, mpz_class expected_ub);

/* NEW */
void compute_ofst_uint128(int& ofst, size_t& pad, const int& bit_bound, int margin = 0);
void mpz_to_uint128_opt(uint128_t& r, mpz_class& x, const int& ofst);
template<class Iterator>
void packhalf_it_opt(Iterator start, Iterator end, std::vector<uint64_t>& highlist, std::vector<uint64_t>& lowlist, const int& size, const int& ofst) {
	highlist.reserve(size);
	lowlist.reserve(size);
	for (auto it = start; it != end; ++it) {
		uint128_t tmp;
		mpz_to_uint128_opt(tmp, it->h, ofst);
		highlist.emplace_back((uint64_t)(tmp >> 64));
		lowlist.emplace_back((uint64_t)tmp);
	}
}

void mpz_to_gs(gs& sc, mpz_class& x);
void mpz_to_gs(sc25519& sc, mpz_class& x);
void mpz_to_gss(gss& sc, mpz_class& x, const int& offset);
void mpz_to_uint128(uint128_t& r, mpz_class& x, const int& offset);
void mpz_to_uint64(uint64_t& rhigh, uint64_t& rlow, mpz_class& x, const int& offset);
void mpz_to_uint96(uint96_t& r, mpz_class& x, const int& offset);
void gs_to_mpz(gs& sc, mpz_class& x);
void gs_to_mpz(sc25519& sc, mpz_class& x);
void uint64v_to_uint128v(std::vector<uint64_t>& highlist, std::vector<uint64_t>& lowlist, std::vector<uint128_t>& list);
void uint128_to_mpz(uint128_t& r, mpz_class& x);
void print_uint128(uint128_t& r);
void uint128v_to_uint96v(std::vector<uint96_t>& rlist, std::vector<uint128_t>& xlist);
void gssv_to_uintv(std::vector<uint64_t>& h0list, std::vector<uint64_t>& h1list, std::vector<gss>& hlist);
void uintv_to_gssv(std::vector<uint64_t>& h0list, std::vector<uint64_t>& h1list, std::vector<gss>& hlist);
void pack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr);
void packhalf(std::vector<gss>& hlist, std::vector<SignatureSimple>* sigsptr, const int& offset);
template<class Iterator>
void packhalf_it(Iterator start, Iterator end, std::vector<gss>& hlist, const int& size, const int& offset) {
	hlist.reserve(size);
	for (auto it = start; it != end; ++it) {
		gss half_h;
		mpz_to_gss(half_h, it->h, offset);
		hlist.emplace_back(half_h);
	}
}
template<class Iterator>
void packhalf_it(Iterator start, Iterator end, std::vector<uint128_t>& hlist, const int& size, const int& offset) {
	hlist.reserve(size);
	for (auto it = start; it != end; ++it) {
		uint128_t half_h;
		mpz_to_uint128(half_h, it->h, offset);
		hlist.emplace_back(half_h);
	}
}
template<class Iterator>
void packhalf_it(Iterator start, Iterator end, std::vector<uint64_t>& highlist, std::vector<uint64_t>& lowlist, const int& size, const int& offset) {
	highlist.reserve(size);
	lowlist.reserve(size);
	for (auto it = start; it != end; ++it) {
		uint64_t high, low;
		mpz_to_uint64(high, low, it->h, offset);
		highlist.emplace_back(high);
		lowlist.emplace_back(low);
	}
}
void unpack(std::vector<gs>& hlist, std::vector<SignatureSimple>* sigsptr);
void packedprint(std::vector<gs>& hlist);
void gsprint(gs& sc);
void gsprint(sc25519& sc);
void gss_add(gss *r, const gss *x, const gss *y);
void gss_sub(gss *r, const gss *x, const gss *y);
int gss_lt(const gss *a, const gss *b);
int gss_lteq(const gss *a, const gss *b);
#endif
