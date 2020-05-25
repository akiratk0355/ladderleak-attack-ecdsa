#ifndef FFT_H
#define FFT_H

class WCandidate {

public:
	uint64_t pos;
	double norm;

	WCandidate(uint64_t pp = 0, double nn = 0) :
		pos(pp), norm(nn) {
	}
	inline bool operator <(const WCandidate& c) const {
		return norm < c.norm;
	}
	inline bool operator >(const WCandidate& c) const {
		return norm > c.norm;
	}
};

double compute_norm(fftw_complex w, uint32_t L);
void compute_norm_table(double* W_norm, fftw_complex* W, uint64_t C, uint32_t L);
void find_peak(fftw_complex* W, uint64_t C, uint32_t L, uint64_t keep, double& peak, uint64_t& peak_at,
		std::priority_queue<WCandidate,
		std::vector<WCandidate>,
		std::greater<WCandidate>>& heap);
double compute_noise_avg(fftw_complex* W, uint64_t C, uint32_t L, uint64_t peak_at);
void compute_bias(fftw_complex* W, uint64_t C, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs);
void compute_bias(fftw_complex* W, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSC25519>& sigs);
void compute_bias_h(fftw_complex* W, uint64_t C, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs);
void compute_bias_mpi(uint64_t C, uint32_t L, uint32_t kb, uint32_t ub, mpz_class sk, mpz_class sk_hi,
		Domain pp, std::vector<SignatureSimple>& sigs, uint64_t fft_outlim, uint32_t batchsize, std::string fname, std::string dir);
#endif
