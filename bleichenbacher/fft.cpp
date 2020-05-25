#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>

#include <boost/mpi.hpp>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <gmpxx.h>

#include "mocksig.h"
#include "fft.h"

double compute_norm(fftw_complex w, uint32_t L) {
	return sqrt(w[0]*w[0] + w[1]*w[1]) / L;
}

void compute_norm_table(double* W_norm, fftw_complex* W, uint64_t C, uint32_t L) {
	for (uint64_t i = 0; i < C; i++)
		W_norm[i] = compute_norm(W[i], L);
}

void find_peak(fftw_complex* W, uint64_t C, uint32_t L, uint64_t keep, double& peak, uint64_t& peak_at,
		std::priority_queue<WCandidate,
		std::vector<WCandidate>,
		std::greater<WCandidate>>& heap) {
	peak_at = 0;
	peak = 0;
	double w_norm;
	for (uint64_t i = 0; i < C; i++) {
		/* compute norm */
		w_norm = compute_norm(W[i], L);

		/* save top largest candidates in heap */
		heap.push(WCandidate(i, w_norm));
		if (i > keep)
			heap.pop();

		//printf("Bn(%u) = %lf = |%+lf %+lf*i|\n", j, W_norm[j], W[j][0], W[j][1]);
		if (w_norm > peak) {
			peak = w_norm;
			peak_at = i;
		}
	}
	return;
}

double compute_noise_avg(fftw_complex* W, uint64_t C, uint32_t L, uint64_t peak_at) {
	double sum = 0;
	for (uint64_t i = 0; i < C; i++) {
		if (i != peak_at) {
			sum += compute_norm(W[i], L);
		}
	}
	return sum / (C - 1);
}

/* Space complexity: |W|= 16 * C bytes since fftw_complex = 16 byte */
void compute_bias(fftw_complex* W, uint64_t C, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs) {
	boost::mpi::environment env;
	fftw_plan p;

	/* init Z_t */
	/* Use Guru interface */
    unsigned long int DOF = 1;

    int rank = 1;
    fftw_iodim64 *dims = (fftw_iodim64*) malloc(rank*sizeof(fftw_iodim64));
    if (dims == NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    dims[0].n = C;
    dims[0].is = DOF;
    dims[0].os = DOF;

    int howmany_rank = 1;
    fftw_iodim64 *howmany_dims = (fftw_iodim64*) malloc(howmany_rank*sizeof(fftw_iodim64));
    if (howmany_dims == NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    howmany_dims[0].n = DOF;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;

    printf("creating the plan\n");
    p = fftw_plan_guru64_dft(rank, dims, howmany_rank, howmany_dims, W, W, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (p == NULL){fprintf(stderr,"plan creation failed\n");exit(1);}
    printf("created the plan\n");

    printf("Initializing FFT table with 0...\n");
#pragma omp parallel shared(W)
    {
#pragma omp for schedule(static)
	for (uint64_t  j = 0; j < C; j++) {
		W[j][0] = 0;
		W[j][1] = 0;
	}
    }
	printf("Preparing FFT table of size %lu...\n", C);
#pragma omp parallel shared(W,sigs)
	{
#pragma omp for schedule(static)
	for (uint64_t  j = 0; j < L; j++) {
		//printf( "%s: thread %d checks j = %lu.\n", env.processor_name().c_str(), omp_get_thread_num (), j);
		mpf_class tmp = (mpf_class) 2 * sigs[j].s / pp.n;
		mpz_class idx = sigs[j].h >> kb;
		W[idx.get_ui()][0] += cos(M_PI * tmp.get_d());
		W[idx.get_ui()][1] += sin(M_PI * tmp.get_d());
	}
	}

	/* doit */
	printf("Computing bias...\n");
	fftw_execute(p);
	printf("Done!\n");

	printf("Cleaning up..\n");
	fftw_destroy_plan(p);
    free(dims);
    free(howmany_dims);

	return;
}


void compute_bias_h(fftw_complex* W, uint64_t C, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSimple>& sigs) {
	boost::mpi::environment env;
	fftw_plan p;

	/* init Z_t */
	/* Use Guru interface */
    unsigned long int DOF = 1;

    int rank = 1;
    fftw_iodim64 *dims = (fftw_iodim64*) malloc(rank*sizeof(fftw_iodim64));
    if (dims == NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    dims[0].n = C;
    dims[0].is = DOF;
    dims[0].os = DOF;

    int howmany_rank = 1;
    fftw_iodim64 *howmany_dims = (fftw_iodim64*) malloc(howmany_rank*sizeof(fftw_iodim64));
    if (howmany_dims == NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    howmany_dims[0].n = DOF;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;

    printf("sizeof fftw complex %ld\n",sizeof(fftw_complex));
    printf("sizeof fftw_iodim64 %ld\n",sizeof(fftw_iodim64));
    printf("creating the plan\n");
    p = fftw_plan_guru64_dft(rank, dims, howmany_rank, howmany_dims, W, W, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (p == NULL){fprintf(stderr,"plan creation failed\n");exit(1);}
    printf("created the plan\n");

    printf("Initializing FFT table with 0...\n");
#pragma omp parallel shared(W)
    {
#pragma omp for schedule(static)
	for (uint64_t j = 0; j < C; j++) {
		W[j][0] = 0;
		W[j][1] = 0;
	}
    }

	printf("Preparing FFT table...\n");
#pragma omp parallel shared(W,sigs)
	{
#pragma omp for schedule(static)
	for (uint64_t j = 0; j < L; j++) {
		mpz_class idx = sigs[j].h >> kb;
		W[idx.get_ui()][0] += 1;
	}
	}

	/* doit */
	printf("Computing bias...\n");
	fftw_execute(p);
	printf("Done!\n");

	printf("Cleaning up..\n");
	fftw_destroy_plan(p);
    free(dims);
    free(howmany_dims);

	return;
}

void compute_bias(fftw_complex* W, uint32_t L, uint32_t kb, Domain pp, std::vector<SignatureSC25519>& sigs) {
	uint32_t j;
	fftw_plan p;

	/* init Z_t */
	p = fftw_plan_dft_1d(L, W, W, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (j = 0; j < L; j++) {
		W[j][0] = 0;
		W[j][1] = 0;
	}
	printf("Preparing FFT table...\n");
	mpf_class tmp;
	for (j = 0; j < L; j++) {
		mpz_class h, s;
		gs_to_mpz(sigs[j].h, h);
		gs_to_mpz(sigs[j].s, s);
		tmp = (mpf_class) 2 * s / pp.n;
		mpz_class idx = h >> kb;
		W[idx.get_ui()][0] += cos(M_PI * tmp.get_d());
		W[idx.get_ui()][1] += sin(M_PI * tmp.get_d());
	}

	/* doit */
	printf("Computing bias...\n");
	fftw_execute(p);
	printf("Done!\n");
	fftw_destroy_plan(p);

	return;
}


/* Space complexity: |W|= 16 * C bytes since fftw_complex = 16 byte */
void compute_bias_mpi(uint64_t C, uint32_t L, uint32_t kb, uint32_t ub, mpz_class sk, mpz_class sk_hi,
			Domain pp, std::vector<SignatureSimple>& sigs, uint64_t fft_outlim, uint32_t batchsize, std::string fname, std::string dir = "") {
	boost::mpi::environment env;
	boost::mpi::communicator world;
	const int master = 0;
	const int myrank = world.rank();
	const int worldsize = world.size();
	if (dir.length())
		dir = dir + "/";

	fftw_plan p;
	fftw_complex* W;
	ptrdiff_t N0 = C;
	ptrdiff_t alloc_local, local_ni, local_i_start, local_no, local_o_start, j;
	fftw_mpi_init();
	alloc_local = fftw_mpi_local_size_1d(N0, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE,
										&local_ni, &local_i_start, &local_no, &local_o_start);
	printf("[%d]/[%d]: C=%lu\n", myrank, worldsize, C);
	printf("[%d]/[%d]: alloc_local=%lu\n", myrank, worldsize, alloc_local);
	printf("[%d]/[%d]: local_ni=%lu\n", myrank, worldsize, local_ni);
	printf("[%d]/[%d]: local_i_start=%lu\n", myrank, worldsize, local_i_start);
	printf("[%d]/[%d]: local_no=%lu\n", myrank, worldsize, local_no);
	printf("[%d]/[%d]: local_o_start=%lu\n", myrank, worldsize, local_o_start);
	W = fftw_alloc_complex(alloc_local);

    printf("[%d]/[%d]: creating the plan\n", myrank, worldsize);
    p = fftw_mpi_plan_dft_1d(N0, W, W, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
    printf("[%d]/[%d]: created the plan\n", myrank, worldsize);
    if (p==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}

    printf("[%d]/[%d]: Initializing FFT table with 0...\n", myrank, worldsize);
#pragma omp parallel shared(W)
    {
#pragma omp for schedule(static)
	for (j = 0; j < local_ni; j++) {
		W[j][0] = 0;
		W[j][1] = 0;
	}
    }

	printf("[%d]/[%d]: Preparing FFT table of size %lu...\n", myrank, worldsize, local_ni);
	uint32_t resume_offset = 0;
	uint32_t sigs_size = 0;
	int batchcount = 0;
	while (resume_offset < L) {
		sigs.reserve(batchsize);
		printf("[%d]/[%d]: Loading batch %d of %s; [%u, %u) \n", myrank, worldsize, batchcount, fname.c_str(), resume_offset, resume_offset+batchsize);
		sigload(sigs, fname, batchsize, false, mpz_class(C) << kb, resume_offset);
		sigs_size = sigs.size();
		printf("[%d]/[%d]: Computing local input...\n", myrank, worldsize);
#pragma omp parallel shared(W,sigs)
		{
#pragma omp for schedule(static)
		for (j = 0; j < sigs_size; j++) {
			//printf( "%s: thread %d checks j = %lu.\n", env.processor_name().c_str(), omp_get_thread_num (), j);
			mpf_class tmp = (mpf_class) 2 * sigs[j].s / pp.n;
			mpz_class idx = sigs[j].h >> kb;
			if ((local_i_start <= idx.get_ui()) &&  (idx.get_ui() < local_i_start + local_ni)) {
				W[idx.get_ui() - local_i_start][0] += cos(M_PI * tmp.get_d());
				W[idx.get_ui() - local_i_start][1] += sin(M_PI * tmp.get_d());
			} else if (idx.get_ui() >= C){
				//printf("WARNING: found index %lu outside FFT table size.\n", idx.get_ui());
				continue;
			}
		}
		}
		printf("[%d]/[%d]: freeing up sigs vector...\n", myrank, worldsize);
		std::vector<SignatureSimple>().swap(sigs);
		malloc_trim(0); // required to clean up the orphaned memory allocated by GMP
		resume_offset += batchsize;
		batchcount++;
	}

	/* doit */
	world.barrier();
	printf("[%d]/[%d]: Computing bias...\n", myrank, worldsize);
	fftw_execute(p);
	printf("[%d]/[%d]: Done!\n", myrank, worldsize);

	/* Find peak and compute noise */
	double local_peak;
	uint64_t local_peak_at;
	std::priority_queue<WCandidate,
			std::vector<WCandidate>,
			std::greater<WCandidate>> heap;
	printf("[%d]/[%d]: Looking for the peak...\n", myrank, worldsize);
	find_peak(W, local_no, L, fft_outlim, local_peak, local_peak_at, heap);
	printf("[%d]/[%d]: local_peak = %lf, local_o_start+local_peak_at = %lu\n", myrank, worldsize,
			local_peak, local_o_start + local_peak_at);
	std::vector<uint64_t> results_peak_at;
	std::vector<double> results_peak;
	world.barrier();
	gather(world, local_o_start+local_peak_at, results_peak_at, master);
	gather(world, local_peak, results_peak, master);

	if (myrank == master) {
		double peak = 0;
		uint64_t peak_at = 0;
		for (size_t i = 0; i < results_peak.size(); i++) {
			if (results_peak[i] > peak) {
				peak = results_peak[i];
				peak_at = results_peak_at[i];
			}
		}

		/* Compute W_m = m*n/C */
		mpz_class w;
		if (kb == 0)
			w = (mpz_class) (peak_at) * pp.n / C;
		else
			w = (mpz_class) (peak_at) * (mpz_class(1) << ub) / C + sk_hi;
		mpf_class w_f = mpf_class(w);
		mpf_class sk_f = mpf_class(sk);
		mpf_class rel_error = abs(w_f - sk_f) / sk_f;
		gmp_printf("	 Estimated secret   w = %Zd \n"
				   "	 Real secret       sk = %Zd \n"
				   "	 Relative error       = %.Ff \n", w.get_mpz_t(), sk.get_mpz_t(), rel_error);

		// Comparison with the actual secret
		std::string skbin = sk.get_str(2);
		std::string wbin = w.get_str(2);
		int count_msb = 0;
		int idx = 0;
		while (1) {
			if (skbin[idx] != wbin[idx])
				break;
			count_msb++;
			idx++;
		}
		std::cout << "	 " << skbin << std::endl;
		std::cout << "	 " << wbin << std::endl;
		printf("	 Recovered %d-MSBs of sk d\n", count_msb);
	}

	/* Saving bias to file */
	std::string bias_fname = "bias-" + std::to_string(myrank) + ".csv";
	printf("[%d]/[%d]: saving top %lu bias candidates to file %s\n",
			myrank, worldsize, fft_outlim, (dir + bias_fname).c_str());
	std::ofstream bias_file((dir + bias_fname).c_str(), std::ofstream::trunc);
	while (!heap.empty()) {
		WCandidate w = heap.top();
		bias_file << w.pos << "," << w.norm << std::endl;
		heap.pop();
	}
	bias_file.close();

	/* Compute noise */
	printf("[%d]/[%d]: Computing average noise...\n", myrank, worldsize);
	double noise = compute_noise_avg(W, local_no, L, local_peak_at);
	printf("	 Average noise             = %lf\n"
		   "	 Estimated noise 1/sqrt(L) = %lf\n", noise, 1 / sqrt(L));

	printf("[%d]/[%d]: Freeing up FFT table...\n", myrank, worldsize);
	fftw_free(W);

	printf("[%d]/[%d]: Cleaning up..\n", myrank, worldsize);
	fftw_destroy_plan(p);
	return;
}
