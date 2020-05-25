#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <omp.h>
#include <malloc.h>
#include <iostream>
#include <algorithm>
#include <queue>
#include <bitset>
#include <chrono>

#include <gmpxx.h>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/sort/spreadsort/spreadsort.hpp>

#include "mocksig.h"
#include "reduction.h"

namespace mpi = boost::mpi;
namespace spsort = boost::sort::spreadsort;

LRComb::LRComb(mpz_class hh, uint32_t i, uint32_t j) {
	hsum = hh;
	idx_L = i;
	idx_R = j;
}

/* serialize user-defined types */
namespace boost {
namespace serialization {
template<class Archive>
void serialize(Archive &ar, gs &scalar, const unsigned int version) {
	ar & scalar.v;
}

template<class Archive>
void serialize(Archive &ar, gss &scalar, const unsigned int version) {
	ar & scalar.v;
}
}
}

void idxsave(std::vector<Index>& is, const std::string& filename) {
	FILE *fp = fopen(filename.c_str(), "wb");

	for (auto& i : is) {
		fwrite(&i.idx_L1, sizeof(uint32_t), 1, fp);
		fwrite(&i.idx_R1, sizeof(uint32_t), 1, fp);
		fwrite(&i.idx_L2, sizeof(uint32_t), 1, fp);
		fwrite(&i.idx_R2, sizeof(uint32_t), 1, fp);
		fwrite(&i.flip, sizeof(bool), 1, fp);
	}
	fclose(fp);
}

void idxload(std::vector<Index>& is, const std::string& filename) {
	FILE *fp = fopen(filename.c_str(), "rb");
	Index i;
	while (fread(&i.idx_L1, sizeof(uint32_t), 1, fp) && fread(&i.idx_R1, sizeof(uint32_t), 1, fp)
	        && fread(&i.idx_L2, sizeof(uint32_t), 1, fp) && fread(&i.idx_R2, sizeof(uint32_t), 1, fp)
	        && fread(&i.flip, sizeof(bool), 1, fp)) {
		is.emplace_back(i);
	}
	fclose(fp);
}

std::string msb(mpz_class m, uint32_t a, uint32_t bit) {
	std::string bin = m.get_str(2);
	uint32_t padding = bit - bin.length();
	if (padding != 0)
		bin = std::string(padding, '0') + bin;

	return bin.substr(0, a);
}

/* helper method that collects linear combinations of two whose a-MSB is equal to A */
void collect_lc_two(std::vector<LRComb64>& combs, const std::vector<uint128_t>& L, const std::vector<uint128_t>& R,
        const uint32_t& A, const uint32_t& a, const uint32_t& current_threshold_bit, const int& ofst, const size_t& pad,
        const uint32_t& lim) {
	uint32_t i = 0, j = 0;
	bool flag_i = false;
	bool flag_j = false;
	uint32_t bad = 0;
	size_t lshift = 128 - (a + 1) - pad;
	size_t rshift = 128 - (a + 1) - 64;
	uint128_t one_128 = 1;
	uint128_t A_128 = (uint128_t) A;
	uint128_t A0 = A_128 << lshift;
	uint128_t A1_low = (one_128 << lshift) - 1;
	uint128_t A1 = A0 + A1_low;
	uint128_t Amid = A0 + A1_low / 2;
#if 0
	print_uint128(A_128);
	print_uint128(A0);
	print_uint128(Amid);
	print_uint128(A1);
#endif

	uint32_t lsize = L.size();
	uint32_t rsize = R.size();

	uint128_t sum = L[i] + R[j];
	while (combs.size() < lim) {
		if (sum < A0) {
			if (j == rsize - 1)
				break;
			j++;
		} else if (A1 < sum) {
			if (i == lsize - 1)
				break;
			i++;
		} else {
			combs.emplace_back((uint64_t)((sum << pad) >> rshift), i, j);
			/* check if indices are at the end */
			if (i == lsize - 1 && j == rsize - 1)
				break;
			else if (i == lsize - 1) {
				j++;
				sum = L[i] + R[j];
				continue;
			} else if (j == rsize - 1) {
				i++;
				sum = L[i] + R[j];
				continue;
			}

			uint128_t peek_i, peek_j;
			peek_i = L[i + 1] + R[j];
			peek_j = L[i]+ R[j + 1];
			flag_i = (A0 <= peek_i) && (peek_i <= A1);
			flag_j = (A0 <= peek_j) && (peek_j <= A1);

			if (flag_i ^ flag_j) {
				if (flag_i)
					i++;
				else
					j++;
			} else {
				if (flag_i && flag_j)
					bad++;
				uint32_t c = 1;
				if (sum < Amid) {
					while (flag_i && i + c < lsize && combs.size() < lim) {
						peek_i = L[i + c] + R[j];
						flag_i = (A0 <= peek_i) && (peek_i <= A1);
						if (flag_i)
							combs.emplace_back((uint64_t)((peek_i << pad) >> rshift), i + c, j);
						c++;
					}
					j++;
				} else {
					while (flag_j && j + c < rsize && combs.size() < lim) {
						peek_j = L[i] + R[j + c];
						flag_j = (A0 <= peek_j) && (peek_j <= A1);
						if (flag_j)
							combs.emplace_back((uint64_t)((peek_j << pad) >> rshift), i, j + c);
						c++;
					}
					i++;
				}
			}
		} // else
		sum = L[i] + R[j];
	} // while (combs.size() < lim)
	//std::cout << "bad: " << bad << std::endl;
	//std::cout << i << "," << j << std::endl;
}


void exhaustive_four_sum(std::vector<SignatureSimple>& sigs, const uint32_t threshold_bit, const uint32_t keep_max, const int log_prec) {
	mpi::environment env;
	mpi::communicator world;
	const int master = 0;
	const int myrank = world.rank();
	const int worldsize = world.size();

	uint32_t S, q1, q2, q3;
	mpz_class threshold_mpz = mpz_class(1) << threshold_bit;
	std::vector<SignatureSimple> result;
	result.reserve(keep_max);
	uint32_t local_keep_max = keep_max / omp_get_max_threads();
	uint32_t n_threads = omp_get_max_threads();

	/* Split sigs into L1 || R1 || L2 || R2 */
	S = sigs.size();
	q1 = S / 4;
	q2 = S / 2;
	q3 = S * 3 / 4;
	gmp_printf("[%d]/[%d]: exhaustively looking for %u sums with threshold = %Zd\n",
			myrank, worldsize, keep_max, threshold_mpz.get_mpz_t());
	printf("[%d]/[%d]: using %u threads\n", myrank, worldsize, n_threads);
#pragma omp parallel shared(sigs,result)
	{
	bool break_flag = false;
	mpz_class h, s, hi, hj, hk, hl, si, sj, sk, sl;
	std::vector<SignatureSimple> local_result;
	local_result.reserve(local_keep_max);

#pragma omp for schedule(static)
	for (uint32_t i = 0; i < q1; i++) {
		if (break_flag) continue;
		hi = sigs[i].h;
		si = sigs[i].s;
		for (uint32_t j = q1; j < q2; j++) {
			if (break_flag) continue;
			hj = sigs[j].h;
			sj = sigs[j].s;
			for (uint32_t k = q2; k < q3; k++) {
				if (break_flag) continue;
				hk = sigs[k].h;
				sk = sigs[k].s;
				for (uint32_t l = q3; l < S; l++) {
					if (break_flag) continue;
					hl = sigs[l].h;
					sl = sigs[l].s;
					h = hi + hj - hk - hl;
					if (h >= 0 && h < threshold_mpz) {
						local_result.emplace_back(SignatureSimple(h, si + sj - sk - sl));
					} else if (h < 0 && h > -threshold_mpz) {
						local_result.emplace_back(SignatureSimple(-h, sk + sl - si - sj));
					}
					if (local_result.size() >= local_keep_max) {
						break_flag = true;
					}
				}
			}
		}
	}
#pragma omp critical
	{
		result.insert(result.end(), local_result.begin(), local_result.end());
	}
	} // pragma omp parallel

	printf("[%d]/[%d]: done!\n", myrank, worldsize);
	printf("[%d]/[%d]: freeing up memory..\n", myrank, worldsize);
	std::vector<SignatureSimple>().swap(sigs);
	malloc_trim(0);
	if (result.size() < keep_max) {
		printf("[%d]/[%d]: WARNING: found only %lu < %u collisions \n", myrank, worldsize, result.size(), keep_max);
	} else {
		printf("[%d]/[%d]: found %lu collisions\n", myrank, worldsize, result.size());
	}
	printf("[%d]/[%d]: copying result..\n", myrank, worldsize);
	sigs.reserve(result.size());
	std::copy(result.begin(), result.end(), std::back_inserter(sigs));
}

void schroeppel_shamir_mpi(std::vector<SignatureSimple>& sigs, mpz_class& bound, const uint32_t n_bit, uint32_t l, uint32_t b,
        const uint32_t filter, uint32_t a, const int log_prec,
        const std::vector<uint32_t>& b_info, const std::vector<uint32_t>& l_info, const size_t iota = 1,
        const std::string out_prefix = "", const bool out_index = false, std::string dir = "", const bool istest = true) {
	mpi::environment env;
	mpi::communicator world;
	const int master = 0;
	const int myrank = world.rank();
	const int worldsize = world.size();

	uint32_t threshold_bit = n_bit - filter; // threshold bit for the next samples
	uint32_t c_threshold_bit = n_bit - filter; // threshold bit for the current samples
	mpz_class threshold_mpz;
	uint32_t keep_min;
	uint32_t keep_max;
	uint32_t l_next;
	uint32_t S, q1, q2, q3;
	const std::string fname_prefix = "ss-mpi";
	if ((b_info.size() != iota) | (l_info.size() != iota + 1)) {
		printf("ERROR: invalid offset_info/b_info/l_info\n");
		return;
	}
	if (dir.length()) {
		dir = dir + "/";
	}

	for (size_t round = 0; round < iota; round++) {
		if (myrank == master)
			printf("-------------------------------------------------- Round %lu begins --------------------------------------------------\n", round);
		/* Compute how many samples to be kept */
		if (l_info.empty()) {
			keep_min = 1 << l;
			keep_max = keep_min * 2; // tentative
			a = l - 2;
		} else {
			l_next = l_info[round + 1];
			keep_min = 1 << (l_next);
			if (round == 0) {
				keep_max = keep_min * 1.01; // tentative
			} else if (round == 1) {
				keep_max = keep_min * 1.01;
			} else {
				keep_max = keep_min * 1.01;
			}

			a = l_info[round] - 2;
		}

		/* Compute the next threshold value */
		c_threshold_bit = threshold_bit;
		if (b_info.empty())
			threshold_bit -= b;
		else
			threshold_bit -= b_info[round];

		threshold_mpz = mpz_class(1) << threshold_bit;
		// Note that b has to be larger than a
		// Low word has to be strictly smaller than this bound
		uint64_t threshold_64 = 1ULL << (64 - (b_info[round] - a));
		int ofst;
		size_t pad;
		compute_ofst_uint128(ofst, pad, c_threshold_bit, 1);

		if (myrank == master) {
			printf("master: b=%u, a=%u\n", b_info[round], a);
			printf("master: ofst = %d, pad = %lu\n", ofst, pad);
			printf("master: c_threshold_bit = %u\n", c_threshold_bit);
			printf("master: threshold_bit = %u\n", threshold_bit);
			printf("master: threshold_64 = %lu\n", threshold_64);
			gmp_printf("master: threshold_mpz = %Zd\n", threshold_mpz.get_mpz_t());
		}

		// File format: ss-mpi_round-i_{L1,R1,L2,R2}.bin
		std::string L1_fname, R1_fname, L2_fname, R2_fname;
		L1_fname = fname_prefix + "_round-" + std::to_string(round) + "_L1.bin";
		R1_fname = fname_prefix + "_round-" + std::to_string(round) + "_R1.bin";
		L2_fname = fname_prefix + "_round-" + std::to_string(round) + "_L2.bin";
		R2_fname = fname_prefix + "_round-" + std::to_string(round) + "_R2.bin";

		std::vector<uint128_t> L1_h, R1_h, L2_h, R2_h;
        std::vector<uint64_t> L1_hhigh, R1_hhigh, L2_hhigh, R2_hhigh, L1_hlow, R1_hlow, L2_hlow, R2_hlow;
		if (myrank == master) {
			/* Split sigs into L1 || R1 || L2 || R2 */
			S = sigs.size();
			q1 = S / 4;
			q2 = S / 2;
			q3 = S * 3 / 4;
#if 1
			std::cout << "master: sorting L1 in descending order..." << std::endl;
			std::sort(sigs.begin(), sigs.begin() + q1, std::greater<SignatureSimple>());
			std::cout << "master: sorting R1 in ascending order..." << std::endl;
			std::sort(sigs.begin() + q1, sigs.begin() + q2);
			std::cout << "master: sorting L2 in descending order..." << std::endl;
			std::sort(sigs.begin() + q2, sigs.begin() + q3, std::greater<SignatureSimple>());
			std::cout << "master: sorting R2 in ascending order..." << std::endl;
			std::sort(sigs.begin() + q3, sigs.end());
			std::cout << "master: sorting done" << std::endl;
#endif
			//std::cin.ignore();

			/* Save split sigs */
			if (!istest) {
				std::cout << "master: saving split lists..." << std::endl;
				sigsave_it(sigs.begin(), sigs.begin() + q1, dir + L1_fname);
				sigsave_it(sigs.begin() + q1, sigs.begin() + q2, dir + R1_fname);
				sigsave_it(sigs.begin() + q2, sigs.begin() + q3, dir  + L2_fname);
				sigsave_it(sigs.begin() + q3, sigs.end(), dir + R2_fname);
				//std::cin.ignore();
			}

			/* Conversion: mpz_t -> uint64_t */
			std::cout << "master: converting mpz_t to uint64_t..." << std::endl;
			packhalf_it_opt(sigs.begin(), sigs.begin() + q1, L1_hhigh, L1_hlow, q1, ofst);
			packhalf_it_opt(sigs.begin() + q1, sigs.begin() + q2, R1_hhigh, R1_hlow, q2 - q1, ofst);
			packhalf_it_opt(sigs.begin() + q2, sigs.begin() + q3, L2_hhigh, L2_hlow, q3 - q2, ofst);
			packhalf_it_opt(sigs.begin() + q3, sigs.end(), R2_hhigh, R2_hlow, S - q3, ofst);
			//std::cin.ignore();
#if 1
			std::cout << "master: cleaning up the original GMP vector..." << std::endl;
			std::vector<SignatureSimple>().swap(sigs);
			malloc_trim(0); // required to clean up the orphaned memory allocated by GMP
#endif
		}
		world.barrier();
		//std::cin.ignore();

		/* Broadcast sigs */
		//TODO: how to force broadcast to preallocate the exact memory for vector?
		printf("[%d]/[%d]: distributing signature data..\n", myrank, worldsize);
		broadcast(world, L1_hhigh, master);;
		broadcast(world, L1_hlow, master);;
		broadcast(world, R1_hhigh, master);
		broadcast(world, R1_hlow, master);
		broadcast(world, L2_hhigh, master);
		broadcast(world, L2_hlow, master);
		broadcast(world, R2_hhigh, master);
		broadcast(world, R2_hlow, master);
		printf("[%d]/[%d]: successfully distributed signature data!\n", myrank, worldsize);

		/* Conversion: uint64_t -> uint128_t */
		printf("[%d]/[%d]: converting uint64 -> uint128...\n", myrank, worldsize);
		uint64v_to_uint128v(L1_hhigh, L1_hlow, L1_h);
		uint64v_to_uint128v(R1_hhigh, R1_hlow, R1_h);
		uint64v_to_uint128v(L2_hhigh, L2_hlow, L2_h);
		uint64v_to_uint128v(R2_hhigh, R2_hlow, R2_h);
		printf("[%d]/[%d]: filled L1 with %lu signatures \n", myrank, worldsize, L1_h.size());
		printf("[%d]/[%d]: filled R1 with %lu signatures \n", myrank, worldsize, R1_h.size());
		printf("[%d]/[%d]: filled L2 with %lu signatures \n", myrank, worldsize, L2_h.size());
		printf("[%d]/[%d]: filled R2 with %lu signatures \n", myrank, worldsize, R2_h.size());

		/* Cleanup */
		printf("[%d]/[%d]: cleaning up uint64_t vectors...\n", myrank, worldsize);
        std::vector<uint64_t>().swap(L1_hhigh);
        std::vector<uint64_t>().swap(L1_hlow);
        std::vector<uint64_t>().swap(R1_hhigh);
        std::vector<uint64_t>().swap(R1_hlow);
        std::vector<uint64_t>().swap(L2_hhigh);
        std::vector<uint64_t>().swap(L2_hlow);
        std::vector<uint64_t>().swap(R2_hhigh);
        std::vector<uint64_t>().swap(R2_hlow);
		//std::cin.ignore();

		/* Job scheduling */
		const uint32_t A_lim = 1U << a;
		const uint32_t range = (A_lim  + worldsize - 1) / worldsize;

		if (myrank == master) {
			printf("master: trying to get 2^%u values less than 2^%u\n"
					"[keep_min, keep_max] = [%u, %u]\n", l_next, threshold_bit, keep_min, keep_max);
			printf("master: looking for collisions on top %u-bits, range is %u\n", a, range);
		}

		/* Start measuring time */
		printf("[%d]/[%d]: reduction started\n", myrank, worldsize);
		auto start = std::chrono::high_resolution_clock::now();

		float percent = 0;
		std::vector<Index> subresult;
		const uint32_t keep_min_sub = keep_min / worldsize;
		const uint32_t keep_max_sub = keep_max / worldsize;
		// Don't change this! The scaling factor below affects total memory a lot with many threads
		const uint32_t keep_max_combs = (1 << a) * 1.05;
		/*
		 * Decide how many collisions to be kept per round based on the tradeoff formula.
		 * 		m' = a' + 2 =  3*a + v - n
		 * Hence for v = 0 (i.e., one iteration of HGJ) m' = 3*a - n collisions are expected.
		 * If n = 3*a - 2 then m' = 4 collisions are expected per each iteration.
		 * If this is repeated for each [0, 2^a) then we get 2^m' = 2^a * 4 samples again.
		 * In practice we get even more and that's why the bound below is multiplied a bit.
		 * */
		uint64_t keep_max_combs_A = (1UL << (3*a - b_info[round])) * 16;
		uint32_t n_threads = omp_get_max_threads();
		if (keep_max_combs_A > keep_max_sub / n_threads)
			keep_max_combs_A = keep_max_sub / n_threads;
		printf("[%d]/[%d]: using %u threads\n", myrank, worldsize, n_threads);
		if (myrank == master)
			printf("[%d]/[%d]: keep_max_combs_A = %lu \n", myrank, worldsize, keep_max_combs_A);
		const uint32_t log_mod = keep_min_sub / log_prec;
		size_t log_th = log_mod;
		const float percent_delta = 100.0 / log_prec;
		subresult.reserve(keep_max_sub * 1.1); // we need some margin because multiple threads may push partial results at a time
#pragma omp parallel shared(subresult,L1_h,R1_h,L2_h,R2_h,threshold_64,c_threshold_bit,threshold_bit,ofst,pad,log_th,percent)
		{
			std::vector<LRComb64> combs1;
			std::vector<LRComb64> combs2;
			std::vector<Index> combs_A;
			combs1.reserve(keep_max_combs);
			combs2.reserve(keep_max_combs);
			combs_A.reserve(keep_max_combs_A);
			uint64_t before = 0, after = 0, ubefore = 0, uafter = 0;
#pragma omp for schedule(dynamic)
			for (uint32_t k = 0; k < range; k++) {
				for (int rev = 0; rev < 2; rev++) {
#if 0
					ubefore = rdtsc();
#endif
					uint32_t A = worldsize * k + myrank;
					if (A >= A_lim || subresult.size() >= keep_max_sub)
						continue;

					if (rev == 0) {
						A = A_lim + A;
					} else {
						A = A_lim - 1 - A;
					}
					//printf("[%d]/[%d]: Finding collisions with A = %u\n", myrank, worldsize, A);

					/* Look for collisions on top a-bits */
#if 0
					before = rdtsc();
					collect_lc_two(combs1, L1_h, R1_h, A, a, c_threshold_bit, ofst, pad, keep_max_combs);
					after = rdtsc();
					printf("%lu cycles @collect_lc_two-1:\n", after-before);
#else
					collect_lc_two(combs1, L1_h, R1_h, A, a, c_threshold_bit, ofst, pad, keep_max_combs);
#endif
					if (combs1.size() == 0)
						continue;
#if 0
					before = rdtsc();
					//std::sort(combs1.begin(), combs1.end());
					spsort::integer_sort(combs1.begin(), combs1.end());
					after = rdtsc();
					printf("%lu cycles @sort-1\n", after-before);
#else
					spsort::integer_sort(combs1.begin(), combs1.end());
#endif

#if 0
					before = rdtsc();
					collect_lc_two(combs2, L2_h, R2_h, A, a, c_threshold_bit, ofst, pad, keep_max_combs);
					after = rdtsc();
					printf("%lu cycles @collect_lc_two-2:\n", after-before);
#else
					collect_lc_two(combs2, L2_h, R2_h, A, a, c_threshold_bit, ofst, pad, keep_max_combs);
#endif
					if (combs2.size() == 0)
						continue;
#if 0
					before = rdtsc();
					//std::sort(combs2.begin(), combs2.end());
					spsort::integer_sort(combs2.begin(), combs2.end());
					after = rdtsc();
					printf("%lu cycles @sort-2\n", after-before);
#else
					spsort::integer_sort(combs2.begin(), combs2.end());
#endif
#if 0
					printf("[%d]/[%d]: A = %u: Found (%lu, %lu) partial collisions\n",
							myrank, worldsize, A, combs1.size(), combs2.size());
#endif
#if 0
					printf( "[%d]/[%d]@%s: thread %d checks A = %u.\n",
							myrank, worldsize, env.processor_name().c_str(), omp_get_thread_num (), A);
#endif
					uint32_t i = 0;
					uint32_t j = 0;
					bool flip;

#if 0
					before = rdtsc();
#endif
					while (i < combs1.size() && j < combs2.size()) {
						uint64_t lsum = combs1[i].hsum;
						uint64_t rsum = combs2[j].hsum;
						flip = lsum < rsum;
						uint64_t h_diff;
						if (flip)
							h_diff = rsum - lsum;
						else
							h_diff = lsum - rsum;

						if (h_diff < threshold_64 && combs_A.size() < keep_max_combs_A) {
							combs_A.emplace_back(combs1[i].idx_L, combs1[i].idx_R, combs2[j].idx_L, combs2[j].idx_R,
							        flip);
						}
						if (flip)
							i++;
						else
							j++;
					}
#if 0
					after = rdtsc();
					printf("%lu cycles @while-loop\n", after-before);
#endif
#if 0
					printf("[%d]/[%d]: A = %u: Found %lu total collisions\n", myrank, worldsize, A, combs_A.size());
#endif
#pragma omp critical
					{
						subresult.insert(subresult.end(), combs_A.begin(), combs_A.end());
						if (subresult.size() >= log_th) {
							uint32_t incr = (subresult.size() - log_th) / log_mod + 1;
							percent += percent_delta * incr;
							log_th += log_mod * incr;
							printf("[%d]/[%d]: %.2f %% done\n", myrank, worldsize, percent);
						}
					}
					combs1.resize(0);
					combs2.resize(0);
					combs_A.resize(0);
					combs1.reserve(keep_max_combs);
					combs2.reserve(keep_max_combs);
					combs_A.reserve(keep_max_combs_A);
#if 0
					uafter = rdtsc();
					printf("%lu cycles @one-loop\n", uafter-ubefore);
#endif
				} // for rev
			} // for k
			printf( "[%d]/[%d]@%s: thread %d is freeing up memory...\n",
					myrank, worldsize, env.processor_name().c_str(), omp_get_thread_num ());
			std::vector<LRComb64>().swap(combs1);
			std::vector<LRComb64>().swap(combs2);
			std::vector<Index>().swap(combs_A);
		} // #pragma omp parallel shared
		if (subresult.size() < keep_min_sub) {
			printf("[%d]/[%d]: WARNING: got %lu subresult < keep_min_sub = %u; failed to get sufficiently many collisions!\n",
					myrank, worldsize, subresult.size(), keep_min_sub);
		}
		if (subresult.size() > keep_max_sub) {
			printf("[%d]/[%d]: got %lu subresult > keep_max_sub = %u; truncating.\n",
					myrank, worldsize,	subresult.size(), keep_max_sub);
			subresult.resize(keep_max_sub);
		}
		printf("[%d]/[%d]: got %lu subresult\n", myrank, worldsize, subresult.size());


		/* End measuring time */
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		printf("[%d]/[%d]: Elapsed time for round %lu = %1.f seconds\n", myrank, worldsize, round, elapsed.count());

		/*
		 * Save interim report (only indices).
		 * file format: index_round-i_rank-j.bin
		 * Turn this on when launching the large attack.
		 * In case of MPI communication failure in the following steps we can still manually
		 * recover the linear combinations from index files.
	 	 */
		if (out_index == true) {
			std::string outidx = "index_round-" + std::to_string(round) + "_rank-" \
					+ std::to_string(myrank) + ".bin";
			printf("[%d]/[%d]: Saving Index subresult to %s...; leading to h < 2^%u \n", myrank, worldsize, outidx.c_str(), threshold_bit);
			idxsave(subresult,  dir + outidx);
			/* debug logging */
			/*
			 std::vector<Index>().swap(subresult);
			 idxload(subresult, outidx);
			 printf("(%u, %u, %u, %u), flip=%d \n", subresult[0].idx_L1, subresult[0].idx_R1, subresult[0].idx_L2, subresult[0].idx_R2, subresult[0].flip);
			 */
		}

		/* Cleanup */
		printf("[%d]/[%d]: cleaning up the original lists L1_h, R1_h, L2_h, R2_h... \n", myrank, worldsize);
		std::vector<uint128_t>().swap(L1_h);
		std::vector<uint128_t>().swap(R1_h);
		std::vector<uint128_t>().swap(L2_h);
		std::vector<uint128_t>().swap(R2_h);
		printf("[%d]/[%d]: waiting for other nodes to finish... \n", myrank, worldsize);
		world.barrier();

		/* Prepare subres_sizes for gatherv */
		printf("[%d]/[%d]: preparing subres_sizes for gatherv... \n", myrank, worldsize);
		std::vector<int> subres_sizes;
		gather(world, static_cast<int>(subresult.size()), subres_sizes, master);
		// casting is required since MPI gatherv does not allow size specification with unsigned int...

		/* Convert Index to vector */
		printf("[%d]/[%d]: converting Index struct to separate vectors... \n", myrank, worldsize);
		std::vector<uint32_t> subresult_L1, subresult_R1, subresult_L2, subresult_R2;
		std::vector<uint8_t> subresult_fl;
		subresult_L1.reserve(subresult.size());
		subresult_R1.reserve(subresult.size());
		subresult_L2.reserve(subresult.size());
		subresult_R2.reserve(subresult.size());
		subresult_fl.reserve(subresult.size());
		for (auto& I : subresult) {
			subresult_L1.emplace_back(I.idx_L1);
			subresult_R1.emplace_back(I.idx_R1);
			subresult_L2.emplace_back(I.idx_L2);
			subresult_R2.emplace_back(I.idx_R2);
			subresult_fl.emplace_back(static_cast<uint8_t>(I.flip));
		}
		printf("[%d]/[%d]: cleaning up subresult... \n", myrank, worldsize);
		std::vector<Index>().swap(subresult);

		/* gatherv subresult into result*/
		std::vector<uint32_t> result_L1, result_R1, result_L2, result_R2;
		std::vector<uint8_t> result_fl;
		if (myrank == master) {
			printf("master: start gathering results...\n");
			uint32_t res_size = 0;

			for (auto& subsize : subres_sizes)
				res_size += subsize;

			result_L1.resize(res_size);
			result_R1.resize(res_size);
			result_L2.resize(res_size);
			result_R2.resize(res_size);
			result_fl.resize(res_size);

			printf("master: gathering subresults...\n");
			gatherv(world, subresult_L1, &result_L1[0], subres_sizes, master);
			printf("master: gathered %lu elements in result_L1\n", result_L1.size());
			gatherv(world, subresult_R1, &result_R1[0], subres_sizes, master);
			printf("master: gathered %lu elements in result_R1\n", result_R1.size());
			gatherv(world, subresult_L2, &result_L2[0], subres_sizes, master);
			printf("master: gathered %lu elements in result_L2\n", result_L2.size());
			gatherv(world, subresult_R2, &result_R2[0], subres_sizes, master);
			printf("master: gathered %lu elements in result_R2\n", result_R2.size());
			gatherv(world, subresult_fl, &result_fl[0], subres_sizes, master);
			printf("master: gathered %lu elements in result_fl\n", result_fl.size());
			printf("master: gatherv completed successfully!\n");
		} else {
			printf("[%d]/[%d]: sending subresult... \n", myrank, worldsize);
			gatherv(world, subresult_L1, master);
			gatherv(world, subresult_R1, master);
			gatherv(world, subresult_L2, master);
			gatherv(world, subresult_R2, master);
			gatherv(world, subresult_fl, master);
		}
		printf("[%d]/[%d]: cleaning up subresult vectors... \n", myrank, worldsize);
		std::vector<int>().swap(subres_sizes);

		std::vector<uint32_t>().swap(subresult_L1);
		std::vector<uint32_t>().swap(subresult_R1);
		std::vector<uint32_t>().swap(subresult_L2);
		std::vector<uint32_t>().swap(subresult_R2);
		std::vector<uint8_t>().swap(subresult_fl);

		if (istest) {
			printf("master: reduction ended; exiting test mode.\n");
			return;
		}

		/* Fukugen */
		std::vector<Index> result;
		if (myrank == master) {
			printf("master: reconstructing Index struct...\n");
			result.reserve(result_L1.size());
			for (size_t i = 0; i < result_L1.size(); i++) {
				result.emplace_back(result_L1[i], result_R1[i], result_L2[i], result_R2[i],
				        static_cast<bool>(result_fl[i]));
			}
			printf("master: gathered %lu elements in result\n", result.size());
			printf("master: cleaning up separate result vectors...\n");
			std::vector<uint32_t>().swap(result_L1);
			std::vector<uint32_t>().swap(result_R1);
			std::vector<uint32_t>().swap(result_L2);
			std::vector<uint32_t>().swap(result_R2);
			std::vector<uint8_t>().swap(result_fl);

			std::vector<SignatureSimple> L1, R1, L2, R2;
			L1.reserve(q1);
			R1.reserve(q2 - q1);
			L2.reserve(q3 - q2);
			R2.reserve(S - q3);
			printf("master: reloading the original samples from %s\n", (dir + L1_fname).c_str());
			sigload(L1, dir + L1_fname, q1);
			printf("master: reloading the original samples from %s\n", (dir + R1_fname).c_str());
			sigload(R1, dir + R1_fname, q2 - q1);
			printf("master: reloading the original samples from %s\n", (dir + L2_fname).c_str());
			sigload(L2, dir + L2_fname, q3 - q2);
			printf("master: reloading the original samples from %s\n", (dir + R2_fname).c_str());
			sigload(R2, dir + R2_fname, S - q3);
			sigs.reserve(result.size());
			printf("master: computing linear combinations from indices...\n");
			restore_from_idx(sigs, result, L1, R1, L2, R2, threshold_mpz);

			printf("master: got %lu result \n", sigs.size());
			if (sigs.size() < keep_min)
				puts("WARNING: failed to get expected amount of reduced values!");

			/* Save the result after each round */
			if (out_prefix.length()) {
				// file format: prefix_round-i.bin
				std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
				printf("master: saving signatures of h < 2^%u to %s... \n", threshold_bit, (dir + outsig).c_str());
				sigsave(sigs, dir + outsig);
			}

			/* cleanups */
			printf("master: cleaning up result, L1, R1, L2, R2...\n");
			std::vector<Index>().swap(result);
			std::vector<SignatureSimple>().swap(L1);
			std::vector<SignatureSimple>().swap(R1);
			std::vector<SignatureSimple>().swap(L2);
			std::vector<SignatureSimple>().swap(R2);
			malloc_trim(0);
			printf("-------------------------------------------------- Round %lu finished --------------------------------------------------\n", round);
		}
	} // for round

	if (myrank == master) {
		bound = threshold_mpz;
		gmp_printf("master: reduction ended after %u rounds, h < %Zd\n", iota, bound.get_mpz_t());
	}
}

void restore_from_idx(std::vector<SignatureSimple>& sigs, const std::vector<Index>& idxlist,
        const std::vector<SignatureSimple>& L1, const std::vector<SignatureSimple>& R1,
        const std::vector<SignatureSimple>& L2, const std::vector<SignatureSimple>& R2,
        const mpz_class& threshold_mpz) {
	for (auto& r : idxlist) {
		mpz_class h = L1[r.idx_L1].h + R1[r.idx_R1].h - L2[r.idx_L2].h - R2[r.idx_R2].h;
		mpz_class s = L1[r.idx_L1].s + R1[r.idx_R1].s - L2[r.idx_L2].s - R2[r.idx_R2].s;
		if (abs(h) >= threshold_mpz) {
			printf("WARNING: found h of %lu-bit at (%u, %u, %u, %u), skipping\n",
			mpz_sizeinbase(h.get_mpz_t(), 2), r.idx_L1, r.idx_R1, r.idx_L2, r.idx_R2);
			continue;
		}
		if (h < 0) {
			if (!r.flip) {
				printf("WARNING: flip info was incorrect; recovering\n");
				mpz_class leftsum = L1[r.idx_L1].h + R1[r.idx_R1].h;
				mpz_class rightsum = L2[r.idx_L2].h + R2[r.idx_R2].h;
				gmp_printf("lsum = %Zd\n", leftsum.get_mpz_t());
				gmp_printf("rsum = %Zd\n", rightsum.get_mpz_t());
			}
			sigs.emplace_back(-h, -s);
		} else {
			if (r.flip) {
				printf("WARNING: flip info was incorrect; recovering\n");
				mpz_class leftsum = L1[r.idx_L1].h + R1[r.idx_R1].h;
				mpz_class rightsum = L2[r.idx_L2].h + R2[r.idx_R2].h;
				gmp_printf("lsum = %Zd\n", leftsum.get_mpz_t());
				gmp_printf("rsum = %Zd\n", rightsum.get_mpz_t());
			}
			sigs.emplace_back(h, s);
		}
	}
}

/* Old SS code without parallelization. */
void schroeppel_shamir(std::vector<SignatureSimple> * sigsptr, const uint32_t n_bit, const uint32_t l, const uint32_t b,
        const uint32_t filter, const size_t iota = 1) {
	uint32_t threshold_bit = n_bit - filter;
	uint32_t S = sigsptr->size();
	uint32_t keep = 1 << l;
	std::vector<SignatureSimple>* L1;
	std::vector<SignatureSimple>* R1;
	std::vector<SignatureSimple>* L2;
	std::vector<SignatureSimple>* R2;

	printf("schroeppel_shamir: got %u sigs\n", S);

	for (size_t round = 0; round < iota; round++) {
		S = sigsptr->size();

		/* split sigs into L1 || R1 || L2 || R2 */
		uint32_t q1 = S / 4;
		uint32_t q2 = S / 2;
		uint32_t q3 = S * 3 / 4;

		std::cout << "Splitting into 4 lists..." << std::endl;
		L1 = new std::vector<SignatureSimple>(sigsptr->begin(), sigsptr->begin() + q1);
		R1 = new std::vector<SignatureSimple>(sigsptr->begin() + q1, sigsptr->begin() + q2);
		L2 = new std::vector<SignatureSimple>(sigsptr->begin() + q2, sigsptr->begin() + q3);
		R2 = new std::vector<SignatureSimple>(sigsptr->begin() + q3, sigsptr->end());
		delete sigsptr;
		sigsptr = new std::vector<SignatureSimple>();

		threshold_bit -= b;

		mpz_class threshold;
		mpz_ui_pow_ui(threshold.get_mpz_t(), 2, threshold_bit);

		/* Init heaps */
		std::cout << "Sorting R1 and R2..." << std::endl;
		std::sort(R1->begin(), R1->end());
		std::sort(R2->begin(), R2->end());

		std::cout << "Initializing heaps..." << std::endl;
		std::priority_queue<LRComb, std::vector<LRComb>> heap1;
		std::priority_queue<LRComb, std::vector<LRComb>> heap2;
		for (uint32_t i = 0; i != L1->size(); i++) {
			heap1.push(LRComb((*L1)[i].h + (*R1)[0].h, i, 0));
			heap2.push(LRComb((*L2)[i].h + (*R2)[0].h, i, 0));
		}

		/* Start looking for collisions */
		std::cout << "Trying to get 2^" << l << " values less than 2^" << threshold_bit << std::endl;
		uint32_t counter = 0;
		uint32_t num_col = 0;
		mpz_class hsum1, hsum2, h_diff, s_diff;
		uint32_t i1, j1, i2, j2;
		while (num_col < keep && heap1.empty() == false && heap2.empty() == false) {
			counter++;
			hsum1 = heap1.top().hsum;
			i1 = heap1.top().idx_L;
			j1 = heap1.top().idx_R;

			hsum2 = heap2.top().hsum;
			i2 = heap2.top().idx_L;
			j2 = heap2.top().idx_R;

			h_diff = hsum1 - hsum2;
			s_diff = (*L1)[i1].s + (*R1)[j1].s - (*L2)[i2].s - (*R2)[j2].s;

			if (h_diff > 0) {
				if (j2 < S / 4 - 1) {
					heap2.pop();
					heap2.push(LRComb((*L2)[i2].h + (*R2)[j2 + 1].h, i2, j2 + 1));
				} else {
					heap2.pop();
				}
			} else {
				h_diff = -h_diff;
				s_diff = -s_diff;
				if (j1 < S / 4 - 1) {
					heap1.pop();
					heap1.push(LRComb((*L1)[i1].h + (*R1)[j1 + 1].h, i1, j1 + 1));
				} else {
					heap1.pop();
				}
			}

			if (h_diff < threshold) {
				sigsptr->emplace_back(SignatureSimple(h_diff, s_diff));
				num_col++;
				if (num_col % 100 == 0)
					printf("%u/%u collisions found, %.2lf %% done \n", num_col, counter, num_col * 100.0 / keep);
			}
		}
		if (num_col < keep)
			puts("WARNING: failed to get expected amount of reduced values!");
		std::cout << "completed after " << counter << " loops" << std::endl;
		puts("-------------------------");
		delete L1;
		delete R1;
		delete L2;
		delete R2;
	}
}


/* Sort-and-difference code */
void sort_and_difference(std::vector<SignatureSC25519>& sigs, const uint32_t n_bit, const uint32_t l, const uint32_t b,
        const uint32_t filter, const uint32_t a, const int log_prec, const size_t iota =
                1, const std::string out_prefix = "", const bool istest = true) {
	uint32_t threshold_bit = n_bit - filter;
	uint32_t S;

	for (size_t round = 0; round < iota; round++) {
		S = sigs.size();
		threshold_bit -= b;
		mpz_class threshold_mpz = mpz_class(1) << threshold_bit;
		sc25519 threshold_sc;
		mpz_to_gs(threshold_sc, threshold_mpz);

		std::sort(sigs.begin(), sigs.end());

		sc25519 newh, news;
		uint32_t j = 0;
		for (uint32_t i = 0; i < S - 1; i++) {
			sc25519_sub(&newh, &(sigs[i + 1].h), &(sigs[i].h));
			if (sc25519_lt(&newh, &threshold_sc)) {
				sc25519_sub(&news, &(sigs[i + 1].s), &(sigs[i].s));
				sigs[j] = SignatureSC25519(newh, news);
				j++;
			}
		}
		sigs.resize(j);
		printf("sorting done; %lu elements with h < 2^%u obtained\n", sigs.size(), threshold_bit);
		if (out_prefix.length()) {
			// file format: redsigs_sd_round-i.bin
			std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
			printf("saving signatures to %s... \n", outsig.c_str());
			sigsave(sigs, outsig);
		}
	}
	return;
}

void sort_and_difference(std::vector<SignatureSimple>& sigs, mpz_class& bound, const uint32_t n_bit, const uint32_t l,
        const uint32_t b, const uint32_t filter, const uint32_t a, const int log_prec,
        const size_t iota = 1, const std::string out_prefix = "", std::string dir = "", const bool istest = true) {
	uint32_t threshold_bit = n_bit - filter;
	mpz_class threshold_mpz;
	uint32_t S;
	if (dir.length())
		dir += "/";

	for (size_t round = 0; round < iota; round++) {
		S = sigs.size();
		threshold_bit -= b;
		threshold_mpz = mpz_class(1) << threshold_bit;
		std::sort(sigs.begin(), sigs.end());

		mpz_class newh, news;
		uint32_t j = 0;
		for (uint32_t i = 0; i < S - 1; i++) {
			newh = sigs[i + 1].h - sigs[i].h;
			if (newh < threshold_mpz) {
				news = sigs[i + 1].s - sigs[i].s;
				sigs[j] = SignatureSimple(newh, news);
				j++;
			}
		}
		sigs.resize(j);
		printf("sorting done; %lu elements with h < 2^%u obtained\n", sigs.size(), threshold_bit);
		bound = threshold_mpz.get_d();
		if (out_prefix.length()) {
			// file format: redsigs_sd_round-i.bin
			std::string outsig = out_prefix + "_round-" + std::to_string(round) + ".bin";
			printf("saving signatures to %s... \n", (dir + outsig).c_str());
			sigsave(sigs, dir + outsig);
		}
	}
	return;
}
