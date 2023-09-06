// Copyright 2023 mimicry AG
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ==============================================================================

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

// Sort2mu - sort 2 uint32_t values pointed to by buffer
void Sort2mu(uint32_t* buffer) {
	asm volatile(
		"mov (%0), %%eax          \n"
		"mov 4(%0), %%edx         \n"
		"cmp %%eax, %%edx         \n"
		"sbb %%rbx, %%rbx         \n"
		"mov %%edx, 4(%0,%%rbx,4) \n"
		"neg %%rbx                \n"
		"mov %%eax, (%0,%%rbx,4)  \n"
		: "+r"(buffer)
		:
		: "eax","edx","rbx","memory");
}

// Sort3mu - sort 3 uint32_t values pointed to by buffer
uint8_t shuf3mu0[] = {
	8, 4, 8, 0, 4, 0
};

uint8_t shuf3mu1[] = {
	4, 8, 0, 8, 0, 4
};

void Sort3mu(uint32_t* buffer) {
	asm volatile(
		"mov (%0), %%eax                        \n"
		"mov 4(%0), %%edx                       \n"
		"cmp %%eax, %%edx                       \n"
		"sbb %%rbx, %%rbx                       \n"
		"mov 8(%0), %%ecx                       \n"
		"cmp %%eax, %%ecx                       \n"
		"sbb %%r8, %%r8                         \n"
		"cmp %%ecx, %%edx                       \n"
		"adc $0, %%r8                           \n"
		"movzb shuf3mu0+3(%%rbx,%%r8,2), %%r9   \n"
		"mov %%ecx, 4(%0,%%r8,4)                \n"
		"movzb shuf3mu1+3(%%rbx,%%r8,2), %%rcx  \n"
		"mov %%eax, (%0,%%r9)                   \n"
		"mov %%edx, (%0,%%rcx)                  \n"
		: "+r"(buffer)
		:
		: "eax","edx","rbx","rcx","r8","r9","memory");
}

// Sort3ms - sort 3 int32_t values pointed to by buffer;
// reuses the write shuffle vectors from Sort3mu
void Sort3ms(int32_t* buffer) {
	asm volatile(
		"mov (%0), %%eax \n"
		"mov 4(%0), %%edx \n"
		"lea 0x80000000(%%eax), %%r10d \n"
		"lea 0x80000000(%%edx), %%r11d \n"
		"cmp %%r10d, %%r11d \n"
		"sbb %%rbx, %%rbx \n"		// range -1..0
		"mov 8(%0), %%ecx \n"
		"lea 0x80000000(%%ecx), %%r12d \n"
		"cmp %%r10d, %%r12d \n"
		"sbb %%r8, %%r8 \n"
		"cmp %%r12d, %%r11d \n"
		"adc $0, %%r8 \n"	// range -1..1
		"movzb shuf3mu0+3(%%rbx,%%r8,2), %%r9 \n"
		"mov %%ecx, 4(%0,%%r8,4) \n"
		"movzb shuf3mu1+3(%%rbx,%%r8,2), %%rcx \n"
		"mov %%eax, (%0,%%r9) \n"
		"mov %%edx, (%0,%%rcx) \n"
		: "+r"(buffer)
		:
		: "eax","edx","rbx","rcx","r8","r9","r10","r11","r12","memory");
}

// Sort3mv - sort 3 int32_t values pointed to by buffer, using SIMD
// instructions

// first three elements mask for masked move (vpmaskmovd) instruction
int32_t mask3[] = { -1, -1, -1, 0 }; 

// byte shuffles representing 32-bit word positions
const uint32_t p0 = 0x03020100;
const uint32_t p1 = 0x07060504;
const uint32_t p2 = 0x0b0a0908;
const uint32_t p3 = 0x0f0e0d0c;
const uint32_t invp = 0xffffffff;

/* Sort3 shuffles
0b000 -> 0 1 2
0b001 -> 1 2 0
0b010 -> 2 0 1
0b011 -> 2 1 0
0b100 -> 0 1 2
0b101 -> 1 0 2
0b110 -> 0 2 1
*/

// In the shuffle definition, if we want to retain the last (unsorted) word,
// the last element needs to be set to p3, if we want to clear it to 0, the
// last element must be invp
uint32_t vshuf3[7*4] = {
	p0, p1, p2, p3,
	p1, p2, p0, p3,
	p2, p0, p1, p3,
	p2, p1, p0, p3,
	p0, p1, p2, p3,
	p1, p0, p2, p3,
	p0, p2, p1, p3
};

void Sort3mv(int32_t* buffer) {
	asm volatile(
		"vmovdqa mask3, %%xmm1                 \n"
		"vpmaskmovd (%0), %%xmm1, %%xmm2       \n"
		"vpshufd $0xc9, %%xmm2, %%xmm3         \n"	// shuffle 0b1100_1001 -> 1 2 0 3
		"vpcmpgtd %%xmm3, %%xmm2, %%xmm3       \n"
		"movmskps %%xmm3, %%eax                \n"
		"shl $4, %%eax                         \n"
		"vpshufb vshuf3(%%eax), %%xmm2, %%xmm2 \n"
		"vpextrq $0, %%xmm2, (%0)              \n"
		"vpextrd $2, %%xmm2, 8(%0)             \n"
		: "+r"(buffer)
		:
		: "eax","xmm1","xmm2","xmm3","memory");
}

// Sort4mu - sort 4 uint32_t values pointed to by buffer
int8_t shuf4mu0[] = {	// src position 0
	2, 3, 3, 1, 1, 2, -1, -1,		// last elem to pos 0
	2, 3, 3, 0, 0, 2, -1, -1,		// last elem to pos 1
	1, 3, 3, 0, 0, 1, -1, -1,		// last elem to pos 2
	1, 2, 2, 0, 0, 1, -1, -1		// last elem to pos 3
};

int8_t shuf4mu1[] = {	// src position 1
	1, 2, 1, 3, 2, 3, -1, -1,		// last elem to pos 0
	0, 2, 0, 3, 2, 3, -1, -1,		// last elem to pos 1
	0, 1, 0, 3, 1, 3, -1, -1,		// last elem to pos 2
	0, 1, 0, 2, 1, 2, -1, -1		// last elem to pos 3
};

int8_t shuf4mu2[] = {	// src position 2
	3, 1, 2, 2, 3, 1, -1, -1,		// last elem to pos 0
	3, 0, 2, 2, 3, 0, -1, -1,		// last elem to pos 1
	3, 0, 1, 1, 3, 0, -1, -1,		// last elem to pos 2
	2, 0, 1, 1, 2, 0, -1, -1		// last elem to pos 3
};

void Sort4mu(uint32_t* buffer) {
	asm volatile(
		"mov (%0), %%eax                        \n"
		"mov 4(%0), %%edx                       \n"
		"cmp %%eax, %%edx                       \n"
		"sbb %%rbx, %%rbx                       \n"
		"mov 8(%0), %%ecx                       \n"
		"mov 12(%0), %%r10d                     \n"
		"cmp %%eax, %%ecx                       \n"
		"adc %%rbx, %%rbx                       \n"
		"cmp %%eax, %%r10d                      \n"
		"sbb %%r11, %%r11                       \n"
		"cmp %%r10d, %%edx                      \n"
		"adc $0, %%r11                          \n"
		"cmp %%ecx, %%edx                       \n"
		"adc %%rbx, %%rbx                       \n"
		"cmp %%r10d, %%ecx                      \n"
		"adc $0, %%r11                          \n"
		"movzb shuf4mu0+11(%%rbx,%%r11,8), %%r8 \n"
		"mov %%r10d,4(%0,%%r11,4)               \n"
		"movzb shuf4mu1+11(%%rbx,%%r11,8), %%r9 \n"
		"mov %%eax, (%0,%%r8,4)                 \n"
		"movzb shuf4mu2+11(%%rbx,%%r11,8), %%r8 \n"
		"mov %%edx, (%0,%%r9,4)                 \n"
		"mov %%ecx, (%0,%%r8,4)                 \n"
		: "+r"(buffer)
		:
		: "eax","rbx","ecx","edx","r8","r9","r10","r11","memory");
}

// Sort4ms - sort 4 int32_t values pointed to by buffer;
// reuses the write shuffle vectors from Sort4mu
void Sort4ms(int32_t* buffer) {
	asm volatile(
		"mov (%0), %%eax                        \n"
		"mov 4(%0), %%edx                       \n"
		"lea 0x80000000(%%eax), %%r12d          \n"
		"lea 0x80000000(%%edx), %%r13d          \n"
		"cmp %%r12d, %%r13d                     \n"
		"sbb %%rbx, %%rbx                       \n"
		"mov 8(%0), %%ecx                       \n"
		"mov 12(%0), %%r10d                     \n"
		"lea 0x80000000(%%ecx), %%r14d          \n"
		"lea 0x80000000(%%r10d), %%r9d          \n"
		"cmp %%r12d, %%r14d                     \n"
		"adc %%rbx, %%rbx                       \n"
		"cmp %%r12d, %%r9d                      \n"
		"sbb %%r11, %%r11                       \n"
		"cmp %%r9d, %%r13d                      \n"
		"adc $0, %%r11                          \n"
		"cmp %%r14d, %%r13d                     \n"
		"adc %%rbx, %%rbx                       \n"
		"cmp %%r9d, %%r14d                      \n"
		"adc $0, %%r11                          \n"
		"movzb shuf4mu0+11(%%rbx,%%r11,8), %%r8 \n"
		"mov %%r10d,4(%0,%%r11,4)               \n"
		"movzb shuf4mu1+11(%%rbx,%%r11,8), %%r9 \n"
		"mov %%eax, (%0,%%r8,4)                 \n"
		"movzb shuf4mu2+11(%%rbx,%%r11,8), %%r8 \n"
		"mov %%edx, (%0,%%r9,4)                 \n"
		"mov %%ecx, (%0,%%r8,4)                 \n"
		: "+r"(buffer)
		:
		: "eax","rbx","ecx","edx","r8","r9","r10","r11","r12","r13","r14","memory");
}

// Sort4mv - sort 4 int32_t values pointed to by buffer, using SIMD
// instructions

// Explanation: Sort4mv performs the six comparisons
//      v0 > v1, v1 > v2, v2 > v3, v3 > v0, v0 > v2, v1 > v3
// and packs the comparison result into a 6-bit vector, used for selecting
// a read shuffle vector. Note that not all indices are possible (filled with
// entries indicating invalid below), and that other combinations of six
// comparisons could be used to select the shuffle, with function finder
// having selected one that can be executed with a minimal number of shuffle
// operations in preparation for comparing values.
uint32_t shuffles[64*4] = {
	p0, p1, p2, p3,	        // index 0
	invp, invp, invp, invp,	// index 1
	invp, invp, invp, invp,	// index 2
	invp, invp, invp, invp,	// index 3
	p1, p2, p3, p0,	        // index 4
	p1, p2, p3, p0,	        // index 5
	invp, invp, invp, invp,	// index 6
	invp, invp, invp, invp,	// index 7
	invp, invp, invp, invp,	// index 8
	p2, p3, p0, p1,        	// index 9
	p0, p2, p3, p1,	        // index 10
	p2, p3, p0, p1,	        // index 11
	invp, invp, invp, invp,	// index 12
	p2, p1, p3, p0,	        // index 13
	invp, invp, invp, invp,	// index 14
	p2, p3, p1, p0,        	// index 15
	p0, p3, p1, p2,	        // index 16
	invp, invp, invp, invp,	// index 17
	p3, p0, p1, p2,	        // index 18
	invp, invp, invp, invp,	// index 19
	p1, p3, p0, p2,	        // index 20
	p1, p3, p2, p0,	        // index 21
	p3, p1, p0, p2,	        // index 22
	p3, p1, p2, p0,	        // index 23
	invp, invp, invp, invp,	// index 24
	invp, invp, invp, invp,	// index 25
	p3, p0, p2, p1,        	// index 26
	p3, p2, p0, p1,	        // index 27
	invp, invp, invp, invp,	// index 28
	invp, invp, invp, invp,	// index 29
	invp, invp, invp, invp,	// index 30
	p3, p2, p1, p0,	        // index 31
	p0, p1, p2, p3,	        // index 32
	invp, invp, invp, invp,	// index 33
	invp, invp, invp, invp,	// index 34
	invp, invp, invp, invp,	// index 35
	p1, p0, p2, p3,	        // index 36
	p1, p2, p0, p3,	        // index 37
	invp, invp, invp, invp,	// index 38
	invp, invp, invp, invp,	// index 39
	p0, p2, p1, p3,	        // index 40
	p2, p0, p1, p3,	        // index 41
	p0, p2, p3, p1,	        // index 42
	p2, p0, p3, p1,	        // index 43
	invp, invp, invp, invp,	// index 44
	p2, p1, p0, p3,	        // index 45
	invp, invp, invp, invp,	// index 46
	invp, invp, invp, invp,	// index 47
	p0, p1, p3, p2,       	// index 48
	invp, invp, invp, invp,	// index 49
	p0, p3, p1, p2,	        // index 50
	invp, invp, invp, invp,	// index 51
	p1, p0, p3, p2,	        // index 52
	invp, invp, invp, invp,	// index 53
	invp, invp, invp, invp,	// index 54
	invp, invp, invp, invp,	// index 55
	invp, invp, invp, invp,	// index 56
	invp, invp, invp, invp,	// index 57
	p0, p3, p2, p1,	        // index 58
	invp, invp, invp, invp,	// index 59
	invp, invp, invp, invp,	// index 60
	invp, invp, invp, invp,	// index 61
	invp, invp, invp, invp,	// index 62
	invp, invp, invp, invp	// index 63
};

void Sort4mv(int32_t* buffer) {
	asm volatile(
		"vmovdqu (%0), %%xmm1            \n"
		"vpshufd $0x39, %%xmm1, %%xmm3   \n"	// shuffle 0b00111001 -> 1, 2, 3, 0
		"vpshufd $0xee, %%xmm1, %%xmm4   \n"  // shuffle 0b11101110 -> 2, 3, 2, 3
		"vpcmpgtd %%xmm3, %%xmm1, %%xmm3 \n"
		"vpcmpgtd %%xmm4, %%xmm1, %%xmm4 \n"
		"movmskps %%xmm3, %%eax          \n"
		"movmskps	%%xmm4, %%ecx          \n"
		"lea 0(%%ecx,%%eax,4), %%eax     \n"
		"shl $4, %%eax                   \n"
		"movdqa shuffles(%%eax), %%xmm2  \n"
		"vpshufb %%xmm2, %%xmm1, %%xmm1  \n"
		"vmovdqu %%xmm1, (%0)            \n"
		: "+r"(buffer)
		:
		: "xmm1","xmm2","xmm3","xmm4","eax","ecx");
}

// StdSort is a wrapper around std::sort to make it fit the signature of our
// sorting functions
template<typename T, size_t len>
void StdSort(T* buffer) {
	std::sort(buffer, buffer+len);
}

// equalValues - are two arrays pointed to by p and q equal for len elements
template<typename T, size_t len>
bool equalValues(T* p, T* q) {
	for (auto i = len; i > 0; i--) {
		if (p[i-1] != q[i-1]) {
			return false;
		}
	}
	return true;
}

// isSorted - return whether the array pointed to by p of size len is sorted
template<typename T, size_t len>
bool isSorted(T* p) {
	T sp[len];		// temporary buffer for a sorted copy
	std::copy(p, p+len, sp);
	std::sort(std::begin(sp), std::end(sp));
	return equalValues<T, len>(p, sp);
}

// slice is a helper for outputting slices
template<typename T>
struct slice {
	T* p;
	size_t l;
	slice(T* buf, size_t len) { p = buf; l = len; }
};

template<typename T>
std::ostream& operator<<(std::ostream& o, slice<T> sl) {
	o << "[";
	for (size_t i = 0; i < sl.l-1; i++) {
		o << sl.p[i] << ",";
	}
	o << sl.p[sl.l-1] << "]";
	return o;
}

// powUint calculates x^y, without overflow check
uint64_t powUint(uint64_t x, uint64_t y) {
	uint64_t res = 1, pot = x;
	while (y > 0) {
		if ((y&1) == 1) res *= pot;
		pot *= pot;
		y >>= 1;
	}
	return res;
}

// rshuffle implements a read shuffle; shuf includes n positions like a
// digit in a base-n number
template<typename T, size_t n>
void rshuffle(T src[], T dest[], uint64_t shuf) {
	for (size_t i = 0; i < n; i++) {
		dest[i] = src[shuf%n];
		shuf /= n;
	}
}

// checkSortFunction determines all permutations of array v including
// repetitions of individual elements (of which there are len^len), sorts
// them, and checks whether the array is the sorted copy of the original
// permutation; outputs failures and returns whether it was able to sort
// all permutations
template<typename T, size_t len>
bool checkSortFunction(void f(T*), T v[]) {
	bool failed = false;
	auto nofPerm = powUint(len, len);		// number of permutations with repetitions
	for (size_t shuf = 0; shuf < nofPerm; shuf++) {
		T org[len], work[len];
		rshuffle<T, len>(v, org, shuf);
		std::copy(org, org+len, work);
		f(work);
		if (!isSorted<T, len>(work)) {
			failed = true;
			std::cout << "failed to sort " << slice<T>(org, len) << ", result " << slice<T>(work, len) << "\n";
		}
	}
	return !failed;
}

// workBody iterates sorting function F over the allocated array, changing
// the offset to obtain different alignments
template<typename T, void F(T*)>
void workBody(T work[], size_t offsetLen, size_t nofArrays) {
	for (size_t offset = 0; offset < offsetLen; offset++) {
		for (size_t i = 0; i < nofArrays; i++) {
			F(&work[offset+i*offsetLen]);
		}
	}
}

// benchSortFunction runs provided sort function on a provided array of
// random values. It checks that the values would be sorted in a first
// run (just to be on the safe side), and then sorts nofArrays subarrays,
// shifted by all offsets from 0 upto less than offsetLen, which provides
// for different alignments of the subarrays to be benchmarked. It executes
// the benchmark runs times, re-initializing the arrays to its initial state
// after each run. It returns the total number of nanoseconds performing
// sorting.
template<typename T, size_t len>
std::chrono::nanoseconds benchSortFunction(void F(T*, size_t, size_t),
																					T v[], size_t nofArrays, size_t offsetLen, size_t runs) {
	T work[(nofArrays+1)*offsetLen];
	// measure execution time as elapsed time (to be improved if needed)
	std::chrono::nanoseconds res(0);
	auto t0 = std::chrono::steady_clock::now();
	for (size_t r = 0; r < runs; r++) {
		std::copy(v, v+(nofArrays+1)*offsetLen, work);	// reset arrays

		// body to measure
		F(work, offsetLen, nofArrays);

	}
	auto t1 = std::chrono::steady_clock::now();
	return std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0);
}

int main() {

	// correctness checking
	int32_t vs[] = { -10, 0, 10, 20 };
	uint32_t vu[] = { 0, 5, 10, 15 };
	std::cout << "checking sort functions for correct sorting on all permutations\n";
	auto success = checkSortFunction<int32_t,2>(StdSort<int32_t,2>, vs);
	std::cout << "StdSort2: " << (success? "success\n":"failed\n");
	success = checkSortFunction<uint32_t,2>(Sort2mu, vu);
	std::cout << "Sort2mu: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,3>(StdSort<int32_t,3>, vs);
	std::cout << "StdSort3: " << (success? "success\n":"failed\n");
	success = checkSortFunction<uint32_t,3>(Sort3mu, vu);
	std::cout << "Sort3mu: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,3>(Sort3ms, vs);
	std::cout << "Sort3ms: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,3>(Sort3mv, vs);
	std::cout << "Sort3mv: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,4>(StdSort<int32_t,4>, vs);
	std::cout << "StdSort4: " << (success? "success\n":"failed\n");
	success = checkSortFunction<uint32_t,4>(Sort4mu, vu);
	std::cout << "Sort4mu: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,4>(Sort4ms, vs);
	std::cout << "Sort4ms: " << (success? "success\n":"failed\n");
	success = checkSortFunction<int32_t,4>(Sort4mv, vs);
	std::cout << "Sort4mv: " << (success? "success\n":"failed\n");
	
	// benchmarking
	std::minstd_rand generator;
	size_t const n = 501*4;
	// 500 subarrays of length 4, one extra so we can have different alignment
	// by adding an offset
	std::vector<uint32_t> au;
	std::uniform_int_distribution<> udist(0, 10000);
	for (size_t i = 0; i < n; i++) {
		au.push_back((uint32_t)udist(generator));
	}
	std::vector<int32_t> as;
	std::uniform_int_distribution<> sdist(-10000, 10000);
	for (size_t i = 0; i < n; i++) {
		as.push_back(sdist(generator));
	}

	size_t const nofRuns = 1000000;
	double const nofSorts = 500.0*4.0*nofRuns;
	std::cout << "benchmarking sort functions (" << nofSorts << " sorts)\n";
	auto ns = benchSortFunction<int32_t, 2>(workBody<int32_t, StdSort<int32_t, 2>>, &as[0], 500, 4, nofRuns);
	std::cout << "StdSort2: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<uint32_t, 2>(workBody<uint32_t, Sort2mu>, &au[0], 500, 4, nofRuns);
	std::cout << "Sort2mu: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 3>(workBody<int32_t, StdSort<int32_t, 3>>, &as[0], 500, 4, nofRuns);
	std::cout << "StdSort3: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<uint32_t, 3>(workBody<uint32_t, Sort3mu>, &au[0], 500, 4, nofRuns);
	std::cout << "Sort3mu: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 3>(workBody<int32_t, Sort3ms>, &as[0], 500, 4, nofRuns);
	std::cout << "Sort3ms: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 3>(workBody<int32_t, Sort3mv>, &as[0], 500, 4, nofRuns);
	std::cout << "Sort3mv: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 4>(workBody<int32_t, StdSort<int32_t, 4>>, &as[0], 500, 4, nofRuns);
	std::cout << "StdSort4: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<uint32_t, 4>(workBody<uint32_t, Sort4mu>, &au[0], 500, 4, nofRuns);
	std::cout << "Sort4mu: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 4>(workBody<int32_t, Sort4ms>, &as[0], 500, 4, nofRuns);
	std::cout << "Sort4ms: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";
	ns = benchSortFunction<int32_t, 4>(workBody<int32_t, Sort4mv>, &as[0], 500, 4, nofRuns);
	std::cout << "Sort4mv: " << ns.count() << " ns, " << double(ns.count())/nofSorts << " ns/sort\n";

	return 0;
}
