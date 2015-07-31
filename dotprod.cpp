
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gtest/gtest.h>

#include <random>

#include "accurate_math.hpp"

template <typename fptype>
void genVector(fptype *vec, unsigned vecSize,
               std::mt19937_64 rgen,
               std::uniform_real_distribution<fptype> dist) {
  for(unsigned i = 0; i < vecSize; i++) {
    vec[i] = dist(rgen);
  }
}

template <typename fptype>
fptype dotProd(const fptype *v1, const fptype *v2, unsigned len) {
  fptype total = 0.0;
  for(unsigned i = 0; i < len; i++) {
    total += v1[i] * v2[i];
  }
  return total;
}

template <typename fptype>
double correctDotProd(const fptype *v1, const fptype *v2,
                      unsigned len) {
  double total = 0.0;
  for(unsigned i = 0; i < len; i++) {
    total += v1[i] * v2[i];
  }
  return total;
}

struct timespec subtractTimes(struct timespec start,
                              struct timespec end) {
  struct timespec delta;
  delta.tv_sec = end.tv_sec - start.tv_sec;
  delta.tv_nsec = end.tv_nsec - start.tv_nsec;
  if(delta.tv_nsec < 0) {
    delta.tv_sec--;
    delta.tv_nsec += 1e9;
  }
  return delta;
}

struct timespec addTimes(struct timespec t1,
                         struct timespec t2) {
  struct timespec sum;
  sum.tv_nsec = t1.tv_nsec + t2.tv_nsec;
  sum.tv_sec = t1.tv_sec + t2.tv_sec;
  const int maxNSec = 1000000000;
  if(sum.tv_nsec > maxNSec) {
     sum.tv_nsec -= maxNSec;
     sum.tv_sec++;
  }
  return sum;
}

int main(int argc, char **argv) {
  int testSize = 1024;
  int numTests = 65536;
  float *vec1, *vec2;
  vec1 = (float *)malloc(sizeof(float[testSize]));
  vec2 = (float *)malloc(sizeof(float[testSize]));
  assert(vec1 != NULL);
  assert(vec2 != NULL);
  const float maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag, maxMag);
  struct timespec runningTimes[3] = {{0, 0}, {0, 0}, {0, 0}};
  double totalErr[2] = {0.0, 0.0};
  for(int i = 0; i < numTests; i++) {
    genVector(vec1, testSize, engine, rgenf);

    struct timespec start;
    int error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    assert(!error);
    double correct = correctDotProd(vec1, vec2, testSize);
    struct timespec end;
    error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    assert(!error);
    struct timespec delta = subtractTimes(start, end);
    runningTimes[0] = addTimes(runningTimes[0], delta);
    
    error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    assert(!error);
    double bad1 = dotProd(vec1, vec2, testSize);
    error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    assert(!error);
    delta = subtractTimes(start, end);
    runningTimes[1] = addTimes(runningTimes[1], delta);
    totalErr[0] += fabs(correct - bad1);

    error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    assert(!error);
    double bad2 = compensatedDotProd(vec1, vec2, testSize);
    error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    assert(!error);
    delta = subtractTimes(start, end);
    runningTimes[2] = addTimes(runningTimes[2], delta);
    totalErr[1] += fabs(correct - bad2);
  }
  printf("Ran %d tests of size %d\n"
         "Correct Running Time: %ld.%09ld s\n"
         "Naive Time: %ld.%09ld s; Error %e\n"
         "Compensated Time: %ld.%09ld s; Error %e\n",
         numTests, testSize,
         runningTimes[0].tv_sec, runningTimes[0].tv_nsec,
         runningTimes[1].tv_sec, runningTimes[1].tv_nsec, totalErr[0],
         runningTimes[2].tv_sec, runningTimes[2].tv_nsec, totalErr[1]);
  free(vec2);
  free(vec1);
  return 0;
}
