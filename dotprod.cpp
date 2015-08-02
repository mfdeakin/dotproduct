
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gtest/gtest.h>
#include <unistd.h>
#include <getopt.h>

#include <random>

#include "accurate_math.hpp"

template <typename fptype>
void genVector(
    fptype *vec, unsigned vecSize, std::mt19937_64 rgen,
    std::uniform_real_distribution<fptype> dist) {
  for (unsigned i = 0; i < vecSize; i++) {
    vec[i] = dist(rgen);
  }
}

template <typename fptype>
fptype dotProd(const fptype *v1, const fptype *v2,
               unsigned len) {
  fptype total = 0.0;
  for (unsigned i = 0; i < len; i++) {
    total += v1[i] * v2[i];
  }
  return total;
}

template <typename fptype>
fptype kahanDotProd(const fptype *v1, const fptype *v2,
                    unsigned len) {
  fptype total = 0.0;
  fptype c = 0.0;
  for(unsigned i = 0; i < len; i++) {
    fptype mod = fma(v1[i], v2[i], -c);
    fptype tmp = total + mod;
    c = (tmp - total) - mod;
    total = tmp;
  }
  return total;
}

template <typename fptype>
long double correctDotProd(const fptype *v1, const fptype *v2,
                           unsigned len) {
  long double total = 0.0;
  for (unsigned i = 0; i < len; i++) {
    total = fmal(v1[i], v2[i], total);
  }
  return total;
}

struct timespec subtractTimes(struct timespec start,
                              struct timespec end) {
  struct timespec delta;
  delta.tv_sec = end.tv_sec - start.tv_sec;
  delta.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (delta.tv_nsec < 0) {
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
  if (sum.tv_nsec > maxNSec) {
    sum.tv_nsec -= maxNSec;
    sum.tv_sec++;
  }
  return sum;
}

template <typename fptype>
struct testResult {
  struct timespec elapsedTime;
  fptype result;
};

template <typename fptype, typename retType,
          retType (*dp)(const fptype *, const fptype *,
                        unsigned len)>
testResult<retType> testFunction(fptype *vec1, fptype *vec2,
                                 unsigned len) {
  struct timespec start;
  int error =
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  assert(!error);
  retType result = dp(vec1, vec2, len);
  struct timespec end;
  error = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  assert(!error);
  struct timespec delta = subtractTimes(start, end);
  struct testResult<retType> ret = {
    delta, result
  };
  return ret;
}

void parseOptions(int argc, char **argv, int &testSize,
                  int &numTests) {
  int ret = 0;
  do {
    ret = getopt(argc, argv, "d:t:");
    switch (ret) {
      case 'd':
        testSize = atoi(optarg);
        break;
      case 't':
        numTests = atoi(optarg);
        break;
    }
  } while (ret != -1);
}

int main(int argc, char **argv) {
  int testSize = 1024;
  int numTests = 65536;
  parseOptions(argc, argv, testSize, numTests);

  float *vec1, *vec2;
  vec1 = (float *)malloc(sizeof(float[testSize]));
  vec2 = (float *)malloc(sizeof(float[testSize]));
  assert(vec1 != NULL);
  assert(vec2 != NULL);
  const float maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag,
                                              maxMag);
  struct timespec runningTimes[] = {
      {0, 0}, {0, 0}, {0, 0}, {0, 0}};
  long double totalErr[] = {0.0, 0.0, 0.0};
  for (int i = 0; i < numTests; i++) {
    genVector(vec1, testSize, engine, rgenf);
    genVector(vec2, testSize, engine, rgenf);
    struct testResult<long double> correctResult =
        testFunction<float, long double,
                     correctDotProd<float> >(vec1, vec2,
                                             testSize);
    runningTimes[0] = addTimes(correctResult.elapsedTime,
                               runningTimes[0]);

    struct testResult<float> spResult =
        testFunction<float, float, dotProd<float> >(
            vec1, vec2, testSize);
    runningTimes[1] =
        addTimes(spResult.elapsedTime, runningTimes[1]);
    double err1 =
        std::fabs(spResult.result - correctResult.result) /
        correctResult.result;
    totalErr[0] += logl(err1) / logl(10.0);

    struct testResult<float> compensatedResult =
        testFunction<float, float,
                     compensatedDotProd<float> >(vec1, vec2,
                                                 testSize);
    runningTimes[2] = addTimes(
        compensatedResult.elapsedTime, runningTimes[2]);
    double err2 = std::fabs(compensatedResult.result -
                            correctResult.result) /
                  correctResult.result;
    totalErr[1] += logl(err2) / logl(10.0);

    struct testResult<float> kahanResult =
      testFunction<float, float,
                   kahanDotProd<float> >(vec1, vec2,
                                         testSize);
    runningTimes[3] = addTimes(
        kahanResult.elapsedTime, runningTimes[3]);
    double err3 = std::fabs(kahanResult.result -
                            correctResult.result) /
                  correctResult.result;
    totalErr[2] += logl(err3) / logl(10.0);
  }
  printf(
      "Ran %d tests of size %d\n"
      "Correct Running Time: %ld.%09ld s\n"
      "Naive Time: %ld.%09ld s; Error %Le\n"
      "Compensated Time: %ld.%09ld s; Error %Le, "
      "Error-Naive: %Le\n"
      "Kahan Time: %ld.%09ld s; Error %Le, Error-Naive: "
      "%Le\n",
      numTests, testSize, runningTimes[0].tv_sec,
      runningTimes[0].tv_nsec, runningTimes[1].tv_sec,
      runningTimes[1].tv_nsec, totalErr[0] / numTests,
      runningTimes[2].tv_sec, runningTimes[2].tv_nsec,
      totalErr[1] / numTests,
      (totalErr[1] - totalErr[0]) / numTests,
      runningTimes[3].tv_sec, runningTimes[3].tv_nsec,
      totalErr[2] / numTests,
      (totalErr[2] - totalErr[0]) / numTests);
  free(vec2);
  free(vec1);
  return 0;
}
