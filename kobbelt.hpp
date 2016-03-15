
#ifndef _KOBBELT_HPP_
#define _KOBBELT_HPP_

#include <cmath>
#include <map>
#include <array>

#include "accurate_math.hpp"
#include "genericfp.hpp"

template <typename T>
constexpr int sign(const T &val) {
  return (T(0) < val) - (val < T(0));
}

template <typename fptype>
int computeGenus(fptype val) {
  fpconvert<fptype> hwFloatFields = gfFPStruct(val);
  return hwFloatFields.exponent * 2 +
         (hwFloatFields.mantissa & 1);
}

template <typename fptype>
void tableInsert(std::map<int, fptype> &table, fptype val) {
  /* First determine where in the table the value is to go */
  int genus = computeGenus(val);
  if(table.count(genus) == 0) {
    /* There is no value here
     * We might be able to exactly add this number
     * and the one with the same exponent
     * but different final bit in the mantissa,
     * though, so check for that one.
     * If it exists and is of opposite sign,
     * then it will work
     */
    int otherGenus = genus ^ 1;
    if(table.count(otherGenus) > 0 &&
       sign(val) != sign(table[otherGenus])) {
      /* The other value exists and is of the correct sign,
       * so add them together, remove the other value,
       * and insert their sum
       */
      val += table[otherGenus];
      table.erase(otherGenus);
      tableInsert(table, val);
    } else {
      /* Nothing else to do, just insert it */
      table.insert({genus, val});
    }
  } else {
    /* There is already a value of the same genus,
     * so we can add them exactly,
     * remove the other value from the table,
     * and then insert their sum
     */
    val += table[genus];
    table.erase(genus);
    tableInsert(table, val);
  }
}

template <typename fptype, typename rettype>
rettype kobbeltDotProd(const fptype *v1, const fptype *v2,
                       const unsigned int size) {
  /* Start by inserting the exact products
   * of the values into a table ordered by their genus
   */
  std::map<int, fptype> table;
  for(unsigned int i = 0; i < size; i++) {
    std::array<fptype, 2> prod = twoProd(v1[i], v2[i]);
    tableInsert(table, prod[0]);
    tableInsert(table, prod[1]);
  }
  /* Now add them together in the order
   * from least genus to greatest
   */
  rettype ret = 0.0;
  for(auto kvpair : table) {
    ret += kvpair.second;
  }
  return ret;
}

#endif
