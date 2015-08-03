
#ifndef _KOBBELT_HPP_
#define _KOBBELT_HPP_

#include <cmath>
#include <map>
#include <array>

#include "accurate_math.hpp"
#include "genericfp.hpp"

template <typename fptype>
int fpGenus(fptype val) {
  fpconvert<fptype> bitwiseVal = gfFPStruct(val);
  int genus =
      bitwiseVal.exponent * 2 + (bitwiseVal.mantissa & 1);
  return genus;
}

template <typename fptype>
bool genusEqual(fptype lhs, fptype rhs) {
  int gLhs = fpGenus(lhs);
  int gRhs = fpGenus(rhs);
  return gLhs == gRhs;
}

template <typename fptype>
typename std::map<int, fptype>::iterator tableInsert(
    fptype value, std::map<int, fptype> table,
    typename std::map<int, fptype>::iterator hint) {
  int genus = fpGenus(value);
  int otherGenus = genus ^ 1;
  if(table.count(genus) > 0) {
    value += table.at(genus);
    table.erase(genus);
    return tableInsert(value, table, hint);
  } else if(table.count(otherGenus) > 0 &&
            std::signbit(value) !=
                std::signbit(table.at(otherGenus))) {
    /* Maintain the invariant that two values of opposite
     * sign
     * of the same genus are always added together.
     * This can be done exactly, and reduces complexity
     * later
     */
    value += table.at(otherGenus);
    table.erase(otherGenus);
    return tableInsert(value, table, hint);
  } else {
    /* No more values are in our way */
    return table.insert(hint, {genus, value});
  }
}

template <typename fptype>
void tableInsert(fptype value,
                 std::map<int, fptype> &table) {
  int genus = fpGenus(value);
  auto hint = table.find(genus);
  tableInsert(value, table, hint);
}

template <typename fptype>
fptype kobbeltDotProd(const fptype *vec1,
                      const fptype *vec2, unsigned len) {
  if(len == 0) return 0.0;
  std::map<int, fptype> fpTable;
  /* First exactly store the values of the products */
  for(unsigned i = 0; i < len; i++) {
    /* A minor simplification of Kobbelt's algorithm -
     * store only two values for each product,
     * instead of four.
     */
    std::array<fptype, 2> prod = twoProd(vec1[i], vec2[i]);
    tableInsert(prod[0], fpTable);
    tableInsert(prod[1], fpTable);
  }
  /* Now eliminate all values of opposite sign to the
   * highest value.
   * This prevents catastrophic losses of precision in the
   * last step, at the expense of potentially exponential
   * runtime.
   */
  auto itr = fpTable.end();
  /* We already know there must be at least one value in the
   * table, as the vectors are gauranteed to have at least
   * one value by this point.
   * Thus, itr cannot equal fpTable.begin(),
   * so it is safe to decrement
   */
  itr--;
  auto prevPos = itr;
  fptype prevVal = (*prevPos).second;
  while(itr != fpTable.begin()) {
    itr--;
    fptype curVal = (*itr).second;
    fpTable.erase(itr);
    if(std::signbit(prevVal) != std::signbit(curVal)) {
      /* Fix curVal's sign! */
      int curGenus = fpGenus(curVal);
      prevVal /= 2;
      int prevGenus = fpGenus(prevVal);
      auto erasePos = prevPos;
      while((prevGenus & ~1) >= (curGenus & ~1)) {
        prevPos = tableInsert(prevVal, fpTable, prevPos);
        prevVal /= 2;
        prevGenus -= 2;
      }
      fpTable.erase(erasePos);
      itr = prevPos;
      tableInsert(curVal + prevVal * 2, fpTable, prevPos);
    } else {
      prevPos = itr;
      prevVal = curVal;
    }
  }
  /* Finally, sum values in the table from the bottom up */
  fptype dotProd = 0.0;
  for(auto itr = fpTable.begin(); itr != fpTable.end();
      itr++)
    dotProd += (*itr).second;
  return dotProd;
}

#endif
