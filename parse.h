/*
 * parse.h
 *
 *  Created on: 26 окт. 2021 г.
 *      Author: User
 */

#ifndef PARSE_H_
#define PARSE_H_
#pragma once
#include "iterator_range.h"
#include <string_view>
#include <sstream>
#include <vector>
using namespace std;

template <typename Container>
string Join(char c, const Container& cont) {
  ostringstream os;
  for (const auto& item : Head(cont, cont.size() - 1)) {
    os << item << c;
  }
  os << *rbegin(cont);
  return os.str();
}

string_view Strip(string_view s);
vector<string_view> SplitBy(string_view s, char sep);
#endif /* PARSE_H_ */
