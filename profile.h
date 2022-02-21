/*
 * profile.h
 *
 *  Created on: 26 окт. 2021 г.
 *      Author: User
 */
//#pragma once;
#ifndef PROFILE_H_
#define PROFILE_H_
#pragma once
#include <chrono>
#include <iostream>
#include <string>
#include<sstream>
using namespace std;
using namespace std::chrono;
class LogDuration {
public:
  explicit LogDuration(const string& msg = "")
    : message(msg + ": ")
    , start(steady_clock::now())
  {
  }

  ~LogDuration() {
    auto finish = steady_clock::now();
    auto dur = finish - start;
    cerr << message
       << duration_cast<milliseconds>(dur).count()
       << " ms" << endl;
  }
private:
  string message;
  steady_clock::time_point start;
};

#define UNIQ_ID_IMPL(lineno) _a_local_var_##lineno
#define UNIQ_ID(lineno) UNIQ_ID_IMPL(lineno)

#define LOG_DURATION(message) \
  LogDuration UNIQ_ID(__LINE__){message};
struct TotalDuration {
string message;
steady_clock::duration value;
explicit TotalDuration(const string& msg = "")
: message(msg + ": ")
, value(0)
{
}
~TotalDuration() {
ostringstream os;
os << message
<< duration_cast<milliseconds>(value).count()
<< " ms" << endl;
cerr << os.str();
}
};
class AddDuration {
public:
explicit AddDuration(steady_clock::duration& dest)
: add_to(dest)
, start(steady_clock::now())
{
}
explicit AddDuration(TotalDuration& dest)
: AddDuration(dest.value)
{
}
~AddDuration() {
add_to += steady_clock::now() - start;
}
private:
steady_clock::duration& add_to;
steady_clock::time_point start;
};
#define ADD_DURATION(value) \
AddDuration UNIQ_ID(__LINE__){value};
#endif /* PROFILE_H_ */
