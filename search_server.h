/*
 * search_server.h
 *
 *  Created on: 26 окт. 2021 г.
 *      Author: User
 */

#ifndef SEARCH_SERVER_H_
#define SEARCH_SERVER_H_
#pragma once
#include "profile.h"
#include <istream>
#include <ostream>
#include <set>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include<string_view>
#include<mutex>
#include<thread>
#include<future>
using namespace std;
template <typename T>
class Synchronized {
public:
  explicit Synchronized(T initial = T())
    : value(move(initial))
  {}

  struct Access {
    T& ref_to_value;
    lock_guard<mutex> g;
  };

  Access GetAccess() {
    return {value, lock_guard(m)};
  }

private:
  T value;
  mutex m;
};
class InvertedIndex {
public:
  void Add(const string& document);
  vector<pair<size_t,size_t>> Lookup(const string& word) const;
  size_t Get_count()
  {
	  return last_docid;
  }
private:
  unordered_map<string, vector<pair<size_t,size_t>>> index;
  size_t last_docid = 0;
};

class SearchServer {
public:
  SearchServer() = default;
  explicit SearchServer(istream& document_input);//, TotalDuration& dur);
  void UpdateDocumentBase(istream& document_input);//, TotalDuration& dur);
  void AddQueriesStream(istream& query_input, ostream& search_results_output);
  //void UpdateDocumentBase(istream& document_input, TotalDuration& dur);
  //void AddQueriesStream(istream& query_input, ostream& search_results_output, TotalDuration& dur);

private:
  Synchronized<InvertedIndex> index;
  mutex m;
  vector<future<void>> f;
};





#endif /* SEARCH_SERVER_H_ */
