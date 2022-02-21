/*
 * search_server.cpp
 *
 *  Created on: 26 окт. 2021 г.
 *      Author: User
 */
#include "search_server.h"
#include "iterator_range.h"
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iostream>
using namespace std;
void Print(const array<size_t, 5>& docid)
{
	for(int i = 0; i < 5; i++)
		cout<< docid[i]<< " ";
	cout<<endl;
}
void Print_(const map<string, vector<pair<size_t, size_t>>>& index)
{
	for(const auto& [word, vec] : index)
	{
		cout<<word<<" : "<<endl;
		for(const auto& p: vec)
		{
			cout<<p.first<<" "<<p.second<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}
vector<string> SplitIntoWords(const string& line) {
  istringstream words_input(line);
  return {istream_iterator<string>(words_input), istream_iterator<string>()};
}

SearchServer::SearchServer(istream& document_input) {
  UpdateDocumentBase(document_input);
}

void SearchServer::UpdateDocumentBase(istream& document_input) {
  InvertedIndex new_index;
  for (string current_document; getline(document_input, current_document); ) {
    new_index.Add(move(current_document));
  }
  auto access = index.GetAccess();
  	access.ref_to_value = move(new_index);

}
void AddQueriesStreamAsync(istream& query_input, ostream& search_results_output, Synchronized<InvertedIndex>& index) {
	size_t count = index.GetAccess().ref_to_value.Get_count();
	for (string current_query; getline(query_input, current_query); ) {
		    const auto words = SplitIntoWords(current_query);
		    vector<size_t> docid_count(count, 0);
		    for (const auto& word : words) {
		    	auto access = index.GetAccess();
		    	for (const auto &p : access.ref_to_value.Lookup(word)) {
		        docid_count[p.first] += p.second;
		      }
		    }
		    array<size_t, 5> search_results = {0, 0, 0, 0, 0};
		    array<size_t, 5> index_results = {0, 0, 0, 0, 0};
		    for(size_t i = 0; i < count; i++)
		    {
		    	if(docid_count[i] > search_results[4]){
		    		search_results[4] = docid_count[i];
		    		index_results[4] = i;
		    		int j = 3;
		    		while( search_results[j+1] > search_results[j])
		    		{

		    			swap(search_results[j+1], search_results[j]);
		    			swap(index_results[j+1], index_results[j]);
		    			if(j == 0) break;
		    			j--;
		    		}
		    	}
		    }
		    //Print(search_results);
		    search_results_output << current_query << ':';
		    for (int i = 0; i < 5; i++) {
		    	if(search_results[i] > 0){
		    	search_results_output << " {"
		        << "docid: " << index_results[i] << ", "
		        << "hitcount: " << search_results[i] << '}';
		    	}
		    }
		    search_results_output<<endl;
		  }
}
void SearchServer::AddQueriesStream(istream& query_input, ostream& search_results_output) {
	f.push_back(async(launch::async, AddQueriesStreamAsync, ref(query_input), ref(search_results_output), ref(index)));

}
void InvertedIndex::Add(const string& document) {
  for (auto& word : SplitIntoWords(document)) {
	  auto it = index.find(word);
	  if(it != index.end())
	  {
		  auto it1 = prev(it->second.end());
		  if(it1->first != last_docid)
			it->second.push_back(make_pair(last_docid, 1));
		  else
			it1->second++;
		  }
	  else
	  {
		  index[move(word)].push_back(make_pair(last_docid, 1));
	  }
  }
  last_docid++;
}

vector<pair<size_t, size_t>> InvertedIndex::Lookup(const string& word) const {
  if (auto it = index.find(word); it != index.end()) {
    return it->second;
  } else {
    return {};
  }
}



