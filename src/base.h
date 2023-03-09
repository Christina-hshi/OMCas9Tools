/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

// For storing basic settings for the entire project

#pragma once

#ifndef BASE_H
#define BASE_H

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<unordered_map>
#include<set>
#include<list>
#include<stack>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<ctime>
#include<cstdlib>
//#include<cmath>
#include<cstring>
#include<thread>
#include<chrono>
#include<mutex>
#include<condition_variable>
#include<algorithm>
#include<limits>
#include<numeric>
#include<assert.h>
//#include<malloc.h>
//#include<malloc/malloc.h>
//#include<stdlib.h>
#include<queue>
#include<utility>
//#include<atomic>
#include<sys/resource.h>
#include<sys/time.h>
// #include <exception>

#include<boost/config.hpp>
#include<boost/tokenizer.hpp>
#include<boost/bimap.hpp>
#include<boost/program_options.hpp>
#include<boost/iostreams/filter/gzip.hpp>
#include<boost/iostreams/copy.hpp>
#include<boost/iostreams/filtering_streambuf.hpp>
#include<boost/atomic.hpp>
#include<boost/algorithm/string.hpp>

// Modules used as default in this project
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::vector;
using std::map;
using std::unordered_map;
using std::set;
using std::list;
using std::pair;
using std::bitset;
using std::array;
using std::stack;
using std::to_string;
using std::ios;

const std::string currentDateTime();

namespace MSG{
  inline void message(string str){
    cout<<"[Message] "<<str<<endl;
  }
  inline void warning(string str){
    cout<<"[Warning] "<<str<<endl;
  }
  // Show error message and exit the program.
  inline void error(string str){
    cout<<"[Error] "<<str<<endl;
    exit(1);
  }
};

#endif
