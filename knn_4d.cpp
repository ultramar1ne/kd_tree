#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Epick_d.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>   // instead of traits_3
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Dimension.h>

#define CGAL_EIGEN3_ENABLED

typedef CGAL::Epick_d<CGAL::Dimension_tag<4> >              Kernel;
typedef Kernel::Point_d Point_d;
typedef boost::tuple<Point_d,int>                           Point_and_int;
typedef CGAL::Search_traits_d<Kernel, CGAL::Dimension_tag<4> >    Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
  CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
  Traits_base>                                              Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;


using namespace std;
using namespace boost;


int getLabel(string s){
  if( !(s.compare("\"Setosa\"") )){
    return 0;
  }else if ( !s.compare("\"Versicolor\"")){
    return 1;
  }else if (!s.compare("\"Virginica\"")){
    return 2;
  }else{
    return -1;
  }
}

const unsigned int K = 3;

int main() {

    std::vector<Point_d> points;
    std::vector<int> labels;

    std::ifstream myFile("./iris.csv");
    if (myFile) {  
      std::vector<double> tmp(4);
      std::string buf;
      getline(myFile, buf); //ignore lines[0]
      while (getline(myFile, buf)){
        char_separator<char> sep(", ");
        tokenizer<char_separator<char>> tokens(buf, sep);
        int i = 0;
        for (const auto& t : tokens) {
          int label = getLabel(t);
          if (label ==-1){
            tmp[i++]=stod(t);
          }else{
            labels.push_back(label);
          }
        }
        points.push_back( Point_d(4, tmp.begin(),tmp.end()) );
      }
    }
    
  // Insert number_of_data_points in the tree
  Tree tree(boost::make_zip_iterator(boost::make_tuple( points.begin(),labels.begin())),boost::make_zip_iterator(boost::make_tuple( points.end(),labels.end())));
  
  // search K nearest neighbours
  Point_d query(0.0, 0.0, 0.0,0.0);
  Distance tr_dist;

  K_neighbor_search search(tree, query, K);

  return 0;
}