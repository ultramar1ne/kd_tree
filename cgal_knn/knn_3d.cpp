#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef boost::tuple<Point_3,int>                           Point_and_int;
typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
  CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
  Traits_base>                                              Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;

#include <boost/tokenizer.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace boost;


int main()
{
  const unsigned int K = 5;
  std::vector<Point_3> points;
  std::vector<int>     indices;

    std::ifstream myFile("./survival.data");
    if (myFile) {  
      std::string buf;
      char_separator<char> sep(",");
      while (getline(myFile, buf)){
        tokenizer<char_separator<char>> tokens(buf, sep);
        int i = 0;
        std::vector<double> tmp(3);
        for (const auto& t : tokens) {
          if(i<3){
            if(i==2){
              tmp[i++] = 5*stod(t); //normarlize
            }else{
              tmp[i++] = stod(t);
            }
          }else{
            indices.push_back(stod(t));
          }  
        }
        points.push_back( Point_3(tmp[0],tmp[1],tmp[2]) );
      }
    }

  // Insert number_of_data_points in the tree
  Tree tree(boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin())),
            boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end())));
  // search K nearest neighbours
  Point_3 query(62, 59, 13);
  Distance tr_dist;
  K_neighbor_search search(tree, query, K);
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << "distance: "
              << tr_dist.inverse_of_transformed_distance(it->second) << " posistion: "
              << boost::get<0>(it->first)<< " tag:" << boost::get<1>(it->first) << std::endl;
  }
  return 0;
}