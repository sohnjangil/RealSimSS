#include <string>
#include <iostream>
#include <fstream>
#include <map>


#include "realsimss.hpp"

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "Usage:  calibrator  genome.fasta  selected_variation.vcf  SV.vcf  outfile_prefix\n" ;
    return 0;
  }
  
  std::string genome_file(argv[1]);
  std::string selected_file(argv[2]);
  std::string sv_file(argv[3]);
  std::string prefix(argv[4]);

  GENOME genome(genome_file.c_str());
  
  VCF_CLASS variation1;
  VCF_CLASS variation2;

  variation1.read_vcf_unbal(genome,selected_file);
  variation2.read_vcf(genome,sv_file);

  std::string id;
  std::size_t size;

  std::map < std::string , std::size_t > genome_size_map;
  std::map < std::string , std::vector < int > > shift_vector;
  std::map < std::string , std::map < std::size_t , int > > shift_map;

  for ( auto & i : genome.sequence){
    id = i.first;
    size = i.second.size();
    genome_size_map[id]=size;
    shift_vector[id].resize(size+1,0);
  }
  
  VCF vcf;
  std::size_t position;
  int diff;
  for ( auto & i : variation1 ){
    id = i.first;
    size = genome_size_map[id];
    for ( auto & j : i.second ){
      vcf = j.second;
      position = vcf.pos;
      diff = vcf.ref.size() - vcf.alt.size();
      if ( diff ){
	shift_map[id][position]=diff;
      }
    }
    shift_map[id][size]=0;
  }
    
  std::size_t from ;
  std::size_t to;
  
  for ( auto & i : shift_map ){
    id = i.first;
    from = 0 ;
    to = 0 ;
    for ( auto & j : i.second ){
      from = to ;
      to = j.first + 1 ;
      diff = j.second ;
      for ( std::size_t pos = from ; pos < to ; pos ++ ){
	shift_vector[id][pos] = diff;
      }
    }
  }

  
  ////////////////////////////////////////
  //
  // TO-DO
  //
  for ( auto & i : variation2 ){
    id = i.first ;
    for ( auto & j : i.second ){
      j.
    }
  }


  return 0;
}
