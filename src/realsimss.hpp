/*
  Coded by Jang-il Sohn
  
  This program was disigned to simulate SNPs or small variations
  implementing dbSNP vcf file.

  This program can be redistributed under GPLv3.
 */

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <assert.h>
#include <iterator>
#include <time.h>

#ifndef MY_VCF
#define MY_VCF


const std::string currentDate();
const std::string currentDateTime();


class VCF{
public:
  VCF(){};
  VCF(std::string);
  ~VCF(){};

  void parse(std::string);

  std::string chrom;
  std::size_t pos;
  std::string id;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
  
  std::string string();
};

void select_alt ( std::string & , std::vector< std::string > & );
void select_alt ( std::string & , std::mt19937 & generator);

inline std::size_t end_pos( VCF );

using VCF_SUB_MAP=std::map < std::size_t , VCF >;
using VCF_MAP=std::map < std::string, VCF_SUB_MAP >;

VCF_MAP merge_VCF_MAP(VCF_MAP & , VCF_MAP & );




class VCF_CLASS{
private:
  void select_alt ( std::string & , std::mt19937 & );
  void select_alt ( std::string & , std::vector< std::string > & );
  void make_header();
public:
  VCF_CLASS(){
    make_header();
  };
  VCF_CLASS(std::string);
  ~VCF_CLASS(){};

  std::map < std::string , std::size_t >  genome_info;

  std::string vcf_file;
  std::string reference;

  std::string header;

  VCF_MAP vcf_map;
  
  void read_info(std::string);

  void read_vcf();
  void read_vcf_unbal();
  
  void select_randomly(std::size_t, std::mt19937 &);

  void clear();

  VCF_MAP::iterator begin(){return vcf_map.begin();};
  VCF_MAP::iterator end(){return vcf_map.end();};

  VCF_CLASS operator+(VCF_CLASS & input){
    VCF_CLASS output;
    std::string chr;
    std::size_t pos;
    output.vcf_map = this->vcf_map;
    for ( auto & i : input.vcf_map ){
      chr = i.first;
      for ( auto & j : i.second ){
	pos = j.first;
	output.vcf_map[chr][pos]=j.second;
      }
    }
    output.reference=this->reference;
    return output;
  }

  void operator += (VCF_CLASS & input ){
    std::string chr;
    std::size_t pos;
    for ( auto & i : input.vcf_map ){
      chr = i.first;
      for ( auto & j : i.second ){
	pos = j.first;
	this->vcf_map[chr][pos]=j.second;
      }
    }
  }
  
  void write(std::string);
  std::size_t size();
};

#endif




#ifndef MY_GENOME
#define MY_GENOME

class SVTYPE{
public:
  SVTYPE(){};
  ~SVTYPE(){};

  std::size_t number;
  std::size_t min_length;
  std::size_t max_length;

  std::map < std::string , std::size_t > dist;

};


class GENOME{
private:
  template < class T >
  T read_conf(std::string conf, std::string key);
  void read_genome();
  std::map < std::string , std::size_t >  genome_info;
  std::map < std::string, std::size_t > sequence_size_sum;
  std::size_t total_length = 0 ;
  
public:
  GENOME(){};
  GENOME(std::string);
  ~GENOME(){};

  std::string ref_file;

  std::size_t SNP_number;

  SVTYPE germline_ins;
  SVTYPE germline_del;
  SVTYPE germline_dup;
  SVTYPE germline_inv;
  SVTYPE germline_tra;

  SVTYPE somatic_ins;
  SVTYPE somatic_del;
  SVTYPE somatic_dup;
  SVTYPE somatic_inv;
  SVTYPE somatic_tra;

  std::map < std::string, std::string > sequence;

  // std::map < std::string, std::string > novel_insertion;
  // std::map < std::string, std::map < std::size_t , std::string > > ins_sequence;

  std::map < std::string , std::vector < bool > > masking;

  VCF_MAP sim_ins ( SVTYPE & , std::mt19937 & );
  VCF_MAP sim_del ( SVTYPE & , std::mt19937 & );
  VCF_MAP sim_dup ( SVTYPE & , std::mt19937 & );
  VCF_MAP sim_inv ( SVTYPE & , std::mt19937 & );
  VCF_MAP sim_tra ( SVTYPE & , std::mt19937 & );

  VCF_MAP germline_sim_SV( std::mt19937 & );
  VCF_MAP somatic_sim_SV( std::mt19937 & );
  VCF_MAP add_dbSNP( VCF_MAP & );
  void modify_genome( VCF_MAP &, std::mt19937 & );
  void write(std::string); 
  
  bool somatic_check();
};



#endif
