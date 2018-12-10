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

#include "realsimss.hpp"


#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>



int main ( int argc , char ** argv ){
  if ( argc != 3 ){
    std::cout << "Usage: RealSimSS  configure_file  outfile_prefix " << std::endl;
    return 0;
  }

  std::cout << "[Program start]\n" << currentDateTime() << std::endl;


  unsigned seed = 0;
  std::string key;
  std::ifstream fin ( argv[1] );
  while ( fin >> key ){
    if ( key == "seed" ){
      fin >> seed;
      break;
    }
  }
  fin.close();

  std::string pref=argv[2];
  std::string outfile;
  
  bool somatic_check=1;
  

  if ( seed == 0 ) seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "(SEED number for random number)=" << seed << std::endl;
  std::mt19937 generator (seed);


  std::cout << std::endl
	    << "[Reading reference genome file]" << std::endl ;  
  std::string conf_file=argv[1];
  GENOME germline(conf_file);
  GENOME somatic;
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;

  std::cout << std::endl 
	    << "[Reading dbSNP vcf and selecting variations]" << std::endl;
  VCF_CLASS dbSNP_vcf(conf_file);
  dbSNP_vcf.select_randomly(germline.SNP_number,generator);
  outfile = pref + ".selected_dbSNP.vcf";
  dbSNP_vcf.write(outfile);
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;


  std::cout << std::endl
	    << "[Simulation of SV]" << std::endl;
  VCF_CLASS germline_SV(conf_file); // germline SV
  VCF_CLASS somatic_SV(conf_file); // somatic SV

  std::cout << "Germline:\n" ;
  germline_SV.vcf_map = germline.germline_sim_SV(generator);
  outfile = pref + ".germline_SV.vcf";
  std::cout << "vcf file of germline SV: " << outfile << std::endl;
  germline_SV.write(outfile);

  somatic_check = germline.somatic_check();
  if ( somatic_check ){
    // copy germline genome to somatic genome
    std::cout << "Somatic:\n" ;
    somatic=germline;
    somatic_SV.vcf_map  = somatic.somatic_sim_SV(generator);
    outfile = pref + ".somatic_SV.vcf";
    std::cout << "vcf file of somatic SV: " << outfile << std::endl;
    somatic_SV.write(outfile);
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;
  
  
  
  std::cout << std::endl
	    << "[Add small variations]" << std::endl;
  VCF_CLASS germline_SNP(conf_file);
  VCF_CLASS somatic_SNP(conf_file);

  std::cout << "Germline:" << std::endl;
  germline_SNP.vcf_map = germline.add_dbSNP( dbSNP_vcf.vcf_map);
  outfile = pref + ".germline_dbSNP.vcf";
  germline_SNP.write(outfile);

  if ( somatic_check ){
    std::cout << "Somatic:" << std::endl;
    somatic_SNP.vcf_map  = somatic.add_dbSNP( dbSNP_vcf.vcf_map);
    outfile = pref + ".somatic_dbSNP.vcf";
    somatic_SNP.write(outfile);
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;


  std::cout << std::endl
	    << "[Merge SV and small variations]" << std::endl;
  VCF_CLASS germline_vcf(conf_file);
  VCF_CLASS somatic_vcf (conf_file);
  germline_vcf = germline_SV + germline_SNP;
  if ( somatic_check ){
    somatic_vcf = germline_SV + somatic_SV + somatic_SNP;
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;


  std::cout << std::endl
	    << "[Modifying genome]" << std::endl;
  std::cout << "Germline: germline_SV + dbSNP" << std::endl;
  germline.modify_genome(germline_vcf.vcf_map , generator );
  if ( somatic_check ){
    std::cout << "Somatic: germline_SV + somatic_SV + dbSNP" << std::endl;
      somatic.modify_genome(somatic_vcf.vcf_map , generator );
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;


  std::cout << std::endl 
	    << "[Writing genome file]" << std::endl;
  outfile = pref + ".germline.fasta";
  std::cout << "Germline: " << outfile << std::endl;
  germline.write(outfile);
  
  if ( somatic_check ){
    outfile = pref + ".somatic.fasta";
    std::cout << "Somatic: " << outfile << std::endl;
    somatic.write(outfile);
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;


  std::cout << std::endl
	    << "[Writing variation file]" << std::endl;
  outfile = pref + ".germline.vcf";
  std::cout << "Germline: germline_SV + dbSNP:\n" << outfile << std::endl;
  germline_vcf.write(outfile);
  if ( somatic_check ){
    outfile = pref + ".somatic.vcf";
    std::cout << "Somatic: germline_SV + somatic_SV + dbSNP:\n" << outfile << std::endl;
    somatic_vcf.write(outfile);
  }
  std::cout << currentDateTime() << std::endl << "[Done]" << std::endl;
  
  std::cout << std::endl << "[Program finished]" << std::endl;
  return 0;
}


