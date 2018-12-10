#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <omp.h>

std::string modify_info(std::string info, std::map < std::string , std::vector < long > > & reori){
  std::string output;
  std::size_t found;
  std::string key;
  std::string tmp;

  std::string first;
  std::string second;

  std::string CHR2;
  std::string END;
  long end;
  
  key = "CHR2=";
  found = info.find(key.c_str());
  first = info.substr(0,found+key.size());
  tmp = info.substr(found+key.size());
  key = ";END=";
  found = tmp.find(key.c_str());
  CHR2 = tmp.substr(0,found);
  tmp = tmp.substr(found);
  second = tmp.substr(0,key.size());
  tmp = tmp.substr(key.size());
  key = ";";
  found = tmp.find(key.c_str());
  END = tmp.substr(0,found);
  tmp = tmp.substr(found);

  end = strtol(END.c_str(),NULL,0);
  std::cerr << CHR2 << "\t" << end << "\t" << reori[CHR2][end] << std::endl;
  end += reori[CHR2][end];
  
  return first + CHR2 + second + std::to_string(end) + tmp ;
}


int main ( int argc , char ** argv ){
  if ( argc == 1 ) {
    std::cout << "Usage: calibrate  hg19.fa  dbSNP.vcf  SURVIVOR.vcf  outfile" << std::endl;
    return 0;
  }

  std::map < std::string , std::vector < long > > shift;
  
  //std::map < std::string , long > chr_size;
  long size;
  std::string id;
  std::string seq;
  std::string tmp;
  long pos;
  std::string ref;
  std::string alt;
  long diff;
  
  std::ifstream fin (argv[1]);
  getline ( fin , tmp );
  id = tmp.substr(1);
  while ( getline ( fin , tmp ) ){
    if ( tmp[0] == '>' ){
      size = seq.size();
      shift[id].resize(size,0);
      std::cerr << id << "\t" << size << std::endl;
      seq.clear();
      id = tmp.substr(1);
    }
    else{
      seq += tmp;
    }
  }
  size = seq.size();
  shift[id].resize(size,0);
  std::cerr << id << "\t" << size << std::endl;
  seq.clear();
  id = tmp.substr(0);
  fin.close();

  std::cerr << "reading genome complete\n" ;

  
  fin.open(argv[2]);
  while ( fin >> id >> pos >> tmp >> ref >> alt ){
    getline ( fin , tmp );
    if ( ref.size() != alt.size() ){
      diff = ref.size() - alt.size();
      shift[id][pos+1] = diff;
      //std::cout << "diff=" << diff << std::endl;
    }
  }
  fin.close();

  std::cerr << "reading dbSNP.vcf complete\n" ;

  std::map < std::string , std::vector < long > > reori;

  std::vector < std::string > id_vec;
  for ( auto & i : shift ){
    id_vec.push_back(i.first);
    reori[i.first].resize(size,0);
  }
  
  omp_set_num_threads(id_vec.size());
  {
#pragma omp parallel for schedule (static)
    for ( std::size_t i= 0 ; i < id_vec.size() ; i ++ ){
      int tid = omp_get_thread_num();
      std::string chr = id_vec[i];
      std::size_t size = shift[chr].size();
      std::cerr << "calibrating in " << chr << std::endl;
      for ( std::size_t j = 1 ; j < size ; j ++ ){
	reori[chr][j] = reori[chr][j-1] + shift[chr][j-1];
      }
    }
  }
  
  std::string chr;
  std::string header;
  std::string qual;
  std::string filter;
  std::string info;
  std::size_t found;
  fin.open(argv[3]);
  std::ofstream fout ( argv[4]);
  while ( getline ( fin , tmp ) ){
    if ( tmp.substr(0,2)=="##" ){
      header += tmp + "\n";
    }
    else {
      found = tmp.find("INFO");
      tmp = tmp.substr(0,found+4);
      header += tmp;
      break;
    }
  }
  fout << header << std::endl;
  while ( fin >> chr >> pos >> id >> ref >> alt >> qual >> filter >> info ){
    getline ( fin , tmp );
    pos += reori[chr][pos];
    info=modify_info(info,reori);
    //std::cout << pos << "\n";
    fout << chr << "\t"
	 << pos << "\t"
	 << id << "\t"
	 << ref << "\t"
	 << alt << "\t"
	 << qual << "\t"
	 << filter << "\t"
	 << info << "\n";
  }
  fin.close();
  fout.close();
  std::cerr << "shift complet\n";
  
  return 0;
}

