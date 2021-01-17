/*
  Coded by Jang-il Sohn
  
  This program was disigned to simulate SNPs or small variations
  implementing dbSNP vcf file.

  This program can be redistributed under GPLv3.
 */

#include "realsimss.hpp" 


inline VCF_MAP merge_VCF_MAP(VCF_MAP & input1, VCF_MAP & input2){
  VCF_MAP output;
  output = input1;
  for ( auto & i : input2 ){
    for ( auto & j : i.second){
      output[i.first][j.first]=j.second;
    }
  }
  return output;
}



inline std::string cut_info_in_line(std::string line){
  std::size_t found;
  std::string output;
  found = line.find('\t');
  output = line.substr(0,found);
  line = line.substr(found+1);
  
  found = line.find('\t');
  output += "\t" + line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  output += "\t" + line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  output += "\t" + line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  output += "\t" + line.substr(0,found) + "\t.\t.\t";
  line = line.substr(found+1);

  output += "SVMTHOD=dbSNP;";

  std::string key="dbSNPBuildID=";
  found = line.find(key);
  std::string tmp = line.substr(found+key.size());
  found =tmp.find(";");
  tmp = tmp.substr(0,found);
  output += key + tmp;

  key="VC=";
  found = line.find(key);
  tmp = line.substr(found+key.size());
  found = tmp.find(";");
  tmp = tmp.substr(0,found);
  output += ";" + key + tmp;


  return output;
}


VCF::VCF(std::string line){
  std::size_t found;

  found = line.find('\t');
  chrom = line.substr(0,found);
  line = line.substr(found+1);
  
  found = line.find('\t');
  pos =  atoi(line.substr(0,found).c_str());
  line = line.substr(found+1);

  found = line.find('\t');
  id =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  ref =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  alt =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  qual =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  filter =  line.substr(0,found);
  line = line.substr(found+1);

  info =  line;
}


void VCF::parse(std::string line){
  std::size_t found;

  found = line.find('\t');
  chrom = line.substr(0,found);
  line = line.substr(found+1);
  
  found = line.find('\t');
  pos =  atoi(line.substr(0,found).c_str());
  line = line.substr(found+1);

  found = line.find('\t');
  id =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  ref =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  alt =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  qual =  line.substr(0,found);
  line = line.substr(found+1);

  found = line.find('\t');
  filter =  line.substr(0,found);
  line = line.substr(found+1);

  info =  line;
}



std::string VCF::string(){
  std::string tmp;
  tmp = chrom + "\t" + std::to_string(pos) + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info ;
  return tmp;
}





void VCF_CLASS::write(std::string outfile){
  std::ofstream fout (outfile.c_str());
  fout << header << std::endl;
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      fout << j.second.string() << std::endl;
    }
  }

  fout.close();
}



std::size_t VCF_CLASS::size(){
  std::size_t sum = 0;
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      sum ++;
    }
  }
  return sum;
}



const std::string currentDate() {
    time_t     now = time(0); 
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct); // YYYY-MM-DD

    return buf;
}


const std::string currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct); // YYYY-MM-DD.HH:mm:ss

  return buf;
}


void VCF_CLASS::make_header(){
  header = "##fileformat=VCFv4.0\n";
  header = header
    + "##fileDate="+currentDate()+"\n"
    + "##source=RealSimSS\n"
    + "##reference="+reference+"\n";
  for ( auto & i : genome_info ){
    header +="##contig=<ID="+i.first+",length="+std::to_string(i.second)+">\n";
  }
  header = header
    + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
    + "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n"
    + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n"
    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n"
    + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n"
    + "##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description=\"First dbSNP Build for RS\">\n"
    + "##INFO=<ID=VC,Number=1,Type=String,Description=\"Variation Class of small variations from dbSNP\">\n"
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ;
}




void GENOME::read_genome(){
  std::ifstream fin ( ref_file.c_str() );
  std::string chr;
  std::string seq;
  std::string tmp;
  getline ( fin , chr );
  chr = chr.substr(1);
  while ( getline ( fin , tmp ) ){
    if ( tmp[0] == '>' ){
      sequence[chr]=seq;
      chr=tmp.substr(1);
      seq.clear();
    }
    else{
      seq += tmp;
    }
  }
  sequence[chr]=seq;
  fin.close();
  
}


GENOME::GENOME(std::string conf_file){
  ref_file=read_conf<std::string>(conf_file,"Reference");
  SNP_number=read_conf<std::size_t>(conf_file,"Number_of_SNP");
  // pref=read_conf<std::string>(conf_file,"Outfile_prefix");

  read_genome();

  std::ifstream fin(conf_file.c_str());
  std::string input;
  std::size_t value;

  germline_ins.number = 0 ;
  germline_del.number = 0 ;
  germline_inv.number = 0 ;
  germline_dup.number = 0 ;
  germline_tra.number = 0 ;

  germline_ins.min_length = 1000 ;
  germline_del.min_length = 1000 ;
  germline_inv.min_length = 1000 ;
  germline_dup.min_length = 1000 ;
  germline_tra.min_length = 1000 ;

  germline_ins.max_length = 10000 ;
  germline_del.max_length = 10000 ;
  germline_inv.max_length = 10000 ;
  germline_dup.max_length = 10000 ;
  germline_tra.max_length = 10000 ;

  somatic_ins.number = 0 ;
  somatic_del.number = 0 ;
  somatic_inv.number = 0 ;
  somatic_dup.number = 0 ;
  somatic_tra.number = 0 ;

  somatic_ins.min_length = 1000 ;
  somatic_del.min_length = 1000 ;
  somatic_inv.min_length = 1000 ;
  somatic_dup.min_length = 1000 ;
  somatic_tra.min_length = 1000 ;

  somatic_ins.max_length = 10000 ;
  somatic_del.max_length = 10000 ;
  somatic_inv.max_length = 10000 ;
  somatic_dup.max_length = 10000 ;
  somatic_tra.max_length = 10000 ;


  while ( fin >> input ){
    if(input=="INS"){
      fin >> input ;
      if ( input == "germline" ){
	fin >> value; germline_ins.number = value;
	fin >> value; germline_ins.min_length=value;
	fin >> value; germline_ins.max_length=value;
	
      }
      else if (input == "somatic" ){
	fin >> value; somatic_ins.number = value;
	fin >> value; somatic_ins.min_length=value;
	fin >> value; somatic_ins.max_length=value;
      }
      else{
	std::cerr << "!!!ERROR!!!" << std::endl;
	std::cerr << "Check conf file. It must has been modified." << std::endl;
	assert(0);
      }
    }
    else if(input=="DEL"){
      fin >> input ;
      if ( input == "germline" ){
	fin >> value; germline_del.number = value;
	fin >> value; germline_del.min_length=value;
	fin >> value; germline_del.max_length=value;
      }
      else if (input == "somatic" ){
	fin >> value; somatic_del.number = value;
	fin >> value; somatic_del.min_length=value;
	fin >> value; somatic_del.max_length=value;
      }
      else{
	std::cerr << "!!!ERROR!!!" << std::endl;
	std::cerr << "Check conf file. It must has been modified." << std::endl;
	assert(0);
      }
    }
    else if(input=="DUP"){
      fin >> input ;
      if ( input == "germline" ){
	fin >> value; germline_dup.number = value;
	fin >> value; germline_dup.min_length=value;
	fin >> value; germline_dup.max_length=value;
      }
      else if (input == "somatic" ){
	fin >> value; somatic_dup.number = value;
	fin >> value; somatic_dup.min_length=value;
	fin >> value; somatic_dup.max_length=value;
      }
      else{
	std::cerr << "!!!ERROR!!!" << std::endl;
	std::cerr << "Check conf file. It must has been modified." << std::endl;
	assert(0);
      }
    }
    else if(input=="INV"){
      fin >> input ;
      if ( input == "germline" ){
	fin >> value; germline_inv.number = value;
	fin >> value; germline_inv.min_length=value;
	fin >> value; germline_inv.max_length=value;
      }
      else if (input == "somatic" ){
	fin >> value; somatic_inv.number = value;
	fin >> value; somatic_inv.min_length=value;
	fin >> value; somatic_inv.max_length=value;
      }
      else{
	std::cerr << "!!!ERROR!!!" << std::endl;
	std::cerr << "Check conf file. It must has been modified." << std::endl;
	assert(0);
      }
    }
    else if(input=="TRA"){
      fin >> input ;
      if ( input == "germline" ){
	fin >> value; germline_tra.number = value;
	fin >> value; germline_tra.min_length=value;
	fin >> value; germline_tra.max_length=value;
      }
      else if (input == "somatic" ){
	fin >> value; somatic_tra.number = value;
	fin >> value; somatic_tra.min_length=value;
	fin >> value; somatic_tra.max_length=value;
      }
      else{
	std::cerr << "!!!ERROR!!!" << std::endl;
	std::cerr << "Check conf file. It must has been modified." << std::endl;
	assert(0);
      }
    }
  }

  fin.close();


  std::size_t seq_length;
  for ( auto & i : sequence ){
    
    seq_length = i.second.size();
    genome_info[i.first]=seq_length;

    total_length += seq_length;
    sequence_size_sum[i.first]=total_length;
    masking[i.first].resize(seq_length+1,1);
  }
  
}


template < class T >
T GENOME::read_conf(std::string conf, std::string key){
  std::ifstream fin(conf.c_str());
  std::string input;
  T value;
    while ( fin >> input ){
      if ( key == input ){
	fin >> value;
	break;
      }
    }
  fin.close();
  return value;
}


//void GENOME::modify_genome( VCF_MAP & vcf_map, std::mt19937 & generator){
void GENOME::modify_genome( VCF_MAP & vcf_map){

  // VCF_MAP output_vcf_map ;


  std::string chr;
  std::string seq;
  std::string output;

  for ( VCF_MAP :: iterator it = vcf_map.begin() ; it != vcf_map.end() ; it ++ ){
    chr = it->first;
    std::size_t variation_number = 0 ;
    if ( sequence.find(chr) != sequence.end() ) {
      seq = sequence[chr];
      std::vector < std::string > inter;
      std::vector < std::string > substitute;
  
      std::string front;
      std::string back;
      std::size_t refsize;
      std::size_t from;
      std::size_t to;
      std::size_t cut_size;

      VCF tmp_vcf;

      std::cout << "chr" << chr << "\tsize=" << seq.size() ;


      std::map < std::size_t , VCF > :: iterator jt = it->second.begin();
      tmp_vcf = jt->second;
      substitute.push_back(tmp_vcf.alt);
      variation_number++;
      front = sequence[chr].substr(0,tmp_vcf.pos-1);

      jt ++;

      std::map < std::size_t , VCF > :: iterator kt = it->second.begin();
      while ( jt != it->second.end() ){
	tmp_vcf = kt->second;
	refsize = tmp_vcf.ref.size();
	from = tmp_vcf.pos - 1 + refsize - 1;

	tmp_vcf = jt->second;
	to = tmp_vcf.pos - 1 ;

	if ( to <= from ){
	  std::cout << std::endl;
	  std::cout << kt->second.chrom << "\t"  << kt->second.pos << "\t" << kt->second.ref.size() << std::endl;
	  std::cout << jt->second.chrom << "\t"  << jt->second.pos << "\t" << jt->second.ref.size() << std::endl;
	}
	assert ( to > from );
	
	variation_number++;
	substitute.push_back(tmp_vcf.alt);
	
	cut_size = to - from ;
	output = seq.substr(from,cut_size);
	inter.push_back(output);

	jt ++ ;
	kt ++;
      }
      tmp_vcf = kt->second;
      refsize = tmp_vcf.ref.size();
      from = tmp_vcf.pos - 1 + refsize ;
      back  = sequence[chr].substr(from);


      output = front ;
      std::size_t order ;
      for ( order = 0 ; order < inter.size() ; order ++ ){
	output += substitute[order] + inter[order] ;
      }
      output += substitute[order] + back;
    
      sequence[chr]=output;

      // std::cout << " --> " << sequence[chr].size() << "\tNo_of_variations=" << output_vcf_map[chr].size() << "\n";
      std::cout << " --> " << sequence[chr].size() << "\tNo_of_variations=" << variation_number << "\n";
    }
  }
}

void GENOME::write(std::string outfile){
  
  // std::string outfile = pref;
  // outfile += ".fasta";
  std::ofstream fout(outfile.c_str());

  for ( auto & i : sequence ){
    fout << ">" << i.first << "\n" << i.second << "\n" ;
  }
  fout.close();
}



VCF_MAP GENOME::germline_sim_SV( std::mt19937 & generator ){
  VCF_MAP ins = sim_ins (germline_ins,generator);
  VCF_MAP del = sim_del (germline_del,generator);
  VCF_MAP dup = sim_dup (germline_dup,generator);
  VCF_MAP inv = sim_inv (germline_inv,generator);
  VCF_MAP tra = sim_tra (germline_tra,generator);

  VCF_MAP output = merge_VCF_MAP(ins,del);
  output = merge_VCF_MAP ( output, dup );
  output = merge_VCF_MAP ( output, inv );
  output = merge_VCF_MAP ( output, tra );
  return output;
}

VCF_MAP GENOME::somatic_sim_SV( std::mt19937 & generator ){
  VCF_MAP ins = sim_ins (somatic_ins,generator);
  VCF_MAP del = sim_del (somatic_del,generator);
  VCF_MAP dup = sim_dup (somatic_dup,generator);
  VCF_MAP inv = sim_inv (somatic_inv,generator);
  VCF_MAP tra = sim_tra (somatic_tra,generator);

  VCF_MAP output = merge_VCF_MAP(ins,del);
  output = merge_VCF_MAP ( output, dup );
  output = merge_VCF_MAP ( output, inv );
  output = merge_VCF_MAP ( output, tra );
  return output;
}




////////////////////////////////////////////////////////////////////////////////////////////
//
// small variations (SNPs and indels)
//

VCF_MAP GENOME::add_dbSNP( VCF_MAP & SNP_vcf ){
  // std::cout << "Add of small variations using dbSNP" << std::endl;

  VCF_MAP output ;
  std::size_t num_of_snp=0;
    
  for ( auto & i : SNP_vcf ){
    std::string chr=i.first;
    for ( auto & j : i.second ){
      VCF vcf=j.second;
      std::size_t from=j.first;
      std::size_t to=from+vcf.ref.size();
      if (to >= genome_info[chr]) to = genome_info[chr] - 1; 
      if ( from < to ){
	bool masked=1;
	for ( std::size_t k = from ; k < to ; k ++ ){
	  masked = masked && masking[chr][k];
	}
	if ( masked ){
	  output[chr][from]=vcf;
	  for ( std::size_t k = from ; k < to ; k ++ ){
	    masking[chr][k]=0;
	  }
	  num_of_snp ++;
	}
      }
    }
  }
  std::cout << "Number of added small variations: " << num_of_snp << std::endl;
  return output;
}


inline std::size_t end_pos( VCF vcf){
  std::string info = vcf.info;
  std::size_t found = info.find("END=");
  info = info.substr(found+4);
  found = info.find(";");
  info = info.substr(0,found);
  return atoi(info.c_str());
}



bool GENOME::somatic_check(){
  return somatic_ins.number + somatic_del.number + somatic_dup.number + somatic_inv.number + somatic_tra.number;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

VCF_CLASS::VCF_CLASS(std::string conf){
  read_info(conf);
  make_header();
}


void VCF_CLASS::read_info( std::string conf){
  std::ifstream fin(conf.c_str());
  std::string input;
  std::string tmp;
  std::string id;

  while ( fin >> input ){
    if ( input == "Reference"){
      fin >> input;
      reference = input;
    }
    else if ( input == "SNP_vcf_file" ){
      fin >> input;
      vcf_file = input;
    }
  }
  fin.close();

  fin.open( reference.c_str());
  getline ( fin , input );
  id = input.substr(1);
  while ( fin >> input ){
    if ( input[0] == '>' ){
      genome_info[id] = tmp.size();
      tmp.clear();
      id = input.substr(1);
    }
    else {
      tmp += input;
    }
  }
  genome_info[id] = tmp.size();
  tmp.clear();
  id.clear();
  fin.close();

  
}



void VCF_CLASS::read_vcf(){


  std::ifstream fin(vcf_file.c_str());
  std::string tmp;
  VCF vcf;


  while ( getline ( fin , tmp ) ){
    if ( tmp[0] == '#' ){
      header += tmp + "\n";
    }
    else{
      vcf.parse(tmp);
      if ( genome_info.find(vcf.chrom) != genome_info.end() ){
	vcf_map[vcf.chrom][vcf.pos]=vcf;
	break;
      }
    }
  }



  while ( getline ( fin , tmp ) ){
    vcf.parse(tmp);
    if ( genome_info.find(vcf.chrom) != genome_info.end() ){
      vcf_map[vcf.chrom][vcf.pos]=vcf;
    }
  }


  fin.close();
}






void VCF_CLASS::read_vcf_unbal(){


  std::ifstream fin( vcf_file.c_str());
  std::string tmp;
  VCF vcf;

  while ( getline ( fin , tmp ) ){
    if ( tmp[0] == '#' ){
      header += tmp + "\n";
    }
    else{
      vcf.parse(tmp);
      if ( genome_info.find(vcf.chrom) != genome_info.end() && vcf.ref.size() != vcf.alt.size() ){
	vcf_map[vcf.chrom][vcf.pos]=vcf;
	break;
      }
    }
  }


  while ( getline ( fin , tmp ) ){
    vcf.parse(tmp);
    if ( genome_info.find(vcf.chrom) != genome_info.end() && vcf.ref.size() != vcf.alt.size() ){
      vcf_map[vcf.chrom][vcf.pos]=vcf;
    }
  }


  fin.close();

}



void VCF_CLASS::select_alt ( std::string & alt , std::mt19937 & generator){
  std::size_t found = alt.find(",");
  if ( found != std::string::npos ){
    std::string input=alt;
    std::vector< std::string > input_vec;
    select_alt ( input, input_vec);
    std::uniform_int_distribution<std::size_t> distribution(0,input_vec.size()-1);
    alt = input_vec[distribution(generator)];
  }
}

void VCF_CLASS::select_alt ( std::string & input, std::vector< std::string > & input_vec){
  std::size_t found = input.find(",");
  if ( found == std::string::npos ){
    input_vec.push_back(input);
  }
  else{
    input_vec.push_back(input.substr(0,found));
    input=input.substr(found+1);
    select_alt(input,input_vec);
  }
}



void VCF_CLASS::select_randomly(std::size_t SNP_number, std::mt19937 & generator){
  std::size_t count = 0 ;
  std::string tmp;
  std::ifstream fin ( vcf_file.c_str() );
  std::vector<VCF> vcf_vec;
  VCF vcf;

  while ( getline ( fin , tmp ) ){
    if ( tmp[0]!='#'){
      count ++;
      break;
    }
  }
  while ( getline ( fin , tmp ) ){
    count ++;
  }
  fin.close();
  
  std::cout << count << " lines in input vcf" << std::endl;
  if ( SNP_number > count ) SNP_number = count;


  std::set<std::size_t> selected;
  std::size_t pick;

  std::uniform_int_distribution<std::size_t> distribution(1,count);

  fin.open(vcf_file.c_str());
  count = 0 ;
  while ( getline ( fin , tmp ) ){
    if ( tmp[0]!='#'){
      pick = distribution(generator);
      tmp=cut_info_in_line(tmp);
      vcf.parse(tmp);
      select_alt(vcf.alt,generator);
      if ( pick <= SNP_number && genome_info.find(vcf.chrom) != genome_info.end() ) {
	vcf_map[vcf.chrom][vcf.pos]=vcf;
	count ++;
      }
      break;
    }
  }
  while ( getline ( fin , tmp ) ){
    pick = distribution(generator);
    tmp=cut_info_in_line(tmp);
    vcf.parse(tmp);
    select_alt(vcf.alt,generator);
    if ( pick <= SNP_number && genome_info.find(vcf.chrom) != genome_info.end() ){
      vcf_map[vcf.chrom][vcf.pos]=vcf;
      count ++;
    }
  }
  fin.close();
  std::cout << count << " lines selected" << std::endl;
}



void VCF_CLASS::clear(){
  vcf_map.clear();
  genome_info.clear();
  vcf_file.clear();
  reference.clear();
  header.clear();
}



