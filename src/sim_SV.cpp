#include "realsimss.hpp"
#include "flip_strand.cpp"
 
////////////////////////////////////////////////////////////////////////////////
//
// simulation of insertion
//

VCF_MAP GENOME::sim_ins( SVTYPE & ins, std::mt19937 & generator){
  VCF vcf;
  std::size_t pick;
  std::string selected_id;


  std::uniform_int_distribution<std::size_t> whole_genome_dist(1,total_length);

  for ( std::size_t i = 0 ; i < ins.number ; i ++ ){
    pick = whole_genome_dist(generator);
    for ( auto & j : sequence_size_sum ){
      selected_id = j.first;
      if ( pick < j.second ){
	break;
      }
    }
    ins.dist[selected_id] ++ ;
  }

  // std::cout << "Simulation of insertion" << std::endl;

  VCF_MAP ins_map;

  char nucl[4];
  nucl[0]='A';
  nucl[1]='C';
  nucl[2]='G';
  nucl[3]='T';
  std::uniform_int_distribution<std::size_t> random_nucl(0,3);

  std::size_t ins_count=0;

  for ( auto & i : ins.dist ){
    std::size_t size = genome_info[i.first];
    std::uniform_int_distribution<std::size_t> position(1,size);
    std::uniform_int_distribution<std::size_t> length(ins.min_length,ins.max_length);
    std::size_t ins_position;
    std::size_t ins_length;
    for ( std::size_t j = 0 ; j < i.second ; ){
      ins_position=position(generator);
      ins_length=length(generator);

      bool masked=masking[i.first][ins_position];
      masked = masked && masking[i.first][ins_position+1];

      if (masked && sequence[i.first][ins_position-1]!='N' && sequence[i.first][ins_position-1]!='n' &&
	  ins_position > 1000 && ins_position < size - 1000 ){
	vcf.chrom=i.first;
	vcf.pos=ins_position;
	vcf.id="INS_"+std::to_string(ins_count++);
	vcf.ref=sequence[vcf.chrom][ins_position-1];

	vcf.alt=vcf.ref;
	for ( std::size_t j = 0 ; j < ins_length ; j ++ ){
          vcf.alt += nucl[random_nucl(generator)];
        }

	vcf.qual=".";
	vcf.filter=".";
	vcf.info="SVTYPE=INS;SVMETHOD=RealSimSS;CHR2=";
	vcf.info+=vcf.chrom+";END="+std::to_string(ins_position+1) +";SVLEN=" + std::to_string(ins_length);


	ins_map[vcf.chrom][vcf.pos]=vcf;

	masking[vcf.chrom][vcf.pos]=0;
	
	j ++;
      }
    }
  }
  std::cout << "Number of simulated INS: " << ins_count << std::endl;
  return ins_map;
}





////////////////////////////////////////////////////////////////////////////////
//
// Simulation of deletion
//

VCF_MAP GENOME::sim_del( SVTYPE & del, std::mt19937 & generator){
  VCF vcf;

  std::size_t pick;
  std::string selected_id;


  std::uniform_int_distribution<std::size_t> whole_genome_dist(1,total_length);


  for ( std::size_t i = 0 ; i < del.number ; i ++ ){
    pick = whole_genome_dist(generator);
    for ( auto & j : sequence_size_sum ){
      selected_id = j.first;
      if ( pick < j.second ){
	break;
      }
    }
    del.dist[selected_id] ++ ;
  }

  // std::cout << "Simulation of deletion" << std::endl;
  
  VCF_MAP del_map;

  std::size_t del_count = 0 ;
  for ( auto & i : del.dist ){
    std::size_t size = genome_info[i.first];
    std::uniform_int_distribution<std::size_t> position(1,size);
    std::uniform_int_distribution<std::size_t> length(del.min_length,del.max_length);
    std::size_t del_position;
    std::size_t del_length;
    for ( std::size_t j = 0 ; j < i.second ; ){
      del_position=position(generator);
      del_length=length(generator);
      
      int Ncount = 0 ;
      for ( std::size_t k = del_position-1 ; k < del_position+del_length ; k ++ ){
	if ( sequence[i.first][k]=='N' ||  sequence[i.first][k]=='n' ) Ncount ++ ;
      }
      Ncount *= 100;
      Ncount /= (int) del_length ;

      bool masked=1;
      for ( std::size_t k = del_position ; k < del_position + del_length + 1 ; k ++ ){
	masked = masked && masking[i.first][k];
      }

      if (sequence[i.first][del_position-1]!='N' && sequence[i.first][del_position-1+del_length] != 'N' &&
	  sequence[i.first][del_position-1]!='n' && sequence[i.first][del_position-1+del_length] != 'n' &&
	  Ncount < 10 && masked && del_position > 1000 && del_position < size - del_length - 1000 ){

	vcf.chrom=i.first;
	vcf.pos=del_position;
	vcf.id="DEL_"+std::to_string(del_count++);
	vcf.ref=sequence[vcf.chrom].substr(vcf.pos-1,del_length+1);
	vcf.alt=sequence[vcf.chrom][vcf.pos-1];
	vcf.qual=".";
	vcf.filter=".";
	vcf.info="SVTYPE=DEL;SVMETHOD=RealSimSS;CHR2=";
	vcf.info+=vcf.chrom+";END="+std::to_string(del_position+del_length) +";SVLEN=" + std::to_string(del_length);
	
	del_map[vcf.chrom][del_position]=vcf;
	
	for ( std::size_t k = del_position ; k < del_position + del_length + 1 ; k ++ ){
	  masking[vcf.chrom][k] = 0 ;
	}
	
	j ++ ;
      }
    }
  }
  std::cout << "Number of simulated DEL: " << del_count << std::endl;
  return del_map;
}








////////////////////////////////////////////////////////////////////////////////
//
// duplication
//

VCF_MAP GENOME::sim_dup( SVTYPE & dup, std::mt19937 & generator){
  VCF vcf;

  std::size_t pick;
  std::string selected_id;


  std::uniform_int_distribution<std::size_t> whole_genome_dist(1,total_length);


  for ( std::size_t i = 0 ; i < dup.number ; i ++ ){
    pick = whole_genome_dist(generator);
    for ( auto & j : sequence_size_sum ){
      selected_id = j.first;
      if ( pick < j.second ){
	break;
      }
    }
    dup.dist[selected_id] ++ ;
  }

  // std::cout << "Simulation of duplication" << std::endl;

  VCF_MAP dup_map;

  std::size_t dup_count=0;
  for ( auto & i : dup.dist ){
    std::size_t size = genome_info[i.first];
    std::uniform_int_distribution<std::size_t> position(1,size);
    std::uniform_int_distribution<std::size_t> length(dup.min_length,dup.max_length);
    std::size_t dup_position;
    std::size_t dup_length;
    for ( std::size_t j = 0 ; j < i.second ; ){
      dup_position=position(generator);
      dup_length=length(generator);

      int Ncount = 0 ;
      for ( std::size_t k = dup_position-1 ; k < dup_position+dup_length ; k ++ ){
	if ( sequence[i.first][k]=='N' || sequence[i.first][k]=='n' ) Ncount ++ ;
      }
      Ncount *=100;
      Ncount /=(int) dup_length ;

      bool masked=1;
      for ( std::size_t k = dup_position ; k < dup_position + dup_length + 1; k ++ ){
        masked = masked && masking[i.first][k];
      }

      if (sequence[i.first][dup_position-1]!='N' && sequence[i.first][dup_position-1+dup_length] !='N' &&
          sequence[i.first][dup_position-1]!='n' && sequence[i.first][dup_position-1+dup_length] !='n' &&
          Ncount < 10 && masked && dup_position > 1000 && dup_position < size - 1000 ){
	vcf.chrom=i.first;
	vcf.pos=dup_position;
	vcf.id="DUP_"+std::to_string(dup_count++);
	vcf.ref=sequence[i.first].substr(vcf.pos-1,dup_length);
	vcf.alt=vcf.ref;
	vcf.alt+=vcf.ref;
	vcf.qual=".";
	vcf.filter=".";
	vcf.info="SVTYPE=DUP;SVMETHOD=RealSimSS;CHR2=";
	vcf.info+=i.first+";END="+std::to_string(dup_position+dup_length-1) +";SVLEN=" + std::to_string(dup_length);

	dup_map[vcf.chrom][vcf.pos]=vcf;

	for ( std::size_t k = vcf.pos ; k < vcf.pos + dup_length + 1 ; k ++ ){
	  masking[vcf.chrom][k]=0;
	}

	j ++;
      }
    }
  }
  std::cout << "Number of simulated DUP: " << dup_count << std::endl;

  return dup_map;
}








////////////////////////////////////////////////////////////////////////////////
//
// inversion
//

VCF_MAP GENOME::sim_inv( SVTYPE & inv, std::mt19937 & generator){
  VCF vcf;

  std::size_t pick;
  std::string selected_id;


  std::uniform_int_distribution<std::size_t> whole_genome_dist(1,total_length);


  for ( std::size_t i = 0 ; i < inv.number ; i ++ ){
    pick = whole_genome_dist(generator);
    for ( auto & j : sequence_size_sum ){
      selected_id = j.first;
      if ( pick < j.second ){
	break;
      }
    }
    inv.dist[selected_id] ++ ;
  }

  // std::cout << "Simulation of inversion" << std::endl;

  VCF_MAP inv_map;

  std::size_t inv_count=0;
  for ( auto & i : inv.dist ){
    std::size_t size = genome_info[i.first];
    std::uniform_int_distribution<std::size_t> position(1,size);
    std::uniform_int_distribution<std::size_t> length(inv.min_length,inv.max_length);
    std::size_t inv_position;
    std::size_t inv_length;
    for ( std::size_t j = 0 ; j < i.second ; ){
      inv_position=position(generator);
      inv_length=length(generator);

      int Ncount = 0 ;
      for ( std::size_t k = inv_position-1 ; k < inv_position+inv_length ; k ++ ){
	if ( sequence[i.first][k]=='N' || sequence[i.first][k]=='n' ) Ncount ++ ;
      }
      Ncount *=100;
      Ncount /=(int) inv_length ;

      bool masked=1;
      for ( std::size_t k = inv_position ; k < inv_position + inv_length + 1; k ++ ){
        masked = masked && masking[i.first][k];
      }

      if (sequence[i.first][inv_position-1]!='N' && sequence[i.first][inv_position-1+inv_length] !='N' &&
          sequence[i.first][inv_position-1]!='n' && sequence[i.first][inv_position-1+inv_length] !='n' &&
          Ncount < 10 && masked && inv_position > 1000 && inv_position < size - 1000 ){
	vcf.chrom=i.first;
	vcf.pos=inv_position;
	vcf.id="INV_"+std::to_string(inv_count++);
	vcf.ref=sequence[i.first].substr(vcf.pos-1,inv_length);
	vcf.alt=vcf.ref;
	flip_seq(vcf.alt);
	vcf.qual=".";
	vcf.filter=".";
	vcf.info="SVTYPE=INV;SVMETHOD=RealSimSS;CHR2=";
	vcf.info+=i.first+";END="+std::to_string(inv_position+inv_length) +";SVLEN=" + std::to_string(inv_length);

        inv_map[vcf.chrom][vcf.pos]=vcf;

	for ( std::size_t k = vcf.pos ; k < vcf.pos + inv_length + 1 ; k ++ ){
          masking[vcf.chrom][k]=0;
        }

	j ++;
      }
    }
  }

  std::cout << "Number of simulated INV: " << inv_count << std::endl;

  return inv_map;
}








////////////////////////////////////////////////////////////////////////////////
//
// translocation
//

VCF_MAP GENOME::sim_tra( SVTYPE & tra, std::mt19937 & generator){
  VCF vcf;

  std::size_t pick;
  std::string selected_id;


  std::uniform_int_distribution<std::size_t> whole_genome_dist(1,total_length);


  for ( std::size_t i = 0 ; i < tra.number ; i ++ ){
    pick = whole_genome_dist(generator);
    for ( auto & j : sequence_size_sum ){
      selected_id = j.first;
      if ( pick < j.second ){
	break;
      }
    }
    tra.dist[selected_id] ++ ;
  }

  // std::cout << "Simulation of translocation" << std::endl;

  VCF_MAP tra_map;

  std::size_t tra_count=0;
  for ( auto & i : tra.dist ){
    std::size_t size = genome_info[i.first];
    std::uniform_int_distribution<std::size_t> position(1,size);
    std::uniform_int_distribution<std::size_t> flip(0,1);
    std::uniform_int_distribution<std::size_t> length(tra.min_length,tra.max_length);
    std::size_t tra_position_1;
    std::size_t tra_position_2;
    std::size_t tra_length;
    std::string seq;
    bool strand;
    for ( std::size_t j = 0 ; j < i.second ; ){
      tra_position_1=position(generator);
      tra_length=length(generator);

      int Ncount = 0 ;
      for ( std::size_t k = tra_position_1-1 ; k < tra_position_1+tra_length ; k ++ ){
        if ( sequence[i.first][k]=='N' || sequence[i.first][k]=='n' ) Ncount ++ ;
      }
      Ncount *=100;
      Ncount /=(int) tra_length ;

      bool masked=1;
      for ( std::size_t k = tra_position_1 ; k < tra_position_1 + tra_length + 1; k ++ ){
        masked = masked && masking[i.first][k];
      }

      if (sequence[i.first][tra_position_1-1]!='N' && sequence[i.first][tra_position_1-1+tra_length] !='N' &&
          sequence[i.first][tra_position_1-1]!='n' && sequence[i.first][tra_position_1-1+tra_length] !='n' &&
          Ncount < 10 && masked && tra_position_1 > 1000 && tra_position_1 < size - tra_length - 1000 ){
	
	pick = whole_genome_dist(generator);
	for ( auto & j : sequence_size_sum ){
	  if ( pick < j.second ){
	    selected_id = j.first;
	  }
	  else{
	    break;
	  }
	}
	size = genome_info[selected_id];
	std::uniform_int_distribution<std::size_t> position_2(1,size);
	tra_position_2=position_2(generator);
	
	if (sequence[selected_id][tra_position_2-1]!='N' && sequence[selected_id][tra_position_2-1+tra_length] !='N' &&
	    sequence[selected_id][tra_position_2-1]!='n' && sequence[selected_id][tra_position_2-1+tra_length] !='n' &&
	    masking[selected_id][tra_position_2] && masking[selected_id][tra_position_2+1] &&
	    tra_position_2 > 1000 && tra_position_2 < size - 1000 ){

	  seq=sequence[i.first].substr(tra_position_1,tra_length);
	  strand = flip(generator) ;

	  vcf.chrom=i.first;
	  vcf.pos=tra_position_1;
	  vcf.id="TRA_"+std::to_string(tra_count);
	  vcf.ref=sequence[vcf.chrom].substr(vcf.pos-1,tra_length+1);
	  vcf.alt=sequence[vcf.chrom][vcf.pos-1];
	  vcf.qual=".";
	  vcf.filter=".";
	  vcf.info="SVTYPE=TRA;SVMETHOD=RealSimSS;CHR2=";
	  vcf.info+=selected_id+";END="+std::to_string(tra_position_2) +";SVLEN=" + std::to_string(tra_length);
	  
	  tra_map[vcf.chrom][vcf.pos]=vcf;

	  if (strand){
	    flip_seq(seq);
	  }
	  
	  vcf.chrom=selected_id;
	  vcf.pos=tra_position_2;
	  vcf.id="TRA_"+std::to_string(tra_count);
	  vcf.ref=sequence[i.first][tra_position_2-1];
	  vcf.alt=vcf.ref;
	  vcf.alt+=seq;
	  vcf.qual=".";
	  vcf.filter=".";
	  vcf.info="SVTYPE=TRA;SVMETHOD=RealSimSS;CHR2=";
	  vcf.info+=i.first+";END="+std::to_string(tra_position_1) +";SVLEN=" + std::to_string(tra_length);

	  tra_map[vcf.chrom][vcf.pos]=vcf;
	   
	  for ( std::size_t k = tra_position_1 ; k < tra_position_1 + tra_length + 1 ; k ++ ){
	    masking[i.first][k]=0;
	  }
	  
	  masking[selected_id][tra_position_2]=0;

	  j ++ ;
	  tra_count ++;
	}
      }
    }
  }

  std::cout << "Number of simulated TRA: " << tra_count << std::endl;

  return tra_map;
}



