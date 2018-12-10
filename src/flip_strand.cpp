#ifndef FLIP_STRAND
#define FLIP_STRAND

#include <string>
#include <string.h>
#include <algorithm>


char flip_base(char a){
  switch ( a ){
  case 'A' : return 'T' ;
  case 'T' : return 'A' ;
  case 'U' : return 'A' ;
    
  case 'C' : return 'G' ;
  case 'G' : return 'C' ;
    
  case 'R' : return 'Y' ; //A, G -->  T(U), C
  case 'Y' : return 'R' ; //C, T(U) ---> G, A

  case 'K' : return 'M' ; //G, T(U) ---> C, A
  case 'M' : return 'K' ; //C, A ---> G, T(U)

  case 'S' : return 'S' ; //G, C ---> G, C
  case 'W' : return 'W' ; //A, T(U) ---> A, T(U)
    
  case 'B' : return 'V' ; //not A --> not T(U)
  case 'V' : return 'B' ; //not T(U) --> not A

  case 'D' : return 'H' ; //not C --> not G
  case 'H' : return 'D' ; //not G --> not C

  case 'N' : return 'N' ;

  default : return 'i';
  }
    
};



void flip_base(std::string::iterator & a){
  switch ( *a ){
  case 'A' : *a = 'T' ; break ;
  case 'T' : *a = 'A' ; break ;
  case 'U' : *a = 'A' ; break ;
    
  case 'C' : *a = 'G' ; break ;
  case 'G' : *a = 'C' ; break ;
    
  case 'a' : *a = 't' ; break ;
  case 't' : *a = 'a' ; break ;
  case 'u' : *a = 'a' ; break ;
    
  case 'c' : *a = 'g' ; break ;
  case 'g' : *a = 'c' ; break ;
    
  case 'R' : *a = 'Y' ; break ;  //A, G -->  C, T(U)
  case 'Y' : *a = 'R' ; break ;  //C, T(U) ---> G, A

  case 'K' : *a = 'M' ; break ;  //G, T(U) ---> C, A
  case 'M' : *a = 'K' ; break ;  //C, A ---> G, T(U)

  case 'S' : *a = 'S' ; break ;  //G, C ---> G, C
  case 'W' : *a = 'W' ; break ;  //A, T(U) ---> A, T(U)
    
  case 'B' : *a = 'V' ; break ;  //not A --> not T(U)
  case 'V' : *a = 'B' ; break ;  //not T(U) --> not A

  case 'D' : *a = 'H' ; break ;  //not C --> not G
  case 'H' : *a = 'D' ; break ;  //not G --> not C

  case 'N' : *a = 'N' ; break ;

  default :  *a = 'N' ; break ;
  }
};


    //flip_base(i) ;

void flip_seq(std::string & seq){
  reverse(seq.begin(),seq.end());
  for ( std::string::iterator a = seq.begin() ; a != seq.end() ; a ++ ){
    switch ( *a ){
    case 'A' : *a = 'T' ; break ;
    case 'T' : *a = 'A' ; break ;
    case 'U' : *a = 'A' ; break ;
      
    case 'C' : *a = 'G' ; break ;
    case 'G' : *a = 'C' ; break ;
      
    case 'a' : *a = 't' ; break ;
    case 't' : *a = 'a' ; break ;
    case 'u' : *a = 'a' ; break ;
      
    case 'c' : *a = 'g' ; break ;
    case 'g' : *a = 'c' ; break ;
    
    case 'R' : *a = 'Y' ; break ;  //A, G -->  C, T(U)
    case 'Y' : *a = 'R' ; break ;  //C, T(U) ---> G, A
      
    case 'K' : *a = 'M' ; break ;  //G, T(U) ---> C, A
    case 'M' : *a = 'K' ; break ;  //C, A ---> G, T(U)
      
    case 'S' : *a = 'S' ; break ;  //G, C ---> G, C
    case 'W' : *a = 'W' ; break ;  //A, T(U) --->  A, T(U)
      
    case 'B' : *a = 'V' ; break ;  //not A --> not T(U)
    case 'V' : *a = 'B' ; break ;  //not T(U) --> not A
      
    case 'D' : *a = 'H' ; break ;  //not C --> not G
    case 'H' : *a = 'D' ; break ;  //not G --> not C
      
    case 'N' : *a = 'N' ; break ;
      
    default :  *a = 'N' ; break ;
    }
  }
}



#endif
