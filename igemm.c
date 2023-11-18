#include <assert.h>

static inline int madd(int x, int a, int b) {
  return x + (a*b);
}

void igemm(int *Cg, int *Ag, int *Bg, int n, int m, int _k, int tid, int nthr) {
  static const int VL = 1;
  static const int II_BLK = 4;
  static const int JJ_VL_BLK = 6;
  static const int JJ_BLK = VL * JJ_VL_BLK;

  static const int CACHE_BLK = 120;
  assert(n % CACHE_BLK == 0);
  assert(n % II_BLK == 0);

  int vC0_0, vC0_1, vC0_2, vC0_3, vC0_4,  vC0_5;
  int vC1_0, vC1_1, vC1_2, vC1_3, vC1_4,  vC1_5;
  int vC2_0, vC2_1, vC2_2, vC2_3, vC2_4,  vC2_5;
  int vC3_0, vC3_1, vC3_2, vC3_3, vC3_4,  vC3_5;  

  int vA, vB0, vB1, vB2, vB3, vB4, vB5;

  int slice = n/nthr;
  int start = slice*tid;
  int stop = (tid==(n-1)) ? n : slice*(tid+1);
  
  //the whole matrix
  for(int ii = start; ii < stop; ii += CACHE_BLK) {
    for(int jj = 0; jj < m; jj+=CACHE_BLK) {
      for(int kk = 0; kk < _k; kk += CACHE_BLK) {      
	for(int i = ii; i < (ii+CACHE_BLK); i+=II_BLK) {
	  for(int j = jj; j < (jj+CACHE_BLK); j+=JJ_BLK) {
	    
	    vC0_0 = *( &Cg[(i+0)*n+j+(0*VL)] );
	    vC0_1 = *( &Cg[(i+0)*n+j+(1*VL)] );
	    vC0_2 = *( &Cg[(i+0)*n+j+(2*VL)] );
	    vC0_3 = *( &Cg[(i+0)*n+j+(3*VL)] );
	    vC0_4 = *( &Cg[(i+0)*n+j+(4*VL)] );
	    vC0_5 = *( &Cg[(i+0)*n+j+(5*VL)] );                  
	    
	    vC1_0 = *( &Cg[(i+1)*n+j+(0*VL)] );
	    vC1_1 = *( &Cg[(i+1)*n+j+(1*VL)] );
	    vC1_2 = *( &Cg[(i+1)*n+j+(2*VL)] );
	    vC1_3 = *( &Cg[(i+1)*n+j+(3*VL)] );
	    vC1_4 = *( &Cg[(i+1)*n+j+(4*VL)] );
	    vC1_5 = *( &Cg[(i+1)*n+j+(5*VL)] );                  
	    
	    vC2_0 = *( &Cg[(i+2)*n+j+(0*VL)] );
	    vC2_1 = *( &Cg[(i+2)*n+j+(1*VL)] );
	    vC2_2 = *( &Cg[(i+2)*n+j+(2*VL)] );
	    vC2_3 = *( &Cg[(i+2)*n+j+(3*VL)] );
	    vC2_4 = *( &Cg[(i+2)*n+j+(4*VL)] );
	    vC2_5 = *( &Cg[(i+2)*n+j+(5*VL)] );                  
	    
	    vC3_0 = *( &Cg[(i+3)*n+j+(0*VL)] );
	    vC3_1 = *( &Cg[(i+3)*n+j+(1*VL)] );
	    vC3_2 = *( &Cg[(i+3)*n+j+(2*VL)] );
	    vC3_3 = *( &Cg[(i+3)*n+j+(3*VL)] );
	    vC3_4 = *( &Cg[(i+3)*n+j+(4*VL)] );
	    vC3_5 = *( &Cg[(i+3)*n+j+(5*VL)] );                  

	    for(int k = kk; k < (kk+CACHE_BLK); k++) {
	      vB0 = (Bg[(k)*m+j+(0*VL)]);
	      vB1 = (Bg[(k)*m+j+(1*VL)]);
	      vB2 = (Bg[(k)*m+j+(2*VL)]);
	      vB3 = (Bg[(k)*m+j+(3*VL)]);
	      vB4 = (Bg[(k)*m+j+(4*VL)]);
	      vB5 = (Bg[(k)*m+j+(5*VL)]);

	      //i = 0
	      vA = (Ag[(i+0)*_k+k]);
	      vC0_0 = madd(vC0_0, vA, vB0);
	      vC0_1 = madd(vC0_1, vA, vB1);
	      vC0_2 = madd(vC0_2, vA, vB2);
	      vC0_3 = madd(vC0_3, vA, vB3);
	      vC0_4 = madd(vC0_4, vA, vB4);
	      vC0_5 = madd(vC0_5, vA, vB5);	  	  
	      
	      //i = 1
	      vA = (Ag[(i+1)*_k+k]);
	      vC1_0 = madd(vC1_0, vA, vB0);
	      vC1_1 = madd(vC1_1, vA, vB1);
	      vC1_2 = madd(vC1_2, vA, vB2);
	      vC1_3 = madd(vC1_3, vA, vB3);
	      vC1_4 = madd(vC1_4, vA, vB4);
	      vC1_5 = madd(vC1_5, vA, vB5);	  	  
	      
	      //i = 2
	      vA = (Ag[(i+2)*_k+k]);
	      vC2_0 = madd(vC2_0, vA, vB0);
	      vC2_1 = madd(vC2_1, vA, vB1);
	      vC2_2 = madd(vC2_2, vA, vB2);
	      vC2_3 = madd(vC2_3, vA, vB3);
	      vC2_4 = madd(vC2_4, vA, vB4);
	      vC2_5 = madd(vC2_5, vA, vB5);	  	  
	      
	      //i = 3
	      vA = (Ag[(i+3)*_k+k]);
	      vC3_0 = madd(vC3_0, vA, vB0);
	      vC3_1 = madd(vC3_1, vA, vB1);
	      vC3_2 = madd(vC3_2, vA, vB2);
	      vC3_3 = madd(vC3_3, vA, vB3);
	      vC3_4 = madd(vC3_4, vA, vB4);
	      vC3_5 = madd(vC3_5, vA, vB5);	  	  
	    }
	    ( Cg[(i+0)*n+j+(0*VL)]= vC0_0 );
	    ( Cg[(i+0)*n+j+(1*VL)]= vC0_1 );
	    ( Cg[(i+0)*n+j+(2*VL)]= vC0_2 );
	    ( Cg[(i+0)*n+j+(3*VL)]= vC0_3 );
	    ( Cg[(i+0)*n+j+(4*VL)]= vC0_4 );
	    ( Cg[(i+0)*n+j+(5*VL)]= vC0_5 );    
	    
	    ( Cg[(i+1)*n+j+(0*VL)]= vC1_0 );
	    ( Cg[(i+1)*n+j+(1*VL)]= vC1_1 );
	    ( Cg[(i+1)*n+j+(2*VL)]= vC1_2 );
	    ( Cg[(i+1)*n+j+(3*VL)]= vC1_3 );
	    ( Cg[(i+1)*n+j+(4*VL)]= vC1_4 );
	    ( Cg[(i+1)*n+j+(5*VL)]= vC1_5 );    
	    
	    ( Cg[(i+2)*n+j+(0*VL)]= vC2_0 );
	    ( Cg[(i+2)*n+j+(1*VL)]= vC2_1 );
	    ( Cg[(i+2)*n+j+(2*VL)]= vC2_2 );
	    ( Cg[(i+2)*n+j+(3*VL)]= vC2_3 );
	    ( Cg[(i+2)*n+j+(4*VL)]= vC2_4 );
	    ( Cg[(i+2)*n+j+(5*VL)]= vC2_5 );    
	    
	    ( Cg[(i+3)*n+j+(0*VL)]= vC3_0 );
	    ( Cg[(i+3)*n+j+(1*VL)]= vC3_1 );
	    ( Cg[(i+3)*n+j+(2*VL)]= vC3_2 );
	    ( Cg[(i+3)*n+j+(3*VL)]= vC3_3 );
	    ( Cg[(i+3)*n+j+(4*VL)]= vC3_4 );
	    ( Cg[(i+3)*n+j+(5*VL)]= vC3_5 );    
	  }
	}
      }
    }
  }
}
