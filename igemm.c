#include <assert.h>

static inline int madd(int x, int a, int b) {
  return x + (a*b);
}

void igemm(int *Cg, int *Ag, int *Bg, int n, int m, int _k) {
  static const int VL = 1;
  static const int II_BLK = 2;
  static const int JJ_VL_BLK = 4;
  static const int JJ_BLK = VL * JJ_VL_BLK;

  static const int CACHE_BLK = 120;
  assert(n % CACHE_BLK == 0);
  assert(n % II_BLK == 0);

  int vC0_0, vC0_1, vC0_2, vC0_3;
  int vC1_0, vC1_1, vC1_2, vC1_3;
  int vC2_0, vC2_1, vC2_2, vC2_3;
  int vC3_0, vC3_1, vC3_2, vC3_3;

  int vA, vB0, vB1, vB2, vB3;

  
  //the whole matrix
  for(int ii = 0; ii < n; ii += CACHE_BLK) {
    for(int jj = 0; jj < m; jj+=CACHE_BLK) {
      for(int kk = 0; kk < _k; kk += CACHE_BLK) {      
	for(int i = ii; i < (ii+CACHE_BLK); i+=II_BLK) {
	  for(int j = jj; j < (jj+CACHE_BLK); j+=JJ_BLK) {
	    
	    vC0_0 = *( &Cg[(i+0)*n+j+(0*VL)] );
	    vC0_1 = *( &Cg[(i+0)*n+j+(1*VL)] );
	    vC0_2 = *( &Cg[(i+0)*n+j+(2*VL)] );
	    vC0_3 = *( &Cg[(i+0)*n+j+(3*VL)] );
	    
	    vC1_0 = *( &Cg[(i+1)*n+j+(0*VL)] );
	    vC1_1 = *( &Cg[(i+1)*n+j+(1*VL)] );
	    vC1_2 = *( &Cg[(i+1)*n+j+(2*VL)] );
	    vC1_3 = *( &Cg[(i+1)*n+j+(3*VL)] );
	    
	    for(int k = kk; k < (kk+CACHE_BLK); k++) {
	      vB0 = (Bg[(k)*m+j+(0*VL)]);
	      vB1 = (Bg[(k)*m+j+(1*VL)]);
	      vB2 = (Bg[(k)*m+j+(2*VL)]);
	      vB3 = (Bg[(k)*m+j+(3*VL)]);

	      //i = 0
	      vA = (Ag[(i+0)*_k+k]);
	      vC0_0 = madd(vC0_0, vA, vB0);
	      vC0_1 = madd(vC0_1, vA, vB1);
	      vC0_2 = madd(vC0_2, vA, vB2);
	      vC0_3 = madd(vC0_3, vA, vB3);
	      
	      //i = 1
	      vA = (Ag[(i+1)*_k+k]);
	      vC1_0 = madd(vC1_0, vA, vB0);
	      vC1_1 = madd(vC1_1, vA, vB1);
	      vC1_2 = madd(vC1_2, vA, vB2);
	      vC1_3 = madd(vC1_3, vA, vB3);

	    }
	    
	    Cg[(i+0)*n+j+(0*VL)]= vC0_0;
	    Cg[(i+0)*n+j+(1*VL)]= vC0_1;
	    Cg[(i+0)*n+j+(2*VL)]= vC0_2;
	    Cg[(i+0)*n+j+(3*VL)]= vC0_3;
	    
	    Cg[(i+1)*n+j+(0*VL)]= vC1_0;
	    Cg[(i+1)*n+j+(1*VL)]= vC1_1;
	    Cg[(i+1)*n+j+(2*VL)]= vC1_2;
	    Cg[(i+1)*n+j+(3*VL)]= vC1_3;
	  }
	}
      }
    }
  }
}
