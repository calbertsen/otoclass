#include "../inst/include/efd.hpp"

template<class Type>
Matrix<Type, Dynamic, Dynamic> efd2coord_work(Matrix<Type, Dynamic, Dynamic> efd, int N){
  Matrix<Type, Dynamic, Dynamic> coord(2,N);
  coord.setZero();

  for(int tt = 0; tt < N; ++tt){
    Type t = (Type)tt / (Type)(N-1);
    for(int i = 0; i < efd.cols(); ++i){
      coord(0,tt) += efd(0,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(0,tt) += efd(1,i) * sin(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(2,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(3,i) * sin(2.0 * Type(i+1) * M_PI * t);
    }
  }
  return coord;
}


extern "C" {
  SEXP efd2coordSEXP(SEXP efd, SEXP N, SEXP A0, SEXP C0){
    if(!Rf_isMatrix(efd) || Rf_ncols(efd) != 4)
      Rf_error("efd must be a matrix with four columns");
    if(!Rf_isInteger(N) || Rf_length(N) != 1)
      Rf_error("N must be a length one integer vector");
    if(!Rf_isNumeric(A0) || Rf_length(A0) != 1)
      Rf_error("A0 must be a length one numeric vector");
    if(!Rf_isNumeric(C0) || Rf_length(C0) != 1)
      Rf_error("C0 must be a length one numeric vector");

    int NN = INTEGER(N)[0];
    SEXP coord;
    PROTECT(coord = Rf_allocMatrix(REALSXP, NN, 2));
    double* rcoord = REAL(coord);
    double* refd = REAL(efd);
    int NE = Rf_nrows(efd);
    
    for(int i = 0; i < NN; ++i){
      rcoord[i + NN * 0] = REAL(A0)[0];
      rcoord[i + NN * 1] = REAL(C0)[0];
      for(int j = 0; j < NE; ++j){
	double tmp = 2.0 * ((double)j + 1.0) * M_PI * (double)i / (double)(NN-1);
	rcoord[i + NN * 0] += refd[j + NE * 0] * cos(tmp);
	rcoord[i + NN * 0] += refd[j + NE * 1] * sin(tmp);
	rcoord[i + NN * 1] += refd[j + NE * 2] * cos(tmp);
	rcoord[i + NN * 1] += refd[j + NE * 3] * sin(tmp);
      }
    }
    UNPROTECT(1);
    return coord;
  }

  
  SEXP efd(SEXP dat, SEXP N, SEXP normalize){
    if(!Rf_isMatrix(dat) || Rf_ncols(dat) != 2)
      Rf_error("dat must be a matrix with two columns");
    if(!Rf_isInteger(N) || Rf_length(N) != 1)
      Rf_error("N must be a length one integer vector");
    if(!Rf_isLogical(normalize) || Rf_length(normalize) != 1)
      Rf_error("normalize must be a length one logical vector");

  int NN = INTEGER(N)[0];
  int nr = Rf_nrows(dat);
  //int nc = Rf_ncols(dat);
  int norma = LOGICAL(normalize)[0];

  SEXP A, B, C, D, A0, C0;
  A = PROTECT(Rf_allocVector(REALSXP,NN));
  B = PROTECT(Rf_allocVector(REALSXP,NN));
  C = PROTECT(Rf_allocVector(REALSXP,NN));
  D = PROTECT(Rf_allocVector(REALSXP,NN));
  A0 = PROTECT(Rf_allocVector(REALSXP,1));
  C0 = PROTECT(Rf_allocVector(REALSXP,1));

  double *rA = REAL(A);
  double *rB = REAL(B);
  double *rC = REAL(C);
  double *rD = REAL(D);
  double *rDat = REAL(dat);
  
  REAL(A0)[0] = 0.0;
  REAL(C0)[0] = 0.0;
  for(int i = 0; i < NN; ++i) {
    rA[i] = 0.0;
    rB[i] = 0.0;
    rC[i] = 0.0;
    rD[i] = 0.0;
  }

  double M_PI_SQ_2 = M_PI * M_PI * 2.0;
  double t = 0.0;
  double tm1 = 0.0;
  double T = 0.0;
  for(int j = 1; j < nr; ++j){    
    // double dx = pd[j + nr * 0] - pd[j-1 + nr * 0];
    double dx = rDat[j + nr * 0] - rDat[j-1 + nr * 0];
    //double dy = pd[j + nr * 1] - pd[j-1 + nr * 1];
    double dy = rDat[j + nr * 1] - rDat[j-1 + nr * 1];
    double dt = sqrt(dx*dx + dy*dy);
    T += dt;
  }
  
  for(int j = 1; j < nr; ++j){
    // [r+nr*c] r is row c is col
    // double dx = pd[j + nr * 0] - pd[j-1 + nr * 0];
    double dx = rDat[j + nr * 0] - rDat[j-1 + nr * 0];
    //double dy = pd[j + nr * 1] - pd[j-1 + nr * 1];
    double dy = rDat[j + nr * 1] - rDat[j-1 + nr * 1];
    double dt = sqrt(dx*dx + dy*dy);
    //if(dt > 0){
    t += dt;
    // REAL(A0)[0] += 0.5 * dx / dt * (t*t - tm1*tm1) + (pd[j-1 + nr * 0] - dx / dt * tm1) * dt;
    // REAL(C0)[0] += 0.5 * dy / dt * (t*t - tm1*tm1) + (pd[j-1 + nr * 1] - dy / dt * tm1) * dt;
    REAL(A0)[0] += 0.5 * dx / dt * (t*t - tm1*tm1) + (rDat[j-1 + nr * 0] - dx / dt * tm1) * dt;
    REAL(C0)[0] += 0.5 * dy / dt * (t*t - tm1*tm1) + (rDat[j-1 + nr * 1] - dy / dt * tm1) * dt;
    for(int i = 1; i <= NN; ++i){
      double f = T / ((double)i*(double)i * M_PI_SQ_2);
      double di = i;
      double pi2i1 = 2.0 * di * M_PI * t / T;
      double pi2i2 = 2.0 * di * M_PI * tm1 / T;
      rA[i-1] += f * dx / dt * ( cos( pi2i1 ) - cos( pi2i2 ) );
      rB[i-1] += f * dx / dt * ( sin( pi2i1 ) - sin( pi2i2 ) );
      rC[i-1] += f * dy / dt * ( cos( pi2i1 ) - cos( pi2i2 ) );
      rD[i-1] += f * dy / dt * ( sin( pi2i1 ) - sin( pi2i2 ) );
    }
    tm1 += dt;
    //}
  }
  REAL(A0)[0] /= T;
  REAL(C0)[0] /= T;


  if(norma){
    double theta = 0.5 * atan2( 2.0 * ( rA[0] * rB[0] + rC[0] * rD[0] ),
				(rA[0]*rA[0] - rB[0]*rB[0] + rC[0]*rC[0] - rD[0]*rD[0]));
    double as = rA[0] * cos(theta) + rB[0] * sin(theta);
    double cs = rC[0] * cos(theta) + rD[0] * sin(theta);
    double scale = 1.0 / sqrt(as*as + cs*cs);
    double phi = atan2(cs,as);

    double cp = cos(phi);
    double sp = sin(phi);
    
    for(int i = 1; i <= NN; ++i){
      double di = i;
      double cit = cos(theta * di);
      double sit = sin(theta * di);
      double Ai = rA[i-1];
      double Bi = rB[i-1];
      double Ci = rC[i-1];
      double Di = rD[i-1];
      rA[i-1] = scale * ((cp * Ai + sp * Ci) * cit + (cp * Bi + sp * Di) * sit);
      rB[i-1] = scale * ((cp * Ai + sp * Ci) * (-sit) + (cp * Bi + sp * Di) * cit);
      rC[i-1] = scale * ((-sp * Ai + cp * Ci) * cit + (-sp * Bi + cp * Di) * sit);
      rD[i-1] = scale * ((-sp * Ai + cp * Ci) * (-sit) + (-sp * Bi + cp * Di) * cit);
    }

    REAL(A0)[0] = 0.0;
    REAL(C0)[0] = 0.0;
  }

  SEXP res = PROTECT(Rf_allocVector(VECSXP,6));
  SET_VECTOR_ELT(res, 0, A);
  SET_VECTOR_ELT(res, 1, B);
  SET_VECTOR_ELT(res, 2, C);
  SET_VECTOR_ELT(res, 3, D);
  SET_VECTOR_ELT(res, 4, A0);
  SET_VECTOR_ELT(res, 5, C0);

  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nms, 0, Rf_mkChar("A"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("B"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("C"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("D"));
  SET_STRING_ELT(nms, 4, Rf_mkChar("A0"));
  SET_STRING_ELT(nms, 5, Rf_mkChar("C0"));

  Rf_setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(8);
  return res;

}

}
