#include "../inst/include/floodfill.hpp"


struct Pixel {
  int x;
  int y;
  Pixel(int x_, int y_): x(x_), y(y_) {};
};

struct Image {
  double *rPic;
  double dtol;
  double dInitCol;
  int h;
  int w;

  Image(double *rPic_, double dtol_, int h_, int w_): dtol(dtol_), h(h_), w(w_), dInitCol(0.0) {
    rPic = rPic_;
  };

  void setInitCol(double col){
    dInitCol = col;
  }

  double getInitCol(){
    return dInitCol;
  }

  double getPixelColor(Pixel pp){
    if(pp.x >= 0 && pp.x < w  && pp.y >= 0 && pp.y < h)
      return rPic[pp.x * h + pp.y];
    return -1.0;
  }
  
  bool isReplaceablePixelColor(Pixel pp){
    if(pp.x >= 0 && pp.x < w  && pp.y >= 0 && pp.y < h)
      return fabs(rPic[pp.x * h + pp.y]-getInitCol()) < dtol && rPic[pp.x * h + pp.y] > -1.0;
    return false;
  }

  void setPixelColor(Pixel pp){
    if(pp.x >= 0 && pp.x < w  && pp.y >= 0 && pp.y < h)
      rPic[pp.x * h + pp.y] = -1.0;
  }

};
    

// From doi: 10.1007/s11554-017-0732-1
extern "C" {
  SEXP scanlineFill(SEXP pic, SEXP tol, SEXP x0, SEXP y0, SEXP setInitCol){
    SEXP picOut;
    int h = Rf_nrows(pic);
    int w = Rf_ncols(pic);
    int xx = INTEGER(x0)[0];
    int yy = INTEGER(y0)[0];
    bool sic = LOGICAL(setInitCol)[0];
    PROTECT(picOut = Rf_allocMatrix(REALSXP,h,w));
    Rf_copyMatrix(picOut,pic,FALSE);
    Image img(REAL(picOut), REAL(tol)[0], h, w);
    std::vector<Pixel> pixels;
    Pixel start(xx,yy);
    if(sic)
      img.setInitCol(img.getPixelColor(start));
    pixels.push_back(start);
    while(!pixels.empty()) {
      Pixel pp = pixels.back();
      pixels.pop_back();
      int x1 = pp.x;      
      while(x1 >= 0 && img.isReplaceablePixelColor(Pixel(x1,pp.y))){
	x1--;
      }
      x1++;
      bool upFlag = false;
      bool downFlag = false;

      while(x1 < w && img.isReplaceablePixelColor(Pixel(x1,pp.y))){
	img.setPixelColor(Pixel(x1,pp.y));
	if(!upFlag && pp.y > 0){
	  if(img.isReplaceablePixelColor(Pixel(x1,pp.y-1))){
	    pixels.push_back(Pixel(x1,pp.y-1));
	    upFlag = true;
	  }else{
	    upFlag = false;
	  }
	}else if(upFlag && pp.y > 0 && x1 > 0){
	  if(!img.isReplaceablePixelColor(Pixel(x1-1,pp.y-1)) &&
	     img.isReplaceablePixelColor(Pixel(x1,pp.y-1))){
	    pixels.push_back(Pixel(x1,pp.y-1));	    
	  }
	}
	if(!downFlag && pp.y < h-1){
	  if(img.isReplaceablePixelColor(Pixel(x1,pp.y+1))){
	    pixels.push_back(Pixel(x1,pp.y+1));
	    downFlag = true;
	  }else{
	      downFlag = false;
	  }
	}else if(downFlag && pp.y < h-1 && x1 > 0){
	  if(!img.isReplaceablePixelColor(Pixel(x1-1,pp.y+1)) &&
	     img.isReplaceablePixelColor(Pixel(x1,pp.y+1))){
	    pixels.push_back(Pixel(x1,pp.y+1));	    
	  }
	}
	x1++;
      }
    }
    UNPROTECT(1);
    return picOut;
  }
}
