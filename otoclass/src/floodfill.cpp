#include "../inst/include/floodfill.hpp"


struct Pixel {
  int x;
  int y;
  Pixel(int x_, int y_): x(x_), y(y_) {};
};

struct Image {
  double *rPic;
  double dtol;
  int h;
  int w;

  Image(double *rPic_, double dtol_, int h_, int w_): dtol(dtol_), h(h_), w(w_) {
    rPic = rPic_;
  };

  bool isReplaceablePixelColor(Pixel pp){
    if(pp.x >= 0 && pp.x < w  && pp.y >= 0 && pp.y < h)
      return rPic[pp.x * h + pp.y] < dtol && rPic[pp.x * h + pp.y] > -1.0;
    return false;
  }

  void setPixelColor(Pixel pp){
    if(pp.x >= 0 && pp.x < w  && pp.y >= 0 && pp.y < h)
      rPic[pp.x * h + pp.y] = -1.0;
  }

};
    

// From doi: 10.1007/s11554-017-0732-1
extern "C" {
  SEXP scanlineFill(SEXP pic, SEXP tol){
    SEXP picOut;
    int h = Rf_nrows(pic);
    int w = Rf_ncols(pic);
    PROTECT(picOut = Rf_allocMatrix(REALSXP,h,w));
    Rf_copyMatrix(picOut,pic,FALSE);
    Image img(REAL(picOut), REAL(tol)[0], h, w);
    std::vector<Pixel> pixels;
    pixels.push_back(Pixel(0, 0));
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
