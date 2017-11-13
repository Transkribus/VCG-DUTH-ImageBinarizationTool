#pragma once

#include <math.h>
#include <vector>
#include <stdlib.h>
#include <malloc.h>

using namespace std;


struct ALLT_struct {
  int  x,y;
};


class GrayScaleImage
        {
        public:
                //Constructors
                GrayScaleImage();
                //Copy Constructor
                GrayScaleImage(const GrayScaleImage &Eikona);
                //Destructors
                ~GrayScaleImage();
                //Alles Sinartiseis
                //void Load(AnsiString FileName);
                void Clear();
                int Get_Image_Width();
                int Get_Image_Height();
                //AnsiString Get_ImageFileName();
                unsigned char Get_Bitmap(int x,int y);
                int Binarize(bool oldOtsu,double Kaugment);
                int Binarize(std::string FileName,bool oldOtsu,double Kaugment);
								int Create_KBitmap(int Ix,int Iy);
                BinaryImage ImageB;
				unsigned char *Kbitmap;
                                
        private:
                std::string ImageFileName;
                unsigned char **bitmap;
                int Image_Width;
                int Image_Height;
                //void StoreImage(TImagXpress7_ *ImagXpress7_1);
                //unsigned char* CreateIM(HANDLE Dib,int Ix,int mIy);
                //unsigned char ypo(int x,int y,int Iy,int Dx,unsigned char *IMAGE);
                int Otsu(bool doubleChar,double Kaugment);
                int Sauvola();
                int max(int a,int b);
                int min(int a,int b);
                int GPP(unsigned char *in,int Ix,int Iy);
                int ALLT_Ypologismos_Thres(unsigned char *EIK,int SW,double Glob_a,int x_t, int y_t,int wind,int Ix,int Iy);
                int ALLT_Ypologismos_Pj(unsigned char *EIK,int SW,int x_j, int y_j,int wind,int Ix,int Iy);
                int ALLT(unsigned char *EIK,int SW,int Ix,int Iy);
                int ALLT_SW(unsigned char *IM,int Ix,int Iy);
                
        };


