#pragma once

#include <math.h>
#include <vector>
#include <stdlib.h>
#include <malloc.h>

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/core.hpp"

using namespace std;


#define PLATOS_POU_THELO        -1
#define IPSOS_POU_THELO         -1


struct cc
{
        int xmin,xmax,ymin,ymax,label,xg,yg,points;
};

class BinaryImage
        {
        public:
                //Constructors
                BinaryImage();

                //Copy Constructor
                //BinaryImage(const BinaryImage &Eikona);
                //Destructors
                ~BinaryImage();
                //Alles Sinartiseis
                void Clear();
                int Calculate_RLSA(int hthr,int vthr); //Calculate RLSA
                int Connected_Component_Labeling();
                int Skeleton();
                int Normalize(int width,int height);
                int NormalizeBGAT(int norm_width, int norm_height);
                //Getters
                int Get_Average_Character_Height();
                int Get_Image_Width();
                int Get_Image_Height();
                int Get_Labelmap(int x,int y);
                unsigned char Get_Bitmap(int x,int y);
                vector <cc> Get_Connected_Components();
                unsigned char Get_Bitmap_RLSA(int x,int y);
                unsigned char Get_Bitmap_Skeleton(int x,int y);
                unsigned char Get_Bitmap_Normalized(int x,int y);
                //Setters
                int Set_Bitmap(int x,int y,unsigned char value);
                int Create_Bitmap(int Ix,int Iy);
								unsigned char **bitmap;
								int **labelmap;
								unsigned char **bitmap_Skeleton;
								vector <cc> Connected_Components;
        private:
                unsigned char **bitmap_RLSA;
                unsigned char **bitmap_Normalized;
                int Image_Width;
                int Image_Height;
                int Normalized_Height;
                int Normalized_Width;
                int Average_Character_Height;
                unsigned char ypo(int x,int y,int Iy,int Dx,unsigned char *IMAGE);
                vector <cc> connectedcomponents(unsigned char **bitmap,int **labelmap,int Ix,int Iy,bool *MemStatus);
                void Tracer_8(int *cy, int *cx, int *tracingdirection,int **labelmap,unsigned char **bitmap);
                void ContourTracing_8(int cy, int cx, int labelindex, int tracingdirection,int **labelmap,unsigned char**bitmap);
                int ConnectedComponentLabeling_8(unsigned char **bitmap,int **labelmap,int *plithos,int Ix,int Iy);
                void Calculate_Average_Character_Height();
                bool Ench(unsigned char *IMA);
                int ZNZT(int x,int y,unsigned char *IMA);
                int NZN(int x,int y,unsigned char *IMA);
        };
