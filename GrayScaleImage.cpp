#include "BinaryImage.h"
#include "GrayScaleImage.h"
#include <assert.h>


//YES
GrayScaleImage::GrayScaleImage()
{
        bitmap=NULL;
		Kbitmap=NULL;
        Image_Height=0;                 
        Image_Width=0;  
        ImageFileName="";

}


GrayScaleImage::GrayScaleImage(const GrayScaleImage &Eikona)
{
   Image_Width=Eikona.Image_Width;
   Image_Height=Eikona.Image_Height;
   if (Eikona.bitmap!=NULL)
   {
   bitmap=new (nothrow) unsigned char* [Image_Height];
   for (int y=0;y<Image_Height;y++)
      {
       bitmap[y]=new (nothrow) unsigned char [Image_Width];
       for (int x=0;x<Image_Width;x++)
          bitmap[y][x]=Eikona.bitmap[y][x];
       }
   }
   ImageFileName=Eikona.ImageFileName;
}

//YES
GrayScaleImage::~GrayScaleImage()
{
    if (bitmap!=NULL)
    {
      for (int y=0;y<Image_Height;y++)
        delete bitmap[y];
      delete bitmap;
      bitmap=NULL;
    }
    Image_Height=0;
    Image_Width=0;
    ImageFileName="";
		if(Kbitmap!=NULL)
		{
			free(Kbitmap);
			Kbitmap=NULL;
		}
    
}

int GrayScaleImage::Create_KBitmap(int Ix,int Iy)
{
  //Destroying old data after multiple calls of Load
  this->Clear();
  //returns 0 if success
  //-1 if Ix and or Iy <0
  //-2 if new returns error
  if ( (Ix>0) && (Iy>0))
  {
     Kbitmap=(unsigned char *)calloc(Ix*Iy,sizeof(unsigned char));
     if (Kbitmap==NULL)
        return -2;

     Image_Height=Iy;
     Image_Width=Ix;
     return 0;
  }
  else
     return -1;
}

/*void GrayScaleImage::StoreImage(TImagXpress7_ *ImagXpress7_1)
{
//To IM einai mia global metavliti
//i opoia den xreiazetai na dilothei gia na epireastei apo ti sinartisi
 unsigned char *IM;
 IM=NULL;
 if (ImagXpress7_1->IBPP==1)
  {
    if (IM!=NULL)
    {
    free(IM);
    IM=NULL;
    }
    IM = CreateIM( (HANDLE)ImagXpress7_1->hDIB,Image_Width,Image_Height);

  }
  else if (ImagXpress7_1->IBPP==8)
  {
    ImagXpress7_1->SaveToBuffer = true;
    ImagXpress7_1->SaveFileType =  FT_TIFF;
    ImagXpress7_1->SaveFile ();
    HANDLE hIM = (HANDLE)ImagXpress7_1->SaveBufferHandle;

    IM = (unsigned char *)GlobalLock(hIM);
    if (IM==NULL)
       ShowMessage("Problem with IM");
    long l = GlobalSize(hIM);
    long offs=l-Image_Width*Image_Height;
    GlobalUnlock(hIM);
    for (int y=0;y<Image_Height;y++)
       for (int x=0;x<Image_Width;x++)
                  *(IM+y*Image_Width+x) = *(IM+y*Image_Width+x+offs);



    }


  if (IM!=NULL)
  {
      bitmap = new (nothrow) unsigned char*[Image_Height];
      if (bitmap==NULL)
      {
         ShowMessage("Problem with bitmap");
      }
      for(int y = 0; y <Image_Height; y++)
            {
            bitmap[y] = new (nothrow) unsigned char [Image_Width];
            if (bitmap[y]==NULL)
              ShowMessage("Problem with bitmap");
             }

     for (int y=0;y<Image_Height;y++)
           for (int x=0;x<Image_Width;x++)
                   {
                      bitmap[y][x] = *(IM+y*Image_Width+x);
                      }*/
     /* Kostas - Do not Erase Image Frame
     for (int x=0;x<Image_Width;x++)
            {
             bitmap[0][x]=255;
             bitmap[Image_Height-1][x]=255;

            }
     for (int y=0;y<Image_Height;y++)
             {
              bitmap[y][0]=255;
              bitmap[y][Image_Width-1]=255;
             } */


     /*if (ImagXpress7_1->IBPP==1)
     {
     free(IM);
     IM=NULL;
     }
     else
     {
     GlobalFree(IM);
     IM=NULL;
     }
 }
}*/

/*
unsigned char *GrayScaleImage::CreateIM(HANDLE Dib,int Ix,int mIy)
{

  unsigned char *IM;
  int x,y,i,j;
  HANDLE hIM = Dib;

  int Iy = mIy;

  unsigned char *IMAGE = (unsigned char *)GlobalLock(hIM);IMAGE+=48;
  int IBx = Ix/8;if (IBx*8!=Ix) IBx++;

  int IBxx = IBx,IBxx4;
  int dx = 0;
  IBxx4 = IBxx/4;
  while (IBxx4*4!=IBxx)
  {
     IBxx++;IBxx4 = IBxx/4;dx++;
  }

  int Dx = IBx+dx;

  if ((IM = (unsigned char *)calloc ((Ix)*(Iy),1)) == NULL)  return NULL;


  for (x=0;x<Ix;x++)
     for (y=0;y<Iy;y++)
       *(IM+y*Ix+x) = ypo(x+1,y+1,Iy,Dx,IMAGE);

  GlobalUnlock(hIM);

  return IM;
}*/

/*
unsigned char GrayScaleImage::ypo(int x,int y,int Iy,int Dx,unsigned char *IMAGE)
{

  unsigned int xx2=(x-1) >> 3;
  long se=Dx*(Iy-y)+xx2;
  unsigned char xx=*(IMAGE+se);
	if (((xx & (1 << (8-x+8*xx2)))) != 0 ) return(0);
		else return(1);
}*/


/*void GrayScaleImage::Load(AnsiString FileName)
{
 //Destroy old data and then Load
 this->Clear();
 //Load new data
 TImagXpress7_ *ImagXpress7_1=new TImagXpress7_(Application);
 ImagXpress7_1->FileName=FileName;

 //Kostas - to Handle Color Images
 if(ImagXpress7_1->IBPP>8)
        ImagXpress7_1->ColorDepth(8,IPAL_Gray,DI_None);

 ImageFileName=FileName;
 Image_Width=ImagXpress7_1->IWidth;
 Image_Height=ImagXpress7_1->IHeight;
 StoreImage(ImagXpress7_1);
 delete ImagXpress7_1;
}*/


int GrayScaleImage::Get_Image_Width()
{
        return this->Image_Width;
}


int GrayScaleImage::Get_Image_Height()
{
        return this->Image_Height;
}

unsigned char GrayScaleImage::Get_Bitmap(int x,int y)
{
  //returns 0 or 1 if the coordinates are fine else returns 255
  //Kostas - added x>=0 kai y>=0
  if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height))
        return bitmap[y][x];
  else
        return 255;
}


/*AnsiString GrayScaleImage::Get_ImageFileName()
{
  return ImageFileName;
}*/

//YES
void GrayScaleImage::Clear()
{
 if (bitmap!=NULL)
 {
    for (int y=0;y<Image_Height;y++)
       delete bitmap[y];
    delete bitmap;
    bitmap=NULL;
 }

 if(Kbitmap!=NULL)
 {
	 free(Kbitmap);
	 Kbitmap=NULL;
 }

 ImageB.Clear();
 this->ImageFileName="";
 this->Image_Width=0;
 this->Image_Height=0;

}

//Noise Parameter
int GrayScaleImage::Otsu(bool doubleChar,double Kaugment)
{

      const int MAXVAL=256;
      int    hist[MAXVAL];
      double pdf[MAXVAL]; //probability distribution
      double cdf[MAXVAL]; //cumulative probability distribution
      double myu[MAXVAL];   // mean value for separation
      double max_sigma, sigma[MAXVAL]; // inter-class variance

      /* Histogram generation */
      for(int i=0; i<MAXVAL; i++)
      {
          hist[i] = 0;
      }
      
			if(doubleChar!=false)
			{
				for(int x=0; x<Image_Width; x++)
				{
          for(int y=0; y<Image_Height; y++)
          {
              hist[bitmap[y][x]]++;
          }
				}
			}
			else
			{
				for(int x=0; x<Image_Width; x++)
				{
          for(int y=0; y<Image_Height; y++)
          {
              hist[(*(Kbitmap+y*Image_Width+x))]++;
          }
				}
			}

      /* calculation of probability density */
      for(int i=0; i<MAXVAL; i++)
      {
          pdf[i] = (double)hist[i] / (Image_Width * Image_Height);
      }

      /* cdf & myu generation */
      cdf[0] = pdf[0];
      myu[0] = 0.0;       /* 0.0 times prob[0] equals zero */
      for(int i=1; i<MAXVAL; i++)
      {
          cdf[i] = cdf[i-1] + pdf[i];
          myu[i] = myu[i-1] + i*pdf[i];
      }

      /* sigma maximization
         sigma stands for inter-class variance
         and determines optimal threshold value */
      int threshold = 0;
      max_sigma = 0.0;
      for(int i=0; i<MAXVAL-1; i++)
      {
          if(cdf[i] != 0.0 && cdf[i] != 1.0)
          {
              double p1p2 = cdf[i]*(1.0 - cdf[i]);
              double mu1mu2diff = myu[MAXVAL-1]*cdf[i]-myu[i];
              sigma[i] = mu1mu2diff * mu1mu2diff / p1p2;
          }
          else
              sigma[i] = 0.0;
          if(sigma[i] > max_sigma)
          {
              max_sigma = sigma[i];
              threshold = i;
          }
      }

      ImageB.Create_Bitmap(Image_Width,Image_Height);

      threshold=(int)((double)threshold*(Kaugment)); //Kostas

      if(doubleChar!=false)
			{
				for(int x=0; x<Image_Width; x++)
				{
          for(int y=0; y<Image_Height; y++)
          {
               if (bitmap[y][x] > threshold)
                  ImageB.Set_Bitmap(x,y,0);
              else
                  ImageB.Set_Bitmap(x,y,1);
          }
				}
			}
			else
			{
				for(int x=0; x<Image_Width; x++)
				{
          for(int y=0; y<Image_Height; y++)
          {
               if ((*(Kbitmap+y*Image_Width+x))> threshold)
                  ImageB.Set_Bitmap(x,y,0);
              else
                  ImageB.Set_Bitmap(x,y,1);
          }
				}
			}
      
  return 0;
}


int GrayScaleImage::Sauvola()
{
  float k=0.3;
  //int w=40;
  int whalf=20;

  // Calculate the integral image, and integral of the squared image
  //Desmevo mnimi gia tous 4 pinakes  xrisimopoiontas long int

  long int **integral_image=new (nothrow) long int* [Image_Height];
  long int **rowsum_image=new (nothrow) long int* [Image_Height];
  long int **integral_sqimg=new (nothrow) long int* [Image_Height];
  long int **rowsum_sqimg=new (nothrow) long int* [Image_Height];

  for (int y=0;y<Image_Height;y++)
    {
       integral_image[y]=new (nothrow) long int [Image_Width];
       rowsum_image[y]=new (nothrow) long int [Image_Width];
       integral_sqimg[y]=new (nothrow) long int [Image_Width];
       rowsum_sqimg[y]=new (nothrow) long int [Image_Width];
    }

  int xmin,ymin,xmax,ymax;
  double diagsum,idiagsum,diff,sqdiagsum,sqidiagsum,sqdiff,area;
  double mean,std,threshold;

  for(int j=0; j<Image_Height; j++)
  {
      rowsum_image[j][0] = bitmap[j][0];
      rowsum_sqimg[j][0] = bitmap[j][0]*bitmap[j][0];
  }
  for(int i=1; i<Image_Width; i++)
  {
      for(int j=0; j<Image_Height; j++)
      {
          rowsum_image[j][i] = rowsum_image[j][i-1] + bitmap[j][i];
          rowsum_sqimg[j][i] = rowsum_sqimg[j][i-1] + bitmap[j][i]*bitmap[j][i];
      }
  }

  for(int i=0; i<Image_Width; i++)
  {
      integral_image[0][i] = rowsum_image[0][i];
      integral_sqimg[0][i] = rowsum_sqimg[0][i];
  }
  for(int i=0; i<Image_Width; i++)
  {
      for(int j=1; j<Image_Height; j++)
      {
          integral_image[j][i] = integral_image[j-1][i] + rowsum_image[j][i];
          integral_sqimg[j][i] = integral_sqimg[j-1][i] + rowsum_sqimg[j][i];
      }
  }

  //Calculate the mean and standard deviation using the integral image

  for(int i=0; i<Image_Width; i++)
  {
      for(int j=0; j<Image_Height; j++)
      {
          xmin = max(0,i-whalf);
          ymin = max(0,j-whalf);
          xmax = min(Image_Width-1,i+whalf);
          ymax = min(Image_Height-1,j+whalf);
          area = (xmax-xmin+1)*(ymax-ymin+1);
          // area can't be 0 here
          // proof (assuming whalf >= 0):
          // we'll prove that (xmax-xmin+1) > 0,
          // (ymax-ymin+1) is analogous
          // It's the same as to prove: xmax >= xmin
          // Image_Width - 1 >= 0         since Image_Width > i >= 0
          // i + whalf >= 0               since i >= 0, whalf >= 0
          // i + whalf >= i - whalf       since whalf >= 0
          // Image_Width - 1 >= i - whalf since Image_Width > i
          // --IM
          assert(area);
          if(!xmin && !ymin){ // Point at origin
              diff   = integral_image[ymax][xmax];
              sqdiff = integral_sqimg[ymax][xmax];
          }
          else if(!xmin && ymin){ // first column
              diff   = integral_image[ymax][xmax] - integral_image[ymin-1][xmax];
              sqdiff = integral_sqimg[ymax][xmax] - integral_sqimg[ymin-1][xmax];
          }
          else if(xmin && !ymin){ // first row
              diff   = integral_image[ymax][xmax] - integral_image[ymax][xmin-1];
              sqdiff = integral_sqimg[ymax][xmax] - integral_sqimg[ymax][xmin-1];
          }
          else
          { // rest of the image
              diagsum    = integral_image[ymax][xmax] + integral_image[ymin-1][xmin-1];
              idiagsum   = integral_image[ymin-1][xmax] + integral_image[ymax][xmin-1];
              diff       = diagsum - idiagsum;
              sqdiagsum  = integral_sqimg[ymax][xmax] + integral_sqimg[ymin-1][xmin-1];
              sqidiagsum = integral_sqimg[ymin-1][xmax] + integral_sqimg[ymax][xmin-1];
              sqdiff     = sqdiagsum - sqidiagsum;
          }

          mean = diff/area;
          std  = sqrt((sqdiff - diff*diff/area)/(area-1));
          threshold = mean*(1+k*((std/128)-1));
          ImageB.Create_Bitmap(Image_Width,Image_Height);
          if(bitmap[j][i] < threshold)
              ImageB.Set_Bitmap(i,j,1);
          else
              ImageB.Set_Bitmap(i,j,0);
      }
  }

  for (int y=0;y<Image_Height;y++)
  {
       delete integral_image[y];
       delete rowsum_image[y];
       delete integral_sqimg[y];
       delete rowsum_sqimg[y];
  }
  delete integral_image;
  delete rowsum_image;
  delete integral_sqimg;
  delete rowsum_sqimg;

  integral_image=NULL;
  rowsum_image=NULL;
  integral_sqimg=NULL;
  rowsum_sqimg=NULL;

  return 0;
}

int GrayScaleImage::max(int a,int b)
{
  if (a>b)
     return a;
  else
     return b;
}

int GrayScaleImage::min(int a,int b)
{
  if (a<b)
     return a;
  else
     return b;
}

int GrayScaleImage::Binarize(bool oldOtsu,double Kaugment)
{

	if(oldOtsu)
	{
		if (this->bitmap!=NULL)
		{
			int a=Otsu(oldOtsu,Kaugment);
			return a;
		}
		else
     return -1;
	}
	else
	{
		if (this->Kbitmap!=NULL)
		{
			int a=Otsu(oldOtsu,Kaugment);
			return a;
		}
		else
     return -1;
	}
}


int GrayScaleImage::Binarize(std::string FileName,bool oldOtsu,double Kaugment)
{
   int a;
   if (bitmap!=NULL)
   {
      if (strcmp(FileName.c_str(),"OTSU") || strcmp(FileName.c_str(),"Otsu") || strcmp(FileName.c_str(),"otsu"))
         {
           a=Otsu(oldOtsu,Kaugment);
           return a;
         }
      else if (strcmp(FileName.c_str(),"SAUVOLA") || strcmp(FileName.c_str(),"Sauvola") || strcmp(FileName.c_str(),"sauvola"))
         {
           a=Sauvola();
           return a;
         }
      else if (strcmp(FileName.c_str(),"GPP") || strcmp(FileName.c_str(),"Gpp") || strcmp(FileName.c_str(),"gpp"))
         {
           //Edo prepei na ftiakso mnimi gia input kai meta na tin free
           unsigned char *IM=new unsigned char [Image_Height*Image_Width];
           for (int y=0;y<Image_Height;y++)
                for (int x=0;x<Image_Width;x++)
                        *(IM+y*Image_Width+x)=Get_Bitmap(x,y);
           a=GPP(IM,Image_Width,Image_Height);
           delete IM;
           IM=NULL;
           return a;
         }
      else if (strcmp(FileName.c_str(),"ALLT") || strcmp(FileName.c_str(),"Allt") || strcmp(FileName.c_str(),"allt"))
         {
          int Ix=Image_Width;
          int Iy=Image_Height;
          unsigned char *IM=new unsigned char [Image_Height*Image_Width];
          for (int y=0;y<Image_Height;y++)
               for (int x=0;x<Image_Width;x++)
                      *(IM+y*Image_Width+x)=Get_Bitmap(x,y);
          int ALLT_swParam;
          ALLT_swParam = ALLT_SW(IM,Image_Width,Image_Height);

          unsigned char *EIK;
          int DX,DY,L,wind;

          wind=2*ALLT_swParam+1;

          L=ALLT_swParam+wind/2;

          //Ayksanoyme tis diastaseis symfwna me to Stroke width
          DX=Ix-2+2*L;
          DY=Iy-2+2*L;

          //Xreiazomaste megalyteri eikona dioti tha kanoyme ypologismoys se
          //parathyro SWxSW poy to kentro toy vrisketai se apostasi SW
          EIK = (unsigned char *) calloc((DX+3)*(DY+3),sizeof(unsigned char)) ;

          for (int y=0;y<Iy;y++)
          for (int x=0;x<Ix;x++)
                  *(EIK+DX*(L+y)+L+x) = *(IM+y*Ix+x);

          a=ALLT(EIK,ALLT_swParam,Image_Width,Image_Height);
          delete IM;
          IM=NULL;
          free(EIK);
          EIK=NULL;
          return a;
         }
      else
         return -1;
   }
   else
   {
      return -1;
   }
}


int GrayScaleImage::GPP(unsigned char *in,int Ix,int Iy)
{
  //params
  int dW = 15;
  float k = 1;k=k/10;
  int R = 128;
  float q = (float)1/1.6;
  float p1 = 0.5;
  float p2 = 0.8;
  int  upsampling = 1;
  int dW1 = 20;

  //int RH[256],GH[256],BH[256];
  unsigned char *IMAGE2;
  unsigned char *IMAGE11;
  unsigned char I,I1;
  int pix,gray,gray2,TH,graygray,pix1,gray1;
  int ydWIx,ydW_1Ix,yIx,ydW1Ix,ydW1_1Ix;
  float m,s;
  int d;
  int d2;

  //float a;
  //float b;

  float aa,bb,cc,dd;

  unsigned char *IMAGE;
  int offs;

  //for (int i=0;i<256;i++) {RH[i]=0;GH[i]=0;BH[i]=0;}

  int *IX_pix, *IX_gray, *IX_graygray, *IX_pix1, *IX_gray1, *IX_graygray1;

  IMAGE=(unsigned char*) malloc((Ix+1)*(Iy+1)*sizeof(unsigned char));
  IMAGE2=(unsigned char*) malloc((Ix+1)*(Iy+1)*sizeof(unsigned char));
  IMAGE11=(unsigned char*) malloc((Ix+1)*(Iy+1)*sizeof(unsigned char));

  for(int y=0;y<Iy;y++)
  for(int x=0;x<Ix;x++)
  {
        *(IMAGE+y*Ix+x)=*(in+y*Ix+x);
        *(IMAGE11+y*Ix+x)=*(in+y*Ix+x);
  }

  offs=0;

  IX_pix=(int*) malloc((Ix+1)*sizeof(int));
  IX_gray=(int*) malloc((Ix+1)*sizeof(int));
  IX_graygray=(int*) malloc((Ix+1)*sizeof(int));
  IX_pix1=(int*) malloc((Ix+1)*sizeof(int));
  IX_gray1=(int*) malloc((Ix+1)*sizeof(int));
  IX_graygray1=(int*) malloc((Ix+1)*sizeof(int));


  //1st step
   for (int y= 0;y<Iy;y++)
  {
    ydWIx = offs+(y+dW)*Ix;
    ydW_1Ix = offs+(y-dW-1)*Ix;
    yIx = offs+y*Ix;

    if (y==0)
    {
      for (int x= 0;x<Ix;x++)
      {
        pix=0;gray=0;gray2=0;
        for (int iy=0;iy<=dW;iy++)
        {
          I = *(IMAGE+offs+iy*Ix+x);
          gray+=I;gray2+=I*I;
        }
        pix+=dW+1;
        IX_pix[x]=pix;IX_gray[x]=gray;IX_graygray[x]=gray2;
        IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];IX_graygray1[x]=IX_graygray[x];
      }
    }
    else
    if (y<=dW)
    {
      for (int x= 0;x<Ix;x++)
      {
         I = *(IMAGE+ydWIx+x);
         IX_pix[x] = IX_pix1[x] + 1;
         IX_gray[x] = IX_gray1[x] + I;
         IX_graygray[x] = IX_graygray1[x] + I*I;
         IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];IX_graygray1[x]=IX_graygray[x];
      }
    }
    else
    if (y<=Iy-1-dW)
    {
      for (int x= 0;x<Ix;x++)
      {
         I = *(IMAGE+ydWIx+x);I1 = *(IMAGE+ydW_1Ix+x);
         IX_pix[x] = IX_pix1[x] ;
         IX_gray[x] = IX_gray1[x]+ I-I1;;
         IX_graygray[x] = IX_graygray1[x]+ I*I-I1*I1;
         IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];IX_graygray1[x]=IX_graygray[x];
      }
    }
    else
    {
      for (int x= 0;x<Ix;x++)
      {
         I = *(IMAGE+ydW_1Ix+x);
         IX_pix[x] = IX_pix1[x]-1;
         IX_gray[x] = IX_gray1[x]-I;
         IX_graygray[x] = IX_graygray1[x]-I*I;
         IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];IX_graygray1[x]=IX_graygray[x];
      }
    }

    int pix1,gray1,graygray1;
    for (int x= 0;x<Ix;x++)
    {
      if (x==0)
      {
        pix=0;gray=0;graygray=0;
        for (int ix=0;ix<=dW;ix++)
        {
          pix+=IX_pix[ix];gray+=IX_gray[ix];graygray+=IX_graygray[ix];
        }
      }
      else
      if (x<=dW)
      {
         pix = pix1 + IX_pix[x+dW];
         gray = gray1 + IX_gray[x+dW];
         graygray = graygray1 + IX_graygray[x+dW];
      }
      else
      if (x<=Ix-1-dW)
      {
         pix = pix1 + IX_pix[x+dW] - IX_pix[x-dW-1];
         gray = gray1 + IX_gray[x+dW] - IX_gray[x-dW-1];
         graygray = graygray1 + IX_graygray[x+dW] - IX_graygray[x-dW-1];
      }
      else
      {
         pix = pix1 - IX_pix[x-dW-1];
         gray = gray1 - IX_gray[x-dW-1];
         graygray = graygray1 - IX_graygray[x-dW-1];
      }

      pix1 = pix;gray1=gray;graygray1=graygray;

      m = gray/pix;
      s = sqrt(fabs((float)(pix*graygray-gray*gray)/(pix*(pix-1))));

      TH = m*(1-k*(1-(float)s/R));

      if ((*(IMAGE+yIx+x)>=TH)) *(IMAGE2+yIx+x)=255; else *(IMAGE2+yIx+x)=0;
    }
  }

  //Interpolation

  for (int y= 0;y<Iy;y++)
  {
    ydW1Ix = offs+(y+dW1)*Ix;
    ydW1_1Ix = offs+(y-dW1-1)*Ix;
    yIx = offs+y*Ix;

    if (y==0)
    {
      for (int x= 0;x<Ix;x++)
      {
        int pix=0;int gray=0;
        for (int iy=0;iy<=dW1;iy++)
        {
          if (*(IMAGE2+offs+iy*Ix+x)==255)
          {
            pix++;gray+=*(IMAGE+offs+iy*Ix+x);
          }
        }
        IX_pix[x]=pix;IX_gray[x]=gray;
        IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];
      }
    }
    else
    if (y<=dW1)
    {
      for (int x= 0;x<Ix;x++)
      {
         if (*(IMAGE2+ydW1Ix+x)==255)
         {
           IX_pix[x] = IX_pix1[x] + 1;
           IX_gray[x] = IX_gray1[x] + *(IMAGE+ydW1Ix+x);
         }
         IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];
      }
    }
    else
    if (y<=Iy-1-dW1)
    {
      for (int x= 0;x<Ix;x++)
      {
         IX_pix[x] = IX_pix1[x];
         IX_gray[x] = IX_gray1[x];

         if (*(IMAGE2+ydW1Ix+x)==255)
         {
           IX_pix[x] ++;
           IX_gray[x] += *(IMAGE+ydW1Ix+x);
         }
         if (*(IMAGE2+ydW1_1Ix+x)==255)
         {
           IX_pix[x] --;
           IX_gray[x] -= (*(IMAGE+ydW1_1Ix+x));
         }

         IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];
      }
    }
    else
    {
      for (int x= 0;x<Ix;x++)
      {
        if (*(IMAGE2+ydW1_1Ix+x)==255)
        {
          IX_pix[x] = IX_pix1[x]-1;
          IX_gray[x] = IX_gray1[x]-(*(IMAGE+ydW1_1Ix+x));
        }
        IX_pix1[x]=IX_pix[x];IX_gray1[x]=IX_gray[x];
      }
    }

    for (int x= 0;x<Ix;x++)
    {

      pix=pix1;gray=gray1;
      if (x==0)
      {
        pix=0;gray=0;
        for (int ix=0;ix<=dW1;ix++)
        {
           pix+=IX_pix[ix];gray+=IX_gray[ix];
        }
      }
      else
      if (x<=dW1)
      {
         pix = pix1 + IX_pix[x+dW1];
         gray = gray1 + IX_gray[x+dW1];
      }
      else
      if (x<=Ix-1-dW1)
      {
          pix = pix1 + IX_pix[x+dW1] - IX_pix[x-dW1-1];
          gray = gray1 + IX_gray[x+dW1]- IX_gray[x-dW1-1];
      }
      else
      {
         pix = pix1 - IX_pix[x-dW1-1];
         gray = gray1 - IX_gray[x-dW1-1];
      }

      pix1 = pix;gray1=gray;

      //NEW
      if (pix==0) *(IMAGE11+yIx+x) = 255;
      else
      if (*(IMAGE2+yIx+x)==0)
        *(IMAGE11+yIx+x)=(float)gray/pix;
      else  *(IMAGE11+yIx+x) = *(IMAGE+yIx+x);

    }
  }

  //Thresholding

  long PixFor = 0;
  long D=0,D1=0,D2=0;
  int Hist[256];for (int i=0;i<256;i++) Hist[i]=0;

  for (int y= 0;y<Iy;y++)
  {
    int yIx = offs+y*Ix;
    for (int x= 0;x<Ix;x++)
    {

      I = *(IMAGE+yIx+x); I1 = *(IMAGE11+yIx+x);
      if ((*(IMAGE2+yIx+x)==0) && (I1>I))
      {
        PixFor++;D=D+I1-I;D1=D1+I;D2=D2+I1;
        Hist[I]++;
      }
    }
  }
  D = D/PixFor;D1 = D1/PixFor;D2 = D2/PixFor;

  long HistMax=0;int H;

  for (int i=0;i<256;i++)
  {
    if (Hist[i]>HistMax) {HistMax = Hist[i];H=i;}
  }
  H+=5;

  d = D;
  d2 = D2;

  //a = ((float)d*(p2-1))/((float)d2*q*(p1-1));
  //b = ((float)d*(p1-p2))/((float)q*(p1-1));
  aa = q*d*(1-p2);
  bb = -(2*(1+p1))/(1-p1);
  cc = p2*d*q;
  dd = -4/(d2*(1-p1));

  TH = (float)D*0.4;
  if (upsampling==1)
  {
    for (int y= 0;y<Iy;y++)
    {
      yIx = offs+y*Ix;
      for (int x= 0;x<Ix;x++)
      {
        if (*(IMAGE2+yIx+x)==0)
        {
          TH = aa/(1+exp(dd*(*(IMAGE+yIx+x))-bb))+cc;

          if (*(IMAGE11+yIx+x)-(*(IMAGE+yIx+x))>TH)
          {
            *(IMAGE2+yIx+x)=0;     //cc
          }
          else
          {
            *(IMAGE2+yIx+x)=255;     //cc
          }
        }
        else
        {
          *(IMAGE2+yIx+x)=255;       //cc
        }
      }
    }

  }

  ImageB.Create_Bitmap(Ix,Iy);

    for (int x=0;x<Ix;x++)
    for (int y=0;y<Iy;y++)
    {
        ImageB.Set_Bitmap(x,y,(*(IMAGE2+y*Ix+x))%254);
        //*(out+y*Ix+x)=(*(IMAGE2+y*Ix+x))%254;
    }



  free(IMAGE);IMAGE=NULL;
  free(IMAGE2);IMAGE2=NULL;
  free(IMAGE11);IMAGE11=NULL;
  free(IX_pix);free(IX_gray);free(IX_graygray);free(IX_pix1);free(IX_gray1);free(IX_graygray1);
  IX_pix=NULL;
  IX_gray=NULL;
  IX_graygray=NULL;
  IX_pix1=NULL;
  IX_gray1=NULL;
  IX_graygray1=NULL;
  return 0;
}




int GrayScaleImage::ALLT(unsigned char *EIK,int SW,int Ix,int Iy)
{
        ImageB.Create_Bitmap(Ix,Iy);
        int i,j,LP_j,x_j,y_j,a,b,telos;
        int nai=0;
        struct ALLT_struct camel[8];

        int Thres;
        int wind=2*SW+1;

        int L=SW+wind/2;
        int DX=Ix-2+2*L;

        camel[0].x=-SW;
        camel[0].y=0;

        camel[1].x=-SW;
        camel[1].y=-SW;

        camel[2].x=0;
        camel[2].y=-SW;

        camel[3].x=SW;
        camel[3].y=-SW;

        camel[4].x=SW;
        camel[4].y=0;

        camel[5].x=SW;
        camel[5].y=SW;

        camel[6].x=0;
        camel[6].y=SW;

        camel[7].x=-SW;
        camel[7].y=SW;



        for(b=0;b<Iy;b++)
{
        for(a=0;a<Ix;a++)
        {
                Thres=ALLT_Ypologismos_Thres(EIK,SW,0.35,a,b,wind,Image_Width,Image_Height);
                telos=4;
        for(i=0;i<telos;i++)
        {
                for(j=i;j<(i+2);j++)
                {
                        x_j=a+camel[j].x;
                        y_j=b+camel[j].y;

                        LP_j=ALLT_Ypologismos_Pj(EIK,SW,x_j,y_j,wind,Image_Width,Image_Height)-(*(EIK+(L+b)*DX+L+a));

                        if(LP_j>Thres)
                                nai=1;
                        else
                        {
                                nai=0;
                                break;
                        }
                        ///////////////////////////////

                        x_j=a+camel[(j+4)%8].x;
                        y_j=b+camel[(j+4)%8].y;

                        LP_j=ALLT_Ypologismos_Pj(EIK,SW,x_j,y_j,wind,Image_Width,Image_Height)-(*(EIK+(L+b)*DX+L+a));
                        if(LP_j>Thres)
                                nai=1;
                        else
                        {
                                nai=0;
                                break;
                        }
                }

                if(nai==1)
                        break;
                else
                {
                        if(j==0)
                           telos--;
                }

        }

        if(nai==1)
                ImageB.Set_Bitmap(a,b,0);//*(EIKout+b*Ix+a)=0;
        else
                ImageB.Set_Bitmap(a,b,1);//*(EIKout+b*Ix+a)=1;
        }
}

return 0;
}

int GrayScaleImage::ALLT_Ypologismos_Pj(unsigned char *EIK,int SW,int x_j, int y_j,int wind,int Ix,int Iy)
{
  int konta;
  int sum=0;
  int k,n,average;

  int L,DX;

  konta=wind/2;
  L=SW+konta;
  DX=Ix-2+2*L;

  for(n=y_j-konta;n<y_j+konta+1;n++)
  {
        for(k=x_j-konta;k<x_j+konta+1;k++)
        {
                //sum+=(*(IM+n*Ix+k));
                sum+=(*(EIK+(L+n)*DX+L+k));
        }
  }

  average=sum/(wind*wind);

       return average;
}

int GrayScaleImage::ALLT_Ypologismos_Thres(unsigned char *EIK,int SW,double Glob_a,int x_t, int y_t,int wind,int Ix,int Iy)
{
int L,DX;

  int Xam,Nim,absolute_1,absolute_2;
  int sum=0;
  int strok;
  double loc_b;
  int k,n,ave,yes,temp,i,plin;
  bool bg=false;
  bool fg=false;
  temp=wind;

  L=SW+wind/2;
  DX=Ix-2+2*L;
  Xam=0;
  Nim=255;
  plin=0;
  strok=0;
  for(i=0;;i++)
  {

  for(n=y_t-SW-strok;n<=y_t+SW+strok;n++)
  {
        for(k=x_t-SW-strok;k<=x_t+SW+strok;k++)
        {
                //sum+=(*(IM+n*Ix+k));
                if(n<0 || k<0 || n>Iy-1 || k>Ix-1)
                {
                        plin++;
                        continue;
                }
                if(*(EIK+(L+n)*DX+L+k)>Xam)
                        Xam=*(EIK+(L+n)*DX+L+k);
               if(*(EIK+(L+n)*DX+L+k)<Nim)
                        Nim=*(EIK+(L+n)*DX+L+k);

                sum+=(*(EIK+(L+n)*DX+L+k));
        }
  }

  if((temp*temp)<plin+1)
        plin=0;
  ave=sum/((temp*temp)-plin);

  absolute_1=abs(ave-Nim); //diafora fwtinotitas, skotadi-meso
  absolute_2=abs(Xam-ave); //diafora fwtinotitas, fws-meso

  loc_b=1.0;
  //An eimaste pros to skotadi
  if(absolute_1<absolute_2){
        //yes=Glob_a*(2*ave+Nim)/3;
        //if(CheckBox2->Checked)
        //        loc_b=(double)((double)(Nim*Nim)/(double)(ave*ave)); //1
        //if((Xam-Nim)!=0 && (ave-Nim)!=0)
        //        loc_b=(double)((Xam-ave)/(double)(Xam-Nim)); //ave-Nim
        //else
        //        loc_b=(double)(Nim/(double)ave);
        //if(absolute_1!=0 && absolute_2!=0)
        //loc_b=(double)absolute_1/(double)absolute_2;
        yes=loc_b*Glob_a*(2*Nim+ave)/3; //yes
        //if(absolute_2-absolute_1<10 && yes<15)
        //yes+=25;
        break;}

  //An eimaste pros to fws
  if(absolute_2<absolute_1){
        //yes=Glob_a*(2*Nim+ave)/3;
        //if(CheckBox2->Checked)
       // loc_b=(double)((double)(ave*ave)/(double)(Xam*Xam)); //1
        //if((Xam-Nim)!=0 && (Xam-ave)!=0)
       //         loc_b=(double)((ave-Nim)/(double)(Xam-Nim));  //Xam-ave
        //else
        //        loc_b=(double)(ave/(double)Xam);
        //if(absolute_1!=0 && absolute_2!=0)
        //loc_b=(double)absolute_2/(double)absolute_1;
        yes=loc_b*Glob_a*(2*ave+Nim)/3; //yes
        //if(ave<96)
        //yes+=25;
        break;}
  if(absolute_1==absolute_2)
  {
        //if(i==0)
        //{
        if(!bg)
        {
           if(Nim==Xam)
           {
                temp+=2;
                sum=0;
                Xam=-1;
                Nim=256;
                strok++;
                bg=true;
           }
           else
           {
                if(fg)
                {
                        yes=Glob_a*ave;
                        break;
                }
                temp+=2;
                sum=0;
                Xam=-1;
                Nim=256;
                strok++;
                fg=true;
           }
        }
        else
        {
                yes=Glob_a*ave;
                break;
        }
        //}
        //if(i==1){
        //yes=Glob_a*ave; break;}
  }

  }


  return yes;
}




int GrayScaleImage::ALLT_SW(unsigned char *IM,int Ix,int Iy)
{
        int sum,sum_1,ena;
        int i,j,k,l,p,q,arxi_x,arxi_y,telos_x,telos_y,a,b;
        int dyo,trexon,max_1,max_2,thesi_1,thesi_2,thesi,min;
        int Ni,Co,SW;
        int *pinax;
        int *pinax_2;
        int *xanip;
        int *pinakaki;
        int Hist[255];
        bool vges;
        int Run_thres;


        //Ni = plithos parathirwn
        Ni = 6;

        //Timi katoflioy gia ton ypologismo SW
        Co = 129;

        ena=(Ix-1)/Ni;
        dyo=(Iy-1)/Ni;
        telos_x=telos_y=-1;
        sum=sum_1=0;

        pinax=(int *)(malloc((2*Ni)*sizeof(int)));
        pinax_2=(int *)(malloc((2*Ni)*sizeof(int)));
        pinakaki=(int *)(malloc((ena+1)*sizeof(int)));
        xanip=(int *)(malloc((ena+1)*sizeof(int)));

        for(k=0;k<ena+1;k++)
        {
                xanip[k]=0;
                pinakaki[k]=0;
        }

  //Gia na kathorisw to xrwma twn grammatwn
  for (int x=0;x<256;x++) Hist[x]=0;

        for (int a=0;a<Ix;a++)
        for (int b=0;b<Iy;b++)
        Hist[*(IM+b*Ix+a)]++;

        Run_thres=0;

        //Koryfh 1 ews to Co
        for(k=1;k<Co;k++)
                if(Hist[k]>Run_thres)
                {
                        Run_thres=Hist[k];
                        thesi_1=k;
                }

        //An i prwti korifi ews to Co einai poly konta sto Co
        //tote kane Co th deyteri korifi
        if(thesi_1>Co-7)
                thesi_2=Co;

        //Vathos meta thn koryfh 1
        thesi=thesi_1+10;
        min=Hist[thesi_1];
         for(l=thesi_1;l<Co-thesi_1;l++)
                if(Hist[l]<min && Hist[l]!=0)
                {
                        min=Hist[l];
                        thesi=l;
                }

        //Koryfh 2
        Run_thres=0;
        thesi_2=thesi+5;
        for(k=thesi;k<Co;k++)
                if(Hist[k]>Run_thres)
                {
                        Run_thres=Hist[k];
                        thesi_2=k;
                }

         //Deytero Vathos ws 128 peripoy
         min=Hist[thesi_2];
         for(l=thesi_2-1;l<Co;l++)
                if(Hist[l]<min && Hist[l]!=0)
                {
                        min=Hist[l];
                        thesi=l;
                }

        if(thesi_2<Co-7)
        {

                Run_thres=thesi;
        }
        else
                Run_thres=thesi_1;



       /*if(Run_thres<26)
                Glob_a = 70;
                else
        if(Run_thres<52)
                Glob_a = 60;
                else
        if(Run_thres<78)
                Glob_a = 40;
                else
        if(Run_thres<104)
                Glob_a = 30;
                else
        if(Run_thres<130)
                Glob_a = 20;
                else
                Glob_a = 10;

        Glob_a/=100;*/

  //H prwth diagwnios (6x6)
   for(i=0;i<Ni;i++)
   {
        arxi_x=telos_x+1;
        arxi_y=telos_y+1;
        telos_x+=ena;
        telos_y+=dyo;

        //Kataskeyi Istogrammatos sto pinakaki
        //me vasi to Run_thres gia to parathyro
        for(p=arxi_y;p<=telos_y;p++)
        {
                trexon=0;
                vges=true;
        for(q=arxi_x;q<=telos_x;q++)
        {
                if(*(IM+p*Ix+q)<=Run_thres)
                        trexon+=1;
                else
                {
                        if(vges)
                        {
                                vges=false;
                                trexon=0;
                                continue;
                        }
                        if(trexon>ena/2)
                        {
                                trexon=0;
                                continue;
                        }
                        pinakaki[trexon]+=1;
                        trexon=0;
                }
        }
                pinakaki[0]=0;
        }


        //eyresi sto a to megisto SW
        max_1=2;
        a=0;
        for(k=2;k<ena;k++)
        {
            if(pinakaki[k]>=max_1)
            {
                max_1=pinakaki[k];
                a=k;
            }
        }

        //eyresi to deytero megalytero SW sto b
        max_2=2;
        pinakaki[a]=0;
        b=0;
        for(l=2;l<ena;l++)
        {
            if(pinakaki[l]>=max_2)
            {
                max_2=pinakaki[l];
                b=l;
            }
        }

        //Ston pinax_2 vazoyme to b
        //enw ston pinax vazoyme to a
        pinax_2[i]=b;

        //Synthiki me deyterh koryfh
        pinax[i]=a;
        if( (b-a)==1 || (b-a)==2)
                if(10*max_2>=8*max_1)
                       pinax[i]=b;

    for(k=0;k<ena+1;k++)
                pinakaki[k]=0;

   }

   max_1=-1;
   for(l=0;l<Ni;l++)
                if(pinax[l]>max_1)
                        max_1=pinax[l];
  for(k=0;k<Ni;k++)
                if(pinax[k]<2)
                        pinax[k]=max_1;

  max_2=-1;
  for(l=0;l<Ni;l++)
                if(pinax_2[l]>max_2)
                        max_2=pinax_2[l];
  for(k=0;k<Ni;k++)
                if(pinax_2[k]<2)
                        pinax_2[k]=max_2;

   arxi_x=Ix;
   telos_y=-1;

   //H deyterh diagwnios (6x6)
   for(j=0;j<Ni;j++)
   {
        telos_x=arxi_x-1;
        arxi_x-=ena;
        arxi_y=telos_y+1;
        telos_y+=dyo;

        for(p=arxi_y;p<=telos_y;p++)
        {
                trexon=0;
                vges=true;
        for(q=arxi_x;q<=telos_x;q++)
        {
                if(*(IM+p*Ix+q)<=Run_thres)
                        trexon+=1;
                else
                {
                        if(vges)
                        {
                                vges=false;
                                trexon=0;
                                continue;
                        }
                        if(trexon>ena/2)
                        {
                                trexon=0;
                                continue;
                        }
                        pinakaki[trexon]+=1;
                        trexon=0;
                }
        }
                pinakaki[0]=0;
        }

        max_1=2;
        a=0;
        for(k=2;k<ena;k++)
        {
            if(pinakaki[k]>max_1)
            {
                max_1=pinakaki[k];
                a=k;
            }
        }

        max_2=2;
        pinakaki[a]=0;
        b=0;
        for(l=2;l<ena;l++)
        {
            if(pinakaki[l]>max_2)
            {
                max_2=pinakaki[l];
                b=l;
            }
        }

        pinax_2[j+Ni]=b;

        //Symthiki me deyterh koryfh
        pinax[j+Ni]=a;
        if( (b-a)==1 || (b-a)==2)
                if(10*max_2>=8*max_1)
                        pinax[j+Ni]=b;

    for(k=0;k<ena+1;k++)
                pinakaki[k]=0;

   }

   max_1=-1;
   for(l=Ni;l<2*Ni;l++)
                if(pinax[l]>max_1)
                        max_1=pinax[l];
  for(k=Ni;k<2*Ni;k++)
                if(pinax[k]<2)
                        pinax[k]=max_1;
  max_2=-1;
  for(l=Ni;l<2*Ni;l++)
                if(pinax_2[l]>max_2)
                        max_2=pinax_2[l];
  for(k=Ni;k<2*Ni;k++)
                if(pinax_2[k]<2)
                        pinax_2[k]=max_2;

   max_1=-1;
   for(trexon=0;trexon<(2*Ni);trexon++)
   {
        if(pinax[trexon]>=max_1)
                max_1=pinax[trexon];
        xanip[pinax_2[trexon]]+=1;
        sum_1+=pinax[trexon];
   }

   max_2=0;
   for(k=0;k<ena;k++)
   {
        if(xanip[k]>max_2)
        {
                max_2=xanip[k];
                thesi=k;
        }
   }

   if(sum_1%(2*Ni)==0)
        sum_1/=(2*Ni);
   else
   {
        sum_1/=(2*Ni);
        sum_1++;
   }

   sum=(2*sum_1+/*max_1+*/thesi)/3;

   if(sum%2==0)
        sum++;

   SW=sum;

   //int wind=2*SW+1;


 free(pinax);
 free(pinax_2);
 free(xanip);
 free(pinakaki);
 pinax=NULL;
  pinax_2=NULL;
   xanip=NULL;
   pinakaki=NULL;


   //
   /*int DY,DX,L;

   L=SW+wind/2;
   DX=Ix-2+2*L;
   DY=Iy-2+2*L;*/

   return SW;


}

