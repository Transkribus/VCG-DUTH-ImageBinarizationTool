//============================================================================
// Name        : binarization.cpp
// Author      : Kntir_Louloud_Nstam_AlexPap_Bgat
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <stdio.h>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <string>
#include <stdlib.h>
#include <malloc.h>

#include "BinaryImage.h"
#include "GrayScaleImage.h"

#include <iostream>

typedef unsigned char uchar;

using namespace std;
using namespace cv;

//Niblack with/out fix window
bool NibSW(float k0, int s0, unsigned char *inputIm, unsigned char *outputIm, int Ix, int Iy);
void NibSW(float k0, int s0, unsigned char *inputIm, unsigned char *outputIm, int *inputSWIm, int Ix, int Iy);
void Dilationx3(unsigned char *inputImB, unsigned char *outputImB, int Ix, int Iy);
bool Inpaintx4Passes(unsigned char *inputIm, unsigned char *inputMaskB, unsigned char *outputIm, unsigned char *outputImStats, int noizLev, int Ix, int Iy);
bool InpaintAssistPass(int whoPass, unsigned char *ImOriginal, unsigned char *InpMaskB, unsigned char *EstBG, int Ix, int Iy);
int InpaintAssistCheck4(int x, int y, unsigned char *im, unsigned char *out, int Ix, int Iy);
void CalcStatsInpaint(unsigned char Fst, unsigned char Snd, unsigned char Trd, unsigned char Frth, int *minn, int *max, int *avee);
bool ImageNormalization(unsigned char *inputIm, unsigned char *ImDevisor, unsigned char *outputIm, int Ix, int Iy);
int PostProcOtsu(int **inputIm, vector<cc> rCCs, unsigned char **outputIm, int Ix, int Iy);
void PostProcOtsuAssist(int who, int plithos, vector<cc> rCCs, int **UnB, double *cccnt, double *ccall, double *pixelscnt, double *pixelsall);
int CalcSWAssistNearCont(int who, vector<cc> BCCs, unsigned char *Contours, int x, int y, int Ix, int Iy);
void CombineBins(unsigned char **inputImB1, vector<cc> BCCs2, int **Blabels2, double LContrast, unsigned char *Excluding, int noisL);
void check_NIB(int xx, int yy, unsigned char *ioImB, int Ix, int Iy);
void check_EDG(int xx, int yy, unsigned char *ioImB, cv::Mat inMat, int Ix, int Iy);
int Binarization_Function(int Ix, int Iy, unsigned char *IMor, int NoiseParam_Norm, double NoiseParam_Otsu, int NoiseParam_OtsPost, int NoiseParam_Exclude, int sw_window);

Mat SW(Mat& bw8u_orig)
{
	// bw8u : we want to calculate the SWT of this. NOTE: Its background pixels are 0 and foreground pixels are 1 (not 255!)
    Mat bw32f, swt32f, kernel;
    Mat bw8u;
    threshold(bw8u_orig, bw8u, -1 , 1, THRESH_BINARY_INV | THRESH_OTSU);
	
	double min, max;
	int strokeRadius;

	bw8u.convertTo(bw32f, CV_32F);  // format conversion for multiplication
	distanceTransform(bw8u, swt32f, CV_DIST_L2, 5); // distance transform
	minMaxLoc(swt32f, NULL, &max);  // find max
	strokeRadius = (int)ceil(max);  // half the max stroke width
	kernel = getStructuringElement(MORPH_RECT, Size(3, 3)); // 3x3 kernel used to select 8-connected neighbors

	for (int j = 0; j < strokeRadius; j++)
	{
		dilate(swt32f, swt32f, kernel); // assign the max in 3x3 neighborhood to each center pixel
		swt32f = swt32f.mul(bw32f); // apply mask to restore original shape and to avoid unnecessary max propogation
	}
    return swt32f;
}

int main(int argc, char **argv)
{

    //Put Input/Ouptut Filename into strings
    if (argc != 3)
    {
	printf("\nUsage: binarization [infile] [outfile]\n");
	return -1;
    }
    string NameIn = argv[1];
    string NameOut = argv[2];

    //Read Image
    cv::Mat ImageIn;
    ImageIn.data = NULL;
    ImageIn = cv::imread(NameIn, CV_LOAD_IMAGE_GRAYSCALE);
	// Mat halfSWImg = SW(ImageIn);
	// double sw_max;
    // minMaxLoc(halfSWImg, nullptr, &sw_max);




    //Check for Success/Failure
    int nothing = -1;
    if (ImageIn.data == NULL)
    {
	printf("\nError[%d] Loading Image: %s\n", nothing, argv[1]);
	return nothing;
    }

    //Perform Operations
    int Ix = ImageIn.cols;
    int Iy = ImageIn.rows;

    unsigned char *IM = NULL;
    nothing = -3;
    IM = (unsigned char *)malloc(Ix * Iy * sizeof(unsigned char));

    if (IM == NULL)
    {
	printf("\nError Allocating pointer\n");
	ImageIn.release();
	return nothing;
    }

    //load image to imagePointer
    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	    *(IM + y * Ix + x) = ImageIn.at<uchar>(y, x);

    //check apo do kai kato!
    int BinRes;
    BinRes = Binarization_Function(Ix, Iy, IM, 0, 1.0, 0, 1, 20);

    //endtime=clock();
    //printf("\tTime: %6.2lf\n", double((endtime-starttime)));

    nothing = -4;
    if (BinRes != 0)
    {
	printf("\nError[%d] Binarizing Image: %s\n", nothing, argv[1]);
	free(IM);
	ImageIn.release();
	return nothing;
    }

    //Write image to file (hard disc)
    for (int yy = 0; yy < Iy; yy++)
	for (int xx = 0; xx < Ix; xx++)
	    ImageIn.at<uchar>(yy, xx) = (*(IM + yy * Ix + xx));

    //Save Image to disk
    bool checkWrite = cv::imwrite(NameOut, ImageIn);

    //Check for Success/Failure
    nothing = -2;
    if (!checkWrite)
    {
	printf("\nError[%d] Saving Image: %s\n", nothing, argv[2]);
	free(IM);
	ImageIn.release();
	return nothing;
    }

    //Success
    nothing = 0;
    free(IM);

    //Release OpenCV Image
    ImageIn.release();

    return nothing;
}

int Binarization_Function(int Ix, int Iy, unsigned char *IMor, int NoiseParam_Norm, double NoiseParam_Otsu, int NoiseParam_OtsPost, int NoiseParam_Exclude, int sw_window)
{
    unsigned char *NibRes1 = NULL;
    NibRes1 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (NibRes1 == NULL)
	return -1;

    bool checkMe = false;
    int checkthis = 0;

    //1.2 (avoid crashing-- adjust window for small images)
    int nibWind = 30;
    long minDim = (Iy < Ix) ? Iy : Ix;
    if (minDim < 60)
	nibWind = minDim / 2;

    //1.1 Perform Niblack
    NibSW(-2, nibWind, IMor, NibRes1, Ix, Iy);

    //2.3 Get mem for Otsu
    GrayScaleImage UnOts;
    UnOts.Create_KBitmap(Ix, Iy);

    //2.2 Preparing data
    for (int yy = 0; yy < Iy; yy++)
	for (int xx = 0; xx < Ix; xx++)
	    *(UnOts.Kbitmap + yy * Ix + xx) = (*(IMor + yy * Ix + xx));

    //2.1 Perform Otsu
    UnOts.Binarize(false, 1.0);

    //3. (avoid crashing-- detect black borders)
    //Combine Otsu(borders) and Niblack
    int cmb_y = Iy / 11;
    int cmb_x = Ix / 11;

    for (int cty = 0; cty < cmb_y; cty++)
	for (int ctx = 0; ctx < Ix; ctx++)
	    if (UnOts.ImageB.Get_Bitmap(ctx, cty) != 0)
		*(NibRes1 + cty * Ix + ctx) = 1;

    for (int cty = Iy - cmb_y; cty < Iy; cty++)
	for (int ctx = 0; ctx < Ix; ctx++)
	    if (UnOts.ImageB.Get_Bitmap(ctx, cty) != 0)
		*(NibRes1 + cty * Ix + ctx) = 1;

    for (int ctx = 0; ctx < cmb_x; ctx++)
	for (int cty = 0; cty < Iy; cty++)
	    if (UnOts.ImageB.Get_Bitmap(ctx, cty) != 0)
		*(NibRes1 + cty * Ix + ctx) = 1;

    for (int ctx = Ix - cmb_x; ctx < Ix; ctx++)
	for (int cty = 0; cty < Iy; cty++)
	    if (UnOts.ImageB.Get_Bitmap(ctx, cty) != 0)
		*(NibRes1 + cty * Ix + ctx) = 1;

    //4.3 Get mem for InpaintMask
    unsigned char *Dil = NULL;
    Dil = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));

    //4.2 Perform Dilation - Inpaint Mask
    Dilationx3(NibRes1, Dil, Ix, Iy);

    //4.1 (avoid crashing-- total black/white pixels belong to inpainting mask)
    for (int yb = 0; yb < Iy; yb++)
	for (int xb = 0; xb < Ix; xb++)
	    if (*(IMor + yb * Ix + xb) < 1 || *(IMor + yb * Ix + xb) > 254)
		*(Dil + yb * Ix + xb) = 1;

    //Release Mem not needed anymore
    UnOts.Clear();

    //5. Get mem for Statistics BG image
    unsigned char *Stats = NULL;
    Stats = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));

    //6. Perform BG Estimation to NibRes1 and Stats
    Inpaintx4Passes(IMor, Dil, NibRes1, Stats, NoiseParam_Norm, Ix, Iy);

    //MemFree - Release Dil
    if (Dil != NULL)
	free(Dil);

    //7.2 Preparing Image Normalization
    GrayScaleImage UnNorm;
    UnNorm.Create_KBitmap(Ix, Iy);
    //7.1 Normlized image to Kbitmap
    ImageNormalization(IMor, NibRes1, UnNorm.Kbitmap, Ix, Iy);
    //8.3 Perform Otsu on Normalized image
    UnNorm.Binarize(false, NoiseParam_Otsu);

    //8.2 Preparing CC detection
    BinaryImage UnB1;
    UnB1.Create_Bitmap(Ix + 2, Iy + 2);

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	    UnB1.Set_Bitmap(x + 1, y + 1, UnNorm.ImageB.Get_Bitmap(x, y));

    //8.1 Otsu CC detection
    UnB1.Connected_Component_Labeling();

    vector<cc> CCs1;
    CCs1 = UnB1.Get_Connected_Components();

    //9.2 Preparing Intermediate Post-Proc
    BinaryImage UnB1P1;
    UnB1P1.Create_Bitmap(Ix + 2, Iy + 2);

    //9.1 Intermediate postproc
    int removHei = PostProcOtsu(UnB1.labelmap, CCs1, UnB1P1.bitmap, Ix, Iy);

    //10.3 Otsu-post CC detection
    UnB1P1.Connected_Component_Labeling();

    vector<cc> CCs1P;
    CCs1P = UnB1P1.Get_Connected_Components();

    //Also 11.3 BG Stats
    long double BGave = 0;
    long double cntBG = Ix * Iy;
    for (int yy = 0; yy < Iy; yy++)
	for (int xx = 0; xx < Ix; xx++)
	{
	    //UnSkel0.Set_Bitmap(xx+1,yy+1,UnB1P1.Get_Bitmap_Skeleton(xx+1,yy+1));
	    BGave = BGave + *(Stats + yy * Ix + xx);
	}

    long double BGaverage = BGave / cntBG;

    //11.1 Calculating FG/BG features/stats
    long double FGave = 0;
    long double cntFG = 0;
    //bool vrika;
    int Ex = 1;
    for (unsigned int i = 0; i < CCs1P.size(); i++)
    {
	if (UnB1P1.bitmap[CCs1P[i].yg][CCs1P[i].xg])
	{
	    //vrika=true;
	    FGave = FGave + *(IMor + (CCs1P[i].yg - Ex) * Ix + (CCs1P[i].xg - Ex));
	    cntFG = cntFG + 1;
	}
    }

    long double FGaverage = FGave / cntFG;

    double Fsd = 0;
    double Bsd = 0;
    for (int b = 0; b < Iy; b++)
	for (int a = 0; a < Ix; a++)
	{
	    double difB = BGaverage - *(Stats + b * Ix + a);
	    Bsd = Bsd + pow(difB, 2.0);
	}

    for (unsigned int j = 0; j < CCs1P.size(); j++)
    {

	if (UnB1P1.bitmap[CCs1P[j].yg][CCs1P[j].xg])
	{
	    double difF = FGaverage - *(IMor + (CCs1P[j].yg - Ex) * Ix + (CCs1P[j].xg - Ex));
	    Fsd = Fsd + pow(difF, 2.0);
	}
    }

    //MemFree - Release Stats
    if (Stats != NULL)
	free(Stats);

    double BGstd = Bsd / cntBG;
    BGstd = pow(BGstd, 0.5);
    double FGstd = Fsd / cntFG;
    FGstd = pow(FGstd, 0.5);

    
    //Calculate k and s(SW) for noblack
    double ContrastRatio;
    if ((BGaverage - BGstd) < 1)
	ContrastRatio = (FGaverage + FGstd) / BGaverage;
    else if ((FGaverage + FGstd) / (BGaverage - BGstd) > 0.99)
	ContrastRatio = (FGaverage) / (BGaverage);
    else
	ContrastRatio = (FGaverage + FGstd) / (BGaverage - BGstd);

    double LogContrast = -50 * log10(ContrastRatio);

    float ks = (-2.0) + (-1) * (floor((LogContrast - 5) / 10.0));

    //12 Perform Niblack
    NibSW(ks, sw_window,UnNorm.Kbitmap,NibRes1,Ix,Iy);
    //checkMe = GPP_no_dll(UnNorm.Kbitmap, Ix, Iy);
    //NibRes1 = UnNorm.Kbitmap;
  
  
    //Get CCs of Niblack No.2
    BinaryImage UnB2;
    checkthis = UnB2.Create_Bitmap(Ix + 2 * Ex, Iy + 2 * Ex);


    for (int yy = 0; yy < Iy; yy++)
	for (int xx = 0; xx < Ix; xx++)
	    UnB2.Set_Bitmap(xx + Ex, yy + Ex, *(NibRes1 + yy * Ix + xx));

    UnB2.Connected_Component_Labeling();
  

    vector<cc> CCs2;
    CCs2 = UnB2.Get_Connected_Components();

    //Combine Bins

    unsigned char *Excl = NULL;
    Excl = (unsigned char *)calloc(CCs2.size(), sizeof(unsigned char));


    
    //FreeMem - Release UnB1 and CCs1
    if (CCs1.size() > 0)
    {
	CCs1.clear();
	vector<cc>().swap(CCs1);
    }
    else
	vector<cc>().swap(CCs1);

    UnB1.Clear();
    //-//

    CombineBins(UnB1P1.bitmap, CCs2, UnB2.labelmap, LogContrast, Excl, NoiseParam_Exclude);

    //free(IM_S);

    for (unsigned int i = 0; i < CCs2.size(); i++)
    {
	if (Excl[i] < 1)
	    continue;

	for (int y = CCs2[i].ymin; y < CCs2[i].ymax + 1; y++)
	    for (int x = CCs2[i].xmin; x < CCs2[i].xmax + 1; x++)
		if (CCs2[i].label == (UnB2.Get_Labelmap(x, y)))
		    *(NibRes1 + (y - 1) * Ix + x - 1) = 0;
    }

    //FreeMem - Release UnB2 and CCs2
    if (CCs2.size() > 0)
    {
	CCs2.clear();
	vector<cc>().swap(CCs2);
    }
    else
	vector<cc>().swap(CCs2);

    UnB2.Clear();

    //Store data to IMor that is both for input/output
    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    unsigned char val = (*(NibRes1 + y * Ix + x)) % 2;
	    *(IMor + y * Ix + x) = (val + 255) % 256;
		
	}

    
    
    UnB1P1.Clear();
    
    return 0;
}

bool NibSW(float k0, int s0, unsigned char *inputIm, unsigned char *outputImB, int Ix, int Iy)
{
    int dW;
    float k;

    dW = s0; //30
    k = k0;  //-2

    //float m,s;
    //long TH,pix,gray,gray2;

    bool ret = true;

    long *IX_pix, *IX_gray, *IX_graygray, *IX_pix1, *IX_gray1, *IX_graygray1;

    IX_pix = NULL;
    IX_gray = NULL;
    IX_graygray = NULL;
    IX_pix1 = NULL;
    IX_gray1 = NULL;
    IX_graygray1 = NULL;

    IX_pix = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_pix == NULL)
    {
	return ret;
    }

    IX_gray = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_gray == NULL)
    {
	if (IX_pix != NULL)
	    free(IX_pix);
	return ret;
    }

    IX_graygray = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_graygray == NULL)
    {
	if (IX_pix != NULL)
	    free(IX_pix);
	if (IX_gray != NULL)
	    free(IX_gray);
	return ret;
    }

    IX_gray1 = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_gray1 == NULL)
    {
	if (IX_pix != NULL)
	    free(IX_pix);
	if (IX_gray != NULL)
	    free(IX_gray);
	if (IX_graygray != NULL)
	    free(IX_graygray);
	return ret;
    }

    IX_graygray1 = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_graygray1 == NULL)
    {
	if (IX_pix != NULL)
	    free(IX_pix);
	if (IX_gray != NULL)
	    free(IX_gray);
	if (IX_graygray != NULL)
	    free(IX_graygray);
	if (IX_gray1 != NULL)
	    free(IX_gray1);
	return ret;
    }

    IX_pix1 = (long *)calloc((Ix + 2), sizeof(long));
    ret = false;
    if (IX_pix1 == NULL)
    {
	if (IX_pix != NULL)
	    free(IX_pix);
	if (IX_gray != NULL)
	    free(IX_gray);
	if (IX_graygray != NULL)
	    free(IX_graygray);
	if (IX_gray1 != NULL)
	    free(IX_gray1);
	if (IX_graygray1 != NULL)
	    free(IX_graygray1);
	return ret;
    }

    for (int y = 0; y < Iy; y++)
    {
	if (y == 0)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		long pix = 0;
		long gray = 0;
		long gray2 = 0;
		for (int iy = 0; iy <= dW; iy++)
		{
		    pix++;
		    {
			gray += *(inputIm + iy * Ix + x);
			gray2 += (*(inputIm + iy * Ix + x)) * (*(inputIm + iy * Ix + x));
		    }
		}
		IX_pix[x] = pix;
		IX_gray[x] = gray;
		IX_graygray[x] = gray2;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else if (y <= dW)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		IX_pix[x] = IX_pix1[x] + 1;
		{
		    IX_gray[x] = IX_gray1[x] + *(inputIm + (y + dW) * Ix + x);
		    IX_graygray[x] = IX_graygray1[x] + (*(inputIm + (y + dW) * Ix + x)) * (*(inputIm + (y + dW) * Ix + x));
		}
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else if (y <= Iy - 1 - dW)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		IX_pix[x] = IX_pix1[x];
		{
		    IX_gray[x] = IX_gray1[x] + *(inputIm + (y + dW) * Ix + x) - *(inputIm + (y - dW - 1) * Ix + x);
		    IX_graygray[x] = IX_graygray1[x] + (*(inputIm + (y + dW) * Ix + x)) * (*(inputIm + (y + dW) * Ix + x)) - (*(inputIm + (y - dW - 1) * Ix + x)) * (*(inputIm + (y - dW - 1) * Ix + x));
		}
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else
	{
	    for (int x = 0; x < Ix; x++)
	    {
		IX_pix[x] = IX_pix1[x] - 1;
		{
		    IX_gray[x] = IX_gray1[x] - (*(inputIm + (y - dW - 1) * Ix + x));
		    IX_graygray[x] = IX_graygray1[x] - (*(inputIm + (y - dW - 1) * Ix + x)) * (*(inputIm + (y - dW - 1) * Ix + x));
		}
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}

	long pix1, gray1, graygray1;
	for (int x = 0; x < Ix; x++)
	{
	    float m, s;
	    long TH, pix, gray, graygray;

	    if (x == 0)
	    {
		pix = 0;
		gray = 0;
		graygray = 0;
		for (int ix = 0; ix <= dW; ix++)
		{
		    pix += IX_pix[ix];
		    gray += IX_gray[ix];
		    graygray += IX_graygray[ix];
		}
	    }
	    else if (x <= dW)
	    {
		pix = pix1 + IX_pix[x + dW];
		gray = gray1 + IX_gray[x + dW];
		graygray = graygray1 + IX_graygray[x + dW];
	    }
	    else if (x <= Ix - 1 - dW)
	    {
		pix = pix1 + IX_pix[x + dW] - IX_pix[x - dW - 1];
		gray = gray1 + IX_gray[x + dW] - IX_gray[x - dW - 1];
		graygray = graygray1 + IX_graygray[x + dW] - IX_graygray[x - dW - 1];
	    }
	    else
	    {
		pix = pix1 - IX_pix[x - dW - 1];
		gray = gray1 - IX_gray[x - dW - 1];
		graygray = graygray1 - IX_graygray[x - dW - 1];
	    }

	    pix1 = pix;
	    gray1 = gray;
	    graygray1 = graygray;

	    m = gray / pix;
	    s = sqrt(fabs((float)(pix * graygray - gray * gray) / (pix * (pix - 1))));

	    TH = m + (float)(k * s) / 10;
	    {
		if ((*(inputIm + y * Ix + x) >= TH))
		    *(outputImB + y * Ix + x) = 0;
		else
		    *(outputImB + y * Ix + x) = 1;
	    }
	}
    }

    //Release Mem
    if (IX_pix != NULL)
	free(IX_pix);
    if (IX_gray != NULL)
	free(IX_gray);
    if (IX_graygray != NULL)
	free(IX_graygray);
    if (IX_gray1 != NULL)
	free(IX_gray1);
    if (IX_graygray1 != NULL)
	free(IX_graygray1);
    if (IX_pix1 != NULL)
	free(IX_pix1);

    return true;
}

//Niblack adaptive window (slow) - For Printed with various sizes
void NibSW(float k0, int s0, unsigned char *inputIm, unsigned char *outputImB, int *inputSWIm, int Ix, int Iy)
{
    int dW;
    float k;

    dW = s0; //30
    k = k0;  //-2

    float m, s;
    long TH, pix, gray, gray2;

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    dW = *(inputSWIm + y * Ix + x);
	    pix = 0;
	    gray = 0;
	    gray2 = 0;
	    for (int ix = x - dW; ix <= x + dW; ix++)
		for (int iy = y - dW; iy <= y + dW; iy++)
		{
		    if ((ix >= 0) && (iy >= 0) && (ix < Ix) && (iy < Iy))
		    {
			pix++;
			gray += *(inputIm + iy * Ix + ix);
			gray2 += (*(inputIm + iy * Ix + ix)) * (*(inputIm + iy * Ix + ix));
		    }
		}
	    m = gray / pix;
	    s = sqrt(fabs((float)(pix * gray2 - gray * gray) / (pix * (pix - 1))));

	    TH = m + (float)(k * s) / 10;

	    if ((*(inputIm + y * Ix + x)) >= TH)
		*(outputImB + y * Ix + x) = 0;
	    else
		*(outputImB + y * Ix + x) = 1;
	}
}

void Dilationx3(unsigned char *inputImB, unsigned char *outputImB, int Ix, int Iy)
{

    int w = 1;

    unsigned char boxy, chk;

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    boxy = *(inputImB + y * Ix + x);
	    //If BG check for FG in the 8-neighborhood
	    if (boxy < 1)
	    {
		chk = 0;
		for (int b = y - w; b < y + w + 1; b++)
		{
		    for (int a = x - w; a < x + w + 1; a++)
		    {
			if (a < 0 || b < 0 || a > Ix - 1 || b > Iy - 1)
			    continue;
			chk = *(inputImB + b * Ix + a);
			if (chk == 1)
			    break;
		    }
		    if (chk == 1)
			break;
		}
		*(outputImB + y * Ix + x) = chk;
	    }
	    else
		*(outputImB + y * Ix + x) = 1;
	}
}

//Background Estimation using Inpainting
//Noise Parameter
bool Inpaintx4Passes(unsigned char *inputIm, unsigned char *inputMaskB, unsigned char *outputIm, unsigned char *outputImStats, int noizLev, int Ix, int Iy)
{
    unsigned char *InPass1 = NULL;
    unsigned char *InPass2 = NULL;
    unsigned char *InPass3 = NULL;
    unsigned char *InPass4 = NULL;

    InPass1 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (InPass1 == NULL)
	return false;

    InPass2 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (InPass2 == NULL)
    {
	if (InPass1 != NULL)
	    free(InPass1);
	return false;
    }

    InPass3 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (InPass3 == NULL)
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	return false;
    }

    InPass4 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (InPass4 == NULL)
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	if (InPass3 != NULL)
	    free(InPass3);
	return false;
    }

    //Perform 4 Inpainting passed
    if (!InpaintAssistPass(0, inputIm, inputMaskB, InPass1, Ix, Iy))
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	if (InPass3 != NULL)
	    free(InPass3);
	if (InPass4 != NULL)
	    free(InPass4);
	return false;
    }
    if (!InpaintAssistPass(1, inputIm, inputMaskB, InPass2, Ix, Iy))
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	if (InPass3 != NULL)
	    free(InPass3);
	if (InPass4 != NULL)
	    free(InPass4);
	return false;
    }
    if (!InpaintAssistPass(2, inputIm, inputMaskB, InPass3, Ix, Iy))
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	if (InPass3 != NULL)
	    free(InPass3);
	if (InPass4 != NULL)
	    free(InPass4);
	return false;
    }
    if (!InpaintAssistPass(3, inputIm, inputMaskB, InPass4, Ix, Iy))
    {
	if (InPass1 != NULL)
	    free(InPass1);
	if (InPass2 != NULL)
	    free(InPass2);
	if (InPass3 != NULL)
	    free(InPass3);
	if (InPass4 != NULL)
	    free(InPass4);
	return false;
    }

    //Keep the min or average from the inpaint passes
    int minval = 0;
    int aveval = 0;
    int maxval = 0;
    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    CalcStatsInpaint(*(InPass1 + y * Ix + x), *(InPass2 + y * Ix + x), *(InPass3 + y * Ix + x), *(InPass4 + y * Ix + x), &minval, &maxval, &aveval);
	    if (noizLev < 1)
		*(outputIm + y * Ix + x) = minval;
	    else if (noizLev < 2)
		*(outputIm + y * Ix + x) = aveval;
	    else
		*(outputIm + y * Ix + x) = maxval;

	    *(outputImStats + y * Ix + x) = aveval;
	}

    //Release Mem
    if (InPass1 != NULL)
	free(InPass1);
    if (InPass2 != NULL)
	free(InPass2);
    if (InPass3 != NULL)
	free(InPass3);
    if (InPass4 != NULL)
	free(InPass4);

    return true;
}

bool InpaintAssistPass(int whoPass, unsigned char *ImOriginal, unsigned char *InpMaskB, unsigned char *EstBG, int Ix, int Iy)
{
    unsigned char *IM22 = NULL;
    IM22 = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));

    if (IM22 == NULL)
	return false;

    memcpy(EstBG, ImOriginal, Ix * Iy * sizeof(unsigned char));
    memcpy(IM22, InpMaskB, Ix * Iy * sizeof(unsigned char));

    /*for(int y=0;y<Iy;y++)
					for(int x=0;x<Ix;x++)
					{
						*(EstBG+y*Ix+x)=*(ImOriginal+y*Ix+x);
						*(v+y*Ix+x)=*(InpMaskB+y*Ix+x);
					}*/

    if (whoPass == 0)
    {
	for (int i = 0; i < Ix; i++)
	{
	    for (int j = 0; j < Iy; j++)
		if (*(IM22 + j * Ix + i) != 0)
		{
		    *(EstBG + j * Ix + i) = InpaintAssistCheck4(i, j, IM22, EstBG, Ix, Iy);
		    *(IM22 + j * Ix + i) = 0;
		}
	}
    }
    else if (whoPass == 1)
    {
	for (int i = 0; i < Ix; i++)
	{
	    for (int j = Iy - 1; j > -1; j--)
		if (*(IM22 + j * Ix + i) != 0)
		{
		    *(EstBG + j * Ix + i) = InpaintAssistCheck4(i, j, IM22, EstBG, Ix, Iy);
		    *(IM22 + j * Ix + i) = 0;
		}
	}
    }
    else if (whoPass == 2)
    {
	for (int i = Ix - 1; i > -1; i--)
	{
	    for (int j = 0; j < Iy; j++)
		if (*(IM22 + j * Ix + i) != 0)
		{
		    *(EstBG + j * Ix + i) = InpaintAssistCheck4(i, j, IM22, EstBG, Ix, Iy);
		    *(IM22 + j * Ix + i) = 0;
		}
	}
    }
    else if (whoPass == 3)
    {
	for (int i = Ix - 1; i > -1; i--)
	{
	    for (int j = Iy - 1; j > -1; j--)
		if (*(IM22 + j * Ix + i) != 0)
		{
		    *(EstBG + j * Ix + i) = InpaintAssistCheck4(i, j, IM22, EstBG, Ix, Iy);
		    *(IM22 + j * Ix + i) = 0;
		}
	}
    }

    //Release Mem
    if (IM22 != NULL)
	free(IM22);

    return true;
}

int InpaintAssistCheck4(int x, int y, unsigned char *imB, unsigned char *out, int Ix, int Iy)
{
    int val[4];

    val[1] = val[2] = val[3] = val[0] = -1;

    if (x > 0)
	if (*(imB + y * Ix + x - 1) != 1)
	    val[1] = *(out + y * Ix + x - 1);

    if (y > 0)
	if (*(imB + (y - 1) * Ix + x) != 1)
	    val[0] = *(out + (y - 1) * Ix + x);

    if (x < Ix - 1)
	if (*(imB + y * Ix + x + 1) != 1)
	    val[3] = *(out + y * Ix + x + 1);

    if (y < Iy - 1)
	if (*(imB + (y + 1) * Ix + x) != 1)
	    val[2] = *(out + (y + 1) * Ix + x);

    int cnt = 0;

    double value = 0.0;
    for (int k = 0; k < 4; k++)
	if (val[k] > 0) //-1
	{
	    cnt = cnt + 1;
	    value = value + val[k];
	}

    //Mesi timi apo tis 4 (egkyres) times
    if (cnt > 0)
	value = value / double(cnt);
    else
	value = 255;

    return (int)value;
}

void CalcStatsInpaint(unsigned char Fst, unsigned char Snd, unsigned char Trd, unsigned char Frth, int *minn, int *maxx, int *avee)
{
    unsigned char val[4];
    val[0] = Fst;
    val[1] = Snd;
    val[2] = Trd;
    val[3] = Frth;

    int min, ave, max, cnt;
    min = 1000;
    max = -1;
    cnt = 0;
    ave = 0;
    for (int i = 0; i < 4; i++)
    {
	if (val[i] < min)
	    min = val[i];
	if (val[i] > max && val[i] < 255)
	    max = val[i];
	if (val[i] < 255)
	{
	    ave += val[i];
	    cnt++;
	}
    }

    (*minn) = min;
    (*maxx) = max;

    double retave;
    if (cnt != 0)
	retave = (double)ave / (double)cnt;
    retave += (double)max;
    retave /= 2;

    *(avee) = (int)retave;
}

bool ImageNormalization(unsigned char *inputIm, unsigned char *ImDevisor, unsigned char *outputIm, int Ix, int Iy)
{
    //Ypologismos min kai max tis Original eikonas
    double MnOr = 1000;
    double MxOr = -1;
    double dif;
    int valOr;
    for (int b = 0; b < Iy; b++)
	for (int a = 0; a < Ix; a++)
	{
	    valOr = *(inputIm + b * Ix + a);
	    if (valOr < MnOr)
		MnOr = valOr;
	    if (valOr > MxOr)
		MxOr = valOr;
	}

    dif = MxOr - MnOr;

    //Normalization (Original/BGEstimation)
    //and min max calculation
    double *NormD = NULL;
    NormD = (double *)calloc(Ix * Iy, sizeof(double));
    if (NormD == NULL)
	return false;

    double MnNorm = 1000;
    double MxNorm = -1;
    double valNormD;
    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    int vohelp = *(inputIm + y * Ix + x);
	    int vbhelp = *(ImDevisor + y * Ix + x);
	    if (vbhelp > 0)
	    {
		valNormD = (double)((double)vohelp / (double)vbhelp);
		if (valNormD > 1.1 && (vbhelp < (MnOr + dif / 5)))
		    valNormD = 1.1;
	    }
	    else
	    {
		valNormD = 0;
	    }
	    *(NormD + y * Ix + x) = valNormD;
	    if (valNormD < MnNorm)
		MnNorm = valNormD;
	    if (valNormD > MxNorm)
		MxNorm = valNormD;
	}

    double valNorm;
    for (int d = 0; d < Iy; d++)
	for (int c = 0; c < Ix; c++)
	{
	    valNormD = *(NormD + d * Ix + c);
	    valNorm = ((dif * (valNormD - MnNorm)) / (MxNorm - MnNorm)) + MnOr;
	    *(outputIm + d * Ix + c) = (int)(valNorm + 0.5);
	}

    //Release Mem
    if (NormD != NULL)
	free(NormD);

    return true;
}

int PostProcOtsu(int **inputIm, vector<cc> rCCs, unsigned char **outputIm, int Ix, int Iy)
{
    int *HistP1 = NULL;
    HistP1 = (int *)calloc(Iy, sizeof(int));

    if (HistP1 == NULL)
	return -1;

    int plithos = rCCs.size();

    bool vrika = false;
    int ena = 0;
    int dyo = 0;
    int metrw = 0;
    int afairw = 0;
    int summ = 0;
    bool vrika2 = false;
    bool vrika3 = false;
    int wheresum = -1;

    for (int who = 0; (who < Iy && !vrika3); who++) // /supHe
    {
	double allpixels = 0;
	double pixelsccs = 0;
	double ccs = 0;
	double allccs = 0;
	PostProcOtsuAssist(who, plithos, rCCs, inputIm, &ccs, &allccs, &pixelsccs, &allpixels);

	double pososto;
	double pososto2, d2;

	pososto = pixelsccs / allpixels;
	pososto *= 100;

	pososto2 = ccs / allccs;
	pososto2 *= 100;
	if (ccs < 1)
	    d2 = 0;
	else
	    d2 = pososto / pososto2;

	int num2 = 0;
	num2 = (int)(d2 * 100);

	HistP1[who] = num2;

	int j = who;

	if (HistP1[j] > 0)
	{
	    if (!vrika)
	    {
		ena = j;
		vrika = true;
	    }
	    dyo = j;
	    afairw = dyo - ena - 1;
	    metrw += afairw;
	    ena = dyo;
	    summ += HistP1[j];
	    if (summ > 100 && !vrika3)
	    {
		wheresum = j;
		vrika3 = true;
	    }
	    if (!vrika2)
		if (HistP1[j] >= 100)
		{
		    vrika2 = true;
		}
	}
    }

    if (wheresum > (Iy / 33))
	wheresum = Iy / 50;

    //Exclude CCs (Keep the Correct CCs)
    for (int i = 0; i < plithos; i++)
    {
	int hei = abs(rCCs[i].ymax - rCCs[i].ymin);
	if (hei <= wheresum)
	    continue;
	for (int y = rCCs[i].ymin; y < rCCs[i].ymax + 1; y++)
	    for (int x = rCCs[i].xmin; x < rCCs[i].xmax + 1; x++)
		if (rCCs[i].label == inputIm[y][x])
		    outputIm[y][x] = 1;
    }

    if (HistP1 != NULL)
    {
	free(HistP1);
	HistP1 = NULL;
    }

    return wheresum;
}

void PostProcOtsuAssist(int who, int plithos, vector<cc> rCCs, int **UnB, double *cccnt, double *ccall, double *pixelscnt, double *pixelsall)
{
    double allpixels = 0;
    double pixels;
    double pixelsccs = 0;
    double ccs = 0;
    double ccsall = 0;
    int i, x, y, h;
    for (i = 0; i < plithos; i++)
    {
	pixels = 0;
	ccsall += 1;
	h = abs(rCCs[i].ymax - rCCs[i].ymin);
	for (y = rCCs[i].ymin; y < rCCs[i].ymax + 1; y++)
	{
	    for (x = rCCs[i].xmin; x < rCCs[i].xmax + 1; x++)
	    {
		if (rCCs[i].label == UnB[y][x])
		{
		    allpixels += 1;
		    pixels += 1;
		}
	    }
	}
	if (h == who)
	{
	    ccs += 1;
	    pixelsccs += pixels;
	}
    }

    *pixelscnt = pixelsccs;
    *cccnt = ccs;
    *pixelsall = allpixels;
    *ccall = ccsall;
}

void boundaries(int Rccid, vector<cc> rCCs, int **labels, unsigned char *outputContour, int Ix, int Iy)
{
    int x, y, a, b;
    bool on;

    for (y = rCCs[Rccid].ymin; y < rCCs[Rccid].ymax + 1; y++)
    {
	on = false;
	int yy;
	int xx;
	for (x = rCCs[Rccid].xmin; x < rCCs[Rccid].xmax + 1; x++)
	{
	    yy = y - 1;
	    xx = x - 1;
	    if (rCCs[Rccid].label == labels[y][x])
	    {
		if (!on)
		{
		    *(outputContour + yy * Ix + xx) = 1;
		    on = true;
		}
	    }
	    else
	    {
		if (on)
		{
		    *(outputContour + yy * Ix + xx - 1) = 1;
		}
		on = false;
	    }
	}
	xx = x - 1;
	if (on)
	{
	    *(outputContour + yy * Ix + xx - 1) = 1;
	}
    }

    for (a = rCCs[Rccid].xmin; a < rCCs[Rccid].xmax + 1; a++)
    {
	on = false;
	int aa, bb;
	for (b = rCCs[Rccid].ymin; b < rCCs[Rccid].ymax + 1; b++)
	{
	    aa = a - 1;
	    bb = b - 1;
	    if (rCCs[Rccid].label == labels[b][a])
	    {
		if (!on)
		{
		    *(outputContour + bb * Ix + aa) = 1;
		    on = true;
		}
	    }
	    else
	    {
		if (on)
		{
		    *(outputContour + (bb - 1) * Ix + aa) = 1;
		}
		on = false;
	    }
	}
	bb = b - 1;
	if (on)
	{
	    *(outputContour + (bb - 1) * Ix + aa) = 1;
	}
    }
}

int CalcSW(vector<cc> SkCCs, int **Sklabels, unsigned char *Contours, int *InSkelWidth, int Ix, int Iy)
{
    int swH;
    double aveSW = 0.0;
    double cntSW = 0;
    int maxSW;
    int Ex = 1;

    for (unsigned int j = 0; j < SkCCs.size(); j++)
    {
	maxSW = -1;
	for (int n = SkCCs[j].ymin; n < SkCCs[j].ymax + 1; n++)
	    for (int m = SkCCs[j].xmin; m < SkCCs[j].xmax + 1; m++)
	    {
		if (SkCCs[j].label == Sklabels[n][m])
		    if (*(InSkelWidth + (n - Ex) * Ix + m - Ex) > 0 && *(InSkelWidth + (n - Ex) * Ix + m - Ex) > maxSW)
			maxSW = *(InSkelWidth + (n - Ex) * Ix + m - Ex);
	    }
	aveSW += maxSW;
	cntSW += 1;
    }

    swH = (int)((aveSW / cntSW));

    return swH;
}

void CalcSWAssist(vector<cc> BCCs, int **Blabels, vector<cc> SkCCs, int **Sklabels, unsigned char *Contours, int *OutSkelWidth, int Ix, int Iy)
{
    for (unsigned int j = 0; j < SkCCs.size(); j++)
	for (int n = SkCCs[j].ymin; n < SkCCs[j].ymax + 1; n++)
	    for (int m = SkCCs[j].xmin; m < SkCCs[j].xmax + 1; m++)
		if (SkCCs[j].label == Sklabels[n][m])
		    *(OutSkelWidth + (n - 1) * Ix + m - 1) = CalcSWAssistNearCont((Blabels[n][m]) - 1, BCCs, Contours, m - 1, n - 1, Ix, Iy);
}

int CalcSWAssistNearCont(int who, vector<cc> BCCs, unsigned char *Contours, int x, int y, int Ix, int Iy)
{
    int ret, a, b, c, d, k, l;
    int cnt = 0;
    bool fnd;
    int m, n;
    int xxmin, yymin, xxmax, yymax;
    int Ex = 1;

    xxmin = BCCs[who].xmin - Ex;
    xxmax = BCCs[who].xmax - Ex;
    yymin = BCCs[who].ymin - Ex;
    yymax = BCCs[who].ymax - Ex;

    if (*(Contours + y * Ix + x) < 1)
    {
	do
	{
	    cnt++;
	    fnd = false;
	    a = x - cnt;
	    b = x + cnt;
	    c = y - cnt;
	    d = y + cnt;
	    if (a < xxmin)
		a = xxmin;
	    if (b > xxmax)
		b = xxmax;
	    if (c < yymin)
		c = yymin;
	    if (d > yymax)
		d = yymax;

	    k = a;
	    for (l = c + 1; l < d; l++)
	    {
		if (*(Contours + l * Ix + k) != 0)
		{
		    fnd = true;
		}
	    }

	    m = b;
	    for (n = c + 1; n < d; n++)
	    {
		if (*(Contours + n * Ix + m) != 0)
		{
		    fnd = true;
		}
	    }

	    k = c;
	    for (l = a + 1; l < b; l++)
	    {
		if (*(Contours + k * Ix + l) != 0)
		{
		    fnd = true;
		}
	    }

	    m = d;
	    for (n = a + 1; n < b; n++)
	    {
		if (*(Contours + m * Ix + n) != 0)
		{
		    fnd = true;
		}
	    }

	    if (*(Contours + c * Ix + a) != 0)
	    {
		fnd = true;
	    }

	    if (*(Contours + d * Ix + a) != 0)
	    {
		fnd = true;
	    }

	    if (*(Contours + c * Ix + b) != 0)
	    {
		fnd = true;
	    }

	    if (*(Contours + d * Ix + b) != 0)
	    {
		fnd = true;
	    }
	} while (!fnd);
    }
    else
    {
	ret = 1;
	return ret;
    }

    ret = 2 * cnt + 1;

    return ret;
}

void CombineBins(unsigned char **inputImB1, vector<cc> BCCs2, int **Blabels2, double LContrast, unsigned char *Excluding, int noisL)
{
    //Exclude (1) CCs that do not have a match and (2) big CCs

    if (noisL < 1)
    {
	bool foundi = false;
	long koina, ola;
	for (unsigned int i = 0; i < BCCs2.size(); i++)
	{
	    foundi = false;
	    koina = 0;
	    ola = 0;
	    for (int y = BCCs2[i].ymin; y < BCCs2[i].ymax + 1; y++)
		for (int x = BCCs2[i].xmin; x < BCCs2[i].xmax + 1; x++)
		    if (BCCs2[i].label == (Blabels2[y][x]))
		    {
			ola++;
			if ((inputImB1[y][x]) != 0)
			{
			    foundi = true;
			    koina++;
			}
		    }

	    if (!foundi)
		Excluding[i] = 1;
	    else
	    {
		double epitis100 = double(koina) / double(ola);
		epitis100 = 100.0 * epitis100;
		if (epitis100 < LContrast)
		    Excluding[i] = 2;
	    }
	}
    }
    else
    {
	bool foundi = false;
	for (unsigned int i = 0; i < BCCs2.size(); i++)
	{
	    foundi = false;
	    for (int y = BCCs2[i].ymin; y < BCCs2[i].ymax + 1 && !foundi; y++)
		for (int x = BCCs2[i].xmin; x < BCCs2[i].xmax + 1; x++)
		    if (BCCs2[i].label == (Blabels2[y][x]))
			if ((inputImB1[y][x]) != 0)
			    foundi = true;

	    if (!foundi)
		Excluding[i] = 1;
	}
    }
}

bool Enhance(unsigned char *inputIm, unsigned char *ioImB, unsigned char **inputImB1P, unsigned char *Contours, unsigned char *Excluding, double LContrast, int Ix, int Iy)
{
    int Ex = 1;

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	    if (*(Contours + y * Ix + x) != 0)
		if (inputImB1P[y + Ex][x + Ex] != 0)
		    check_NIB(x, y, ioImB, Ix, Iy);

    /*for(int yy=0;yy<Iy;yy++)
		for(int xx=0;xx<Ix;xx++)
			if(*(ioImB+yy*Ix+xx)==1)
				check_OTS(xx,yy,ioImB,inputImB1P);*/

    Mat cvOriginal, Edg;
    cvOriginal.data = NULL;
    Edg.data = NULL;

    cvOriginal.create(Size(Ix, Iy), CV_8UC1);
    if (cvOriginal.data == NULL)
	return false;

    Edg = cvOriginal.clone();
    if (Edg.data == NULL)
    {
	cvOriginal.release();
	return false;
    }

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	    cvOriginal.at<uchar>(y, x) = *(inputIm + y * Ix + x);

    float th = LContrast * 10;
    float p = 2;

    cv::Canny(cvOriginal, Edg, th / p, th, 3, false);
    //Canny(cvOriginal,Edg,th/p,th,3,false);

    for (int yy = 0; yy < Iy; yy++)
	for (int xx = 0; xx < Ix; xx++)
	    if ((*(ioImB + yy * Ix + xx)) != 0 && (*(ioImB + yy * Ix + xx)) < 4)
		check_EDG(xx, yy, ioImB, Edg, Ix, Iy);

    cvOriginal.release();
    Edg.release();

    return true;
}

void check_NIB(int xx, int yy, unsigned char *ioImB, int Ix, int Iy)
{
    int x, y, sx, sy, ex, ey;
    int count = 0;

    sx = xx - 1;
    sy = yy - 1;
    ex = xx + 1;
    ey = yy + 1;

    for (y = sy; y < ey + 1; y++)
	for (x = sx; x < ex + 1; x++)
	    if (y < 0 || x < 0 || y > Iy - 1 || x > Ix - 1)
		;
	    else
	    {
		if (*(ioImB + y * Ix + x) != 0 && *(ioImB + y * Ix + x) < 2)
		{
		    count++;
		    if (count > 2)
		    {
			*(ioImB + y * Ix + x) = 3;
			return;
		    }
		}
		else
		    ;
	    }
}

void check_OTS(int xx, int yy, unsigned char *ioImB, unsigned char **inputImB1P, int Ix, int Iy)
{
    int x, y, sx, sy, ex, ey;
    sx = xx - 1;
    sy = yy - 1;
    ex = xx + 1;
    ey = yy + 1;
    int Ex = 1;

    for (y = sy; y < ey + 1; y++)
	for (x = sx; x < ex + 1; x++)
	    if (y < 0 || x < 0 || y > Iy - 1 || x > Ix - 1)
		;
	    else
	    {
		if (inputImB1P[y + Ex][x + Ex] != 0)
		    *(ioImB + y * Ix + x) = 1;
		else
		    ;
	    }
}

void check_EDG(int xx, int yy, unsigned char *ioImB, cv::Mat inMat, int Ix, int Iy)
{
    int x, y, sx, sy, ex, ey;
    sx = xx - 1;
    sy = yy - 1;
    ex = xx + 1;
    ey = yy + 1;

    int dat;
    int pou = 0;

    dat = (int)(inMat.at<uchar>(yy, xx));
    if (dat > 0)
	return;

    for (y = sy; y < ey + 1; y++)
	for (x = sx; x < ex + 1; x++)
	{
	    pou++;
	    if (pou % 2 != 1)
	    {
		if (y < 0 || x < 0 || y > Iy - 1 || x > Ix - 1)
		    ;
		else
		{
		    dat = (int)(inMat.at<uchar>(y, x));
		    if (dat > 0)
			*(ioImB + y * Ix + x) = 5;
		    else
			;
		}
	    }
	}
}

bool Filter(unsigned char *IM_S, int Ix, int Iy)
{
    unsigned char *flt = NULL;
    int d = 1; //3x3
    flt = (unsigned char *)calloc(Ix * Iy, sizeof(unsigned char));
    if (flt == NULL)
	return false;

    float filt[25];
    if (d == 1)
    {
	filt[0] = 0.0751;
	filt[1] = 0.1238;
	filt[2] = 0.0751;
	filt[3] = 0.1238;
	filt[4] = 0.2042;
	filt[5] = 0.1238;
	filt[6] = 0.0751;
	filt[7] = 0.1238;
	filt[8] = 0.0751;
    }
    else
    {
	filt[0] = 0.0030;
	filt[1] = 0.0133;
	filt[2] = 0.0219;
	filt[3] = 0.0133;
	filt[4] = 0.0030;
	filt[5] = 0.0133;
	filt[6] = 0.0596;
	filt[7] = 0.0983;
	filt[8] = 0.0596;
	filt[9] = 0.0133;
	filt[10] = 0.0219;
	filt[11] = 0.0983;
	filt[12] = 0.1621;
	filt[13] = 0.0983;
	filt[14] = 0.0219;
	filt[15] = 0.0133;
	filt[16] = 0.0596;
	filt[17] = 0.0983;
	filt[18] = 0.0596;
	filt[19] = 0.0133;
	filt[20] = 0.0030;
	filt[21] = 0.0133;
	filt[22] = 0.0219;
	filt[23] = 0.0133;
	filt[24] = 0.0030;
    }
    for (int ky = d; ky < Iy - d; ky++)
	for (int kx = d; kx < Ix - d; kx++)
	{
	    float sum = 0;
	    int cnt = 0;
	    for (int i = -d; i <= d; i++)
		for (int j = -d; j <= d; j++)
		{
		    int kk = *(IM_S + (ky + j) * Ix + kx + i);
		    sum += ((float)kk * filt[cnt]);
		    cnt++;
		}
	    *(flt + ky * Ix + kx) = sum;
	}

    for (int y = d; y < Iy - d; y++)
	for (int x = d; x < Ix - d; x++)
	    *(IM_S + y * Ix + x) = *(flt + y * Ix + x);

    //        Sleep(10);
    free(flt);
    flt = NULL;

    return true;
}

bool GPP_no_dll(unsigned char *IM_S, int Ix, int Iy)
{
    float inp_k = 5;
    float inp_q = 1;
    float inp_p2 = 10;

    //float inp_q_f = (float)inp_q/10;
    float inp_p2_f = (float)inp_p2 / 10;

    //params
    int dW = 20;
    float k = inp_k; //= 1;
    k = k / 100;
    int R = 128;
    float q = (float)1 / inp_q; //(float)1/1.6;
    float p1 = 0.5;
    float p2 = inp_p2_f; //0.8;
    int upsampling = 1;
    int dW1 = 20;

    //int RH[256],GH[256],BH[256];
    unsigned char *IMAGE11;
    unsigned char I, I1;
    long pix, gray, gray2, TH, graygray, pix1, gray1;
    long ydWIx, ydW_1Ix, yIx, ydW1Ix, ydW1_1Ix;
    float m, s;
    long d, d2;

    float a;
    float b;
    // int xx;
    float aa, bb, cc, dd;

    unsigned char *IMAGE;

    //for (int i=0;i<256;i++) {RH[i]=0;GH[i]=0;BH[i]=0;}

    long *IX_pix, *IX_gray, *IX_graygray, *IX_pix1, *IX_gray1, *IX_graygray1;

    IMAGE = IMAGE11 = NULL;
    IX_pix = IX_gray = IX_graygray = IX_pix1 = IX_gray1 = IX_graygray1 = NULL;

    IMAGE = (unsigned char *)malloc((Ix + 1) * (Iy + 1) * sizeof(unsigned char));

    IMAGE11 = (unsigned char *)malloc((Ix + 1) * (Iy + 1) * sizeof(unsigned char));

    for (int y = 0; y < Iy; y++)
	for (int x = 0; x < Ix; x++)
	{
	    IMAGE[y * Ix + x] = *(IM_S + y * Ix + x);
	    IMAGE11[y * Ix + x] = *(IM_S + y * Ix + x);
	}

    IX_pix = (long *)malloc((Ix + 1) * sizeof(long));
    IX_gray = (long *)malloc((Ix + 1) * sizeof(long));
    IX_graygray = (long *)malloc((Ix + 1) * sizeof(long));
    IX_pix1 = (long *)malloc((Ix + 1) * sizeof(long));
    IX_gray1 = (long *)malloc((Ix + 1) * sizeof(long));
    IX_graygray1 = (long *)malloc((Ix + 1) * sizeof(long));

    //1st step
    for (int y = 0; y < Iy; y++)
    {
	ydWIx = (y + dW) * Ix;
	ydW_1Ix = (y - dW - 1) * Ix;
	yIx = y * Ix;

	if (y == 0)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		pix = 0;
		gray = 0;
		gray2 = 0;
		for (int iy = 0; iy <= dW; iy++)
		{
		    I = *(IMAGE + iy * Ix + x);
		    gray += I;
		    gray2 += I * I;
		}
		pix += dW + 1;
		IX_pix[x] = pix;
		IX_gray[x] = gray;
		IX_graygray[x] = gray2;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else if (y <= dW)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		I = *(IMAGE + ydWIx + x);
		IX_pix[x] = IX_pix1[x] + 1;
		IX_gray[x] = IX_gray1[x] + I;
		IX_graygray[x] = IX_graygray1[x] + I * I;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else if (y <= Iy - 1 - dW)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		I = *(IMAGE + ydWIx + x);
		I1 = *(IMAGE + ydW_1Ix + x);
		IX_pix[x] = IX_pix1[x];
		IX_gray[x] = IX_gray1[x] + I - I1;
		;
		IX_graygray[x] = IX_graygray1[x] + I * I - I1 * I1;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}
	else
	{
	    for (int x = 0; x < Ix; x++)
	    {
		I = *(IMAGE + ydW_1Ix + x);
		IX_pix[x] = IX_pix1[x] - 1;
		IX_gray[x] = IX_gray1[x] - I;
		IX_graygray[x] = IX_graygray1[x] - I * I;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
		IX_graygray1[x] = IX_graygray[x];
	    }
	}

	long graygray1;
	for (int x = 0; x < Ix; x++)
	{
	    if (x == 0)
	    {
		pix = 0;
		gray = 0;
		graygray = 0;
		for (int ix = 0; ix <= dW; ix++)
		{
		    pix += IX_pix[ix];
		    gray += IX_gray[ix];
		    graygray += IX_graygray[ix];
		}
	    }
	    else if (x <= dW)
	    {
		pix = pix1 + IX_pix[x + dW];
		gray = gray1 + IX_gray[x + dW];
		graygray = graygray1 + IX_graygray[x + dW];
	    }
	    else if (x <= Ix - 1 - dW)
	    {
		pix = pix1 + IX_pix[x + dW] - IX_pix[x - dW - 1];
		gray = gray1 + IX_gray[x + dW] - IX_gray[x - dW - 1];
		graygray = graygray1 + IX_graygray[x + dW] - IX_graygray[x - dW - 1];
	    }
	    else
	    {
		pix = pix1 - IX_pix[x - dW - 1];
		gray = gray1 - IX_gray[x - dW - 1];
		graygray = graygray1 - IX_graygray[x - dW - 1];
	    }

	    pix1 = pix;
	    gray1 = gray;
	    graygray1 = graygray;

	    m = gray / pix;
	    s = sqrt(fabs((float)(pix * graygray - gray * gray) / (pix * (pix - 1))));

	    //if (CheckBox2->Checked)   TH = m+k*s;
	    //else
	    TH = m * (1 - k * (1 - (float)s / R));

	    //if (TH<ThrMin) TH=ThrMin;
	    if ((*(IMAGE + yIx + x) >= TH))
		*(IM_S + yIx + x) = 255;
	    else
		*(IM_S + yIx + x) = 0;
	}
    }
    //if (CheckBox1->Checked) return;
    //Interpolation

    for (int y = 0; y < Iy; y++)
    {
	ydW1Ix = (y + dW1) * Ix;
	ydW1_1Ix = (y - dW1 - 1) * Ix;
	yIx = y * Ix;

	if (y == 0)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		int pix = 0;
		int gray = 0;
		for (int iy = 0; iy <= dW1; iy++)
		{
		    if (*(IM_S + iy * Ix + x) == 255)
		    {
			pix++;
			gray += *(IMAGE + iy * Ix + x);
		    }
		}
		IX_pix[x] = pix;
		IX_gray[x] = gray;
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
	    }
	}
	else if (y <= dW1)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		if (*(IM_S + ydW1Ix + x) == 255)
		{
		    IX_pix[x] = IX_pix1[x] + 1;
		    IX_gray[x] = IX_gray1[x] + *(IMAGE + ydW1Ix + x);
		}
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
	    }
	}
	else if (y <= Iy - 1 - dW1)
	{
	    for (int x = 0; x < Ix; x++)
	    {
		IX_pix[x] = IX_pix1[x];
		IX_gray[x] = IX_gray1[x];

		if (*(IM_S + ydW1Ix + x) == 255)
		{
		    IX_pix[x]++;
		    IX_gray[x] += *(IMAGE + ydW1Ix + x);
		}
		if (*(IM_S + ydW1_1Ix + x) == 255)
		{
		    IX_pix[x]--;
		    IX_gray[x] -= (*(IMAGE + ydW1_1Ix + x));
		}

		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
	    }
	}
	else
	{
	    for (int x = 0; x < Ix; x++)
	    {
		if (*(IM_S + ydW1_1Ix + x) == 255)
		{
		    IX_pix[x] = IX_pix1[x] - 1;
		    IX_gray[x] = IX_gray1[x] - (*(IMAGE + ydW1_1Ix + x));
		}
		IX_pix1[x] = IX_pix[x];
		IX_gray1[x] = IX_gray[x];
	    }
	}

	for (int x = 0; x < Ix; x++)
	{

	    pix = pix1;
	    gray = gray1;
	    if (x == 0)
	    {
		pix = 0;
		gray = 0;
		for (int ix = 0; ix <= dW1; ix++)
		{
		    pix += IX_pix[ix];
		    gray += IX_gray[ix];
		}
	    }
	    else if (x <= dW1)
	    {
		pix = pix1 + IX_pix[x + dW1];
		gray = gray1 + IX_gray[x + dW1];
	    }
	    else if (x <= Ix - 1 - dW1)
	    {
		pix = pix1 + IX_pix[x + dW1] - IX_pix[x - dW1 - 1];
		gray = gray1 + IX_gray[x + dW1] - IX_gray[x - dW1 - 1];
	    }
	    else
	    {
		pix = pix1 - IX_pix[x - dW1 - 1];
		gray = gray1 - IX_gray[x - dW1 - 1];
	    }

	    pix1 = pix;
	    gray1 = gray;

	    //NEW
	    if (pix == 0)
		*(IMAGE11 + yIx + x) = 255;
	    else if (*(IM_S + yIx + x) == 0)
		*(IMAGE11 + yIx + x) = (float)gray / pix;
	    else
		*(IMAGE11 + yIx + x) = *(IMAGE + yIx + x);
	}
    }

    //Thresholding

    long PixFor = 0;
    long D = 0, D1 = 0, D2 = 0;
    int Hist[256];
    for (int i = 0; i < 256; i++)
	Hist[i] = 0;

    for (int y = 0; y < Iy; y++)
    {
	int yIx = y * Ix;
	for (int x = 0; x < Ix; x++)
	{

	    I = *(IMAGE + yIx + x);
	    I1 = *(IMAGE11 + yIx + x);
	    if ((*(IM_S + yIx + x) == 0) && (I1 > I))
	    {
		PixFor++;
		D = D + I1 - I;
		D1 = D1 + I;
		D2 = D2 + I1;
		Hist[I]++;
	    }
	}
    }
    D = D / PixFor;
    D1 = D1 / PixFor;
    D2 = D2 / PixFor;

    long HistMax = 0;
    int H;

    for (int i = 0; i < 256; i++)
    {
	if (Hist[i] > HistMax)
	{
	    HistMax = Hist[i];
	    H = i;
	}
    }
    H += 5;

    d = D;
    d2 = D2;

    a = ((float)d * (p2 - 1)) / ((float)d2 * q * (p1 - 1));
    b = ((float)d * (p1 - p2)) / ((float)q * (p1 - 1));
    aa = q * d * (1 - p2);
    bb = -(2 * (1 + p1)) / (1 - p1);
    cc = p2 * d * q;
    dd = -4 / (d2 * (1 - p1));

    TH = (float)D * 0.4;
    if (upsampling == 1)
    {
	for (int y = 0; y < Iy; y++)
	{
	    yIx = y * Ix;
	    for (int x = 0; x < Ix; x++)
	    {
		//if ((*(IMAGE+yIx+x)<=H) && (*(IMAGE2+yIx+x)==0))   //new
		//  *(IMAGE3+yIx+x)=0;
		//else
		if (*(IM_S + yIx + x) == 0)
		{
		    TH = aa / (1 + exp(dd * (*(IMAGE + yIx + x)) - bb)) + cc;

		    if (*(IMAGE11 + yIx + x) - (*(IMAGE + yIx + x)) > TH)
		    {
			*(IM_S + yIx + x) = 0; //cc
		    }
		    else
		    {
			*(IM_S + yIx + x) = 255; //cc
		    }
		}
		else
		{
		    *(IM_S + yIx + x) = 255; //cc
		}
	    }
	}
    }

    free(IMAGE);
    IMAGE = NULL;
    free(IMAGE11);
    IMAGE11 = NULL;
    free(IX_pix);
    free(IX_gray);
    free(IX_graygray);
    free(IX_pix1);
    free(IX_gray1);
    free(IX_graygray1);

    return true;
}

bool Closing(unsigned char *IM_S, int Ix, int Iy)
{
    //Dilation
    BinaryImage UnD;
    int che = 0;
    che = UnD.Create_Bitmap(Ix, Iy);
    if (che < 0)
	return false;

    for (int d2 = 2; d2 < Iy - 2; d2++)
	for (int d1 = 2; d1 < Ix - 2; d1++)
	{
	    if (*(IM_S + d2 * Ix + d1) > 0)
		continue;
	    for (int f2 = -1; f2 < 2; f2++)
		for (int f1 = -1; f1 < 2; f1++)
		    UnD.Set_Bitmap(d1 + f1, d2 + f2, 1);
	}

    //Erosion
    for (int e2 = 1; e2 < Iy - 1; e2++)
	for (int e1 = 1; e1 < Ix - 1; e1++)
	{
	    if (UnD.Get_Bitmap(e1, e2) < 1)
		continue;

	    int sm = 0;
	    for (int g2 = -1; g2 < 2; g2++)
		for (int g1 = -1; g1 < 2; g1++)
		{
		    sm += UnD.Get_Bitmap(e1 + g1, e2 + g2);
		}

	    if (sm < 9)
		*(IM_S + e2 * Ix + e1) = 255;
	    else
		*(IM_S + e2 * Ix + e1) = 0;
	}
    UnD.Clear();

    return true;
}
