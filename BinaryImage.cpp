#include "BinaryImage.h"
#include <assert.h>

int SearchDirection_8[8][2] = {{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{-1,1}};
int NextTracingDirection_8[8] = {7,7,1,1,3,3,5,5};
int SearchDirection_4[4][2] = {{0,1},{1,0},{0,-1},{-1,0}};
int NextTracingDirection_4[4] = {3,0,1,2};

static unsigned char lc_look_up[16][16]=
/* 0 1 2 3 4 5 6 7 8 9 A B C D E F */
{{255, 255, 255, 01, 255, 00, 01, 01, 255, 255, 255, 255, 01, 00, 00, 01}, /* 0 */
 {255, 255, 255, 255, 00, 255, 00, 00, 01, 255, 255, 255, 01, 00, 01, 01}, /* 1 */
 {255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255}, /* 2 */ 
 {01, 255, 255, 255, 00, 00, 00, 00, 01, 255, 255, 255, 01, 00, 01, 01}, /* 3 */
 {255, 00, 255, 00, 255, 0, 255, 00, 255, 255, 255, 255, 255, 00, 255, 00}, /* 4 */
 {00, 255, 255, 00, 0, 255, 00, 255, 00, 00, 255, 00, 00, 255, 00, 255}, /* 5 */
 {01, 00, 255, 00, 255, 00, 255, 00, 255, 255, 255, 255, 255, 00, 255, 00}, /* 6 */ 
 {01, 00, 255, 00, 00, 255, 00, 255, 01, 00, 255, 00, 01, 255, 01, 255}, /* 7 */ 
 {255, 01, 255, 01, 255, 00, 255, 01, 255, 255, 255, 255, 255, 00, 255, 01}, /* 8 */ 
 {255, 255, 255, 255, 255, 00, 255, 00, 255, 255, 255, 255, 255, 00, 255, 01}, /* 9 */
 {255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255}, /* A */ 
 {255, 255, 255, 255, 255, 00, 255, 00, 255, 255, 255, 255, 255, 00, 255, 00}, /* B */
 {01, 01, 255, 01, 255, 00, 255, 01, 255, 255, 255, 255, 255, 00, 255, 01}, /* C */
 {00, 00, 255, 00, 00, 255, 00, 255, 00, 00, 255, 00, 00, 255, 00, 255}, /* D */
 {00, 01, 255, 01, 255, 00, 255, 01, 255, 255, 255, 255, 255, 00, 255, 00}, /* E */
 {01, 01, 255, 01, 00, 255, 00, 255, 01, 01, 255, 00, 01, 255, 00, 255}};/* F */


BinaryImage::BinaryImage()
{
        bitmap= nullptr;
        bitmap_RLSA= nullptr;
        bitmap_Skeleton= nullptr;
        bitmap_Normalized= nullptr;
        labelmap= nullptr;
        Average_Character_Height=0;
        Image_Height=0;
        Image_Width=0;
        Normalized_Height=0;
        Normalized_Width=0;
        Connected_Components.clear();
}


BinaryImage::~BinaryImage()
{
    if (bitmap!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
        delete bitmap[y];
      delete bitmap;
      bitmap= nullptr;
    }
    if (labelmap!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
        delete labelmap[y];
      delete labelmap;
      labelmap= nullptr;
    }

    if (bitmap_RLSA!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
        delete bitmap_RLSA[y];
      delete bitmap_RLSA;
      bitmap_RLSA= nullptr;
    }

    if (bitmap_Skeleton!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
        delete bitmap_Skeleton[y];
      delete bitmap_Skeleton;
      bitmap_Skeleton= nullptr;
    }

    if (bitmap_Normalized!= nullptr)
    {
      for (int y=0;y<Normalized_Height;y++)
        delete bitmap_Normalized[y];
      delete bitmap_Normalized;
      bitmap_Normalized= nullptr;
    }

    Average_Character_Height=0;
    Image_Height=0;
    Image_Width=0;
    Normalized_Height=0;
    Normalized_Width=0;
    if(Connected_Components.size()>0)
		{
			Connected_Components.clear(); 
			vector <cc>().swap(Connected_Components);
		}
		else
			vector <cc>().swap(Connected_Components);
}



unsigned char BinaryImage::ypo(int x,int y,int Iy,int Dx,unsigned char *IMAGE)
{

  unsigned int xx2=(x-1) >> 3;
  long se=Dx*(Iy-y)+xx2;
  unsigned char xx=*(IMAGE+se);
	if (((xx & (1 << (8-x+8*xx2)))) != 0 ) return(0);
		else return(1);
}


//YES
vector <cc> BinaryImage::connectedcomponents(unsigned char **bitmap,int **labelmap,int Ix,int Iy,bool *MemStatus)
{
  //Efarmozo ton klasiko algorithmo ccs apo to paper
  //To plithos periexei ton arithmo ton ccs pou brethikan sto bitmap
  int plithos;
  vector <cc> pinakas;
  
 ConnectedComponentLabeling_8(bitmap,labelmap,&plithos,Ix,Iy);
  //Bazo 0 ston pinaka labelmap giati o algorithmos kapou exei afisei -1
  for (int y=0;y<Iy;y++)
       for (int x=0;x<Ix;x++)
           {
                if (labelmap[y][x]==-1)
                        labelmap[y][x]=0;
           }
   //Edo na kano orthogonia giro apo ta connected components
   //I domi cc krataei gia kathe cc tis eksis plirofories
   //1)sintetagmenes tou bounding box tou cc
   //2)to id tou connected component ston pinaka labelmap
   //3)tis sintetagmenes tou kentrou mazas tou
   //4)ena flag analoga an exei xrisimopoiithei sti diadikasia i oxi
   //5)ton arithmo ton simeion an afto exei spasei se block(ton arithmo ton blocks)
      struct cc* pinakas1 = new (nothrow) struct cc [plithos+1];
			if (pinakas1== nullptr)
			{
				(*MemStatus)=false;
				return pinakas;
			}
//Ston pinakas apothikevo se kathe grammi [ymin,ymax,xmin,xmax]
        for (int y=1;y<plithos+1;y++)
                {
                        pinakas1[y].ymin=Iy+1;
                        pinakas1[y].ymax=-1;
                        pinakas1[y].xmin=Ix+1;
                        pinakas1[y].xmax=-1;
                }
//Tha saroso tin eikona kai gia kathe simeio tha elesgkso kathe stoixeio tou pinakas
       for (int y=0;y<Iy;y++)
                for (int x=0;x<Ix;x++)
                {
                       if (labelmap[y][x]>0)
                       {
                                if (pinakas1[labelmap[y][x]].ymin>=y)
                                    pinakas1[labelmap[y][x]].ymin=y;
                                if (pinakas1[labelmap[y][x]].ymax<=y)
                                    pinakas1[labelmap[y][x]].ymax=y;
                                if (pinakas1[labelmap[y][x]].xmin>=x)
                                    pinakas1[labelmap[y][x]].xmin=x;
                                if (pinakas1[labelmap[y][x]].xmax<=x)
                                    pinakas1[labelmap[y][x]].xmax=x;
                         }
                   }

        int check=0;
        //Edo petao afta ta ccs pou exoun platos i ipsos mikrotero i iso ton 2 pixels
        for (int y=1;y<plithos+1;y++)
                if ((pinakas1[y].ymax-pinakas1[y].ymin>IPSOS_POU_THELO) || (pinakas1[y].xmax-pinakas1[y].xmin>PLATOS_POU_THELO) )//&& (pinakas1[y].ymax-pinakas1[y].ymin<Iy/3)
                    check+=1;
        //struct cc* pinakas=new (nothrow) struct cc [check+1];
       pinakas.clear();
       pinakas.begin();

       cc help;

        for (int y=1;y<plithos+1;y++)
                {
                //Exei mpei kai periorismos gia poli megala ccs
                  if ((pinakas1[y].ymax-pinakas1[y].ymin>IPSOS_POU_THELO) || (pinakas1[y].xmax-pinakas1[y].xmin>PLATOS_POU_THELO))//&& (pinakas1[y].ymax-pinakas1[y].ymin<Iy/3)
                  {

                        help.ymin=pinakas1[y].ymin;
                        help.ymax=pinakas1[y].ymax;
                        help.xmin=pinakas1[y].xmin;
                        help.xmax=pinakas1[y].xmax;
                        help.label=y;
                        help.points=0;
                        long double m00,m01,m10;
                        m00=m01=m10=0;
                        int arxix,arxiy,telosx,telosy,label;
                        arxix=help.xmin;
                        telosx=help.xmax;
                        arxiy=help.ymin;
                        telosy=help.ymax;
                        label=help.label;
                        for (int a=arxiy;a<=telosy;a++)
                           for (int b=arxix;b<=telosx;b++)
                                 if (labelmap[a][b]==label)
                                 {
                                    m00+=1;
                                    m10+=b;
                                    m01+=a;
                                    help.points++;
                                  }
                        help.yg=m01/m00;
                        help.xg=m10/m00;
                        pinakas.push_back(help);

                }

           }

 //Na apodesmefso ton boithitiko pinaka
	delete [] pinakas1;

	(*MemStatus)=true;
 plithos=check;
 return pinakas;

}


void BinaryImage::Tracer_8(int *cy, int *cx, int *tracingdirection,int **labelmap,unsigned char**bitmap)
{
	int i, y, x;

	for(i = 0; i < 7; i++)
	{
		y = *cy + SearchDirection_8[*tracingdirection][0];
		x = *cx + SearchDirection_8[*tracingdirection][1];

		if(bitmap[y][x] == 0)
		{
			labelmap[y][x] = -1;
			*tracingdirection = (*tracingdirection + 1) % 8;
		}
		else
		{
			*cy = y;
			*cx = x;
			break;
		}
	}
}



void BinaryImage::ContourTracing_8(int cy, int cx, int labelindex, int tracingdirection,int **labelmap,unsigned char **bitmap)
{
	char tracingstopflag = 0, SearchAgain = 1;
	int fx, fy, sx = cx, sy = cy;

	Tracer_8(&cy, &cx, &tracingdirection,labelmap,bitmap);


	if(cx != sx || cy != sy)
	{
		fx = cx;
		fy = cy;

		while(SearchAgain)
		{
			tracingdirection = NextTracingDirection_8[tracingdirection];
			labelmap[cy][cx] = labelindex;
			Tracer_8(&cy, &cx, &tracingdirection,labelmap,bitmap);


			if(cx == sx && cy == sy)
			{
				tracingstopflag = 1;
			}
			else if(tracingstopflag)
			{
				if(cx == fx && cy == fy)
				{
					SearchAgain = 0;
				}
				else
				{
					tracingstopflag = 0;
				}
			}
		}
	}
}


int BinaryImage::ConnectedComponentLabeling_8(unsigned char **bitmap,int **labelmap,int *plithos,int Ix,int Iy)
{
	int height, width, cx, cy, tracingdirection, ConnectedComponentsCount = 0, labelindex;
	height=Iy;
        width=Ix;

	for(cy = 1; cy < height - 1; cy++)
	{
		for(cx = 1, labelindex = 0; cx < width - 1; cx++)
		{
			if(bitmap[cy][cx] == 1)// black pixel
			{
				if(labelindex != 0)// use pre-pixel label
				{

					labelmap[cy][cx] = labelindex;
				}
				else
				{
					labelindex = labelmap[cy][cx];

					if(labelindex == 0)
					{
						labelindex = ++ConnectedComponentsCount;
						tracingdirection = 0;
						ContourTracing_8(cy, cx, labelindex, tracingdirection,labelmap,bitmap);// external contour
						labelmap[cy][cx] = labelindex;
					}
				}
			}
			else if(labelindex != 0)// white pixel & pre-pixel has been labeled
			{
				if(labelmap[cy][cx] == 0)
				{
					tracingdirection = 1;
					ContourTracing_8(cy, cx - 1, labelindex, tracingdirection,labelmap,bitmap);// internal contour
				}

				labelindex = 0;
			}
		}
	}

  *plithos=ConnectedComponentsCount;
  return 1;
}

//YES
int BinaryImage::Connected_Component_Labeling()
{
	bool memstat;
 Connected_Components.clear();
 Connected_Components.begin();
 if (bitmap!= nullptr)
  {
    labelmap=new (nothrow) int* [Image_Height];
    if (labelmap== nullptr)
       return -2;
    for (int y=0;y<Image_Height;y++)
      {
      labelmap[y]=new (nothrow) int [Image_Width];
      if (labelmap[y]== nullptr)
			{
				for (int yx=0;yx<y;yx++)
					delete labelmap[yx];
				delete labelmap;
				labelmap= nullptr;
        return -3;
			}
      for (int x=0;x<Image_Width;x++)
         labelmap[y][x]=0;
      }
      Connected_Components=connectedcomponents(bitmap,labelmap,Image_Width,Image_Height,&memstat);
      if(memstat!=true) return -4;
			else return 0;        
  }
 else
   return -1;

}

//YES
vector <cc> BinaryImage::Get_Connected_Components()
{
 /*int a=Connected_Components.size();
 if (a==0)
 {
    Connected_Component_Labeling();
 }*/
 return Connected_Components;
}



void BinaryImage::Calculate_Average_Character_Height()
{
        int max=0;
        int plithos;
        plithos=Connected_Components.size();
        for (int y=0;y<plithos;y++)
            if (max<Connected_Components[y].ymax-Connected_Components[y].ymin+1)
                max=Connected_Components[y].ymax-Connected_Components[y].ymin+1;
        int *hist;
        //Desmevo mnimi gia to istogramma
        hist=new (nothrow) int [max+1];
        if (hist== nullptr)
           printf("Problem with hist");
        for (int y=0;y<max+1;y++)
            hist[y]=0;
        for (int y=0;y<plithos;y++)
           {
                int ipsos=Connected_Components[y].ymax-Connected_Components[y].ymin+1;
                //int platos=Connected_Components[y].xmax-Connected_Components[y].xmin;
                hist[ipsos]+=2;
                if (ipsos!=0)
                 hist[ipsos-1]+=1;
                if (ipsos!=max)
                    hist[ipsos+1]+=1;

	    }

        int korifi=0;
        int poukorifi;
        for (int y=6;y<max+1;y++)
            if (hist[y]>=korifi)
            {
              korifi=hist[y];
              poukorifi=y;
              }

        //Apodesmefsi tis mnimis hist
        delete [] hist;
        hist= nullptr;
     //    return poukorifi;
     int *dianisma;
     dianisma=new (nothrow) int [plithos+1];
     if (dianisma== nullptr)
        printf("Problem with dianisma");
     for (int y=0;y<plithos;y++)
        dianisma[y]=Connected_Components[y].ymax-Connected_Components[y].ymin+1;
     //Taksinomo kai pairno to mesaio
     for (int i=2;i<plithos+1;i++)
        for (int j=plithos;j>=i;j--)
           if (dianisma[j]<dianisma[j-1])
           {
             int help=dianisma[j];
             dianisma[j]=dianisma[j-1];
             dianisma[j-1]=help;
           }
    int mmmmm=dianisma[plithos/2];
    delete dianisma;
    dianisma= nullptr;
    int gianado=(poukorifi+mmmmm)/2;
    Average_Character_Height=gianado;
}


int BinaryImage::Get_Average_Character_Height()
{
  //returns the average character height or -1 if it can not be calculated      
  if (Connected_Components.empty()!=true)
  {
        Calculate_Average_Character_Height();
        return Average_Character_Height;
  }
  else
        return -1;
}


int BinaryImage::Get_Image_Width()
{
        return this->Image_Width;
}


int BinaryImage::Get_Image_Height()
{
        return this->Image_Height;
}

unsigned char BinaryImage::Get_Bitmap(int x,int y)
{
  //returns 0 or 1 if the coordinates are fine else returns 255
  if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height) && (bitmap!= nullptr))
        return bitmap[y][x];
  else
        return 255;
}



unsigned char BinaryImage::Get_Bitmap_RLSA(int x,int y)
{
  //returns 0 or 1 if the coordinates are fine else returns 255
  if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height) && (bitmap_RLSA!= nullptr))
        return bitmap_RLSA[y][x];
  else
        return 255;
}


unsigned char BinaryImage::Get_Bitmap_Skeleton(int x,int y)
{
  //returns 0 or 1 if the coordinates are fine else returns 255
  if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height) && (bitmap_Skeleton!= nullptr))
        return bitmap_Skeleton[y][x];
  else
        return 255;
}


unsigned char BinaryImage::Get_Bitmap_Normalized(int x,int y)
{
  //returns 0 or 1 if the coordinates are fine else returns 255
  if ( (x>=0) && (x<Normalized_Width) && (y>=0) && (y<Normalized_Height) && (bitmap_Normalized!= nullptr))
        return bitmap_Normalized[y][x];
  else
        return 255;
}

int BinaryImage::Get_Labelmap(int x,int y)
{
  //returns the label if the coordinates are fine else returns -1
  if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height) && (labelmap!= nullptr))
        return labelmap[y][x];
  else
        return -1;
}


int BinaryImage::Set_Bitmap(int x,int y,unsigned char value)
{
   if (bitmap!= nullptr)
   {
         if ( (x>=0) && (x<Image_Width) && (y>=0) && (y<Image_Height))
         {
            bitmap[y][x]=value;
            return 0;
         }
         else
            return -1;
   }
   else
         return -1;
}


//YES
int BinaryImage::Create_Bitmap(int Ix,int Iy)
{
  //Destroying old data after multiple calls of Load
  this->Clear();
  //returns 0 if success
  //-1 if Ix and or Iy <0
  //-2 if new returns error
  if ( (Ix>0) && (Iy>0))
  {
     bitmap=new (nothrow) unsigned char* [Iy];
     if (bitmap== nullptr)
        return -2;
     for (int y=0;y<Iy;y++)
        {
           bitmap[y]=new (nothrow) unsigned char [Ix];
           if (bitmap[y]== nullptr)
					 {
						 for (int yx=0;yx<y;yx++)
								delete bitmap[yx];
						delete bitmap;
						bitmap= nullptr;
              return -2;
					 }
        }
     for (int y=0;y<Iy;y++)
        for (int x=0;x<Ix;x++)
           bitmap[y][x]=0;

     Image_Height=Iy;
     Image_Width=Ix;
     return 0;

  }
  else
     return -1;
}


//YES
void BinaryImage::Clear()
{
 //Destroying old data after multiple calls of Load
 if (bitmap!= nullptr)
 {
    for (int y=0;y<Image_Height;y++)
       delete bitmap[y];
    delete bitmap;
    bitmap= nullptr;
 }
 if (labelmap!= nullptr)
 {
    for (int y=0;y<Image_Height;y++)
       delete labelmap[y];
    delete labelmap;
    labelmap= nullptr;
 }
 if (bitmap_RLSA!= nullptr)
 {
    for (int y=0;y<Image_Height;y++)
       delete bitmap_RLSA[y];
    delete bitmap_RLSA;
    bitmap_RLSA= nullptr;
 }

 if (bitmap_Skeleton!= nullptr)
 {
    for (int y=0;y<Image_Height;y++)
       delete bitmap_Skeleton[y];
    delete bitmap_Skeleton;
    bitmap_Skeleton= nullptr;
 }

 if (bitmap_Normalized!= nullptr)
 {
    for (int y=0;y<Normalized_Height;y++)
       delete bitmap_Normalized[y];
    delete bitmap_Normalized;
    bitmap_Normalized= nullptr;
 }

 if (Connected_Components.size()>0)
 {
    Connected_Components.clear();
		vector <cc>().swap(Connected_Components);
 }
 else
	 vector <cc>().swap(Connected_Components);

 Image_Height=0;
 Image_Width=0;
 Normalized_Height=0;
 Normalized_Width=0;
 Average_Character_Height=0;

}



int BinaryImage::Calculate_RLSA(int hthr,int vthr)
{
     if (bitmap== nullptr)
        return -1;

     if (bitmap_RLSA!= nullptr)
       {
          for (int y=0;y<Image_Height;y++)
             delete bitmap_RLSA[y];
          delete bitmap_RLSA;
          bitmap_RLSA= nullptr;   
       }
       
     bitmap_RLSA=new (nothrow) unsigned char* [Image_Height];
     if (bitmap_RLSA== nullptr)
        return -2;
     for (int y=0;y<Image_Height;y++)
        {
           bitmap_RLSA[y]=new (nothrow) unsigned char [Image_Width];
           if (bitmap_RLSA[y]== nullptr)
              return -3;
           for (int x=0;x<Image_Width;x++)
              bitmap_RLSA[y][x]=bitmap[y][x];   
        }
    
	//CCA
     int ThHor;
     ThHor=hthr;

	 //Hor RLSA

     for (int y=0;y<Image_Height;y++)
     {
        int x=0;
        while (x<Image_Width-1)
        {
          while ((bitmap[y][x]==1) && (x<Image_Width-1)) x++;
          int X=x;
          while ((bitmap[y][x]!=1) && (x<Image_Width-1)) x++;
          if ((x-1-X<ThHor) && (X!=0) && (x<Image_Width-1))
          {
            for (int ix=X;ix<x;ix++)
            {
              bitmap_RLSA[y][ix]=1;
            }
          }
        }
     }

    int ThVer;
    ThVer=vthr;


  //Ver RLSA
   for (int x=0;x<Image_Width;x++)
   {
      int y=0;
      while (y<Image_Height-1)
      {
        while ((bitmap[y][x]==1) && (y<Image_Height-1))
                y++;
        int Y=y;
        while ((bitmap[y][x]!=1) && (y<Image_Height-1))
                y++;
        if ((y-1-Y<ThVer) && (Y!=0) && (y<Image_Height-1))
        {
          for (int iy=Y;iy<y;iy++)
          {
            bitmap_RLSA[iy][x]=1;
          }
        }
      }
   }

  return 0;


}

//YES
int BinaryImage::Skeleton()
{
  //IMA kai IMA1 unsigned char
  //IMA2x int
  unsigned char *IMA,*IMA1;
  unsigned int *IMA2x;
	IMA= nullptr; IMA1= nullptr; IMA2x= nullptr;

  long Yq;


  bool Change = true;

  // Elegxos gia apodesmefsi mnimis tou bitmap_Skeleton
  if (bitmap_Skeleton!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
         delete bitmap_Skeleton[y];
      delete bitmap_Skeleton;
      bitmap_Skeleton= nullptr;
    }

  // Desmefsi mnimis gia ton pinaka bitmap_Skeleton
  bitmap_Skeleton=new (nothrow) unsigned char *[Image_Height];
	if(bitmap_Skeleton== nullptr) return -1;

  for (int y=0;y<Image_Height;y++)
     {
        bitmap_Skeleton[y]=new (nothrow) unsigned char [Image_Width];
        if(bitmap_Skeleton[y]== nullptr)
				{
					for (int yx=0;yx<y;yx++)
						delete bitmap_Skeleton[yx];
						
					delete bitmap_Skeleton;
					bitmap_Skeleton= nullptr;
					return -1;
				}

				for (int x=0;x<Image_Width;x++)
           bitmap_Skeleton[y][x]=0;
     }

  IMA = (unsigned char *)malloc ((Image_Width)*(Image_Height)*sizeof(unsigned char));
	if(IMA== nullptr)
	{
		if (bitmap_Skeleton!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
         delete bitmap_Skeleton[y];
      delete bitmap_Skeleton;
      bitmap_Skeleton= nullptr;
    }
		return -1;
	}

  IMA1 = (unsigned char *)calloc ((Image_Width)*(Image_Height),sizeof(unsigned char));
	if(IMA1== nullptr)
	{
		if (bitmap_Skeleton!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
         delete bitmap_Skeleton[y];
      delete bitmap_Skeleton;
      bitmap_Skeleton= nullptr;
    }
		if(IMA!= nullptr) free(IMA);
		return -1;
	}

  IMA2x = (unsigned int *)malloc ((Image_Width)*(Image_Height)*sizeof(unsigned int));
	if(IMA2x== nullptr)
	{
		if (bitmap_Skeleton!= nullptr)
    {
      for (int y=0;y<Image_Height;y++)
         delete bitmap_Skeleton[y];
      delete bitmap_Skeleton;
      bitmap_Skeleton= nullptr;
    }
		if(IMA!= nullptr) free(IMA);
		if(IMA1!= nullptr) free(IMA1);
		return -1;
	}
  
  
  //Orizo 4 neous pinakes apo pano
  //O IMA kai IMA2x antigrafo ton bitmap
  for (int y=0;y<Image_Height;y++)
	for (int x=0;x<Image_Width;x++)
		{
		*(IMA+y*Image_Width+x) = Get_Bitmap(x,y);
		*(IMA2x+y*Image_Width+x) = Get_Bitmap(x,y);
		}

  bool Ep=false;
  while (Change)
  {
    Change = false;
    if (Ep==true)
        Ep=false;
    else
        Ep=true;
    for (int y=1;y<Image_Height-1;y++)
    {
      Yq = y*Image_Width;
      for (int x=1;x<Image_Width-1;x++)
      {
        if (*(IMA+Yq+x)==1)
        //An to pixel einai mavro
        {
            if (Ep)
            {
              if ((ZNZT(x,y,IMA)==1) &&  (NZN(x,y,IMA)>=2)  &&  (NZN(x,y,IMA)<=6)  &&  ((*(IMA+(y-1)*Image_Width+x)!=1) || (*(IMA+Yq+x+1)!=1) || (*(IMA+(y+1)*Image_Width+x)!=1)) &&   ((*(IMA+Yq+x+1)!=1) || (*(IMA+(y+1)*Image_Width+x)!=1) || (*(IMA+Yq+x-1)!=1)))
                {
                   *(IMA1+Yq+x)=1;
                   Change = true;
                }
            }
            else
            {
                if ((ZNZT(x,y,IMA)==1) && (NZN(x,y,IMA)>=2)  &&    (NZN(x,y,IMA)<=6)  &&  ((*(IMA+(y-1)*Image_Width+x)!=1) || (*(IMA+Yq+x+1)!=1) || (*(IMA+Yq+x-1)!=1)) && ((*(IMA+(y-1)*Image_Width+x)!=1) || (*(IMA+(y+1)*Image_Width+x)!=1) || (*(IMA+Yq+x-1)!=1)))
                {
                   *(IMA1+Yq+x)=1;
                   Change = true;
                }
            }
        }
      }
    }

    if (Change)//An exei anapsei to Change tote sbino to pixel ston pinaka IMA 
    {
      for (int y=1;y<Image_Height-1;y++)
      {
        Yq = y*Image_Width;
        for (int x=1;x<Image_Width-1;x++)
        {
          if (*(IMA1+Yq+x)==1)
          {
             *(IMA+Yq+x)=0;
          }
        }
      }
    }

  }   //Edo telionei to while

  //Ench //Afto einai mallon kati parapano sti methodo
  Change=true;
  while(Change) 
	Change = Ench(IMA);


  for (int y=0;y<Image_Height;y++)
  for (int x=0;x<Image_Width;x++)
  {
    if ((*(IMA+y*Image_Width+x)==1) || (*(IMA+y*Image_Width+x)==2) || (*(IMA+y*Image_Width+x)==3))
      bitmap_Skeleton[y][x]=1; //Edo na mpei afto pou prepei!!
    else
      bitmap_Skeleton[y][x]=0; //Kai edo!!!
  }

  if(IMA!= nullptr) free(IMA);
  if(IMA1!= nullptr) free(IMA1);
  if(IMA2x!= nullptr) free(IMA2x);
  IMA= nullptr;
  IMA1= nullptr;
  IMA2x= nullptr; 
  return 0; 
}




int BinaryImage::NZN(int x,int y,unsigned char *IMA)
{
   int mNZN = 0;

   if (*(IMA+(y-1)*Image_Width+x-1)==1) mNZN++;
   if (*(IMA+(y-1)*Image_Width+x)==1)  mNZN++;
   if (*(IMA+(y-1)*Image_Width+x+1)==1)  mNZN++;
   if (*(IMA+y*Image_Width+x+1)==1)  mNZN++;
   if (*(IMA+(y+1)*Image_Width+x+1)==1)  mNZN++;
   if (*(IMA+(y+1)*Image_Width+x)==1)  mNZN++;
   if (*(IMA+(y+1)*Image_Width+x-1)==1)  mNZN++;
   if (*(IMA+y*Image_Width+x-1)==1)  mNZN++;

   return mNZN;
}

int BinaryImage::ZNZT(int x,int y,unsigned char *IMA)
{
   int mZNZT = 0;
   if  (*(IMA+y*Image_Width+x)==0) return 0;
   if  (((*(IMA+(y-1)*Image_Width+x-1)==0)) && ((*(IMA+(y-1)*Image_Width+x)==1))) mZNZT++;
   if  (((*(IMA+(y-1)*Image_Width+x)==0)) && ((*(IMA+(y-1)*Image_Width+x+1)==1))) mZNZT++;
   if  (((*(IMA+(y-1)*Image_Width+x+1)==0)) && ((*(IMA+y*Image_Width+x+1)==1))) mZNZT++;
   if  (((*(IMA+y*Image_Width+x+1)==0)) && ((*(IMA+(y+1)*Image_Width+x+1)==1))) mZNZT++;
   if  (((*(IMA+(y+1)*Image_Width+x+1)==0)) && ((*(IMA+(y+1)*Image_Width+x)==1))) mZNZT++;
   if  (((*(IMA+(y+1)*Image_Width+x)==0)) && ((*(IMA+(y+1)*Image_Width+x-1)==1))) mZNZT++;
   if  (((*(IMA+(y+1)*Image_Width+x-1)==0)) && ((*(IMA+y*Image_Width+x-1)==1))) mZNZT++;
   if  (((*(IMA+y*Image_Width+x-1)==0)) && ((*(IMA+(y-1)*Image_Width+x-1)==1))) mZNZT++;
   return mZNZT;
}

bool BinaryImage::Ench(unsigned char *IMA)
{
  long Yq;
  bool changed=false;
  int P2, P3, P4, P5;
  int P6, P7, P8, P9,i,j;//,cnt=0
  for (int y=1; y<Image_Height-1; y++)
  {
    Yq = y*Image_Width;
    for (int x=1; x<Image_Width-1; x++)
    {
      if (*(IMA+Yq+x) == 1)
      {
	P2 = (*(IMA+(y-1)*Image_Width+x)==1) << 0;
	P3 = (*(IMA+(y-1)*Image_Width+x+1)==1) << 1;
	P4 = (*(IMA+Yq+x+1)==1) << 2;
	P5 = (*(IMA+(y+1)*Image_Width+x+1)==1) << 3;
	P6 = (*(IMA+(y+1)*Image_Width+x)==1) << 0;
	P7 = (*(IMA+(y+1)*Image_Width+x-1)==1) << 1;
	P8 = (*(IMA+Yq+x-1)==1) << 2;
	P9 = (*(IMA+(y-1)*Image_Width+x-1)==1) << 3;

	i = P9+P8+P7+P6;
	j = P5+P4+P3+P2;

      	if (lc_look_up[i][j] == 00)
        {
          *(IMA+Yq+x)=0;
          changed=true;
        }
      }

    }

  }
  return changed;
}





int BinaryImage::Normalize(int norm_width, int norm_height)
{
        float logos_x, logos_y;
        float decimal, decimal2;
        double ip;
        float heightOrwidth;
        int diafora;
        int log_x,log_y;

        Normalized_Width=norm_width;
        Normalized_Height=norm_height;
        bitmap_Normalized=new (nothrow) unsigned char *[norm_height];
        if (bitmap_Normalized== nullptr)
           return -1;
        for (int y=0;y<norm_height;y++)
           {
             bitmap_Normalized[y]=new (nothrow) unsigned char [norm_width];
             if (bitmap_Normalized[y]== nullptr)
                return -2;
             for (int x=0;x<norm_width;x++)
                bitmap_Normalized[y][x]=0;
           }

        unsigned char **temp=(unsigned char **)malloc(norm_width*sizeof(int));
        for (int i=0; i<norm_width; i++)
                temp[i]=(unsigned char *)malloc(norm_height*sizeof(int));

        for (int i=0; i<norm_width; i++)
        for (int j=0; j<norm_height; j++)
                temp[i][j]=0;

        float a=(float)norm_width/(float)Image_Width;
        float b=(float)norm_height/(float)Image_Height;


        for (int x=0; x<norm_width; x++)
        for (int y=0; y<norm_height; y++)
        {
                if((a/b>=1))
                {
                        logos_x=((float)x)/b;
                        logos_y=((float)y)/b;
                }
                else
                {
                        logos_x=((float)x)/a;
                        logos_y=((float)y)/a;
                }

                decimal=modf(logos_x, &ip);
                decimal2=modf(logos_y, &ip);
                if (decimal<0.5)
                        logos_x=floor(logos_x);
                else
                        logos_x=ceil(logos_x);
                if (decimal2<0.5)
                        logos_y=floor(logos_y);
                else
                        logos_y=ceil(logos_y);

                log_x=logos_x;
                log_y=logos_y;

                if  (log_x>=Image_Width || log_y>=Image_Height )
                        break;
                if (bitmap[log_y][log_x]==1)
                        temp[x][y]=1;
        }

     
        if((a/b>=1))
                heightOrwidth=b*Image_Width;
        else
                heightOrwidth=a*Image_Height;

        heightOrwidth=heightOrwidth/2;
        decimal=modf(heightOrwidth, &ip);
        if (decimal<0.5)
                heightOrwidth=floor(heightOrwidth);
        else
                heightOrwidth=ceil(heightOrwidth);

        if((a/b>=1))
                diafora=(norm_width/2)-heightOrwidth;
        else
                diafora=(norm_height/2)-heightOrwidth;

        diafora=abs(diafora);

        for (int x=0; x<norm_width; x++)
        {
                for (int y=0; y<norm_height; y++)
                {
                        if (temp[x][y]==1)
                        {
                                if((a/b>=1))
                                        bitmap_Normalized[y][x+diafora]=1;
                                else
                                        bitmap_Normalized[y+diafora][x]=1;
                        }
                }
        }


        for (int x=0; x<norm_width; x++)
            free(temp[x]);
        free(temp);
        temp= nullptr;    
        return 0;
}




int BinaryImage::NormalizeBGAT(int norm_width, int norm_height)
{
     int Hist[1000];
     int Iy=Image_Height;
     int Ix=Image_Width;

     for (int y=0;y<Iy;y++)
     {
        Hist[y]=0;
        for (int x=0;x<Ix;x++)
        if (bitmap[y][x]==1)
                Hist[y]++;
     }

     int y1=Iy/2; int Hm=Hist[y1];
     while ((Hist[y1]>Hm/2) && (y1>0)) y1--;

     int y2=Iy/2;
     while ((Hist[y2]>Hm/2) && (y2<Iy-1)) y2++;


     int x1=0;
     bool IsOK=false;
     while ((!IsOK) && (x1<Ix-1))
     {
       for (int y=0;y<Iy;y++)
         if (bitmap[y][x1]==1) {IsOK=true;y=Iy;}
       if (!IsOK) x1++;
     }

     int x2=Ix-1;
     IsOK=false;
     while ((!IsOK) && (x2>0))
     {
       for (int y=0;y<Iy;y++)
         if (bitmap[y][x2]==1) {IsOK=true;y=Iy;}
       if (!IsOK) x2--;
     }

     float a1=1;if (y2!=y1) a1 = (float)30/(y2-y1);
     float a2=1;if (x2!=x1) a2 = (float)300/(x2-x1);
     float b1=0;if (y2!=y1) b1 = (float)(30*y2-60*y1)/(y2-y1);
     float b2=0;if (x1!=x2) b2 = (float)(300*x1)/(x1-x2);
     
    Normalized_Width=norm_width;
    Normalized_Height=norm_height;
    bitmap_Normalized=new (nothrow) unsigned char *[norm_height];
    if (bitmap_Normalized== nullptr)
       return -1;
    for (int y=0;y<norm_height;y++)
       {
         bitmap_Normalized[y]=new (nothrow) unsigned char [norm_width];
         if (bitmap_Normalized[y]== nullptr)
            return -2;
         for (int x=0;x<norm_width;x++)
            bitmap_Normalized[y][x]=0;
       }

   for (int x=0;x<300;x++)
     for (int y=0;y<90;y++)
     {
      int xx = (x-b2)/a2;
      int yy = (y-b1)/a1;
      if ((xx>=0) && (xx<Ix) &&  (yy>=0) && (yy<Iy))
      if  (bitmap[yy][xx]==1)
        bitmap_Normalized[y][x]=1;

     }
     return 0;
}
