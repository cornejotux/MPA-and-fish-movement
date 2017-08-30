
#include "tga.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


TGAImage::TGAImage()
{
	mWidth = mHeight = mBitsPerPixel = 0;
	mData = NULL;
}

TGAImage::~TGAImage()
{
	if(mData!=NULL)
		free(mData);
}



bool TGAImage::resize(int inBitsPerPixel, int inWidth, int inHeight) //For making a tga from scratch
{
	mWidth = inWidth;
	mHeight = inHeight;
	mBitsPerPixel = inBitsPerPixel;
	
	int bytesPerPixel = mBitsPerPixel/8;
	
	int dataSize = mWidth*mHeight*bytesPerPixel;
	
	//reserve memory to hold the raw image mData
	if(mData==NULL)
		mData = (unsigned char *)malloc(dataSize);
	else
		mData = (unsigned char *)realloc(mData, dataSize);
	
	if( !mData )
	{
		printf( "memory failed to allocate\n" );
		return false;
	}
	
	return true;
}




bool TGAImage::read(const char* filename)
{
    // Uncompressed TGA header
    unsigned char TGAheader[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    unsigned char TGAcompare[12];   // Used to compare TGA header
    unsigned char header[6];        // First 6 useful bytes from the header
	
    int identSize;		// Size of an optional user-defined string header to occur after the real header
    int bytesPerPixel;	        // Recall the member variable is BITS per pixel
    int dataSize;		// Used to store the image size when setting aside RAM
    int temp, i;
	
	
	
    FILE *file = fopen( filename, "rb" );   // Open the TGA file
	
	// Check if the file successfully opened, we can read its header, its
	// header matches what we want, and we can read the following 6 bytes
	if( file==NULL )
	{
		printf( "tga file %s didn't even open\n", filename );
		return false;
	}
	
	
	if( fread( TGAcompare, 1, sizeof( TGAcompare ), file ) != sizeof( TGAcompare ) )
	{
		fclose( file );
		printf( "tga file %s TGAheader didn't load.\n", filename );
		return false;
	}
	
	identSize = TGAcompare[0];
	
	if(	memcmp( TGAheader+1, TGAcompare+1, sizeof( TGAheader-1 ) ) != 0 )
	{
		fclose( file );
		printf( "tga file %s header is incorrect.  Is:\n", filename );
		for(int i=0; i<12; i++)
			printf( "%2.x ", TGAheader[i] );
		printf( "\n" );
		for(int i=0; i<12; i++)
			printf( "%2.x ", TGAcompare[i] );
		printf( "\n" );
		return false;
	}
		
	if( fread(header, 1, sizeof( header ), file) != sizeof(header) )
    {
		fclose( file );
		printf( "tga file %s header (first 6 _useful_ bytes) didn't load.\n", filename );
		return false;
    }
	
	if( fseek(file, identSize, SEEK_CUR)!=0 )
	{
		fclose(file);
		printf( "tga file %s failed to contain its optional user-defined header\n", filename );
		return false;
	}

	
	mWidth = header[1]*256 + header[0];
	mHeight =  header[3]*256 + header[2];
	
	
	// Make sure width and height are nonnegative, and TGA is 24 or 32 bit
	if(mWidth <= 0 || mHeight <= 0) 
	{
		fclose( file );
		printf( "tga file %s has bad width and height: %d %d\n", filename, mWidth, mHeight );
		return false;
	}
	
	mBitsPerPixel = header[4];
	
	if( header[4] != 24 && header[4] != 32 )
	{
		fclose( file );
		printf( "tga file %s has bad number of bits-per-pixel: %d\n", filename, mBitsPerPixel );
		return false;
	}

	
	bytesPerPixel = mBitsPerPixel/8;
	dataSize = mWidth * mHeight * bytesPerPixel;
	
	//reserve memory to hold the raw image mData
	if(mData==NULL)
		mData = (unsigned char *)malloc(dataSize); 
	else
		mData = (unsigned char *)realloc(mData, dataSize);
	
	
	if( mData == NULL )
	{
		fclose( file );
		printf( "tga file %s required %d bytes for image mData.  Unable to allocate.\n", filename, dataSize );
		return false;
	}
	
	if( fread(mData, 1, dataSize, file) != dataSize )
	{
		fclose( file );
		printf( "tga file %s mData size is supposed to be %d.  Unable to load that many bytes from file.\n", filename, dataSize );
		return false;
	}
	
	#if 0
	for( i=0; i < dataSize; i += bytesPerPixel )
	{
		temp = mData[i];  // Temporarily store the value at i
		// Set the first byte to the value of the third byte
		mData[i] = mData[i+2];
		// Set the third byte to the value in temp (first byte value)
		mData[i+2] = temp;
	}
	#endif
	
	fclose( file );

	return true;
}


bool TGAImage::write(const char* filename)
{
	unsigned char TGAheader[12] = {0,0,2,0,0,0,0,0,0,0,0,0};
	unsigned char header[6] = {0,0,0,0,0,0};
	
	FILE* file = fopen(filename, "wb");
	if( file==NULL )
	{
		printf( "file %s failed to open for writing.\n", filename );
		return false;
	}
	
	int bytesPerPixel = mBitsPerPixel/8;
	int dataSize = mWidth*mHeight*bytesPerPixel;
	
	
	header[0] = mWidth&0xff;
	header[1] = mWidth>>8;
	
	header[2] = mHeight&0xff;
	header[3] = mHeight>>8;
	
	header[4] = mBitsPerPixel;
	
	mWidth = header[1] * 256 + header[0];
	mHeight = header[3] * 256 + header[2];
	mBitsPerPixel = header[4];
	
	fwrite(TGAheader, 1, sizeof(TGAheader), file);
	fwrite(header, 1, sizeof(header), file);
	fwrite(mData, 1, dataSize, file);
	fclose(file);
	
	return true;
}





void TGAImage::setPixel(int x, int y, double inRed, double inGreen, double inBlue, double inAlpha)
{
	int r=(int)(inRed*255.0), g=(int)(inGreen*255.0), b=(int)(inBlue*255.0), a=(int)(inAlpha*255.0);
	
	r = (r>255)?255:r;
	g = (g>255)?255:g;
	b = (b>255)?255:b;
	a = (a>255)?255:a;
	
	r = (r<0)?0:r;
	g = (g<0)?0:g;
	b = (b<0)?0:b;
	a = (a<0)?0:a;
    
	setPixelExactly(x,y, r,g,b,a);
}


void TGAImage::setPixelExactly(int x, int y, unsigned char inRed, unsigned char inGreen, unsigned char inBlue, unsigned char inAlpha)
{
	int bytesPerPixel = mBitsPerPixel/8;
	int index = bytesPerPixel * (y*mWidth+x);
	
	if( bytesPerPixel>=4 )
		mData[index+3] = inAlpha;
	
	mData[index+2] = inRed;
	mData[index+1] = inGreen;
	mData[index+0] = inBlue;
}




