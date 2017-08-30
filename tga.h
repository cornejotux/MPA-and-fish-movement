/*
	A quick and dirty class for reading and writing tga files.
	Can read and write tga files which are 24 or 32 bit color.
	Can set individual pixels of an image.
	Endian independent.
	Cannot read or write tga files with RLE image data.
	Cannot read or write tga files with bytes-per-pixel other than 24 or 32.
*/


#ifndef _TGA_
#define _TGA_

class TGAImage
{
public:
	TGAImage();
	~TGAImage();
	
	bool read(const char* filename);  //reads in a file replacing current data
	bool write(const char* filename); //writes current data as a tga file
	
	bool resize(int inBitsPerPixel, int inWidth, int inHeight); //For making a tga from scratch.  Calls realloc
	
	void setPixel(int x, int y, double inRed, double inGreen, double inBlue, double inAlpha=1.0);
	void setPixelExactly(int x, int y, unsigned char inRed, unsigned char inGreen, unsigned char inBlue, unsigned char inAlpha = 255);
	
	// Accessors:
	int bitsPerPixel() {return mBitsPerPixel;}
	int width() {return mWidth;}
	int height() {return mHeight;}
	unsigned char* data() {return mData;}
	
private:
	int mBitsPerPixel;
	int mWidth;    
	int mHeight;   // width and height in pixels
	unsigned char* mData;   // Raw image data
};

#endif

