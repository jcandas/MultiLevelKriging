/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef BYTESWAPPING_H
#define BYTESWAPPING_H

#include <assert.h>

static inline bool isBigEndian();
static inline bool isLittleEndian();

static inline void swapByteOrder(float& input);
static inline void swapByteOrder(double& input);
static inline void swapByteOrder(unsigned int& input);
static inline void swapByteOrder(int& input);
static inline void swapByteOrder(unsigned short& input);

void swapByteOrder(float* input, unsigned int num);
void swapByteOrder(double* input, unsigned int num);
void swapByteOrder(unsigned int* input, unsigned int num);
void swapByteOrder(int* input, unsigned int num);
void swapByteOrder(unsigned short* input, unsigned int num);

static inline float flip(float input);
static inline double flip(double input);
static inline unsigned int flip(unsigned int input);
static inline int flip(int input);
static inline unsigned short flip(unsigned short input);


static inline bool isBigEndian()
{
	assert(sizeof(unsigned int)==4);
	unsigned int intVersion= 0x89ABCDEF;
	unsigned char* charVersion = (unsigned char*)&intVersion;
	return charVersion[0]==0x89 && 
		charVersion[1]==0xAB && 
		charVersion[2]==0xCD &&
		charVersion[3]==0xEF;
}

static inline bool isLittleEndian()
{
	assert(sizeof(unsigned int)==4);
	unsigned int intVersion= 0x89ABCDEF;
	unsigned char* charVersion = (unsigned char*)&intVersion;
	return charVersion[3]==0x89 && 
		charVersion[2]==0xAB && 
		charVersion[1]==0xCD &&
		charVersion[0]==0xEF;
}

static inline void swapByteOrder(float& input)
{
	assert(sizeof(float)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
}

static inline void swapByteOrder(double& input)
{
	assert(sizeof(double)==8);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<4; c++) {
		temp = bytes[c];
		bytes[c] = bytes[7-c];
		bytes[7-c] = temp;
	}
}

static inline void swapByteOrder(unsigned int& input)
{
	assert(sizeof(unsigned int)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
}

static inline void swapByteOrder( int& input)
{
	assert(sizeof(int)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
}

static inline void swapByteOrder(unsigned short& input)
{
	assert(sizeof(unsigned short)==2);
	unsigned char * bytes = (unsigned char*)&input;
	unsigned char temp;
	temp = bytes[0];
	bytes[0] = bytes[1];
	bytes[1] = temp;
}

static inline float flip(float input)
{
	assert(sizeof(float)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
	return input;
}

static inline double flip(double input)
{
	assert(sizeof(double)==8);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<4; c++) {
		temp = bytes[c];
		bytes[c] = bytes[7-c];
		bytes[7-c] = temp;
	}
	return input;
}

static inline unsigned int flip(unsigned int input)
{
	assert(sizeof(unsigned int)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
	return input;
}

static inline int flip(int input)
{
	assert(sizeof(int)==4);
	unsigned char * bytes = (unsigned char*)&input;
	int c;
	unsigned char temp;
	for (c=0; c<2; c++) {
		temp = bytes[c];
		bytes[c] = bytes[3-c];
		bytes[3-c] = temp;
	}
	return input;
}

static inline unsigned short flip(unsigned short input)
{
	assert(sizeof(unsigned short)==2);
	unsigned char * bytes = (unsigned char*)&input;
	unsigned char temp;
	temp = bytes[0];
	bytes[0] = bytes[1];
	bytes[1] = temp;
	return input;
}


#endif


