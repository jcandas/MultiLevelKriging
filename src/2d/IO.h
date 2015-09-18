#ifndef CCV_IO_H
#define CCV_IO_H

#include <stdio.h>
#include <sys/time.h>
#include "diag_blk.h"
#include "hb.h"

void readRawiv( const char* fileName, fftw_complex* data, int size );
void read_diag_blk(int argc, char* argv[], Diag_Blk* db);
void write_diag_blk(int argc, char* argv[], HB* hb, Diag_Blk* db);

// extern "C" double second();   /* Use for IBM processors. */

#endif
