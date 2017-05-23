/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Markus P. Fuelscher                              *
***********************************************************************/
/************************************************************************/
/*                                                                      */
/*     purpose:                                                         */
/*     -  pack/unpack a buffer of double precision floating             */
/*        point numbers with a given absolute accuracy                  */
/*     -  pack/unpack a buffer of of integer numbers                    */
/*                                                                      */
/*     entry points:                                                    */
/*     rzip    : pack double precision numbers                          */
/*     runzip  : unpack double precision numbers                        */
/*     rziplen : return length of packed numbers in units of bytes      */
/*     izip    : pack integer numbers                                   */
/*     iunzip  : unpack integer numbers                                 */
/*     iziplen : return length of packed numbers in units of bytes      */
/*                                                                      */
/*     calling parameters:                                              */
/*     nData : length of input vector in units of double precison       */
/*             words (rzip/runzip) or integer words (izip/iunzip)       */
/*     nByte : length of pack array in units of bytes                   */
/*     thrs  : absolute accuracy to be retained in packing/upacking     */
/*             double precision numbers                                 */
/*     InBuf : vector of input data                                     */
/*     OutBuf: vector of packed data ([ir]zip/[irunzip] or              */
/*             length of packed word in units of bytes ([ir]ziplen)     */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*     written by:                                                      */
/*     M.P. Fuelscher                                                   */
/*     University of Lund, Sweden, 1999                                 */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*     history:                                                         */
/*     - replacement of packing routines used in MOLCAS 4.1             */
/*       M.P. Fuelscher                                                 */
/*       University of Lund, Sweden, 1999                               */
/*                                                                      */
/************************************************************************/
#ifndef _WIN32_
#include <unistd.h>
#endif
#include <stdio.h>
#include <string.h>
#include "molcastype.h"

#ifdef _CAPITALS_
#define rzip RZIP
#define runzip RUNZIP
#define rziplen RZIPLEN
#define izip IZIP
#define iunzip IUNZIP
#define iziplen IZIPLEN
#else
#ifndef ADD_
#define rzip rzip_
#define runzip runzip_
#define rziplen rziplen_
#define izip izip_
#define iunzip iunzip_
#define iziplen iziplen_
#endif
#endif

void rzip(INT       *p_N,
          double    *p_Thrs,
          INT       *p_Bytes,
          double    *p_X,
          short int *p_Y)
{

  INT        N;
  INT        i;
  INT        inci;
  INT        j;
  INT        jmax;
  unsigned INT        zip;
  unsigned INT        mask;
  unsigned INT        shift;
  double     value;
  double     abs_value;
  int        i_value;
  short int  si_value;
  short int *p_zip;
  short int *p_data;
  double     thrs1;
  double     thrs2;
  double     thrs3;
  double     scale;

  N = *p_N;
  inci = 4*sizeof(INT     );
  p_data = p_Y;
  thrs1 = 0.5*(*p_Thrs);
  thrs2 = (double)(0x7FFC)*thrs1;
  thrs3 = (double)(0x7FFFFFFC)*thrs1;
  scale = (double)(1)/thrs1;

  for ( i=0; i<N; i+=inci )
  {
    p_zip = p_data;
    p_data += sizeof(INT     )/sizeof(short int);
    zip = 0;
    shift = 1;
    jmax = inci-1;
    if ( (i+jmax) >= N ) jmax = N-i-1;
    for ( j=0; j<=jmax; j++ )
    {
      value = *p_X; p_X++;
      abs_value = value;
      if ( value < (double)(0) ) abs_value = -value;
      if ( abs_value < thrs1 )
      {
        mask = 0;
        zip = zip+mask*shift; shift = shift*4;
        continue;
      }
      else
      {
        if ( abs_value < thrs2 )
        {
          mask = 1;
          zip = zip+mask*shift; shift = shift*4;
          si_value = (short int)(value*scale);
          *p_data = si_value;
          p_data++;
          continue;
        }
        else
        {
          if ( abs_value < thrs3 )
          {
            mask = 2;
            zip = zip+mask*shift; shift = shift*4;
            i_value = (int)(value*scale);
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
            (void)memcpy(p_data,&i_value,sizeof(double));
#else
            (*(int *)(p_data)) = i_value;
#endif
            p_data += (sizeof(int)/sizeof(short int));
            continue;
          }
          else
          {
            mask = 3;
            zip = zip+mask*shift; shift = shift*4;
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
            (void)memcpy(p_data,&value,sizeof(double));
#else
            (*(double *)(p_data)) = value;
#endif
            p_data += (sizeof(double)/sizeof(short int));
            continue;
          }
        }
      }
    }
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
    (void)memcpy(p_zip,&zip,sizeof(INT     ));
#else
    (*(INT      *)(p_zip)) = zip;
#endif
  }
  *p_Bytes = (INT     )(sizeof(short int)*(p_data-p_Y));
}

void runzip(INT       *p_N,
            double    *p_Thrs,
            INT       *p_Bytes,
            short int *p_X,
            double    *p_Y)
{

  INT        N;
  INT        i;
  INT        inci;
  INT        j;
  INT        jmax;
  unsigned INT        zip;
  unsigned INT        mask;
  double     value;
  int        i_value;
  short int  si_value;
  short int *p_zip;
  short int *p_data;
  double     thrs;

  N = *p_N;
  inci = 4*sizeof(INT     );
  p_data = p_X;
  thrs = 0.5*(*p_Thrs);

  for ( i=0; i<N; i+=inci )
  {
    p_zip = p_data;
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
    (void)memcpy(&zip,p_zip,sizeof(INT     ));
#else
    zip = (*(INT      *)(p_zip));
#endif
    p_data += sizeof(INT     )/sizeof(short int);
    jmax = inci-1;
    if ( (i+jmax) >= N ) jmax = N-i-1;
    for ( j=0; j<=jmax; j++ )
    {
      mask = 3;
      mask = zip&mask;
      zip = zip/4;
      if ( mask == 0 )
      {
        value = (double)(0);
        *p_Y = value; p_Y++;
        continue;
      }
      if ( mask == 1 )
      {
        si_value = *p_data;
        value = thrs*((double)si_value);
        p_data++;
        *p_Y = value; p_Y++;
        continue;
      }
      if ( mask == 2 )
      {
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
        (void)memcpy(&i_value,p_data,sizeof(int));
#else
        i_value = (*(int *)(p_data));
#endif
        p_data += (sizeof(int)/sizeof(short int));
        value = thrs*((double)i_value);
        *p_Y = value; p_Y++;
        continue;
      }
      if ( mask == 3 )
      {
#if defined(_DECAXP_) || defined(_IRIX64_) || defined(_SOLARIS_) || defined(_HP_UX_) || defined(_PRIMEPOWER_)
        (void)memcpy(&value,p_data,sizeof(double));
#else
        value = (*(double *)(p_data));
#endif
        p_data += (sizeof(double)/sizeof(short int));
        *p_Y = value; p_Y++;
        continue;
      }
    }
  }
  *p_Bytes = (INT     )(sizeof(short int)*(p_data-p_X));
}

void rziplen(INT       *p_N,
             double    *p_Thrs,
             double    *p_X,
             INT       *p_Y)
{

  INT        N;
  INT        i;
  INT        inci;
  INT        j;
  INT        jmax;
  double     value;
  double     abs_value;
  double     thrs1;
  double     thrs2;
  double     thrs3;

  N = *p_N;
  inci = 4*sizeof(INT     );
  thrs1 = 0.5*(*p_Thrs);
  thrs2 = (double)(0x7FFC)*thrs1;
  thrs3 = (double)(0x7FFFFFFC)*thrs1;

  for ( i=0; i<N; i+=inci )
  {
    jmax = inci-1;
    if ( (i+jmax) >= N ) jmax = N-i-1;
    for ( j=0; j<=jmax; j++ )
    {
      value = *p_X; p_X++;
      abs_value = value;
      if ( value < (double)(0) ) abs_value = -value;
      if ( abs_value < thrs1 )
      {
        *p_Y = (INT     )0;
      }
      else
      {
        if ( abs_value < thrs2 )
        {
          *p_Y = (INT     )sizeof(short int);
        }
        else
        {
          if ( abs_value < thrs3 )
          {
            *p_Y = (INT     )sizeof(int);
          }
          else
          {
            *p_Y = (INT     )sizeof(double);
          }
        }
      }
      if ( j == 0 ) *p_Y = *p_Y + sizeof(INT     );
      p_Y++;
    }
  }
}

void izip(INT      *OpCode,
          INT      *nData,
          INT      *nBytes,
          INT      *InBuf,
          INT      *OutBuf)

{

  INT       byte_1,byte_2,byte_3,byte_4,byte_5,byte_6,byte_7,byte_8,byte_9,byte_10;
  char      *p_byte_1,*p_byte_2,*p_byte_3,*p_byte_4,*p_byte_5,*p_byte_6,*p_byte_7,*p_byte_8,*p_byte_9,*p_byte_10;
  INT       prev_number;
  INT       number;
  INT       stop_bits;
  size_t    l_char;
  char     *p_zip;
  int       i;
  static char      magic[8] = { 0x7f , 'M' , 'O' , 'L' , 'C' , 'A' , 'S' , 0x7f };
  int      *MAGIC;

  l_char = sizeof(char);
  p_zip = (char *)OutBuf;
  *nBytes = 0;

  MAGIC = (int *)(&magic[0]);
  if ( *MAGIC == 1280265599 )
  {
    p_byte_1 = (char *)(&byte_1);
    p_byte_2 = (char *)(&byte_2);
    p_byte_3 = (char *)(&byte_3);
    p_byte_4 = (char *)(&byte_4);
    p_byte_5 = (char *)(&byte_5);
    p_byte_6 = (char *)(&byte_6);
    p_byte_7 = (char *)(&byte_7);
    p_byte_8 = (char *)(&byte_8);
    p_byte_9 = (char *)(&byte_9);
    p_byte_10 = (char *)(&byte_10);
  }
  else if ( *MAGIC == 2135773004 )
  {
    p_byte_1 = (char *)(&byte_1)+sizeof(INT     )-sizeof(char);
    p_byte_2 = (char *)(&byte_2)+sizeof(INT     )-sizeof(char);
    p_byte_3 = (char *)(&byte_3)+sizeof(INT     )-sizeof(char);
    p_byte_4 = (char *)(&byte_4)+sizeof(INT     )-sizeof(char);
    p_byte_5 = (char *)(&byte_5)+sizeof(INT     )-sizeof(char);
    p_byte_6 = (char *)(&byte_6)+sizeof(INT     )-sizeof(char);
    p_byte_7 = (char *)(&byte_7)+sizeof(INT     )-sizeof(char);
    p_byte_8 = (char *)(&byte_8)+sizeof(INT     )-sizeof(char);
    p_byte_9 = (char *)(&byte_9)+sizeof(INT     )-sizeof(char);
    p_byte_10 = (char *)(&byte_10)+sizeof(INT     )-sizeof(char);
  }
  else
  {
    (void)printf("\n *** Error in pocedure IZIP *** \n");
    (void)printf(" You attempt to run the program on an unknow architecture\n");
    (void)printf("MAGIC = %d\n\n",*MAGIC);
    _exit (16);
  }

  prev_number = 0;
  for ( i=1; i<=*nData; i++ )
  {
    number = *InBuf;
    if ( *OpCode == 1 ) { number = number-prev_number; }
    prev_number = *InBuf; InBuf++;
    stop_bits = 0x80;
    if ( number < 0 )
    {
      number = -number;
      stop_bits = 0xC0;
    }
    byte_1 = number&0x3F;
    byte_1 = byte_1|stop_bits;
    if ( number <= 0x3F )
    {
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+1;
      continue;
    }
    number = number>>0x06;
    byte_2 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+2;
      continue;
    }
    number = number>>0x07;
    byte_3 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+3;
      continue;
    }
    number = number>>0x07;
    byte_4 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+4;
      continue;
    }
    number = number>>0x07;
    byte_5 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+5;
      continue;
    }
    number = number>>0x07;
    byte_6 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_6,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+6;
      continue;
    }
    number = number>>0x07;
    byte_7 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_7,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_6,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+7;
      continue;
    }
    number = number>>0x07;
    byte_8 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_8,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_7,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_6,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+8;
      continue;
    }
    number = number>>0x07;
    byte_9 = number&0x7F;
    if ( number <= 0x7F )
    {
      (void)memcpy(p_zip,p_byte_9,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_8,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_7,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_6,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+9;
      continue;
    }
    else
    {
      byte_10 = number>>0x07;
      (void)memcpy(p_zip,p_byte_10,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_9,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_8,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_7,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_6,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_5,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_4,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_3,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_2,l_char); p_zip++;
      (void)memcpy(p_zip,p_byte_1,l_char); p_zip++;
      *nBytes = *nBytes+10;
      continue;
    }
  }
}

void iunzip(INT      *OpCode,
            INT      *nData,
            INT      *nBytes,
            INT      *InBuf,
            INT      *OutBuf)

{

  INT      byte_1,byte_2,byte_3,byte_4,byte_5,byte_6,byte_7,byte_8,byte_9,byte_10;
  char      *p_byte_1,*p_byte_2,*p_byte_3,*p_byte_4,*p_byte_5,*p_byte_6,*p_byte_7,*p_byte_8,*p_byte_9,*p_byte_10;
  INT       stop_bit;
  INT       sign_bit;
  INT       data_bits;
  size_t    l_char;
  INT      *prev_number;
  INT       number;
  char     *p_unzip;
  int       i;
  static char      magic[8] = { 0x7f , 'M' , 'O' , 'L' , 'C' , 'A' , 'S' , 0x7f };
  int      *MAGIC;

  l_char = sizeof(char);
  p_unzip = (char *)InBuf;
  *nBytes = 0;

  byte_1 = 0;
  byte_2 = 0;
  byte_3 = 0;
  byte_4 = 0;
  byte_5 = 0;
  byte_6 = 0;
  byte_7 = 0;
  byte_8 = 0;
  byte_9 = 0;
  byte_10 = 0;

  MAGIC = (int *)(&magic[0]);
  if ( *MAGIC == 1280265599 )
  {
    p_byte_1 = (char *)(&byte_1);
    p_byte_2 = (char *)(&byte_2);
    p_byte_3 = (char *)(&byte_3);
    p_byte_4 = (char *)(&byte_4);
    p_byte_5 = (char *)(&byte_5);
    p_byte_6 = (char *)(&byte_6);
    p_byte_7 = (char *)(&byte_7);
    p_byte_8 = (char *)(&byte_8);
    p_byte_9 = (char *)(&byte_9);
    p_byte_10 = (char *)(&byte_10);
  }
  else if ( *MAGIC == 2135773004 )
  {
    p_byte_1 = (char *)(&byte_1)+sizeof(INT     )-sizeof(char);
    p_byte_2 = (char *)(&byte_2)+sizeof(INT     )-sizeof(char);
    p_byte_3 = (char *)(&byte_3)+sizeof(INT     )-sizeof(char);
    p_byte_4 = (char *)(&byte_4)+sizeof(INT     )-sizeof(char);
    p_byte_5 = (char *)(&byte_5)+sizeof(INT     )-sizeof(char);
    p_byte_6 = (char *)(&byte_6)+sizeof(INT     )-sizeof(char);
    p_byte_7 = (char *)(&byte_7)+sizeof(INT     )-sizeof(char);
    p_byte_8 = (char *)(&byte_8)+sizeof(INT     )-sizeof(char);
    p_byte_9 = (char *)(&byte_9)+sizeof(INT     )-sizeof(char);
    p_byte_10 = (char *)(&byte_10)+sizeof(INT     )-sizeof(char);
  }
  else
  {
    (void)printf("\n *** Error in pocedure IUNZIP *** \n");
    (void)printf(" You attempt to run the program on an unknow architecture\n");
    (void)printf("MAGIC = %d\n\n",*MAGIC);
    _exit (16);
  }

  prev_number = OutBuf;
  for ( i=1; i<=*nData; i++ )
  {
    number = 0;
    (void)memcpy(p_byte_1,p_unzip,l_char); p_unzip++;
    stop_bit = byte_1&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_1&0x40;
      data_bits = byte_1&0x3F;
      number = data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+1;
      continue;
    }
    number = byte_1;
    (void)memcpy(p_byte_2,p_unzip,l_char); p_unzip++;
    stop_bit = byte_2&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_2&0x40;
      data_bits = byte_2&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+2;
      continue;
    }
    number = number<<0x07|byte_2;
    (void)memcpy(p_byte_3,p_unzip,l_char); p_unzip++;
    stop_bit = byte_3&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_3&0x40;
      data_bits = byte_3&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+3;
      continue;
    }
    number = number<<0x07|byte_3;
    (void)memcpy(p_byte_4,p_unzip,l_char); p_unzip++;
    stop_bit = byte_4&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_4&0x40;
      data_bits = byte_4&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+4;
      continue;
    }
    number = number<<0x07|byte_4;
    (void)memcpy(p_byte_5,p_unzip,l_char); p_unzip++;
    stop_bit = byte_5&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_5&0x40;
      data_bits = byte_5&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+5;
      continue;
    }
    else
    {
      number = number<<0x07|byte_5;
    }
    (void)memcpy(p_byte_6,p_unzip,l_char); p_unzip++;
    stop_bit = byte_6&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_6&0x40;
      data_bits = byte_6&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+6;
      continue;
    }
    else
    {
      number = number<<0x07|byte_6;
    }
    (void)memcpy(p_byte_7,p_unzip,l_char); p_unzip++;
    stop_bit = byte_7&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_7&0x40;
      data_bits = byte_7&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+7;
      continue;
    }
    else
    {
      number = number<<0x07|byte_7;
    }
    (void)memcpy(p_byte_8,p_unzip,l_char); p_unzip++;
    stop_bit = byte_8&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_8&0x40;
      data_bits = byte_8&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+8;
      continue;
    }
    else
    {
      number = number<<0x07|byte_8;
    }
    (void)memcpy(p_byte_9,p_unzip,l_char); p_unzip++;
    stop_bit = byte_9&0x80;
    if ( stop_bit != 0 )
    {
      sign_bit = byte_9&0x40;
      data_bits = byte_9&0x3F;
      number = number<<0x06|data_bits;
      if ( sign_bit != 0 ) number = -number;
      *OutBuf = number; OutBuf++;
      *nBytes = *nBytes+9;
      continue;
    }
    else
    {
      number = number<<0x07|byte_9;
    }
    (void)memcpy(p_byte_10,p_unzip,l_char); p_unzip++;
    sign_bit = byte_10&0x40;
    data_bits = byte_10&0x3F;
    number = number<<0x06|data_bits;
    if ( sign_bit != 0 ) number = -number;
    *OutBuf = number; OutBuf++;
    *nBytes = *nBytes+10;
  }
  if ( *OpCode == 1 )
  {
    OutBuf = prev_number; OutBuf++;
    for ( i=2; i<=*nData; i++ )
    {
      *OutBuf = *OutBuf+*prev_number;
      OutBuf++; prev_number++;
    }
  }

}

void iziplen(INT      *OpCode,
             INT      *nData,
             INT      *InBuf,
             INT      *OutBuf)

{

  INT       prev_number;
  INT       number;
/*char     *p_zip;*/
  int       i;

/*p_zip = (char *)OutBuf;*/

  prev_number = 0;
  for ( i=1; i<=*nData; i++ )
  {
    number = *InBuf;
    if ( *OpCode == 1 ) { number = number-prev_number; }
    prev_number = *InBuf; InBuf++;
    if ( number < 0 )
    {
      number = -number;
    }
    if ( number <= 0x3F )
    {
      *OutBuf = 1; OutBuf++;
      continue;
    }
    number = number>>0x06;
    if ( number <= 0x7F )
    {
      *OutBuf = 2; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 3; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 4; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 5; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 6; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 7; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 8; OutBuf++;
      continue;
    }
    number = number>>0x07;
    if ( number <= 0x7F )
    {
      *OutBuf = 9; OutBuf++;
      continue;
    }
    else
    {
      *OutBuf = 10; OutBuf++;
      continue;
    }
  }
}
