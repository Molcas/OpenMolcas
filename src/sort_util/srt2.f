************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*---------------------------------------------------------------------*
*                                                                     *
*     Information pertinent to the bin sort algorithm                 *
*                                                                     *
*     *--------------------------------------------------------*      *
*     *                                                        *      *
*     *  Note:                                                 *      *
*     *  These definitions depend on the record structure of   *      *
*     *       2el integral file which are given in the         *      *
*     *         common /TWODEF/ of the RDORD utility           *      *
*     *                                                        *      *
*     *--------------------------------------------------------*      *
*                                                                     *
*     Parameter definitions:                                          *
*     mxBin  : maximum number of bin allowed. The number should       *
*              not be smaller than:                                   *
*              If nSyOp=1 mxBin=  1 and If Square=.true. mxBin=  1    *
*                 nSyOp=2 mxBin=  4                      mxBin=  5    *
*                 nSyOp=4 mxBin= 19                      mxBin= 28    *
*                 nSyOp=8 mxBin=106                      mxBin=176    *
*     lBin   : buffer length of a bin                                 *
*              Note: if lBin is not a multiple of 4*lDaRec            *
*              the buffers used for intermediate storing of           *
*              labels on LuTmp are used very inefficiently            *
*     mxSrtA : maximum sorting area to be used                        *
*                                                                     *
*     Entries to common SRT2:                                         *
*     lwIBin : array used to store index Bins                         *
*     lwVBin : array  used to store value Bins                        *
*     nRec   : number of record per slice                             *
*     iDIBin : disk adresses of index bins                            *
*     iDVBin : disk adresses of value bins                            *
*     nByte  : packed size of value bins                              *
*     nInt   : number of integrals in a bin                           *
*     mInt   : number of integrals and bytes in a slice               *
*     nOff1  : memory allocation offset for index bins                *
*     nOff2  : memory allocation offset for value bins                *
*     LuTwo  : logical unit number of ordered 2el file                *
*     LuTmp  : logical unit number of temporary file                  *
*     iDaTw0 : first disk adress after header of ordered 2el file     *
*     iDaTwo : current disk position of LuTwo                         *
*     iDaTmp : current disk position of LuTmp                         *
*     mDaTwo : highest accessed disk position of LuTwo                *
*     mDaTwo : highest accessed disk position of LuTmp                *
*     MxOrd  : total number of records on LuTwo                       *
*                                                                     *
*---------------------------------------------------------------------*
*
      Module Srt2
#include "TwoDef.fh"
      Parameter ( mxBin = 2048 )
      Parameter ( lBin_tce   = 4*lDaRec          )
      Parameter ( lBin_rle   =32*lDaRec          )
*
      Integer iDIBin(3,mxBin),iDVBin(4,mxBin)
      Integer nRec(mxBin),nByte(mxBin),nInt(mxBin),mInt(3,mxBin)
      Integer nOffV(mxBin),nOffI(mxBin)
*
      Integer LuTwo,LuTmp,iDaTw0,iDaTwo,iDaTmp,mDaTwo,mDaTmp,MxOrd,lbin
*
      Real*8, Allocatable:: ValBin(:)
      Integer, Allocatable:: IndBin(:), lIndx(:), lInts(:)
      Integer, Allocatable:: lwIBin(:)
      Real*8,  Allocatable:: lwVBin(:)
      End Module Srt2
