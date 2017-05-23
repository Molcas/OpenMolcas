************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Per-Olof Widmark                                 *
*               1993, Markus P. Fuelscher                              *
*               1996, Luis Serrano-Andres                              *
*               2012, Victor P. Vysotskiy                              *
************************************************************************
      Subroutine FIOInit()
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     set initial values of all arguments in the common blocks         *
*     /FIO1/, /FIO2/ and /FIO3/                                        *
*     /PFIO1/,/PFIO2/,/PFIO3/,/PFIO4/                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
#include "fio.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#else
#include "pfio.fh"
#endif
#include "ComDat.fh"
*
      Do i = 1, MxFile
         isOpen(i)=0
         FSCB(i)=0
         Addr(i)=0
         Multi_File(i)=.False.
#ifdef _OLD_IO_STAT_
         MxAddr(i)=0
         Do j = 1, 4
            Count(j,i)=0
         End Do
#else
         Do j=1,8
            PRofData(j,i)=0.d0
         End Do
#endif
         LuName(i)(1:2)='FT'
         Write (LuName(i)(3:4),'(I2.2)') i
         LuName(i)(5:8)='F001'
         LuMark(i)=-1
         Do j = 0, MaxSplitFile-1
            MPUnit(j,i)=0
         End Do
         MBL(i)=0
      End Do
*
#ifndef _OLD_IO_STAT_
      NProfFiles=0
#endif
      MaxFileSize=0
      FirstCall=.true.
      Trace=.false.
      Query=.false.
*
*---- Initiate the COMFILE as closed
*
      AuxCom(pOpen)=NaN
*
      Return
      End
