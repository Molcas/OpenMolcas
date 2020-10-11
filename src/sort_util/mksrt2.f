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
* Copyright (C) 1991, Markus P. Fuelscher                              *
*               1991, Per Ake Malmqvist                                *
************************************************************************
      Subroutine MkSrt2
************************************************************************
*                                                                      *
*     Purpose: Inizialize counters and offsets required                *
*              for bin sorting algorithm                               *
*                                                                      *
*     Called from: Sort1                                               *
*                                                                      *
*     Calls to : none                                                  *
*                                                                      *
*     Calling parameters: none                                         *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*                                                                      *
*     Local data declarations: none                                    *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      use srt2
      Implicit Integer (A-Z)
*
#include "srt0.fh"
#include "srt1.fh"
#include "print.fh"
*
      iRout = 80
      iPrint = nPrint(iRout)
      If ( iPrint.gt.10) Write(6,*) ' >>> Enter MKSRT2 <<<'
*
*----------------------------------------------------------------------*
*     initialize various pointers, counters and disk adresses          *
*----------------------------------------------------------------------*
*
      iBin=0
      Do 30 iSyBlk=1,mSyBlk
         nSlice=nSln(iSyBlk)
         If ( nSlice.ne.0 ) then
            Do 40 iSlice=1,nSlice
               iBin=iBin+1
               iDIBin(2,iBin)=-1
               iDVBin(2,iBin)=-1
               iDVBin(3,iBin)=-1
               iDVBin(4,iBin)=-1
               nInt(iBin)=0
               nRec(iBin)=0
40          Continue
         End If
30    Continue
*
      Return
      End
