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
* Copyright (C) 1993,1996, Markus P. Fuelscher                         *
*               1993, Per Ake Malmqvist                                *
************************************************************************
      Subroutine SORT1B
************************************************************************
*                                                                      *
*     Purpose: Phase 1 of the bin sorting algorithm                    *
*              The integral generation is completed and the remaining  *
*              2el integrals stored in the bins are dumped to disk.    *
*                                                                      *
*     Called from: Seward_main                                         *
*                                                                      *
*     Calls to : PKI4,PKR8,SetVec,ErrTra,ISORTX,I4Len,R8Len            *
*                                                                      *
*     Calling parameters: none                                         *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     TwoDat : definitions of sorting flags and address tables         *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*     WSpc    : dynamic work space                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.-AA. Malmqvist                             *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - modified to use a virtual disk                                 *
*       M. P. Fuelscher, University of Lund, Sweden, 1996              *
*                                                                      *
************************************************************************
*
      use srt2
      Implicit Real*8 (A-H,O-Z)
*
#include "Molcas.fh"
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "stdalloc.fh"
#include "print.fh"
*
*----------------------------------------------------------------------*
*     pick up print level                                              *
*----------------------------------------------------------------------*
*
      iRout = 82
      iPrint = nPrint(iRout)
      If ( iPrint.ge.99 ) Write(6,*) ' >>> Enter SORT1B <<<'
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
*
*----------------------------------------------------------------------*
*     dump remaining integrals to disk                                 *
*----------------------------------------------------------------------*
*
      iOpt=0 ! Always tight!
      Do iBin=1,nBin
        Do while ( nInt(iBin).gt.0 )
          Call SaveBin(iBin,iOpt)
        End Do
      End Do
*
*----------------------------------------------------------------------*
*     release the work space used to store bins                        *
*----------------------------------------------------------------------*
*
      call mma_deallocate(lwVBin)
      call mma_deallocate(lwIBin)
*
      call mma_deallocate(lIndx)
      call mma_deallocate(lInts)
*
*----------------------------------------------------------------------*
*     Turn timing OFF and exit                                         *
*----------------------------------------------------------------------*
*
      Return
      End
