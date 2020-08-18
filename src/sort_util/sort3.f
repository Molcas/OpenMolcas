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
* Copyright (C) 1993,1996,1998, Markus P. Fuelscher                    *
*               1993, Per Ake Malmqvist                                *
*               1998, Roland Lindh                                     *
************************************************************************
      Subroutine SORT3(MaxDax)
************************************************************************
*                                                                      *
*     Purpose: Up to this point the two electron integral file         *
*              has been used like a random stack of records. To        *
*              save hunting for records when reading the records       *
*              bring them into sequential order.                       *
*                                                                      *
*     Called from: Seward_main                                         *
*                                                                      *
*     Calls to : DaFile,DCopy,MkOrd,ClsOrd,ErrOrd                      *
*                                                                      *
*     Calling parameters:                                              *
*     MaxDax  : Higest disk adress of the final 2el integral file      *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     TwoDat  : definitions of sorting flags and address tables        *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*                                                                      *
*     local data declarations:                                         *
*     Buf1   : I/O buffer contains packed integral values              *
*     Buf2   : I/O buffer contains packed integral values              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.-AA. Malmqvist                             *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     - modified to use a virtual disk                                 *
*       M. P. Fuelscher, University of Lund, Sweden, 1996              *
*     - modified to get rid of all disk address calculations           *
*       M. P. Fuelscher, University of Lund, Sweden, 1998              *
*     - modified to eliminate copy statement                           *
*       R. Lindh, University of Lund, Sweden, 1998                     *
*     - modified to eliminate free records                             *
*       R. Lindh, University of Lund, Sweden, 1998                     *
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

#include "SysDef.fh"
#include "print.fh"
#include "stdalloc.fh"

      Dimension Buf(2*lStRec)
      Integer, Allocatable:: SrtKey(:), SrtAdr(:)

*
*----------------------------------------------------------------------*
*     pick up the print level                                          *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
      iRout = 88
      iPrint = nPrint(iRout)
      If ( iPrint.gt.5 ) Write(6,*) ' >>> Enter SORT3 <<<'
#endif
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
      Call qEnter('Sort3')
*
*----------------------------------------------------------------------*
*     Scan once the two-electron integral file a pick up the sort      *
*     key as well as disk addresses                                    *
*----------------------------------------------------------------------*
*
      Call mma_allocate(SrtKey,mxOrd,Label='SrtKey')
      Call mma_allocate(SrtAdr,mxOrd,Label='SrtAdr')
      iOpt = 2
      iDisk = iDaTw0
      MaxDax=0
      Do iOrd = 1,mxOrd
         SrtAdr(iOrd) = iDisk
         MaxDax=Max(iDisk,MaxDax)
         Call dDAFILE(LuTwo,iOpt,Buf,lStRec,iDisk)
         SrtKey(iOrd) = Int(Buf(2))
      End Do
      MaxDax=iDisk
#ifdef _DEBUG_
      If ( iPrint.ge.10 ) then
        Call iVcPrt('Sort keys',' ',SrtKey,mxOrd)
        Call iVcPrt('Disk addresses',' ',SRtAdr,mxOrd)
      End If
#endif
*
*----------------------------------------------------------------------*
*     Sort records in ascending order of the sort key                  *
*----------------------------------------------------------------------*
*
      iWr=1
      iRd=2
*
*---- Loop over all records
*
      Do i = 1,mxOrd
        j1 = i
        j2 = SrtKey(i)
        If ( j2.ne.i ) then
          iB1 = 1
          iB2 = lStRec+1
          iDisk = SrtAdr(j1)
          Call dDAFILE(LuTwo,iRd,Buf(iB1),lStRec,iDisk)
          Do while ( j2.ne.i )
            iDisk = SrtAdr(j2)
            Call dDAFILE(LuTwo,iRd,Buf(iB2),lStRec,iDisk)
            iDisk = SrtAdr(j2)
            Call dDAFILE(LuTwo,iWr,Buf(iB1),lStRec,iDisk)
            iTmp=iB2
            iB2 = iB1
            iB1 = iTmp
            j1 = j2
            j2 = SrtKey(j2)
            SrtKey(j1) = j1
          End Do
          iDisk = SrtAdr(j2)
          Call dDAFILE(LuTwo,iWr,Buf(iB1),lStRec,iDisk)
          SrtKey(j2) = j2
        End If
      End Do
#ifdef _DEBUG_
      If ( iPrint.ge.10 ) then
        Call iVcPrt('Sort keys',' ',SrtKey,mxOrd)
      End If
#endif
*
*----------------------------------------------------------------------*
*     Update the disk start adressed of each slice                     *
*----------------------------------------------------------------------*
*
      j = 1
      Do iBin=1,nBin
        iDVBin(2,iBin) = SrtAdr(j)
        j = j+nRec(iBin)
      End Do
      Call mma_deallocate(SrtKey)
      Call mma_deallocate(SrtAdr)
*
*----------------------------------------------------------------------*
*     Write the final table of content to disk                         *
*----------------------------------------------------------------------*
*
      Call MkOrd(iDummy)
*
*----------------------------------------------------------------------*
*     Close the ordered 2el integral file                              *
*----------------------------------------------------------------------*
*
      iRc=-1
      iOpt=0
      Call ClsOrd(iRc,iOpt)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'SORT3: Error closing ORDINT'
         Call Abend()
      End If
      Call DaClos(LuTmp)
*
*----------------------------------------------------------------------*
*     Release RAMD                                                     *
*----------------------------------------------------------------------*
*
      If (RAMD) Call GetMem('RAMD','Free','Real',ip_RAMD,RAMD_Size)
*
*----------------------------------------------------------------------*
*     Turn timing OFF and exit                                         *
*----------------------------------------------------------------------*
*
      Call qExit('Sort3')
      Return
      End
