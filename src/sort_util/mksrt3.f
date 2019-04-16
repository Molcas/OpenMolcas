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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine MkSrt3(iRc,iSquar,nIrrep,nBas,nSkip)
************************************************************************
*                                                                      *
*    Purpose: Prepare the address table for the virtual disk and       *
*             initialize the memory                                    *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat : table of contents and auxiliary information              *
*    TowRc  : Table of return code                                     *
*                                                                      *
*    calling arguments:                                                *
*    iRc    : return code                                              *
*             A value of 0 (zero) is returned upon successful          *
*             completion of the request. A nonzero value indi-         *
*             cates an error.                                          *
*    iSquar  : ordering flag                                           *
*    nIrrep  : number of irreducible representations                   *
*    nBas    : number of basis functions per irred. rep.               *
*    nSkip   : flag to exlude symmetry combinations                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
      Dimension nBas(nIrrep),nSkip(nIrrep)
*
#include "srt1.fh"
#include "Molcas.fh"
#include "TwoDat.fh"
#include "print.fh"
      Integer nSln_Temp(mSyBlk), lSll_Temp(mSyBlk)
      Logical Square

      iRout = 80
      iPrint = nPrint(iRout)
      If ( iPrint.gt.5 ) Write(6,*) ' >>> Enter MKSRT3 <<<'

      Square = .true.
      If( iSquar.eq.0 )  Square=.false.
*
      Call iCopy(mSyBlk,nSln,1,nSln_Temp,1)
      Call iCopy(mSyBlk,lSll,1,lSll_Temp,1)
      Do iSyBlk=1,mSyBlk
        nSln(iSyBlk)=0
        lSll(iSyBlk)=0
      End Do
*
*     loop over all symmetry blocks
*
*     Write(*,*)
*     Write(*,'(2X,A,I3.3,A)') '*** (I)-level message MKSRT3 ***'
*     Write(*,'(2X,A,I8    )') 'RAMD_size   =',RAMD_size
*     Write(*,'(2X,A,I8    )') 'RAMD_anchor =',RAMD_anchor
*     Write(*,*)
*     Write (*,'(2X,84A1)') ('*',i=1,84)
*
      iRc   = 0
      iOff  = RAMD_anchor
      nInts = 0
      nSym=TocTwo(isSym)
      mxSyP=nSym*(nSym+1)/2
      Do iSymi=1,nSym
        ib=nBas(iSymi)
        iSkip=nSkip(iSymi)
        Do jSymj=1,iSymi
          iSymj=1+ieor(iSymi-1,jSymj-1)
          iSyblj=jSymj+iSymi*(iSymi-1)/2
          jb=nBas(jSymj)
          ibj=ib*jb
          If( iSymi.eq.jSymj ) ibj=ib*(ib+1)/2
          jSkip=nSkip(jSymj)
          kSymMx=iSymi
          If( Square ) kSymMx=nSym
          Do kSymk=1,kSymMx
            kb=nBas(kSymk)
            kSkip=nSkip(kSymk)
            lSymMx=jSymj
            If( kSymk.ne.iSymi .or. Square ) lSymMx=kSymk
            Do lSyml=1,lSymMx
              kSyml=1+ieor(kSymk-1,lSyml-1)
              kSybll=lSyml+kSymk*(kSymk-1)/2
              If( ieor(iSymj-1,kSyml-1).eq.0 ) then
                lb=nBas(lSyml)
                kbl=kb*lb
                If( kSymk.eq.lSyml ) kbl=kb*(kb+1)/2
                lSkip=nSkip(lSyml)
                iSyBlk=kSybll+mxSyP*(iSyblj-1)
                If ( (iSkip+jSkip+kSkip+lSkip).eq.0 .and.
     &               (ibj*kbl).ne.0 ) then
C                  Write (6,*) 'iSymi,jSymj,kSymk,lSyml=',
C    &                          iSymi,jSymj,kSymk,lSyml
C                  Write (6,*) 'iSyBlk=',iSyBlk
C                  Write (6,*) 'nBatch(iSyBlk)=',nBatch(iSyBlk)
*     check if there is enough space available
*
                   lBuf = ibj*kbl
                   nInts = nInts + lBuf
                   If ( nInts.ge.RAMD_size ) then
                      iRc = 001
                      Write(6,*)
                      Write(6,'(2X,A,I3.3,A)')
     &                '*** (W)-level message MKSRT3',iRc,' ***'
                      Write(6,'(2X,A)')
     &                'There is not enough space on the RAM disk'
                      Write(6,'(2X,A)')
     &                'The program will resume normal activity'
                      Write(6,*)
                      Call iCopy(mSyBlk,nSln_Temp,1,nSln,1)
                      Call iCopy(mSyBlk,lSll_Temp,1,lSll,1)
                      Return
                   End If
*
*     save start address of this symmetry block
*
                   iBatch = nBatch(iSyBlk)
                   RAMD_adr(iBatch) = iOff
                   nSln(iSyBlk)=1
                   lSll(iSyBlk)=lBuf
*                  Write (*,'(2X,A,4I2,A,2I8,A,2I8)')
*    &                   ' iSym,jSym,kSym,lSym',
*    &                     iSymi,jSymj,kSymk,lSyml,
*    &                   ' lBuf,nInts',
*    &                     lBuf,nInts,
*    &                   ' iBatch,iOff',
*    &                     iBatch,iOff
*
*     Init integrals
*
                   Call dCopy_(lBuf,[0.0d0],0,RAMD_ints(iOff),1)
*
*     update pointers
*
                   iOff = iOff + lBuf

                End If
              End If
            End Do
          End Do
        End Do
      End Do
*     Write (*,'(2X,84A1)') ('*',i=1,84)
*
*----------------------------------------------------------------------*
*     compute offsets                                                  *
*----------------------------------------------------------------------*
*
*     Note that we fake the information so that the iSyBlk index is
*     really the iSyBlk index as it will be needed in Sort1c!
*
      Do iSyBlk=1,mSyBlk
         iStBin(iSyBlk)=iSyBlk
      End Do
*
*     define initial load point
*
*     RAMD_next=RAMD_adr(1)
*
      If ( iPrint.gt.5 ) Write(6,*) ' >>> Exit MKSRT3 <<<'
      Return
      End
