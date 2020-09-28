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
      Subroutine MkSrt1()
************************************************************************
*                                                                      *
*     Purpose: Determine the splitting of the integrals into           *
*              submatrices (slices) and hence determine the            *
*              number of bins. Also construct the offsets              *
*              which allow to reduce an integral symmetry block        *
*              number into a relative bin number with respect to       *
*              the bin number of the first integral in an given        *
*              symmetry block.                                         *
*                                                                      *
*     Called from: Sort0                                               *
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

#include "SysDef.fh"
#include "print.fh"
#include "warnings.fh"

      iRout = 80
      iPrint = nPrint(iRout)
      If ( iPrint.gt.10) Write (6,*) ' >>> Enter MKSRT1 <<<'
*----------------------------------------------------------------------*
*                                                                      *
*     grab memory                                                      *
*                                                                      *
*     The available memory should not be smaller than the number       *
*     of unique symmetry blocks times the size of a bin.               *
*     If nSyOp=1 mxBin=  1 and If Square=.true. mxBin=  1              *
*        nSyOp=2 mxBin=  4                      mxBin=  5              *
*        nSyOp=4 mxBin= 19                      mxBin= 28              *
*        nSyOp=8 mxBin=106                      mxBin=176              *
*                                                                      *
*----------------------------------------------------------------------*
*
      lSrtA=1
      If ( nSyOp.eq.2 ) then
        lSrtA = 4
        If ( Square ) lSrtA = 5
      Else If ( nSyOp.eq.4 ) then
        lSrtA = 19
        If ( Square ) lSrtA = 28
      Else If ( nSyOp.eq.8 ) then
        lSrtA = 106
        If ( Square ) lSrtA = 176
      End If
      lSrtA = lSrtA*lBin
      lSrtA = ((1+RtoI)*lSrtA)/RtoI
      Call mma_MaxDBLE(MaxMem)
      lSrtA=Max(lSrtA,MaxMem/2)
*
*----------------------------------------------------------------------*
*     determine the partitioning of the 2el integrals into             *
*     submatrices by dividing the symmetry blocks into slices that     *
*     fit into the available memory.                                   *
*----------------------------------------------------------------------*
*
      Do iSyBlk=1,mSyBlk
        nSln(iSyBlk)=0
        lSll(iSyBlk)=0
      End Do
*
#ifndef _I8_
      lim_32=2**30
#endif
      nSym=nSyOp
      Do iSymi=1,nSym
        ib=nBs(iSymi)
        iSkip=nSkip(iSymi)
        If (ib.eq.0.or.iSkip.eq.1) Go To 100
        Do jSymj=1,iSymi
          iSymj=1+ieor(iSymi-1,jSymj-1)
          iSyblj=jSymj+iSymi*(iSymi-1)/2
          jb=nBs(jSymj)
          ibj=ib*jb
          If( iSymi.eq.jSymj ) ibj=ib*(ib+1)/2
          jSkip=nSkip(jSymj)
          If (jb.eq.0.or.jSkip.eq.1) Go To 200
          kSymMx=iSymi
          If( Square ) kSymMx=nSym
          Do kSymk=1,kSymMx
            kb=nBs(kSymk)
            kSkip=nSkip(kSymk)
            If (kb.eq.0.or.kSkip.eq.1) Go To 300
            lSymMx=jSymj
            If( kSymk.ne.iSymi .or. Square ) lSymMx=kSymk
            Do lSyml=1,lSymMx
               kSyml=1+ieor(kSymk-1,lSyml-1)
               kSybll=lSyml+kSymk*(kSymk-1)/2
               If( ieor(iSymj-1,kSyml-1).ne.0 ) Go To 400
*
C              Write (*,*) 'i,j,k,l=',iSymi,jSymj,kSymk,lSyml
               lb=nBs(lSyml)
               lSkip=nSkip(lSyml)
               If (lb.eq.0.or.lSkip.eq.1) Go To 400
               kbl=kb*lb
#ifndef _I8_
               ix=lim_32/kbl
#endif
               If( kSymk.eq.lSyml ) kbl=kb*(kb+1)/2
               iSyBlk=kSybll+mxSyP*(iSyblj-1)
                nij = 0
                Do ij = 1,ibj
#ifdef _I8_
                   If ( (ij*kbl).lt.lSrtA ) nij = ij
#else
                   If (ij.le.ix) Then
                      If ( (ij*kbl).lt.lSrtA ) nij = ij
                   End If
#endif
                End Do
                If ( nij.eq.0 ) then
                    rc=001
                    Write(6,*)
                    Write(6,'(2X,A,I3.3,A)')
     &              '*** Error in MKSRT1 ***'
                    Write(6,'(2X,A)') 'Insufficient memory'
                    Write(6,'(2X,A)') 'Increase the value of the'//
     &                                'MOLCAS_MEM environement variable'
                    Write(6,*)
                    Call qTrace
                    Call Quit(_RC_MEMORY_ERROR_)
                End If
                nSlice = 1+(ibj-1)/nij
                Do ij = nij,1,-1
                   If ( (ij*nSlice).ge.ibj ) nij=ij
                End Do
                nSln(iSyBlk)=nSlice
                lSll(iSyBlk)=nij*kbl
C               Write (*,*) 'iSyBlk=',iSyBlk
C               Write (*,*) 'nSln(iSyBlk)=',nSln(iSyBlk)
C               Write (*,*) 'lSll(iSyBlk)=',lSll(iSyBlk)
 400          Continue
            End Do
 300        Continue
          End Do
 200      Continue
        End Do
 100    Continue
      End Do
*
*----------------------------------------------------------------------*
*     Determine the final amount of memory and bins that will be used. *
*     Check again for consistency                                      *
*----------------------------------------------------------------------*
*
      nBin = 0
      MxSrtA2 = 0
      Do iSyBlk = 1,mSyBlk
        nSlice  = nSln(iSyBlk)
        lSlice  = lSll(iSyBlk)
        nBin    = nBin+nSlice
        MxSrtA2 = Max(MxSrtA2,lSlice)
      End Do
      MxSrtA1 = ((1+RtoI)*nBin*lBin)/RtoI
*
      If (iPrint.gt. 5) Then
      Write (6,*)
      Write(6,'(A,I12,A,I4,A)') '  SEWARD will use a sorting area of',
     &     MxSrtA1,' Words(Real*8) in the first phase (=',nBin,' bins).'
      Write(6,'(A,I12,A)') '  SEWARD will use a sorting area of',
     &     MxSrtA2,' Words(Real*8) in the second phase.'
      Write (6,*)
      End If
*
      If ( nBin.gt.mxBin ) Then
        rc=003
        Write(6,*)
        Write(6,'(2X,A,I3.3,A)')
     &  '*** Error in MKSRT1 ***'
        Write(6,'(2X,A)') 'Insufficient memory'
        Write(6,'(2X,A)') 'nBin exceeds limits (nBin>mxBin)'
        Write(6,*)
        Write(6,*) 'nBin=',nBin
        Write(6,*) 'mxBin=',mxBin
        Write(6,*)
        Write (6,*) 'Increase MOLCAS_MEM and try again!'
        Write(6,*)
        Call qTrace
        Call Quit(_RC_MEMORY_ERROR_)
      End If
*
      If ( MxSrtA1.gt.lSrtA ) Then
        rc=004
        Write(6,*)
        Write(6,'(2X,A,I3.3,A)')
     &  '*** Error in MKSRT1 ***'
        Write(6,'(2X,A)') 'Insufficient memory'
        Write(6,'(2X,A)') 'MxSrtA1>lSrtA'
        Write(6,*)
        Write(6,*) 'MxSrtA1=',MxSrtA1
        Write(6,*) 'lSrtA=',lSrtA
        Write(6,*)
        Write (6,*) 'Increase MOLCAS_MEM and try again!'
        Write(6,*)
        Call qTrace
        Call Quit(_RC_MEMORY_ERROR_)
      End If
*
*----------------------------------------------------------------------*
*     compute offsets                                                  *
*----------------------------------------------------------------------*
*
      iOff=1
      Do iSyBlk=1,mSyBlk
        nSlice=nSln(iSyBlk)
        IstBin(iSyBlk)=iOff
        iOff=iOff+nSlice
      End Do
*
*----------------------------------------------------------------------*
*     initialize counter for integrals per slice                       *
*----------------------------------------------------------------------*
*
      Call ICopy(3*nBin,[0],0,mInt,1)
*
      Return
      End
