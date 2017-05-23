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
      Subroutine RdOrd_(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
************************************************************************
*                                                                      *
*    Purpose: Read a buffer of ordered integrals of length lBuf        *
*             The subroutine returns the number of submatrices         *
*             it has picked up.                                        *
*                                                                      *
*    Calling parameters:                                               *
*    iSym   : irred. representation of first symmetry label            *
*    jSym   : irred. representation of second symmetry label           *
*    kSym   : irred. representation of third symmetry label            *
*    lSym   : irred. representation of fourth symmetry label           *
*    Buf    : containts on output the integrals                        *
*    lBuf   : length of the integral buffer                            *
*    nMat   : number of submatrices read in                            *
*    iOpt   : option code (iOpt=1:start reading at first integral)     *
*                         (iOpt=2:continue reading)                    *
*    rc     : return code                                              *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat : table of contents and auxiliary information              *
*    TowRc  : Table of return code                                     *
*                                                                      *
*    Local data declarations:                                          *
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
#include "Molcas.fh"
#include "TwoDat.fh"
#include "TwoRc.fh"
*
      Real*8 Buf(*)
      Logical Square
*----------------------------------------------------------------------*
*     Start the procedure                                              *
*----------------------------------------------------------------------*
*     Call qEnter('RdOrd')
      rc=rc0000
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
      LuTwo=AuxTwo(isUnit)
      Open=AuxTwo(isStat)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      If ( Open.ne.1 ) Then
        rc=rcRD10
        Write (6,*) 'RdOrd: ORDINT not opened yet!'
        Call QTrace()
        Call Abend()
      End If
*----------------------------------------------------------------------*
*     Check if the packing table has been loaded                       *
*----------------------------------------------------------------------*
      If ( TocTwo(isPkPa).lt.0 .or. TocTwo(isPkPa).gt.1 .or.
     &     TocTwo(isPkAs).lt.0 .or. TocTwo(isPkAs).gt.1      ) then
        rc=rcRD11
        Write (6,*) 'RdOrd: the packing flags are spoiled'
        Call QTrace()
        Call Abend()
      End If
*---------------------------------------------------------------------*
*     Check the symmetry labels                                       *
*---------------------------------------------------------------------*
      Square=TocTwo(isOrd).eq.1
      If ( ieor(iSym-1,jSym-1).ne.ieor(kSym-1,lSym-1) ) then
        rc=rcRD01
        Write (6,*) 'RdOrd: Wrong symmetry labels, direct product',
     &              ' is not total symmetric'
        Call QTrace()
        Call Abend()
      End If
      If ( iSym.lt.jSym .or. kSym.lt.lSym ) then
        rc=rcRD02
        Write (6,*) 'RdOrd: invalid order of symmetry labels'
        Call QTrace()
        Call Abend()
      End If
      ijS=jSym+iSym*(iSym-1)/2
      klS=lSym+kSym*(kSym-1)/2
      If ( ijS.lt.klS .and. .not.Square ) then
        rc=rcRD03
        Write (6,*) 'RdOrd: invalid combination of symmetry labels'
        Call QTrace()
        Call Abend()
      End If
      nSym=TocTwo(isSym)
      nPairs=nSym*(nSym+1)/2
      iSyBlk=(ijS-1)*nPairs+klS
      iBatch=nBatch(iSyBlk)
*---------------------------------------------------------------------*
*     Check the skipping flags                                        *
*---------------------------------------------------------------------*
      iSkip=TocTwo(isSkip+iSym-1)
      jSkip=TocTwo(isSkip+jSym-1)
      kSkip=TocTwo(isSkip+kSym-1)
      lSkip=TocTwo(isSkip+lSym-1)
      If ( (iSkip+jSkip+kSkip+lSkip).ne.0 ) then
        rc=rcRD07
        Write (6,*) 'RdOrd: Requested symmetry block has not been',
     &              ' computed'
        Call QTrace()
        Call Abend()
      End If
*---------------------------------------------------------------------*
*     Check options                                                   *
*---------------------------------------------------------------------*
      If ( iOpt.ne.1 .and. iOpt.ne.2 ) then
        rc=rcRD06
        Write (6,*) 'RdOrd: Invalid option'
        Write (6,*) 'iOpt=',iOpt
        Call QTrace()
        Call Abend()
      End If
*---------------------------------------------------------------------*
*     Check the buffer size                                           *
*---------------------------------------------------------------------*
      If ( lBuf.le.0  ) then
        rc=rcRD04
        Write (6,*) 'RdOrd: invalid buffer size'
        Write (6,*) 'lbuf=',lBuf
        Call QTrace()
        Call Abend()
      End If
*---------------------------------------------------------------------*
*     Compute matrix dimensions                                       *
*---------------------------------------------------------------------*
      iB=TocTwo(isBas+iSym-1)
      jB=TocTwo(isBas+jSym-1)
      kB=TocTwo(isBas+kSym-1)
      lB=TocTwo(isBas+lSym-1)
      ijB=iB*jB
      if ( iSym.eq.jSym ) ijB=jB*(iB+1)/2
      klB=kB*lB
      if ( kSym.eq.lSym ) klB=lB*(kB+1)/2
*---------------------------------------------------------------------*
*     Compute submatrix dimensions                                    *
*---------------------------------------------------------------------*
      If ( lBuf.le.0  ) then
        rc=rcRD04
        Write (6,*) 'RdOrd: invalid buffer size'
        Write (6,*) 'lbuf=',lBuf
        Call QTrace()
        Call Abend()
      End If
      If (klB.le.0) Then
         nMat=0
*        Call qExit('RdOrd')
         Return
      Else
         nMat=(lBuf-1)/klB
      End If
      If ( nMat.gt.ijB ) nMat=ijB
      If ( nMat.eq.0  ) then
        rc=rcRD05
        Write (6,*) 'RdOrd: too small buffer'
        Write (6,*) 'Buffer size is lBuf  =',lBuf
        Write (6,*) 'Size of submatrix klB=',klB
        Write (6,*) 'Call parameters to rdord are:'
        Write (6,*) 'iOpt=',iOpt
        Write (6,*) 'iSym=',iSym
        Write (6,*) 'jSym=',jSym
        Write (6,*) 'kSym=',kSym
        Write (6,*) 'lSym=',lSym
        Write (6,*) 'lBuf=',lBuf
        Write (6,*) 'nMat=',nMat
        Write (6,*) 'Symmetry block iSyBlk=',iSyBlk
        Write (6,*) 'Batch nr       iBatch=',iBatch
        Write (6,*) 'iB=TocTwo(isBas+iSym-1), etc:'
        Write (6,*) 'iB=',iB
        Write (6,*) 'jB=',jB
        Write (6,*) 'kB=',kB
        Write (6,*) 'lB=',lB
        Call QTrace()
        Call Abend()
      End If
*---------------------------------------------------------------------*
*     Check that the number of submatrices do not run beyond          *
*     the last integral of a symmetry allowed batch                   *
*---------------------------------------------------------------------*
      Leftpq=AuxTwo(isNpq)
      If ( iOpt.eq.1 ) then
        Leftpq=ijB
      Else
        nMat=min(Leftpq,nMat)
      End If
      Leftpq=Leftpq-nMat
      AuxTwo(isNpq)=Leftpq
      nInts=nMat*klB
*---------------------------------------------------------------------*
*     Transfer integrals                                              *
*---------------------------------------------------------------------*
      If ( RAMD ) then
         Call ORDIN2(iOpt,Buf,nInts,iBatch)
      Else
         Call ORDIN1(iOpt,Buf,nInts,iBatch)
      End If
*---------------------------------------------------------------------*
*     exit                                                            *
*---------------------------------------------------------------------*
*     Call qExit('RdOrd')
      Return
      End
