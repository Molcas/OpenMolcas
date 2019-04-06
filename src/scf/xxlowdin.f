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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
      SubRoutine xxlowdin(S,CMO,nBas,nSym)
************************************************************************
*                                                                      *
* This routine computes S^(-1/2), i.e. Lowdin orthogonalization.       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         Lund University, Sweden                                      *
*                                                                      *
************************************************************************
      Implicit none
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8 S(*)
      Real*8 CMO(*)
      Integer nSym
      Integer nBas(nSym)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Real*8  t
      Integer iSym
      Integer MaxBas
      Integer MaxTri
      Integer MaxSqr
      Integer nSqr
      Integer nTri
      Integer iSqrOff
      Integer iTriOff

      Integer npSmat, npVect, npEige
      Integer i,j,k, ki, kj, ktemp
      Real*8, Dimension(:), Allocatable:: SMat, Vect, Eige
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
      MaxBas=0
      Do iSym=1,nSym
         MaxBas=Max(MaxBas,nBas(iSym))
      End Do
      MaxTri=MaxBas*(MaxBas+1)/2
      MaxSqr=MaxBas*MaxBas
*----------------------------------------------------------------------*
* Allocate                                                             *
*----------------------------------------------------------------------*
      npSmat=MaxTri
      Call mma_allocate(Smat,npSmat,Label='Smat')
      npVect=MaxSqr
      Call mma_allocate(Vect,npVect,Label='Vect')
      npEige=MaxBas
      Call mma_allocate(Eige,npEige,Label='Eige')
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      iSqrOff=0
      iTriOff=0
      Do iSym=1,nSym
*        Write(6,'(a,2i5)') 'iSqrOff,iTriOff',iSqrOff,iTriOff
         nSqr=nBas(iSym)*nBas(iSym)
         nTri=nBas(iSym)*(nBas(iSym)+1)/2
         If(nBas(iSym).le.0) Goto 101
         Call dCopy_(nTri,S(1+iTriOff),1,Smat,1)
         Call FZero(Vect,nSqr)
         Call dCopy_(nBas(iSym),[1.0D0],0,Vect,nBas(iSym)+1)
         Call NIdiag(Smat,Vect,nBas(iSym),nBas(iSym),0)
         Do i=1,nBas(iSym)
*           Eige(i)=Smat(+i*(i+1)/2)
            Eige(i)=1.0d0/Sqrt(Smat(+i*(i+1)/2))
         End Do
         Do i=1,nBas(iSym)
            Do j=1,nBas(iSym)
               t=0.0d0
               Do k=1,nBas(iSym)
                  ktemp = (k-1)*nBas(iSym)
                  ki = ktemp + i
                  kj = ktemp + j
                  t=t+Vect(ki) * Eige(k) * Vect(kj)
               End Do
               CMO(i+(j-1)*nBas(iSym)+iSqrOff)=t
            End Do
         End Do
101      Continue
         iSqrOff=iSqrOff+nSqr
         iTriOff=iTriOff+nTri
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call mma_deallocate(Eige)
      Call mma_deallocate(Vect)
      Call mma_deallocate(SMat)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
