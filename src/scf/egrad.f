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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2017, Roland Lindh                                     *
************************************************************************
      SubRoutine EGrad(O,T,V,S,D,nOTSD,C,nC,G,nG,nD,CMO)
************************************************************************
*                                                                      *
*     purpose: This routine calculates the gradient of the SCF energy  *
*              with respect to the rotation parameters.                *
*              Grad(E) = C(t)(FDS-SDF)C                                *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       O       : one-electron hamiltonian of length nOTSD             *
*       T       : two-electron hamiltonian of length nOTSD             *
*       V       : external potential       of length nOTSD             *
*       S       : overlap in AO basis of length nOTSD                  *
*       D       : density matrix of length nOTSD                       *
*       C       : matrix transforming to the set of orthonormal        *
*                 (and spherical, if needed) functions of length nC    *
*                                                                      *
*     output:                                                          *
*       G       : gradient of the SCF energy with respect to the       *
*                 rotation parameters of length nG                     *
*                                                                      *
*     called from: GrdClc                                              *
*                                                                      *
*     calls to: ModFck, Asym                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*                                                                      *
*     Corrected for an error of 2.0                                    *
*     R. Lindh, Harvard University, 2017                               *
*                                                                      *
************************************************************************
      Use Orb_Type
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 O(nOTSD),T(nOTSD,nD),S(nOTSD),D(nOTSD,nD),C(nC,nD),
     &       V(nOTSD,nD), G(nG,nD), CMO(nC,nD)
*
      Real*8, Dimension(:,:), Allocatable:: FckM
      Real*8, Dimension(:), Allocatable:: Aux1, Aux2, Aux3
*
*----------------------------------------------------------------------*
*     Start
*define _DEBUG_
#ifdef _DEBUG_
      Call NrmClc(O,nOTSD   ,'EGrad','O')
      Call NrmClc(S,nOTSD   ,'EGrad','S')
      Call NrmClc(D,nOTSD*nD,'EGrad','D')
      Call NrmClc(T,nOTSD*nD,'EGrad','T')
      Call NrmClc(V,nOTSD*nD,'EGrad','V')
      Call NrmClc(C,nC   *nD,'EGrad','C')
#endif
*----------------------------------------------------------------------*
*
*---- Allocate memory for modified fock matrix
      Call mma_allocate(FckM,nBT,nD,Label='FckM')
      Call FZero(FckM,nBT*nD)
      Call FZero(G,nG*nD)
*
*---- Allocate memory for auxiliary matrices
      Call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
      Call mma_allocate(Aux2,MaxBas**2,Label='Aux2')
      Call mma_allocate(Aux3,MaxBas**2,Label='Aux3')
*
      Do iD = 1, nD
*
         Call DZAXPY(nBT,1.0D0,O,1,T(1,iD),1,FckM(1,iD),1)
#ifdef _DEBUG_
         Call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
#endif
         If (nnFr.gt.0)
     &      Call ModFck(FckM(1,iD),S,nBT,CMO(1,iD),nBO,nOcc(1,1))
*
         Call DaXpY_(nBT,1.0D0,V(1,iD),1,FckM(1,iD),1)
#ifdef _DEBUG_
         Call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
#endif
*
         iOff = 0
         ij = 1
         it = 1
         ig = 1
         Do iSym = 1, nSym
            nBs = nBas(iSym)
            nOr = nOrb(iSym)
            nOrbmF = nOrb(iSym)-nFro(iSym)
*
            If (nOrb(iSym).gt.0) Then
               lth = nBs*(nBs + 1)/2
*define _ALTERNATIVE_CODE_
#ifdef _ALTERNATIVE_CODE_
*
*              This is an alternative section to compute the gradients.
*              In this section you have an early difference followed
*              by a late purification to guarantee G is stricktly
*              anti-symmetric.
*
*              Observations so far is only minute changes on very small
*              gradients.
*
*----------    Square Fock matrix
               Call Square(FckM(ij,iD),Aux1,1,nBs,nBs)
*----------    Square density matrix
               Call DSq(D(ij,iD),Aux2,1,nBs,nBs)
*----------    Perform FD
               Call DGEMM_('N','N',
     &                     nBs,nBs,nBs,
     &                     1.0D0,Aux1,nBs,
     &                           Aux2,nBs,
     &                     0.0D0,Aux3,nBs)
*----------    Square overlap matrix and perform FDS
               Call Square(S(ij),Aux1,1,nBs,nBs)
               Call DGEMM_('N','N',
     &                     nBs,nBs,nBs,
     &                     1.0D0,Aux3,nBs,
     &                           Aux1,nBs,
     &                     0.0D0,Aux2,nBs)
*----------    Form FDS-SDF
               Call Asym(Aux2,Aux3,nBs)
*----------    Perform C(T)(FDS-SDF)
               Call DGEMM_('T','N',
     &                     nOr,nBs,nBs,
     &                     1.0d0,C(it,iD),nBs,
     &                           Aux3,nBs,
     &                     0.0d0,Aux1,nOr)
*----------    Perform C(T)(FDS-SDF)C
               Call DGEMM_('N','N',
     &                     nOr,nOr,nBs,
     &                     1.0d0,Aux1,nOr,
     &                           C(it,iD),nBs,
     &                     0.0d0,G(ig,iD),nOr)
               Call Purify(G(ig,iD),nOr)
#else
*
*----------    Square Fock matrix and perform C(T)F
               Call Square(FckM(ij,iD),Aux2,1,nBs,nBs)
               Call DGEMM_('T','N',
     &                     nOr,nBs,nBs,
     &                     1.0d0,C(it,iD),nBs,
     &                           Aux2,nBs,
     &                     0.0d0,Aux1,nOr)
*
*----------    Square density matrix and perform C(T)FD
               Call DSq(D(ij,iD),Aux2,1,nBs,nBs)
               Call DGEMM_('N','N',
     &                     nOr,nBs,nBs,
     &                     1.0d0,Aux1,nOr,
     &                           Aux2,nBs,
     &                     0.0d0,Aux3,nOr)
*
*----------    Square overlap matrix and perform C(T)FDS
               Call Square(S(ij),Aux2,1,nBs,nBs)
               Call DGEMM_('N','N',
     &                     nOr,nBs,nBs,
     &                     1.0d0,Aux3,nOr,
     &                           Aux2,nBs,
     &                     0.0d0,Aux1,nOr)
*----------    C(T)FDSC
               Call DGEMM_('N','N',
     &                     nOr,nOr,nBs,
     &                     1.0d0,Aux1,nOr,
     &                           C(it,iD),nBs,
     &                     0.0d0,Aux2,nOr)
*
               Call Asym(Aux2,G(ig,iD),nOr)
#endif
*
*              At this point enforce that the gradient is exactly zero
*              for elements corresponding to orbitals of different
*              fermion types.
*
               Do i = 1, nOr
                  If (i.le.nFro(iSym)) Then
                     k=-1
                  Else
                     k=OrbType(iOff+i-nFro(iSym),iD)
                  End If
*
                  Do j = 1, nOr
                     If (i.le.nFro(iSym)) Then
                        l=-1
                     Else
                        l=OrbType(iOff+j-nFro(iSym),iD)
                     End If
*
                     ih = ig + (i-1)*nOr + j - 1
                     If (k.lt.0 .or. l.lt.0 .or. k.ne.l) G(ih,iD)=0.0D0
*
                  End Do
               End Do
*
            End If
            ij = ij + nBs*(nBs + 1)/2
            it = it + nBs*nOr
            ig = ig + nOr*nOr
            iOff = iOff + nOrbmF
*
         End Do ! iSym
*
      End Do ! iD
*
*---- Deallocate memory
      Call mma_deallocate(Aux3)
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(FckM)
*
      Call DScal_(nG*nD,2.0D0,G,1)
*
#ifdef _DEBUG_
*     Call RecPrt('EGrad: G',' ',G,nG,nD)
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
      SubRoutine Asym(H,A,n)
************************************************************************
*                                                                      *
*     purpose: Antisymmetrize matrix                                   *
*                                                                      *
*     input:                                                           *
*       H       : input matrix                                         *
*       n       : dimension                                            *
*                                                                      *
*     output:                                                          *
*       A       : antisymmetrized matrix                               *
*                                                                      *
*     called from: EGrad                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 H(n,n),A(n,n)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
      Do j = 1, n
         Do i = 1, n
            A(i,j) = H(i,j) - H(j,i)
         End Do
      End Do
*
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
#ifdef _NEW_CODE_
      Subroutine Purify(A,n)
      Implicit Real*8 (a-h,o-z)
      Real*8 A(n,n)
*
      Do i = 1, n
         Do j = 1, n
            tmp = A(i,j) - A(j,i)
            A(i,j) = 0.5D0*tmp
            A(j,i) = -A(i,j)
         End Do
      End Do
      Return
      End
#endif
