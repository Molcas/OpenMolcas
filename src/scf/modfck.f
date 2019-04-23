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
************************************************************************
      SubRoutine ModFck(Fock,Ovlp,nFO,CMO,nCMO,mynOcc)
************************************************************************
*                                                                      *
*     purpose: Modify Fock matrix taking into account frozen orbitals  *
*              F <- [1 - S*D(f)]*F  (symmetrized).                     *
*                                                                      *
*     input:                                                           *
*       Fock    : Fock matrix of length nFO                            *
*       Ovlp    : overlap matrix of length nFO                         *
*                                                                      *
*     output:                                                          *
*       Fock    : modified Fock matrix                                 *
*                                                                      *
*     called from: NewOrb, DCore, EGrad                                *
*                                                                      *
*     calls to: DFroz, Sym                                             *
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 Fock(nFO),Ovlp(nFO),CMO(nCMO)
      Integer mynOcc(*)
      Real*8, Dimension(:), Allocatable:: DFro, DFSq, OvSq, Aux1
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Allocate memory for density contribution
      Call mma_allocate(DFro,nBT,Label='DFro')
*
*---- Allocate memory for squared density (then for squared Fock matrix)
      Call mma_allocate(DFSq,MaxBas**2,Label='DFSq')
*
*---- Allocate memory for squared overlap and later for 1 - S*D(f)
      Call mma_allocate(OvSq,MaxBas**2,Label='OvSq')
*
*---- Allocate memory for S*D(f)
      Call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
*
*---- Generate contribution to the density mat. from frozen orbitals
      Call DFroz(DFro,nBT,CMO,nBO,mynOcc)
*
      ij   = 1
      Do iSym = 1, nSym
*
         If (nBas(iSym).gt.0) Then
*---------- Square density and overlap
            Call DSq(DFro(ij),DFSq,1,nBas(iSym),nBas(iSym))
            Call Square(Ovlp(ij),OvSq,1,nBas(iSym),nBas(iSym))
*
*---------- Compute S*D(f)
            Call DGEMM_('N','N',
     &                  nBas(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,OvSq,nBas(iSym),
     &                        DFSq,nBas(iSym),
     &                  0.0d0,Aux1,nBas(iSym))
*
*---------- Put unit matrix to memory ascribed for overlap and compute 1 - S*D(f)
            call dcopy_(nBas(iSym)**2,[Zero],0,OvSq,1)
            call dcopy_(nBas(iSym)   ,[One] ,0,OvSq,nBas(iSym) + 1)
            call daxpy_(nBas(iSym)**2,-One,Aux1,1,OvSq,1)
*
*---------- Square Fock matrix and store in memory ascribed for density
            Call Square(Fock(ij),DFSq,1,nBas(iSym),nBas(iSym))
*
*---------- Compute F <- [1 - S*D(f)]*F
            Call DGEMM_('N','N',
     &                  nBas(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,OvSq,nBas(iSym),
     &                        DFSq,nBas(iSym),
     &                  0.0d0,Aux1,nBas(iSym))
*
*---------- Symmetrize modified Fock matrix
            Call Sym(Aux1,Fock(ij),nBas(iSym))
         End If
*
*------- Update pointers
         ij = ij + nBas(iSym)*(nBas(iSym) + 1)/2
*
      End Do
*
*---- Deallocate memory
      Call mma_deallocate(Aux1)
      Call mma_deallocate(OvSq)
      Call mma_deallocate(DFSq)
      Call mma_deallocate(DFro)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
      SubRoutine Sym(A,S,n)
************************************************************************
*                                                                      *
*     purpose: Symmetrize matrix                                       *
*                                                                      *
*     input:                                                           *
*       A       : input matrix (square)                                *
*       n       : dimension                                            *
*                                                                      *
*     output:                                                          *
*       S       : Symmetrized matrix (triangular)                      *
*                                                                      *
*     called from: ModFck                                              *
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
      Real*8 A(n,n),S(n*(n+1)/2)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      ij = 0
      Do i = 1, n
         Do j = 1, i
            ij = ij + 1
            S(ij) = (A(i,j) + A(j,i))/2.0d+00
         End Do
      End Do
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
