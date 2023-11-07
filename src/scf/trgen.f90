!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine TrGen(TrMat,nTrMat,Ovrlp,OneHam,mBT)
!***********************************************************************
!                                                                      *
!     purpose: Generate transformation matrix from AO's in which       *
!              integrals are computed to orthogonal, symmetry adapted  *
!              (and spherical, if desired) functions - near linear     *
!              dependencies are also removed.                          *
!                                                                      *
!     input:                                                           *
!       Ovrlp    : overlap matrix in AO basis of length mBT            *
!                                                                      *
!     output:                                                          *
!       TrMat   : Transformation matrix of length nTrMat               *
!                                                                      *
!***********************************************************************
!
      use InfSCF, only: DelThr, nBO, nBT, nnFr, nSym, nBas
      use Constants, only: Zero, One
      Implicit None
!
      Integer nTrMat, mBT
      Real*8 TrMat(nTrMat),Ovrlp(mBT),OneHam(mBT)

      Integer i, j, iSym, ind
!
      ind=0
      Do iSym = 1, nSym
         Do i = 1, nBas(iSym)
            Do j = 1, nBas(iSym)
               ind = ind + 1
               TrMat(ind) = Zero
               If (I.eq.J) TrMat(ind) = One
            End Do
         End Do
      End Do
!
!---- Set up certain parameters (nOrb(i) may be changed)
      Call SetUp_SCF()
!
!---- Move frozen atomic orbitals to the begining
      If (nnFr.gt.0) Then
         Call Freeze(TrMat,nBO,OneHam,mBT)
         Call SetUp_SCF()
      End If
!
!---- Remove near linear dependencies from basis set
      If (DelThr.ne.Zero) Then
         Call OvlDel(Ovrlp,nBT,TrMat,nBO)
         Call SetUp_SCF()
      End If
!
!---- Orthogonalize final orbitals
      Call Ortho(TrMat,nBO,Ovrlp,nBT)
!
      End subroutine TrGen
