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
! Copyright (C) 2000, Markus P. Fuelscher                              *
!***********************************************************************
      Subroutine ClnMO(CMO)

!***********************************************************************
!                                                                      *
!     In order to preserve symmetry of orbitals which belong           *
!     to a point group of higher order than those allowed in           *
!     MOLCAS, it is sometimes necessary to enforce the MO coefficients *
!     to remain zero. This subroutine does that job in every           *
!     iteration.                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     This subroutine replaces the one written earlier by L. Serrano   *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!***********************************************************************
      use Constants, only: Zero
      use general_data, only: CleanMask, NSYM, NBAS

      Implicit None
      Real*8 CMO(*)
      Integer ij, iSym, mBas, i, j

! Prelude


! Body

      ij = 0
      Do iSym = 1,nSym
        mBas = nBas(iSym)
        Do i = 1,mBas
          Do j = 1,mBas
            ij = ij+1
            If ( CleanMask(ij).eq.1 ) CMO(ij) = Zero
          End Do
        End Do
      End Do

! Epilogue

      End Subroutine ClnMO
