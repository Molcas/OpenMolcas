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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
      Subroutine Get_Ntri(nTri)
      use input_mclr, only: nSym,nBas
      Implicit None

      INTEGER nTri,kSym
      nTri=0
      DO kSym=1,nSym
       nTri=nTri+nBas(kSym)*(nBas(kSym)+1)/2
      END DO
      END Subroutine Get_Ntri
