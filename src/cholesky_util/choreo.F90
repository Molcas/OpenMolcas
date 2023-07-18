!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! Data for Cholesky vector reordering.
!
Module ChoReO
      INTEGER NNBST(8)
      INTEGER NABPK(8,8)
      INTEGER LUFV(8,8)
      CHARACTER(LEN=4), PARAMETER:: REONAM = 'CHFV'
End Module ChoReO
