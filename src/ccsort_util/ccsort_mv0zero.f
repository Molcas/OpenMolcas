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
       SUBROUTINE ccsort_mv0zero                                        &
     & (DD,LENGTH,MAT)
!
       INTEGER           DD
       INTEGER           LENGTH
       real*8  MAT(1:DD)
       INTEGER           INIT
       real*8  ZERO
!
       DATA              ZERO/0.0D+00/
!
!     ...loop over all elements
!
       DO 10 INIT=1,LENGTH
       MAT(INIT) = ZERO
 10     CONTINUE
!
       RETURN
       END
