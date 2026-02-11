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
      SUBROUTINE ORTHO2(S,U,V,N)
!
!     Purpose: normalize vector U and calculate V=S*U.
!
!     Called from: ORTHO1.
!
!          ****** IBM 3090 MOLCAS Release: 90 02 22 ******
!
      use output_ras, only: LF
      IMPLICIT None
      Integer N
#include "warnings.h"
      REAL*8 S(*),U(*),V(*)

      REAL*8 THR,SUM,X
      REAl*8, External:: DDot_
      INTEGER I

      THR=1.D-10
      IF ( N.EQ.0 ) Return
      CALL DGEMM_('N','N',                                              &
     &             N,1,N,                                               &
     &             1.0d0,S,N,                                           &
     &             U,N,                                                 &
     &             0.0d0,V,N)
      SUM=DDOT_(N,U,1,V,1)
      IF (SUM.LT.THR) THEN
        Write(LF,*)' TEST IN ORTHO2: N=',N
        Write(LF,'(1X,5G16.8)') (U(I),I=1,N)
        Write(LF,'(1X,5G16.8)') (V(I),I=1,N)
        Write(LF,*)' Error in ORTHO2. Norm=', SUM
        Write(LF,*)' RASSCF tried to orthonormalize orbitals, but'
        Write(LF,*)' failed due to a condition that should not be'
        Write(LF,*)' possible in a low-level subroutine. Either'
        Write(LF,*)' some extremely strange orbitals have been'
        Write(LF,*)' produced, or something is seriously wrong'
        Write(LF,*)' with the program. Please check, and consider'
        Write(LF,*)' issuing a bug report.'
        CALL QUIT(_RC_GENERAL_ERROR_)
      ENDIF
      X=1.0D0/SQRT(SUM)
      DO I=1,N
        U(I)=X*U(I)
        V(I)=X*V(I)
      END DO
      END SUBROUTINE ORTHO2
