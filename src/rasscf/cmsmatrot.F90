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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine CMSMatRot(Mat,A,I,J,N)
      Implicit None
      INTEGER I,J,N
      Real*8 A
      Real*8,DIMENSION(N,N)::Mat,TM
      INTEGER K
      DO K=1,N
       TM(I,K)=Mat(I,K)
       TM(J,K)=Mat(J,K)
      END DO
      DO K=1,N
       Mat(J,K)= cos(A)*TM(J,K)+sin(A)*TM(I,K)
       Mat(I,K)=-sin(A)*TM(J,K)+cos(A)*TM(I,K)
      END DO
      END SUBROUTINE CMSMatRot
