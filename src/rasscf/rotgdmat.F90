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
      Subroutine RotGDMat(R,GD)
      use rasscf_global, only: lRoots, NAC
      Implicit None


#include "warnings.h"

      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD,GD2
      Real*8,DIMENSION(lRoots,lRoots)::R

      INTEGER I,J,K,L,p,q,iIJ,iKL,ip,iq

      DO p=1,nac
      DO q=1,nac
       Do I=1,lRoots
       Do J=1,I
        iIJ=(I-1)*I/2+J
        GD2(iIJ,p,q)=0.0d0
        do K=1,lRoots
        do L=1,lRoots
         IF(K.gt.L) THEN
          iKL=(K-1)*K/2+L
          ip=p
          iq=q
         ELSE
          iKL=(L-1)*L/2+K
          ip=q
          iq=p
         END IF
         GD2(iIJ,p,q)=GD2(iIJ,p,q)+GD(iKL,ip,iq)*R(I,K)*R(J,L)
        end do
        end do
       End Do
       End Do
      END DO
      END DO

      DO p=1,nac
      DO q=1,nac
       Do I=1,lRoots
       Do J=1,I
        iIJ=(I-1)*I/2+J
        GD(iIJ,p,q)=GD2(iIJ,p,q)
       End Do
       End Do
      END DO
      END DO
      END SUBROUTINE RotGDMat
