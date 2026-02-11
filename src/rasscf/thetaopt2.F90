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
      Subroutine ThetaOpt2(R,theta,deltaQ,SPair,NP,GD,Vee,G)
      use rasscf_global, only: lRoots, NAC
      Implicit None


#include "warnings.h"
      INTEGER NP
      Real*8,DIMENSION(NP)::theta
      Real*8 Change,deltaQ
      INTEGER,DIMENSION(NP,2)::SPair
      Real*8,DIMENSION(lroots,lroots)::R
      Real*8,DIMENSION(lroots)::Vee
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G

      INTEGER IP,I,J
      deltaQ=0.0d0
      DO IP=1,NP
       I=SPair(IP,1)
       J=SPair(IP,2)
       CALL OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
       deltaQ=deltaQ+change
      END DO

      DO IP=NP-1,1,-1
       I=SPair(IP,1)
       J=SPair(IP,2)
       CALL OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
       deltaQ=deltaQ+change
      END DO

      END SUBROUTINE ThetaOpt2
