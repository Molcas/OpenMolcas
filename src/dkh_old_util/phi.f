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
! Copyright (C) 1984,1986, Bernd Artur Hess                            *
!***********************************************************************
      REAL*8 FUNCTION PHI(M,N)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "crelop.fh"
      IF (MOD(N,2).EQ.1.OR.MOD(M,2).EQ.1) GOTO 10
      PHI=2.D0*GA(M+1)*GA(N+1)/GA(M+N+2)
      RETURN
10    PHI=0.D0
      RETURN
      END
