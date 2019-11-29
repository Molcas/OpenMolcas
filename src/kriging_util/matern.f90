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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
SUBROUTINE matern(dh, m, d1, d2)
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
  integer d1,d2,i
  REAL*8 a,d,dh(d1,d2),m(d1,d2)
  REAL*8, Allocatable :: d0(:,:)
  INTEGER*8 c
!
  Call mma_Allocate(d0,d1,d2,label="d0")
!
! For this expresion you can check https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
! and equations (11) and (12) on ref.
  d0(:,:) = sqrt(dh)
  c = idint(pAI)
  a = Gamma(pAI+1)/Gamma(2*pAI+1)
  m = 0.0D0
  do i=0, c
    d=DBLE(i)
    m = m + (Gamma(pAI+1.0D0+d)/(Gamma(d+1.D0)*Gamma(pAI+1.0D0-d)))*(2.0D0*Sqrt(2.0D0*pAI+1.0D0)*d0)**(pAI-i)
  enddo
  m = a*m*exp(-sqrt(2.0D0*pAI+1)*d0)
!
  Call mma_deallocate(d0)
!
END SUBROUTINE matern
