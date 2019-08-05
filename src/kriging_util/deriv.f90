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
subroutine deriv(finder,foutder,nd,d1,d2)
    use globvar
#include "stdalloc.fh"
    integer d1,d2,nd
    real*8 finder(d1,d2),foutder(d1,d2),a,kr,nr
    real*8, Allocatable :: b(:,:)
!
    Call mma_allocate(b,d1,d2,label="b")
!
    nr=DBLE(nd)
    a=Gamma(nr+1.0d0)/h**nd
    b=0.0
    do k=0,nd
        kr=DBLE(k)
        finder=finder+kr*h
        b=b+DBLE((-1)**(k+1))/(Gamma(nr-kr+1.0D0)*Gamma(kr+1.0D0))*finder
    enddo
    foutder=a*b*DBLE((-1)**(nr+1))
!
    Call mma_deallocate(b)
!
end
