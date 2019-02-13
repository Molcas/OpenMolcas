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
subroutine deriv(finder,foutder,nd,d1,d2)
    use globvar
    integer d1,d2,nd
    real*8 finder(d1,d2),foutder(d1,d2),a,b(d1,d2),kr,nr
    nr=real(nd)
    a=Gamma(nr+1.0)/h**nd
    b=0.0
    do k=0,nd
        kr=real(k)
        finder=finder+kr*h
        b=b+(-1)**(k+1)/(Gamma(nr-kr+1.0)*Gamma(kr+1.0))*finder
    enddo
    foutder=a*b*(-1)**(nr+1)
end
