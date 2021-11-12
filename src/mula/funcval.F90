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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine funcval(x,coef,ipow,fval,nterm,nvar)
!  Purpose:
!    Return function value at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One

real*8 x(nvar)
real*8 coef(nterm)
integer ipow(nterm,nvar)
real*8 sum, prod, fval

sum = Zero
do iterm=1,nterm
  prod = One
  do ivar=1,nvar
    nsum = ipow(iterm,ivar)
    prod = prod*(x(ivar)**(nsum))
  end do
  sum = sum+coef(iterm)*prod
end do
fval = sum

end subroutine funcval
!####
subroutine gradient(x,coef,ipow,grad,nterm,nvar)
!  Purpose:
!    Return gradient at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One

real*8 x(nvar)
real*8 coef(nvar)
integer ipow(nterm,nvar)
real*8 grad(nvar)
real*8 sum, prod
logical ijEq
real*8 rfactor

do ivar=1,nvar
  sum = Zero
  do iterm=1,nterm
    prod = One
    do jvar=1,nvar
      ijEq = (ivar == jvar)
      if ((ipow(iterm,jvar) >= 1) .and. ijEq) then
        nder = 1
        nsum = ipow(iterm,jvar)-nder
      else if ((ipow(iterm,jvar) >= 0) .and. (.not. ijEq)) then
        nder = 0
        nsum = ipow(iterm,jvar)
      else
        nder = -1
        nsum = 0
      end if
      call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(jvar)**(nsum))
    end do
    sum = sum+coef(iterm)*prod
  end do
  grad(ivar) = sum
end do

end subroutine gradient
!####
subroutine Hessian(x,coef,ipow,Hess,nterm,nvar)
!  Purpose:
!    Return Hessian at x.
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One

real*8 x(nvar)
real*8 coef(nvar)
integer ipow(nterm,nvar)
real*8 Hess(nvar,nvar)
real*8 sum, prod
logical ijEq, ikEq, jkEq
real*8 rfactor

do ivar=1,nvar
  do jvar=ivar,nvar
    ijEq = (ivar == jvar)
    sum = Zero
    do iterm=1,nterm
      prod = One
      do kvar=1,nvar
        ikEq = (ivar == kvar)
        jkEq = (jvar == kvar)
        if ((ipow(iterm,kvar) >= 2) .and. ikEq .and. jkEq) then
          nder = 2
          nsum = ipow(iterm,kvar)-nder
        else if ((ipow(iterm,kvar) >= 1) .and. (ikEq .or. jkEq) .and. (.not. ijEq)) then
          nder = 1
          nsum = ipow(iterm,kvar)-nder
        else if ((ipow(iterm,kvar) >= 0) .and. (.not. ikEq) .and. (.not. jkEq)) then
          nder = 0
          nsum = ipow(iterm,kvar)
        else
          nder = -1
          nsum = 0
        end if
        call factor(nsum,nder,rfactor)
        prod = prod*rfactor*(x(kvar)**(nsum))
      end do
      sum = sum+coef(iterm)*prod
    end do
    Hess(ivar,jvar) = sum
    Hess(jvar,ivar) = sum
  end do
end do

end subroutine Hessian
!####
subroutine thirdDer(x,coef,ipow,D3,nterm,nvar)
!  Purpose:
!    Return third derivatives at x.
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One

real*8 x(nvar)
real*8 coef(nvar)
integer ipow(nterm,nvar)
real*8 D3(nvar,nvar,nvar)
real*8 sum, prod
logical ilEq, jlEq, klEq
real*8 rfactor

D3(:,:,:) = Zero
do ivar=1,nvar
  do jvar=ivar,nvar
    do kvar=jvar,nvar
      sum = Zero
      do iterm=1,nterm
        prod = One
        do lvar=1,nvar
          ilEq = (ivar == lvar)
          jlEq = (jvar == lvar)
          klEq = (kvar == lvar)
          if ((ipow(iterm,lvar) >= 3) .and. ilEq .and. jlEq .and. klEq) then
            nder = 3
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 2) .and. ((ilEq .and. jlEq .and. (.not. klEq)) .or. &
                                                  (jlEq .and. klEq .and. (.not. ilEq)))) then
            nder = 2
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 1) .and. ((ilEq .and. (.not. jlEq) .and. (.not. klEq)) .or. &
                                                  (jlEq .and. (.not. ilEq) .and. (.not. klEq)) .or. &
                                                  (klEq .and. (.not. ilEq) .and. (.not. jlEq)))) then
            nder = 1
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 0) .and. (.not. ilEq) .and. (.not. jlEq) .and. (.not. klEq)) then
            nder = 0
            nsum = ipow(iterm,lvar)
          else
            nder = -1
            nsum = 0
          end if
          call factor(nsum,nder,rfactor)
          prod = prod*rfactor*(x(lvar)**(nsum))
        end do
        sum = sum+coef(iterm)*prod
      end do
      D3(ivar,jvar,kvar) = sum
      D3(kvar,ivar,jvar) = sum
      D3(jvar,kvar,ivar) = sum
      D3(jvar,ivar,kvar) = sum
      D3(kvar,jvar,ivar) = sum
      D3(ivar,kvar,jvar) = sum
    end do
  end do
end do

end subroutine thirdDer
!####
subroutine fourthDer(x,coef,ipow,D4,nterm,nvar)
!  Purpose:
!    Return fourth derivatives at x.
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One

real*8 x(nvar)
real*8 coef(nvar)
integer ipow(nterm,nvar)
real*8 D4(nvar,nvar,nvar,nvar)
real*8 sum, prod
logical imEq, jmEq, kmEq, lmEq
real*8 rfactor

D4(:,:,:,:) = Zero
do ivar=1,nvar
  do jvar=ivar,nvar
    do kvar=jvar,nvar
      do lvar=kvar,nvar
        sum = Zero
        do iterm=1,nterm
          prod = One
          do mvar=1,nvar
            imEq = (ivar == mvar)
            jmEq = (jvar == mvar)
            kmEq = (kvar == mvar)
            lmEq = (lvar == mvar)
            if ((ipow(iterm,mvar) >= 4) .and. imEq .and. jmEq .and. kmEq .and. lmEq) then
              nder = 4
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 3) .and. ((imEq .and. jmEq .and. kmEq .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. kmEq .and. lmEq .and. (.not. imEq)))) then
              nder = 3
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 2) .and. ((imEq .and. jmEq .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. kmEq .and. (.not. imEq) .and. (.not. lmEq)) .or. &
                                                    (kmEq .and. lmEq .and. (.not. imEq) .and. (.not. jmEq)))) then
              nder = 2
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 1) .and. ((imEq .and. (.not. jmEq) .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. (.not. imEq) .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (kmEq .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. lmEq)) .or. &
                                                    (lmEq .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. kmEq)))) then
              nder = 1
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 0) .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. kmEq) .and. (.not. lmEq)) then
              nder = 0
              nsum = ipow(iterm,mvar)
            else
              nder = -1
              nsum = 0
            end if
            call factor(nsum,nder,rfactor)
            prod = prod*rfactor*(x(mvar)**(nsum))
          end do
          sum = sum+coef(iterm)*prod
        end do
        D4(ivar,jvar,kvar,lvar) = sum
        D4(ivar,kvar,jvar,lvar) = sum
        D4(ivar,jvar,lvar,kvar) = sum
        D4(ivar,lvar,kvar,jvar) = sum
        D4(ivar,lvar,jvar,kvar) = sum
        D4(ivar,kvar,lvar,jvar) = sum

        D4(jvar,ivar,kvar,lvar) = sum
        D4(jvar,kvar,ivar,lvar) = sum
        D4(jvar,ivar,lvar,kvar) = sum
        D4(jvar,lvar,kvar,ivar) = sum
        D4(jvar,lvar,ivar,kvar) = sum
        D4(jvar,kvar,lvar,ivar) = sum

        D4(kvar,ivar,jvar,lvar) = sum
        D4(kvar,jvar,ivar,lvar) = sum
        D4(kvar,ivar,lvar,jvar) = sum
        D4(kvar,lvar,jvar,ivar) = sum
        D4(kvar,lvar,ivar,jvar) = sum
        D4(kvar,jvar,lvar,ivar) = sum

        D4(lvar,ivar,jvar,kvar) = sum
        D4(lvar,jvar,ivar,kvar) = sum
        D4(lvar,kvar,jvar,ivar) = sum
        D4(lvar,ivar,kvar,jvar) = sum
        D4(lvar,kvar,ivar,jvar) = sum
        D4(lvar,jvar,kvar,ivar) = sum
      end do
    end do
  end do
end do

end subroutine fourthDer
