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

subroutine qpasy(alpha,dfac,npi,l,lambu,lmahi,lmbhi,ltot1,xka,xkb,prd,dum,qsum)
! compute type 2 radial integrals, scaled by exp(-arc2)/sqrt(pi),
! using the partially asymptotic method,
! for lama=l to lmahi, lamb=l to lmbhi, n=lama+lamb-l-l

use Constants, only: Zero, One, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: npi, l, lambu, lmahi, lmbhi, ltot1
real(kind=wp), intent(in) :: alpha, dfac(*), xka, xkb, prd, dum
real(kind=wp), intent(inout) :: qsum(ltot1,lambu,lmahi)
integer(kind=iwp) :: j, la, lama, lamb, lb, n, nprime
real(kind=wp) :: alf1, f1, f2, prde, prefac, qnew, qold1, qold2, sqalpi, dsum, t, tk, xk
real(kind=wp), parameter :: accrcy = 1.0e-13_wp
real(kind=wp), external :: qcomp

sqalpi = One/sqrt(alpha)
alf1 = One

if (xka <= xkb) then

  ! xka is smaller: set up parameters for qcomp using xkb

  xk = xkb*sqalpi
  t = Quart*xk*xk
  prde = prd*exp(t-dum)*sqalpi**(npi+l)
  if (l >= 2) prde = prde*xka**(l-1)
  tk = xka*xka/(alpha+alpha)
  do lama=l,lmahi
    la = lama-1
    prefac = prde
    do lamb=l,lmbhi
      lb = lamb-1
      n = ((1-l-l)+lama)+lamb
      ! run power series using xka, obtaining initial
      ! q(n,l) values from qcomp, then recurring upwards
      ! j=0 term in sum
      nprime = npi+n+la-1
      qold1 = qcomp(alf1,dfac,nprime,lb,t,xk)/dfac(la+la+3)
      dsum = qold1
      if (tk /= Zero) then
        ! j=1 term in sum
        nprime = nprime+2
        qnew = qcomp(alf1,dfac,nprime,lb,t,xk)/dfac(la+la+3)
        f1 = (la+la+3)
        qold2 = (tk/f1)*qold1
        qold1 = (tk/f1)*qnew
        dsum = dsum+qold1
        j = 1
        do
          ! increment j for next term
          j = j+1
          nprime = nprime+2
          f1 = (nprime+nprime-5)
          f2 = ((lb-nprime+4)*(lb+nprime-3))
          qnew = (t+Half*f1)*qold1+Quart*f2*qold2
          f1 = (j*(la+la+j+j+1))
          qold2 = (tk/f1)*qold1
          qold1 = (tk/f1)*qnew
          dsum = dsum+qold1
          if (qold1 <= accrcy*dsum) exit
        end do
      end if
      qsum(n,lamb,lama) = qsum(n,lamb,lama)+prefac*dsum
      prefac = prefac*sqalpi
    end do
    prde = prde*(xka/alpha)
  end do

else

  ! xkb is smaller: set up parameters for qcomp using xka

  xk = xka*sqalpi
  t = Quart*xk*xk
  prde = prd*exp(t-dum)*sqalpi**(npi+l)
  if (l >= 2) prde = prde*xkb**(l-1)
  tk = xkb*xkb/(alpha+alpha)
  do lama=l,lmahi
    la = lama-1
    prefac = prde
    do lamb=l,lmbhi
      lb = lamb-1
      n = ((1-l-l)+lama)+lamb
      ! run power series using xkb, obtaining initial
      ! q(n,l) values from qcomp, then recurring upwards
      ! j=0 term in sum
      nprime = npi+n+lb-1
      qold1 = qcomp(alf1,dfac,nprime,la,t,xk)/dfac(lb+lb+3)
      dsum = qold1
      if (tk /= Zero) then
        ! j=1 term in sum
        nprime = nprime+2
        qnew = qcomp(alf1,dfac,nprime,la,t,xk)/dfac(lb+lb+3)
        f1 = (lb+lb+3)
        qold2 = (tk/f1)*qold1
        qold1 = (tk/f1)*qnew
        dsum = dsum+qold1
        j = 1
        do
          ! increment j for next term
          j = j+1
          nprime = nprime+2
          f1 = (nprime+nprime-5)
          f2 = ((la-nprime+4)*(la+nprime-3))
          qnew = (t+Half*f1)*qold1+Quart*f2*qold2
          f1 = (j*(lb+lb+j+j+1))
          qold2 = (tk/f1)*qold1
          qold1 = (tk/f1)*qnew
          dsum = dsum+qold1
          if (qold1 <= accrcy*dsum) exit
        end do
      end if
      qsum(n,lamb,lama) = qsum(n,lamb,lama)+prefac*dsum
      prefac = prefac*(xkb/alpha)
    end do
    prde = prde*sqalpi
  end do

end if

return

end subroutine qpasy
