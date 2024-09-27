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

subroutine RRINT(K,ALFA,A,BETA,R0,GRINT,lmax)

use welcom, only: binom, fac, fiint, kmax
use Constants, only: Zero, One, Two, Three, Four, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K, lMax
real(kind=wp), intent(in) :: Alfa, A, Beta, R0
real(kind=wp), intent(out) :: grint(0:lmax,lmax)
integer(kind=iwp) :: i, Ind, kk, l, ll, M, mm, mMax, n
real(kind=wp) :: AA, AA2, AA3, AA4, AA5, AExp, Al, bExp, bExp1, bExp2, Bi, Exp1, ExpA, FiIntM, ggg, Pi4, rri(0:kmax+2), Test, Tmp, &
                 Tmp1, Tmp2, Tmp3, Tmp4
real(kind=wp), external :: QRint

M = K+1
EXPA = -A*A*ALFA-BETA*R0*R0
AEXP = ALFA+BETA
BEXP1 = -Two*(BETA*R0+A*ALFA)/AEXP
BEXP2 = -Two*(BETA*R0-A*ALFA)/AEXP
if (A == Zero) then
  Test = Zero
else
  if (R0 == Zero) then
    Test = One
  else
    TEST = A*ALFA
  end if
end if

if (TEST >= 0.005_wp) then

  !write(u6,*) ' Large A'
  ! K=0 ONE CONTRIBUTION SS-INTEGRAL
  do i=0,k
    rri(i) = qrint(i+1,aexp,bexp1,expa)*real((-1)**i,kind=wp)-qrint(i+1,aexp,bexp2,expa)
  end do
  !call RecPrt(' In RRInt: rri',' ',rri,k+1,1)
  AL = One/(Two*ALFA*A)
  do i=0,k
    mmax = i/2
    do m=1,mmax+1
      ! calculate integral ri(i,m)
      fiintm = fiint(m-1,0)
      grint(i,m) = Zero
      do n=1,m
        Bi = binom(m-1,n-1)*real((-1)**(n+1),kind=wp)
        ind = i-(m-n)*2
        do kk=0,ind
          ggg = fac(ind)/fac(ind-kk)*al**(kk+1)
          grint(i,m) = grint(i,m)+bi*ggg*rri(i-kk)*fiintm
        end do
      end do
    end do
  end do

else

  !write(u6,*) ' SERIES EXPANSION FOR SMALL A'

  ! SERIES EXPANSION FOR SMALL A

  ! K=0 FIRST
  l = (k+1)/2
  EXP1 = -(ALFA*A*A+BETA*R0*R0)
  BEXP = -Two*BETA*R0/(ALFA+BETA)
  do i=0,l+2
    rri(i) = qrint(2*(i+1),AExp,BExp,Exp1)
  end do
  !call RecPrt(' rri',' ',rri,l+3,1)
  pi4 = pi*Four
  AA = Two*(A*Alfa)
  AA2 = Two*(A*Alfa)**2
  AA3 = Four*(A*Alfa)**3
  AA4 = Two*(A*Alfa)**4
  AA5 = Four*(A*Alfa)**5
  GRINT(0,1) = pi4*(rri(0)+AA2/Three*rri(1)+AA4/15.0_wp*rri(2))
  if (K == 0) return
  do ll=1,l
    do kk=1,ll+1

      tmp1 = fiint(kk-1,0)/fiint(0,0)
      tmp2 = real(2-2*kk,kind=wp)
      tmp = Zero
      do mm=0,kk-1
        tmp3 = tmp1*binom(kk-1,mm)*(-One)**mm
        tmp4 = tmp2+real(2*mm+1,kind=wp)
        tmp = tmp+tmp3*(One/(real(2*ll,kind=wp)+tmp4)*rri(ll)+AA2/(real(2*ll+2,kind=wp)+tmp4)*rri(ll+1)+ &
                        AA4/(Three*(real(2*ll+4,kind=wp)+tmp4))*rri(ll+2))
      end do
      Grint(ll*2,kk) = pi4*tmp

      if (kk == 1) cycle
      tmp1 = fiint(kk-2,0)/fiint(0,0)
      tmp = Zero
      do mm=1,kk-1
        tmp3 = tmp1*binom(kk-2,mm-1)*(-One)**(mm+1)
        tmp4 = tmp2+real(2*mm+1)
        tmp = tmp-tmp3*(AA/(real(2*ll,kind=wp)+tmp4)*rri(ll)+AA3/(Three*(real(2*ll+2,kind=wp)+tmp4))*rri(ll+1)+ &
                        AA5/(15.0_wp*(real(2*ll+4,kind=wp)+tmp4))*rri(ll+2))
      end do
      Grint(ll*2-1,kk-1) = pi4*tmp
      Grint(ll*2-1,kk) = Zero

    end do
  end do

end if

!call RecPrt(' In RRint:grint',' ',grint,1+lmax,lmax)

end subroutine RRINT
