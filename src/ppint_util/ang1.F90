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

subroutine ang1(ang,dfac,nanb,lalb,mamb,lamu,lmf,lml,lmx,lmy,lmz,ltot1,xab,yab,zab,xk,yk,zk,zlm)
! compute type 1 angular integrals

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ltot1, nanb, lalb, mamb, lamu, lmf(*), lml(*), lmx(*), lmy(*), lmz(*)
real(kind=wp), intent(out) :: ang(ltot1,lamu)
real(kind=wp), intent(in) :: dfac(*), xab(*), yab(*), zab(*), xk, yk, zk, zlm(*)
integer(kind=iwp) :: i, iend, indx, indy, indz, istart, l, l2, lam, lamhi, lamlo, loc, m, mu1, n, nlm
real(kind=wp) :: a_int, angt, pre, xkp, ykp, zkp

ang(:,:) = Zero
do n=1,nanb
  if (xab(n) == Zero) cycle
  do l=1,lalb
    if (yab(l) == Zero) cycle
    do m=1,mamb
      if (zab(m) == Zero) cycle
      nlm = ((n-2)+l)+m
      lamlo = mod(nlm-1,2)+1
      lamhi = min(nlm,lamu)
      if (lamlo > lamhi) cycle
      do lam=lamlo,lamhi,2
        l2 = lam+lam-1
        angt = Zero
        loc = (lam-1)**2
        do mu1=1,l2
          istart = lmf(loc+mu1)
          if ((mod(n,2) == mod(lmx(istart),2)) .or. &
              (mod(l,2) == mod(lmy(istart),2)) .or. &
              (mod(m,2) == mod(lmz(istart),2))) cycle
          pre = Zero
          a_int = Zero
          iend = lml(loc+mu1)
          do i=istart,iend
            indx = lmx(i)
            indy = lmy(i)
            indz = lmz(i)
            if (indx == 0) then
              xkp = One
            else
              xkp = xk**indx
            end if
            if (indy == 0) then
              ykp = One
            else
              ykp = yk**indy
            end if
            if (indz == 0) then
              zkp = One
            else
              zkp = zk**indz
            end if
            pre = pre+zlm(i)*xkp*ykp*zkp
            a_int = a_int+zlm(i)*dfac(n+indx)*dfac(l+indy)*dfac(m+indz)/dfac((n+indx)+(l+indy)+(m+indz))
          end do
          angt = angt+pre*a_int
        end do
        ang(nlm,lam) = ang(nlm,lam)+((xab(n)*yab(l))*zab(m))*angt
      end do
    end do
  end do
end do

return

end subroutine ang1
