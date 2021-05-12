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

subroutine ang2(ang,binom,crda,dfac,it,l,lit,lmlo,lmhi,lmf,lml,lmx,lmy,lmz,mproju,xk,yk,zk,zlm)
! compute type 2 angular integrals

use ppint_arrays, only: lmnv
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: it, l, lit, lmlo, lmhi, lmf(*), lml(*), lmx(*), lmy(*), lmz(*), mproju
real(kind=wp), intent(out) :: ang(lit,mproju,lmhi)
real(kind=wp), intent(in) :: binom(*), crda(lit,3), dfac(*), xk, yk, zk, zlm(*)
integer(kind=iwp) :: i, ia, ib, ic, iend, indx, indy, indz, istart, j, l2, la1, laind, lam, lamhi, lamlo, loc1, loc2, m, ma1, &
                     maind, mend, mhi, mndx, mndy, mndz, mstart, mu, n, na1, naind
real(kind=wp) :: a_int, angt, pab1, pab2, pab3, pre, xkp, ykp, zkp

ang(:,:,:) = Zero
na1 = lmnv(1,it)+1
la1 = lmnv(2,it)+1
ma1 = lmnv(3,it)+1
naind = (na1*(na1-1))/2
laind = (la1*(la1-1))/2
maind = (ma1*(ma1-1))/2
loc1 = (l-1)**2
mhi = l+l-1
do ia=1,na1
  pab1 = binom(naind+ia)*crda((na1+1)-ia,1)
  if (pab1 == Zero) cycle
  do ib=1,la1
    pab2 = pab1*binom(laind+ib)*crda((la1+1)-ib,2)
    if (pab2 == Zero) cycle
    do ic=1,ma1
      pab3 = pab2*binom(maind+ic)*crda((ma1+1)-ic,3)
      if (pab3 == Zero) cycle
      n = ((ia-3)+ib)+ic
      lamlo = max(l-n,lmlo+mod(l+n+lmlo,2))
      lamhi = min(l+n,lmhi-mod(l+n+lmhi,2))
      if (lamlo > lamhi) cycle
      do m=1,mhi
        mstart = lmf(loc1+m)
        mend = lml(loc1+m)
        do lam=lamlo,lamhi,2
          l2 = lam+lam-1
          angt = Zero
          loc2 = (lam-1)**2
          do mu=1,l2
            istart = lmf(loc2+mu)
            if ((mod(ia+lmx(mstart)+lmx(istart),2) /= 1) .or. &
                (mod(ib+lmy(mstart)+lmy(istart),2) /= 1) .or. &
                (mod(ic+lmz(mstart)+lmz(istart),2) /= 1)) cycle
            pre = Zero
            iend = lml(loc2+mu)
            a_int = Zero
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
              do j=mstart,mend
                mndx = lmx(j)
                mndy = lmy(j)
                mndz = lmz(j)
                a_int = a_int+zlm(i)*zlm(j)*dfac(ia+indx+mndx)*dfac(ib+indy+mndy)*dfac(ic+indz+mndz)/ &
                        dfac(ia+indx+mndx+ib+indy+mndy+ic+indz+mndz)
              end do
            end do
            angt = angt+pre*a_int
          end do
          ang(n+1,m,lam) = ang(n+1,m,lam)+angt*pab3
        end do
      end do
    end do
  end do
end do

return

end subroutine ang2
