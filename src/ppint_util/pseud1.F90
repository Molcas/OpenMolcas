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

subroutine pseud1(ccr,gout,ltot1,ncr,nkcrl,nkcru,zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lproju1,crda,crdb)
! compute type 1 core potential integrals

use ppint_arrays, only: binom, dfac, lmf, lml, lmnv, lmx, lmy, lmz, zlm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ltot1, ncr(*), lit, ljt, kcrs, lproju1, nkcrl(lproju1,*), nkcru(lproju1,*)
real(kind=wp), intent(in) :: ccr(*), zcr(*), ai, aj, xi, yi, zi, xj, yj, zj, xc, yc, zc
real(kind=wp), intent(inout) :: gout(*)
real(kind=wp), intent(out) :: crda(lit,3), crdb(ljt,3)
integer(kind=iwp) :: i, ijt, it, itl, itu, j, jt, jtl, jtu, kcrl, kcru, la1, lam, lamu, lb1, ma1, mb1, n, na1, nb1, nhi
real(kind=wp) :: aa, aaa, aarr1, alpt, arp2, fctr2, rk, rp, rp2, rr, s, taa, xij, xijm, xk, xka, xkb, yij, yijm, yk, yka, ykb, &
                 zij, zijm, zk, zka, zkb
integer(kind=iwp), parameter :: llt_(7,2) = reshape([1,2,5,11,21,36,57,1,4,10,20,35,56,84],shape(llt_))
real(kind=wp), allocatable :: ang(:,:), qsum(:,:), xab(:), yab(:), zab(:)
real(kind=wp), parameter :: tol = 20.0_wp*log(Ten)

call mma_allocate(ang,ltot1,ltot1,label='ang')
call mma_allocate(qsum,ltot1,ltot1,label='qsum')
call mma_allocate(xab,ltot1,label='xab')
call mma_allocate(yab,ltot1,label='yab')
call mma_allocate(zab,ltot1,label='zab')

fctr2 = Four
itl = llt_(lit,1)
itu = llt_(lit,2)
jtl = llt_(ljt,1)
jtu = llt_(ljt,2)
aa = ai+aj
do i=1,3
  crda(1,i) = One
  crdb(1,i) = One
end do
xka = xc-xi
yka = yc-yi
zka = zc-zi
!ca = sqrt(xka*xka+yka*yka+zka*zka)
if (lit /= 1) then
  crda(2,1) = xka
  crda(2,2) = yka
  crda(2,3) = zka
  if (lit /= 2) then
    do i=1,3
      do j=3,lit
        crda(j,i) = crda(2,i)*crda(j-1,i)
      end do
    end do
  end if
end if
xkb = xc-xj
ykb = yc-yj
zkb = zc-zj
!cb = sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
if (ljt /= 1) then
  crdb(2,1) = xkb
  crdb(2,2) = ykb
  crdb(2,3) = zkb
  if (ljt /= 2) then
    do i=1,3
      do j=3,ljt
        crdb(j,i) = crdb(2,i)*crdb(j-1,i)
      end do
    end do
  end if
end if
xij = Half*(xi+xj)
yij = Half*(yi+yj)
zij = Half*(zi+zj)
xijm = Half*(xi-xj)
yijm = Half*(yi-yj)
zijm = Half*(zi-zj)
rr = (xi-xj)**2+(yi-yj)**2+(zi-zj)**2
aaa = (ai-aj)/aa
xk = xij+aaa*xijm-xc
yk = yij+aaa*yijm-yc
zk = zij+aaa*zijm-zc
aarr1 = (ai*aj/aa)*rr
taa = aa+aa

rp2 = xk*xk+yk*yk+zk*zk
if (rp2 == Zero) then
  rp = Zero
  arp2 = Zero
  alpt = Zero
  rk = Zero
  lamu = 1
else
  rp = sqrt(rp2)
  xk = xk/rp
  yk = yk/rp
  zk = zk/rp
  arp2 = aa*rp2
  alpt = aa*arp2
  rk = taa*rp
  lamu = ltot1
end if

! compute radial integrals and sum over potential terms

qsum(:,:) = Zero
kcrl = nkcrl(1,kcrs)
kcru = nkcru(1,kcrs)
call rad1(aa,aarr1,alpt,arp2,ccr,dfac,fctr2,kcrl,kcru,lamu,ltot1,ncr,qsum,rk,tol,zcr)
ijt = 0
do it=itl,itu
  na1 = lmnv(1,it)+1
  la1 = lmnv(2,it)+1
  ma1 = lmnv(3,it)+1
  do jt=jtl,jtu
    ijt = ijt+1
    s = Zero
    nb1 = lmnv(1,jt)+1
    lb1 = lmnv(2,jt)+1
    mb1 = lmnv(3,jt)+1

    ! compute angular integrals

    call facab(binom,na1,nb1,crda(1,1),crdb(1,1),xab)
    call facab(binom,la1,lb1,crda(1,2),crdb(1,2),yab)
    call facab(binom,ma1,mb1,crda(1,3),crdb(1,3),zab)
    call ang1(ang,dfac,na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,lmf,lml,lmx,lmy,lmz,ltot1,xab,yab,zab,xk,yk,zk,zlm)

    ! combine angular and radial integrals

    do lam=1,lamu
      nhi = ltot1-mod(ltot1-lam,2)
      do n=lam,nhi,2
        s = s+ang(n,lam)*qsum(n,lam)
      end do
    end do
    gout(ijt) = gout(ijt)+s
  end do
end do

call mma_deallocate(ang)
call mma_deallocate(qsum)
call mma_deallocate(xab)
call mma_deallocate(yab)
call mma_deallocate(zab)

return

end subroutine pseud1
