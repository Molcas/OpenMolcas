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

subroutine pseud2(ccr,gout,lambu,ltot1,mproju,ncr,nkcrl,nkcru,zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lcru,lproju1,crda, &
                  crdb)
! compute type 2 core potential integrals

use ppint_arrays, only: binom, dfac, lmf, lml, lmx, lmy, lmz, zlm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lambu, ltot1, mproju, ncr(*), lit, ljt, kcrs, lcru, lproju1, nkcrl(lproju1,*), nkcru(lproju1,*)
real(kind=wp), intent(inout) :: gout(*)
real(kind=wp), intent(in) :: ccr(*), zcr(*), ai, aj, xi, yi, zi, xj, yj, zj, xc, yc, zc
real(kind=wp), intent(out) :: crda(lit,3), crdb(ljt,3)
integer(kind=iwp) :: i, ijt, inc, it, itl, itu, j, jt, jtl, jtu, kcrl, kcru, l, lama, lamb, ldifa1, ldifb, lhi, llo, lmahi, lmalo, &
                     lmau, lmbhi, lmblo, lmbu, m, mhi, n, nhi, nlma, nlmahi, nlmalo, nlmau, nlmbu, nlo
real(kind=wp) :: aa, aarr2, angp, ca, cb, fctr2, rka, rkb, s, tai, taj, xka, xkb, yka, ykb, zka, zkb
real(kind=wp), allocatable :: anga(:,:,:), angb(:,:,:), qsum(:,:,:)
integer(kind=iwp), parameter :: llt_(7,2) = reshape([1,2,5,11,21,36,57,1,4,10,20,35,56,84],shape(llt_))

fctr2 = Four
itl = llt_(lit,1)
itu = llt_(lit,2)
jtl = llt_(ljt,1)
jtu = llt_(ljt,2)
tai = ai+ai
taj = aj+aj
aa = ai+aj
do i=1,3
  crda(1,i) = One
  crdb(1,i) = One
end do
xka = xc-xi
yka = yc-yi
zka = zc-zi
ca = sqrt(xka*xka+yka*yka+zka*zka)
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
cb = sqrt(xkb*xkb+ykb*ykb+zkb*zkb)
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
aarr2 = (ai*aj/aa)*(ca-cb)**2

if (ca == Zero) then
  rka = Zero
  lmau = 1
else
  xka = -xka/ca
  yka = -yka/ca
  zka = -zka/ca
  rka = tai*ca
  lmau = lcru+(lit-1)
end if
if (cb == Zero) then
  rkb = Zero
  lmbu = 1
else
  xkb = -xkb/cb
  ykb = -ykb/cb
  zkb = -zkb/cb
  rkb = taj*cb
  lmbu = lcru+(ljt-1)
end if
if ((ca == Zero) .and. (cb == Zero)) then
  lhi = min(lcru,lit,ljt)
  llo = mod((lit-1),2)+1
  if ((llo /= mod((ljt-1),2)+1) .or. (llo > lhi)) return
  inc = 2
else if (ca == Zero) then
  lhi = min(lcru,lit)
  llo = mod((lit-1),2)+1
  if (llo > lhi) return
  inc = 2
else if (cb == Zero) then
  lhi = min(lcru,ljt)
  llo = mod((ljt-1),2)+1
  if (llo > lhi) return
  inc = 2
else
  lhi = lcru
  llo = 1
  inc = 1
end if

call mma_allocate(anga,lit,mproju,lcru+lit-1,label='anga')
call mma_allocate(angb,ljt,mproju,lcru+ljt-1,label='angb')
call mma_allocate(qsum,ltot1,lambu,lcru+lit-1,label='qsum')

do l=llo,lhi,inc
  mhi = l+l-1
  lmalo = max(l-(lit-1),1)
  lmahi = min(lmau,l+(lit-1))
  lmblo = max(l-(ljt-1),1)
  lmbhi = min(lmbu,l+(ljt-1))

  ! compute radial integrals

  kcrl = nkcrl(l+1,kcrs)
  kcru = nkcru(l+1,kcrs)
  call rad2(ccr,kcrl,kcru,l,lambu,lmahi,lmalo,lmbhi,lmblo,ltot1,ncr,qsum,rka,rkb,zcr,lit,ljt,ca,cb,tai,taj,aa,aarr2,fctr2)

  ! compute angular integrals and combine with radial integrals

  ijt = 0
  do it=itl,itu
    call ang2(anga,binom,crda,dfac,it,l,lit,lmalo,lmahi,lmf,lml,lmx,lmy,lmz,mproju,xka,yka,zka,zlm)
    do jt=jtl,jtu
      ijt = ijt+1
      s = Zero
      call ang2(angb,binom,crdb,dfac,jt,l,ljt,lmblo,lmbhi,lmf,lml,lmx,lmy,lmz,mproju,xkb,ykb,zkb,zlm)
      do lama=lmalo,lmahi
        ldifa1 = abs(l-lama)+1
        nlmau = lit-mod(lit-ldifa1,2)
        do lamb=lmblo,lmbhi
          ldifb = abs(l-lamb)
          nlmbu = (ljt-1)-mod((ljt-1)-ldifb,2)
          nlo = ldifa1+ldifb
          nhi = nlmau+nlmbu
          do n=nlo,nhi,2
            nlmalo = max(ldifa1,n-nlmbu)
            nlmahi = min(nlmau,n-ldifb)
            angp = Zero
            do m=1,mhi
              do nlma=nlmalo,nlmahi,2
                angp = angp+anga(nlma,m,lama)*angb((n+1)-nlma,m,lamb)
              end do
            end do
            s = s+angp*qsum(n,lamb,lama)
          end do
        end do
      end do
      gout(ijt) = gout(ijt)+s
    end do
  end do
end do

call mma_deallocate(anga)
call mma_deallocate(angb)
call mma_deallocate(qsum)

return

end subroutine pseud2
