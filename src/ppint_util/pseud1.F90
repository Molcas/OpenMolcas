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

subroutine pseud1(ang,ccr,gout,ipt,lmnv,ltot1,ncr,nkcrl,nkcru,qsum,xab,yab,zab,zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs, &
                  lproju1,crda,crdb)
! compute type 1 core potential integrals

use Constants, only: Zero, One, Four, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipt(*), lmnv(3,*), ltot1, ncr(*), lit, ljt, kcrs, lproju1, nkcrl(lproju1,*), nkcru(lproju1,*)
real(kind=wp), intent(in) :: ccr(*), zcr(*), ai, aj, xi, yi, zi, xj, yj, zj, xc, yc, zc
real(kind=wp), intent(out) :: ang(ltot1,ltot1), qsum(ltot1,ltot1), xab(ltot1), yab(ltot1), zab(ltot1)
real(kind=wp), intent(inout) :: gout(*)
real(kind=wp), intent(out) :: crda(lit,3), crdb(ljt,3)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ijt, it, itl, itu, j, jt, jtl, jtu, kcrl, kcru, la1, lam, lamu, lb1, ma1, mb1, n, na1, nb1, nhi
real(kind=wp) :: aa, aaa, aarr1, alpt, arp2, fctr2, rk, rp, rp2, rr, s, taa, xij, xijm, xk, xka, xkb, yij, yijm, yk, yka, ykb, &
                 zij, zijm, zk, zka, zkb
integer(kind=iwp), parameter :: llt_(7,2) = reshape([1,2,5,11,21,36,57,1,4,10,20,35,56,84],shape(llt_))
real(kind=wp), parameter :: tol = 20.0_wp*log(Ten)

call pseud1_internal(Work)

! This is to allow type punning without an explicit interface
contains

subroutine pseud1_internal(a)

  use iso_c_binding
  real(kind=wp), intent(in), target :: a(*)
  integer(kind=iwp), pointer :: ia13(:), ia14(:), ia15(:), ia16(:), ia17(:)

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
  call rad1(aa,aarr1,alpt,arp2,ccr,a(ipt(11)),fctr2,kcrl,kcru,lamu,ltot1,ncr,qsum,rk,tol,zcr)
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

      call facab(a(ipt(12)),na1,nb1,crda(1,1),crdb(1,1),xab)
      call facab(a(ipt(12)),la1,lb1,crda(1,2),crdb(1,2),yab)
      call facab(a(ipt(12)),ma1,mb1,crda(1,3),crdb(1,3),zab)
      call c_f_pointer(c_loc(a(ipt(13))),ia13,[1])
      call c_f_pointer(c_loc(a(ipt(14))),ia14,[1])
      call c_f_pointer(c_loc(a(ipt(15))),ia15,[1])
      call c_f_pointer(c_loc(a(ipt(16))),ia16,[1])
      call c_f_pointer(c_loc(a(ipt(17))),ia17,[1])
      call ang1(ang,a(ipt(11)),na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,ia13,ia14,ia15,ia16,ia17,ltot1,xab,yab,zab,xk,yk,zk,a(ipt(18)))
      nullify(ia13,ia14,ia15,ia16,ia17)

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

  return

end subroutine pseud1_internal

end subroutine pseud1
