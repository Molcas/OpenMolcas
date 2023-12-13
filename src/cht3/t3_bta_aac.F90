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

subroutine t3_bta_aac(nuga,nugc,kab,kca,kac,kc,la,lxa,lxc,mi,mij,adim,cdim,N,noab_a,noab_b,lu,iasblock,nga,ngc,oehi,oehk,oepa, &
                      oepc,enx,vab,vca,t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nuga, nugc, adim, cdim, N, noab_a, noab_b, lu(6), iasblock(5), nga, ngc
real(kind=wp), intent(_OUT_) :: kab(nTri_Elem(adim-1),N,*), kca(adim*cdim,N,*), kac(adim*cdim,N,*), kc(*), la(N*adim,*), &
                                lxa(N*adim,*), lxc(N*cdim,*), mi(cdim*nTri_Elem(adim-1),*), mij(*), vab(nTri_Elem(adim-1),*), &
                                vca(adim*cdim,*), t3a(*), t3b(*)
real(kind=wp), intent(in) :: oehi(*), oehk(*), oepa(*), oepc(*)
real(kind=wp), intent(inout) :: enx, t1aa(noab_a,*), t1ba(noab_a,*), t1ac(noab_b,*), t1bc(noab_b,*)
logical(kind=iwp), intent(in) :: ifvo
integer(kind=iwp) :: a, ab, abb, b, bab, c, i, ias, iasabi, iasack, iascai, ij, ik, j, jk, k, ki, kj, nadim, ncdim, ngab_offset, &
                     ngac_offset, ngca_offset, nno_a, nnoab, nuga_offset, nugc_offset
real(kind=wp) :: den, dena, denb, denc, xx, yy

! iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
if (adim == 1) return
nno_a = nTri_Elem(noab_a-1)
nnoab = noab_a*noab_b
nadim = nTri_Elem(adim-1)
ncdim = adim*cdim
nuga_offset = iasblock(1)*nTri_Elem(nuga)
nugc_offset = iasblock(1)*nuga*nugc
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nno_a*adim*N,lu(2),ias)
ias = iasblock(3)*(nga-1)+1
call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
ias = iasblock(3)*(ngc-1)+1
call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)
! vvoo ints reading
ngab_offset = iasblock(4)*(nTri_Elem(nga-1)+nga-1)+1
ias = iasblock(2)*nuga+ngab_offset
call multi_readir(vab,nno_a*nadim,lu(2),ias)
ngca_offset = iasblock(5)*(nugc*(nga-1)+ngc-1)+1
ias = iasblock(2)*nuga+iasblock(4)*nTri_Elem(nuga)+ngca_offset
call multi_readir(vca,nnoab*adim*cdim,lu(2),ias)
!!ngac_offset = iasblock(5)*(nuga*(ngc-1)+nga-1)+1
! end readin vvoo ints
ngab_offset = iasblock(1)*(nTri_Elem(nga-1)+nga-1)+1
ngac_offset = iasblock(1)*(nuga*(ngc-1)+nga-1)+1
ngca_offset = iasblock(1)*(nugc*(nga-1)+ngc-1)+1

do i=1,noab_a
  iasabi = (i-1)*nuga_offset+ngab_offset
  call multi_readir(kab(:,:,i),N*nadim,lu(1),iasabi)
end do
do i=1,noab_a
  iascai = (i-1)*nugc_offset+ngca_offset
  call multi_readir(kca(:,:,i),N*ncdim,lu(3),iascai)
end do
do k=1,noab_b
  do i=1,noab_a
    ik = (k-1)*noab_a+i
    ki = (i-1)*noab_b+k
    ! K_ab^ir x L_rc^ik     cba
    call DGEMM_('T','T',cdim,nadim,N,one,lxc(:,ik),N,kab(:,:,i),nadim,zero,mi(:,i),cdim)

    ! K_ac^ir x L_rb^ki     cab
    call DGEMM_('N','N',ncdim,adim,N,one,kca(:,:,i),ncdim,lxa(:,ki),N,zero,t3b,ncdim)
    ab = 1
    do a=2,adim
      abb = (a-1)*cdim+1
      bab = (a-1)*ncdim+1
      do b=1,a-1
        mi(ab:ab+cdim-1,i) = mi(ab:ab+cdim-1,i)-t3b(abb:abb+cdim-1)+t3b(bab:bab+cdim-1)
        ab = ab+cdim
        abb = abb+ncdim
        bab = bab+cdim
      end do
    end do
  end do    ! i
  ! end prefactors
  iasack = (k-1)*nugc_offset+ngac_offset
  call multi_readir(kac,N*ncdim,lu(4),iasack)
  ij = 0
  do i=2,noab_a
    ki = (i-1)*noab_b+k
    ik = (k-1)*noab_a+i
    kj = k-noab_b
    jk = (k-1)*noab_a
    do j=1,i-1
      ij = ij+1
      kj = kj+noab_b
      jk = jk+1
      ! K_bc^kr x L_ra^ij
      call DGEMM_('N','N',ncdim,adim,N,one,kac,ncdim,la(:,ij),N,zero,t3a,ncdim)
      ! transpose the first two inicesd
      ab = 1
      do a=1,adim
        call map2_21_t3(t3a(ab),t3b(ab),adim,cdim)
        ab = ab+ncdim
      end do
      ! K_ab^ir x L_rc^jk -K_ab^jr x L_rc^ik
      kc(1:N*nadim) = pack(kab(:,:,i)-kab(:,:,j),.true.)
      mij(1:N*cdim) = lxc(:,jk)+lxc(:,ik)
      call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,zero,t3a,cdim)
      ! K_ab^ir x L_rc^jk
      !!call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,zero,t3a,cdim)
      ! -K_ab^jr x L_rc^ik
      !!call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,one,t3a,cdim)

      ! K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki
      kc(1:N*ncdim) = pack(kca(:,:,i)-kca(:,:,j),.true.)
      mij(1:N*adim) = lxa(:,kj)+lxa(:,ki)
      call DGEMM_('N','N',ncdim,adim,N,one,kc,ncdim,mij,N,one,t3b,ncdim)
      ! K_bc^ir x L_ra^kj
      !!call DGEMM_('N','N',ncdim,adim,N,one,ka,ncdim,lxa(1,kj),N,one,t3b,ncdim)
      ! -K_bc^jr x L_ra^ki
      !!call DGEMM_('N','N',ncdim,adim,N,-one,kb,ncdim,lxa(1,ki),N,one,t3b,ncdim)
      ab = 1
      do a=2,adim
        abb = (a-1)*cdim+1
        bab = (a-1)*ncdim+1
        do b=1,a-1
          t3a(ab:ab+cdim-1) = t3a(ab:ab+cdim-1)-t3b(abb:abb+cdim-1)+t3b(bab:bab+cdim-1)
          ab = ab+cdim
          abb = abb+ncdim
          bab = bab+cdim
        end do
      end do

      t3a(1:nadim*cdim) = t3a(1:nadim*cdim)-mi(:,i)+mi(:,j)
      den = oehi(i)+oehi(j)+oehk(k)
      ab = 0
      do a=2,adim
        dena = den-oepa(a)
        do b=1,a-1
          denb = dena-oepa(b)
          do c=1,cdim
            denc = denb-oepc(c)
            ab = ab+1
            xx = t3a(ab)
            yy = xx/denc
            enx = enx+yy*xx
            t3a(ab) = yy
            !! t1aa(j,a) = t1aa(j,a)-yy*vca((a-1)*cdim+c,ki)
          end do
        end do
      end do
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
      call DGEMM_('N','T',1,cdim,nadim,one,vab(:,ij),1,t3a,cdim,one,t1ac(k,1),noab_b)
      call DGEMM_('N','N',1,adim,ncdim,one,vca(:,kj),1,t3b,ncdim,one,t1aa(i,1),noab_a)
      call DGEMM_('N','N',1,adim,ncdim,-one,vca(:,ki),1,t3b,ncdim,one,t1aa(j,1),noab_a)
      ! ccsd(T) part t2*t3
      if (ifvo) then
        call DGEMM_('N','T',1,cdim,nadim,one,kab(:,i,j),1,t3a,cdim,one,t1bc(k,1),noab_b)
        call DGEMM_('N','N',1,adim,ncdim,-one,kca(:,k,j),1,t3b,ncdim,one,t1ba(i,1),noab_a)
        call DGEMM_('N','N',1,adim,ncdim,one,kca(:,k,i),1,t3b,ncdim,one,t1ba(j,1),noab_a)
      end if
    end do !j
  end do !i
end do !k

return

end subroutine t3_bta_aac
