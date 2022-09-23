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

subroutine t3_bta_abc(nuga,nugc,kab,kcb,kca,kac,kbc,kc,la,lb,lxa,lxb,lxc,mi,mij,adim,bdim,cdim,N,noab_a,noab_b,lu,iasblock,nga, &
                      ngb,ngc,oehi,oehk,oepa,oepb,oepc,enx,vab,vcb,vca,t1aa,t1ba,t1ab,t1bb,t1ac,t1bc,t3a,t3b,ifvo)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nuga, nugc, adim, bdim, cdim,N, noab_a, noab_b, lu(6), iasblock(5), nga, ngb, ngc
real(kind=wp), intent(_OUT_) :: kab(adim*bdim,N,*), kcb(cdim*bdim,N,*), kca(cdim*adim,N,*), kc(*), la(N*adim,*), lb(N*bdim,*), &
                                lxa(N*adim,*), lxb(N*bdim,*), lxc(N*cdim,*), mi(adim*bdim*cdim,*), mij(*), vab(adim*bdim,*), &
                                vcb(bdim*cdim,*), vca(adim*cdim,*), t3a(*), t3b(*)
real(kind=wp), intent(out) :: kac(cdim*adim,N), kbc(cdim*bdim,N)
real(kind=wp), intent(in) :: oehi(*), oehk(*), oepa(*), oepb(*), oepc(*)
real(kind=wp), intent(inout) :: enx, t1aa(noab_a,*), t1ba(noab_a,*), t1ab(noab_a,*), t1bb(noab_a,*), t1ac(noab_b,*), t1bc(noab_b,*)
logical(kind=iwp), intent(in) :: ifvo
integer(kind=iwp) :: a, ab, b, ba, c, i, ias, iasabi, iasack, iasbck, iascai, iascbi, ij, ik, j, jk, k, ki, kj, nadim, nbdim, &
                     ncdim, ngab_offset, ngac_offset, ngbc_offset, ngca_offset, ngcb_offset, nno_a, nnoab, nuga_offset, nugc_offset
real(kind=wp) :: den, dena, denb, denc, xx, yy

! iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
!!sumt3 = Zero
nno_a = nTri_Elem(noab_a-1)
nnoab = noab_a*noab_b
nadim = adim*bdim
nbdim = bdim*cdim
ncdim = adim*cdim
nuga_offset = iasblock(1)*nTri_Elem(nuga)
nugc_offset = iasblock(1)*nuga*nugc
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nno_a*adim*N,lu(2),ias)
ias = iasblock(2)*(ngb-1)+1
call multi_readir(lb,nno_a*bdim*N,lu(2),ias)
ias = iasblock(3)*(nga-1)+1
call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
ias = iasblock(3)*(ngb-1)+1
call multi_readir(lxb,nnoab*bdim*N,lu(5),ias)
ias = iasblock(3)*(ngc-1)+1
call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)

! vvoo ints reading
ngab_offset = iasblock(4)*(nTri_Elem(nga-1)+ngb-1)+1
ias = iasblock(2)*nuga+ngab_offset
call multi_readir(vab,nno_a*nadim,lu(2),ias)
ngca_offset = iasblock(5)*(nugc*(nga-1)+ngc-1)+1
ias = iasblock(2)*nuga+iasblock(4)*nTri_Elem(nuga)+ngca_offset
call multi_readir(vca,nnoab*ncdim,lu(2),ias)
ngcb_offset = iasblock(5)*(nugc*(ngb-1)+ngc-1)+1
ias = iasblock(2)*nuga+iasblock(4)*nTri_Elem(nuga)+ngcb_offset
call multi_readir(vcb,nnoab*nbdim,lu(2),ias)
! end readin vvoo ints

ngab_offset = iasblock(1)*(nTri_Elem(nga-1)+ngb-1)+1
ngac_offset = iasblock(1)*(nuga*(ngc-1)+nga-1)+1
ngbc_offset = iasblock(1)*(nuga*(ngc-1)+ngb-1)+1
ngca_offset = iasblock(1)*(nugc*(nga-1)+ngc-1)+1
ngcb_offset = iasblock(1)*(nugc*(ngb-1)+ngc-1)+1
! saves reading:
do i=1,noab_a
  iasabi = (i-1)*nuga_offset+ngab_offset
  call multi_readir(kab(:,:,i),N*nadim,lu(1),iasabi)
end do
do i=1,noab_a
  iascai = (i-1)*nugc_offset+ngca_offset
  call multi_readir(kca(:,:,i),N*ncdim,lu(3),iascai)
  !mp
  !mp write(u6,*) 'ze tak teraz ju dam',i
  !mp call check_mat(kca(1,1,i),1,N*ncdim)
  !mp
end do
do i=1,noab_a
  iascbi = (i-1)*nugc_offset+ngcb_offset
  call multi_readir(kcb(:,:,i),N*nbdim,lu(3),iascbi)
end do
do k=1,noab_b
  iasbck = (k-1)*nugc_offset+ngbc_offset
  call multi_readir(kbc,N*nbdim,lu(4),iasbck)
  iasack = (k-1)*nugc_offset+ngac_offset
  call multi_readir(kac,N*ncdim,lu(4),iasack)
  ! start calculating prefactors:
  !mp
  do i=1,noab_a
    ik = (k-1)*noab_a+i
    ki = (i-1)*noab_b+k
    !mp
    ! K_ab^ir x L_rc^ik     cba
    call DGEMM_('T','T',cdim,nadim,N,one,lxc(:,ik),N,kab(:,:,i),nadim,zero,mi(:,i),cdim)
    ! K_bc^ir x L_ra^ki          cba
    call DGEMM_('N','N',nbdim,adim,N,one,kcb(:,:,i),nbdim,lxa(:,ki),N,one,mi(:,i),nbdim)
    !mp
    !mp do imm=0,adim*nbdim-1
    !mp   if (abs(mi(1+imm,i)) > 1.0e5_wp) then
    !mp     write(u6,*) 'uz mi dojebane 2',imm+1,i,mi(1+imm,i)
    !mp     stop
    !mp   end if
    !mp end do
    !mp

    ! K_ac^ir x L_rb^ki     cab
    call DGEMM_('N','N',ncdim,bdim,N,one,kca(:,:,i),ncdim,lxb(:,ki),N,zero,t3b,ncdim)
    !mp
    ab = 1
    do a=1,adim
      ba = (a-1)*cdim+1
      do b=1,bdim
        mi(ab:ab+cdim-1,i) = mi(ab:ab+cdim-1,i)-t3b(ba:ba+cdim-1)
        ab = ab+cdim
        ba = ba+ncdim
      end do
    end do
  end do    ! i
  !mp
  ! end prefactors
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
      !mp
      ! K_ac^kr x L_rb^ij        bac
      call DGEMM_('T','T',bdim,ncdim,N,-one,lb(:,ij),N,kac,ncdim,zero,t3b,bdim)
      ! K_bc^kr x L_ra^ij        abc
      call DGEMM_('T','T',adim,nbdim,N,one,la(:,ij),N,kbc,nbdim,zero,t3a,adim)
      ! transpose the first two indices
      !mp
      ab = 1
      do c=1,cdim
        ba = (c-1)*nadim+1
        do b=1,bdim
          call daxpy_(adim,One,t3a(ab),1,t3b(ba),bdim)
          ba = ba+1
          ab = ab+adim
        end do
      end do
      ! t3b  bac
      call map2_21_t3(t3b,t3a,nadim,cdim)
      ! cba in t3a
      ! K_ab^ir x L_rc^jk  -K_ab^jr x L_rc^ik
      !call ka(1:N*nadim) = ka(1:N*nadim)-kb(1:N*nadim)
      kc(1:N*nadim) = pack(kab(:,:,i)-kab(:,:,j),.true.)
      mij(1:N*cdim) = lxc(:,jk)+lxc(:,ik)
      call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,one,t3a,cdim)
      ! K_ab^ir x L_rc^jk
      !!call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,one,t3a,cdim)
      ! -K_ab^jr x L_rc^ik
      !!call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,one,t3a,cdim)

      ! K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki         cba
      kc(1:N*nbdim) = pack(kcb(:,:,i)-kcb(:,:,j),.true.)
      mij(1:N*adim) = lxa(:,kj)+lxa(:,ki)
      call DGEMM_('N','N',nbdim,adim,N,one,kc,nbdim,mij,N,one,t3a,nbdim)
      ! K_bc^ir x L_ra^kj          cba
      !!call DGEMM_('N','N',nbdim,adim,N,one,ka,nbdim,lxa(1,kj),N,one,t3a,nbdim)
      ! -K_bc^jr x L_ra^ki         cba
      !!call DGEMM_('N','N',nbdim,adim,N,-one,kb,nbdim,lxa(1,ki),N,one,t3a,nbdim)
      ! K_ac^ir x L_rb^kj  -K_ac^jr x L_rb^ki      cab
      kc(1:N*ncdim) = pack(kca(:,:,i)-kca(:,:,j),.true.)
      mij(1:N*bdim) = lxb(:,kj)+lxb(:,ki)
      call DGEMM_('N','N',ncdim,bdim,N,one,kc,ncdim,mij,N,zero,t3b,ncdim)
      ! K_ac^ir x L_rb^kj      cab
      !!call DGEMM_('N','N',ncdim,bdim,N,one,ka,ncdim,lxb(1,kj),N,zero,t3b,ncdim)
      ! -K_ac^jr x L_rb^ki     cab
      !!call DGEMM_('N','N',ncdim,bdim,N,-one,kb,ncdim,lxb(1,ki),N,one,t3b,ncdim)
      !mp
      ab = 1
      do a=1,adim
        ba = (a-1)*cdim+1
        do b=1,bdim
          t3a(ab:ab+cdim-1) = t3a(ab:ab+cdim-1)-t3b(ba:ba+cdim-1)
          ab = ab+cdim
          ba = ba+ncdim
        end do
      end do

      !mp do a=1,adim
      !mp   do b=1,bdim
      !mp     do c=1,cdim
      !mp       ab = c+(b-1)*cdim+(a-1)*bdim*cdim
      !mp       ba = c+(a-1)*cdim+(b-1)*adim*cdim
      !mp       t3a(ab) = t3a(ab)-t3b(ba)
      !mp     end do
      !mp   end do
      !mp end do

      t3a(1:nadim*cdim) = t3a(1:nadim*cdim)-mi(:,i)+mi(:,j)
      den = oehi(i)+oehi(j)+oehk(k)
      ab = 0
      do a=1,adim
        ba = ab+1
        dena = den-oepa(a)
        do b=1,bdim
          denb = dena-oepb(b)
          do c=1,cdim
            denc = denb-oepc(c)
            ab = ab+1
            !mp if ((i == 14) .and. (j == 1) .and. (k == 1)) then
            !mp   write(u6,'(A,3(i5,x),3(f18.10,x))') 'a,b,c, oepa(a), oepb(b), oepc(c)',a,b,c,oepa(a),oepb(b),oepc(c)
            !mp   write(u6,*) 'denc = ',denc
            !mp   write(u6,*) 'ab, t3a(ab) ',ab,t3a(ab)
            !mp end if
            xx = t3a(ab)
            yy = xx/denc
            !!sumt3 = sumt3+xx
            enx = enx+yy*xx
            !mp if ((i == 14) .and. (j == 1) .and. (k == 1)) write(u6,'(A,3(f18.10,x))') 'xx,yy,enx = ',xx,yy,enx
            t3a(ab) = yy
          end do
        end do
        call map2_21_t3(t3a(ba),t3b(ba),cdim,bdim)
      end do
      call DGEMM_('N','T',1,cdim,nadim,one,vab(:,ij),1,t3a,cdim,one,t1ac(k,1),noab_b)
      call DGEMM_('N','N',1,adim,nbdim,-one,vcb(:,kj),1,t3a,nbdim,one,t1aa(i,1),noab_a)
      call DGEMM_('N','N',1,adim,nbdim,one,vcb(:,ki),1,t3a,nbdim,one,t1aa(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,-one,vca(:,ki),1,t3b,bdim,one,t1ab(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,one,vca(:,kj),1,t3b,bdim,one,t1ab(i,1),noab_a)
      if (ifvo) then
        call DGEMM_('N','T',1,cdim,nadim,one,kab(:,i,j),1,t3a,cdim,one,t1bc(k,1),noab_b)
        call DGEMM_('N','N',1,adim,nbdim,one,kcb(:,k,j),1,t3a,nbdim,one,t1ba(i,1),noab_a)
        call DGEMM_('N','N',1,adim,nbdim,-one,kcb(:,k,i),1,t3a,nbdim,one,t1ba(j,1),noab_a)
        call DGEMM_('N','T',1,bdim,ncdim,one,kca(:,k,i),1,t3b,bdim,one,t1bb(j,1),noab_a)
        call DGEMM_('N','T',1,bdim,ncdim,-one,kca(:,k,j),1,t3b,bdim,one,t1bb(i,1),noab_a)
      end if
      !mp write(u6,'(3(i5,x),f18.10)') i,j,k,enx
    end do !j
  end do !i
end do !k

return

end subroutine t3_bta_abc
