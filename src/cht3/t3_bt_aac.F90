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

subroutine t3_bt_aac(nug,ka,kc,la,lc,adim,cdim,N,noab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,t1aa,t1ba,t1ac,t1bc,t3a, &
                     t3b,ifvo)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nug, adim, cdim, N, noab, nnoab, lu(2), iasblock(3), nga, ngc
real(kind=wp), intent(_OUT_) :: ka(nTri_Elem(adim-1),N,*), kc(adim*cdim,N,*), voa(nTri_Elem(adim-1),*), voc(adim*cdim,*), &
                                t3a(cdim,*), t3b(cdim,adim,*)
real(kind=wp), intent(out) :: la(N*adim,nnoab), lc(N*cdim,nnoab)
real(kind=wp), intent(in) :: oeh(noab), oepa(adim), oepc(cdim)
real(kind=wp), intent(inout) :: enx, t1aa(noab,*), t1ba(noab,*), t1ac(noab,*), t1bc(noab,*)
logical(kind=iwp), intent(in) :: ifvo
integer(kind=iwp) :: a, ab, b, c, i, ias, iasai, iasci, ij, ik, j, jk, k, nadim, ncdim, nga_offset, ngc_offset, nug_offset
real(kind=wp) :: den, dena, denb, denc, sumt3, xx, yy

sumt3 = Zero
if (adim == 1) return
nadim = nTri_Elem(adim-1)
ncdim = adim*cdim
nug_offset = iasblock(1)*nTri_Elem(nug)
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nnoab*adim*N,lu(2),ias)
ias = iasblock(2)*(ngc-1)+1
call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
! reads vvoo
nga_offset = iasblock(3)*(nTri_Elem(nga-1)+nga-1)+1
ngc_offset = iasblock(3)*(nTri_Elem(nga-1)+ngc-1)+1
ias = iasblock(2)*nug+nga_offset
call multi_readir(voa,nnoab*nadim,lu(2),ias)
ias = iasblock(2)*nug+ngc_offset
call multi_readir(voc,nnoab*ncdim,lu(2),ias)

nga_offset = iasblock(1)*(nTri_Elem(nga-1)+nga-1)+1
ngc_offset = iasblock(1)*(nTri_Elem(nga-1)+ngc-1)+1
do i=1,noab
  iasci = (i-1)*nug_offset+ngc_offset
  call multi_readir(kc(:,:,i),N*ncdim,lu(1),iasci)
  iasai = (i-1)*nug_offset+nga_offset
  call multi_readir(ka(:,:,i),N*nadim,lu(1),iasai)
end do
do i=3,noab
  jk = 0
  do j=2,i-1
    ij = nTri_Elem(i-2)+j
    ik = nTri_Elem(i-2)
    do k=1,j-1
      jk = jk+1
      ik = ik+1
      !!write(u6,'(9I5)') i,j,k,iasai,iasaj,iasak,iasci,iascj,iasck
      ! K_ab^ir x L_rc^jk
      call DGEMM_('T','T',cdim,nadim,N,one,lc(:,jk),N,ka(:,:,i),nadim,zero,t3a,cdim)
      ! K_ab^kr x L_rc^ij
      call DGEMM_('T','T',cdim,nadim,N,one,lc(:,ij),N,ka(:,:,k),nadim,one,t3a,cdim)
      ! -K_ab^jr x L_rc^ik
      call DGEMM_('T','T',cdim,nadim,N,-one,lc(:,ik),N,ka(:,:,j),nadim,one,t3a,cdim)
      ! K_bc^ir x L_ra^jk
      call DGEMM_('N','N',ncdim,adim,N,one,kc(:,:,i),ncdim,la(:,jk),N,zero,t3b,ncdim)
      ! K_bc^kr x L_ra^ij
      call DGEMM_('N','N',ncdim,adim,N,one,kc(:,:,k),ncdim,la(:,ij),N,one,t3b,ncdim)
      ! -K_bc^jr x L_ra^ik
      call DGEMM_('N','N',ncdim,adim,N,-one,kc(:,:,j),ncdim,la(:,ik),N,one,t3b,ncdim)

      ab = 0
      do a=2,adim
        do b=1,a-1
          ab = ab+1
          t3a(:,ab) = t3a(:,ab)-t3b(:,a,b)+t3b(:,b,a)
        end do
      end do
      den = oeh(i)+oeh(j)+oeh(k)
      ab = 0
      !!kaka = 0
      do a=2,adim
        dena = den-oepa(a)
        do b=1,a-1
          ab = ab+1
          denb = dena-oepa(b)
          do c=1,cdim
            denc = denb-oepc(c)
            !!ab = ab+1
            !!kaka = kaka+1
            !!xx = t3a(ab,1)
            xx = t3a(c,ab)
            yy = xx/denc
            sumt3 = sumt3+xx
            enx = enx+yy*xx
            !!t3a(ab,1) = yy
            t3a(c,ab) = yy
            !!write(u6,*) ab,t3a(ab,1)
            !!write(u6,*) kaka,t3a(c,ab)
          end do
        end do
      end do
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
      ! ccsd(T) part vvoo*t3
      call DGEMM_('N','T',1,cdim,nadim,one,voa(:,ij),1,t3a,cdim,one,t1ac(k,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,one,voa(:,jk),1,t3a,cdim,one,t1ac(i,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,-one,voa(:,ik),1,t3a,cdim,one,t1ac(j,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,-one,voc(:,ij),1,t3b,ncdim,one,t1aa(k,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,-one,voc(:,jk),1,t3b,ncdim,one,t1aa(i,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,one,voc(:,ik),1,t3b,ncdim,one,t1aa(j,1),noab)
      ! ccsd(T) part t2*t3
      if (ifvo) then
        call DGEMM_('N','T',1,cdim,nadim,one,ka(:,i,j),1,t3a,cdim,one,t1bc(k,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,one,ka(:,j,k),1,t3a,cdim,one,t1bc(i,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,-one,ka(:,i,k),1,t3a,cdim,one,t1bc(j,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,-one,kc(:,i,j),1,t3b,ncdim,one,t1ba(k,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,-one,kc(:,j,k),1,t3b,ncdim,one,t1ba(i,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,one,kc(:,i,k),1,t3b,ncdim,one,t1ba(j,1),noab)
      end if
    end do !k
  end do !j
end do !i

return

end subroutine t3_bt_aac
