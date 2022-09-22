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

subroutine t3_bt_aaa(nug,ka,la,adim,N,noab,nnoab,lu,iasblock,nga,oeh,oep,enx,voa,t1a,t1b,t3a,t3b,ifvo)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nug, adim, N, noab, nnoab, lu(2), iasblock(3), nga
real(kind=wp) :: ka(nTri_Elem(adim-1),N,*), la(N*adim,nnoab), oeh(noab), oep(adim), enx, voa(nTri_Elem(adim-1),nnoab), &
                 t1a(noab,*), t1b(noab,*), t3a(nTri_Elem(adim-1),adim), t3b(nTri_Elem(adim-1),adim)
logical(kind=iwp) :: ifvo
integer(kind=iwp) :: a, aa, ab, ac, b, bc, c, i, ias, ij, ik, j, jk, k, ndim, nga_offset, nug_offset
real(kind=wp) :: den, dena, denb, denc, XX, YY

if (adim == 1) return
!mp write(u6,*) 'enter_aaa enx,nga,',nga,enx
ndim = nTri_Elem(adim-1)
t3b(:,:) = Zero
nug_offset = iasblock(1)*nTri_Elem(nug)
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nnoab*adim*N,lu(2),ias)
ias = iasblock(2)*nug+iasblock(3)*(nTri_Elem(nga)-1)+1
call multi_readir(voa,nnoab*ndim,lu(2),ias)

nga_offset = iasblock(1)*(nTri_Elem(nga)-1)+1
do i=1,noab
  ias = (i-1)*nug_offset+nga_offset
  call multi_readir(ka(1,1,i),N*ndim,lu(1),ias)
end do
do i=3,noab
  jk = 0
  do j=2,i-1
    ij = nTri_Elem(i-2)+j
    ik = nTri_Elem(i-2)
    do k=1,j-1
      jk = jk+1
      ik = ik+1
      ! K_ab^ir x L_rc^jk
      call DGEMM_('N','N',ndim,adim,N,one,ka(1,1,i),ndim,la(1,jk),N,zero,t3a,ndim)
      ! K_ab^kr x L_rc^ij
      call DGEMM_('N','N',ndim,adim,N,one,ka(1,1,k),ndim,la(1,ij),N,one,t3a,ndim)
      ! -K_ab^jr x L_rc^ik
      call DGEMM_('N','N',ndim,adim,N,-one,ka(1,1,j),ndim,la(1,ik),N,one,t3a,ndim)
      den = oeh(i)+oeh(j)+oeh(k)
      do a=3,adim
        aa = nTri_Elem(a-2)
        dena = den-oep(a)
        bc = 0
        do b=2,a-1
          denb = dena-oep(b)
          ab = aa+b
          do c=1,b-1
            bc = bc+1
            denc = denb-oep(c)
            ac = aa+c
            xx = t3a(ab,c)+t3a(bc,a)-t3a(ac,b)
            yy = xx/denc
            enx = enx+yy*xx
            t3b(ab,c) = yy
            t3b(bc,a) = yy
            t3b(ac,b) = -yy
          end do
        end do
      end do
      ! ccsd(T) part vvoo*t3
      call DGEMM_('N','N',1,adim,ndim,one,voa(1,ij),1,t3b,ndim,one,t1a(k,1),noab)
      call DGEMM_('N','N',1,adim,ndim,one,voa(1,jk),1,t3b,ndim,one,t1a(i,1),noab)
      call DGEMM_('N','N',1,adim,ndim,-one,voa(1,ik),1,t3b,ndim,one,t1a(j,1),noab)
      ! ccsd(T) part t2*t3
      if (ifvo) then
        call DGEMM_('N','N',1,adim,ndim,one,ka(1,i,j),1,t3b,ndim,one,t1b(k,1),noab)
        call DGEMM_('N','N',1,adim,ndim,one,ka(1,j,k),1,t3b,ndim,one,t1b(i,1),noab)
        call DGEMM_('N','N',1,adim,ndim,-one,ka(1,i,k),1,t3b,ndim,one,t1b(j,1),noab)
      end if
    end do !k
  end do !j
end do !i

return

end subroutine t3_bt_aaa
