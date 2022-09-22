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

subroutine t3_bt_acc(nug,ka,kc,la,lc,adim,cdim,N,noab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,t1aa,t1ba,t1ac,t1bc,t3a, &
                     t3b,ifvo)

use Constants, only: Zero, One
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nug, adim, cdim, N, noab, nnoab, lu(2), iasblock(3), nga, ngc
real(kind=wp) :: ka(adim*cdim,N,*), kc(nTri_Elem(cdim-1),N,*), la(N*adim,nnoab), lc(N*cdim,nnoab), oeh(noab), oepa(adim), &
                 oepc(cdim), enx, voa(adim*cdim,*), voc(nTri_Elem(cdim-1),*), t1aa(noab,*), t1ba(noab,*), t1ac(noab,*), &
                 t1bc(noab,*), t3a(cdim,cdim,*), t3b(*)
logical(kind=iwp) :: ifvo
integer(kind=iwp) :: a, abc, b, bc, c, i, ias, iasai, iasci, ij, ik, j, jk, k, nadim, ncdim, nga_offset, ngc_offset, nug_offset
real(kind=wp) :: den, dena, denb, denc, sumt3, xx, yy

sumt3 = Zero
if (cdim == 1) return
ncdim = nTri_Elem(cdim-1)
nadim = adim*cdim
nug_offset = iasblock(1)*nTri_Elem(nug)
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nnoab*adim*N,lu(2),ias)
ias = iasblock(2)*(ngc-1)+1
call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
! reads vvoo
nga_offset = iasblock(3)*(nTri_Elem(nga-1)+ngc-1)+1
ngc_offset = iasblock(3)*(nTri_Elem(ngc-1)+ngc-1)+1
ias = iasblock(2)*nug+nga_offset
call multi_readir(voa,nnoab*nadim,lu(2),ias)
ias = iasblock(2)*nug+ngc_offset
call multi_readir(voc,nnoab*ncdim,lu(2),ias)

nga_offset = iasblock(1)*(nTri_Elem(nga-1)+ngc-1)+1
ngc_offset = iasblock(1)*(nTri_Elem(ngc-1)+ngc-1)+1
do i=1,noab
  iasci = (i-1)*nug_offset+ngc_offset
  call multi_readir(kc(1,1,i),N*ncdim,lu(1),iasci)
  iasai = (i-1)*nug_offset+nga_offset
  call multi_readir(ka(1,1,i),N*nadim,lu(1),iasai)
end do
do i=3,noab
  jk = 0
  do j=2,i-1
    ij = nTri_Elem(i-2)+j
    ik = nTri_Elem(i-2)
    do k=1,j-1
      jk = jk+1
      ik = ik+1
      ! K_ab^ir x L_rc^jk   (stored as   c,b,a)
      call DGEMM_('T','T',cdim,nadim,N,one,lc(1,jk),N,ka(1,1,i),nadim,zero,t3a,cdim)
      ! K_ab^kr x L_rc^ij
      call DGEMM_('T','T',cdim,nadim,N,one,lc(1,ij),N,ka(1,1,k),nadim,one,t3a,cdim)
      ! -K_ab^jr x L_rc^ik
      call DGEMM_('T','T',cdim,nadim,N,-one,lc(1,ik),N,ka(1,1,j),nadim,one,t3a,cdim)

      ! K_bc^ir x L_ra^jk
      call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,i),ncdim,la(1,jk),N,zero,t3b,ncdim)
      ! K_bc^kr x L_ra^ij
      call DGEMM_('N','N',ncdim,adim,N,one,kc(1,1,k),ncdim,la(1,ij),N,one,t3b,ncdim)
      ! -K_bc^jr x L_ra^ik
      call DGEMM_('N','N',ncdim,adim,N,-one,kc(1,1,j),ncdim,la(1,ik),N,one,t3b,ncdim)

      den = oeh(i)+oeh(j)+oeh(k)
      abc = 0
      do a=1,adim
        dena = den-oepa(a)
        bc = 0
        do b=2,cdim
          denb = dena-oepc(b)
          do c=1,b-1
            denc = denb-oepc(c)
            bc = bc+1
            abc = abc+1
            xx = t3b(abc)-t3a(b,c,a)+t3a(c,b,a)
            yy = xx/denc
            t3b(abc) = yy
            sumt3 = sumt3+xx
            enx = enx+yy*xx
          end do
        end do
      end do
      call expa1_uhf(t3b,adim,cdim,-1,t3a)
      ! t3a  bac
      ! blok1   voa ka
      call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,ij),1,t3a,cdim,one,t1ac(k,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,jk),1,t3a,cdim,one,t1ac(i,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,one,voa(1,ik),1,t3a,cdim,one,t1ac(j,1),noab)
      ! blok 2 t3b bca
      call DGEMM_('N','N',1,adim,ncdim,one,voc(1,ij),1,t3b,ncdim,one,t1aa(k,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,one,voc(1,jk),1,t3b,ncdim,one,t1aa(i,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,ik),1,t3b,ncdim,one,t1aa(j,1),noab)
      if (ifvo) then
        call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,i,j),1,t3a,cdim,one,t1bc(k,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,j,k),1,t3a,cdim,one,t1bc(i,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,one,ka(1,i,k),1,t3a,cdim,one,t1bc(j,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,one,kc(1,i,j),1,t3b,ncdim,one,t1ba(k,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,one,kc(1,j,k),1,t3b,ncdim,one,t1ba(i,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,i,k),1,t3b,ncdim,one,t1ba(j,1),noab)
      end if
    end do !k
  end do !j
end do !i

return

end subroutine t3_bt_acc
