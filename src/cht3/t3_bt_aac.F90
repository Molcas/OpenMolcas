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

subroutine t3_bt_aac(nug,ka,kb,kc,la,lc,mi,mij,adim,cdim,N,noab,nuab,nnoab,lu,iasblock,nga,ngc,oeh,oepa,oepc,enx,voa,voc,t1aa, &
                     t1ba,t1ac,t1bc,t3a,t3b,ifvo)

implicit none
real*8 one, zero, sumt3
logical ifvo
parameter(one=1.d0,zero=0.d0)
integer nadim, adim, noab, nuab, i, j, k, nga, iasblock(3), lu(2), N
integer ias, nga_offset, ngc_offset, nug_offset, jk, ij, ik, ncdim, cdim
integer a, b, c, ab, nug, nnoab, ngc, iasci, iasai !, kaka
real*8 dena, denb, denc, xx, yy
real*8 ka(adim*(adim-1)/2,n,*), kc(adim*cdim,n,*), kb(*), la(N*adim,nnoab), lc(N*cdim,nnoab)
!!real*8 lb(N,adim,nnoab),lc(N,adim,nnoab)
real*8 t3a(cdim,*), t3b(cdim,adim,*), t1ac(noab,*), t1bc(noab,*)
real*8 mi(*), mij(*), voa(adim*(adim-1)/2,*), voc(adim*cdim,*)
real*8 t1aa(noab,*), t1ba(noab,*), enx, oeh(noab), oepa(adim), den
real*8 oepc(cdim)

sumt3 = 0.d0
if (adim == 1) return
nadim = adim*(adim-1)/2
ncdim = adim*cdim
nug_offset = iasblock(1)*nug*(nug+1)/2
ias = iasblock(2)*(nga-1)+1
call multi_readir(la,nnoab*adim*N,lu(2),ias)
ias = iasblock(2)*(ngc-1)+1
call multi_readir(lc,nnoab*cdim*N,lu(2),ias)
! reads vvoo
nga_offset = iasblock(3)*(nga*(nga-1)/2+nga-1)+1
ngc_offset = iasblock(3)*(nga*(nga-1)/2+ngc-1)+1
ias = iasblock(2)*nug+nga_offset
call multi_readir(voa,nnoab*nadim,lu(2),ias)
ias = iasblock(2)*nug+ngc_offset
call multi_readir(voc,nnoab*ncdim,lu(2),ias)

nga_offset = iasblock(1)*(nga*(nga-1)/2+nga-1)+1
ngc_offset = iasblock(1)*(nga*(nga-1)/2+ngc-1)+1
do i=1,noab
  iasci = (i-1)*nug_offset+ngc_offset
  call multi_readir(kc(1,1,i),N*ncdim,lu(1),iasci)
  iasai = (i-1)*nug_offset+nga_offset
  call multi_readir(ka(1,1,i),N*nadim,lu(1),iasai)
end do
do i=3,noab
  jk = 0
  do j=2,i-1
    ij = (i-1)*(i-2)/2+j
    ik = (i-1)*(i-2)/2
    do k=1,j-1
      jk = jk+1
      ik = ik+1
      !!write(6,'(9I5)') i,j,k,iasai,iasaj,iasak,iasci,iascj,iasck
      ! K_ab^ir x L_rc^jk
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

      ab = 0
      do a=2,adim
        do b=1,a-1
          ab = ab+1
          call daxpy_(cdim,-1.d0,t3b(1,a,b),1,t3a(1,ab),1)
          call daxpy_(cdim,1.d0,t3b(1,b,a),1,t3a(1,ab),1)
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
            !!write(6,*) ab,t3a(ab,1)
            !!write(6,*) kaka,t3a(c,ab)
          end do
        end do
      end do
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
      ! ccsd(T) part vvoo*t3
      call DGEMM_('N','T',1,cdim,nadim,one,voa(1,ij),1,t3a,cdim,one,t1ac(k,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,one,voa(1,jk),1,t3a,cdim,one,t1ac(i,1),noab)
      call DGEMM_('N','T',1,cdim,nadim,-one,voa(1,ik),1,t3a,cdim,one,t1ac(j,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,ij),1,t3b,ncdim,one,t1aa(k,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,-one,voc(1,jk),1,t3b,ncdim,one,t1aa(i,1),noab)
      call DGEMM_('N','N',1,adim,ncdim,one,voc(1,ik),1,t3b,ncdim,one,t1aa(j,1),noab)
      ! ccsd(T) part t2*t3
      if (ifvo) then
        call DGEMM_('N','T',1,cdim,nadim,one,ka(1,i,j),1,t3a,cdim,one,t1bc(k,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,one,ka(1,j,k),1,t3a,cdim,one,t1bc(i,1),noab)
        call DGEMM_('N','T',1,cdim,nadim,-one,ka(1,i,k),1,t3a,cdim,one,t1bc(j,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,i,j),1,t3b,ncdim,one,t1ba(k,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,-one,kc(1,j,k),1,t3b,ncdim,one,t1ba(i,1),noab)
        call DGEMM_('N','N',1,adim,ncdim,one,kc(1,i,k),1,t3b,ncdim,one,t1ba(j,1),noab)
      end if
    end do !k
  end do !j
end do !i

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(kb)
  call Unused_real_array(mi)
  call Unused_real_array(mij)
  call Unused_integer(nuab)
end if

end subroutine t3_bt_aac
