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

subroutine cidenmat()

implicit none
integer :: i, ncount1, ncount2, neigen
real*8 :: sechc !, x1e(50000)
!character(len=256) :: filename
!logical :: logic_mulroot
real*8, parameter :: htoklm = 627.50956d0, zero = 0.0d0
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "grad_h.fh"
#include "files_gugaci.fh"
#include "scratch.fh"
#include "lgrn.fh"
#include "iaib.fh"
#include "vect.fh"
#include "grad_xyz.fh"
#include "ncprhf.fh"

!=======================================================================
! the main subroutine for ci gradient calculations.
! ican_a and ican_b save canonical order.
! e(*) save the scf orbital energies.
! xlgrn(*) save the lagrangian matrix.

!-----------------------------------------------------------------------
!  initialize the integral index nxo by canonical order
call init_canonical()

! calculate the frozen mo contribution to density matrix
!call density_scf_frz()
!-----------------------------------------------------------------------
! calculate the ci reduced one and two electron density matrix
ncount1 = ican_a(norb_all)+norb_all
ncount2 = ican_b(ncount1)+ncount1

!logic_mulroot = .false.
!ndim = mroot*nci_dim
!if (ndim <= max_vector) logic_mulroot = .true.

neigen = mroot
vector1(1:nci_dim) = zero

! calculated ci density matrix
! need debug here, density matrix is saved in vector2, it is too large t
! keep it.
! read ci vector
!if (logic_mulroot) then
!  call read_ml(lucivec,1,vector1,neigen*nci_dim,1)
!end if
!write(6,*) 'density matrix',neigen,ncount2
! should we calculated two-electronic density matrix ?
do i=1,neigen
  call read_ml(lucivec,1,vector1,nci_dim,i)
  vector2 = zero
  dm1tmp = 0.d0
  call memcidiag_alloc()
  call diagonal_loop_wyb_g()
  call memcidiag_dealloc()
  call matrix_vector_multi_parallel_drt_g(sechc)
  call ci_density_label_sm(i,ncount2)
  !call ci_dentest(i)
end do

return

end subroutine cidenmat

subroutine ci_dentest(iroot)

implicit none
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
#include "grad_h.fh"
integer :: iroot
integer :: i, idisk, idisk0, idisk_array(max_root+1), idx, iout, itratoc(ntratoc), kbuf, lenrd, li, lj, lk, ll, lri, lrj, nbpq, &
           nbrs, nc, nc1, nc2, nidx, nintb, nintone, nism, nmob, noidx(8), nop, noq, nor, norb(8), nos, nsmint, nsp, nspq, nspqr, &
           nsq, nsr, nss, nssm, ntj, ntk, numax, numin
real*8 :: buff(ntrabuf), cienergy, ecor, val, vnuc, xfock(max_orb*(max_orb+1)/2) !, x(1024*1024)
real*8, pointer :: x(:)
integer, external :: ipair

kbuf = ntrabuf

#ifndef MOLPRO
norb = nlsm_all
nintone = 0
nmob = 0
nidx = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  nsmint = nism*(nism+1)/2
  itratoc(i) = nism*nism
  nmob = nmob+nlsm_bas(i)*nlsm_bas(i)
  noidx(i) = nidx
  nidx = nidx+nism
  nintone = nintone+nsmint
end do

allocate(x(nmob))
call daname(luonemo,fnonemo)
!call copen_molcas(nft,filename,lenstr)
call readtraonehead(luonemo,ecor,idisk)
vnuc = ecor

! read mo coeff, need debug, if frozen and delete orbitals are not zero
!call ddatard(nft,x,nmob,idisk)
call ddafile(luonemo,2,x,nmob,idisk)

deallocate(x)

! read one electron fock matrix
call ddafile(luonemo,2,xfock,nintone,idisk)
!call ddatard(nft,xfock,nintone,idisk)
call daclos(luonemo)
! read one elctron kenetic intergrals
!call ddatard(nft,x1e,nintone,idisk)
!call cclose_molcas(nft)
! write one electron fock matrix into voint
cienergy = 0.d0
nidx = 0
nc = 0
do i=1,ng_sm
  nism = nlsm_all(i)
  idx = noidx(i)
  nsmint = nism*(nism+1)/2
  do lri=1,nism
    do lrj=1,lri
      nc = nc+1
      val = xfock(nc)*denm1(nc)
      if (lri /= lrj) val = 2.d0*val
      cienergy = cienergy+val
      write(6,'(1x,i3,1x,i3,2(1x,f18.9))') lri+idx,lrj+idx,xfock(nc),denm1(nc)
    end do
  end do
  nidx = nidx+nsmint
end do

idisk0 = 0
call idafile(luciden,2,idisk_array,max_root+1,idisk0)
idisk0 = idisk_array(iroot)
call ddafile(luciden,2,denm1,nc,idisk0)

call daname(lutwomo,fntwomo)

denm2 = 0.d0
idisk = 0
lenrd = ntratoc*lenintegral
write(6,*) lenrd
call idafile(lutwomo,2,itratoc,ntratoc,idisk)
write(6,2000)
2000 format(/7x,'symmetry',6x,' orbitals',8x,'integrals')
do nsp=1,ng_sm
  nop = norb(nsp)
  do nsq=1,nsp
    noq = norb(nsq)
    nspq = mul_tab(nsp,nsq)
    do nsr=1,nsp
      nor = norb(nsr)
      nspqr = mul_tab(nspq,nsr)
      nssm = nsr
      if (nsr == nsp) nssm = nsq
      do nss=1,nssm
        if (nspqr /= nss) cycle
        nos = norb(nss)

        if (nsr == nss) then
          nbpq = (nop+nop**2)/2
          nbrs = (nos+nos**2)/2
          if (nsp == nsr) then
            ! (ii|ii) type 1 int
            nintb = (nbpq+nbpq**2)/2
          else
            ! (ii|jj) type 3 int
            nintb = nbpq*nbrs
          end if
        else
          nbpq = nop*noq
          nbrs = nor*nos
          if (nsp == nsr) then
            ! (ij|ij) type 2 int
            nintb = (nbpq+nbpq**2)/2
          else
            ! (ij|kl) type 4 int
            nintb = nbpq*nbrs
          end if
        end if

        if (nintb == 0) cycle
        !write(6,2100) nsp,nsq,nsr,nss,nop,noq,nor,nos,nintb
        !2100 format(7x,4i2,1x,4i4,2x,3x,i9)

        iout = 0
        call ddafile(lutwomo,2,buff,kbuf,idisk)
        call ddafile(luciden,2,denm2,nintb,idisk0)
        idx = 0

        do li=1,nor
          ntj = nos
          if (nsr == nss) ntj = li
          do lj=1,ntj
            ntk = 1
            if (nsp == nsr) ntk = li
            do lk=ntk,nop
              numin = 1
              if ((nsp == nsr) .and. (lk == li)) numin = lj
              numax = noq
              if (nsp == nsq) numax = lk
              do ll=numin,numax
                iout = iout+1
                if (iout > kbuf) then
                  !call ddatard(nft,buff,kbuf,idisk)
                  call ddafile(lutwomo,2,buff,kbuf,idisk)
                  iout = 1
                end if
                idx = idx+1
                val = 0.5d0*buff(iout)*denm2(idx)
                if (li /= lj) val = 2.d0*val
                if (lk /= ll) val = 2.d0*val
                nc1 = ipair(li,lj)
                nc2 = ipair(lk,ll)
                if (nc1 /= nc2) val = 2.d0*val
                cienergy = cienergy+val
                !write(6,'(5(1x,i4),2(1x,f18.9))') li,lj,lk,ll,iout,buff(iout),denm2(idx)
              end do
            end do
          end do
        end do
        if (idx /= nintb) then
          write(6,*) 'in ci_dentest,count error'
          call abend()
        end if
      end do
    end do
  end do
end do

call daclos(lutwomo)
write(6,'(a11,3(2x,f18.9))') 'ci energy=',cienergy,vnuc,cienergy+vnuc

#endif

return

end subroutine ci_dentest

subroutine init_canonical()

implicit none
integer :: i, l1, l2
#include "drt_h.fh"
#include "iaib.fh"

!=======================================================================
! calculate the canonical order for index transform

l1 = max_orb
l2 = max_orb*(max_orb+1)/2
do i=1,l1
  ican_a(i) = (i*i-i)/2
end do
do i=1,l2
  ican_b(i) = (i*i-i)/2
end do

end subroutine init_canonical

subroutine density_scf_frz()

implicit none
integer :: i, j, k
real*8 :: val
real*8, parameter :: four = 4.0d0, one = 1.0d0, two = 2.0d0, zero = 0.0d0
#include "drt_h.fh"
#include "vect.fh"
#include "density.fh"

do i=1,naorbs
  do j=1,naorbs
    val = zero
    do k=1,norb_frz
      val = val+cf(i,k)*cf(j,k)
    end do
    p(i,j) = val
    !write(2,'(2i8,f18.10)') i,j,p(i,j)
  end do
end do

end subroutine density_scf_frz

subroutine matrix_vector_multi_parallel_drt_g(sechc)

implicit none
real*8 :: sechc
real*8 :: sc1, sc2
real*8, external :: c_time
#include "drt_h.fh"

write(6,*)
sc1 = c_time()
call ext_space_loop_g()
call inner_space_loop_g()
call vd_drt_ci_new()
call dv_drt_ci_new()
call dd_drt_ci_new()
call dt_drt_ci_new()
call ds_drt_ci_new()
call tv_drt_ci_new()
call td_drt_ci_new()
call tt_drt_ci_new()
call ts_drt_ci_new()
call sv_drt_ci_new()
call sd_drt_ci_new_den()
call st_drt_ci_new()
call ss_drt_ci_new()
sc2 = c_time()
sechc = sc2-sc1
!write(6,'(a42,2x,f10.2,2x,a1)') 'End of calculating density',' matrix, takes',sechc,'s'

return

end subroutine matrix_vector_multi_parallel_drt_g

subroutine matrix_vector_multi_parallel_prt_g(sechc)

implicit none
real*8 :: sechc
real*8 :: sc1, sc10, sc11, sc12, sc13, sc14, sc15, sc16, sc2, sc3, sc4, sc5, sc6, sc7, sc8, sc9
real*8, external :: c_time
#include "drt_h.fh"

write(6,*)
sc1 = c_time()
call ext_space_loop_g()
sc2 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of ext, takes',sc2-sc1,'s'
call inner_space_loop_g()
sc3 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of inn, takes',sc3-sc2,'s'
call vd_drt_ci_new()
sc4 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of vd , takes',sc4-sc3,'s'
call dv_drt_ci_new()
sc5 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of dv , takes',sc5-sc4,'s'
call dd_drt_ci_new()
sc6 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of dd , takes',sc6-sc5,'s'
call dt_drt_ci_new()
sc7 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of dt , takes',sc7-sc6,'s'
call ds_drt_ci_new()
sc8 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of ds , takes',sc8-sc7,'s'
call tv_drt_ci_new()
sc9 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of tv , takes',sc9-sc8,'s'
call td_drt_ci_new()
sc10 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of td , takes',sc10-sc9,'s'
call tt_drt_ci_new()
sc11 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of tt , takes',sc11-sc10,'s'
call ts_drt_ci_new()
sc12 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of ts , takes',sc12-sc11,'s'
call sv_drt_ci_new()
sc13 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of sv , takes',sc13-sc12,'s'
call sd_drt_ci_new_den()
sc14 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of sd , takes',sc14-sc13,'s'
call st_drt_ci_new()
sc15 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of st , takes',sc15-sc14,'s'
call ss_drt_ci_new()
sc16 = c_time()
write(6,'(a16,2x,f10.2,2x,a1)') 'end of ss , takes',sc16-sc15,'s'
sc2 = c_time()
sechc = sc2-sc1
write(6,'(a16,2x,f10.2,2x,a1)') 'end of run, takes',sechc,'s'

return

end subroutine matrix_vector_multi_parallel_prt_g

subroutine ci_density_label_sm(iroot,ncount2)

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "files_gugaci.fh"
#include "iaib.fh"
integer :: iroot, ncount2
integer :: i, ic, idisk, ii, ij, ijm, im, indx_m(maxgdm), idisk_array(max_root+1), j, jc, je, jj, jm, k, kk, kl, klm, km, l, lc, &
           le, ll, lm, nc0, nc1
real*8 :: val
real*8, parameter :: half = 0.5d0, two = 2.0d0, zero = 0.0d0
integer, external :: ipair

!=======================================================================
! transfer the ci density matrices (dm1 and dm2)
! based on the gamess save rule.
!vector1(1:ncount2) = zero
!vector1(1:ncount2) = vector2(1:ncount2)
!vector2(1:ncount2) = zero

denm1 = 0.0d0
nc0 = 1
nc1 = 0
do im=1,ng_sm
  if (nlsm_all(im) == 0) cycle
  do i=1,nlsm_all(im)
    do j=1,i
      ii = map_orb_order(i+nc1)
      jj = map_orb_order(j+nc1)
      ij = ipair(ii,jj)
      denm1(nc0) = dm1tmp(ij)
      nc0 = nc0+1
    end do
  end do
  indx_m(im) = nc1
  nc1 = nc1+nlsm_all(im)
end do
idisk = 0
if (iroot == 1) then
  idisk_array = 0
  call idafile(luciden,1,idisk_array,max_root+1,idisk)
  idisk_array(1) = idisk
else
  idisk = 0
  call idafile(luciden,2,idisk_array,max_root+1,idisk)
  idisk = idisk_array(iroot)
end if
call ddafile(luciden,1,denm1,nc0,idisk)

denm2 = 0.d0
! label 2-elc den matrix
! cycle on symmetry
do im=1,ng_sm
  if (nlsm_all(im) == 0) cycle
  do jm=1,im
    if (nlsm_all(jm) == 0) cycle
    ijm = mul_tab(im,jm)
    do km=1,im
      if (nlsm_all(km) == 0) cycle
      le = km
      if (km == im) le = jm
      do lm=1,le
        if (nlsm_all(lm) == 0) cycle
        klm = mul_tab(km,lm)
        if (ijm /= klm) cycle
        ! ityp 1 (ii|jj) 2 (ii|jj) 3 (ij|ij) 4 (ij|kl)

        ! cycle on mo index
        ! i>=j, k>=l, ij>=kl
        nc0 = 0
        do k=1,nlsm_all(km)
          kk = map_orb_order(k+indx_m(km))
          lc = nlsm_all(lm)
          if (km == lm) lc = k
          do l=1,lc
            ll = map_orb_order(l+indx_m(lm))
            kl = ipair(kk,ll)
            ic = 1
            if (im == km) ic = k
            do i=ic,nlsm_all(im)
              ii = map_orb_order(i+indx_m(im))
              jc = 1
              if ((im == km) .and. (i == k)) jc = l
              je = nlsm_all(jm)
              if (im == jm) je = i
              do j=jc,je
                jj = map_orb_order(j+indx_m(jm))
                nc0 = nc0+1
                ij = ipair(ii,jj)
                nc1 = ipair(ij,kl)

                val = vector2(nc1)*half
                if ((ii /= jj) .and. (kk /= ll)) then
                  val = val*half
                else if ((ii == jj) .and. (kk == ll) .and. (ii == kk)) then
                  val = val*two
                end if
                denm2(nc0) = val
              end do
            end do
          end do
        end do
        call ddafile(luciden,1,denm2,nc0,idisk)
        ! end cycle on mo index
      end do
    end do
  end do
end do

idisk_array(iroot+1) = idisk
idisk = 0
call idafile(luciden,1,idisk_array,max_root+1,idisk)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(ncount2)

end subroutine ci_density_label_sm
