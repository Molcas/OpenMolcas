!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine reprt2_cvb(orbs,cvb,civec,civb,civbs,civbh,citmp,sstruc,sstruc2,orbinv,sorbs,owrk,gjorb,gjorb2,gjorb3,cvbstot,cvbsspn, &
                      cvbdet,dvbdet,evbdet)

use casvb_global, only: formE, formroot, formSymW, formVBWnorm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), civec(ndet), civb(ndet), civbs(ndet), civbh(ndet), citmp(ndet), sstruc(nvb,nvb), &
                 sstruc2(nvb,nvb), orbinv(norb,norb), sorbs(norb,norb), owrk(norb,norb), gjorb(*), gjorb2(*), gjorb3(*), &
                 cvbstot(nvb), cvbsspn(nvb), cvbdet(ndetvb), dvbdet(ndetvb), evbdet(ndetvb)
#include "files_cvb.fh"
#include "print_cvb.fh"
integer(kind=iwp) :: i, ii, iimx, iorb, iroot, k, l, nr_print
real(kind=wp) :: cnrm, fac, occ_nel, rsum, sum1, sum2, swp
logical(kind=iwp) :: make_sstruc
real(kind=wp), allocatable :: dmat(:,:), occ(:)
real(kind=wp), external :: ddot_
logical(kind=iwp), external :: valid_cvb, ifcasci_cvb, ifhamil_cvb ! ... Files/Hamiltonian available ...

if (ip(5) >= 1) then
  call report_cvb(orbs,norb)
  write(u6,'(/,a)') ' Structure coefficients :'
  write(u6,'(a)') ' ------------------------'
  call vecprint_cvb(cvb,nvb)
end if

! First save CI vector
call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb)
call gaussj_cvb(orbs,gjorb)
call applyt_cvb(civb,gjorb)
call proj_cvb(civb)
call cinorm2_cvb(civb,cnrm)
call ciscale_cvb(civb,one/cnrm)
call symweight_cvb(civb,civb,wsym)

! Save VB wavefunction
call setsavvb_cvb(savvb)
if (valid_cvb(savvb)) then
  if (ip(5) >= 1) then
    write(u6,'(a)') ' '
    call prtfid_cvb(' Saving VB wavefunction to ',savvb)
  end if
  call putguess_cvb(orbs,cvb,savvb)
end if
call putci_cvb(civb)
call mma_allocate(dmat,norb,norb,label='dmat')
call fzero(dmat,norb*norb)
call onedens_cvb(civb,civb,dmat,.true.,0)
! Before overwriting CIVB/CIVEC get SVB:
if (lcalcsvb) then
  call cidot_cvb(civec,civb,svb)
  call untouch_cvb('SVB')
end if
if (lcalcevb) then
  call cicopy_cvb(civb,civbh)
  call applyh_cvb(civbh)
  call cidot_cvb(civb,civbh,evb)
  evb = evb+corenrg
  call untouch_cvb('EVB')
  ! ESYM needed for variational calculation
  call symweight_cvb(civb,civbh,esym)
  do i=1,nirrep
    if (wsym(i) > 1.0e-10_wp) then
      esym(i) = esym(i)/wsym(i)+corenrg
    else
      esym(i) = Zero
    end if
  end do
end if
if (lcalcsvb .and. (ip(5) >= 1)) then
  write(u6,'(a)') ' '
  write(u6,formE) ' Svb :      ',svb
  if ((lcalcevb .and. (ip(5) >= 1)) .or. ((icrit == 2) .and. (ip(3) < 0) .and. (ip(5) >= 1))) write(6,formE) ' Evb :      ',evb
else
  if ((lcalcevb .and. (ip(5) >= 1)) .or. ((icrit == 2) .and. (ip(3) < 0) .and. (ip(5) >= 1))) then
    write(u6,'(a)') ' '
    write(u6,formE) ' Evb :      ',evb
  end if
end if

if ((.not. variat) .or. endvar) then

  call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)
  call makecivbs_cvb(civbs,orbs,gjorb,gjorb2,gjorb3,cvbdet)

  call ciscale_cvb(civbs,one/cnrm)
  call dscal_(nvb,one/cnrm,cvb,1)
  call dscal_(ndetvb,one/cnrm,cvbdet,1)
  call finalresult_cvb()
  ! -- Analysis
  if (.not. lcalccivbs) then
    ! CIVBS has been evaluated previously (NB. based on unnormalized CVB):
    call ciscale_cvb(civbs,one/cnrm)
    call ci2vbg_cvb(civbs,dvbdet)
  else
    call cicopy_cvb(civb,civbs)
    call transp_cvb(orbs,owrk,norb,norb)
    call gaussj_cvb(owrk,gjorb2)
    call applyt_cvb(civbs,gjorb2)
    call ci2vbg_cvb(civbs,dvbdet)
  end if
  call vb2strg_cvb(dvbdet,cvbstot)
  sum1 = ddot_(nvb,cvb,1,cvbstot,1)
  call vb2strg_cvb(cvbdet,cvbsspn)
  sum2 = ddot_(nvb,cvb,1,cvbsspn,1)
  if (((ip(5) >= 1) .and. (ivbweights < 0)) .or. ((ivbweights >= 0) .and. (mod(ivbweights,2) == 1))) then
    write(u6,'(/,a)') ' Chirgwin-Coulson weights of structures :'
    write(u6,'(a)') ' ----------------------------------------'
    write(u6,formVBWnorm) ' VB spin+space (norm ',sum1,') :'
    do i=1,nvb
      dvbdet(i) = cvb(i)*cvbstot(i)/sum1
    end do
    call vecprint_cvb(dvbdet,nvb)
    if (sum2 < 1.0e3_wp) then
      write(u6,formVBWnorm) ' VB spin only  (norm ',sum2,') :'
    else
      formVBWnorm(4:4) = 'e'
      write(u6,formVBWnorm) ' VB spin only  (norm ',sum2,') :'
      formVBWnorm(4:4) = 'f'
    end if
    do i=1,nvb
      dvbdet(i) = cvb(i)*cvbsspn(i)/sum2
    end do
    call vecprint_cvb(dvbdet,nvb)
  end if
  make_sstruc = ((ivbweights > 1) .or. (ishstruc == 1))
  if (make_sstruc) then
    if (.not. proj) then
      call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
      call gaussj_cvb(sorbs,gjorb3)
    else
      call gaussj_cvb(orbs,gjorb)
      call transp_cvb(orbs,owrk,norb,norb)
      call gaussj_cvb(owrk,gjorb2)
    end if
    do k=1,nvb
      call fzero(sstruc(1,k),nvb)
      sstruc(k,k) = One
      call str2vbc_cvb(sstruc(1,k),dvbdet)
      call vb2cif_cvb(dvbdet,civb)
      if (.not. proj) then
        call applyt_cvb(civb,gjorb3)
      else
        call applyt_cvb(civb,gjorb)
        call proj_cvb(civb)
        call applyt_cvb(civb,gjorb2)
      end if
      call ci2vbg_cvb(civb,dvbdet)
      call vb2strg_cvb(dvbdet,sstruc(1,k))
    end do
    do k=1,nvb
      do l=k+1,nvb
        sstruc(k,l) = Half*(sstruc(k,l)+sstruc(l,k))
        sstruc(l,k) = sstruc(k,l)
      end do
    end do
  end if

  if (ivbweights > 1) then
    ! --  VB analysis Lowdin and Inverse  - begin  --
    if (mod(ivbweights,8) > 3) then
      write(u6,'(/,a)') ' Inverse-overlap weights of structures :'
      write(u6,'(a)') ' ---------------------------------------'
      call fmove_cvb(sstruc,sstruc2,nvb*nvb)
      call mxinv_cvb(sstruc2,nvb)
      ! Use DVBDET for weights:
      rsum = Zero
      do k=1,nvb
        dvbdet(k) = cvb(k)*cvb(k)/sstruc2(k,k)
        rsum = rsum+dvbdet(k)
      end do
      write(u6,formVBWnorm) ' VB spin+space (norm ',rsum,') :'
      call dscal_(nvb,one/rsum,dvbdet,1)
      call vecprint_cvb(dvbdet,nvb)
    end if

    if (mod(ivbweights,4) > 1) then
      write(u6,'(/,a)') ' Weights of Lowdin-orthogonalized structures :'
      write(u6,'(a)') ' ---------------------------------------------'
      call fmove_cvb(sstruc,sstruc2,nvb*nvb)
      ! Normalise overlap matrix before square root:
      ! Use CVBDET for normalized structure coefficients:
      call fmove_cvb(cvb,cvbdet,nvb)
      do k=1,nvb
        fac = sqrt(sstruc2(k,k))
        cvbdet(k) = fac*cvbdet(k)
        fac = one/fac
        call dscal_(nvb,fac,sstruc2(1,k),1)
        call dscal_(nvb,fac,sstruc2(k,1),nvb)
        sstruc2(k,k) = one
      end do

      call mxsqrt_cvb(sstruc2,nvb,1)
      ! Use DVBDET for weights:
      call mxatb_cvb(sstruc2,cvbdet,nvb,nvb,1,dvbdet)
      rsum = Zero
      do k=1,nvb
        dvbdet(k) = dvbdet(k)*dvbdet(k)
        rsum = rsum+dvbdet(k)
      end do
      write(u6,formVBWnorm) ' VB spin+space (norm ',rsum,') :'
      call dscal_(nvb,one/rsum,dvbdet,1)
      call vecprint_cvb(dvbdet,nvb)
    end if
  end if
  ! --  VB analysis Lowdin and Inverse  - end    --

  ! Spin correlation analysis
  ! Transform spin -> det to cater for non-orthogonal spin basis
  if (sij .and. sc) then
    call str2vbc_cvb(cvb,cvbdet)
    call str2vbg_cvb(cvbstot,dvbdet)
    call str2vbg_cvb(cvbsspn,evbdet)
    call scorr_cvb(cvbdet,dvbdet,evbdet)
  end if

  ! Weights of CASSCF vector in VB basis
  if (lciweights) then
    if (.not. ifcasci_cvb()) then
      write(u6,'(a)') ' Warning - no CIWEIGHTS without CASSCF vector.'
    else
      call str2vbc_cvb(cvb,cvbdet)
      call vb2cic_cvb(cvbdet,civb)
      call ciweight_cvb(civec,civbs,civb,citmp,civbh,orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3)
    end if
  end if

  if (ishstruc == 1) then
    write(u6,'(/,a)') ' Overlap matrix between structures :'
    write(u6,'(a)') ' -----------------------------------'
    call mxprintd_cvb(sstruc,nvb,nvb,0)
    if (ifhamil_cvb()) then
      call gaussj_cvb(orbs,gjorb)
      call transp_cvb(orbs,owrk,norb,norb)
      call gaussj_cvb(owrk,gjorb2)
      do k=1,nvb
        call fzero(sstruc2(1,k),nvb)
        sstruc2(k,k) = One
        call str2vbc_cvb(sstruc2(1,k),dvbdet)
        call vb2cif_cvb(dvbdet,civb)
        call applyt_cvb(civb,gjorb)
        call proj_cvb(civb)
        call applyh_cvb(civb)
        call proj_cvb(civb)
        call applyt_cvb(civb,gjorb2)
        call ci2vbg_cvb(civb,dvbdet)
        call vb2strg_cvb(dvbdet,sstruc2(1,k))
      end do
      do k=1,nvb
        do l=k+1,nvb
          sstruc2(k,l) = Half*(sstruc2(k,l)+sstruc2(l,k))
          sstruc2(l,k) = sstruc2(k,l)
        end do
      end do
      write(u6,'(/,a)') ' Hamiltonian matrix between structures :'
      write(u6,'(a)') ' ---------------------------------------'
      call mxprintd_cvb(sstruc2,nvb,nvb,0)
      call mxgendiag_cvb(sstruc2,sstruc,dvbdet,nvb)
      nr_print = min(20,nvb)
      if (nr_print < nvb) write(u6,'(/,a,i4,a)') ' Printing',nr_print,' lowest roots:'
      do iroot=1,nr_print
        write(u6,formroot) ' Root no.',iroot,' energy=',dvbdet(iroot),' :'
        call vecprint_cvb(sstruc2(1,iroot),nvb)
      end do
    end if
  end if

  call dscal_(nvb,cnrm,cvb,1)

  if ((ip(5) >= 1) .and. (nirrep > 1)) then
    write(u6,'(/,a)') ' Symmetry contributions to total VB wavefunction :'
    write(u6,'(a)') ' -------------------------------------------------'
    iimx = min(4,nirrep)
    write(u6,formSymW) ' Irreps 1 to',iimx,' :',(wsym(ii),ii=1,iimx)
    if (nirrep > 4) then
      iimx = min(8,nirrep)
      write(u6,formSymW) ' Irreps 5 to',iimx,' :',(wsym(ii),ii=5,iimx)
    end if
    if (lcalcevb) then
      write(u6,'(/,a)') ' Energies for components > 1e-10 :'
      write(u6,'(a)') ' ---------------------------------'
      iimx = min(4,nirrep)
      write(u6,formSymW) ' Irreps 1 to',iimx,' :',(esym(ii),ii=1,iimx)
      if (nirrep > 4) then
        iimx = min(8,nirrep)
        write(u6,formSymW) ' Irreps 5 to',iimx,' :',(esym(ii),ii=5,iimx)
      end if
    end if
  end if
  if (ip(5) >= 1) then
    write(u6,'(/,a)') ' One-electron density :'
    write(u6,'(a)') ' ----------------------'
    call mxprint_cvb(dmat,norb,norb,0)
    call mma_allocate(occ,norb,label='occ')
    call mxdiag_cvb(dmat,occ,norb)
    ! Sort NOs in order of increasing occ. numbers:
    do iorb=1,norb/2
      call dswap_(norb,dmat(1,iorb),1,dmat(1,norb+1-iorb),1)
      swp = occ(iorb)
      occ(iorb) = occ(norb+1-iorb)
      occ(norb+1-iorb) = swp
    end do
    write(u6,'(/,a)') ' Natural orbitals :'
    write(u6,'(a)') ' ------------------'
    call mxprint_cvb(dmat,norb,norb,0)
    write(u6,'(/,a)') ' Occupation numbers :'
    write(u6,'(a)') ' --------------------'
    call mxprint_cvb(occ,1,norb,0)
    occ_nel = Zero
    do i=1,norb
      occ_nel = occ_nel+occ(i)
    end do
    if (abs(occ_nel-real(nel,kind=wp)) > 0.1_wp) then
      write(u6,*) ' Error, sum of occupation numbers not equal to number of electrons :',occ_nel,nel
      call abend_cvb()
    end if
    call mma_deallocate(occ)
  end if
  call molden_cvb()
end if
call mma_deallocate(dmat)

return

end subroutine reprt2_cvb
