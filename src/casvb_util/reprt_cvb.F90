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

subroutine reprt_cvb()

use casvb_global, only: civb1, civb2, civb3, civb4, civb5, corenrg, cvb, cvbdet, cvbsspn, cvbstot, dvbdet, endvar, esym, evb, &
                        evbdet, formE, formroot, formSymW, formVBWnorm, gjorb, gjorb2, gjorb3, icrit, ifhamil, ipr, ishstruc, &
                        ivbweights, lcalccivbs, lcalcevb, lcalcsvb, lciweights, mxirrep, nel, nirrep, norb, nvb, orbinv, orbs, &
                        owrk2, proj, savvb, sc, sij, sorbs, sstruc, sstruc2, svb, variat
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, ii, iimx, iorb, iroot, k, l, nr_print
real(kind=wp) :: cnrm, fac, occ_nel, rsum, sum1, sum2, swp, wsym(mxirrep)
logical(kind=iwp) :: make_sstruc
real(kind=wp), allocatable :: dmat(:,:), occ(:)
real(kind=wp), external :: ddot_
logical(kind=iwp), external :: valid_cvb, ifcasci_cvb ! ... Files available ...

if (ipr(5) >= 1) then
  call report_cvb(orbs,norb)
  write(u6,'(/,a)') ' Structure coefficients :'
  write(u6,'(a)') ' ------------------------'
  call vecprint_cvb(cvb,nvb)
end if

! First save CI vector
call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb2)
call gaussj_cvb(orbs,gjorb)
call applyt_cvb(civb2,gjorb)
call proj_cvb(civb2)
call cinorm2_cvb(civb2,cnrm)
call ciscale_cvb(civb2,One/cnrm)
call symweight_cvb(civb2,civb2,wsym)

! Save VB wavefunction
call setsavvb_cvb(savvb)
if (valid_cvb(savvb)) then
  if (ipr(5) >= 1) then
    write(u6,'(a)') ' '
    call prtfid_cvb(' Saving VB wavefunction to ',savvb)
  end if
  call putguess_cvb(orbs,cvb,savvb)
end if
call putci_cvb(civb2)
call mma_allocate(dmat,norb,norb,label='dmat')
dmat(:,:) = Zero
call onedens_cvb(civb2,civb2,dmat,.true.,0)
! Before overwriting CIVB/CIVEC get SVB:
if (lcalcsvb) then
  call cidot_cvb(civb1,civb2,svb)
  call untouch_cvb('SVB')
end if
if (lcalcevb) then
  call cicopy_cvb(civb2,civb4)
  call applyh_cvb(civb4)
  call cidot_cvb(civb2,civb4,evb)
  evb = evb+corenrg
  call untouch_cvb('EVB')
  ! ESYM needed for variational calculation
  call symweight_cvb(civb2,civb4,esym)
  do i=1,nirrep
    if (wsym(i) > 1.0e-10_wp) then
      esym(i) = esym(i)/wsym(i)+corenrg
    else
      esym(i) = Zero
    end if
  end do
end if
if (lcalcsvb .and. (ipr(5) >= 1)) then
  write(u6,'(a)') ' '
  write(u6,formE) ' Svb :      ',svb
  if ((lcalcevb .and. (ipr(5) >= 1)) .or. ((icrit == 2) .and. (ipr(3) < 0) .and. (ipr(5) >= 1))) write(6,formE) ' Evb :      ',evb
else
  if ((lcalcevb .and. (ipr(5) >= 1)) .or. ((icrit == 2) .and. (ipr(3) < 0) .and. (ipr(5) >= 1))) then
    write(u6,'(a)') ' '
    write(u6,formE) ' Evb :      ',evb
  end if
end if

if ((.not. variat) .or. endvar) then

  call makegjorbs_cvb(orbs)
  call makecivbs_cvb(civb3,orbs,cvbdet)

  call ciscale_cvb(civb3,One/cnrm)
  cvb(1:nvb) = cvb(1:nvb)/cnrm
  cvbdet(:) = cvbdet(:)/cnrm
  call finalresult_cvb()
  ! -- Analysis
  if (.not. lcalccivbs) then
    ! CIVBS has been evaluated previously (NB. based on unnormalized CVB):
    call ciscale_cvb(civb3,One/cnrm)
    call ci2vbg_cvb(civb3,dvbdet)
  else
    call cicopy_cvb(civb2,civb3)
    call trnsps(norb,norb,orbs,owrk2)
    call gaussj_cvb(owrk2,gjorb2)
    call applyt_cvb(civb3,gjorb2)
    call ci2vbg_cvb(civb3,dvbdet)
  end if
  call vb2strg_cvb(dvbdet,cvbstot)
  sum1 = ddot_(nvb,cvb,1,cvbstot,1)
  call vb2strg_cvb(cvbdet,cvbsspn)
  sum2 = ddot_(nvb,cvb,1,cvbsspn,1)
  if (((ipr(5) >= 1) .and. (ivbweights < 0)) .or. ((ivbweights >= 0) .and. (mod(ivbweights,2) == 1))) then
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
      call trnsps(norb,norb,orbs,owrk2)
      call gaussj_cvb(owrk2,gjorb2)
    end if
    do k=1,nvb
      sstruc(:,k) = Zero
      sstruc(k,k) = One
      call str2vbc_cvb(sstruc(1,k),dvbdet)
      call vb2cif_cvb(dvbdet,civb2)
      if (.not. proj) then
        call applyt_cvb(civb2,gjorb3)
      else
        call applyt_cvb(civb2,gjorb)
        call proj_cvb(civb2)
        call applyt_cvb(civb2,gjorb2)
      end if
      call ci2vbg_cvb(civb2,dvbdet)
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
      sstruc2(:,:) = sstruc(:,:)
      call mxinv_cvb(sstruc2,nvb)
      ! Use DVBDET for weights:
      rsum = Zero
      do k=1,nvb
        dvbdet(k) = cvb(k)*cvb(k)/sstruc2(k,k)
        rsum = rsum+dvbdet(k)
      end do
      write(u6,formVBWnorm) ' VB spin+space (norm ',rsum,') :'
      dvbdet(1:nvb) = dvbdet(1:nvb)/rsum
      call vecprint_cvb(dvbdet,nvb)
    end if

    if (mod(ivbweights,4) > 1) then
      write(u6,'(/,a)') ' Weights of Lowdin-orthogonalized structures :'
      write(u6,'(a)') ' ---------------------------------------------'
      sstruc2(:,:) = sstruc(:,:)
      ! Normalise overlap matrix before square root:
      ! Use CVBDET for normalized structure coefficients:
      cvbdet(1:nvb) = cvb(1:nvb)
      do k=1,nvb
        fac = sqrt(sstruc2(k,k))
        cvbdet(k) = fac*cvbdet(k)
        sstruc2(:,k) = sstruc2(:,k)/fac
        sstruc2(k,:) = sstruc2(k,:)/fac
        sstruc2(k,k) = One
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
      dvbdet(1:nvb) = dvbdet(1:nvb)/rsum
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
      call vb2cic_cvb(cvbdet,civb2)
      call ciweight_cvb(civb1,civb3,civb2,civb5,civb4,orbs,sorbs,orbinv,owrk2)
    end if
  end if

  if (ishstruc == 1) then
    write(u6,'(/,a)') ' Overlap matrix between structures :'
    write(u6,'(a)') ' -----------------------------------'
    call mxprintd_cvb(sstruc,nvb,nvb,0)
    if (ifhamil) then
      call gaussj_cvb(orbs,gjorb)
      call trnsps(norb,norb,orbs,owrk2)
      call gaussj_cvb(owrk2,gjorb2)
      do k=1,nvb
        sstruc2(:,k) = Zero
        sstruc2(k,k) = One
        call str2vbc_cvb(sstruc2(1,k),dvbdet)
        call vb2cif_cvb(dvbdet,civb2)
        call applyt_cvb(civb2,gjorb)
        call proj_cvb(civb2)
        call applyh_cvb(civb2)
        call proj_cvb(civb2)
        call applyt_cvb(civb2,gjorb2)
        call ci2vbg_cvb(civb2,dvbdet)
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

  cvb(1:nvb) = cnrm*cvb(1:nvb)

  if ((ipr(5) >= 1) .and. (nirrep > 1)) then
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
  if (ipr(5) >= 1) then
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
    occ_nel = sum(occ(:))
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

end subroutine reprt_cvb
