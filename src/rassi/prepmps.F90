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
! Copyright (C) 2016-2017, Stefan Knecht                               *
!***********************************************************************

!ifdef _DEBUGPRINT_
subroutine prepMPS(trorb,istate,lsym,mplet,mspro,nacte,tra,ntra,nish,nash,nosh,nsym,lupri,job,ist)
!-------------------------------------------------------------------------------
!
!    driver routine for the MPS rotation wrt the orbital transformation matrix
!    given in TRA.
!    --> rotate MPS state ISTATE
!    NOTE: TRA contains square matrices, one per symmetry
!
!-------------------------------------------------------------------------------

! module dependencies
#ifdef _DMRG_
use qcmaquis_info, only: qcm_group_names, qcm_prefixes
use qcmaquis_interface_cfg, only: dmrg_orbital_space, dmrg_state, dmrg_symmetry
use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_rotate
use fortran_strings, only: str
use Constants, only: Zero, One
use Definitions, only: u6
#endif
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: trorb
integer(kind=iwp), intent(in) :: istate, lsym, mplet, mspro, nacte, ntra, nsym, nish(nsym), nash(nsym), nosh(nsym), lupri, job, ist
real(kind=wp), intent(inout) :: tra(ntra)
#ifdef _DMRG_
integer(kind=iwp) :: i, ii, ista, isym, jorb, ni, no
real(kind=wp) :: fac
real(kind=wp), allocatable :: tmat(:,:)

! tmat: active-active rotation matrix
! Leon 8/12/2016 -- istatereal: "real" state index across all JobIphs, needed when we read checkpoint file names

if (.not. trorb) then
  write(lupri,'(a,a)') ' prepMPS: no MPS rotation requested for state ', &
                       trim(qcm_group_names(job)%states(ist))//' jobiph: '//str(job)//', root: '//str(istate)
  return
else
  write(lupri,'(a,a)') ' prepMPS:    MPS rotation requested for state ', &
                       trim(qcm_group_names(job)%states(ist))//' jobiph: '//str(job)//', root: '//str(istate)
end if

dmrg_orbital_space%nash(1:nsym) = nash(1:nsym)
dmrg_symmetry%nirrep = nsym
dmrg_state%nactel = nacte
dmrg_state%ms2 = mplet
dmrg_state%irefsm = lsym

#ifdef _DEBUGPRINT_
write(lupri,*) ' Entering prepMPS. TRA='
write(lupri,'(1x,5f16.8)') (TRA(I),I=1,NTRA)
#endif

!> find first the scaling factor to transform wrt the inactive orbitals
fac = One
ista = 1
do isym=1,nsym
  no = nosh(isym)
  do i=1,nish(isym)
    ii = ista+(no+1)*(i-1)
    fac = fac*tra(ii)
  end do
  ista = ista+no**2
end do
fac = fac**2
#ifdef _DEBUGPRINT_
write(lupri,*) ' scaling factor for MPS (inactive orbital rotations)',fac
#endif

dmrg_state%ms2 = mspro

if (nsym > 1) then
  write(u6,*) 'MPS rotation not supported with symmetry'
  call abend()
end if

!call mma_allocate(tmat,nash(1),nash(1))
! Leon: I get a maybe-uninitialized error if I use mma_allocate on tmat
allocate(tmat(nash(1),nash(1)))
tmat = Zero
ista = 1
jorb = 0
do isym=1,nsym
  ni = nish(isym)
  no = nosh(isym)
  do i=1,nash(isym)
    ! copy the active-active part of the rotation matrix into tmat
    tmat(i,:) = tra(1+(no+i)*ni+nash(isym)*(i-1):(no+i)*ni+nash(isym)*i)
  end do
  ista = ista+no**2
end do

! rotate MPS
call qcmaquis_mpssi_rotate(qcm_prefixes(job),istate,tmat,nash(1)**2,fac,mspro)

if (allocated(tmat)) deallocate(tmat)
!call mma_deallocate(tmat)

#else
write(lupri,*) ' calling prepMPS w/o DMRG interface - foolish!'
write(lupri,*) ' ... no actual task is performed.'
! Avoid unused variable warnings if DMRG is disabled
#include "macros.fh"
unused_var(trorb)
unused_var(istate)
unused_var(lsym)
unused_var(mplet)
unused_var(mspro)
unused_var(nacte)
unused_var(tra)
unused_var(nish)
unused_var(nash)
unused_var(nosh)
unused_var(job)
unused_var(ist)
#endif

end subroutine prepMPS
