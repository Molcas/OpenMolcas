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

subroutine DO_SONATORB(NSS,USOR,USOI)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: EIGVEC, JBNUM
use rassi_data, only: NBTRI
use Cntrl, only: IfCurd, MLTPLT, NOSO, NSTATE, SODIAG, SODIAGNSTATE, SONAT, SONATNSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSS
real(kind=wp), intent(in) :: USOR(NSS,NSS), USOI(NSS,NSS)
integer(kind=iwp) :: I, IC, INATSTATE, IOPT, ISS, ISTATE, JOB1, JOB2, JSS, JSTATE, MPLET1, MPLET2, MSPROJ1, MSPROJ2
real(kind=wp) :: DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, IDENTMAT(3,3)
character(len=8) :: LAB
real(kind=wp), allocatable :: DMATTMP(:), UMATI(:), UMATR(:), VMAT(:,:)

! Calculates natural orbitals, including spinorbit effects
write(u6,*)
write(u6,*)
write(u6,*) '*****************************************'
write(u6,*) '* RUNNING SONATORB CODE *****************'
write(u6,*) '*****************************************'
write(u6,*)

call unitmat(IDENTMAT,3)

call mma_allocate(UMATR,NSS**2,Label='UMATR')
call mma_allocate(UMATI,NSS**2,Label='UMATI')
call mma_allocate(VMAT,NSS,NSS,Label='VMAT')

VMAT(:,:) = Zero

! transform V matrix in SF basis to spin basis
! This was taken from smmat and modified slightly
ISS = 0
do ISTATE=1,NSTATE
  JOB1 = JBNUM(ISTATE)
  MPLET1 = MLTPLT(JOB1)

  do MSPROJ1=-MPLET1+1,MPLET1-1,2
    ISS = ISS+1
    JSS = 0

    do JSTATE=1,NSTATE
      JOB2 = JBNUM(JSTATE)
      MPLET2 = MLTPLT(JOB2)

      do MSPROJ2=-MPLET2+1,MPLET2-1,2
        JSS = JSS+1

        if ((MPLET1 == MPLET2) .and. (MSPROJ1 == MSPROJ2)) VMAT(ISS,JSS) = EIGVEC(JSTATE,ISTATE)
      end do
    end do
  end do
end do

if (.not. NOSO) then
  ! combine this matrix with the SO eigenvector matrices
  call DGEMM_('N','N',NSS,NSS,NSS,One,VMAT,NSS,USOR,NSS,Zero,UMATR,NSS)
  call DGEMM_('N','N',NSS,NSS,NSS,One,VMAT,NSS,USOI,NSS,Zero,UMATI,NSS)
else
  ! Spinorbit contributions to this are disabled
  UMATR(:) = pack(VMAT(:,:),.true.)
  UMATI(:) = Zero
end if

! Holds the density matrices for all three directions
call mma_allocate(DMATTMP,6*NBTRI,Label='DMATTMP')

! SONATNSTATE = number of states to calculate.
! These states are stored beginning in SONAT
do I=1,SONATNSTATE
  INATSTATE = SONAT(I)

  write(u6,*)
  write(u6,*) 'CALCULATING NAT ORBITALS FOR SSTATE: ',INATSTATE
  if ((INATSTATE > NSS) .or. (INATSTATE <= 0)) then
    write(u6,*) '...WHICH DOES NOT EXIST!'
    call ABEND()
  end if
  write(u6,*)

  ! Calculate overall density, store in DMATTMP
  iOpt = 0
  call SONATORBM('HERMSING',UMATR,UMATI,INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

  ! Integrate for the expectation value
  if (IPGLOB >= 3) then
    IC = 1
    iOpt = 0
    LAB = 'MLTPL  0'
    call SONATORBM_INT(DMATTMP,LAB,IC,'HERMSING',INATSTATE,INATSTATE,iOpt,IDENTMAT,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

    !call ADD_INFO('MLTPL0SING_INT3',DUM3,1,6)
    !call ADD_INFO('MLTPL0SING_INT6',DUM6,1,6)
  end if

  ! Create SONATTDENS total density orbital file for this (I,I) state
  call SONATORB_PLOT(DMATTMP,'SONATTDENS','HERMSING',INATSTATE,INATSTATE)

  if (IPGLOB >= 4) call SONATORB_CPLOT(DMATTMP,'TDENSTESTX','HERMSING',INATSTATE,INATSTATE)

  ! Calculate spin density, store in LDMATTMP
  iOpt = 0
  call SONATORBM('HERMTRIP',UMATR,UMATI,INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

  ! Integrate for the expectation value
  if (IPGLOB >= 3) then
    IC = 1
    iOpt = 0
    LAB = 'MLTPL  0'
    call SONATORBM_INT(DMATTMP,LAB,IC,'HERMTRIP',INATSTATE,INATSTATE,iOpt,IDENTMAT,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
    !call ADD_INFO('MLTPL0TRIP_INT3',DUM3,1,6)
    !call ADD_INFO('MLTPL0TRIP_INT6',DUM6,1,6)
  end if

  ! Create SONATSDENS spin density orbital file for this (I,I) state
  call SONATORB_PLOT(DMATTMP,'SONATSDENS','HERMTRIP',INATSTATE,INATSTATE)

  if (IPGLOB >= 4) call SONATORB_CPLOT(DMATTMP,'SDENSTESTX','HERMTRIP',INATSTATE,INATSTATE)

  ! Type 2 - current density
  if (IFCURD) then
    iOpt = 0
    call SONATORBM('ANTISING',UMATR,UMATI,INATSTATE,INATSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

    if (IPGLOB >= 3) then
      IC = 1
      iOpt = 0
      LAB = 'ANGMOM'
      call SONATORBM_INT(DMATTMP,LAB,IC,'ANTISING',INATSTATE,INATSTATE,iOpt,IDENTMAT,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
      !call ADD_INFO('CURD1_INT3',DUM3,1,6)
      !call ADD_INFO('CURD1_INT6',DUM6,1,6)

      IC = 2
      iOpt = 0
      LAB = 'ANGMOM'
      call SONATORBM_INT(DMATTMP,LAB,IC,'ANTISING',INATSTATE,INATSTATE,iOpt,IDENTMAT,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)

      !call ADD_INFO('CURD2_INT3',DUM3,1,6)
      !call ADD_INFO('CURD2_INT6',DUM6,1,6)

      IC = 3
      iOpt = 0
      LAB = 'ANGMOM'
      call SONATORBM_INT(DMATTMP,LAB,IC,'ANTISING',INATSTATE,INATSTATE,iOpt,IDENTMAT,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6)
      !call ADD_INFO('CURD3_INT3',DUM3,1,6)
      !call ADD_INFO('CURD3_INT6',DUM6,1,6)

    end if

    call SONATORB_CPLOT(DMATTMP,'SONATLDENS','ANTISING',INATSTATE,INATSTATE)
  end if

end do

call mma_deallocate(DMATTMP)
call mma_deallocate(SONAT)

! perform the state diagonalization similar to
! what is done in single_aniso
if (SODIAGNSTATE > 0) then

  ! This actually does all the work
  call mkSODIAG(UMATR,UMATI,NSS)

  ! This is only allocated if SODIAGNSTATE > 0
  call mma_deallocate(SODIAG)
end if

call mma_deallocate(UMATR)
call mma_deallocate(UMATI)
call mma_deallocate(VMAT)

end subroutine DO_SONATORB
