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
! Copyright (C) 2021, Rulin Feng                                       *
!***********************************************************************

!****************************************************
!                    Do SO-NTO
!****************************************************
! This routine is made to prepare the transition density
! This routine is modified from do_sonatorb.
! matrices, transition spin density matrices.
! In the future, this may include antisymmetric, transition
! densities. Namely the keyword antisin or antitrip
!
!                                               -RF 8/18,2021
subroutine DO_SONTO(NSS,USOR,USOI)

use rassi_global_arrays, only: JBNUM, EIGVEC
use cntrl, only: SONTOSTATES, SONTO
use Cntrl, only: NSTATE, NOSO, MLTPLT
use rassi_data, only: NBST
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer NSS
real*8 USOR(NSS,NSS), USOI(NSS,NSS)
real*8 IDENTMAT(3,3)
real*8, allocatable :: UMATR(:), UMATI(:), VMAT(:,:)
real*8, allocatable :: TDMAO(:), TSDMAO(:)
real*8, allocatable :: ANTSIN(:)
integer I, ISS, ISTATE, JOB1, MPLET1, MSPROJ1, JSS, JSTATE, JOB2, MPLET2, MSPROJ2, INTOSTATE, JNTOSTATE, IOPT

! Calculates natural orbitals, including spinorbit effects
write(u6,*)
write(u6,*)
write(u6,*) '*****************************************'
write(u6,*) '* RUNNING SONTO CODE ********************'
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
  !S1 = Half*real(MPLET1-1,kind=wp)

  do MSPROJ1=-MPLET1+1,MPLET1-1,2
    !SM1 = Half*real(MSPROJ1,kind=wp)
    ISS = ISS+1
    JSS = 0

    do JSTATE=1,NSTATE
      JOB2 = JBNUM(JSTATE)
      MPLET2 = MLTPLT(JOB2)
      !S2 = Half*real(MPLET2-1,kind=wp)

      do MSPROJ2=-MPLET2+1,MPLET2-1,2
        !SM2 = Half*real(MSPROJ2,kind=wp)
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
  call DCOPY_(NSS,VMAT,1,UMATR,1)
  call DCOPY_(NSS,[Zero],0,UMATI,1)
end if

! SONTONSTATE = number of state pairs to calculate.
! These states are stored as pairs beginning in SONTO
do I=1,SONTOSTATES
  INTOSTATE = SONTO(1,I)
  JNTOSTATE = SONTO(2,I)
  write(u6,*)
  write(u6,*) 'CALCULATING SO-NTOs BETWEEN SO STATES: ',INTOSTATE,JNTOSTATE
  if ((INTOSTATE > NSS) .or. (INTOSTATE <= 0) .or. (JNTOSTATE > NSS) .or. (JNTOSTATE <= 0)) then
    write(u6,*) '...WHICH DOES NOT EXIST!'
    call ABEND()
  end if
  write(u6,*)
  iOpt = 0
  ! Currently only HERMISING TDMs are dealt with here
  call mma_allocate(TDMAO,6*NBST**2,Label='TDMAO')
  call mma_allocate(TSDMAO,6*NBST**2,Label='TSDMAO')
  call mma_allocate(ANTSIN,6*NBST**2,Label='ANTSIN')
  ! Initialization is important
  TDMAO(:) = Zero
  TSDMAO(:) = Zero
  ANTSIN(:) = Zero

  call MAKETDMAO('HERMSING',UMATR,UMATI,INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,TDMAO)
  ! Following codes are left for other types of SO-TDMs
  !call print_matrixt('TDM after MAKETDMAO 1',nbst,nbst**2,1,TDMAO)
  !call MAKETDMAO('HERMTRIP',UMATR,UMATI,INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,TSDMAO,NBST)
  !call print_matrixt('TSDM after MAKETDMAO 1',nbst,nbst**2,1,TSDMAO)
  !call MAKETDMAO('ANTISING',UMATR,UMATI,INTOSTATE,JNTOSTATE,NSS,iOpt,IDENTMAT,ANTSIN,NBST)
  !call print_matrixt('ANTITDM after MAKETDMAO 1',nbst,nbst**2,1,ANTSIN)
  call DO_AOTDMNTO(TDMAO,TSDMAO,ANTSIN,INTOSTATE,JNTOSTATE,NBST,NBST**2)
  call mma_deallocate(TDMAO)
  call mma_deallocate(TSDMAO)
  call mma_deallocate(ANTSIN)
end do
call mma_deallocate(UMATR)
call mma_deallocate(UMATI)
call mma_deallocate(VMAT)
call mma_deallocate(SONTO)

end subroutine DO_SONTO
