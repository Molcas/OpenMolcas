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

subroutine HFCSD(LABEL,IC,BUFF,NBUFF,NSIZ,ISCHK)
!***********************************************************************
!     Objective: to compute the 'spin-dependent' part of the hyperfine *
!                from the magnetic integrals contributions             *
!     Output: BUFF                                                     *
!***********************************************************************

use hfc_logical, only: MAG_X2C
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, Four
use Definitions, only: wp, iwp, u6

implicit none
character(len=8), intent(inout) :: LABEL
integer(kind=iwp), intent(in) :: IC, NBUFF, NSIZ
real(kind=wp), intent(inout) :: BUFF(NBUFF)
integer(kind=iwp), intent(inout) :: ISCHK
integer(kind=iwp) :: ICM, IOPT, IRC
real(kind=wp) :: DA
real(kind=wp), allocatable :: TA(:)

! Set MAG_X2C to avoid add_info in hfcts
MAG_X2C = .true.
IOPT = 0
call mma_allocate(TA,NBUFF,Label='TA')
! BUFF needs to be initialized
BUFF(:) = Zero
! end of initialization
DA = Two

select case (IC)

  case (1)
    ! EF2(1) = (2*MAG(1)-MAG(5)-MAG(9))*(2/3)
    ICM = 1
    DA = Four/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

    ICM = 5
    DA = -Two/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

    ICM = 9
    DA = -Two/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

  case (2)
    ! EF2(2) = MAG(2)*2
    ICM = 2
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

  case (3)
    ! EF2(3) = MAG(3)*2
    ICM = 3
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

  case (4)
    ! EF2(4) = (2*MAG(5)-MAG(1)-MAG(9))*(2/3)
    ICM = 5
    DA = Four/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

    ICM = 1
    DA = -Two/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

    ICM = 9
    DA = -Two/Three
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

  case (5)
    ! EF2(5) = MAG(6)*2
    ICM = 6
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

  case (6)
    ! EF2(6) = (MAG(1)+MAG(5)+MAG(9))*2
    ICM = 1
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)
    BUFF(NSIZ+1:NBUFF) = TA(NSIZ+1:NBUFF)

    ICM = 5
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

    ICM = 9
    call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
    if (IRC /= 0) call ErrStop()
    call DAXPY_(NSIZ,DA,TA,1,BUFF,1)

  case default
    write(u6,'(6X,A)') '*** ERROR IN SUBROUTINE HFCSD ***'
    call Abend()

end select

call mma_deallocate(TA)

contains

subroutine ErrStop()

  write(u6,*)
  write(u6,'(6X,A)') '*** ERROR IN SUBROUTINE HFCSD ***'
  write(u6,'(6X,A)') '  FAILED IN READING FROM  ONEINT'
  write(u6,'(6X,A)') ' PLEASE MAKE SURE THE MAGNETIC'
  write(u6,'(6X,A)') ' HYPERFINE INTEGRALS ARE AVAILABLE'
  write(u6,'(6X,A,A)') '  LABEL     = ',LABEL
  write(u6,'(6X,A,I2)') '  COMPONENT = ',ICM
  write(u6,*)
  call ABEND()

end subroutine ErrSTop

end subroutine HFCSD
