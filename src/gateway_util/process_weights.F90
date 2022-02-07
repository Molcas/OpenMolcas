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
! Copyright (C) 2013, Ignacio Fdez. Galvan                             *
!               2016, Roland Lindh                                     *
!***********************************************************************
!  Process_Weights
!
!> @brief
!>   Process and store the weights used for alignment
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Process the weights according to the `WEIG` option in `GATEWAY`.
!> These weights are used in alignment and in the "sphere" constraint.
!> The list of weights is stored in the runfile, first the symmetry-unique
!> atoms and then the symmetric images, in the manner of ::Expand_Coor.
!>
!> @param[in] iPrint Print level
!***********************************************************************

subroutine Process_Weights(iPrint)

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, UtoAU
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint
integer(kind=iwp) :: i, iAt, iCnt, iErr, iSymAt, j, k, nAt, ndc, nSymAt
real(kind=wp) :: wTot
logical(kind=iwp) :: Small
character(len=512) :: Align_Weights
real(kind=wp), allocatable :: W(:)
real(kind=wp), parameter :: thr = 1.0e-6_wp

call Get_cArray('Align_Weights',Align_Weights,512)

! Count the total and symmetry-unique number of atoms
nAt = 0
nSymAt = 0
ndc = 0
do i=1,nCnttp
  do j=1,dbsc(i)%nCntr
    ndc = ndc+1
    if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
      nAt = nAt+nIrrep/dc(ndc)%nStab
      nSymAt = nSymAt+1
    end if
  end do
end do
call mma_allocate(W,nAt,label='W')

! By default, all weights are 1
W(:) = One

if (Align_Weights(1:4) == 'MASS') then
  ! Set the weights to the mass of each atom
  j = 1
  do i=1,nCnttp
    if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
      do iCnt=1,dbsc(i)%nCntr
        W(j) = dbsc(i)%CntMass/UtoAU
        j = j+1
      end do
    end if
  end do
else if (Align_Weights(1:5) == 'HEAVY') then
  ! Set the the weight to 1 for heavy atoms, 0 for hydrogens
  j = 1
  do i=1,nCnttp
    if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
      do iCnt=1,dbsc(i)%nCntr
        if (dbsc(i)%AtmNr <= 1) W(j) = Zero
        j = j+1
      end do
    end if
  end do
else if (Align_Weights(1:5) == 'EQUAL') then
  ! EQUAL is already the default: 1 for all
  !continue
else
  ! Read the weights from the input line
  read(Align_Weights,*,iostat=iErr) (W(i),i=1,nAt)
  if (iErr > 0) then
    call WarningMessage(2,'Unable to read data from WEIG')
    call Quit_OnUserError()
  end if
end if

! Unfold the symmetry
iSymAt = 1
iAt = 1+nSymAt
ndc = 0
do i=1,nCnttp
  do j=1,dbsc(i)%ncntr
    ndc = ndc+1
    if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
      do k=1,nIrrep/dc(ndc)%nStab-1
        W(iAt) = W(iSymAt)
        iAt = iAt+1
      end do
      iSymAt = iSymAt+1
    end if
  end do
end do

! Check for zero total weight
wTot = Zero
do i=1,nAt
  wTot = wTot+W(i)
end do
if (wTot < thr) then
  call WarningMessage(1,'Total weight too small. Setting equal weights.')
  do i=1,nAt
    W(i) = One
  end do
end if
! Prevent zero weights, it could break the "sphere" constraint
! (a value between 1e-6 and 1e-1 can still be used)
Small = .false.
do i=1,nAt
  if (W(i) < thr) then
    W(i) = 1.0e-1_wp
    Small = .true.
  end if
end do
if (iPrint >= 6) then
  if (Small) then
    call WarningMessage(1,'Small weights were increased to avoid problems with constraints.')
  end if
  call RecPrt('Weights used for alignment and distance',' ',W,nAt,1)
  write(u6,*)
end if

! Store weights in the runfile too
call Put_dArray('Weights',W,nAt)
call mma_deallocate(W)

end subroutine Process_Weights
