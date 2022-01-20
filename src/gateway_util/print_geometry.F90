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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine Print_Geometry(iOpt)
!***********************************************************************
!                                                                      *
! Object: to print the molecular coordinates, bonds, angles and        *
!         torsional angles.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
!             September 2006                                           *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Period, only: AdCell, Cell_l, lthCell
use Gateway_global, only: Expert
use Sizes_of_Seward, only: S
use Gateway_Info, only: Rtrnc
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iOpt
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: i, iPrint, iRout, j, jCnt, jCnttp, mCnt, nc, ndc
character(len=16) :: frmt
character :: help_c
real(kind=wp), allocatable :: Centr(:,:)
character(len=LenIn), allocatable :: Lblxxx(:)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Centr,3,S%mCentr,Label='Centr')
call mma_allocate(Lblxxx,S%mCentr,Label='Lblxxx')
!                                                                      *
!***********************************************************************
!                                                                      *

write(u6,*)
call CollapseOutput(1,'   Molecular structure info:')
write(u6,'(3X,A)') '   -------------------------'
write(u6,*)
if (iOpt == 0) frmt = '(19X,A)'
if (iOpt == 1) frmt = '(11X,A)'

write(u6,frmt) ' ************************************************ '
if (iOpt == 0) then
  write(u6,frmt) ' **** Cartesian Coordinates / Bohr, Angstrom **** '
else
  write(u6,frmt) ' **** Cartesian Coordinates / Angstrom       **** '
end if
write(u6,frmt) ' ************************************************ '
write(u6,*)
if (iOpt == 0) then
  write(u6,'(A)') '     Center  Label                x              y              z'// &
                               '                     x              y              z'
else
  write(u6,'(A)') '     Center  Label                x              y              z'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
ndc = 0
nc = 1
do jCnttp=1,nCnttp
  mCnt = dbsc(jCnttp)%nCntr
  if (dbsc(jCnttp)%Aux .or. dbsc(jCnttp)%Frag) then
    ndc = ndc+mCnt
    cycle
  end if
  do jCnt=1,mCnt
    ndc = ndc+1
    do i=0,nIrrep/dc(ndc)%nStab-1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Centr(1:3,nc))
      if (Show) then
        help_c = ' '
        if (Cell_l) then
          do j=1,lthCell
            if (AdCell(j) == ndc) help_c = '*'
          end do
        end if
        if (iOpt == 0) then
          write(u6,'(6X,I3,A1,5X,A,3F15.6,7X,3F15.6)') nc,help_c,dc(ndc)%LblCnt,Centr(1:3,nc),Centr(1:3,nc)*Angstrom
        else
          write(u6,'(6X,I3,A1,5X,A,3F15.6)') nc,help_c,dc(ndc)%LblCnt,Centr(1:3,nc)*Angstrom
        end if
      end if
      if (nc > 8*MxAtom) then
        call WarningMessage(2,'lblxxx too small')
        call Abend()
      end if
      lblxxx(nc) = dc(ndc)%LblCnt(1:LenIn)
      nc = nc+1
    end do
  end do
end do
nc = nc-1
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute distances

if (S%mCentr > 2) then
  call Dstncs(lblxxx,Centr,nc,Angstrom,S%Max_Center,6)
  if (.not. Expert) call DstChk(Centr,lblxxx,nc)

  if (iPrint >= 5) then
    ! Compute valence bond angles

    if ((S%mCentr >= 3) .and. (iOpt /= 1)) then
      call Angles(lblxxx,Centr,nc,rtrnc,S%Max_Center)
    end if

    ! Compute dihedral angles

    if (S%mCentr >= 4) then
      call Dihedr(lblxxx,Centr,nc,rtrnc,S%Max_Center)
    end if
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *

call mma_deallocate(Lblxxx)
call mma_deallocate(Centr)
!                                                                      *
!***********************************************************************
!                                                                      *
call CollapseOutput(0,'   Molecular structure info:')
write(u6,*)

return

end subroutine Print_Geometry
