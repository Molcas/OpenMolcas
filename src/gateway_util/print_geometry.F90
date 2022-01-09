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

use Basis_Info
use Center_Info
use Period
use Temporary_Parameters, only: Expert
use Sizes_of_Seward, only: S
use Real_Info, only: Rtrnc
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
character(LEN=1) help_c
character(LEN=16) FMT
character(LEN=LENIN), allocatable :: Lblxxx(:)
real*8, dimension(:,:), allocatable :: Centr
#include "angstr.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
LuWr = 6
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Centr,3,S%mCentr,Label='Centr')
call mma_allocate(Lblxxx,S%mCentr,Label='Lblxxx')
!                                                                      *
!***********************************************************************
!                                                                      *

write(LuWr,*)
call CollapseOutput(1,'   Molecular structure info:')
write(LuWr,'(3X,A)') '   -------------------------'
write(LuWr,*)
if (iOpt == 0) FMT = '(19X,A)'
if (iOpt == 1) FMT = '(11X,A)'

write(LuWr,FMT) ' ************************************************ '
if (iOpt == 0) then
  write(LuWr,FMT) ' **** Cartesian Coordinates / Bohr, Angstrom **** '
else
  write(LuWr,FMT) ' **** Cartesian Coordinates / Angstrom       **** '
end if
write(LuWr,FMT) ' ************************************************ '
write(LuWr,*)
if (iOpt == 0) then
  write(LuWr,'(A)') '     Center  Label                x              y              z'// &
                                 '                     x              y              z'
else
  write(LuWr,'(A)') '     Center  Label                x              y              z'
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
    Go To 32
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
          write(LuWr,'(6X,I3,A1,5X,A,3F15.6,7X,3F15.6)') nc,help_c,dc(ndc)%LblCnt,Centr(1:3,nc),Centr(1:3,nc)*angstr
        else
          write(LuWr,'(6X,I3,A1,5X,A,3F15.6)') nc,help_c,dc(ndc)%LblCnt,Centr(1:3,nc)*angstr
        end if
      end if
      if (nc > 8*MxAtom) then
        call WarningMessage(2,'lblxxx too small')
        call Abend()
      end if
      lblxxx(nc) = dc(ndc)%LblCnt(1:LENIN)
      nc = nc+1
    end do
  end do
32 continue
end do
nc = nc-1
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute distances

if (S%mCentr <= 2) Go To 55
call Dstncs(lblxxx,Centr,nc,angstr,S%Max_Center,6)
if (.not. Expert) call DstChk(Centr,lblxxx,nc)

! Compute valence bond angles

if ((iPrint < 5) .or. (S%mCentr < 3) .or. (iOpt == 1)) Go To 55
call Angles(lblxxx,Centr,nc,rtrnc,S%Max_Center)

! Compute dihedral angles

if ((iPrint < 5) .or. (S%mCentr < 4)) Go To 55
call Dihedr(lblxxx,Centr,nc,rtrnc,S%Max_Center)
!                                                                      *
!***********************************************************************
!                                                                      *
55 continue

call mma_deallocate(Lblxxx)
call mma_deallocate(Centr)
!                                                                      *
!***********************************************************************
!                                                                      *
call CollapseOutput(0,'   Molecular structure info:')
write(LuWr,*)

return

end subroutine Print_Geometry
