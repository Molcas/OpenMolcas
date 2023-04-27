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

subroutine Print_T_Values(T_Values,iT_Sets,iANr,EC,Bond_Threshold,nAtoms,nij,Standard,iWarnings,Num_Warnings,iPrint)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nij, iT_Sets(nij), nAtoms, iANr(nAtoms), iWarnings(nij), Num_Warnings, iPrint
real(kind=wp), intent(in) :: T_Values(nij), EC(3,nij), Bond_Threshold
logical(kind=iwp), intent(in) :: Standard
#include "Molcas.fh"
integer(kind=iwp) :: i, iAtom, ii, ij, j, jAtom, jj, Last_NonBlank
real(kind=wp) :: Bond_Length, Bond_Max, bs_t, Factor, Radius_i, Radius_j
integer(kind=iwp), parameter :: iLength = 25
character(len=iLength) :: Warning
character(len=17) :: BondLbl
character(len=LenIn), allocatable :: AtomLbl(:)
character(len=LenIn4), allocatable :: AtomLbl4(:)
real(kind=wp), external :: Bragg_Slater

! Print header

call mma_allocate(AtomLbl,nAtoms,label='AtomLbl')
call mma_allocate(AtomLbl4,nAtoms,label='AtomLbl4')

call Get_cArray('LP_L',AtomLbl4,LenIn4*nAtoms)
do i=1,nAtoms
  AtomLbl(i) = AtomLbl4(i)(1:LenIn)
end do
call mma_deallocate(AtomLbl4)
ij = 0
write(u6,*)
if (Num_Warnings > 0) then
  write(u6,'(A,I3,A)') 'During optimization of the expansion centers ',Num_Warnings,' warnings were encountered.'
  write(u6,*)
  write(u6,'(A)') ' iAtom   jAtom   Atom(s)          Factor   Bragg-Slater      t      Warning'
else
  write(u6,'(A)') ' iAtom   jAtom   Atom(s)          Factor   Bragg-Slater      t'
end if

! Print informations for the atoms

do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2
  if ((iT_sets(ii) == 1) .or. (iPrint >= 2)) then
    write(BondLbl,'(A)') AtomLbl(iAtom)
    Last_NonBlank = 0
    j = 0
    do i=1,17
      if (BondLbl(i:i) /= ' ') then
        j = j+1
        BondLbl(j:j) = BondLbl(i:i)
        Last_NonBlank = j
      end if
    end do
    do i=Last_NonBlank+1,17
      BondLbl(i:i) = ' '
    end do

    if (Num_Warnings > 0) then
      call Warnings_lp(iWarnings(ii),Warning,iLength)
      if (iT_sets(ii) == 1) then
        write(u6,'(1X,I3,5X,8X,A17,24X,F7.4,3X,A)') iAtom,BondLbl,T_Values(ij),Warning
      else if (Standard) then
        write(u6,'(1X,I3,5X,8X,A17,24X,A,3X,A)') iAtom,BondLbl,'Standard',Warning
      else
        write(u6,'(1X,I3,5X,8X,A17,24X,A,3X,A)') iAtom,BondLbl,'Skipped',Warning
      end if
    else
      if (iT_sets(ii) == 1) then
        write(u6,'(1X,I3,5X,8X,A17,24X,F7.4)') iAtom,BondLbl,T_Values(ij)
      else if (Standard) then
        write(u6,'(1X,I3,5X,8X,A17,24X,A)') iAtom,BondLbl,'Standard'
      else
        write(u6,'(1X,I3,5X,8X,A17,24X,A)') iAtom,BondLbl,'Skipped'
      end if
    end if
  end if
end do

write(u6,'(A)') repeat('-',79)

! Print informations for the bonds

do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=1,iAtom-1
    ij = iAtom*(iAtom-1)/2+jAtom
    if ((iT_sets(ij) == 1) .or. (iPrint >= 2)) then
      jj = jAtom*(jAtom+1)/2
      Factor = Zero
      write(BondLbl,'(3A)') AtomLbl(iAtom),'-',AtomLbl(jAtom)
      Radius_i = Bragg_Slater(iANr(iAtom))
      Radius_j = Bragg_Slater(iANr(jAtom))
      Bond_Length = sqrt((EC(1,ii)-EC(1,jj))**2+(EC(2,ii)-EC(2,jj))**2+(EC(3,ii)-EC(3,jj))**2)
      Bond_Max = Bond_Threshold*(Radius_i+Radius_j)
      Factor = Bond_Length/Bond_Max
      bs_t = Radius_i/(Radius_i+Radius_j)-Half

      Last_NonBlank = 0
      j = 0
      do i=1,17
        if (BondLbl(i:i) /= ' ') then
          j = j+1
          BondLbl(j:j) = BondLbl(i:i)
          Last_NonBlank = j
        end if
      end do
      do i=Last_NonBlank+1,17
        BondLbl(i:i) = ' '
      end do

      if (Num_Warnings > 0) then
        call Warnings_lp(iWarnings(ij),Warning,iLength)
        if (iT_sets(ij) == 1) then
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,F7.4,3X,A)') iAtom,jAtom,BondLbl,Factor,BS_t,T_Values(ij),Warning
        else if (Standard) then
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A,3X,A)') iAtom,jAtom,BondLbl,Factor,BS_t,'Standard',Warning
        else
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A,3X,A)') iAtom,jAtom,BondLbl,Factor,BS_t,'Skipped',Warning
        end if
      else
        if (iT_sets(ij) == 1) then
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,F7.4)') iAtom,jAtom,BondLbl,Factor,BS_t,T_Values(ij)
        else if (Standard) then
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A)') iAtom,jAtom,BondLbl,Factor,BS_t,'Standard'
        else
          write(u6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A)') iAtom,jAtom,BondLbl,Factor,BS_t,'Skipped'
        end if
      end if
    end if
  end do
end do

call mma_deallocate(AtomLbl)

return

end subroutine Print_T_Values
