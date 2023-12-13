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

subroutine Write_QMMM(Coord,nAtIn,Iter)

use espf_global, only: MMI, MMO, QM
use Isotopes, only: PTab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Angstrom
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Iter, nAtIn
real(kind=wp), intent(in) :: Coord(3*nAtIn,Iter)
#include "LenIn.fh"
integer(kind=iwp) :: i, iAtIn, iAtNum, iAtOut, iAtTot, iFirst, k, Lu_XYZ, nAtOut, nAtTot
logical(kind=iwp) :: Found, isMMI, isMMO, isQM
character(len=LenIn) :: Symbol
character(len=16) :: FileName
character(len=LenIn), allocatable :: LabMMO(:)
integer(kind=iwp), allocatable :: AT(:)
real(kind=wp), allocatable :: Charge(:), CoordMMO(:,:)

call Qpg_dArray('MMO Coords',Found,nAtOut)
if (Found) then
  nAtOut = nAtOut/3
  nAtTot = nAtIn+nAtOut
  call mma_allocate(Charge,nAtIn)
  call mma_allocate(CoordMMO,3,nAtOut)
  call mma_allocate(LabMMO,nAtOut)
  call mma_allocate(AT,nAtTot)
  call Get_dArray('Nuclear charge',Charge,nAtIn)
  call Get_dArray('MMO Coords',CoordMMO,3*nAtOut)
  call Get_cArray('MMO Labels',LabMMO,LenIn*nAtOut)
  call Get_iArray('Atom Types',AT,nAtTot)

  ! Write one file per iteration and one for the current geometry
  do k=1,2
    if (k == 1) then
      write(FileName,'(I16.4)') Iter
      FileName = 'QMMMITXYZ.'//trim(adjustl(FileName))
    else
      FileName = 'QMMMENDXYZ'
    end if
    Lu_XYZ = 1
    call molcas_open(Lu_XYZ,FileName)
    write(Lu_XYZ,'(I6)') nAtTot
    write(Lu_XYZ,*)
    iAtIn = 1
    iAtOut = 1
    do iAtTot=1,nAtTot
      isQM = AT(iAtTot) == QM
      isMMI = AT(iAtTot) == MMI
      isMMO = AT(iAtTot) == MMO
      if (isQM .or. isMMI) then
        iAtNum = nint(Charge(iAtIn))
        iFirst = index(PTab(iAtNum),' ')+1
        Symbol = PTab(iAtNum)(iFirst:2)
        write(Lu_XYZ,100) Symbol,(Coord(3*(iAtIn-1)+i,Iter)*Angstrom,i=1,3)
        iAtIn = iAtIn+1
      else if (isMMO) then
        Symbol = LabMMO(iAtOut)
        if (index(Symbol,'_') > 1) then
          Symbol = Symbol(1:index(Symbol,'_')-1)
        end if
        write(Lu_XYZ,100) Symbol,(CoordMMO(i,iAtOut)*Angstrom,i=1,3)
        iAtOut = iAtOut+1
      end if
    end do
    close(Lu_XYZ)
  end do
  call mma_deallocate(Charge)
  call mma_deallocate(CoordMMO)
  call mma_deallocate(LabMMO)
  call mma_deallocate(AT)
end if

return

100 format(A6,3(1X,F12.6))

end subroutine Write_QMMM
