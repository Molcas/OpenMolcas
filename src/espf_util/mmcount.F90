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

subroutine MMCount(natom,nAtMM,IsMM)
! Count the number of MM atoms

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom
integer(kind=iwp), intent(out) :: nAtMM, IsMM(nAtom)
integer(kind=iwp) :: iAt, iAtom, iPL, nBla
logical(kind=iwp) :: Exists
integer(kind=iwp), allocatable :: IsMM1(:), NTC(:)
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

call Qpg_iArray('IsMM',Exists,nBla)
if (.not. Exists) then
  write(u6,'(A)') 'MMCount: IsMM not on the runfile'
  call Abend()
end if
if (nBla <= 0) then
  write(u6,'(A,I5)') 'MMCount: IsMM bad length:',nBla
  call Abend()
end if
call mma_allocate(IsMM1,nBla,Label='IsMM1')
call Get_iArray('IsMM',IsMM1,nBla)

call mma_allocate(NTC,natom,Label='NTC')
call Get_iArray('Atom -> Basis',NTC,natom)

do iAtom=1,natom
  IsMM(iAtom) = IsMM1(NTC(iAtom))
end do

call mma_deallocate(NTC)
call mma_deallocate(IsMM1)

nAtMM = 0
do iAt=1,natom
  if (IsMM(iAt) == 1) nAtMM = nAtMM+1
end do

if (nAtMM < 0) then
  write(u6,'(A)') 'Error in MMCount: nAtMM < 0!'
  call Quit_OnUserError()
else if (nAtMM > natom) then
  write(u6,'(A)') 'Error in MMCount: nAtMM >= natom!'
  call Quit_OnUserError()
else if ((nAtMM /= 0) .and. (iPL >= 3)) then
  write(u6,'(A,I5,A)') ' QM/MM: found ',nAtMM,' MM atoms'
end if

return

end subroutine MMCount
