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

function ExtNuc(Ext,natom)
! Compute Z - ExtPot interactions

use espf_global, only: MxExtPotComp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(in) :: Ext(MxExtPotComp,natom)
integer(kind=iwp) :: INuc, iPL, Len
real(kind=wp) :: E, ExtNuc
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Charge(:)
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

call qpg_dArray('Effective nuclear Charge',Found,Len)

if (Found) then
  call mma_allocate(Charge,Len,Label='Charge')
  if (Len /= nAtom) then
    write(u6,*) 'ExtNuc: Len /= nAtom'
    call Abend()
  end if
else
  write(u6,*) 'ExtNuc: Effective nuclear Charges not found.'
  call Abend()
end if
call Get_dArray('Effective nuclear Charge',Charge,nAtom)

E = Zero
do INuc=1,natom
  E = E+Charge(iNuc)*Ext(1,INuc)
end do
if ((E /= zero) .and. (iPL >= 3)) then
  write(u6,*) ' '
  write(u6,1000) E
end if
ExtNuc = E

call mma_deallocate(Charge)

return

1000 format(' Ext Pot/(QM nuclei and MM charges) energy =',F16.10,' hartrees')

end function ExtNuc
