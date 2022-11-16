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

function ExtNuc(ipExt,natom)
! Compute Z - ExtPot interactions

implicit real*8(A-H,O-Z)
#include "espf.fh"
#include "stdalloc.fh"
real*8 E, ExtNuc
real*8, allocatable :: Charge(:)
logical Found
integer Len

iPL = iPL_espf()

call qpg_dArray('Effective nuclear Charge',Found,Len)

if (Found) then
  call mma_allocate(Charge,Len,Label='Charge')
  if (Len /= nAtom) then
    write(6,*) 'ExtNuc: Len /= nAtom'
    call Abend()
  end if
else
  write(6,*) 'ExtNuc: Effective nuclear Charges not found.'
  call Abend()
end if
call Get_dArray('Effective nuclear Charge',Charge,nAtom)

E = Zero
do INuc=1,natom
  E = E+Charge(iNuc)*Work(ipExt+(INuc-1)*10)
end do
if ((E /= zero) .and. (iPL >= 3)) then
  write(6,*) ' '
  write(6,1000) E
end if
ExtNuc = E

call mma_deallocate(Charge)

return

1000 format(' Ext Pot/(QM nuclei and MM charges) energy =',F16.10,' hartrees')

end function ExtNuc
