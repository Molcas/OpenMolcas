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

subroutine outmo(imo,ipower,cmo,clincomb,cout,nbas,nmo)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: imo, ipower, nbas, nmo
real(kind=wp), intent(in) :: cmo(nbas,nmo), clincomb(nmo)
real(kind=wp), intent(out) :: cout(nbas)
integer(kind=iwp) :: i
real(kind=wp), allocatable :: TmpMo(:)

if (imo /= 0) then
  cout(:) = cmo(:,imo)
  call power(cout,nbas,ipower)
else
  cout(:) = Zero
  call mma_allocate(TmpMo,nbas,label='TmpMo')
  do i=1,nmo
    if (clincomb(i) /= Zero) then
      TmpMo(:) = cmo(:,i)
      call power(TmpMo,nbas,ipower)
      cout(:) = cout(:)+clincomb(i)*TmpMo(:)
    end if
  end do
  call mma_deallocate(TmpMo)
end if

return

end subroutine outmo
