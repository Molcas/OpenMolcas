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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine asc2ab2_cvb(detvec,nvec,nel,nalf,nbet,ndet)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec, nel, nalf, nbet, ndet
real(kind=wp), intent(inout) :: detvec(ndet,nvec)
integer(kind=iwp) :: inddet, iorb, rc
integer(kind=iwp), allocatable :: locc(:), maxdet(:), mindet(:), nkdet(:), xdet(:,:)
real(kind=wp), external :: party_cvb

call mma_allocate(mindet,[0,nel],label='mindet')
call mma_allocate(maxdet,[0,nel],label='maxdet')
call mma_allocate(nkdet,[0,nel],label='nkdet')
call mma_allocate(xdet,[0,nel],[0,nalf],label='xdet')
call mma_allocate(locc,nel,label='locc')

do iorb=0,nel
  mindet(iorb) = max(iorb-nbet,0)
  maxdet(iorb) = min(iorb,nalf)
end do
call weight_cvb(xdet,mindet,maxdet,nalf,nel)
nkdet(:) = maxdet(:)
call occupy_cvb(nkdet,nel,locc,locc(nalf+1))
inddet = 1
do
  detvec(inddet,:) = party_cvb(locc,nel)*detvec(inddet,:)
  call loind_cvb(nel,nalf,nkdet,mindet,maxdet,locc,locc(nalf+1),inddet,xdet,rc)
  if (rc == 0) exit
end do

call mma_deallocate(mindet)
call mma_deallocate(maxdet)
call mma_deallocate(nkdet)
call mma_deallocate(xdet)
call mma_deallocate(locc)

return

end subroutine asc2ab2_cvb
