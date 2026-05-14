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

subroutine mol2vb2_cvb(vecvb,vecmol,isyml,fac,iwr,nsa,nsb)

use Symmetry_Info, only: Mul
use casvb_global, only: mxirrep, nda, ndet
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vecvb(ndet), vecmol(*)
integer(kind=iwp), intent(in) :: isyml, iwr, nsa, nsb
real(kind=wp), intent(in) :: fac
integer(kind=iwp) :: idet, indbet, indx, ioffsa, ioffsb, isa, isb, isyma, isymb, nnsa, nnsb, nstra(mxirrep), nstrb(mxirrep)
integer(kind=iwp), allocatable :: indxa(:), indxb(:)

call mma_allocate(indxa,nsa,label='indxa')
call mma_allocate(indxb,nsb,label='indxb')

call indxab_cvb(indxa,indxb,nstra,nstrb,nsa,nsb)

! Now loop casvb -> molcas
idet = 0
do isyma=1,mxirrep
  isymb = Mul(isyma,isyml)
  nnsa = nstra(isyma)
  nnsb = nstrb(isymb)
  if ((nnsa <= 0) .or. (nnsb <= 0)) cycle

  ioffsa = sum(nstra(1:isyma-1))
  ioffsb = sum(nstrb(1:isymb-1))

  do isb=1,nnsb
    indbet = indxb(isb+ioffsb)
    do isa=1,nnsa
      indx = indxa(isa+ioffsa)+(indbet-1)*nda
      idet = idet+1
      if (iwr == 0) then
        vecmol(idet) = vecvb(indx)
      else if (iwr == 1) then
        vecvb(indx) = vecmol(idet)
      else if (iwr == 2) then
        vecvb(indx) = vecvb(indx)+fac*vecmol(idet)
      end if
    end do
  end do
end do

call mma_deallocate(indxa)
call mma_deallocate(indxb)

return

end subroutine mol2vb2_cvb
