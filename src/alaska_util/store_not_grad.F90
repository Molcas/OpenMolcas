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
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Store_Not_Grad
!
!> @brief Mark a vector as non-computable in the gradients file
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Marks a gradient or coupling vector in the gradients file (GRADS)
!> as non-computable. This is useful, for example, to allow SLAPAF
!> to proceed with an approximate coupling vector if it cannot be
!> computed with the current method.
!> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
!> \p iRoot = 0.
!>
!> @param[in] iRoot Root number to which the gradient belongs
!> @param[in] iNAC  First root of a coupling vector
!> @param[in] jNAC  Second root of a coupling vector
!***********************************************************************

subroutine Store_Not_Grad(iRoot,iNAC,jNAC)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iRoot, iNAC, jNAC
integer(kind=iwp) :: iAd, idx, iSt, jSt, Length(1), LuGrad, nCoup, nGrad, nRoots, TOC(5)
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: i_grad(:), i_nac(:)
character(len=5), parameter :: Filename = 'GRADS'

! Create GRADS file if it does not exist

call Get_iScalar('Number of roots',nRoots)
call Get_iScalar('Unique atoms',nGrad)
nGrad = 3*nGrad
LuGrad = 20
call f_Inquire(Filename,Found)
if (.not. Found) call Create_Grads(Filename,nRoots,nGrad)

! Read the header

call DaName(LuGrad,Filename)
iAd = 0
call iDaFile(LuGrad,2,TOC,size(TOC),iAd)
call iDaFile(LuGrad,2,Length,1,iAd)
if (Length(1) /= nRoots) then
  call WarningMessage(2,'Bad number of roots in GRADS file')
  call Abend()
end if
call iDaFile(LuGrad,2,Length,1,iAd)
if (Length(1) /= nGrad) then
  call WarningMessage(2,'Bad length in GRADS file')
  call Abend()
end if
nCoup = max(1,nRoots*(nRoots-1)/2)
call mma_Allocate(i_grad,nRoots)
call mma_Allocate(i_nac,nCoup)
call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
call iDaFile(LuGrad,2,i_nac,nCoup,iAd)

! Write the negative index that marks it as non-computable

if (iRoot == 0) then
  if ((iNAC /= 0) .and. (jNAC /= 0)) then
    iSt = max(iNAC,jNAC)-1
    jSt = min(iNAC,jNAC)
    idx = iSt*(iSt-1)/2+jSt
    i_nac(idx) = -1
    iAd = TOC(4)
    call iDaFile(LuGrad,1,i_nac,nCoup,iAd)
  end if
else
  idx = iRoot
  i_grad(idx) = -1
  iAd = TOC(3)
  call iDaFile(LuGrad,1,i_grad,nRoots,iAd)
end if

call DaClos(LuGrad)
call mma_Deallocate(i_grad)
call mma_Deallocate(i_nac)

end subroutine Store_Not_Grad
