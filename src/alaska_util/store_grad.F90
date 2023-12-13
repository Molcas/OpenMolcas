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
!  Store_Grad
!
!> @brief Store a vector in the gradients file
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Stores a gradient or coupling vector in the gradients file (GRADS).
!> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
!> \p iRoot = 0.
!>
!> @param[in] Grad  Gradient vector to store
!> @param[in] nGrad Length of \p Grad
!> @param[in] iRoot Root number to which the gradient belongs
!> @param[in] iNAC  First root of a coupling vector
!> @param[in] jNAC  Second root of a coupling vector
!***********************************************************************

subroutine Store_Grad(Grad,nGrad,iRoot,iNAC,jNAC)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nGrad, iRoot, iNAC, jNAC
real(kind=wp), intent(_IN_) :: Grad(nGrad)
integer(kind=iwp) :: iAd, idx, iSt, jSt, Length(1), LuGrad, nCoup, nRoots, TOC(5)
logical(kind=iwp) :: Found, BadFile
integer(kind=iwp), allocatable :: i_grad(:), i_nac(:)
integer(kind=iwp), external :: AixRm
character(len=5), parameter :: Filename = 'GRADS'

! Create GRADS file if it does not exist

call Get_iScalar('Number of roots',nRoots)
LuGrad = 20
call f_Inquire(Filename,Found)
if (.not. Found) call Create_Grads(Filename,nRoots,nGrad)

! Read the header

BadFile = .false.
call DaName(LuGrad,Filename)
iAd = 0
call iDaFile(LuGrad,2,TOC,size(TOC),iAd)
iAd = TOC(1)
call iDaFile(LuGrad,2,Length,1,iAd)
if (Length(1) /= nRoots) BadFile = .true.
iAd = TOC(2)
call iDaFile(LuGrad,2,Length,1,iAd)
if (Length(1) /= nGrad) BadFile = .true.
if (BadFile) then
  call DaClos(LuGrad)
  if (AixRm(Filename) /= 0) call Abend()
  call WarningMessage(1,'Number of roots and/or length of gradients do not match, re-creating GRADS file')
  call Create_Grads(Filename,nRoots,nGrad)
  call DaName(LuGrad,Filename)
  iAd = 0
  call iDaFile(LuGrad,1,TOC,size(TOC),iAd)
end if
nCoup = max(1,nRoots*(nRoots-1)/2)
call mma_Allocate(i_grad,nRoots)
call mma_Allocate(i_nac,nCoup)
iAd = TOC(3)
call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
iAd = TOC(4)
call iDaFile(LuGrad,2,i_nac,nCoup,iAd)

! Write the gradient or NAC vector

if (iRoot == 0) then
  if ((iNAC /= 0) .and. (jNAC /= 0)) then
    iSt = max(iNAC,jNAC)-1
    jSt = min(iNAC,jNAC)
    idx = iSt*(iSt-1)/2+jSt
    if (i_nac(idx) == 0) then
      i_nac(idx) = TOC(5)
      call dDaFile(LuGrad,1,Grad,nGrad,TOC(5))
      iAd = 0
      call iDaFile(LuGrad,1,TOC,size(TOC),iAd)
      iAd = TOC(4)
      call iDaFile(LuGrad,1,i_nac,nCoup,iAd)
    else
      iAd = i_nac(idx)
      call dDaFile(LuGrad,1,Grad,nGrad,iAd)
    end if
  end if
else
  idx = iRoot
  if (i_grad(idx) == 0) then
    i_grad(idx) = TOC(5)
    call dDaFile(LuGrad,1,Grad,nGrad,TOC(5))
    iAd = 0
    call iDaFile(LuGrad,1,TOC,size(TOC),iAd)
    iAd = TOC(3)
    call iDaFile(LuGrad,1,i_grad,nRoots,iAd)
  else
    iAd = i_grad(idx)
    call dDaFile(LuGrad,1,Grad,nGrad,iAd)
  end if
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'iRoot, iNAC, jNAC:',iRoot,iNAC,jNAC
write(u6,*) 'grads:',i_grad
write(u6,*) 'nacs:',i_nac
#endif
call DaClos(LuGrad)
call mma_Deallocate(i_grad)
call mma_Deallocate(i_nac)

end subroutine Store_Grad
