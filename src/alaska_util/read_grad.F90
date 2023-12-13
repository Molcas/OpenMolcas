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
!  Read_Grad
!
!> @brief Read a vector from the gradients file
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Reads a gradient or coupling vector from the gradients file (GRADS).
!> For a gradient use \p iNAC, \p jNAC = 0. For a coupling vector use
!> \p iRoot = 0.
!>
!> @param[out] Grad  Gradient vector to read
!> @param[in]  nGrad Length of \p Grad
!> @param[in]  iRoot Root number to which the gradient belongs
!> @param[in]  iNAC  First root of a coupling vector
!> @param[in]  jNAC  Second root of a coupling vector
!>
!> @return \p 0 if the vector was not found, \p 1 if the vector was
!>         found, \p -1 if the vector was marked as non-computable
!***********************************************************************

function Read_Grad(Grad,nGrad,iRoot,iNAC,jNAC)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: Read_Grad
integer(kind=iwp), intent(in) :: nGrad, iRoot, iNAC, jNAC
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: iAd, iDum(1), idx, iSt, jSt, LuGrad, nCoup, nRoots, TOC(5)
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: i_grad(:), i_nac(:)
character(len=5), parameter :: Filename = 'GRADS'

! If the GRADS file does not exist, there is no gradient

call f_Inquire(Filename,Found)
if (.not. Found) then
  Read_Grad = 0
else

  ! Read the header

  LuGrad = 20
  call DaName(LuGrad,Filename)
  iAd = 0
  call iDaFile(LuGrad,2,TOC,size(TOC),iAd)
  call iDaFile(LuGrad,2,iDum,1,iAd)
  nRoots = iDum(1)
  if (max(iRoot,iNAC,jNAC) > nRoots) then
    call WarningMessage(2,'Bad number of roots in GRADS file')
    call Abend()
  end if
  call iDaFile(LuGrad,2,iDum,1,iAd)
  if (iDum(1) /= nGrad) then
    call WarningMessage(2,'Bad length in GRADS file')
    call Abend()
  end if
  nCoup = max(1,nRoots*(nRoots-1)/2)
  call mma_allocate(i_grad,nRoots)
  call mma_allocate(i_nac,nCoup)
  call iDaFile(LuGrad,2,i_grad,nRoots,iAd)
  call iDaFile(LuGrad,2,i_nac,nCoup,iAd)

  ! Read the gradient or NAC vector

  if (iRoot == 0) then
    if ((iNAC /= 0) .and. (jNAC /= 0)) then
      iSt = max(iNAC,jNAC)-1
      jSt = min(iNAC,jNAC)
      idx = iSt*(iSt-1)/2+jSt
      iAd = i_nac(idx)
    else
      iAd = -1
    end if
  else
    idx = iRoot
    iAd = i_grad(idx)
  end if

  if (iAd == 0) then
    Read_Grad = 0
  else if (iAd < 0) then
    Read_Grad = -1
  else
    call dDaFile(LuGrad,2,Grad,nGrad,iAd)
    Read_Grad = 1
  end if

  call DaClos(LuGrad)
  call mma_deallocate(i_grad)
  call mma_deallocate(i_nac)

end if
if (Read_Grad <= 0) call FZero(Grad,nGrad)

end function Read_Grad
