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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine SphAve()

use Genano_globals, only: MxLqn, nPrim, tDsym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, ij, iLqn, ind, iPtr, nDim
real(kind=wp) :: t

iPtr = 0
do iLqn=0,MxLqn
  nDim = nPrim(iLqn)*(nPrim(iLqn)+1)/2
  do ij=1,nDim
    ind = iPtr+ij
    t = Zero
    do i=0,2*iLqn
      t = t+tDsym(ind+i*nDim)
    end do
    t = t/(2*iLqn+1)
    do i=0,2*iLqn
      tDsym(ind+i*nDim) = t
    end do
  end do
  iPtr = iPtr+(2*iLqn+1)*nDim
end do
!write(u6,*) '*** Density matrix in SphAve ***'
!iBlk = 0
!do iLqn=0,MxLqn
!  do iShell=-iLqn,iLqn
!    iBlk = iBlk+1
!    if (nPrim(iLqn) > 0) then
!      write(u6,'(a,2i5)') ' Block',iLqn,iShell
!      call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
!    end if
!  end do
!end do

return

end subroutine SphAve
