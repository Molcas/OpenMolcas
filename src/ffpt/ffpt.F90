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
! Copyright (C) 1992, Markus P. Fuelscher                              *
!               1999, Roland Lindh                                     *
!               2005, Anders Ohrn                                      *
!***********************************************************************

subroutine FFPT(ireturn)
!***********************************************************************
!                                                                      *
!                    ######  ######  #####    #####                    *
!                    #       #       #    #     #                      *
!                    #####   #####   #    #     #                      *
!                    #       #       #####      #                      *
!                    #       #       #          #                      *
!                    #       #       #          #                      *
!                                                                      *
!     A utility to perform finite field perturbation calculations      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!     Modified for dynamic memory allocation by                        *
!     R. Lindh March 1999                                              *
!                                                                      *
!     Added local perturbation by                                      *
!     A. Ohrn October 2005                                             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nBas, nSym, Cleanup
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, nSize, nTemp
real(kind=wp), allocatable :: H0(:), RR(:), Temp(:)

!----------------------------------------------------------------------*
call MkCom()
call Rd1Int_FFPT()
!----------------------------------------------------------------------*
nSize = 0
nTemp = 0
do i=1,nSym
  nSize = nSize+nBas(i)*(nBas(i)+1)/2
  nTemp = nTemp+nBas(i)
end do
nTemp = nTemp**2+4
nSize = nSize+4
call mma_allocate(H0,nSize,label='H0')
call mma_allocate(RR,nSize,label='RR')
call mma_allocate(Temp,nTemp,label='Temp')
!----------------------------------------------------------------------*
call RdInp_FFPT()
call PrInp_FFPT()
call PtAdd(H0,RR,nSize,Temp,nTEmp)
!----------------------------------------------------------------------*
call mma_deallocate(H0)
call mma_deallocate(RR)
call mma_deallocate(Temp)
!----------------------------------------------------------------------*
call FastIO('STATUS')
!----------------------------------------------------------------------*
call Cleanup()

ireturn = 0

return

end subroutine FFPT
