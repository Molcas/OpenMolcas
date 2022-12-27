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
! Copyright (C) 1993,1996, Markus P. Fuelscher                         *
!               1993, Per Ake Malmqvist                                *
!***********************************************************************

subroutine SORT1B()
!***********************************************************************
!                                                                      *
!     Purpose: Phase 1 of the bin sorting algorithm                    *
!              The integral generation is completed and the remaining  *
!              2el integrals stored in the bins are dumped to disk.    *
!                                                                      *
!     Called from: Seward_main                                         *
!                                                                      *
!     Calls to : PKI4,PKR8,SetVec,ISORTX,I4Len,R8Len                   *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.-AA. Malmqvist                             *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     - modified to use a virtual disk                                 *
!       M. P. Fuelscher, University of Lund, Sweden, 1996              *
!                                                                      *
!***********************************************************************

use sort_data, only: lIndx, lInts, lwIBin, lwVBin, n_Int, nBin
use stdalloc, only: mma_deallocate
use Definitions, only: iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iBin, iOpt, iPrint, iRout

!----------------------------------------------------------------------*
!     pick up print level                                              *
!----------------------------------------------------------------------*

iRout = 82
iPrint = nPrint(iRout)
if (iPrint >= 99) write(u6,*) ' >>> Enter SORT1B <<<'

!----------------------------------------------------------------------*
!     dump remaining integrals to disk                                 *
!----------------------------------------------------------------------*

iOpt = 0 ! Always tight!
do iBin=1,nBin
  do while (n_Int(iBin) > 0)
    call SaveBin(iBin,iOpt)
  end do
end do

!----------------------------------------------------------------------*
!     release the work space used to store bins                        *
!----------------------------------------------------------------------*

call mma_deallocate(lwVBin)
call mma_deallocate(lwIBin)

call mma_deallocate(lIndx)
call mma_deallocate(lInts)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT1B
