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

use FFPT_Global, only: nBas, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ipH0, ipOvlp, ipRR, ipTemp, nSize, nTemp

!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!call Hello()
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
call GetMem('H0','Allo','Real',ipH0,nSize)
call GetMem('Ovlp','Allo','Real',ipOvlp,nSize)
call GetMem('RR','Allo','Real',ipRR,nSize)
call GetMem('Temp','Allo','Real',ipTemp,nTemp)
!----------------------------------------------------------------------*
call RdInp_FFPT()
call PrInp_FFPT()
call PtAdd(Work(ipH0),Work(ipOvlp),Work(ipRR),nSize,Work(ipTemp),nTemp)
!----------------------------------------------------------------------*
call GetMem('Temp','Free','Real',ipTemp,nTemp)
call GetMem('RR','Free','Real',ipRR,nSize)
call GetMem('Ovlp','Free','Real',ipOvlp,nSize)
call GetMem('H0','Free','Real',ipH0,nSize)
!----------------------------------------------------------------------*
call FastIO('STATUS')

ireturn = 0

return

end subroutine FFPT
