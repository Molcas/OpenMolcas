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
! Copyright (C) Francesco Aquilante                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine swap_full2rs(irc,iLoc,nRS,nDen,JSYM,XLT,Xab,add)

use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nnBstR
use Index_Functions, only: iTri
use Data_Structures, only: DSBA_Type
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iLoc, nRS, nDen, JSYM
type(DSBA_Type), intent(in) :: XLT(nDen)
real(kind=wp), intent(out) :: Xab(nRS,nDen)
logical(kind=iwp), intent(in) :: add
integer(kind=iwp) :: iab, iag, ias, ibg, ibs, iRab, iSyma, jDen, jRab, kRab
integer(kind=iwp), external :: cho_isao

#include "macros.fh"
unused_var(add)

!                                                                      *
!***********************************************************************
!                                                                      *
if (JSYM == 1) then ! TOTAL SYMMETRIC

  do jRab=1,nnBstR(jSym,iLoc)

    kRab = iiBstr(jSym,iLoc)+jRab
    iRab = IndRed(kRab,iLoc)

    iag = iRS2F(1,iRab) ! global address
    ibg = iRS2F(2,iRab)

    iSyma = cho_isao(iag) ! symmetry block; Sym(b)=Sym(a)

    ias = iag-ibas(iSyma) ! address within that symm block
    ibs = ibg-ibas(iSyma)
    iab = iTri(ias,ibs)

    do jDen=1,nDen

      Xab(jRab,jDen) = XLT(jDen)%SB(iSyma)%A1(iab)

    end do

  end do ! jRab loop

else

  write(u6,*) 'Wrong input parameters. JSYM = ',JSYM
  irc = 66
  call abend()

end if

irc = 0

return

end subroutine swap_full2rs
