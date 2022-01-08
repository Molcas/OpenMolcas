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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine GetUmat_Localisation(U,C,S,X,Scr,nBas,nOrb)
! Author: T.B. Pedersen
!
! Purpose: compute transformation matrix U=C^TSX.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, nOrb
real(kind=wp), intent(out) :: U(nOrb,nOrb), Scr(nBas,nOrb)
real(kind=wp), intent(in) :: C(nBas,nOrb), S(nBas,nBas), X(nBas,nOrb)
character(len=*), parameter :: SecNam = 'GetUmat_Localisation'

if ((nOrb < 1) .or. (nBas < 1)) return

call DGEMM_('N','N',nBas,nOrb,nBas,One,S,nBas,X,nBas,Zero,Scr,nBas)
call DGEMM_('T','N',nOrb,nOrb,nBas,One,C,nBas,Scr,nBas,Zero,U,nOrb)

end subroutine GetUmat_Localisation
