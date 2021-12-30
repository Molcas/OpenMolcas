!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GetUmat_T1(U,C,S,X,Scr,lScr,nBas,nOrb1,nOrb2)
! Purpose: compute transformation matrix U=C^TSX.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: lScr, nBas, nOrb1, nOrb2
real(kind=wp), intent(_OUT_) :: U(*)
real(kind=wp), intent(in) :: C(*), S(*), X(*)
real(kind=wp), intent(out) :: Scr(lScr)
integer(kind=iwp) :: Need
character(len=80) :: Txt
character(len=*), parameter :: SecNam = 'GetUmat_T1'

if ((nOrb1*nOrb2 < 1) .or. (nBas < 1)) return

Need = nBas*nOrb2
if (lScr < Need) then
  write(Txt,'(A,I9,A,I9)') 'lScr =',lScr,'     Need =',Need
  call SysAbendMsg(SecNam,'Insufficient dimension of scratch array!',Txt)
end if

call DGEMM_('N','N',nBas,nOrb2,nBas,One,S,nBas,X,nBas,Zero,Scr,nBas)
call DGEMM_('T','N',nOrb1,nOrb2,nBas,One,C,nBas,Scr,nBas,Zero,U,nOrb1)

end subroutine GetUmat_T1
