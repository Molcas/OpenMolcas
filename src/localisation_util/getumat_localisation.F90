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

subroutine GetUmat_Localisation(U,C,S,X,Scr,lScr,nBas,nOrb)
! Author: T.B. Pedersen
!
! Purpose: compute transformation matrix U=C^TSX.

implicit none
real*8 U(*), C(*), S(*), X(*)
integer lScr
real*8 Scr(lScr)
integer nBas, nOrb
character*80 Txt
character*20 SecNam
parameter(SecNam='GetUmat_Localisation')
real*8 d0, d1
parameter(d0=0.0d0,d1=1.0d0)
integer Need

if ((nOrb < 1) .or. (nBas < 1)) return

Need = nBas*nOrb
if (lScr < Need) then
  write(Txt,'(A,I9,A,I9)') 'lScr =',lScr,'     Need =',Need
  call SysAbendMsg(SecNam,'Insufficient dimension of scratch array!',Txt)
end if

call DGEMM_('N','N',nBas,nOrb,nBas,d1,S,nBas,X,nBas,d0,Scr,nBas)
call DGEMM_('T','N',nOrb,nOrb,nBas,d1,C,nBas,Scr,nBas,d0,U,nOrb)

end subroutine GetUmat_Localisation
