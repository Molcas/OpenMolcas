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

subroutine GenBMp_Localisation(D,C,X,nShell,iSym,ColD,ColC,ColX,PreFix)

implicit real*8(a-h,o-z)
#include "Molcas.fh"
real*8 D(nShell,nShell), C(nShell,*), X(nShell,*)
character*1 ColD, ColC, ColX
character*2 PreFix
#include "inflocal.fh"

character*12 FilNam

! Generate bitmap for density.
! ----------------------------

write(FilNam,'(A2,A5,I1,A4)') PreFix,'Dnsty',iSym,'.bmp'
call GenBMp_Loc(D,nShell,nShell,FilNam,ColD)

! Generate bitmap for original MOs.
! ---------------------------------

write(FilNam,'(A2,A5,I1,A4)') PreFix,'MOrig',iSym,'.bmp'
call GenBMp_Loc(C,nShell,nOrb2Loc(iSym),FilNam,ColC)

! Generate bitmap for localised MOs.
! ----------------------------------

write(FilNam,'(A2,A5,I1,A4)') PreFix,'MOloc',iSym,'.bmp'
call GenBMp_Loc(X,nShell,nOrb2Loc(iSym),FilNam,ColX)

end subroutine GenBMp_Localisation
