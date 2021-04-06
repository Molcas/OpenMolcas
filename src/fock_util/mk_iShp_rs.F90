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
Subroutine Mk_iShp_rs(iShp_rs,nShell)
Implicit None
Integer nShell
Integer iShp_rs( nShell*(nShell+1)/2 )

Integer iaSh, ibSh, iShp
Integer, External:: Cho_F2SP

! *** Mapping shell pairs from the full to the reduced set

Do iaSh=1,nShell
   Do ibSh=1,iaSh
      iShp = iaSh*(iaSh-1)/2 + ibSh
      iShp_rs(iShp) = Cho_F2SP(iShp)
  End Do
End Do

End Subroutine Mk_iShp_rs





