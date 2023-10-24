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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mkvbinfo_cvb()

use casvb_global, only: iapr, ibpr, iconfs, idetvb, ixapr, ixbpr, nalf, nbet, nconf, nconfion_fr, nda, ndb, ndetvb, nel, nfrag, &
                        noe, norb

implicit none

if (nfrag > 1) then
  call dpgendet_cvb()
else
  call vbgendet_cvb(iapr,ixapr,ibpr,ixbpr,iconfs,idetvb,nconf,nconfion_fr(0,1),nda,ndb,ndetvb,nel,noe,nalf,nbet,norb)
end if

return

end subroutine mkvbinfo_cvb
