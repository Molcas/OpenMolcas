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

subroutine setifnss_cvb()

use casvb_global, only: ifnss1, ifnss2, ndetvbs, nel
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iretval1, iretval2, n, nalfa, nbeta

ifnss1(:,:) = 0
ifnss2(:,:) = 0
ndetvbs(:,:) = 0

do n=0,nel
  do nalfa=(n+1)/2,n
    nbeta = n-nalfa
    if (nbeta > nalfa) cycle
    call icomb_cvb(n,nbeta,iretval1)
    call icomb_cvb(n,nbeta-1,iretval2)
    ifnss1(n,nalfa-nbeta) = iretval1-iretval2
    call icomb_cvb(n,nalfa,ifnss2(n,nalfa-nbeta))
    if (nalfa == nbeta) ifnss2(n,nalfa-nbeta) = (ifnss2(n,nalfa-nbeta)+1)/2
    call icomb_cvb(n,nalfa,ndetvbs(n,nalfa))
  end do
end do

return

end subroutine setifnss_cvb
