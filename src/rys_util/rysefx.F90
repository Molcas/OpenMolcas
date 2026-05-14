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
! Copyright (C) 1990,1994, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine RysEFX(xyz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,iye,iyf,ixye,ixyf,nzeMin,nzeMax, &
                  nzfMin,nzfMax)
!***********************************************************************
!                                                                      *
!     Object: kernel routine to assemble the integrals from the Ixy    *
!             and Iz integrals.                                        *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             August '90                                               *
!                                                                      *
!             Modified for decreased memory access January '94.        *
!***********************************************************************

use Index_Functions, only: C3_Ind
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, mArg, nRys, neMax, nfMax, meMin, meMax, mfMin, mfMax, ixe, ixf, iye, iyf, ixye, ixyf, &
                                 nzeMin, nzeMax, nzfMin, nzfMax
real(kind=wp), intent(in) :: xyz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg)
real(kind=wp), intent(inout) :: EFInt(nArg,meMin:meMax,mfMin:mfMax)
integer(kind=iwp) :: Inde, Indf, iRys, ize, izf

! General code

do izf=nzfMin,nzfMax
  Indf = C3_Ind(ixyf+izf,ixf,izf)-1
  do ize=nzeMin,nzeMax
    Inde = C3_Ind(ixye+ize,ixe,ize)-1
    EFInt(1:mArg,Inde,Indf) = xyz2D(1,:,1,ixe,ixf)*xyz2D(1,:,2,iye,iyf)*xyz2D(1,:,3,ize,izf)
    do iRys=2,nRys
      EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)+xyz2D(iRys,:,1,ixe,ixf)*xyz2D(iRys,:,2,iye,iyf)*xyz2D(iRys,:,3,ize,izf)
    end do
    EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)*PreFct(:)
  end do
end do

return

end subroutine RysEFX
