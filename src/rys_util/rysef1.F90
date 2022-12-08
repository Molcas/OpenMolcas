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

subroutine RysEF1(Iz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin,nzeMax,nzfMin, &
                  nzfMax)
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
integer(kind=iwp), intent(in) :: nArg, mArg, nRys, neMax, nfMax, meMin, meMax, mfMin, mfMax, ixe, ixf, ixye, ixyf, nzeMin, nzeMax, &
                                 nzfMin, nzfMax
real(kind=wp), intent(in) :: Iz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg)
real(kind=wp), intent(inout) :: EFInt(nArg,meMin:meMax,mfMin:mfMax)
integer(kind=iwp) :: Inde, Indf, iRys, ize, izf

select case (nRys)
  case (1)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = PreFct(:)*Iz2D(1,:,3,ize,izf)
      end do
    end do
  case (2)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = PreFct(:)*(Iz2D(1,:,3,ize,izf)+Iz2D(2,:,3,ize,izf))
      end do
    end do
  case (3)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = PreFct(:)*(Iz2D(1,:,3,ize,izf)+Iz2D(2,:,3,ize,izf)+Iz2D(3,:,3,ize,izf))
      end do
    end do
  case (4)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = PreFct(:)*(Iz2D(1,:,3,ize,izf)+Iz2D(2,:,3,ize,izf)+Iz2D(3,:,3,ize,izf)+Iz2D(4,:,3,ize,izf))
      end do
    end do
  case (5)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = PreFct(:)*(Iz2D(1,:,3,ize,izf)+Iz2D(2,:,3,ize,izf)+Iz2D(3,:,3,ize,izf)+Iz2D(4,:,3,ize,izf)+ &
                                             Iz2D(5,:,3,ize,izf))
      end do
    end do
  case default

    ! General code

    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        EFInt(1:mArg,Inde,Indf) = Iz2D(1,:,3,ize,izf)
        do iRys=2,nRys
          EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)+Iz2D(iRys,:,3,ize,izf)
        end do
        EFInt(1:mArg,Inde,Indf) = EFInt(1:mArg,Inde,Indf)*PreFct(:)
      end do
    end do
end select

return

end subroutine RysEF1
