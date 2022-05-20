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

subroutine RysEF0(Ixy4D,Iz2D,nArg,mArg,nRys,neMax,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin, &
                  nzeMax,nzfMin,nzfMax)
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
integer(kind=iwp) :: nArg, mArg, nRys, neMax, nfMax, meMin, meMax, mfMin, mfMax, ixe, ixf, ixye, ixyf, nzeMin, nzeMax, nzfMin, &
                     nzfMax
real(kind=wp) :: Ixy4D(nRys,mArg), Iz2D(nRys,mArg,3,0:neMax,0:nfMax), EFInt(nArg,meMin:meMax,mfMin:mfMax), PreFct(mArg)
integer(kind=iwp) :: iArg, Inde, Indf, iRys, ize, izf

select case (nRys)
  case (1)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = PreFct(iArg)*Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)
        end do
      end do
    end do
  case (2)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)+Ixy4D(2,iArg)*Iz2D(2,iArg,3,ize,izf))
        end do
      end do
    end do
  case (3)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)+Ixy4D(2,iArg)*Iz2D(2,iArg,3,ize,izf)+ &
                                  Ixy4D(3,iArg)*Iz2D(3,iArg,3,ize,izf))
        end do
      end do
    end do
  case (4)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)+Ixy4D(2,iArg)*Iz2D(2,iArg,3,ize,izf)+ &
                                  Ixy4D(3,iArg)*Iz2D(3,iArg,3,ize,izf)+Ixy4D(4,iArg)*Iz2D(4,iArg,3,ize,izf))
        end do
      end do
    end do
  case (5)
    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)+Ixy4D(2,iArg)*Iz2D(2,iArg,3,ize,izf)+ &
                                  Ixy4D(3,iArg)*Iz2D(3,iArg,3,ize,izf)+Ixy4D(4,iArg)*Iz2D(4,iArg,3,ize,izf)+ &
                                  Ixy4D(5,iArg)*Iz2D(5,iArg,3,ize,izf))
        end do
      end do
    end do
  case default

    ! General code

    do izf=nzfMin,nzfMax
      Indf = C3_Ind(ixyf+izf,ixf,izf)-1
      do ize=nzeMin,nzeMax
        Inde = C3_Ind(ixye+ize,ixe,ize)-1
        do iArg=1,mArg
          EFInt(iArg,Inde,Indf) = Ixy4D(1,iArg)*Iz2D(1,iArg,3,ize,izf)
          do iRys=2,nRys
            EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf)+Ixy4D(iRys,iArg)*Iz2D(iRys,iArg,3,ize,izf)
          end do
          EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf)*PreFct(iArg)
        end do
      end do
    end do
end select

return

end subroutine RysEF0
