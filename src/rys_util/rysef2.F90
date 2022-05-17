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

subroutine RysEF2(Iz2D,nArg,mArg,nRys,neMin,neMax,nfMin,nfMax,EFInt,meMin,meMax,mfMin,mfMax,PreFct,ixe,ixf,ixye,ixyf,nzeMin, &
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

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 Iz2D(nRys,mArg,3,0:neMax,0:nfMax), PreFct(mArg), EFInt(nArg,meMin:meMax,mfMin:mfMax)
! Statement function to compute canonical index
iCan(ixyz,ix,iz) = ixyz*(ixyz+1)*(ixyz+2)/6+(ixyz-ix)*(ixyz-ix+1)/2+iz

izf = nzfMax
ize = nzeMax
Indf = iCan(ixyf+izf,ixf,izf)
Inde = iCan(ixye+ize,ixe,ize)

select case (nRys)
  case (1)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = PreFct(iArg)*Iz2D(1,iArg,3,ize,izf)
    end do
  case (2)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Iz2D(1,iArg,3,ize,izf)+Iz2D(2,iArg,3,ize,izf))
    end do
  case (3)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Iz2D(1,iArg,3,ize,izf)+Iz2D(2,iArg,3,ize,izf)+Iz2D(3,iArg,3,ize,izf))
    end do
  case (4)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Iz2D(1,iArg,3,ize,izf)+Iz2D(2,iArg,3,ize,izf)+Iz2D(3,iArg,3,ize,izf)+ &
                              Iz2D(4,iArg,3,ize,izf))
    end do
  case (5)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = PreFct(iArg)*(Iz2D(1,iArg,3,ize,izf)+Iz2D(2,iArg,3,ize,izf)+Iz2D(3,iArg,3,ize,izf)+ &
                              Iz2D(4,iArg,3,ize,izf)+Iz2D(5,iArg,3,ize,izf))
    end do
  case default

    ! General code

    Indf = iCan(ixyf+izf,ixf,izf)
    Inde = iCan(ixye+ize,ixe,ize)
    do iArg=1,mArg
      EFInt(iArg,Inde,Indf) = Iz2D(1,iArg,3,ize,izf)
      do iRys=2,nRys
        EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf)+Iz2D(iRys,iArg,3,ize,izf)
      end do
      EFInt(iArg,Inde,Indf) = EFInt(iArg,Inde,Indf)*PreFct(iArg)
    end do
end select

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(neMin)
  call Unused_integer(nfMin)
  call Unused_integer(nzeMin)
  call Unused_integer(nzfMin)
end if

end subroutine RysEF2
