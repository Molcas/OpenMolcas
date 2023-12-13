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

subroutine ReaW3(Ww,Wx,aSGrp,beSGrp,bSGrp,LunAux)
! this routine does:
! define Ww(a",be",b",i) <- (a",be"|b",i)
!
! integrals (a",be"|b",i) are stored in files V3xxyyzz for xx>=yy
!
! Storing of integrals in V3files
! do i=1,no
!   record of V3(a"be",b",_i): a">=be",b" for given i
! end do

use chcc_global, only: DimSGrpa, DimSGrpbe, no
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Ww(*), Wx(*)
integer(kind=iwp), intent(in) :: aSGrp, beSGrp, bSGrp, LunAux
integer(kind=iwp) :: dima, dimb, dimbe
character(len=8) :: LunName

!0 def dimensions

dima = DimSGrpa(aSGrp)
dimb = DimSGrpa(bSGrp)
dimbe = DimSGrpbe(beSGrp)

if (aSGrp > beSGrp) then
  !1 case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
  !1.1 make Name
  call MkNameV3(aSGrp,beSGrp,bSGrp,'W3',LunName)
  !1.2 read (a",be"|b",i) from V3:(a",be"|b",i)
  call ReaW3hlp1(Ww,dima,dimbe,dimb,no,LunName,LunAux)

else if (aSGrp == beSGrp) then
  !2 case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
  !2.1 make Name
  call MkNameV3(aSGrp,beSGrp,bSGrp,'W3',LunName)
  !2.2 read (a",be"|b",i) from V3:(a">=be"|b",i)
  call ReaW3hlp2(Ww,Wx,dima,dimb,no,LunName,LunAux)

else
  !3 case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
  !3.1 make Name
  call MkNameV3(beSGrp,aSGrp,bSGrp,'W3',LunName)
  !3.2 read (a",be"|b",i) from V3:(be",a"|b",i)
  call ReaW3hlp3(Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)

end if

return

end subroutine ReaW3
