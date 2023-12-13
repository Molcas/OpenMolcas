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

subroutine Ext_W4(V2,M1,nc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)
! this routine is a control routine to:
! Extract M1(m,a"b") <- V2(m,a'b')
! for all combinations of aGrp,bGrp,aSGrp,bSGrp

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: V2(*)
real(kind=wp), intent(_OUT_) :: M1(*)
integer(kind=iwp) :: nc, dima, dimb, dimab, dimapp, dimbpp, dimabpp, addapp, addbpp, aGrp, bGrp, aSGrp, bSGrp

if (aGrp == bGrp) then
  ! case V2(m,a'b')
  if (aSGrp == bSGrp) then
    ! subcase M1(m,a"b") <- V2(m,a'b')
    call Ext_W4hlp1(V2,M1,nc,dimab,dimapp,dimabpp,addapp)
  else
    ! subcase M1(m,a",b") <- V2(m,a'b')
    call Ext_W4hlp2(V2,M1,nc,dimab,dimapp,dimbpp,addapp,addbpp)
  end if
else
  ! case M1(m,a",b") <- V2(m,a',b')
  call Ext_W4hlp3(V2,M1,nc,dima,dimb,dimapp,dimbpp,addapp,addbpp)
end if

return

end subroutine Ext_W4
