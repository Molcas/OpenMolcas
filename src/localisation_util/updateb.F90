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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine UpdateB(Col,nOrb2Loc,Lbl,nComp,Gamma_rot,iMO_s,iMO_t,Debug)
! Author: T.B. Pedersen
!
! Purpose: update MO dipole matrices for Boys localisation.
!          (Almost exact copy of UpdateP by Y. Carissan.)

use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOrb2Loc, nComp, iMO_s, iMO_t
real(kind=wp), intent(out) :: Col(nOrb2Loc,2)
real(kind=wp), intent(inout) :: Lbl(nOrb2Loc,nOrb2Loc,nComp)
real(kind=wp), intent(in) :: Gamma_rot
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: iComp
real(kind=wp) :: cos2g, cosg, cosing, Dss, Dst, Dtt, sin2g, sing
character(len=18) :: Label
#ifdef _DEBUGPRINT_
real(kind=wp) :: Dts, Tst
#endif

cosg = cos(Gamma_rot)
sing = sin(Gamma_rot)
cos2g = cosg*cosg
sin2g = sing*sing
cosing = cosg*sing

do iComp=1,nComp

  Dss = Lbl(iMO_s,iMO_s,iComp)
  Dst = Lbl(iMO_s,iMO_t,iComp)
  Dtt = Lbl(iMO_t,iMO_t,iComp)
# ifdef _DEBUGPRINT_
  Dts = Lbl(iMO_t,iMO_s,iComp)
  Tst = Dst-Dts
  if (abs(Tst) > 1.0e-14_wp) then
    write(u6,*) 'Broken symmetry in UpdateB!!'
    write(u6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(u6,*) 'Component  : ',iComp
    write(u6,*) 'Dst  = ',Dst
    write(u6,*) 'Dts  = ',Dts
    write(u6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateB','Broken symmetry!','[1]')
  end if
# endif

  Col(:,1) = Lbl(:,iMO_s,iComp)
  Col(:,2) = Lbl(:,iMO_t,iComp)

  Lbl(:,iMO_s,iComp) = cosg*Col(:,1)+sing*Col(:,2)
  Lbl(:,iMO_t,iComp) = cosg*Col(:,2)-sing*Col(:,1)

  Lbl(iMO_s,iMO_s,iComp) = Dss*cos2g+Dtt*sin2g+Two*Dst*cosing
  Lbl(iMO_t,iMO_s,iComp) = (Dtt-Dss)*cosing+Dst*(cos2g-sin2g)
  Lbl(iMO_s,iMO_t,iComp) = Lbl(iMO_t,iMO_s,iComp)
  Lbl(iMO_t,iMO_t,iComp) = Dtt*cos2g+Dss*sin2g-Two*Dst*cosing

  Lbl(iMO_s,:,iComp) = Lbl(:,iMO_s,iComp)
  Lbl(iMO_t,:,iComp) = Lbl(:,iMO_t,iComp)

# ifdef _DEBUGPRINT_
  Dst = Lbl(iMO_s,iMO_t,iComp)
  Dts = Lbl(iMO_t,iMO_s,iComp)
  Tst = Dst-Dts
  if (abs(Tst) > 1.0e-14_wp) then
    write(u6,*) 'Broken symmetry in UpdateB!!'
    write(u6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(u6,*) 'Component  : ',iComp
    write(u6,*) 'Dst  = ',Dst
    write(u6,*) 'Dts  = ',Dts
    write(u6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateB','Broken symmetry!','[2]')
  end if
# endif

end do

if (Debug) then
  write(u6,*) 'In UpdateB'
  write(u6,*) '----------'
  do iComp=1,nComp
    write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_s
    call RecPrt(Label,' ',Lbl(:,iMO_s,iComp),nOrb2Loc,1)
    write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_t
    call RecPrt(Label,' ',Lbl(:,iMO_t,iComp),nOrb2Loc,1)
  end do
end if

end subroutine UpdateB
