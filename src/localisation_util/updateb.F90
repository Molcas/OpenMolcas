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

subroutine UpdateB(Col,nOrb2Loc,ipLbl,nComp,Gamma_rot,iMO_s,iMO_t,Debug)
! Author: T.B. Pedersen
!
! Purpose: update MO dipole matrices for Boys localisation.
!          (Almost exact copy of UpdateP by Y. Carissan.)

implicit real*8(a-h,o-z)
real*8 Col(nOrb2Loc,2)
integer ipLbl(nComp)
logical Debug
#include "WrkSpc.fh"
character*18 Label

cosg = cos(Gamma_rot)
sing = sin(Gamma_rot)
cos2g = cosg*cosg
sin2g = sing*sing
cosing = cosg*sing

do iComp=1,nComp

  ip = ipLbl(iComp)
  ip0 = ip-1

  kOff_s = ip0+nOrb2Loc*(iMO_s-1)
  kOff_t = ip0+nOrb2Loc*(iMO_t-1)
  kOff_ss = kOff_s+iMO_s
  kOff_st = kOff_t+iMO_s
  kOff_tt = kOff_t+iMO_t
  Dss = Work(kOff_ss)
  Dst = Work(kOff_st)
  Dtt = Work(kOff_tt)
# if defined (_DEBUGPRINT_)
  kOff_ts = kOff_s+iMO_t
  Dts = Work(kOff_ts)
  Tst = Dst-Dts
  if (abs(Tst) > 1.0d-14) then
    write(6,*) 'Broken symmetry in UpdateB!!'
    write(6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(6,*) 'Component  : ',iComp
    write(6,*) 'Dst  = ',Dst
    write(6,*) 'Dts  = ',Dts
    write(6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateB','Broken symmetry!','[1]')
  end if
# endif

  call dCopy_(nOrb2Loc,Work(kOff_s+1),1,Col(1,1),1)
  call dCopy_(nOrb2Loc,Work(kOff_t+1),1,Col(1,2),1)

  call dScal_(nOrb2Loc,cosg,Work(kOff_s+1),1)
  call dAXPY_(nOrb2Loc,sing,Col(1,2),1,Work(kOff_s+1),1)
  call dScal_(nOrb2Loc,cosg,Work(kOff_t+1),1)
  call dAXPY_(nOrb2Loc,-sing,Col(1,1),1,Work(kOff_t+1),1)

  Work(kOff_s+iMO_s) = Dss*cos2g+Dtt*sin2g+2.0d0*Dst*cosing
  Work(kOff_s+iMO_t) = (Dtt-Dss)*cosing+Dst*(cos2g-sin2g)
  Work(kOff_t+iMO_s) = Work(kOff_s+iMO_t)
  Work(kOff_t+iMO_t) = Dtt*cos2g+Dss*sin2g-2.0d0*Dst*cosing

  call dCopy_(nOrb2Loc,Work(kOff_s+1),1,Work(ip0+iMO_s),nOrb2Loc)
  call dCopy_(nOrb2Loc,Work(kOff_t+1),1,Work(ip0+iMO_t),nOrb2Loc)

# if defined (_DEBUGPRINT_)
  Dst = Work(kOff_st)
  Dts = Work(kOff_ts)
  Tst = Dst-Dts
  if (abs(Tst) > 1.0d-14) then
    write(6,*) 'Broken symmetry in UpdateB!!'
    write(6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(6,*) 'Component  : ',iComp
    write(6,*) 'Dst  = ',Dst
    write(6,*) 'Dts  = ',Dts
    write(6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateB','Broken symmetry!','[2]')
  end if
# endif

end do

if (Debug) then
  write(6,*) 'In UpdateB'
  write(6,*) '----------'
  do iComp=1,nComp
    write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_s
    ip = ipLbl(iComp)+nOrb2Loc*(iMO_s-1)
    call RecPrt(Label,' ',Work(ip),nOrb2Loc,1)
    write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_t
    ip = ipLbl(iComp)+nOrb2Loc*(iMO_t-1)
    call RecPrt(Label,' ',Work(ip),nOrb2Loc,1)
  end do
end if

end subroutine UpdateB
