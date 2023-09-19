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

subroutine psym2_cvb(civec1,civec2,isymalf,isymbet,iasyind,ibsyind,osym,ips)

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension civec1(nda,ndb), civec2(nda,ndb)
dimension isymalf(nda), isymbet(ndb)
dimension iasyind(0:mxirrep), ibsyind(0:mxirrep)
dimension osym(mxirrep)

if (ips == 1) then
  do irp=1,nirrep
    if (isympr(irp) == 1) cycle
    do jrpa=1,nirrep
      jrpb = md2h(irp,jrpa)
      do ida=iasyind(jrpa-1)+1,iasyind(jrpa)
        inda = isymalf(ida)
        do idb=ibsyind(jrpb-1)+1,ibsyind(jrpb)
          civec1(inda,isymbet(idb)) = zero
        end do
      end do
    end do
  end do
else if (ips == 2) then
  do irp=1,nirrep
    osym(irp) = zero
    do jrpa=1,nirrep
      jrpb = md2h(irp,jrpa)
      do ida=iasyind(jrpa-1)+1,iasyind(jrpa)
        inda = isymalf(ida)
        do idb=ibsyind(jrpb-1)+1,ibsyind(jrpb)
          osym(irp) = osym(irp)+civec1(inda,isymbet(idb))*civec2(inda,isymbet(idb))
        end do
      end do
    end do
  end do
end if

return

end subroutine psym2_cvb
