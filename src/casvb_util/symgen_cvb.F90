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

subroutine symgen_cvb(nalf1,nbet1,nda1,ndb1,isymalf,isymbet,iasyind,ibsyind,ialfsym,ibetsym,irpdet,irpalf,irpbet,mingrph,maxgrph, &
                      nk,locc,lunocc,xalf,xbet,icount)

implicit real*8(a-h,o-w,y-z),integer(x)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension isymalf(nda1), isymbet(ndb1)
dimension iasyind(0:mxirrep), ibsyind(0:mxirrep)
dimension ialfsym(nda1), ibetsym(ndb1)
dimension irpdet(mxirrep), irpalf(mxirrep), irpbet(mxirrep)
dimension mingrph(0:norb), maxgrph(0:norb)
dimension nk(0:norb), locc(norb+1), lunocc(norb+1)
dimension xalf(0:norb,0:nalf1), xbet(0:norb,0:nbet1)
dimension icount(mxirrep)

! Alpha loop:
call izero(irpalf,mxirrep)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf1,0)
  maxgrph(iorb) = min(iorb,nalf1)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf1,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
index = 1
200 continue
irp = 1
do ia=1,nalf1
  irp = md2h(irp,ityp(locc(ia)))
end do
irpalf(irp) = irpalf(irp)+1
ialfsym(index) = irp
call loind_cvb(norb,nalf1,nk,mingrph,maxgrph,locc,lunocc,index,xalf,*200)
iasyind(0) = 0
do irp=1,mxirrep
  iasyind(irp) = iasyind(irp-1)+irpalf(irp)
end do
call izero(icount,mxirrep)
do ida=1,nda1
  irrep = ialfsym(ida)
  icount(irrep) = icount(irrep)+1
  isymalf(icount(irrep)+iasyind(irrep-1)) = ida
end do

! Beta loop:
call izero(irpbet,mxirrep)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet1,0)
  maxgrph(iorb) = min(iorb,nbet1)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet1,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
index = 1
400 continue
irp = 1
do ib=1,nbet1
  irp = md2h(irp,ityp(locc(ib)))
end do
irpbet(irp) = irpbet(irp)+1
ibetsym(index) = irp
call loind_cvb(norb,nbet1,nk,mingrph,maxgrph,locc,lunocc,index,xbet,*400)
ibsyind(0) = 0
do irp=1,mxirrep
  ibsyind(irp) = ibsyind(irp-1)+irpbet(irp)
end do
call izero(icount,mxirrep)
do idb=1,ndb1
  irrep = ibetsym(idb)
  icount(irrep) = icount(irrep)+1
  isymbet(icount(irrep)+ibsyind(irrep-1)) = idb
end do

do irp=1,mxirrep
  irpdet(irp) = 0
  do jrp=1,mxirrep
    irpdet(irp) = irpdet(irp)+irpalf(jrp)*irpbet(md2h(irp,jrp))
  end do
end do

return

end subroutine symgen_cvb
