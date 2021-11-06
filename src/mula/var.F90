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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine var_to_qvar(var,qvar,ref,qref,alpha,trfName,ndata,nvar)
!  Purpose:
!    Transform coordinates given in input using tranformation
!    specified in input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 var(ndata,nvar)
real*8 par(ndata,nvar)
real*8 qvar(ndata,nvar)
real*8 ref(nvar), qref(nvar)
real*8 alpha(nvar)
character*80 trfName(nvar)
character*32 trfCode
character*32 Inline

do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ix = index(trfcode,'AS IT IS')
  ia = index(trfcode,'-AVG')
  ie = index(trfcode,'EXP')
  ir = index(trfcode,'RAD')
  id = index(trfcode,'DEG')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  angsc = 1.0d0
  if (ir > 0) angsc = 1.0d0
  if (id > 0) angsc = rpi/180.0d0
  do idata=1,ndata
    v = var(idata,ivar)
    if (ic > 0) then
      par(idata,ivar) = cos(angsc*v)
    else if (is > 0) then
      par(idata,ivar) = sin(angsc*v)
    else if ((ix > 0) .or. (ie > 0)) then
      par(idata,ivar) = v
    else
      write(6,*) ' TRFCODE ERROR.'
      call abend()
    end if
  end do
  ! Calculate refrence value.
  sum = 0.0d0
  do idata=1,ndata
    sum = sum+var(idata,ivar)
  end do
  ref(ivar) = sum/ndata
  if (ic > 0) then
    ref(ivar) = cos(angsc*ref(ivar))
  else if (is > 0) then
    ref(ivar) = sin(angsc*ref(ivar))
  end if

  ! Subtract reference value if requested.
  if (ia > 0) then
    do idata=1,ndata
      v = par(idata,ivar)
      if ((ic > 0) .or. (is > 0)) then
        par(idata,ivar) = ref(ivar)-v
      else
        par(idata,ivar) = v-ref(ivar)
      end if
    end do
  end if
end do

! Transform coordinates.
do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ie = index(trfcode,'EXP')
  ifit = index(trfcode,'FIT')
  if ((ie > 0) .and. (ifit == 0)) then
    Inline = trfName(ivar)(1:32)
    istart = index(Inline,'ALPHA=')
    istart = istart+6
    Inline = trfCode(istart:len(trfCode))
    istop = index(Inline,' ')
    istop = istop-1
    read(Inline(1:istop),*) alpha(ivar)
  end if
  do idata=1,ndata
    if (ie > 0) then
      qvar(idata,ivar) = 1.0d0-exp(-alpha(ivar)*par(idata,ivar))
    else
      qvar(idata,ivar) = par(idata,ivar)
    end if
  end do
end do

! Calculate refrence value of transformed coordinates.
do ivar=1,nvar
  sum = 0.0d0
  do idata=1,ndata
    sum = sum+qvar(idata,ivar)
  end do
  qref(ivar) = sum/ndata
end do

! Subtract reference value from transformed coordinates.
do ivar=1,nvar
  do idata=1,ndata
    qvar(idata,ivar) = qvar(idata,ivar)-qref(ivar)
  end do
end do

end subroutine var_to_qvar
!####
subroutine x_to_qvar(x,ref,qref,alpha,trfName,nDimX)
!  Purpose:
!    Transform the coordinates of a given point.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

implicit real*8(a-h,o-z)
real*8 x(nDimX)
real*8 par(nDimX)
real*8 ref(nDimX), qref(nDimX)
real*8 alpha(nDimx)
character*80 trfName(nDimX)
character*32 trfCode
character*32 Inline

do ivar=1,nDimX
  trfcode = trfName(ivar)(1:32)
  ia = index(trfcode,'-AVG')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  if (ic > 0) then
    par(ivar) = cos(x(ivar))
  else if (is > 0) then
    par(ivar) = sin(x(ivar))
  else
    par(ivar) = x(ivar)
  end if
  if (ia > 0) then
    v = par(ivar)
    if ((ic > 0) .or. (is > 0)) then
      par(ivar) = ref(ivar)-v
    else
      par(ivar) = v-ref(ivar)
    end if
  end if
end do

! Transform coordinates.
do ivar=1,nDimX
  trfcode = trfName(ivar)(1:32)
  ie = index(trfcode,'EXP')
  if (ie > 0) then
    Inline = trfName(ivar)(1:32)
    istart = index(Inline,'ALPHA=')
    istart = istart+6
    Inline = trfCode(istart:len(trfCode))
    istop = index(Inline,' ')
    istop = istop-1
    read(Inline(1:istop),*) alpha(ivar)
  end if
  if (ie > 0) then
    x(ivar) = 1.0d0-exp(-alpha(ivar)*par(ivar))
  else
    x(ivar) = par(ivar)
  end if
end do

! Subtract reference value from transformed coordinates.
do ivar=1,nDimX
  x(ivar) = x(ivar)-qref(ivar)
end do

end subroutine x_to_qvar
