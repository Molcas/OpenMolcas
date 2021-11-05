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
!!-----------------------------------------------------------------------!
!!
      Subroutine var_to_qvar(var,qvar,ref,qref,alpha,                   &
     &  trfName,ndata,nvar)
!!
!!  Purpose:
!!    Transform coordinates given in input using tranformation
!!    specified in input.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 var  (ndata, nvar)
      Real*8 par  (ndata, nvar)
      Real*8 qvar  (ndata, nvar)
      Real*8 ref(nvar),qref (nvar)
      Real*8 alpha (nvar)
      Character*80 trfName (nvar)
      Character*32  trfCode
      Character*32  Inline
!!
!!
      Do ivar = 1,nvar
      trfcode = trfName(ivar)(1:32)
      ix = index(trfcode,'AS IT IS')
      ia = index(trfcode,'-AVG')
      ie = index(trfcode,'EXP')
      ir = index(trfcode,'RAD')
      id = index(trfcode,'DEG')
      ic = index(trfcode,'COS')
      is = index(trfcode,'SIN')
      angsc = 1.0d0
      If ( ir.gt.0 ) angsc = 1.0d0
      If ( id.gt.0 ) angsc = rpi/180.0d0
      Do idata = 1,ndata
      v = var(idata,ivar)
      If ( ic.gt.0 ) Then
      par(idata,ivar) = cos(angsc*v)
      Else If ( is.gt.0 ) Then
      par(idata,ivar) = sin(angsc*v)
      Else If (( ix.gt.0 ).or.( ie.gt.0 )) Then
      par(idata,ivar) = v
      Else
      Write(6,*)' TRFCODE ERROR.'
      call abend()
      End if
      End Do
!!
!!---- Calculate refrence value.
      sum = 0.0d0
      Do idata = 1,ndata
      sum = sum+var(idata,ivar)
      End Do
      ref(ivar) = sum/ndata
      If ( ic.gt.0 ) Then
      ref(ivar) = cos(angsc*ref(ivar))
      Else If ( is.gt.0 ) Then
      ref(ivar) = sin(angsc*ref(ivar))
      End If
!!
!!---- Subtract reference value if requested.
      If ( ia.gt.0 ) Then
      Do idata = 1,ndata
      v = par(idata,ivar)
      If (( ic.gt.0 ).or.( is.gt.0 )) Then
      par(idata,ivar) = ref(ivar)-v
      Else
      par(idata,ivar) = v-ref(ivar)
      End If
      End Do
      End If
      End Do
!!
!!---- Transform coordinates.
      Do ivar = 1,nvar
      trfcode = trfName(ivar)(1:32)
      ie   = index(trfcode,'EXP')
      ifit = index(trfcode,'FIT')
      If (( ie.gt.0 ).and.( ifit.eq.0 )) Then
      Inline = trfName(ivar)(1:32)
      istart = index(Inline,'ALPHA=')
      istart = istart+6
      Inline = trfCode(istart:len(trfCode))
      istop = index(Inline,' ')
      istop = istop-1
      Read(Inline(1:istop),*) alpha(ivar)
      End If
      Do idata = 1,ndata
      If ( ie.gt.0 ) Then
      qvar(idata,ivar) = 1.0d0-exp(-alpha(ivar)*                        &
     &          par(idata,ivar))
      Else
      qvar(idata,ivar) = par(idata,ivar)
      End If
      End Do
      End Do
!!
!!---- Calculate refrence value of transformed coordinates.
      Do ivar = 1,nvar
      sum = 0.0d0
      Do idata = 1,ndata
      sum = sum+qvar(idata,ivar)
      End Do
      qref(ivar) = sum/ndata
      End Do
!!
!!---- Subtract reference value from transformed coordinates.
      Do ivar = 1,nvar
      Do idata = 1,ndata
      qvar(idata,ivar) = qvar(idata,ivar)-qref(ivar)
      End Do
      End Do
!!
      End
!!
!!-----------------------------------------------------------------------!
!!-----------------------------------------------------------------------!
!!
      Subroutine x_to_qvar(x,ref,qref,alpha,trfName,nDimX)
!!
!!  Purpose:
!!    Transform the coordinates of a given point.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Implicit Real*8 ( a-h,o-z )
      Real*8  x (nDimX)
      Real*8 par (nDimX)
      Real*8 ref(nDimX),qref (nDimX)
      Real*8 alpha (nDimx)
      Character*80 trfName (nDimX)
      Character*32  trfCode
      Character*32  Inline
!!
!!
      Do ivar = 1,nDimX
      trfcode = trfName(ivar)(1:32)
      ia = index(trfcode,'-AVG')
      ic = index(trfcode,'COS')
      is = index(trfcode,'SIN')
      If ( ic.gt.0 ) Then
      par(ivar) = cos(x(ivar))
      Else If ( is.gt.0 ) Then
      par(ivar) = sin(x(ivar))
      Else
      par(ivar) = x(ivar)
      End If
      If ( ia.gt.0 ) Then
      v = par(ivar)
      If (( ic.gt.0 ).or.( is.gt.0 )) Then
      par(ivar) = ref(ivar)-v
      Else
      par(ivar) = v-ref(ivar)
      End If
      End If
      End Do
!!
!!---- Transform coordinates.
      Do ivar = 1,nDimX
      trfcode = trfName(ivar)(1:32)
      ie = index(trfcode,'EXP')
      If ( ie.gt.0 ) Then
      Inline = trfName(ivar)(1:32)
      istart = index(Inline,'ALPHA=')
      istart = istart+6
      Inline = trfCode(istart:len(trfCode))
      istop = index(Inline,' ')
      istop = istop-1
      Read(Inline(1:istop),*) alpha(ivar)
      End If
      If ( ie.gt.0 ) Then
      x(ivar) = 1.0d0-exp(-alpha(ivar)*par(ivar))
      Else
      x(ivar) = par(ivar)
      End If
      End Do
!!
!!---- Subtract reference value from transformed coordinates.
      Do ivar = 1,nDimX
      x(ivar) = x(ivar)-qref(ivar)
      End Do
!!
      End
