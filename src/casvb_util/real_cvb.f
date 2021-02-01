************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine real_cvb(arr,nmax,nread,ifc)
      implicit real*8 (a-h,o-z)
#include "inpmod_cvb.fh"
      dimension arr(nmax)

      call real_cvb_internal(arr)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine real_cvb_internal(arr)
      use iso_c_binding
      real*8, target :: arr(*)
      integer, pointer :: iarr(:)
      if(inputmode.eq.2)then
        call c_f_pointer(c_loc(arr(1)),iarr,[nread])
        call gethr_cvb(iarr,nread)
        nullify(iarr)
        return
      endif
      nread=0
      if(nmax.le.0)goto 2000

c  Treat first field differently
      ifcuse=mod(ifc,4)
      if(ifcuse.ge.2)ifcuse=2
      call popfield_cvb(ifcuse)
      call rdreal_cvb(arr(1),ierr)
      if(ierr.gt.0)goto 1000
      nread=nread+1

      ifcuse=mod(ifc,2)
      do 100 i=2,nmax
      call popfield_cvb(ifcuse)
      call rdreal_cvb(arr(i),ierr)
      if(ierr.gt.0)goto 1000
      nread=nread+1
100   continue
      goto 2000
1000  continue
c  Crash if invalid field and IFC +4 :
      if(ierr.eq.4.and.ifc.ge.4)then
        write(6,*)' Invalid field found while reading real!'
        call abend_cvb()
      endif
      call pushfield_cvb()
2000  continue
      if(inputmode.eq.1)then
        call c_f_pointer(c_loc(arr(1)),iarr,[nread])
        call sethr_cvb(iarr,nread)
        nullify(iarr)
      endif
      return
      end subroutine real_cvb_internal
*
      end
