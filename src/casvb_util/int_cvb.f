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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine int_cvb(iarr,nmax,nread,ifc)
      implicit real*8 (a-h,o-z)
#include "inpmod_cvb.fh"
      dimension iarr(*)

      if(inputmode.eq.2)then
        call gethi_cvb(iarr,nread)
        return
      endif
      nread=0
      if(nmax.le.0)goto 2000

c  Treat first field differently
      ifcuse=mod(ifc,4)
      if(ifcuse.ge.2)ifcuse=2
      call popfield_cvb(ifcuse)
      call rdint_cvb(iarr(1),ierr)
      if(ierr.gt.0)goto 1000
      nread=nread+1

      ifcuse=mod(ifc,2)
      do 100 i=2,nmax
      call popfield_cvb(ifcuse)
      call rdint_cvb(iarr(i),ierr)
      if(ierr.gt.0)goto 1000
      nread=nread+1
100   continue
      goto 2000
1000  continue
c  Crash if invalid field and IFC +4 :
      if(ierr.eq.4.and.ifc.ge.4)then
        write(6,*)' Invalid field found while reading integer!'
        call abend_cvb()
      endif
      call pushfield_cvb()
2000  continue
      if(inputmode.eq.1)then
        call sethi_cvb(iarr,nread)
      endif
      return
      end
