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
c  **************************************
c  ** Low-level input parsing routines **
c  **************************************
      subroutine popfield_cvb(ifc)
c  IFC is field control :
c  IFC=1 --> read to end of line only, no new line -- DISABLED in MOLCAS
c  IFC=2 --> begin read from next line
      implicit real*8 (a-h,o-z)
c      character*8 string
      common /pop_cvb/ ifield,nfield,nfold
      save initpop
      data initpop/0/
      if(initpop.eq.0) then
      ifield=0
      nfield=0
      nfold=0
      endif
      initpop=1

      if((ifield.eq.nfield).or.ifc.eq.2)then
        nfold=nfield
        call rdline_cvb(nfield)
        ifield=1
      else
c  IFIELD > NFIELD will signify no read
        ifield=min(ifield+1,nfield+1)
      endif
      return
      end

      subroutine pushfield_cvb()
      implicit real*8 (a-h,o-z)
c      character*8 string
      common /pop_cvb/ ifield,nfield,nfold
      if(ifield.eq.1.or.nfield.eq.-1)then
        call pushline_cvb()
        ifield=nfold
        nfield=nfold
      else
        ifield=ifield-1
      endif
      return
      end
      subroutine rdstring_cvb(string,ierr)
c  Check if field is applicable:
      implicit real*8 (a-h,o-z)
      character*8 string
      common /pop_cvb/ ifield,nfield,nfold
      ierr=0
      if(nfield.eq.-1)ierr=1
      if(ifield.gt.nfield)ierr=2
      if(ierr.ne.0)then
        string='        '
        return
      endif
c      call gtstring_cvb(string,ifield)
      call gtany_cvb(string,idi,rdr,1,ifield,ierr)
      return
      end
      subroutine rdint_cvb(intval,ierr)
c  Check if field is applicable:
      implicit real*8 (a-h,o-z)
      character*8 string
      common /pop_cvb/ ifield,nfield,nfold
      ierr=0
      if(nfield.eq.-1)ierr=1
      if(ifield.gt.nfield)ierr=2
      if(ierr.ne.0)return
c      call gtint_cvb(intval,ifield,jerr)
      call gtany_cvb(string,intval,rdr,2,ifield,jerr)
      if(jerr.eq.1)then
        if(ifield.eq.1)ierr=3
        if(ifield.ne.1)ierr=4
      endif
      return
      end
      subroutine rdreal_cvb(realval,ierr)
c  Check if field is applicable:
      implicit real*8 (a-h,o-z)
      character*8 string
      common /pop_cvb/ ifield,nfield,nfold
      ierr=0
      if(nfield.eq.-1)ierr=1
      if(ifield.gt.nfield)ierr=2
      if(ierr.ne.0)return
c      call gtreal_cvb(realval,ifield,jerr)
      call gtany_cvb(string,idi,realval,3,ifield,jerr)
      if(jerr.eq.1)then
        if(ifield.eq.1)ierr=3
        if(ifield.ne.1)ierr=4
      endif
      return
      end
