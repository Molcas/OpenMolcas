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
      subroutine fstring_cvb(strings,nstring,istring,ncmp,ifc)
      implicit real*8 (a-h,o-z)
#include "inpmod_cvb.fh"
      character*(*) strings(nstring)
      character*8 string
      logical debug
      data debug/.false./

      if(inputmode.eq.2)then
        call gethfs_cvb(istring)
        if(debug)then
          write(6,*)' fstring :',istring
          if(istring.ne.0)write(6,*)' fstring :',strings(istring)
        endif
        return
      endif
      call popfield_cvb(ifc)
      call rdstring_cvb(string,ierr)
      do 100 istring=1,nstring
      if(string(1:ncmp).eq.strings(istring)(1:ncmp))then
c  For keywords starting by END we check more letters. This
c  implementation is a bit ungainly, but simpler to code :
        if(string(1:3).eq.'END')then
          if(string(4:ncmp+3).ne.strings(istring)(4:ncmp+3))goto 100
        endif
        goto 200
      endif
100   continue
      istring=0
      call pushfield_cvb()
200   continue
      if(inputmode.eq.1)then
        call sethfs_cvb(istring)
      endif
      if(debug)then
        write(6,*)' fstring :',istring
        if(istring.ne.0)write(6,*)' fstring :',strings(istring)
      endif
      return
      end
