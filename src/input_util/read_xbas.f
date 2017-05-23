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
* Copyright (C) 2017, Valera Veryazov                                  *
************************************************************************
        subroutine read_xbas(LuRd,iglobal,nxbas,xb_label,xb_bas,ierr)
        character *(*) xb_label(*)
        character *(*) xb_bas(*)
        character*128 Line
        icount=0
        iglobal=0
        ierr=0
100        read(LuRd,'(a)',err=300,end=300) Line
        if(Line.eq.' '.or.Line(1:3).eq.'END'.or.Line(1:3).eq.'end'
     &        .or. Line(1:3).eq.'End') goto 200
        i=2
        if(icount.eq.0) then
          i=index(Line,'.')
          if(i.eq.0) then
            iglobal=1
            xb_bas(1)=Line
            goto 200
           endif
         endif
        icount=icount+1
        nxbas=icount
        xb_label(nxbas)=Line(1:i-1)
        xb_bas(nxbas)=Line(i+1:)
        goto 100
200    continue
       return
300    ierr=1
       return
       end
