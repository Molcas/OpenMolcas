************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       Subroutine fcopy(In,Ut,iErr)
       character*(*) In, Ut
       character*1024 myIn,myUt
       Integer c_open, c_copy, c_openw, c_close
       Integer rcIn, rcUt, rc
       iErr=0
       if(len(In).gt.1024.or. len(Ut).gt.1024) then
        Write(6,*) 'Error in fcopy: long filenames'
        iErr=1
        return
       endif
       Call PrgmTranslate(In,myIn,lIn)
       myIn(1+lIn:1+lIn)=Char(0)
       Call PrgmTranslate(Ut,myUt,lUt)
       myUt(1+lUt:1+lUt)=Char(0)
       rcIn=c_open(myIn)
       if(rcIn.lt.0) then
         Write(6,*) 'Can not open file ',myIn(1:lIn)
         iErr=1
         return
       endif
       rcUt=c_openw(myUt)
       if(rcUt.lt.0) then
         Write(6,*) 'Can not open file ',myUt(1:lUt)
         iErr=1
         return
       endif
       rc=c_copy(rcIn,rcUt)
       if(rc.lt.0) then
         Write(6,*) 'Can not copy file ',myIn(1:lIn)
         iErr=1
         return
       endif
       rc=c_close(rcIn)
       if(rc.lt.0) then
         Write(6,*) 'Can not close file ',myIn(1:lIn)
         iErr=1
         return
       endif
       rc=c_close(rcUt)
       if(rc.lt.0) then
         Write(6,*) 'Can not close file ',myUt(1:lUt)
         iErr=1
         return
       endif

       return
       end
