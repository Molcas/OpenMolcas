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
        Subroutine fileorb(filein,fileout)
        Character*(*) filein
        Character*(*) fileout
        Character*256 tmp
        Logical Exist
        if(index(filein,'/').ne.0) then
           fileout=filein
           goto 100
        endif
        tmp=' '
        call getenvf('MOLCAS_SUBMIT_DIR',tmp)
        if(tmp.ne.' ') then
         fileout=tmp(1:mylen(tmp))//'/'//filein
c         print *,'vv',fileout
         call f_inquire(fileout,Exist)
         if(Exist) goto 100
        endif
         fileout=filein
         call f_inquire(fileout,Exist)
         if(.not.Exist)        then
          tmp='file '//fileout(1:mylen(fileout))//' not found'
          Call WarningMessage(2,tmp)
          Call Quit_OnUserError()
         endif
100      continue
c         write(6,*) 'INPORB file=',fileout
c         call fcopy(fileout(1:mylen(fileout)),'INPORB')
         return
         end
