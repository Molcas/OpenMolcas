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
* Copyright (C) 2001-2005, Valera Veryazov                             *
************************************************************************

        subroutine molcas_open_ext2(Lu,f_Name,f_access,
     &   f_form,f_iostat,is_recl,f_recl,f_status,is_error)
        Integer Lu,f_recl,f_iostat
        Character*(*) f_Name
        Character*4096 RealName
        Character*(*) f_access,f_form,f_status
        Logical is_recl,is_error
        is_error=.false.
       Call prgmtranslate(f_Name, RealName, lRealName)
       if(index(RealName,'UNK_VAR').ne.0) then
         write(6,*) '*** attempt to open ',RealName(1:lRealName)
         RealName=f_Name
         lRealName=index(RealName,' ')
       endif

        if(is_recl) then
                   open(Unit=Lu,File=Realname(1:lRealName),
     &              status=f_status,err=100,
     &             access=f_access,  form=f_form, iostat=f_iostat,
     &             recl=f_recl)
        else

                   open(Unit=Lu,File=RealName(1:lRealName),
     &                   status=f_status,err=100,
     &             access=f_access,  form=f_form, iostat=f_iostat)
        endif
c        print *,'DEBUG open ',RealName(1:lRealName)
c        print *,'Unit ', Lu

        return
100     is_error=.true.
        end

        subroutine molcas_binaryopen_vanilla(Lu,f_name)
        integer Lu
        Character*(*) f_name
        Character*4096 RealName
c        RealName=f_Name
        Call PrgmTranslate(f_Name, RealName,lRealName)
c        print *,'DEBUG binopen ',RealName(1:lRealName)
        open(Unit=Lu,File=RealName(1:lRealName), form='unformatted')
        return
        end
