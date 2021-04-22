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
        subroutine SaveX (X,length,Lun,LunName,keyopen,keyclose)
c
c       this routine do
c       1) keyopen = 1 - open LunName file with Lun
c                    2 - rewind Lun file
c                    3 - open LunName file with Lun with ACCESS='append'
c                 else - nothing (i.e) file is opened
c       2) write X  of dimension length
c       3) keyclose= 1 - close Lun file
c                 else - nothing
c
c
        implicit none
        integer length,Lun,keyopen,keyclose
        real*8 X(1)
        character*6 LunName
c
c1
        if (keyopen.eq.1) then
*         open (unit=Lun,file=LunName,form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
        else if (keyopen.eq.2) then
          rewind(Lun)
        else if (keyopen.eq.3) then
cmp!          open (unit=Lun,file=LunName,form='unformatted',
cmp!     c          ACCESS='append')

          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
          call append_file_u(Lun)

        end if
c
c2
        call wri_chcc (Lun,length,X(1))
c
c3
        if (keyclose.eq.1) then
          close (Lun)
        end if
c
        return
        end
