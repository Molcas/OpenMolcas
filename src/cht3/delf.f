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
       Subroutine DELF(FNAM,INUM1,INUM2)
       Implicit None
       Integer I,inum1,inum2
       character FNAM*6,FN*8
       FN(1:6)=FNAM
       do I=inum1,inum2
       write(fn(7:8),'(I2.2)')I
c       write(6,*)'File ',FN,' to be deleted'
       call Molcas_Open(8,fn)
       close(8,status='DELETE')
       enddo
       return
       end
