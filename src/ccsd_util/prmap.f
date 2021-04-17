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
        subroutine prmap (mapd,mapi)
c

       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
       integer ii,i,j
       ii=mapd(0,5)
       do i=0,ii
       write (6,99) i,(mapd(i,j),j=1,6)
99      format (i3,6x,i10,1x,5(i6,1x))
       end do
c
        write (6,*) mapi(1,1,1),mapi(2,1,1)
       return
       end
