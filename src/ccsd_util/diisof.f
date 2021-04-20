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
       subroutine diisof (diispoint,ndiis)
c
c     this routine open 1-ndiis Temp files for DIIS procedure
c     and store lun's in stack (diispoint)
c
c     diispoint - stack of lun numbers (I)
c     ndiis     - size of diis procedure (I)
c
       integer ndiis
       integer diispoint(1:4)
c
c     help variables
c
       integer lun,rc
c
c
       if (ndiis.gt.0) then
       call filemanager (1,lun,rc)
       diispoint(1)=lun
       end if
c
       if (ndiis.gt.1) then
       call filemanager (1,lun,rc)
       diispoint(2)=lun
       end if
c
       if (ndiis.gt.2) then
       call filemanager (1,lun,rc)
       diispoint(3)=lun
       end if
c
       if (ndiis.gt.3) then
       call filemanager (1,lun,rc)
       diispoint(4)=lun
       end if
c
       return
       end
