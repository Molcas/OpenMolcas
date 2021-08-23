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
       subroutine saverest2 (lunrst,energy,niter,iokey,daddr)
c
c     this routine save restart informations:
c     energy, niter
c     to prepaired possition in lunrst
c

#include "SysDef.fh"
       integer lunrst,niter,iokey,daddr,idum(1)
       real*8 energy,dum(1)
c
c1    write energy,niter
       if (iokey.eq.1) then
c      Fortran IO
       write (lunrst) energy,niter
c
       else
c      MOLCAS IO
       dum(1)=energy
       call ddafile (lunrst,1,dum,1,daddr)
       idum(1)=niter
       call idafile (lunrst,1,idum,1,daddr)
       end if
c
       return
       end
