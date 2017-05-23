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
      subroutine initired
      implicit real*8 (a-h,o-z)
cbs   initialize all information for ireducible representations
cbs   later on, it might be useful to have a switch for
cbs    changing to other orders of IREDs like e.g. in TURBOMOLE
c
c
c   HOW2ADD another symmetry:
c
c   1. add it in readbas.f to be accepted. Add the number of IRs
c
c   2. copy one of the symmetry-blocks in this subroutine and
c      edit the multiplication-table for the group
c
c   3. assign the right IRs to L_X, L_Y and L_Z
c
c   that is  all. Good luck!!!
c
#include "para.fh"
#include "ired.fh"
      character*3 symmetry
      symmetry='D2H'  ! MOLCAS-Version
      if (symmetry.eq.'D2H') then
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
      mult(5,1)=5
      mult(6,1)=6
      mult(7,1)=7
      mult(8,1)=8
c
      mult(3,2)=4
      mult(4,2)=3
      mult(5,2)=6
      mult(6,2)=5
      mult(7,2)=8
      mult(8,2)=7
c
      mult(4,3)=2
      mult(5,3)=7
      mult(6,3)=8
      mult(7,3)=5
      mult(8,3)=6
c
      mult(5,4)=8
      mult(6,4)=7
      mult(7,4)=6
      mult(8,4)=5
c
      mult(6,5)=2
      mult(7,5)=3
      mult(8,5)=4
c
      mult(7,6)=4
      mult(8,6)=3
c
      mult(8,7)=2
c
C
      do ired=1,8
      mult(ired,ired)=1
      enddo
      do irun=2,8
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
CBS   write(6,*)
CBS   write(6,*)
CBS  *'multiplicitation table (atkins,child and phillips)'
CBS   write(6,*)
CBS   do ired=1,8
CBS   write(6,'(8I5)') (mult(jred,ired),jred=1,8)
CBS   write(6,*)
CBS   enddo

c
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   assume same order of ireds as Atkins Child and Phillips use..
cbs   would lead to an order with 1 to 1, 2 to 2 ...
cbs   however, this is the molecule/ seward order.
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=6
      iredorder(4)=7
      iredorder(5)=8
      iredorder(6)=5
      iredorder(7)=3
      iredorder(8)=2
      do ired=1,8
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(5)
      ipow2ired(0,0,1)=iredorder(6)
      ipow2ired(0,1,0)=iredorder(7)
      ipow2ired(1,0,0)=iredorder(8)
c     write(6,*) 'interacting IRs '
      do ired=1,8
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2V') then
cbs   1. A1 2. A2 3. B1 4. B2
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c
      mult(3,2)=4
      mult(4,2)=3
c
      mult(4,3)=2
C
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(6,*)
      write(6,*)
     *'multiplicitation table '
      write(6,*)
      do ired=1,4
      write(6,'(4I5)') (mult(jred,ired),jred=1,4)
      write(6,*)
      enddo

c
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   this is the molecule/ seward order.
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=2
      iredorder(4)=3
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(1)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(3)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'D2 ') then
cbs   1. A1 2. B1 3. B2 4. B3
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c
      mult(3,2)=4
      mult(4,2)=3
      mult(4,3)=2
C
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(6,*)
      write(6,*)
     *'multiplicitation table '
      write(6,*)
      do ired=1,4
      write(6,'(4I5)') (mult(jred,ired),jred=1,4)
      write(6,*)
      enddo

c
      IRLX=4
      IRLY=3
      IRLZ=2
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(1)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(3)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2H') then
cbs   assume 1.Ag 2.Au 3.Bg 4.Bu
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c
      mult(3,2)=4
      mult(4,2)=3
c
      mult(4,3)=2
C
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(6,*)
      write(6,*)
     *'multiplicitation table '
      write(6,*)
      do ired=1,4
      write(6,'(4I5)') (mult(jred,ired),jred=1,4)
      write(6,*)
      enddo

c
      IRLX=3
      IRLY=3
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(3)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'CS ') then
      write(6,*) 'CS in initired '
cbs   assume 1.A' 2.A'
      mult(2,1)=2
C
      do ired=1,2
      mult(ired,ired)=1
      enddo
      do irun=2,2
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(6,*)
      write(6,*)
     *'multiplicitation table '
      write(6,*)
      do ired=1,2
      write(6,'(2I5)') (mult(jred,ired),jred=1,2)
      write(6,*)
      enddo

c
      IRLX=2
      IRLY=2
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      do ired=1,2
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(2)
      ipow2ired(0,1,1)=iredorder(2)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(1)
      ipow2ired(1,0,0)=iredorder(1)
c     write(6,*) 'interacting IRs '
      do ired=1,2
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      endif
      return
      end
