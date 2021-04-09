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
        subroutine DefT3par (noa,nsym)
c
c        this routine do:
c        0) Open t3nam file, with lunt3 lun
c        define parameters, required for T3 integral handling, namely
c        1) def T3IndPoss(i)
c          address possitions for all occupied orbitals in t3nam file
c        2) def T3Off(ii,isym)
c           relative shifts of address for ii-th block of R_i(a,bc)
c          for each symmetry
c
c     noa   - array with occupation numbers
c     nsym  - actual number of irreps
c
        implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
c
       integer noa(1:8)
       integer nsym
c
c        help variables
        integer iorb,ii,i,symi,length,posst
c
c
c0        open t3nam file
c        lunt3=1
        call daname (lunt3,t3nam)
c
c1        set address poiter to 0
        daddr(lunt3)=0
c
c2        first record in t3nam file is T3IntPoss
c       (emulate writing of T3IntPoss)
        call idafile (lunt3,0,[0],mbas,daddr(lunt3))
c
        iorb=0
c3        cycle over irreps
        do symi=1,nsym
c
c3.1      make mapd and mapi for  R_i(a,bc)
          call ccsort_t3grc0
     c         (3,8,4,4,4,0,symi,possri0,posst,mapdri,mapiri)
c
c3.2          cycle over occupied orbitals in symi
          do i=1,noa(symi)
c
c3.2.1      save initial addres for this orbital
            iorb=iorb+1
            T3IntPoss(iorb)=daddr(lunt3)
c
c3.2.2      emulate writing of mapd and mapp
            call idafile (lunt3,0,[0],513*6,daddr(lunt3))
            call idafile (lunt3,0,[0],8*8*8,daddr(lunt3))
c
c3.2.3      cycle over all blocks of R_i(a,bc), which will
c           be stored separately
            do ii=1,mapdri(0,5)
c
c3.2.3.1      def T3Off(ii,symi)
c              note, that iorb is always proper one, since only besides
c             first occ. orbital in given irrep T3Off is defined
              if (i.eq.1) then
                T3Off(ii,symi)=daddr(lunt3)-T3IntPoss(iorb)
              end if
c
c3.2.3.2      emulate writing of each block
              length=mapdri(ii,2)
              call ddafile (lunt3,0,[0.0d0],length,daddr(lunt3))
c
            end do
          end do
        end do
c
        return
        end
