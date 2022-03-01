!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine DefT3par (noa,nsym)
!
!        this routine do:
!        0) Open t3nam file, with lunt3 lun
!        define parameters, required for T3 integral handling, namely
!        1) def T3IndPoss(i)
!          address possitions for all occupied orbitals in t3nam file
!        2) def T3Off(ii,isym)
!           relative shifts of address for ii-th block of R_i(a,bc)
!          for each symmetry
!
!     noa   - array with occupation numbers
!     nsym  - actual number of irreps
!
        implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
!
       integer noa(1:8)
       integer nsym
!
!        help variables
        integer iorb,ii,i,symi,length,posst,idum(1)
        real*8 dum(1)
!
!
!0        open t3nam file
!        lunt3=1
        call daname (lunt3,t3nam)
!
!1        set address poiter to 0
        daddr(lunt3)=0
!
!2        first record in t3nam file is T3IntPoss
!       (emulate writing of T3IntPoss)
        idum(1)=0
        dum(1)=0.0d0
        call idafile (lunt3,0,idum,mbas,daddr(lunt3))
!
        iorb=0
!3        cycle over irreps
        do symi=1,nsym
!
!3.1      make mapd and mapi for  R_i(a,bc)
          call ccsort_t3grc0                                            &
     &         (3,8,4,4,4,0,symi,possri0,posst,mapdri,mapiri)
!
!3.2          cycle over occupied orbitals in symi
          do i=1,noa(symi)
!
!3.2.1      save initial addres for this orbital
            iorb=iorb+1
            T3IntPoss(iorb)=daddr(lunt3)
!
!3.2.2      emulate writing of mapd and mapp
            call idafile (lunt3,0,idum,513*6,daddr(lunt3))
            call idafile (lunt3,0,idum,8*8*8,daddr(lunt3))
!
!3.2.3      cycle over all blocks of R_i(a,bc), which will
!           be stored separately
            do ii=1,mapdri(0,5)
!
!3.2.3.1      def T3Off(ii,symi)
!              note, that iorb is always proper one, since only besides
!             first occ. orbital in given irrep T3Off is defined
              if (i.eq.1) then
                T3Off(ii,symi)=daddr(lunt3)-T3IntPoss(iorb)
              end if
!
!3.2.3.2      emulate writing of each block
              length=mapdri(ii,2)
              call ddafile (lunt3,0,dum,length,daddr(lunt3))
!
            end do
          end do
        end do
!
        return
        end
