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
       subroutine t3reorg (wrk,wrksize,                                 &
     & noa,nsym)
!
!     this routine do final reorganization of t3nam file
!     and produce final form of this file
!     as it will be required in T3 and close t3nam file
!
!     noa   - array with occupation numbers
!     nsym  - actual number of irreps
!
       implicit none
#include "wrk.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       integer noa(1:8)
       integer nsym
!
!     help variables
!
       integer length,iri,possri
       integer posst
       integer symi,i,iaddr,iindex,iPossPack
!
!*        def iPossPack
!           iPossPack - possition of (maps+Ri) set in packed
!                        (i.e. final) of T3nam file
        iPossPack=T3IntPoss(1)
!
        iindex=0
        do symi=1,nsym
!
!0      get map's of R_i(a,bc)
        call ccsort_t3grc0                                              &
     &       (3,8,4,4,4,0,symi,possri0,posst,mapdri,mapiri)
!
        do i=1,noa(symi)
        iindex=iindex+1
!
!1        reconstruct R_i(a,bc) per blocks as in is
!         actually written in t3man file
          do iri=1,mapdri(0,5)
!
!1.1        iind address of this R_i block in t3nam file
            iaddr=T3IntPoss(iindex)+T3Off(iri,symi)
!
!1.2        def possition of of this block in R1
            possri=mapdri(iri,1)
!
!1.3        read integrals into proper possition
            length=mapdri(iri,2)
            if (length.gt.0) then
            call ddafile (lunt3,2,wrk(possri),length,iaddr)
            end if
!
          end do
!
!2          write into t3nam file in packed form
!            1) mapdri, mapiri
!           2) R_i
!2.1          def final (packed) address for i-th set (maps+Ri)
          T3intPoss(iindex)=iPossPack
          iaddr=T3intPoss(iindex)
!
!2.2      write maps
          call idafile (lunt3,1,mapdri,3078,iaddr)
          call idafile (lunt3,1,mapiri,512,iaddr)
!
!2.3      def actual length of Ri
          length=0
          do iri=1,mapdri(0,5)
          length=length+mapdri(iri,2)
          end do
!         length=mapdri(iri,1)+mapdri(iri,2)-mapdri(1,1)
!
!2.4          write Ri as one block
          call ddafile (lunt3,1,wrk(possri0),length,iaddr)
!
!2.5          save updated address as a new packed (final) possition
!         for next i
          iPossPack=iaddr
!
        end do
        end do
!
!3        store new packed (final) addreses T3IntPoss in t3nam file
!       (at the beggining)
        iaddr=0
        call idafile (lunt3,1,T3IntPoss,mbas,iaddr)
!
!4        close t3nam file
        call daclos (lunt3)
!
       return
       end
