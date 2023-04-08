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
        subroutine GetChVHlp2 (L2Status,NL2,used,kam)
!
!       this routine do:
!       Define, which one of unused L2 arrays allocated in memory
!       will be used for placing newly needed L2x
!       Criteria of present approach:
!       1) First those possitions that are allocated, but never used yet
!       2) if no 1), take first unused in this step, according L2Status
!          N.B. 2) druhy krok moze byt eventuelne vylepseny, zatial
!               takto odflaknute
!
!       description ov variables:
!       L2Status  - L2 Status array (I)
!       NL2       - number of L2 arrays reserver in memory (I)
!       used      - array indicating, which L2x are already used
!                   in this step (those cannot be touched) (I)
!       kam       - index of possition, where new L2 can be placed (O)
!
        implicit none
        integer L2Status(1:4,1:3)
        integer used(1:4)
        integer NL2,kam
!
!       help variables
        integer i
!
!
!1      search, if there are never used possitions
!
        do i=1,NL2
          if (L2Status(i,1).eq.0) then
          kam=i
          return
          end if
        end do
!
!2      find first unused in this step
!
        do i=1,NL2
          if (used(i).eq.0) then
          kam=i
          return
          end if
        end do
!
!       Jaj nieje dobre ak sme sa dostali az sem
        write (6,*) ' Sorry fish getChVHlp2 '
        call Abend()
!
        return
        end
