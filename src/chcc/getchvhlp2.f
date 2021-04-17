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
        subroutine GetChVHlp2 (L2Status,NL2,used,kam)
c
c       this routine do:
c       Define, which one of unused L2 arrays allocated in memory
c       will be used for placing newly needed L2x
c       Criteria of present approach:
c       1) First those possitions that are allocated, but never used yet
c       2) if no 1), take first unused in this step, according L2Status
c          N.B. 2) druhy krok moze byt eventuelne vylepseny, zatial
c               takto odflaknute
c
c       description ov variables:
c       L2Status  - L2 Status array (I)
c       NL2       - number of L2 arrays reserver in memory (I)
c       used      - array indicating, which L2x are already used
c                   in this step (those cannot be touched) (I)
c       kam       - index of possition, where new L2 can be placed (O)
c
        implicit none
        integer L2Status(1:4,1:3)
        integer used(1:4)
        integer NL2,kam
c
c       help variables
        integer i
c
c
c1      search, if there are never used possitions
c
        do i=1,NL2
          if (L2Status(i,1).eq.0) then
          kam=i
          return
          end if
        end do
c
c2      find first unused in this step
c
        do i=1,NL2
          if (used(i).eq.0) then
          kam=i
          return
          end if
        end do
c
c       Jaj nieje dobre ak sme sa dostali az sem
        write (6,*) ' Sorry fish getChVHlp2 '
        call Abend()
c
        return
        end
