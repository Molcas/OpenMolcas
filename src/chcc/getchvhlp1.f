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
        subroutine GetChVHlp1 (cGrp,deGrp,yes,NL2,L2Status)
c
c       this routine do:
c       check, if Cholesky vector block cGrp,deGrp is actually situated
c       in the memory as one of L2x
c
c       descrition of parameters:
c       cGrp, deGrp - required Group of c, delta in L2(m,c',de') (I)
c       yes    -  0 - not loaded  (O)
c                 x - loaded, value of L2x (x=1-4)
c       NL2         - Number of L2 arrays really reserved in memory (I)
c       L2Status    - L2 status matrix (I)
c                     L2Status(i,1) - cGrp
c                     L2Status(i,2) - deGrp
c                     L2Status(i,3) - Poss
c
        implicit none
        integer cGrp,deGrp,yes,NL2
        integer L2Status(1:4,1:3)
c
c       help variables
        integer i
c
        yes=0
        do i=1,NL2
        if ((cGrp.eq.L2Status(i,1)).and.(deGrp.eq.L2Status(i,2))) then
        yes=i
        end if
        end do
c
        return
        end
