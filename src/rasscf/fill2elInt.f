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
      SUBROUTINE Fill2elInt(TUVX,nacpr2)
*
*     Purpose: Generate 2-electron integral entries and their index in FCIDUMP file for NECI.
*              This should be working with Choleski as well.
*
* Elements (tu|vx) are in triangular order over (t,u), (v,x) and (tu,vx).
* So that index idx = ijkl will lead to (ij,kl) via
*                              (ij-1)*ij/2 < ijkl <= ij*(ij+1)/2
* which will provide solution for
*                              ij = ceiling((-0.5d0+sqrt(0.25+2*ijkl)))
* I removed the 0.25 inside the sqrt to make sure that we do not hit numerical instabilities for diagonal terms
*                              ij = ceiling(-0.5d0+sqrt(2.0d0*ijkl))
*                              kl = ijkl - (ij-1)*ij/2
* From ij and kl indices i, j, k and l are obtained with same argument.
*
#include "fciqmc_global.fh"
#include "trafo_fciqmc.fh"
#include "files_fciqmc.fh"
#include "fciqmc.fh"
* Declaring explicitely the needed ones...
      Real*8 TUVX(nacpr2),CPT,CPE,TIOT,TIOE
      integer i, ijidx, klidx, iorb, jorb, korb, lorb
*
      Call qEnter('Fill2elInt')
*
*     Set time at start of transformation
      CALL SETTIM
*
* loop over the entire set of tuvx integrals...
c      write(6,*) 'writing 2-el integrals in FCIDUMP'
      do i = 1, nacpr2
* compute ij and kl i = ijkl
        ijidx = ceiling(-0.5d0+sqrt(2.0d0*i))
        klidx = i - (ijidx-1)*ijidx/2
*
* now we are going to compute i and j from ij....
        iorb = ceiling(-0.5d0+sqrt(2.0d0*ijidx))
        jorb = ijidx - (iorb-1)*iorb/2
* ...and k and l from kl...
        korb = ceiling(-0.5d0+sqrt(2.0d0*klidx))
        lorb = klidx - (korb-1)*korb/2
        if(abs(TUVX(i)).ge.1.0d-11)then
            write(LuFCI,'(1X,G20.11,4I5)') TUVX(i),
     &                                      iorb,jorb,korb,lorb
c            write(6,'(1X,G20.11,4I5)')      TUVX(i),
c     &                                      iorb,jorb,korb,lorb
        endif
      enddo

*     Set time at end of 2-electron integrals entries process ....
      CALL TIMING(CPT,CPE,TIOT,TIOE)
      If (iPrint.GE.5) WRITE(6,2200) CPT,TIOT
2200  FORMAT(/6X,' TOTAL CPU TIME(SEC)',F8.2,'TOTAL I/O TIME(SEC)',F8.2)
*     exits
      Call qExit('Fill2elInt')
      RETURN
      END
