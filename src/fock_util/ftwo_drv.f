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
      SUBROUTINE FTwo_Drv(nSym,nBas,nAsh,nSkipX,
     &                    DI,D1A,FA,nTot1,
     &                    ExFac,nTot2,nBMX,CMO)


      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "real.fh"

      Integer nBas(8), nAsh(8), nSkipX(8)
      Dimension CMO(*) , D1A(*) , DI(*), FA(*)
      Logical DoCholesky,REORD,DECO
      Integer ALGO
      COMMON /CHORAS/ REORD,DECO,ALGO


      Call DecideOncholesky(DoCholesky)

      IF (DoCholesky.and.ALGO.eq.2)THEN
*
* Building of the Fock matrix directly from Cholesky
* vectors
*
         Call GetMem('LWFSQ','Allo','Real',LWFSQ,nTot2)
         call dcopy_(nTot2,Zero,0,Work(LWFSQ),1)

         Call Allocate_Work(ipTemp,nTot1)
         Call FZero(Work(ipTemp),nTot1)
*
         CALL CHORAS_DRV(nSym,nBas,nAsh,D1A,DI,Work(ipTemp),
     &                   ExFac,LWFSQ,CMO)

         Call DaXpY_(nTot1,One,Work(ipTemp),1,FA,1)
*
         Call Free_Work(ipTemp)
         Call GetMem('LWFSQ','Free','Real',LWFSQ,nTot2)

      ELSE

*
* Standard building of the Fock matrix from Two-el integrals
*
         Call FockTwo_Drv(nSym,nBas,nAsh,nSkipX,
     &                    DI,D1A,FA,nTot1,
     &                    ExFac,nTot2,nBMX)

      ENDIF


      RETURN
      END
