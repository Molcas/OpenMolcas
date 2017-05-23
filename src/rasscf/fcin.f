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
      SUBROUTINE FCIN_rasscf(FLT,nFLT,DLT,FSQ,DSQ,EMY,CMO)
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
*     Include 'global.inc'
#include "fciqmc_global.fh"
#include "trafo_fciqmc.fh"
#include "WrkSpc.fh"
*
      Real*8 CMO(*)
      DIMENSION DLT(*),FLT(nFLT),DSQ(*),FSQ(*)
*
      Call qEnter('FCIN_rasscf')
*
*     Construct the one-electron density matrix for frozen space
*
      CALL DONEI_rasscf(DLT,DSQ,CMO)
*
*     Compute the one-electron energy contribution to EMY
*
      EONE=0.0D0
      DO 105 NPQ=1,NTOT1
         EONE=EONE+DLT(NPQ)*FLT(NPQ)
105   CONTINUE
      EMY=EONE
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
         WRITE(6,'(6X,A,E20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
      ENDIF
*
*     Quit here if there are no frozen orbitals
*
      NTFRO=0
      n_Bas=0
      DO 33 ISYM=1,NSYM
        NTFRO=NTFRO+NFRO(ISYM)
        n_Bas=Max(n_Bas,nBas(iSym))
33    CONTINUE
      IF ( NTFRO.EQ.0 ) THEN
        Call qExit('FCIN_rasscf')
        RETURN
      ENDIF
*
*     Compute the two-electron contribution to the Fock matrix
*
      Call Allocate_Work(ipTemp,nFlt)
      Call FZero(Work(ipTemp),nFlt)
*
*
*
*
      CALL GETMEM('FCIN2','ALLO','REAL',LW2,n_Bas**2)
      Call FZero(Work(LW2),n_Bas**2)
      CALL GETMEM('FCIN1','MAX','REAL',LW1,LBUF)
      LBUF = Max(LBUF-LBUF/10,0)
      CALL GETMEM('FCIN1','ALLO','REAL',LW1,LBUF)
      Call FZero(Work(LW1),LBUF)

      CALL FTWOI_rasscf(DLT,DSQ,Work(ipTemp),nFlt,FSQ,LBUF,WORK(LW1),
     &                                                     WORK(LW2))

      CALL GETMEM('FCIN1','FREE','REAL',LW1,LBUF)
      CALL GETMEM('FCIN2','FREE','REAL',LW2,n_Bas**2)

      Call DaXpY_(nFlt,1.0D+0,Work(ipTemp),1,Flt,1)
      Call Free_Work(ipTemp)
*
*     Add the two-electron contribution to EMY
      ETWO=-EONE
      DO 110 NPQ=1,NTOT1
         ETWO=ETWO+DLT(NPQ)*FLT(NPQ)
110   CONTINUE
      EMY=EONE+0.5D0*ETWO
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
         WRITE(6,'(6X,A,E20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
      ENDIF
*
*     Exit
*
      Call qExit('FCIN_rasscf')
      RETURN
      END
