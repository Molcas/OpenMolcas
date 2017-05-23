************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2011, Jeppe Olsen                                      *
*               2011, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE CHECK_BLOCKS_FOR_BK_APPROX(
     &           IATP,IBTP,JATP,JBTP,
     &           IASM,IBSM,JASM,JBSM,
     &           IOCTPA,IOCTPB,I_DO_EXACT_BLOCK)
*. Check whether block <IATP, IBTP! H! JATP, JBTP> should
* be calculated exactly or by BK approx
*. Input
* ======
* IATP IBTP JATP JBTP : Supergroups (relative numbers)
* IOCTPA, IOBTPB : Offset for type
*. Output
* ======
* I_DO_EXACT_BLOCK = 1 => Do exact block
*                  = 0 => Set block to zero
*                  =-1 => Use diagonal aproximation
* Giovanni +Jeppe Olsen, Sept 2011, on a bench at Torre Normanna, Sicily
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "bk_approx.fh"
*. Local
      INTEGER IOCC(MXPNGAS), JOCC(MXPNGAS)
      NTEST = 000

      IONE = 1
      IOCTPA=IBSPGPFTP(1)
      IOCTPB=IBSPGPFTP(2)
      if(ntest.ge.100) then
        write(6,*) 'IOCTPA, IOCTPB', IOCTPA, IOCTPB
        write(6,*) 'IATP,IBTP,JATP,JBTP'
        write(6,'(5x,I4,5x,I4,5x,I4,5X,I4)') IATP,IBTP,JATP,JBTP
        write(6,*)  'NELFSPGP (IA), (IB), (JA), (JB)'
        write(6,*) (NELFSPGP(igas,IOCTPA-1+IATP),igas=1,NGAS)
        write(6,*) (NELFSPGP(igas,IOCTPB-1+IBTP),igas=1,NGAS)
        write(6,*) (NELFSPGP(igas,IOCTPA-1+JATP),igas=1,NGAS)
        write(6,*) (NELFSPGP(igas,IOCTPB-1+JBTP),igas=1,NGAS)
      end if
      CALL IVCSUM(IOCC,
     &     NELFSPGP(1,IOCTPA-1+IATP),NELFSPGP(1,IOCTPB-1+IBTP),
     &     IONE,IONE,NGAS)
      CALL IVCSUM(JOCC,
     &     NELFSPGP(1,IOCTPA-1+JATP),NELFSPGP(1,IOCTPB-1+JBTP),
     &     IONE,IONE,NGAS)
        if(ntest.ge.100) then
          write(6,*) ' Routine CHECK_BLOCKS_FOR_BK_APPROX is speaking! '
          write(6,*) ' I am doing BK-type of approximation '
          write(6,*) ' Min and Max for subspace with exact Hamiltonian '
          write(6,*) ' =============================================== '
          write(6,*) 'NGASBK : ',NGASBK
          write(6,*) '              Min. Occ.      Max. Occ.           '
          Do IGAS = 1, NGASBK
            write(6,'(A,I2,10X,I3,9X,I3)')
     &      '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
          End do
         end if
*
      IOCC_IN = ICHECK_OCC_IN_ACCSPC(IOCC,IOCCPSPC,NGASBK,20)
      JOCC_IN = ICHECK_OCC_IN_ACCSPC(JOCC,IOCCPSPC,NGASBK,20)
* = 1 if the Occupation is IN
* = 0 of the Occupation is OUT
*. If both occupation are outside of IOCCPSPC, we make approximations
      IF(IOCC_IN.EQ.0.AND.JOCC_IN.EQ.0) THEN
*. If the blocks are identical use diagonal approximation
        IF(IATP.EQ.JATP.AND.IASM.EQ.JASM.AND.
     &     IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*. Compute as diagonal
           I_DO_EXACT_BLOCK = -1
        ELSE
*. Neglect
           I_DO_EXACT_BLOCK =  0
        END IF
      ELSE
*. Atleast one block is in PSPC, so calculate exactly
        I_DO_EXACT_BLOCK = 1
      END IF
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*)  ' CHECK_BLOCKS_FOR_BK_APPROX is speaking '
       WRITE(6,'(A, 4I4)') ' Input blocks IA, IB, JA, JB = ',
     &                     IATP,IBTP, JATP, JBTP
       WRITE(6,'(A, 4I4)') ' Input blocks IASM, IBSM, JASM, JBSM = ',
     &                     IASM,IBSM, JASM, JBSM
      END IF
*
      IF(NTEST.GE.10) THEN
       WRITE(6,'(A,I4)') ' I_DO_EXACT_BLOCK = ', I_DO_EXACT_BLOCK
      END IF
*
      RETURN
      END
