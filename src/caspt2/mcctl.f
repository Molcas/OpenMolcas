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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MCCTL(HEFF)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      REAL*8 HEFF(NSTATE,NSTATE)

      INTEGER ISTATE
      REAL*8 DVALUE

      CALL QENTER('MCCTL')

C The ket state is JSTATE.
C Loop over the bra states
      DO ISTATE=1,NSTATE
        IF(ISTATE.EQ.JSTATE) THEN
          HEFF(ISTATE,JSTATE)=E2TOT
        ELSE
C Compute the effective Hamiltonian:
          CALL HEFVAL(ISTATE,JSTATE,DVALUE)
          HEFF(ISTATE,JSTATE)=DVALUE
        END IF
      END DO

      IF ((IPRGLB.GE.VERBOSE).OR.(NLYROOT.NE.0)) THEN
       WRITE(6,*)
       WRITE(6,*) 'Hamiltonian Effective Couplings'
       WRITE(6,*) '-------------------------------'
       WRITE(6,*)
       WRITE(6,'(10X,6X,A3,I4,A3)') ' | ', MSTATE(JSTATE), ' > '
       DO ISTATE=1,NSTATE
        WRITE(6,'(A3,I4,A3,F16.8)') ' < ',MSTATE(ISTATE),' | ',
     &    HEFF(ISTATE,JSTATE)
       ENDDO
      ENDIF

      CALL QEXIT('MCCTL')
      RETURN
      END
