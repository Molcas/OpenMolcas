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
      SUBROUTINE HEFVAL(IST,JST,DVALUE)
      IMPLICIT NONE
C Apart from input call parameters, we need two vectors stored on
C LUSOLV. Vector nr IVECC (presently=2) contains the contravariant
C elements of the solution to the CASPT2 equations.
C IVECW is the number (presently=6) of the vector on LUSOLV
C where a contravariant representation of the CASPT2 Right-Hand Side
C vector is stored. This depends on the MOs used, but is actually
C the same for all the root states.

#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

      INTEGER IST,JST
      REAL*8 DVALUE

      INTEGER I
      INTEGER NTG1,NTG2,NTG3,LTG1,LTG2,LTG3
      INTEGER IDCI,LCI1,LCI2
      REAL*8 OVL,DUMMY(1)

C We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
C Note: Need proper allocation even if unused, sinced allocated
C arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      CALL GETMEM('TG1','ALLO','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','ALLO','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','ALLO','REAL',LTG3,NTG3)
      WORK(LTG1)=0.0D0
      WORK(LTG2)=0.0D0
      WORK(LTG3)=0.0D0

      CALL GETMEM('MCCI1','ALLO','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','ALLO','REAL',LCI2,MXCI)
      IF(ISCF.EQ.0) THEN
C Read root vectors nr. IST and JST from LUCI.
        IDCI=IDTCEX
        DO I=1,NSTATE
          IF(I.EQ.IST) THEN
            CALL DDAFILE(LUCIEX,2,WORK(LCI1),NCONF,IDCI)
            IF(I.EQ.JST) THEN
              CALL DCOPY_(NCONF,WORK(LCI1),1,WORK(LCI2),1)
            END IF
          ELSE IF(I.EQ.JST) THEN
            CALL DDAFILE(LUCIEX,2,WORK(LCI2),NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          END IF
        END DO
      END IF

      CALL MKTG3(STSYM,STSYM,WORK(LCI1),WORK(LCI2),OVL,
     &           WORK(LTG1),WORK(LTG2),NTG3,WORK(LTG3))
      CALL GETMEM('MCCI1','FREE','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','FREE','REAL',LCI2,MXCI)

      CALL HCOUP(IVECW,IVECC,OVL,WORK(LTG1),WORK(LTG2),
     &           WORK(LTG3),DVALUE)

      CALL GETMEM('TG1','FREE','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','FREE','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','FREE','REAL',LTG3,NTG3)

      RETURN
      END
