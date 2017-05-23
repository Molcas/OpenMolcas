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
* Copyright (C) Chad E. Hoyer                                          *
************************************************************************
      SUBROUTINE DQVDiabat(PROP)
************************************************************
*
*   <DOC>
*     <Name>DQVDiabat</Name>
*     <Syntax>Call DQVDiabat(PROP)</Syntax>
*     <Arguments>
*       \Argument{PROP}{Properties computed in RASSI}
*          {Real*8 array}{in}
*     </Arguments>
*     <Purpose> Compute diabats with DQV</Purpose>
*     <Dependencies>RASSI</Dependencies>
*     <Author>C. E. Hoyer (CEH)</Author>
*     <Modified_by>C.E. Hoyer</Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        This subroutine takes in properties that have been
*       computed by RASSI and uses them to compute diabatic
*       states and thus diabatic energies and couplings.
*       Currently, the user must compute x, y, z, xx,
*       yy, zz, and 1/r, and they must be computed in that
*       order.
*     </Description>
*    </DOC>
*
************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DQVDiabat')
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "tshcntrl.fh"
#include "constants2.fh"
      REAL*8 PROP(NSTATE,NSTATE,NPROP)
      REAL*8 TROT(NSTATE,NSTATE)
      REAL*8 TRQ(NSTATE,NSTATE)
      REAL*8 TROTT(NSTATE,NSTATE)
      REAL*8 HDIA(NSTATE,NSTATE)
      REAL*8 HDIAI(NSTATE,NSTATE)
      REAL*8 HAMT(NSTATE,NSTATE)

*  These are the blocks of parameters
      INTEGER, PARAMETER :: MAX=50
      REAL*8, PARAMETER :: MTE=1.0D-8, MTF=1.0D-14
      REAL*8, PARAMETER :: THRS=1.0D-8
      INTEGER :: PNUM(7)

*


      CALL QENTER(ROUTINE)

*
* Printing some stuff
*
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,'(6X,100A1)') ('*',i=1,100)
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,A,33X,A,34X,A)')
     &     '*',' The DQV Diabatization Section ','*'
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,100A1)') ('*',i=1,100)
      WRITE(6,*)
      WRITE(6,*)

      Call CollapseOutput(1,'DQV Diabatization section')
      WRITE(6,'(3X,A)')     '-------------------------'

* Find the properties we need
      DO IPROP=1,7
        PNUM(IPROP)=0
      END DO
      DO IPROP=1,NPROP
        IF (PNAME(IPROP).EQ.'MLTPL  1') THEN
          IF (ICOMP(IPROP).EQ.1) PNUM(1)=IPROP
          IF (ICOMP(IPROP).EQ.2) PNUM(2)=IPROP
          IF (ICOMP(IPROP).EQ.3) PNUM(3)=IPROP
        END IF
        IF (PNAME(IPROP).EQ.'MLTPL  2') THEN
          IF (ICOMP(IPROP).EQ.1) PNUM(4)=IPROP
          IF (ICOMP(IPROP).EQ.4) PNUM(5)=IPROP
          IF (ICOMP(IPROP).EQ.6) PNUM(6)=IPROP
        END IF
        IF (PNAME(IPROP).EQ.'EF0    1') THEN
          IF (ICOMP(IPROP).EQ.1) PNUM(7)=IPROP
        END IF
      END DO
      DO IPROP=1,7
        IF (PNUM(IPROP).EQ.0) THEN
          WRITE(6,*) 'DQVDiabat: '//
     &               'Some required properties are not available'
          CALL Abend()
        END IF
      END DO

* CEH Now we maximize f_DQV through a Jacobi sweep algorithm


* Initialize rotation matrix to unit matrix (theta=0)
        DO ISTA=1, NSTATE
         DO JSTA=1, NSTATE
          TROT(ISTA,JSTA) = 0.0D0
          IF (ISTA .eq. JSTA) THEN
           TROT(ISTA,JSTA) = 1.0D0
          END IF
         END DO
       END DO

*Now, we have to form the trace of the quadrupole
*The user has to specify dipole then quadrupole
       DO ISTA=1, NSTATE
        DO JSTA=1, NSTATE
         TRQ(ISTA,JSTA) = PROP(ISTA,JSTA,PNUM(4)) +
     &                    PROP(ISTA,JSTA,PNUM(5)) +
     &                    PROP(ISTA,JSTA,PNUM(6))
        END DO
       END DO

       Call RecPrt('The TRQ matrix in DQVDiabat','',TRQ,NSTATE,NSTATE)
       WRITE(6,*) ''

       DO i=1, MAX
        THRSCH = 0.0D0

        DO ISTA=1, NSTATE
         DO JSTA=ISTA+1, NSTATE

* Equations 9 and 10 of Kleir et al. with \mu written
* as a vector in addition to the quadrupole trace

          ATERM = PROP(ISTA,JSTA,PNUM(1))**2
     &           +PROP(ISTA,JSTA,PNUM(2))**2
     &           +PROP(ISTA,JSTA,PNUM(3))**2
     &           +ALPHZ*TRQ(ISTA,JSTA)**2
     &           +BETAE*PROP(ISTA,JSTA,PNUM(7))**2
     &           -0.25D0*(PROP(ISTA,ISTA,PNUM(1))**2
     &           +PROP(ISTA,ISTA,PNUM(2))**2
     &           +PROP(ISTA,ISTA,PNUM(3))**2
     &           +ALPHZ*TRQ(ISTA,ISTA)**2
     &           +BETAE*PROP(ISTA,ISTA,PNUM(7))**2)
     &           -0.25D0*(PROP(JSTA,JSTA,PNUM(1))**2
     &           +PROP(JSTA,JSTA,PNUM(2))**2
     &           +PROP(JSTA,JSTA,PNUM(3))**2
     &           +ALPHZ*TRQ(JSTA,JSTA)**2
     &           +BETAE*PROP(JSTA,JSTA,PNUM(7))**2)
     &           +0.5D0*
     &           (PROP(ISTA,ISTA,PNUM(1))*PROP(JSTA,JSTA,PNUM(1))
     &           +PROP(ISTA,ISTA,PNUM(2))*PROP(JSTA,JSTA,PNUM(2))
     &           +PROP(ISTA,ISTA,PNUM(3))*PROP(JSTA,JSTA,PNUM(3))
     &           +ALPHZ*TRQ(ISTA,ISTA)*TRQ(JSTA,JSTA)
     &           +BETAE*PROP(ISTA,ISTA,PNUM(7))*PROP(JSTA,JSTA,PNUM(7)))


          BTERM = PROP(ISTA,JSTA,PNUM(1))*(PROP(ISTA,ISTA,PNUM(1))-
     &            PROP(JSTA,JSTA,PNUM(1))) +
     &            PROP(ISTA,JSTA,PNUM(2))*(PROP(ISTA,ISTA,PNUM(2))-
     &            PROP(JSTA,JSTA,PNUM(2))) +
     &            PROP(ISTA,JSTA,PNUM(3))*(PROP(ISTA,ISTA,PNUM(3))-
     &            PROP(JSTA,JSTA,PNUM(3))) +
     &            ALPHZ*TRQ(ISTA,JSTA)*(TRQ(ISTA,ISTA)-TRQ(JSTA,JSTA))
     &            +BETAE*
     &            PROP(ISTA,JSTA,PNUM(7))*(PROP(ISTA,ISTA,PNUM(7))-
     &            PROP(JSTA,JSTA,PNUM(7)))

          IF (ABS(BTERM) .gt. MTF) THEN
           IF((ATERM .gt. MTE) .or. (ATERM .lt. -1.0D0*MTE)) THEN
* This part of the code computes the rotation angle.
            CTERM=-1.0D0*BTERM/ATERM
            ROTANGF = Atan(CTERM)

            IF((ROTANGF .lt. 0.0D0) .and. (BTERM .gt. 0.0D0)) THEN
             ROTANGF = ROTANGF + RPI
            END IF

            IF((ROTANGF .gt. 0.0D0) .and. (BTERM .lt. 0.0D0)) THEN
             ROTANGF = ROTANGF - RPI
            END IF

            ROTANGO = ROTANGF/4.0D0

            IF((ROTANGO .gt. RPI/4.0D0) .or.
     &       (ROTANGO .lt. -RPI/4.0D0)) THEN
             WRITE(6,*) 'Quadrants 2 and 3'
            END IF

            COSO = cos(ROTANGO)
            SINO = sin(ROTANGO)

            ELSE

            WRITE(6,*) 'A is pretty close to zero.  Something is fishy.'
             COSO=0.0D0
             SINO=0.0D0

            END IF

*Equation 13 of Kleir et al.
            CHNG = ATERM + SQRT(ATERM**2+BTERM**2)

            IF(CHNG .gt. THRSCH) THEN
             THRSCH = CHNG
            END IF

            IF(CHNG .gt. MTE) THEN


*Update T rotation matrix
            DO k=1, NSTATE
             T1 = TROT(k,ISTA)
             T2 = TROT(k,JSTA)
             TROT(k,ISTA) = COSO*T1 + SINO*T2
             TROT(k,JSTA) = -1.0D0*SINO*T1 + COSO*T2
            END DO

*Rotates transition dipole matrices
*These equations can be found in a book on Jacobi sweeps

            TII = PROP(ISTA,ISTA,PNUM(1))
            TJJ = PROP(JSTA,JSTA,PNUM(1))
            TIJ = PROP(ISTA,JSTA,PNUM(1))
            DO k=1, NSTATE
             T1 = PROP(ISTA,k,PNUM(1))
             T2 = PROP(JSTA,k,PNUM(1))
             PROP(ISTA,k,PNUM(1)) = COSO*T1 + SINO*T2
             PROP(k,ISTA,PNUM(1)) = PROP(ISTA,k,PNUM(1))
             PROP(JSTA,k,PNUM(1)) = COSO*T2 - SINO*T1
             PROP(k,JSTA,PNUM(1)) = PROP(JSTA,k,PNUM(1))
            END DO
            PROP(ISTA,JSTA,PNUM(1)) = (COSO**2-SINO**2)*TIJ +
     &          COSO*SINO*(TJJ-TII)
            PROP(JSTA,ISTA,PNUM(1)) = PROP(ISTA,JSTA,PNUM(1))
            PROP(ISTA,ISTA,PNUM(1)) = COSO**2*TII + SINO**2*TJJ +
     &          2*COSO*SINO*TIJ
            PROP(JSTA,JSTA,PNUM(1)) = SINO**2*TII + COSO**2*TJJ -
     &          2*COSO*SINO*TIJ


            TII = PROP(ISTA,ISTA,PNUM(2))
            TJJ = PROP(JSTA,JSTA,PNUM(2))
            TIJ = PROP(ISTA,JSTA,PNUM(2))
            DO k=1, NSTATE
             T1 = PROP(ISTA,k,PNUM(2))
             T2 = PROP(JSTA,k,PNUM(2))
             PROP(ISTA,k,PNUM(2)) = COSO*T1 + SINO*T2
             PROP(k,ISTA,PNUM(2)) = PROP(ISTA,k,PNUM(2))
             PROP(JSTA,k,PNUM(2)) = COSO*T2 - SINO*T1
             PROP(k,JSTA,PNUM(2)) = PROP(JSTA,k,PNUM(2))
            END DO
            PROP(ISTA,JSTA,PNUM(2)) = (COSO**2-SINO**2)*TIJ +
     &          COSO*SINO*(TJJ-TII)
            PROP(JSTA,ISTA,PNUM(2)) = PROP(ISTA,JSTA,PNUM(2))
            PROP(ISTA,ISTA,PNUM(2)) = COSO**2*TII + SINO**2*TJJ +
     &          2*COSO*SINO*TIJ
            PROP(JSTA,JSTA,PNUM(2)) = SINO**2*TII + COSO**2*TJJ -
     &          2*COSO*SINO*TIJ


            TII = PROP(ISTA,ISTA,PNUM(3))
            TJJ = PROP(JSTA,JSTA,PNUM(3))
            TIJ = PROP(ISTA,JSTA,PNUM(3))
            DO k=1, NSTATE
             T1 = PROP(ISTA,k,PNUM(3))
             T2 = PROP(JSTA,k,PNUM(3))
             PROP(ISTA,k,PNUM(3)) = COSO*T1 + SINO*T2
             PROP(k,ISTA,PNUM(3)) = PROP(ISTA,k,PNUM(3))
             PROP(JSTA,k,PNUM(3)) = COSO*T2 - SINO*T1
             PROP(k,JSTA,PNUM(3)) = PROP(JSTA,k,PNUM(3))
            END DO
            PROP(ISTA,JSTA,PNUM(3)) = (COSO**2-SINO**2)*TIJ +
     &          COSO*SINO*(TJJ-TII)
            PROP(JSTA,ISTA,PNUM(3)) = PROP(ISTA,JSTA,PNUM(3))
            PROP(ISTA,ISTA,PNUM(3)) = COSO**2*TII + SINO**2*TJJ +
     &          2*COSO*SINO*TIJ
            PROP(JSTA,JSTA,PNUM(3)) = SINO**2*TII + COSO**2*TJJ -
     &          2*COSO*SINO*TIJ

*Rotates trace of quadrupole matrix

            TII = TRQ(ISTA,ISTA)
            TJJ = TRQ(JSTA,JSTA)
            TIJ = TRQ(ISTA,JSTA)
            DO k=1, NSTATE
             T1 = TRQ(ISTA,k)
             T2 = TRQ(JSTA,k)
             TRQ(ISTA,k) = COSO*T1 + SINO*T2
             TRQ(k,ISTA) = TRQ(ISTA,k)
             TRQ(JSTA,k) = COSO*T2 - SINO*T1
             TRQ(k,JSTA) = TRQ(JSTA,k)
            END DO
            TRQ(ISTA,JSTA) = (COSO**2-SINO**2)*TIJ + COSO*SINO*(TJJ-TII)
            TRQ(JSTA,ISTA) = TRQ(ISTA,JSTA)
            TRQ(ISTA,ISTA) = COSO**2*TII + SINO**2*TJJ + 2*COSO*SINO*TIJ
            TRQ(JSTA,JSTA) = SINO**2*TII + COSO**2*TJJ - 2*COSO*SINO*TIJ


*Rotates electric potential
            TII = PROP(ISTA,ISTA,PNUM(7))
            TJJ = PROP(JSTA,JSTA,PNUM(7))
            TIJ = PROP(ISTA,JSTA,PNUM(7))
            DO k=1, NSTATE
             T1 = PROP(ISTA,k,PNUM(7))
             T2 = PROP(JSTA,k,PNUM(7))
             PROP(ISTA,k,PNUM(7)) = COSO*T1 + SINO*T2
             PROP(k,ISTA,PNUM(7)) = PROP(ISTA,k,PNUM(7))
             PROP(JSTA,k,PNUM(7)) = COSO*T2 - SINO*T1
             PROP(k,JSTA,PNUM(7)) = PROP(JSTA,k,PNUM(7))
            END DO
            PROP(ISTA,JSTA,PNUM(7)) = (COSO**2-SINO**2)*TIJ +
     &          COSO*SINO*(TJJ-TII)
            PROP(JSTA,ISTA,PNUM(7)) = PROP(ISTA,JSTA,PNUM(7))
            PROP(ISTA,ISTA,PNUM(7)) = COSO**2*TII + SINO**2*TJJ +
     &          2*COSO*SINO*TIJ
            PROP(JSTA,JSTA,PNUM(7)) = SINO**2*TII + COSO**2*TJJ -
     &          2*COSO*SINO*TIJ

            END IF
           END IF
          END DO
         END DO

         IF(THRSCH .lt. THRS) THEN
          WRITE(6,*) 'Converged in this many iterations:'
          WRITE(6,*) i
          WRITE(6,*) ''
          EXIT
         END IF

         IF(i .eq. MAX) THEN
          write(6,*) 'Max iterations'
         END IF
        END DO

        CALL RecPrt('Diabatic Coefficients','',TROT,NSTATE,NSTATE)

        DO ISTA = 1, NSTATE
         DO JSTA = 1, NSTATE
          TROTT(JSTA,ISTA)=TROT(JSTA,ISTA)**2
         END DO
        END DO
        CALL RecPrt('Weights of adiabatic states','',
     &              TROTT,NSTATE,NSTATE)

        DO ISTA = 1, NSTATE
         DO JSTA = 1, NSTATE
          HAMT(ISTA,JSTA)=HAM(ISTA,JSTA)
         END DO
        END DO

        TROTT = TROT


        CALL DGEMM_('T','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &              TROTT, NSTATE, HAMT, NSTATE, 0.0D0,HDIAI,NSTATE)

        CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &              HDIAI, NSTATE, TROTT, NSTATE, 0.0D0,HDIA,NSTATE)

        WRITE(6,*) ''
        WRITE(6,*) 'The eigenvectors may no longer be in'
        WRITE(6,*) 'energetic order'
        WRITE(6,*) ''
        CALL RecPrt('Diabatic Hamiltonian','',HDIA,NSTATE,NSTATE)

        Call CollapseOutput(0,'DQV Diabatization section')
        WRITE(6,*) ''

*Molcas verify calls
       CALL Add_Info('Adiabat1',HAM(1,1),1,4)
       CALL Add_Info('Adiabat2',HAM(2,2),1,4)
       CALL Add_Info('Adiabat3',HAM(3,3),1,4)
       CALL Add_Info('DQVHam11',HDIA(1,1),1,4)
       CALL Add_Info('DQVHam12',HDIA(1,2),1,4)
       CALL Add_Info('DQVHam13',HDIA(1,3),1,4)
       CALL Add_Info('DQVHam22',HDIA(2,2),1,4)
       CALL Add_Info('DQVHam23',HDIA(2,3),1,4)
       CALL Add_Info('DQVHam33',HDIA(3,3),1,4)
*End Molcas verify calls

       CALL QEXIT(ROUTINE)
       RETURN
       END
