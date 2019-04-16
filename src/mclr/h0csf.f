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
* Copyright (C) 1993, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE H0CSF(H0,IPQCSF,IPQCNF,MXP1DM,MXP2DM,MXQDM,
     &                  DTOC,IPRODT,ICONF,
     &                  IREFSM,ECORE,NINOB,NACTOB,
     &                  SCR,ISCR,NCONF,NEL,NAEL,NBEL,IPWAY,
     &                  NP1CSF,NP1CNF,NP2CSF,NP2CNF,
     &                  NQCSF,NQCNF,NPQCSF,NPQCNF,
     &                  DIAG,DIAGCN,
     &                  NTEST,INTSPC,ICOMBI,PSSIGN)
*
* Obtain H0 subspace defined by the three parameters
* MXP1DM,MXP2DM,MXQDM and obtain
* explicit representation of hamilton matrix in subspace
*
*
* The H0 space consist of three subspaces
* P1,P2 AND Q
* The H0 matrix can be pictorized as
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* The exact Hamiltonian is therefore calculated in subspace P1 and P2
* but only the interaction between a larger space Q and subspace P1
* is calculated exactly
*
* Lucia Version, September 1993, Jeppe Olsen
*
* ===========
* ARGUMENTS :
* ===========
*   H0 : combined block to contain :
*      PHP : hamilton matrix in subspace P1+P2 output)
*      PHQ : Hamiltonian matrix block with row in P1 and column in Q
*      QHQ : Diagonal approximation to matrix in Q-Q space
*   IPQCSF : CSFs defining subspace (output)
*   IPQCNF : Configurations defining subspace (output )
*   MXP1DM : Largest allowed dimension of subspace P1 (Input)
*   MXP2DM : Largest allowed dimension of subspace P2 (Input)
*   MXQDM  : Largest allowed dimension of subspace Q  (Input)
*   DTOC : Transformation matrix between CSFs and DETs (input)
*   IPRODT : Prototype determinants (input)
*   ICONF : List of configurations  (input)
*   IREFSM : symmetry of considered CI space (input)
*   ECORE : Core energy (input)
*   NINOB : Number of inactive orbitals(input)
*   NACTOB : Number of active orbitals (input)
*   SCR    : Scratch array of length ????
*   NCONF : Number of configurations of symmetry IREFSM
*   IPWAY : Defines way of choosing Primary space
*    = 1  : use the first configurations (untill atmost
*           MXPDIM CSFs have been included )
*   NP1CNF : Number of primary configurations obtained in P1(output)
*   NP1CSF : Number of primary CSFs obtained in P1 (OUTPUT)
*   NP2CNF : Number of primary configurations obtained in P2 (output)
*   NP2CSF : Number of primary CSFs obtained in P2 (OUTPUT)
*   NQCNF  : Number of primary configurations obtained in Q (output)
*   NQCSF  : Number of primary CSFs obtained in Q (OUTPUT)
*
*   DIAG   : Hamilton diagonal over CSFs ( INPUT )
*   DIAGCN : space for diagonal over configurations
*   INTSPC : Internal space number of actual expansion
*   ICOMBI : = 0 => no spin combinations
*            = 1 =>    spin combinations
*   PSSIGN : Spin combination sign
*
* =========================================
* Jeppe Olsen , Spring of 90 , from PHPCSF
* =========================================
*  Lucia Version, September 1993
* =========================================
      IMPLICIT REAL*8 (A-H,O-Z)
*.Output
      DIMENSION H0(*),IPQCSF(*),IPQCNF(*)
*.Input
      DIMENSION DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
*.Scratch space
      DIMENSION SCR(*),ISCR(*),DIAGCN(*)
*. SCR and ISCR are supposed to refer to the same array
*
#include "detdim.fh"

#include "SysDef.fh"
#include "spinfo_mclr.fh"
*
      CALL H0CSF_INTERNAL(SCR,DIAGCN)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE H0CSF_INTERNAL(SCR,DIAGCN)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*),DIAGCN(*)
      INTEGER, POINTER :: iPTR(:)
*
** 1 : Obtain primary subspace
*
      NP1CSF = 0
      NP1CNF = 0
      NP2CSF = 0
      NP2CNF = 0
      NPCSF  = 0
      NQCSF  = 0
      NQCNF  = 0
      NPQCSF = 0
      NPQCNF = 0
      ICSFMN = 0 ! dummy initialize
      MXPQDM = MXP1DM + MXP2DM + MXQDM
*
*
* ===================================================================== *
* 1 :                 Generate initial subspace
* ===================================================================== *
*
      IF(IPWAY .EQ. 1 ) THEN
*
* Just use the first CSFs as subspace
*
        ICNF = 0
        ICSF = 0
        DO 100 ITYP = 1, NTYP
          NJCNF = NCNATS(ITYP,IREFSM)
          NIRREP = NCPCNT(ITYP)
          DO 90 IICNF = 1, NJCNF
            ICNF = ICNF + 1
            IF(NP1CSF+NIRREP .LE. MXP1DM) THEN
              NPQCNF = NPQCNF + 1
              NP1CNF = NP1CNF + 1
              IPQCNF(NPQCNF) = ICNF
              DO 80 IICSF = 1, NIRREP
                NPQCSF = NPQCSF + 1
                NP1CSF = NP1CSF + 1
                IPQCSF(NPQCSF) = NPQCSF
   80         CONTINUE
            ELSE IF( NP2CSF+NIRREP .LE. MXP2DM ) THEN
              NPQCNF = NPQCNF + 1
              NP2CNF = NP2CNF + 1
              IPQCNF(NPQCNF) = ICNF
              DO 180 IICSF = 1, NIRREP
                NPQCSF = NPQCSF + 1
                NP2CSF = NP2CSF + 1
                IPQCSF(NPQCSF) = NPQCSF
  180         CONTINUE
            ELSE IF( NQCSF+NIRREP .LE. MXQDM ) THEN
              NPQCNF = NPQCNF + 1
              NQCNF = NQCNF + 1
              IPQCNF(NPQCNF) = ICNF
              DO 280 IICSF = 1, NIRREP
                NPQCSF = NPQCSF + 1
                NQCSF = NQCSF + 1
                IPQCSF(NPQCSF) = NPQCSF
  280         CONTINUE
            ELSE
              GOTO 101
            END IF
   90     CONTINUE
  100   CONTINUE
  101   CONTINUE
*
      ELSE IF( IPWAY .EQ. 2 ) THEN
*          Obtain lowest CSFs
*
* 1 :  local Diagonal elements over configurations
*
        IICSF = 1
        IICNF = 1
*
        klFREE=1
        KLDIPQ = 1
        KLFREE = KLFREE + MXPQDM
*
        KLCONF = KLFREE
        KLFREE = KLCONF + NEL
*

        DO 300 ITYP = 1, NTYP
          NJCNF = NCNATS(ITYP,IREFSM)
          NIRREP = NCPCNT(ITYP)
          DO 290 ICNF = 1, NJCNF
            DIAGCN(IICNF) = DIAG(IICSF)
            IICNF = IICNF + 1
            IICSF = IICSF + NIRREP
  290     CONTINUE
  300   CONTINUE
* Largest element
        XMAX = FNDMNX(DIAGCN,NCONF,2)
* loop over lowest configurations
*
        IFINIT = 0
  400   CONTINUE
*
        XMIN = XMAX + 1.0D0
        IMIN = 0
*
        IICNF = 1
        ICSFOF = 1
        DO 500 ITYP = 1, NTYP
          NIRREP = NCPCNT(ITYP)
          DO 450 ICNF = 1,NCNATS(ITYP,IREFSM)
            IF(DIAGCN(     IICNF).LT.XMIN) THEN
              XMIN = DIAGCN(     IICNF)
              IMIN = IICNF
              ICSFMN = ICSFOF
              NCSFMN = NIRREP
            END IF
*
            IICNF = IICNF + 1
            ICSFOF = ICSFOF + NIRREP
  450     CONTINUE
  500   CONTINUE
*
** Next lowest element has been found
*
        IF(NPQCSF + NCSFMN .LE. MXPQDM ) THEN
*.1       add new configuration
          NPQCNF = NPQCNF + 1
          IPQCNF(NPQCNF) = IMIN
          SCR(KLDIPQ-1+NPQCNF) = XMIN
          CALL ISTVC2(IPQCSF(NPQCSF+1),ICSFMN-1,1,NCSFMN)
          NPQCSF = NPQCSF + NCSFMN

*. Mask
          DIAGCN(     IMIN) = XMAX + 1.0D0
        ELSE
          IFINIT = 1
*.2       No space for this configuration , remove previous
*         configurations with the same diagonal value
          IICNF = NPQCNF+1
  600     CONTINUE
            IICNF = IICNF - 1
            DIAVAL = DIAGCN(     IPQCNF(IICNF) )
            IF(ABS(DIAVAL-XMIN) .LE. 1.0D-10) THEN
              NPQCNF = NPQCNF -1
              CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iPTR,[1])
              CALL GETCNF(iPTR,ITYP,IPQCNF(IICNF),
     &        ICONF,IREFSM,NEL,NTEST)
              NULLIFY(iPTR)
              NPQCSF = NPQCSF - NCPCNT(ITYP)
              GOTO 600
            END IF
        END IF
        IF(IFINIT.EQ.0.AND.NPQCNF.LT.NCONF) GOTO 400
*. NPQCSF has now been collected, obtain P1,P2 and Q space
*. so that degenerate configurations are not  in
*. different  subspaces
*
*. Arrange selected configurations in degenerate pairs
*
        KLFREI = rtoi*(KLFREE-1) + 1
*
        KLIDEG = KLFREI
        KLFREE = KLFREE + NPQCNF
*
        KLCONF = KLFREE
        KLFREE = KLFREE + NEL
*
        CALL DEGVEC(SCR(KLDIPQ),NPQCNF,NDGVL,ISCR(KLIDEG))
*. Number of configurations in P1,P2,Q
        ICNF = 0
        DO 800 IDEG = 1, NDGVL
*.Number of CSFs in this group of degenerate values
          IDGVL = ISCR(KLIDEG-1+IDEG)
          IDGCSF = 0
          DO 780 IDGCNF = 1, IDGVL
            ICNF = ICNF + 1
              CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iPTR,[1])
              CALL GETCNF(iPTR,ITYP,IPQCNF(ICNF),
     &        ICONF,IREFSM,NEL,NTEST)
              NULLIFY(iPTR)
              IDGCSF = IDGCSF + NCPCNT(ITYP)
  780    CONTINUE
         IF(NP1CSF+IDGCSF .LE. MXP1DM. AND. NP2CSF+NQCSF.EQ.0) THEN
*. Add to P1
           NP1CSF = NP1CSF + IDGCSF
           NP1CNF = NP1CNF + IDGVL
         ELSE IF(NP2CSF+IDGCSF .LE. MXP2DM .AND. NQCSF .EQ. 0 ) THEN
*. Add to P2
           NP2CSF = NP2CSF + IDGCSF
           NP2CNF = NP2CNF + IDGVL
         ELSE IF( NQCSF+IDGCSF .LE. MXQDM ) THEN
*. Add to Q
           NQCSF = NQCSF + IDGCSF
           NQCNF = NQCNF + IDGVL
         ELSE
*. No space for configuration so
           GOTO 801
         END IF
  800  CONTINUE
  801  CONTINUE
*
       NPCSF = NP1CSF + NP2CSF
       NPCNF = NP1CNF + NP2CNF
*
       NPQCSF = NPCSF + NQCSF
       NPQCNF = NPCNF + NQCNF
*
      END IF
*. End if for IWAY = 2
*
*. This is not beautiful, but neccessary
      MXP1DM = NP1CSF
      MXP2DM = NP2CSF
      MXQDM = NQCSF
*
*
* ============================================================
* 2          Construct Hamiltonian matrix in subspace
* ============================================================
*
*. Do not add core energy to subspace Hamiltonian, add to eigenvalues
*. Pointers in H0
      KLPHP  = 1
      KLPHQ = KLPHP + NPCSF*(NPCSF+1)/2
      KLQHQ = KLPHQ + NP1CSF*NQCSF
*
*.PHP matrix
*
*     CALL C_F_POINTER(C_LOC(DIAGCN(1)),iPTR,[1])
      CALL CNHCNM(H0(KLPHP),1,IPQCNF,NPCNF,IPQCNF,NPCNF,NPCSF,NPCSF,
     &     DIAGCN,ICONF,NEL,IREFSM,NAEL,NBEL,NINOB,NACTOB,ECORE,
     &     IPRODT,DTOC,INTSPC,ICOMBI,PSSIGN,NTEST)
*
*. PHQ matrix
*
      CALL CNHCNM(H0(KLPHQ),0,IPQCNF,NP1CNF,IPQCNF(1+NPCNF),NQCNF,
     &   NP1CSF,NQCSF,DIAGCN,ICONF,NEL,IREFSM,NAEL,NBEL,NINOB,NACTOB,
     &   ECORE,IPRODT,DTOC,INTSPC,ICOMBI,PSSIGN,NTEST)
*     NULLIFY(iPTR)
*
      RETURN
      END SUBROUTINE H0CSF_INTERNAL
*
      END
