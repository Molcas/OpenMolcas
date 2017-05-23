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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ADVICE_SIGMA(  IAOCC,  IBOCC,  JAOCC,  JBOCC, LADVICE)
*
* Advice Sigma routine about best route to take
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 :
*           LADVICE = 1 : Business as usual, no transpose of matrix
*                         (resolution on alpha strings, direct exc on beta)
*           LADVICE = 2 = Transpose matrices
*                         (resolution on beta strings, direct exc on alpha)
* (SVC: one call to this routine and ITERM is one, so I removed that
* argument and skipped the checking. Also, the arguments are all scalar,
* so that has been hard-coded now too)
*
* Jeppe Olsen, Tirstrup Airport, Jan 12, 98
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.fh"
#include "gasstr.fh"
#include "orbinp.fh"
#include "cgas.fh"
#include "crun.fh"
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Local Scratch
       DIMENSION ITP(16),JTP(16),KTP(16),LTP(16)
*
      NTEST = 00
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
        CALL SXTYP2_GAS(  NIJTYP,     ITP,     JTP,    NGAS,   IAOCC,
     &                     JAOCC,  IPHGAS)
        CALL SXTYP2_GAS(  NKLTYP,     KTP,     LTP,    NGAS,   IBOCC,
     &                     JBOCC,  IPHGAS)
C?      WRITE(6,*) 'NIJTYP, NKLTYP', NIJTYP,NKLTYP
*. P-h modifications ( I cannot predict these at the moment
        IF(NIJTYP.GE.1.AND.NKLTYP.GE.1) THEN
*
        IF((IPHGAS(ITP(1)).EQ.2.AND.IPHGAS(JTP(1)).EQ.2).OR.
     &     (IPHGAS(KTP(1)).EQ.2.AND.IPHGAS(LTP(1)).EQ.2)     ) THEN
           IPHMODI = 1
         ELSE
           IPHMODI = 0
         END IF
        ELSE
           IPHMODI = 0
        END IF

*
        IF(IPHMODI.EQ.1.OR.NIJTYP.NE.1.OR.NKLTYP.NE.1
     &     .OR.IADVICE.EQ.0) THEN
*. Several connections, i.e. the alpha or the beta blocks are identical,
*. or ph modifications
*. just continue
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
*. going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
          XNJOB = dble(NOBPT(JTP(1)))
          XNJEL = dble(JAOCC(JTP(1)))
          XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1.0D0)
*. Number of kl excitations per beta string :
          XNKLSX = dble((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)))
*. Number of ops (relative to dim of C)
          XNIOB = dble(NOBPT(ITP(1)))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = dble(NOBPT(LTP(1)))
          XNLEL = dble(JBOCC(LTP(1)))
          XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1.0D0)
*. Number of ij excitations per alpha string :
          XNIJSX = dble((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)))
*. Number of ops (relative to dim of C)
          XNKOB = dble(NOBPT(KTP(1)))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8D0*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(LHOLEA+LHOLEB.EQ.MNHL.AND.
     &       (IPHGAS(JTP(1)).EQ.2.OR.IPHGAS(LTP(1)).EQ.2))  THEN
*
             IF(IPHGAS(JTP(1)).EQ.2) THEN
              KHOLEA = LHOLEA-1
              KHOLEB = LHOLEB
             ELSE
              KHOLEA = LHOLEA
              KHOLEB = LHOLEB - 1
             END IF
*
             IF(KHOLEA.EQ.KHOLEB) THEN
               LLADVICE = LADVICE
             ELSE IF(KHOLEA.LT.KHOLEB) THEN
               LLADVICE= 1
             ELSE
               LLADVICE = 2
             END IF
             IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
               WRITE(6,*) ' Advice changed by hole considetions'
               WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
             END IF
             LADVICE = LLADVICE
          END IF
*
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' ADVICE active '
            WRITE(6,*) ' IAOCC JAOCC IBOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
            WRITE(6,*) ' XFLOPA,XFLOPB', XFLOPA,XFLOPB
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
        END IF
*       ^ End if several types/ph modi
C     WRITE(6,*) ' MEMCHECK at end of ADVICE'
C     CALL MEMCHK
C     WRITE(6,*) ' MEMCHECK passed '
      RETURN
      END
