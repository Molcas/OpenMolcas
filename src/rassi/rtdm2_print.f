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
* Copyright (C) 2020, Bruno Tenorio                                    *
************************************************************************

C Print the reduced 2-e TDM in ASCII format.
C Code adapted from trd_print.f written by P. A. Malmqvist.
      SUBROUTINE RTDM2_PRINT(ISTATE, JSTATE, EIJ, NDYSAB, DYSAB,
     &                 NRT2MAB , RT2M , CMO1, CMO2, AUGSPIN)

      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "stdalloc.fh"
      INTEGER ISTATE, JSTATE, SYM12
      INTEGER NDYSAB,NRT2MAB,AUGSPIN
      Real*8  DYSAB(*), RT2M(*), CMO1(*), CMO2(*)
      Real*8  EIJ
C -------------------------------------------------------------
C The spin coupling matrix elements have the following index-code:
C             SPIN=1 means  K2V (AAB+BBB)
C             SPIN=-1 means SDA (AAA+BBA)
C Notice, SPIN here has nothing to do with the spin quantum number. It
C is just a printing code.
C ------------------------------------------------------------
C Other variables
      CHARACTER*3 NUM1,NUM2
      CHARACTER*16 FNM
      DIMENSION IOFFA(8), IOFFO(8)
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      END DO
C IOFFO=NR OF OCC ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFO(1)=0
      DO J=1,NSYM-1
        IOFFO(J+1)=IOFFO(J)+NOSH(J)
      END DO
C Subroutine starts
      LU=51
      LU=IsFreeUnit(LU)
      WRITE(NUM1,'(I3.3)') ISTATE
      WRITE(NUM2,'(I3.3)') JSTATE
C AUGSPIN
      IF(AUGSPIN.EQ.1) THEN
       FNM='r2TM_K2V_'//NUM1//'_'//NUM2
      ELSE IF(AUGSPIN.EQ.-1) THEN
       FNM='r2TM_SDA_'//NUM1//'_'//NUM2
C if AAB (all spin) true
C      ELSE IF(aab.and.AUGSPIN.EQ.2) THEN
C       FNM='r2TM_BBB_'//NUM1//'_'//NUM2
C      ELSE IF(aab.and.AUGSPIN.EQ.3) THEN
C       FNM='r2TM_AAA_'//NUM1//'_'//NUM2
C      ELSE IF(aab.and.AUGSPIN.EQ.4) THEN
C       FNM='r2TM_AAB_'//NUM1//'_'//NUM2
C      ELSE IF(aab.and.AUGSPIN.EQ.5) THEN
C       FNM='r2TM_BBA_'//NUM1//'_'//NUM2
C      ELSE IF(aab.and.AUGSPIN.EQ.6) THEN
C       FNM='r2TM_ABA_'//NUM1//'_'//NUM2
C      ELSE IF(aab.and.AUGSPIN.EQ.7) THEN
C       FNM='r2TM_BAB_'//NUM1//'_'//NUM2
      END IF
      CALL Molcas_Open(LU,FNM)
      WRITE(LU,*)'# Auger Densities: CMO1, CMO2, 1-e Dyson, 2-e Dyson'
      IF(AUGSPIN.EQ.1) THEN
      WRITE(LU,*)'# Spin Matrix (K-2V) for RAES and NAES'
      ELSE IF(AUGSPIN.EQ.-1) THEN
      WRITE(LU,*)'# Spin matrix SDA for NAES.'
C     if AAB (all spin) true
C      ELSE IF(aab.and.AUGSPIN.EQ.2) THEN
C      WRITE(LU,*)'# Spin BBB 2-el rTDM density matrix.'
C      ELSE IF(aab.and.AUGSPIN.EQ.3) THEN
C      WRITE(LU,*)'# Spin AAA 2-el rTDM density matrix.'
C      ELSE IF(aab.and.AUGSPIN.EQ.4) THEN
C      WRITE(LU,*)'# Spin AAB 2-el rTDM density matrix.'
C      ELSE IF(aab.and.AUGSPIN.EQ.5) THEN
C      WRITE(LU,*)'# Spin BBA 2-el rTDM density matrix.'
C      ELSE IF(aab.and.AUGSPIN.EQ.6) THEN
C      WRITE(LU,*)'# Spin ABA 2-el rTDM density matrix.'
C      ELSE IF(aab.and.AUGSPIN.EQ.7) THEN
C      WRITE(LU,*)'# Spin BAB 2-el rTDM density matrix.'
      END IF

      WRITE(LU,*)'# OCA for scattering atom:'
      WRITE(LU,*) OCAN
      DO I=1,OCAN
      WRITE(LU,*) OCAA(I)
      END DO
      WRITE(LU,*)'# Binding energy (eV)'
      WRITE(LU,'(5ES19.12)') EIJ

      SYM12=MUL(LSYM1,LSYM2)
      WRITE(LU,*)'# Total Symmetry of the WF product (<N-1,N>):'
      WRITE(LU,*) SYM12
      WRITE(LU,*)'# States:',ISTATE, JSTATE
      WRITE(LU,*)'# Nr of irreps:',NSYM
      WRITE(LU,*) NSYM
      WRITE(LU,'(A17,8I5)')' # Basis func   :',(NBASF(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(8I5)') (NBASF(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(A17,8I5)')' # Frozen  orb   :',(NFRO(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(A17,8I5)')' # Inactive orb  :',(NISH(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(A17,8I5)')' # Active orb    :',(NASH(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(8I5)') (NASH(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(A17,8I5)')' # Total num orb :',(NOSH(ISYM),ISYM=1,NSYM)
      WRITE(LU,'(8I5)') (NOSH(ISYM),ISYM=1,NSYM)
      WRITE(LU,*)'# CMO1 Molecular orbitals'
      LPOS=1
      DO ISYM=1,NSYM
        NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        NB=NBASF(ISYM)
        DO IO=1,NO
          WRITE(LU,*)'# Symm ',ISYM,'   Orbital ',IO
          DO i=0,NB-1
          WRITE(LU,'(5ES19.12)') CMO1(LPOS+NB*(IO-1)+i)
          END DO
        END DO
        LPOS=LPOS+NB*NO
      END DO
      WRITE(LU,*)'# CMO2 Molecular orbitals'
      LPOS=1
      DO ISYM=1,NSYM
        NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        NB=NBASF(ISYM)
        DO IO=1,NO
          WRITE(LU,*)'# Symm ',ISYM,'   Orbital ',IO
          DO i=0,NB-1
          WRITE(LU,'(5ES19.12)') CMO2(LPOS+NB*(IO-1)+i)
          END DO
        END DO
        LPOS=LPOS+NB*NO
      END DO
C Write Dyson orbitals in CI basis
      WRITE(LU,*)'# 1-e Dyson orbital for CI coeff. in MO biorth. basis'
      WRITE(LU,*)'# Symmetry Block elements:',NDYSAB
      WRITE(LU,*)'# sub-Block info:Sym(I), NumOrb in SymmBlock'
      IOFFTD=0
      DO ISYI=1,NSYM
       NOI=NOSH(ISYI)
       NII=NISH(ISYI)
       IF(NOI.EQ.0) GOTO 400
       WRITE(LU,'(A10,8I7)')' # sub-Block:',ISYI,NOI
       DO I=1,NOI
         IA=I+IOFFTD
C        eliminate small numbers
         IF(ABS(DYSAB(IA)).LT.1.0D-29) THEN
            DYSAB(IA)=0.0D0
         END IF
         write(LU,'(I7,ES22.12)') IA,DYSAB(IA)
       END DO
400    CONTINUE
       NORBSYM=NOI
       IOFFTD=IOFFTD+NORBSYM
      END DO
C Write reduced 2-e TDM in CI basis.
      IOFFTD=0
      WRITE(LU,*)'# 2-e reduced TDM for CI coeff. in MO biorth. basis'
      WRITE(LU,*)'# Symmetry Block elements:',NRT2MAB
      WRITE(LU,*) NRT2MAB
      WRITE(LU,*)'# sub-Block info:Sym(I,J,L), NumOrb in SymmBlock'
      DO ISYI=1,NSYM
       NOI=NOSH(ISYI)
       NAI=NASH(ISYI)
       NII=NISH(ISYI)
       IF(NOI.EQ.0) GOTO 270
       DO ISYJ=1,NSYM
        NOJ=NOSH(ISYJ)
        NAJ=NASH(ISYJ)
        NIJ=NISH(ISYJ)
        IF(NOJ.EQ.0) GOTO 370
        DO ISYL=1,NSYM
         NOL=NOSH(ISYL)
         NAL=NASH(ISYL)
         NIL=NISH(ISYL)
         IF(NOL.EQ.0) GOTO 470
          IF(MUL(ISYI,MUL(ISYJ,ISYL)).EQ.SYM12) THEN
            IF(NAI.EQ.0) GOTO 670
            IF(NAJ.EQ.0) GOTO 670
            IF(NAL.EQ.0) GOTO 670
            WRITE(LU,'(A10,8I7,8I7,8I7,8I7)')' # sub-Block:',ISYI,ISYJ,
     &       ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFA(ISYI)+I-NII
             IO=IOFFO(ISYI)+I
             DO J=1,NOJ
              JA=IOFFA(ISYJ)+J-NIJ
              JO=IOFFO(ISYJ)+J
              DO L=1,NOL
               LA=IOFFA(ISYL)+L-NIL
               LO=IOFFO(ISYL)+L
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
               write(LU,'(I7,I7,I7,ES26.12)') IO,JO,LO,0.0D0
               ELSE
                KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
                IF(ABS(RT2M(KPOS)).LT.1.0D-39) THEN
                 RT2M(KPOS) = 0.0D0
                END IF
                write(LU,'(I7,I7,I7,ES26.12)') IO,JO,LO,
     &           RT2M(KPOS)
               END IF
              END DO
             END DO
            END DO
670         CONTINUE
          END IF
470     CONTINUE
        END DO
370    CONTINUE
       END DO
270   CONTINUE
      END DO
      CLOSE (LU)
      END SUBROUTINE
