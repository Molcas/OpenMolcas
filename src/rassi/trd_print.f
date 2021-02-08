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
* Copyright (C) 2011, Per Ake Malmqvist                                *
************************************************************************

! Print the transition density matrices in ASCII format.
!   Code written by P. A. Malmqvist.
! This code was moved from the main gtdmctl.f file for clarity.
! - F. Plasser
      SUBROUTINE TRD_PRINT(ISTATE, JSTATE, DO22, TDMAB, TDM2,
     &                     CMO1, CMO2, SIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='TRD_PRINT')
#include "rasdim.fh"
!#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "stdalloc.fh"
! Variables passed
      INTEGER ISTATE, JSTATE
      REAL*8 TDMAB(*), TDM2(*), CMO1(*), CMO2(*), SIJ
      LOGICAL DO22
! Other variables
      CHARACTER*3 NUM1,NUM2
      CHARACTER*12 FNM
      DIMENSION WBUF(5)
! Subroutine starts
      LU=50
      LU=IsFreeUnit(LU)
      WRITE(NUM1,'(I3.3)') ISTATE
      WRITE(NUM2,'(I3.3)') JSTATE
      FNM='TRD2_'//NUM1//'_'//NUM2
      CALL Molcas_Open(LU,FNM)
      WRITE(LU,*)'#Transition density file from RASSI.'
      WRITE(LU,*)'#  States:'
      WRITE(LU,*) ISTATE, JSTATE
      WRITE(LU,*)'#  Nr of irreps:'
      WRITE(LU,*) NSYM
      WRITE(LU,*)'#  Basis functions:'
      WRITE(LU,'(8I5)') (NBASF(ISYM),ISYM=1,NSYM)
      WRITE(LU,*)'#  Frozen orbitals:'
      WRITE(LU,'(8I5)') (NFRO(ISYM),ISYM=1,NSYM)
      WRITE(LU,*)'#  Inactive orbitals:'
      WRITE(LU,'(8I5)') (NISH(ISYM),ISYM=1,NSYM)
      WRITE(LU,*)'#  Active orbitals:'
      WRITE(LU,'(8I5)') (NASH(ISYM),ISYM=1,NSYM)
      WRITE(LU,*)'#  State ',ISTATE,'    CMO coefficients:'
      LPOS=1
      DO ISYM=1,NSYM
        NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        NB=NBASF(ISYM)
        DO IO=1,NO
          WRITE(LU,*)'#  Symm ',ISYM,'   Orbital ',IO
          WRITE(LU,'(5D19.12)')(CMO1(LPOS+NB*(IO-1)+i),i=0,NB-1)
        END DO
        LPOS=LPOS+NB*NO
      END DO
      WRITE(LU,*)'#  State ',JSTATE,'    CMO coefficients:'
      LPOS=1
      DO ISYM=1,NSYM
        NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        NB=NBASF(ISYM)
        DO IO=1,NO
          WRITE(LU,*)'#  Symm ',ISYM,'   Orbital ',IO
          WRITE(LU,'(5D19.12)')(CMO2(LPOS+NB*(IO-1)+i),i=0,NB-1)
        END DO
        LPOS=LPOS+NB*NO
      END DO
      WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Overlap:'
      WRITE(LU,'(5D19.12)') SIJ
      WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Active TRD1:'
      LSYM12=MUL(LSYM1,LSYM2)
      LPOS=1
      DO ISYM1=1,NSYM
        NO1=NOSH(ISYM1)
        ISYM2=MUL(ISYM1,LSYM12)
        NO2=NOSH(ISYM2)
        IF (NO1*NO2 .gt. 0) THEN
          NA1=NASH(ISYM1)
          NA2=NASH(ISYM2)
          IF (NA1*NA2 .gt. 0) THEN
            NI1=NISH(ISYM1)
            NI2=NISH(ISYM2)
            WRITE(LU,*)'#  Symmetries ',ISYM1,ISYM2
            WRITE(LU,'(5D19.12)')((TDMAB(LPOS-1+II+NO1*(JJ-1)),
     &                                  JJ=NI2+1,NO2),II=NI1+1,NO1)
          END IF
          LPOS=LPOS+NO1*NO2
        END IF
      END DO

      IF (DO22) THEN
        WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Active TRD2:'
        DO ISYT=1,NSYM
          DO ISYU=1,NSYM
            DO ISYV=1,ISYT
              LIMX=ISYV
              IF(ISYV.EQ.ISYT) LIMX=ISYU
              DO ISYX=1,LIMX
!               > Write out one symmetry block (4 indices!) of two-electron
!               > transition density matrix elements.
!               > Write a full 'rectangular' array, even if it could be made
!               > smaller by permutation symmetry.
                WRITE(LU,*)'#  Orbital symm:',ISYT,ISYU,ISYV,ISYX
                IWBUF=0
                DO IT=1,NASH(ISYT)
                  ITABS=NAES(ISYT)+IT
                  DO IU=1,NASH(ISYU)
                    IUABS=NAES(ISYU)+IU
                    ITU=ITABS+NASHT*(IUABS-1)
                    DO IV=1,NASH(ISYV)
                      IVABS=NAES(ISYV)+IV
                      DO IX=1,NASH(ISYX)
                        IXABS=NAES(ISYX)+IX
                        IVX=IVABS+NASHT*(IXABS-1)
                        IF(ITU.GE.IVX) THEN
                          ITUVX=(ITU*(ITU-1))/2+IVX
                        ELSE
                          ITUVX=(IVX*(IVX-1))/2+ITU
                        END IF
                        IWBUF=IWBUF+1
                        WBUF(IWBUF)=TDM2(ITUVX)
                        IF(IWBUF.EQ.5) THEN
                          WRITE(LU,'(5D19.12)')(WBUF(I),I=1,IWBUF)
                          IWBUF=0
                        END IF
                      END DO
                    END DO
                  END DO
                END DO
                IF(IWBUF.GT.0) THEN
                  WRITE(LU,'(5D19.12)')(WBUF(I),I=1,IWBUF)
                  IWBUF=0
                END IF
* End of writing a symmetry block.
              END DO
            END DO
          END DO
        END DO
      END IF
      CLOSE (LU)
      END SUBROUTINE
