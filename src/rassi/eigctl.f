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
      SUBROUTINE EIGCTL(PROP,OVLP,DYSAMPS,HAM,EIGVEC,ENERGY)
      USE kVectors
#include "compiler_features.h"
#ifndef POINTER_REMAP
      USE ISO_C_Binding
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='EIGCTL')
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
* constants.fh defines values of physical constants.
#include "constants.fh"
#include "stdalloc.fh"
#include "rassiwfn.fh"

      character*100 line
      REAL*8 PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &       HAM(NSTATE,NSTATE),EIGVEC(NSTATE,NSTATE),ENERGY(NSTATE),
     &       DYSAMPS(NSTATE,NSTATE)
      REAL*8, ALLOCATABLE :: ESFS(:)
      REAL*8, Dimension(:,:), Allocatable :: TDS
      Logical Get_TDS, Diagonal
      Integer, Dimension(:), Allocatable :: IndexE,TMOgrp1,TMOgrp2
* Short array, just for putting transition dipole values
* into Add_Info, for generating check numbers:
      Real*8 TDIPARR(3)
      Integer  cho_x_gettol
      External cho_x_gettol
      INTEGER IOFF(8), SECORD(4), IPRTMOM(12)
      CHARACTER*8 LABEL
      Real*8 TM_R(3), TM_I(3), TM_C(3)
      Character*60 FMTLINE
      Real*8 Wavevector(3), UK(3)
      Real*8, Allocatable :: pol_Vector(:,:)
#ifdef _HDF5_
      Real*8, Allocatable, Target :: Storage(:,:,:,:)
      Real*8, Pointer :: flatStorage(:)
#endif

#ifdef _DEBUG_RASSI_
      logical :: debug_dmrg_rassi_code = .true.
#else
      logical :: debug_dmrg_rassi_code = .false.
#endif
      logical TMOgroup
      REAL*8 COMPARE



      CALL QENTER(ROUTINE)
C CONSTANTS:
      AU2EV=CONV_AU_TO_EV_
      AU2CM=CONV_AU_TO_CM1_
      DEBYE=CONV_AU_TO_DEBYE_
      AU2REDR=2.0D2*DEBYE
      SPEED_OF_LIGHT=CONST_C_IN_AU_
      HALF=0.5D0
      PI= CONST_PI_
*
      DIAGONAL=.TRUE.

#ifdef _DEBUG_RASSI_
      write(6,*) 'BLUBB start of eigctl: debug print of property matrix'
        do istate = 1, nstate
        do jstate = 1, nstate
        DO IPROP=1,NPROP
          if(abs(prop(istate,jstate,iprop)) > 1.0d-14)
     &    write(6,*) 'prop(',istate,',',jstate,',',iprop,') = ',
     &                prop(istate,jstate,iprop)
        end do
        end do
        end do
#endif

C DIAGONALIZE SCALAR HAMILTONIAN.

C Initialize eigenvector array.
      DO J=1,NSTATE
        DO I=1,NSTATE
          EIGVEC(I,J)=0.0D0
        END DO
      END DO
C NOTE: It is imperative that we do not mix, or change order, of
C states with different nr of electrons, spin, or symmetry. It is
C assumed in some subsequent parts of the program that these
C remain good properties of the wave functions, with values as
C listed in various tables in common /CNTRL/.
C So it is worth the extra inconvenience to construct an outer
C loop over sets of interacting wave functions.
C Make a list of interacting sets of states:
      CALL GETMEM('LIST','ALLO','INTE',LLIST,NSTATE)
      DO I=1,NSTATE
        IWORK(LLIST-1+I)=0
      END DO
      ISET=0
      DO I=1,NSTATE
        IF(IWORK(LLIST-1+I).GT.0) GOTO 20
        ISET=ISET+1
        IWORK(LLIST-1+I)=ISET
        JOB1=iWork(lJBNUM+I-1)
        NACTE1=NACTE(JOB1)
        MPLET1=MLTPLT(JOB1)
        LSYM1=IRREP(JOB1)
        DO J=I+1,NSTATE
          IF(IWORK(LLIST-1+J).GT.0) GOTO 10
          JOB2=iWork(lJBNUM+J-1)
          NACTE2=NACTE(JOB2)
          IF(NACTE2.NE.NACTE1) GOTO 10
          MPLET2=MLTPLT(JOB2)
          IF(MPLET2.NE.MPLET1) GOTO 10
          LSYM2=IRREP(JOB2)
          IF(LSYM2.NE.LSYM1) GOTO 10
          IWORK(LLIST-1+J)=ISET
  10      CONTINUE
        END DO
  20    CONTINUE
      END DO
      NSETS=ISET
CTEST      write(*,*)' EIGCTL. There are NSETS sets of interacting states.'
CTEST      write(*,'(1x,a,i3)')' where NSETS=',NSETS
CTEST      write(*,*)' The LIST array:'
CTEST      write(*,'(1x,20i3)')(IWORK(LLIST-1+I),I=1,NSTATE)

      NHH=(NSTATE*(NSTATE+1))/2
      CALL GETMEM('HH','ALLO','REAL',LHH,NHH)
      CALL GETMEM('HSQ','ALLO','REAL',LHSQ,NSTATE**2)
      CALL GETMEM('SS','ALLO','REAL',LSS,NHH)
      CALL GETMEM('UU','ALLO','REAL',LUU,NSTATE**2)
      CALL GETMEM('SCR','ALLO','REAL',LSCR,NSTATE**2)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LSCR),1)
      CALL GETMEM('STACK','ALLO','INTE',LSTK,NSTATE)
C Loop over the sets:
      DO ISET=1,NSETS
C Stack up the states belonging to this set:
       MSTATE=0
       DO I=1,NSTATE
        JSET=IWORK(LLIST-1+I)
        IF(JSET.EQ.ISET) THEN
         MSTATE=MSTATE+1
         IWORK(LSTK-1+MSTATE)=I
        END IF
       END DO

       if(debug_dmrg_rassi_code)then
         write(6,*) 'BLUBB DEBUG print of Hamiltonian and overlap'
       end if

C 1. PUT UNIT MATRIX INTO UU
      CALL DCOPY_(MSTATE**2,[0.0D0],0,WORK(LUU),1)
      CALL DCOPY_(MSTATE   ,[1.0D0],0,WORK(LUU),MSTATE+1)
C 2. COPY OVERLAP MATRIX INTO TRIANGULAR STORAGE,
C    and Hamiltonian into square storage:
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LHSQ),1)
      IJ=0
      DO II=1,MSTATE
        I=IWORK(LSTK-1+II)
        DO JJ=1,II
          J=IWORK(LSTK-1+JJ)
          IJ=IJ+1
          if(debug_dmrg_rassi_code)then
            write(6,*) 'overlap     for i,j',i,j,ovlp(i,j)
            write(6,*) 'Hamiltonian for i,j',i,j,HAM(i,j)
          end if
          WORK(LSS-1+IJ)=OVLP(I,J)
          WORK(LHSQ-1+II+MSTATE*(JJ-1))=HAM(I,J)
          WORK(LHSQ-1+JJ+MSTATE*(II-1))=HAM(I,J)
        END DO
      END DO
C 3. SPECTRAL DECOMPOSITION OF OVERLAP MATRIX:
      CALL Jacob(WORK(LSS),WORK(LUU),MSTATE,MSTATE)
      II=0
      DO I=1,MSTATE
        II=II+I
        X=1.0D00/SQRT(MAX(0.5D-14,WORK(LSS-1+II)))
        DO K=1,MSTATE
          LPOS=LUU-1+K+MSTATE*(I-1)
          WORK(LPOS)=X*WORK(LPOS)
        END DO
      END DO
      IJ=0
      DO I=1,MSTATE
        DO J=1,I
          IJ=IJ+1
          If (I.NE.J .AND.
     &    ABS(WORK(LHSQ-1+I+MSTATE*(J-1))).gt.1.0D-9) Then
            DIAGONAL=.FALSE.
            go to 11
          EndIf
        End Do
      End Do
 11   Continue
      If (.not.diagonal) Then
C 4. TRANSFORM HAMILTON MATRIX.
*        CALL MXMA(WORK(LHSQ),1,MSTATE,
*     &            WORK(LUU),1,MSTATE,
*     &            WORK(LSCR),1,MSTATE,
*     &            MSTATE,MSTATE,MSTATE)
        CALL DGEMM_('N','N',MSTATE,MSTATE,MSTATE,1.0D0,
     &             WORK(LHSQ),MSTATE,WORK(LUU),MSTATE,
     &             0.0D0,WORK(LSCR),MSTATE)
*        CALL MXMA(WORK(LUU),MSTATE,1,
*     &            WORK(LSCR),1,MSTATE,
*     &            WORK(LHSQ),1,MSTATE,
*     &            MSTATE,MSTATE,MSTATE)
        CALL DGEMM_('T','N',MSTATE,MSTATE,MSTATE,1.0D0,
     &             WORK(LUU),MSTATE,WORK(LSCR),MSTATE,
     &             0.0D0,WORK(LHSQ),MSTATE)

C 5. DIAGONALIZE HAMILTONIAN.
      IJ=0
      DO I=1,MSTATE
        DO J=1,I
          IJ=IJ+1
          WORK(LHH-1+IJ)=WORK(LHSQ-1+I+MSTATE*(J-1))
        END DO
      END DO

      CALL Jacob (WORK(LHH),WORK(LUU),MSTATE,MSTATE)
      CALL Jacord(WORK(LHH),WORK(LUU),MSTATE,MSTATE)

      IDIAG=0
      DO II=1,MSTATE
        IDIAG=IDIAG+II
        I=IWORK(LSTK-1+II)
        ENERGY(I)=WORK(LHH-1+IDIAG)
        DO JJ=1,MSTATE
          J=IWORK(LSTK-1+JJ)
          EIGVEC(I,J)=WORK(LUU-1+II+MSTATE*(JJ-1))
        END DO
      END DO
*if diagonal
      Else
        DO II=1,MSTATE
          I=IWORK(LSTK-1+II)
          ENERGY(I)=HAM(I,I)
          Do JJ=1,MSTATE
            J=IWORK(LSTK-1+JJ)
            EIGVEC(I,J)=0.0d0
          End Do
          EIGVEC(I,I)=1.0d0
        End Do
      EndIf


CUNGUR
c   Correct for diagonal energies in case of orbital degeneracy:
c   Convention: two energies are considered degenerate if their energy difference is
c               lower than 1.0D-4 cm-1
      TMP=0.d0
      DLT=0.d0
      IDIAG=0
      DO II=1,MSTATE
        I=IWORK(LSTK-1+II)
        TMP=ENERGY(I)
        JDIAG=0
        Do JJ=1,MSTATE
          J=IWORK(LSTK-1+JJ)
          IF(I==J) CYCLE
          DLT=ABS(ENERGY(J)-TMP)*AU2CM
          If(DLT<1.0D-4) THEN
            ENERGY(J)=TMP
          End If
        End Do
      End Do

      IDIAG=0
      DO II=1,MSTATE
        IDIAG=IDIAG+II
        I=IWORK(LSTK-1+II)
        WORK(LHH-1+IDIAG)=ENERGY(I)
      END DO
C End of loop over sets.
      END DO
C Morgane Vacher 02/17 - Fix the "arbitrary" sign of
C the eigenvectors such that the largest coefficient
C is positive. This is to avoid spurious changes of
C sign of the SFS with respect to the original ones,
C especially for already diagonal Hamiltonian matrix.
      do i=1,nstate
         j = maxloc(abs(eigvec(:,i)),1)
         if (eigvec(i,j) .lt. 0.0d0) then
           eigvec(:,i) = -eigvec(:,i)
         endif
      enddo
      CALL GETMEM('HH','FREE','REAL',LHH,NHH)
      CALL GETMEM('SS','FREE','REAL',LSS,NHH)
      CALL GETMEM('UU','FREE','REAL',LUU,NSTATE**2)
      CALL GETMEM('HSQ','FREE','REAL',LHSQ,NSTATE**2)
      CALL GETMEM('STACK','FREE','INTE',LSTK,NSTATE)
      CALL GETMEM('LIST','FREE','INTE',LLIST,NSTATE)

      IF(IPGLOB.GE.TERSE) THEN
       DO ISTATE=1,NSTATE
        Call PrintResult(6,'(6x,A,I5,5X,A,F16.8)',
     &    'RASSI State',ISTATE,'Total energy:',ENERGY(ISTATE),1)
       END DO
      END IF
#ifdef _HDF5_
      call mh5_put_dset(wfn_sfs_energy, ENERGY)
#endif
C To handle extreme cases of large energies/small energy differences
C all TOTAL energies will undergo a universal constant shift:
      EMIN=1.0D12
      DO ISTATE=1,NSTATE
cvv NAG compiler overoptimize this!
c       EMIN=MIN(EMIN,ENERGY(ISTATE))
       if(ENERGY(ISTATE).lt.EMIN) EMIN=ENERGY(ISTATE)
      END DO
      KAU= INT(EMIN/1000.0D0)
      EVAC=1000.0D0*DBLE(KAU)
c      EMIN=EVAC
      DO ISTATE=1,NSTATE
c        ENERGY(ISTATE)=ENERGY(ISTATE)-EVAC
        ENERGY(ISTATE)=ENERGY(ISTATE)-EMIN
      END DO

C Put energies onto info file for automatic verification runs:
CPAM06 Added error estimate, based on independent errors for all
C components of H and S in original RASSCF wave function basis:
      EPSH=MAX(5.0D-10,ABS(ENERGY(1))*5.0D-11)
      EPSS=5.0D-11
      IDX=100
      DO I=1,NSTATE
       EI=ENERGY(I)
       V2SUM=0.0D0
       DO J=1,NSTATE
        V2SUM=V2SUM+EIGVEC(J,I)**2
       END DO
       ERMS=SQRT(EPSH**2+EI**2*EPSS**2)*V2SUM
       IDX=MIN(IDX,INT(-LOG10(ERMS)))
      END DO
      iTol=cho_x_gettol(IDX) ! reset thr iff Cholesky
      Call Add_Info('E_RASSI',ENERGY+EMIN-EVAC,NSTATE,iTol)

C Experimental addition: Effective L and/or M quantum numbers.

C Identify which properties are angular moment matrix elements:
      IAMX=0
      IAMY=0
      IAMZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO
*                                                                      *
************************************************************************
*                                                                      *
C The matrix elements are actually for i*Lx, etc.
C Now form matrix elements of L**2 assuming closure:
C i.e. assume L**2 = -(iLx)*(iLx)-(iLy)*(iLy)-(iLz)*(iLz)
C within the basis formed by the states.
      IAMXYZ=0
      IF (IAMZ.GT.0) THEN
         CALL GETMEM('L2','ALLO','REAL',LL2,NSTATE**2)
         CALL GETMEM('M2DIA','ALLO','REAL',LM2DIA,NSTATE)
         CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-1.0D0,
     &             PROP(1,1,IAMZ),NSTATE,PROP(1,1,IAMZ),NSTATE,
     &             0.0D0,WORK(LL2),NSTATE)
         CALL DCOPY_(NSTATE,WORK(LL2),(NSTATE+1),WORK(LM2DIA),1)
         IF (IAMX.GT.0.and.IAMY.gt.0) THEN
            IAMXYZ=1
            CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-1.0D0,
     &               PROP(1,1,IAMX),NSTATE,PROP(1,1,IAMX),NSTATE,
     &               1.0D0,WORK(LL2),NSTATE)
            CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,-1.0D0,
     &               PROP(1,1,IAMY),NSTATE,PROP(1,1,IAMY),NSTATE,
     &               1.0D0,WORK(LL2),NSTATE)
            CALL GETMEM('L2DIA','ALLO','REAL',LL2DIA,NSTATE)
            CALL DCOPY_(NSTATE,WORK(LL2),(NSTATE+1),WORK(LL2DIA),1)
         END IF
         CALL GETMEM('L2','FREE','REAL',LL2,NSTATE**2)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
*      Sort the states energywise
*
       Call mma_Allocate(IndexE,nState,Label='IndexE')
       Do iState = 1, nState
          IndexE(iState)=iState
       End Do
       Do iState = 1, nState-1
          EX=ENERGY(IndexE(iState))
*
          kState=iState
          Do jState = iState+1, nState
             If (ENERGY(IndexE(jState)).lt.EX) Then
                kState=jState
                EX=ENERGY(IndexE(jState))
             End If
          End Do
          If (kState.ne.iState) Then
             lState=IndexE(iState)
             IndexE(iState)=IndexE(kState)
             IndexE(kState)=lState
          End If
       End Do
*                                                                      *
************************************************************************
*                                                                      *
C REPORT ON SECULAR EQUATION RESULT:
      CALL MMA_ALLOCATE(ESFS,NSTATE)
      IF(IPGLOB.ge.TERSE) THEN
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A,34X,A)')
     &     '*','       Spin-free section      ','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,17X,A,17X,A)')
     &     '*','Note: index according to input order, order according'
     &     //' to energy.','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)' SPIN-FREE ENERGIES:'
       WRITE(6,'(1X,A,F22.10,A1)')' (Shifted by EVAC (a.u.) =',EMIN,')'
       WRITE(6,*)
       IF(IFJ2.ne.0 .and. IAMXYZ.gt.0) THEN
        IF(IFJZ.ne.0 .and. IAMZ.gt.0) THEN
        WRITE(6,*)'SF State       Relative EVAC(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)      L_eff   Abs_M'
        ELSE
        WRITE(6,*)'SF State       Relative EVAC(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)      L_eff'
        END IF
       ELSE
        IF(IFJZ.ne.0 .and. IAMZ.gt.0) THEN
        WRITE(6,*)'SF State       Relative EVAC(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)      Abs_M'
        ELSE
        WRITE(6,*)'SF State       Relative EVAC(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)'
        END IF
       END IF
       WRITE(6,*)
*
       E0=ENERGY(IndexE(1))
       Do kSTATE=1,NSTATE
          iState=IndexE(kState)
          E1=ENERGY(ISTATE)
          E2=AU2EV*(E1-E0)
          E3=AU2CM*(E1-E0)
*
         IF(IFJ2.ne.0 .and. IAMXYZ.gt.0) THEN
          IF(IFJZ.ne.0 .and. IAMZ.gt.0) THEN
           EFFL=SQRT(MAX(0.5D-12,0.25D0+WORK(LL2DIA-1+ISTATE)))-Half
           EFFM=SQRT(MAX(0.5D-12,WORK(LM2DIA-1+ISTATE)))
           FMTLINE='(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1,2X,'//
     &             'F6.1)'
           WRITE(6,FMTLINE) ISTATE,E1,E2,E3,EFFL,EFFM
          ELSE
           EFFL=SQRT(MAX(0.5D-12,0.25D0+WORK(LL2DIA-1+ISTATE)))-Half
           FMTLINE='(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1)'
           WRITE(6,FMTLINE) ISTATE,E1,E2,E3,EFFL
          END IF
         ELSE
          IF(IFJZ.ne.0 .and. IAMZ.gt.0) THEN
           EFFM=SQRT(MAX(0.5D-12,WORK(LM2DIA-1+ISTATE)))
           FMTLINE='(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4,6X,F6.1)'
           WRITE(6,FMTLINE) ISTATE,E1,E2,E3,EFFM
          ELSE
           FMTLINE='(1X,I5,7X,F18.10,2X,F18.10,2X,F18.4)'
           WRITE(6,FMTLINE) ISTATE,E1,E2,E3
          END IF
         END IF
       ESFS(ISTATE)=E3
*
       End Do
*
       IF(IAMZ.GT.0) THEN
        CALL GETMEM('M2DIA','FREE','REAL',LM2DIA,NSTATE)
       END IF
       IF(IAMXYZ.GT.0) THEN
        CALL GETMEM('L2DIA','FREE','REAL',LL2DIA,NSTATE)
       END IF
      END IF
c LU: save esfs array
       CALL Put_dArray( 'ESFS_SINGLE',ESFS,NSTATE)
       CALL MMA_DEALLOCATE(ESFS)
c

      IF(IPGLOB.ge.VERBOSE) THEN
       WRITE(6,*)
       WRITE(6,*)'  Spin-free eigenstates in basis of input states:'
       WRITE(6,*)'  -----------------------------------------------'
       WRITE(6,*)
       DO L=1,NSTATE
          I=IndexE(L)
         Write(6,'(5X,A,I5,A,F18.10)')'Eigenstate No.',I,
     &         ' energy=',ENERGY(I)
         WRITE(6,'(5X,5F15.7)')(EIGVEC(K,I),K=1,NSTATE)
       END DO
       CALL GETMEM('ILST','ALLO','INTE',LILST,NSTATE)
       CALL GETMEM('VLST','ALLO','REAL',LVLST,NSTATE)
       DO L=1,NSTATE
          I=IndexE(L)
          Write(6,'(5X,A,I5,A,F18.10)')'Eigenstate No.',I,
     &          ' energy=',ENERGY(I)
        EVMAX=0.0D0
        DO K=1,NSTATE
         EVMAX=MAX(EVMAX,ABS(EIGVEC(IndexE(K),I)))
        END DO
        EVLIM=0.10D0*EVMAX
        NLST=0
        DO K=1,NSTATE
         EV=EIGVEC(IndexE(K),I)
         IF(ABS(EV).GE.EVLIM) THEN
           NLST=NLST+1
           WORK(LVLST-1+NLST)=EV
           IWORK(LILST-1+NLST)=IndexE(K)
         END IF
        END DO
         DO KSTA=1,NLST,6
          KEND=MIN(NLST,KSTA+4)
          WRITE(Line,'(5X,5(I5,F12.6))')
     &     (IWORK(LILST-1+K),WORK(LVLST-1+K),K=KSTA,KEND)
          CALL NORMAL(Line)
          WRITE(6,*) Line
         END DO
         WRITE(6,*)
        END DO
        CALL GETMEM('ILST','FREE','INTE',LILST,NSTATE)
        CALL GETMEM('VLST','FREE','REAL',LVLST,NSTATE)
        WRITE(6,*)
        WRITE(6,*)' THE INPUT RASSCF STATES REEXPRESSED IN EIGENSTATES:'
        WRITE(6,*)
        DO L=1,NSTATE
           I=IndexE(L)
         CALL DGEMM_('T','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &             EIGVEC,NSTATE,OVLP,NSTATE,
     &             0.0D0,WORK(LSCR),NSTATE)
         WRITE(6,'(A,I5)')' INPUT STATE NR.:',I
         WRITE(6,*)' OVERLAP WITH THE EIGENSTATES:'
         WRITE(6,'(5(1X,F15.7))')(WORK(LSCR-1+IndexE(K)+NSTATE*(I-1)),
     &         K=1,NSTATE)
         WRITE(6,*)
       END DO
      END IF

C TRANSFORM AND PRINT OUT PROPERTY MATRICES:
      DO IP=1,NPROP
        CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &             PROP(1,1,IP),NSTATE,EIGVEC,NSTATE,
     &             0.0D0,WORK(LSCR),NSTATE)
        CALL DGEMM_('T','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &             EIGVEC,NSTATE,WORK(LSCR),NSTATE,
     &             0.0D0,PROP(1,1,IP),NSTATE)
      END DO

! +++ J. Norell 12/7 - 2018
C And the same for the Dyson amplitudes
        CALL DGEMM_('N','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &             DYSAMPS,NSTATE,EIGVEC,NSTATE,
     &             0.0D0,WORK(LSCR),NSTATE)
        CALL DGEMM_('T','N',NSTATE,NSTATE,NSTATE,1.0D0,
     &             EIGVEC,NSTATE,WORK(LSCR),NSTATE,
     &             0.0D0,DYSAMPS,NSTATE)
! +++ J. Norell

      CALL GETMEM('SCR','FREE','REAL',LSCR,NSTATE**2)
*
* Initial setup for both dipole, quadrupole etc. and exact operator
*
!
! There are debug statements thoughout - look for ZVAL
! If you want to debug in length gauge then comment out velocity dipole
!
!     ZVAL(1) = 1.0D0 ! Simulation of moving the origin along Z
!     ZVAL(2) = 2.0D0
!     ZVAL(3) = 3.0D0
!     ZVAL(4) = 5.0D0
!     ZVAL(5) = 7.0D0
!     ZVAL(6) = 10.0D0
!     ZVAL(7) = 15.0D0
!     ZVAL(8) = 20.0D0
!     ZVAL(9) = 25.0D0

      OSTHR=1.0D-5
      IF(DIPR) OSTHR = OSTHR_DIPR
      IF(DIPR) WRITE(6,*) ' Dipole threshold changed to ',OSTHR
! this is to ensure that the total transistion strength is non-zero
! Negative transitions strengths can occur for quadrupole transistions
! due to the truncation of the Taylor expansion.
      IF(QIPR) OSTHR = OSTHR_QIPR
      IF(QIPR) WRITE(6,*) ' Dipole threshold changed to ',OSTHR,
     &                       ' since quadrupole threshold is given '
      OSTHR2=1.0D-5
      IF(QIPR) OSTHR2 = OSTHR_QIPR
      IF(QIPR) WRITE(6,*) ' Quadrupole threshold changed to ',OSTHR2
      IF(QIALL) WRITE(6,*) ' Will write all quadrupole contributions '
!
!     Reducing the loop over states - good for X-rays
!     At the moment memory is not reduced
!
      IF(REDUCELOOP) THEN
        IEND = LOOPDIVIDE
        JSTART = LOOPDIVIDE+1
      ELSE
        IEND = NSTATE
        JSTART = 1
      END IF
*
      ! AFACTOR = 2*pi*e^2*E_h^2 / eps_0*m_e*c^3*h^2
      ! 1/c^3 (in a.u. of time ^ -1)
      AFACTOR = 2.0D0/CONST_C_IN_AU_**3
     &          /CONST_AU_TIME_IN_SI_
*
      IF(IPGLOB.le.SILENT) GOTO 900
!
* CALCULATION OF THE DIPOLE TRANSITION STRENGTHS
!
!     Initialize arrays for indentifying problematic transitions
!     These stores all dipole oscillator strengths in
!     length and velocity gauge for a later comparison.
!
      CALL GETMEM('DL   ','ALLO','REAL',LDL,NSTATE**2)
      CALL GETMEM('DV   ','ALLO','REAL',LDV,NSTATE**2)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LDL),1)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LDV),1)
      I_HAVE_DL = 0
      I_HAVE_DV = 0
!
      IF(IPGLOB.ge.TERSE) THEN

       IPRDX=0
       IPRDY=0
       IPRDZ=0
       IFANYD=0
       DO IPROP=1,NPROP
        IF(IPUSED(IPROP).NE.0) THEN
         IF(PNAME(IPROP).EQ.'MLTPL  1') THEN
          IFANYD=1
          IF(ICOMP(IPROP).EQ.1) IPRDX=IPROP
          IF(ICOMP(IPROP).EQ.2) IPRDY=IPROP
          IF(ICOMP(IPROP).EQ.3) IPRDZ=IPROP
         END IF
        END IF
       END DO

       IF(IFANYD.NE.0) THEN
        WRITE(6,*)
        CALL CollapseOutput(1,'Dipole transition strengths '//
     &                        '(spin-free states):')
        WRITE(6,'(3X,A)')     '----------------------------'//
     &                        '-------------------'
        IF(OSTHR.GT.0.0D0) THEN
         WRITE(6,30) 'for osc. strength at least',OSTHR
        END IF
        WRITE(6,*)

! Printout the osc. strength in 3 dimensions into a file
       ! Should be if something happen
        losc_strength=20
        losc_strength=isFreeUnit(losc_strength)
        Call Molcas_Open(losc_strength,'osc_strength.au')

        If (Do_SK) Then
           nVec = nk_Vector
        Else
           nVec = 1
        End If
*
        Do iVec = 1, nVec
*
        LNCNT=0
        FMAX=0.0D0
        Two3rds=2.0D0/3.0D0
        DO K_=1,IEND
           I=IndexE(K_)
         DO L_=JSTART,NSTATE
            J=IndexE(L_)
          IJ=I+NSTATE*(J-1)
          EDIFF=ENERGY(J)-ENERGY(I)
*
          IF(EDIFF.GT.0.0D0) THEN
           DX2=0.0D0
           DY2=0.0D0
           DZ2=0.0D0
           If (Do_SK) Then
              tmp=0.0D0
              IF(IPRDX.GT.0) tmp=tmp+PROP(J,I,IPRDX)*k_vector(1,iVec)
              IF(IPRDY.GT.0) tmp=tmp+PROP(J,I,IPRDY)*k_vector(2,iVec)
              IF(IPRDZ.GT.0) tmp=tmp+PROP(J,I,IPRDZ)*k_vector(2,iVec)
              IF(IPRDX.GT.0)
     &           DX2=(PROP(J,I,IPRDX)-tmp*k_vector(1,iVec))**2
              IF(IPRDY.GT.0)
     &           DY2=(PROP(J,I,IPRDY)-tmp*k_vector(2,iVec))**2
              IF(IPRDZ.GT.0)
     &           DZ2=(PROP(J,I,IPRDZ)-tmp*k_vector(3,iVec))**2
           Else
              IF(IPRDX.GT.0) DX2=PROP(J,I,IPRDX)**2
              IF(IPRDY.GT.0) DY2=PROP(J,I,IPRDY)**2
              IF(IPRDZ.GT.0) DZ2=PROP(J,I,IPRDZ)**2
           End If
           FX=Two3rds*EDIFF*(DX2)
           FY=Two3rds*EDIFF*(DY2)
           FZ=Two3rds*EDIFF*(DZ2)
           F =FX+FY+FZ
           FMAX=MAX(F,FMAX)
           AX=(AFACTOR*EDIFF**2)*FX
           AY=(AFACTOR*EDIFF**2)*FY
           AZ=(AFACTOR*EDIFF**2)*FZ
           A =(AFACTOR*EDIFF**2)*F
           IF (F.ge.OSTHR) THEN
              IF (LNCNT.EQ.0) THEN
                 If (Do_SK) Then
                    WRITE(6,*)
                    WRITE(6,'(4x,a,3F8.4)')
     &                 'Direction of the k-vector: ',
     &                  (k_vector(k,iVec),k=1,3)
                    WRITE(6,'(4x,a)')
     &                 'The light is assumed to be unpolarized.'
                    WRITE(6,*)
                 End If
                 WRITE(6,31) 'From','To','Osc. strength',
     &                   'Einstein coefficients Ax, Ay, Az (sec-1)   ',
     &                   'Total A (sec-1)'
                 WRITE(6,32)
                 WRITE(losc_strength,34) 'From','To','Osc. strength',
     &                   'Fx','Fy','Fz','(a.u.)'
                 WRITE(losc_strength,32)
              END IF
              LNCNT=LNCNT+1
              WRITE(6,33) I,J,F,AX,AY,AZ,A
              write(losc_strength,33) I,J,F,Fx,Fy,Fz
           END IF
! Store dipole oscillator strength
            WORK(LDL-1+IJ) = F
*
           If (F.gt.1.0D0) Then
              k=INT(LOG10(F))+1
              F=F/(10.0D0)**k
           End If
           Call Add_Info('TMS(SF,Len)',[F],1,6)
          END IF
         END DO
        END DO
        IF (LNCNT.EQ.0) THEN
           WRITE(6,*)' ( Max oscillator strength is only ',FMAX,')'
        ELSE
           WRITE(6,32)
           WRITE(losc_strength,32)
        END IF
*
        End Do ! iVec
*
        close(losc_strength)
        CALL CollapseOutput(0,'Dipole transition strengths '//
     &                        '(spin-free states):')
        I_HAVE_DL = 1
       END IF

* Key words for printing transition dipole vectors
* PRDIPVEC TDIPMIN
       IF(PRDIPVEC .AND. (NSTATE.gt.1) .and. (IFANYD.NE.0)) THEN
        WRITE(6,*)
        CALL CollapseOutput(1,'Dipole transition vectors '//
     &                        '(spin-free states):')
        WRITE(6,'(3X,A)')     '--------------------------'//
     &                        '-------------------'
        IF(TDIPMIN.GT.0.0D0) THEN
         WRITE(6,30) 'for vector sizes at least',TDIPMIN
        END IF
        WRITE(6,*)

        If (Do_SK) Then
           nVec = nk_Vector
        Else
           nVec = 1
        End If
*
        Do iVec = 1, nVec
*
        LNCNT=0
        DMAX=0.0D0
        DO K_=1,NSTATE-1
           I=IndexE(K_)
         DO L_=I+1,NSTATE
            J=IndexE(L_)
           DX=0.0D0
           DY=0.0D0
           DZ=0.0D0
           If (Do_SK) Then
              tmp=0.0D0
              IF(IPRDX.GT.0) tmp=tmp+PROP(J,I,IPRDX)*k_vector(1,iVec)
              IF(IPRDY.GT.0) tmp=tmp+PROP(J,I,IPRDY)*k_vector(2,iVec)
              IF(IPRDZ.GT.0) tmp=tmp+PROP(J,I,IPRDZ)*k_vector(3,iVec)
              IF(IPRDX.GT.0) DX=PROP(J,I,IPRDX)-tmp*k_vector(1,iVec)
              IF(IPRDY.GT.0) DY=PROP(J,I,IPRDY)-tmp*k_vector(2,iVec)
              IF(IPRDZ.GT.0) DZ=PROP(J,I,IPRDZ)-tmp*k_vector(3,iVec)
           Else
              IF(IPRDX.GT.0) DX2=PROP(J,I,IPRDX)**2
              IF(IPRDY.GT.0) DY2=PROP(J,I,IPRDY)**2
              IF(IPRDZ.GT.0) DZ2=PROP(J,I,IPRDZ)**2
           End If
           IF(IPRDX.GT.0) DX=PROP(J,I,IPRDX)
           IF(IPRDY.GT.0) DY=PROP(J,I,IPRDY)
           IF(IPRDZ.GT.0) DZ=PROP(J,I,IPRDZ)
           DSZ = SQRT(DX**2+DY**2+DZ**2)
           DMAX=MAX(DSZ,DMAX)
           IF(DSZ.ge.TDIPMIN) THEN
            IF(LNCNT.EQ.0) THEN
             WRITE(6,34) 'From','To','Dx','Dy','Dz','Total D (a.u.)'
             WRITE(6,42)
            END IF
            LNCNT=LNCNT+1
            WRITE(6,33) I,J,DX,DY,DZ,DSZ
* Put values into array for add_info:
            TDIPARR(1)=DX
            TDIPARR(2)=DY
            TDIPARR(3)=DZ
            Call Add_Info('TRDIP',TDIPARR,3,6)
           END IF
         END DO
        END DO

        IF(LNCNT.EQ.0) THEN
         WRITE(6,*)' ( Max transition dipole is only ',DMAX,')'
        ELSE
         WRITE(6,32)
        END IF
        CALL CollapseOutput(0,'Dipole transition vectors '//
     &                        '(spin-free states):')
*
*
      End Do ! iVec
*
      End If
      END IF
*
*     Transition moments computed with the velocity operator.
*
      IPRDX=0
      IPRDY=0
      IPRDZ=0
      IFANYD=0
      DO IPROP=1,NPROP
         IF (IPUSED(IPROP).NE.0) THEN
            IF (PNAME(IPROP).EQ.'VELOCITY') THEN
               IFANYD=1
               IF(ICOMP(IPROP).EQ.1) IPRDX=IPROP
               IF(ICOMP(IPROP).EQ.2) IPRDY=IPROP
               IF(ICOMP(IPROP).EQ.3) IPRDZ=IPROP
            END IF
         END IF
      END DO

      IF (IFANYD.NE.0) THEN
         WRITE(6,*)
         CALL CollapseOutput(1,'Velocity transition strengths '//
     &                         '(spin-free states):')
         WRITE(6,'(3X,A)')     '------------------------------'//
     &                         '-------------------'
         IF (OSTHR.GT.0.0D0) THEN
            WRITE(6,30) 'for osc. strength at least',OSTHR
         END IF
         WRITE(6,*)
*
         If (Do_SK) Then
            nVec = nk_Vector
         Else
            nVec = 1
         End If
*
         Do iVec = 1, nVec
*
         LNCNT=0
         FMAX=0.0D0
         Two3rds=2.0D0/3.0D0
         DO K_=1,IEND
            I=IndexE(K_)
            DO L_=JSTART,NSTATE
               J=IndexE(L_)
               IJ=I+NSTATE*(J-1)
               EDIFF=ENERGY(J)-ENERGY(I)
               IF(EDIFF.GT.0.0D0) THEN
               DX2=0.0D0
               DY2=0.0D0
               DZ2=0.0D0
               If (Do_SK) Then
                 tmp=0.0D0
                 IF(IPRDX.GT.0) tmp=tmp+PROP(J,I,IPRDX)*k_vector(1,iVec)
                 IF(IPRDY.GT.0) tmp=tmp+PROP(J,I,IPRDY)*k_vector(2,iVec)
                 IF(IPRDZ.GT.0) tmp=tmp+PROP(J,I,IPRDZ)*k_vector(3,iVec)
                 IF(IPRDX.GT.0)
     &              DX2=(PROP(J,I,IPRDX)-tmp*k_vector(1,iVec))**2
                 IF(IPRDY.GT.0)
     &              DY2=(PROP(J,I,IPRDY)-tmp*k_vector(2,iVec))**2
                 IF(IPRDZ.GT.0)
     &              DZ2=(PROP(J,I,IPRDZ)-tmp*k_vector(3,iVec))**2
               Else
                  IF(IPRDX.GT.0) DX2=PROP(J,I,IPRDX)**2
                  IF(IPRDY.GT.0) DY2=PROP(J,I,IPRDY)**2
                  IF(IPRDZ.GT.0) DZ2=PROP(J,I,IPRDZ)**2
               End If
               FX=Two3rds*(DX2)/EDIFF
               FY=Two3rds*(DY2)/EDIFF
               FZ=Two3rds*(DZ2)/EDIFF
               F =FX+FY+FZ
               FMAX=MAX(F,FMAX)
               AX=(AFACTOR*EDIFF**2)*FX
               AY=(AFACTOR*EDIFF**2)*FY
               AZ=(AFACTOR*EDIFF**2)*FZ
               A =(AFACTOR*EDIFF**2)*F
               IF (F.ge.OSTHR) THEN
                  IF (LNCNT.EQ.0) THEN
                     If (Do_SK) Then
                        WRITE(6,*)
                        WRITE(6,'(4x,a,3F8.4)')
     &                        'Direction of the k-vector: ',
     &                         (k_vector(k,ivec),k=1,3)
                        WRITE(6,'(4x,a)')
     &                        'The light is assumed to be unpolarized.'
                        WRITE(6,*)
                     End If
                     WRITE(6,31) 'From','To','Osc. strength',
     &                   'Einstein coefficients Ax, Ay, Az (sec-1)   ',
     &                   'Total A (sec-1)'
                     WRITE(6,32)
                  END IF
                  LNCNT=1
                  WRITE(6,33) I,J,F,AX,AY,AZ,A
! Store dipole oscillator strength
                  WORK(LDV-1+IJ) = F
               END IF
               Call Add_Info('TMS(SF,Vel)',[F],1,6)
               END IF
            END DO
         END DO
         IF (LNCNT.EQ.0) THEN
            WRITE(6,*)' ( Max oscillator strength is only ',FMAX,')'
         ELSE
            WRITE(6,32)
         END IF
*
         End Do ! iVec
         CALL CollapseOutput(0,'Velocity transition strengths '//
     &                         '(spin-free states):')
         WRITE(6,*)
         I_HAVE_DV = 1
      END IF
!
!      Compare oscillator strengths in length and velocity gauge
!      All differences in oscillator strengths above the tolerance
!      of 0.1 (10 percent) will be printed.
!
       IF(I_HAVE_DL.EQ.1.AND.I_HAVE_DV.EQ.1) THEN
         CALL CollapseOutput(1,'Length and velocity gauge comparison '//
     &                         '(spin-free states):')
         WRITE(6,'(3X,A)')     '-------------------------------------'//
     &                         '-------------------'
!
! I guess that I have to explain it when I print a warning
!
         WRITE(6,*)
         WRITE(6,*) "--------------------------------------------------"
         WRITE(6,*) "A comparison between the dipole oscillator "//
     &              "strengths in "
         WRITE(6,*) "length and velocity gauge "//
     &              "will be performed"
         WRITE(6,*)
         WRITE(6,*) "All dipole oscillator differences above the "//
     &              "tolerance of ",TOLERANCE," will be printed "
         WRITE(6,*)
         WRITE(6,*) "Due to basis set deficiency these oscillator "//
     &              "may be problematic "
         WRITE(6,*)
         WRITE(6,*) "The tolerance is defined as ABS(1-O_l/O_v) "
         WRITE(6,*) "O_l : dipole oscillator strength in length gauge"
         WRITE(6,*) "O_p : dipole oscillator strength in velocity gauge"
         WRITE(6,*) "--------------------------------------------------"
!
          I_PRINT_HEADER = 0
          DO K_=1,IEND
             I=IndexE(K_)
            DO L_=JSTART,NSTATE
               J=IndexE(L_)
               IJ=I+NSTATE*(J-1)
               EDIFF=ENERGY(J)-ENERGY(I)
               IF(JSTART.EQ.1.AND.EDIFF.LT.0.0D0) CYCLE
*
             COMPARE=0.0D0
             dlt=1.0D-18 ! Add small value to avoid zero divide.
             IF(WORK(LDL-1+IJ).GE.OSTHR+dlt .AND.
     &          WORK(LDV-1+IJ).GE.OSTHR+dlt) THEN
               COMPARE = ABS(1-WORK(LDL-1+IJ)/WORK(LDV-1+IJ))
             ELSE IF((WORK(LDL-1+IJ).GE.OSTHR).AND.
     &               (WORK(LDL-1+IJ).GT.0.0D0)) THEN
               COMPARE = -1.5D0
             ELSE IF((WORK(LDV-1+IJ).GE.OSTHR).AND.
     &               (WORK(LDV-1+IJ).GT.0.0D0)) THEN
               COMPARE = -2.5D0
             END IF

             IF(ABS(COMPARE).GE.TOLERANCE) THEN
               I_PRINT_HEADER = I_PRINT_HEADER + 1
               IF(I_PRINT_HEADER.EQ.1) THEN
                 WRITE(6,*)
                 WRITE(6,*) " Problematic transitions have been found"
                 WRITE(6,*)
                 WRITE(6,39) "From","To","Difference (%)",
     &                       "Osc. st. (len.)","Osc. st. (vel.)"
                 WRITE(6,40)
               END IF
               IF (COMPARE.GE.0.0D0) THEN
                 WRITE(6,33) I,J,COMPARE*100D0,
     &                      WORK(LDL-1+IJ),WORK(LDV-1+IJ)
               ELSE IF (COMPARE.GE.-2.0D0) THEN
                 WRITE(6,36) I,J,WORK(LDL-1+IJ),"below threshold"
               ELSE
                 WRITE(6,37) I,J,"below threshold",WORK(LDV-1+IJ)
               END IF
             END IF
*
            END DO ! L_
          END DO ! K_
          IF(I_PRINT_HEADER.EQ.0) THEN
            WRITE(6,*)
            WRITE(6,*) "No problematic oscillator strengths above "//
     &                 "the tolerance ", TOLERANCE," have been found"
            WRITE(6,*)
          ELSE
            WRITE(6,40)
            WRITE(6,*)
            WRITE(6,*) "Number of problematic transitions = ",
     &                  I_PRINT_HEADER
            WRITE(6,*)
          END IF
         CALL CollapseOutput(0,'Length and velocity gauge comparison '//
     &                         '(spin-free states):')
         WRITE(6,*)
        END IF
*
* Free the memory
*
      CALL GETMEM('DV   ','FREE','REAL',LDV,NSTATE**2)
      CALL GETMEM('DL   ','FREE','REAL',LDL,NSTATE**2)
*
* CALCULATION OF THE QUADRUPOLE TRANSITION STRENGTHS
*
      SECORD = 0
!
! Lazy mans version
!
      NSS = NSTATE
!
! We will first allocate a matrix for the total of the second order wave vector
!
      CALL GETMEM('TOT2K','ALLO','REAL',LTOT2K,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LTOT2K),1)

* Magnetic-Dipole - Magnetic-Dipole transitions
!
! Magnetic-Dipole
!
! DEBUG
!
        IPRDX_TEMP=IPRDX
        IPRDY_TEMP=IPRDY
        IPRDZ_TEMP=IPRDZ
! BEBUG END
        IPRDX=0
        IPRDY=0
        IPRDZ=0

        IFANYD=0
        DO IPROP=1,NPROP
          IF(PNAME(IPROP).EQ.'ANGMOM  ') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZ=IPROP
          END IF
        END DO

        IF(IFANYD.NE.0) THEN
!
! Only print the part calculated
!
        IF(QIALL) THEN
          WRITE(6,*)
          Call CollapseOutput(1,
     &                  'Magnetic-Dipole - Magnetic-Dipole '//
     &                  'transition strengths (spin-free states):')
          WRITE(6,'(3X,A)')
     &                  '----------------------------------'//
     &                  '----------------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30) 'for osc. strength at least',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
        END IF

         ONEOVER6C2=1.0D0/(6.0D0*CONST_C_IN_AU_**2)

         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN
            IJSS=JSS+NSS*(ISS-1)

            DX2=0.0D0
            DY2=0.0D0
            DZ2=0.0D0

            IF(IPRDX.GT.0) DX2=PROP(JSS,ISS,IPRDX)**2
            IF(IPRDY.GT.0) DY2=PROP(JSS,ISS,IPRDY)**2
            IF(IPRDZ.GT.0) DZ2=PROP(JSS,ISS,IPRDZ)**2

            F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F
            IF(ABS(F).GE.OSTHR2) THEN
!            WRITE(6,*) ' value at distance '
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
!
! Debug to move along z. Change DX and DY (1 Aangstrom)
!
!     DO I = 1, 9
!     AA = 1.889726D0*ZVAL(I)*EDIFF!/(2.0D0*CONST_C_IN_AU_)
! z-direction
!     DX2=(PROP(JSS,ISS,IPRDX)
!    &    -AA*PROP(JSS,ISS,IPRDY_TEMP))**2
!     DY2=(PROP(JSS,ISS,IPRDY)
!    &    +AA*PROP(JSS,ISS,IPRDX_TEMP))**2
!           F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
! x-direction
!     DZ2=(PROP(JSS,ISS,IPRDZ)
!    &    +AA*PROP(JSS,ISS,IPRDY_TEMP))**2
!     DY2=(PROP(JSS,ISS,IPRDY)
!    &    -AA*PROP(JSS,ISS,IPRDZ_TEMP))**2
!           F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
! y-direction
!     DX2=(PROP(JSS,ISS,IPRDX)
!    &    +AA*PROP(JSS,ISS,IPRDZ_TEMP))**2
!     DZ2=(PROP(JSS,ISS,IPRDZ)
!    &    -AA*PROP(JSS,ISS,IPRDX_TEMP))**2
!           F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
!           IF(ABS(F).GE.OSTHR2) THEN
!            WRITE(6,*) ' moved value ',ZVAL(I)
!            WRITE(6,'(5X,2I5,5X,G16.8)') ISS,JSS,F
!           END IF
!     END DO
           END IF
          END DO
         END DO
        IF(QIALL) THEN
         WRITE(6,35)

         Call CollapseOutput(0,
     &                  'Magnetic-Dipole - Magnetic-Dipole '//
     &                  'transition strengths (spin-free states):')
        END IF
! Magnetic-dipole - Magnetic-dipole calculated
          SECORD(1) = 1
        END IF

*Electric-Quadrupole Electric-Quadrupole transitions

        IPRDXX=0
        IPRDXY=0
        IPRDXZ=0
        IPRDYY=0
        IPRDYZ=0
        IPRDZZ=0

        IFANYD=0
        DO IPROP=1,NPROP
          IF(PNAME(IPROP).EQ.'MLTPL  2') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDXX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDXY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDXZ=IPROP
           IF(ICOMP(IPROP).EQ.4) IPRDYY=IPROP
           IF(ICOMP(IPROP).EQ.5) IPRDYZ=IPROP
           IF(ICOMP(IPROP).EQ.6) IPRDZZ=IPROP
          END IF
        END DO

        IF(IFANYD.NE.0) THEN
        IF(QIALL) THEN
         WRITE(6,*)
         Call CollapseOutput(1,
     &            'Quadrupole transition strengths (spin-free states):')
         WRITE(6,'(3X,A)')
     &            '---------------------------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30) 'for osc. strength at least',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
        END IF

         ONEOVER10C=1.0D0/(10.0D0*CONST_C_IN_AU_**2)
         ONEOVER30C=ONEOVER10C/3.0D0

         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF3=EDIFF**3
            IJSS=JSS+NSS*(ISS-1)

            DXX=0.0D0
            DYY=0.0D0
            DZZ=0.0D0
            DXY=0.0D0
            DXZ=0.0D0
            DYZ=0.0D0
            IF(IPRDXX.GT.0) DXX=PROP(JSS,ISS,IPRDXX)
            IF(IPRDYY.GT.0) DYY=PROP(JSS,ISS,IPRDYY)
            IF(IPRDZZ.GT.0) DZZ=PROP(JSS,ISS,IPRDZZ)
            IF(IPRDXY.GT.0) DXY=PROP(JSS,ISS,IPRDXY)
            IF(IPRDXZ.GT.0) DXZ=PROP(JSS,ISS,IPRDXZ)
            IF(IPRDYZ.GT.0) DYZ=PROP(JSS,ISS,IPRDYZ)

            DXX2=DXX**2
            DYY2=DYY**2
            DZZ2=DZZ**2
            FXX=ONEOVER30C*EDIFF3*(DXX2)
            FYY=ONEOVER30C*EDIFF3*(DYY2)
            FZZ=ONEOVER30C*EDIFF3*(DZZ2)

            DXY2=DXY**2
            DXZ2=DXZ**2
            DYZ2=DYZ**2
            FXY=ONEOVER10C*EDIFF3*(DXY2)
            FXZ=ONEOVER10C*EDIFF3*(DXZ2)
            FYZ=ONEOVER10C*EDIFF3*(DYZ2)

            DXXDYY=DXX*DYY
            DXXDZZ=DXX*DZZ
            DYYDZZ=DYY*DZZ
            FXXFYY=-ONEOVER30C*EDIFF3*(DXXDYY)
            FXXFZZ=-ONEOVER30C*EDIFF3*(DXXDZZ)
            FYYFZZ=-ONEOVER30C*EDIFF3*(DYYDZZ)

            F =FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
!
! Debug to move along z. Change DZZ, DXZ and DYZ
!
!     DO I = 1, 9
!           DZZ2=(PROP(JSS,ISS,IPRDZZ)+
!    &           1.889726D0*ZVAL(I)*2.0D0*PROP(JSS,ISS,IPRDZ_TEMP))**2
!           DXZ2=(PROP(JSS,ISS,IPRDXZ)+
!    &           1.889726D0*ZVAL(I)*  PROP(JSS,ISS,IPRDX_TEMP))**2
!           DYZ2=(PROP(JSS,ISS,IPRDYZ)+
!    &           1.889726D0*ZVAL(I)*  PROP(JSS,ISS,IPRDY_TEMP))**2
!           FZZ=ONEOVER30C*EDIFF3*(DZZ2)
!           FXZ=ONEOVER10C*EDIFF3*(DXZ2)
!           FYZ=ONEOVER10C*EDIFF3*(DYZ2)
!           DXXDZZ=PROP(JSS,ISS,IPRDXX)*(PROP(JSS,ISS,IPRDZZ)+
!    &             1.889726D0*ZVAL(I)*2.0D0*PROP(JSS,ISS,IPRDZ_TEMP))
!           DYYDZZ=PROP(JSS,ISS,IPRDYY)*(PROP(JSS,ISS,IPRDZZ)+
!    &             1.889726D0*ZVAL(I)*2.0D0*PROP(JSS,ISS,IPRDZ_TEMP))
!           FXXFZZ=-ONEOVER30C*EDIFF3*(DXXDZZ)
!           FYYFZZ=-ONEOVER30C*EDIFF3*(DYYDZZ)
!           F =FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
!           IF(ABS(F).GE.OSTHR2) THEN
!            WRITE(6,*) ' moved value ',ZVAL(I)
!            WRITE(6,'(5X,2I5,5X,G16.8)') ISS,JSS,F
!           END IF
!     END DO
           END IF
          END DO
         END DO
        IF(QIALL) THEN
         WRITE(6,35)

         Call CollapseOutput(0,
     &            'Quadrupole transition strengths (spin-free states):')
        END IF
          SECORD(2) = 1
        END IF

*Electric-Dipole Electric-Octupole transitions

! Octupole
! This is a real symmetric rank 3 tensor so only 10 and not 27 is needed
! The order which comes in
!
! DEBUG
!
       IPRDXX_TEMP=IPRDXX
       IPRDXY_TEMP=IPRDXY
       IPRDXZ_TEMP=IPRDXZ
       IPRDYY_TEMP=IPRDYY
       IPRDYZ_TEMP=IPRDYZ
       IPRDZZ_TEMP=IPRDZZ
! DEBUG END
        IPRDXXX=0 !
        IPRDXXY=0 !
        IPRDXXZ=0 !

!       IPRDXYX=0
!       IPRDXYY=0 ! YYX These are the same due to symmetry
        IPRDXYZ=0 ! Not present

!       IPRDXZX=0
!       IPRDXZY=0
!       IPRDXZZ=0 ! ZZX

!       IPRDYXX=0
!       IPRDYXY=0
!       IPRDYXZ=0

        IPRDYYX=0 ! Taking the XYY order
        IPRDYYY=0 !
        IPRDYYZ=0 !

!       IPRDYZX=0
!       IPRDYZY=0
!       IPRDYZZ=0 ! ZZY

!       IPRDZXX=0
!       IPRDZXY=0
!       IPRDZXZ=0

!       IPRDZYX=0
!       IPRDZYY=0
!       IPRDZYZ=0

        IPRDZZX=0 ! Taking order from XZZ
        IPRDZZY=0 ! Taking order from YZZ
        IPRDZZZ=0 !
! Dipole
        IPRDX=0
        IPRDY=0
        IPRDZ=0


        IFANYD=0
        DO IPROP=1,NPROP
          IF(PNAME(IPROP).EQ.'MLTPL  1') THEN
           IF(ICOMP(IPROP).EQ.1) IPRDX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZ=IPROP
          ELSE IF(PNAME(IPROP).EQ.'MLTPL  3') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDXXX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDXXY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDXXZ=IPROP
           IF(ICOMP(IPROP).EQ.4) IPRDYYX=IPROP ! Changed from XYY
           IF(ICOMP(IPROP).EQ.5) IPRDXYZ=IPROP
           IF(ICOMP(IPROP).EQ.6) IPRDZZX=IPROP ! Changed from XZZ
           IF(ICOMP(IPROP).EQ.7) IPRDYYY=IPROP
           IF(ICOMP(IPROP).EQ.8) IPRDYYZ=IPROP
           IF(ICOMP(IPROP).EQ.9) IPRDZZY=IPROP ! Changed from YZZ
           IF(ICOMP(IPROP).EQ.10) IPRDZZZ=IPROP
          END IF
        END DO
! Sanity check. Only check that dipole are there
! since it will give problems the other way when
! only calculating dipole transitions
        IF(((IPRDXXX.GT.0.OR.IPRDYYX.GT.0.OR.IPRDZZX.GT.0)
     &   .AND.IPRDX.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Octupole'
         CALL ABEND()
        END IF
        IF(((IPRDXXY.GT.0.OR.IPRDYYY.GT.0.OR.IPRDZZY.GT.0)
     &   .AND.IPRDY.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Octupole'
         CALL ABEND()
        END IF
        IF(((IPRDXXZ.GT.0.OR.IPRDYYZ.GT.0.OR.IPRDZZZ.GT.0)
     &   .AND.IPRDZ.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Octupole'
         CALL ABEND()
        END IF

        IF(IFANYD.NE.0) THEN
        IF(QIALL) THEN
         WRITE(6,*)
         Call CollapseOutput(1,
     &                     'Electric-Dipole - Electric-Octupole '//
     &                     'transition strengths (spin-free states):')
         WRITE(6,'(3X,A)') '------------------------------------'//
     &                     '----------------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30) 'for osc. strength at least',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
        END IF

         TWOOVERM45C=-2.0D0/(45.0D0*CONST_C_IN_AU_**2)
         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF3=EDIFF**3
            IJSS=JSS+NSS*(ISS-1)

            DXXXDX=0.0D0
            DYYXDX=0.0D0
            DZZXDX=0.0D0
            IF(IPRDXXX.GT.0) DXXXDX=PROP(JSS,ISS,IPRDXXX)
     &                             *PROP(JSS,ISS,IPRDX)
            IF(IPRDYYX.GT.0) DYYXDX=PROP(JSS,ISS,IPRDYYX)
     &                             *PROP(JSS,ISS,IPRDX)
            IF(IPRDZZX.GT.0) DZZXDX=PROP(JSS,ISS,IPRDZZX)
     &                             *PROP(JSS,ISS,IPRDX)
            FXXX=TWOOVERM45C*EDIFF3*(DXXXDX)
            FYYX=TWOOVERM45C*EDIFF3*(DYYXDX)
            FZZX=TWOOVERM45C*EDIFF3*(DZZXDX)

            DXXYDY=0.0D0
            DYYYDY=0.0D0
            DZZYDY=0.0D0
            IF(IPRDXXY.GT.0) DXXYDY=PROP(JSS,ISS,IPRDXXY)
     &                             *PROP(JSS,ISS,IPRDY)
            IF(IPRDYYY.GT.0) DYYYDY=PROP(JSS,ISS,IPRDYYY)
     &                             *PROP(JSS,ISS,IPRDY)
            IF(IPRDZZY.GT.0) DZZYDY=PROP(JSS,ISS,IPRDZZY)
     &                             *PROP(JSS,ISS,IPRDY)
            FXXY=TWOOVERM45C*EDIFF3*(DXXYDY)
            FYYY=TWOOVERM45C*EDIFF3*(DYYYDY)
            FZZY=TWOOVERM45C*EDIFF3*(DZZYDY)

            DXXZDZ=0.0D0
            DYYZDZ=0.0D0
            DZZZDZ=0.0D0
            IF(IPRDXXZ.GT.0) DXXZDZ=PROP(JSS,ISS,IPRDXXZ)
     &                             *PROP(JSS,ISS,IPRDZ)
            IF(IPRDYYZ.GT.0) DYYZDZ=PROP(JSS,ISS,IPRDYYZ)
     &                             *PROP(JSS,ISS,IPRDZ)
            IF(IPRDZZZ.GT.0) DZZZDZ=PROP(JSS,ISS,IPRDZZZ)
     &                             *PROP(JSS,ISS,IPRDZ)
            FXXZ=TWOOVERM45C*EDIFF3*(DXXZDZ)
            FYYZ=TWOOVERM45C*EDIFF3*(DYYZDZ)
            FZZZ=TWOOVERM45C*EDIFF3*(DZZZDZ)

            F =FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
!
! Debug to move along z. Change DZZX,DZZY,DXXZ,DYYZ and DZZZ
!
!     DO I = 1, 9
!     DZZXDX=(PROP(JSS,ISS,IPRDZZX)+
!    &        1.889726D0*ZVAL(I)*2*PROP(JSS,ISS,IPRDXZ_TEMP)+
!    &       (1.889726D0*ZVAL(I))**2*PROP(JSS,ISS,IPRDX_TEMP))*
!    &        PROP(JSS,ISS,IPRDX)
!           FZZX=TWOOVERM45C*EDIFF3*(DZZXDX)

!     DZZYDY=(PROP(JSS,ISS,IPRDZZY)+
!    &        1.889726D0*ZVAL(I)*2*PROP(JSS,ISS,IPRDYZ_TEMP)+
!    &       (1.889726D0*ZVAL(I))**2*PROP(JSS,ISS,IPRDY_TEMP))*
!    &        PROP(JSS,ISS,IPRDY)
!           FZZY=TWOOVERM45C*EDIFF3*(DZZYDY)

!     DXXZDZ=(PROP(JSS,ISS,IPRDXXZ)+
!    &        1.889726D0*ZVAL(I)*PROP(JSS,ISS,IPRDXX_TEMP))*
!    &        PROP(JSS,ISS,IPRDZ)
!     DYYZDZ=(PROP(JSS,ISS,IPRDYYZ)+
!    &        1.889726D0*ZVAL(I)*PROP(JSS,ISS,IPRDYY_TEMP))*
!    &        PROP(JSS,ISS,IPRDZ)
!     DZZZDZ=(PROP(JSS,ISS,IPRDZZZ)+
!    &        1.889726D0*ZVAL(I)*3*PROP(JSS,ISS,IPRDZZ_TEMP)+
!    &       (1.889726D0*ZVAL(I))**2*3*PROP(JSS,ISS,IPRDZ_TEMP))*
!    &        PROP(JSS,ISS,IPRDZ)
!           FXXZ=TWOOVERM45C*EDIFF3*(DXXZDZ)
!           FYYZ=TWOOVERM45C*EDIFF3*(DYYZDZ)
!           FZZZ=TWOOVERM45C*EDIFF3*(DZZZDZ)

!           F =FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ

!           IF(ABS(F).GE.OSTHR2) THEN
!            WRITE(6,*) ' moved value ',ZVAL(I)
!            WRITE(6,'(5X,2I5,5X,G16.8)') ISS,JSS,F
!           END IF
!     END DO
           END IF
          END DO
         END DO
        IF(QIALL) THEN
         WRITE(6,35)

         Call CollapseOutput(0,
     &                     'Electric-Dipole - Electric-Octupole '//
     &                     'transition strengths (spin-free states):')
        END IF
          SECORD(3) = 1
        END IF
*
*Electric-Dipole - Magnetic-Quadrupole transitions
!
! Magnetic-Quadrupole
        IPRDXX=0
        IPRDXY=0
        IPRDXZ=0

        IPRDYX=0
        IPRDYY=0
        IPRDYZ=0

        IPRDZX=0
        IPRDZY=0
        IPRDZZ=0
! Electric-Dipole
        IPRDX=0
        IPRDY=0
        IPRDZ=0

        IFANYD=0
        DO IPROP=1,NPROP
          IF(PNAME(IPROP).EQ.'MLTPL  1') THEN
           IF(ICOMP(IPROP).EQ.1) IPRDX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZ=IPROP
          ELSE IF(PNAME(IPROP).EQ.'OMQ') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDXX=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDXY=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDXZ=IPROP

           IF(ICOMP(IPROP).EQ.4) IPRDYX=IPROP
           IF(ICOMP(IPROP).EQ.5) IPRDYY=IPROP
           IF(ICOMP(IPROP).EQ.6) IPRDYZ=IPROP

           IF(ICOMP(IPROP).EQ.7) IPRDZX=IPROP
           IF(ICOMP(IPROP).EQ.8) IPRDZY=IPROP
           IF(ICOMP(IPROP).EQ.9) IPRDZZ=IPROP
          END IF
        END DO
! Sanity check. Only check that dipole are there
! since it will give problems the other way when
! only calculating dipole transitions
        IF(((IPRDYZ.GT.0.OR.IPRDZY.GT.0)
     &   .AND.IPRDX.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Quadrupole'
         CALL ABEND()
        END IF
        IF(((IPRDZX.GT.0.OR.IPRDXZ.GT.0)
     &   .AND.IPRDY.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Quadrupole'
         CALL ABEND()
        END IF
        IF(((IPRDXY.GT.0.OR.IPRDYX.GT.0)
     &   .AND.IPRDZ.LE.0)) THEN
         WRITE(6,*) ' Remember to include both Dipole and Quadrupole'
         CALL ABEND()
        END IF

        IF(IFANYD.NE.0) THEN
        IF(QIALL) THEN
          WRITE(6,*)
          Call CollapseOutput(1,
     &                  'Electric-Dipole - Magnetic-Quadrupole '//
     &                  'transition strengths (spin-free states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------'//
     &                  '---------------------------------------'

         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30) 'for osc. strength at least',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         END IF

         ONEOVER9C2=1.0D0/(9.0D0*CONST_C_IN_AU_**2)
         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF2=EDIFF**2
            IJSS=JSS+NSS*(ISS-1)
!
            DXYDZ=0.0D0
            DYXDZ=0.0D0
            IF(IPRDXY.GT.0) DXYDZ=PROP(JSS,ISS,IPRDXY)
     &                           *PROP(JSS,ISS,IPRDZ)
            IF(IPRDXY.GT.0) DYXDZ=PROP(JSS,ISS,IPRDYX)
     &                           *PROP(JSS,ISS,IPRDZ)
            FXY=ONEOVER9C2*EDIFF2*(DXYDZ)
            FYX=-ONEOVER9C2*EDIFF2*(DYXDZ)

            DZXDY=0.0D0
            DXZDY=0.0D0
            IF(IPRDZX.GT.0) DZXDY=PROP(JSS,ISS,IPRDZX)
     &                           *PROP(JSS,ISS,IPRDY)
            IF(IPRDXZ.GT.0) DXZDY=PROP(JSS,ISS,IPRDXZ)
     &                           *PROP(JSS,ISS,IPRDY)
            FZX=ONEOVER9C2*EDIFF2*(DZXDY)
            FXZ=-ONEOVER9C2*EDIFF2*(DXZDY)

            DYZDX=0.0D0
            DZYDX=0.0D0
            IF(IPRDYZ.GT.0) DYZDX=PROP(JSS,ISS,IPRDYZ)
     &                           *PROP(JSS,ISS,IPRDX)
            IF(IPRDZY.GT.0) DZYDX=PROP(JSS,ISS,IPRDZY)
     &                           *PROP(JSS,ISS,IPRDX)
            FYZ=ONEOVER9C2*EDIFF2*(DYZDX)
            FZY=-ONEOVER9C2*EDIFF2*(DZYDX)

            F =FYX+FXY+FZX+FXZ+FYZ+FZY
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
!
! Debug to move along z.
!
!           ONEOVER3C = 1.0D0/(3.0D0*CONST_C_IN_AU_)
!           DYXDZ=(PROP(JSS,ISS,IPRDYX)*ONEOVER3C
!    &           + PROP(JSS,ISS,IPRDXX_TEMP)*ONEOVER3C*EDIFF*1.889726D0)
!    &           * PROP(JSS,ISS,IPRDZ)
!           DXYDZ=(PROP(JSS,ISS,IPRDXY)*ONEOVER3C
!    &           - PROP(JSS,ISS,IPRDYY_TEMP)*ONEOVER3C*EDIFF*1.889726D0)
!    &           * PROP(JSS,ISS,IPRDZ)
!     print*,'YX,XY moved',PROP(JSS,ISS,IPRDYX)
!    &                   + PROP(JSS,ISS,IPRDXX_TEMP)*EDIFF*1.889726D0,
!    &                     PROP(JSS,ISS,IPRDXY)
!    &                   - PROP(JSS,ISS,IPRDYY_TEMP)*EDIFF*1.889726D0
!           FXY=ONEOVER3C*EDIFF2*(DXYDZ)
!           FYX=-ONEOVER3C*EDIFF2*(DYXDZ)

!           DZXDY=PROP(JSS,ISS,IPRDZX)*ONEOVER3C*PROP(JSS,ISS,IPRDY) !independent
!           DXZDY=(PROP(JSS,ISS,IPRDXZ)*ONEOVER3C
!    &           - PROP(JSS,ISS,IPRDYZ_TEMP)*ONEOVER3C*EDIFF*1.889726D0  ! changed from IPRDZY_TEMP to IPRDYZ_TEMP
!    &     + PROP(JSS,ISS,IPRDY_TEMP)*2.0D0*ONEOVER3C*EDIFF*1.889726D0**2)
!    &           * PROP(JSS,ISS,IPRDY) ! skipped magnetic dipole
!     print*,'ZX,XZ moved',PROP(JSS,ISS,IPRDZX),
!    &                     PROP(JSS,ISS,IPRDXZ)
!    &                   - PROP(JSS,ISS,IPRDYZ_TEMP)*EDIFF*1.889726D0
!    &                + PROP(JSS,ISS,IPRDY_TEMP)*2.0D0*EDIFF*1.889726D0**2
!           FZX=ONEOVER3C*EDIFF2*(DZXDY)
!           FXZ=-ONEOVER3C*EDIFF2*(DXZDY)

!           DYZDX=(PROP(JSS,ISS,IPRDYZ)*ONEOVER3C
!    &           + PROP(JSS,ISS,IPRDXZ_TEMP)*ONEOVER3C*EDIFF*1.889726D0 ! changed from IPRDZX_TEMP to IPRDXZ_TEMP
!    &     - PROP(JSS,ISS,IPRDX_TEMP)*2.0D0*ONEOVER3C*EDIFF*1.889726D0**2)
!    &           * PROP(JSS,ISS,IPRDX)
!           DZYDX=PROP(JSS,ISS,IPRDZY)*ONEOVER3C*PROP(JSS,ISS,IPRDX)
!     print*,'YZ,ZY moved',PROP(JSS,ISS,IPRDYZ)
!    &                   + PROP(JSS,ISS,IPRDXZ_TEMP)*EDIFF*1.889726D0
!    &            - PROP(JSS,ISS,IPRDX_TEMP)*2.0D0*EDIFF*1.889726D0**2,
!    &              PROP(JSS,ISS,IPRDZY)
!           FYZ=ONEOVER3C*EDIFF2*(DYZDX)
!           FZY=-ONEOVER3C*EDIFF2*(DZYDX)
! The new diagonal ones?
!           F =FYX+FXY+FZX+FXZ+FYZ+FZY

!           IF(ABS(F).GE.OSTHR2) THEN
!            WRITE(6,*) ' The moved value '
!            WRITE(6,'(5X,2I5,5X,G16.8)') ISS,JSS,F
!           END IF
! End debug
           END IF
          END DO
         END DO
        IF(QIALL) THEN
         WRITE(6,35)

         Call CollapseOutput(0,
     &                  'Electric-Dipole - Magnetic-Quadrupole '//
     &                  'transition strengths (spin-free states):')
        END IF
          SECORD(4) = 1
        END IF
!
! Now write out the total
!
! Add it to the total
!
      I2TOT = 0
      DO I = 1, 4
        IF(SECORD(I).EQ.1) THEN
          I2TOT = I2TOT + 1
        END IF
      END DO
       IF(I2TOT.GE.1) THEN
         IF(SECORD(1).EQ.0)
     &   WRITE(6,*) 'Magnetic-dipole - magnetic-dipole not included'
         IF(SECORD(2).EQ.0)
     &   WRITE(6,*) 'Electric-quadrupole - electric-quadrupole not '//
     &              'included'
         IF(SECORD(3).EQ.0)
     &   WRITE(6,*) 'Electric-dipole - electric-octupole not included'
         IF(SECORD(4).EQ.0)
     &   WRITE(6,*) 'Electric-dipole - magnetic-quadrupole not included'
         iPrint=0
         DO ISS_=1,IEND
          ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
           JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            IJSS=JSS+NSS*(ISS-1)
            F = WORK(LTOT2K-1+IJSS)
            IF(ABS(F).GE.OSTHR2) THEN
            If (iPrint.eq.0) Then
         WRITE(6,*)
         Call CollapseOutput(1,
     &                'Total transition strengths ' //
     &                'for the second-order expansion of the wave ' //
     &                'vector (spin-free states):')
         WRITE(6,'(3X,A)')
     &                '---------------------------'//
     &                '-------------------------------------------'//
     &                '--------------------------'
!
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30) 'for osc. strength at least',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         iPrint=1
             End If
             WRITE(6,33) ISS,JSS,F
            END IF
           END IF
          END DO
         END DO
         If (iPrint.eq.1) Then
         WRITE(6,35)
         Call CollapseOutput(0,
     &                'Total transition strengths ' //
     &                'for the second-order expansion of the wave ' //
     &                'vector (spin-free states):')
         End If
       END IF
! release the memory again
       CALL GETMEM('TOT2K','FREE','REAL',LTOT2K,NSS**2)
!
!
      IF(DOCD) THEN
* Lasse 2019
* New CD here with electric dipole and magnetic-dipole - velocity gauge
        IPRDXD=0
        IPRDYD=0
        IPRDZD=0
        IPRDXM=0
        IPRDYM=0
        IPRDZM=0

        IFANYD=0
        IFANYM=0
        DO IPROP=1,NPROP
          IF (PNAME(IPROP).EQ.'VELOCITY') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDXD=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDYD=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZD=IPROP
          END IF
          IF(PNAME(IPROP).EQ.'ANGMOM  ') THEN
           IFANYM=1
           IF(ICOMP(IPROP).EQ.1) IPRDXM=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDYM=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZM=IPROP
          END IF
        END DO

        IF((IFANYD.NE.0).AND.(IFANYM.NE.0)) THEN
!
! Only print the part calculated
!
          WRITE(6,*)
          Call CollapseOutput(1,
     &                  'Circular Dichroism - velocity gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (spin-free states):')
          WRITE(6,'(3X,A)')
     &                  '------------------------------------'//
     &                  '----------------------------------'//
     &                  '--------------------------------------'
          WRITE(6,*)
          WRITE(6,31) 'From','To','Red. rot. str.'
          WRITE(6,35)
!
         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN

*           R = e^2*hbar/(2*m^2*E) <J|p|I>.<I|l|J>
*             = e^2*hbar/(2*m^2*E) -i*hbar*<J|nabla|I>.-i*hbar*<I|r x nabla|J>
*             = e^2*hbar^3/(2*m^2*E) <J|nabla|I>.<J|r x nabla|I>
            R=0.0D0

            IF((IPRDXM.GT.0).AND.(IPRDXD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDXD)*PROP(JSS,ISS,IPRDXM)
            END IF
            IF((IPRDYM.GT.0).AND.(IPRDYD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDYD)*PROP(JSS,ISS,IPRDYM)
            END IF
            IF((IPRDZM.GT.0).AND.(IPRDZD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDZD)*PROP(JSS,ISS,IPRDZM)
            END IF

            R = R*Half/EDIFF*AU2REDR

            WRITE(6,33) ISS,JSS,R
!
            Call Add_Info('CD_V(SF)',[R],1,6)
           END IF
          END DO
         END DO
         WRITE(6,35)

         Call CollapseOutput(0,
     &                  'Circular Dichroism - velocity gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (spin-free states):')
        END IF
!
* Lasse 2019
* New CD here with electric dipole and magnetic-dipole - mixed gauge
* Usually refered to as the length gauge
        IPRDXD=0
        IPRDYD=0
        IPRDZD=0
        IPRDXM=0
        IPRDYM=0
        IPRDZM=0

        IFANYD=0
        IFANYM=0
        DO IPROP=1,NPROP
          IF (PNAME(IPROP).EQ.'MLTPL  1') THEN
           IFANYD=1
           IF(ICOMP(IPROP).EQ.1) IPRDXD=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDYD=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZD=IPROP
          END IF
          IF(PNAME(IPROP).EQ.'ANGMOM  ') THEN
           IFANYM=1
           IF(ICOMP(IPROP).EQ.1) IPRDXM=IPROP
           IF(ICOMP(IPROP).EQ.2) IPRDYM=IPROP
           IF(ICOMP(IPROP).EQ.3) IPRDZM=IPROP
          END IF
        END DO

        IF((IFANYD.NE.0).AND.(IFANYM.NE.0)) THEN
!
! Only print the part calculated
!
          WRITE(6,*)
          Call CollapseOutput(1,
     &                  'Circular Dichroism - mixed gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (spin-free states):')
          WRITE(6,'(3X,A)')
     &                  '---------------------------------'//
     &                  '----------------------------------'//
     &                  '--------------------------------------'
          WRITE(6,*)
          WRITE(6,*) ' WARNING WARNING WARNING !!! '
          WRITE(6,*)
          WRITE(6,*) ' Circular Dichroism in the mixed gauge '
          WRITE(6,*) ' is NOT origin independent - check your results '
          WRITE(6,*)
          WRITE(6,31) 'From','To','Red. rot. str.'
          WRITE(6,35)
!
         DO ISS_=1,IEND
            ISS=IndexE(ISS_)
          DO JSS_=JSTART,NSS
             JSS=IndexE(JSS_)
           EDIFF=ENERGY(JSS)-ENERGY(ISS)
           IF(EDIFF.GT.0.0D0) THEN

*           R = -i*e^2/(2*m) <J|r|I>.<I|l|J>
*             = -i*e^2/(2*m) <J|r|I>.-i*hbar*<I|r x nabla|J>
*             = e^2*hbar/(2*m) <J|r|I>.<J|r x nabla|I>
            R=0.0D0

            IF((IPRDXM.GT.0).AND.(IPRDXD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDXD)*PROP(JSS,ISS,IPRDXM)
            END IF
            IF((IPRDYM.GT.0).AND.(IPRDYD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDYD)*PROP(JSS,ISS,IPRDYM)
            END IF
            IF((IPRDZM.GT.0).AND.(IPRDZD.GT.0)) THEN
              R=R+PROP(JSS,ISS,IPRDZD)*PROP(JSS,ISS,IPRDZM)
            END IF

            R = R*Half*AU2REDR

            WRITE(6,33) ISS,JSS,R
!
            Call Add_Info('CD_M(SF)',[R],1,6)
           END IF
          END DO
         END DO
         WRITE(6,35)

         Call CollapseOutput(0,
     &                  'Circular Dichroism - mixed gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (spin-free states):')
        END IF
      END IF
* CD end

*
! +++ J. Norell 12/7 - 2018
! Dyson amplitudes for (1-electron) ionization transitions
       IF (DYSO) THEN
        DYSTHR=1.0D-5
        WRITE(6,*)
        CALL CollapseOutput(1,'Dyson amplitudes '//
     &                        '(spin-free states):')
        WRITE(6,'(3X,A)')     '----------------------------'//
     &                        '-------------------'
        IF (DYSTHR.GT.0.0D0) THEN
           WRITE(6,30) 'for Dyson intensities at least',DYSTHR
           WRITE(6,30)
        END IF
        WRITE(6,*) '       From      To        '//
     &   'BE (eV)       Dyson intensity'
        WRITE(6,32)
        FMAX=0.0D0
        DO I_=1,NSTATE
           I=IndexE(I_)
         DO J_=1,NSTATE
            J=IndexE(J_)
          F=DYSAMPS(I,J)*DYSAMPS(I,J)
          EDIFF=AU2EV*(ENERGY(J)-ENERGY(I))
          IF (F.GT.0.00001) THEN
           IF (EDIFF.GT.0.0D0) THEN
            WRITE(6,'(A,I8,I8,F15.3,E22.5)') '    ',I,J,EDIFF,F
           END IF
          END IF
         END DO ! J
        END DO ! I
       END IF
! +++ J. Norell


************************************************************************
*                                                                      *
*     Start of section for transition moments                          *
*                                                                      *
************************************************************************
*
*     Find the section of transition moments in the property list.
*
*     The operator is split in 4 different component, each with three
*     elements corresponding to differentiation in the x, y, and z
*     direction.
*     The B.s term is split into 2 component.
*
************************************************************************
*                                                                      *
*     Computation of the isotropic oscillator strength.                *
*                                                                      *
************************************************************************
*
      If (.Not.Do_TMOM) Go To 900
#define _TIME_TMOM_
#ifdef _TIME_TMOM_
      Call CWTime(TCpu1,TWall1)
#endif
*
*     Here we will use a Lebedev grid to integrate over all possible
*     directions of the wave vector, k. The property integrals will be
*     computed on the fly and traced with the density to generate the
*     corresponding values in the PROP matrix.
*
*     Find the slot on the one-electron file where we will store the
*     on-the-fly generated property integrals.
*
      IPRTMOM(:)=-1
      DO IPROP=1,NPROP
         IF (PNAME(IPROP).EQ.'TMOM  RS') THEN
            IF (IPRTMOM(0+ICOMP(IPROP)).EQ.-1)
     &          IPRTMOM(0+ICOMP(IPROP))=IPROP
         END IF
         IF (PNAME(IPROP).EQ.'TMOM  IS') THEN
            IF (IPRTMOM(3+ICOMP(IPROP)).EQ.-1)
     &          IPRTMOM(3+ICOMP(IPROP))=IPROP
         END IF
         IF (PNAME(IPROP).EQ.'TMOM  RA') THEN
            IF (IPRTMOM(6+ICOMP(IPROP)).EQ.-1)
     &          IPRTMOM(6+ICOMP(IPROP))=IPROP
         END IF
         IF (PNAME(IPROP).EQ.'TMOM  IA') THEN
            IF (IPRTMOM(9+ICOMP(IPROP)).EQ.-1)
     &          IPRTMOM(9+ICOMP(IPROP))=IPROP
         END IF
      ENDDO
      IF (ANY(IPRTMOM.EQ.-1)) RETURN
*
*     Initiate the Seward environment
*
      nDiff=0
      Call IniSew(Info,.FALSE.,nDiff)
*
*     Generate the quadrature points.
*
      If (Do_SK) Then
         nQuad=1
         Call GetMem('SK','ALLO','REAL',ipR,4*nQuad)
         nVec = nk_Vector
      Else
         Call Setup_O()
*        In the spin-free case, oscillator and rotatory strengths for k and -k
*        are equal, so we compute only half the quadrature points and multiply
*        the weights by 2
         Call Do_Lebedev_Sym(L_Eff,nQuad,ipR)
         Call DScal_(nQuad,2.0D0,Work(ipR+3),4)
         nVec = 1
      End If
      If (Do_Pol) Call mma_allocate(pol_Vector,3,nVec*nQuad,Label='POL')
*
*     Get table of content for density matrices.
*
      Call DaName(LuToM,FnToM)
      iDisk=0
      Call iDaFile(LuToM,2,iWork(liTocM),nState*(nState+1)/2,iDisk)
*
      NIP=4+(NBST*(NBST+1))/2
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
      NSCR=(NBST*(NBST+1))/2
      CALL GETMEM('TDMSCR','Allo','Real',LSCR,4*NSCR)
*                                                                      *
************************************************************************
*                                                                      *
*     Transform the transition densities to the new basis.
*
*     This procedure is redundant if the Hamiltonian matrix was
*     diagonal and the energies where in increasing order.
*
      If (.NOT.DIAGONAL) Then
*
      Call mma_Allocate(TDS,4*nSCR,nState*(nState+1)/2,Label='TDS')
      Call FZero(TDS,4*nSCR*nState*(nState+1)/2)
*
*     Loop over all unique TDs and distribute their contributions to the
*     TDs in the new basis.
*
      DO I=1, NSTATE
         DO J= 1, NSTATE
            ISTATE=MAX(i,j)
            JSTATE=MIN(i,j)
            ij=ISTATE*(ISTATE-1)/2+JSTATE
            iDisk=iWork(liTocM+ij-1)
            Get_TDS=.True.
*
*           The reading is postponed until it is really needed
*           (see below).
*           Call dDaFile(LuToM,2,Work(LSCR),4*NSCR,iDisk)
*
*           Loop over the TDs which will be used in the subsequent part
*           of the code.
*
            Do K = JSTART, NSTATE
               C_ik=EigVec(i,K)
               Do L = 1, Min(IEND,K)
                  C_ikjl=C_ik*EigVec(j,L)
                  If (Abs(C_ikjl).gt.1.0D-14) Then
                     If (Get_TDS) Then
                        If (iDisk.gt.0) Then
                           Call dDaFile(LuToM,2,Work(LSCR),4*NSCR,iDisk)
                        Else
                           Call FZero(Work(LSCR),4*NSCR)
                        End If
                        Get_TDS=.False.
                     End If
                     kl=K*(K-1)/2+L
                     Call DaXpY_(4*NSCR,C_ikjl,Work(LSCR),1,
     &                                         TDS(1,kl),1)
                  End If
               End Do
            End Do
*
         END DO
      END DO
*
*     Replace the old TDs with the new ones on the file.
*
      DO I=1, IEND
         DO J=JSTART, NSTATE
            ISTATE=MAX(i,j)
            JSTATE=MIN(i,j)
            ij=ISTATE*(ISTATE-1)/2+JSTATE
            iDisk=iWork(liTocM+ij-1)
            If (iDisk.ge.0) Then
               Call dDaFile(LuToM,1,TDS(1,ij),4*NSCR,iDisk)
            Else
*              This should probably never happen
               dNorm=dDot_(4*NSCR,TDS(1,ij),1,TDS(1,ij),1)
               If (dNorm.gt.0.0d0) Then
                  Write(6,*) 'A transition density matrix should be '//
     &                       'zero, but it is not.'
                  Call AbEnd()
               End If
            End If
         END DO
      END DO
      Call mma_DeAllocate(TDS)
*
      END IF
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _HDF5_
*
*     Allocate vector to store all individual transition moments.
*     We do this for
*     all unique pairs ISO-JSO, iSO=/=JSO (NSS*(NSS-1)/2)
*         all k-vectors (nQuad or nVec)
*             we store:
*                 the weight (1)
*                 the k-vector (3)
*                 the projected transition vector (real and imaginary parts) (2*3)
*
      nIJ=nState*(nState-1)/2
      ip_w       = 1
      ip_kvector = ip_w + 1
      ip_TMR     = ip_kvector + 3
      ip_TMI     = ip_TMR + 3
      nData      = ip_TMI + 3 - 1
      Call mma_allocate(Storage,nData,nQuad,nIJ,nVec,label='Storage')
      Call dCopy_(Size(Storage),[0.0D0],0,Storage,1)
#endif
*MGD create the groups of indices
*Only with reduce loop to make things easier
      TMOgroup=.false.
      ngroup1=IEND
      ngroup2=NSTATE-JSTART+1
      nmax2=1
      IF(REDUCELOOP.and.TMGr_thrs.ge.0.0d0) THEN
        TMOgroup=.true.
        THRS=TMGr_thrs
        i=IndexE(IEND)
        RefEne=0
        TAU=-1
        ngroup2=1
        Do j_=JSTART,NSTATE
           j=IndexE(j_)
           if (ENERGY(J)-Refene.gt.TAU) then
              NGROUP2=NGROUP2+1
              Refene=ENERGY(J)
              ediff=Refene-ENERGY(I)
              TAU=ediff*THRS
           EndIf
        End Do
        Call mma_Allocate(TMOgrp2,NGROUP2,Label='TMOgrp2')
        ngroup2=0
        TAU=-1
        RefEne=0
        Do j_=JSTART,NSTATE
           j=IndexE(j_)
           if (ENERGY(J)-Refene.gt.TAU) then
              NGROUP2=NGROUP2+1
              TMOgrp2(NGROUP2)=J_
              Refene=ENERGY(J)
              ediff=Refene-ENERGY(I)
              TAU=ediff*THRS
           EndIf
        End Do
        TMOgrp2(ngroup2+1)=NSTATE+1
*
        j=IndexE(JSTART)
        Refene=ENERGY(j)
        TAU=-1
        ngroup1=1
        Do i_=IEND,1,-1
          i=IndexE(i_)
          if (Refene-ENERGY(i).gt.TAU) then
            ngroup1=ngroup1+1
            Refene=energy(i)
            ediff=energy(j)-Refene
            Tau=ediff*THRS
          EndIf
        End Do
        Call mma_Allocate(TMOgrp1,NGROUP1,Label='TMOgrp1')
        Ntmp=Ngroup1
        Ngroup1=Ngroup1-1
        Refene=ENERGY(j)
        TAU=-1
        Do i_=IEND,1,-1
          i=IndexE(i_)
          if (Refene-ENERGY(i).gt.TAU) then
            TMOgrp1(ntmp)=i_+1
            ntmp=ntmp-1
            Refene=energy(i)
            ediff=energy(j)-Refene
            Tau=ediff*THRS
          EndIf
        End Do
        TMOgrp1(1)=1
*       write(6,*) (TMOgrp1(i),i=1,ngroup1+1)
*       write(6,*) (TMOgrp2(i),i=1,ngroup2+1)
        maxgrp1=0
        Do i=1,ngroup1
          maxgrp1=max(maxgrp1,TMOgrp1(i+1)-TMOgrp1(i))
        End Do
        maxgrp2=0
        Do i=1,ngroup2
          maxgrp2=max(maxgrp2,TMOgrp2(i+1)-TMOgrp2(i))
        End Do
        nmax2=maxgrp1*maxgrp2
      EndIF
*
*     Array for printing contributions from different directions
*
      CALL GETMEM('RAW   ','ALLO','REAL',LRAW,NQUAD*6*nmax2)
      CALL GETMEM('OSCSTR','ALLO','REAL',LF,2*nmax2)
      CALL GETMEM('MAXMIN','ALLO','REAL',LMAX,8*nmax2)
      LMAX_=0
*
      Do iVec = 1, nVec
         If (Do_SK) Then
            Work(ipR  )=k_Vector(1,iVec)
            Work(ipR+1)=k_Vector(2,iVec)
            Work(ipR+2)=k_Vector(3,iVec)
            Work(ipR+3)=1.0D0   ! Dummy weight
         End If
*
      iPrint=0
      IJSO=0
      Do igrp=1,ngroup1
         Do jgrp=1,ngroup2
*
            If (TMOgroup) Then
              istart_=TMOgrp1(igrp)
              iend_=TMOgrp1(igrp+1)-1
              jstart_=TMOgrp2(jgrp)
              jend_=TMOgrp2(jgrp+1)-1
              EDIFF_=(ENERGY(IndexE(jstart_))+ENERGY(IndexE(jend_))
     &              -ENERGY(IndexE(istart_))-ENERGY(IndexE(iend_)))
     &              *Half
            Else
              istart_=igrp
              iend_=igrp
              jstart_=jgrp+jstart-1
              jend_=jgrp+jstart-1
              EDIFF_=ENERGY(IndexE(jstart_))-ENERGY(IndexE(istart_))
            EndIf
            n12=(iend_-istart_+1)*(jend_-jstart_+1)
            If (ABS(EDIFF_).le.1.0D-8) CYCLE
*
            If (JSTART.eq.1 .AND.  EDIFF_.LT.0.0D0) CYCLE
*
            IJSO=IJSO+1
*
*           The energy difference is used to define the norm of the
*           wave vector.
*
            rkNorm=ABS(EDIFF_)/SPEED_OF_LIGHT
*
*           Iterate over the quadrature points.
*
            CALL DCOPY_(2*n12,[0.0D0],0,WORK(LF),1)
*
*           Initialize output arrays
*
            CALL DCOPY_(NQUAD*6*n12,[0.0D0],0,WORK(LRAW),1)
            CALL DCOPY_(8*n12,[0.0D0],0,WORK(LMAX),1)
*
            Do iQuad = 1, nQuad
               iVec_=(iVec-1)*nQuad+iQuad
*
*              Read or generate the wavevector
*
*              Generate the wavevector associated with this quadrature
*              point and pick up the associated quadrature weight.
*
               UK(1)=Work((iQuad-1)*4  +ipR)
               UK(2)=Work((iQuad-1)*4+1+ipR)
               UK(3)=Work((iQuad-1)*4+2+ipR)
               Wavevector(:)=rkNorm*UK(:)
*
*              Note that the weights are normalized to integrate to
*              4*pi over the solid angles.
*
               Weight=Work((iQuad-1)*4+3+ipR)
               If (.Not.Do_SK) Weight=Weight/(4.0D0*PI)
*
*              Generate the polarization vector
*
               If (Do_Pol) Then
                  pol_Vector(:,iVec_)=
     &               e_Vector-DDot_(3,UK,1,e_Vector,1)*UK
                  rNorm=DDot_(3,pol_Vector(:,iVec_),1,
     &                         pol_Vector(:,iVec_),1)
                  If (rNorm.gt.1.0D-12) Then
                     pol_Vector(:,iVec_)=pol_Vector(:,iVec_)/Sqrt(rNorm)
                  Else
                     pol_Vector(:,iVec_)=0.0D0
                  End If
               End If
*
*              Generate the property integrals associated with this
*              direction of the wave vector k.
*
               iOpt=1
               Call TMOMInt(Wavevector,iOpt)
*
*              Compute the transition property of the property
*              integrals between the two states.
*
               ij_=0
               Do i_=istart_,iend_
                  I=IndexE(I_)
                  MPLET_I=MLTPLT(iWork(lJBNUM+I-1))
                  Do j_=jstart_,jend_
                     J=IndexE(J_)
                     MPLET_J=MLTPLT(iWork(lJBNUM+J-1))
                     EDIFF=ENERGY(J)-ENERGY(I)
                     ij_=ij_+1
                     LFIJ=LF+(ij_-1)*2
C COMBINED SYMMETRY OF STATES:
                     JOB1=iWork(lJBNUM+I-1)
                     JOB2=iWork(lJBNUM+J-1)
                     LSYM1=IRREP(JOB1)
                     LSYM2=IRREP(JOB2)
                     ISY12=MUL(LSYM1,LSYM2)
C THE SYMMETRY CHECK MASK:
                     MASK=2**(ISY12-1)
C ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
C FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
                     IOF=0
                     Call IZERO(IOFF,8)
                     DO ISY1=1,NSYM
                       ISY2=MUL(ISY1,ISY12)
                       IF (ISY1.LT.ISY2) CYCLE
                       IOFF(ISY1)=IOF
                       IOFF(ISY2)=IOF
                       NB1=NBASF(ISY1)
                       NB2=NBASF(ISY2)
                       NB12=NB1*NB2
                       IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
                       IOF=IOF+NB12
                     END DO ! ISY1
C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
*
*           Pick up the transition density between the two states from
*           disc. Generated in PROPER.
*
                     ISTATE=MAX(i,j)
                     JSTATE=MIN(i,j)
                     ij=ISTATE*(ISTATE-1)/2+JSTATE
                     iDisk=iWork(liTocM+ij-1)
                     If (iDisk.gt.0) Then
                        Call dDaFile(LuToM,2,Work(LSCR),4*NSCR,iDisk)
                     Else
                        Call FZero(Work(LSCR),4*NSCR)
                     End If
               DO IPRP = 1,12
                  IPROP=IPRTMOM(IPRP)
                  ITYPE=0
                  IF (PTYPE(IPROP).EQ.'HERMSING') ITYPE=1
                  IF (PTYPE(IPROP).EQ.'ANTISING') ITYPE=2
                  IF (PTYPE(IPROP).EQ.'HERMTRIP') ITYPE=3
                  IF (PTYPE(IPROP).EQ.'ANTITRIP') ITYPE=4
                  LABEL=PNAME(IPROP)
                  Call MK_PROP(PROP,IPROP,ISTATE,JSTATE,LABEL,ITYPE,
     &                         WORK(LIP),NIP,WORK(LSCR),NSCR,
     &                         MASK,ISY12,IOFF)
               END DO ! IPRP
*
*              (1) the oam part
*
*              The contribution to the generalized momentum operator.
*              Note that the integrals contain nabla, but we need p,
*              so we multiply by -i
*
               DO iCar=1,3
                  TM_R(iCar)=+PROP(I,J,IPRTMOM(3+iCar)) ! IS
     &                       +PROP(I,J,IPRTMOM(9+iCar)) ! IA
                  TM_I(iCar)=-PROP(I,J,IPRTMOM(0+iCar)) ! RS
     &                       -PROP(I,J,IPRTMOM(6+iCar)) ! RA
               END DO
*
*              (2) the magnetic-spin part
*
C                 Well the B.S term is overkill, get rid of it.
C                 Why do it when we don't do the L.S-term!
*
*              Finally, evaluate the transition moment from the two
*              different contributions.
*
#ifdef _HDF5_
*              Fix the triangular index because we are not storing the diagonal
               IJSF=IJ-ISTATE+1
               Storage(ip_w,iQuad,IJSF,iVec)=Weight
               Call DCopy_(3,Wavevector,1,
     &                       Storage(ip_kvector,iQuad,IJSF,iVec),1)
               Call DCopy_(3,TM_R,1,Storage(ip_TMR,iQuad,IJSF,iVec),1)
               Call DCopy_(3,TM_I,1,Storage(ip_TMI,iQuad,IJSF,iVec),1)
#endif
*
*              Project out the k direction from the real and imaginary components
*
               Call DaXpY_(3,-DDot_(3,TM_R,1,UK,1),UK,1,TM_R,1)
               Call DaXpY_(3,-DDot_(3,TM_I,1,UK,1),UK,1,TM_I,1)
*
*              Implicitly integrate over all directions of the
*              polarization vector to get the average value.
*
               TM1 = DDot_(3,TM_R,1,TM_R,1)
               TM2 = DDot_(3,TM_I,1,TM_I,1)
               TM_2 = Half*(TM1+TM2)
*
*              Compute maximum and minimum oscillator strengths
*              and the corresponding polarization vectors
*
               If (Do_SK) Then
                  LMAX_ = LMAX+8*(ij_-1)
                  TM3 = DDot_(3,TM_R,1,TM_I,1)
                  Rng = Sqrt((TM1-TM2)**2+4.0D0*TM3**2)
                  Work(LMAX_+0) = TM_2+Half*Rng
                  Work(LMAX_+4) = TM_2-Half*Rng
*                 The direction for the maximum
                  Ang = Half*Atan2(2.0D0*TM3,TM1-TM2)
                  Call daXpY_(3, Cos(Ang),TM_R,1,Work(LMAX_+1),1)
                  Call daXpY_(3, Sin(Ang),TM_I,1,Work(LMAX_+1),1)
*                 Normalize and compute the direction for the minimum
*                 as a cross product with k
                  rNorm = DDot_(3,Work(LMAX_+1),1,Work(LMAX_+1),1)
                  If (rNorm.gt.1.0D-12) Then
                     Call dScal_(3,1.0/Sqrt(rNorm),Work(LMAX_+1),1)
                     Work(LMAX_+5)=Work(LMAX_+2)*UK(3)-
     &                             Work(LMAX_+3)*UK(2)
                     Work(LMAX_+6)=Work(LMAX_+3)*UK(1)-
     &                             Work(LMAX_+1)*UK(3)
                     Work(LMAX_+7)=Work(LMAX_+1)*UK(2)-
     &                             Work(LMAX_+2)*UK(1)
                     rNorm = DDot_(3,Work(LMAX_+5),1,Work(LMAX_+5),1)
                     Call dScal_(3,1.0/Sqrt(rNorm),Work(LMAX_+5),1)
                  Else
                     Call dCopy_(3,[0.0D0],0,Work(LMAX_+1),1)
                     Call dCopy_(3,[0.0D0],0,Work(LMAX_+5),1)
                  End If
               End If
*
*              Oscillator strength for a specific polarization vector
*
               If (Do_Pol) Then
                  TM1 = DDot_(3,TM_R,1,pol_Vector(1,iVec_),1)
                  TM2 = DDot_(3,TM_I,1,pol_Vector(1,iVec_),1)
                  TM_2 = TM1*TM1+TM2*TM2
               End If
*
*              Compute the oscillator strength
*
               F_Temp = 2.0D0*TM_2/EDIFF
               If (Do_SK) Then
                  Work(LMAX_+0) = 2.0D0*Work(LMAX_+0)/EDIFF
                  Work(LMAX_+4) = 2.0D0*Work(LMAX_+4)/EDIFF
               End If
*
*              Compute the rotatory strength
*
               TM_C(1) = TM_R(2)*TM_I(3)-TM_R(3)*TM_I(2)
               TM_C(2) = TM_R(3)*TM_I(1)-TM_R(1)*TM_I(3)
               TM_C(3) = TM_R(1)*TM_I(2)-TM_R(2)*TM_I(1)
               TM_2 = 2.0D0*DDot_(3,TM_C,1,UK,1)
*
*              R = 3/4 * c*hbar^2/DeltaE^2 * (|T^L|^2 - |T^R|^2)
*
               R_Temp=0.75D0*SPEED_OF_LIGHT/EDIFF**2*TM_2
*
*              Now let's convert this to reduced rotational strength
*              (units of 1e-2 debye*Bohr_magneton)
*
               R_Temp=R_Temp*AU2REDR
*
*              Save the raw oscillator and rotatory strengths in a given direction
*
               LRAW_=LRAW+6*NQUAD*(ij_-1)
               WORK(LRAW_+(IQUAD-1)+0*NQUAD) = F_Temp
               WORK(LRAW_+(IQUAD-1)+1*NQUAD) = R_Temp
*
*              Save the direction and weight too
*
               WORK(LRAW_+(IQUAD-1)+2*NQUAD) = UK(1)
               WORK(LRAW_+(IQUAD-1)+3*NQUAD) = UK(2)
               WORK(LRAW_+(IQUAD-1)+4*NQUAD) = UK(3)
               WORK(LRAW_+(IQUAD-1)+5*NQUAD) = Weight
*
*              Compute the oscillator and rotatory strength
*
               Work(LFIJ  ) = Work(LFIJ  ) + Weight * F_Temp
               Work(LFIJ+1) = Work(LFIJ+1) + Weight * R_Temp
                  End Do ! j_
               End Do ! i_
*
            End Do ! iQuad
*
            ij_=0
            Do i_=istart_,iend_
              I=IndexE(I_)
              Do j_=jstart_,jend_
                J=IndexE(J_)
                ij_=ij_+1
                LFIJ=LF+(ij_-1)*2
*
                F=Work(LFIJ)
                R=Work(LFIJ+1)
*
                Call Add_Info('ITMS(SF)',[F],1,6)
                Call Add_Info('ROTS(SF)',[R],1,6)
*
                IF (Do_Pol) THEN
                   LMAX_=LMAX+8*(ij_-1)
                   F_CHECK=ABS(WORK(LMAX_+0))
                ELSE
                   F_CHECK=ABS(F)
                END IF
                IF (F_CHECK.LT.OSTHR) CYCLE
                A =(AFACTOR*EDIFF**2)*F
*
                If (iPrint.eq.0) Then
                   WRITE(6,*)
                   If (Do_SK) Then
                      CALL CollapseOutput(1,
     &                  'Transition moment strengths '//
     &                  '(spin-free states):')
                      WRITE(6,'(3X,A)')
     &                  '----------------------------'//
     &                  '-------------------'
                      If (Do_Pol) Then
                         iVec_=(iVec-1)*nQuad+1
                         WRITE(6,'(4x,a,3F8.4)')
     &                     'Direction of the polarization: ',
     &                     (pol_vector(k,iVec_),k=1,3)
                      Else
                         WRITE(6,'(4x,a)')
     &                     'The oscillator strength is integrated '//
     &                     'over all directions of the polarization '//
     &                     'vector'
                      End If
                      WRITE(6,'(4x,a,3F8.4)')
     &                  'Direction of the k-vector: ',
     &                  (k_vector(k,iVec),k=1,3)
                   Else
                      CALL CollapseOutput(1,
     &                  'Isotropic transition moment strengths '//
     &                  '(spin-free states):')
                       WRITE(6,'(3X,A)')
     &                  '--------------------------------------'//
     &                  '-------------------'
                   End If
                   IF (OSTHR.GT.0.0D0) THEN
                      WRITE(6,30) 'for osc. strength at least',OSTHR
                   END IF
                   WRITE(6,*)
                   If (.NOT.Do_SK) Then
                      WRITE(6,'(4x,a,I4,a)')
     &                  'Integrated over ',nQuad,' directions of the '//
     &                  'wave vector'
                      WRITE(6,'(4x,a)')
     &                  'The oscillator strength is '//
     &                  'integrated over all directions of the polar'//
     &                  'ization vector'
                      WRITE(6,*)
                   End If
                   WRITE(6,39) 'From','To','Osc. strength',
     &                         'Red. rot. str.','Total A (sec-1)'
                   WRITE(6,40)
                   iPrint=1
                END IF
*
*     Regular print
*
                WRITE(6,33) I,J,F,R,A
*
                IF (Do_SK) THEN
                   WRITE(6,50) 'maximum',WORK(LMAX_+0),
     &                'for polarization direction:',
     &                WORK(LMAX_+1),WORK(LMAX_+2),WORK(LMAX_+3)
                   WRITE(6,50) 'minimum',WORK(LMAX_+4),
     &                'for polarization direction:',
     &                WORK(LMAX_+5),WORK(LMAX_+6),WORK(LMAX_+7)
                END IF
*
*     Printing raw (unweighted) and direction for every transition
*
                IF(PRRAW) THEN
                  WRITE(6,*)
                  WRITE(6,*)
                  WRITE(6,41) 'From', 'To', 'Raw osc. str.',
     &                        'Rot. str.','kx','ky','kz'
                  WRITE(6,32)
                  LRAW_=LRAW+6*NQUAD*(ij_-1)
                  DO IQUAD = 1, NQUAD
                    WRITE(6,33) I,J,
     &              WORK(LRAW_+(IQUAD-1)+0*NQUAD),
     &              WORK(LRAW_+(IQUAD-1)+1*NQUAD),
     &              WORK(LRAW_+(IQUAD-1)+2*NQUAD),
     &              WORK(LRAW_+(IQUAD-1)+3*NQUAD),
     &              WORK(LRAW_+(IQUAD-1)+4*NQUAD)
                  END DO
                  WRITE(6,32)
                  WRITE(6,*)
                END IF
*
*     Printing weighted and direction for every transition
*
                IF(PRWEIGHT) THEN
                  WRITE(6,*)
                  WRITE(6,*)
                  WRITE(6,41) 'From', 'To', 'Weig. osc. str.',
     &                        'Rot. str.','kx','ky','kz'
                  WRITE(6,32)
                  LRAW_=LRAW+5*NQUAD*(ij_-1)
                  DO IQUAD = 1, NQUAD
                    Weight=WORK(LRAW+(IQUAD-1)+5*NQUAD)
                    WRITE(6,33) I,J,
     &                WORK(LRAW+(IQUAD-1)+0*NQUAD)*Weight,
     &                WORK(LRAW+(IQUAD-1)+1*NQUAD)*Weight,
     &                WORK(LRAW+(IQUAD-1)+2*NQUAD),
     &                WORK(LRAW+(IQUAD-1)+3*NQUAD),
     &                WORK(LRAW+(IQUAD-1)+4*NQUAD)
                  END DO
                  WRITE(6,32)
                  WRITE(6,*)
                END IF
              End Do
            End Do
*
         END DO
      END DO
*
      If (iPrint.EQ.1) THEN
         WRITE(6,40)
         If (Do_SK) Then
            CALL CollapseOutput(0,
     &                'Transition moment strengths (spin-free states):')
         Else
         CALL CollapseOutput(0,
     &                'Isotropic transition moment strengths '//
     &                '(spin-free states):')
         End If
      END IF
*
      End Do ! iVec
*
#ifdef _HDF5_
#ifdef POINTER_REMAP
      flatStorage(1:SIZE(Storage)) => Storage
#else
      Call C_F_Pointer(C_Loc(Storage), flatStorage, [SIZE(Storage)])
#endif
      Call mh5_put_dset(wfn_sfs_tm,flatStorage)
      Nullify(flatStorage)
      Call mma_deallocate(Storage)
#endif
*
*     Deallocate some arrays.
*
      CALL GETMEM('RAW   ','FREE','REAL',LRAW,NQUAD*5*nmax2)
      CALL GETMEM('TDMSCR','FREE','Real',LSCR,4*NSCR)
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      CALL GETMEM('OSCSTR','FREE','REAL',LF,2*nmax2)
      CALL GETMEM('MAXMIN','FREE','REAL',LMAX,8*nmax2)
      if (TMOgroup) Then
        Call mma_DeAllocate(TMOgrp1)
        Call mma_DeAllocate(TMOgrp2)
      EndIf
      If (Do_Pol) Call mma_deallocate(pol_Vector)
#ifdef _TIME_TMOM_
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for TMOM : ',TCpu2-TCpu1,TWall2-TWall1
#endif
*
*     Do some cleanup
*
      Call DaClos(LuToM)
      If (.NOT.Do_SK) Call Free_O()
      Call Free_Work(ipR)
      Call ClsSew()
*
************************************************************************
*                                                                      *
*     End of section for transition moments                            *
*                                                                      *
************************************************************************
*

 900  CONTINUE
      if(debug_dmrg_rassi_code)then
        write(6,*) 'end of eigctl: BLUBB debug print of property matrix'
        do istate = 1, nstate
        do jstate = 1, nstate
        DO IPROP=1,NPROP
          if(abs(prop(istate,jstate,iprop)) > 1.0d-14)
     &    write(6,*) 'prop(',istate,',',jstate,',',iprop,') = ',
     &                prop(istate,jstate,iprop)
        end do
        end do
        end do
      end if
      Call mma_DeAllocate(IndexE)

      CALL QEXIT(ROUTINE)
      RETURN
30    FORMAT (5X,A,1X,ES15.8)
31    FORMAT (5X,2(1X,A4),6X,A15,1X,A47,1X,A15)
32    FORMAT (5X,95('-'))
33    FORMAT (5X,2(1X,I4),5X,5(1X,ES15.8))
34    FORMAT (5X,2(1X,A4),5X,4(1X,A15),1X,A)
35    FORMAT (5X,31('-'))
36    FORMAT (5X,2(1X,I4),6X,15('-'),1X,ES15.8,1X,A15)
37    FORMAT (5X,2(1X,I4),6X,15('-'),1X,A15,1X,ES15.8)
39    FORMAT (5X,2(1X,A4),5X,3(1X,A15))
40    FORMAT (5X,63('-'))
41    FORMAT (5X,2(1X,A4),5X,5(1X,A15))
42    FORMAT (5X,79('-'))
50    FORMAT (10X,A7,3X,1(1X,ES15.8),5X,A27,3(1X,F7.4))
      END Subroutine EigCtl

