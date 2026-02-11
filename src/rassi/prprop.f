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
#include "macros.fh"
      SUBROUTINE PRPROP(PROP,USOR,USOI,ENSOR,NSS,OVLP,ENERGY,JBNUM,
     &                  EigVec)

      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: SODYSAMPS
      USE kVectors
#ifdef _HDF5_
      USE mh5, ONLY: mh5_put_dset
      use RASSIWfn, only: wfn_sfs_angmom, wfn_sos_angmomr,
     &                    wfn_sos_angmomi, wfn_sos_spinr, wfn_sos_spini,
     &                    wfn_sfs_edipmom, wfn_sfs_amfi,
     &                    wfn_sos_edipmomr, wfn_sos_edipmomi,
     &                    wfn_sos_dys
      use Cntrl, only: RhoDyn
#endif
      use Constants, only: Pi, auTocm, auToeV, auTofs, auTokJ, auToT,
     &                     c_in_au, Debye, gElectron, kBoltzmann, mBohr,
     &                     rNAVO, Half, Two, Three, Zero, One, Six, Ten,
     &                     Nine
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, NPROP, NTS, PRXVE, PRMEE, LPRPR,
     &                 PRMES, IFSO, NSOPR, DIPR, OSThr_DIPR, QIPR,
     &                 OSthr_QIPR, QIALL, RSPR, RSThr, ReduceLoop,
     &                 LoopDivide, LoopMax, Do_SK, PRDIPCOM, Tolerance,
     &                 DoCD, DYSO, Do_TMom, IfGCAL, EPrThr,
     &                 IfvanVleck, TMins, TMaxs, IfGTCALSA, IfGTSHSA,
     &                 MULTIP, IfACAL, IfXCal, BStart, NBSTep, BIncre,
     &                 TStart, nTStep, TIncre, IfMCal, BAngRes, iComp,
     &                 IPUSED, ISOCMP, MLTPLT, PNAME, PNUC, PORIG,
     &                 PTYPE, SOPRNM, SOPRTP

      IMPLICIT NONE
      Integer NSS, JBNUM(NSTATE)
      Real*8 USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
      Real*8 PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &       ENERGY(NSTATE), EigVec(NSTATE,NSTATE)

      Real*8, parameter:: THRSH=1.0D-10
      Character(LEN=1), Parameter :: xyzchr(3)=['x','y','z']
      Integer IPAMFI(3),IPAM(3)
      Real*8 DTENS(3,3),GTENS(3,3),GSTENS(3,3),SOSTERM(9)
      Real*8 TMPMAT(3,3),TMPVEC(3,3),EVR(3),EVI(3)
      COMPLEX*16 ZEKL(2,2,3,NSTATE),GCONT(9,NSTATE)
      COMPLEX*16 DIPSOm(3,NSS,NSS),Z(NSS,NSS),DIPSOn(3,NSS,NSS)
      COMPLEX*16 SPNSFS(3,NSS,NSS)
      REAL*8 GTOTAL(9),ANGMOME(3,NSTATE,NSTATE),ESO(NSS)
      REAL*8 EDIP1MOM(3,NSTATE,NSTATE),AMFIINT(3,NSTATE,NSTATE)
      Real*8 TMPm(NTS)!,TMPf(NTP)
      Real*8 c_1(3,3),c_2(3,3)!,Zstat1m(NTS),Zstat1f(NTP)
      Real*8 curit(3,3),paramt(3,3)
      Real*8 chiT_tens(NTS,3,3)!,PNMRT(NTP,3,3),PNMR(NTP,3,3)
      Real*8 chicuriT_tens(NTS,3,3),chiparamT_tens(NTS,3,3)
      REAL*8 DLTT,Zstat,p_Boltz!,DLTTA
      LOGICAL ISGS(NSS),IFANGM,IFDIP1,IFAMFI
      Real*8 RMAGM(3),Chi(3)
      INTEGER IFUNCT, SECORD(4)
      Complex*16 T0(3), TM1
      REAL*8 COMPARE
      REAL*8 Rtensor(6)
      REAL*8, Allocatable:: SOPRR(:,:), SOPRI(:,:)
#ifdef _HDF5_
      REAL*8, Allocatable:: TMP(:,:,:)
#endif
      Integer, allocatable:: PMAP(:)
! Dipole
      Real*8, allocatable:: DXR(:,:), DXI(:,:),
     &                      DYR(:,:), DYI(:,:),
     &                      DZR(:,:), DZI(:,:)
! Magnetic-Dipole
      Real*8, allocatable:: MDXR(:,:), MDXI(:,:),
     &                      MDYR(:,:), MDYI(:,:),
     &                      MDZR(:,:), MDZI(:,:)
! Electric-Quadrupole
      Real*8, allocatable:: QXXR(:,:), QXXI(:,:),
     &                      QXYR(:,:), QXYI(:,:),
     &                      QXZR(:,:), QXZI(:,:),
     &                      QYYR(:,:), QYYI(:,:),
     &                      QYZR(:,:), QYZI(:,:),
     &                      QZZR(:,:), QZZI(:,:)
! Magnetic-Quadrupole
      Real*8, allocatable:: MQZXR(:,:), MQZXI(:,:),
     &                      MQXZR(:,:), MQXZI(:,:),
     &                      MQXYR(:,:), MQXYI(:,:),
     &                      MQYXR(:,:), MQYXI(:,:),
     &                      MQYZR(:,:), MQYZI(:,:),
     &                      MQZYR(:,:), MQZYI(:,:)
! Octupole
      Real*8, allocatable:: DXXXR(:,:),DXXXI(:,:),
     &                      DXXYR(:,:),DXXYI(:,:),
     &                      DXXZR(:,:),DXXZI(:,:),
     &                      DYYXR(:,:),DYYXI(:,:),
     &                      DYYYR(:,:),DYYYI(:,:),
     &                      DYYZR(:,:),DYYZI(:,:),
     &                      DZZXR(:,:),DZZXI(:,:),
     &                      DZZYR(:,:),DZZYI(:,:),
     &                      DZZZR(:,:),DZZZI(:,:)
! Spin-Magnetic-Dipole
      Real*8, allocatable:: SXR(:,:), SXI(:,:),
     &                      SYR(:,:), SYI(:,:),
     &                      SZR(:,:), SZI(:,:)
! Spin-Magnetic-Quadrupole
      Real*8, allocatable:: SXYR(:,:), SXYI(:,:), SYXR(:,:), SYXI(:,:),
     &                      SYZR(:,:), SYZI(:,:), SZYR(:,:), SZYI(:,:),
     &                      SZXR(:,:), SZXI(:,:), SXZR(:,:), SXZI(:,:)

      Real*8, allocatable:: DV(:,:), DL(:,:), TOT2K(:,:)
      Real*8, Allocatable:: LXI(:,:), LYI(:,:), LZI(:,:)
      Real*8, Allocatable:: ZR(:,:), ZI(:,:)
      Type A2_Array
           Real*8, Pointer:: A2(:,:)
      End Type A2_Array
      Type (A2_array):: pZMR(3), pZMI(3)
      Type (A2_array):: pMR(3), pMI(3)
      Real*8, allocatable, Target:: ZXR(:,:), ZXI(:,:)
      Real*8, allocatable, Target:: ZYR(:,:), ZYI(:,:)
      Real*8, allocatable, Target:: ZZR(:,:), ZZI(:,:)
      Real*8, allocatable, Target:: MXR(:,:), MXI(:,:)
      Real*8, allocatable, Target:: MYR(:,:), MYI(:,:)
      Real*8, allocatable, Target:: MZR(:,:), MZI(:,:)
      Real*8, allocatable, Target:: UZR(:,:), UZI(:,:)

      Real*8, allocatable:: MAGM(:)

      Real*8, Parameter:: AU2J=auTokJ*1.0D3
      Real*8, Parameter:: J2CM=auTocm/AU2J
      Real*8, Parameter:: AU2JTM=(AU2J/auToT)*rNAVO
      Real*8, Parameter:: AU2REDR=2.0D2*Debye

      Real*8, Parameter:: BOLTZ_K=kBoltzmann*J2CM
      Real*8, Parameter:: coeff_chi=0.1D0*rNAVO/kBoltzmann*mBohr**2
      Real*8, Parameter:: FEGVAL=-gElectron
      Real*8, Parameter:: BOLTZ=kBoltzmann/AU2J
      Real*8, Parameter:: Rmu0=4.0D-7*Pi
      Real*8, Parameter:: Two3rds=Two/Three
      Real*8, Parameter:: ONEOVER6C2=One/(Six*c_in_au**2)
      Real*8, Parameter:: ONEOVER10C=One/(Ten*c_in_au**2)
      Real*8, Parameter:: ONEOVER30C=ONEOVER10C/Three
      Real*8, Parameter:: TWOOVERM45C=-Two/(45.0D0*c_in_au**2)
      REAL*8, Parameter:: ONEOVER9C2=One/(Nine*c_in_au**2)

      Integer nCol, iProp,               I, ISTA, IEND, J, ICMP, NPMSIZ,
     &        nMiss, iSOPr, JSTART, JEND, I_Have_DL, I_Have_DV, nVec,
     &      i_print,ISS, JSS, K, I_Print_Header, IfAnyM, IfAnyS, IfAnyQ,
     &        IfAnyO, I2Tot, iAMFIx, iAMFIy, iAMFIz, iXYZ, jXYZ, iState,
     &        MPLET1, jState, MPLET2, iAMx, iAMy, iAMz, MPLET, kXYZ,
     &        iERR, iMLTPL, iStart, iFinal, ijXYZ, IT, IC, JC, KDGN,
     &        ISO, JSO, LMStep, IBStep, ITStep, nPhiStep, nTheStep,
     &        NORIENT, ITHE, iPhiStep, iPhi, IfAnyD, iPrDXY, iPrDXZ,
     &        iPrDYZ, iVec
      Real*8 AFactor, OSthr, OSThr2, EDiff, DX2, DY2, DZ2, FX, FY, FZ,
     &       A, DLT, EDIFF3, DXX2, DYY2, DZZ2, FXX, FYY, FYZ, DXXDYY,
     &       DXXDZZ, DYYDZZ, FXXFYY, FXXFZZ, FYYFZZ, G, DXXXDX, DYYXDX,
     &       DZZXDX, FXXX, FYYY, FZZZ, DXXYDY, DYYYDY, DZZYDY, FXXY,
     &       FZZY, DXXZDZ, DYYZDZ, DZZZDZ, EDIFF2, DXYDZ, DYXDZ,
     &       FYX, DZXDY, DXZDY, FZX, DYZDX, DZYDX, FZY, D_XR, D_YR,
     &       D_ZR, D_XI, D_YI, D_ZI, D_MXR, D_MYR, D_MZR, D_MXI, D_MYI,
     &       D_MZI, RXX, RYY, RZZ, R, PLIMIT, PMAX,
     &       Q_XXR, Q_XYR, Q_XZR, Q_YYR, Q_YZR, Q_ZZR,
     &       Q_XXI, Q_XYI, Q_XZI, Q_YYI, Q_YZI, Q_ZZI,
     &       RXY, RXZ, RYX, RYZ, RZX, RZY,
     &       RXXY, RXXZ, RXYX, RXYZ, RXZX, RXZY, RXYY, RYYX, RYYZ,
     &       RYZX, RYZY, RXZZ, RZZX, RZZY, DysThr, S1, FACT0, FACTP,
     &       FACTM, DTIJ, S2, DELTA, CONTRIB, S, Factor, GTij, GSEnergy,
     &       DLT_E, BFinal, TFinal, B, HZer, T, RKT, RPart, Fact,
     &       rMagm2, rMagMO, GTR, bPhiRes, Phi, Bx, By, Bz, DIPSOM_SA,
     &       EEX, EEY, EEZ, ThreEJ, AX, AY, AZ, F, FZZ, DXY2, DXZ2,
     &       FXY, FXZ, FYYX, FZZX, FXXZ, FYYZ, DYZ2, RYZZ, THE

******************************************************
* printout of properties over the spin-free states
******************************************************

      IF(IPGLOB.LE.0) GOTO 400

      IF( PRXVE.OR.PRMEE ) THEN
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,'(6X,A)') repeat('*',100)
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,A,34X,A,34X,A)')
     &     '*',' Spin-free properties section ','*'
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,A)') repeat('*',100)
      WRITE(6,*)
      WRITE(6,*)
      END IF

* Did the user want printed expectation values?
      IF( PRXVE ) THEN
       Call CollapseOutput(1,'Expectation values')
       WRITE(6,*)
       WRITE(6,*)' ============================================'
       WRITE(6,*)'  EXPECTATION VALUES OF 1-ELECTRON OPERATORS'
       WRITE(6,*)'  FOR THE SPIN-FREE EIGENSTATES:'
       WRITE(6,*)' ============================================'
       WRITE(6,*)' (note: negative sign used for electronic multipoles)'
       WRITE(6,*)
       NCOL=4
       DO IPROP=1,NPROP
        IF(IPUSED(IPROP).EQ.0) GOTO 100

* Skip printing if all the diagonal values are very small
*  (presumed zero for reasons of selection rules)
        PLIMIT=1.0D-10
        PMAX=ZERO


        DO I=1,NSTATE
         PMAX=MAX(PMAX,ABS(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
        END DO
        IF(PMAX.LT.PLIMIT) GOTO 100

        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     *  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
          WRITE(6,'(1X,A,3(1X,ES16.9))')
     *'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,4I17)')
     *'STATE     :',(I,I=ISTA,IEND)
          WRITE(6,*)
          WRITE(6,'(1X,A,4(1X,ES16.9))')
     *'ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1X,ES16.9))')
     *'NUCLEAR   :',(PNUC(IPROP),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1X,ES16.9))')
     *'TOTAL     :',(PROP(I,I,IPROP)+PNUC(IPROP),I=ISTA,IEND)
          WRITE(6,*)
        END DO
 100    CONTINUE
       END DO
       Call CollapseOutput(0,'Expectation values')
       WRITE(6,*)
      END IF

* include nuclear contribution
      DO IPROP=1,NPROP
        DO I=1,NSTATE
          PROP(I,I,IPROP)=PROP(I,I,IPROP)+PNUC(IPROP)
        END DO
      END DO

* Did the user want printed matrix elements?
      IF( PRMEE ) THEN
       Call CollapseOutput(1,'Matrix elements')
       WRITE(6,*)
       WRITE(6,*)' ========================================='
       WRITE(6,*)'  MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
       WRITE(6,*)'  FOR THE SPIN-FREE EIGENSTATES:'
       WRITE(6,*)' ========================================='
       WRITE(6,*)' (including nuclear contrib.)'
       WRITE(6,*)
       WRITE(6,*)' SELECTED PROPERTIES:'
       DO I=1,NPROP,5
         WRITE(6,'(1X,5(A8,1X,I2,4X))')
     &          (PNAME(IPROP),ICOMP(IPROP),IPROP=I,MIN(NPROP,I+4))
       END DO

       NCOL=4
       DO IPROP=1,NPROP
         IF(IPUSED(IPROP).EQ.0) GOTO 200
         WRITE(6,*)
         WRITE(6,'(1X,A,A8,A,I4)')
     *   'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
         WRITE(6,'(1X,A,3(1X,ES16.9))')
     *   'ORIGIN: ',(PORIG(I,IPROP),I=1,3)
         DO ISTA=1,NSTATE,NCOL
           IEND=MIN(NSTATE,ISTA+NCOL-1)
           WRITE(6,*)
           WRITE(6,'(1X,A,I8,3I17)')
     *     ' STATE   ',(I,I=ISTA,IEND)
           WRITE(6,*)
           DO J=1,NSTATE
            WRITE(6,'(1X,I4,6X,4(1X,ES16.9))')
     *      J,(PROP(J,I,IPROP),I=ISTA,IEND)
           END DO
         END DO
 200     CONTINUE
       END DO

       Call CollapseOutput(0,'Matrix elements')
       WRITE(6,*)
      END IF
C Added by Ungur Liviu on 04.11.2009.
C Addition of ANGMOM to Runfile.

      IFANGM=.FALSE.
      IFDIP1=.FALSE.
      IFAMFI=.FALSE.
      ANGMOME(:,:,:) = 0.0D0
      EDIP1MOM(:,:,:) = 0.0D0
      AMFIINT(:,:,:) = 0.0D0
      DO IPROP=1,NPROP
         IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
            IFANGM=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  ANGMOME(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
               ENDDO
            ENDDO
         ENDIF
c add dipole moment integrals:
         IF(PNAME(IPROP).EQ.'MLTPL  1'.AND.
     &      PTYPE(IPROP).EQ.'HERMSING') THEN
            IFDIP1=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  EDIP1MOM(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
               ENDDO
            ENDDO
         ENDIF
c add spin-orbit AMFI integrals:
         IF(PNAME(IPROP)(1:8).EQ.'AMFI    ') THEN
            IFAMFI=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  AMFIINT(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      IF(IFANGM) CALL Put_dArray('ANGM_SINGLE',ANGMOME,3*NSTATE*NSTATE)
      IF(IFDIP1) CALL Put_dArray('DIP1_SINGLE',EDIP1MOM,3*NSTATE*NSTATE)
      IF(IFAMFI) CALL Put_dArray('AMFI_SINGLE',AMFIINT,3*NSTATE*NSTATE)
#ifdef _HDF5_
      call mma_allocate(TMP,NSTATE,NSTATE,3,Label='TMP')
      if (IFANGM) then
        do i=1,3
          TMP(:,:,i) = ANGMOME(i,:,:)
        end do
        call mh5_put_dset(wfn_sfs_angmom,TMP,[NSTATE,NSTATE,3],[0,0,0])
      end if
      if (IFDIP1) then
        do i=1,3
          TMP(:,:,i) = EDIP1MOM(i,:,:)
        end do
        call mh5_put_dset(wfn_sfs_edipmom,TMP,[NSTATE,NSTATE,3],[0,0,0])
      end if
      if (IFAMFI) then
        do i=1,3
          TMP(:,:,i) = AMFIINT(i,:,:)
        end do
        call mh5_put_dset(wfn_sfs_amfi,TMP,[NSTATE,NSTATE,3],[0,0,0])
      end if
      call mma_deallocate(TMP)
#endif
*******************************************************
* printout of properties over the spin-orbit states
*******************************************************
c If PRPR requested, print the spin matrices
#ifdef _HDF5_
      IF (LPRPR.OR.PRMES) THEN
#else
      IF (LPRPR) THEN
#endif
         Call mma_Allocate(SOPRR,NSS,NSS,Label='SOPRR')
         Call mma_Allocate(SOPRI,NSS,NSS,Label='SOPRI')
         DO ICMP=1,3
            SOPRR(:,:)=0.0D0
            SOPRI(:,:)=0.0D0
            IF (ICMP.EQ.2) THEN
              CALL SMMAT(PROP,SOPRI,NSS,0,ICMP)
            ELSE
              CALL SMMAT(PROP,SOPRR,NSS,0,ICMP)
            END IF
            CALL ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
#ifdef _HDF5_
            Call mh5_put_dset(wfn_sos_spinr,
     &                        SOPRR,[NSS,NSS,1],[0,0,ICMP-1])
            Call mh5_put_dset(wfn_sos_spini,
     &                        SOPRI,[NSS,NSS,1],[0,0,ICMP-1])
#endif
            IF (LPRPR) CALL PRCMAT3(NSS,SOPRR,SOPRI,ICMP)
         END DO
         Call mma_deallocate(SOPRR)
         Call mma_deallocate(SOPRI)
      END IF

      IF(.not.IFSO) GOTO 300
      NPMSIZ=NSOPR
      IF(NSOPR.EQ.0) GOTO 300

      IF( PRMES ) THEN
* match the SO property list to the SF property list
       CALL mma_allocate(PMAP,NPMSIZ,Label='PMap')
       NMISS=0
       DO ISOPR=1,NSOPR
        PMAP(ISOPR)=0
        DO IPROP=1,NPROP
         IF(PNAME(IPROP).EQ.SOPRNM(ISOPR).AND.
     &     ICOMP(IPROP).EQ.ISOCMP(ISOPR)) THEN
          PMAP(ISOPR)=IPROP
          GOTO 10
         END IF
        END DO
        NMISS=NMISS+1
 10     CONTINUE
       END DO

c check for inconsistencies
       IF(NMISS.GT.0) THEN
         Call WarningMessage(1,'Missing data integrals.')
         WRITE(6,*)'WARNING: You have requested matrix elements'
         WRITE(6,*)'over spin states of some operators. The present'
         WRITE(6,*)'code uses matrix elements computed over spin-free'
         WRITE(6,*)'states to compute those over spin states.'
         WRITE(6,*)'Matrix elements of the following operator(s)'
         WRITE(6,*)'were never computed and must be skipped.'
         WRITE(6,*)'   (If you need these properties, change the'
         WRITE(6,*)'    input to SEWARD and recompute.)'
         DO ISOPR=1,NSOPR
          IF(PMAP(ISOPR).EQ.0)
     &       WRITE(6,*)'Property:',SOPRNM(ISOPR),
     &                 '      Component:',ISOCMP(ISOPR)
         END DO
        WRITE(6,*)
       END IF
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A,34X,A)')
     &     '*','Spin-orbit properties section ','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,*)
       WRITE(6,*)

       Call CollapseOutput(1,'Matrix elements over SO states')
       WRITE(6,*)
       WRITE(6,*)' ========================================='
       WRITE(6,*)'  MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
       WRITE(6,*)'  FOR THE SPIN-ORBIT EIGENSTATES:'
       WRITE(6,*)' ========================================='
       WRITE(6,*)
       WRITE(6,*)' SELECTED PROPERTIES:'
       DO I=1,NPROP,5
         WRITE(6,'(1X,5(A8,1X,I2,4X))')
     &          (SOPRNM(ISOPR),ISOCMP(ISOPR),ISOPR=I,MIN(NSOPR,I+4))
       END DO

C Remove zeroes to make SOPRNM and ISOCMP lists contiguous. New NSOPR.
       ISOPR=0
       DO I=1,NSOPR
        IPROP=PMAP(I)
        IF(IPROP.GT.0) THEN
         ISOPR=ISOPR+1
         SOPRNM(ISOPR)=SOPRNM(I)
         ISOCMP(ISOPR)=ISOCMP(I)
        END IF
       END DO
       CALL mma_deallocate(PMAP)
       NSOPR=ISOPR

       Call mma_Allocate(SOPRR,NSS,NSS,Label='SOPRR')
       Call mma_Allocate(SOPRI,NSS,NSS,Label='SOPRI')
C Print out the matrix elements:
       NCOL=4
       DO ISOPR=1,NSOPR
        WRITE(6,*)
        WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',SOPRNM(ISOPR),'   COMPONENT:',ISOCMP(ISOPR)
CIFG  should print the origin, but where is it stored (for SO properties)?
        SOPRR(:,:)=0.0D0
        SOPRI(:,:)=0.0D0

        CALL SMMAT(PROP,SOPRR,NSS,ISOPR,0)
        CALL ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
        CALL PRCMAT(NSS,SOPRR,SOPRI)
C prpr keyword: Print selected spin-orbit properties to ext. data files
        IF((LPRPR).AND.(SOPRNM(ISOPR)(1:5).EQ.'MLTPL')) THEN
          IF(SOPRTP(ISOPR).EQ.'HERMSING') THEN
            CALL PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
          ENDIF
        ELSE IF((LPRPR).AND.(SOPRNM(ISOPR)(1:6).EQ.'ANGMOM')) THEN
          IF(SOPRTP(ISOPR).EQ.'ANTISING') THEN
            CALL PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
          ENDIF
        ELSE IF((LPRPR).AND.(SOPRNM(ISOPR)(1:8).EQ.'VELOCITY')) THEN
          IF(SOPRTP(ISOPR).EQ.'ANTISING') THEN
            CALL PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
          ENDIF
        ELSE IF((LPRPR).AND.(SOPRNM(ISOPR)(1:5).EQ.'MLTPV')) THEN
          IF(SOPRTP(ISOPR).EQ.'ANTISING') THEN
            CALL PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
          ENDIF
! prpr end
        ENDIF

#ifdef _HDF5_
        IF( SOPRNM(ISOPR)(1:6) .EQ.'ANGMOM') THEN
           Call mh5_put_dset(wfn_sos_angmomr,
     $                SOPRR,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
           Call mh5_put_dset(wfn_sos_angmomi,
     $                SOPRI,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
        ENDIF

        IF( (SOPRNM(ISOPR)(1:8) .EQ.'MLTPL  1').AND.
     &      (SOPRTP(ISOPR).EQ.'HERMSING') ) THEN
           Call mh5_put_dset(wfn_sos_edipmomr,
     $                SOPRR,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
           Call mh5_put_dset(wfn_sos_edipmomi,
     $                SOPRI,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
        ENDIF
#endif

       END DO
       Call mma_deallocate(SOPRR)
       Call mma_deallocate(SOPRI)
       Call CollapseOutput(0,'Matrix elements over SO states')
       WRITE(6,*)

      END IF

 300  CONTINUE

 400  CONTINUE

******************************************************
* printout of special properties
******************************************************

       ! AFACTOR = 2*pi*e^2*E_h^2 / eps_0*m_e*c^3*h^2
       ! numerically: 2/c^3 (in a.u. of time ^ -1)
       AFACTOR = 2.0D0/c_in_au**3/(auTofs*1.0D-15)

      IF (IPGLOB.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(6X,A)') repeat('*',100)
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A,34X,A,34X,A)')
     &       '*','  Special properties section  ','*'
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A)') repeat('*',100)
        WRITE(6,*)
        WRITE(6,*)
      END IF

C Compute transition strengths for spin-orbit states:
      IF(.not.IFSO) GOTO 500
*
* Initial setup for both dipole, quadrupole etc. and exact operator
*
C printing threshold
!     IF(IPGLOB.eq.2) OSTHR=1.0D-8 ! first order
!     IF(IPGLOB.eq.2) OSTHR2=1.0D-12 ! second order (weaker)
!     IF(IPGLOB.gt.2) OSTHR=0.0D0
!     IF(IPGLOB.gt.2) OSTHR2=0.0D0
      OSTHR=1.0D-5
      OSTHR2=1.0D-5
      IF(DIPR) OSTHR = OSTHR_DIPR
      IF(DIPR) WRITE(6,30) 'Dipole printing threshold changed to ',OSTHR
! Again to avoid total negative transition strengths
      IF(QIPR) OSTHR = OSTHR_QIPR
      IF(QIPR) THEN
        WRITE(6,49)  'Printing threshold changed to ',OSTHR,
     &              'since quadrupole threshold is given '
      END IF
      IF(QIPR) OSTHR2 = OSTHR_QIPR
      IF(QIPR) THEN
      WRITE(6,30) 'Quadrupole printing threshold changed to ',OSTHR2
      END IF
      IF(QIALL) WRITE(6,*) ' Will write all quadrupole contributions '

! Rotatory strength threshold
      IF(RSPR) THEN
        WRITE(6,30) 'Rotatory strength printing threshold changed '//
     &             'to ',RSTHR
      ELSE
        RSTHR = 1.0D-07 !Default
      END IF
!
!     Reducing the loop over states - good for X-rays
!     At the moment memory is not reduced
!
      JEND = NSS
      IF(REDUCELOOP) THEN
        IEND = LOOPDIVIDE
        JSTART = LOOPDIVIDE+1
        IF(LOOPMAX.GT.0) JEND=MIN(NSS, LOOPDIVIDE + LOOPMAX)
      ELSE
        IEND = NSS
        JSTART = 1
      END IF
!
      IF (IPGLOB.GE.1) THEN
!
!     Initialize arrays for indentifying problematic transitions
!     These stores all dipole oscillator strengths in
!     length and velocity gauge for a later comparison.
!
      CALL mma_allocate(DL,NSS,NSS,Label='DL')
      CALL mma_allocate(DV,NSS,NSS,Label='DV')
      DL(:,:)=0.0D0
      DV(:,:)=0.0D0
      I_HAVE_DL = 0
      I_HAVE_DV = 0

*Electric-Dipole Electric-Dipole transitions


        If (Do_SK) Then
           nVec = nk_Vector
        Else
           nVec = 1
        End If
*
        Call Allocate_and_Load_electric_dipoles(IFANYD)

        IF(IFANYD.NE.0) THEN
*
           Do iVec = 1, nVec
*
         i_Print=0

         DO ISS=1,IEND
          DO JSS=JSTART,JEND
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF).LE.1.0D-8) CYCLE
           IF(EDIFF.GT.0.0D0) THEN
            T0(1)=CMPLX(DXR(JSS,ISS),DXI(JSS,ISS),kind=8)
            T0(2)=CMPLX(DYR(JSS,ISS),DYI(JSS,ISS),kind=8)
            T0(3)=CMPLX(DZR(JSS,ISS),DZI(JSS,ISS),kind=8)
            If (Do_SK) Then
               TM1=k_vector(1,iVec)*T0(1)+
     &             k_vector(2,iVec)*T0(2)+
     &             k_vector(3,iVec)*T0(3)
               T0(1) = T0(1) - TM1 * k_vector(1,iVec)
               T0(2) = T0(2) - TM1 * k_vector(2,iVec)
               T0(3) = T0(3) - TM1 * k_vector(3,iVec)
            End If
            DX2=ABS(CONJG(T0(1))*T0(1))
            DY2=ABS(CONJG(T0(2))*T0(2))
            DZ2=ABS(CONJG(T0(3))*T0(3))
            FX=Two3rds*EDIFF*(DX2)
            FY=Two3rds*EDIFF*(DY2)
            FZ=Two3rds*EDIFF*(DZ2)
            F =FX+FY+FZ
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
! Store dipole oscillator strength
            DL(JSS,ISS) = F
            IF(ABS(F).GE.OSTHR) THEN
              If (i_Print.eq.0) Then
                 i_Print=1

! Print full COMPLEX transition dipole moment vectors?
! J. Norell 7/5 2020
         IF(PRDIPCOM) THEN

          Call CollapseOutput(1,
     & 'Complex transition dipole vectors (SO states):')
          WRITE(6,'(3X,A)') '----------------------------------------'
          IF(OSTHR.GT.0.0D0) THEN
           WRITE(6,30)'   for osc. strength at least ',OSTHR
           WRITE(6,*)
          END IF
          WRITE(6,*) '     From   To',
     &'       Re(Dx)       Im(Dx)',
     &'       Re(Dy)       Im(Dy)',
     &'       Re(Dz)       Im(Dz)'
          WRITE(6,32)
          GOTO 137 ! Skip past "regular" print
         END IF
! END print COMPLEX vectors

         Call CollapseOutput(1,
     &                     'Dipole transition strengths (SO states):')
         WRITE(6,'(3X,A)') '----------------------------------------'
         IF(OSTHR.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR
          WRITE(6,*)
         END IF
         If (Do_SK) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F10.6,a)')
     &            'Direction of the k-vector: ',
     &             (k_vector(k,iVec),k=1,3),' (a.u.)'
            WRITE(6,'(4x,a)')
     &            'The light is assumed to be unpolarized.'
            WRITE(6,*)
         End If
         WRITE(6,31) 'From','To','Osc. strength',
     &               'Einstein coefficients Ax, Ay, Az (sec-1)   ',
     &               'Total A (sec-1)'
         WRITE(6,32)
               End If

! Print full COMPLEX transition dipole moment vectors?
137      IF(PRDIPCOM) THEN
             WRITE(6,'(5X,I5,I5,A,A,ES12.3,A,ES12.3,A,ES12.3,A,ES12.3,A,
     &       ES12.3,A,ES12.3)')
     &       ISS,JSS,'   ',
     &       ' ',REAL(T0(1)),' ',AIMAG(T0(1)),
     &       ' ',REAL(T0(2)),' ',AIMAG(T0(2)),
     &       ' ',REAL(T0(3)),' ',AIMAG(T0(3))
        ELSE
             WRITE(6,33) ISS,JSS,F,AX,AY,AZ,A ! "Regular" print instead
         END IF
! END print COMPLEX vectors

            END IF
            Call Add_Info('TMS(SO,Len)',[F],1,6)
            IF(PRDIPCOM) Call Add_Info('TVC(SO,Len)',[DX2+DY2+DZ2],1,6)

           END IF
          END DO
         END DO

         If (i_Print.eq.1) THEN
           WRITE(6,32)
           Call CollapseOutput(0,
     &                     'Dipole transition strengths (SO states):')
           WRITE(6,*)
         END IF
*
         End Do ! iVec
*
         I_HAVE_DL = 1
        END IF

        Call Deallocate_electric_dipoles()

*       Now the same in velocity representation

*
        If (Do_SK) Then
           nVec = nk_Vector
        Else
           nVec = 1
        End If
*
        Call Allocate_and_Load_velocities(IFANYD)

        IF(IFANYD.NE.0) THEN
*
        Do iVec = 1, nVec
*
         i_Print=0

         DO ISS=1,IEND
          DO JSS=JSTART,JEND
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF).LE.1.0D-8) CYCLE
           IF(EDIFF.GT.0.0D0) THEN
            T0(1)=CMPLX(DXR(JSS,ISS),DXI(JSS,ISS),kind=8)
            T0(2)=CMPLX(DYR(JSS,ISS),DYI(JSS,ISS),kind=8)
            T0(3)=CMPLX(DZR(JSS,ISS),DZI(JSS,ISS),kind=8)
            If (Do_SK) Then
               TM1=k_vector(1,iVec)*T0(1)+
     &             k_vector(2,iVec)*T0(2)+
     &             k_vector(3,iVec)*T0(3)
               T0(1) = T0(1) - TM1 * k_vector(1,iVec)
               T0(2) = T0(2) - TM1 * k_vector(2,iVec)
               T0(3) = T0(3) - TM1 * k_vector(3,iVec)
            End If
            DX2=ABS(CONJG(T0(1))*T0(1))
            DY2=ABS(CONJG(T0(2))*T0(2))
            DZ2=ABS(CONJG(T0(3))*T0(3))
            FX=Two3rds*(DX2)/EDIFF
            FY=Two3rds*(DY2)/EDIFF
            FZ=Two3rds*(DZ2)/EDIFF
            F =FX+FY+FZ
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
! Store dipole oscillator strength
            DV(JSS,ISS) = F
            IF(ABS(F).GE.OSTHR) THEN
              If (i_Print.eq.0) Then
                 i_Print=1
         Call CollapseOutput(1,
     &                     'Velocity transition strengths (SO states):')
         WRITE(6,'(3X,A)') '------------------------------------------'
         IF(OSTHR.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR
          WRITE(6,*)
         END IF
         If (Do_SK) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F10.6,a)')
     &            'Direction of the k-vector: ',
     &             (k_vector(k,iVec),k=1,3),' (a.u.)'
            WRITE(6,'(4x,a)')
     &            'The light is assumed to be unpolarized.'
            WRITE(6,*)
         End If
         WRITE(6,31) 'From','To','Osc. strength',
     &               'Einstein coefficients Ax, Ay, Az (sec-1)   ',
     &               'Total A (sec-1)'
         WRITE(6,32)
               END IF
             WRITE(6,33) ISS,JSS,F,AX,AY,AZ,A
            END IF
            Call Add_Info('TMS(SO,Vel)',[F],1,6)
           END IF
          END DO
         END DO

         If (i_Print.eq.1) THEN
           WRITE(6,32)
           Call CollapseOutput(0,
     &                     'Velocity transition strengths (SO states):')
           WRITE(6,*)
         END IF
*
         End Do ! iVec
*
         I_HAVE_DV = 1
        END IF

        Call Deallocate_electric_dipoles()

!
!      Compare oscillator strengths in length and velocity gauge
!      All differences in oscillator strengths above the tolerance
!      of 0.1 (10 percent) will be printed.
!
       IF(I_HAVE_DL.EQ.1.AND.I_HAVE_DV.EQ.1) THEN
         CALL CollapseOutput(1,'Length and velocity gauge comparison '//
     &                         '(SO states):')
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
         WRITE(6,49) "All dipole oscillator differences above the "//
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
          DO I=1,IEND
            DO J=JSTART,JEND
               EDIFF=ENSOR(J)-ENSOR(I)
               IF(EDIFF.LT.0.0D0) CYCLE
             COMPARE=0.0D0
             dlt=1.0D-18 ! Add small value to avoid zero divide.
             IF(DL(J,I).GE.OSTHR+dlt .AND.
     &          DV(J,I).GE.OSTHR+dlt) THEN
               COMPARE = ABS(1-DL(J,I)/DV(J,I))
             ELSE IF((DL(J,I).GE.OSTHR+dlt).AND.
     &               (DL(J,I).GT.0.0D0)) THEN
               COMPARE = -1.5D0
             ELSE IF((DV(J,I).GE.OSTHR+dlt).AND.
     &               (DV(J,I).GT.0.0D0)) THEN
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
                 WRITE(6,38) I,J,COMPARE*100D0,
     &                    DL(J,I),DV(J,I)
               ELSE IF (COMPARE.GE.-2.0D0) THEN
                 WRITE(6,36) I,J,DL(J,I),"below threshold"
               ELSE
                 WRITE(6,37) I,J,"below threshold",DV(J,I)
               END IF
             END IF
            END DO
          END DO
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
     &                         '(SO states):')
         WRITE(6,*)
        END IF
*
* Free the memory
*
      CALL mma_deallocate(DL)
      CALL mma_deallocate(DV)
!
! We will first allocate a matrix for the total of the second order wave vector
!
        CALL mma_allocate(TOT2K,NSS,NSS,Label='TOT2K')
        TOT2K(:,:)=0.0D0
!
! Checking if all are in
        SECORD = 0

* Magnetic-Dipole - Magnetic-Dipole transitions and
* Spin-Magnetic-Dipole - Spin-Magnetic-Dipole transitions
!
! I will not separate these for SO states since there would then be
! M^2 + Ms^2 + 2*MMs (three terms to be programmed)
! M^2 and Ms^2 can be calculated separately but the cross term not directly
!
! Magnetic-Dipole
        Call Allocate_and_Load_Magnetic_Dipoles(IFANYM)
! Spin-Magnetic-Dipole ---- notice the S
        Call Allocate_and_Load_Spin_Magnetic_dipoles(IFANYS)

        IF(IFANYM.NE.0.OR.IFANYS.NE.0) THEN
!
! Only print the part calculated
!
         IF(QIALL) THEN
         IF(IFANYM.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Magnetic-dipole - magnetic-dipole and '//
     &                  'spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------'//
     &                  '--------------------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYM.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(1,
     &                  'Magnetic-dipole - magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '----------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYM.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------------'//
     &                  '---------------------------------'
         END IF
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         END IF


! Magnetic-Dipole

! Spin-Magnetic-Dipole

         g = FEGVAL
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN

            DX2=(MDXI(JSS,ISS)+g*SXR(JSS,ISS))**2
     &         +(MDXR(JSS,ISS)-g*SXI(JSS,ISS))**2
            DY2=(MDYI(JSS,ISS)+g*SYR(JSS,ISS))**2
     &         +(MDYR(JSS,ISS)-g*SYI(JSS,ISS))**2
            DZ2=(MDZI(JSS,ISS)+g*SZR(JSS,ISS))**2
     &         +(MDZR(JSS,ISS)-g*SZI(JSS,ISS))**2

            F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
! Add it to the total
            TOT2K(JSS,ISS) = TOT2K(JSS,ISS) + F
            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

       IF(QIALL) THEN
         WRITE(6,35)
         IF(IFANYM.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Magnetic-dipole - magnetic-dipole and '//
     &                  'spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         ELSE IF(IFANYM.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(0,
     &                  'Magnetic-dipole - magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         ELSE IF(IFANYM.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Spin-magnetic-dipole - Spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         END IF
        END IF
        SECORD(1) = 1
        END IF

! Magnetic-Dipole
        Call Deallocate_Magnetic_dipoles()

! Spin-Magnetic-Dipole
        Call Deallocate_Spin_Magnetic_dipoles()

*Electric-Quadrupole Electric-Quadrupole transitions

        Call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)

        IF(IFANYQ.NE.0) THEN
        IF(QIALL) THEN
         Call CollapseOutput(1,
     &                 'Quadrupole transition strengths (SO states):')
         WRITE(6,'(3X,A)')
     &                 '--------------------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         END IF


         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN
!
! D should be purely real since D is a real symmetric matrix
!
            EDIFF3=EDIFF**3

            DXX2=QXXR(JSS,ISS)**2+QXXI(JSS,ISS)**2
            DYY2=QYYR(JSS,ISS)**2+QYYI(JSS,ISS)**2
            DZZ2=QZZR(JSS,ISS)**2+QZZI(JSS,ISS)**2
            FXX=ONEOVER30C*EDIFF3*(DXX2)
            FYY=ONEOVER30C*EDIFF3*(DYY2)
            FZZ=ONEOVER30C*EDIFF3*(DZZ2)

            DXY2=QXYR(JSS,ISS)**2+QXYI(JSS,ISS)**2
            DXZ2=QXZR(JSS,ISS)**2+QXZI(JSS,ISS)**2
            DYZ2=QYZR(JSS,ISS)**2+QYZI(JSS,ISS)**2
            FXY=ONEOVER10C*EDIFF3*(DXY2)
            FXZ=ONEOVER10C*EDIFF3*(DXZ2)
            FYZ=ONEOVER10C*EDIFF3*(DYZ2)

            DXXDYY=QXXR(JSS,ISS)*QYYR(JSS,ISS)
     &            +QXXI(JSS,ISS)*QYYI(JSS,ISS)
            DXXDZZ=QXXR(JSS,ISS)*QZZR(JSS,ISS)
     &            +QXXI(JSS,ISS)*QZZI(JSS,ISS)
            DYYDZZ=QYYR(JSS,ISS)*QZZR(JSS,ISS)
     &            +QYYI(JSS,ISS)*QZZI(JSS,ISS)
            FXXFYY=-ONEOVER30C*EDIFF3*(DXXDYY)
            FXXFZZ=-ONEOVER30C*EDIFF3*(DXXDZZ)
            FYYFZZ=-ONEOVER30C*EDIFF3*(DYYDZZ)

            F =FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
! Add it to the total
            TOT2K(JSS,ISS) = TOT2K(JSS,ISS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

        IF(QIALL) THEN
         WRITE(6,35)
         Call CollapseOutput(0,
     &                 'Quadrupole transition strengths (SO states):')
         WRITE(6,*)
        END IF
        SECORD(2) = 1
        END IF

        Call Deallocate_Electric_Quadrupoles()

*Electric-Dipole Electric-Octupole transitions

! Octupole
         Call Allocate_and_Load_Octupoles(IFANYO)
! Dipole
         Call Allocate_and_Load_electric_dipoles(IFANYD)

        IF(IFANYD.NE.0.AND.IFANYO.NE.0) THEN
        IF(QIALL) THEN
         Call CollapseOutput(1,
     &                     'Electric-dipole - electric-octupole '//
     &                     'transition strengths (SO states):')
         WRITE(6,'(3X,A)') '------------------------------------'//
     &                     '---------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         END IF

         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF3=EDIFF**3

            DXXXDX=DXXXR(JSS,ISS)*DXR(JSS,ISS)
     &            +DXXXI(JSS,ISS)*DXI(JSS,ISS)
            DYYXDX=DYYXR(JSS,ISS)*DXR(JSS,ISS)
     &            +DYYXI(JSS,ISS)*DXI(JSS,ISS)
            DZZXDX=DZZXR(JSS,ISS)*DXR(JSS,ISS)
     &            +DZZXI(JSS,ISS)*DXI(JSS,ISS)
            FXXX=TWOOVERM45C*EDIFF3*(DXXXDX)
            FYYX=TWOOVERM45C*EDIFF3*(DYYXDX)
            FZZX=TWOOVERM45C*EDIFF3*(DZZXDX)

            DXXYDY=DXXYR(JSS,ISS)*DYR(JSS,ISS)
     &            +DXXYI(JSS,ISS)*DYI(JSS,ISS)
            DYYYDY=DYYYR(JSS,ISS)*DYR(JSS,ISS)
     &            +DYYYI(JSS,ISS)*DYI(JSS,ISS)
            DZZYDY=DZZYR(JSS,ISS)*DYR(JSS,ISS)
     &            +DZZYI(JSS,ISS)*DYI(JSS,ISS)
            FXXY=TWOOVERM45C*EDIFF3*(DXXYDY)
            FYYY=TWOOVERM45C*EDIFF3*(DYYYDY)
            FZZY=TWOOVERM45C*EDIFF3*(DZZYDY)

            DXXZDZ=DXXZR(JSS,ISS)*DZR(JSS,ISS)
     &            +DXXZI(JSS,ISS)*DZI(JSS,ISS)
            DYYZDZ=DYYZR(JSS,ISS)*DZR(JSS,ISS)
     &            +DYYZI(JSS,ISS)*DZI(JSS,ISS)
            DZZZDZ=DZZZR(JSS,ISS)*DZR(JSS,ISS)
     &            +DZZZI(JSS,ISS)*DZI(JSS,ISS)
            FXXZ=TWOOVERM45C*EDIFF3*(DXXZDZ)
            FYYZ=TWOOVERM45C*EDIFF3*(DYYZDZ)
            FZZZ=TWOOVERM45C*EDIFF3*(DZZZDZ)

            F =FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
! Add it to the total
            TOT2K(JSS,ISS) = TOT2K(JSS,ISS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

        IF(QIALL) THEN
         WRITE(6,35)
         Call CollapseOutput(0,
     &                     'Electric-dipole - electric-octupole '//
     &                     'transition strengths (SO states):')
         WRITE(6,*)
        END IF
        SECORD(3) = 1
        END IF

        Call Deallocate_Octupoles()

        Call Deallocate_electric_dipoles()

*Electric-Dipole - Magnetic-Quadrupole transitions and
*Electric-Dipole - Spin-Magnetic-Quadrupole transitions
!
! Again I will just include the spin-term so both terms are calculated
! (Can also be done separately)
! DM + DMs
!
! Magnetic-Quadrupole
! Spin-Magnetic-Quadrupole
! Spin-Magnetic-Quadrupole = M^s_ab = r_b * s_a

! Electric-Dipole
        Call Allocate_and_Load_electric_dipoles(IFANYD)
        IF(IFANYD.NE.0) THEN
! Magnetic-Quadrupole
        Call Allocate_and_Load_Magnetic_Quadrupoles(IFANYD)
! Spin-Magnetic-Quadrupole
        Call Allocate_and_Load_Spin_Magnetic_Quadrupoles(IFANYS)

        IF(IFANYD.NE.0.OR.IFANYS.NE.0) THEN
        IF(QIALL) THEN
         IF(IFANYD.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Electric-dipole - magnetic-quadrupole and '//
     &                  'electric-dipole - spin-magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '------------------------------------------'//
     &                  '-------------------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYD.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(1,
     &                  'Electric-dipole - magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYD.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Electric-dipole - spin-magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '-------------------------------------------'//
     &                  '---------------------------------'
         END IF

         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,30)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,31) 'From','To','Osc. strength'
         WRITE(6,35)
         END IF

         g = FEGVAL*Three/Two ! To remove the 2/3 factor in ONEOVER9C2
         g = g*2.0d0 ! Seem to be needed to agree with the exact term,
                     ! needs to be looked further into!
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF2=EDIFF**2
!
! Since the Spin-Magnetic-Quadrupole is made from the multiplication of two complex integrals we have
! M^s = (a+ib)(c+id) = ac-bd + i(ad+bc) hence the long expressions below
! Also, since the magnetic quadrupole terms are real and the electric dipole are imaginary
! we multiply the real components of MQ with the imaginary of the dipole term, and vice versa.
! However, the spin y component is imaginary
!
!                  Magnetic-Quadrupole   Spin-Magnetic-Quadrupole
!                  Electric-Dipole
            DXYDZ=-((MQXYI(JSS,ISS) + g*SXYI(JSS,ISS)) * DZI(JSS,ISS))
     &           +(( MQXYR(JSS,ISS) + g*SXYR(JSS,ISS)) * DZR(JSS,ISS))
            DYXDZ=-((MQYXI(JSS,ISS) + g*SYXR(JSS,ISS)) * DZI(JSS,ISS))
     &           +(( MQYXR(JSS,ISS) + g*SYXI(JSS,ISS)) * DZR(JSS,ISS))
            FXY=ONEOVER9C2*EDIFF2*(DXYDZ)
            FYX=-ONEOVER9C2*EDIFF2*(DYXDZ)

            DZXDY=-((MQZXI(JSS,ISS) + g*SZXR(JSS,ISS)) * DYI(JSS,ISS))
     &            +((MQZXR(JSS,ISS) + g*SZXI(JSS,ISS)) * DYR(JSS,ISS))
            DXZDY=-((MQXZI(JSS,ISS) + g*SXZR(JSS,ISS)) * DYI(JSS,ISS))
     &            +((MQXZR(JSS,ISS) + g*SXZI(JSS,ISS)) * DYR(JSS,ISS))
            FZX=ONEOVER9C2*EDIFF2*(DZXDY)
            FXZ=-ONEOVER9C2*EDIFF2*(DXZDY)

            DYZDX=-((MQYZI(JSS,ISS) + g*SYZR(JSS,ISS)) * DXI(JSS,ISS))
     &            +((MQYZR(JSS,ISS) + g*SYZI(JSS,ISS)) * DXR(JSS,ISS))
            DZYDX=-((MQZYI(JSS,ISS) + g*SZYI(JSS,ISS)) * DXI(JSS,ISS))
     &           + ((MQZYR(JSS,ISS) + g*SZYR(JSS,ISS)) * DXR(JSS,ISS))
            FYZ=ONEOVER9C2*EDIFF2*(DYZDX)
            FZY=-ONEOVER9C2*EDIFF2*(DZYDX)

            F =FYX+FXY+FZX+FXZ+FYZ+FZY
! Add it to the total
            TOT2K(JSS,ISS) = TOT2K(JSS,ISS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,33) ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

        IF(QIALL) THEN
         WRITE(6,35)
         IF(IFANYD.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Electric-dipole - magnetic-quadrupole and '//
     &                  'electric-dipole - spin-magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
         WRITE(6,*)
         ELSE IF(IFANYD.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(0,
     &                  'Electric-dipole - magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
         WRITE(6,*)
         ELSE IF(IFANYD.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Electric-dipole - spin-magnetic-quadrupole '//
     &                  'transition strengths (SO states):')
         WRITE(6,*)
         END IF
        END IF
        SECORD(4) = 1
        END IF

! Magnetic-Quadrupole
        Call Deallocate_Magnetic_Quadrupoles()
! Spin-Magnetic-Quadrupole
        Call Deallocate_Spin_Magnetic_Quadrupoles()
        END IF
! Electric-Dipole
        Call Deallocate_electric_dipoles()
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
         i_Print=0
         DO ISS=1,IEND
          DO JSS=JSTART,JEND
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN
!
            F = TOT2K(JSS,ISS)
            IF(ABS(F).GE.OSTHR2) THEN
             IF(i_Print.eq.0) THEN
              i_Print=1
              Call CollapseOutput(1,
     &                  'Second-order contribution to the '//
     &                  'transition strengths (SO states):')
              WRITE(6,'(3X,A)')
     &                  '---------------------------------'//
     &                  '---------------------------------'
!
              IF(OSTHR2.GT.0.0D0) THEN
               WRITE(6,30)'   for osc. strength at least ',OSTHR2
               WRITE(6,*)
              END IF
              WRITE(6,31) 'From','To','Osc. strength'
              WRITE(6,35)
             END IF
             WRITE(6,33) ISS,JSS,F
             Call Add_Info('TMS(SO,2nd)',[F],1,6)
            END IF
           END IF
          END DO
         END DO
         If (i_Print.eq.1) THEN
           WRITE(6,35)
           Call CollapseOutput(0,
     &                  'Second-order contribution to the '//
     &                  'transition strengths (SO states):')
           WRITE(6,*)
         END IF
       END IF
! release the memory again
       CALL mma_deallocate(TOT2K)

      END IF
!
!
      IF(DOCD) THEN
* Lasse 2019
* New CD here with electric dipole and magnetic-dipole - velocity gauge


! Electric dipole (linear momentum, p)
        Call Allocate_and_Load_velocities(IFANYD)

! Magnetic-Dipole (angular momentum, l = r x p)
        Call Allocate_and_Load_Magnetic_Dipoles(IFANYM)

        IF((IFANYD.NE.0).AND.(IFANYM.NE.0)) THEN

! Spin-Magnetic-Dipole
         Call Allocate_and_Load_Spin_Magnetic_Dipoles(IFANYS)

! Electric quadrupole (r:p+p:r)

         Call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)
!
! Only print the part calculated
!
         WRITE(6,*)
         Call CollapseOutput(1,
     &                 'Circular Dichroism - velocity gauge '//
     &                 'Electric-Dipole - Magnetic-Dipole '//
     &                 'rotatory strengths (SO states):')
         WRITE(6,'(3X,A)')
     &                 '------------------------------------'//
     &                 '----------------------------------'//
     &                 '-------------------------------'
         IF (DO_SK) THEN
           WRITE(6,30) 'For red. rot. strength at least',RSTHR
         ELSE
           WRITE(6,30) 'For isotropic red. rot. strength at least',RSTHR
         END IF
         WRITE(6,*)
*
         If (Do_SK.AND.(IFANYQ.NE.0)) Then
            nVec = nk_Vector
         Else
            nVec = 1
         End If
*
         Do iVec = 1, nVec
*
         If (Do_SK.AND.(IFANYQ.NE.0)) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F10.6)')
     &         'Direction of the k-vector: ',
     &          (k_vector(k,iVec),k=1,3)
            WRITE(6,*)
            WRITE(6,31) 'From','To','Red. rot. str.'
         Else
            WRITE(6,31) 'From','To','Red. rot. str.'
            IF (IFANYQ.NE.0)
     &         WRITE(6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
         End If
         WRITE(6,35)
!
         g = FEGVAL
         DO ISS=1,IEND
          DO JSS=JSTART,JEND
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN

! These are all complex quantities, and their products are complex too,
! but eventually every piece will be:
!   <I|a|J> <J|b|I> + <I|b|J> <J|b|I>
! for a and b Hermitian operators, so it will be reduced to:
!   2*(Re(A_ij)*Re(B_ij)+Im(A_ij)*Im(B_ij))
! and the imaginary parts of the products can be ignored.

!           Note p = -i*hbar*nabla
            D_XR= DXI(JSS,ISS)
            D_YR= DYI(JSS,ISS)
            D_ZR= DZI(JSS,ISS)
            D_XI=-DXR(JSS,ISS)
            D_YI=-DYR(JSS,ISS)
            D_ZI=-DZR(JSS,ISS)
!           Note r x p = -i*hbar * (r x nabla)
            D_MXR= MDXI(JSS,ISS)+g*SXR(JSS,ISS)
            D_MYR= MDYI(JSS,ISS)+g*SYR(JSS,ISS)
            D_MZR= MDZI(JSS,ISS)+g*SZR(JSS,ISS)
            D_MXI=-MDXR(JSS,ISS)+g*SXI(JSS,ISS)
            D_MYI=-MDYR(JSS,ISS)+g*SYI(JSS,ISS)
            D_MZI=-MDZR(JSS,ISS)+g*SZI(JSS,ISS)

*           R = 1/3 tr(Rtensor)
            RXX=D_XR*D_MXR+D_XI*D_MXI
            RYY=D_YR*D_MYR+D_YI*D_MYI
            RZZ=D_ZR*D_MZR+D_ZI*D_MZI
            IF (ABS(EDIFF)>1.0D-8) THEN
               R = Half/EDIFF*AU2REDR*(RXX+RYY+RZZ)
            ELSE
               R = ZERO
            END IF
            WRITE(6,43) '1/3 Tr(RTensor): ',R
*
* Compute full rotatory strength tensor
* (see Hansen and Bak, 10.1021/jp001899+)
*
            IF (IFANYQ.NE.0) THEN
!            Note r:p+p:r = -i*hbar * (r:nabla+nabla:r)
             Q_XXR= QXXI(JSS,ISS)
             Q_XYR= QXYI(JSS,ISS)
             Q_XZR= QXZI(JSS,ISS)
             Q_YYR= QYYI(JSS,ISS)
             Q_YZR= QYZI(JSS,ISS)
             Q_ZZR= QZZI(JSS,ISS)
             Q_XXI=-QXXR(JSS,ISS)
             Q_XYI=-QXYR(JSS,ISS)
             Q_XZI=-QXZR(JSS,ISS)
             Q_YYI=-QYYR(JSS,ISS)
             Q_YZI=-QYZR(JSS,ISS)
             Q_ZZI=-QZZR(JSS,ISS)
             RXY=D_XR*D_MYR+D_XI*D_MYI
             RXZ=D_XR*D_MZR+D_XI*D_MZI
             RYX=D_YR*D_MXR+D_YI*D_MXI
             RYZ=D_YR*D_MZR+D_YI*D_MZI
             RZX=D_ZR*D_MXR+D_ZI*D_MXI
             RZY=D_ZR*D_MYR+D_ZI*D_MYI
             RXXY=Q_XXR*D_YR+Q_XXI*D_YI
             RXXZ=Q_XXR*D_ZR+Q_XXI*D_ZI
             RXYX=Q_XYR*D_XR+Q_XYI*D_XI
             RXYZ=Q_XYR*D_ZR+Q_XYI*D_ZI
             RXZX=Q_XZR*D_XR+Q_XZI*D_XI
             RXZY=Q_XZR*D_YR+Q_XZI*D_YI
             RXYY=Q_XYR*D_YR+Q_XYI*D_YI
             RYYX=Q_YYR*D_XR+Q_YYI*D_XI
             RYYZ=Q_YYR*D_ZR+Q_YYI*D_ZI
             RYZX=Q_YZR*D_XR+Q_YZI*D_XI
             RYZY=Q_YZR*D_YR+Q_YZI*D_YI
             RXZZ=Q_XZR*D_ZR+Q_XZI*D_ZI
             RYZZ=Q_YZR*D_ZR+Q_YZI*D_ZI
             RZZX=Q_ZZR*D_XR+Q_ZZI*D_XI
             RZZY=Q_ZZR*D_YR+Q_ZZI*D_YI
             ! xx, xy, xz, yy, yz, zz
             Rtensor(1) =  0.75D0 *(RYY+RZZ + (RXYZ-RXZY))
             Rtensor(2) = -0.375D0*(RXY+RYX + (RXXZ+RYZY-RXZX-RYYZ))
             Rtensor(3) = -0.375D0*(RXZ+RZX + (RXYX+RZZY-RXXY-RYZZ))
             Rtensor(4) =  0.75D0 *(RXX+RZZ + (RYZX-RXYZ))
             Rtensor(5) = -0.375D0*(RYZ+RZY + (RYYX+RXZZ-RXYY-RZZX))
             Rtensor(6) =  0.75D0 *(RXX+RYY + (RXZY-RYZX))
             If (ABS(EDIFF)>1.0D-8) Then
                CALL DSCAL_(6,AU2REDR/EDIFF,Rtensor,1)
             ELSE
                Rtensor(:)=ZERO
             END IF
             IF (Do_SK) THEN
              ! k^T R k
              R = k_vector(1,iVec)**2*Rtensor(1)+
     &            k_vector(2,iVec)**2*Rtensor(4)+
     &            k_vector(3,iVec)**2*Rtensor(6)+
     &            2.0D0*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+
     &            2.0D0*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+
     &            2.0D0*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
             ELSE
                WRITE(6,43) 'tensor: ',Rtensor(:)
             END IF
            END IF
*
            IF (ABS(R).GT.RSTHR) THEN
              WRITE(6,33) ISS,JSS,R
            END IF
!
            Call Add_Info('CD_V(SO)',[R],1,6)
           END IF
          END DO
         END DO

         WRITE(6,35)
         End Do

         Call Deallocate_Spin_Magnetic_Dipoles()

         Call Deallocate_Electric_Quadrupoles()

         Call CollapseOutput(0,
     &                  'Circular Dichroism - velocity gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (SO states):')
        END IF

        Call Deallocate_electric_dipoles()

        Call Deallocate_magnetic_dipoles()

* Lasse 2019
* New CD here with electric dipole and magnetic-dipole - mixed gauge

! Electric dipole (r)
        Call Allocate_and_Load_electric_dipoles(IFANYD)
! Magnetic-Dipole (angular momentum, l = r x p)
        Call Allocate_and_Load_Magnetic_dipoles(IFANYM)

        IF((IFANYD.NE.0).AND.(IFANYM.NE.0)) THEN

! Spin-Magnetic-Dipole
         Call Allocate_and_Load_Spin_Magnetic_Dipoles(IFANYS)

! Electric quadrupole (r:r)
         Call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)
!
! Only print the part calculated
!
         WRITE(6,*)
         Call CollapseOutput(1,
     &                 'Circular Dichroism - mixed gauge '//
     &                 'Electric-Dipole - Magnetic-Dipole '//
     &                 'rotatory strengths (SO states):')
         WRITE(6,'(3X,A)')
     &                 '---------------------------------'//
     &                 '----------------------------------'//
     &                 '-------------------------------'
         WRITE(6,*)
         WRITE(6,*) ' WARNING WARNING WARNING !!! '
         WRITE(6,*)
         WRITE(6,*) ' Circular Dichroism in the mixed gauge '
         WRITE(6,*) ' is NOT origin independent - check your results '
         IF (DO_SK) THEN
           WRITE(6,30) 'For red. rot. strength at least',RSTHR
         ELSE
           WRITE(6,30) 'For isotropic red. rot. strength at least',RSTHR
         END IF
         WRITE(6,*)
*
         If (Do_SK.AND.(IFANYQ.NE.0)) Then
            nVec = nk_Vector
         Else
            nVec = 1
         End If
*
         Do iVec = 1, nVec
*
         If (Do_SK.AND.(IFANYQ.NE.0)) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F10.6)')
     &         'Direction of the k-vector: ',
     &          (k_vector(k,iVec),k=1,3)
            WRITE(6,*)
            WRITE(6,31) 'From','To','Red. rot. str.'
         Else
            WRITE(6,31) 'From','To','Red. rot. str.'
            IF (IFANYQ.NE.0)
     &         WRITE(6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
         End If
         WRITE(6,35)
!
         g = FEGVAL
         DO ISS=1,IEND
          DO JSS=JSTART,JEND
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF)<1.0D-8) Cycle
           IF(EDIFF.GT.0.0D0) THEN

! These are all complex quantities, and their products are complex too,
! but eventually every piece will be:
!   <I|a|J> <J|b|I> + <I|b|J> <J|b|I>
! for a and b Hermitian operators, so it will be reduced to:
!   2*(Re(A_ij)*Re(B_ij)+Im(A_ij)*Im(B_ij))
! and the imaginary parts of the products can be ignored.

            D_XR=DXR(JSS,ISS)
            D_YR=DYR(JSS,ISS)
            D_ZR=DZR(JSS,ISS)
            D_XI=DXI(JSS,ISS)
            D_YI=DYI(JSS,ISS)
            D_ZI=DZI(JSS,ISS)
!           Note r x p = -i*hbar * (r x nabla),
!           but we will need i * (r x p) = hbar * r x nabla
            D_MXI=MDXI(JSS,ISS)+g*SXR(JSS,ISS)
            D_MYI=MDYI(JSS,ISS)+g*SYR(JSS,ISS)
            D_MZI=MDZI(JSS,ISS)+g*SZR(JSS,ISS)
            D_MXR=MDXR(JSS,ISS)+g*SXI(JSS,ISS)
            D_MYR=MDYR(JSS,ISS)+g*SYI(JSS,ISS)
            D_MZR=MDZR(JSS,ISS)+g*SZI(JSS,ISS)

*           R = 1/3 tr(Rtensor)
            RXX=D_XR*D_MXR+D_XI*D_MXI
            RYY=D_YR*D_MYR+D_YI*D_MYI
            RZZ=D_ZR*D_MZR+D_ZI*D_MZI
            R = Half*AU2REDR*(RXX+RYY+RZZ)
*
* Compute full rotatory strength tensor
* (see Hansen and Bak, 10.1021/jp001899+)
*
            IF (IFANYQ.NE.0) THEN
             Q_XXR=QXXR(JSS,ISS)
             Q_XYR=QXYR(JSS,ISS)
             Q_XZR=QXZR(JSS,ISS)
             Q_YYR=QYYR(JSS,ISS)
             Q_YZR=QYZR(JSS,ISS)
             Q_ZZR=QZZR(JSS,ISS)
             Q_XXI=QXXI(JSS,ISS)
             Q_XYI=QXYI(JSS,ISS)
             Q_XZI=QXZI(JSS,ISS)
             Q_YYI=QYYI(JSS,ISS)
             Q_YZI=QYZI(JSS,ISS)
             Q_ZZI=QZZI(JSS,ISS)
             RXY=D_XR*D_MYR+D_XI*D_MYI
             RXZ=D_XR*D_MZR+D_XI*D_MZI
             RYX=D_YR*D_MXR+D_YI*D_MXI
             RYZ=D_YR*D_MZR+D_YI*D_MZI
             RZX=D_ZR*D_MXR+D_ZI*D_MXI
             RZY=D_ZR*D_MYR+D_ZI*D_MYI
             RXXY=Q_XXR*D_YR+Q_XXI*D_YI
             RXXZ=Q_XXR*D_ZR+Q_XXI*D_ZI
             RXYX=Q_XYR*D_XR+Q_XYI*D_XI
             RXYZ=Q_XYR*D_ZR+Q_XYI*D_ZI
             RXZX=Q_XZR*D_XR+Q_XZI*D_XI
             RXZY=Q_XZR*D_YR+Q_XZI*D_YI
             RXYY=Q_XYR*D_YR+Q_XYI*D_YI
             RYYX=Q_YYR*D_XR+Q_YYI*D_XI
             RYYZ=Q_YYR*D_ZR+Q_YYI*D_ZI
             RYZX=Q_YZR*D_XR+Q_YZI*D_XI
             RYZY=Q_YZR*D_YR+Q_YZI*D_YI
             RXZZ=Q_XZR*D_ZR+Q_XZI*D_ZI
             RYZZ=Q_YZR*D_ZR+Q_YZI*D_ZI
             RZZX=Q_ZZR*D_XR+Q_ZZI*D_XI
             RZZY=Q_ZZR*D_YR+Q_ZZI*D_YI
             ! xx, xy, xz, yy, yz, zz
             Rtensor(1) =  0.75D0 *(RYY+RZZ+EDIFF*(RXYZ-RXZY))
             Rtensor(2) = -0.375D0*(RXY+RYX+EDIFF*(RXXZ+RYZY-RXZX-RYYZ))
             Rtensor(3) = -0.375D0*(RXZ+RZX+EDIFF*(RXYX+RZZY-RXXY-RYZZ))
             Rtensor(4) =  0.75D0 *(RXX+RZZ+EDIFF*(RYZX-RXYZ))
             Rtensor(5) = -0.375D0*(RYZ+RZY+EDIFF*(RYYX+RXZZ-RXYY-RZZX))
             Rtensor(6) =  0.75D0 *(RXX+RYY+EDIFF*(RXZY-RYZX))
             CALL DSCAL_(6,AU2REDR,Rtensor,1)
             IF (Do_SK) THEN
              ! k^T R k
              R = k_vector(1,iVec)**2*Rtensor(1)+
     &            k_vector(2,iVec)**2*Rtensor(4)+
     &            k_vector(3,iVec)**2*Rtensor(6)+
     &            2.0D0*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+
     &            2.0D0*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+
     &            2.0D0*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
             ELSE
                WRITE(6,43) 'tensor: ',Rtensor(:)
             END IF
            END IF
*
            IF (ABS(R).GT.RSTHR) THEN
              WRITE(6,33) ISS,JSS,R
            END IF
!
            Call Add_Info('CD_M(SO)',[R],1,6)
           END IF
          END DO
         END DO
         WRITE(6,35)
         End Do

         Call Deallocate_Spin_Magnetic_Dipoles()

         Call Deallocate_Electric_Quadrupoles()

         Call CollapseOutput(0,
     &                  'Circular Dichroism - mixed gauge '//
     &                  'Electric-Dipole - Magnetic-Dipole '//
     &                  'rotatory strengths (SO states):')
        END IF

        Call Deallocate_electric_dipoles()
        Call Deallocate_Magnetic_Dipoles()
      END IF
* CD end

! +++ J. Norell 19/7 - 2018
! Dyson amplitudes for (1-electron) ionization transitions
      IF (DYSO) THEN
        Call Add_Info('SODYSAMPS',SODYSAMPS,NSS*NSS,4)
        DYSTHR=1.0D-5
        WRITE(6,*)
        CALL CollapseOutput(1,'Dyson amplitudes '//
     &                        '(SO states):')
        WRITE(6,'(3X,A)')     '----------------------------'//
     &                        '-------------------'
        IF (DYSTHR.GT.0.0D0) THEN
           WRITE(6,*) 'for Dyson intensities at least',DYSTHR
           WRITE(6,*)
        END IF
        WRITE(6,*) '       From      To        '//
     &   'BE (eV)       Dyson intensity'
              WRITE(6,'(3X,A)')
     &                  '---------------------------' //
     &                  '-------------------------------------------'//
     &                  '-------------------'
        DO I=1,NSS
         DO J=1,NSS
          F=SODYSAMPS(I,J)*SODYSAMPS(I,J)
          EDIFF=auToeV*(ENSOR(J)-ENSOR(I))
          IF (ABS(EDIFF)<1.0D-8) Cycle
          IF (F.GT.0.00001) THEN
           IF (EDIFF.GT.0.0D0) THEN
            WRITE(6,'(A,I8,I8,F15.3,ES22.5)') '    ',I,J,EDIFF,F
           END IF
          END IF
         END DO ! J
        END DO ! I
        WRITE(6,*)
        WRITE(6,*)
        CALL CollapseOutput(0,'Dyson amplitudes '//
     &                        '(SO states):')
        WRITE(6,*)
! VKochetov 2021 put SO Dyson amplitudes to hdf5
#ifdef _HDF5_
        if (rhodyn) then
            call mh5_put_dset(wfn_sos_dys, SODYSAMPS)
        endif
#endif
       END IF
! +++ J. Norell

*
************************************************************************
*                                                                      *
*     Start of section for transition moments using the exact operator *
*     for the vector potential.                                        *
*                                                                      *
************************************************************************
*
      If (Do_TMOM)
     &   Call PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,JBNUM,EigVec)
*
 500  CONTINUE


C CALCULATION OF THE D-TENSOR (experimental)
C IFDCAL to implement keyword that will activate computation
C of d-tensor
*     IF(.NOT.IFDCAL) GOTO 600
      GOTO 600

      WRITE(6,*)
      WRITE(6,*) '  D-Matrix'
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  calculated using 2nd order perturbation'
      WRITE(6,*) '  > any spin degeneracy, no spatial degeneracy'
      WRITE(6,*) '  > weak spin-orbit coupling'
      WRITE(6,*)

* SVC 2006: Compute D-tensor through second order perturbation theory.
* no orbitally-degenerate groundstates!
      IAMFIX=0
      IAMFIY=0
      IAMFIZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:4).EQ.'AMFI') THEN
         IF(ICOMP(IPROP).EQ.1) IAMFIX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMFIY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMFIZ=IPROP
       END IF
      END DO
      IPAMFI(1)=IAMFIX
      IPAMFI(2)=IAMFIY
      IPAMFI(3)=IAMFIZ
* initialisations
      DO IXYZ=1,3
       DO JXYZ=1,3
        DTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

* loop over all excited states, different factors will arise depending
* on the difference in spin between ground and excited states.
      ISTATE=1
      MPLET1=MLTPLT(JBNUM(ISTATE))
      S1=0.5D0*DBLE(MPLET1-1)
      FACT0=THREEJ(S1,1.0D0,S1,S1,0.0D0,-S1)*
     &       THREEJ(S1,1.0D0,S1,S1,0.0D0,-S1)/
     &       (S1*S1)
      FACTP=THREEJ(S1+1.0D0,1.0D0,S1,S1+1.0D0,-1.0D0,-S1)*
     &       THREEJ(S1,1.0D0,S1+1.0D0,S1,1.0D0,-(S1+1.0D0))/
     &       ((S1+1.0D0)*(2.0D0*S1+1.0D0))
      FACTM=THREEJ(S1-1.0D0,1.0D0,S1,S1-1.0D0,1.0D0,-S1)*
     &       THREEJ(S1,1.0D0,S1-1.0D0,S1,-1.0D0,-(S1-1.0D0))/
     &       (S1*(2.0D0*S1-1.0D0))
C     WRITE(6,*)
C     WRITE(6,*)'S1 ', S1
C     WRITE(6,*)'FACT0 ', FACT0
C     WRITE(6,*)'FACTP ', FACTP
C     WRITE(6,*)'FACTM ', FACTM
C     WRITE(6,*)
      DO IXYZ=1,3
       DO JXYZ=1,3
        DTIJ=0.0D0
         DO JSTATE=2,NSTATE
          MPLET2=MLTPLT(JBNUM(JSTATE))
          S2=0.5D0*DBLE(MPLET2-1)
          DELTA=ENERGY(JSTATE)-ENERGY(ISTATE)
          IF(DELTA.LT.1.0D-05) GOTO 600
          CONTRIB=0.0D0
          CONTRIB=PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*
     &         PROP(JSTATE,ISTATE,IPAMFI(JXYZ))/DELTA
C         WRITE(6,*) 'ISTATE, JSTATE, IXYZ, JXYZ, DELTA ',
C    &                 ISTATE, JSTATE, IXYZ, JXYZ, DELTA
C         WRITE(6,*) 'CONTRIB ', CONTRIB
C         WRITE(6,*) 'PROP(ISTATE,JSTATE,IXYZ) ',
C    &                 PROP(ISTATE,JSTATE,IPAMFI(IXYZ))
C         WRITE(6,*) 'PROP(JSTATE,ISTATE,JXYZ) ',
C    &                 PROP(JSTATE,ISTATE,IPAMFI(JXYZ))
          IF(S2.EQ.S1) THEN
           DTIJ=DTIJ+FACT0*CONTRIB
          ELSE IF(S2.EQ.S1+1.0D0) THEN
           DTIJ=DTIJ+FACTP*CONTRIB
          ELSE IF(S2.EQ.S1-1.0D0) THEN
           DTIJ=DTIJ+FACTM*CONTRIB
          END IF
         END DO
         DTENS(IXYZ,JXYZ)=DTIJ
       END DO
      END DO

* diagonalisation of the D-tensor matrix
      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
        TMPMAT(IXYZ,JXYZ)=DTENS(IXYZ,JXYZ)
        IF(IXYZ.EQ.JXYZ) THEN
         TMPVEC(IXYZ,JXYZ)=1.0D0
        ELSE
         TMPVEC(IXYZ,JXYZ)=0.0D0
        END IF
       END DO
      END DO
      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

* D-factor printout
* D = D_zz - 1/2 * (D_xx + D_yy)
* E = 1/2 * (D_xx - D_yy)
      WRITE(6,*)
      WRITE(6,*) 'The D matrix and eigenvalues:'
      WRITE(6,*)
      WRITE(6,'(2x,2x,2x,3(5x,a2,5x),'//
     & '4x,4x,2x,10x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3),
     & ('D_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.8,1x),'//
     & '4x,a2,i1,a1,2x,f10.8,'//
     & '2x,a2,2x,3(1x,f8.4,1x),'//
     & '3x,a2,i1,a1,2x,f8.3,2x,a5)')
     & xyzchr(IXYZ), (DTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'D_',IXYZ,':',EVR(IXYZ),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3),
     & 'D_',IXYZ,':',EVR(IXYZ)*auTocm,'cm^-1'
      ENDDO


 600  CONTINUE

C CALCULATION OF THE G-TENSOR
C IFGCAL is set by keyword EPRG
      IF(.NOT.IFGCAL) GOTO 800
* PAM 2005 Experimental: Compute g-tensor through 2-nd order
* perturbation approach, mixed ang mom / spin orbit terms
* Declarations for gtensor(3,3) and some other odds and ends
* for nice output have been added to declaration head above.

      IF(.not.IFSO) THEN
          WRITE(6,*) 'keyword SPIN needed together with EPRG'
          WRITE(6,*)
          GOTO 800
      ENDIF

      WRITE(6,*)
      WRITE(6,*) '  g-Matrix Approach I                               '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  calculated using 2nd order perturbation'
      WRITE(6,*) '  > any spin degeneracy, no spatial degeneracy'
      WRITE(6,*) '  > weak spin-orbit coupling'
      WRITE(6,*)

      IAMFIX=0
      IAMFIY=0
      IAMFIZ=0
      IAMX=0
      IAMY=0
      IAMZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:4).EQ.'AMFI') THEN
         IF(ICOMP(IPROP).EQ.1) IAMFIX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMFIY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMFIZ=IPROP
       ELSE IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO
      IPAMFI(1)=IAMFIX
      IPAMFI(2)=IAMFIY
      IPAMFI(3)=IAMFIZ
      IPAM(1)=IAMX
      IPAM(2)=IAMY
      IPAM(3)=IAMZ

C start loop over the states ISTATE:
      ISTATE=1
      DO WHILE (
     &          ((ENERGY(MIN(ISTATE,NSTATE))-ENERGY(1)).LE.EPRTHR)
     &  .AND. (ISTATE.LE.NSTATE)
     &          )

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      MPLET=MLTPLT(JBNUM(ISTATE))
      S=0.5D0*DBLE(MPLET-1)

      WRITE(6,*)
      WRITE(6,'(3x,A6,I4,3x,A4,F4.1,3x,A4,F18.8)')
     & 'STATE ',ISTATE,'S = ',S,'E = ',ENERGY(ISTATE)
      WRITE(6,'(3x,A46)')
     & '----------------------------------------------'

      IF(MPLET.NE.1) THEN
              FACTOR=1.0D0/SQRT(S*(S+1.0D0)*(2.0D0*S+1.0D0))
      ELSE
              GOTO 690
      END IF

C print separate contributions if verbose
      IF (IPGLOB.GE.3) THEN
       WRITE(6,*)
       WRITE(6,*) 'contributions from the SOS expansion'//
     &            ' to delta(g_pq) in *ppt* (p,q=x,y,z)'
       WRITE(6,*)
       WRITE(6,'(2x,a8,2x,9(4x,a2,3x))')
     &  ' states ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO JSTATE=1,NSTATE
       IF (JSTATE.NE.ISTATE) THEN
       DELTA=ENERGY(JSTATE)-ENERGY(ISTATE)
        IF(ABS(DELTA).LT.1.0D-04) THEN
            WRITE(6,'(1x,i3,2x,i3,3x,A20,1x,A20,F18.8)')
     &       ISTATE,JSTATE,'possible degeneracy,',
     &       'energy difference = ',DELTA
             GOTO 610
        ENDIF
        DO IXYZ=1,3
         DO JXYZ=1,3
          CONTRIB=PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*
     &             PROP(ISTATE,JSTATE,IPAM(JXYZ))
          CONTRIB=CONTRIB/DELTA
          SOSTERM(3*(IXYZ-1)+JXYZ)=-2.0D0*FACTOR*CONTRIB
         ENDDO
        ENDDO
        WRITE(6,'(1x,i3,2x,i3,3x,9(f8.3,1x))')
     &   ISTATE,JSTATE,(SOSTERM(I)*1.0D3,I=1,9)
       ENDIF
 610   CONTINUE
       ENDDO
      END IF

C calculate sum-over-states for each g_pq (p,q = x,y,z)
      DO IXYZ=1,3
       DO JXYZ=1,3
         GTIJ=0.0D0
         DELTA=0.0D0
         CONTRIB=0.0D0
         DO JSTATE=1,NSTATE
         IF (JSTATE.NE.ISTATE) THEN
          DELTA=ENERGY(JSTATE)-ENERGY(ISTATE)
C SVC 2008: no good criterium for spatial degeneracy, use rasscf
C energies ?
          IF(ABS(DELTA).LT.1.0D-04) THEN
              WRITE(6,*)
              WRITE(6,*) 'SPATIALLY DEGENERATE STATE: '//
     &                   'sum-over-states not applicable'
*             WRITE(6,*) '> lower the degeneracy treshold if this '//
*    &                   'is not a spatially degenerate state'
              WRITE(6,*)
              GOTO 690
          ENDIF
          CONTRIB=PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*
     &             PROP(ISTATE,JSTATE,IPAM(JXYZ))
          CONTRIB=CONTRIB/DELTA
          GTIJ=GTIJ+CONTRIB
         ENDIF
         ENDDO
        GTENS(IXYZ,JXYZ)=-2.0D0*FACTOR*GTIJ
       END DO
      END DO

C put g_e on the diagonal
      DO IXYZ=1,3
       DO JXYZ=IXYZ,3
       IF (IXYZ.EQ.JXYZ) GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+FEGVAL
       END DO
      END DO

C determine symmetric G = gg+ tensor, this is what can be measured
C experimentally, and store as GSTENS
      DO IXYZ=1,3
       DO JXYZ=1,3
        GSTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        GSTENS(IXYZ,JXYZ)=GSTENS(IXYZ,JXYZ)+
     &                     GTENS(IXYZ,KXYZ)*GTENS(JXYZ,KXYZ)
        END DO
       END DO
      END DO

C determine the eigenvalues of the g matrix
      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
C XEIGEN alters the input matrix! copy GTENS to TMPMAT
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      IERR=0
      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)
      IF (IERR.NE.0) THEN
          WRITE(6,*) 'Error: xEigen returned IERR = ', IERR
          RETURN
      END IF

      WRITE(6,*)
      WRITE(6,*) 'The g matrix and eigenvalues:'
      WRITE(6,*)
      WRITE(6,'(6x,3(5x,a2,5x))')
     & (xyzchr(IXYZ),IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.8,1x),'//
     & '4x,a2,i1,a1,2x,f8.4,3x,a8,i1,a2,2x,f10.3,2x,a3)')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'g_',IXYZ,':',EVR(IXYZ),
     & 'delta(g_',IXYZ,'):',(EVR(IXYZ)-FEGVAL)*1.0D3,'ppt'
      ENDDO

*     WRITE(6,'(6x,3(5x,i1,4x))') (IXYZ,IXYZ=1,3)

C determine the eigenvalues of the G = gg* matrix
      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      ENDDO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GSTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      IERR=0
      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)
      IF (IERR.NE.0) THEN
          WRITE(6,*) 'Error: xEigen returned IERR = ', IERR
          RETURN
      END IF

C reconstruct g_s from the square root of the eigenvalues
C and the eigenvectors of G = gg+ by back transformation
      DO IXYZ=1,3
       DO JXYZ=1,3
       GSTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        GSTENS(IXYZ,JXYZ)=GSTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,*)
      WRITE(6,*) 'The symmetric g matrix is the '//
     &           'actual experimentally determined g matrix.'
      WRITE(6,*) 'The sign of the eigenvalues is undetermined '//
     &           '(assumed positive).'
      WRITE(6,*)
      WRITE(6,'(2x,2x,2x,3(5x,a2,5x),'//
     & '4x,4x,2x,10x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3),
     & ('g_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.8,1x),'//
     & '4x,a2,i1,a1,2x,f10.6,'//
     & '2x,a2,2x,3(1x,f8.4,1x),'//
     & '3x,a8,i1,a2,2x,f8.3,2x,a3)')
     & xyzchr(IXYZ), (GSTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'g_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3),
     & 'delta(g_',IXYZ,'):',(SQRT(EVR(IXYZ))-FEGVAL)*1.0D3,'ppt'
      ENDDO
      DO I=1,3
      EVR(I)=SQRT(EVR(I))-FEGVAL
      ENDDO
      Call Add_Info('EPRGVAL',EVR,3,6)

 690  CONTINUE

      ISTATE=ISTATE+1

* end long loop over states ISTATE
      ENDDO

* SVC alternative approach for the g-tensor:
* using first order degenerate perturbation theory

      IF(IFVANVLECK) THEN
      WRITE(6,*)
      WRITE(6,*) '  VAN VLECK Tensor and g-Matrix Approach II        '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ELSE
      WRITE(6,*)
      WRITE(6,*) '  g-Matrix Approach II                     '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ENDIF

      IAMX=0
      IAMY=0
      IAMZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
        !write(6,*)"3****ANGMOM rassi/prprop.f "
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO

      Call mma_allocate(LXI,NSS,NSS,Label='LXI')
      LXI(:,:)=0.0D0
      Call mma_allocate(LYI,NSS,NSS,Label='LYI')
      LYI(:,:)=0.0D0
      Call mma_allocate(LZI,NSS,NSS,Label='LZI')
      LZI(:,:)=0.0D0

      IF(IAMX.GT.0) CALL SMMAT(PROP,LXI,NSS,IAMX,0)
      IF(IAMY.GT.0) CALL SMMAT(PROP,LYI,NSS,IAMY,0)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,LZI,NSS,IAMZ,0)


      Call mma_allocate(ZXR,NSS,NSS,Label='ZXR')
      Call mma_allocate(ZXI,NSS,NSS,Label='ZXI')
      ZXR(:,:)=0.0D0
      ZXI(:,:)=0.0D0
      pZMR(1)%A2=>ZXR(:,:)
      pZMI(1)%A2=>ZXI(:,:)
      Call mma_allocate(ZYR,NSS,NSS,Label='ZYR')
      Call mma_allocate(ZYI,NSS,NSS,Label='ZYI')
      ZYR(:,:)=0.0D0
      ZYI(:,:)=0.0D0
      pZMR(2)%A2=>ZYR(:,:)
      pZMI(2)%A2=>ZYI(:,:)
      Call mma_allocate(ZZR,NSS,NSS,Label='ZZR')
      Call mma_allocate(ZZI,NSS,NSS,Label='ZZI')
      ZZR(:,:)=0.0D0
      ZZI(:,:)=0.0D0
      pZMR(3)%A2=>ZZR(:,:)
      pZMI(3)%A2=>ZZI(:,:)

      CALL SMMAT(PROP,ZXR,NSS,0,1)
      CALL SMMAT(PROP,ZYI,NSS,0,2)
      CALL SMMAT(PROP,ZZR,NSS,0,3)

      CALL DSCAL_(NSS**2,FEGVAL,ZXR,1)
      CALL DSCAL_(NSS**2,FEGVAL,ZYI,1)
      CALL DSCAL_(NSS**2,FEGVAL,ZZR,1)

      CALL DAXPY_(NSS**2,1.0D0,LXI,1,ZXI,1)
      CALL DAXPY_(NSS**2,1.0D0,LYI,1,ZYI,1)
      CALL DAXPY_(NSS**2,1.0D0,LZI,1,ZZI,1)

      Call mma_deallocate(LXI)
      Call mma_deallocate(LYI)
      Call mma_deallocate(LZI)

*     SVC 20090926 Experimental
*     Add analysis of different contributions

*     Establish which spin components of SFS belong to the ground state
      DO I=1,NSS
       ISGS(I)=.FALSE.
      ENDDO

      GSENERGY=ENERGY(1)
      DO ISTATE=2,NSTATE
       IF (ENERGY(ISTATE).LT.GSENERGY) THEN
        GSENERGY = ENERGY(ISTATE)
       ENDIF
      ENDDO

      IMLTPL=1
      DO ISTATE=1,NSTATE
       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.TRUE.
        ENDDO
       ELSE
        DO I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
         ISGS(I)=.FALSE.
        ENDDO
       ENDIF
       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

*     Analyze the different contributions to the GS Kramers doublet
*     Zeeman matrix elements.  There are 4 different ME's: <1|Ze|1>,
*     <1|Ze|2>, <2|Ze|1>, and <2|Ze|2>, stored in ZEKL.  Contributions
*     of SFS i,j to SOS k,l (k,l=1,2): <k|Ze|l> = Sum(i,j) U(i,k)*
*     <i|Ze|j> U(j,l).  This sum is decomposed into parts belonging to
*     each SFS state i as follows:
*     -> GS's contain only MEs with themselves and other GS's
*     -> ES's contain MEs with themselves, the GS's (2x) and other ES's
*        The ME's with the GS's are counted twice as they do not belong to
*        any GS's (they contain only ME's within their own GS group)
*        The contributions with other ES's are split between the ES's,
*        counting them double (<i|Ze|j> and <j|Ze|i>) and divide by two later.

      IMLTPL=1
      ZEKL(:,:,:,:)=CMPLX(0.0d0,0.0d0,kind=8)

      DO ISTATE=1,NSTATE

       ISTART=IMLTPL
       IFINAL=IMLTPL-1+MLTPLT(JBNUM(ISTATE))

       IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contribution of the GS spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                    ZEKL(:,:,IXYZ,ISTATE)
           ENDIF
          ENDDO
         ENDDO
        ENDDO

       ELSE

* Contributions of the ES spin components
        DO IXYZ=1,3
         DO ISS=ISTART,IFINAL
          DO JSS=1,NSS
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           pZMR(IXYZ)%A2,pZMI(IXYZ)%A2,
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

C       DO IXYZ=1,3
C        WRITE(6,FMT=720) 'ZEKL', IXYZ, ISTATE,
C     &      ZEKL(:,:,IXYZ,ISTATE)
C       ENDDO

C 710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
C 720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO

*     We now have decomposed the <k|Ze|l> into terms belonging to either
*     a GS or an ES for each k,l=1,2 and p=x,y,z stored in ZEKL(k,l,p,SFS)
*     Now, these new decomposed terms of the ZEKL ME's are combined to
*     form the G tensor.  Consider e.g. that <k|Ze|l> is decomposed into
*     <k|GS|l> + <k|ES1|l> + <k|ES2|l>, then the contributions to G are given as:
*     -> G_pq/2 = <k|Ze_p|l> <l|Ze_q|k>
*             = (<k|GS_p|l> + <k|ES1_p|l> + <k|ES2_p|l>)
*             * (<l|GS_q|k> + <l|ES1_q|k> + <l|ES2_q|k>)
*
*     from GS: (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2
*     from ES1: 2*((<k|ES1_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_q|l>/2 * <l|ES1_p|k>/2)/2)
*               + (<k|ES1_p|l>/2 * <l|ES1_q|k>/2)/2 + (<k|ES1_q|l>/2 * <l|ES2_p|k>/2)/2
*               + (<k|ES1_p|l>/2 * <l|ES2_q|k>/2)/2 + (<k|ES2_q|l>/2 * <l|ES1_p|k>/2)/2
*     In the end, the outer division by 2 cancels on both sides, and the
*     inner divisions by two combine to a division by 4.

      GCONT(:,:)=CMPLX(0.0d0,0.0d0,kind=8)
      GTOTAL(:)=0.0d0

      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO JXYZ=1,3
         IJXYZ=3*(IXYZ-1)+JXYZ

         IF (ABS(ENERGY(ISTATE)-GSENERGY).LT.1.0d-6) THEN

* Contributions for the GS's
          DO JSTATE=1,NSTATE
           IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
            DO I=1,2
             DO J=1,2
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDDO
            ENDDO
           ENDIF
          ENDDO

         ELSE

* Contributions for the ES's
          DO JSTATE=1,NSTATE
           DO I=1,2
            DO J=1,2
             GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &            +(ZEKL(I,J,IXYZ,ISTATE)*
     &            ZEKL(J,I,JXYZ,JSTATE))/4
     &            +(ZEKL(I,J,IXYZ,JSTATE)*
     &            ZEKL(J,I,JXYZ,ISTATE))/4
             IF (ABS(ENERGY(JSTATE)-GSENERGY).LT.1.0d-6) THEN
              GCONT(IJXYZ,ISTATE)=GCONT(IJXYZ,ISTATE)
     &             +(ZEKL(I,J,IXYZ,ISTATE)*
     &             ZEKL(J,I,JXYZ,JSTATE))/4
     &             +(ZEKL(I,J,IXYZ,JSTATE)*
     &             ZEKL(J,I,JXYZ,ISTATE))/4
             ENDIF
            ENDDO
           ENDDO
          ENDDO

         ENDIF

        ENDDO
       ENDDO

       DO IJXYZ=1,9
        GTOTAL(IJXYZ)=GTOTAL(IJXYZ)+DBLE(GCONT(IJXYZ,ISTATE))
       ENDDO
      ENDDO

       DIPSOm(:,:,:)=(0.0d0,0.0d0)
       DIPSOn(:,:,:)=(0.0d0,0.0d0)


*     Continue original calculation of G tensor (=gg^*)
      CALL get_dArray( 'ESO_SINGLE',ESO,NSS)
      CALL ZTRNSF(NSS,USOR,USOI,ZXR,ZXI)
      CALL MULMAT(NSS,ZXR,ZXI,eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOm(1,ISS,JSS)=0.5d0*Z(ISS,JSS)
      DIPSOn(1,ISS,JSS)=-Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,ZYR,ZYI)
      CALL MULMAT(NSS,ZYR,ZYI,eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOm(2,ISS,JSS)=0.5d0*Z(ISS,JSS)
      DIPSOn(2,ISS,JSS)=-Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,ZZR,ZZI)
      CALL MULMAT(NSS,ZZR,ZZI,eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOm(3,ISS,JSS)=0.5d0*Z(ISS,JSS)
      DIPSOn(3,ISS,JSS)=-Z(ISS,JSS)
      enddo
      enddo
      WRITE(6,*)''

      IF(IFVANVLECK) THEN

      iT=0
      do iT=1,NTS
      do ic=1,3
      do jc=1,3
      chiT_tens(iT,ic,jc)=0.d0
      chicuriT_tens(iT,ic,jc)=0.d0
      chiparamT_tens(iT,ic,jc)=0.d0
      enddo
      enddo
      enddo
      iT=0
      do iT=1,NTS
      if(iT.eq.1) then
      TMPm(iT)=TMINS+0.0001d0
      ELSE
      DLTT=(TMAXS-TMINS)/(dble(NTS-1))
      TMPm(iT)=TMINS+DLTT*dble(iT-1)
      ENDIF
      Zstat=0.d0
      do Iss=1,Nss
      p_Boltz=EXP(-ESO(Iss)/Boltz_k/TMPm(iT))
      Zstat=Zstat+p_Boltz
      do IC=1,3
      do JC=1,3
      c_2(IC,JC)   =0.d0
      curit(IC,JC) =0.d0
      paramt(IC,JC)=0.d0
      enddo
      enddo
      do Jss=1,Nss
      dlt_E=Eso(Iss)-Eso(Jss)
      do IC=1,3
      do JC=1,3
      c_1(IC,JC)=0.d0
      enddo
      enddo
      do ic=1,3
      do jc=1,3
      c_1(ic,jc)=DBLE(DIPSOn(ic,Iss,Jss)*CONJG(DIPSOn(jc,Iss,Jss)))
      if(ABS(dlt_E).LT.10.97d0) then
      c_2(ic,jc)=    c_2(ic,jc)  +  c_1(ic,jc)
      curit(ic,jc)= curit(ic,jc) +  c_1(ic,jc)
      paramt(ic,jc)=paramt(ic,jc)+  0.d0*c_1(ic,jc)
      else
      c_2(ic,jc)= c_2(ic,jc)-2.d0*Boltz_k*TMPm(iT)* c_1(ic,jc)/dlt_E
      curit(ic,jc)= curit(ic,jc)-0.d0*(2.d0*Boltz_k*TMPm(iT)*
     &c_1(ic,jc)/dlt_E)
      paramt(ic,jc)= paramt(ic,jc)-2.d0*Boltz_k*TMPm(iT)*
     &c_1(ic,jc)/dlt_E
      endif
      enddo
      enddo
      enddo !Jss
      do ic=1,3
      do jc=1,3
      chiT_tens(iT,ic,jc)=    chiT_tens(iT,ic,jc)+p_Boltz*
     & c_2(ic,jc)
      chicuriT_tens(iT,ic,jc)= chicuriT_tens(iT,ic,jc) +
     & p_Boltz*curit(ic,jc)
      chiparamT_tens(iT,ic,jc)= chiparamT_tens(iT,ic,jc) +
     & p_Boltz*paramt(ic,jc)
      enddo
      enddo
      enddo !Iss
      !Zstat1m(iT)=Zstat
      do ic=1,3
      do jc=1,3
      chiT_tens(iT,ic,jc)= coeff_chi*(chiT_tens(iT,ic,jc)/Zstat)
      chicuriT_tens(iT,ic,jc)=coeff_chi*(chicuriT_tens(iT,ic,jc)/
     & Zstat)
      chiparamT_tens(iT,ic,jc)=coeff_chi*(chiparamT_tens(iT,ic,jc)/
     & Zstat)
      enddo
      enddo
      enddo ! iT


      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR'//
     & '  (cm3*K/mol)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTS
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),
     & ((chiT_tens(iT,ic,jc),jc=1,3),ic=1,3)
      enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Curie contrib. to VAN VLECK TENSOR'//
     & '  (cm3*K/mol)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTS
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),
     & ((chicuriT_tens(iT,ic,jc),jc=1,3),ic=1,3)
      enddo
      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Parama. contrib. to VAN VLECK TENSOR'//
     & '  (cm3*K/mol)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTS
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),
     & ((chiparamT_tens(iT,ic,jc),jc=1,3),ic=1,3)
      enddo
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '  g-Matrix'
      WRITE(6,*) '  =========='
      ENDIF ! IFVANVLECK
      !do I=1,3
      !do J=1,3
      !do iT=1,NT
      !chiT_tens(iT,I,J)=0.d0
      !enddo
      !enddo
      !enddo

      ISS=1
      DO WHILE ((ISS.LE.NSS).AND.
     &          (ENSOR(MIN(ISS,NSS))-ENSOR(1).LE.EPRTHR))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFGTCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
      !WRITE(6,*) 'KDGN=',KDGN
      ELSE IF (ABS(EDIFF).LT.1.0D-06) THEN
      KDGN=KDGN+1
      ENDIF
      ENDDO

      WRITE(6,*)
      DO I=1,KDGN
      WRITE(6,'(3x,A9,I4,3x,A4,F18.8)')
     & 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      ENDDO
      WRITE(6,'(3x,A46)')
     &'----------------------------------------------'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(.NOT.IFGTCALSA) GOTO 450
      IF(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOn,SPNSFS,DIPSOm_SA)
      IFUNCT=IFUNCT+KDGN
  450 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.3)) THEN
       WRITE(6,*) 'Experimental: SFS contributions to G=gg+'
       WRITE(6,*)
       WRITE(6,'(a6,9(5x,a2,5x))')
     &      'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
       WRITE(6,*)
       DO ISTATE=1,NSTATE
        WRITE(6,'(2x,I2,2x,9(F12.6))')
     &       ISTATE, (DBLE(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
       ENDDO

       WRITE(6,*)
       WRITE(6,'(A6,9(F12.6))')
     &      'total ', (GTOTAL(IJXYZ),IJXYZ=1,9)
      ENDIF


      JSS=ISS+1

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTIJ=0.0D0
       CONTRIB=0.0D0
       DO ISO=ISS,JSS
        DO JSO=ISS,JSS
        CONTRIB= pZMR(IXYZ)%A2(ISO,JSO)*pZMR(JXYZ)%A2(JSO,ISO)
     &          -pZMI(IXYZ)%A2(ISO,JSO)*pZMI(JXYZ)%A2(JSO,ISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      IF(IPGLOB.GT.3) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
       END IF

      DO I=1,3
      EVR(I)=0.0D0
      EVI(I)=0.0D0
      END DO
      DO IXYZ=1,3
       DO JXYZ=1,3
       TMPMAT(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)
       IF(IXYZ.EQ.JXYZ) THEN
        TMPVEC(IXYZ,JXYZ)=1.0D0
       ELSE
        TMPVEC(IXYZ,JXYZ)=0.0D0
       END IF
       END DO
      ENDDO

      CALL XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

C construct g_s matrix from G by back-transormation of the
C square root of the G eigenvalues
      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('g_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f8.4,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'g_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

 780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      Call mma_deallocate(ZXR)
      Call mma_deallocate(ZXI)
      Call mma_deallocate(ZYR)
      Call mma_deallocate(ZYI)
      Call mma_deallocate(ZZR)
      Call mma_deallocate(ZZI)
      nullify(pZMR(1)%A2,pZMI(1)%A2,pZMR(2)%A2,pZMI(2)%A2,pZMR(3)%A2,
     &        pZMI(3)%A2)

 800  CONTINUE


******************************************************
******************************************************
******************************************************
** Experimental hyperfine tensor stuff starts here
******************************************************
******************************************************
******************************************************

* Skip if not a hyperfine calculation
      IF(.NOT.IFACAL) GOTO 1801
        CALL HFCTS(PROP,USOR,USOI,ENSOR,NSS,ENERGY,JBNUM,
     &             DIPSOM,ESO,XYZCHR,BOLTZ_K)

1801  CONTINUE
******************************************************
******************************************************
******************************************************
** Experimental hyperfine tensor stuff ends here
******************************************************
******************************************************
******************************************************



C SVC20080312 calculation of magnetization

      IF(.NOT.IFXCAL) GOTO 900

      IF(.not.IFSO) THEN
          WRITE(6,*) 'keyword SPIN needed with MAGN'
          WRITE(6,*)
          GOTO 900
      ENDIF

      WRITE(6,*)
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  Magnetization and Magnetic Susceptibility'
      WRITE(6,*) '  ========================================='
      WRITE(6,*)

C initialization same as G-tensor, construct L+gS matrix elements
      IAMX=0
      IAMY=0
      IAMZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
        !write(6,*)"4****ANGMOM rassi/prprop.f "
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO

      Call mma_allocate(LXI,NSS,NSS,Label='LXI')
      LXI(:,:)=0.0D0
      Call mma_allocate(LYI,NSS,NSS,Label='LYI')
      LYI(:,:)=0.0D0
      Call mma_allocate(LZI,NSS,NSS,Label='LZI')
      LZI(:,:)=0.0D0

      IF(IAMX.GT.0) CALL SMMAT(PROP,LXI,NSS,IAMX,0)
      IF(IAMY.GT.0) CALL SMMAT(PROP,LYI,NSS,IAMY,0)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,LZI,NSS,IAMZ,0)

      Call mma_allocate(MXR,NSS,NSS,Label='MXR')
      Call mma_allocate(MXI,NSS,NSS,Label='MXI')
      MXR(:,:)=0.0D0
      MXI(:,:)=0.0D0
      Call mma_allocate(MYR,NSS,NSS,Label='MYR')
      Call mma_allocate(MYI,NSS,NSS,Label='MYI')
      MYR(:,:)=0.0D0
      MYI(:,:)=0.0D0
      Call mma_allocate(MZR,NSS,NSS,Label='MZR')
      Call mma_allocate(MZI,NSS,NSS,Label='MZI')
      MZR(:,:)=0.0D0
      MZI(:,:)=0.0D0

      pMR(1)%A2=>MXR(:,:)
      pMI(1)%A2=>MXI(:,:)
      pMR(2)%A2=>MYR(:,:)
      pMI(2)%A2=>MYI(:,:)
      pMR(3)%A2=>MZR(:,:)
      pMI(3)%A2=>MZI(:,:)

      CALL SMMAT(PROP,MXR,NSS,0,1)
      CALL SMMAT(PROP,MYI,NSS,0,2)
      CALL SMMAT(PROP,MZR,NSS,0,3)

      CALL DSCAL_(NSS**2,FEGVAL,MXR,1)
      CALL DSCAL_(NSS**2,FEGVAL,MYI,1)
      CALL DSCAL_(NSS**2,FEGVAL,MZR,1)

      CALL DAXPY_(NSS**2,1.0D0,LXI,1,MXI,1)
      CALL DAXPY_(NSS**2,1.0D0,LYI,1,MYI,1)
      CALL DAXPY_(NSS**2,1.0D0,LZI,1,MZI,1)

      Call mma_deallocate(LXI)
      Call mma_deallocate(LYI)
      Call mma_deallocate(LZI)

      CALL ZTRNSF(NSS,USOR,USOI,MXR,MXI)
      CALL ZTRNSF(NSS,USOR,USOI,MYR,MYI)
      CALL ZTRNSF(NSS,USOR,USOI,MZR,MZI)

      Call mma_allocate(ZXR,NSS,NSS,Label='ZXR')
      Call mma_allocate(ZXI,NSS,NSS,Label='ZXI')
      ZXR(:,:)=0.0D0
      ZXR(:,:)=0.0D0
      pZMR(1)%A2=>ZXR(:,:)
      pZMI(1)%A2=>ZXI(:,:)
      Call mma_allocate(ZYR,NSS,NSS,Label='ZYR')
      Call mma_allocate(ZYI,NSS,NSS,Label='ZYI')
      ZYR(:,:)=0.0D0
      ZYR(:,:)=0.0D0
      pZMR(2)%A2=>ZYR(:,:)
      pZMI(2)%A2=>ZYI(:,:)
      Call mma_allocate(ZZR,NSS,NSS,Label='ZZR')
      Call mma_allocate(ZZI,NSS,NSS,Label='ZZI')
      ZZR(:,:)=0.0D0
      ZZR(:,:)=0.0D0
      pZMR(3)%A2=>ZZR(:,:)
      pZMI(3)%A2=>ZZI(:,:)

      Call mma_allocate(ZR,NSS,NSS,Label='ZR')
      Call mma_allocate(ZI,NSS,NSS,Label='ZI')
      CALL mma_allocate(UZR,NSS,NSS,Label='UZR')
      CALL mma_allocate(UZI,NSS,NSS,Label='UZI')

      BFINAL=BSTART+(NBSTEP-1)*BINCRE
      TFINAL=TSTART+(NTSTEP-1)*TINCRE

      WRITE(6,*) "Magnetic flux density range (T): "
      WRITE(6,'(2x,f6.2,a3,f6.2,a4,i4,a6)')
     & BSTART," - ",BFINAL," in ",NBSTEP," steps"
      WRITE(6,*)
      WRITE(6,*) "Temperature range (K): "
      WRITE(6,'(2x,f6.2,a3,f6.2,a4,i4,a6)')
     & TSTART," - ",TFINAL," in ",NTSTEP," steps"

      CALL mma_allocate(MAGM,9*NBSTEP*NTSTEP,Label='MAGM')

      LMSTEP=0

      DO IXYZ=1,3

      WRITE(6,*)
      WRITE(6,'(3x,a1,3x,8(1x,a12,1x))') "T",
     & "    B"//xyzchr(IXYZ)//" (T)  ","   M (J/T)  ",
     & "  Mx (J/T)  ","  My (J/T)  ","  Mz (J/T)  ",
     & "Xx"//xyzchr(IXYZ)//" (m3/mol)",
     & "Xy"//xyzchr(IXYZ)//" (m3/mol)",
     & "Xz"//xyzchr(IXYZ)//" (m3/mol)"
      WRITE(6,*)

       DO IBSTEP=1,NBSTEP
        B=BSTART+BINCRE*(IBSTEP-1)
        ZR(:,:)=0.0D0
        ZI(:,:)=0.0D0
        CALL DAXPY_(NSS**2,0.5D0*B/auToT,pMR(IXYZ)%A2,1,ZR,1)
        CALL DAXPY_(NSS**2,0.5D0*B/auToT,pMI(IXYZ)%A2,1,ZI,1)
        DO ISS=1,NSS
         HZER=ZR(ISS,ISS)
         ZR(ISS,ISS)=HZER+ENSOR(ISS)
        END DO
        CALL DCOPY_(NSS**2,[0.0D0],0,UZR,1)
        CALL DCOPY_(NSS   ,[1.0D0],0,UZR,NSS+1)
        CALL DCOPY_(NSS**2,[0.0D0],0,UZI,1)
        CALL ZJAC(NSS,ZR,ZI,NSS,UZR,UZI)
        DO JXYZ=1,3
         CALL DCOPY_(NSS**2,pMR(JXYZ)%A2,1,pZMR(JXYZ)%A2,1)
         CALL DCOPY_(NSS**2,pMI(JXYZ)%A2,1,pZMI(JXYZ)%A2,1)
         CALL DSCAL_(NSS**2,-0.5D0,pZMR(JXYZ)%A2,1)
         CALL DSCAL_(NSS**2,-0.5D0,pZMI(JXYZ)%A2,1)
         CALL ZTRNSF(NSS,UZR,UZI,pZMR(JXYZ)%A2,pZMI(JXYZ)%A2)
        ENDDO
        DO ITSTEP=1,NTSTEP
         T=TSTART+TINCRE*(ITSTEP-1)
         RkT=T*BOLTZ
         RMAGM(1)=0.0D0
         RMAGM(2)=0.0D0
         RMAGM(3)=0.0D0
         RPART=0.0D0
         IF(IPGLOB.GT.2) THEN
          WRITE(6,*)
          WRITE(6,'(2x,a14,3(4x,a4,4x),2x,a6)') "Energy (cm^-1)",
     &     "mu_x", "mu_y", "mu_z","weight"
          WRITE(6,*)
         ENDIF
         DO ISS=1,NSS
          DELTA=ZR(ISS,ISS)-ZR(1,1)
          FACT=EXP(-DELTA/RkT)
          RMAGM(1)=RMAGM(1)+pZMR(1)%A2(ISS,ISS)*FACT
          RMAGM(2)=RMAGM(2)+pZMR(2)%A2(ISS,ISS)*FACT
          RMAGM(3)=RMAGM(3)+pZMR(3)%A2(ISS,ISS)*FACT
          RPART=RPART+FACT
          IF(IPGLOB.GT.2) THEN
           WRITE(6,'(2x,f14.3,3(1x,f10.6,1x),2x,f6.3)')
     &      (ZR(ISS,ISS)-ZR(1,1))*auTocm,
     &      pZMR(1)%A2(ISS,ISS),
     &      pZMR(2)%A2(ISS,ISS),
     &      pZMR(3)%A2(ISS,ISS),FACT
          ENDIF
         ENDDO
         IF(IPGLOB.GT.2) THEN
          WRITE(6,*)
         ENDIF
         RMAGM(1)=(RMAGM(1)/RPART)*AU2JTM
         RMAGM(2)=(RMAGM(2)/RPART)*AU2JTM
         RMAGM(3)=(RMAGM(3)/RPART)*AU2JTM
         RMAGM2=RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
         RMAGMO=SQRT(RMAGM2)
         DO JXYZ=1,3
          LMSTEP=LMSTEP+1
          MAGM(LMSTEP)=RMAGM(JXYZ)
          IF(IBSTEP.GT.1) THEN
              Chi(JXYZ)=RMAGM(JXYZ)-MAGM(LMSTEP-3*NTSTEP)
              Chi(JXYZ)=Chi(JXYZ)*Rmu0/BINCRE
          ENDIF
         ENDDO
          IF(IBSTEP.EQ.1) THEN
         WRITE(6,'(1x,f6.2,5(1x,es12.5,1x))')
     &    T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
          ELSE
         WRITE(6,'(1x,f6.2,8(1x,es12.5,1x))')
     &    T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3),Chi(1),Chi(2),Chi(3)
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      CALL mma_deallocate(MAGM)

      WRITE(6,*)

C powder magnetization, useful in nonlinear cases

      IF(.NOT.IFMCAL) GOTO 810

      WRITE(6,*)
      WRITE(6,*) "Powder Magnetization"
      WRITE(6,*)
      WRITE(6,'(3x,a1,3x,5(1x,a12,1x))') "T",
     & "    B  (T)  ","   M (J/T)  ",
     & "  Mx (J/T)  ","  My (J/T)  ","  Mz (J/T)  "

      CALL mma_allocate(MAGM,3*NBSTEP*NTSTEP,Label='MAGM')
      MAGM(:)=0.0D0

      NPHISTEP=INT(360.0D0/BANGRES)
      NTHESTEP=INT(180.0D0/BANGRES)

      GTR=ACOS(-1.0D0)/180

C scale number of points on phi via sin(theta)
      NORIENT=0
      DO ITHE=1,NTHESTEP+1
       THE=BANGRES*(ITHE-1)*GTR
       IPHISTEP=INT((NPHISTEP-1)*SIN(THE)+1)
       BPHIRES=360/IPHISTEP
       DO IPHI=1,IPHISTEP
       PHI=BPHIRES*(IPHI-1)*GTR

       NORIENT=NORIENT+1

       LMSTEP=0
*      WRITE(6,*)
*      WRITE(6,'(1x,2(A6,I4))') ' ITHE ',ITHE,' IPHI',IPHI
*      WRITE(6,'(6(5x,A4,5x))')
*    &  ' B  ','THE ','PHI ',' Mx ',' My ',' Mz '
       DO IBSTEP=1,NBSTEP
        B=BSTART+BINCRE*(IBSTEP-1)
        BX=B*SIN(THE)*COS(PHI)
        BY=B*SIN(THE)*SIN(PHI)
        BZ=B*COS(THE)
        CALL DCOPY_(NSS**2,[0.0D0],0,ZR,1)
        CALL DAXPY_(NSS**2,0.5D0*BX/auToT,MXR,1,ZR,1)
        CALL DAXPY_(NSS**2,0.5D0*BY/auToT,MYR,1,ZR,1)
        CALL DAXPY_(NSS**2,0.5D0*BZ/auToT,MZR,1,ZR,1)
        CALL DCOPY_(NSS**2,[0.0D0],0,ZI,1)
        CALL DAXPY_(NSS**2,0.5D0*BX/auToT,MXI,1,ZI,1)
        CALL DAXPY_(NSS**2,0.5D0*BY/auToT,MYI,1,ZI,1)
        CALL DAXPY_(NSS**2,0.5D0*BZ/auToT,MZI,1,ZI,1)
        DO ISS=1,NSS
         HZER=ZR(ISS,ISS)
         ZR(ISS,ISS)=HZER+ENSOR(ISS)
        END DO
        CALL DCOPY_(NSS**2,[0.0D0],0,UZR,1)
        CALL DCOPY_(NSS   ,[1.0D0],0,UZR,NSS+1)
        CALL DCOPY_(NSS**2,[0.0D0],0,UZI,1)
        CALL ZJAC(NSS,ZR,ZI,NSS,UZR,UZI)
        DO IXYZ=1,3
         CALL DCOPY_(NSS**2,pMR(IXYZ)%A2,1,pZMR(IXYZ)%A2,1)
         CALL DCOPY_(NSS**2,pMI(IXYZ)%A2,1,pZMI(IXYZ)%A2,1)
         CALL DSCAL_(NSS**2,-0.5D0,pZMR(IXYZ)%A2,1)
         CALL DSCAL_(NSS**2,-0.5D0,pZMI(IXYZ)%A2,1)
         CALL ZTRNSF(NSS,UZR,UZI,pZMR(IXYZ)%A2,pZMI(IXYZ)%A2)
        ENDDO
        DO ITSTEP=1,NTSTEP
         T=TSTART+TINCRE*(ITSTEP-1)
         RkT=T*BOLTZ
         RMAGM(1)=0.0D0
         RMAGM(2)=0.0D0
         RMAGM(3)=0.0D0
         RPART=0.0D0
         DO ISS=1,NSS
          DELTA=ZR(ISS,ISS)-ZR(1,1)
          FACT=EXP(-DELTA/RkT)
          RMAGM(1)=RMAGM(1)+pZMR(1)%A2(ISS,ISS)*FACT
          RMAGM(2)=RMAGM(2)+pZMR(2)%A2(ISS,ISS)*FACT
          RMAGM(3)=RMAGM(3)+pZMR(3)%A2(ISS,ISS)*FACT
          RPART=RPART+FACT
         ENDDO
         RMAGM(1)=(RMAGM(1)/RPART)*AU2JTM
         RMAGM(2)=(RMAGM(2)/RPART)*AU2JTM
         RMAGM(3)=(RMAGM(3)/RPART)*AU2JTM
*        WRITE(6,'(6(1x,es12.5,1x))')
*    &    B,THE,PHI,RMAGM(1),RMAGM(2),RMAGM(3)
C backtransformation in two steps, -phi and -theta
         A=RMAGM(1)
         B=RMAGM(2)
         RMAGM(1)=A*COS(PHI)+B*SIN(PHI)
         RMAGM(2)=B*COS(PHI)-A*SIN(PHI)
         A=RMAGM(1)
         B=RMAGM(3)
         RMAGM(1)=A*COS(THE)-B*SIN(THE)
         RMAGM(3)=B*COS(THE)+A*SIN(THE)
         DO IXYZ=1,3
          LMSTEP=LMSTEP+1
          MAGM(LMSTEP)=MAGM(LMSTEP)+RMAGM(IXYZ)
         ENDDO
        ENDDO
       ENDDO

       ENDDO
      ENDDO

      WRITE(6,*)
      LMSTEP=0
      DO IBSTEP=1,NBSTEP
       B=BSTART+BINCRE*(IBSTEP-1)
       DO ITSTEP=1,NTSTEP
        T=TSTART+TINCRE*(ITSTEP-1)
        DO IXYZ=1,3
        LMSTEP=LMSTEP+1
          RMAGM(IXYZ)=MAGM(LMSTEP)/NORIENT
        ENDDO
        RMAGM2=RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
        RMAGMO=SQRT(RMAGM2)
        WRITE(6,'(1x,f6.2,5(1x,es12.5,1x))')
     &   T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
       ENDDO
      ENDDO

      CALL mma_deallocate(MAGM)

 810  CONTINUE

      WRITE(6,*)

      Call mma_deallocate(ZR)
      Call mma_deallocate(ZI)
      Call mma_deallocate(UZR)
      Call mma_deallocate(UZI)

      Call mma_deallocate(ZXR)
      Call mma_deallocate(ZXI)
      Call mma_deallocate(ZYR)
      Call mma_deallocate(ZYI)
      Call mma_deallocate(ZZR)
      Call mma_deallocate(ZZI)

      Call mma_deallocate(MXR)
      Call mma_deallocate(MXI)
      Call mma_deallocate(MYR)
      Call mma_deallocate(MYI)
      Call mma_deallocate(MZR)
      Call mma_deallocate(MZI)
      nullify(pMR(1)%A2,pMI(1)%A2,pMR(2)%A2,pMI(2)%A2,pMR(3)%A2,
     &        pMI(3)%A2)

 900  CONTINUE

30    FORMAT (5X,A,1X,ES15.8)
31    FORMAT (5X,2(1X,A4),6X,A15,1X,A47,1X,A15)
32    FORMAT (5X,95('-'))
33    FORMAT (5X,2(1X,I4),5X,5(1X,ES15.8))
35    FORMAT (5X,31('-'))
36    FORMAT (5X,2(1X,I4),6X,15('-'),1X,ES15.8,1X,A15)
37    FORMAT (5X,2(1X,I4),6X,15('-'),1X,A15,1X,ES15.8)
38    FORMAT (5X,2(1X,I4),6X,F15.6,4(1X,ES15.8))
39    FORMAT (5X,2(1X,A4),6X,A15,1X,A15,1X,A15)
40    FORMAT (5X,63('-'))
43    FORMAT (12X,A8,6(1X,ES15.6))
44    FORMAT (20X,6(1X,A15))
49    FORMAT (5X,A,1X,ES15.8,1X,A)

      Contains

      Subroutine Allocate_and_Load_electric_dipoles(IFANY)
      Integer ISOPR
      Integer IPRDX, IPRDY, IPRDZ, IFANY
         IPRDX=0
         IPRDY=0
         IPRDZ=0
         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'MLTPL  1'.AND.
     &        SOPRTP(ISOPR).EQ.'HERMSING') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
           END IF
         END DO
         CALL mma_allocate(DXR,NSS,NSS,Label='DXR')
         CALL mma_allocate(DXI,NSS,NSS,Label='DXI')
         CALL mma_allocate(DYR,NSS,NSS,Label='DYR')
         CALL mma_allocate(DYI,NSS,NSS,Label='DYI')
         CALL mma_allocate(DZR,NSS,NSS,Label='DZR')
         CALL mma_allocate(DZI,NSS,NSS,Label='DZI')
         DXR(:,:)=0.0D0
         DXI(:,:)=0.0D0
         DYR(:,:)=0.0D0
         DYI(:,:)=0.0D0
         DZR(:,:)=0.0D0
         DZI(:,:)=0.0D0
         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,DXR,NSS,IPRDX,0)
          CALL ZTRNSF(NSS,USOR,USOI,DXR,DXI)
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,DYR,NSS,IPRDY,0)
          CALL ZTRNSF(NSS,USOR,USOI,DYR,DYI)
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,DZR,NSS,IPRDZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,DZR,DZI)
         END If
      End Subroutine Allocate_and_Load_electric_dipoles

      Subroutine Allocate_and_Load_velocities(IFANY)
      Integer ISOPR
      Integer IPRDX, IPRDY, IPRDZ, IFANY
         IPRDX=0
         IPRDY=0
         IPRDZ=0
         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'VELOCITY') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
           END IF
         END DO
         CALL mma_allocate(DXR,NSS,NSS,Label='DXR')
         CALL mma_allocate(DXI,NSS,NSS,Label='DXI')
         CALL mma_allocate(DYR,NSS,NSS,Label='DYR')
         CALL mma_allocate(DYI,NSS,NSS,Label='DYI')
         CALL mma_allocate(DZR,NSS,NSS,Label='DZR')
         CALL mma_allocate(DZI,NSS,NSS,Label='DZI')
         DXR(:,:)=0.0D0
         DXI(:,:)=0.0D0
         DYR(:,:)=0.0D0
         DYI(:,:)=0.0D0
         DZR(:,:)=0.0D0
         DZI(:,:)=0.0D0
         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,DXR,NSS,IPRDX,0)
          CALL ZTRNSF(NSS,USOR,USOI,DXR,DXI)
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,DYR,NSS,IPRDY,0)
          CALL ZTRNSF(NSS,USOR,USOI,DYR,DYI)
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,DZR,NSS,IPRDZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,DZR,DZI)
         END If
      End Subroutine Allocate_and_Load_velocities

      Subroutine Deallocate_electric_dipoles()
         CALL mma_deallocate(DXR)
         CALL mma_deallocate(DXI)
         CALL mma_deallocate(DYR)
         CALL mma_deallocate(DYI)
         CALL mma_deallocate(DZR)
         CALL mma_deallocate(DZI)
      End Subroutine Deallocate_electric_dipoles

      Subroutine Allocate_and_Load_magnetic_dipoles(IFANY)
      Integer ISOPR
      Integer IPRMDX, IPRMDY, IPRMDZ, IFANY
         IPRMDX=0
         IPRMDY=0
         IPRMDZ=0

         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'ANGMOM  ') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRMDX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRMDY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRMDZ=ISOPR
           END IF
         END DO
         CALL mma_allocate(MDXR,NSS,NSS,Label='MDXR')
         CALL mma_allocate(MDXI,NSS,NSS,Label='MDXI')
         CALL mma_allocate(MDYR,NSS,NSS,Label='MDYR')
         CALL mma_allocate(MDYI,NSS,NSS,Label='MDYI')
         CALL mma_allocate(MDZR,NSS,NSS,Label='MDZR')
         CALL mma_allocate(MDZI,NSS,NSS,Label='MDZI')
         MDXR(:,:)=0.0D0
         MDXI(:,:)=0.0D0
         MDYR(:,:)=0.0D0
         MDYI(:,:)=0.0D0
         MDZR(:,:)=0.0D0
         MDZI(:,:)=0.0D0
         IF(IPRMDX.GT.0) THEN
          CALL SMMAT(PROP,MDXR,NSS,IPRMDX,0)
          CALL ZTRNSF(NSS,USOR,USOI,MDXR,MDXI)
         END IF
         IF(IPRMDY.GT.0) THEN
          CALL SMMAT(PROP,MDYR,NSS,IPRMDY,0)
          CALL ZTRNSF(NSS,USOR,USOI,MDYR,MDYI)
         END IF
         IF(IPRMDZ.GT.0) THEN
          CALL SMMAT(PROP,MDZR,NSS,IPRMDZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,MDZR,MDZI)
         END If
      End Subroutine Allocate_and_Load_magnetic_dipoles

      Subroutine Deallocate_magnetic_dipoles()
         CALL mma_deallocate(MDXR)
         CALL mma_deallocate(MDXI)
         CALL mma_deallocate(MDYR)
         CALL mma_deallocate(MDYI)
         CALL mma_deallocate(MDZR)
         CALL mma_deallocate(MDZI)
      End Subroutine Deallocate_magnetic_dipoles

      Subroutine Allocate_and_Load_Spin_Magnetic_dipoles(IFANY)
      Integer ISOPR
      Integer IPRSX, IPRSY, IPRSZ, IFANY
         IPRSX=0
         IPRSY=0
         IPRSZ=0

         IFANY=0
         DO ISOPR=1,NSOPR
            IF(SOPRNM(ISOPR).EQ.'MLTPL  0'.AND.
     &         SOPRTP(ISOPR).EQ.'ANTITRIP') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRSX=ISOPR
            IF(ISOCMP(ISOPR).EQ.1) IPRSY=ISOPR
            IF(ISOCMP(ISOPR).EQ.1) IPRSZ=ISOPR
           END IF
         END DO
         CALL mma_allocate(SXR,NSS,NSS,Label='SXR')
         CALL mma_allocate(SXI,NSS,NSS,Label='SXI')
         CALL mma_allocate(SYR,NSS,NSS,Label='SYR')
         CALL mma_allocate(SYI,NSS,NSS,Label='SYI')
         CALL mma_allocate(SZR,NSS,NSS,Label='SZR')
         CALL mma_allocate(SZI,NSS,NSS,Label='SZI')
         SXR(:,:)=0.0D0
         SXI(:,:)=0.0D0
         SYR(:,:)=0.0D0
         SYI(:,:)=0.0D0
         SZR(:,:)=0.0D0
         SZI(:,:)=0.0D0
         IF(IPRSX.GT.0) THEN
          CALL SMMAT(PROP,SXR,NSS,IPRSX,1)
          CALL ZTRNSF(NSS,USOR,USOI,SXR,SXI)
         END IF
         IF(IPRSY.GT.0) THEN
          CALL SMMAT(PROP,SYR,NSS,IPRSY,2)
          CALL ZTRNSF(NSS,USOR,USOI,SYR,SYI)
         END IF
         IF(IPRSZ.GT.0) THEN
          CALL SMMAT(PROP,SZR,NSS,IPRSZ,3)
          CALL ZTRNSF(NSS,USOR,USOI,SZR,SZI)
         END IF
      End Subroutine Allocate_and_Load_Spin_Magnetic_dipoles

      Subroutine Deallocate_Spin_Magnetic_dipoles()
         CALL mma_deallocate(SXR)
         CALL mma_deallocate(SXI)
         CALL mma_deallocate(SYR)
         CALL mma_deallocate(SYI)
         CALL mma_deallocate(SZR)
         CALL mma_deallocate(SZI)
      End Subroutine Deallocate_Spin_Magnetic_dipoles

      Subroutine Allocate_and_Load_Spin_Magnetic_Quadrupoles(IFANY)
      Integer ISOPR
      Integer IPRSXY, IPRSXZ, IPRSYX, IPRSYZ, IPRSZX, IPRSZY, IFANY
         IPRSXY=0
         IPRSXZ=0

         IPRSYX=0
         IPRSYZ=0

         IPRSZX=0
         IPRSZY=0
         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'MLTPL  1'.AND.
     &             SOPRTP(ISOPR).EQ.'ANTITRIP') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRSXY=ISOPR
            IF(ISOCMP(ISOPR).EQ.1) IPRSXZ=ISOPR

            IF(ISOCMP(ISOPR).EQ.2) IPRSYX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRSYZ=ISOPR

            IF(ISOCMP(ISOPR).EQ.3) IPRSZX=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRSZY=ISOPR

           END IF
         END DO
         CALL mma_allocate(SZXR,NSS,NSS,Label='SZXR')
         CALL mma_allocate(SZXI,NSS,NSS,Label='SZXI')
         SZXR(:,:)=0.0D0
         SZXI(:,:)=0.0D0
         CALL mma_allocate(SXZR,NSS,NSS,Label='SXZR')
         CALL mma_allocate(SXZI,NSS,NSS,Label='SXZI')
         SXZR(:,:)=0.0D0
         SXZI(:,:)=0.0D0

         CALL mma_allocate(SXYR,NSS,NSS,Label='SXYR')
         CALL mma_allocate(SXYI,NSS,NSS,Label='SXYI')
         SXYR(:,:)=0.0D0
         SXYI(:,:)=0.0D0
         CALL mma_allocate(SYXR,NSS,NSS,Label='SYXR')
         CALL mma_allocate(SYXI,NSS,NSS,Label='SYXI')
         SYXR(:,:)=0.0D0
         SYXI(:,:)=0.0D0

         CALL mma_allocate(SYZR,NSS,NSS,Label='SYZR')
         CALL mma_allocate(SYZI,NSS,NSS,Label='SYZI')
         SYZR(:,:)=0.0D0
         SYZI(:,:)=0.0D0
         CALL mma_allocate(SZYR,NSS,NSS,Label='SZYR')
         CALL mma_allocate(SZYI,NSS,NSS,Label='SZYI')
         SZYR(:,:)=0.0D0
         SZYI(:,:)=0.0D0
         IF(IPRSXY.GT.0) THEN
          CALL SMMAT(PROP,SXYR,NSS,IPRSXY,2)
          CALL ZTRNSF(NSS,USOR,USOI,SXYR,SXYI)
         END IF
         IF(IPRSYX.GT.0) THEN
          CALL SMMAT(PROP,SYXR,NSS,IPRSYX,1)
          CALL ZTRNSF(NSS,USOR,USOI,SYXR,SYXI)
         END IF

         IF(IPRSXZ.GT.0) THEN
          CALL SMMAT(PROP,SXZR,NSS,IPRSXZ,3)
          CALL ZTRNSF(NSS,USOR,USOI,SXZR,SXZI)
         END IF
         IF(IPRSZX.GT.0) THEN
          CALL SMMAT(PROP,SZXR,NSS,IPRSZX,1)
          CALL ZTRNSF(NSS,USOR,USOI,SZXR,SZXI)
         END IF

         IF(IPRSYZ.GT.0) THEN
          CALL SMMAT(PROP,SYZR,NSS,IPRSYZ,3)
          CALL ZTRNSF(NSS,USOR,USOI,SYZR,SYZI)
         END IF
         IF(IPRSZY.GT.0) THEN
          CALL SMMAT(PROP,SZYR,NSS,IPRSZY,2)
          CALL ZTRNSF(NSS,USOR,USOI,SZYR,SZYI)
         END IF
      End Subroutine Allocate_and_Load_Spin_Magnetic_Quadrupoles

      Subroutine Deallocate_Spin_Magnetic_Quadrupoles()
         Call mma_deallocate(SXYR)
         Call mma_deallocate(SXYI)
         Call mma_deallocate(SYXR)
         Call mma_deallocate(SYXI)

         Call mma_deallocate(SYZR)
         Call mma_deallocate(SYZI)
         Call mma_deallocate(SZYR)
         Call mma_deallocate(SZYI)

         Call mma_deallocate(SZXR)
         Call mma_deallocate(SZXI)
         Call mma_deallocate(SXZR)
         Call mma_deallocate(SXZI)
      End Subroutine Deallocate_Spin_Magnetic_Quadrupoles

      Subroutine Allocate_and_Load_Electric_Quadrupoles(IFANY)
      Integer ISOPR
      Integer IPRDXX, IPRDXY, IPRDXZ, IPRDYY, IPRDYZ, IPRDZZ, IFANY
         IPRDXX=0
         IPRDXY=0
         IPRDXZ=0
         IPRDYY=0
         IPRDYZ=0
         IPRDZZ=0

         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'MLTPL  2') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRDXX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRDXY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRDXZ=ISOPR
            IF(ISOCMP(ISOPR).EQ.4) IPRDYY=ISOPR
            IF(ISOCMP(ISOPR).EQ.5) IPRDYZ=ISOPR
            IF(ISOCMP(ISOPR).EQ.6) IPRDZZ=ISOPR
           END IF
         END DO
         CALL mma_allocate(QXXR,NSS,NSS,Label='QXXR')
         CALL mma_allocate(QXXI,NSS,NSS,Label='DXXI')
         CALL mma_allocate(QXYR,NSS,NSS,Label='DXYR')
         CALL mma_allocate(QXYI,NSS,NSS,Label='DXYI')
         CALL mma_allocate(QXZR,NSS,NSS,Label='DXZR')
         CALL mma_allocate(QXZI,NSS,NSS,Label='DXZI')
         CALL mma_allocate(QYYR,NSS,NSS,Label='DYYR')
         CALL mma_allocate(QYYI,NSS,NSS,Label='DYYI')
         CALL mma_allocate(QYZR,NSS,NSS,Label='DYZR')
         CALL mma_allocate(QYZI,NSS,NSS,Label='DYZI')
         CALL mma_allocate(QZZR,NSS,NSS,Label='DZZR')
         CALL mma_allocate(QZZI,NSS,NSS,Label='DZZI')
         QXXR(:,:)=0.0D0
         QXXI(:,:)=0.0D0
         QXYR(:,:)=0.0D0
         QXYI(:,:)=0.0D0
         QXZR(:,:)=0.0D0
         QXZI(:,:)=0.0D0
         QYYR(:,:)=0.0D0
         QYYI(:,:)=0.0D0
         QYZR(:,:)=0.0D0
         QYZI(:,:)=0.0D0
         QZZR(:,:)=0.0D0
         QZZI(:,:)=0.0D0
         IF(IPRDXX.GT.0) THEN
          CALL SMMAT(PROP,QXXR,NSS,IPRDXX,0)
          CALL ZTRNSF(NSS,USOR,USOI,QXXR,QXXI)
         END IF
         IF(IPRDXY.GT.0) THEN
          CALL SMMAT(PROP,QXYR,NSS,IPRDXY,0)
          CALL ZTRNSF(NSS,USOR,USOI,QXYR,QXYI)
         END IF
         IF(IPRDXZ.GT.0) THEN
          CALL SMMAT(PROP,QXZR,NSS,IPRDXZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,QXZR,QXZI)
         END IF
         IF(IPRDYY.GT.0) THEN
          CALL SMMAT(PROP,QYYR,NSS,IPRDYY,0)
          CALL ZTRNSF(NSS,USOR,USOI,QYYR,QYYI)
         END IF
         IF(IPRDYZ.GT.0) THEN
          CALL SMMAT(PROP,QYZR,NSS,IPRDYZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,QYZR,QYZI)
         END IF
         IF(IPRDZZ.GT.0) THEN
          CALL SMMAT(PROP,QZZR,NSS,IPRDZZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,QZZR,QZZI)
         END IF
      End Subroutine Allocate_and_Load_Electric_Quadrupoles

      Subroutine Deallocate_Electric_Quadrupoles()
         CALL mma_deallocate(QXXR)
         CALL mma_deallocate(QXXI)
         CALL mma_deallocate(QXYR)
         CALL mma_deallocate(QXYI)
         CALL mma_deallocate(QXZR)
         CALL mma_deallocate(QXZI)
         CALL mma_deallocate(QYYR)
         CALL mma_deallocate(QYYI)
         CALL mma_deallocate(QYZR)
         CALL mma_deallocate(QYZI)
         CALL mma_deallocate(QZZR)
         CALL mma_deallocate(QZZI)
      End Subroutine Deallocate_Electric_Quadrupoles

      Subroutine Allocate_and_Load_Magnetic_Quadrupoles(IFANY)
      Integer ISOPR
      Integer IPRDZX, IPRDYX, IPRDZY, IFANY
         IPRDXY=0
         IPRDXZ=0
         IPRDYX=0
         IPRDYZ=0
         IPRDZX=0
         IPRDZY=0

         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'OMQ') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.2) IPRDXY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRDXZ=ISOPR

            IF(ISOCMP(ISOPR).EQ.4) IPRDYX=ISOPR
            IF(ISOCMP(ISOPR).EQ.6) IPRDYZ=ISOPR

            IF(ISOCMP(ISOPR).EQ.7) IPRDZX=ISOPR
            IF(ISOCMP(ISOPR).EQ.8) IPRDZY=ISOPR
           END IF
         END DO

         CALL mma_allocate(MQXYR,NSS,NSS,Label='MQXYR')
         CALL mma_allocate(MQXYI,NSS,NSS,Label='MQXYI')
         CALL mma_allocate(MQYXR,NSS,NSS,Label='MQYXR')
         CALL mma_allocate(MQYXI,NSS,NSS,Label='MQYXI')
         CALL mma_allocate(MQXZR,NSS,NSS,Label='MQXZR')
         CALL mma_allocate(MQXZI,NSS,NSS,Label='MQXZI')
         CALL mma_allocate(MQZXR,NSS,NSS,Label='MQZXR')
         CALL mma_allocate(MQZXI,NSS,NSS,Label='MQZXI')
         CALL mma_allocate(MQYZR,NSS,NSS,Label='MQYZR')
         CALL mma_allocate(MQYZI,NSS,NSS,Label='MQYZI')
         CALL mma_allocate(MQZYR,NSS,NSS,Label='MQZYR')
         CALL mma_allocate(MQZYI,NSS,NSS,Label='MQZYI')
         MQXYR(:,:)=0.0D0
         MQXYI(:,:)=0.0D0
         MQYXR(:,:)=0.0D0
         MQYXI(:,:)=0.0D0
         MQXZR(:,:)=0.0D0
         MQXZI(:,:)=0.0D0
         MQZXR(:,:)=0.0D0
         MQZXI(:,:)=0.0D0
         MQYZR(:,:)=0.0D0
         MQYZI(:,:)=0.0D0
         MQZYR(:,:)=0.0D0
         MQZYI(:,:)=0.0D0
         IF(IPRDXY.GT.0) THEN
          CALL SMMAT(PROP,MQXYR,NSS,IPRDXY,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQXYR,MQXYI)
         END IF
         IF(IPRDYX.GT.0) THEN
          CALL SMMAT(PROP,MQYXR,NSS,IPRDYX,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQYXR,MQYXI)
         END IF

         IF(IPRDXZ.GT.0) THEN
          CALL SMMAT(PROP,MQXZR,NSS,IPRDXZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQXZR,MQXZI)
         END IF
         IF(IPRDZX.GT.0) THEN
          CALL SMMAT(PROP,MQZXR,NSS,IPRDZX,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQZXR,MQZXI)
         END IF

         IF(IPRDYZ.GT.0) THEN
          CALL SMMAT(PROP,MQYZR,NSS,IPRDYZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQYZR,MQYZI)
         END IF
         IF(IPRDZY.GT.0) THEN
          CALL SMMAT(PROP,MQZYR,NSS,IPRDZY,0)
          CALL ZTRNSF(NSS,USOR,USOI,MQZYR,MQZYI)
         END IF
      End Subroutine Allocate_and_Load_Magnetic_Quadrupoles

      Subroutine Deallocate_Magnetic_Quadrupoles()
         CALL mma_deallocate(MQXYR)
         CALL mma_deallocate(MQXYI)
         CALL mma_deallocate(MQYXR)
         CALL mma_deallocate(MQYXI)
         CALL mma_deallocate(MQXZR)
         CALL mma_deallocate(MQXZI)
         CALL mma_deallocate(MQZXR)
         CALL mma_deallocate(MQZXI)
         CALL mma_deallocate(MQYZR)
         CALL mma_deallocate(MQYZI)
         CALL mma_deallocate(MQZYR)
         CALL mma_deallocate(MQZYI)
      End Subroutine Deallocate_Magnetic_Quadrupoles

      Subroutine Allocate_and_Load_Octupoles(IFANY)
      Integer ISOPR
      Integer IPRDZZX, IPRDZZY, IPRDZZZ, IPRDXXX, IPRDXXY, IPRDXXZ,
     &        IPRDYYX, IPRDYYY, IPRDYYZ, IFANY
! This is a real symmetric rank 3 tensor so only 10 and not 27 is needed
! The order which comes in
         IPRDXXX=0 !
         IPRDXXY=0 !
         IPRDXXZ=0 !

!        IPRDXYX=0
!        IPRDXYY=0 ! YYX These are the same due to symmetry
!        IPRDXYZ=0 ! Not present

!        IPRDXZX=0
!        IPRDXZY=0
!        IPRDXZZ=0 ! ZZX

!        IPRDYXX=0
!        IPRDYXY=0
!        IPRDYXZ=0

         IPRDYYX=0 ! Taking the XYY order
         IPRDYYY=0 !
         IPRDYYZ=0 !

!        IPRDYZX=0
!        IPRDYZY=0
!        IPRDYZZ=0 ! ZZY

!        IPRDZXX=0
!        IPRDZXY=0
!        IPRDZXZ=0

!        IPRDZYX=0
!        IPRDZYY=0
!        IPRDZYZ=0

         IPRDZZX=0 ! Taking order from XZZ
         IPRDZZY=0 ! Taking order from YZZ
         IPRDZZZ=0 !

         IFANY=0
         DO ISOPR=1,NSOPR
           IF(SOPRNM(ISOPR).EQ.'MLTPL  3') THEN
            IFANY=1
            IF(ISOCMP(ISOPR).EQ.1) IPRDXXX=ISOPR
            IF(ISOCMP(ISOPR).EQ.2) IPRDXXY=ISOPR
            IF(ISOCMP(ISOPR).EQ.3) IPRDXXZ=ISOPR
            IF(ISOCMP(ISOPR).EQ.4) IPRDYYX=ISOPR ! Changed from XYY
            !IF(ISOCMP(ISOPR).EQ.5) IPRDXYZ=ISOPR
            IF(ISOCMP(ISOPR).EQ.6) IPRDZZX=ISOPR ! Changed from XZZ
            IF(ISOCMP(ISOPR).EQ.7) IPRDYYY=ISOPR
            IF(ISOCMP(ISOPR).EQ.8) IPRDYYZ=ISOPR
            IF(ISOCMP(ISOPR).EQ.9) IPRDZZY=ISOPR ! Changed from YZZ
            IF(ISOCMP(ISOPR).EQ.10) IPRDZZZ=ISOPR

           END IF
         END DO
         CALL mma_allocate(DXXXR,NSS,NSS,Label='DXXXR')
         CALL mma_allocate(DXXXI,NSS,NSS,Label='DXXXI')
         CALL mma_allocate(DXXYR,NSS,NSS,Label='DXXYR')
         CALL mma_allocate(DXXYI,NSS,NSS,Label='DXXYI')
         CALL mma_allocate(DXXZR,NSS,NSS,Label='DXXZR')
         CALL mma_allocate(DXXZI,NSS,NSS,Label='DXXZI')
         DXXXR(:,:)=0.0D0
         DXXXI(:,:)=0.0D0
         DXXYR(:,:)=0.0D0
         DXXYI(:,:)=0.0D0
         DXXZR(:,:)=0.0D0
         DXXZI(:,:)=0.0D0
         CALL mma_allocate(DYYXR,NSS,NSS,Label='DYYXR')
         CALL mma_allocate(DYYXI,NSS,NSS,Label='DYYXI')
         CALL mma_allocate(DYYYR,NSS,NSS,Label='DYYYR')
         CALL mma_allocate(DYYYI,NSS,NSS,Label='DYYYI')
         CALL mma_allocate(DYYZR,NSS,NSS,Label='DYYZR')
         CALL mma_allocate(DYYZI,NSS,NSS,Label='DYYZI')
         DYYXR(:,:)=0.0D0
         DYYXI(:,:)=0.0D0
         DYYYR(:,:)=0.0D0
         DYYYI(:,:)=0.0D0
         DYYZR(:,:)=0.0D0
         DYYZI(:,:)=0.0D0
         CALL mma_allocate(DZZXR,NSS,NSS,Label='DZZXR')
         CALL mma_allocate(DZZXI,NSS,NSS,Label='DZZXI')
         CALL mma_allocate(DZZYR,NSS,NSS,Label='DZZYR')
         CALL mma_allocate(DZZYI,NSS,NSS,Label='DZZYI')
         CALL mma_allocate(DZZZR,NSS,NSS,Label='DZZZR')
         CALL mma_allocate(DZZZI,NSS,NSS,Label='DZZZI')
         DZZXR(:,:)=0.0D0
         DZZXI(:,:)=0.0D0
         DZZYR(:,:)=0.0D0
         DZZYI(:,:)=0.0D0
         DZZZR(:,:)=0.0D0
         DZZZI(:,:)=0.0D0
         IF(IPRDXXX.GT.0) THEN
          CALL SMMAT(PROP,DXXXR,NSS,IPRDXXX,0)
          CALL ZTRNSF(NSS,USOR,USOI,DXXXR,DXXXI)
         END IF
         IF(IPRDXXY.GT.0) THEN
          CALL SMMAT(PROP,DXXYR,NSS,IPRDXXY,0)
          CALL ZTRNSF(NSS,USOR,USOI,DXXYR,DXXYI)
         END IF
         IF(IPRDXXZ.GT.0) THEN
          CALL SMMAT(PROP,DXXZR,NSS,IPRDXXZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,DXXZR,DXXZI)
         END IF

         IF(IPRDYYX.GT.0) THEN
          CALL SMMAT(PROP,DYYXR,NSS,IPRDYYX,0)
          CALL ZTRNSF(NSS,USOR,USOI,DYYXR,DYYXI)
         END IF
         IF(IPRDYYY.GT.0) THEN
          CALL SMMAT(PROP,DYYYR,NSS,IPRDYYY,0)
          CALL ZTRNSF(NSS,USOR,USOI,DYYYR,DYYYI)
         END IF
         IF(IPRDYYZ.GT.0) THEN
          CALL SMMAT(PROP,DYYZR,NSS,IPRDYYZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,DYYZR,DYYZI)
         END IF

         IF(IPRDZZX.GT.0) THEN
          CALL SMMAT(PROP,DZZXR,NSS,IPRDZZX,0)
          CALL ZTRNSF(NSS,USOR,USOI,DZZXR,DZZXI)
         END IF
         IF(IPRDZZY.GT.0) THEN
          CALL SMMAT(PROP,DZZYR,NSS,IPRDZZY,0)
          CALL ZTRNSF(NSS,USOR,USOI,DZZYR,DZZYI)
         END IF
         IF(IPRDZZZ.GT.0) THEN
          CALL SMMAT(PROP,DZZZR,NSS,IPRDZZZ,0)
          CALL ZTRNSF(NSS,USOR,USOI,DZZZR,DZZZI)
         END IF
      End Subroutine Allocate_and_Load_Octupoles

      Subroutine deallocate_Octupoles()
         CALL mma_deallocate(DXXXR)
         CALL mma_deallocate(DXXXI)
         CALL mma_deallocate(DXXYR)
         CALL mma_deallocate(DXXYI)
         CALL mma_deallocate(DXXZR)
         CALL mma_deallocate(DXXZI)
         CALL mma_deallocate(DYYXR)
         CALL mma_deallocate(DYYXI)
         CALL mma_deallocate(DYYYR)
         CALL mma_deallocate(DYYYI)
         CALL mma_deallocate(DYYZR)
         CALL mma_deallocate(DYYZI)
         CALL mma_deallocate(DZZXR)
         CALL mma_deallocate(DZZXI)
         CALL mma_deallocate(DZZYR)
         CALL mma_deallocate(DZZYI)
         CALL mma_deallocate(DZZZR)
         CALL mma_deallocate(DZZZI)
      End Subroutine deallocate_Octupoles

      END SUBROUTINE PRPROP

      SUBROUTINE SINANI(KDGN,IFUNCT,NSS,DIPSOm,SPNSFS,DIPSOm_SA)
      use definitions, only: iwp, wp, u6
!      IMPLICIT NONE
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER(kind=iwp), intent(in):: KDGN,IFUNCT,NSS
      COMPLEX(kind=wp), intent(in):: DIPSOm(3,NSS,NSS)
      COMPLEX(kind=wp), intent(in):: SPNSFS(3,NSS,NSS)
      real(kind=wp), intent(in):: DIPSOm_SA

      COMPLEX(kind=wp) DIPSOmSA(3,KDGN,KDGN)
      INTEGER(kind=iwp) l,Iso1,Jso2,Ico1,i,j
      COMPLEX(kind=wp) SPNSOSA(3,KDGN,KDGN)
      COMPLEX(kind=wp) Z(NSS,NSS),MATL(NSS,NSS),FINL(NSS,NSS)
      COMPLEX(kind=wp) SPNSO(3,NSS,NSS)
      real(kind=wp) UMATR(NSS,NSS),UMATI(NSS,NSS),gtens(3),maxes(3,3)
      CHARACTER(LEN=1) angm

      if(.FALSE.) then
      write(u6,'(/)')
      write(u6,'(10A)') (('############'),J=1,10)
      write(u6,'(25X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT IN '//
     &'THE BASIS OF SPIN ORBIT STATES'
      write(u6,'(10A)') (('############'),J=1,10)

      do l=1,3
      if(l.eq.1)  angm='X'
      if(l.eq.2)  angm='Y'
      if(l.eq.3)  angm='Z'
      write(u6,'(/)')
      write(u6,'(4X,A12,A2)') 'PROJECTION: ', angm
      write(u6,'(/)')
       do Iso1=1,NSS
      write(u6,'(20(2X,2F10.6))') (DIPSOm(l,Iso1,Jso2), Jso2=1,NSS)
       enddo
      enddo
      write(6,'(/)')

       endif !if(IPGLOB.GE.4)

         do l=1,3
         do Ico1=1,KDGN
         do Jco1=1,KDGN
        DIPSOmSA(l,Ico1,Jco1)=(0.d0,0.d0)
      !S_SOM( L,I,J)=(0.d0,0.d0)
         enddo
         enddo
         enddo


         do Iso1=1,KDGN
         do Jso2=1,KDGN
        Ico1=Iso1+IFUNCT
        Jco1=Jso2+IFUNCT
         do l=1,3
      !write(6,*)'DIPSOm',DIPSOm(l,Ico1,Jco1)
       DIPSOmSA(l,Iso1,Jso2)=-DIPSOm(l,Ico1,Jco1)
      !S_SOM( l,i,j)=S_SO(l,ic1,ic2)
        enddo
        enddo
        enddo

         if(.False.) then
      write(u6,*)
      write(u6,'(10X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT '//
     & 'IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
      do l=1,3
      write(u6,'(/)')
      write(u6,'(5X,A6,I3)') 'AXIS= ',l
      write(u6,*)
       do Ico1=1,KDGN
      write(u6,'(16(2F12.8,2x))') (DIPSOmSA(l,Ico1,Jco1), Jco1=1,KDGN)
       enddo
      enddo

       endif ! if(IPGLOB.GE.4)

          CALL ATENS_RASSI(DIPSOmSA, KDGN, gtens, maxes, 3)

           if(.False.) then

          call get_dArray('UMATR_SINGLE',UMATR,NSS**2)
          call get_dArray('UMATI_SINGLE',UMATI,NSS**2)
       write(u6,'(/)')
       write(u6,'(5x,a)') 'umatr and umati'
       write(u6,'(/)')
       do i=1,NSS
       write(u6,'(5x,10(2f14.10,2x))') (UMATR(i,j),UMATI(i,j),j=1,NSS)
       enddo

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       SPNSO(L,I,J)=(0.0d0,0.0d0)
          enddo
       Z(I,J)=(0.0d0,0.0d0)
         enddo
       enddo


        do i=1,NSS
        do j=1,NSS
        Z(i,j)=Z(i,j)+CMPLX(UMATR(i,j),UMATI(i,j),kind=8)
        enddo
        enddo

       do l=1,3
         do i=1,NSS
            do j=1,NSS
      MATL(i,j)=(0.0d0,0.0d0)
      MATL(i,j)=SPNSFS(L,i,j)
            enddo
         enddo


         do i=1,NSS
            do j=1,NSS
      FINL(i,j)=(0.0d0,0.0d0)
            enddo
         enddo

       call ADARASSI(NSS,Z,MATL,FINL)

         do i=1,NSS
            do j=1,NSS
      SPNSO(L,i,j) = FINL(i,j)
            enddo
         enddo
      enddo !l

      write(u6,'(/)')
      write(u6,'(10A)') (('############'),J=1,10)
      write(u6,'(30X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT IN '//
     & 'THE BASIS OF SPIN ORBIT STATES'
      write(u6,'(10A)') (('############'),J=1,10)
      write(u6,'(/)')
      do l=1,3
      if(l.eq.1)  angm='X'
      if(l.eq.2)  angm='Y'
      if(l.eq.3)  angm='Z'
      write(u6,'(/)')
      write(u6,'(4X,A,A)') 'PROJECTION: ', angm
      write(u6,'(/)')
       do Iso1=1,NSS
      write(u6,'(20(2F10.6,2X))') (SPNSO(l,Iso1,Jso2), Jso2=1,NSS)
       enddo
      enddo

         do l=1,3
         do Ico1=1,KDGN
         do Jco1=1,KDGN
        SPNSOSA(l,Ico1,Jco1)=(0.d0,0.d0)
         enddo
         enddo
         enddo

         do Iso1=1,KDGN
         do Jso2=1,KDGN
        Ico1=Iso1+IFUNCT
        Jco1=Jso2+IFUNCT
         do l=1,3
       SPNSOSA(l,Iso1,Jso2)=SPNSO(l,Ico1,Jco1)
        enddo
        enddo
        enddo

      write(u6,*)
      write(u6,'(10X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT '//
     & 'IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
      do l=1,3
      write(u6,'(/)')
      write(u6,'(5X,A6,I3)') 'AXIS= ',l
      write(u6,*)
       do Ico1=1,KDGN
      write(u6,'(16(2F12.8,2x))') (SPNSOSA(l,Ico1,Jco1), Jco1=1,KDGN)
       enddo
      enddo

       endif! if(IPGLOB.GT.3)

       !!do l=1,3
       !!do Iso1=1,KDGN
       !!do Jso2=1,KDGN
      !!Ico1=Iso1+IFUNCT
      !!Jco1=Jso2+IFUNCT
       ! do l=1,3
       !!write(6,*)'DIPSOm',DIPSOm(l,Ico1,Jco1)
      !DIPSOm_SA(l,Iso1,Jso2)=-DIPSOm(l,Ico1,Jco1)
      !S_SOM( l,i,j)=S_SO(l,ic1,ic2)
        !!enddo
       !!enddo
      !!enddo

c Avoid unused argument warnings
      unused_var(DIPSOm_SA)
      END SUBROUTINE SINANI

      SUBROUTINE ADARASSI(N,A,D,DROT)
      use definitions, only: iwp, wp
      use constants, only: Zero,One

      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: N
      COMPLEX(kind=wp), intent(in)::  A(N,N), D(N,N)
      COMPLEX(kind=wp), intent(out)::  DROT(N,N)

      INTEGER(kind=iwp) I, J
      COMPLEX(kind=wp)  TEMP(N,N)

C initialization
      do I=1,N
       do J=1,N
      DROT(I,J)=(Zero,Zero)
      TEMP(I,J)=(Zero,Zero)
       enddo
      enddo

C actual multiplication
      call ZGEMM('C','N',N,N,N,(One,Zero),A,N,D,N,(Zero,Zero),
     &TEMP,N)
      call ZGEMM('N','N',N,N,N,(One,Zero),TEMP,N,A,N,(Zero,Zero),
     &DROT,N)

      END SUBROUTINE ADARASSI

      SUBROUTINE ZECON(NSTATE,N,UR,UI,AR,AI,ZEKL,IXYZ,ISTATE,ISS,JSS)
      use definitions, only: iwp, wp
      IMPLICIT None
      integer(kind=iwp), intent(in):: NSTATE,N,IXYZ,ISTATE,ISS,JSS
      real(kind=wp), intent(in):: UR(N,N),UI(N,N)
      real(kind=wp), intent(in):: AR(N,N),AI(N,N)
      COMPLEX(kind=wp), intent(inout):: ZEKL(2,2,3,NSTATE)

      real(kind=wp) TMPR1,TMPR2,TMPI1,TMPI2

      TMPR1=AR(ISS,JSS)*UR(JSS,1)-AI(ISS,JSS)*UI(JSS,1)
      TMPR2=AR(ISS,JSS)*UR(JSS,2)-AI(ISS,JSS)*UI(JSS,2)
      TMPI1=AI(ISS,JSS)*UR(JSS,1)+AR(ISS,JSS)*UI(JSS,1)
      TMPI2=AI(ISS,JSS)*UR(JSS,2)+AR(ISS,JSS)*UI(JSS,2)

      ZEKL(1,1,IXYZ,ISTATE)=ZEKL(1,1,IXYZ,ISTATE)+
     $     CMPLX(UR(ISS,1)*TMPR1+UI(ISS,1)*TMPI1,
     $     UR(ISS,1)*TMPI1-UI(ISS,1)*TMPR1,kind=8)
      ZEKL(1,2,IXYZ,ISTATE)=ZEKL(1,2,IXYZ,ISTATE)+
     $     CMPLX(UR(ISS,1)*TMPR2+UI(ISS,1)*TMPI2,
     $     UR(ISS,1)*TMPI2-UI(ISS,1)*TMPR2,kind=8)
      ZEKL(2,1,IXYZ,ISTATE)=ZEKL(2,1,IXYZ,ISTATE)+
     $     CMPLX(UR(ISS,2)*TMPR1+UI(ISS,2)*TMPI1,
     $     UR(ISS,2)*TMPI1-UI(ISS,2)*TMPR1,kind=8)
      ZEKL(2,2,IXYZ,ISTATE)=ZEKL(2,2,IXYZ,ISTATE)+
     $     CMPLX(UR(ISS,2)*TMPR2+UI(ISS,2)*TMPI2,
     $     UR(ISS,2)*TMPI2-UI(ISS,2)*TMPR2,kind=8)

      END SUBROUTINE ZECON
