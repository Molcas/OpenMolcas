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
      SUBROUTINE PRPROP(PROP,USOR,USOI,ENSOR,NSS,OVLP,ENERGY,JBNUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='PRPROP')
      parameter (THRSH=1.0D-10)
      parameter (ZERO=0.0D0)
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "constants.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &          ENERGY(NSTATE),JBNUM(NSTATE)
#include "SysDef.fh"
#include "rassiwfn.fh"
      Character*1 xyzchr(3)
      Character*3 ASDLAB
      Character*8 EFPROP
      Character*8 PSOPROP
      Character*8 DMPPROP
*     Character*8 OVRPROP
      Dimension IPAMFI(3),IPAM(3),IZMR(3),IZMI(3)
      Dimension DTENS(3,3),GTENS(3,3),GSTENS(3,3),SOSTERM(9)
      Dimension TMPMAT(3,3),TMPVEC(3,3),EVR(3),EVI(3)
      COMPLEX*16 ZEKL(2,2,3,NSTATE),GCONT(9,NSTATE)
      COMPLEX*16 DIPSOm(3,NSS,NSS),Z(NSS,NSS),DIPSOn(3,NSS,NSS)
      COMPLEX*16 SPNSFS(3,NSS,NSS)
      COMPLEX*16 DIPSOf(3,NSS,NSS),DIMSO(3,3,NSS,NSS)
      COMPLEX*16 DIPSOfc(3,NSS,NSS),DIPSOfsd(3,NSS,NSS)
      COMPLEX*16 DIPSOfcsd(3,NSS,NSS),DIPSOfpso(3,NSS,NSS)
       !REAL*8  DIMSOIJ(3,3,NSS)
      REAL*8 GTOTAL(9),ANGMOME(3,NSTATE,NSTATE),ESO(NSS)
      REAL*8 EDIP1MOM(3,NSTATE,NSTATE),AMFIINT(3,NSTATE,NSTATE)
      REAL*8 TMPL(NSTATE,NSTATE,3),TMPE(NSTATE,NSTATE,3)
      REAL*8 TMPA(NSTATE,NSTATE,3)
      Dimension TMPm(NTS),TMPf(NTP)
*     Dimension TMPm(NTS),TMPf(NTP),TMFC(NTF)
      Dimension c_1(3,3),c_2(3,3)!,Zstat1m(NTS),Zstat1f(NTP)
      Dimension curit(3,3),paramt(3,3)
      Dimension HFC_1(3,3),HFC_2(3,3),HFC_3(3,3)
      Dimension CurieT(3,3),DiamT(3,3),PNMRCPS(NTP,NSS,3,3)
      Dimension chiT_tens(NTS,3,3),PNMRT(NTP,3,3),PNMR(NTP,3,3)
      Dimension chicuriT_tens(NTS,3,3),chiparamT_tens(NTS,3,3)
      Dimension PNMRC(NTP,3,3),PNMRD(NTP,3,3)
*     Dimension PNMRC(NTP,3,3),PNMRD(NTP,3,3),PNMRFCC(NTP,3,3)
*     Dimension NMRFT(NTF,3,3),NMRFP(NTF,3,3),NMRFC(NTF,3,3)
*     Dimension NMRFD(NTF,3,3)
      REAL*8 DLTTA,DLTT,Zstat,p_Boltz,Boltz_k,coeff_chi
      LOGICAL ISGS(NSS),IFANGM,IFDIP1,IFAMFI, Sparse_I,Sparse_J
      Dimension IMR(3),IMI(3),RMAGM(3),Chi(3)
      INTEGER IFUNCT, SECORD(4)
      REAL*8 J2CM
      Real*8 P1(3), P2(3), kxe1(3), kxe2(3)
      INTEGER IOFF(8)
      CHARACTER*8 LABEL
      Complex*16 T0(3), T1(3), TM1, TM2, E1A, E2A, E1B, E2B,
     &           IMAGINARY
      REAL*8 COMPARE


      CALL QENTER(ROUTINE)

      AVOGADRO=CONST_AVOGADRO_
      AU2EV=CONV_AU_TO_EV_
      AU2CM=CONV_AU_TO_CM1_
      AU2T=CONV_AU_TO_T_
      AU2J=CONV_AU_TO_KJ_*1.0D3
      J2CM=AU2CM/AU2J
      AU2JTM=(AU2J/AU2T)*AVOGADRO
      ALPHA=CONST_AU_VELOCITY_IN_SI_/CONST_C_IN_SI_
      ALPHA2= ALPHA*ALPHA
      IMAGINARY=DCMPLX(0.0D0,1.0D0)

      BOLTZ_K=CONST_BOLTZMANN_*J2CM
      coeff_chi=0.1D0*AVOGADRO/CONST_BOLTZMANN_*
     &          CONST_BOHR_MAGNETON_IN_SI_**2
      FEGVAL=-(CONST_ELECTRON_G_FACTOR_)
      BOLTZ=CONST_BOLTZMANN_/AU2J
      Rmu0=4.0D-7*CONST_PI_

      xyzchr(1)='x'
      xyzchr(2)='y'
      xyzchr(3)='z'

******************************************************
* printout of properties over the spin-free states
******************************************************

      IF(IPGLOB.LE.SILENT) GOTO 400

      IF( PRXVE.OR.PRMEE ) THEN
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,'(6X,100A1)') ('*',i=1,100)
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,A,34X,A,34X,A)')
     &     '*',' Spin-free properties section ','*'
      WRITE(6,'(6X,A,98X,A)') '*','*'
      WRITE(6,'(6X,100A1)') ('*',i=1,100)
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
        PMAX=0.0D0


        DO I=1,NSTATE
         PMAX=MAX(PMAX,ABS(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
        END DO
        IF(PMAX.LT.PLIMIT) GOTO 100

        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     *  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
          WRITE(6,'(1X,A,3D17.8)')
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
         WRITE(6,'(1X,A,3D17.8)')
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
      DO IPROP=1,NPROP
         IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
            IFANGM=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  ANGMOME(ICOMP(IPROP),I,J)=0.0D0
                  ANGMOME(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
                  TMPL(I,J,ICOMP(IPROP))=0.0D0
                  TMPL(I,J,ICOMP(IPROP))=PROP(I,J,IPROP)
               ENDDO
            ENDDO
c#ifdef _HDF5_
c            call mh5_put_dset_array_real(wfn_sfs_angmom,
c     $      PROP(1:NSTATE,1:NSTATE,ICOMP(IPROP)),
c     $      [NSTATE,NSTATE,1], [0,0,ICOMP(IPROP)-1])
c#endif
         ENDIF
c add dipole moment integrals:
         IF(PNAME(IPROP)(1:8).EQ.'MLTPL  1') THEN
            IFDIP1=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  EDIP1MOM(ICOMP(IPROP),I,J)=0.0D0
                  EDIP1MOM(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
                  TMPE(I,J,ICOMP(IPROP))=0.0D0
                  TMPE(I,J,ICOMP(IPROP))=PROP(I,J,IPROP)
               ENDDO
            ENDDO
c#ifdef _HDF5_
c            call mh5_put_dset_array_real(wfn_sfs_edipmom,
c     $      PROP(1:NSTATE,1:NSTATE,ICOMP(IPROP)),
c     $      [NSTATE,NSTATE,1], [0,0,ICOMP(IPROP)-1])
c#endif
         ENDIF
c add spin-orbit AMFI integrals:
         IF(PNAME(IPROP)(1:8).EQ.'AMFI    ') THEN
            IFAMFI=.TRUE.
            DO I=1,NSTATE
               DO J=1,NSTATE
                  AMFIINT(ICOMP(IPROP),I,J)=0.0D0
                  AMFIINT(ICOMP(IPROP),I,J)=PROP(I,J,IPROP)
                  TMPA(I,J,ICOMP(IPROP))=0.0D0
                  TMPA(I,J,ICOMP(IPROP))=PROP(I,J,IPROP)
               ENDDO
            ENDDO
c#ifdef _HDF5_
c            call mh5_put_dset_array_real(wfn_sfs_amfi,
c     $      PROP(1:NSTATE,1:NSTATE,ICOMP(IPROP)),
c     $      [NSTATE,NSTATE,1], [0,0,ICOMP(IPROP)-1])
c#endif
         ENDIF
      ENDDO
      IF(IFANGM.EQV..TRUE.) THEN
       CALL Put_dArray('ANGM_SINGLE',ANGMOME,3*NSTATE*NSTATE)
#ifdef _HDF5_
            call mh5_put_dset_array_real(wfn_sfs_angmom,
     $      TMPL(:,:,:),
     $      [NSTATE,NSTATE,3], [0,0,0])
#endif
      ENDIF
      IF(IFDIP1.EQV..TRUE.) THEN
       CALL Put_dArray('DIP1_SINGLE',EDIP1MOM,3*NSTATE*NSTATE)
#ifdef _HDF5_
            call mh5_put_dset_array_real(wfn_sfs_edipmom,
     $      TMPE(:,:,:),
     $      [NSTATE,NSTATE,3], [0,0,0])
#endif
      ENDIF
      IF(IFAMFI.EQV..TRUE.) THEN
       CALL Put_dArray('AMFI_SINGLE',AMFIINT,3*NSTATE*NSTATE)
#ifdef _HDF5_
            call mh5_put_dset_array_real(wfn_sfs_amfi,
     $      TMPA(:,:,:),
     $      [NSTATE,NSTATE,3], [0,0,0])
#endif
      ENDIF

*******************************************************
* printout of properties over the spin-orbit states
*******************************************************

      IF(.not.IFSO) GOTO 300
      NPMSIZ=NSOPR
      IF(NSOPR.EQ.0) GOTO 300

      IF( PRMES ) THEN
* match the SO property list to the SF property list
       CALL GETMEM('PMAP','ALLO','INTE',LPMAP,NPMSIZ)
       NMISS=0
       DO ISOPR=1,NSOPR
        IWORK(LPMAP-1+ISOPR)=0
        DO IPROP=1,NPROP
         IF(PNAME(IPROP).EQ.SOPRNM(ISOPR).AND.
     &     ICOMP(IPROP).EQ.ISOCMP(ISOPR)) THEN
          IWORK(LPMAP-1+ISOPR)=IPROP
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
          IF(IWORK(LPMAP-1+ISOPR).EQ.0)
     &       WRITE(6,*)'Property:',SOPRNM(ISOPR),
     &                 '      Component:',ISOCMP(ISOPR)
         END DO
        WRITE(6,*)
       END IF
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A,34X,A)')
     &     '*','Spin-orbit properties section ','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
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
        IPROP=IWORK(LPMAP-1+I)
        IF(IPROP.GT.0) THEN
         ISOPR=ISOPR+1
         SOPRNM(ISOPR)=SOPRNM(I)
         ISOCMP(ISOPR)=ISOCMP(I)
        END IF
       END DO
       CALL GETMEM('PMAP','FREE','INTE',LPMAP,NPMSIZ)
       NSOPR=ISOPR

C Print out the matrix elements:
       NCOL=4
       DO ISOPR=1,NSOPR
        WRITE(6,*)
        WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',SOPRNM(ISOPR),'   COMPONENT:',ISOCMP(ISOPR)
CIFG  should print the origin, but where is it stored (for SO properties)?
        CALL GETMEM('SOPROPR','ALLO','REAL',LSOPRR,NSS**2*NSOPR)
        CALL GETMEM('SOPROPI','ALLO','REAL',LSOPRI,NSS**2*NSOPR)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSOPRR),1)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSOPRI),1)

        CALL SMMAT(PROP,WORK(LSOPRR),NSS,SOPRNM(ISOPR),ISOCMP(ISOPR))
        CALL ZTRNSF(NSS,USOR,USOI,WORK(LSOPRR),WORK(LSOPRI))
        CALL PRCMAT(NSS,WORK(LSOPRR),WORK(LSOPRI))
#ifdef _HDF5_
        IF( SOPRNM(ISOPR)(1:6) .EQ.'ANGMOM') THEN
         call mh5_put_dset_array_real(wfn_sos_angmomr,
     $   WORK(LSOPRR),[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
         call mh5_put_dset_array_real(wfn_sos_angmomi,
     $   WORK(LSOPRI),[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
        ENDIF
#endif
        CALL GETMEM('SOPROPR','FREE','REAL',LSOPRR,NSS**2*NSOPR)
        CALL GETMEM('SOPROPI','FREE','REAL',LSOPRI,NSS**2*NSOPR)
       END DO
       Call CollapseOutput(0,'Matrix elements over SO states')
       WRITE(6,*)

      END IF

 300  CONTINUE

 400  CONTINUE

******************************************************
* printout of special properties
******************************************************

      IF (IPGLOB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(6X,100A1)') ('*',i=1,100)
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,A,34X,A,34X,A)')
     &       '*','  Special properties section  ','*'
        WRITE(6,'(6X,A,98X,A)') '*','*'
        WRITE(6,'(6X,100A1)') ('*',i=1,100)
        WRITE(6,*)
        WRITE(6,*)
      END IF

C Compute transition strengths for spin-orbit states:
      IF(.not.IFSO) GOTO 500
*
* Initial setup for both dipole, quadrupole etc. and exact operator
*
C printing threshold
!     IF(IPGLOB.eq.USUAL) OSTHR=1.0D-8 ! first order
!     IF(IPGLOB.eq.USUAL) OSTHR2=1.0D-12 ! second order (weaker)
!     IF(IPGLOB.gt.USUAL) OSTHR=0.0D0
!     IF(IPGLOB.gt.USUAL) OSTHR2=0.0D0
      OSTHR=1.0D-5
      OSTHR2=1.0D-5
      IF(DIPR) OSTHR = OSTHR_DIPR
      IF(DIPR) WRITE(6,*) ' Dipole threshold changed to ',OSTHR
! Again to avoid total negative transition strengths
      IF(QIPR) OSTHR = OSTHR_QIPR
      IF(QIPR) WRITE(6,*) ' Dipole threshold changed to ',OSTHR,
     &                    ' since quadrupole threshold is given '
      IF(QIPR) OSTHR2 = OSTHR_QIPR
      IF(QIPR) WRITE(6,*) ' Quadrupole threshold changed to ',OSTHR2

      IF(QIALL) WRITE(6,*) ' Will write all quadrupole contributions '
!
!     Reducing the loop over states - good for X-rays
!     At the moment memory is not reduced
!
      IF(REDUCELOOP) THEN
        EX=ENSOR(1)
        L=1
        LD=1
        DO ISO = 2, NSS
           If (ABS(ENSOR(ISO)-EX).gt.1.0D-8) Then
              LD = LD + 1
              EX = ENSOR(ISO)
           Else
              L = L + 1
           End If
           If (LD.gt.LOOPDIVIDE) Exit
        End Do
        IEND = L
        JSTART = L+1
      ELSE
        IEND = NSS
        JSTART = 1
      END IF
!
      IF (IPGLOB.GE.TERSE) THEN
!
!     Initialize arrays for indentifying problematic transitions
!     These stores all dipole oscillator strengths in
!     length and velocity gauge for a later comparison.
!
      CALL GETMEM('DL   ','ALLO','REAL',LDL,NSS**2)
      CALL GETMEM('DV   ','ALLO','REAL',LDV,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDL),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDV),1)
      I_HAVE_DL = 0
      I_HAVE_DV = 0

*Electric-Dipole Electric-Dipole transitions

        IPRDX=0
        IPRDY=0
        IPRDZ=0
        IFANYD=0
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'MLTPL  1') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
          END IF
        END DO

        i_Print=0
        IF(IFANYD.NE.0) THEN
         AFACTOR=32.1299D09

         CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
         CALL GETMEM('DYR','ALLO','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','ALLO','REAL',LDYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYI),1)
         CALL GETMEM('DZR','ALLO','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','ALLO','REAL',LDZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZI),1)

         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXR),NSS,SOPRNM(IPRDX),ISOCMP(IPRDX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYR),NSS,SOPRNM(IPRDY),ISOCMP(IPRDY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYR),WORK(LDYI))
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZR),NSS,SOPRNM(IPRDZ),ISOCMP(IPRDZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZR),WORK(LDZI))
         END IF

         Two3rds=2.0D0/3.0D0
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF).LE.1.0D-8) CYCLE
           IF(EDIFF.GT.0.0D0) THEN
            IJSS=ISS+NSS*(JSS-1)
            If (Do_SK) Then
               T0(1)=DCMPLX(WORK(LDXR-1+IJSS),WORK(LDXI-1+IJSS))
               T0(2)=DCMPLX(WORK(LDYR-1+IJSS),WORK(LDYI-1+IJSS))
               T0(3)=DCMPLX(WORK(LDZR-1+IJSS),WORK(LDZI-1+IJSS))
               TM1=k_vector(1)*T0(1)+k_vector(2)*T0(2)
     &                              +k_vector(3)*T0(3)
               T0(1) = T0(1) - TM1 * k_vector(1)
               T0(2) = T0(2) - TM1 * k_vector(2)
               T0(3) = T0(3) - TM1 * k_vector(3)
               DX2=ABS(DCONJG(T0(1))*T0(1))
               DY2=ABS(DCONJG(T0(2))*T0(2))
               DZ2=ABS(DCONJG(T0(3))*T0(3))
            Else
               DX2=WORK(LDXR-1+IJSS)**2+WORK(LDXI-1+IJSS)**2
               DY2=WORK(LDYR-1+IJSS)**2+WORK(LDYI-1+IJSS)**2
               DZ2=WORK(LDZR-1+IJSS)**2+WORK(LDZI-1+IJSS)**2
            End If
            FX=Two3rds*EDIFF*(DX2)
            FY=Two3rds*EDIFF*(DY2)
            FZ=Two3rds*EDIFF*(DZ2)
            F =FX+FY+FZ
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
! Store dipole oscillator strength
            WORK(LDL-1+IJSS) = F
            IF(ABS(F).GE.OSTHR) THEN
              If (i_Print.eq.0) Then
                 i_Print=1
         Call CollapseOutput(1,
     &                     'Dipole transition strengths (SO states):')
         WRITE(6,'(3X,A)') '----------------------------------------'
         IF(OSTHR.GT.0.0D0) THEN
          WRITE(6,*)'   for osc. strength at least ',OSTHR
          WRITE(6,*)
         END IF
         If (Do_SK) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F8.4,a)')
     &            'Direction of the k-vector: ',
     &             (k_vector(k),k=1,3),' (au)'
            WRITE(6,'(4x,a)')
     &            'The light is assumed to be unpolarized.'
            WRITE(6,*)
         End If
         WRITE(6,*)"        To  From     Osc. strength"//
     &    "   Einstein coefficients Ax, Ay, Az (sec-1) "//
     &    "      Total A (sec-1)  "
         WRITE(6,*)'        -----------------------------------------'//
     &    '--------------------------------------------------'
               End If
             WRITE(6,'(5X,2I5,5X,5ES16.8)') ISS,JSS,F,AX,AY,AZ,A
            END IF
            Call Add_Info('TMS(SO,Len)',F,1,6)
           END IF
          END DO
         END DO

         CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LDYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LDZI,NSS**2)

         If (i_Print.eq.1) THEN
           WRITE(6,*)'        -----------------------------------------'
     &      //'--------------------------------------------------'
           Call CollapseOutput(0,
     &                     'Dipole transition strengths (SO states):')
           WRITE(6,*)
         END IF
         I_HAVE_DL = 1
        END IF

*       Now the same in velocity representation

        IPRDX=0
        IPRDY=0
        IPRDZ=0
        IFANYD=0
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'VELOCITY') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
          END IF
        END DO

        i_Print=0
        IF(IFANYD.NE.0) THEN
         AFACTOR=32.1299D09

         CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
*
         CALL GETMEM('DYR','ALLO','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','ALLO','REAL',LDYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYI),1)
*
         CALL GETMEM('DZR','ALLO','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','ALLO','REAL',LDZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZI),1)

         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXR),NSS,SOPRNM(IPRDX),ISOCMP(IPRDX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYR),NSS,SOPRNM(IPRDY),ISOCMP(IPRDY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYR),WORK(LDYI))
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZR),NSS,SOPRNM(IPRDZ),ISOCMP(IPRDZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZR),WORK(LDZI))
         END IF

         Two3rds=2.0D0/3.0D0
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF (ABS(EDIFF).LE.1.0D-8) CYCLE
           IF(EDIFF.GT.0.0D0) THEN
            IJSS=ISS+NSS*(JSS-1)
            If (Do_SK) Then
               T0(1)=DCMPLX(WORK(LDXR-1+IJSS),WORK(LDXI-1+IJSS))
               T0(2)=DCMPLX(WORK(LDYR-1+IJSS),WORK(LDYI-1+IJSS))
               T0(3)=DCMPLX(WORK(LDZR-1+IJSS),WORK(LDZI-1+IJSS))
               TM1=k_vector(1)*T0(1)+k_vector(2)*T0(2)
     &                              +k_vector(3)*T0(3)
               T0(1) = T0(1) - TM1 * k_vector(1)
               T0(2) = T0(2) - TM1 * k_vector(2)
               T0(3) = T0(3) - TM1 * k_vector(3)
               DX2=ABS(DCONJG(T0(1))*T0(1))
               DY2=ABS(DCONJG(T0(2))*T0(2))
               DZ2=ABS(DCONJG(T0(3))*T0(3))
            Else
               DX2=WORK(LDXR-1+IJSS)**2+WORK(LDXI-1+IJSS)**2
               DY2=WORK(LDYR-1+IJSS)**2+WORK(LDYI-1+IJSS)**2
               DZ2=WORK(LDZR-1+IJSS)**2+WORK(LDZI-1+IJSS)**2
            End If
            FX=Two3rds*(DX2)/EDIFF
            FY=Two3rds*(DY2)/EDIFF
            FZ=Two3rds*(DZ2)/EDIFF
            F =FX+FY+FZ
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
! Store dipole oscillator strength
            WORK(LDV-1+IJSS) = F
            IF(ABS(F).GE.OSTHR) THEN
              If (i_Print.eq.0) Then
                 i_Print=1
         Call CollapseOutput(1,
     &                     'Velocity transition strengths (SO states):')
         WRITE(6,'(3X,A)') '------------------------------------------'
         IF(OSTHR.GT.0.0D0) THEN
          WRITE(6,*)'   for osc. strength at least ',OSTHR
          WRITE(6,*)
         END IF
         If (Do_SK) Then
            WRITE(6,*)
            WRITE(6,'(4x,a,3F8.4,a)')
     &            'Direction of the k-vector: ',
     &             (k_vector(k),k=1,3),' (au)'
            WRITE(6,'(4x,a)')
     &            'The light is assumed to be unpolarized.'
            WRITE(6,*)
         End If
         WRITE(6,*)"        To  From     Osc. strength"//
     &    "   Einstein coefficients Ax, Ay, Az (sec-1) "//
     &    "      Total A (sec-1)  "
         WRITE(6,*)'        -----------------------------------------'//
     &    '--------------------------------------------------'
               END IF
             WRITE(6,'(5X,2I5,5X,5ES16.8)') ISS,JSS,F,AX,AY,AZ,A
            END IF
            Call Add_Info('TMS(SO,Vel)',F,1,6)
           END IF
          END DO
         END DO

         CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LDYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LDZI,NSS**2)

         If (i_Print.eq.1) THEN
           WRITE(6,*)'        -----------------------------------------'
     &     //'--------------------------------------------------'
           Call CollapseOutput(0,
     &                     'Velocity transition strengths (SO states):')
           WRITE(6,*)
         END IF
         I_HAVE_DV = 1
        END IF

!
!      Compare oscillator strengths in length and velocity gauge
!      All differences in oscillator strengths above the tolerance
!      of 0.1 (10 percent) will be printed.
!
       IF(I_HAVE_DL.EQ.1.AND.I_HAVE_DV.EQ.1) THEN
!
! I guess that I have to explain it when I print a warning
!
         WRITE(6,*)
         WRITE(6,*) "--------------------------------------------------"
         WRITE(6,*)
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
         WRITE(6,*)
         WRITE(6,*) "--------------------------------------------------"
!
          I_PRINT_HEADER = 0
          DO I=1,IEND
            DO J=JSTART,NSS
               IJ=I+NSS*(J-1)
               EDIFF=ENSOR(J)-ENSOR(I)
               IF(EDIFF.LT.0.0D0) CYCLE
               COMPARE=0.0D0
             IF(WORK(LDL-1+IJ).GE.OSTHR.AND.WORK(LDV-1+IJ).GE.OSTHR)
     &          THEN
               COMPARE = ABS(1-WORK(LDL-1+IJ)/WORK(LDV-1+IJ))
             ELSE IF(WORK(LDL-1+IJ).GE.OSTHR) THEN
               COMPARE = -1.5D0
             ELSE IF(WORK(LDV-1+IJ).GE.OSTHR) THEN
               COMPARE = -2.5D0
             END IF
             IF(ABS(COMPARE).GE.TOLERANCE) THEN
               I_PRINT_HEADER = I_PRINT_HEADER + 1
               IF(I_PRINT_HEADER.EQ.1) THEN
                 WRITE(6,*)
                 WRITE(6,*) " Problematic transitions have been found"
                 WRITE(6,*) "     From   To      Difference (%)  "//
     &                      "Osc. st. (len.) Osc. st. (vel.)"
                 WRITE(6,*) "     -------------------------------"//
     &                      "-------------------------------"
                 WRITE(6,*)
               END IF
               IF (COMPARE.GE.0.0D0) THEN
                 WRITE(6,33) I,J,COMPARE*100D0,
     &                    WORK(LDL-1+IJ),WORK(LDV-1+IJ)
               ELSE IF (COMPARE.GE.2.0D0) THEN
                 WRITE(6,36) I,J,WORK(LDL-1+IJ),"below threshold"
               ELSE
                 WRITE(6,37) I,J,"below threshold",WORK(LDV-1+IJ)
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
            WRITE(6,*) "     -------------------------------"//
     &                 "-------------------------------"
            WRITE(6,*)
            WRITE(6,*) "Number of problematic transitions = ",
     &                  I_PRINT_HEADER
            WRITE(6,*)
          END IF
        END IF
33    FORMAT (5X,2(1X,I4),5X,5(1X,ES15.8))
36    FORMAT (5X,2(1X,I4),6X,15('-'),1X,ES15.8,1X,A15)
37    FORMAT (5X,2(1X,I4),6X,15('-'),1X,A15,1X,ES15.8)
*
* Free the memory
*
      CALL GETMEM('DL   ','FREE','REAL',LDL,NSTATE**2)
      CALL GETMEM('DV   ','FREE','REAL',LDV,NSTATE**2)
!
! We will first allocate a matrix for the total of the second order wave vector
!
        CALL GETMEM('TOT2K','ALLO','REAL',LTOT2K,NSS**2)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LTOT2K),1)
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
        IPRDX=0
        IPRDY=0
        IPRDZ=0
! Spin-Magnetic-Dipole ---- notice the S
        IPRSX=0
        IPRSY=0
        IPRSZ=0

        IFANYD=0
        IFANYS=0
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'ANGMOM  ') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
          ELSE IF(SOPRNM(ISOPR).EQ.'AMFI') THEN
           IFANYS=1
           IF(ISOCMP(ISOPR).EQ.1) IPRSX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRSY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRSZ=ISOPR
          END IF
        END DO

        IF(IFANYD.NE.0.OR.IFANYS.NE.0) THEN
!
! Only print the part calculated
!
         IF(QIALL) THEN
         IF(IFANYD.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Magnetic-dipole - magnetic-dipole and '//
     &                  'spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------'//
     &                  '--------------------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYD.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(1,
     &                  'Magnetic-dipole - magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '----------------------------------'//
     &                  '---------------------------------'
         ELSE IF(IFANYD.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(1,
     &                  'Spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,'(3X,A)')
     &                  '--------------------------------------------'//
     &                  '---------------------------------'
         END IF
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,*)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,*)"        To  From     Osc. strength"
         WRITE(6,*)'        --------------------------'
         END IF
! Magnetic-Dipole
         CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
         CALL GETMEM('DYR','ALLO','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','ALLO','REAL',LDYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYI),1)
         CALL GETMEM('DZR','ALLO','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','ALLO','REAL',LDZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZI),1)

! Spin-Magnetic-Dipole
         CALL GETMEM('SXR','ALLO','REAL',LSXR,NSS**2)
         CALL GETMEM('SXI','ALLO','REAL',LSXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSXI),1)
         CALL GETMEM('SYR','ALLO','REAL',LSYR,NSS**2)
         CALL GETMEM('SYI','ALLO','REAL',LSYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSYI),1)
         CALL GETMEM('DSR','ALLO','REAL',LSZR,NSS**2)
         CALL GETMEM('DSI','ALLO','REAL',LSZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSZI),1)

! Magnetic-Dipole
         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXR),NSS,SOPRNM(IPRDX),ISOCMP(IPRDX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYR),NSS,SOPRNM(IPRDY),ISOCMP(IPRDY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYR),WORK(LDYI))
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZR),NSS,SOPRNM(IPRDZ),ISOCMP(IPRDZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZR),WORK(LDZI))
         END IF

! Spin-Magnetic-Dipole
         IF(IPRSX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSXR),NSS,SOPRNM(IPRSX),ISOCMP(IPRSX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSXR),WORK(LSXI))
         END IF
         IF(IPRSY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSYR),NSS,SOPRNM(IPRSY),ISOCMP(IPRSY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSYR),WORK(LSYI))
         END IF
         IF(IPRSZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSZR),NSS,SOPRNM(IPRSZ),ISOCMP(IPRSZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSZR),WORK(LSZI))
         END IF

         ONEOVER6C2=1.0D0/(6.0D0*CONST_C_IN_AU_**2)
         g = 2.0D0 ! Well I assume that the g-factor for the electron is 2
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF(EDIFF.GT.0.0D0) THEN
            IJSS=ISS+NSS*(JSS-1)

            DX2=(WORK(LDXR-1+IJSS)+g*WORK(LSXR-1+IJSS))**2
     &         +(WORK(LDXI-1+IJSS)+g*WORK(LSXI-1+IJSS))**2
            DY2=(WORK(LDYR-1+IJSS)+g*WORK(LSYR-1+IJSS))**2
     &         +(WORK(LDYI-1+IJSS)+g*WORK(LSYI-1+IJSS))**2
            DZ2=(WORK(LDZR-1+IJSS)+g*WORK(LSZR-1+IJSS))**2
     &         +(WORK(LDZI-1+IJSS)+g*WORK(LSZI-1+IJSS))**2

            F = (DX2 + DY2 + DZ2)*EDIFF*ONEOVER6C2
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F
            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,'(5X,2I5,5X,ES16.8)') ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

! Magnetic-Dipole
         CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LDYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LDZI,NSS**2)

! Spin-Magnetic-Dipole
         CALL GETMEM('DXR','FREE','REAL',LSXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LSXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LSYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LSYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LSZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LSZI,NSS**2)

       IF(QIALL) THEN
         WRITE(6,*)'        --------------------------'
         IF(IFANYD.NE.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Magnetic-dipole - magnetic-dipole and '//
     &                  'spin-magnetic-dipole - spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         ELSE IF(IFANYD.NE.0.AND.IFANYS.EQ.0) THEN
          Call CollapseOutput(0,
     &                  'Magnetic-dipole - magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         ELSE IF(IFANYD.EQ.0.AND.IFANYS.NE.0) THEN
          Call CollapseOutput(0,
     &                  'Spin-magnetic-dipole - Spin-magnetic-dipole '//
     &                  'transition strengths (SO states):')
          WRITE(6,*)
         END IF
        END IF
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
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'MLTPL  2') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDXX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDXY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDXZ=ISOPR
           IF(ISOCMP(ISOPR).EQ.4) IPRDYY=ISOPR
           IF(ISOCMP(ISOPR).EQ.5) IPRDYZ=ISOPR
           IF(ISOCMP(ISOPR).EQ.6) IPRDZZ=ISOPR
          END IF
        END DO

        IF(IFANYD.NE.0) THEN
        IF(QIALL) THEN
         Call CollapseOutput(1,
     &                 'Quadrupole transition strengths (SO states):')
         WRITE(6,'(3X,A)')
     &                 '--------------------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,*)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,*)"        To  From     Osc. strength"
         WRITE(6,*)'        ----------------------------'
         END IF

         CALL GETMEM('DXXR','ALLO','REAL',LDXXR,NSS**2)
         CALL GETMEM('DXXI','ALLO','REAL',LDXXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXI),1)
         CALL GETMEM('DXYR','ALLO','REAL',LDXYR,NSS**2)
         CALL GETMEM('DXYI','ALLO','REAL',LDXYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXYI),1)
         CALL GETMEM('DXZR','ALLO','REAL',LDXZR,NSS**2)
         CALL GETMEM('DXZI','ALLO','REAL',LDXZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXZI),1)
         CALL GETMEM('DYYR','ALLO','REAL',LDYYR,NSS**2)
         CALL GETMEM('DYYI','ALLO','REAL',LDYYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYI),1)
         CALL GETMEM('DYZR','ALLO','REAL',LDYZR,NSS**2)
         CALL GETMEM('DYZI','ALLO','REAL',LDYZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYZI),1)
         CALL GETMEM('DZZR','ALLO','REAL',LDZZR,NSS**2)
         CALL GETMEM('DZZI','ALLO','REAL',LDZZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZI),1)

         IF(IPRDXX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXXR),NSS,SOPRNM(IPRDXX),ISOCMP(IPRDXX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXXR),WORK(LDXXI))
         END IF
         IF(IPRDXY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXYR),NSS,SOPRNM(IPRDXY),ISOCMP(IPRDXY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXYR),WORK(LDXYI))
         END IF
         IF(IPRDXZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXZR),NSS,SOPRNM(IPRDXZ),ISOCMP(IPRDXZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXZR),WORK(LDXZI))
         END IF
         IF(IPRDYY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYYR),NSS,SOPRNM(IPRDYY),ISOCMP(IPRDYY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYYR),WORK(LDYYI))
         END IF
         IF(IPRDYZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYZR),NSS,SOPRNM(IPRDYZ),ISOCMP(IPRDYZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYZR),WORK(LDYZI))
         END IF
         IF(IPRDZZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZZR),NSS,SOPRNM(IPRDZZ),ISOCMP(IPRDZZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZZR),WORK(LDZZI))
         END IF

         ONEOVER10C=1.0D0/(10.0D0*CONST_C_IN_AU_**2)
         ONEOVER30C=ONEOVER10C/3.0D0

         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
! D should be purely real since D is a real symmetric matrix
!
            EDIFF3=EDIFF**3
            IJSS=ISS+NSS*(JSS-1)

            DXX2=WORK(LDXXR-1+IJSS)**2+WORK(LDXXI-1+IJSS)**2
            DYY2=WORK(LDYYR-1+IJSS)**2+WORK(LDYYI-1+IJSS)**2
            DZZ2=WORK(LDZZR-1+IJSS)**2+WORK(LDZZI-1+IJSS)**2
            FXX=ONEOVER30C*EDIFF3*(DXX2)
            FYY=ONEOVER30C*EDIFF3*(DYY2)
            FZZ=ONEOVER30C*EDIFF3*(DZZ2)

            DXY2=WORK(LDXYR-1+IJSS)**2+WORK(LDXYI-1+IJSS)**2
            DXZ2=WORK(LDXZR-1+IJSS)**2+WORK(LDXZI-1+IJSS)**2
            DYZ2=WORK(LDYZR-1+IJSS)**2+WORK(LDYZI-1+IJSS)**2
            FXY=ONEOVER10C*EDIFF3*(DXY2)
            FXZ=ONEOVER10C*EDIFF3*(DXZ2)
            FYZ=ONEOVER10C*EDIFF3*(DYZ2)

            DXXDYY=WORK(LDXXR-1+IJSS)*WORK(LDYYR-1+IJSS)
     &            +WORK(LDXXI-1+IJSS)*WORK(LDYYI-1+IJSS)
            DXXDZZ=WORK(LDXXR-1+IJSS)*WORK(LDZZR-1+IJSS)
     &            +WORK(LDXXI-1+IJSS)*WORK(LDZZI-1+IJSS)
            DYYDZZ=WORK(LDYYR-1+IJSS)*WORK(LDZZR-1+IJSS)
     &            +WORK(LDYYI-1+IJSS)*WORK(LDZZI-1+IJSS)
            FXXFYY=-ONEOVER30C*EDIFF3*(DXXDYY)
            FXXFZZ=-ONEOVER30C*EDIFF3*(DXXDZZ)
            FYYFZZ=-ONEOVER30C*EDIFF3*(DYYDZZ)

            F =FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,'(5X,2I5,5X,ES16.8)') ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

         CALL GETMEM('DXXR','FREE','REAL',LDXXR,NSS**2)
         CALL GETMEM('DXXI','FREE','REAL',LDXXI,NSS**2)
         CALL GETMEM('DXYR','FREE','REAL',LDXYR,NSS**2)
         CALL GETMEM('DXYI','FREE','REAL',LDXYI,NSS**2)
         CALL GETMEM('DXZR','FREE','REAL',LDXZR,NSS**2)
         CALL GETMEM('DXZI','FREE','REAL',LDXZI,NSS**2)
         CALL GETMEM('DYYR','FREE','REAL',LDYYR,NSS**2)
         CALL GETMEM('DYYI','FREE','REAL',LDYYI,NSS**2)
         CALL GETMEM('DYZR','FREE','REAL',LDYZR,NSS**2)
         CALL GETMEM('DYZI','FREE','REAL',LDYZI,NSS**2)
         CALL GETMEM('DZZR','FREE','REAL',LDZZR,NSS**2)
         CALL GETMEM('DZZI','FREE','REAL',LDZZI,NSS**2)

        IF(QIALL) THEN
         WRITE(6,*)'        ----------------------------'
         Call CollapseOutput(0,
     &                 'Quadrupole transition strengths (SO states):')
         WRITE(6,*)
        END IF
        SECORD(2) = 1
        END IF

*Electric-Dipole Electric-Octupole transitions

! Octupole
! This is a real symmetric rank 3 tensor so only 10 and not 27 is needed
! The order which comes in
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
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'MLTPL  1') THEN
           IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
          ELSE IF(SOPRNM(ISOPR).EQ.'MLTPL  3') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDXXX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDXXY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDXXZ=ISOPR
           IF(ISOCMP(ISOPR).EQ.4) IPRDYYX=ISOPR ! Changed from XYY
           IF(ISOCMP(ISOPR).EQ.5) IPRDXYZ=ISOPR
           IF(ISOCMP(ISOPR).EQ.6) IPRDZZX=ISOPR ! Changed from XZZ
           IF(ISOCMP(ISOPR).EQ.7) IPRDYYY=ISOPR
           IF(ISOCMP(ISOPR).EQ.8) IPRDYYZ=ISOPR
           IF(ISOCMP(ISOPR).EQ.9) IPRDZZY=ISOPR ! Changed from YZZ
           IF(ISOCMP(ISOPR).EQ.10) IPRDZZZ=ISOPR

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
         Call CollapseOutput(1,
     &                     'Electric-dipole - electric-octupole '//
     &                     'transition strengths (SO states):')
         WRITE(6,'(3X,A)') '------------------------------------'//
     &                     '---------------------------------'
         IF(OSTHR2.GT.0.0D0) THEN
          WRITE(6,*)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,*)"        To  From     Osc. strength"
         WRITE(6,*)'        ----------------------------'
         END IF
! Octupole
         CALL GETMEM('DXXXR','ALLO','REAL',LDXXXR,NSS**2)
         CALL GETMEM('DXXXI','ALLO','REAL',LDXXXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXXI),1)
         CALL GETMEM('DXXYR','ALLO','REAL',LDXXYR,NSS**2)
         CALL GETMEM('DXXYI','ALLO','REAL',LDXXYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXYI),1)
         CALL GETMEM('DXXZR','ALLO','REAL',LDXXZR,NSS**2)
         CALL GETMEM('DXXZI','ALLO','REAL',LDXXZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXXZI),1)

         CALL GETMEM('DYYXR','ALLO','REAL',LDYYXR,NSS**2)
         CALL GETMEM('DYYXI','ALLO','REAL',LDYYXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYXI),1)
         CALL GETMEM('DYYYR','ALLO','REAL',LDYYYR,NSS**2)
         CALL GETMEM('DYYYI','ALLO','REAL',LDYYYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYYI),1)
         CALL GETMEM('DYYZR','ALLO','REAL',LDYYZR,NSS**2)
         CALL GETMEM('DYYZI','ALLO','REAL',LDYYZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYYZI),1)

         CALL GETMEM('DZZXR','ALLO','REAL',LDZZXR,NSS**2)
         CALL GETMEM('DZZXI','ALLO','REAL',LDZZXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZXI),1)
         CALL GETMEM('DZZYR','ALLO','REAL',LDZZYR,NSS**2)
         CALL GETMEM('DZZYI','ALLO','REAL',LDZZYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZYI),1)
         CALL GETMEM('DZZZR','ALLO','REAL',LDZZZR,NSS**2)
         CALL GETMEM('DZZZI','ALLO','REAL',LDZZZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZZZI),1)
! Dipole
         CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
         CALL GETMEM('DYR','ALLO','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','ALLO','REAL',LDYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYI),1)
         CALL GETMEM('DZR','ALLO','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','ALLO','REAL',LDZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZI),1)
! Octupole
         IF(IPRDXXX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXXXR),NSS,SOPRNM(IPRDXXX),
     &               ISOCMP(IPRDXXX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXXXR),WORK(LDXXXI))
         END IF
         IF(IPRDXXY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXXYR),NSS,SOPRNM(IPRDXXY),
     &               ISOCMP(IPRDXXY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXXYR),WORK(LDXXYI))
         END IF
         IF(IPRDXXZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXXZR),NSS,SOPRNM(IPRDXXZ),
     &               ISOCMP(IPRDXXZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXXZR),WORK(LDXXZI))
         END IF

         IF(IPRDYYX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYYXR),NSS,SOPRNM(IPRDYYX),
     &               ISOCMP(IPRDYYX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYYXR),WORK(LDYYXI))
         END IF
         IF(IPRDYYY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYYYR),NSS,SOPRNM(IPRDYYY),
     &               ISOCMP(IPRDYYY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYYYR),WORK(LDYYYI))
         END IF
         IF(IPRDYYZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYYZR),NSS,SOPRNM(IPRDYYZ),
     &               ISOCMP(IPRDYYZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYYZR),WORK(LDYYZI))
         END IF

         IF(IPRDZZX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZZXR),NSS,SOPRNM(IPRDZZX),
     &               ISOCMP(IPRDZZX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZZXR),WORK(LDZZXI))
         END IF
         IF(IPRDZZY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZZYR),NSS,SOPRNM(IPRDZZY),
     &               ISOCMP(IPRDZZY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZZYR),WORK(LDZZYI))
         END IF
         IF(IPRDZZZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZZZR),NSS,SOPRNM(IPRDZZZ),
     &               ISOCMP(IPRDZZZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZZZR),WORK(LDZZZI))
         END IF
! Dipole
         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXR),NSS,SOPRNM(IPRDX),
     &               ISOCMP(IPRDX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYR),NSS,SOPRNM(IPRDY),
     &               ISOCMP(IPRDY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYR),WORK(LDYI))
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZR),NSS,SOPRNM(IPRDZ),
     &               ISOCMP(IPRDZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZR),WORK(LDZI))
         END IF

         TWOOVERM45C=-2.0D0/(45.0D0*CONST_C_IN_AU_**2)
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF3=EDIFF**3
            IJSS=ISS+NSS*(JSS-1)

            DXXXDX=WORK(LDXXXR-1+IJSS)*WORK(LDXR-1+IJSS)
     &            +WORK(LDXXXI-1+IJSS)*WORK(LDXI-1+IJSS)
            DYYXDX=WORK(LDYYXR-1+IJSS)*WORK(LDXR-1+IJSS)
     &            +WORK(LDYYXI-1+IJSS)*WORK(LDXI-1+IJSS)
            DZZXDX=WORK(LDZZXR-1+IJSS)*WORK(LDXR-1+IJSS)
     &            +WORK(LDZZXI-1+IJSS)*WORK(LDXI-1+IJSS)
            FXXX=TWOOVERM45C*EDIFF3*(DXXXDX)
            FYYX=TWOOVERM45C*EDIFF3*(DYYXDX)
            FZZX=TWOOVERM45C*EDIFF3*(DZZXDX)

            DXXYDY=WORK(LDXXYR-1+IJSS)*WORK(LDYR-1+IJSS)
     &            +WORK(LDXXYI-1+IJSS)*WORK(LDYI-1+IJSS)
            DYYYDY=WORK(LDYYYR-1+IJSS)*WORK(LDYR-1+IJSS)
     &            +WORK(LDYYYI-1+IJSS)*WORK(LDYI-1+IJSS)
            DZZYDY=WORK(LDZZYR-1+IJSS)*WORK(LDYR-1+IJSS)
     &            +WORK(LDZZYI-1+IJSS)*WORK(LDYI-1+IJSS)
            FXXY=TWOOVERM45C*EDIFF3*(DXXYDY)
            FYYY=TWOOVERM45C*EDIFF3*(DYYYDY)
            FZZY=TWOOVERM45C*EDIFF3*(DZZYDY)

            DXXZDZ=WORK(LDXXZR-1+IJSS)*WORK(LDZR-1+IJSS)
     &            +WORK(LDXXZI-1+IJSS)*WORK(LDZI-1+IJSS)
            DYYZDZ=WORK(LDYYZR-1+IJSS)*WORK(LDZR-1+IJSS)
     &            +WORK(LDYYZI-1+IJSS)*WORK(LDZI-1+IJSS)
            DZZZDZ=WORK(LDZZZR-1+IJSS)*WORK(LDZR-1+IJSS)
     &            +WORK(LDZZZI-1+IJSS)*WORK(LDZI-1+IJSS)
            FXXZ=TWOOVERM45C*EDIFF3*(DXXZDZ)
            FYYZ=TWOOVERM45C*EDIFF3*(DYYZDZ)
            FZZZ=TWOOVERM45C*EDIFF3*(DZZZDZ)

            F =FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,'(5X,2I5,5X,ES16.8)') ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

         CALL GETMEM('DXXXR','FREE','REAL',LDXXXR,NSS**2)
         CALL GETMEM('DXXXI','FREE','REAL',LDXXXI,NSS**2)
         CALL GETMEM('DXXYR','FREE','REAL',LDXXYR,NSS**2)
         CALL GETMEM('DXXYI','FREE','REAL',LDXXYI,NSS**2)
         CALL GETMEM('DXXZR','FREE','REAL',LDXXZR,NSS**2)
         CALL GETMEM('DXXZI','FREE','REAL',LDXXZI,NSS**2)

         CALL GETMEM('DYYXR','FREE','REAL',LDYYXR,NSS**2)
         CALL GETMEM('DYYXI','FREE','REAL',LDYYXI,NSS**2)
         CALL GETMEM('DYYYR','FREE','REAL',LDYYYR,NSS**2)
         CALL GETMEM('DYYYI','FREE','REAL',LDYYYI,NSS**2)
         CALL GETMEM('DYYZR','FREE','REAL',LDYYZR,NSS**2)
         CALL GETMEM('DYYZI','FREE','REAL',LDYYZI,NSS**2)

         CALL GETMEM('DZZXR','FREE','REAL',LDZZXR,NSS**2)
         CALL GETMEM('DZZXI','FREE','REAL',LDZZXI,NSS**2)
         CALL GETMEM('DZZYR','FREE','REAL',LDZZYR,NSS**2)
         CALL GETMEM('DZZYI','FREE','REAL',LDZZYI,NSS**2)
         CALL GETMEM('DZZZR','FREE','REAL',LDZZZR,NSS**2)
         CALL GETMEM('DZZZI','FREE','REAL',LDZZZI,NSS**2)

         CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LDYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LDZI,NSS**2)

        IF(QIALL) THEN
         WRITE(6,*)'        ----------------------------'
         Call CollapseOutput(0,
     &                     'Electric-dipole - electric-octupole '//
     &                     'transition strengths (SO states):')
         WRITE(6,*)
        END IF
        SECORD(3) = 1
        END IF

*Electric-Dipole - Magnetic-Quadrupole transitions and
*Electric-Dipole - Spin-Magnetic-Quadrupole transitions
!
! Again I will just include the spin-term so both terms are calculated
! (Can also be done separately)
! DM + DMs
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
! Spin-Magnetic-Quadrupole
        IPRSXX=0
        IPRSXY=0
        IPRSXZ=0

        IPRSYX=0
        IPRSYY=0
        IPRSYZ=0

        IPRSZX=0
        IPRSZY=0
        IPRSZZ=0
! Electric-Dipole
        IPRDX=0
        IPRDY=0
        IPRDZ=0
! AMFI
! Spin-Magnetic-Quadrupole = M^s_ab = r_b * s_a
! Will create this integral as the product of
! the electric-dipole and amfi integral (of course in the length representation)
! Will still keep SMQ as keyword since this term could otherwise
! be calculated with other inputs!!!
        IPRSX=0
        IPRSY=0
        IPRSZ=0

        IFANYD=0
        IFANYS=0
        DO ISOPR=1,NSOPR
          IF(SOPRNM(ISOPR).EQ.'MLTPL  1') THEN
           IF(ISOCMP(ISOPR).EQ.1) IPRDX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDZ=ISOPR
          ELSE IF(SOPRNM(ISOPR).EQ.'AMFI') THEN
           IF(ISOCMP(ISOPR).EQ.1) IPRSX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRSY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRSZ=ISOPR
          ELSE IF(SOPRNM(ISOPR).EQ.'OMQ') THEN
           IFANYD=1
           IF(ISOCMP(ISOPR).EQ.1) IPRDXX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRDXY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRDXZ=ISOPR

           IF(ISOCMP(ISOPR).EQ.4) IPRDYX=ISOPR
           IF(ISOCMP(ISOPR).EQ.5) IPRDYY=ISOPR
           IF(ISOCMP(ISOPR).EQ.6) IPRDYZ=ISOPR

           IF(ISOCMP(ISOPR).EQ.7) IPRDZX=ISOPR
           IF(ISOCMP(ISOPR).EQ.8) IPRDZY=ISOPR
           IF(ISOCMP(ISOPR).EQ.9) IPRDZZ=ISOPR

          ELSE IF(SOPRNM(ISOPR).EQ.'SMQ') THEN
           IFANYS=1
           IF(ISOCMP(ISOPR).EQ.1) IPRSXX=ISOPR
           IF(ISOCMP(ISOPR).EQ.2) IPRSXY=ISOPR
           IF(ISOCMP(ISOPR).EQ.3) IPRSXZ=ISOPR

           IF(ISOCMP(ISOPR).EQ.4) IPRSYX=ISOPR
           IF(ISOCMP(ISOPR).EQ.5) IPRSYY=ISOPR
           IF(ISOCMP(ISOPR).EQ.6) IPRSYZ=ISOPR

           IF(ISOCMP(ISOPR).EQ.7) IPRSZX=ISOPR
           IF(ISOCMP(ISOPR).EQ.8) IPRSZY=ISOPR
           IF(ISOCMP(ISOPR).EQ.9) IPRSZZ=ISOPR

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

        IF(((IPRSYZ.GT.0.OR.IPRSZY.GT.0)
     &   .AND.(IPRDX.LE.0.OR.IPRSX.LE.0))) THEN
         WRITE(6,*) ' Remember to include Dipole, AMFI and Quadrupole'
         CALL ABEND()
        END IF
        IF(((IPRSZX.GT.0.OR.IPRSXZ.GT.0)
     &   .AND.(IPRDY.LE.0.OR.IPRSY.LE.0))) THEN
         WRITE(6,*) ' Remember to include Dipole, AMFI and Quadrupole'
         CALL ABEND()
        END IF
        IF(((IPRSXY.GT.0.OR.IPRSYX.GT.0)
     &   .AND.(IPRDZ.LE.0.OR.IPRSZ.LE.0))) THEN
         WRITE(6,*) ' Remember to include Dipole, AMFI and Quadrupole'
         CALL ABEND()
        END IF

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
          WRITE(6,*)'   for osc. strength at least ',OSTHR2
          WRITE(6,*)
         END IF
         WRITE(6,*)"        To  From     Osc. strength"
         WRITE(6,*)'        ----------------------------'
         END IF
! Magnetic-Quadrupole
         CALL GETMEM('DZXR','ALLO','REAL',LDZXR,NSS**2)
         CALL GETMEM('DZXI','ALLO','REAL',LDZXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZXI),1)
         CALL GETMEM('DXZR','ALLO','REAL',LDXZR,NSS**2)
         CALL GETMEM('DXZI','ALLO','REAL',LDXZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXZI),1)

         CALL GETMEM('DXYR','ALLO','REAL',LDXYR,NSS**2)
         CALL GETMEM('DXYI','ALLO','REAL',LDXYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXYI),1)
         CALL GETMEM('DYXR','ALLO','REAL',LDYXR,NSS**2)
         CALL GETMEM('DYXI','ALLO','REAL',LDYXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYXI),1)

         CALL GETMEM('DYZR','ALLO','REAL',LDYZR,NSS**2)
         CALL GETMEM('DYZI','ALLO','REAL',LDYZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYZI),1)
         CALL GETMEM('DZYR','ALLO','REAL',LDZYR,NSS**2)
         CALL GETMEM('DZYI','ALLO','REAL',LDZYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZYI),1)
! Spin-Magnetic-Quadrupole
! Only AMFI integrals will be needed
         CALL GETMEM('SXR','ALLO','REAL',LSXR,NSS**2)
         CALL GETMEM('SXI','ALLO','REAL',LSXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSXI),1)
         CALL GETMEM('SYR','ALLO','REAL',LSYR,NSS**2)
         CALL GETMEM('SYI','ALLO','REAL',LSYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSYI),1)
         CALL GETMEM('SZR','ALLO','REAL',LSZR,NSS**2)
         CALL GETMEM('SZI','ALLO','REAL',LSZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSZI),1)
! Electric-Dipole
         CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
         CALL GETMEM('DYR','ALLO','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','ALLO','REAL',LDYI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDYI),1)
         CALL GETMEM('DZR','ALLO','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','ALLO','REAL',LDZI,NSS**2)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZR),1)
         CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDZI),1)
! Magnetic-Quadrupole
         IF(IPRDXY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXYR),NSS,SOPRNM(IPRDXY),
     &               ISOCMP(IPRDXY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXYR),WORK(LDXYI))
         END IF
         IF(IPRDYX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYXR),NSS,SOPRNM(IPRDYX),
     &               ISOCMP(IPRDYX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYXR),WORK(LDYXI))
         END IF

         IF(IPRDXZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXZR),NSS,SOPRNM(IPRDXZ),
     &               ISOCMP(IPRDXZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXZR),WORK(LDXZI))
         END IF
         IF(IPRDZX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZXR),NSS,SOPRNM(IPRDZX),
     &               ISOCMP(IPRDZX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZXR),WORK(LDZXI))
         END IF

         IF(IPRDYZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYZR),NSS,SOPRNM(IPRDYZ),
     &               ISOCMP(IPRDYZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYZR),WORK(LDYZI))
         END IF
         IF(IPRDZY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZYR),NSS,SOPRNM(IPRDZY),
     &               ISOCMP(IPRDZY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZYR),WORK(LDZYI))
         END IF
! Spin-Magnetic-Quadrupole
! Again only AMFI integrals
! And only check first index in M^s
         IF(IPRSXX.GT.0.OR.IPRSXY.GT.0.OR.IPRSXZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSXR),NSS,SOPRNM(IPRSX),
     &               ISOCMP(IPRSX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSXR),WORK(LSXI))
         END IF
         IF(IPRSYX.GT.0.OR.IPRSYY.GT.0.OR.IPRSYZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSYR),NSS,SOPRNM(IPRSY),
     &               ISOCMP(IPRSY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSYR),WORK(LSYI))
         END IF
         IF(IPRSZX.GT.0.OR.IPRSZY.GT.0.OR.IPRSZZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LSZR),NSS,SOPRNM(IPRSZ),
     &               ISOCMP(IPRSZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LSZR),WORK(LSZI))
         END IF
! Electric-Dipole
         IF(IPRDX.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDXR),NSS,SOPRNM(IPRDX),
     &               ISOCMP(IPRDX))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
         END IF
         IF(IPRDY.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDYR),NSS,SOPRNM(IPRDY),
     &               ISOCMP(IPRDY))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDYR),WORK(LDYI))
         END IF
         IF(IPRDZ.GT.0) THEN
          CALL SMMAT(PROP,WORK(LDZR),NSS,SOPRNM(IPRDZ),
     &               ISOCMP(IPRDZ))
          CALL ZTRNSF(NSS,USOR,USOI,WORK(LDZR),WORK(LDZI))
         END IF

         ONEOVER9C2=1.0D0/(9.0D0*CONST_C_IN_AU_**2)
         g = 2
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            EDIFF2=EDIFF**2
            IJSS=ISS+NSS*(JSS-1)
!
! Since the Spin-Magnetic-Quadrupole is made from the multiplication of two complex integrals we have
! M^s = (a+ib)(c+id) = ac-bd + i(ad+bc) hence the long expressions below
!
            DYXDZ=((WORK(LDYXR-1+IJSS) ! Magnetic-Quadrupole
     &         +g*(WORK(LSYR-1+IJSS)*WORK(LDXR-1+IJSS)
     &            -WORK(LSYI-1+IJSS)*WORK(LDXI-1+IJSS))) ! Spin-Magnetic-Quadrupole
     &           *WORK(LDZR-1+IJSS)) ! Electric-Dipole
     &           +((WORK(LDYXI-1+IJSS)
     &         +g*(WORK(LSYR-1+IJSS)*WORK(LDXI-1+IJSS)
     &            -WORK(LSYI-1+IJSS)*WORK(LDXR-1+IJSS)))
     &           *WORK(LDZI-1+IJSS))
            DXYDZ=((WORK(LDXYR-1+IJSS)
     &         +g*(WORK(LSXR-1+IJSS)*WORK(LDYR-1+IJSS)
     &            -WORK(LSXI-1+IJSS)*WORK(LDYI-1+IJSS)))
     &           *WORK(LDZR-1+IJSS))
     &           +((WORK(LDXYI-1+IJSS)
     &         +g*(WORK(LSXR-1+IJSS)*WORK(LDYI-1+IJSS)
     &            -WORK(LSXI-1+IJSS)*WORK(LDYR-1+IJSS)))
     &           *WORK(LDZI-1+IJSS))
            FXY=ONEOVER9C2*EDIFF2*(DXYDZ)
            FYX=-ONEOVER9C2*EDIFF2*(DYXDZ)

            DZXDY=((WORK(LDZXR-1+IJSS)
     &         +g*(WORK(LSZR-1+IJSS)*WORK(LDXR-1+IJSS)
     &            -WORK(LSZI-1+IJSS)*WORK(LDXI-1+IJSS)))
     &           *WORK(LDYR-1+IJSS))
     &           +((WORK(LDZXI-1+IJSS)
     &         +g*(WORK(LSZR-1+IJSS)*WORK(LDXI-1+IJSS)
     &            -WORK(LSZI-1+IJSS)*WORK(LDXR-1+IJSS)))
     &           *WORK(LDYI-1+IJSS))
            DXZDY=((WORK(LDXZR-1+IJSS)
     &         +g*(WORK(LSXR-1+IJSS)*WORK(LDZR-1+IJSS)
     &            -WORK(LSXI-1+IJSS)*WORK(LDZI-1+IJSS)))
     &           *WORK(LDYR-1+IJSS))
     &           +((WORK(LDXZI-1+IJSS)
     &         +g*(WORK(LSXR-1+IJSS)*WORK(LDZI-1+IJSS)
     &            -WORK(LSXI-1+IJSS)*WORK(LDZR-1+IJSS)))
     &           *WORK(LDYI-1+IJSS))
            FZX=ONEOVER9C2*EDIFF2*(DZXDY)
            FXZ=-ONEOVER9C2*EDIFF2*(DXZDY)

            DYZDX=((WORK(LDYZR-1+IJSS)
     &         +g*(WORK(LSYR-1+IJSS)*WORK(LDZR-1+IJSS)
     &            -WORK(LSYI-1+IJSS)*WORK(LDZI-1+IJSS)))
     &           *WORK(LDXR-1+IJSS))
     &           +((WORK(LDYZI-1+IJSS)
     &         +g*(WORK(LSYR-1+IJSS)*WORK(LDZI-1+IJSS)
     &            -WORK(LSYI-1+IJSS)*WORK(LDZR-1+IJSS)))
     &           *WORK(LDXI-1+IJSS))
            DZYDX=((WORK(LDZYR-1+IJSS)
     &         +g*(WORK(LSZR-1+IJSS)*WORK(LDYR-1+IJSS)
     &            -WORK(LSZI-1+IJSS)*WORK(LDYI-1+IJSS)))
     &           *WORK(LDXR-1+IJSS))
     &           +((WORK(LDZYI-1+IJSS)
     &         +g*(WORK(LSZR-1+IJSS)*WORK(LDYI-1+IJSS)
     &            -WORK(LSZI-1+IJSS)*WORK(LDYR-1+IJSS)))
     &           *WORK(LDXI-1+IJSS))
            FYZ=ONEOVER9C2*EDIFF2*(DYZDX)
            FZY=-ONEOVER9C2*EDIFF2*(DZYDX)

            F =FYX+FXY+FZX+FXZ+FYZ+FZY
            F = -F
! Add it to the total
            WORK(LTOT2K-1+IJSS) = WORK(LTOT2K-1+IJSS) + F

            IF(ABS(F).GE.OSTHR2) THEN
             IF(QIALL) WRITE(6,'(5X,2I5,5X,ES16.8)') ISS,JSS,F
            END IF
           END IF
          END DO
         END DO

! Magnetic-Quadrupole
         CALL GETMEM('DXYR','FREE','REAL',LDXYR,NSS**2)
         CALL GETMEM('DXYI','FREE','REAL',LDXYI,NSS**2)
         CALL GETMEM('DXZR','FREE','REAL',LDXZR,NSS**2)
         CALL GETMEM('DXZI','FREE','REAL',LDXZI,NSS**2)

         CALL GETMEM('DYXR','FREE','REAL',LDYXR,NSS**2)
         CALL GETMEM('DYXI','FREE','REAL',LDYXI,NSS**2)
         CALL GETMEM('DYZR','FREE','REAL',LDYZR,NSS**2)
         CALL GETMEM('DYZI','FREE','REAL',LDYZI,NSS**2)

         CALL GETMEM('DZXR','FREE','REAL',LDZXR,NSS**2)
         CALL GETMEM('DZXI','FREE','REAL',LDZXI,NSS**2)
         CALL GETMEM('DZYR','FREE','REAL',LDZYR,NSS**2)
         CALL GETMEM('DZYI','FREE','REAL',LDZYI,NSS**2)
! Spin-Magnetic-Quadrupole
         CALL GETMEM('SXR','FREE','REAL',LSXR,NSS**2)
         CALL GETMEM('SXI','FREE','REAL',LSXI,NSS**2)
         CALL GETMEM('SYR','FREE','REAL',LSYR,NSS**2)
         CALL GETMEM('SYI','FREE','REAL',LSYI,NSS**2)
         CALL GETMEM('SZR','FREE','REAL',LSZR,NSS**2)
         CALL GETMEM('SZI','FREE','REAL',LSZI,NSS**2)
! Electric-Dipole
         CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
         CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
         CALL GETMEM('DYR','FREE','REAL',LDYR,NSS**2)
         CALL GETMEM('DYI','FREE','REAL',LDYI,NSS**2)
         CALL GETMEM('DZR','FREE','REAL',LDZR,NSS**2)
         CALL GETMEM('DZI','FREE','REAL',LDZI,NSS**2)

        IF(QIALL) THEN
         WRITE(6,*)'        ----------------------------'
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
     &   WRITE(6,*) 'Magnetic-dipole - Magnetic-dipole not included'
         IF(SECORD(2).EQ.0)
     &   WRITE(6,*) 'Electric-Quadrupole - Electric-Quadrupole not in'
         IF(SECORD(3).EQ.0)
     &   WRITE(6,*) 'Electric-Dipole - Electric-Octupole not included'
         IF(SECORD(4).EQ.0)
     &   WRITE(6,*) 'Electric-Dipole - Magnetic-Quadrupole not included'
         i_Print=0
         DO ISS=1,IEND
          DO JSS=JSTART,NSS
           EDIFF=ENSOR(JSS)-ENSOR(ISS)
           IF(EDIFF.GT.0.0D0) THEN
!
            IJSS=ISS+NSS*(JSS-1)
            F = WORK(LTOT2K-1+IJSS)
            IF(ABS(F).GE.OSTHR2) THEN
             IF(i_Print.eq.0) THEN
              i_Print=1
              Call CollapseOutput(1,
     &                  'Total transition strengths ' //
     &                  'for the second-order expansion of the wave '//
     &                  'vector (SO states):')
              WRITE(6,'(3X,A)')
     &                  '---------------------------' //
     &                  '-------------------------------------------'//
     &                  '-------------------'
!
              IF(OSTHR2.GT.0.0D0) THEN
               WRITE(6,*)'   for osc. strength at least ',OSTHR2
               WRITE(6,*)
              END IF
              WRITE(6,*)"        To  From     Osc. strength"
              WRITE(6,*)'        ----------------------------'
             END IF
             WRITE(6,'(5X,2I5,5X,ES16.8)') ISS,JSS,F
            END IF
           END IF
          END DO
         END DO
         If (i_Print.eq.1) THEN
           Call CollapseOutput(0,
     &                  'Total transition strengths ' //
     &                  'for the second-order expansion of the wave '//
     &                  'vector (SO states):')
           WRITE(6,*)
         END IF
       END IF
! release the memory again
         CALL GETMEM('TOT2K','FREE','REAL',LTOT2K,NSS**2)

      END IF
*
************************************************************************
*                                                                      *
*     Start of section for transition moments                          *
*                                                                      *
*     This section has two parts. (1) for matrix elements computed by  *
*     Seward, i.e. for a specific wavevector, (2) for the computation  *
*     of the isotropic oscillator strength.                            *
*                                                                      *
************************************************************************
*
*     Find the section of transition moments in the property list.
*
*     The operator is split in 4 different component, each with three
*     elements corresponding to differentiation in the x, y, and z
*     direction. The four parts are labels as:
*     EMFR  RS: The symmetric part of the real comp. of the op.
*     EMFR  RA: The asymmetric part of the real comp. of the op.
*     EMFR  IS: The symmetric part of the imaginary comp. of the op.
*     EMFR  IA: The asymmetric part of the imaginary comp. of the op.
*
************************************************************************
*                                                                      *
*     Section (2): Computation of the isotropic oscillator strength.
*                                                                      *
************************************************************************
*
      If (.Not.Do_TMOS) Go To 500
*
*     Here we will use a Lebedev grid to integrate over all possible
*     directions of the wave vector, k. The property integrals will be
*     computed on the fly and traced with the density to generate the
*     corresponding values in the PROP matrix.
*
*     Find the slot on the one-electron file where we will store the
*     on-the-fly generated property integrals.
*
      IPREMFR_RS=-1
      IPORIG=-1
      DO IPROP=1,NPROP
         IF (PNAME(IPROP).EQ.'TMOS  RS'.AND.IPREMFR_RS.EQ.-1) THEN
            IPREMFR_RS=IPROP
            IPORIG=IPROP
         END IF
      ENDDO
      IPREMFR_0R=IPREMFR_RS-6
      IPREMFR_0I=IPREMFR_RS-3
      IPREMFR_RA=IPREMFR_RS+3
      IPREMFR_IS=IPREMFR_RS+6
      IPREMFR_IA=IPREMFR_RS+9
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
         Call GetMem('SK','ALLO','REAL',ipR,4)
         Work(ipR  )=k_Vector(1)
         Work(ipR+1)=k_Vector(2)
         Work(ipR+2)=k_Vector(3)
         Work(ipR+3)=1.0D0   ! Dummy weight
      Else
         Call Setup_O()
         Call Do_Lebedev(L_Eff,nQuad,ipR)
      End If
*
      IF(DO_KVEC) THEN
*
*     Specific directions specified by user
*     The arrays from the Lebedev procedure is still used
*
        NQUAD = NKVEC
      END IF
*
*     Get table of content for density matrices.
*
      Call DaName(LuToM,FnToM)
      iDisk=0
      Call iDaFile(LuToM,2,iWork(liTocM),nState*(nState+1)/2,iDisk)
*
*     Allocate some temporary arrays for handling the
*     properties of the spin-orbit states.
*
      CALL GETMEM('DXR','ALLO','REAL',LDXR,NSS**2)
      CALL GETMEM('DXI','ALLO','REAL',LDXI,NSS**2)
      CALL GETMEM('TMP','ALLO','REAL',LTMP,NSS**2)
*
C     ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
      NIP=4+(NBST*(NBST+1))/2
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
      NSCR=(NBST*(NBST+1))/2
      CALL GETMEM('TDMSCR','Allo','Real',LSCR,4*NSCR)
*
*     Array for printing contributions from different directions
*
      CALL GETMEM('RAW   ','ALLO','REAL',LRAW,NQUAD*5)
      CALL GETMEM('WEIGHT','ALLO','REAL',LWEIGH,NQUAD*5)
*
*     Allocate vector to store all individual transition moments.
*     We do this for
*     all unique pairs ISO-JSO, iSO=/=JSO (NSS*(NSS-1)/2)
*         all k-vectors (3*nQuad)
*             all polarization directions (2*3)
*                 we store the transition moment (a complex number) (2*2)
*
      nIJ=nSS*(nSS-1)/2
      nData= 1 + 3 + 2*3 + 2*2
      nStorage = nIJ * nQuad * nData
      Call GetMem('STORAGE','Allo','Real',ipStorage,nStorage)
      ip_w      = 1
      ip_kvector= ip_w + 1
      ip_e1     = ip_kvector + 3
      ip_e2     = ip_e1 + 3
      ip_TM1R   = ip_e2 + 3
      ip_TM1I   = ip_TM1R + 1
      ip_TM2R   = ip_TM1I + 1
      ip_TM2I   = ip_TM2R + 1
*
      AFACTOR=32.1299D09
      HALF=0.5D0
      PI= CONST_PI_
      HBAR=1.0D0 ! in a.u.
      SPEED_OF_LIGHT=CONST_C_IN_AU_
      G_Elec=CONST_ELECTRON_G_FACTOR_
      iPrint=0
      IJSO=0
*
      Sparse_Limit=0.25D0
      ThrSparse=1.0D-6
      DO ISO=1, IEND
*
*        Check the sparseness of the coefficient array
*
         mSS=0
         Do i = 1, nSS
            temp = USOR(i,ISO)**2 + USOI(i,ISO)**2
            If (temp.gt.ThrSparse) mSS=mSS+1
         End Do
*
         Sparse_I = DBLE(MSS)/DBLE(NSS) .le. Sparse_Limit
*
         DO JSO=JSTART, NSS
*
*           Check the sparseness of the coefficient array
*
            mSS=0
            Do i = 1, nSS
               temp = USOR(i,JSO)**2 + USOI(i,JSO)**2
               If (temp.gt.ThrSparse) mSS=mSS+1
            End Do
*
            Sparse_J = DBLE(MSS)/DBLE(NSS) .le. Sparse_Limit
*
            EDIFF=ENSOR(JSO)-ENSOR(ISO)
            IF (ABS(EDIFF).LE.1.0D-8) CYCLE
            IF(EDIFF.LT.0.0D0) CYCLE
            IJSO=IJSO+1
            iOff_= (IJSO-1)*nQuad*nData
*
*           The energy difference is used to define the norm of the
*           wave vector.
*
            rkNorm=EDIFF/(HBAR*SPEED_OF_LIGHT)
*
*           Iterate over the quadrature points.
*
            FX=0.0D0
            FY=0.0D0
            FZ=0.0D0
            F =0.0D0
            Fm=0.0D0
*
*           Initialize output arrays
*
            CALL DCOPY_(NQUAD*5,0.0D0,0,WORK(LRAW),1)
            CALL DCOPY_(NQUAD*5,0.0D0,0,WORK(LWEIGH),1)

            Do iQuad = 1, nQuad
               iStorage = iOff_+ (iQuad-1)*nData + ipStorage - 1
*
*              Read or generate the wavevector
*
               IF(DO_KVEC) THEN
*
*              Get wavevector from input
*
               XCOOR = WORK(PKVEC+IQUAD-1)
               YCOOR = WORK(PKVEC+IQUAD-1+NQUAD)
               ZCOOR = WORK(PKVEC+IQUAD-1+2*NQUAD)
               ELSE
*
*              Generate the wavevector associated with this quadrature
*              point and pick up the associated quadrature weight.
*
               xcoor=Work((iQuad-1)*4  +ipR)
               ycoor=Work((iQuad-1)*4+1+ipR)
               zcoor=Work((iQuad-1)*4+2+ipR)
               END IF

               PORIG(1,IPREMFR_RS)=rkNorm*xcoor
               PORIG(2,IPREMFR_RS)=rkNorm*ycoor
               PORIG(3,IPREMFR_RS)=rkNorm*zcoor
               Call DCopy_(3,PORIG(1,IPREMFR_RS),1,
     &                       Work(iStorage+ip_kvector),1)
*
               Weight=Work((iQuad-1)*4+3+ipR)
               Work(iStorage+ip_w)=Weight
*
*              Generate the associated polarization vectors.
*
               IF (PORIG(1,IPREMFR_RS).EQ.0.0D0 .and.
     &             PORIG(2,IPREMFR_RS).EQ.0.0D0) Then
                  P1(1)=1.0D0
                  P1(2)=0.0D0
                  P1(3)=0.0D0
               ELSE
                  P1(1)= PORIG(2,IPREMFR_RS)
                  P1(2)=-PORIG(1,IPREMFR_RS)
                  P1(3)= 0.0D0
               END IF
               Tmp=1.0D0/SQRT(P1(1)**2+P1(2)**2+P1(3)**2)
               P1(1)=P1(1)*Tmp
               P1(2)=P1(2)*Tmp
               P1(3)=P1(3)*Tmp
               P2(1)=PORIG(2,IPREMFR_RS)*P1(3)-P1(2)*PORIG(3,IPREMFR_RS)
               P2(2)=PORIG(3,IPREMFR_RS)*P1(1)-P1(3)*PORIG(1,IPREMFR_RS)
               P2(3)=PORIG(1,IPREMFR_RS)*P1(2)-P1(1)*PORIG(2,IPREMFR_RS)
               Tmp=1.0D0/SQRT(P2(1)**2+P2(2)**2+P2(3)**2)
               P2(1)=P2(1)*Tmp
               P2(2)=P2(2)*Tmp
               P2(3)=P2(3)*Tmp
               Call DCopy_(3,P1,1,Work(iStorage+ip_e1),1)
               Call DCopy_(3,P2,1,Work(iStorage+ip_e2),1)
*
*              Compute the vectors (k x e1) and  (k x e2).
*
               kxe1(1)=PORIG(2,IPORIG)*P1(3)-PORIG(3,IPORIG)*P1(2)
               kxe1(2)=PORIG(3,IPORIG)*P1(1)-PORIG(1,IPORIG)*P1(3)
               kxe1(3)=PORIG(1,IPORIG)*P1(2)-PORIG(2,IPORIG)*P1(1)
               kxe2(1)=PORIG(2,IPORIG)*P2(3)-PORIG(3,IPORIG)*P2(2)
               kxe2(2)=PORIG(3,IPORIG)*P2(1)-PORIG(1,IPORIG)*P2(3)
               kxe2(3)=PORIG(1,IPORIG)*P2(2)-PORIG(2,IPORIG)*P2(1)
*
*              Generate the property integrals associated with this
*              direction of the wave vector k.
*
               Call TMOSInt(PORIG(1,IPREMFR_RS))
*
************************************************************************
*                                                                      *
*              Recompute the needed properties for all the spin-free   *
*              states.                                                 *
*                                                                      *
************************************************************************
*
               DO IPROP = IPREMFR_RS-6, IPREMFR_RS+11
                  Call FZero(PROP(1,1,IPROP),NSTATE**2)
               End Do
               ISS = 0
               DO I=1, NSTATE
*
*                 Does this spin-free state contribute to any of the
*                 two spin states? Check the corresponding coefficients.
*
                  JOB1=iWork(lJBNUM+I-1)
                  MPLET1 = MLTPLT(JOB1)
*
                  temp=0.0D0
                  Do MS1 = 1, MPLET1
                     ISS = ISS + 1
                     temp=Max(temp,USOR(ISS,ISO)**2 + USOI(ISS,ISO)**2)
                     temp=Max(temp,USOR(ISS,JSO)**2 + USOI(ISS,JSO)**2)
                  End Do
                  If (temp.le.ThrSparse)  Cycle
*
                  JSS=0
                  DO J=1, I
*
                     JOB2=iWork(lJBNUM+J-1)
                     MPLET2 = MLTPLT(JOB2)
*
                     temp=0.0D0
                     Do MS2 = 1, MPLET2
                        JSS = JSS + 1
                        temp=
     &                     Max(temp,USOR(JSS,ISO)**2 + USOI(JSS,ISO)**2)
                        temp=
     &                     Max(temp,USOR(JSS,JSO)**2 + USOI(JSS,JSO)**2)
                     End Do
                     If (temp.le.ThrSparse)  Cycle
*
*                    COMBINED SYMMETRY OF STATES:
                     JOB1=JBNUM(I)
                     JOB2=JBNUM(J)
                     LSYM1=IRREP(JOB1)
                     LSYM2=IRREP(JOB2)
                     ISY12=MUL(LSYM1,LSYM2)
*                    THE SYMMETRY CHECK MASK:
                     MASK=2**(ISY12-1)
*                    FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF
*                    TDMSCR
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
                     END DO
*
*                    Pick up the transition density between the two
*                    states from disc. Generated in PROPER.
*
                     ij=I*(I-1)/2+J
                     iDisk=iWork(liTocM+ij-1)
                     Call dDaFile(LuToM,2,Work(LSCR),4*NSCR,iDisk)
*
*                    Compute the transition property of the property
*                    integrals between the two states.
*
                     DO IPROP = IPREMFR_RS-6, IPREMFR_RS+11
                        ITYPE=0
                        IF (PTYPE(IPROP).EQ.'HERMSING') ITYPE=1
                        IF (PTYPE(IPROP).EQ.'ANTISING') ITYPE=2
                        IF (PTYPE(IPROP).EQ.'HERMTRIP') ITYPE=3
                        IF (PTYPE(IPROP).EQ.'ANTITRIP') ITYPE=4
                        LABEL=PNAME(IPROP)
                        Call MK_PROP(PROP,IPROP,I,J,LABEL,ITYPE,
     &                               WORK(LIP),NIP,WORK(LSCR),NSCR,
     &                               MASK,ISY12,IOFF)
                     END DO
*
                  END DO ! J
               END DO ! I
*
*              (1) the oam part
*
*              The contribution to the generalized momentum operator.
*
C              IJ=(JSO-1)*NSS + ISO - 1
               Do iCar = 1, 3
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LTMP),1)
*                 the real symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXR),NSS,'TMOS  RS',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*                 the real anti-symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LTMP),NSS,'TMOS  RA',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
                  CALL DAXPY_(NSS**2,1.0D0,WORK(LTMP),1,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LTMP),1)
*                 the imaginary symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXI),NSS,'TMOS  IS',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*                 the imaginary anti-symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LTMP),NSS,'TMOS  IA',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
                  CALL DAXPY_(NSS**2,1.0D0,WORK(LTMP),1,WORK(LDXI),1)
*                 Transform properties to the spin-orbit basis
*                 and pick up correct element
                  CALL ZTRNSF_IJ(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI),
     &                           WORK(LTMP),T0(iCar),ISO,JSO)
               End Do
*
               E1A = P1(1)*T0(1) + P1(2)*T0(2) + P1(3)*T0(3)
               E2A = P2(1)*T0(1) + P2(2)*T0(2) + P2(3)*T0(3)
*
*              (2) the magnetic-spin part
*
*              iCar=1: 1/2(S(+)+S(-))
*              iCar=2: 1/2i(S(+)-S(-))
*              iCar=3: Sz
*
               Do iCar = 1, 3
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
                  If (iCar.eq.1.or.iCar.eq.3) Then
*                    pick up the real component
                     CALL SMMAT_CHECK(PROP,WORK(LDXR),NSS,'TMOS0  R',
     &                                iCar,
     &                                USOR,USOI,ISO,JSO,ThrSparse)
*                 Else
*                    pick up the imaginary component
                     CALL SMMAT_CHECK(PROP,WORK(LDXI),NSS,'TMOS0  I',
     &                                iCar,
     &                                USOR,USOI,ISO,JSO,ThrSparse)
                  End If
                  CALL ZTRNSF_IJ(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI),
     &                           WORK(LTMP),T1(iCar),ISO,JSO)
               End Do
*
               E1B=kxe1(1)*T1(1)+kxe1(2)*T1(2)+kxe1(3)*T1(3)
               E2B=kxe2(1)*T1(1)+kxe2(2)*T1(2)+kxe2(3)*T1(3)
               TM1 = E1A + IMAGINARY*(g_Elec/2.0D0)*E1B
               TM2 = E2A + IMAGINARY*(g_Elec/2.0D0)*E2B
               Work(iStorage+ip_TM1R)=DBLE(TM1)
               Work(iStorage+ip_TM1I)=AIMAG(TM1)
               Work(iStorage+ip_TM2R)=DBLE(TM2)
               Work(iStorage+ip_TM2I)=AIMAG(TM2)
*
*              Integrate over all directions of the polarization
*              vector and divide with the "distance", 2*pi, to get
*              the average value.
*
               TM_2 = Half*DBLE(DCONJG(TM1)*TM1 + DCONJG(TM2)*TM2)
*
*              Compute the oscillator strength
*
               F_Temp = 2.0D0*TM_2/EDIFF
*
*              Do only the magnetic part.
*
               TM1 = IMAGINARY*(g_Elec/2.0D0)*E1B
               TM2 = IMAGINARY*(g_Elec/2.0D0)*E2B
               TM_2 = Half*DBLE(DCONJG(TM1)*TM1 + DCONJG(TM2)*TM2)
               F_Tempm = 2.0D0*TM_2/EDIFF
*
*              Save the raw oscillator strengths in a given direction
*
               WORK(LRAW+(IQUAD-1)+0*NQUAD) = F_TEMP
               WORK(LRAW+(IQUAD-1)+1*NQUAD) = F_TEMPM
               WORK(LRAW+(IQUAD-1)+2*NQUAD) = XCOOR
               WORK(LRAW+(IQUAD-1)+3*NQUAD) = YCOOR
               WORK(LRAW+(IQUAD-1)+4*NQUAD) = ZCOOR
*
*              Compute the oscillator strength
*
               F = F + Weight * F_Temp
*
*              Do only the magnetic part.
               Fm= Fm+ Weight * F_Tempm
*
*              Save the weighted oscillator strengths in a given direction
*
               WORK(LWEIGH+(IQUAD-1)+0*NQUAD) = F_TEMP*WEIGHT
               WORK(LWEIGH+(IQUAD-1)+1*NQUAD) = F_TEMPM*WEIGHT
               WORK(LWEIGH+(IQUAD-1)+2*NQUAD) = XCOOR
               WORK(LWEIGH+(IQUAD-1)+3*NQUAD) = YCOOR
               WORK(LWEIGH+(IQUAD-1)+4*NQUAD) = ZCOOR

            End Do ! iQuad
*
*           Note that the weights are normalized to integrate to
*           4*pi over the solid angles.
*
            F = F / (4.0D0*PI)
            Fm= Fm/ (4.0D0*PI)
            IF (ABS(F).LT.OSTHR) CYCLE
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
*
            IF(.NOT.DO_KVEC) THEN
              If (iPrint.eq.0) Then
                 WRITE(6,*)
                 If (Do_SK) Then
                    CALL CollapseOutput(1,
     &                  'Transition moment strengths:')
                 Else
                    CALL CollapseOutput(1,
     &                  'Isotropic transition moment strengths:')
                 End If
                 WRITE(6,'(3X,A)')
     &                  '--------------------------------------'
                 IF (OSTHR.GT.0.0D0) THEN
                    WRITE(6,'(1x,a,ES16.8)')
     &                    '   for osc. strength at least ',OSTHR
                 END IF
                 WRITE(6,*)
                 If (Do_SK) Then
                    WRITE(6,'(4x,a,3F8.4,a)')
     &                    'Direction of the k-vector: ',
     &                     (Work(ipR+k),k=0,2),' (au)'
                    WRITE(6,'(4x,a)')
     &                    'The light is assumed to be unpolarized.'
                 Else
                    WRITE(6,'(1x,a,I4,a)')
     &                '   Integrated over ',nQuad,' directions of the'//
     &                ' wave vector'
                    WRITE(6,'(1x,a)')
     &                '   Integrated over all directions of the polar'//
     &                'ization vector'
                 End If
                 WRITE(6,*)
                 WRITE(6,*)"        To  From     Osc. strength"//
     &             "   Einstein coefficients Ax, Ay, Az (sec-1) "//
     &             "      Total A (sec-1)  "
                 WRITE(6,*)
     &  '        -------------------------------------------'//
     &  '------------------------------------------------'
                iPrint=1
              END IF
*
              WRITE(6,'(5X,2I5,5X,5G16.8)') ISO,JSO,F,AX,AY,AZ,A
              WRITE(6,'(A,6X,G16.8)') ' Magnetic only',Fm
            END IF
*
*     Printing raw (unweighted) and direction for every transition
*
            IF(PRRAW) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)"        To  From     Raw Osc. str."//
     &          "   Mag. cont.       "//
     &          "   kx,            ky,            kz "
              WRITE(6,*)
     &  '        -------------------------------------------'//
     &  '------------------------------------------------'
              DO IQUAD = 1, NQUAD
                WRITE(6,'(5X,2I5,5X,5G16.8)') ISO,JSO,
     &          WORK(LRAW+(IQUAD-1)+0*NQUAD),
     &          WORK(LRAW+(IQUAD-1)+1*NQUAD),
     &          WORK(LRAW+(IQUAD-1)+2*NQUAD),
     &          WORK(LRAW+(IQUAD-1)+3*NQUAD),
     &          WORK(LRAW+(IQUAD-1)+4*NQUAD)
              END DO
              WRITE(6,*)
              WRITE(6,*)
            END IF
*
*     Printing weighted and direction for every transition
*
            IF(PRWEIGHT) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)"        To  From     Wei Osc. str."//
     &          "   Mag. cont.       "//
     &          "   kx,            ky,            kz "
              WRITE(6,*)
     &  '        -------------------------------------------'//
     &  '------------------------------------------------'
              DO IQUAD = 1, NQUAD
                WRITE(6,'(5X,2I5,5X,5G16.8)') ISO,JSO,
     &          WORK(LWEIGH+(IQUAD-1)+0*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH+(IQUAD-1)+1*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH+(IQUAD-1)+2*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH+(IQUAD-1)+3*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH+(IQUAD-1)+4*NQUAD)/ (4.0D0*PI)
              END DO
              WRITE(6,*)
              WRITE(6,*)
            END IF
*
            Call Add_Info('ITMS(SO)',F,1,6)
*
         END DO
      END DO
      If (IPGLOB.EQ.1) THEN
         If (Do_SK) Then
            CALL CollapseOutput(0,
     &                'Transition moment strengths:')
         Else
            CALL CollapseOutput(0,
     &                'Isotropic transition moment strengths:')
         End If
      END IF
*
#ifdef _HDF5_
      Call mh5_put_dset(wfn_sos_tm,Work(ipStorage))
#endif
*
*     Do some cleanup
*
      Call GetMem('STORAGE','FREE','Real',ipStorage,nStorage)
      CALL GETMEM('RAW   ','FREE','REAL',LRAW,NQUAD*5)
      CALL GETMEM('WEIGHT','FREE','REAL',LWEIGH,NQUAD*5)
      CALL GETMEM('TDMSCR','FREE','Real',LSCR,4*NSCR)
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
      CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
      CALL GETMEM('TMP','FREE','REAL',LTMP,NSS**2)
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
     & 'D_',IXYZ,':',EVR(IXYZ)*AU2CM,'cm^-1'
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

C print seperate contributions if verbose
      IF (IPGLOB.GE.VERBOSE) THEN
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

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZI),1)

      IF(IAMX.GT.0)
     &     CALL SMMAT(PROP,WORK(LLXI),NSS,PNAME(IAMX),ICOMP(IAMX))
      IF(IAMY.GT.0)
     &    CALL SMMAT(PROP,WORK(LLYI),NSS,PNAME(IAMY),ICOMP(IAMY))
      IF(IAMZ.GT.0)
     &    CALL SMMAT(PROP,WORK(LLZI),NSS,PNAME(IAMZ),ICOMP(IAMZ))

* PAM09 -- This code appears to be unused:
*      CALL GETMEM('LXR','ALLO','REAL',LLXR,NSS**2)
*      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXR),1)
*      CALL GETMEM('LYR','ALLO','REAL',LLYR,NSS**2)
*      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYR),1)
*      CALL GETMEM('LZR','ALLO','REAL',LLZR,NSS**2)
*      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZR),1)
*------------------------

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)
      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)
      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)


      CALL SMMAT(PROP,WORK(LZXR),NSS,'SPIN    ',1)
      CALL SMMAT(PROP,WORK(LZYI),NSS,'SPIN    ',2)
      CALL SMMAT(PROP,WORK(LZZR),NSS,'SPIN    ',3)

      CALL DSCAL_(NSS**2,FEGVAL,WORK(LZXR),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LZYI),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LZZR),1)

      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LZXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LZYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LZZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

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
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOm(L,I,J)=(0.0d0,0.0d0)
       DIPSOn(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo


*     Continue original calculation of G tensor (=gg^*)
      CALL get_dArray( 'ESO_SINGLE',ESO,NSS)
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOm(1,ISS,JSS)=0.5d0*Z(ISS,JSS)
      DIPSOn(1,ISS,JSS)=-Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOm(2,ISS,JSS)=0.5d0*Z(ISS,JSS)
      DIPSOn(2,ISS,JSS)=-Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
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
      c_1(ic,jc)=DBLE(DIPSOn(ic,Iss,Jss)*DCONJG(DIPSOn(jc,Iss,Jss)))
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

      !IF (KDGN.NE.2) THEN
      !    WRITE(6,*) 'no twofold degeneracy'
      !    GOTO 780
      !ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       IF(.NOT.IFGTCALSA) GOTO 450
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOn,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  450 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      IF(IPGLOB.GT.VERBOSE) THEN
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

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

 800  CONTINUE



* Skip if not a hyperfine calculation
      IF(.NOT.IFACAL) GOTO 1801


*****************************************************
*****************************************************
*****************************************************
* Experimental hyperfine tensor stuff starts here
*****************************************************
*****************************************************
*****************************************************
      IF(IFSONCINI) THEN
      WRITE(6,*)
      WRITE(6,*) '  Soncini pNMR Tensor and A-Matrix Approach II '
      WRITE(6,*) '  ============================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ELSE
      WRITE(6,*)
      WRITE(6,*) '  A-Matrix Approach II                     '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  1st order degenerate perturbation theory '
      WRITE(6,*) '  within isolated kramers doublets.        '
      WRITE(6,*) '  > spatial degeneracy'
      WRITE(6,*) '  > strong spin-orbit coupling'
      WRITE(6,*)
      ENDIF

C Mapping from spin states to spin-free state and to spin:
      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)

      ISS=0
      DO ISTATE=1,NSTATE
       JOB=JBNUM(ISTATE)
       MPLET=MLTPLT(JOB)
       DO MSPROJ=-MPLET+1,MPLET-1,2
        ISS=ISS+1
        IWORK(LMAPST-1+ISS)=ISTATE
        IWORK(LMAPSP-1+ISS)=MPLET
        IWORK(LMAPMS-1+ISS)=MSPROJ
       END DO
      END DO


      IAMX=0

      DO IPROP=1,NPROP

        IF(PNAME(IPROP)(1:3).EQ.'ASD'.AND.ICOMP(IPROP).EQ.1) THEN

* Get the center number
         Read (PNAME(IPROP),'(a4,i4)') ASDLAB,ICEN

      WRITE(6,*) '  ========================================='
      write(6,*) '  A(Total)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='



C Identify which properties are ASD matrix elements:
c Labeled AMFI for now
c 1,2,3,4,5,6 -> xx,xy,xz,yy,yz,zz
      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
C Identify which properties are Orbital Paramagnetic (PSOP) matrix elements:
      IAMX=0
      IAMY=0
      IAMZ=0

      WRITE(EFPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",EFPROP
      WRITE(PSOPROP,'(a4,i4)') 'PSOP',ICEN
      WRITE(6,*) "Looking for ",PSOPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP).EQ.EFPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       ELSE IF(PNAME(KPROP).EQ.PSOPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMX=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMY=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMZ=KPROP
       END IF
      END DO

cccccccccccccccccccccccccccccccccccccccc
c Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc
c      EFPROP = 'MLTPL  0'
c      WRITE(6,*) "Looking for overlap matrix ",EFPROP
c      DO KPROP=1,NPROP
c       IF(PNAME(KPROP).EQ.EFPROP) THEN
c         IAMFI1=KPROP
c         IAMFI2=0
c         IAMFI3=0
c         IAMFI4=0
c         IAMFI5=0
c         IAMFI6=0
c       END IF
c      END DO
cccccccccccccccccccccccccccccccccccccccc
c End testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc


* These will hold the entire expression
* For g: L+2S
* For A: ?

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZI),1)

      IF(IAMX.GT.0)
     &     CALL SMMAT(PROP,WORK(LLXI),NSS,PNAME(IAMX),ICOMP(IAMX))
      IF(IAMY.GT.0)
     &    CALL SMMAT(PROP,WORK(LLYI),NSS,PNAME(IAMY),ICOMP(IAMY))
      IF(IAMZ.GT.0)
     &    CALL SMMAT(PROP,WORK(LLZI),NSS,PNAME(IAMZ),ICOMP(IAMZ))


      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)


      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)


* Note that the 6th element contains r*r
* This is 4*pi * the contact term we want
* We need 8*pi/3 * < delta > so we just want
* a factor of 2/3
          ACNT = 0d0

          IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
           !write(6,*)"ACNT",ACNT
          ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
             WRITE(6,*) "********************"
             WRITE(6,*) "* Skipping FC Part *"
             WRITE(6,*) "********************"
          END IF

* Since this is a traceless tensor, generate the 6th
* element (3z^2-r^2)
          AMFI6=-AMFI1-AMFI4

          IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          ELSE
             IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
                WRITE(6,*) "********************"
                WRITE(6,*) "* Skipping SD Part *"
                WRITE(6,*) "********************"
             END IF
             AMFI1=0.0d0
             AMFI2=0.0d0
             AMFI3=0.0d0
             AMFI4=0.0d0
             AMFI5=0.0d0
             AMFI6=0.0d0
          END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          CALL ADD_INFO("ASDFC1",AMFI1,1,5)
          CALL ADD_INFO("ASDFC2",AMFI2,1,5)
          CALL ADD_INFO("ASDFC3",AMFI3,1,5)
          CALL ADD_INFO("ASDFC4",AMFI4,1,5)
          CALL ADD_INFO("ASDFC5",AMFI5,1,5)
          CALL ADD_INFO("ASDFC6",AMFI6,1,5)

cccccccccccccccccccccccccccccccccccccccc
c Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc
c          AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
c          AMFI2=0.0d0
c          AMFI3=0.0d0
c          AMFI4=AMFI1
c          AMFI5=0.0d0
c          AMFI6=AMFI1
cccccccccccccccccccccccccccccccccccccccc
c END Testing - use overlap matrix
cccccccccccccccccccccccccccccccccccccccc

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO


      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LZXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LZYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LZZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)


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
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOf(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo

c
c*     Continue original calculation of G tensor (=gg^*)
c

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL PRCMAT(NSS,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOf(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      IF(IFSONCINI) THEN

      do ic=1,3
      do jc=1,3
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIMSO(ic,jc,ISS,JSS)=0.0d0
      enddo
      enddo
      enddo
      enddo
      WRITE(DMPPROP,'(a4,i4)') 'DMP ',ICEN
      WRITE(6,*) "Looking for ",DMPPROP
      DO KPROP=1,NPROP
      IF(PNAME(KPROP).EQ.DMPPROP) THEN
      CALL GETMEM('SOPROPR','ALLO','REAL',LSOPRR,NSS**2)
      CALL GETMEM('SOPROPI','ALLO','REAL',LSOPRI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSOPRR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LSOPRI),1)

      CALL SMMAT(PROP,WORK(LSOPRR),NSS,PNAME(KPROP),ICOMP(KPROP))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LSOPRR),WORK(LSOPRI))
      !CALL PRCMAT(NSS,WORK(LSOPRR),WORK(LSOPRI))
      CALL MULMAT(NSS,WORK(LSOPRR),WORK(LSOPRI),eex,Z)
      ic=0
      jc=0
      !write(6,*)"ICOMP(KPROP)",ICOMP(KPROP)
      if (ICOMP(KPROP).EQ.1) then
      ic=1
      jc=1
      elseif (ICOMP(KPROP).EQ.2) then
      ic=1
      jc=2
      elseif (ICOMP(KPROP).EQ.3) then
      ic=1
      jc=3
      elseif (ICOMP(KPROP).EQ.4) then
      ic=2
      jc=1
      elseif (ICOMP(KPROP).EQ.5) then
      ic=2
      jc=2
      elseif (ICOMP(KPROP).EQ.6) then
      ic=2
      jc=3
      elseif (ICOMP(KPROP).EQ.7) then
      ic=3
      jc=1
      elseif (ICOMP(KPROP).EQ.8) then
      ic=3
      jc=2
      else
      ic=3
      jc=3
      endif
      DO ISS=1,NSS
      DO JSS=1,NSS

      if (ISS.eq.JSS) then
      !DIMSO(ICOMP(KPROP),ISS,JSS)=Z(ISS,JSS)
      DIMSO(ic,jc,ISS,JSS)=Z(ISS,JSS)
      ELSE
      DIMSO(ic,jc,ISS,JSS)=0
      ENDIF
      enddo
      enddo
      CALL GETMEM('SOPROPR','FREE','REAL',LSOPRR,NSS**2)
      CALL GETMEM('SOPROPI','FREE','REAL',LSOPRI,NSS**2)
      ENDIF

      END DO
      iT=0
      do iT=1,NTP
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=0.d0
      PNMR(iT,ic,jc)=0.d0
      PNMRC(iT,ic,jc)=0.d0
      PNMRD(iT,ic,jc)=0.d0
      enddo
      enddo
      enddo

      do iT=1,NTP
      do Iss=1,Nss
      do ic=1,3
      do jc=1,3
      PNMRCPS(iT,Iss,ic,jc)=0.d0
      enddo
      enddo
      enddo
      enddo

      iT=0
      do iT=1,NTP
      if(iT.eq.1) then
      TMPf(iT)=TMINP+0.0001d0
      ELSE
      DLTTA=(TMAXP-TMINP)/(dble(NTP-1))
      TMPf(iT)=TMINP+DLTTA*dble(iT-1)
      ENDIF
      Zstat=0.d0
      do Iss=1,Nss
      p_Boltz=EXP(-ESO(Iss)/Boltz_k/TMPf(iT))
      !write(6,*)'p_Boltz',p_Boltz
      Zstat=Zstat+p_Boltz
      do IC=1,3
      do JC=1,3
      DiamT(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      CurieT(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      HFC_2(IC,JC)=0.d0
      enddo
      enddo
      do IC=1,3
      do JC=1,3
      HFC_3(IC,JC)=0.d0
      enddo
      enddo
      do Jss=1,Nss
      dlt_E=Eso(Iss)-Eso(Jss)
      do IC=1,3
      do JC=1,3
      HFC_1(IC,JC)=0.d0
      enddo
      enddo
      do ic=1,3
      do jc=1,3

      HFC_1(ic,jc)=DBLE(DIPSOm(ic,Iss,Jss)*
     & DCONJG(DIPSOf(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))

!!      HFC_1(ic,jc)=DBLE(DIPSOf(ic,Iss,Jss)*
!!     & DCONJG(DIPSOm(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))


!      HFC_3(ic,jc)= 1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))/(ALPHA2*AU2CM)

      !if(ABS(dlt_E).LT.1.D-3) then
      if(ABS(dlt_E).LT.10.97D0) then
      CurieT(ic,jc)=   CurieT(ic,jc)+HFC_1(ic,jc)
      HFC_2(ic,jc)=    HFC_2(ic,jc)+ HFC_1(ic,jc)
      HFC_3(ic,jc)=    HFC_3(ic,jc)+ 0.d0*HFC_1(ic,jc)
      DiamT(ic,jc)=    DiamT(ic,jc)+ 0.d0*HFC_1(ic,jc)
      else

      HFC_2(ic,jc)= HFC_2(ic,jc)-(DBLE(DCONJG(DIPSOm(ic,Iss,Jss))*
     & DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & DCONJG(DIPSOf(jc,Iss,Jss))))/dlt_E
!      HFC_2(ic,jc)= HFC_2(ic,jc)-2.d0*(DBLE(DCONJG(DIPSOm(ic,Iss,Jss))*
!     & DIPSOf(jc,Iss,Jss)))/dlt_E

      HFC_3(ic,jc)= HFC_3(ic,jc)-(DBLE(DCONJG(DIPSOm(ic,Iss,Jss))*
     & DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & DCONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      CurieT(ic,jc)= CurieT(ic,jc)-0.d0*
     &(DBLE(DCONJG(DIPSOm(ic,Iss,Jss))*
     &DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & DCONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      DiamT(ic,jc)=DiamT(ic,jc)-0.d0*
     &(DBLE(DCONJG(DIPSOm(ic,Iss,Jss))*
     &DIPSOf(jc,Iss,Jss)) + DBLE(DIPSOm(ic,Iss,Jss)*
     & DCONJG(DIPSOf(jc,Iss,Jss))))/dlt_E

      endif

      HFC_2(ic,jc)= HFC_2(ic,jc) + (1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))
     & /(ALPHA2*AU2CM))
      DiamT(ic,jc)=DiamT(ic,jc)+ (1.D-6*DBLE(DIMSO(ic,jc,Iss,Jss))
     & /(ALPHA2*AU2CM))
      enddo
      enddo
      enddo !Jss
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=    PNMRT(iT,ic,jc)+ p_Boltz*
     & HFC_2(ic,jc)

      PNMR(iT,ic,jc)=    PNMR(iT,ic,jc)+ p_Boltz*
     & HFC_3(ic,jc)

      PNMRC(iT,ic,jc)=    PNMRC(iT,ic,jc)+ p_Boltz*
     & CurieT(ic,jc)

      PNMRCPS(iT,Iss,ic,jc)= PNMRCPS(iT,Iss,ic,jc)+CurieT(ic,jc)
      PNMRCPS(iT,Iss,ic,jc)=1.D6*AU2CM*ALPHA2*PNMRCPS(iT,Iss,ic,jc)

      PNMRD(iT,ic,jc)=    PNMRD(iT,ic,jc)+ p_Boltz*
     & DiamT(ic,jc)

      enddo
      enddo
      enddo !Iss
      do ic=1,3
      do jc=1,3
      PNMRT(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRT(iT,ic,jc)/Zstat)
      PNMR(iT,ic,jc) =1.D6*AU2CM*ALPHA2*(PNMR(iT,ic,jc)/Zstat)
      PNMRC(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRC(iT,ic,jc)/Zstat)
      PNMRD(iT,ic,jc)=1.D6*AU2CM*ALPHA2*(PNMRD(iT,ic,jc)/Zstat)
      !write(6,*) PNMRT(iT,ic,jc)
      !AU2CM*ALPHA2*
      enddo
      enddo

      enddo ! iT

      !IT=0
      !IC=0
      !JC=0
      !Do IT=1,NTP
      !Do IC=1,3
      !Do JC=1,3
      !PNMRT(IT,IC,JC)=1.d6*PNMRT(IT,IC,JC)
      !enddo
      !enddo
      !enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') ' TOTAL SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),
     & ((PNMRT(iT,ic,jc),jc=1,3),ic=1,3)
      enddo
      !write(6,*)'*******'
      !do iT=1,NTP
      !do ic=1,3
      !do jc=1,3
      !write(6,*)PNMRT_tens(iT,ic,jc)
      !enddo
      !enddo
      !enddo
      !write(6,*)'*******'

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'paramagnetic SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),
     & ((PNMR(iT,ic,jc),jc=1,3),ic=1,3)
      enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Diamagnetic  SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),
     & ((PNMRD(iT,ic,jc),jc=1,3),ic=1,3)
      enddo

      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Curie term contrib. to SONCINI PNMR TENSOR'//
     & ' in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,9X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      do iT=1,NTP
      write(6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),
     & ((PNMRC(iT,ic,jc),jc=1,3),ic=1,3)
!       IF (.EQV..TRUE.) THEN
      write(6,'(/)')
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,'(30X,A)') 'Curie per state in (ppm)'
      write(6,'(10A)') (('------------'), K=1,10)
      write(6,*)
C      write(6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)',
C     & '(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(6,'(6X,A,8X,9(A,9X))') 'NSS','xx','xy','xz','yx','yy',
     & 'yz','zx','zy','zz'
      write(6,*)
      !write(6,*)""
      !write(6,*) PNMRC(iT,3,3)
      !write(6,*)""
      do Iss=1,NSS
      write(6,'(4X,I5,3X,11(F10.4,2X),F8.4)') Iss,
     & ((PNMRCPS(iT,Iss,ic,jc),jc=1,3),ic=1,3)
       enddo

!       ENDIF
       write(6,'(/)')
       write(6,'(10A)') (('============='), K=1,10)
       write(6,'(/)')
      enddo
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '  A-Matrix'
      WRITE(6,*) '  =========='
      ENDIF !IFSONCINI

! 1111 Continue

      ISS=1
      DO WHILE ( (ISS.LE.NSS) .AND.
     &          ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRTHR) )

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
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
       IF(.NOT.IFATCALSA) GOTO 451
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOf,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  451 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 1780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.VERBOSE) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

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
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

      CALL ADD_INFO("ATENS",GTENS,9,5)
      CALL ADD_INFO("ATENS2",EVR,3,5)

 1780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

      IF(.NOT.IFACALFCON) GOTO 1901

      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (FC)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(EFPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",EFPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP).EQ.EFPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       END IF
      END DO

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZI),1)

      IF(IAMX.GT.0)
     &     CALL SMMAT(PROP,WORK(LLXI),NSS,PNAME(IAMX),ICOMP(IAMX))
      IF(IAMY.GT.0)
     &    CALL SMMAT(PROP,WORK(LLYI),NSS,PNAME(IAMY),ICOMP(IAMY))
      IF(IAMZ.GT.0)
     &    CALL SMMAT(PROP,WORK(LLZI),NSS,PNAME(IAMZ),ICOMP(IAMZ))

      CALL GETMEM('MXR','ALLO','REAL',LMXR,NSS**2)
      CALL GETMEM('MXI','ALLO','REAL',LMXI,NSS**2)
      CALL GETMEM('MYR','ALLO','REAL',LMYR,NSS**2)
      CALL GETMEM('MYI','ALLO','REAL',LMYI,NSS**2)
      CALL GETMEM('MZR','ALLO','REAL',LMZR,NSS**2)
      CALL GETMEM('MZI','ALLO','REAL',LMZI,NSS**2)

      IMR(1)=LMXR
      IMI(1)=LMXI
      IMR(2)=LMYR
      IMI(2)=LMYI
      IMR(3)=LMZR
      IMI(3)=LMZI

      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMXI),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMYI),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMZI),1)

      CALL SMMAT(PROP,WORK(LMXR),NSS,'SPIN    ',1)
      CALL SMMAT(PROP,WORK(LMYI),NSS,'SPIN    ',2)
      CALL SMMAT(PROP,WORK(LMZR),NSS,'SPIN    ',3)

      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMXR),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMYI),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMZR),1)

      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LMXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LMYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LMZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMXR),WORK(LMXI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMYR),WORK(LMYI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMZR),WORK(LMZI))

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0
* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
         ! END IF
          !AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
           !   AMFI1=-1.0*AMFI1
           !   AMFI2=-1.0*AMFI2
           !   AMFI3=-1.0*AMFI3
           !   AMFI4=-1.0*AMFI4
           !   AMFI5=-1.0*AMFI5
           !   AMFI6=-1.0*AMFI6
          !ELSE
          !   IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !      WRITE(6,*) "********************"
          !      WRITE(6,*) "* Skipping SD Part *"
          !      WRITE(6,*) "********************"
          !   END IF
             AMFI1=0.0d0
             AMFI2=0.0d0
             AMFI3=0.0d0
             AMFI4=0.0d0
             AMFI5=0.0d0
             AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO

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
      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfc(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfc(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRTHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
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
       IF(.NOT.IFATCALSA) GOTO 452
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfc,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  452 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 2780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.VERBOSE) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

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

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
       GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &  TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO


 2780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1901  CONTINUE

      IF(.NOT.IFACALSDON) GOTO 1902

      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (SD)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(EFPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",EFPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP).EQ.EFPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
       END IF
      END DO

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
          !   ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
          !END IF
          AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          !ELSE
           !  IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
            !    WRITE(6,*) "********************"
             !   WRITE(6,*) "* Skipping SD Part *"
              !  WRITE(6,*) "********************"
             !END IF
             !AMFI1=0.0d0
             !AMFI2=0.0d0
             !AMFI3=0.0d0
             !AMFI4=0.0d0
             !AMFI5=0.0d0
             !AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)

          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5
        END DO
      END DO

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

      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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

      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfsd(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfsd(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo

      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRTHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
      KDGN=MULTIP
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
       IF(.NOT.IFATCALSA) GOTO 453
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfsd,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  453 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 3780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.VERBOSE) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

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

      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
       GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO

      !DO IXYZ=1,3
      !   write(6,*) 'IXYZ',IXYZ
      !   write(6,*) ''
      ! DO JXYZ=1,3
      !   write(6,*) 'JXYZ',JXYZ
      !   write(6,*) ''
      ! GTENS(IXYZ,JXYZ)=0.0D0
      !  DO KXYZ=1,3
      !   write(6,*) 'KXYZ',KXYZ
      !   write(6,*) ''
      !   write(6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
      !   write(6,*) ''
      !   write(6,*) 'TMPVEC(IXYZ,KXYZ)',TMPVEC(IXYZ,KXYZ),IXYZ,KXYZ
      !   write(6,*) ''
      !   write(6,*) 'SQRT(EVR(KXYZ))',SQRT(EVR(KXYZ)),KXYZ
      !   write(6,*) ''
      !   write(6,*) 'TMPVEC(JXYZ,KXYZ)',TMPVEC(JXYZ,KXYZ),JXYZ,KXYZ
      !   write(6,*) '*************************'
      !  GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
      !&   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
      !  write(6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
      !   write(6,*) '*************************'
      !  END DO
      ! END DO
      !ENDDO

      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

 3780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1902  CONTINUE

* Skip if not a hyperfine calculation
      IF(.NOT.IFACALFCSDON) GOTO 1903


      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (FC+SD)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

      IAMFI1=0
      IAMFI2=0
      IAMFI3=0
      IAMFI4=0
      IAMFI5=0
      IAMFI6=0
      WRITE(EFPROP,'(a4,i4)') 'ASD ',ICEN
      WRITE(6,*) "Looking for ",EFPROP
      DO KPROP=1,NPROP
       IF(PNAME(KPROP).EQ.EFPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMFI1=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMFI2=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMFI3=KPROP
         IF(ICOMP(KPROP).EQ.4) IAMFI4=KPROP
         IF(ICOMP(KPROP).EQ.5) IAMFI5=KPROP
         IF(ICOMP(KPROP).EQ.6) IAMFI6=KPROP
      END IF
      END DO

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)

      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFI1=0.0D0
          AMFI2=0.0D0
          AMFI3=0.0D0
          AMFI4=0.0D0
          AMFI5=0.0D0
          AMFI6=0.0D0

* SD
          IF(IAMFI1.NE.0) AMFI1=PROP(ISTATE,JSTATE,IAMFI1)
          IF(IAMFI2.NE.0) AMFI2=PROP(ISTATE,JSTATE,IAMFI2)
          IF(IAMFI3.NE.0) AMFI3=PROP(ISTATE,JSTATE,IAMFI3)
          IF(IAMFI4.NE.0) AMFI4=PROP(ISTATE,JSTATE,IAMFI4)
          IF(IAMFI5.NE.0) AMFI5=PROP(ISTATE,JSTATE,IAMFI5)
          IF(IAMFI6.NE.0) AMFI6=PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = 0d0

          !IF(IFACALFC) THEN
             ACNT=AMFI6*2.0d0/3.0d0
          !ELSE IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
          !   WRITE(6,*) "********************"
          !   WRITE(6,*) "* Skipping FC Part *"
          !   WRITE(6,*) "********************"
          !END IF

          AMFI6=-AMFI1-AMFI4

          !IF(IFACALSD) THEN
              AMFI1=-1.0*AMFI1
              AMFI2=-1.0*AMFI2
              AMFI3=-1.0*AMFI3
              AMFI4=-1.0*AMFI4
              AMFI5=-1.0*AMFI5
              AMFI6=-1.0*AMFI6
          !ELSE
           !  IF(ISS.EQ.1.AND.JSS.EQ.1) THEN
           !     WRITE(6,*) "********************"
           !     WRITE(6,*) "* Skipping SD Part *"
           !     WRITE(6,*) "********************"
           !  END IF
           !  AMFI1=0.0d0
           !  AMFI2=0.0d0
           !  AMFI3=0.0d0
           !  AMFI4=0.0d0
           !  AMFI5=0.0d0
           !  AMFI6=0.0d0
          !END IF

          AMFI1=AMFI1+ACNT
          AMFI4=AMFI4+ACNT
          AMFI6=AMFI6+ACNT

          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)


          WORK(LZXR-1+IJSS)=CGX*AMFI1+CG0*AMFI3
          WORK(LZXI-1+IJSS)=CGY*AMFI2
          WORK(LZYR-1+IJSS)=CGX*AMFI2+CG0*AMFI5
          WORK(LZYI-1+IJSS)=CGY*AMFI4
          WORK(LZZR-1+IJSS)=CGX*AMFI3+CG0*AMFI6
          WORK(LZZI-1+IJSS)=CGY*AMFI5

        END DO
      END DO

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

      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,JSS,ISS)
           ENDIF
C     WRITE(6,FMT=710) 'ZEKL', ISTATE, IXYZ, ISS, JSS,
C     &                 ZEKL(:,:,IXYZ,ISTATE)
          ENDDO
         ENDDO
        ENDDO

       ENDIF

       IMLTPL=IMLTPL+MLTPLT(JBNUM(ISTATE))
      ENDDO
      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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


        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfcsd(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfcsd(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo


      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRTHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
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
       IF(.NOT.IFATCALSA) GOTO 454
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfcsd,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  454 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 4780
      ENDIF

      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)
        GTIJ=GTIJ+CONTRIB
        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      !IF(IPGLOB.GT.VERBOSE) THEN
       WRITE(6,*) 'G tensor = gg+'
       WRITE(6,*)
       WRITE(6,'(6x,3(6x,a2,4x))')
     &  (xyzchr(IXYZ),IXYZ=1,3)
       DO IXYZ=1,3
       WRITE(6,'(2x,a2,2x,3(1x,f18.8,1x))')
     &  xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3)
       ENDDO
      !END IF

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
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO

 4780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

 1903  CONTINUE

* Skip if not a hyperfine calculation
      IF(.NOT.IFACALPSO) GOTO 1904


      WRITE(6,*) ' '
      WRITE(6,*) '  ========================================='
      WRITE(6,*) '  A (PSO)-Matrix for center:',ICEN
      WRITE(6,*) '  ========================================='

       IAMX=0
       IAMY=0
       IAMZ=0

      WRITE(PSOPROP,'(a4,i4)') 'PSOP',ICEN
      WRITE(6,*) "Looking for ",PSOPROP
      DO KPROP=1,NPROP
         IF(PNAME(KPROP).EQ.PSOPROP) THEN
         IF(ICOMP(KPROP).EQ.1) IAMX=KPROP
         IF(ICOMP(KPROP).EQ.2) IAMY=KPROP
         IF(ICOMP(KPROP).EQ.3) IAMZ=KPROP
       END IF
      END DO

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZI),1)

      IF(IAMX.GT.0)
     &     CALL SMMAT(PROP,WORK(LLXI),NSS,PNAME(IAMX),ICOMP(IAMX))
      IF(IAMY.GT.0)
     &    CALL SMMAT(PROP,WORK(LLYI),NSS,PNAME(IAMY),ICOMP(IAMY))
      IF(IAMZ.GT.0)
     &    CALL SMMAT(PROP,WORK(LLZI),NSS,PNAME(IAMZ),ICOMP(IAMZ))


      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      IZMR(1)=LZXR
      IZMI(1)=LZXI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZXI),1)

      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZYI),1)

      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)
      IZMR(3)=LZZR
      IZMI(3)=LZZI
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZZI),1)


      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LZXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LZYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LZZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

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
******************************************************
      IMLTPL=1
      DO ISTATE=1,NSTATE
       DO IXYZ=1,3
        DO J=1,2
         DO I=1,2
          ZEKL(I,J,IXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

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
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,ISS,JSS)
           CALL ZECON(NSTATE,NSS,USOR,USOI,
     &          WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &          ZEKL,IXYZ,ISTATE,JSS,ISS)
           IF (ISGS(JSS)) THEN
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
     &           ZEKL,IXYZ,ISTATE,ISS,JSS)
            CALL ZECON(NSTATE,NSS,USOR,USOI,
     &           WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)),
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


      DO ISTATE=1,NSTATE
       DO IJXYZ=1,9
        GCONT(IJXYZ,ISTATE)=DCMPLX(0.0d0,0.0d0)
       ENDDO
      ENDDO
      DO IJXYZ=1,9
       GTOTAL(IJXYZ)=0.0d0
      ENDDO

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

        do I=1,NSS
         do J=1,NSS
          do L=1,3
       DIPSOfpso(L,I,J)=(0.0d0,0.0d0)
          enddo
         enddo
       enddo
c
c*     Continue original calculation of G tensor (=gg^*)
c
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZXR),WORK(LZXI))
      CALL MULMAT(NSS,WORK(LZXR),WORK(LZXI),eex,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(1,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZYR),WORK(LZYI))
      CALL MULMAT(NSS,WORK(LZYR),WORK(LZYI),eey,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(2,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LZZR),WORK(LZZI))
      CALL MULMAT(NSS,WORK(LZZR),WORK(LZZI),eez,Z)
      DO ISS=1,NSS
      DO JSS=1,NSS
      DIPSOfpso(3,ISS,JSS)=Z(ISS,JSS)
      enddo
      enddo


      ISS=1
      DO WHILE ((ENSOR(MIN(ISS,NSS))-ENSOR(1)).LE.EPRTHR
     &  .AND.(ISS.LE.NSS))

      DO IXYZ=1,3
       DO JXYZ=1,3
        GTENS(IXYZ,JXYZ)=0.0D0
       END DO
      END DO

      KDGN=1
      DO JSS=ISS+1,NSS
      EDIFF=ENSOR(JSS)-ENSOR(ISS)
      IF (IFATCALSA.AND.IFGTSHSA) THEN
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
       IF(.NOT.IFATCALSA) GOTO 455
       if(ISS.EQ.1) IFUNCT=0
      call SINANI(KDGN,IFUNCT,NSS,DIPSOfpso,SPNSFS,DIPSOm_SA)
       IFUNCT=IFUNCT+KDGN
  455 CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (KDGN.NE.2) THEN
          WRITE(6,*) 'no twofold degeneracy'
          GOTO 5780
      ENDIF


      IF ((ISS.EQ.1).AND.(KDGN.EQ.2).AND.(IPGLOB.GE.VERBOSE)) THEN
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
        IJSO=ISO+NSS*(JSO-1)
        JISO=JSO+NSS*(ISO-1)
        CONTRIB=WORK(IZMR(IXYZ)-1+IJSO)*WORK(IZMR(JXYZ)-1+JISO)
     &          -WORK(IZMI(IXYZ)-1+IJSO)*WORK(IZMI(JXYZ)-1+JISO)


        GTIJ=GTIJ+CONTRIB


        END DO
       END DO
       GTENS(IXYZ,JXYZ)=2.0D0*GTIJ
       END DO
      END DO

      IF(IPGLOB.GT.VERBOSE) THEN
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

      WRITE(6,*)" "

      ICOUNT = 0
      DO IXYZ=1,3
       DO JXYZ=1,3
       GTENS(IXYZ,JXYZ)=0.0D0
        DO KXYZ=1,3
        If (EVR(KXYZ).LT.THRSH) THEN ! EVR(KXYZ) should probably be ZERO
          IF (EVR(KXYZ).LT.-THRSH) THEN
           write(6,*) "WARNNING : negative real eigenvalue of G"
           write(6,*) "WARNNING : this may cause numerical problem"
           write(6,*) "WARNNING : and may give 'NAN' in G tensor"
          ELSE
           EVR(KXYZ)=ZERO
           If (ICOUNT.EQ.0) then
           write(6,*) " The eigenvalue 'EVR(KXYZ)' will be set to ZERO"
           write(6,*) " each time it is very small"
           write(6,*) " to avoid the roots of negative eigenvalues"
           endif
           ICOUNT=1
          ENDIF
        ENDIF
        GTENS(IXYZ,JXYZ)=GTENS(IXYZ,JXYZ)+
     &   TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
        END DO
       END DO
      ENDDO
      WRITE(6,*) " "
      WRITE(6,'(6x,3(5x,a2,5x),'//
     & '4x,4x,2x,8x,'//
     & '2x,2x,2x,3(4x,a2,i1,3x))')
     & (xyzchr(IXYZ),IXYZ=1,3), ('a_',IXYZ,IXYZ=1,3)
      WRITE(6,*)
      DO IXYZ=1,3
      WRITE(6,'(2x,a2,2x,3(1x,f10.6,1x),'//
     & '4x,a2,i1,a1,2x,f12.9,'//
     & '2x,a2,2x,3(1x,f8.4,1x))')
     & xyzchr(IXYZ), (GTENS(IXYZ,JXYZ),JXYZ=1,3),
     & 'a_',IXYZ,':',SQRT(EVR(IXYZ)),
     & xyzchr(IXYZ), (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
      ENDDO
       WRITE(6,*)
 5780  CONTINUE

      ISS=ISS+KDGN

      ENDDO

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)
 1904  CONTINUE

* End loop over CNT properties
      END IF
      END DO

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

 1801  CONTINUE
*****************************************************
*****************************************************
*****************************************************
* Experimental hyperfine tensor stuff ends here
*****************************************************
*****************************************************
*****************************************************



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

      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LLZI),1)

      IF(IAMX.GT.0)
     &     CALL SMMAT(PROP,WORK(LLXI),NSS,PNAME(IAMX),ICOMP(IAMX))
      IF(IAMY.GT.0)
     &    CALL SMMAT(PROP,WORK(LLYI),NSS,PNAME(IAMY),ICOMP(IAMY))
      IF(IAMZ.GT.0)
     &    CALL SMMAT(PROP,WORK(LLZI),NSS,PNAME(IAMZ),ICOMP(IAMZ))

      CALL GETMEM('MXR','ALLO','REAL',LMXR,NSS**2)
      CALL GETMEM('MXI','ALLO','REAL',LMXI,NSS**2)
      CALL GETMEM('MYR','ALLO','REAL',LMYR,NSS**2)
      CALL GETMEM('MYI','ALLO','REAL',LMYI,NSS**2)
      CALL GETMEM('MZR','ALLO','REAL',LMZR,NSS**2)
      CALL GETMEM('MZI','ALLO','REAL',LMZI,NSS**2)

      IMR(1)=LMXR
      IMI(1)=LMXI
      IMR(2)=LMYR
      IMI(2)=LMYI
      IMR(3)=LMZR
      IMI(3)=LMZI

      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMXR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMXI),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMYR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMYI),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMZR),1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LMZI),1)

      CALL SMMAT(PROP,WORK(LMXR),NSS,'SPIN    ',1)
      CALL SMMAT(PROP,WORK(LMYI),NSS,'SPIN    ',2)
      CALL SMMAT(PROP,WORK(LMZR),NSS,'SPIN    ',3)

      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMXR),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMYI),1)
      CALL DSCAL_(NSS**2,FEGVAL,WORK(LMZR),1)

      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LMXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LMYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LMZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMXR),WORK(LMXI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMYR),WORK(LMYI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LMZR),WORK(LMZI))

      CALL GETMEM('ZXR','ALLO','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','ALLO','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','ALLO','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','ALLO','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','ALLO','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','ALLO','REAL',LZZI,NSS**2)

      IZMR(1)=LZXR
      IZMI(1)=LZXI
      IZMR(2)=LZYR
      IZMI(2)=LZYI
      IZMR(3)=LZZR
      IZMI(3)=LZZI

      CALL GETMEM('LZR','ALLO','REAL',LZR,NSS**2)
      CALL GETMEM('LZI','ALLO','REAL',LZI,NSS**2)
      CALL GETMEM('UZR','ALLO','REAL',LUZR,NSS**2)
      CALL GETMEM('UZI','ALLO','REAL',LUZI,NSS**2)

      BFINAL=BSTART+(NBSTEP-1)*BINCRE
      TFINAL=TSTART+(NTSTEP-1)*TINCRE

      WRITE(6,*) "Magnetic flux density range (T): "
      WRITE(6,'(2x,f6.2,a3,f6.2,a4,i4,a6)')
     & BSTART," - ",BFINAL," in ",NBSTEP," steps"
      WRITE(6,*)
      WRITE(6,*) "Temperature range (K): "
      WRITE(6,'(2x,f6.2,a3,f6.2,a4,i4,a6)')
     & TSTART," - ",TFINAL," in ",NTSTEP," steps"

      CALL GETMEM('MAGM','ALLO','REAL',LMAGM,9*NBSTEP*NTSTEP)

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
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZR),1)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZI),1)
        CALL DAXPY_(NSS**2,0.5D0*B/AU2T,WORK(IMR(IXYZ)),1,WORK(LZR),1)
        CALL DAXPY_(NSS**2,0.5D0*B/AU2T,WORK(IMI(IXYZ)),1,WORK(LZI),1)
        DO ISS=1,NSS
         IISS=ISS+NSS*(ISS-1)
         HZER=WORK(LZR-1+IISS)
         WORK(LZR-1+IISS)=HZER+ENSOR(ISS)
        END DO
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUZR),1)
        CALL DCOPY_(NSS   ,1.0D0,0,WORK(LUZR),NSS+1)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUZI),1)
        CALL ZJAC(NSS,WORK(LZR),WORK(LZI),NSS,WORK(LUZR),WORK(LUZI))
        DO JXYZ=1,3
         CALL DCOPY_(NSS**2,WORK(IMR(JXYZ)),1,WORK(IZMR(JXYZ)),1)
         CALL DCOPY_(NSS**2,WORK(IMI(JXYZ)),1,WORK(IZMI(JXYZ)),1)
         CALL DSCAL_(NSS**2,-0.5D0,WORK(IZMR(JXYZ)),1)
         CALL DSCAL_(NSS**2,-0.5D0,WORK(IZMI(JXYZ)),1)
         CALL ZTRNSF(NSS,WORK(LUZR),WORK(LUZI),
     &    WORK(IZMR(JXYZ)),WORK(IZMI(JXYZ)))
        ENDDO
        DO ITSTEP=1,NTSTEP
         T=TSTART+TINCRE*(ITSTEP-1)
         RkT=T*BOLTZ
         RMAGM(1)=0.0D0
         RMAGM(2)=0.0D0
         RMAGM(3)=0.0D0
         RPART=0.0D0
         IF(IPGLOB.GT.USUAL) THEN
          WRITE(6,*)
          WRITE(6,'(2x,a14,3(4x,a4,4x),2x,a6)') "Energy (cm^-1)",
     &     "mu_x", "mu_y", "mu_z","weight"
          WRITE(6,*)
         ENDIF
         DO ISS=1,NSS
          IISS=ISS+NSS*(ISS-1)
          DELTA=WORK(LZR-1+IISS)-WORK(LZR)
          FACT=EXP(-DELTA/RkT)
          RMAGM(1)=RMAGM(1)+WORK(IZMR(1)-1+IISS)*FACT
          RMAGM(2)=RMAGM(2)+WORK(IZMR(2)-1+IISS)*FACT
          RMAGM(3)=RMAGM(3)+WORK(IZMR(3)-1+IISS)*FACT
          RPART=RPART+FACT
          IF(IPGLOB.GT.USUAL) THEN
           WRITE(6,'(2x,f14.3,3(1x,f10.6,1x),2x,f6.3)')
     &      (WORK(LZR-1+IISS)-WORK(LZR))*AU2CM,
     &      WORK(IZMR(1)-1+IISS),WORK(IZMR(2)-1+IISS),
     &      WORK(IZMR(3)-1+IISS),FACT
          ENDIF
         ENDDO
         IF(IPGLOB.GT.USUAL) THEN
          WRITE(6,*)
         ENDIF
         RMAGM(1)=(RMAGM(1)/RPART)*AU2JTM
         RMAGM(2)=(RMAGM(2)/RPART)*AU2JTM
         RMAGM(3)=(RMAGM(3)/RPART)*AU2JTM
         RMAGM2=RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
         RMAGMO=SQRT(RMAGM2)
         DO JXYZ=1,3
          LMSTEP=LMSTEP+1
          WORK(LMAGM-1+LMSTEP)=RMAGM(JXYZ)
          IF(IBSTEP.GT.1) THEN
              Chi(JXYZ)=RMAGM(JXYZ)-WORK(LMAGM-1+LMSTEP-3*NTSTEP)
              Chi(JXYZ)=Chi(JXYZ)*Rmu0/BINCRE
          ENDIF
         ENDDO
          IF(IBSTEP.EQ.1) THEN
         WRITE(6,'(1x,f6.2,5(1x,e12.5,1x))')
     &    T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
          ELSE
         WRITE(6,'(1x,f6.2,8(1x,e12.5,1x))')
     &    T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3),Chi(1),Chi(2),Chi(3)
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      CALL GETMEM('MAGM','FREE','REAL',LMAGM,9*NBSTEP*NTSTEP)

      WRITE(6,*)

C powder magnetization, useful in nonlinear cases

      IF(.NOT.IFMCAL) GOTO 810

      WRITE(6,*)
      WRITE(6,*) "Powder Magnetization"
      WRITE(6,*)
      WRITE(6,'(3x,a1,3x,5(1x,a12,1x))') "T",
     & "    B  (T)  ","   M (J/T)  ",
     & "  Mx (J/T)  ","  My (J/T)  ","  Mz (J/T)  "

      CALL GETMEM('MAGM','ALLO','REAL',LMAGM,3*NBSTEP*NTSTEP)
      CALL DCOPY_(3*NBSTEP*NTSTEP,0.0D0,0,WORK(LMAGM),1)

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
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZR),1)
        CALL DAXPY_(NSS**2,0.5D0*BX/AU2T,WORK(LMXR),1,WORK(LZR),1)
        CALL DAXPY_(NSS**2,0.5D0*BY/AU2T,WORK(LMYR),1,WORK(LZR),1)
        CALL DAXPY_(NSS**2,0.5D0*BZ/AU2T,WORK(LMZR),1,WORK(LZR),1)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LZI),1)
        CALL DAXPY_(NSS**2,0.5D0*BX/AU2T,WORK(LMXI),1,WORK(LZI),1)
        CALL DAXPY_(NSS**2,0.5D0*BY/AU2T,WORK(LMYI),1,WORK(LZI),1)
        CALL DAXPY_(NSS**2,0.5D0*BZ/AU2T,WORK(LMZI),1,WORK(LZI),1)
        DO ISS=1,NSS
         IISS=ISS+NSS*(ISS-1)
         HZER=WORK(LZR-1+IISS)
         WORK(LZR-1+IISS)=HZER+ENSOR(ISS)
        END DO
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUZR),1)
        CALL DCOPY_(NSS   ,1.0D0,0,WORK(LUZR),NSS+1)
        CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUZI),1)
        CALL ZJAC(NSS,WORK(LZR),WORK(LZI),NSS,WORK(LUZR),WORK(LUZI))
        DO IXYZ=1,3
         CALL DCOPY_(NSS**2,WORK(IMR(IXYZ)),1,WORK(IZMR(IXYZ)),1)
         CALL DCOPY_(NSS**2,WORK(IMI(IXYZ)),1,WORK(IZMI(IXYZ)),1)
         CALL DSCAL_(NSS**2,-0.5D0,WORK(IZMR(IXYZ)),1)
         CALL DSCAL_(NSS**2,-0.5D0,WORK(IZMI(IXYZ)),1)
         CALL ZTRNSF(NSS,WORK(LUZR),WORK(LUZI),
     &    WORK(IZMR(IXYZ)),WORK(IZMI(IXYZ)))
        ENDDO
        DO ITSTEP=1,NTSTEP
         T=TSTART+TINCRE*(ITSTEP-1)
         RkT=T*BOLTZ
         RMAGM(1)=0.0D0
         RMAGM(2)=0.0D0
         RMAGM(3)=0.0D0
         RPART=0.0D0
         DO ISS=1,NSS
          IISS=ISS+NSS*(ISS-1)
          DELTA=WORK(LZR-1+IISS)-WORK(LZR)
          FACT=EXP(-DELTA/RkT)
          RMAGM(1)=RMAGM(1)+WORK(IZMR(1)-1+IISS)*FACT
          RMAGM(2)=RMAGM(2)+WORK(IZMR(2)-1+IISS)*FACT
          RMAGM(3)=RMAGM(3)+WORK(IZMR(3)-1+IISS)*FACT
          RPART=RPART+FACT
         ENDDO
         RMAGM(1)=(RMAGM(1)/RPART)*AU2JTM
         RMAGM(2)=(RMAGM(2)/RPART)*AU2JTM
         RMAGM(3)=(RMAGM(3)/RPART)*AU2JTM
*        WRITE(6,'(6(1x,e12.5,1x))')
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
          WORK(LMAGM-1+LMSTEP)=WORK(LMAGM-1+LMSTEP)+RMAGM(IXYZ)
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
          RMAGM(IXYZ)=WORK(LMAGM-1+LMSTEP)/NORIENT
        ENDDO
        RMAGM2=RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
        RMAGMO=SQRT(RMAGM2)
        WRITE(6,'(1x,f6.2,5(1x,e12.5,1x))')
     &   T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
       ENDDO
      ENDDO

      CALL GETMEM('MAGM','FREE','REAL',LMAGM,3*NBSTEP*NTSTEP)

 810  CONTINUE

      WRITE(6,*)

      CALL GETMEM('LZR','FREE','REAL',LZR,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LZI,NSS**2)
      CALL GETMEM('UZR','FREE','REAL',LUZR,NSS**2)
      CALL GETMEM('UZI','FREE','REAL',LUZI,NSS**2)

      CALL GETMEM('ZXR','FREE','REAL',LZXR,NSS**2)
      CALL GETMEM('ZXI','FREE','REAL',LZXI,NSS**2)
      CALL GETMEM('ZYR','FREE','REAL',LZYR,NSS**2)
      CALL GETMEM('ZYI','FREE','REAL',LZYI,NSS**2)
      CALL GETMEM('ZZR','FREE','REAL',LZZR,NSS**2)
      CALL GETMEM('ZZI','FREE','REAL',LZZI,NSS**2)

      CALL GETMEM('MXR','FREE','REAL',LMXR,NSS**2)
      CALL GETMEM('MXI','FREE','REAL',LMXI,NSS**2)
      CALL GETMEM('MYR','FREE','REAL',LMYR,NSS**2)
      CALL GETMEM('MYI','FREE','REAL',LMYI,NSS**2)
      CALL GETMEM('MZR','FREE','REAL',LMZR,NSS**2)
      CALL GETMEM('MZI','FREE','REAL',LMZI,NSS**2)

 900  CONTINUE

      RETURN
      END

      SUBROUTINE SINANI(KDGN,IFUNCT,NSS,DIPSOm,SPNSFS,DIPSOm_SA)
!      IMPLICIT NONE
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER KDGN,IFUNCT,NSS,l,Iso1,Jso2,Ico1,i,j
      COMPLEX*16 DIPSOm(3,NSS,NSS),DIPSOmSA(3,KDGN,KDGN)
      COMPLEX*16 SPNSOSA(3,KDGN,KDGN)
      COMPLEX*16 Z(NSS,NSS),MATL(NSS,NSS),FINL(NSS,NSS)
      COMPLEX*16 SPNSO(3,NSS,NSS),SPNSFS(3,NSS,NSS)
      real*8 UMATR(NSS,NSS),UMATI(NSS,NSS),gtens(3),maxes(3,3)
      CHARACTER*1 angm

      if(.FALSE.) then
      write(6,'(/)')
      write(6,'(10A)') (('############'),J=1,10)
      write(6,'(25X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT IN '//
     &'THE BASIS OF SPIN ORBIT STATES'
      write(6,'(10A)') (('############'),J=1,10)

      do l=1,3
      if(l.eq.1)  angm='X'
      if(l.eq.2)  angm='Y'
      if(l.eq.3)  angm='Z'
      write(6,'(/)')
      write(6,'(4X,A12,A2)') 'PROJECTION: ', angm
      write(6,'(/)')
       do Iso1=1,NSS
      write(6,'(20(2X,2F10.6))') (DIPSOm(l,Iso1,Jso2), Jso2=1,NSS)
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
      write(6,*)
      write(6,'(10X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT '//
     & 'IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
      do l=1,3
      write(6,'(/)')
      write(6,'(5X,A6,I3)') 'AXIS= ',l
      write(6,*)
       do Ico1=1,KDGN
      write(6,'(16(2F12.8,2x))') (DIPSOmSA(l,Ico1,Jco1), Jco1=1,KDGN)
       enddo
      enddo

       endif ! if(IPGLOB.GE.4)

          CALL ATENS(DIPSOmSA, KDGN, gtens, maxes, 2)

           if(.False.) then

          call get_dArray('UMATR_SINGLE',UMATR,NSS**2)
          call get_dArray('UMATI_SINGLE',UMATI,NSS**2)
       write(6,'(/)')
       write(6,'(5x,a)') 'umatr and umati'
       write(6,'(/)')
       do i=1,NSS
       write(6,'(5x,10(2f14.10,2x))') (UMATR(i,j),UMATI(i,j),j=1,NSS)
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
        Z(i,j)=Z(i,j)+DCMPLX(UMATR(i,j),UMATI(i,j))
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

      write(6,'(/)')
      write(6,'(10A)') (('############'),J=1,10)
      write(6,'(30X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT IN '//
     & 'THE BASIS OF SPIN ORBIT STATES'
      write(6,'(10A)') (('############'),J=1,10)
      write(6,'(/)')
      do l=1,3
      if(l.eq.1)  angm='X'
      if(l.eq.2)  angm='Y'
      if(l.eq.3)  angm='Z'
      write(6,'(/)')
      write(6,'(4X,A,A)') 'PROJECTION: ', angm
      write(6,'(/)')
       do Iso1=1,NSS
      write(6,'(20(2F10.6,2X))') (SPNSO(l,Iso1,Jso2), Jso2=1,NSS)
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

      write(6,*)
      write(6,'(10X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT '//
     & 'IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
      do l=1,3
      write(6,'(/)')
      write(6,'(5X,A6,I3)') 'AXIS= ',l
      write(6,*)
       do Ico1=1,KDGN
      write(6,'(16(2F12.8,2x))') (SPNSOSA(l,Ico1,Jco1), Jco1=1,KDGN)
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


      RETURN
c Avoid unused argument warnings
      IF (.FALSE.)  CALL Unused_real(DIPSOm_SA)
      END

      SUBROUTINE ADARASSI(N,A,D,DROT)

      IMPLICIT NONE
      INTEGER I, J,  N
      COMPLEX*16  A(N,N), D(N,N), DROT(N,N), TEMP(N,N)

C initialization
      do I=1,N
       do J=1,N
      DROT(I,J)=(0.0D0,0.0D0)
      TEMP(I,J)=(0.0D0,0.0D0)
       enddo
      enddo

C actual multiplication
      call ZGEMM('C','N',N,N,N,(1.0D0,0.0D0),A,N,D,N,(0.0D0,0.0D0),
     &TEMP,N)
      call ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),TEMP,N,A,N,(0.0D0,0.0D0),
     &DROT,N)

      RETURN
      END

      SUBROUTINE ZECON(NSTATE,N,UR,UI,AR,AI,ZEKL,IXYZ,ISTATE,ISS,JSS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UR(N,N),UI(N,N)
      DIMENSION AR(N,N),AI(N,N)
      COMPLEX*16 ZEKL(2,2,3,NSTATE)
#include "WrkSpc.fh"

      TMPR1=AR(ISS,JSS)*UR(JSS,1)-AI(ISS,JSS)*UI(JSS,1)
      TMPR2=AR(ISS,JSS)*UR(JSS,2)-AI(ISS,JSS)*UI(JSS,2)
      TMPI1=AI(ISS,JSS)*UR(JSS,1)+AR(ISS,JSS)*UI(JSS,1)
      TMPI2=AI(ISS,JSS)*UR(JSS,2)+AR(ISS,JSS)*UI(JSS,2)
      ZEKL(1,1,IXYZ,ISTATE)=ZEKL(1,1,IXYZ,ISTATE)+
     $     DCMPLX(UR(ISS,1)*TMPR1+UI(ISS,1)*TMPI1,
     $     UR(ISS,1)*TMPI1-UI(ISS,1)*TMPR1)
      ZEKL(1,2,IXYZ,ISTATE)=ZEKL(1,2,IXYZ,ISTATE)+
     $     DCMPLX(UR(ISS,1)*TMPR2+UI(ISS,1)*TMPI2,
     $     UR(ISS,1)*TMPI2-UI(ISS,1)*TMPR2)
      ZEKL(2,1,IXYZ,ISTATE)=ZEKL(2,1,IXYZ,ISTATE)+
     $     DCMPLX(UR(ISS,2)*TMPR1+UI(ISS,2)*TMPI1,
     $     UR(ISS,2)*TMPI1-UI(ISS,2)*TMPR1)
      ZEKL(2,2,IXYZ,ISTATE)=ZEKL(2,2,IXYZ,ISTATE)+
     $     DCMPLX(UR(ISS,2)*TMPR2+UI(ISS,2)*TMPI2,
     $     UR(ISS,2)*TMPI2-UI(ISS,2)*TMPR2)

      RETURN
      END
