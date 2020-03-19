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
* Copyright (C) 2018,2019, Roland Lindh                                *
*               2019, Mickael G. Delcey                                *
*               2019, Ignacio Fdez. Galvan                             *
************************************************************************
      SUBROUTINE PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,JBNUM,EigVec)
      USE RASSI_AUX
      USE kVectors
#include "compiler_features.h"
#ifndef POINTER_REMAP
      USE ISO_C_Binding
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='PRPROP_TM')
      parameter (THRSH=1.0D-10)
      parameter (ZERO=0.0D0)
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "constants.fh"
#include "stdalloc.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP),JBNUM(NSTATE),
     &          EigVec(NSTATE,NSTATE)
#include "SysDef.fh"
#include "rassiwfn.fh"
      LOGICAL TMOgroup
      INTEGER IOFF(8),IJSS(4),IPRTMOM(14)
      CHARACTER*8 LABEL
      CHARACTER*6 STLNE1
      CHARACTER*52 STLNE2
      Integer, Dimension(:), Allocatable :: TMOgrp1,TMOgrp2,ISS_INDEX,
     &   iMask,jMask,iSSMask,jSSMask
      Real*8 TM_R(3), TM_I(3), TM_C(3)
      Real*8 wavevector(3), UK(3)
      Real*8 kPhase(2)
      Real*8, Allocatable :: pol_Vector(:,:)
#ifdef _HDF5_
      Real*8, Allocatable, Target :: Storage(:,:,:,:)
      Real*8, Pointer :: flatStorage(:)
#endif
      Real*8, Allocatable:: TDMZZ(:),TSDMZZ(:),WDMZZ(:), SCR(:,:)
      Real*8, Allocatable:: VSOR(:,:), VSOI(:,:), TMP(:)

      CALL QENTER(ROUTINE)
#define _TIME_TMOM_
#ifdef _TIME_TMOM_
      Call CWTime(TCpu1,TWall1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Before we start we need to backtransform the coefficients of the
*     SO states from the basis of the SF states which diagonalize the
*     SF Hamiltonian to the basis of the original (input) SF states.
*     This since all transition moments, whether retrieved from disk or
*     recomputed, are in the basis of the original SF states.
*
      Call mma_allocate(VSOR,NSS,NSS,Label='VSOR')
      Call mma_allocate(VSOI,NSS,NSS,Label='VSOI')
*
      Call USOTRANS(USOR,USOI,NSS,
     &              EigVec,NSTATE,
     &              VSOR,VSOI)
*
*                                                                      *
************************************************************************
*                                                                      *
      DEBYE=CONV_AU_TO_DEBYE_
      AU2REDR=2.0D2*DEBYE
      ! AFACTOR = 2*pi*e^2*E_h^2 / eps_0*m_e*c^3*h^2
      ! 1/c^3 (in a.u. of time ^ -1)
      AFACTOR = 2.0D0/CONST_C_IN_AU_**3
     &          /CONST_AU_TIME_IN_SI_
      HALF=0.5D0
      PI= CONST_PI_
      SPEED_OF_LIGHT=CONST_C_IN_AU_
      G_Elec=CONST_ELECTRON_G_FACTOR_

*define _TIME_TMOM_
#ifdef _TIME_TMOM_
      Call CWTime(TCpu1,TWall1)
#endif
C Compute transition strengths for spin-orbit states:
*
* Initial setup for exact operator
*
C printing threshold
      OSTHR=1.0D-5
      IF(DIPR) OSTHR = OSTHR_DIPR
      IF(DIPR) WRITE(6,*) ' Threshold changed to ',OSTHR
! Again to avoid total negative transition strengths
      IF(QIPR) OSTHR = OSTHR_QIPR
      IF(QIPR) WRITE(6,*) ' Threshold changed to ',OSTHR,
     &                    ' since quadrupole threshold is given '
!
!     Reducing the loop over states - good for X-rays
!     At the moment memory is not reduced
!
      IF(REDUCELOOP) THEN
        IEND = LOOPDIVIDE
        JSTART = LOOPDIVIDE+1
      ELSE
        IEND = NSS
        JSTART = 1
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
*     TMOM  RS: The symmetric part of the real comp. of the op.
*     TMOM  RA: The asymmetric part of the real comp. of the op.
*     TMOM  IS: The symmetric part of the imaginary comp. of the op.
*     TMOM  IA: The asymmetric part of the imaginary comp. of the op.
*
************************************************************************
*                                                                      *
*     Section (1): Computation of k specific oscillator strength.      *
*     Section (2): Computation of the isotropic oscillator strength.   *
*                                                                      *
************************************************************************
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
         IF (PNAME(IPROP).EQ.'TMOM0  R') THEN
            IF (IPRTMOM(13).EQ.-1) IPRTMOM(13)=IPROP
         END IF
         IF (PNAME(IPROP).EQ.'TMOM0  I') THEN
            IF (IPRTMOM(14).EQ.-1) IPRTMOM(14)=IPROP
         END IF
      ENDDO
      IF (ANY(IPRTMOM.EQ.-1)) Go To 666
*
*     Initiate the Seward environment
*
      nDiff=0
      Call IniSew(Info,.FALSE.,nDiff)
*
*     Generate the quadrature points.
*
*     In the spin-coupled case, wave functions are complex and there is
*     not a simple relation between oscillator and rotatory strengths for
*     k and -k vectors, but there is between the integrals computed in
*     the spin-free basis, so we compute them only once and obtain separately
*     the results for k and -k, given by kPhase
*
      kPhase = [1.0D0, -1.0D0]
      If (Do_SK) Then
         nQuad = 1
         nVec=nk_Vector
         Call GetMem('SK','ALLO','REAL',ipR,4*nQuad)
         If (.Not.(PRRAW.Or.PRWEIGHT)) kPhase(2) = 0.0D0
      Else
         Call Setup_O()
         Call Do_Lebedev_Sym(L_Eff,nQuad,ipR)
         nVec = 1
      End If
      If (Do_Pol) Call mma_allocate(pol_Vector,3,nVec*nQuad,Label='POL')
*
*     Initialize for density matrices.
*
      Call mma_Allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
      Call mma_Allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
      Call mma_Allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
      nSCR=(NBST*(NBST+1))/2
      Call mma_allocate(SCR,nSCR,4,LABEL='SCR')

*
*     Allocate some temporary arrays for handling the
*     properties of the spin-orbit states.
*
      CALL GETMEM('DXR','ALLO','REAL',LDXR,3*NSS**2)
      CALL GETMEM('DXI','ALLO','REAL',LDXI,3*NSS**2)
      CALL GETMEM('DXRM','ALLO','REAL',LDXRM,3*NSS**2)
      CALL GETMEM('DXIM','ALLO','REAL',LDXIM,3*NSS**2)
      Call mma_Allocate(TMP,NSS**2,Label='TMP')
      CALL GETMEM('TMR','ALLO','REAL',LTMR,3*NSS**2)
      CALL GETMEM('TMI','ALLO','REAL',LTMI,3*NSS**2)
*
C     ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
      NIP=4+(NBST*(NBST+1))/2
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
#ifdef _HDF5_
*
*     Allocate vector to store all individual transition moments.
*     We do this for
*     all unique pairs ISO-JSO, iSO=/=JSO (NSS*(NSS-1)/2)
*         all k-vectors (2*nQuad or 2*nVec)
*             we store:
*                 the weight (1)
*                 the k-vector (3)
*                 we projected transition vector (real and imaginary parts) (2*3)
*
      nIJ=nSS*(nSS-1)/2
      ip_w       = 1
      ip_kvector = ip_w + 1
      ip_TMR     = ip_kvector + 3
      ip_TMI     = ip_TMR + 3
      nData      = ip_TMI + 3 - 1
      Call mma_allocate(Storage,nData,2*nQuad,nIJ,nVec,label='Storage')
      Call dCopy_(Size(Storage),[0.0D0],0,Storage,1)
#endif
*MGD group transitions
      TMOgroup=.false.
      ngroup1=IEND
      ngroup2=NSS-JSTART+1
      nmax2=1
      IF(REDUCELOOP.and.TMGr_thrs.ge.0.0d0) THEN
        TMOgroup=.true.
        THRS=TMGr_thrs
        i=IEND
        RefEne=0
        TAU=-1
        ngroup2=1
        Do j=JSTART,NSS
          if (ENSOR(J)-Refene.gt.TAU) then
              NGROUP2=NGROUP2+1
              Refene=ENSOR(J)
              ediff=Refene-ENSOR(I)
              TAU=ediff*THRS
           EndIf
        End Do
        Call mma_Allocate(TMOgrp2,NGROUP2,Label='TMOgrp2')
        ngroup2=0
        TAU=-1
        RefEne=0
        Do j=JSTART,NSS
          if (ENSOR(J)-Refene.gt.TAU) then
              NGROUP2=NGROUP2+1
              TMOgrp2(NGROUP2)=J
              Refene=ENSOR(J)
              ediff=Refene-ENSOR(I)
              TAU=ediff*THRS
           EndIf
        End Do
        TMOgrp2(ngroup2+1)=NSS+1
*
        j=JSTART
        Refene=ENSOR(j)
        TAU=-1
        ngroup1=1
        Do i=IEND,1,-1
          if (Refene-ENSOR(i).gt.TAU) then
            ngroup1=ngroup1+1
            Refene=ENSOR(i)
            ediff=ENSOR(j)-Refene
            Tau=ediff*THRS
          EndIf
        End Do
        Call mma_Allocate(TMOgrp1,NGROUP1,Label='TMOgrp1')
        Ntmp=Ngroup1
        Ngroup1=Ngroup1-1
        Refene=ENSOR(j)
        TAU=-1
        Do i=IEND,1,-1
          if (Refene-ENSOR(i).gt.TAU) then
            TMOgrp1(ntmp)=i+1
            ntmp=ntmp-1
            Refene=ENSOR(i)
            ediff=ENSOR(j)-Refene
            Tau=ediff*THRS
          EndIf
        End Do
        TMOgrp1(1)=1
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
      CALL GETMEM('RAW   ','ALLO','REAL',LRAW,2*NQUAD*6*nmax2)
      CALL GETMEM('OSCSTR','ALLO','REAL',LF,2*nmax2)
      CALL GETMEM('MAXMIN','ALLO','REAL',LMAX,8*nmax2)
      LMAX_=0
*
      Do iVec = 1, nVec
*
         If (Do_SK) Then
            Work(ipR  )=k_Vector(1,iVec)
            Work(ipR+1)=k_Vector(2,iVec)
            Work(ipR+2)=k_Vector(3,iVec)
            Work(ipR+3)=1.0D0   ! Dummy weight
         End If
*
      iPrint=0
      IJSO=0
*
      ThrSparse=1.0D-12
*
      Call mma_allocate(iMask,NState,LABEL='iMask')
      Call mma_allocate(jMask,NState,LABEL='jMask')
      Call mma_allocate(iSSMask,NSS,LABEL='iSSMask')
      Call mma_allocate(jSSMask,NSS,LABEL='jSSMask')
*
      CALL mma_allocate(ISS_INDEX,NState+1,LABEL='ISS_INDEX')
      ISS_INDEX(1)=0
      Do iState=1,NState
         Job=JBNUM(iState)
         ISS_INDEX(iState+1)=ISS_INDEX(IState)+MLTPLT(Job)
      End Do
*
      Do igrp=1,ngroup1
*
         If (TMOgroup) Then
            istart_=TMOgrp1(igrp)
            iend_=TMOgrp1(igrp+1)-1
            ENSOR1=0.5D0*(ENSOR(istart_)+ENSOR(iend_))
         Else
            istart_=igrp
            iend_=igrp
            ENSOR1=ENSOR(istart_)
         EndIf
*Screening
         Call iCopy(NState,[0],0,iMask,1)
         Call iCopy(NSS,[0],0,iSSMask,1)
         ISM=0
         ISSM=0
         ISFLoop: Do ISF=1,NState
           Do ISS=ISS_INDEX(ISF)+1,ISS_INDEX(ISF+1)
             Do ISO=istart_,iend_
               Temp=VSOR(ISS,ISO)**2 + VSOI(ISS,ISO)**2
               If (Temp.gt.ThrSparse) Then
                 ISM=ISM+1
                 iMask(ISM)=ISF
                 Do IMSS=ISS_INDEX(ISF)+1,ISS_INDEX(ISF+1)
                   ISSM=ISSM+1
                   iSSMask(ISSM)=IMSS
                 End Do
                 Cycle ISFLoop
               End If
             End Do
           End Do
         End Do ISFLoop
*
         Do jgrp=1,ngroup2
*
            If (TMOgroup) Then
              jstart_=TMOgrp2(jgrp)
              jend_=TMOgrp2(jgrp+1)-1
              ENSOR2=0.5D0*(ENSOR(jstart_)+ENSOR(jend_))
            Else
              jstart_=jgrp+jstart-1
              jend_=jgrp+jstart-1
              ENSOR2=ENSOR(jstart_)
            EndIf
            EDIFF_=ENSOR2-ENSOR1
            n12=(iend_-istart_+1)*(jend_-jstart_+1)
*Screening
            Call iCopy(NState,[0],0,jMask,1)
            Call iCopy(NSS,[0],0,jSSMask,1)
            JSM=0
            JSSM=0
            JSFLoop: Do JSF=1,NState
              Do JSS=ISS_INDEX(JSF)+1,ISS_INDEX(JSF+1)
                Do JSO=jstart_,jend_
                  Temp=VSOR(JSS,JSO)**2 + VSOI(JSS,JSO)**2
                  If (Temp.gt.ThrSparse) Then
                    JSM=JSM+1
                    jMask(JSM)=JSF
                    Do JMSS=ISS_INDEX(JSF)+1,ISS_INDEX(JSF+1)
                      JSSM=JSSM+1
                      jSSMask(JSSM)=JMSS
                    End Do
                    Cycle JSFLoop
                  End If
                End Do
              End Do
            End Do JSFLoop
*
            IJSS(1)=istart_
            IJSS(2)=iend_
            IJSS(3)=jstart_
            IJSS(4)=jend_
*
            WRITE(STLNE1,'(A6)') 'RASSI:'
            WRITE(STLNE2,'(A33,I5,A5,I5)')
     &         'Trans. intensities for SO groups ',igrp,' and ',jgrp
            Call StatusLine(STLNE1,STLNE2)
*
            IF (ABS(EDIFF_).LE.1.0D-8) CYCLE
            IF(EDIFF_.LT.0.0D0) CYCLE
            IJSO=IJSO+1
*
*           The energy difference is used to define the norm of the
*           wave vector.
*
            rkNorm=ABS(EDIFF_)/SPEED_OF_LIGHT
*
*           For the case the energy difference is negative we
*           need to change the sign of the B.s term to get
*           consistency.
*
            cst=SIGN(Half,EDIFF_)*g_Elec
*
*           Iterate over the quadrature points.
*
*           Initialize output arrays
*
            CALL DCOPY_(2*n12,[0.0D0],0,WORK(LF),1)
            CALL DCOPY_(2*NQUAD*6*n12,[0.0D0],0,WORK(LRAW),1)
            CALL DCOPY_(8*n12,[0.0D0],0,WORK(LMAX),1)
*
            Do iQuad = 1, nQuad
               iVec_=(iVec-1)*nQuad+iQuad
*
*              Generate the wavevector associated with this quadrature
*              point and pick up the associated quadrature weight.
*
               UK(1)=Work((iQuad-1)*4  +ipR)
               UK(2)=Work((iQuad-1)*4+1+ipR)
               UK(3)=Work((iQuad-1)*4+2+ipR)
               wavevector(:)=rkNorm*UK(:)
*
*              Note that the weights are normalized to integrate to
*              4*pi over the solid angles.
*
               Weight=Work((iQuad-1)*4+3+ipR)
               If (.Not.Do_SK) Weight = Weight/(4.0D0*PI)
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
               iOpt=2
               Call TMOMInt(wavevector,iOpt)
*
************************************************************************
*                                                                      *
*              Recompute the needed properties for all the spin-free   *
*              states.                                                 *
*                                                                      *
************************************************************************
*
               Do IPRP = 1,14
                  IPROP = IPRTMOM(IPRP)
                  Call FZero(PROP(1,1,IPROP),NSTATE**2)
               End Do

               DO ISS=1,ISM
*
*                 Does this spin-free state contribute to any of the
*                 two spin states? Check the corresponding coefficients.
*
                  i=iMask(ISS)
                  DO JSS=1,JSM
                     j=jMask(JSS)
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
                     Call mk_IOFF(IOFF,nSYM,NBASF,ISY12)
*
*                    Pick up the transition density between the two
*                    states from disc. Generated in GTDMCTL.
*
                     IDISK=iDisk_TDM(I,J,1)
                     iEmpty=iDisk_TDM(I,J,2)
                     iOpt=2
                     iGO=5
                     CALL dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,
     &                              LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
                     Call MK_TWDM(nSym,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,
     &                            IOFF,NBASF,ISY12)
*
*                    Compute the transition property of the property
*                    integrals between the two states.
*
                     DO IPRP = 1,14
                        IPROP = IPRTMOM(IPRP)
                        ITYPE=0
                        IF (PTYPE(IPROP).EQ.'HERMSING') ITYPE=1
                        IF (PTYPE(IPROP).EQ.'ANTISING') ITYPE=2
                        IF (PTYPE(IPROP).EQ.'HERMTRIP') ITYPE=3
                        IF (PTYPE(IPROP).EQ.'ANTITRIP') ITYPE=4
                        LABEL=PNAME(IPROP)
                        Call MK_PROP(PROP,IPROP,I,J,
     &                               LABEL,ITYPE,
     &                               WORK(LIP),NIP,SCR,nSCR,
     &                               MASK,ISY12,IOFF)
                     END DO
*
                  END DO ! J
               END DO ! I
*
*              LDXR & LDXI hold the component that does not change from k to -k
*              LDXRM & LDXIM hold the component that changes sign
*
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LDXR),1)
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LDXI),1)
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LDXRM),1)
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LDXIM),1)
               DO iCar = 1, 3
*
*              (1) the spin-free part.
*                 Note that the integrals are not computed for the
*                 momentum operator but rather for nabla. That is
*                 as we assemble to transition momentum we have to
*                 remember to put in a factor of -i.
*
                  LI_=LDXI+(iCar-1)*NSS**2
                  LRM_=LDXRM+(iCar-1)*NSS**2
*
*                 The real part (symmetric and anti-symmetric) becomes imaginary
*
                  CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(0+iCar),
     &                        0,ISS_INDEX,iMask,ISM,jMask,JSM)
                  CALL DAXPY_(NSS**2,-1.0D0,TMP,1,WORK(LI_),1)
                  CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(6+iCar),
     &                        0,ISS_INDEX,iMask,ISM,jMask,JSM)
                  CALL DAXPY_(NSS**2,-1.0D0,TMP,1,WORK(LI_),1)
*
*                 The imaginary part (symmetric and anti-symmetric) becomes real
*
                  CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(3+iCar),
     &                        0,ISS_INDEX,iMask,ISM,jMask,JSM)
                  CALL DAXPY_(NSS**2, 1.0D0,TMP,1,WORK(LRM_),1)
                  CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(9+iCar),
     &                        0,ISS_INDEX,iMask,ISM,jMask,JSM)
                  CALL DAXPY_(NSS**2, 1.0D0,TMP,1,WORK(LRM_),1)
*
*              (2) the spin-dependent part, magnetic
*
*                 iCar=1: 1/2(S(+)+S(-))
*                 iCar=2: 1/2i(S(+)-S(-))
*                 iCar=3: Sz
*
*                 Here we need to be very careful. The operator is
*                 similar to the L.S operator, the spin-orbit coupling.
*                 However, here we have that the operator is B.S. Thus
*                 we can do this in a similar fashion as the spin-orbit
*                 coupling, but with some important difference.
*
*                 For the spin-orbit coupling the L operator is divided
*                 into x-, y-, and z-components, that is the operator
*                 will be represented by three different integrals.
*                 For the integrals over B it is similar but the x-, y-,
*                 and z-components are constants outside of the integral.
*                 For example, the x-component is expressed as
*                 (k x e_l)_x <0|e^(i k.r)|n>. In this section we will
*                 handle the (k x e_l)_x part outside the loop over the
*                 Cartesian components.
*
*                 We have to note one further difference, the integrals
*                 are complex in this case.
*
*                 Let us now compute the contributions T(i), i=x,y,z
*
*                 In the equations we find that the y-component is imaginary,
*                 see the equations on page 234 in Per Ake's paper. Hence,
*                 the y-component is treated slightly differently.
*
*                 Actually, here we compute (<0|e^(i k.r)|n> x k), which
*                 will be later dotted with e_l.
*
*                 i*g/2*(s_y*k_z-s_z*k_y) -> T_x
*                 i*g/2*(s_z*k_x-s_x*k_z) -> T_y
*                 i*g/2*(s_x*k_y-s_y*k_x) -> T_z
*
*                 Note also that the "I" parts should change sign with kPhase,
*                 but so does wavevector, so the net result is the opposite.
*
                  IF (iCar.EQ.1) THEN
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),
     &                           iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2, wavevector(2)*cst,
     &                           TMP,1,WORK(LDXIM+2*NSS**2),1)
                     CALL DAXPY_(NSS**2,-wavevector(3)*cst,
     &                           TMP,1,WORK(LDXIM+1*NSS**2),1)
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),
     &                           iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2,-wavevector(2)*cst,
     &                           TMP,1,WORK(LDXR+2*NSS**2),1)
                     CALL DAXPY_(NSS**2, wavevector(3)*cst,
     &                           TMP,1,WORK(LDXR+1*NSS**2),1)
                  ELSE IF (iCar.EQ.2) THEN
*                    For the y-component we have to interchange the real and
*                    the imaginary components. The real component gets a
*                    minus sign due to the product i*i=-1
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),
     &                               iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2,-wavevector(3)*cst,
     &                           TMP,1,WORK(LDXI+0*NSS**2),1)
                     CALL DAXPY_(NSS**2, wavevector(1)*cst,
     &                           TMP,1,WORK(LDXI+2*NSS**2),1)
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),
     &                           iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2,-wavevector(3)*cst,
     &                           TMP,1,WORK(LDXRM+0*NSS**2),1)
                     CALL DAXPY_(NSS**2, wavevector(1)*cst,
     &                           TMP,1,WORK(LDXRM+2*NSS**2),1)
                  ELSE IF (iCar.EQ.3) THEN
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),
     &                           iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2, wavevector(1)*cst,
     &                           TMP,1,WORK(LDXIM+1*NSS**2),1)
                     CALL DAXPY_(NSS**2,-wavevector(2)*cst,
     &                           TMP,1,WORK(LDXIM+0*NSS**2),1)
                     CALL SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),
     &                           iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
                     CALL DAXPY_(NSS**2,-wavevector(1)*cst,
     &                           TMP,1,WORK(LDXR+1*NSS**2),1)
                     CALL DAXPY_(NSS**2, wavevector(2)*cst,
     &                           TMP,1,WORK(LDXR+0*NSS**2),1)
                  END IF
               END DO
*
*              Now we can compute the transition moments for k and -k
*
               DO kp = 1, 2
*
               IF (ABS(kPhase(kp)).LT.0.5d0) CYCLE
               CALL DCOPY_(3*NSS**2,WORK(LDXR),1,WORK(LTMR),1)
               CALL DCOPY_(3*NSS**2,WORK(LDXI),1,WORK(LTMI),1)
               CALL DAXPY_(3*NSS**2,kPhase(kp),
     &                     WORK(LDXRM),1,WORK(LTMR),1)
               CALL DAXPY_(3*NSS**2,kPhase(kp),
     &                     WORK(LDXIM),1,WORK(LTMI),1)
               DO iCar=1, 3
                  LR_=LTMR+(iCar-1)*NSS**2
                  LI_=LTMI+(iCar-1)*NSS**2
                  CALL ZTRNSF_MASKED(NSS,VSOR,VSOI,
     &                               WORK(LR_),WORK(LI_),
     &                               IJSS,iSSMask,ISSM,jSSMask,JSSM)
               END DO
*
               IJ_=0
               Do ISO=istart_,iend_
                 Do JSO=jstart_,jend_
                   IJ=(JSO-1)*NSS+ISO-1
                   IJ_=IJ_+1
                   LFIJ=LF+(ij_-1)*2
                   EDIFF=ENSOR(JSO)-ENSOR(ISO)
*
*              Store the vectors for this direction
*
                   Call DCopy_(3,Work(LTMR+IJ),NSS**2,TM_R,1)
                   Call DCopy_(3,Work(LTMI+IJ),NSS**2,TM_I,1)
#ifdef _HDF5_
*              Get proper triangular index
                   IJSO_=(JSO-1)*(JSO-2)/2+ISO
                   iQuad_=2*(iQuad-1)+kp
                   Storage(ip_w,iQuad_,IJSO_,iVec)=Weight
                   Call DCopy_(3,kPhase(kp)*Wavevector,1,
     &                  Storage(ip_kvector,iQuad_,IJSO_,iVec),1)
                   Call DCopy_(3,TM_R,1,
     &                  Storage(ip_TMR,iQuad_,IJSO_,iVec),1)
                   Call DCopy_(3,TM_I,1,
     &                  Storage(ip_TMI,iQuad_,IJSO_,iVec),1)
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
*                     The direction for the maximum
                      Ang = Half*Atan2(2.0D0*TM3,TM1-TM2)
                      Call daXpY_(3, Cos(Ang),TM_R,1,Work(LMAX_+1),1)
                      Call daXpY_(3, Sin(Ang),TM_I,1,Work(LMAX_+1),1)
*                     Normalize and compute the direction for the minimum
*                     as a cross product with k
                      rNorm = DDot_(3,Work(LMAX_+1),1,Work(LMAX_+1),1)
                      If (rNorm.gt.1.0D-12) Then
                         Call dScal_(3,1.0/Sqrt(rNorm),Work(LMAX_+1),1)
                         Work(LMAX_+5)=Work(LMAX_+2)*UK(3)-
     &                                 Work(LMAX_+3)*UK(2)
                         Work(LMAX_+6)=Work(LMAX_+3)*UK(1)-
     &                                 Work(LMAX_+1)*UK(3)
                         Work(LMAX_+7)=Work(LMAX_+1)*UK(2)-
     &                                 Work(LMAX_+2)*UK(1)
                         rNorm=DDot_(3,Work(LMAX_+5),1,Work(LMAX_+5),1)
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
*              Compute the rotatory strength, note that it depends on kPhase
*
                   TM_C(1) = TM_R(2)*TM_I(3)-TM_R(3)*TM_I(2)
                   TM_C(2) = TM_R(3)*TM_I(1)-TM_R(1)*TM_I(3)
                   TM_C(3) = TM_R(1)*TM_I(2)-TM_R(2)*TM_I(1)
                   TM_2 = 2.0D0*kPhase(kp)*DDot_(3,TM_C,1,UK,1)
                   R_Temp=0.75D0*SPEED_OF_LIGHT/EDIFF**2*TM_2
                   R_Temp=R_Temp*AU2REDR
*
*              Save the raw oscillator and rotatory strengths in a given direction
*
                   NQUAD_=2*NQUAD
                   LRAW_=LRAW+6*NQUAD_*(ij_-1)
                   IQUAD_=2*(IQUAD-1)+(kp-1)
                   WORK(LRAW_+IQUAD_+0*NQUAD_) = F_Temp
                   WORK(LRAW_+IQUAD_+1*NQUAD_) = R_Temp
*
*              Save the direction and weight too
*
                   WORK(LRAW_+IQUAD_+2*NQUAD_) = UK(1)*kPhase(kp)
                   WORK(LRAW_+IQUAD_+3*NQUAD_) = UK(2)*kPhase(kp)
                   WORK(LRAW_+IQUAD_+4*NQUAD_) = UK(3)*kPhase(kp)
                   WORK(LRAW_+IQUAD_+5*NQUAD_) = Weight
*
*              Do not accumulate if not doing an isotropic integration
*
                   If (Do_SK.And.(kp.gt.1)) Cycle
*
*              Accumulate to the isotropic oscillator strength
*
                   Work(LFIJ  )=Work(LFIJ  ) + Weight * F_Temp
*
*              Accumulate to the isotropic rotatory strength
*
                   Work(LFIJ+1)=Work(LFIJ+1) + Weight * R_Temp
                 End Do
               End Do
*
               End Do ! kp
*
            End Do ! iQuad
*
            IJ_=0
            Do ISO=istart_,iend_
              Do JSO=jstart_,jend_
                 IJ=(ISO-1)*NSS+JSO-1
                 IJ_=IJ_+1
                 LFIJ=LF+(ij_-1)*2
                 EDIFF=ENSOR(JSO)-ENSOR(ISO)
                 F=Work(LFIJ)
                 R=Work(LFIJ+1)
*
                 Call Add_Info('ITMS(SO)',[F],1,6)
                 Call Add_Info('ROTS(SO)',[R],1,4)
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
     &              'Transition moment strengths (SO states):')
                  WRITE(6,'(3X,A)')
     &              '----------------------------------------'
                  If (Do_Pol) Then
                      iVec_=(iVec-1)*nQuad+1
                      WRITE(6,'(4x,a,3F8.4)')
     &                  'Direction of the polarization: ',
     &                  (pol_vector(k,iVec),k=1,3)
                  Else
                     WRITE(6,'(4x,a)')
     &                 'The oscillator strength is integrated '//
     &                 'over all directions of the polarization '//
     &                 'vector'
                  End If
                  WRITE(6,'(4x,a,3F8.4)')
     &                  'Direction of the k-vector: ',
     &                   (Work(ipR+k),k=0,2)
               Else
                  CALL CollapseOutput(1,
     &              'Isotropic transition moment strengths '//
     &              '(SO states):')
                  WRITE(6,'(3X,A)')
     &              '--------------------------------------'//
     &              '------------'
               End If
               IF (OSTHR.GT.0.0D0) THEN
                  WRITE(6,30) 'for osc. strength at least ',OSTHR
               END IF
               WRITE(6,*)
               If (.NOT.Do_SK) Then
                 WRITE(6,'(4x,a,I4,a)')
     &             'Integrated over ',2*nQuad,' directions of the '//
     &             'wave vector'
                 WRITE(6,'(4x,a)')
     &             'The oscillator strength is '//
     &             'integrated over all directions of the polar'//
     &             'ization vector'
                  WRITE(6,*)
               End If
               WRITE(6,31) 'From', 'To', 'Osc. strength',
     &                     'Red. rot. str.', 'Total A (sec-1)'
               WRITE(6,32)
              iPrint=1
            END IF
*
*     Regular print
*
            WRITE(6,33) ISO,JSO,F,R,A
*
            IF (Do_SK) THEN
               WRITE(6,50) 'maximum',WORK(LMAX_+0),
     &            'for polarization direction:',
     &            WORK(LMAX_+1),WORK(LMAX_+2),WORK(LMAX_+3)
               WRITE(6,50) 'minimum',WORK(LMAX_+4),
     &            'for polarization direction:',
     &            WORK(LMAX_+5),WORK(LMAX_+6),WORK(LMAX_+7)
            END IF
*
*
*     Printing raw (unweighted) and direction for every transition
*
            IF(PRRAW) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,34) 'From', 'To', 'Raw osc. str.',
     &                    'Red. rot. str.','kx','ky','kz'
              WRITE(6,35)
              NQUAD_=2*NQUAD
              LRAW_=LRAW+6*NQUAD_*(ij_-1)
              DO IQUAD = 1, NQUAD
                DO kp=1,2
                  IF (ABS(kPhase(kp)).LT.0.5D0) CYCLE
                  IQUAD_=2*(IQUAD-1)+(kp-1)
                  WRITE(6,33) ISO,JSO,
     &            WORK(LRAW_+IQUAD_+0*NQUAD_),
     &            WORK(LRAW_+IQUAD_+1*NQUAD_),
     &            WORK(LRAW_+IQUAD_+2*NQUAD_),
     &            WORK(LRAW_+IQUAD_+3*NQUAD_),
     &            WORK(LRAW_+IQUAD_+4*NQUAD_)
                END DO
              END DO
              WRITE(6,35)
              WRITE(6,*)
            END IF
*
*     Printing weighted and direction for every transition
*
            IF(PRWEIGHT) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,34) 'From', 'To', 'Weig. osc. str.',
     &                    'Red. rot. str.','kx','ky','kz'
              WRITE(6,35)
              NQUAD_=2*NQUAD
              LRAW_=LRAW+6*NQUAD_*(ij_-1)
              DO IQUAD = 1, NQUAD
                DO kp=1,2
                  IF (ABS(kPhase(kp)).LT.0.5D0) CYCLE
                  IQUAD_=2*(IQUAD-1)+(kp-1)
                  Weight=WORK(LRAW_+IQUAD_+5*NQUAD_)
                  WRITE(6,33) ISO,JSO,
     &            WORK(LRAW_+IQUAD_+0*NQUAD_)*Weight,
     &            WORK(LRAW_+IQUAD_+1*NQUAD_)*Weight,
     &            WORK(LRAW_+IQUAD_+2*NQUAD_),
     &            WORK(LRAW_+IQUAD_+3*NQUAD_),
     &            WORK(LRAW_+IQUAD_+4*NQUAD_)
                END DO
              END DO
              WRITE(6,35)
              WRITE(6,*)
            END IF
*
               End Do
            End Do
*
         END DO
      END DO
      Call mma_deallocate(iMask)
      Call mma_deallocate(jMask)
      Call mma_deallocate(iSSMask)
      Call mma_deallocate(jSSMask)
      Call mma_deallocate(ISS_INDEX)
*
      If (iPrint.EQ.1) THEN
         WRITE(6,32)
         If (Do_SK) Then
            CALL CollapseOutput(0,
     &            'Transition moment strengths (SO states):')
         Else
            CALL CollapseOutput(0,
     &            'Isotropic transition moment strengths (SO states):')
         End If
      END IF
*
      End Do ! iVec
#ifdef _TIME_TMOM_
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for TMOM : ',TCpu2-TCpu1,TWall2-TWall1
#endif
*
#ifdef _HDF5_
#ifdef POINTER_REMAP
      flatStorage(1:SIZE(Storage)) => Storage
#else
      Call C_F_Pointer(C_Loc(Storage), flatStorage, [SIZE(Storage)])
#endif
      Call mh5_put_dset(wfn_sos_tm,flatStorage)
      Nullify(flatStorage)
      Call mma_deallocate(Storage)
#endif
*
*     Do some cleanup
*
      CALL GETMEM('RAW   ','FREE','REAL',LRAW,2*NQUAD*6*nmax2)
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
      CALL GETMEM('DXI','FREE','REAL',LDXI,3*NSS**2)
      CALL GETMEM('DXRM','FREE','REAL',LDXRM,NSS**2)
      CALL GETMEM('DXIM','FREE','REAL',LDXIM,3*NSS**2)
      call mma_deallocate(TMP)
      if (TMOgroup) Then
        Call mma_DeAllocate(TMOgrp1)
        Call mma_DeAllocate(TMOgrp2)
      EndIf
      If (Do_Pol) Call mma_deallocate(pol_Vector)
      CALL GETMEM('TMR','FREE','REAL',LTMR,3*NSS**2)
      CALL GETMEM('TMI','FREE','REAL',LTMI,3*NSS**2)
      CALL GETMEM('OSCSTR','FREE','REAL',LF,2*nmax2)
      CALL GETMEM('MAXMIN','FREE','REAL',LMAX,8*nmax2)
*
      Call mma_deallocate(SCR)
      Call mma_deAllocate(TDMZZ)
      Call mma_deAllocate(TSDMZZ)
      Call mma_deAllocate(WDMZZ)
      If (.NOT.Do_SK) Call Free_O()
      Call Free_Work(ipR)
      Call ClsSew()

 666  Continue
      Call mma_deallocate(VSOR)
      Call mma_deallocate(VSOI)
*
30    FORMAT (5X,A,1X,ES15.8)
31    FORMAT (5X,2(1X,A4),4X,3(1X,A15))
32    FORMAT (5X,63('-'))
33    FORMAT (5X,2(1X,I4),5X,5(1X,ES15.8))
34    FORMAT (5X,2(1X,A4),5X,5(1X,A15))
35    FORMAT (5X,95('-'))
50    FORMAT (10X,A7,3X,1(1X,ES15.8),5X,A27,3(1X,F7.4))
*
************************************************************************
*                                                                      *
*     End of section for transition moments                            *
*                                                                      *
************************************************************************
*
#ifdef _TIME_TMOM_
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for TMOM(SO) : ',TCpu2-TCpu1,TWall2-TWall1
#endif
      RETURN
      END Subroutine PRPROP_TM_Exact

