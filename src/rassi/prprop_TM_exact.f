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
      SUBROUTINE PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,OVLP,ENERGY,
     &                           JBNUM)
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
      DIMENSION PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &          ENERGY(NSTATE),JBNUM(NSTATE)
#include "SysDef.fh"
#include "rassiwfn.fh"
      REAL*8 Boltz_k,coeff_chi
      LOGICAL Sparse_I,Sparse_J
      REAL*8 J2CM
      Real*8 E1(3), E2(3), kxe1(3), kxe2(3)
      INTEGER IOFF(8)
      CHARACTER*8 LABEL
      Complex*16 T0_e(3), T0_m(3), T1(3), TM1, TM2, PE1_e, PE1_m,
*    &                                    TMR, TML, PE2_e, PE2_m,
     &                                              PE2_e, PE2_m,
     &           E1B, E2B,
     &           IMAGINARY

      CALL QENTER(ROUTINE)

      Dummy=Energy(1)
      Dummy=OVLP(1,1)
*
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

C Compute transition strengths for spin-orbit states:
*
* Initial setup for both dipole, quadrupole etc. and exact operator
*
C printing threshold
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
*     Section (1): Computation of k specific oscillator strength.
*     Section (2): Computation of the isotropic oscillator strength.
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
         nQuad = 1
         nVec=nK_Vector
         Call GetMem('SK','ALLO','REAL',ipR,4*nQuad)
      Else
         Call Setup_O()
         Call Do_Lebedev(L_Eff,nQuad,ipR)
         nVec = 1
      End If
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
      mStorage = nStorage * nVec
      Call GetMem('STORAGE','Allo','Real',ipStorage,mStorage)
*
      Do iVec = 1, nVec
*
         ip_w      = 1
         ip_kvector= ip_w + 1
         ip_e1     = ip_kvector + 3
         ip_e2     = ip_e1 + 3
         ip_TM1R   = ip_e2 + 3
         ip_TM1I   = ip_TM1R + 1
         ip_TM2R   = ip_TM1I + 1
         ip_TM2I   = ip_TM2R + 1
*
         If (Do_SK) Then
            Work(ipR  )=k_Vector(1,iVec)
            Work(ipR+1)=k_Vector(2,iVec)
            Work(ipR+2)=k_Vector(3,iVec)
            Work(ipR+3)=1.0D0   ! Dummy weight
         End If
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
            R =0.0D0
*
*           Initialize output arrays
*
            CALL DCOPY_(NQUAD*5,0.0D0,0,WORK(LRAW),1)
            CALL DCOPY_(NQUAD*5,0.0D0,0,WORK(LWEIGH),1)

            Do iQuad = 1, nQuad
               iStorage = iOff_+ (iQuad-1)*nData + ipStorage - 1
     &                  + (iVec-1)*nStorage
*
*              Generate the wavevector associated with this quadrature
*              point and pick up the associated quadrature weight.
*
               xcoor=Work((iQuad-1)*4  +ipR)
               ycoor=Work((iQuad-1)*4+1+ipR)
               zcoor=Work((iQuad-1)*4+2+ipR)

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
                  E1(1)=1.0D0
                  E1(2)=0.0D0
                  E1(3)=0.0D0
               ELSE
                  E1(1)= PORIG(2,IPREMFR_RS)
                  E1(2)=-PORIG(1,IPREMFR_RS)
                  E1(3)= 0.0D0
               END IF
               Tmp=1.0D0/SQRT(E1(1)**2+E1(2)**2+E1(3)**2)
               E1(1)=E1(1)*Tmp
               E1(2)=E1(2)*Tmp
               E1(3)=E1(3)*Tmp
               E2(1)=PORIG(2,IPREMFR_RS)*E1(3)-E1(2)*PORIG(3,IPREMFR_RS)
               E2(2)=PORIG(3,IPREMFR_RS)*E1(1)-E1(3)*PORIG(1,IPREMFR_RS)
               E2(3)=PORIG(1,IPREMFR_RS)*E1(2)-E1(1)*PORIG(2,IPREMFR_RS)
               Tmp=1.0D0/SQRT(E2(1)**2+E2(2)**2+E2(3)**2)
               E2(1)=E2(1)*Tmp
               E2(2)=E2(2)*Tmp
               E2(3)=E2(3)*Tmp
               Call DCopy_(3,E1,1,Work(iStorage+ip_e1),1)
               Call DCopy_(3,E2,1,Work(iStorage+ip_e2),1)
*
*              Compute the vectors (k x e1) and  (k x e2).
*
               kxe1(1)=PORIG(2,IPORIG)*E1(3)-PORIG(3,IPORIG)*E1(2)
               kxe1(2)=PORIG(3,IPORIG)*E1(1)-PORIG(1,IPORIG)*E1(3)
               kxe1(3)=PORIG(1,IPORIG)*E1(2)-PORIG(2,IPORIG)*E1(1)
               kxe2(1)=PORIG(2,IPORIG)*E2(3)-PORIG(3,IPORIG)*E2(2)
               kxe2(2)=PORIG(3,IPORIG)*E2(1)-PORIG(1,IPORIG)*E2(3)
               kxe2(3)=PORIG(1,IPORIG)*E2(2)-PORIG(2,IPORIG)*E2(1)
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
*              (1) the spin-free part. This part is split into a
*                  symmetric and an antisymmetric part, being electric
*                  and magnetic, respectively. We will treat them
*                  separate for now to facilitate the calculation of
*                  the rotatory strength.
*
*              Note that the integrals are not computed for the
*              momentum opertor but rather for nabla. That is
*              as we assemble to transition momentum we have to
*              remember to put in a factor of -i.
*
               Do iCar = 1, 3
*
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LTMP),1)
*
*                 The electric (symmetric) part
*
*                 the real symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXR),NSS,'TMOS  RS',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*
*                 the imaginary symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXI),NSS,'TMOS  IS',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*
*                 Transform properties to the spin-orbit basis
*                 and pick up correct element
                  CALL ZTRNSF_IJ(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI),
     &                           WORK(LTMP),T0_e(iCar),ISO,JSO)
*
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LDXI),1)
                  CALL DCOPY_(NSS**2,0.0D0,0,WORK(LTMP),1)
*
*                 The magnetic (antisymmetric) part
*
*                 the real anti-symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXR),NSS,'TMOS  RA',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*
*                 the imaginary anti-symmetric part
                  CALL SMMAT_CHECK(PROP,WORK(LDXI),NSS,'TMOS  IA',iCar,
     &                             USOR,USOI,ISO,JSO,ThrSparse)
*
*                 Transform properties to the spin-orbit basis
*                 and pick up correct element
                  CALL ZTRNSF_IJ(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI),
     &                           WORK(LTMP),T0_m(iCar),ISO,JSO)
*
               End Do
*
*              Assemble the electric and magnetic contribution of the
*              transition moment for the respective polarization
*              directions from the spin-independent part of the
*              Hamiltonian.
*
               PE1_e = E1(1)*T0_e(1) + E1(2)*T0_e(2) + E1(3)*T0_e(3)
               PE1_m = E1(1)*T0_m(1) + E1(2)*T0_m(2) + E1(3)*T0_m(3)
*
               PE2_e = E2(1)*T0_e(1) + E2(2)*T0_e(2) + E2(3)*T0_e(3)
               PE2_m = E2(1)*T0_m(1) + E2(2)*T0_m(2) + E2(3)*T0_m(3)
*
*              Multiply with the factor -i to have the transition
*              moment corresponding to the momentum operator and not
*              as now nabla!
*
               PE1_e = - Imaginary * PE1_e
               PE1_m = - Imaginary * PE1_m
               PE2_e = - Imaginary * PE2_e
               PE2_m = - Imaginary * PE2_m
*
*              (2) the spin-dependent part, magnetic
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
*              Assemble the spin-dependent part of the transition
*              moment.
*
               E1B=kxe1(1)*T1(1)+kxe1(2)*T1(2)+kxe1(3)*T1(3)
               E2B=kxe2(1)*T1(1)+kxe2(2)*T1(2)+kxe2(3)*T1(3)
*
*              Assemble to total transition moment for the two
*              polarization directions.
*
               TM1 =(PE1_e + PE1_m)+ IMAGINARY*(g_Elec/2.0D0)*E1B
               TM2 =(PE2_e + PE2_m)+ IMAGINARY*(g_Elec/2.0D0)*E2B
*
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
*              Compute the rotatory strength
*
*              TMR = (TM1 + IMAGINARY*TM2)/Sqrt(2.0D0)
*              TML = (TM1 - IMAGINARY*TM2)/Sqrt(2.0D0)
*
*              TM_2 = DBLE(DCONJG(TMR)*TMR - DCONJG(TML)*TML)
               TM_2 = - 2.0D0*(
     &                           DBLE(TM1)*AIMAG(TM2)
     &                          -DBLE(TM2)*AIMAG(TM1)
     &                           )
               R_Temp=2.0D0*TM_2/ABS(EDIFF)
*
*              Save the raw oscillator strengths in a given direction
*
               WORK(LRAW+(IQUAD-1)+0*NQUAD) = F_TEMP
               WORK(LRAW+(IQUAD-1)+1*NQUAD) = R_TEMP
               WORK(LRAW+(IQUAD-1)+2*NQUAD) = XCOOR
               WORK(LRAW+(IQUAD-1)+3*NQUAD) = YCOOR
               WORK(LRAW+(IQUAD-1)+4*NQUAD) = ZCOOR
*
*              Accumulate to the isotropic oscillator strength
*
               F = F + Weight * F_Temp
*
*              Accumulate to the isotropic rotatory strength
*
               R = R + Weight * R_Temp
*
*              Save the weighted oscillator and rotatory strengths in a
*              given direction k.
*
               WORK(LWEIGH+(IQUAD-1)+0*NQUAD) = F_TEMP*WEIGHT
               WORK(LWEIGH+(IQUAD-1)+1*NQUAD) = R_TEMP*WEIGHT
               WORK(LWEIGH+(IQUAD-1)+2*NQUAD) = XCOOR
               WORK(LWEIGH+(IQUAD-1)+3*NQUAD) = YCOOR
               WORK(LWEIGH+(IQUAD-1)+4*NQUAD) = ZCOOR

            End Do ! iQuad
*
*           Note that the weights are normalized to integrate to
*           4*pi over the solid angles.
*
            F = F / (4.0D0*PI)
            R = R / (4.0D0*PI)
            IF (ABS(F).LT.OSTHR) CYCLE
            AX=(AFACTOR*EDIFF**2)*FX
            AY=(AFACTOR*EDIFF**2)*FY
            AZ=(AFACTOR*EDIFF**2)*FZ
            A =(AFACTOR*EDIFF**2)*F
*
            If (iPrint.eq.0) Then
               WRITE(6,*)
               If (Do_SK) Then
                  CALL CollapseOutput(1,
     &                'Transition moment strengths:')
               WRITE(6,'(1x,a)')
     &           '   The oscillator strength is'//
     &           ' integrated over all directions of the polar'//
     &           'ization vector'
                  WRITE(6,'(4x,a,3F8.4,a)')
     &                  'Direction of the k-vector: ',
     &                   (Work(ipR+k),k=0,2),' (au)'
               Else
                  CALL CollapseOutput(1,
     &                'Isotropic transition moment strengths:')
               End If
               WRITE(6,'(3X,A)')
     &                '--------------------------------------'
               IF (OSTHR.GT.0.0D0) THEN
                  WRITE(6,'(1x,a,ES16.8)')
     &                  '   for osc. strength at least ',OSTHR
               END IF
               WRITE(6,*)
               If (.NOT.Do_SK) Then
               WRITE(6,'(1x,a,I4,a)')
     &           '   Integrated over ',nQuad,' directions of the'//
     &           ' wave vector'
               WRITE(6,'(1x,a)')
     &           '   The oscillator strength is'//
     &           ' integrated over all directions of the polar'//
     &           'ization vector'
               WRITE(6,*)
               End If
               WRITE(6,*)"        To  From     Osc. strength"//
     &           "    Rot. strength",
     &           "   Einstein coefficients Ax, Ay, Az (sec-1) "//
     &           "      Total A (sec-1)  "
               WRITE(6,*)
     &  '      -------------------------------------------'//
     &  '------------------------------------------------'
              iPrint=1
            END IF
            WRITE(6,'(5X,2I5,5X,2(F8.6,8X),4ES16.8)')
     &           ISO,JSO,F,R,AX,AY,AZ,A
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
            Call Add_Info('ROTS(SO)',R,1,6)
*
         END DO
      END DO
*
      If (iPrint.EQ.1) THEN
         If (Do_SK) Then
            CALL CollapseOutput(0,
     &                'Transition moment strengths:')
         Else
            CALL CollapseOutput(0,
     &                'Isotropic transition moment strengths:')
         End If
      END IF
*
      End Do ! iVec
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
      RETURN
      END

