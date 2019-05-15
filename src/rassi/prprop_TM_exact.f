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
      SUBROUTINE PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,JBNUM)
      USE kVectors
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
      DIMENSION PROP(NSTATE,NSTATE,NPROP),JBNUM(NSTATE)
#include "SysDef.fh"
#include "rassiwfn.fh"
      LOGICAL TMOgroup
      INTEGER IOFF(8),IJSF(4)
      CHARACTER*8 LABEL
      Integer, Dimension(:), Allocatable :: TMOgrp1,TMOgrp2
      Real*8 TM_R(3), TM_I(3), TM_C(3)
      Real*8 wavevector(3), UK(3)

      CALL QENTER(ROUTINE)

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

      Call CWTime(TCpu1,TWall1)
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
      IPRTMOM_RS=-1
      DO IPROP=1,NPROP
         IF (PNAME(IPROP).EQ.'TMOM  RS'.AND.IPRTMOM_RS.EQ.-1) THEN
            IPRTMOM_RS=IPROP
         END IF
      ENDDO
      IF (IPRTMOM_RS.EQ.-1) RETURN
      IPRTMOM_0R=IPRTMOM_RS-2
      IF (PNAME(IPRTMOM_0R).NE.'TMOM0  R') RETURN
      IPRTMOM_0I=IPRTMOM_RS-1
      IF (PNAME(IPRTMOM_0I).NE.'TMOM0  I') RETURN
      IPRTMOM_RA=IPRTMOM_RS+3
      IF (PNAME(IPRTMOM_RA).NE.'TMOM  RA') RETURN
      IPRTMOM_IS=IPRTMOM_RS+6
      IF (PNAME(IPRTMOM_IS).NE.'TMOM  IS') RETURN
      IPRTMOM_IA=IPRTMOM_RS+9
      IF (PNAME(IPRTMOM_IA).NE.'TMOM  IA') RETURN
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
         nVec=nk_Vector
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
      CALL GETMEM('TMR','ALLO','REAL',LTMR,3*NSS**2)
      CALL GETMEM('TMI','ALLO','REAL',LTMI,3*NSS**2)
*
C     ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
      NIP=4+(NBST*(NBST+1))/2
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
      NSCR=(NBST*(NBST+1))/2
      CALL GETMEM('TDMSCR','Allo','Real',LSCR,4*NSCR)
*
*     Allocate vector to store all individual transition moments.
*     We do this for
*     all unique pairs ISO-JSO, iSO=/=JSO (NSS*(NSS-1)/2)
*         all k-vectors (nQuad or nVec)
*             we store:
*                 the weight (1)
*                 the k-vector (3)
*                 we projected transition vector (real and imaginary parts) (2*3)
*
      nIJ=nSS*(nSS-1)/2
      nData= 1 + 3 + 2*3
      nStorage = nIJ * nQuad * nData
      mStorage = nStorage * nVec
*MGD Storage is not well handled with groups yet
      Call GetMem('STORAGE','Allo','Real',ipStorage,mStorage)
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
        write(6,*) (TMOgrp1(i),i=1,ngroup1+1)
        write(6,*) (TMOgrp2(i),i=1,ngroup2+1)
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
      CALL GETMEM('RAW   ','ALLO','REAL',LRAW,NQUAD*5*nmax2)
      CALL GETMEM('WEIGHT','ALLO','REAL',LWEIGH,NQUAD*5*nmax2)
      CALL GETMEM('OSCSTR','ALLO','REAL',LF,5*nmax2)
*
      ip_w       = 1
      ip_kvector = ip_w + 1
      ip_TMR     = ip_kvector + 3
      ip_TMI     = ip_TMR + 3
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
      ThrSparse=1.0D-10

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
         ISFstart=Nstate
         ISSstart=NSS
         Do ISO=istart_,iend_
           ISS=0
           Do ISF=1,ISFstart-1
              JOB1=iWork(lJBNUM+ISF-1)
              MPLET1=MLTPLT(JOB1)
              Do IMS=1,MPLET1
                ISS=ISS+1
                Temp=USOR(ISS,ISO)**2 + USOI(ISS,ISO)**2
                If (Temp.gt.ThrSparse) Then
                   ISFstart=ISF
                   ISSstart=ISS
                   Go to 21
                EndIf
              EndDo
           End Do
 21        Continue
         End Do
         ISFend=ISFstart
         ISSend=ISSstart
         Do ISO=istart_,iend_
           ISS=NSS+1
           Do ISF=Nstate,ISFend+1,-1
             JOB1=iWork(lJBNUM+ISF-1)
             MPLET1=MLTPLT(JOB1)
             Do IMS=MPLET1,1,-1
               ISS=ISS-1
               Temp=USOR(ISS,ISO)**2 + USOI(ISS,ISO)**2
               If (Temp.gt.ThrSparse) Then
                  ISFend=ISF
                  ISSend=ISS
                  Go to 22
               EndIf
             EndDo
           End Do
 22        Continue
         End Do
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
            JSFstart=NState
            JSSstart=NSS
            Do JSO=jstart_,jend_
              JSS=0
              Do JSF=1,JSFstart-1
                 JOB1=iWork(lJBNUM+JSF-1)
                 MPLET1=MLTPLT(JOB1)
                 Do JMS=1,MPLET1
                   JSS=JSS+1
                   Temp=USOR(JSS,JSO)**2 + USOI(JSS,JSO)**2
                   If (Temp.gt.ThrSparse) Then
                      JSFstart=JSF
                      JSSstart=JSS
                      Go to 23
                   EndIf
                 End Do
              End Do
 23           Continue
            End Do
            JSFend=JSFstart
            JSSend=JSSstart
            Do JSO=jstart_,jend_
              JSS=NSS+1
              Do JSF=NState,JSFend+1,-1
                JOB1=iWork(lJBNUM+JSF-1)
                MPLET1=MLTPLT(JOB1)
                Do JMS=MPLET1,1,-1
                  JSS=JSS-1
                  Temp=USOR(JSS,JSO)**2 + USOI(JSS,JSO)**2
                  If (Temp.gt.ThrSparse) Then
                      JSFend=JSF
                      JSSend=JSS
                      Go to 24
                   EndIf
                End Do
              End Do
 24           Continue
            End Do
            IJSF(1)=JSFstart
            IJSF(2)=JSFend
            IJSF(3)=ISFstart
            IJSF(4)=ISFend
*
            IF (ABS(EDIFF_).LE.1.0D-8) CYCLE
            IF(EDIFF_.LT.0.0D0) CYCLE
            IJSO=IJSO+1
            iOff_= (IJSO-1)*nQuad*nData
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
            CALL DCOPY_(5*n12,[0.0D0],0,WORK(LF),1)
            CALL DCOPY_(NQUAD*5*n12,[0.0D0],0,WORK(LRAW),1)
            CALL DCOPY_(NQUAD*5*n12,[0.0D0],0,WORK(LWEIGH),1)
*
            Do iQuad = 1, nQuad
               iStorage = iOff_+ (iQuad-1)*nData + ipStorage - 1
     &                  + (iVec-1)*nStorage
*
*              Generate the wavevector associated with this quadrature
*              point and pick up the associated quadrature weight.
*
               UK(1)=Work((iQuad-1)*4  +ipR)
               UK(2)=Work((iQuad-1)*4+1+ipR)
               UK(3)=Work((iQuad-1)*4+2+ipR)

               wavevector(:)=rkNorm*UK(:)
               Call DCopy_(3,wavevector,1,Work(iStorage+ip_kvector),1)
*
               Weight=Work((iQuad-1)*4+3+ipR)
               Work(iStorage+ip_w)=Weight
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
               DO IPROP = IPRTMOM_RS-2, IPRTMOM_RS+11
                  Call FZero(PROP(1,1,IPROP),NSTATE**2)
               End Do

               DO ISS=ISFstart, ISFend
*
*                 Does this spin-free state contribute to any of the
*                 two spin states? Check the corresponding coefficients.
*
                  DO JSS=JSFstart, JSFend
                     j=JSS
                     i=ISS
                     JOB1=iWork(lJBNUM+I-1)
                     MPLET1 = MLTPLT(JOB1)
*
                     JOB2=iWork(lJBNUM+J-1)
                     MPLET2 = MLTPLT(JOB2)
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
                     ISTATE=max(ISS,JSS)
                     JSTATE=min(ISS,JSS)
                     ij=ISTATE*(ISTATE-1)/2+JSTATE
                     iDisk=iWork(liTocM+ij-1)
                     Call dDaFile(LuToM,2,Work(LSCR),4*NSCR,iDisk)
*
*                    Compute the transition property of the property
*                    integrals between the two states.
*
                     DO IPROP = IPRTMOM_RS-2, IPRTMOM_RS+11
                        ITYPE=0
                        IF (PTYPE(IPROP).EQ.'HERMSING') ITYPE=1
                        IF (PTYPE(IPROP).EQ.'ANTISING') ITYPE=2
                        IF (PTYPE(IPROP).EQ.'HERMTRIP') ITYPE=3
                        IF (PTYPE(IPROP).EQ.'ANTITRIP') ITYPE=4
                        LABEL=PNAME(IPROP)
                        Call MK_PROP(PROP,IPROP,ISTATE,JSTATE,
     &                               LABEL,ITYPE,
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
*              momentum operator but rather for nabla. That is
*              as we assemble to transition momentum we have to
*              remember to put in a factor of -i.
*
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LTMR),1)
               CALL DCOPY_(3*NSS**2,[0.0D0],0,WORK(LTMI),1)
               Do iCar = 1, 3
*
*                 The electric (symmetric) part
*
                  CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LDXI),1)
*
*                 the real symmetric part
                  CALL SMMAT2(PROP,WORK(LDXR),NSS,'TMOM  RS',iCar,IJSF)
*
*                 the imaginary symmetric part
                  CALL SMMAT2(PROP,WORK(LDXI),NSS,'TMOM  IS',iCar,IJSF)
*
*                 Transform properties to the spin-orbit basis
*                 and pick up correct element
                  Call ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
*
                  Call DAXPY_(NSS**2, 1.0D0,WORK(LDXI),1,
     &                                      Work(LTMR+iCar-1),3)
                  Call DAXPY_(NSS**2,-1.0D0,WORK(LDXR),1,
     &                                      Work(LTMI+iCar-1),3)
*
*                 The magnetic (antisymmetric) part
*
                  CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LDXR),1)
                  CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LDXI),1)
*
*                 the real anti-symmetric part
                  CALL SMMAT2(PROP,WORK(LDXR),NSS,'TMOM  RA',iCar,IJSF)
*
*                 the imaginary anti-symmetric part
                  CALL SMMAT2(PROP,WORK(LDXI),NSS,'TMOM  IA',iCar,IJSF)
*
*                 Transform properties to the spin-orbit basis
*                 and pick up correct element
                  CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
*
                  Call DAXPY_(NSS**2, 1.0D0,WORK(LDXI),1,
     &                                      Work(LTMR+iCar-1),3)
                  Call DAXPY_(NSS**2,-1.0D0,WORK(LDXR),1,
     &                                      Work(LTMI+iCar-1),3)
*
*
*              (2) the spin-dependent part, magnetic
*
*              iCar=1: 1/2(S(+)+S(-))
*              iCar=2: 1/2i(S(+)-S(-))
*              iCar=3: Sz
*
*              Here we need to be very careful. The operator is
*              similar to the L.S operator, the spin-orbit coupling.
*              However, here we have that the operator is B.S. Thus
*              we can do this in a similar fashion as the spin-orbit
*              coupling, but with some important difference.
*
*              For the spin-orbit coupling the L operator is divided
*              into x-, y-, and z-components, that is the operator
*              will be represented by three different integrals.
*              For the integrals over B it is similar but the x-, y-,
*              and z-components are constants outside of the integral.
*              For example, the x-component is expressed as
*              (k x e_l)_x <0|e^(i k.r)|n>. In this section we will
*              handle the (k x e_l)_x part outside the loop over the
*              Cartesian components.
*
*              We have to note one further difference, the integrals
*              are complex in this case.
*
*              Let us now compute the contributions T(i), i=x,y,z
*
*              In the equations we find that the y-component is imaginary,
*              see the equations on page 234 in Per Ake's paper. Hence,
*              the y-component is treated slightly differently.
*
                  If (iCar.eq.1.or.iCar.eq.3) Then
*
*                    pick up the real component
                   CALL SMMAT2(PROP,WORK(LDXR),NSS,'TMOM0  R',iCar,IJSF)
*                    pick up the imaginary component
                   CALL SMMAT2(PROP,WORK(LDXI),NSS,'TMOM0  I',iCar,IJSF)
                  Else
*                    For the y-component we have to interchange the real and
*                    the imaginary components. The real component gets a
*                    minus sign due to the product ixi=-1
*
*                    pick up the real component
                   CALL SMMAT2(PROP,WORK(LDXI),NSS,'TMOM0  R',iCar,IJSF)
*                    pick up the imaginary component
                   CALL SMMAT2(PROP,WORK(LDXR),NSS,'TMOM0  I',iCar,IJSF)
                   Call DScal_(NSS**2,-1.0D0,WORK(LDXR),1)
                  End If
                  CALL ZTRNSF(NSS,USOR,USOI,WORK(LDXR),WORK(LDXI))
*
*                 i*g/2*(s_y*k_z-s_z*k_y) -> T_x
*                 i*g/2*(s_z*k_x-s_x*k_z) -> T_y
*                 i*g/2*(s_x*k_y-s_y*k_x) -> T_z
*
                  If (iCar.eq.1) Then
                     Call DAXPY_(NSS**2,-wavevector(2)*cst,
     &                           WORK(LDXI),1,Work(LTMR+2),3)
                     Call DAXPY_(NSS**2, wavevector(2)*cst,
     &                           WORK(LDXR),1,Work(LTMI+2),3)
                     Call DAXPY_(NSS**2, wavevector(3)*cst,
     &                           WORK(LDXI),1,Work(LTMR+1),3)
                     Call DAXPY_(NSS**2,-wavevector(3)*cst,
     &                           WORK(LDXR),1,Work(LTMI+1),3)
                  Else If (iCar.eq.2) Then
                     Call DAXPY_(NSS**2,-wavevector(3)*cst,
     &                           WORK(LDXI),1,Work(LTMR+0),3)
                     Call DAXPY_(NSS**2, wavevector(3)*cst,
     &                           WORK(LDXR),1,Work(LTMI+0),3)
                     Call DAXPY_(NSS**2, wavevector(1)*cst,
     &                           WORK(LDXI),1,Work(LTMR+2),3)
                     Call DAXPY_(NSS**2,-wavevector(1)*cst,
     &                           WORK(LDXR),1,Work(LTMI+2),3)
                  Else If (iCar.eq.3) Then
                     Call DAXPY_(NSS**2,-wavevector(1)*cst,
     &                           WORK(LDXI),1,Work(LTMR+1),3)
                     Call DAXPY_(NSS**2, wavevector(1)*cst,
     &                           WORK(LDXR),1,Work(LTMI+1),3)
                     Call DAXPY_(NSS**2, wavevector(2)*cst,
     &                           WORK(LDXI),1,Work(LTMR+0),3)
                     Call DAXPY_(NSS**2,-wavevector(2)*cst,
     &                           WORK(LDXR),1,Work(LTMI+0),3)
                  End If
               End Do
*
               IJ_=0
               Do ISO=istart_,iend_
                 Do JSO=jstart_,jend_
                   IJ=(JSO-1)*NSS+ISO-1
                   IJ_=IJ_+1
                   LFIJ=LF+(ij_-1)*5
                   EDIFF=ENSOR(JSO)-ENSOR(ISO)
*
*                  Project out the k direction from the real and imaginary components
*
                   Call DCopy_(3,Work(LTMR+IJ*3),1,TM_R,1)
                   Call DCopy_(3,Work(LTMI+IJ*3),1,TM_I,1)
                   Call DaXpY_(3,-DDot_(3,TM_R,1,UK,1),UK,1,TM_R,1)
                   Call DaXpY_(3,-DDot_(3,TM_I,1,UK,1),UK,1,TM_I,1)

*MGD not stored yet
                   Call DCopy_(3,TM_R,1,Work(iStorage+ip_TMR),1)
                   Call DCopy_(3,TM_I,1,Work(iStorage+ip_TMI),1)
*
*              Integrate over all directions of the polarization
*              vector and divide with the "distance", 2*pi, to get
*              the average value.
*
                   TM1 = DDot_(3,TM_R,1,TM_R,1)
                   TM2 = DDot_(3,TM_I,1,TM_I,1)
                   TM_2 = Half*(TM1+TM2)
*
*              Compute the oscillator strength
*
                   F_Temp = 2.0D0*TM_2/EDIFF
*
*              Compute the rotatory strength
*
                   TM_C(1) = TM_R(2)*TM_I(3)-TM_R(3)*TM_I(2)
                   TM_C(2) = TM_R(3)*TM_I(1)-TM_R(1)*TM_I(3)
                   TM_C(3) = TM_R(1)*TM_I(2)-TM_R(2)*TM_I(1)
                   TM_2 = -2.0D0*SQRT(DDot_(3,TM_C,1,TM_C,1))*
     &                           SIGN(1.0D0,DDot_(3,TM_C,1,UK,1))
                   R_Temp=0.75D0*SPEED_OF_LIGHT/EDIFF**2*TM_2
                   R_Temp=R_Temp*AU2REDR
*
*              Save the raw oscillator strengths in a given direction
*
                   LRAW_=LRAW+5*NQUAD*(ij_-1)
                   WORK(LRAW_+(IQUAD-1)+0*NQUAD) = F_Temp
                   WORK(LRAW_+(IQUAD-1)+1*NQUAD) = R_Temp
                   WORK(LRAW_+(IQUAD-1)+2*NQUAD) = UK(1)
                   WORK(LRAW_+(IQUAD-1)+3*NQUAD) = UK(2)
                   WORK(LRAW_+(IQUAD-1)+4*NQUAD) = UK(3)
*
*              Accumulate to the isotropic oscillator strength
*
                   Work(LFIJ+3)=Work(LFIJ+3) + Weight * F_Temp
*
*              Accumulate to the isotropic rotatory strength
*
                   Work(LFIJ+4)=Work(LFIJ+4) + Weight * R_Temp
*
*              Save the weighted oscillator and rotatory strengths in a
*              given direction k.
*
                   LWEIGH_=LWEIGH+5*NQUAD*(ij_-1)
                   WORK(LWEIGH_+(IQUAD-1)+0*NQUAD) = F_Temp*WEIGHT
                   WORK(LWEIGH_+(IQUAD-1)+1*NQUAD) = R_Temp*WEIGHT
                   WORK(LWEIGH_+(IQUAD-1)+2*NQUAD) = UK(1)
                   WORK(LWEIGH_+(IQUAD-1)+3*NQUAD) = UK(2)
                   WORK(LWEIGH_+(IQUAD-1)+4*NQUAD) = UK(3)
                 End Do
               End Do

            End Do ! iQuad
*
*           Note that the weights are normalized to integrate to
*           4*pi over the solid angles.
*
            IJ_=0
            Do ISO=istart_,iend_
              Do JSO=jstart_,jend_
                 IJ=(ISO-1)*NSS+JSO-1
                 IJ_=IJ_+1
                 LFIJ=LF+(ij_-1)*5
                 EDIFF=ENSOR(JSO)-ENSOR(ISO)
                 F=Work(LFIJ+3)
                 R=Work(LFIJ+4)
*
            If (.NOT.Do_SK) Then
               F = F / (4.0D0*PI)
               R = R / (4.0D0*PI)
            End If
            IF (ABS(F).LT.OSTHR) CYCLE
            A =(AFACTOR*EDIFF**2)*F
*
            If (iPrint.eq.0) Then
               WRITE(6,*)
               If (Do_SK) Then
                  CALL CollapseOutput(1,
     &            'Transition moment strengths (SO states):')
                  WRITE(6,'(3X,A)')
     &            '----------------------------------------'
                  WRITE(6,'(4x,a)')
     &            'The oscillator strength is '//
     &            'integrated over all directions of the polar'//
     &            'ization vector'
                  WRITE(6,'(4x,a,3F8.4,a)')
     &                  'Direction of the k-vector: ',
     &                   (Work(ipR+k),k=0,2),' (au)'
               Else
                  CALL CollapseOutput(1,
     &            'Isotropic transition moment strengths (SO states):')
                  WRITE(6,'(3X,A)')
     &            '--------------------------------------------------'
               End If
               IF (OSTHR.GT.0.0D0) THEN
                  WRITE(6,'(4x,a,ES16.8)')
     &                  'for osc. strength at least ',OSTHR
               END IF
               WRITE(6,*)
               If (.NOT.Do_SK) Then
                 WRITE(6,'(4x,a,I4,a)')
     &             'Integrated over ',nQuad,' directions of the '//
     &             'wave vector'
                 WRITE(6,'(4x,a)')
     &             'The oscillator and strength is '//
     &             'integrated over all directions of the polar'//
     &             'ization vector'
                 WRITE(6,*)
               End If

               WRITE(6,31) 'From', 'To', 'Osc. strength',
     &                     'Red. rot. str.', 'Total A (sec-1)'
               WRITE(6,32)
              iPrint=1
            END IF
            WRITE(6,33) ISO,JSO,F,R,A
*
*     Printing raw (unweighted) and direction for every transition
*
            IF(PRRAW) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,34) 'From', 'To', 'Raw osc. str.',
     &                    'Mag. cont.','kx','ky','kz'
              WRITE(6,35)
              LRAW_=LRAW+5*NQUAD*(ij_-1)
              DO IQUAD = 1, NQUAD
                WRITE(6,33) ISO,JSO,
     &          WORK(LRAW_+(IQUAD-1)+0*NQUAD),
     &          WORK(LRAW_+(IQUAD-1)+1*NQUAD),
     &          WORK(LRAW_+(IQUAD-1)+2*NQUAD),
     &          WORK(LRAW_+(IQUAD-1)+3*NQUAD),
     &          WORK(LRAW_+(IQUAD-1)+4*NQUAD)
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
     &                    'Mag. cont.','kx','ky','kz'
              WRITE(6,35)
              LWEIGH_=LWEIGH+5*NQUAD*(ij_-1)
              DO IQUAD = 1, NQUAD
                WRITE(6,33) ISO,JSO,
     &          WORK(LWEIGH_+(IQUAD-1)+0*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH_+(IQUAD-1)+1*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH_+(IQUAD-1)+2*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH_+(IQUAD-1)+3*NQUAD)/ (4.0D0*PI),
     &          WORK(LWEIGH_+(IQUAD-1)+4*NQUAD)/ (4.0D0*PI)
              END DO
              WRITE(6,35)
              WRITE(6,*)
            END IF
*
            Call Add_Info('ITMS(SO)',[F],1,6)
            Call Add_Info('ROTS(SO)',[R],1,6)
*
               End Do
            End Do
*
         END DO
      END DO
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
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for TMOM : ',TCpu2-TCpu1,TWall2-TWall1
*
#ifdef _HDF5_
      Call mh5_put_dset(wfn_sos_tm,Work(ipStorage))
#endif
*
*     Do some cleanup
*
      Call GetMem('STORAGE','FREE','Real',ipStorage,nStorage)
      CALL GETMEM('RAW   ','FREE','REAL',LRAW,NQUAD*5*nmax2)
      CALL GETMEM('WEIGHT','FREE','REAL',LWEIGH,NQUAD*5*nmax2)
      CALL GETMEM('TDMSCR','FREE','Real',LSCR,4*NSCR)
      CALL GETMEM('IP    ','FREE','REAL',LIP,NIP)
      CALL GETMEM('DXR','FREE','REAL',LDXR,NSS**2)
      CALL GETMEM('DXI','FREE','REAL',LDXI,NSS**2)
      if (TMOgroup) Then
        Call mma_DeAllocate(TMOgrp1)
        Call mma_DeAllocate(TMOgrp2)
      EndIf
      CALL GETMEM('TMR','FREE','REAL',LTMR,3*NSS**2)
      CALL GETMEM('TMI','FREE','REAL',LTMI,3*NSS**2)
      CALL GETMEM('OSCSTR','FREE','REAL',LF,5*nmax2)
*
      Call DaClos(LuToM)
      If (.NOT.Do_SK) Call Free_O()
      Call Free_Work(ipR)
      Call ClsSew()
*
31    FORMAT (5X,2(1X,A4),6X,A15,1X,A15,1X,A15)
32    FORMAT (5X,63('-'))
33    FORMAT (5X,2(1X,I4),5X,5(1X,ES15.8))
34    FORMAT (5X,2(1X,A4),6X,A15,1X,A15,1X,A15,1X,A15,1X,A15)
35    FORMAT (5X,95('-'))
*
************************************************************************
*                                                                      *
*     End of section for transition moments                            *
*                                                                      *
************************************************************************
*
      RETURN
      END Subroutine PRPROP_TM_Exact

