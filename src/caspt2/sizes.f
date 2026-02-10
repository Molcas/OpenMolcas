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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE SIZES()
      use definitions, only: iwp, wp
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: usual
      USE SUPERINDEX, only: SUPINI
      use stdalloc, only: mma_MaxDBLE
      use caspt2_global, only: NCMO
      use caspt2_global, only: do_csf, do_grad, do_nac, if_invar,
     &                         if_invaria, if_SSDM, ipea_shift
      use Cholesky, only: NumCho
      use EQSOLV, only: NLSTOT, NLIST, IFCOUP
      use caspt2_module, only: IfChol, IfDW, IfMSCoup, IfProp, IfXMS,
     &                         IfRMS, nActEl, nAshT, nBasT, nBSqT,
     &                         nBTri, nCases, nConf, nFroT, nIMx,
     &                         nIshT, nOSqT, nOTri, nSMx, nState, nSym,
     &                         Zeta, Mul, nTU, nTGTU, nAsh, nISh, nSsh,
     &                         nOrb, nInDep, nASup, nISup
      use pt2_guga, only: MxCI, nG1, nG2, nG3Tot, nPlBuf

      IMPLICIT none

#include "warnings.h"

      integer(kind=iwp) I, ICASE, ICASE1, ICASE2, ILIST, ISL1, ISL2,
     &                  ISL3, ISYM, ISYM1, ISYM2, M, M11, M12, M21, M22,
     &                  MAXAIS, MC1S1DER, MC2DER, MEMBASE, MINBUF,
     &                  MMX, MXLeft, MXLFT, N, NA, NAS, NAS1, NAS2,
     &                  NBOTTOM, NCH, nCLag, NCX, NDD, NEED, NEED0,
     &                  NFIA, NFIT, NFTA, NG02, NG10, NG12, NG20, NG30,
     &                  NGARR, ngrad, ngrad1, ngrad10, ngrad11,
     &                  ngrad11_1, ngrad11_2, ngrad2, ngrad3, ngrad4,
     &                  ngrad5, ngrad6, ngrad6_1, ngrad6_2, ngrad7,
     &                  ngrad7_1, ngrad7_2, ngrad8, ngrad9, NI, NIN,
     &                  NIN1, NIS, NIS1, NLISTS, NMKRHS, NMX, nOLag,
     &                  nOMax, NPoly, nPrp, nPrp1, nPrp2, nRHSP, NS,
     &                  nSgm1, nSgm2, nSigma, nSigma_inner,
     &                  nSigma_outer, NSLag, nTG1, nTG2, nTG3, NumChT,
     &                  nV, nVCUtil, nWLag, NX, M31, NIS2
#ifdef _DEBUGPRINT_
      integer(kind=iwp) NAMX
#endif
      integer(kind=iwp), external :: iParDiv
      real(kind=wp) XX, YY

C Available workspace right now:
      CALL mma_MaxDBLE(MXLEFT)
C 250000 words margin, for various purposes.
      NBOTTOM=250000
C POLY package:
C MINBUF=Min acceptable nr of buffers in POLY3 for efficiency.
      MINBUF=1
      NPOLY=NBOTTOM
C Sizes of G10,G20,G02,G12,G30: Faked, since no more used.
      NG10=1
      NG20=1
      NG02=1
      NG30=1
      NG12=1

      IF(NACTEL.GT.0) THEN
        NGARR=2*(NG10+NG20+NG02+NG30+NG12)
        NPOLY=NPOLY+NGARR+5*MXCI
        MXLFT=MXLEFT-NPOLY
        XX=0.5D0*DBLE(3*MXCI+3)
        YY=DBLE(MXLFT)
C Allowed by memory:
        NPLBUF=INT(SQRT(YY+XX**2)-XX)
C Preferred size for efficiency:
        NPLBUF=MAX(MINBUF,NPLBUF)
C Actual max needed in-core:
        NPLBUF=MIN(NPLBUF,(NASHT*(NASHT+1))/2)
        NPOLY=NPOLY+NPLBUF*(NPLBUF+3+3*MXCI)
*        write(6,*)' Memory requirements for POLY3 (Above SGUGA):'
*        write(6,*)'   CI vector                      :',MXCI
*        write(6,*)'   Extra margin for small scratch :',NBOTTOM
*        write(6,*)'   Ten arrays G10,G20..F12        :',NGARR
*        write(6,*)'   H0,SGM0,C1,C2                  :',4*MXCI
*        write(6,*)'   VXYZB,TUBUF,DTUB               :',3*NPLBUF
*        write(6,*)'   A,B1,B2                        :',3*NPLBUF*MXCI
*        write(6,*)'   SCR                            :',NPLBUF**2
*        write(6,*)
*        write(6,*)' Total, for POLY3                 :',NPOLY
*        write(6,*)
      END IF

C Precompute sizes, offsets etc.
      CALL SUPINI

C SBMAT need:
*     NG3C=0
*     DO ISYM=1,NSYM
*       N=NTUV(ISYM)
*       NG3C=NG3C+(N*(N+1))/2
*     END DO
*     NG3C=iPARDIV(NG3TOT,NG2)

C Sizes and addresses to lists:
      DO ISL1=1,NSYM
       DO ISL3=1,NSYM
        ISL2=MUL(ISL1,ISL3)
        NLIST(ISL1,ISL3,1)= NASH(ISL2)*NTU(ISL3)*2
        NLIST(ISL1,ISL3,2)= NLIST(ISL1,ISL3,1)
        NLIST(ISL1,ISL3,3)= NASH(ISL2)*NTU(ISL3)
        NLIST(ISL1,ISL3,4)= NASH(ISL2)*NTGTU(ISL3)*2
        NLIST(ISL1,ISL3,5)= NLIST(ISL1,ISL3,3)
        NLIST(ISL1,ISL3,6)= NLIST(ISL1,ISL3,4)
        NLIST(ISL1,ISL3,7)= NASH(ISL2)*NASH(ISL3)*2
        NLIST(ISL1,ISL3,8)= NLIST(ISL1,ISL3,7)
        NLIST(ISL1,ISL3,9)= NASH(ISL2)*NASH(ISL3)
        NLIST(ISL1,ISL3,10)= NLIST(ISL1,ISL3,9)
        IF(ISL1.EQ.1) NLIST(ISL1,ISL3,10)
     &   =NLIST(ISL1,ISL3,9)-NASH(ISL2)
        NLIST(ISL1,ISL3,12)= NASH(ISL1)*NASH(ISL2)
        NLIST(ISL1,ISL3,13)= NLIST(ISL1,ISL3,12)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,13)
     &   =NLIST(ISL1,ISL3,12)-NASH(ISL1)
        NLIST(ISL1,ISL3,11)=NLIST(ISL1,ISL3,12)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,11)
     &   =NLIST(ISL1,ISL3,11)+NASH(ISL1)*NASHT
        NLIST(ISL1,ISL3,14)= NISH(ISL1)*NISH(ISL2)
        NLIST(ISL1,ISL3,15)= NLIST(ISL1,ISL3,14)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,15)
     &   =NLIST(ISL1,ISL3,14)-NISH(ISL1)
        NLIST(ISL1,ISL3,16)= NSSH(ISL1)*NSSH(ISL2)
        NLIST(ISL1,ISL3,17)= NLIST(ISL1,ISL3,16)
        IF(ISL3.EQ.1) NLIST(ISL1,ISL3,17)
     &   =NLIST(ISL1,ISL3,16)-NSSH(ISL1)
       END DO
      END  DO
      NLISTS=0
      DO ILIST=1,17
       DO ISL1=1,NSYM
        DO ISL3=1,NSYM
         NLISTS=NLISTS+NLIST(ISL1,ISL3,ILIST)
        END DO
       END DO
      END DO
      NLSTOT=4*NLISTS

C maximum orbitals in an irrep
      NOMAX=0
      DO ISYM=1,NSYM
        NOMAX=MAX(NOMAX,NORB(ISYM))
      END DO
      IF (.NOT.IFCHOL) THEN
C MKRHS needs:
        NMKRHS=2*NOMAX**2
        NMX=0
        DO ICASE=1,13
          DO ISYM=1,NSYM
            NIN=NINDEP(ISYM,ICASE)
            IF(NIN.EQ.0) GOTO 10
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NV=NAS*NIS
            IF(NV.EQ.0) GOTO 10
            IF(ICASE.GT.11) THEN
              N=NV
            ELSE
              N=2*NV+(NAS*(NAS+1))/2
            END IF
            NMX=MAX(N,NMX)
  10        CONTINUE
          END DO
        END DO
        NMKRHS=NMKRHS+NMX
*       IF(MAXIT.EQ.0) THEN
*         NSIGMA=0
*       ELSE
*         NSIGMA=NBOTTOM+NLSTOT+MMX
*       END IF
      ELSE
C RHSALL2 and ADDRHS needs:
        NMKRHS=0
        DO ICASE=1,13
          DO ISYM=1,NSYM
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NMKRHS=MAX(NMKRHS,NAS*NIS)
          END DO
        END DO
        NMKRHS = iPARDIV(NMKRHS,2*NOMAX**2)
        NMKRHS = 2*NMKRHS
      END IF
      NMKRHS=NMKRHS+NBOTTOM

C PCG/(new)SIGMA routine needs: twice a global array RHS size for
C transformations (overrides NSIGMA computed above)
      NSIGMA_OUTER=0
      NSIGMA_INNER=0
      NVCUTIL=0
      DO ICASE=1,13
      DO ISYM=1,NSYM
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
            NRHSP = iPARDIV(NAS*NIS,0)
      IF (ICASE.GT.11) THEN
        NSIGMA_INNER=MAX(NSIGMA_INNER,NRHSP)
      ELSE
        NSIGMA_INNER=MAX(NSIGMA_INNER,NAS*NIS+NRHSP)
        NSIGMA_OUTER=MAX(NSIGMA_OUTER,NAS*NIS+NRHSP)
      END IF
      NVCUTIL=MAX(NVCUTIL,2*NRHSP)
      END DO
      END DO
      NSIGMA=MAX(NSIGMA_INNER+NSIGMA_OUTER,NVCUTIL)
      NSIGMA=NSIGMA+NBOTTOM

C PRPCTL needs:
C In DIADNS alone, NDD words are needed:
#ifdef _DEBUGPRINT_
      WRITE(6,*)' Memory requirements for PRPCTL (Above SGUGA).'
      WRITE(6,*)
      WRITE(6,*)' PRP1) First phase of PRPCTL.'
      WRITE(6,*)
      WRITE(6,*)'    A) DIADNS alone, broken up in case/symm:'
#endif
      MMX=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 11
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NV=NAS*NIS
          IF(NV.EQ.0) GOTO 11
          IF(ICASE.GT.11) THEN
            M=2*NV
          ELSE
            M=3*NV+(NAS*(NAS+1))/2
          END IF
          MMX=MAX(M,MMX)
  11      CONTINUE
        END DO
      END DO

      NDD=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NX=0
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            NIS=NISUP(ISYM,ICASE)
            IF(NIS.GT.0) THEN
              NAS=NASUP(ISYM,ICASE)
#ifdef _DEBUGPRINT_
              WRITE(6,*)' Case, Symm:',ICASE,ISYM
              WRITE(6,*)' NIN,NAS,NIS:',NIN,NAS,NIS
              WRITE(6,*)' NIMX,NSMX:',NIMX,NAMX
#endif
              IF(ICASE.EQ.2 .OR. ICASE.EQ.3) NX=NIN*NIMX**2
              IF(ICASE.EQ.6 .OR. ICASE.EQ.7) NX=NIN*NSMX*NIMX**2
              IF(ICASE.EQ.8 .OR. ICASE.EQ.9) NX=NIN*NSMX**2
              IF(ICASE.EQ.10 .OR. ICASE.EQ.11) NX=NIN*NIMX*NSMX**2
              IF(ICASE.EQ.12 .OR. ICASE.EQ.13) THEN
                NX=MAX(NAS*NIMX**2,NIS*NSMX**2)
              END IF
            END IF
          END IF
#ifdef _DEBUGPRINT_
      WRITE(6,'(1x,a,2i4,5x,i8)')'      Case, symm:',ICASE,ISYM,2*NX
#endif
          NDD=MAX(2*NX,NDD)
        END DO
      END DO

#ifdef _DEBUGPRINT_
      WRITE(6,*)'       DIADNS needs the maximum, or NDD=',NDD
#endif
      NCMO=NBSQT
      NPRP1=NBOTTOM+NCMO+notri+NLSTOT+2*NOSQT+MMX+NDD

#ifdef _DEBUGPRINT_
      WRITE(6,*)
      WRITE(6,*)'    B) Also needed, for 1st phase of PRPCTL:'
      WRITE(6,'(1x,a,i8)')'       NBOTTOM:',NBOTTOM
      WRITE(6,'(1x,a,i8)')'       NCMO   :',NCMO
      WRITE(6,'(1x,a,i8)')'       notri :',notri
      WRITE(6,'(1x,a,i8)')'       NLSTOT :',NLSTOT
      WRITE(6,'(1x,a,i8)')'     2*NOSQT  :',2*NOSQT
      WRITE(6,'(1x,a,i8)')'       MMX    :',MMX
#endif
      NPRP2=NBOTTOM+2*NCMO+notri+NBAST*(NBAST+1)

#ifdef _DEBUGPRINT_
      WRITE(6,*)' PRP2) Second phase of PRPCTL.'
      WRITE(6,*)
      WRITE(6,'(1x,a,i8)')'       NBOTTOM:',NBOTTOM
      WRITE(6,'(1x,a,i8)')'     2*NCMO   :',2*NCMO
      WRITE(6,'(1x,a,i8)')'       notri :',notri
      WRITE(6,'(1x,a,i8)')'NBAST*(NBAST+1)',NBAST*(NBAST+1)
#endif

      ! (rough) memory estimation for gradients
      ! Cholesky vectors try to use all the available memory, so the
      ! memory used for them are not counted (as in other steps)
      ngrad = 0
      if (do_grad) then
        ! memory in caspt2_grad.f
        nCLag = nConf*nState
        nOLag = NBSQT
        NSLag = nState*nState
        nWLag = NBTRI
        ngrad1 = NBSQT*4 + nCLag*2 + nOLag*2 + nSLag*2 + nWLag
     *         + NBSQT*2 + nAshT**2
        if (IFXMS .or. IFRMS) ngrad1 = ngrad1 + NBSQT
        if (IFDW .and. zeta >= 0.0d+00) ngrad1 = ngrad1 + nState
        if (do_nac) ngrad1 = ngrad1 + NBSQT
        if (nFroT /= 0) ngrad1 = ngrad1 + nFroT**2

        ! "base" memory in dens.f
        nch=0
!       If (IfChol) nch=nvloc_chobatch(1)
        ngrad2 = NOSQT*2 + NBSQT*9 + 2*max(nbast**2,nch) + 2*nAshT**2
        If (nFroT.ne.0 .or. .not.if_invaria) ngrad2 = ngrad2 + 2*NBSQT
        if (do_csf) ngrad2 = ngrad2 + NBSQT
        ngrad2 = ngrad2 + nAshT**2
        if (.not.if_invaria) ngrad2 = ngrad2 + nAshT**2

        membase = MXLEFT - ngrad1 - ngrad2 - NBOTTOM

        ! gradient: residual
        ngrad3 = NSIGMA-NBOTTOM

        ! gradient: density (trdns2o)
        ngrad4 = NDD*3/2

        ! gradient: off-diagonal derivatives (sigder)
        NFIT=0
        NFIA=0
        NFTA=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NS=NSSH(ISYM)
          NFIT=NFIT+NA*NI
          NFIA=NFIA+NS*NI
          NFTA=NFTA+NS*NA
        END DO
        NFIT=NFIT+1
        NFIA=NFIA+1
        NFTA=NFTA+1

        MMX = 0
        DO ICASE1=1,11
          DO ISYM1=1,NSYM
            NIS1=NISUP(ISYM1,ICASE1)
            NAS1=NASUP(ISYM1,ICASE1)
            NIN1=NINDEP(ISYM1,ICASE1)
            NSGM2=NIS1*NAS1
            NSGM1=1
            IF(ICASE1.EQ.1) THEN
              NSGM1=NASH(ISYM1)*NISH(ISYM1)
            ELSE IF(ICASE1.EQ.4) THEN
              NSGM1=NASH(ISYM1)*NSSH(ISYM1)
            ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
              NSGM1=NIS1
            END IF
            M11 = NSGM2 + NSGM1
            M12 = 2*iPARDIV(NSGM2,0) + NSGM2
            IF (ICASE1.LE.11) THEN
              M12 = M12 + iPARDIV(MAX(NSGM2,NAS1*(NAS1+1)/2),0)
            ELSE
              M12 = M12 + iPARDIV(NSGM2,0)
            END IF
            DO ICASE2=ICASE1+1,NCASES
              IF (IFCOUP(ICASE2,ICASE1)==0) CYCLE
              DO ISYM2=1,NSYM
                NIS2=NISUP(ISYM2,ICASE2)
                NAS2=NASUP(ISYM2,ICASE2)
                NCX=NIS2*NAS2
                ! first loop
                M21 = M11 + 2*iPARDIV(NCX,0)
                M31 = M11 - NSGM1 + iPARDIV(NCX,0)
                IF (iCASE1.LE.11) THEN
                  ! C1S1DER
                  MC1S1DER = iPARDIV(NAS1*NAS1+NIN1*NIS1+NAS1*NIN1,0)
                  M31 = M31 + iPARDIV(NSGM2,0)
     *                + MAX(MC1S1DER+NAS1*NAS1,iPARDIV(NSGM2,0))
                END IF

                ! second loop
                M22 = NSGM1 + iPARDIV(NCX,0)
                IF (ICASE2.LE.11) THEN
                  ! C2DER
                  MC2DER = iPARDIV(NAS2*(2*NIS2+MAX(NAS2,NIS2)),0)
                  M22 = M22 + NAS2*NAS2 + MC2DER
                END IF
                M = MAX(M21,M31,M12,M22)
                MMX=MAX(M,MMX)
              END DO
            END DO
          END DO
        END DO
        ngrad5 = NFIT*2 + NFIA*2 + NFTA*2 + MMX

        ! gradient: CI derivatives (clagx.f)
        ngrad6 = NG1*4 + NG2*4 + NG3TOT*3
        ! CLagD (hmm, not for the A and C subspaces in parallel)
        MMX = 0
        DO ICASE=1,11
          DO ISYM=1,NSYM
            NIS=NISUP(ISYM,ICASE)
            NAS=NASUP(ISYM,ICASE)
            NIN=NINDEP(ISYM,ICASE)
            M = 2*NAS**2 + 3*NIN*NIS + NAS*NIS
            IF (IFMSCOUP) M = M + NIN*NIS
            ! inside CLagDX
            M = M + 2*NAS**2 + NAS*MIN(NAS,NIS) + NAS*NIN + NIN
            IF (ICASE.EQ.1) THEN
              M = M + NAS*(NAS+1)/2 + NG3TOT
            ELSE IF (ICASE.EQ.2 .OR. ICASE.EQ.3) THEN
              M = M + 2*NASHT**4
              IF (IPEA_SHIFT /= 0.0D+00) M = M + NAS*(NAS+1)/2
            ELSE IF (ICASE.EQ.4) THEN
              M = M + NAS*(NAS+1)/2 + NG3TOT
            ELSE IF (ICASE.EQ.5) THEN
              IF (IPEA_SHIFT /= 0.0D+00) M = M + NAS*(NAS+1)/2
            ELSE IF (ICASE.EQ.6 .OR. ICASE.EQ.7) THEN
              IF (IPEA_SHIFT /= 0.0D+00) M = M + NAS*(NAS+1)/2
            ELSE IF (ICASE.EQ.8 .OR. ICASE.EQ.9) THEN
              IF (IPEA_SHIFT /= 0.0D+00) M = M + NAS*(NAS+1)/2
            ELSE IF (ICASE.EQ.10 .OR. ICASE.EQ.11) THEN
              IF (IPEA_SHIFT /= 0.0D+00) M = M + NAS*(NAS+1)/2
            END IF
            MMX = MAX(M,MMX)
          END DO
        END DO
        ngrad6_1 = MMX
        ! CnstCLag (derfg3)
        ngrad6_2 = NG3TOT + NCONF
        MXLFT = membase - ngrad6 - ngrad6_2
        MXLFT = MIN(MXCI*(3*NASHT**2+NASHT),MXLFT)
        ngrad6_2 = ngrad6_2 + MXLFT
        ngrad6 = ngrad6 + MAX(ngrad6_1,ngrad6_2)

        ! gradient: effective Hamiltonian (DerHEff)
        NTG1=NASHT**2
        NTG2=NASHT**4
        NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
        ngrad7 = NTG1 + NTG2 + NTG3
        ! DerHeffX
        MAXAIS = 0
        DO ICASE=1,13
          DO ISYM=1,NSYM
            NIS=NISUP(ISYM,ICASE)
            NAS=NASUP(ISYM,ICASE)
            NIN=NINDEP(ISYM,ICASE)
            M = 2*iPARDIV(NAS*NIS,0)
            MAXAIS = MAX(M,MAXAIS)
          END DO
        END DO
        ngrad7_1 = MAXAIS
        ! DERTG3
        ngrad7_2 = MXCI*3 + 2*NASHT**2
        MXLFT = membase - ngrad7 - ngrad7_2
        MXLFT = MIN(MXCI*(2*NASHT**2+4),MXLFT)
        ngrad7_2 = ngrad7_2 + MXLFT
        ngrad7 = ngrad7 + MAX(ngrad7_1,ngrad7_2)

        ! gradient: iterative CI derivatives (DEPSAOffC)
        ngrad8 = 0
        if (.not.if_invar) then
          ngrad8 = nConf*nState*5 + nConf + nState**3 + nAshT**2
     &           + nAshT**4
        end if

        ! gradient: kappa (orbital) derivatives
        !           (OLagNS_RI, OLagNS2, OLagVVVO)
        NumChT = 0
        if (IFCHOL) then
          call Get_iArray('NumCho',NumCho,nSym)
          NumChT = sum(NumCho(1:nSym))
          ngrad9 = NumChT**2 + 1
        else
          ngrad9 = (NISHT+NASHT)**2*NBAST**2 + 1
        end if

        ! gradient: state-specific density matrix
        ngrad10 = 0
        if (if_SSDM) then
          ! IFCHOL is always true
          ngrad10 = NCONF + NumChT**2 + 2*nBasT**2
     &           + 2*NumChT + NBSQT
        end if

        ! memory in XMS (XMS_Grad)
        ngrad11 = 0
        IF (IFMSCOUP) THEN
          ngrad11_1 = NCONF*4 + 2*NASHT**2 + 2*NASHT**4 + NASHT**6
          ngrad11_2 = 3*NASHT**2 + NBSQT*4 + NCONF
          ngrad11 = MAX(ngrad11_1,ngrad11_2)
        END IF

        !! NPRP1 corresponds to trdns2d
        ngrad = NBOTTOM + ngrad1 + ngrad2 + max(ngrad3,ngrad4,ngrad5,
     *          ngrad6,ngrad7,ngrad8,ngrad9,ngrad10,ngrad11,
     *          NPRP1-NBOTTOM)

#ifdef _DEBUGPRINT_
        WRITE(6,*)' PRP3) Gradient.'
        WRITE(6,*)
        WRITE(6,'(1x,a,i12,x,"(essential)")')'NBOTTOM        :',NBOTTOM
        WRITE(6,'(1x,a,i12,x,"(essential)")')'caspt2_grad.f  :',ngrad1
        WRITE(6,'(1x,a,i12,x,"(essential)")')'dens.f         :',ngrad2
        WRITE(6,'(1x,a,i12)')'(Part of) PRP1 :',NPRP1-NBOTTOM
        WRITE(6,'(1x,a,i12)')'caspt2_res.f   :',ngrad3
        WRITE(6,'(1x,a,i12)')'trnds2o.f      :',ngrad4
        WRITE(6,'(1x,a,i12)')'sigder.f       :',ngrad5
        WRITE(6,'(1x,a,i12)')'clagx.f        :',ngrad6
        WRITE(6,'(1x,a,i12)')'derheff.f      :',ngrad7
        WRITE(6,'(1x,a,i12)')'DEPSAOffC      :',ngrad8
        WRITE(6,'(1x,a,i12)')'OLagNS/OLagVVVO:',ngrad9
        WRITE(6,'(1x,a,i12)')'CnstAB_SSDM    :',ngrad10
        WRITE(6,'(1x,a,i12)')'XMS_Grad       :',ngrad11
        WRITE(6,'(1x,a,a )')'----------------','------------'
        WRITE(6,'(1x,a,i12)')'Total Estimate :',ngrad
#endif
      end if

      NPRP=0
      IF (IFPROP .OR. do_grad) NPRP=MAX(NPRP1,NPRP2,ngrad)
      IF ( IPRGLB.GE.USUAL) THEN
        WRITE(6,'(20A4)')('----',I=1,20)
        WRITE(6,*)'Estimated memory requirements:'
        WRITE(6,'(a,i12)')'  POLY3 :             ',NPOLY
        WRITE(6,'(a,i12)')'  RHS:                ',NMKRHS
        WRITE(6,'(a,i12)')'  SIGMA :             ',NSIGMA
        WRITE(6,'(a,i12)')'  PRPCTL:             ',NPRP
*SVC: NPRP includes NDD if it is needed, so this is confusing
*       WRITE(6,'(a,i12)')'  DIADNS:             ',NDD
        WRITE(6,'(a,i12)')' Available workspace: ',MXLEFT
        WRITE(6,*)
      ENDIF

      NEED0=MAX(NPOLY,NMKRHS,NSIGMA)
      NEED=MAX(NEED0,NPRP)
      IF(NEED.GT.MXLEFT) THEN
       IF(NEED0.LE.MXLEFT) THEN
        WRITE(6,'(5X,A)') repeat('*',26)
        WRITE(6,'(5X,A)')' Memory problem!! The memory is insufficient '
        WRITE(6,'(5X,A)')' for the property section.'
        WRITE(6,'(5X,A,I5,A)')'* Need at least ',2+NEED/119000,' MB *'
        IF (do_grad) THEN
          WRITE(6,'(5X,A,A)')' The property (including gradient)',
     *                       ' section will be skipped.'
          WRITE(6,'(5X,A,A)')' If you cannot increase the memory,',
     *                       ' consider using numerical gradients'
          IF (IFPROP) WRITE(6,'(5X,A)')' (without PROPerty keyword)'
        ELSE
          WRITE(6,'(5X,A)')' The property section will be skipped.'
        END IF
        WRITE(6,'(5X,A)') repeat('*',26)
        IF (do_grad) CALL Quit(_RC_MEMORY_ERROR_)
        NEED=NEED0
        IFPROP=.False.
       ELSE
        IF (.NOT.IFCHOL) THEN
C not a Cholesky calculation, keep old memory requirements
         WRITE(6,'(5X,A)') repeat('*',26)
         WRITE(6,'(5X,A)')'* Insufficient memory !! *'
         WRITE(6,'(5X,A,I5,A)')'* Need at least ',2+NEED0/119000,' MB *'
         WRITE(6,'(5X,A)') repeat('*',26)
         WRITE(6,*)
         WRITE(6,*)' Program execution stops -- sorry!'
         WRITE(6,*)
         CALL Quit(_RC_MEMORY_ERROR_)
        ELSE
C This is a Cholesky calculation, only give recommended amount
         IF (MAX(NPOLY,NSIGMA).LE.MXLEFT) THEN
          WRITE(6,'(5X,A)') repeat('*',40)
          WRITE(6,'(5X,A,I5,A)')'* Below comfortable memory of ',
     &                           2+NEED0/119000,' MB *'
          WRITE(6,'(5X,A)') repeat('*',40)
          WRITE(6,*)
          WRITE(6,*)' Program will try to adapt'
          WRITE(6,*)' If it fails, please raise the memory to at least',
     &            2+MAX(NPOLY,INT(0.55D0*NMKRHS),NSIGMA)/119000,' MB'
          WRITE(6,*)' (Maybe more, this aint rocket science)'
          WRITE(6,*)
          WRITE(6,*)' With print level DEBUG you will get some more'
          WRITE(6,*)' informative memory estimates during computation'
          WRITE(6,*)
         ELSE
          WRITE(6,'(5X,A)') repeat('*',26)
          WRITE(6,'(5X,A)')'* Insufficient memory... *'
          WRITE(6,'(5X,A)') repeat('*',26)
          WRITE(6,*)
          WRITE(6,'(A,I6,A)')' If possible, I would like to have   ',
     &                NEED0/119000, ' MB'
          WRITE(6,'(A,I6,A)')' Please raise the memory to at least ',
     &            MAX(NPOLY,INT(0.55D0*NMKRHS),NSIGMA)/119000,' MB'
          WRITE(6,*)' (Maybe more, this aint rocket science)'
          WRITE(6,*)
          CALL Quit(_RC_MEMORY_ERROR_)
         END IF
        END IF
       END IF
      END IF

      END SUBROUTINE SIZES
