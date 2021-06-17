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
      SUBROUTINE DENS(IVEC,JVEC,DMAT)
C
      USE CHOVEC_IO
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"
#include "para_info.fh"

#include "caspt2_grad.fh"
#include "csfbas.fh"
C
#include "pt2_guga.fh"
C
#include "chocaspt2.fh"
      DIMENSION DMAT(*)
      Logical   INVAR
      Character*4096 RealName

      CALL QENTER('DENS')
        CALL GETMEM('G1','ALLO','REAL',LG1,NG1)
        CALL GETMEM('G2','ALLO','REAL',LG2,NG2)
        CALL GETMEM('F1','ALLO','REAL',LF1,NG1)
        CALL PT2_GET(NG1,' GAMMA1',WORK(LG1))
        CALL PT2_GET(NG2,' GAMMA2',WORK(LG2))
        CALL PT2_GET(NG1,' DELTA1',WORK(LF1))
      !! CASPT2 is invariant with respect to rotations in active?
      INVAR=.TRUE.
      If (BSHIFT.NE.0.0D+00.OR.NRAS1T+NRAS3T.NE.0) INVAR=.FALSE.
      !! ???
      invar=.true.
      if (bshift.ne.0.0d+00) invar=.false.
      IF (.not.IFINVAR) INVAR = .FALSE.
C     invar=.false.
C
      If (.not.INVAR .and. IPRGLB.GE.USUAL) Then
        Write (6,*)
        Write (6,'(3X,"This is a non-invariant CASPT2 calculation")')
        If (BSHIFT.NE.0.0D+00)
     *    Write (6,'(3X,"- IPEA shift is employed")')
        If (NRAS1T+NRAS3T.ne.0)
     *    Write (6,'(3X,"- RAS reference is employed")')
        Write (6,'(3X,"A linear equation will be solved to obtain ",
     *                "off-diagonal active density")')
        Write (6,*)
      End If
C       do it = 1, nasht
C       do iu = 1, nasht
C       vvv = 0.0d+00
C       do iv = 1, nasht
C         vvv = vvv + work(lg2+it-1+nasht*(iu-1)
C    *  +nasht**2*(iv-1)+nasht**3*(iv-1))*epsa(iv)
C       end do
C       write(6,'(2i3,3f20.10)') it,iu,vvv,work(lf1+it-1+nasht*(iu-1)),
C    *   work(lg1+it-1+nasht*(iu-1))
C       end do
C       end do

c       easum = 0.0d+00
c       do it = 1, nasht
c         easum = easum + epsa(it)*work(lg1+it-1+nasht*(it-1))
c       end do
c       write(6,*) "easum = ", easum
c       do it = 1, nasht
c       do iu = 1, nasht
c         bvalue = 0.0d+00
c         do iv = 1, nasht
c           bvalue = bvalue + work(lg2+it-1+nasht*(iu-1)
c    *    +nasht**2*(iv-1)+nasht**3*(iv-1))*epsa(iv)
c         end do
c         bvalue = bvalue - easum*work(lg1+it-1+nasht*(iu-1))
c
c         write(6,'(2i3,2f20.10)')
c    *      it,iu,work(lg1+it-1+nasht*(iu-1)),bvalue
c       end do
c       end do
        CALL GETMEM('G1','FREE','REAL',LG1,NG1)
        CALL GETMEM('G2','FREE','REAL',LG2,NG2)
        CALL GETMEM('F1','FREE','REAL',LF1,NG1)
C     do icase = 1, 13
C       write(6,*) "-------------"
C       write(6,*) "icase  = ",icase
C       write(6,*) "nasup  = ", nasup(1,icase)
C       write(6,*) "nisup  = ", nisup(1,icase)
C       write(6,*) "nindep = ", nindep(1,icase)
C       write(6,*) "nexc   = ", nexc(1,icase)
C       write(6,*) "nexces = ", nexces(1,icase)
C     end do
C Compute total density matrix as symmetry-blocked array of
C triangular matrices in DMAT. Size of a triangular submatrix is
C  (NORB(ISYM)*(NORB(ISYM)+1))/2.
      NDMAT=0
      NDPT=0
      nDPTAO=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        nAO = nBas(iSym)
        NDPT=NDPT+NO**2
        NDMAT=NDMAT+(NO*(NO+1))/2
        nDPTAO = nDPTAO + nAO**2
      END DO
      CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
C First, put in the reference density matrix.
      IDMOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO II=1,NI
          IDM=IDMOFF+(II*(II+1))/2
          DMAT(IDM)=2.0D0
        END DO
        DO IT=1,NA
          ITABS=NAES(ISYM)+IT
          ITTOT=NI+IT
          DO IU=1,IT
            IUABS=NAES(ISYM)+IU
            IUTOT=NI+IU
            IDRF=(ITABS*(ITABS-1))/2+IUABS
            IDM=IDMOFF+((ITTOT*(ITTOT-1))/2+IUTOT)
            DMAT(IDM)=WORK(LDREF-1+IDRF)
          END DO
        END DO
         IDMOFF=IDMOFF+(NO*(NO+1))/2
      END DO
*     write(6,*)' DENS. Initial DMAT:'
*     WRITE(*,'(1x,8f16.8)')(dmat(i),i=1,ndmat)
C Add the 1st and 2nd order density matrices:
      CALL GETMEM('DPT','ALLO','REAL',LDPT,NDPT)
      CALL GETMEM('DSUM','ALLO','REAL',LDSUM,NDPT)
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDSUM),1)

C The 1st order contribution to the density matrix
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
      !! working on H only now
      !  CALL TRDNS1(IVEC,WORK(LDPT))
      !  CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*     write(6,*)' DPT after TRDNS1.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
C
C
C
      !! Modify the solution (T; amplitude), if the real- or imaginary-
      !! shift is utilized. We need both the unmodified (T) and modified
      !! (T+\lambda) amplitudes. \lambda can be obtained by solving the
      !! CASPT2 equation, but it can alternatively obtained by a direct
      !! summation only if CASPT2-D.
      !! iVecX remains unchanged (iVecX = T)
      !! iVecR will be 2\lambda
      Call CASPT2_Res
C
C
C
      CALL TRDNS2D(IVEC,JVEC,WORK(LDPT),NDPT)
      IF(IFDENS) THEN
C The exact density matrix evaluation:
C       CALL TRDTMP(WORK(LDPT))
      ELSE
C The approximate density matrix evaluation:
C       CALL TRDNS2A(IVEC,IVEC,WORK(LDPT))
      END IF
      CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
*     write(6,*)' DPT after TRDNS2D.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
      IF (MAXIT.NE.0) THEN
        !! off-diagonal are ignored for CASPT2-D
        CALL DCOPY_(NDPT,[0.0D0],0,WORK(LDPT),1)
C       CALL TRDNS2O(IVEC,IVEC,WORK(LDPT))
        CALL TRDNS2O(iVecX,iVecR,WORK(LDPT))
        CALL DAXPY_(NDPT,1.0D00,WORK(LDPT),1,WORK(LDSUM),1)
      END IF
*     write(6,*)' DPT after TRDNS2O.'
*     WRITE(*,'(1x,8f16.8)')(work(ldpt-1+i),i=1,ndpt)
C
      !! for analytic gradient
      IF (IFDENS) THEN
        IF (.not.IFSADREF.and.nState.ge.2) Then
          write(6,*)
     *      "Please add SADREF keyword in CASPT2 section",
     *      "This keyword is recommended with state-averaged reference"
C         call abend
        End If
        IF (.not.IFDORTHO.and.BSHIFT.ne.0.0D+00) Then
          write(6,*)
     *      "It seems that DORT keyword is not used, ",
     *      "even though this calculation uses the IPEA shift"
          write(6,*)
     *      "Sometimes, analytic gradients do not agree ",
     *      "with numerical gradients"
          write(6,*)
     *      "(which are correct?)"
        End If
C
        !! D^PT in MO
        CALL GETMEM('DPT   ','ALLO','REAL',ipDPT   ,nDPTAO)
        !! D^PT(C) in MO
        CALL GETMEM('DPTC  ','ALLO','REAL',ipDPTC  ,nDPTAO)
        !! DPTAO1 (D^PT in AO, but not DPTA-01) couples with
        !! the CASSCF density (assume state-averaged) through ERIs.
        !! This density corresponds to the eigenvalue derivative.
        !! This is sometimes referred to as DPT2(AO) else where.
        CALL GETMEM('DPTAO ','ALLO','REAL',ipDPTAO ,nDPTAO)
        !! DPTAO2 couples with the inactive density.
        !! This density comes from derivative of the generalized
        !! Fock matrix (see for instance Eq. (24) in the 1990 paper).
        !! This is sometimes referred to as DPT2C(AO) else where.
        CALL GETMEM('DPTCAO','ALLO','REAL',ipDPTCAO,nDPTAO)
        !! DPTAO,DPTCAO,FPTAO,FPTCAO are in a block-squared form
        CALL GETMEM('FPT   ','ALLO','REAL',ipFPT   ,nDPTAO)
        CALL GETMEM('FPTC  ','ALLO','REAL',ipFPTC  ,nDPTAO)
        CALL GETMEM('FPTAO ','ALLO','REAL',ipFPTAO ,nDPTAO)
        CALL GETMEM('FPTCAO','ALLO','REAL',ipFPTCAO,nDPTAO)
        !! Transformation matrix
        CALL GETMEM('TRFMAT','ALLO','REAL',ipTrf   ,nBsqT)
        nch=0
        If (IfChol) nch=nvloc_chobatch(1)
        CALL GETMEM('WRK1  ','ALLO','REAL',ipWRK1  ,Max(nBasT**2,nch))
        CALL GETMEM('WRK2  ','ALLO','REAL',ipWRK2  ,Max(nBasT**2,nch))
        !! FIFA and FIMO (due to frozen orbitals)
        CALL GETMEM('FIFA  ','ALLO','REAL',ipFIFA  ,nBsqT)
        CALL GETMEM('FIMO  ','ALLO','REAL',ipFIMO  ,nBsqT)
        !! state-averaged density
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA ,nAshT*nAshT)
        !! Derivative of state-averaged density
        CALL GETMEM('RDMEIG','ALLO','REAL',ipRDMEIG,nAshT*nAshT)
C       write(6,*) "olag before"
C       call sqprt(Work(ipolag),nbast)
C
        Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPT),1)
        Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPTC),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipDPTAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipDPTCAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPT),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTC),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTAO),1)
        Call DCopy_(nDPTAO,[0.0D+00],0,Work(ipFPTCAO),1)
        Call DCopy_(nBsqT ,[0.0D+00],0,Work(ipFIFA),1)
        Call DCopy_(nBsqT ,[0.0D+00],0,Work(ipFIMO),1)
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipRDMSA),1)
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipRDMEIG),1)
C
        If (nFroT.ne.0) Then
          CALL GETMEM('DIA   ','ALLO','REAL',ipDIA ,nBsqT)
          CALL GETMEM('DI    ','ALLO','REAL',ipDI  ,nBsqT)
        End If
C
        !! Work(LDPT) -> Work(ipDPT2)
        !! Note that Work(ipDPT2) has the index of frozen orbitals.
        !! Note also that unrelaxed (w/o Z-vector) dipole moments with
        !! frozen orbitals must be wrong.
C       call dcopy_(ndpt,[0.0d+00],0,work(ldpt),1)
        If (nFroT.eq.0) Then
          Call DCopy_(nOsqT,Work(LDSUM),1,Work(ipDPT),1)
        Else
          Call OLagFro0(Work(LDSUM),Work(ipDPT))
        End If
C
        !! Construct the transformation matrix
        !! It seems that we have to transform quasi-canonical
        !! to CASSCF orbitals. The forward transformation has been
        !! done in ORBCTL.
        !!   C(PT2) = C(CAS)*X    ->    C(CAS) = C(PT2)*X^T
        !!   -> L(CAS) = X*L(PT2)*X^T
        !! inactive and virtual orbitals are not affected.
        Call DCopy_(nBsqT,[0.0D+0],0,Work(ipTrf),1)
        Call CnstTrf(Work(LTOrb),Work(ipTrf))
C       call sqprt(work(iptrf),nbast)
C
        !! Construct state-averaged density matrix
        Call DCopy_(nDRef,[0.0D+00],0,Work(ipWRK1),1)
        Do iState = 1, nState
          If (IFSADREF) Then
C           Wgt  = Work(LDWgt+iState-1+nState*(iState-1))
            Wgt  = 1.0D+00/nState
            Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
     *                  Work(ipWRK1),1)
          Else If (iState.eq.jState) Then
            Call DaXpY_(nDRef,1.0D+00,Work(LDMix+nDRef*(iState-1)),1,
     *                  Work(ipWRK1),1)
          End If
        End Do
        Call SQUARE(Work(ipWRK1),Work(ipRDMSA),1,nAshT,nAshT)
C       write(6,*) "state-averaged density matrix"
C       call sqprt(work(iprdmsa),nasht)
C
        !! For CI coefficient derivatives (CLag)
        !! Calculate the configuration Lagrangian
        !! Already transformed to natural (CASSCF) orbital basis
        CALL GETMEM('DEPSA ','ALLO','REAL',ipDEPSA,nAshT*nAshT)
        Call DCopy_(nAshT*nAshT,[0.0D+00],0,Work(ipDEPSA),1)
        !! Derivative of off-diagonal H0 of <Psi1|H0|Psi1>
        IF (MAXIT.NE.0) Call SIGDER(iVecX,iVecR)
        Call CLagX(1,Work(ipCLag),Work(ipTRF),Work(ipDEPSA))
C       call test3_dens(work(ipclag))
        write(6,*) "original depsa"
        call sqprt(work(ipdepsa),nasht)
        write(6,*) "original depsa (sym)"
          do i = 1, nasht
          do j = 1, i-1
            val =(work(ipdepsa+i-1+nasht*(j-1))
     *           +work(ipdepsa+j-1+nasht*(i-1)))*0.5d+00
            work(ipdepsa+i-1+nasht*(j-1)) = val
            work(ipdepsa+j-1+nasht*(i-1)) = val
          end do
          end do
        call sqprt(work(ipdepsa),nasht)
C
        If (NRAS1T+NRAS3T.NE.0) Then
          !! The density of the independent pairs (off-diagonal blocks)
          !! should be determined by solving Z-vector, so these blocks
          !! should be removed...?
        write(6,*) "removing DEPSA of off-diagonal blocks"
        write(6,*) "before"
        call sqprt(work(ipdepsa),nasht)
          Do II = 1, nRAS1T
            Do JJ = nRAS1T+1, nAshT
              Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
              Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
            End Do
          End Do
          Do II = nRAS1T+1, nRAS1T+nRAS2T
            Do JJ = nRAS1T+nRAS2T+1, nAshT
              Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
              Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
            End Do
          End Do
        write(6,*) "after"
        call sqprt(work(ipdepsa),nasht)
C       work(ipdepsa+1-1+nasht*(1-1)) = -0.0010235084d+00
C       work(ipdepsa+2-1+nasht*(2-1)) = -0.0005207568d+00
C       work(ipdepsa+3-1+nasht*(3-1)) = -0.0001867875d+00
C       work(ipdepsa+4-1+nasht*(4-1)) = -0.0000110703d+00
C       work(ipclag+2-1) = 0.0002565628d+00
C       work(ipclag+5-1) = 0.0004206981d+00
        End If
C
        !! If CASPT2 energy is not invariant to rotations in active
        !! orbitals, off-diagonal elements of the density obtained
        !! as DEPSA is incorrect, so remove them. The true density
        !! is computed after everything.
        If (.not.INVAR) Then
          !! But, save the diagonal elements
          CALL GETMEM('DEPSA ','ALLO','REAL',ipDEPSAD,nAshT)
          Call DCopy_(nAshT,Work(ipDEPSA),nAshT+1,Work(ipDEPSAD),1)
          !! Clear
          Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipDEPSA),1)
        End If
C       call dscal_(25,2.0d+00,work(ipdepsa),1)
C       write(6,*) "2*depsa"
C       call sqprt(work(ipdepsa),5)
C       call dscal_(25,0.5d+00,work(ipdepsa),1)

          do i = 1, nasht
          do j = 1, nasht
C          if (i.eq.j) cycle
C          work(ipdepsa+i-1+nasht*(j-1))= 0.0d+00
          end do
          end do
C     work(ipdepsa+ 3-1) =  - 0.0001234d+00
C     work(ipdepsa+ 4-1) =  - 0.0000048d+00
C     work(ipdepsa+10-1) =  + 0.0001001d+00
C     work(ipdepsa+11-1) =  - 0.0001234d+00
C     work(ipdepsa+14-1) =  - 0.0001061d+00
C     work(ipdepsa+16-1) =  - 0.0000048d+00
C     work(ipdepsa+18-1) =  - 0.0001061d+00
C     work(ipdepsa+22-1) =  + 0.0001001d+00
        !! asdf test
        !***************************
C       Call CLagFinal(Work(ipCLag),Work(ipSLag))
        !***************************
C
        !! Transform the quasi-variational amplitude (T+\lambda/2?)
        !! in SR (iVecX) to C (iVecC2)
        If (SHIFT.ne.0.0D+00.or.SHIFTI.ne.0.0D+00) Then
          CALL PTRTOC(1,iVecX,iVecC2)
        End If
C         ipTrfL = ipTrf+nAshT*nBasT+nAshT
C         Call DGemm_('n','N',nAshT,nAshT,nAshT,
C    *                1.0D+00,Work(ipTrfL),nBasT,Work(ipDEPSA),nAshT,
C    *                0.0D+00,Work(ipdptcao),nAshT)
C         Call DGemm_('N','t',nAshT,nAshT,nAshT,
C    *                1.0D+00,Work(ipdptcao),nAshT,Work(ipTrfL),nBasT,
C    *                0.0D+00,Work(ipDEPSA),nAshT)
C
C     scal= 1.0d+00
C     work(ipdepsa+ 3-1) = work(ipdepsa+ 3-1)*1 + 0.0000612d+00*scal
C     work(ipdepsa+ 4-1) = work(ipdepsa+ 4-1)*1 + 0.0005517d+00*scal
C     work(ipdepsa+10-1) = work(ipdepsa+10-1)*1 + 0.0000000d+00*scal
C     work(ipdepsa+11-1) = work(ipdepsa+11-1)*1 + 0.0000612d+00*scal
C     work(ipdepsa+14-1) = work(ipdepsa+14-1)*1 + 0.0052315d+00*scal
C     work(ipdepsa+16-1) = work(ipdepsa+16-1)*1 + 0.0005517d+00*scal
C     work(ipdepsa+18-1) = work(ipdepsa+18-1)*1 + 0.0052315d+00*scal
C     work(ipdepsa+22-1) = work(ipdepsa+22-1)*1 + 0.0000000d+00*scal
      ! do i = 1,5
      !   do j = 1,5
      !   val = 0.0d+00
      !   if ((i.eq.1.and.j.eq.3) .or. (i.eq.3.and.j.eq.1)) then
      !     val = -0.0001224565d+00
      !     val =  0.0000240d+00
      !     val = -0.0000612d+00 - 0.0000675d+00
      !   else if ((i.eq.1.and.j.eq.4) .or. (i.eq.4.and.j.eq.1)) then
      !     val = -0.0011033251d+00
      !     val =  0.0006138d+00
      !     val = -0.0005517d+00 - 0.0005327d+00
      !   else if ((i.eq.2.and.j.eq.5) .or. (i.eq.5.and.j.eq.2)) then
C     !     val =  0.0006955589d+00
      !     val = -0.0004582d+00
      !     val =  0.0001074d+00
      !   else if ((i.eq.3.and.j.eq.4) .or. (i.eq.4.and.j.eq.3)) then
      !     val = -0.0104630916d+00
      !     val =  0.0056357d+00
      !     val = -0.0052315d+00 - 0.0052882d+00
      !   else
C     !     cycle
      !   end if
C     !   val=val*0.5d+00
C     !   work(ipdepsa+i-1+5*(j-1)) = work(ipdepsa+i-1+5*(j-1))*0
C    *!    + val
      !   end do
      ! end do
C     ! call sqprt(work(ipdepsa),5)
C
C          call dcopy_(144,[0.0d+00],0,work(ipdpt),1)
C       !! Just add DEPSA to DPT2
        Call AddDEPSA(Work(ipDPT),Work(ipDEPSA))
        !! Just transform the density in MO to AO
        CALL DPT2_Trf(Work(LDPT),Work(ipDPTAO),
     *                Work(LCMOPT2),Work(ipDEPSA),
     *                Work(LDSUM))
C       CALL GETMEM('DEPSA ','FREE','REAL',ipDEPSA,nAshT)
        !! Save the AO density
        !! ... write
C
        !! Construct orbital Lagrangian that comes from the derivative
        !! of ERIs. Also, do the Fock transformation of the DPT2 and
        !! DPT2C densities.
C       CALL OLagNS(Work(ipDPTC),Work(ipDPTAO),Work(ipDPTCAO),
C    *              Work(ipFPTAO),Work(ipFPTCAO))
C        write(6,*) "olag"
C       call sqprt(work(ipolag),12)
C     write(6,*) "OLag"
C     do i = 1, 144
C       write(6,'(i3,f20.10)') i,work(ipolag+i-1)
C     end do
C        write(6,*) "dpt2ao ref"
C       call sqprt(work(ipdptao),12)
C        write(6,*) "fpt2ao ref"
C       call sqprt(work(ipfptao),12)

C       call dcopy_(144,[0.0d+00],0,work(ipolag),1)
C       call dcopy_(nbast**2,[0.0d+00],0,work(ipfptao),1)
C       Call DCopy_(nDPTAO,[0.0d+00],0,Work(ipDPTC),1)
C
C
C       write(6,*) "nfrot = ", nfrot
        If (nFroT.ne.0) Then
          !! If frozen orbitals exist, we need to obtain electron-repulsion
          !! integrals with frozen orbitals to construct the orbital
          !! Lagrangian.
          If (.not.IfChol) Call TRAFRO(1)
C
          !! Get density matrix (Work(ipDIA)) and inactive density
          !! matrix (Work(ipDI)) to compute FIFA and FIMO.
          Call OLagFroD(Work(ipDIA),Work(ipDI),Work(ipRDMSA),
     *                  Work(ipTrf))
C         write(6,*) "density matrix"
C         call sqprt(work(ipdia),12)
C         call sqprt(work(ipdi),12)
        End If
C
C
C
        If (.false.) then
C
C
C
C         call dcopy_(144,[0.0d+00],0,work(ipdpt),1)
        do i = 6, 10
          do j = 6, 10
          val = 0.0d+00
          if ((i.eq.6.and.j.eq.8) .or. (i.eq.8.and.j.eq.6)) then
            val = -0.0000264d+00
            val =  0.0004270d+00
            val =  0.0022536d+00
            val = -0.0001392d+00

C           val = -0.0010462839d+00 !! g-T
            val = -0.0002482178d+00 !! g-N
            val = -0.0004752616d+00 !! after all
            val =  0.0000762193d+00 !! after all+depsa
C        val=-0.0001224d+00
         val=-0.0001916426d+00
          else if ((i.eq.6.and.j.eq.9) .or. (i.eq.9.and.j.eq.6)) then
            val = -0.0040663d+00
            val = -0.0040513d+00
            val = -0.0056767d+00
            val = -0.0056977d+00

C           val = -0.0013117791d+00 !! g-T
            val = -0.0013208468d+00 !! g-N
            val = -0.0010672239d+00 !! after all
            val =  0.0000107227d+00 !! after all+depsa
C        val=-0.0001472d+00
         val=-0.0010467525d+00
          else if ((i.eq.7.and.j.eq.10) .or. (i.eq.10.and.j.eq.7)) then
            val = -0.0131419d+00
            val = -0.0139400d+00
            val = -0.0211634d+00
            val = -0.0199518d+00

C           val = -0.0125402532d+00 !! g-T
            val = -0.0128318974d+00 !! g-N
            val = -0.0118058046d+00 !! after all
            val = -0.0000747859d+00 !! after all+depsa
C        val=-0.0125014d+00
         val=-0.0120205513d+00
          else if ((i.eq.8.and.j.eq.9) .or. (i.eq.9.and.j.eq.8)) then
            val =  0.0167787d+00
            val =  0.0160932d+00
            val =  0.0306999d+00
            val =  0.0320075d+00

C           val = -0.0212510757d+00 !! g-T
            val = -0.0212531499d+00 !! g-N
            val = -0.0207823583d+00 !! after all
            val = -0.0000278912d+00 !! after all+depsa
C        val=-0.0100415d+00
         val=-0.0206811594d+00
          else
            cycle
          end if
C         if (i.lt.j) val=-val
        if (.false.) then
           val = -0.25d+00*val
C         work(ipdpt+i-1+12*(j-1)) = work(ipdpt+i-1+12*(j-1))
C    *     + val
C         work(ipdpt+j-1+12*(i-1)) = work(ipdpt+j-1+12*(i-1))
C    *     + val
C         work(ldsum+i-1+12*(j-1)) = work(ldsum+i-1+12*(j-1))
C    *     + val
C         work(ldsum+j-1+12*(i-1)) = work(ldsum+j-1+12*(i-1))
C    *     + val
C      work(iprdmeig+i-6+5*(j-6))=work(iprdmeig+i-6+5*(j-6)) + val
C      work(iprdmeig+j-6+5*(i-6))=work(iprdmeig+j-6+5*(i-6)) + val
        else
           val = -0.5d+00*val
C         work(ipdpt+i-1+12*(j-1)) = work(ipdpt+i-1+12*(j-1))*0
C    *     + val
C         work(ipdpt+j-1+12*(i-1)) = work(ipdpt+j-1+12*(i-1))*0
C    *     + val
C         work(ldsum+i-1+12*(j-1)) = work(ldsum+i-1+12*(j-1))*0
C    *     + val
C         work(ldsum+j-1+12*(i-1)) = work(ldsum+j-1+12*(i-1))*0
C    *     + val
        end if
          end do
        end do
        call dcopy_(25,[0.0d+00],0,Work(ipWRK1),1)
        work(ipwrk1+3-1+5*(1-1)) = 0.0000043d+00
        work(ipwrk1+1-1+5*(3-1)) = 0.0000043d+00
        work(ipwrk1+4-1+5*(1-1)) =-0.0001063d+00
        work(ipwrk1+1-1+5*(4-1)) =-0.0001063d+00
        work(ipwrk1+5-1+5*(2-1)) =-0.0002388d+00
        work(ipwrk1+2-1+5*(5-1)) =-0.0002388d+00
        work(ipwrk1+4-1+5*(3-1)) = 0.0014233d+00
        work(ipwrk1+3-1+5*(4-1)) = 0.0014233d+00
        work(ipwrk1+3-1+5*(1-1)) = 0.0000239518d+00
        work(ipwrk1+1-1+5*(3-1)) = 0.0000239518d+00
        work(ipwrk1+4-1+5*(1-1)) = 0.0006138144d+00
        work(ipwrk1+1-1+5*(4-1)) = 0.0006138144d+00
        work(ipwrk1+5-1+5*(2-1)) =-0.0004581528d+00
        work(ipwrk1+2-1+5*(5-1)) =-0.0004581528d+00
        work(ipwrk1+4-1+5*(3-1)) = 0.0056356566d+00
        work(ipwrk1+3-1+5*(4-1)) = 0.0056356566d+00

        work(ipwrk1+3-1+5*(1-1)) = 0.0000255636d+00
        work(ipwrk1+1-1+5*(3-1)) = 0.0000255636d+00
        work(ipwrk1+4-1+5*(1-1)) = 0.0007603565d+00
        work(ipwrk1+1-1+5*(4-1)) = 0.0007603565d+00
        work(ipwrk1+5-1+5*(2-1)) =-0.0001414500d+00
        work(ipwrk1+2-1+5*(5-1)) =-0.0001414500d+00
        work(ipwrk1+4-1+5*(3-1)) = 0.0055442664d+00
        work(ipwrk1+3-1+5*(4-1)) = 0.0055442664d+00

        !! with diagonal(?)
        work(ipwrk1+3-1+5*(1-1)) = 0.0000270181d+00
        work(ipwrk1+1-1+5*(3-1)) = 0.0000270181d+00
        work(ipwrk1+4-1+5*(1-1)) = 0.0006236086d+00
        work(ipwrk1+1-1+5*(4-1)) = 0.0006236086d+00
        work(ipwrk1+5-1+5*(2-1)) =-0.0004306630d+00
        work(ipwrk1+2-1+5*(5-1)) =-0.0004306630d+00
        work(ipwrk1+4-1+5*(3-1)) = 0.0055647223d+00
        work(ipwrk1+3-1+5*(4-1)) = 0.0055647223d+00
        !! without diagonal
        work(ipwrk1+3-1+5*(1-1)) = 0.0000274314d+00
        work(ipwrk1+1-1+5*(3-1)) =-0.0000274314d+00
        work(ipwrk1+4-1+5*(1-1)) = 0.0006093060d+00
        work(ipwrk1+1-1+5*(4-1)) =-0.0006093060d+00
        work(ipwrk1+5-1+5*(2-1)) =-0.0001112837d+00
        work(ipwrk1+2-1+5*(5-1)) = 0.0001112837d+00
        work(ipwrk1+4-1+5*(3-1)) = 0.0054718711d+00
        work(ipwrk1+3-1+5*(4-1)) =-0.0054718711d+00
        !! without diagonal with IPEA 0.5
C       work(ipwrk1+3-1+5*(1-1)) = 0.0000262871d+00
C       work(ipwrk1+1-1+5*(3-1)) =-0.0000262871d+00
C       work(ipwrk1+4-1+5*(1-1)) = 0.0006068526d+00
C       work(ipwrk1+1-1+5*(4-1)) =-0.0006068526d+00
C       work(ipwrk1+5-1+5*(2-1)) =-0.0001029116d+00
C       work(ipwrk1+2-1+5*(5-1)) = 0.0001029116d+00
C       work(ipwrk1+4-1+5*(3-1)) = 0.0055083426d+00
C       work(ipwrk1+3-1+5*(4-1)) =-0.0055083426d+00
        !! without diagonal with canonical keyword
C       work(ipwrk1+3-1+5*(1-1)) = 0.0000698704d+00
C       work(ipwrk1+1-1+5*(3-1)) = 0.0000698704d+00
C       work(ipwrk1+4-1+5*(1-1)) = 0.0006505549d+00
C       work(ipwrk1+1-1+5*(4-1)) = 0.0006505549d+00
C       work(ipwrk1+5-1+5*(2-1)) =-0.0001112837d+00
C       work(ipwrk1+2-1+5*(5-1)) =-0.0001112837d+00
C       work(ipwrk1+4-1+5*(3-1)) = 0.0054667427d+00
C       work(ipwrk1+3-1+5*(4-1)) = 0.0054667427d+00
        write(6,*) "trf"
        call sqprt(work(iptrf),12)
        write(6,*) "Z in natural"
        call sqprt(work(ipwrk1),5)
        call dgemm_('T','N',5,5,5,
     *             1.0d+00,Work(ipTRF+6-1+12*(6-1)),12,
     *                     Work(ipWRK1),5,
     *             0.0d+00,Work(ipWRK2),5)
        call dgemm_('N','N',5,5,5,
     *             1.0D+00,Work(ipWRK2),5,
     *                     Work(ipTRF+6-1+12*(6-1)),12,
     *             0.0d+00,Work(ipWRK1),5)
        write(6,*) "Z in quasi-canonical"
        call sqprt(work(ipwrk1),5)
C       write(6,*) "D before"
C       call sqprt(work(ipdpt),12)
        do i = 1, 5
C         work(ipwrk1+i-1+5*(i-1)) = 0.0d+00
          ii = i+5
        write(6,*) ii,eps(ii)
          do j = 1, i-1
            if (i.eq.j) cycle
            jj = j+5
            val = work(ipwrk1+i-1+5*(j-1))/(eps(ii)-eps(jj))
            val=val*0.5d+00
             work(ipwrk1+i-1+5*(j-1)) = val
             work(ipwrk1+j-1+5*(i-1)) = val
C     write(6,'(i3,1f20.10)')
C    * indpq,work(ipwrk1+i-1+5*(j-1))/(eps(ii)-eps(jj))
C           work(ipdpt+ii-1+12*(jj-1))
C    *    = work(ipdpt+ii-1+12*(jj-1))*1.0d+00
C    *    + work(ipwrk1+i-1+5*(j-1))*1.0d+00
C           work(ipdptc+ii-1+12*(jj-1))
C    *    = work(ipdptc+ii-1+12*(jj-1))
C    *    - work(ipwrk1+i-1+5*(j-1))*1.0d+00
          end do
        end do
        write(6,*) "Z density"
        call sqprt(work(ipwrk1),5)
C       write(6,*) "D after"
C       call sqprt(work(ipdpt),12)
C
C
C
        end if
C
C
C
C
C
        If (IfChol) Then
          NumChoTot = 0
          Do iSym = 1, nSym
            NumChoTot = NumChoTot + NumCho_PT2(iSym)
          End Do
          Call GetMem('A_PT2 ','ALLO','REAL',ipA_PT2,NumChoTot**2)
          Call dcopy_(NumChoTot**2,[0.0D+00],0,Work(ipA_PT2),1)
        End If
        Do iSym = 1, nSym
          nOcc  = nIsh(iSym)+nAsh(iSym)
          If (.not.IfChol.or.iALGO.ne.1) Then
            lT2AO = nOcc*nOcc*nBasT*nBasT
            Call GetMem('T2AO','Allo','Real',ipT2AO,lT2AO)
            Call DCopy_(lT2AO,[0.0D+00],0,Work(ipT2AO),1)
          End If
C
          !! Orbital Lagrangian that comes from the derivative of ERIs.
          !! OLagNS computes only the particle orbitals.
C         write(6,*) "ialgo = ", ialgo
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          If (IfChol.and.iALGO.eq.1) Then
            !! check the dimension, in particular auxiliary basis
            nChoBuf = (nAsh(iSym)*nIsh(iSym)
     *                +nSsh(iSym)*nIsh(iSym)
     *                +nAsh(iSym)*nAsh(iSym)
     *                +nSsh(iSym)*nAsh(iSym))*NVLOC_CHOBATCH(1)
       write(6,*) "nChoBuf = ", nchobuf,isym
       write(6,*) nish(isym),nash(isym),nssh(isym),NVLOC_CHOBATCH(1)
            Call GetMem('DENBRA','ALLO','REAL',ipDBra,nChoBuf)
            Call dcopy_(nChoBuf,[0.0D+00],0,Work(ipDBra),1)
C           CALL OLagNS_RI(iSym,Work(ipWRK1),Work(ipWRK2),
C    *                     Work(ipDPTC),Work(ipDBra),Work(ipA_PT2),
C    *                     NVLOC_CHOBATCH(1))
            ipAI = ipDBra
            ipSI = ipAI + nAsh(iSym)*nIsh(iSym)*NVLOC_CHOBATCH(1)
            ipAA = ipSI + nSsh(iSym)*nIsh(iSym)*NVLOC_CHOBATCH(1)
            ipSA = ipAA + nAsh(iSym)*nAsh(iSym)*NVLOC_CHOBATCH(1)
            CALL OLagNS_RI(iSym,Work(ipWRK1),Work(ipWRK2),
     *                     Work(ipDPTC),Work(ipAI),Work(ipSI),
     *                     Work(ipAA),Work(ipSA),Work(ipA_PT2),
     *                     NVLOC_CHOBATCH(1))
C           do i = 1, nchobuf
C             write(6,'(i4,f20.10)') i,work(ipdbra+i-1)
C           end do
          Else
            CALL OLagNS2(iSym,Work(ipDPTC),Work(ipT2AO))
          End If
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(6,*) "OLagNS: CPU/WALL TIME=", cput,wallt
C         write(6,*) "DPT2C"
C         call sqprt(work(ipdptc),nbast)
C
          !! MO -> AO transformations for DPT2 and DPT2C
          If ((.not.IfChol.or.iALGO.ne.1).or.nFroT.eq.0) Then
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPT),
     *                   Work(ipDPTAO),Work(ipWRK1))
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPTC),
     *                   Work(ipDPTCAO),Work(ipWRK1))
C           write(6,*) "dpt2"
C           call sqprt(work(ipdpt),nbast)
C           write(6,*) "dpt2ao"
C           call sqprt(work(ipdptao),nbast)
          End If
C
          !! Do some transformations relevant to avoiding (VV|VO)
          !! integrals. Orbital Lagrangian for the hole orbitals are
          !! computed. At the same time, F = G(D) transformations are
          !! also performed for D = DPT2 and DPT2C
          !! The way implemented (what?) is just a shit. I cannot find
          !! FIFA and FIMO for frozen orbitals, so I have to construct
          !! them. Here is the transformation of G(D^inact) and G(D).
          !! Work(ipFIFA) and Work(ipFIMO) computed in this subroutine
          !! is not yet correct. They are just two-electron after this
          !! subroutine.
          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          CALL OLagVVVO(iSym,Work(ipDPTAO),Work(ipDPTCAO),
     *                  Work(ipFPTAO),Work(ipFPTCAO),Work(ipT2AO),
     *                  Work(ipDIA),Work(ipDI),Work(ipFIFA),
     *                  Work(ipFIMO),Work(ipDBra),
     *                  Work(ipA_PT2),NumChoTot)
C       write(6,*) "olag after vvvo"
C       call sqprt(work(ipolag),nbast)
C       call abend
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          write(6,*) "OLagVVVO: CPU/WALL TIME=", cput,wallt
          If (IfChol.and.iALGO.eq.1) Then
            Call GetMem('DENBRA','FREE','REAL',ipDBra,NCHOBUF)
          End If
C     write(6,*) "OLag"
C     do i = 1, 144
C       write(6,'(i3,f20.10)') i,work(ipolag+i-1)
C     end do
C      write(6,*) "fpt2ao"
C     call sqprt(work(ipfptao),12)
C     call abend
C
          !! AO -> MO transformations for FPT2AO and FPT2CAO
          If ((.not.IfChol.or.iALGO.ne.1).or.nFroT.eq.0) Then
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPT),
     *                   Work(ipFPTAO),Work(ipWRK1))
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPTC),
     *                   Work(ipFPTCAO),Work(ipWRK1))
          End If
C
          If (.not.IfChol.or.iALGO.ne.1) Then
            Call GetMem('T2AO','Free','Real',ipT2AO,lT2AO)
          End If
        End Do
        If (IfChol) Then
          Call GetMem('A_PT2 ','FREE','REAL',ipA_PT2,NumChoTot**2)
        End If
C
c
c
C       Call SQUARE(Work(LFIFA),Work(ipFIFA),1,12,12)
C       CALL DGEMM_('N','T',12,12,12,
C    *              2.0D+00,work(ipfifa),12,work(ipdpt),12,
C    *              1.0D+00,Work(ipOLAG),12)
C          call test2_dens(work(ipolag),work(ipdepsa))
C
c
c
C       write(6,*) "fptao after olagns"
C       call sqprt(Work(ipfptao),nbast)
C       write(6,*) "fptcao after olagns"
C       call sqprt(Work(ipfptcao),nbast)
C       write(6,*) "olag after olagns"
C       call sqprt(Work(ipolag),nbast)
C
        !! If frozen orbitals exist, frozen-inactive part of the
        !! unrelaxed PT2 density matrix is computed using the orbital
        !! Lagrangian. Additionally, Fock transformation is also
        !! required.
        If (nFroT.ne.0) Then
          !! Compute DPT2 density for frozen-inactive
C         write(6,*) "dpt before frozen"
C         call sqprt(work(ipdpt),nbast)
          !! Construct FIFA and FIMO
          Call OLagFro3(Work(ipFIFA),Work(ipFIMO),Work(ipWRK1),
     *                  Work(ipWRK2))
          !! Add the FIMO contributions (other eigenvalue derivative
          !! etc. contributions are symmetric in inactive orbitals, so
          !! they do not contribute to frozen density)
          CALL DGEMM_('N','T',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)
          CALL DGEMM_('T','N',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipFIMO),nBasT,Work(ipDPTC),nBasT,
     *                1.0D+00,Work(ipOLAG),nBasT)

          Call OLagFro1(Work(ipDPT),Work(ipOLag),Work(ipTrf))
C         write(6,*) "dpt after frozen"
C         call sqprt(work(ipdpt),nbast)
          !! Fock transformation for frozen-inactive density
          If (IfChol) Then
            iSym=1
            !! MO -> AO transformations for DPT2 and DPT2C
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPT),
     *                   Work(ipDPTAO),Work(ipWRK1))
            Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDPTC),
     *                   Work(ipDPTCAO),Work(ipWRK1))
            !! For DF-CASPT2, Fock transformation of DPT2, DPT2C, DIA,
            !! DA is done here, but not OLagVVVO
            !! It seems that it is not possible to do this
            !! transformation in OLagVVVO, because the frozen-part of
            !! the DPT2 is obtained after OLagVVVO.
            Call OLagFro4(1,1,1,1,1,
     *                    Work(ipDPTAO),Work(ipDPTCAO),Work(ipFPTAO),
     *                    Work(ipFPTCAO),Work(ipWRK1))
            !! AO -> MO transformations for FPT2AO and FPT2CAO
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPT),
     *                   Work(ipFPTAO),Work(ipWRK1))
            Call OLagTrf(2,iSym,Work(LCMOPT2),Work(ipFPTC),
     *                   Work(ipFPTCAO),Work(ipWRK1))
          Else
C           write(6,*) "dpt"
C           call sqprt(work(ipdpt),nbast)
            Call OLagFro2(Work(ipDPT),Work(ipFPT),Work(ipWRK1),
     *                    Work(ipWRK2))
C         write(6,*) "fpt"
C           call sqprt(work(ipdpt),nbast)
          End If
C         write(6,*) "ipfifa"
C         call sqprt(work(ipfifa),12)
C         write(6,*) "ipfimo"
C         call sqprt(work(ipfimo),12)
      !   !! Construct FIFA and FIMO
      !   Call OLagFro3(Work(ipFIFA),Work(ipFIMO),Work(ipWRK1),
     *!                 Work(ipWRK2))
C
          CALL GETMEM('DIA   ','FREE','REAL',ipDIA ,nBsqT)
          CALL GETMEM('DI    ','FREE','REAL',ipDI  ,nBsqT)
        Else
          iSQ = 0
          iTR = 0
          Do iSym = 1, nSym
            nOrbI = nOrb(iSym)
            Call SQUARE(Work(LFIFA+iTR),Work(ipFIFA+iSQ),1,nOrbI,nOrbI)
            Call SQUARE(Work(LFIMO+iTR),Work(ipFIMO+iSQ),1,nOrbI,nOrbI)
            iSQ = iSQ + nOrbI*nOrbI
            iTR = iTR + nOrbI*(nOrbI+1)/2
          End Do
        End If
C         write(6,*) "ipfifa"
C         call sqprt(work(ipfifa),nbast)
C         write(6,*) "ipfimo"
C         call sqprt(work(ipfimo),nbast)
C        write(6,*) "FIFA in quasi-canonical"
C         call sqprt(work(ipfifa),12)
C        write(6,*) "FIFA in natural"
C         Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipTrf),nBasT,Work(ipFIFA),nBasT,
C    *                0.0D+00,Work(ipWRK1),nBasT)
C         Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *                0.0D+00,Work(ipWRK2),nBasT)
C         call sqprt(work(ipwrk2),12)
C
        !! Do some post-process for the contributions that comes from
        !! the above two densities.
C       CALL EigDer(Work(LDPT),Work(ipDPTC),Work(ipFPTAO),
C       write(6,*) "olag before eigder"
C       call sqprt(Work(ipolag),nbast)
C       write(6,*) "fpt2"
C       call sqprt(Work(ipfpt),nbast)
        CALL EigDer(Work(ipDPT),Work(ipDPTC),Work(ipFPTAO),
     *              Work(ipFPTCAO),Work(ipRDMEIG),Work(LCMOPT2),
     *              Work(ipTrf),Work(ipFPT),Work(ipFPTC),
     *              Work(ipFIFA),Work(ipFIMO),Work(ipRDMSA))
C          call test2_dens(work(ipolag),work(ipdepsa))
C       write(6,*) "olag after eigder"
C       call sqprt(Work(ipolag),nbast)
C       write(6,*) "Wlag after eigder"
C       call sqprt(work(ipwlag),nbast)
C       write(6,*) "rdmeig"
C       call sqprt(work(iprdmeig),nasht)
C       call abend
C
        !! Calculate the configuration Lagrangian again.
        !! The contribution comes from the derivative of eigenvalues.
        !! It seems that TRACI_RPT2 uses CI coefficients of RASSCF,
        !! so canonical -> natural transformation is required.
C       ipTrfL = ipTrf+nAshT*nBasT+nAshT
C       Call DGemm_('N','N',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipTrfL),nBasT,Work(ipRDMEIG),nAshT,
C    *              0.0D+00,Work(ipWRK1),nAshT)
C       Call DGemm_('N','T',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipWRK1),nAshT,Work(ipTrfL),nBasT,
C    *              0.0D+00,Work(ipRDMEIG),nAshT)
        If (.not.INVAR) Then !test
          CALL GETMEM('CLagT','ALLO','REAL',ipCLagT,nConf*nState)
          CALL GETMEM('EigT ','ALLO','REAL',ipEigT ,nAshT**2)
          Call DCopy_(nConf*nState,Work(ipCLag),1,Work(ipCLagT),1)
          Call DCopy_(nAshT**2,Work(ipRDMEIG),1,Work(ipEigT),1)
        End If
        !! Use canonical CSFs rather than natural CSFs
        ISAV = IDCIEX
        IDCIEX = IDTCEX
        !! Now, compute the configuration Lagrangian
        Call CLagEig(Work(ipCLag),Work(ipRDMEIG))
        !! Now, compute the state Lagrangian and do some projections
        Call CLagFinal(Work(ipCLag),Work(ipSLag))
C
        !! Now, here is the best place to compute the true off-diagonal
        !! active density for non-invariant CASPT2
        If (.not.INVAR) Then
          !! Add the density that comes from CI Lagrangian
          Call DEPSAOffC(Work(ipCLag),Work(ipDEPSA),Work(ipFIFA),
     *                   Work(ipFIMO),
     *                   Work(ipWRK1),Work(ipWRK2))
          !! Add the density that comes from orbital Lagrangian
          Call DEPSAOffO(Work(ipOLag),Work(ipDEPSA),Work(ipFIFA))
          !! Restore the diagonal elements
          Call DCopy_(nAshT,Work(ipDEPSAD),1,Work(ipDEPSA),nAshT+1)
          CALL GETMEM('DEPSAD','FREE','REAL',ipDEPSAD,nAshT)
          write(6,*) "DEPSA computed again"
          call sqprt(work(ipdepsa),nasht)
          If (NRAS1T+NRAS3T.NE.0) Then
            !! Remove the off-diagonal blocks for RASPT2
            Do II = 1, nRAS1T
              Do JJ = nRAS1T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
            Do II = nRAS1T+1, nRAS1T+nRAS2T
              Do JJ = nRAS1T+nRAS2T+1, nAshT
                Work(ipDEPSA+II-1+nAshT*(JJ-1)) = 0.0D+00
                Work(ipDEPSA+JJ-1+nAshT*(II-1)) = 0.0D+00
              End Do
            End Do
          End If
C         call dcopy_(nasht**2,[0.0d+00],0,work(ipdepsa),1)
C
          !! We have to do many things again...
          !! Just add DEPSA to DPT2
          Call AddDEPSA(Work(ipDPT),Work(ipDEPSA))
          !! Just transform the density in MO to AO
          CALL DPT2_Trf(Work(LDPT),Work(ipDPTAO),
     *                  Work(LCMOPT2),Work(ipDEPSA),
     *                  Work(LDSUM))
          !! Some transformations similar to EigDer
          Call EigDer2(Work(ipRDMEIG),Work(ipTrf),Work(ipFIFA),
     *                 Work(ipRDMSA),Work(ipDEPSA),
     *                 Work(ipWRK1),Work(ipWRK2))
C
          Call DCopy_(nConf*nState,Work(ipCLagT),1,Work(ipCLag),1) !test
          Call DaXpY_(nAshT**2,1.0D+00,Work(ipEigT),1,Work(ipRDMEIG),1) !test
          call DCopy_(nState*(nState-1)/2,[0.0D+00],0,Work(ipSLag),1)
          CALL GETMEM('CLagT','FREE','REAL',ipCLagT,nConf*nState)
          CALL GETMEM('EigT ','FREE','REAL',ipEigT ,nAshT**2)
          !! RDMEIG contributions
          !! Use canonical CSFs rather than natural CSFs
          !! Now, compute the configuration Lagrangian
          Call CLagEig(Work(ipCLag),Work(ipRDMEIG))
          !! Now, compute the state Lagrangian and do some projections
          Call CLagFinal(Work(ipCLag),Work(ipSLag))
        End If
C
        !! Restore integrals without frozen orbitals, although not sure
        !! this operation is required.
        If (nFroT.ne.0.and..not.IfChol) Call TRAFRO(2)
C
        IDCIEX = ISAV
        !! Canonical -> natural transformation
        IF(ORBIN.EQ.'TRANSFOR') Then
          Do iState = 1, nState
            Call CLagX_TrfCI(Work(ipCLag+nConf*(iState-1)))
          End Do
        End If
C       Call CLagFinal(Work(ipCLag),Work(ipSLag))
C
        !! Transformations of DPT2 in quasi-canonical to natural orbital
        !! basis and store the transformed density so that the MCLR
        !! module can use them.
        Call DPT2_TrfStore(1.0D+00,Work(ipDPT),Work(ipDPT2),
     *                     Work(ipTrf),Work(ipWRK1))
        Call DPT2_TrfStore(2.0D+00,Work(ipDPTC),Work(ipDPT2C),
     *                     Work(ipTrf),Work(ipWRK1))
C       !! Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipTrf),nBasT,Work(LDPT),nBasT,
C    *              0.0D+00,Work(ipWRK1),nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *              0.0D+00,Work(ipWRK2),nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C         Call DaXpY_(nSQ,1.0D+00,Work(ipWRK2+iSQ),1,Work(ipDPT2+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C
C       !! Do the same for DPT2C Save MO densities for post MCLR
C       Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipTrf),nBasT,Work(ipDPTC),nBasT,
C    *              0.0D+00,Work(ipWRK1),nBasT)
C       Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *              1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *              0.0D+00,Work(ipWRK2),nBasT)
C       iSQ = 0
C       Do iSym = 1, nSym
C         nOrbI = nBas(iSym)-nDel(iSym)
C         nSQ = nOrbI*nOrbI
C        Call DaXpY_(nSQ,2.0D+00,Work(ipWRK2+iSQ),1,Work(ipDPT2C+iSQ),1)
C         iSQ = iSQ + nSQ
C       End Do
C       call abend()
C       call sqprt(Work(ipRDMEIG),nAshT)
C
        !! square -> triangle so that the MCLR module can use the AO
        !! densities. Do this for DPT2AO and DPT2CAO (defined in
        !! caspt2_grad.f and caspt2_grad.h).
        iBasTr = 1
        iBasSq = 1
        Do iSym = 1, nSym
          nBasI = nBas(iSym)
          liBasTr = iBasTr
          liBasSq = iBasSq
          Do iBasI = 1, nBasI
            Do jBasI = 1, iBasI
              liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
              If (iBasI.eq.jBasI) Then
                Work(ipDPT2AO +liBasTr-1) = Work(ipDPTAO +liBasSq-1)
                Work(ipDPT2CAO+liBasTr-1) = Work(ipDPTCAO+liBasSq-1)
              Else
                Work(ipDPT2AO +liBasTr-1)
     *            = Work(ipDPTAO +liBasSq-1)*2.0D+00
                Work(ipDPT2CAO+liBasTr-1)
     *            = Work(ipDPTCAO+liBasSq-1)*2.0D+00
              End If
              liBasTr = liBasTr + 1
            End Do
          End Do
          iBasTr = iBasTr + nBasI*(nBasI+1)/2
          iBasSq = iBasSq + nBasI*nBasI
        End Do
C
        !! If SS density matrix is used, we need an additional term for
        !! electron-repulsion integral. Here prepares such densities.
        !! The first one is just DPT2AO, while the second one is the
        !! difference between the SS and SA density matrix. because the
        !! SA density-contribution will be added and should be
        !! subtracted
        If (.not.IFSADREF.and.nState.ne.1.and..not.IFXMS) Then
          If (.not.INVAR) Then
            write(6,*) "SS density matrix with BSHIFT is not yet"
            Call abend()
          End If
     *
          !! Subtract the SA-RDM (inactive part is later)
          Call DCopy_(nDRef,[0.0D+00],0,Work(ipWRK1),1)
          Do iState = 1, nState
C           Wgt  = Work(LDWgt+iState-1+nState*(iState-1))
            Wgt  = 1.0D+00/nState
            Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
     *                  Work(ipWRK1),1)
          End Do
          Call SQUARE(Work(ipWRK1),Work(ipWRK2),1,nAshT,nAshT)
          Call DaXpY_(nAshT**2,-1.0D+00,Work(ipWRK2),1,Work(ipRDMSA),1)
          !! Construct the SS density matrix in Work(ipWRK1)
          Call OLagFroD(Work(ipWRK1),Work(ipWRK2),Work(ipRDMSA),
     *                  Work(ipTrf))
          !! Subtract the inactive part
          Call DaXpY_(nBasT**2,-1.0D+00,Work(ipWRK2),1,Work(ipWRK1),1)
          !! Save
          If (IfChol) Then
            Call CnstAB_SSDM(Work(ipDPTAO),Work(ipWRK1))
          Else
            !! Well, it is not working any more. I need to use
            !! Position='APPEND', but it is not possible if I need to
            !! use molcas_open or molcas_open_ext2
            call abend()
            Call PrgmTranslate('CMOPT2',RealName,lRealName)
C           Open (Unit=LuCMOPT2,
C    *            File=RealName(1:lRealName),
C    *            Position='APPEND',
C    *            Status='OLD',
C    *            Form='UNFORMATTED')
            call molcas_Open(LuCMOPT2,RealName(1:lRealName))
            Do iBasI = 1, nBasT
              Do jBasI = 1, iBasI
                Write (LuCMOPT2) Work(ipDPTAO+iBasI-1+nBasT*(jBasI-1)),
     *                           Work(ipWRK1 +iBasI-1+nBasT*(jBasI-1))
              End Do
            End Do
C
            Close (LuCMOPT2)
          End If
        End If
C       write(6,*) "pt2ao"
C       call sqprt(Work(ipDPTAO),12)
C       call prtril(Work(ipDPT2AO),12)
        CALL GETMEM('DEPSA ','FREE','REAL',ipDEPSA,nAshT)
C
        CALL GETMEM('DPT   ','FREE','REAL',ipDPT   ,nDPTAO)
        CALL GETMEM('DPTC  ','FREE','REAL',ipDPTC  ,nDPTAO)
        CALL GETMEM('DPTAO ','FREE','REAL',ipDPTAO ,nDPTAO)
        CALL GETMEM('DPTCAO','FREE','REAL',ipDPTCAO,nDPTAO)
        CALL GETMEM('FPT   ','FREE','REAL',ipFPT   ,nDPTAO)
        CALL GETMEM('FPTC  ','FREE','REAL',ipFPTC  ,nDPTAO)
        CALL GETMEM('FPTAO ','FREE','REAL',ipFPTAO ,nDPTAO)
        CALL GETMEM('FPTCAO','FREE','REAL',ipFPTCAO,nDPTAO)
        CALL GETMEM('FIFA  ','FREE','REAL',ipFIFA  ,nBsqT)
        CALL GETMEM('FIMO  ','FREE','REAL',ipFIMO  ,nBsqT)
C
C       call test_dens(work(ipolag),work(ipclag),work(iptrf),
C    *                 work(ipwrk1),work(ipwrk2))
C
C       CALL GETMEM('WRK1  ','ALLO','REAL',ipWRK1,nBasT*nBasT)
C       CALL GETMEM('WRK2  ','ALLO','REAL',ipWRK2,nBasT*nBasT)
C       call dcopy_(nbast*nbast,[0.0d+00],0,Work(ipWRK1),1)
C       do i = 1, 5
C         Work(ipWRK1+i-1+nBasT*(i-1)) = 1.0D+00
C       end do
C       do i = 1, 5
C         ii = 5+i
C         do j = 1, 5
C           jj = 5+j
C           Work(ipWRK1+ii-1+nBasT*(jj-1)) = Work(LTORB+25+i-1+5*(j-1))
C         End Do
C       end do
C       do i = 11, 12
C         Work(ipWRK1+i-1+nBasT*(i-1)) = 1.0D+00
C       end do
C       write(6,*) "square transformation matrix"
C       call sqprt(work(ipwrk1),12)
C       write(6,*) "OLag before transformation"
C       call sqprt(work(ipolag),12)
C       Do iSym = 1, nSym
C       write(6,*) "olag before"
C       call sqprt(work(ipolag),nbast)
      if (.false.) then
        do i = 1, nbast
        do j = 1, i
          if (i.eq.j) then
            work(ipwlag+i-1+nbast*(j-1)) = work(ipwlag+i-1+nbast*(j-1))
     *        + 0.5d+00*work(ipolag+i-1+nbast*(j-1))
          else
            work(ipwlag+i-1+nbast*(j-1)) = work(ipwlag+i-1+nbast*(j-1))
     *        + 1.0d+00*work(ipolag+j-1+nbast*(i-1))
          end if
          vali = work(ipolag+i-1+nbast*(j-1))
          valj = work(ipolag+j-1+nbast*(i-1))
          work(ipolag+i-1+nbast*(j-1)) = vali-valj
          work(ipolag+j-1+nbast*(i-1)) = 0.0d+00 ! -vali+valj
          if (i.le.nIsh(1).and.j.le.nIsh(1)) then
          work(ipolag+i-1+nbast*(j-1)) = 0.0d+00
          work(ipolag+j-1+nbast*(i-1)) = 0.0d+00
          end if
          if (i.gt.nIsh(1)+nAsh(1).and.j.gt.nIsh(1)+nAsh(1)) then
          work(ipolag+i-1+nbast*(j-1)) = 0.0d+00
          work(ipolag+j-1+nbast*(i-1)) = 0.0d+00
          end if
        end do
        end do
C       write(6,*) "olag after antisymetrization"
C       call sqprt(work(ipolag),nbast)

C       write(6,*) "Wlag square"
C       call sqprt(work(ipwlag),nbast)
        !! It seems Molcas does similar to GAMESS
        call daxpy_(nbast*nbast,0.5d+00,work(ipolag),1,work(ipwlag),1)
C
        !! W(MO) -> W(AO) using the quasi-canonical orbitals
        !! No need to back transform to natural orbital basis
        Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(LCMOPT2),nBasT,Work(ipWLag),nBasT,
     *              0.0D+00,Work(ipWRK1),nBasT)
        Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(ipWRK1),nBasT,Work(LCMOPT2),nBasT,
     *              0.0D+00,Work(ipWLag),nBasT)
C
        !! square -> triangle for WLag(AO)
        Call DCopy_(nBasT*nBasT,Work(ipWLag),1,Work(ipWRK1),1)
        iBasTr = 1
        iBasSq = 1
        Do iSym = 1, nSym
          nBasI = nBas(iSym)
          liBasTr = iBasTr
          liBasSq = iBasSq
          Do iBasI = 1, nBasI
            Do jBasI = 1, iBasI
              liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
              If (iBasI.eq.jBasI) Then
                Work(ipWLag   +liBasTr-1) = Work(ipWRK1  +liBasSq-1)
              Else
              liBasSq2 = iBasSq + jBasI-1 + nBasI*(iBasI-1)
                Work(ipWLag   +liBasTr-1)
     *            = Work(ipWRK1  +liBasSq-1)
     *            + Work(ipWRK1  +liBasSq2-1)
              End If
              liBasTr = liBasTr + 1
            End Do
          End Do
          iBasTr = iBasTr + nBasI*(nBasI+1)/2
          iBasSq = iBasSq + nBasI*nBasI
        End Do
C
C
C
        !! Transform quasi-canonical -> natural MO basis
        !! orbital Lagrangian
          Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipTrf),nBasT,Work(ipOLag),nBasT,
     *                0.0D+00,Work(ipWRK1),nBasT)
          Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
     *                0.0D+00,Work(ipOLag),nBasT)
      else
        Call DaXpY_(nBasSq,0.5D+00,Work(ipOLag),1,Work(ipWLag),1)
        do i = 1, nbast
        do j = 1, i
C         if (i.gt.nIsh(1).and.i.le.nIsh(1)+nAsh(1).and.
C    *        j.gt.nIsh(1).and.j.le.nIsh(1)+nAsh(1).and.i.ne.j) then
C         work(ipolag+i-1+nbast*(j-1))
C    *      = 2.0d+00*work(ipolag+i-1+nbast*(j-1))
C         work(ipolag+j-1+nbast*(i-1))
C    *      = 2.0d+00*work(ipolag+j-1+nbast*(i-1))
C         end if
C         if (i.eq.j) then
C           work(ipwlag+i-1+nbast*(j-1)) = work(ipwlag+i-1+nbast*(j-1))
C    *        + 0.5d+00*work(ipolag+i-1+nbast*(j-1))
C         else
C           work(ipwlag+i-1+nbast*(j-1)) = work(ipwlag+i-1+nbast*(j-1))
C    *        + 1.0d+00*work(ipolag+j-1+nbast*(i-1))
C         end if
       !  vali = work(ipolag+i-1+nbast*(j-1))
       !  valj = work(ipolag+j-1+nbast*(i-1))
       !  work(ipolag+i-1+nbast*(j-1)) = vali-valj
       !  work(ipolag+j-1+nbast*(i-1)) = 0.0d+00 ! -vali+valj
C         if (i.le.nIsh(1).and.j.le.nIsh(1)) then
C         work(ipolag+i-1+nbast*(j-1)) = 0.0d+00
C         work(ipolag+j-1+nbast*(i-1)) = 0.0d+00
C         end if
C         if (i.gt.nIsh(1)+nAsh(1).and.j.gt.nIsh(1)+nAsh(1)) then
C         work(ipolag+i-1+nbast*(j-1)) = 0.0d+00
C         work(ipolag+j-1+nbast*(i-1)) = 0.0d+00
C         end if
C         if (i.gt.nIsh(1).and.i.le.nIsh(1)+nAsh(1).and.
C    *        j.gt.nIsh(1).and.j.le.nIsh(1)+nAsh(1)) then
C         work(ipolag+i-1+nbast*(j-1)) = 0.0d+00
C         work(ipolag+j-1+nbast*(i-1)) = 0.0d+00
C         end if
        end do
        end do
C       write(6,*) "olag after antisymetrization"
C       call sqprt(work(ipolag),nbast)

C       write(6,*) "Wlag square"
C       call sqprt(work(ipwlag),nbast)
        !! It seems Molcas does similar to GAMESS
C       call daxpy_(nbast*nbast,0.5d+00,work(ipolag),1,work(ipwlag),1)
C
        !! W(MO) -> W(AO) using the quasi-canonical orbitals
        !! No need to back transform to natural orbital basis
        Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(LCMOPT2),nBasT,Work(ipWLag),nBasT,
     *              0.0D+00,Work(ipWRK1),nBasT)
        Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *              1.0D+00,Work(ipWRK1),nBasT,Work(LCMOPT2),nBasT,
     *              0.0D+00,Work(ipWLag),nBasT)
C
        !! square -> triangle for WLag(AO)
        Call DCopy_(nBasT*nBasT,Work(ipWLag),1,Work(ipWRK1),1)
        iBasTr = 1
        iBasSq = 1
        Do iSym = 1, nSym
          nBasI = nBas(iSym)
          liBasTr = iBasTr
          liBasSq = iBasSq
          Do iBasI = 1, nBasI
            Do jBasI = 1, iBasI
              liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
              If (iBasI.eq.jBasI) Then
                Work(ipWLag   +liBasTr-1) = Work(ipWRK1  +liBasSq-1)
              Else
              liBasSq2 = iBasSq + jBasI-1 + nBasI*(iBasI-1)
                Work(ipWLag   +liBasTr-1)
     *            = Work(ipWRK1  +liBasSq-1)
     *            + Work(ipWRK1  +liBasSq2-1)
              End If
              liBasTr = liBasTr + 1
            End Do
          End Do
          iBasTr = iBasTr + nBasI*(nBasI+1)/2
          iBasSq = iBasSq + nBasI*nBasI
        End Do
C
C
C
        !! Transform quasi-canonical -> natural MO basis
        !! orbital Lagrangian
C       write(6,*) "OLag before transformation"
C       call sqprt(work(ipolag),12)
C         do i = 1, nbast
C         do j = 1, i
C           vali = work(ipolag+i-1+nbast*(j-1))
C           valj = work(ipolag+j-1+nbast*(i-1))
C           work(ipolag+i-1+nbast*(j-1)) = vali-valj
C           work(ipolag+j-1+nbast*(i-1)) = 0.0d+00 ! -vali+valj
C         end do
C         end do
C         write(6,*) "OLag after antisymmetrization"
C         call sqprt(work(ipolag),12)
C         call abend
          Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipTrf),nBasT,Work(ipOLag),nBasT,
     *                0.0D+00,Work(ipWRK1),nBasT)
          Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *                1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
     *                0.0D+00,Work(ipOLag),nBasT)
C       write(6,*) "OLag after transformation (pre)"
C       call sqprt(work(ipolag),12)
C       write(6,*) "trf"
C       call sqprt(work(ipTrf),12)
C       If (BSHIFT.ne.0.0D+00) Then
          !! sufficient only for active
          nBasI = nBas(1)
        Call DCopy_(nBasI**2,Work(ipOLag),1,Work(ipWRK1),1)
        Call DGeSub(Work(ipWRK1),nBas(1),'N',
     &              Work(ipWRK1),nBas(1),'T',
     &              Work(ipOLag),nBas(1),
     &              nBas(1),nBas(1))
C         do i = 1, nbast
C         do j = 1, i
C           vali = work(ipolag+i-1+nbast*(j-1))
C           valj = work(ipolag+j-1+nbast*(i-1))
C           work(ipolag+i-1+nbast*(j-1)) = vali-valj
C           work(ipolag+j-1+nbast*(i-1)) = 0.0d+00 ! -vali+valj
C         end do
C         end do
C       write(6,*) "OLag after transformation (pre)"
C       call sqprt(work(ipolag),12)
C       End If
C       End Do
      end if
        !! for configuration Lagrangian
        !! Somehow, transformations using TRACI_RPT2 do not work well
C       ipTrfL = ipTrf+nAshT*nBasT+nAshT
C       Call DGemm_('N','N',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipTrfL),nBasT,Work(ipRDMEIG),nAshT,
C    *              0.0D+00,Work(ipWRK1),nAshT)
C       Call DGemm_('N','T',nAshT,nAshT,nAshT,
C    *              1.0D+00,Work(ipWRK1),nAshT,Work(ipTrfL),nBasT,
C    *              0.0D+00,Work(ipRDMEIG),nAshT)
        !! W (renormalization) term
C         Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipTrf),nBasT,Work(ipWLag),nBasT,
C    *                0.0D+00,Work(ipWRK1),nBasT)
C         Call DGemm_('N','T',nBasT,nBasT,nBasT,
C    *                1.0D+00,Work(ipWRK1),nBasT,Work(ipTrf),nBasT,
C    *                0.0D+00,Work(ipWLag),nBasT)
C       write(6,*) "RDMEIG after transformation"
C       call sqprt(work(iprdmeig),5)
C       CALL GETMEM('WRK2  ','FREE','REAL',ipWRK2,nBasT*nBasT)
C


C       call dcopy_(nbast*nbast,[0.0d+00],0,work(ipwlag),1)
C       write(6,*) "Wlag"
C       call sqprt(work(ipwlag),nbast)
C       call test_dens(work(ipolag),work(ipclag),work(iptrf),
C    *                 work(ipwrk1),work(ipwrk2))
        CALL GETMEM('TRFMAT','FREE','REAL',ipTRF   ,nBsqT)
        CALL GETMEM('WRK1  ','FREE','REAL',ipWRK1,nBasT*nBasT)
        CALL GETMEM('WRK2  ','FREE','REAL',ipWRK2,nBasT*nBasT)
C
        !! Calculate the configuration Lagrangian again.
        !! The contribution comes from the derivative of eigenvalues.
C       do i = 1, 5
C         do j = 1, 5
C           val = work(ipolag+5+i-1+12*(5+j-1))
C           work(iprdmeig+i-1+5*(j-1))
C         = work(iprdmeig+i-1+5*(j-1)) + val
C           work(iprdmeig+j-1+5*(i-1))
C         = work(iprdmeig+j-1+5*(i-1)) + val
C         end do
C       end do
C       write(6,*) "clag in det",nclag
C       do i = 1, nclag
C         write(6,'(i3,1f20.10)') i,work(ipclag+i-1)
C       end do
C       write(6,*) "RDMEIG"
C       call sqprt(work(iprdmeig),5)
        !! Here is the original place
C       Call CLagEig(Work(ipCLag),Work(ipRDMEIG),.true.)
C       write(6,*) "clag in det",nclag
C       do i = 1, nclag
C         write(6,'(i3,1f20.10)') i,work(ipclag+i-1)
C       end do
        !! Here is the original place
C       Call CLagFinal(Work(ipCLag),Work(ipSLag))
C       write(6,*) "clag in det",nclag
C       do i = 1, nclag
C         write(6,'(i3,1f20.10)') i,work(ipclag+i-1)
C       end do

C       CALL GETMEM('aaa','allo','REAL',kdtoc,1000)
C       CALL GETMEM('bbb','allo','REAL',kdtoc2,1000)
C       call dcopy_(ndet,[0.0d+00],0,work(ldpt),1)
C       call csdtvc(work(ipclag),work(ldpt),1,
C    *              work(kdtoc),iwork(kicts(1)),1,0)
C       write(6,*) "clag in det",ndet
C       do i = 1, ndet
C         write(6,'(i3,2f20.10)') i,work(ldpt+i-1),work(ipclag+i-1)
C       end do
C       CALL GETMEM('aaa','free','REAL',kdtoc,1000)
C       CALL GETMEM('bbb','free','REAL',kdtoc2,1000)
C
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA ,nAshT*nAshT)
        CALL GETMEM('RDMEIG','FREE','REAL',ipRDMEIG,nAshT*nAshT)
      END IF

      CALL GETMEM('DPT','FREE','REAL',LDPT,NDPT)
      IDMOFF=0
      IDSOFF=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        DO IP=1,NO
          DO IQ=1,IP
            IDM=IDMOFF+(IP*(IP-1))/2+IQ
            IDSUM=IDSOFF+IP+NO*(IQ-1)
            DMAT(IDM)=DMAT(IDM)+WORK(LDSUM-1+IDSUM)
          END DO
        END DO
        IDMOFF=IDMOFF+(NO*(NO+1))/2
        IDSOFF=IDSOFF+NO**2
      END DO
      CALL GETMEM('DSUM','FREE','REAL',LDSUM,NDPT)
C Scale with 1/DENORM to normalize
      X=1.0D0/DENORM
C     write(6,*) "scaling of DM: ",x
        x = 1.0d+00
      CALL DSCAL_(NDMAT,X,DMAT,1)
C       write(6,*) "final dmat in MO"
C       do i = 1, ndmat
C         write(6,'(i4,f20.10)') i,dmat(i)
C       end do
C       call prtril(dmat,nbast)

CSVC: For true parallel calculations, replicate the DMAT array
C so that the slaves have the same density matrix as the master.
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        IF (.NOT.KING()) THEN
          CALL DCOPY_(NDMAT,[0.0D0],0,DMAT,1)
        END IF
        CALL GADSUM(DMAT,NDMAT)
      END IF
#endif

      CALL QEXIT('DENS')
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      Subroutine CnstTrf(Trf0,Trf)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
C#include "output.fh"
C#include "eqsolv.fh"
C#include "WrkSpc.fh"
C#include "sigma.fh"
C#include "para_info.fh"

C#include "caspt2_grad.fh"
C#include "csfbas.fh"
C
C#include "pt2_guga.fh"
C
      Dimension Trf0(*),Trf(*)
C
      iSQ = 0
      iTOrb = 1 ! LTOrb
      ipTrfL = 0
C     write(6,*) "norbt = ",norbt
C     write(6,*) "nosqt = ", nosqt
C     write(6,*) "nbast = ", nbast
      Do iSym = 1, nSym
        nBasI = nBas(iSym)
        nOrbI = nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nSshI = nSsh(iSym)
        nDelI = nDel(iSym)
        NR1   = nRAS1(iSym)
        NR2   = nRAS2(iSym)
        NR3   = nRAS3(iSym)
        write(6,*) "nBasI",nBas(iSym)
        write(6,*) "nOrbI",nOrb(iSym)
        write(6,*) "nFroI",nFro(iSym)
        write(6,*) "nIshI",nIsh(iSym)
        write(6,*) "nAshI",nAsh(iSym)
        write(6,*) "nSshI",nSsh(iSym)
        write(6,*) "nDelI",nDel(iSym)
        nCor  = nFroI + nIshI
        nOcc  = nCor  + nAshI
        nVir  = nSshI + nDelI
        ipTrfL = ipTrfL + iSQ
        !! frozen + inactive
C       Do iIsh = 1, nFroI + nIshI
C         Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = 1.0D+00
C       End Do
        !! frozen
        Do iIsh = 1, nFroI
          Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = 1.0D+00
        End Do
        !! inactive
        Do I = 1, nIshI
          iIsh = nFroI + I
          Do J = 1, nIshI
            jIsh = nFroI + J
            IJ=I-1+nIshI*(J-1)
            Trf(ipTrfL+iIsh+nBasI*(jIsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nIshI*nIshI
        !! RAS1 space
        Do I = 1, NR1
          iAsh = nCor + I
          Do J = 1, NR1
            jAsh = nCor + J
            IJ=I-1+NR1*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR1*NR1
        !! RAS2 space
        Do I = 1, NR2
          iAsh = nCor + NR1 + I
          Do J = 1, NR2
            jAsh = nCor + NR1 + J
            IJ=I-1+NR2*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR2*NR2
        !! RAS3 space
        Do I = 1, NR3
          iAsh = nCor + NR1 + NR2 + I
          Do J = 1, NR3
            jAsh = nCor + NR1 + NR2 + J
            IJ=I-1+NR3*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR3*NR3
C       call sqprt(trf,12)
      ! !! Active
      ! Do iAsh0 = 1, nAshI
      !   iAsh = nCor + iAsh0
      !   Do jAsh0 = 1, nAshI
      !     jAsh = nCor + jAsh0
C     !     Work(ipTrfL+iAsh-1+nBasI*(jAsh-1))
C    *!       = Work(iTOrb+nIshI*nIshI+iAsh0-1+nAshI*(jAsh0-1))
      !     Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     *!       = Trf0(iTOrb+iAsh0-1+nAshI*(jAsh0-1))
      !   End Do
      ! End Do
C       call sqprt(trf,12)
        !! virtual + deleted (deleted is not needed, though)
C       Do iSsh = nOcc+1, nOcc+nVir
C         Trf(ipTrfL+iSsh+nBasI*(iSsh-1)) = 1.0D+00
C       End Do
        Do I = 1, nVir
          iSsh = nCor + nAshI + I
          Do J = 1, nVir
            jSsh = nCor + nAshI + J
            IJ=I-1+nVir*(J-1)
            Trf(ipTrfL+iSsh+nBasI*(jSsh-1))
     *        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nSshI*nSshI
C       call sqprt(trf,12)
        iSQ = iSQ + nBasI*nBasI

C       n123 = nAshI*nAshI !! just for CAS at present
C       iTOrb = iTOrb + n123 + nSshI*nSshI
C     write(6,*) "transformation matrix"
C     call sqprt(trf,nbasi)
      End Do
C
      Return
C
      End Subroutine CnstTrf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AddDEPSA(DPT2,DEPSA)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
C
      DIMENSION DPT2(*),DEPSA(nAshT,nAshT)
C
C     write(6,*) "DPT2MO"
C     call sqprt(dpt2,nbas(1)-ndel(1))
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2.gt.0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nAsh(iSym)
            iOrb1 = nIsh(iSym)+iOrb0
            iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              jOrb1 = nIsh(iSym)+jOrb0
              jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          = DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          + DEPSA(iOrb0,jOrb0)
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nBas(iSym)-nDel(iSym)
            Do jOrb = 1, iOrb-1
              Val =(DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1))
     *             +DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)))*0.5D+00
              DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1)) = Val
              DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
C     write(6,*) "DPT2MO after DEPSA"
C     call sqprt(dpt2,nbas(1)-ndel(1))
C
      End Subroutine AddDEPSA
C
C-----------------------------------------------------------------------
C
      Subroutine TRAFRO(MODE)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      DIMENSION nFroTmp(8),nOshTmp(8),nOrbTmp(8)
C
      If (Mode.eq.1) Then
        Do jSym = 1, 8
          nFroTmp(jSym) = nFro(jSym)
          nOshTmp(jSym) = nOsh(jSym)
          nOrbTmp(jSym) = nOrb(jSym)
          nOsh(jSym) = nFro(jSym)+nIsh(jSym)+nAsh(jSym)
          nOrb(jSym) = nOsh(jSym)+nSsh(jSym)
          nFro(jSym) = 0
        End Do
      End If
C
      Call GetMem('LCMO','ALLO','REAL',LCMO,NCMO)
      Call DCopy_(NCMO,WORK(LCMOPT2),1,WORK(LCMO),1)
      if (IfChol) then
        call TRACHO3(WORK(LCMO))
      else
        call TRACTL(0)
      end if
      Call GetMem('LCMO','FREE','REAL',LCMO,NCMO)
C
      If (Mode.eq.1) Then
        Do jSym = 1, 8
          nFro(jSym) = nFroTmp(jSym)
          nOsh(jSym) = nOshTmp(jSym)
          nOrb(jSym) = nOrbTmp(jSym)
        End Do
      End If
C
      Return
C
      End Subroutine TRAFRO
C
C-----------------------------------------------------------------------
C
      Subroutine DPT2_TrfStore(Scal,DPT2q,DPT2n,Trf,WRK)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
C
      Dimension DPT2q(*),DPT2n(*),Trf(*),WRK(*)
C
      iMO = 1
      Do iSym = 1, nSym
        If (nOrb(iSym).GT.0) Then
          nOrbI = nBas(iSym)-nDel(iSym)
          !! Quasi-canonical -> natural transformation of DPT2
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,
     *                1.0D+00,Trf(iMO),nOrbI,DPT2q(iMO),nOrbI,
     *                0.0D+00,WRK,nOrbI)
          Call DGemm_('N','T',nOrbI,nOrbI,nOrbI,
     *                   Scal,WRK,nOrbI,Trf(iMO),nOrbI,
     *                1.0D+00,DPT2n(iMO),nOrbI)
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
C
      Return
C
      End Subroutine DPT2_TrfStore
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DPT2_Trf(DPT2,DPT2AO,CMO,DEPSA,DSUM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      DIMENSION DPT2(*),DPT2AO(*),CMO(*),DEPSA(nAshT,nAshT),DSUM(*)
C
      !! DPT2 transformation
      !! Just transform DPT2 (in MO, block-squared) to DPT2AO (in AO,
      !! block-squared). Also, for DPT2C which couples with the inactive
C     !! density matrix.
C     write (6,*) "DPT2_Trf"
C     write(6,*) "dpt2"
C     do isym = 1, nsym
C       nbasi = nbas(isym)
C       write(6,*) "for symmetry :", isym,nbasi
C       call sqprt(dpt2,nbasi)
C     end do
      CALL GETMEM('WRK   ','ALLO','REAL',ipWRK,nBSQT)
C     write(6,*) "vec"
C     call sqprt(cmo,12)
C     write(6,*) "DPT2MO"
C     call sqprt(dpt2,norb(1))
C
      !! MO -> AO back transformation
      iCMO =1
      iAO = 1
      iMO = 1
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
C       iOFF = iWTMP + nBas(iSym)*nBas(iSym)
        If (nORB(ISYM).GT.0) Then
          nBasI = nBas(iSym)
          nOrbI = nOrb(iSym)
          !! Add active orbital density
          Do iOrb0 = 1, nAsh(iSym)
            iOrb = nIsh(iSym)+iOrb0
            iOrb2= nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              jOrb = nIsh(iSym)+jOrb0
              jOrb2= nFro(iSym)+nIsh(iSym)+jOrb0
C             if (iorb0.eq.jorb0) then
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + DEPSA(iOrb0,jOrb0)
              DSUM(iMO+iOrb-1+nOrbI*(jOrb-1))
     *          = DSUM(iMO+iOrb2-1+nOrbI*(jOrb2-1)) + DEPSA(iOrb0,jOrb0)
C             end if
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nOrb(iSym)
            Do jOrb = 1, iOrb
              Val =(DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *             +DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*0.5D+00
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
C             Val =(DSUM(iMO+iOrb-1+nOrbI*(jOrb-1))
C    *             +DSUM(iMO+jOrb-1+nOrbI*(iOrb-1)))*0.5D+00
C             DSUM(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
C             DSUM(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
            End Do
          End Do
          !! First, DPT2 -> DPT2AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
     *                 1.0D+00,CMO(iCMO),nBasI,DPT2(iMO),nBasI,
     *                 0.0D+00,Work(ipWRK),nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
     *                 1.0D+00,Work(ipWRK),nBasI,CMO(iCMO),nBasI,
     *                 0.0D+00,DPT2AO(iAO),nBasI)
          !! Second, DPT2C -> DPT2CAO
C         CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
C    *                 1.0D+00,CMO(iCMO),nBasI,DPT2C(iMO),nBasI,
C    *                 0.0D+00,Work(ipWRK),nBasI)
C         CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
C    *                 1.0D+00,Work(ipWRK),nBasI,CMO(iCMO),nBasI,
C    *                 0.0D+00,DPT2CAO(iAO),nBasI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
        iMO  = iMO  + nBasI*nBasI
      End Do
C     write(6,*) "DPT2MO after DEPSA"
C     call sqprt(dpt2,norb(1))
C
C     write(6,*) "dpt2ao"
C     do isym = 1, nsym
C       nbasi = nbas(isym)
C       write(6,*) "for symmetry :", isym,nbasi
C       call sqprt(dpt2ao,nbasi)
C     end do
C
      CALL GETMEM('WRK   ','FREE','REAL',ipWRK,nBSQT)
C
      END SUBROUTINE DPT2_Trf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer(DPT2,DPT2C,FPT2AO,FPT2CAO,RDMEIG,CMO,Trf,
     *                  FPT2,FPT2C,FIFA,FIMO,RDMSA)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      DIMENSION DPT2(*),DPT2C(*),FPT2AO(*),FPT2CAO(*),RDMEIG(*),CMO(*),
     *          Trf(*)
      DIMENSION FPT2(*),FPT2C(*),FIFA(*),FIMO(*),RDMSA(*)
C
C     write (6,*) "here is eigder"
      CALL GETMEM('WRK1 ','ALLO','REAL',ipWRK1 ,nBSQT)
      CALL GETMEM('WRK2 ','ALLO','REAL',ipWRK2 ,nBSQT)
      CALL GETMEM('FPT2 ','ALLO','REAL',ipFPT2 ,nBSQT)
      CALL GETMEM('FPT2C','ALLO','REAL',ipFPT2C,nBSQT)
C
      !! AO -> MO transformation
      iCMO =1
      iAO = 1
C     iMO = 1
C     write(6,*) "fpt2ao"
C     call sqprt(fpt2ao,nbasT)
C     write(6,*) "fpt2cao"
C     call sqprt(fpt2cao,nbasT)
      if (nfrot.ne.0) then
        Call DCopy_(nBsqT,FPT2,1,Work(ipFPT2),1)
        Call DCopy_(nBsqT,FPT2C,1,Work(ipFPT2C),1)
      else
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
C       iOFF = iWTMP + nBas(iSym)*nBas(iSym)
        If (nOrb(iSym).GT.0) Then
          nBasI = nBas(iSym)
          nOrbI = nOrb(iSym)
          !! First, FPT2(AO) -> FPT2(MO)
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                 1.0D+00,CMO(iCMO),nBasI,FPT2AO(iAO),nBasI,
     *                 0.0D+00,Work(ipWRK1),nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                 1.0D+00,Work(ipWRK1),nOrbI,CMO(iCMO),nBasI,
     *                 0.0D+00,Work(ipFPT2+iAO-1),nOrbI)
          !! Second, FPT2C(AO) -> FPT2C(MO)
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                 1.0D+00,CMO(iCMO),nBasI,FPT2CAO(iAO),nBasI,
     *                 0.0D+00,Work(ipWRK1),nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                 1.0D+00,Work(ipWRK1),nOrbI,CMO(iCMO),nBasI,
     *                 0.0D+00,Work(ipFPT2C+iAO-1),nOrbI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
C       iMO  = iMO  + nBasI*nBasI
      End Do
      end if
C     write(6,*) "fpt2"
C     call sqprt(work(ipfpt2),nbasT)
C     write(6,*) "fpt2c"
C     call sqprt(work(ipfpt2c),nbasT)
C
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2) ,1)
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2C),1)
C     write(6,*) "fpt2mo"
C     call sqprt(work(ipfpt2),nbasT)
C
C     write(6,*) "fpt2 in MO"
C     do isym = 1, nsym
C       nbasi = nbas(isym)
C       write(6,*) "for symmetry :", isym,nbasi
C       call sqprt(work(ipfpt2),nbasi)
C     end do
C
      !! construct Fock in MO
C
C     write(6,*) "ndref = ", ndref
C     write(6,*) "nstate = ", state
C     write(6,*) "ldref"
C     do i = 1, ndref
C       write(6,'(i3,f20.10)') i,Work(ldref+i-1)
C     end do
C     write(6,*) "ldmix"
C     do istate = 1, nstate
C     write(6,*) "istate = ", istate
C     do i = 1, ndref
C       write(6,'(i3,f20.10)') i,Work(ldmix+i-1+ndref*(istate-1))
C     end do
C     end do
C     write(6,*) "average"
C     do i = 1, ndref
C       val = 0.0d+00
C       do istate = 1, nstate
C         val = val + Work(ldmix+i-1+ndref*(istate-1))/dble(nstate)
C       end do
C       write(6,'(i3,f20.10)') i,val
C     end do
      iSQ = 1
C     write(6,*) "olag before"
C     call sqprt(work(ipolag),nbast)
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nSshI = nSsh(iSym)
        nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        nOcc  = nCor  + nAshI
        nVir  = nSshI + nDelI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        CALL DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        !! Active orbital contributions: (p,q) = (all,act)
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA,nAshI*nAshI)
        !  Construct the active density of the orbital energy
        !  Assume the state-averaged density (SS- and XMS-CASPT2)
C       nSeq = 0
C       Call DCopy_(nAshI*nAshI,[0.0D+00],0,Work(ipWRK1),1)
C       Do iState = 1, nState
C         Wgt  = Work(LDWgt+iState-1+nState*(iState-1))
C         Wgt  = 1.0D+00/nState
C         Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
C    *                Work(ipWRK1),1)
C       End Do
        !  RDM of CASSCF
        !  This is likely defined by a set of natural orbitals.
        !  Here, we have to transform to a set of quasi-canonical
        !  orbitals, so forward transformation is appropriate.
C       Call SQUARE(Work(ipWRK1),Work(ipRDMSA),1,nAshI,nAshI)
        Call DCopy_(nAshT**2,RDMSA,1,Work(ipRDMSA),1)
        !! nbast?
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *                      Work(ipRDMSA),nAshT,
     *              0.0D+00,Work(ipWRK1),nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipWRK1),nAshT,
     *                      Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipRDMSA),nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              1.0D+00,Work(ipFPT2+iSQ-1+nOrbI*nCor),nOrbI,
     *                      Work(ipRDMSA),nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA,nAshI*nAshI)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
C       If (nFroI.eq.0) Then
C         Call SQUARE(Work(LFIFA+iSQ-1),Work(ipWRK1),1,nOrbI,nOrbI)
C       Else
C         Call OLagFroSq(iSym,Work(LFIFA+iSQ-1),Work(ipWRK1))
C       End If
C       write(6,*) "fock in MO"
C       call sqprt(FIFA(iSQ),norbi)
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
C    *              2.0D+00,Work(ipWRK1),nOrbI,DPT2(iSQ),nOrbI,
     *              2.0D+00,FIFA(iSQ),nOrbI,DPT2(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
C
        !! explicit derivative of the effective Hamiltonian
        !! dfpq/da = d/da(C_{mu p} C_{nu q} f_{mu nu})
        !!         = f_{mu nu}^a + (C_{mu m} U_{mp} C_{nu q} + C_{mu p} C_{nu m} U_{mq}) * f_{mu nu}
        !!         = f_{mu nu}^a + U_{mp} f_{mq} + U_{mq} f_{pm}
        !! U_{pq}  = f_{pm} df_{qm} + f_{mp} df_{mq}
C       Call SQUARE(Work(LFIMO+iSQ-1),Work(ipWRK1),1,nOrbI,nOrbI)
C       write(6,*) "effective fock in MO"
C       call sqprt(FIMO(iSQ),norbi)
        If (nFroT.eq.0) Then
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,
C    *              1.0D+00,Work(ipWRK1),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
        CALL DGEMM_('T','N',nOrbI,nOrbI,nOrbI,
C    *              1.0D+00,Work(ipWRK1),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,
     *              1.0D+00,Work(ipOLAG),nOrbI)
        End If
        !! Implicit derivative of inactive orbitals (DPT2C)
        Call DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2C+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
C     write(6,*) "olag in eigder"
C     call sqprt(work(ipolag),nbasT)
C
C     ----- CASSCF density derivative contribution in active space
C
      iSQ = 1
      iSQA= 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nCor  = nFroI + nIshI
        Do iT = 1, nAshI
          iTabs = nCor + iT
          Do iU = 1, nAshI
            iUabs = nCor + iU
            iTU = iTabs-1 + nOrbI*(iUabs-1)
            iTUA= iT   -1 + nAshI*(iU   -1)
            RDMEIG(iSQA+iTUA)
     *        = RDMEIG(iSQA+iTUA) + Work(ipFPT2+iSq-1+iTU)
C           write(6,'(2i3,f20.10)') it,iu,Work(ipFPT2+iSq-1+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C     write(6,*) "rdmeig"
C     call sqprt(rdmeig,5)
C
      CALL GETMEM('WRK1 ','FREE','REAL',ipWRK1 ,nBSQT)
      CALL GETMEM('WRK2 ','FREE','REAL',ipWRK2 ,nBSQT)
      CALL GETMEM('FPT2 ','FREE','REAL',ipFPT2 ,nBSQT)
      CALL GETMEM('FPT2C','FREE','REAL',ipFPT2C,nBSQT)
C
      END SUBROUTINE EigDer
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EigDer2(RDMEIG,Trf,FIFA,RDMSA,DEPSA,WRK1,WRK2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      DIMENSION RDMEIG(*),Trf(*)
      DIMENSION FIFA(*),RDMSA(*),DEPSA(*),WRK1(*),WRK2(*)
C
      CALL GETMEM('FPT2 ','ALLO','REAL',ipFPT2 ,nBSQT)
C
      !! Compute G(D), where D=DEPSA
      Call DEPSATrf(DEPSA,Work(ipFPT2),WRK1,WRK2)
      Call DScal_(nBSQT,2.0D+00,Work(ipFPT2),1)
C
      iSQ = 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nSshI = nSsh(iSym)
        nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        nOcc  = nCor  + nAshI
        nVir  = nSshI + nDelI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        CALL DaXpY_(nOrbI*nCor,2.0D+00,Work(ipFPT2+iSQ-1),1,
     *              Work(ipOLAG+iSQ-1),1)
        !! Active orbital contributions: (p,q) = (all,act)
        CALL GETMEM('RDMSA ','ALLO','REAL',ipRDMSA,nAshI*nAshI)
        Call DCopy_(nAshT**2,RDMSA,1,Work(ipRDMSA),1)
        Call DGemm_('T','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *                      Work(ipRDMSA),nAshT,
     *              0.0D+00,WRK1,nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,WRK1,nAshT,
     *                      Trf(iSQ+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipRDMSA),nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,
     *              1.0D+00,Work(ipFPT2+iSQ-1+nOrbI*nCor),nOrbI,
     *                      Work(ipRDMSA),nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
        CALL GETMEM('RDMSA ','FREE','REAL',ipRDMSA,nAshI*nAshI)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
        CALL DGEMM_('N','T',nOrbI,nAshI,nAshI,
     *              2.0D+00,FIFA(1+nOrbI*nCor),nOrbI,DEPSA,nAshI,
     *              1.0D+00,Work(ipOLAG+iSQ-1+nOrbI*nCor),nOrbI)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
C
C     ----- CASSCF density derivative contribution in active space
C
      iSQ = 1
      iSQA= 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nCor  = nFroI + nIshI
        Do iT = 1, nAshI
          iTabs = nCor + iT
          Do iU = 1, nAshI
            iUabs = nCor + iU
            iTU = iTabs-1 + nOrbI*(iUabs-1)
            iTUA= iT   -1 + nAshI*(iU   -1)
            RDMEIG(iSQA+iTUA) = Work(ipFPT2+iSq-1+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do
C
      CALL GETMEM('FPT2 ','FREE','REAL',ipFPT2 ,nBSQT)
C
      END SUBROUTINE EigDer2
C
C-----------------------------------------------------------------------
C
      Subroutine DEPSATrf(DEPSA,FPT2,WRK1,WRK2)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C
      Dimension DEPSA(nAshT,nAshT),FPT2(*),WRK1(*),WRk2(*)
C
      Call DCopy_(nBasT**2,[0.0D+00],0,FPT2,1)
C
      iSym = 1
      iSymA= 1
      iSymI= 1
      iSymB= 1
      iSymJ= 1
C
      If (nFroT.ne.0.and.IfChol) Then
        !! DEPSA(MO) -> DEPSA(AO) -> G(D) in AO -> G(D) in MO
        !! The Cholesky vectors do not contain frozen orbitals...
        Call GetMem('DAO ','ALLO','REAL',ipDAO ,nBsqT)
        Call GetMem('DMO ','ALLO','REAL',ipDMO ,nBsqT)
        Call GetMem('WRK1','ALLO','REAL',ipWRK1,nBsqT)
        Call GetMem('WRK2','ALLO','REAL',ipWRK2,nBsqT)
        !! First, MO-> AO transformation of DEPSA
        Do iSym = 1, nSym
          Call DCopy_(nBsqT,[0.0D+00],0,Work(ipDMO),1)
          nCorI = nFro(iSym)+nIsh(iSym)
          nBasI = nBas(iSym)
          Do iAsh = 1, nAsh(iSym)
            Do jAsh = 1, nAsh(iSym)
              Work(ipDMO+nCorI+iAsh-1+nBasI*(nCorI+jAsh-1))
     *          = DEPSA(iAsh,jAsh)
            End Do
          End Do
          Call OLagTrf(1,iSym,Work(LCMOPT2),Work(ipDMO),
     *                 Work(ipDAO),Work(ipWRK1))
        End Do
        !! Compute G(D)
        Call DCopy_(nBsqT,[0.0D+00],0,Work(ipWRK1),1)
        Call DCopy_(nBsqT,[0.0D+00],0,Work(ipDMO),1)
        !! it's very inefficient
        Call OLagFro4(1,1,1,1,1,
     *                Work(ipDAO),Work(ipWRK1),Work(ipDMO),
     *                Work(ipWRK1),Work(ipWRK2))
        !! G(D) in AO -> G(D) in MO
        Do iSym = 1, nSym
          Call OLagTrf(2,iSym,Work(LCMOPT2),FPT2,
     *                 Work(ipDMO),Work(ipWRK1))
        End Do
        Call GetMem('DAO ','FREE','REAL',ipDAO ,nBsqT)
        Call GetMem('DMO ','FREE','REAL',ipDMO ,nBsqT)
        Call GetMem('WRK1','FREE','REAL',ipWRK1,nBsqT)
        Call GetMem('WRK2','FREE','REAL',ipWRK2,nBsqT)
      Else
        nCorI = nFro(iSym)+nIsh(iSym)
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            Call DaXpY_(nBasT**2,DEPSA(iAshI,jAshI),WRK1,1,FPT2,1)
C
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            Call DaXpY_(nBasT**2,-0.5D+00*DEPSA(iAshI,jAshI),
     *                  WRK1,1,FPT2,1)
          End Do
        End Do
      End If
C
      Return
C
      End Subroutine DEPSATrf
C
C-----------------------------------------------------------------------
C
C*MODULE MTHLIB  *DECK PRTRIL
      SUBROUTINE PRTRIL(D,N)
C
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION D(*)
C
      MAX = 5
      MM1 = MAX - 1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(6,9008)
         WRITE(6,9028) (I,I=I0,IL)
         WRITE(6,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            WRITE(6,9048) I,'        ',(D(J),J=J0,JL)
  100    CONTINUE
  120 CONTINUE
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,A8,10F11.6)
      END
C
C-----------------------------------------------------------------------
C
      Subroutine CnstAB_SSDM(DPT2AO,SSDM)
C
      USE CHOVEC_IO
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
#include "warnings.fh"
#include "chocaspt2.fh"
#include "choptr.fh"
#include "choglob.fh"
C
      Dimension DPT2AO(*),SSDM(*)
      Character*4096 RealName
      Integer iSkip(8),ipWRK(8)
      integer nnbstr(8,3)
      Logical is_error
C
      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      call getritrfinfo(nnbstr,maxvec,n2)
      iSym = 1 !! iSym0
      nVec = NVLOC_CHOBATCH(1)
C
      NumChoTot = 0
      Do jSym = 1, nSym
        NumChoTot = NumChoTot + NumCho_PT2(jSym)
      End Do
      NumCho=NumChoTot
      Do jSym = 1, nSym
        iSkip(jSym) = 1
      End Do
C
      nBasI  = nBas(iSym)
C
      Call GetMem('A_PT2 ','ALLO','REAL',ipA_PT2,NumChoTot**2)
      !! Read A_PT2
      Call PrgmTranslate('CMOPT2',RealName,lRealName)
C     Open (Unit=LuCMOPT2,
C    *      File=RealName(1:lRealName),
C    *      Status='OLD',
C    *      Form='UNFORMATTED')
C     call molcas_Open(LuCMOPT2,RealName(1:lRealName))
      Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
     &                      'DIRECT','UNFORMATTED',
     &                      iost,.FALSE.,
     &                        1,'REPLACE',is_error)
      Do i = 1, NumChoTot*NumChoTot
        Read (LuCMOPT2) Work(ipA_PT2+i-1)
      End Do
C
      CALL GETMEM('CHSPC','ALLO','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTVEC','ALLO','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','ALLO','REAL',ipWRK(iSym),nBasT*nBasT)
      !! V(P) = (mu nu|P)*D_{mu nu}
      CALL GETMEM('V1   ','ALLO','REAL',ipV1,NumCho)
      CALL GETMEM('V2   ','ALLO','REAL',ipV2,NumCho)
      !! B_SSDM(mu,nu,P) = D_{mu rho}*D_{nu sigma}*(rho sigma|P)
      Call GetMem('B_SSDM','ALLO','REAL',ipB_SSDM,nBasT**2*NumChoTot)
C
      !! Prepare density matrix
      !! subtract the state-averaged density matrix
C
      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym).EQ.0) Return

      ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      JRED1=iWork(ipnt)
      JRED2=iWork(ipnt-1+NumCho_PT2(iSym))

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED.EQ.0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
* For a reduced set, the structure is known, including
* the mapping between reduced index and basis set pairs.
* The reduced set is divided into suitable batches.
* First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        JEND=JSTART+NVECS_RED-1

* Determine batch length for this reduced set.
* Make sure to use the same formula as in the creation of disk
* address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

* Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
C         write(6,*) "ibatch,nbatch = ", ibatch,nbatch
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1

          JREDC=JRED
* Read a batch of reduced vectors
          CALL CHO_VECRD(WORK(IP_CHSPC),NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)
          IF(NUMV.ne.JNUM) THEN
            write(6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(6,*)' read JNUM vectors. Instead it returned NUMV'
            write(6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC.NE.JRED) THEN
            write(6,*)' Rats! It was assumed that the Cholesky vectors'
            write(6,*)' in HALFTRNSF all belonged to a given reduced'
            write(6,*)' set, but they don''t!'
            write(6,*)' JRED, JREDC:',JRED,JREDC
            write(6,*)' Back to the drawing board?'
            write(6,*)' Let the program continue and see what happens.'
          END IF
C
          ipVecL = ip_CHSPC
          Do iVec = JV1, JV2
C
            !! reduced form -> squared AO vector (mu nu|iVec)
            jVref = 1 !! only for iSwap=1
            lscr  = nBasI*(nBasI+1)/2
            If (l_NDIMRS.LT.1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
            End If
            JVEC1 = 1
            JNUM  = 1
            NUMV  = 1
            iSwap = 2
            Call DCopy_(nBasI**2,[0.0D+00],0,Work(ipWRK(iSym)),1)
            Call Cho_ReOrdr(irc,Work(ipVecL),lscr,jVref,
     *                      JVEC1,JNUM,NUMV,iSym,JREDC,iSwap,ipWRK,
     *                      iSkip)
            ipVecL = ipVecL + lscr
C
            Work(ipV1+iVec-1) = DDot_(nBasI**2,DPT2AO,1,
     *                                         Work(ipWRK(iSym)),1)
            Work(ipV2+iVec-1) = DDot_(nBasI**2,SSDM  ,1,
     *                                         Work(ipWRK(iSym)),1)
C
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  1.0D+00,DPT2AO,nBasI,Work(ipWRK(iSym)),nBasI,
     *                  0.0D+00,Work(ipHTVec),nBasI)
            Call DGemm_('N','N',nBasI,nBasI,nBasI,
     *                  1.0D+00,Work(ipHTVec),nBasI,SSDM,nBasI,
     *                  0.0D+00,Work(ipB_SSDM+nBasT**2*(iVec-1)),nBasI)
            do i = 1, nBasT
              do j = 1, i-1
                Val = (Work(ipB_SSDM+i-1+nBasT*(j-1)+nBasT**2*(iVec-1))
     *                +Work(ipB_SSDM+j-1+nBasT*(i-1)+nBasT**2*(iVec-1)))
     *                *0.5d+00
                Work(ipB_SSDM+i-1+nBasT*(j-1)+nBasT**2*(iVec-1)) = Val
                Work(ipB_SSDM+j-1+nBasT*(i-1)+nBasT**2*(iVec-1)) = Val
              end do
            end do
          End Do
C
          ipVecL = ip_CHSPC
          Do iVec = JV1, JV2
C
            !! reduced form -> squared AO vector (mu nu|iVec)
            jVref = 1 !! only for iSwap=1
            lscr  = nBasI*(nBasI+1)/2
            If (l_NDIMRS.LT.1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
            End If
            JVEC1 = 1
            JNUM  = 1
            NUMV  = 1
            iSwap = 2
            Call DCopy_(nBasI**2,[0.0D+00],0,Work(ipWRK(iSym)),1)
            Call Cho_ReOrdr(irc,Work(ipVecL),lscr,jVref,
     *                      JVEC1,JNUM,NUMV,iSym,JREDC,iSwap,ipWRK,
     *                      iSkip)
            ipVecL = ipVecL + lscr
C
            !! Exchange part of A_PT2
            Do jVec = 1, NumCho
              Work(ipA_PT2+iVec-1+NumCho*(jVec-1))
     *          = Work(ipA_PT2+iVec-1+NumCho*(jVec-1))
     *          - DDot_(nBasT**2,Work(ipWRK(iSym)),1,
     *                  Work(ipB_SSDM+nBasT**2*(jVec-1)),1)
            End Do
          End Do
        End Do
      End Do
C
      !! Coulomb
      Call DGEMM_('N','T',NumCho,NumCho,1,
     *            2.0D+00,Work(ipV1),NumCho,Work(ipV2),NumCho,
     *            1.0D+00,Work(ipA_PT2),NumCho)
      Do i = 1, NumCho
        Do j = 1, i-1
          Val = (Work(ipA_PT2+i-1+NumCho*(j-1))
     *          +Work(ipA_PT2+j-1+NumCho*(i-1)))*0.5d+00
          Work(ipA_PT2+i-1+NumCho*(j-1)) = Val
          Work(ipA_PT2+j-1+NumCho*(i-1)) = Val
        End Do
      End Do
      !! Write A_PT2
      REWIND LuCMOPT2
      Do i = 1, NumCho*NumCho
        Write (LuCMOPT2) Work(ipA_PT2+i-1)
      End Do
      Close (LuCMOPT2)
      Call GetMem('A_PT2 ','FREE','REAL',ipA_PT2,NumChoTot**2)
C
C     Call GetMem('B_PT2 ','ALLO','REAL',ipB_PT2,nBasT**2*NumChoTot)
      !! Read B_PT2
      Call PrgmTranslate('GAMMA',RealName,lRealName)
C     Open (Unit=LuGamma,
C    *      File=RealName(1:lRealName),
C    *      Status='OLD',
C    *      Form='UNFORMATTED',
C    *      Access='DIRECT',
C    *      Recl=nBas(iSym)*nBas(iSym)*8)
C     call molcas_Open(LuGamma,RealName(1:lRealName))
      Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                      'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      nBas(iSym)**2*8,'REPLACE',is_error)
      Do iVec = 1, NumCho
        Read  (Unit=LuGAMMA,Rec=iVec)
     *    (Work(ipWRK(iSym)+i-1),i=1,nBasT**2)
        !! The contributions are doubled, because halved in PGet1_RI3?
        !! Coulomb
        Call DaXpY_(nBasT**2,Work(ipV2+iVec-1)*2.D+00,
     *              DPT2AO,1,Work(ipWRK(iSym)),1)
        Call DaXpY_(nBasT**2,Work(ipV1+iVec-1)*2.D+00,
     *              SSDM  ,1,Work(ipWRK(iSym)),1)
        !! Exchange
        Call DaXpY_(nBasT**2,-2.0D+00,
     *              Work(ipB_SSDM+nBasT**2*(iVec-1)),1,
     *              Work(ipWRK(iSym)),1)
        Write (Unit=LuGAMMA,Rec=iVec)
     *    (Work(ipWRK(iSym)+i-1),i=1,nBasT**2)
      End Do
      Close (LuGamma)
C     Call GetMem('B_PT2 ','FREE','REAL',ipB_PT2,nBasT**2*NumChoTot)
C
      CALL GETMEM('CHSPC','FREE','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTVEC','FREE','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','FREE','REAL',ipWRK(iSym),nBasT*nBasT)
      CALL GETMEM('V1   ','FREE','REAL',ipV1,NCHSPC)
      CALL GETMEM('V2   ','FREE','REAL',ipV2,NCHSPC)
      Call GetMem('B_SSDM','FREE','REAL',ipB_SSDM,nBasT**2*NumChoTot)
C
      End Subroutine CnstAB_SSDM
