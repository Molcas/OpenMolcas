!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      Subroutine DEPSAOffC(NCONF,NSTATE,NASHT,NBAST,CLag,DEPSA,FIFA,    &
     &                     FIMO,WRK1,WRK2,U0)

      use Symmetry_Info, only: Mul
      use caspt2_global, only:IPrGlb
      use PrintLevel, only: VERBOSE
      use caspt2_global, only: ConvInvar,SLag
      use sguga, only: SGS, CIS
      use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp, wp, u6
      use Constants, only: Zero, One, Half
      use caspt2_module, only: IFXMS, IFRMS, NSYM, STSYM, NFRO, NISH,   &
     &                         NASH, NORB, NBAS, ISCF, NBTCH, NBTCHES,  &
     &                         NROOTS
!     use caspt2_module, only: NSSH
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      integer(kind=iwp), intent(in) :: NCONF, NSTATE, NASHT, NBAST
      real(kind=wp), intent(inout) :: CLag(nConf,nState),               &
     &  DEPSA(nAshT,nAshT),WRK1(nBasT,nBasT),WRK2(nBasT**2)
      real(kind=wp), intent(in) :: FIFA(nBasT**2), FIMO(nBasT**2),      &
     &  U0(nState,nState)
      real(kind=wp),allocatable :: VecST(:,:),VecS1(:,:),VecS2(:,:),    &
     &                             VecCID(:,:),VecPre(:),VecFancy(:),   &
     &                             VecCIT(:,:),INT1(:),INT2(:),G2(:)
      real(kind=wp), allocatable ::  Eact(:)

      real(kind=wp) :: Thres, DeltaC, Delta, Delta0, AlphaC, Alpha,     &
     &                 ResCI, Beta, Res
      real(kind=wp), external :: DDot_
      integer(kind=iwp) :: nLev, nMidV, iSym, ID, iState, isyci,        &
     &                     MaxIter, Iter

      real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, TIOE, TIOTF0
      real(kind=wp) :: CPTF1, CPTF2, TIOTF1, TIOTF2

      nLev = SGS%nLev
      nMidV= CIS%nMidV

      Thres = ConvInvar !! 1.0d-07

      If (IPRGLB >= VERBOSE) Then
        Write (u6,*)
        Write (u6,'(3X,"Linear Equation for Non-Invariant CASPT2",      &
     &                 " (threshold =",ES9.2,")")') Thres
        Write (u6,*)
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      End If
!
!     If CASPT2 energy is not invariant with respect to rotations within
!     active space (with IPEA shift and/or with RAS reference), the
!     active density obtained in constructing CI derivative is no longer
!     correct... well, it may be correct, but orbital rotations in the
!     active space cannot be parametrized in Z-vector, so analytic
!     derivatives cannot be computed with the existing module. So,
!     the active density is computed in a differnt way.
!
!     See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
!     Section II C 4 "Non-invariance with respect to active MOs"
!     To be more specific, this subroutine solves the linear equation
!     (Eq. (71)) and computes the second term in Eq. (70) later.
!     CLag corresponds to the RHS in Eq. (71).
!
      !! Some post-processing of CI derivative
      !! Somehow, this has to be done in the XMS basis
      Call CLagFinalOffC(nState,SLag)
!
!     ----- Solve the linear equation -----
!     A_{IS,JR}*X_{JR} = CLag_{IS}, where A_{IS,JR} is the CI-CI Hessian
!     which may be seen in Z-vector
!
      call mma_allocate(VecST,nConf,nState,Label='VecST')
      call mma_allocate(VecS1,nConf,nState,Label='VecS1')
      call mma_allocate(VecS2,nConf,nState,Label='VecS2')
      call mma_allocate(VecCID,nConf,nState,Label='VecCID')
!     call mma_allocate(VecS,nState*(nState-1)/2,Label='VecS')
      call mma_allocate(VecPre,nConf,Label='VecPre')
      call mma_allocate(VecFancy,nState**3,Label='VecFancy')

      call mma_allocate(VecCIT,nConf,nState,Label='VecCIT')
      call mma_allocate(INT1,nAshT**2,Label='INT1')
      call mma_allocate(INT2,nAshT**4,Label='INT2')

      call mma_allocate(Eact,nState,Label='Eact')

      !! We do not have Cholesky vectors for frozen orbitals,
      !! so may be it is not possible to get inactive energies?
      !! It can be computed with TimesE2
      iSym = 1
      Call CnstInt(0,nAshT,INT1,INT2)
      Do iState = 1, nState
        ID = IDTCEX(iState)
        If (ISCF == 0) Then
          !! quasi-canonical, XMS
          Call DDaFile(LUCIEX,2,VecCIT(1,iState),nConf,ID)
        Else
          VecCIT(1,iState) = One
        End If
        !! The second term should be removed
        Eact(iState) = Zero
      End Do
      if (ifxms .or. ifrms) then
        !! Transform the CLag and CI vector from XMS to SCF basis
        !! Maybe, in order to define Eact
        Call DGEMM_('N','T',nConf,nState,nState,                        &
     &              One,CLag,nConf,U0,nState,                           &
     &              Zero,VecST,nConf)
        CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)
        Call DGEMM_('N','T',nConf,nState,nState,                        &
     &              One,VecCIT,nConf,U0,nState,                         &
     &              Zero,VecST,nConf)
        VecCIT(1:nConf,1:nState) = VecST(1:nConf,1:nState)
      end if
      Call TimesE2(0,nConf,nState,nAshT,VecCIT,VecS1,INT1,INT2)
      Do iState = 1, nState
        !! scaling with nState is due to the division in TimesE2
        Eact(iState) = -Half*nState*                                    &
     &    DDot_(nConf,VecS1(1,iState),1,VecCIT(1,iState),1)
      End Do
      isyci = 1

      !! Precondition
      Call CnstInt(2,nAshT,INT1,INT2)
      Call CnstPrec(ISYCI,nConf,nRoots,NLEV,nMidV,VecPre,VecCIT,        &
     &              INT1,INT2,VecFancy)
      Call CnstInt(0,nAshT,INT1,INT2)

      !! Begin!
      VecST(1:nConf,1:nState) = CLag(1:nConf,1:nState)

      !! z0 = M^{-1}*r0
      VecS2(1:nConf,1:nState) = VecST(1:nConf,1:nState)
      Call DoPrec(nConf,nRoots,VecST,VecS2,VecS1,VecPre,VecFancy)
      !! p0 = z0
      VecCId(1:nConf,1:nState) = VecS2(1:nConf,1:nState)
      MaxIter = 100
      Iter    = 1
      iSym    = 1
      ! jspin   = 0
      ! r^T dot z
      ! r (residue) = ipST
      ! z (prec. r) = ipS2
      ! p (...)     = ipCId
      ! x (solution)= ipCIT
      ! Ap          = ipS1
      ! r_{k}z_{k}  = ipST*ipS2 = deltaC
      DeltaC = DDot_(nConf*nState,VecST,1,VecS2,1)
      Delta  = DeltaC
      Delta0 = Delta
!
      If (IPRGLB >= VERBOSE) Write(u6,*)                                &
     &      ' Iteration       Delta           Res(CI)        '//        &
     &      '  DeltaC'
      VecCIT(1:nConf,1:nState) = Zero
      If (Delta0 > Abs(Thres)) then
        Do Iter = 1, MaxIter
          If (nConf == 1) Then
            Do iState = 1, nState
              VecCIT(1,iState) = One
            End Do
            Exit
          End If
          !! Compute Ap
          !! ipS2 is used as a workind array
          Call TimesE2(1,nConf,nState,nAshT,VecCId,VecS1,INT1,INT2)

          !! AlphaC = p^T*A*p
          AlphaC= DDot_(nConf*nState,VecS1,1,VecCId,1)
          !! Alpha = r^T*z / AlphaC
          Alpha = Delta/(AlphaC)
          ! new x of CI
          VecCIT(1:nConf,1:nState) = VecCIT(1:nConf,1:nState)           &
     &      + Alpha*VecCId(1:nConf,1:nState)
          ! new r of CI
          VecST(1:nConf,1:nState) = VecST(1:nConf,1:nState)             &
     &      - Alpha*VecS1(1:nConf,1:nState)
          ResCI = sqrt(DDot_(nConf*nState,VecST,1,VecST,1))
          !! z = M^{-1}*r
          VecS2(1:nConf,1:nState) = VecST(1:nConf,1:nState)
          Call DoPrec(nConf,nRoots,VecST,VecS2,VecS1,VecPre,VecFancy)

          !! Append new vectors
          DeltaC= Ddot_(nConf*nState,VecST,1,VecS2,1)
          Beta  = DeltaC/Delta
          Delta = DeltaC
          VecCId(1:nConf,1:nState)                                      &
     &      = Beta*VecCId(1:nConf,1:nState) + VecS2(1:nConf,1:nState)

          If (IPRGLB >= VERBOSE)                                        &
     &    Write(u6,'(I7,4X,ES17.9,ES17.9,ES17.9)')                      &
     &           iter,delta/delta0,resci,deltac

          Res = ResCI
          If (Res <= Abs(Thres)) Exit
        End Do

        If (Iter == MaxIter+1) Then
          write(u6,*)                                                   &
     &    'CI iteration for non-invariant CASPT2 did not converge...'
          call abend()
        End If
      end if

      If (IPRGLB >= VERBOSE) Then
        CALL TIMING(CPTF1,CPE,TIOTF1,TIOE)
        CPUT =CPTF1-CPTF0
        WALLT=TIOTF1-TIOTF0
        Write (u6,*)
        Write (u6,'(3X,"Linear equation converged in ",I3," steps")')   &
     &         iter-1
        Write (u6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
        Write (u6,*)
      End If

      If (IFXMS .OR. IFRMS) Then
        !! Transform back the CLag from CAS to XMS
        Call DGEMM_('N','N',NConf,nState,nState,                        &
     &              One,CLag,nConf,U0,nState,                           &
     &              Zero,VecST,nConf)
        CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)
      End If

      call mma_deallocate(VecS1)
      call mma_deallocate(VecS2)
      call mma_deallocate(VecCId)
!     call mma_deallocate(VecS)
      call mma_deallocate(VecPre)
      call mma_deallocate(VecFancy)
!
!     ----- Construct (a part of) the true active density -----
!     Compute the second term in Eq. (70) = Eq. (72)
!     The SCF, not XMS, basis is used
!
      Do iState = 1, nState
        ID = IDCIEX(iState) !! idtcex?
        If (ISCF == 0) Then
          If (IFXMS .OR. IFRMS) THen
            !! Use unrotated (SCF) CI vector
           Call LoadCI_XMS('C',1,nConf,nState,VecST(1,iState),iState,U0)
          Else
            Call DDaFile(LUCIEX,2,VecST(1,iState),nConf,ID)
          End If
        Else
          VecST(1,iState) = One
        End If
      End Do
      call mma_allocate(G2,nAshT**4,Label='G2')
      Call CnstInt(1,nAshT,INT1,INT2)
      Call CnstDEPSA(nConf,nState,nAshT,VecST,VecCIT,INT1,G2,INT2)
      call mma_deallocate(G2)

      If (IPRGLB >= VERBOSE) Then
        CALL TIMING(CPTF2,CPE,TIOTF2,TIOE)
        CPUT =CPTF2-CPTF1
        WALLT=TIOTF2-TIOTF1
        Write (u6,'(3X,"Off-diagonal density is constructed")')
        Write (u6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
        Write (u6,*)
      End If

      call mma_deallocate(VecST)
      call mma_deallocate(VecCIT)
      call mma_deallocate(INT1)
      call mma_deallocate(INT2)
      call mma_deallocate(Eact)

      Contains
!
!-----------------------------------------------------------------------
!
      Subroutine CLagFinalOffC(nState,SLag)

      use caspt2_module, only: REFENE

      implicit none

      integer(kind=iwp), intent(in) :: nState
      real(kind=wp), intent(inout) :: SLag(nState**2)

      real(kind=wp), allocatable :: CI1(:), CI2(:)

      integer(kind=iwp) :: ijst, ilStat, jlStat
      real(kind=wp) :: Scal, Ovl
      real(kind=wp), external :: DDot_
!
!     Orthogonalize the partial derivative with respect to CI coeff
!
      call mma_allocate(VecST,nConf,nState,Label='VecST')
      Call DGEMM_('N','T',nConf,nState,nState,                          &
     &            One,CLag,nConf,U0,nState,                             &
     &            Zero,VecST,nConf)
      CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)

      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')

      !! Construct SLag
      ijst = 0
      do ilStat = 1, nState
        If (ISCF == 0) Then
          Call LoadCI_XMS('C',1,nConf,nState,CI1,ilStat,U0)
        Else
          CI1(1) = One
        End If
        Do jlStat = 1, ilStat !! -1
          ijst = ilStat + nState*(jlStat-1)
          If (ilStat == jlStat) Cycle
          If (ISCF == 0) Then
            Call LoadCI_XMS('C',1,nConf,nState,CI2,jlStat,U0)
          Else
            CI2(1) = One
          End If
          Scal = DDOT_(nConf,CI1,1,CLag(1,jlStat),1)                    &
     &         - DDOT_(nConf,CI2,1,CLag(1,ilStat),1)
          Scal = Scal/(REFENE(jlStat)-REFENE(ilStat))
          SLag(ijst) = SLag(ijst) + Scal
          IF (IPRGLB >= VERBOSE) THEN
            write(u6,'(1x,"SLag for State ",i1,"-",i1," = ",f20.10)')   &
     &         ilstat,jlstat,slag(ijst)
          END IF
        end do
      end do

      !! Projection
      Do ilStat = 1, nState
        CI1(1:nConf) = CLag(1:nConf,ilStat)
        Do jlStat = 1, nState
          If (ISCF == 0) Then
            Call LoadCI_XMS('C',1,nConf,nState,CI2,jlStat,U0)
          Else
            CI2(1) = One
          End If
          Ovl = DDot_(nConf,CI1,1,CI2,1)
          CLag(1:nConf,ilStat) = CLag(1:nConf,ilStat) - Ovl*CI2(1:nConf)
        End Do
      End Do

      Call DGEMM_('N','N',nConf,nState,nState,                          &
     &            One,CLag,nConf,U0,nState,                             &
     &            Zero,VecST,nConf)
      CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)

      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
      call mma_deallocate(VecST)

      End Subroutine CLagFinalOffC
!
!-----------------------------------------------------------------------
!
      Subroutine CnstInt(Mode,nAshT,INT1,INT2)

      Use CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_module, only: IfChol, NFRO, NBAS

      implicit none

      integer(kind=iwp), intent(in) :: Mode, nAshT
      real(kind=wp), intent(inout) :: INT1(nAshT,nAshT),                &
     &                                INT2(nAshT,nAshT,nAshT,nAshT)

      integer(kind=iwp),allocatable :: BGRP(:,:)
      real(kind=wp),allocatable :: KET(:)

      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) :: nFroI, nIshI, nCorI, nBasI, iAshI,           &
     &  jAshI, iSymA, iSymI, iSymB, iSymJ, JSYM, IB, IB1, IB2, MXBGRP,  &
     &  IBGRP, NBGRP, NCHOBUF, MXPIQK, NADDBUF, IBSTA, IBEND, NV, nKET, &
     &  kAshI, lAshI, iT, iU, iTU, iV, iX, iVX, iOrb, jOrb
!     integer(kind=iwp) :: nSh(8,3)
      real(kind=wp) :: Val

      INT1(:,:) = Zero
      Int2(:,:,:,:) = Zero

      nFroI = nFro(iSym)
      nIshI = nIsh(iSym)
      nCorI = nFroI+nIshI
      nBasI = nBas(iSym)
!
!     --- One-Electron Integral
!
      !! Read H_{\mu \nu}
!     IRC=-1
!     IOPT=6
!     ICOMP=1
!     ISYLBL=1
!     CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,WRK2,ISYLBL)
!     !! triangular -> square transformation
!     Call Square(WRK2,WRK1,1,nBasT,nBasT)
!     !! AO -> MO transformation
!     Call DGemm_('T','N',nBasT,nBasT,nBasT,
!    *            One,CMOPT2,nBasT,WRK1,nBasT,
!    *            Zero,WRK2,nBasT)
!     Call DGemm_('N','N',nBasT,nBasT,nBasT,
!    *            One,WRK2,nBasT,CMOPT2,nBasT,
!    *            Zero,WRK1,nBasT)
      !! Inactive energy
!     Do iCorI = 1, nFro(iSym)+nIsh(iSym)
!       RIn_Ene = RIn_Ene + Two*WRK1(iCorI,iCorI)
!     End Do
      !! Put in INT1
!     Do iAshI = 1, nAsh(iSym)
!       Do jAshI = 1, nAsh(iSym)
!         Val = WRK1(nCorI+iAshI,nCorI+jAshI)
!         INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
!       End Do
!     End Do
      Do iAshI = 1, nAsh(iSym)
        Do jAshI = 1, nAsh(iSym)
          Val = FIMO(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
          INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
        End Do
      End Do
!
!     --- Two-Electron Integral
!
      iSymA = 1
      iSymI = 1
      iSymB = 1
      iSymJ = 1
!     If (.not.IfChol) Then
!       Do iCorI = 1, nFro(iSym)+nIsh(iSym)
!         iOrb = iCorI
!         jOrb = iCorI
!         Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
!         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
!           RIn_Ene = RIn_Ene + Two*WRK1(jCorI,jCorI)
!         End Do
!         Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
!         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
!           RIn_Ene = RIn_Ene - WRK1(jCorI,jCorI)
!         End Do
!       End Do
!     End If
!
      If (IfChol) Then
!       nSh(1:nSym,Inactive) = NISH(1:nSym)
!       nSh(1:nSym,Active  ) = NASH(1:nSym)
!       nSh(1:nSym,Virtual ) = NSSH(1:nSym)
        DO JSYM=1,NSYM
          IB1=NBTCHES(JSYM)+1
          IB2=NBTCHES(JSYM)+NBTCH(JSYM)

          MXBGRP=IB2-IB1+1
          IF (MXBGRP <= 0) CYCLE
          call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
          IBGRP=1
          DO IB=IB1,IB2
           BGRP(1,IBGRP)=IB
           BGRP(2,IBGRP)=IB
           IBGRP=IBGRP+1
          END DO
          NBGRP=MXBGRP

          CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,                         &
     &                         NCHOBUF,MXPIQK,NADDBUF)
          call mma_allocate(KET,NCHOBUF,Label='KETBUF')
!         write(u6,*) 'nchobuf= ', nchobuf
!         write(u6,*) 'nbgrp= ', nbgrp
!         write(u6,*) 'nbtch= ', nbtch(jsym)
          Do IBGRP=1,NBGRP

            IBSTA=BGRP(1,IBGRP)
            IBEND=BGRP(2,IBGRP)
!           write(u6,*) ibsta,ibend

            NV=0
            DO IB=IBSTA,IBEND
              NV=NV+NVLOC_CHOBATCH(IB)
            END DO

            !! int2(tuvx) = (tu|vx)/2
            !! This can be computed without frozen orbitals
            Call Get_Cholesky_Vectors(Active,Active,JSYM,               &
     &                                KET,SIZE(KET),nKet,               &
     &                                IBSTA,IBEND)

            Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,         &
     &                  Half,KET,NASH(JSYM)**2,                         &
     &                       KET,NASH(JSYM)**2,                         &
     &                  Zero,INT2,NASH(JSYM)**2)
          End Do
          call mma_deallocate(KET)
          call mma_deallocate(BGRP)
        End Do
      Else
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI

            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT1
!           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
!             INT1(iAshI,jAshI) = INT1(iAshI,jAshI)
!    *          + Two*WRK1(iCorI,iCorI)
!           End Do
            !! Put in INT2
            Do kAshI = 1, nAsh(iSym)
              Do lAshI = 1, nAsh(iSym)
                INT2(iAshI,jAshI,kAshI,lAshI)                           &
     &        = INT2(iAshI,jAshI,kAshI,lAshI)                           &
     &        + WRK1(nCorI+kAshI,nCorI+lAshI)*Half
              End Do
            End Do

!           Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT1
!           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
!             INT1(iAshI,jAshI) = INT1(iAshI,jAshI) - WRK1(iCorI,iCorI)
!           End Do
          End Do
        End Do
      End If
#ifdef _MOLCAS_MPP_
      call GADGOP(INT2,nAshT**4,'+')
#endif
!     write(u6,*) 'int2'
!     call sqprt(int2,25)
!     call sqprt(int1,5)
!     call sqprt(fimo,12)
      If (Mode == 0) Then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          iTU = iT + nAshT*(iU-1)
          Do iV = 1, nAshT
            Do iX = 1, nAshT
              iVX = iV + nAshT*(iX-1)
              If (iVX > iTU) Then
               INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX) + INT2(iV,iX,iT,iU)
               INT2(iV,iX,iT,iU) = Zero
              End If
            End Do
          End Do
        End Do
      End Do
      End If

      if (mode == 0 .or. mode == 1) then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          Do iX = 1, nAshT
            INT1(IT,IU) = INT1(IT,IU) - INT2(IT,IX,IX,IU)
          End Do
        End Do
      End Do
      endif

      Return

      End Subroutine CnstInt
!
!-----------------------------------------------------------------------
!
      !! dens2_rpt2.f
      Subroutine TimesE2(Mode,nConf,nState,nAshT,CIin,CIout,INT1,INT2)

      use sguga, only: SGS, L2ACT, CIS
      use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk
      use Constants, only: Two

      implicit none

      integer(kind=iwp), intent(in) :: Mode, nConf, nState, nAshT
      real(kind=wp), intent(in) :: CIin(nConf,nState),                  &
     &  INT1(nAshT,nAshT), INT2(nAshT,nAshT,nAshT,nAshT)
      real(kind=wp), intent(out) :: CIout(nConf,nState)

      real(kind=wp), allocatable :: SGM1(:), SGM2(:)
      integer(kind=iwp), allocatable :: TASK(:,:)

      integer(kind=iwp) :: nLev, nTasks, iTask, LT, LU, kState, IST,    &
     &  IT, ISU, IU, ISTU, ISSG, NSGM, LVX, LV, ISV, IV, LX, ISX, ISVX, &
     &  IX, ilStat, jlStat
      real(kind=wp) :: Ovl

      nLev=SGS%nLev
      ! logical tras,uras,vras,xras
!
!     --- H_{IJ}*P_J
!    <CI1|EtuEvx|CI2>=<CI1|Evx
!
      nTasks = nLev**2
      CALL mma_allocate (Task,nTasks,2,Label='TASK')

      iTask=0
      DO LT=1,nLev
        DO LU=1,nLev
          iTask=iTask+1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        ENDDO
      ENDDO
      IF (iTask /= nTasks) WRITE(u6,*) 'ERROR nTasks'

      call mma_allocate(SGM1,nConf,Label='SGM1')
      call mma_allocate(SGM2,nConf,Label='SGM2')

      CIout(1:nConf,1:nState) = Zero
      Do kState = 1, nState
        !! Start the actual part of dens2_rpt2
        Call Init_Tsk(ID, nTasks)

        do while (Rsv_Tsk(ID,iTask))
          LT=TASK(iTask,1)
          ! tras=.false.
          ! if (lt <= nras1(1)) tras=.true.
          IST=SGS%ISM(LT)
          IT=L2ACT(LT)
          LU=Task(iTask,2)
          ! uras=.false.
          ! if (lu > nras1(1)+nras2(1)) uras=.true.
!         if (tras.and.uras) cycle
          ! LTU=iTask
          ISU=SGS%ISM(LU)
          IU=L2ACT(LU)
          ISTU=Mul(IST,ISU)
          ISSG=Mul(ISTU,STSYM)
          NSGM=CIS%NCSF(ISSG)
          IF(NSGM == 0) cycle
          !! <CIin|Etu
          CALL GETSGM2(LU,LT,STSYM,CIin(1,kState),nConf,SGM1,NSGM)
          IF(ISTU == 1) THEN
            !! <CIin|Etu|CIout>*I1tu
            CIout(1:NSGM,kState) = CIout(1:NSGM,kState)                 &
     &        + INT1(IT,IU)*SGM1(1:NSGM)
          END IF
          LVX=0
          DO LV=1,NLEV
            ISV=SGS%ISM(LV)
            IV=L2ACT(LV)
            ! vras=.false.
            ! if (lv <= nras1(1)) vras=.true.
            DO LX=1,NLEV
              LVX=LVX+1
              ISX=SGS%ISM(LX)
              ISVX=Mul(ISV,ISX)
              ! xras=.false.
              ! if (lx > nras1(1)+nras2(1)) xras=.true.
!             if (vras.and.xras) cycle
              IF (ISVX /= ISTU) cycle
              IX=L2ACT(LX)
              CALL GETSGM2(LX,LV,ISSG,SGM1,NSGM,SGM2,NSGM)
              CIout(1:NSGM,kState) = CIout(1:NSGM,kState)               &
     &          + INT2(IT,IU,IV,IX)*SGM2(1:NSGM)
            END DO
          END DO
        end do
        CALL Free_Tsk(ID)
        !! End the actual part of dens2_rpt2
      End Do

      call mma_deallocate(Task)

#ifdef _MOLCAS_MPP_
      CALL GAdGOP(CIout,nConf*nState,'+')
#endif
!
!     --- -E_{S}*CJ + zL_{KL}
!
      Do kState = 1, nState
        CIout(1:nConf,kState)                                           &
     &    = CIout(1:nConf,kState) + Eact(kState)*CIin(1:nConf,kState)
      End Do

      !! Project out the reference vector, just in case
      If (Mode == 1) Then
        Do ilStat = 1, nState
          SGM1(1:nConf) = CIout(1:nConf,ilStat)
          Do jlStat = 1, nState
            Call LoadCI_XMS('C',1,nConf,nState,SGM2,jlStat,U0)
            Ovl = DDot_(nConf,SGM1,1,SGM2,1)
            CIout(1:nConf,ilStat)                                       &
     &        = CIout(1:nConf,ilStat) - Ovl*SGM2(1:nConf)
          End Do
        End Do
      End If

      call mma_deallocate(SGM1)
      call mma_deallocate(SGM2)

      CIout(1:nConf,1:nState) = Two/nState*CIout(1:nConf,1:nState)

      Return

      End Subroutine TimesE2
!
!-----------------------------------------------------------------------
!
      Subroutine CnstDEPSA(nConf,nState,nAshT,CI,CIT,G1,G2,INT2)

      use sguga, only: SGS
      use caspt2_module, only: MXCI, NG1, NG2

      implicit none

      integer(kind=iwp), intent(in) :: nConf, nState, nAshT
      real(kind=wp), intent(in) :: CI(nConf,nState), CIT(nConf,nState), &
     &  INT2(nAshT,nAshT,nAshT,nAshT)
      real(kind=wp), intent(out) :: G1(nAshT,nAshT),                    &
     &  G2(nAshT,nAshT,nAshT,nAshT)

      real(kind=wp),allocatable :: SGM1(:),SGM2(:),G1T(:),G2T(:),       &
     &                             Fock(:),FockOut(:)
      integer(kind=iwp) :: nLev, kState, ilState, jlState, iS, jS, iA,  &
     &  jA, ip1, ip2, ipS, ijS, kS, lS, kAsh, kAA, lAsh, lAA, iAsh, ipQ,&
     &  jAsh, ipM, imo, iOrb, jOrb
      real(kind=wp) :: Wgt, vSLag, rd, EigI, EigJ, OLagIJ, Tmp

      nLev=SGS%nLev
!
!     This subroutine computes the second term in Eq. (70) or the RHS of
!     Eq. (72) in the CASPT2-IPEA gradient paper
!     CIT corresponds to \overline{Q}, if I remember correctly
!
      G1(:,:) = Zero
      G2(:,:,:,:) = Zero

      !! Construct transition(?) density matrix
      !! (<CI|Etu|CIT>+<CIT|Etu|CI>)/2, where CIT is the solution
      call mma_allocate(SGM1,MXCI,Label='SGM1')
      call mma_allocate(SGM2,MXCI,Label='SGM2')
      call mma_allocate(G1T,NG1,Label='GT1')
      call mma_allocate(G2T,NG2,Label='GT2')
!
!  !! This is for CASSCF orbital Lagrangian, but this may not contribute
!     Call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,jState),CI(1,jState),
!    *                 SGM1,SGM2,G1T,G2T)
!     Call DaXpY_(NG1,-Half,G1T,1,G1,1)
!     Call DaXpY_(NG2,-Half,G2T,1,G2,1)
!
      Do kState = 1, nState
!       Wgt = DWgt(iState,iState)
        Wgt = One/nState

        !! <CI|Etu|CIT>+<CIT|Etu|CI> and the t+ u+ x v variant
        Call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,kState),CIT(1,kState),    &
     &                   SGM1,SGM2,G1T,G2T)
        Call DaXpY_(NG1,WGT,G1T,1,G1,1)
        Call DaXpY_(NG2,WGT,G2T,1,G2,1)

        !! For the orbital contribution of CASSCF Lagrangian
        !! Just add the SLag rotation contributions
        ilState = kState
        Do jlState = 1, ilState-1
          vSLag = -Half*SLag(ilState,jlState)
          If (abs(vSLag) <= 1.0e-08_wp) Cycle
          Call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,ilState),CI(1,jlState), &
     &                     SGM1,SGM2,G1T,G2T)
          Call DaXpY_(NG1,vSLag,G1T,1,G1,1)
          Call DaXpY_(NG2,vSLag,G2T,1,G2,1)
        End Do
      End Do

      call mma_deallocate(SGM1)
      call mma_deallocate(SGM2)
      call mma_deallocate(G1T)
      call mma_deallocate(G2T)

      !! Finally, construct the Fock matrix only for active-active
      !! Should be equivalent to FockGen in MCLR
      call mma_allocate(Fock,nAshT**2,Label='Fock')
      Fock(:) = Zero

      !! 1) FIMO term
      Do iS=1,nSym
        If (nBas(iS) > 0) Then
          jS=iEOr(is-1,iSym-1)+1
          Do iA=1,nAsh(is)
            Do jA=1,nAsh(js)
!             rd=rDens1(iA+nA(iS),jA+nA(js))
!             ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)-1
!             ip2=nBas(iS)*(nIsh(js)+jA-1) +ipmat(is,js)
              rd=G1(iA,jA)
              ip1= 1+nFro(jS)+nIsh(jS)                                  &
     &           + nBas(iS)*(nFro(iS)+nIsh(iS)+iA-1)
              ip2=1+nAsh(iS)*(jA-1)
              Call DaXpY_(nAsh(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
            End Do
          End Do
        End If
      End Do
!     Write(u6,*) 'after 1'
!     call sqprt(fock,nasht)

      !! 2) two-electron term (only CreQADD part)
      Do iS=1,nSym
        ipS=iEOr(is-1,isym-1)+1
        if (norb(ips) /= 0) Then
          Do jS=1,nsym
            ijS=iEOR(is-1,js-1)+1
            Do kS=1,nSym
              ls=iEOr(ijs-1,iEor(ks-1,isym-1))+1
!                                                                      *
!***********************************************************************
!                                                                      *
               Do kAsh=1,nAsh(kS)
                kAA=kAsh+nFro(kS)+nIsh(kS)
                Do lAsh=1,nAsh(lS)
                  lAA=lAsh+nFro(lS)+nIsh(lS)
!
!                 Pick up (pj|kl)
!
                  Call Coul(ipS,jS,kS,lS,kAA,lAA,WRK1,WRK2)
!
                  Do iAsh=1,nAsh(iS)
                    ipQ=nAsh(ipS)*(iAsh-1)
                    Do jAsh=1,nAsh(jS)
                      ipM=nFro(ipS)+nIsh(ipS)                           &
     &                   +(nFro(jS)+nIsh(jS)+jAsh-1)*nBas(ipS)
                      Call DaXpY_(nAsh(ipS),G2(iAsh,jAsh,kAsh,lAsh)*2,  &
     &                            INT2(1,jAsh,kAsh,lAsh),1,             &
     &                            Fock(1+ipQ),1)
                      ipM=ipM+nOrb(ipS)
!
                    End Do
                  End Do
!
                End Do
              End Do
!                                                                      *
!***********************************************************************
!                                                                      *
            End Do  ! kS
          End Do     ! jS
        End If
      End Do           ! iS

      !! 3) anti-symmetrize
      !! 4) Divide by the difference of orbital energies
      call mma_allocate(FockOut,nAshT**2,Label='FockOut')
      Do iS=1,nSym
        jS=iEOR(iS-1,iSym-1)+1
        If (nAsh(is)*nAsh(jS) /= 0) Then
          !! Anti-symmetrize
          Call DGeSub(Fock,nAsh(iS),'N',                                &
     &                Fock,nAsh(jS),'T',                                &
     &                FockOut,nAsh(iS),                                 &
     &                nAsh(iS),nAsh(jS))


          !! Divide
          imo=1
          Do iAsh = 1, nAsh(iSym)
            iOrb = iAsh + nFro(iSym) + nIsh(iSym)
            EigI = FIFA(iMO+iOrb-1+nBas(iSym)*(iOrb-1))
            Do jAsh = 1, iAsh-1
              jOrb = jAsh + nFro(iSym) + nIsh(iSym)
              EigJ = FIFA(iMO+jOrb-1+nBas(iSym)*(jOrb-1))
              OLagIJ = FockOut(iAsh+nAsh(iSym)*(jAsh-1))
              Tmp = OLagIJ/(EigI-EigJ)
              DEPSA(iAsh,jAsh) = DEPSA(iAsh,jAsh) + Tmp
              DEPSA(jAsh,iAsh) = DEPSA(jAsh,iAsh) + Tmp
            End Do
          End Do
        End If
      End Do

      call mma_deallocate(FockOut)
      call mma_deallocate(Fock)

      Return

      End Subroutine CnstDEPSA
!
!-----------------------------------------------------------------------
!
      !! PRWF1_CP2
      SUBROUTINE CnstPrec(ISYCI,NCONF,NROOTS,NLEV,nMidV,PRE,CI,INT1,    &
     &                    INT2,Fancy)
      use molcas, only: MXLEV
      use sguga, only: SGS, CIS
      use Constants, only: Two, Four

      implicit none

#include "intent.fh"

      integer(kind=iwp), INTENT(IN) :: ISYCI, NCONF, NROOTS, nLev, nMidV
      real(kind=wp), intent(in) ::  CI(nConf*nRoots), INT1(NLEV,NLEV),  &
     &  INT2(NLEV,NLEV,NLEV,NLEV)
      real(kind=wp), intent(_OUT_) :: PRE(nConf)
      real(kind=wp), intent(out) :: Fancy(nRoots,nRoots,nRoots)

      real(kind=wp) ::  ICS(MXLEV), val, val2, Ene, dnum
      integer(kind=iwp) :: nIpWlk, LENCSF, ISY, LEV, MV, ISYUP, NCI,    &
     &  NUP, ISYDWN, NDWN, ICONF, IUW0, IDW0, IDWN, IUP, ICDPOS, ICDWN, &
     &  NNN, IC1, ICUP, K, IDWNSV, ICUPOS, LEV2, iSt, jSt, kSt

      nIpWlk = CIS%nIpWlk
!
!     Construct (approximate?) preconditioner for the active linear
!     equation that should be solved for non-invariant CASPT2 methods
!     (with IPEA shift)
!
! -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
! -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
! -- THE   U P P E R   PART OF THE WALK.

! SVC: set up a CSF string length as LENCSF
!     LINE=' '
      LENCSF=0
      ISY=0
      DO LEV=1,NLEV
        IF(ISY /= SGS%ISM(LEV)) THEN
          ISY=SGS%ISM(LEV)
          LENCSF=LENCSF+1
        END IF
        LENCSF=LENCSF+1
      END DO
      LENCSF=MIN(LENCSF,256)
      LENCSF=MAX(LENCSF,10)

!     LINE=' '

! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NCI=CIS%NOCSF(ISYUP,MV,ISYCI)
          IF(NCI == 0) cycle
          NUP=CIS%NOW(1,ISYUP,MV)
          ISYDWN=Mul(ISYUP,ISYCI)
          NDWN=CIS%NOW(2,ISYDWN,MV)
          ICONF=CIS%IOCSF(ISYUP,MV,ISYCI)
          IUW0=1-NIPWLK+CIS%IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+CIS%IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO IDWN=1,NDWN
            DO IUP=1,NUP
              ICONF=ICONF+1
!             COEF=CI(ICONF)
! -- SKIP OR PRINT IT OUT?
!             IF(ABS(COEF) < THR) cycle
              IF(IDWNSV /= IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=CIS%ICASE(ICDPOS)
! -- UNPACK LOWER WALK.
                NNN=0
                DO LEV=1,SGS%MIDLEV
                  NNN=NNN+1
                  IF(NNN == 16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=CIS%ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
                END DO
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=CIS%ICASE(ICUPOS)
! -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=SGS%MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN == 16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=CIS%ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
! -- PRINT IT!
              K=0
              ISY=0
              PRE(ICONF) = Zero
              DO LEV=1,NLEV
                IF(ISY /= SGS%ISM(LEV)) THEN
                  ISY=SGS%ISM(LEV)
                  K=K+1
!                 LINE(K:K)=' '
                END IF
                K=K+1
!               LINE(K:K)=CODE(ICS(LEV))
                IF (ICS(LEV) == 0) THEN
                  VAL = Zero
                ELSE IF (ICS(LEV) == 3) THEN
                  VAL = Two*INT1(LEV,LEV)
!                 L=0
!                 JSY=0
                  DO LEV2=1,NLEV
                    IF (ICS(LEV2) == 0) THEN
                    ELSE IF ((LEV == LEV2 .AND. ICS(LEV2) == 3).OR.     &
     &                       (LEV /= LEV2 .AND. ICS(LEV2) == 1).OR.     &
     &                       (LEV /= LEV2 .AND. ICS(LEV2) == 2)) THEN
                      val2 =  Four*int2(lev,lev ,lev2,lev2)             &
     &                      - Two*int2(lev,lev2,lev ,lev2)
                      val = val + val2
                    ELSE IF (LEV /= LEV2 .AND. ICS(LEV2) == 3) THEN
                      val2 =  Four*int2(lev,lev ,lev2,lev2)             &
     &                      - Two*int2(lev,lev2,lev ,lev2)
                      val = val + val2
                    END IF
                  END DO
                ELSE
                  VAL = INT1(LEV,LEV)
!                 L=0
!                 JSY=0
                  DO LEV2=1,NLEV
                    IF (ICS(LEV2) == 0 .OR. LEV == LEV2) THEN
                    ELSE IF (ICS(LEV2) == 3) THEN
!      val2 =  Two*int2(lev,lev ,lev2,lev2)
!    *       - One*int2(lev,lev2,lev ,lev2)
!      val2 = Zero
!      val = val + val2*Half
                    ELSE
                      val2 = int2(lev,lev ,lev2,lev2)                   &
     &                     + int2(lev,lev2,lev ,lev2)
                      if (ics(lev) == ics(lev2)) then
                      val2 = int2(lev,lev ,lev2,lev2)                   &
     &                     - int2(lev,lev2,lev ,lev2)
                      end if
                      val = val + val2
                    END IF
                  END DO
                END IF
                PRE(ICONF) = PRE(ICONF) + VAL
              END DO
            end do
          end do
        end do
      end do

      !! mclr/sa_prec.f
      !! Prepare so-called fancy preconditioner
      Do iSt = 1, nRoots
        Ene = Eact(iSt)
        Do jSt = 1, nRoots
          Do kSt = 1, nRoots
            Fancy(jSt,kSt,iSt) = Zero
            Do iConf = 1, nConf
              dnum=PRE(iConf)+Ene
              dnum=Sign(Max(Abs(dnum),1.0e-16_wp),dnum)
              Fancy(jSt,kSt,iSt) = Fancy(jSt,kSt,iSt)                   &
     &          + CI(iConf+nConf*(jSt-1))*CI(iConf+nConf*(kSt-1))/dnum
            End Do
          End Do
        End Do
        Call MatInvert(Fancy(1,1,iSt),nRoots)
      End Do

      RETURN

      END SUBROUTINE CnstPrec
!
!-----------------------------------------------------------------------
!
      !! mclr/dminvci_sa.f
      Subroutine DoPrec(nConf,nRoots,VecIN,VecOUT,CI,Pre,Fancy)
!
!     Apply precondition to CI vectors, taken from the MCLR module
!
      implicit none

      integer(kind=iwp), intent(in) :: nConf, nRoots
      real(kind=wp), intent(in) :: VecIN(nConf,nRoots),                 &
     &  Pre(nConf), Fancy(nRoots,nRoots,nRoots)
      real(kind=wp), intent(out) :: VecOUT(nConf,nRoots),               &
     &  CI(nConf,nRoots)

      integer(kind=iwp) :: iRoots, iConf, jRoots, kRoots
      real(kind=wp) :: rcoeff(nRoots), alpha(nRoots)

      !! Standard inverse of the diagonal elements
      Do iRoots = 1, nRoots
        Do iConf = 1, nConf
          VecOUT(iConf,iRoots)                                          &
     &      = VecIN(iConf,iRoots)/(Pre(iConf)+Eact(iRoots))
        End Do
      End Do

      !! Construct reference CI vectors
      Do iRoots = 1, nRoots
        Call LoadCI_XMS('C',1,nConf,nState,CI(1,iRoots),iRoots,U0)
      End Do

      !! The so-called fancy precondioner
      Do iRoots = 1, nRoots
        Do jRoots = 1, nRoots
         rcoeff(jRoots) = DDot_(nconf,VecOUT(1,iRoots),1,CI(1,jRoots),1)
        End Do

        Do jRoots = 1, nRoots
          alpha(jRoots) = Zero
          Do kRoots = 1, nRoots
            alpha(jRoots) = alpha(jRoots)                               &
     &        + Fancy(jRoots,kRoots,iRoots)*rcoeff(kRoots)
          End Do
        End Do

        Do jRoots = 1, nRoots
          Do iConf = 1, nConf
            VecOUT(iConf,iRoots) = VecOUT(iConf,iRoots)                 &
     &        - CI(iConf,jRoots)*alpha(jRoots)/(Pre(iConf)+Eact(iRoots))
          End Do
        End Do
      End Do

      End Subroutine DoPrec

      End Subroutine DEPSAOffC
