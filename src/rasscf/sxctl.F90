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
! Copyright (C) 1989, Bjorn O. Roos                                    *
!               1989, Per Ake Malmqvist                                *
!               1991,1993,1996, Markus P. Fuelscher                    *
!***********************************************************************
      Subroutine SXCtl(CMO,OCC,D,P,PA,FI,FA,D1A,THMAX,IFINAL)
!*******************************************************************
!                                                                  *
! super-CI control section                                         *
!                                                                  *
! calling arguments:                                               *
! CMO     : array of real*8                                        *
!           MO-coefficients                                        *
! OCC     : array of real*8                                        *
!           orbital occupation numbers                             *
! D       : array of real*8                                        *
!           averaged one-body density matrix                       *
! P       : array of real*8                                        *
!           averaged two-body density matrix                       *
! PA      : array of real*8                                        *
!           averaged antisymmetric two-body density matrix         *
! FI      : array of real*8                                        *
!           inactive Fock matrix                                   *
! FA      : array of real*8                                        *
!           active Fock matrix                                     *
! D1A     : array of real*8                                        *
!           active one-body density matrix in AO basis             *
! IFINAL  : integer                                                *
!           termination flag                                       *
!                                                                  *
!------------------------------------------------------------------*
!                                                                  *
! written by:                                                      *
! B.O. Roos and P.Aa. Malmqvist                                    *
! University of Lund, Sweden, 1989                                 *
!                                                                  *
!------------------------------------------------------------------*
!                                                                  *
! history:                                                         *
! - updated for MOLCAS version 2                                   *
!   M.P. Fuelscher, University of Lund, Sweden, 1991               *
! - updated for MOLCAS version 3                                   *
!   M.P. Fuelscher, University of Lund, Sweden, 1993               *
! - updated for integral direct and reaction field calculations    *
!   M.P. Fuelscher, University of Lund, Sweden, 1996               *
!                                                                  *
!*******************************************************************

#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fciqmc, only : DoNECI
#ifdef _HDF5_
      use mh5, only: mh5_put_dset
      use raswfn, only: wfn_mocoef, wfn_occnum
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use Fock_util_global, only: ALGO, DoCholesky
      use Lucia_Interface, only: Lucia_Util
      use wadr, only: DIA, SXN, BM, F1, F2, SXG, SXH, NLX
      use input_ras, only: KeyHEUR
      use rasscf_global, only: ExFac, KSDFT, DoBlockDMRG, DoDMRG, ECAS, &
     &                         ESX, iCIOnly, iPT2, ITER, ITERSX, ITMAX, &
     &                         l_casdft, NAC, nDimSX, nFint, NO2M,      &
     &                         nQune, NROOT, NSXS, NTOT4, QNSTEP,       &
     &                         QNUPDT, SXSEL, TMIN, VIA, ISTORP,        &
     &                         IADR15, EMY
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NACTEL,JOBIPH,LUINTM,LUQUNE,NASH,    &
     &                        NBAS,NDEL,NFRO,NISH,NORB,NRS1,NRS2,NRS3,  &
     &                        NSSH,NTOT,NTOT1,NTOT2


      Implicit None

      Real*8 CMO(*),OCC(*),D(*),P(*),PA(*),FI(*),FA(*),D1A(*)
      Real*8 THMAX
      INTEGER IFINAL

      Character(LEN=16), Parameter :: ROUTINE='SXCTL   '
! PAM 2008 IndType, VecTyp added, see below at call to WrVec
      Integer IndType(56)
      Character(LEN=80) VecTyp
      Integer, Save ::nCall
      Logical TraOnly
      Real*8 P2act(1),CIDUMMY(1)
      Real*8, Allocatable:: SXHD(:)
      Real*8, Allocatable, Target:: SMAT(:)
      Real*8, Allocatable:: PUVX(:), DA(:), STRP(:), P2reo(:), P2Raw(:),&
     &                      Fck(:), QMat(:), EDum(:), CMON(:), FTR(:),  &
     &                      Vec(:), WO(:), SQ(:), CMOX(:), SXDF(:),     &
     &                      SXDD(:), CSX(:), Sigma(:), HH(:), CC(:),    &
     &                      ENER_X(:), SC(:), QQ(:), OVL(:), VT(:),     &
     &                      VL(:), XQN(:), SCR(:), V1(:), V2(:),        &
     &                      XMAT(:), X2(:)
      Real*8 :: Dummy(1), CASDFT_En, CPES, CPTS, P2reo_size, TIOES,     &
     &          TIOS, XSXMAX
      Real*8, External :: DDot_
      Integer :: i, iBas, IC, iDisk, IndT, iOff, iOrb, iPrLev, iRC,     &
     &           iRef, iShift, iSym, iWay, j, kMax, LCSXI, LuvvVec,     &
     &           MIAAE, MNO, NA, NAEAE, NAOAE, NB, NCR, NCR1, NE, NI,   &
     &           NIA, NIAIA, NLCC, NLHH, NLOVL, NLQ, NLX1, NLX2, NO,    &
     &           nP2Act, NQ, NAE
      Integer, External:: IsFreeUnit

! PAM01 The SXCI part has been slightly modified by P-AA M Jan 15, 2001:
! Changes affect several of the subroutines of this part.
! The program should now behave more gracefully when there are (almost)
! redundant orbital rotation modes.
! For individual 2x2 rotations that are to be excluded, the SX hamiltonian
! diagonal, used in the preconditioning step in DAVCRE, is set to a huge
! value. Corresponding elements of trial vectors and sigma vectors are
! zeroed.
! To take care of 'hidden' redundant linear combinations, a small extra
! term, proportional to the square norm of the orbital rotations, is
! added to the optimization problem: this means in practice that a small
! extra term is added in the SIGVEC routine.
! Also, an extra term is added to the overlap calculation. The
! proportionality factor in SIGVEC should be small, but larger than the
! one used in the COVLP routine.
! Presently I try the following values:
! Huge elements in diagonal values: 1.0D32
! (Tested against 1.0D20 in IF-statements)
! Extra term in overlaps (COVLP, SXHAM): 1.0D-14
! Extra term in SIGVEC:                  1.0D-12


! Local print level (if any)
      IPRLEV=IPRLOC(4)
!      write(6,*) 'Entering SXCTL!'
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

! --- Check for Cholesky ---------------

!      Call DecideOnCholesky(DoCholesky)
! --------------------------------------

! compute constants needed for addressing
      NLX1=0
      NLX2=0
      NQ=0
      NSXS=0
      NIAIA=0
      NAEAE=0
      NAOAE=0
      MNO=0
      MIAAE=0
      NLX=0
      DO 10 ISYM=1,NSYM
       NB=NBAS(ISYM)
       NA=NASH(ISYM)
       NI=NISH(ISYM)
       NE=NSSH(ISYM)
       NO=NORB(ISYM)
       NIA=NI+NA
       NAE=NA+NE
       NSXS=NSXS+NIA*NAE
       NIAIA=NIAIA+NIA**2
       NAEAE=NAEAE+NAE**2
       NAOAE=NAOAE+NA*NAE
       NLX1=MAX(NLX1,NO**2)
       NLX2=MAX(NLX2,NB*NO)
       NQ=MAX(NQ,NA*NO)
       MNO=MAX(MNO,NO)
       MIAAE=MAX(MIAAE,NIA*NAE)
       NLX=MAX(NLX,NIA*NA,NIA*NAE)
10    CONTINUE
      IF(NQ.LT.NIAIA) NQ=NIAIA
      NROOT=1
      NDIMSX=NROOT+NSXS
!***********************************************************************
! load back two-electron integrals (pu|vx)
!***********************************************************************
      If (.not.DoCholesky .or. ALGO.eq.1) Then
        If ( nFint.gt.0) then
          iDisk = 0
          Call mma_allocate(PUVX,nFint,Label='PUVX')
          Call DDaFile(LUINTM,2,PUVX,nFint,iDisk)
        Else
          Call mma_allocate(PUVX,1,Label='PUVX')
        EndIf
      End If
      IF(IPRLEV.ge.DEBUG) THEN
       write(6,*) 'PUVX integrals in SXCTL'
       call wrtmat(PUVX,1,nFInt,1,nFInt)
      END IF
!********************************************************************************
! update and transform the Fock matrices FI and FA in MO basis ----> Fmat routine
!********************************************************************************
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call Fmat(CMO,PUVX,D,D1A,FI,FA)

      ElseIf (ALGO.eq.2) Then

!     Inactive-active contribution to ECAS
         Call mma_allocate(DA,nTot1,Label='DA')
         Call Fold(nSym,nBas,D1A,DA) !get the packed DA
         VIA=dDot_(nTot1,FI,1,DA,1)
         ECAS=EMY+VIA
         If ( IPRLEV.ge.DEBUG ) then
           Write(LF,'(A,ES20.10)') ' Total core energy:            ',EMY
           Write(LF,'(A,ES20.10)') ' inactive-active interaction:  ',VIA
           Write(LF,'(A,ES20.10)') ' CAS energy (core+interaction):',   &
     &                             ECAS
         End If
         Call mma_deallocate(DA)

         TraOnly=.true.
         Call CHO_CAS_DRV(irc,CMO,D,FI,D1A,FA,Dummy,TraOnly)

         if (irc.ne.0) then
         Write(LF,*)'SXCTL: Cho_cas_drv non-zero return code! rc= ',irc
         call abend()
         endif

      Else

         Write(LF,*)'SXCTL: Illegal Cholesky parameter ALGO= ',ALGO
         call abend()

      EndIf
!***********************************************************************
! reorder the two-body density matrix P
!***********************************************************************
      IF(.not.l_casdft) then
! ISTORP(NSYM+1) represents the size of the 2-body density matrix,d(vwxy), with vwxy all active.
! the size is computed as NAP*NAQ*NRS (sum over all symmetries). If Sym_R = Sym_S then triangular
! form over NRS... with R.ge.S, rectanguar otherwise.
       IF(ISTORP(NSYM+1).GT.0) THEN
         CALL mma_allocate(STRP,ISTORP(NSYM+1),Label='STRP')

         CALL PMAT_RASSCF(P,STRP)

         If (ExFac.ne.1.0D0 .and.KSDFT(1:3).ne.'SCF') Then
            CALL mma_allocate(P2reo,ISTORP(NSYM+1),Label='P2reo')
            Call Get_Temp('nP2Act  ',P2Act,1)
            nP2Act=Int(P2Act(1))
            CALL mma_allocate(P2RAW,nP2Act,Label='P2Raw')
            Call Get_Temp('P2_RAW  ',P2RAW,nP2Act)
            CALL PMAT_RASSCF(P2RAW,P2reo)
            Call mma_deallocate(P2Raw)
            P2reo_size=DBLE(ISTORP(NSYM+1))
            Call Put_Temp('nP2reo  ',[P2reo_size],1)
            Call Put_Temp('P2_reo  ',P2reo,ISTORP(NSYM+1))
            Call mma_deallocate(P2reo)
         End If
       ELSE
         CALL mma_allocate(STRP,1,Label='STRP')
       END IF
      ELSE ! GLM-CASDFT
! ISTORP(NSYM+1) here represents the size of the Dvw*Dxy array (product of one-body
! density matrix,d(vwxy), with vwxy all active. The size is computed as NAP*NAQ*NRS
! (sum over all symmetries). If Sym_R = Sym_S then triangular form over NRS...
! with R.ge.S, rectanguar otherwise. Basically we will use same simmetry as for dvwxy.
       IF(ISTORP(NSYM+1).GT.0) THEN
!         Write(LF,*)
!         Write(LF,*) ' ---------------------'
         CALL mma_allocate(STRP,ISTORP(NSYM+1),Label='STRP')
         CALL DmatDmat(D,STRP)
       ELSE
         CALL mma_allocate(STRP,1,Label='STRP')
       END IF
      end if
!***********************************************************************
! Compute the MCSCF generalized Fock matrix and Brillouin matrix elements
!***********************************************************************
      Call mma_allocate(FCK,NTOT4,Label='FCK')
      CALL mma_allocate(BM,NSXS,Label='BM')
      CALL mma_allocate(QMat,NQ,Label='QMat') ! q-matrix(1symmblock)
      CALL FOCK(FCK,BM,FI,FA,D,STRP,QMat,PUVX,IFINAL,CMO)
! Now FA = FI + FA. Original FA has been overwritten in FOCK routine.
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)
        Write(LF,*)'FI+FA in MO-basis in sxctl (overwritten on FA)'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      END IF
      Call mma_deallocate(FCK)
      Call mma_deallocate(QMat)
      call mma_deallocate(STRP,safe='*')
      Call mma_deallocate(PUVX,safe='*')

! PAM 2008: Orbital files should be updated each iteration
! for easy access in case of catastrophic failure.
      IF ( IFINAL.NE.1 ) THEN
       iShift=0
       DO ISYM=1,NSYM
         IndT=0
         IndType(1+iShift)= NFRO(ISYM)
         IndT=IndT+NFRO(ISYM)
         IndType(2+iShift)= NISH(ISYM)
         IndT=IndT+NISH(ISYM)
         IndType(3+iShift)= NRS1(ISYM)
         IndT=IndT+NRS1(ISYM)
         IndType(4+iShift)= NRS2(ISYM)
         IndT=IndT+NRS2(ISYM)
         IndType(5+iShift)= NRS3(ISYM)
         IndT=IndT+NRS3(ISYM)
         IndType(7+iShift)= NDEL(ISYM)
         IndT=IndT+NDEL(ISYM)
         IndType(6+iShift)= NBAS(ISYM)-IndT
         iShift=iShift+7
        EndDo
! Note: This is not the final orbitals, and the orbital energies and
! active occupation numbers may be meaningless.
! There is an array with occupation numbers, so use it, even if
! possibly irrelevant. But put zeroes as orbital energies:
        Call mma_allocate(EDUM,NTOT,Label='EDUM')
        EDUM(:)=0.0D0

        Write(VecTyp,'(A)')
        VecTyp='* RASSCF average (pseudo-natural) orbitals (Not final)'
        LuvvVec=50
        LuvvVec=isfreeunit(LuvvVec)
        Call WrVec('RASORB',LuvvVec,'COE',NSYM,NBAS,NBAS,               &
     &           CMO, OCC, EDUM, INDTYPE,VECTYP)
        Call WrVec('RASORB',LuvvVec,'AI',NSYM,NBAS,NBAS,                &
     &           CMO, OCC, EDUM, INDTYPE,VECTYP)
        Call mma_deallocate(EDUM)

#ifdef _HDF5_
        call mh5_put_dset(wfn_mocoef,CMO)
        call mh5_put_dset(wfn_occnum,OCC)
#endif
      END IF

      IF ( IFINAL.EQ.1 ) THEN
! If ifinal=1 , this is the last call to SXCTL for calculating
! the MCSCF Fock matrix for occupied orbitals, and new orbitals).
! First generate final orbitals:
! Diagonalize inactive and secondary part of FP = FI + FA.
!         diagonal blocks of the active density matrix
! for RAS1, RAS2, and RAS3. All this is done in NEWORB.
! Finally generate orbitals and Fock matrices for CAS-PT2 runs.
! in FCKPT2

! Memory allocation and call to NEWORB and FCKPT2
! CMON: New molecular orbitals (NTOT2)
! FTR:  Temporary area for part of the Fock matrix FP (NTOT1)
! VEC:  EIGENVECTORS OF FTR (NO2M)
! SQ and WO: scratch areas

        CALL mma_allocate(CMON,NTOT2,Label='CMON')
        CALL mma_allocate(FTR,NTOT1,Label='FTR')
        CALL mma_allocate(VEC,NTOT2,Label='VEC')
        CALL mma_allocate(WO,NTOT2,Label='WO')
        CALL mma_allocate(SQ,NTOT2,Label='SQ')
        CALL mma_allocate(CMOX,NTOT2,Label='CMOX')
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in SXCTL for IFINAL=1'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=0
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMO(ioff+iBas*(i-1)+j),j=1,iBas)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If

        If ( iCIonly.eq.1 ) then
          IDISK=IADR15(2)
          CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
          CALL DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#ifdef _HDF5_
          call mh5_put_dset(wfn_mocoef,CMO)
          call mh5_put_dset(wfn_occnum,OCC)
#endif
        Else
!          this part (TRACI) need to be changed to "TRAMPS", not yet ! Yingjin
           CALL NEWORB_RASSCF(CMO,CMON,FA,FTR,VEC,WO,SQ,CMOX,D,OCC)
! compute orbital overlap matrix
!           IF (NACTEL.GT.0) THEN
! NN.14 Skip this when DMRG-CASSCF due to CI-vector dependency
           !IF(.NOT.(DoDMRG.or.doBlockDMRG).AND.NACTEL.GT.0) THEN
           IF(NACTEL.GT.0)THEN
             CALL mma_allocate(SMAT,NAC*NAC,Label='SMAT')
             IWAY = 1
             CALL OVLP(IWAY,CMO,CMON,SMAT)

             if(dodmrg)then
#ifdef _DMRG_
#ifdef BLUBB
               call mpsrot(mat,nac,nrs2,nsym)
#endif
#endif
             else if(doBlockDMRG .or. DoNECI)then
             else !CI
               iDisk=IADR15(4)
               CALL LUCIA_UTIL('TRACI',                                 &
     &                         iDisk=iDisk,                             &
     &                         Lu=JOBIPH,                               &
     &                         Array=SMAT(:))
             end if
             Call mma_deallocate(SMAT)
           ELSE
             CIDUMMY=1.0D0
             IDISK=IADR15(4)
             CALL DDAFILE(JOBIPH,1,CIDUMMY,1,IDISK)
           END IF
        End If

        ! IPT2 = 1 for OUTO, CANOnical option...
        IF(IPT2.NE.0)                                                   &
     &  CALL FCKPT2(CMO,CMON,FI,FA,FTR,VEC,WO,SQ,CMOX)

        Call mma_deallocate(CMOX)
        Call mma_deallocate(SQ)
        Call mma_deallocate(WO)
        Call mma_deallocate(VEC)
        Call mma_deallocate(FTR)
        Call mma_deallocate(CMON)

        CALL TIMING(CPTS,CPES,TIOS,TIOES)

        GOTO 9990
      ENDIF

! Memory allocation and calling sequence for SXHAM
! SXN: Normalization constants for super-CI vector
! F1 and F2: parts of the Fock matrix FP
! DIA: Occupied part of the density matrix (squared)
! SXG: The G matrix(used in sigvec)
! SXH: The H matrix( "    "   "   )
! SXHD: The diagonal of the super-CI Hamiltonian
! LDF: The matrix D*FP
! LDDIA: Diagonal of the density matrix (all elements one symmetry)

      CALL mma_allocate(SXN,NSXS,Label='SXN')
      CALL mma_allocate(F1,NIAIA,Label='F1')
      CALL mma_allocate(F2,NAEAE,Label='F2')
      CALL mma_allocate(DIA,NIAIA,Label='DIA')
      CALL mma_allocate(SXG,NIAIA,Label='SXG')
      CALL mma_allocate(SXH,NAOAE,Label='SXH')
      CALL mma_allocate(SXHD,NDIMSX,Label='SXHD')
      CALL mma_allocate(SXDF,NQ,Label='SXDF')
      CALL mma_allocate(SXDD,MNO,Label='SXDD')

!         CALL TRIPRT(' Dmat in MO in SXCTL bf call to SXHAM ',' ',D,NAC)
!         CALL TRIPRT(' Pmat in MO in SXCTL bf call to SXHAM ',
!     &              ' ',P,NACPAR)
!         CALL TRIPRT(' PAmat in MO in SXCTL bf call to SXHAM',
!     &              ' ',PA,NACPAR)
      CALL SXHAM(D,P,PA,FA,SXN,F1,F2,DIA,SXG,SXH,SXHD,SXDF,SXDD)

      Call mma_deallocate(SXDD)
      Call mma_deallocate(SXDF)

! PAM01 Removal of certain rotations from the BLB elements.
! Some additional rotations (besides those listed in IZROT) may
! need to be suppressed. These are rotations that are (very nearly)
! redundant -- they hardly affect the wave function at all.
! All suppressed rotations can be identified because the corresponding
! diagonal elements have been set to a huge number in SXHAM.
! Use this criterion to set some BLB elements exactly =0:
      DO I=1,NSXS
       IF(SXHD(NROOT+I).GT.1.0D20) BM(I)=0.0D0
      END DO

! MEMORY ALLOCATION AND CALLING SEQUENCE FOR SX DIAGONALIZATION

! CSX: The super-CI vectors
! SIGMA: The sigma vectors
! HH:  The Davidson H matrix
! CC:   "     "     egenvectors
! ENER: "     "     energies
! SC:   Scratch area
! QMat:    Davidson update vectors
! QQ:   Norm of update vectors
! OVL:  Overlap matrix

      NCR=NDIMSX*NROOT*ITMAX
      KMAX=ITMAX*NROOT
      NCR1=NDIMSX*NROOT*(ITMAX+1)
      CALL mma_allocate(CSX,NCR1,Label='CSX')
      CALL mma_allocate(SIGMA,NCR,Label='SIGMA')
      NLHH=KMAX**2+KMAX
      NLCC=KMAX**2
      NLQ=NDIMSX*(NROOT+1)
      NLOVL=ITMAX*NROOT**2
      CALL mma_allocate(HH,NLHH,Label='HH')
      CALL mma_allocate(CC,NLCC,Label='CC')
      CALL mma_allocate(ENER_X,KMAX,Label='ENER_X')
      CALL mma_allocate(SC,NDIMSX,Label='SC')
      Call mma_allocate(QMat,NLQ,Label='QMat')
      CALL mma_allocate(QQ,NROOT,Label='QQ')
      CALL mma_allocate(OVL,NLOVL,Label='OVL')

      CALL DAVCRE(CSX,SIGMA,HH,CC,ENER_X,SXHD,SC,                       &
     &            QMat,QQ,OVL,SXSEL,                                    &
     &            NROOT,ITMAX,NDIMSX,ITERSX,NSXS)

      ESX=ENER_X(1)
      Call mma_deallocate(SIGMA)
      Call mma_deallocate(HH)
      Call mma_deallocate(CC)
      Call mma_deallocate(ENER_X)
      Call mma_deallocate(SC)
      Call mma_deallocate(QMat)
      Call mma_deallocate(QQ)
      Call mma_deallocate(OVL)
      Call mma_deallocate(F1)
      Call mma_deallocate(F2)
      Call mma_deallocate(SXG)
      Call mma_deallocate(SXH)
      Call mma_deallocate(SXHD)

! Renormalize the SX-coefficients

      IREF=1
      LCSXI=1+NDIMSX*(IREF-1)
      IC=NROOT+LCSXI-1
      XSXMAX=0.0D0
      DO 54 I=1,NSXS
       CSX(IC+I)=SXN(I)*CSX(IC+I)/CSX(LCSXI)
       XSXMAX=MAX(XSXMAX,abs(CSX(IC+I)))
54    CONTINUE
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*) 'SXCTL after DAVCRE, Renormalized SX coeffs:'
        Write(LF,'(1X,8F14.6)')(CSX(IC+I),I=1,NSXS)
      END IF

! Step size control has been built into qune now.
!C Step length control, just for safety.
!      DO I=1,NSXS
!        VAL=CSX(IC+I)
!        CSX(IC+I)=VAL/(1.0D0+1.7D0*ABS(VAL))
!      END DO

! Intercept XSX and BM, to use (perhaps) Quasi-Newton or Line Search

!      IF(ITER.EQ.1) NCALL=0
      IF(ITER.LE.4) NCALL=0
      IF(KeyHEUR.AND.ITER.GT.10.AND.MOD(ITER,10).LT.4) NCALL=0
      if(doDMRG) then
        IF(ITER.LE.2) NCALL=0  ! YM: change 4 -> 2, for saving time
      end if
      IF(XSXMAX.GT.0.5D0) NCALL=0
      IF(NQUNE.NE.0.AND.XSXMAX.LT.0.5D0) THEN
        CALL mma_allocate(VT,NSXS,Label='VT')
        CALL mma_allocate(VL,NSXS,Label='VL')
        CALL mma_allocate(XQN,NSXS,Label='XQN')
        CALL mma_allocate(SCR,NSXS,Label='SCR')
        CALL mma_allocate(V1,NSXS,Label='V1')
        CALL mma_allocate(V2,NSXS,Label='V2')
        CASDFT_En=0.0d0
        If(KSDFT.ne.'SCF'.and.KSDFT(1:3).ne.'PAM')                      &
     &      Call Get_dScalar('CASDFT energy',CASDFT_En)
        CASDFT_En=ECAS+CASDFT_En
        CALL QUNE(NCALL,CASDFT_En,BM,CSX(NROOT+LCSXI),                  &
     &            VL,VT,XQN,SCR,                                        &
     &            V1,V2,NSXS,LUQUNE,                                    &
     &            TMIN,QNSTEP,QNUPDT,KSDFT)

        Call mma_deallocate(VT)
        Call mma_deallocate(VL)
        Call mma_deallocate(XQN)
        Call mma_deallocate(SCR)
        Call mma_deallocate(V1)
        Call mma_deallocate(V2)
      ENDIF

! Rotation of orbitals with exp(x) where x is obtained from
! the super-CI coefficients, with a Quasi Newton update (NQUNE=1)


! CMO:  before - old MO's           after - new MO's
! CMON: intermediate storage for new MO's (moved to CMO in ORTHO)
! X2:  work area, also in ORTHO (AO overlap matrix)
! Scr: WORK AREA

      CALL mma_allocate(CMON,NTOT2,Label='CMON')
      CALL mma_allocate(XMAT,NO2M,Label='XMAT')
      CALL mma_allocate(X2,NTOT1,Label='X2')
      CALL mma_allocate(Scr,NO2M,Label='SCR')

      CALL ROTORB(CMO,CMON,CSX(LCSXI),XMAT,X2,SCR,THMAX,FA)

      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)
        Write(LF,*)'FI+FA in SXCTL after Unitary transform in ROTORB'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      END IF
      Call mma_deallocate(CMON)
      Call mma_deallocate(XMAT)
      Call mma_deallocate(X2)
      Call mma_deallocate(SCR)
      Call mma_deallocate(SXN)
      Call mma_deallocate(DIA)
      Call mma_deallocate(CSX)

      IDISK=IADR15(2)
      CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
      CALL DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#ifdef _HDF5_
      call mh5_put_dset(wfn_mocoef,CMO)
      call mh5_put_dset(wfn_occnum,OCC)
#endif
      CALL TIMING(CPTS,CPES,TIOS,TIOES)

9990  CONTINUE
      Call mma_deallocate(BM)

      END Subroutine SXCtl
