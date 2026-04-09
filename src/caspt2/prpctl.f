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
      SUBROUTINE PRPCTL(MODE,UEFF,U0)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, Half, One, Five
      USE PT2WFN, only: PT2WFN_DENSSTORE
      use caspt2_global, only:iPrGlb
      use OneDat, only: sNoNuc, sNoOri
      use caspt2_global, only: do_grad,do_nac,iRoot1,iRoot2,SLag,
     &                         DPT2_tot,DPT2C_tot
      use caspt2_global, only: CMO, CMO_Internal, CMOPT2, TORB, NCMO,
     &                       LISTS
      use caspt2_global, only: LUONEM
      use PrintLevel, only: USUAL, VERBOSE
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use EQSOLV, only: IVECX, NLSTOT
      use caspt2_module, only: NSTATE, IFMSCOUP, IFPROP, irlxroot,
     &                         ISCF, JSTATE, BNAME, NASHT, NBAST, NCONF,
     &                         NSYM, OUTFMT, PRORB, THRENE, THROCC,
     &                         NORB, NBAS, NISH, NASH, IAD1M, NFRO,
     &                         NRAS1, NRAS2, NRAS3, MSTATE, NDEL,
     &                         Energy, MSTATE

      IMPLICIT None

      integer(kind=iwp), intent(in):: Mode
      real(kind=wp), intent(in):: UEFF(NSTATE,NSTATE),U0(*)

      Logical(kind=iwp) FullMlk,lSave,Do_ESPF
      Character(Len=8) Label
      Character(Len=128) FILENAME,MDNAME
      Character(Len=80) Note
      integer(kind=iwp) IndType(56)
      real(kind=wp) Dummy(2), DUM(1), SCAL
      integer(kind=iwp) NFROSAV(NSYM),NORBSAV(NSYM)
      real(kind=wp), allocatable :: DMAT(:),CI1(:),CI2(:),SGM(:),
     &                              TG1(:,:),CNAT(:),OCC(:),Scr(:)
      integer(kind=iwp) I, iComp, IDISK, IDMAT, IDMOFF, IERR, II, II2,
     &                  IJ, IJ2, IndT, iOpt, iRc, iShift, ISTATE,
     &                  iSyLbl, ISYM, iUHF, KSTATE, LUTMP, N, nDens,
     &                  NDMAT, NO, NOCC
      integer(kind=iwp), external:: IsFreeUnit

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par() .AND. IPRGLB.GE.USUAL .AND. .not.do_grad) THEN
        WRITE(u6,'(1X,A)') ' ====================================='
        WRITE(u6,'(1X,A)') ' CASPT2 properties were requested, but'
        WRITE(u6,'(1X,A)') ' these are not efficiently implemented'
        WRITE(u6,'(1X,A)') ' in parallel. If you do not need them,'
        WRITE(u6,'(1X,A)') ' not using the PROP keyword could lead'
        WRITE(u6,'(1X,A)') ' to a significant speed up.'
        WRITE(u6,'(1X,A)') ' ====================================='
      END IF
#else
      call unused_logical(do_grad)
#endif

* PAM2008 When this subroutine is called, theMSTATE calculation has been done
* for the (individual) state nr JSTATE in 1,2,..,NSTATE.
* The corresponding CI-root from rasscf, the root state for this PT2,
* is number MSTATE(JSTATE) on the input JOBIPH file.
* JSTATE,NSTATE and MSTATE() are in the caspt2_module.F90 file.
      IERR=0
      IF(NSTATE.GT.1) THEN
       N=MSTATE(JSTATE)
       IF(N.LE.0 .or. N.GT.999) THEN
        WRITE(u6,*)' Subroutine PRPCTL fails -- It seems to get data'
        WRITE(u6,*)' computed for a root nr ',N
        WRITE(u6,*)' which is surely wrong.'
        WRITE(u6,*)' PRPCTL gives up, there will be no calculations'
        WRITE(u6,*)' done of orbitals, properties, etc for this state.'
        WRITE(u6,*)' This was state nr JSTATE=',JSTATE
        WRITE(u6,*)' in the MS-CASPT2 calculation.'
        IERR=1
       END IF
      END IF
      IF(IERR.GT.0) Return

      FullMlk=.True.
      IF (.NOT.PRORB ) FullMlk=.False.

      IF ( IPRGLB.GE.USUAL ) THEN
       WRITE(u6,'(20A4)')('----',I=1,20)
      END IF

C Compute density matrix, output orbitals, and properties.

C Compute density matrix of CASPT2 wave function, in MO basis,
C to produce output orbitals.
C This density matrix may be approximated in several ways, see DENS.
      NDMAT=0
      NOCC=0
      IF (MODE.EQ.0) THEN
        !! Density matrix for each state (single-state)
        DO ISYM=1,NSYM
          NO=NORB(ISYM)
          NDMAT=NDMAT+(NO**2+NO)/2
          NOCC=NOCC+NBAS(ISYM)
        END DO
        call mma_allocate(DMAT,NDMAT,Label='DMAT')
        CALL DCOPY_(NDMAT,[Zero],0,DMAT,1)
        CALL mma_allocate(LISTS,NLSTOT,LABEL='LISTS')
        CALL MKLIST(LISTS)
        CALL DENS(IVECX,DMAT,UEFF,U0)
        CALL mma_deallocate(LISTS)
      ELSE
        !! Density matrix for the target adiabatic state
        !! MODE = 1 is called after gradient stuff only for MS, so
        !! the density computed here is what we call a correct unrelaxed
        !! correlated CASPT2 density
        !! Called after gradient calculations, only, from GrdCls
        DO ISYM=1,NSYM
          NO=NBAS(ISYM)
          NDMAT=NDMAT+(NO**2+NO)/2
          NOCC=NOCC+NBAS(ISYM)
        END DO
        call mma_allocate(DMAT,NDMAT,Label='DMAT')
        CALL DCOPY_(NDMAT,[Zero],0,DMAT,1)
        !! Copy the unrelaxed density matrix to triangular
        !! The basis of DPT2_tot is natural (CASSCF)
        IDMAT = 0
        IDMOFF = 0
        DO ISYM=1,NSYM
          NO=NBAS(ISYM)
          DO II = 1, NO
            DO IJ = 1, II
              !! second-order (DPT2) and first-order (DPT2C)
              DMAT(1+IDMAT) = DPT2_TOT(IDMOFF+II+NO*(IJ-1))
     *                      + DPT2C_TOT(IDMOFF+II+NO*(IJ-1))*0.25d+00
              IF (.NOT.DO_NAC) THEN
                !! Add the reference density matrix (inactive)
                IF (II.EQ.IJ .and. II.LE.NFRO(ISYM)+NISH(ISYM))
     *            DMAT(1+IDMAT) = DMAT(1+IDMAT) + 2.0D+00
              END IF
              IDMAT = IDMAT + 1
            END DO
          END DO
          IDMOFF = IDMOFF + NO*NO
        END DO
        !! Add the reference density matrix (active)
        call mma_allocate(CI1,NCONF,Label='CI1')
        call mma_allocate(CI2,NCONF,Label='CI2')
        call mma_allocate(SGM,NCONF,Label='SGM')
        call mma_allocate(TG1,NASHT,NASHT,Label='TG1')
        DO ISTATE = 1, NSTATE
          IF (ISCF.NE.0) THEN
            CI1(1)=One
          ELSE
            CALL LOADCI_XMS('N',1,CI1,ISTATE,U0)
          END IF
          DO KSTATE = 1, NSTATE
            SCAL = SLag(ISTATE,KSTATE)
            IF (.NOT.DO_NAC) THEN
              IF (ISTATE.EQ.IROOT1.AND.KSTATE.EQ.IROOT2)
     *          SCAL = SCAL + One
            END IF
            IF (ABS(SCAL).LE.1.0D-09) CYCLE
            IF (ISCF.NE.0) THEN
              CI2(1)=One
            ELSE
              CALL LOADCI_XMS('N',1,CI2,KSTATE,U0)
            END IF
            Call Dens1T_RPT2(CI1,CI2,
     *                       SGM,TG1,nAshT)
            CALL DSCAL_(NASHT**2,SCAL,TG1,1)
            DO II = 1, NASH(1)
              II2 = II+NFRO(1)+NISH(1)
              DO IJ = 1, II
                IJ2 = IJ+NFRO(1)+NISH(1)
                DMAT(II2*(II2-1)/2+IJ2)
     *            = DMAT(II2*(II2-1)/2+IJ2)
     *            + TG1(II,IJ)*Half + TG1(IJ,II)*Half
              END DO
            END DO
          END DO
        END DO
        call mma_deallocate(CI1)
        call mma_deallocate(CI2)
        call mma_deallocate(SGM)
        call mma_deallocate(TG1)

        !! The density matrix above has frozen orbitals.
        !! Accordingly, the natural orbitals should be generated
        !! by diagonalizing all orbitals (including frozen), so
        !! make the number of frozen orbitals zero for the moment.
        CALL ICOPY(NSYM,NFRO,1,NFROSAV,1)
        CALL ICOPY(NSYM,NORB,1,NORBSAV,1)
        CALL ICOPY(NSYM,[0],0,NFRO,1)
        CALL ICOPY(NSYM,NBAS,1,NORB,1)
      END IF

C Compute natural orbitals of CASPT2 wave function.
      CALL mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
      CMO=>CMO_Internal
      CMO(:)=CMOPT2(:)
      IF (MODE.EQ.1) THEN
        !! Read the CASSCF orbital (really?)
        IDISK=IAD1M(1)
        call ddafile(LUONEM,2,CMO,NCMO,IDISK)
      END IF
      call mma_allocate(CNAT,NCMO,Label='CNAT')
      call mma_allocate(OCC,NOCC,Label='OCC')
      CALL NATORB_CASPT2(DMAT,CMO,OCC,CNAT)
      CALL mma_deallocate(CMO_Internal)
      nullify(CMO)
C Backtransform density matrix to original MO basis before storing
      CALL TRANSFOCK(TORB,SIZE(TORB),DMAT,NDMAT,-1)
      CALL PT2WFN_DENSSTORE(DMAT,NDMAT)
      call mma_deallocate(DMAT)
      IF (MODE.EQ.1) THEN
        !! Restore with frozen orbitals
        CALL ICOPY(NSYM,NFROSAV,1,NFRO,1)
        CALL ICOPY(NSYM,NORBSAV,1,NORB,1)
      END IF

      IF (IFMSCOUP.AND.MODE.EQ.0.AND..NOT.IFPROP) THEN
        !! Do not show properties, if PROP calculation is not requested
        call mma_deallocate(CNAT)
        call mma_deallocate(OCC)
        RETURN
      END IF

C Write natural orbitals as standard orbital file on PT2ORB
* PAM2008: Separate PT2ORB files for each state:
      FILENAME='PT2ORB'
      MDNAME='MD_PT2'
      IF(NSTATE.GT.1.AND.MODE.EQ.0) THEN
       FILENAME='PT2ORB.x'
       MDNAME='MD_PT2.x'
       N=MSTATE(JSTATE)
       IF(N.LE.9) THEN
        WRITE(FILENAME(8:8),'(I1)') N
        WRITE(MDNAME(8:8),'(I1)') N
       ELSE IF(N.LE.99) THEN
        WRITE(FILENAME(8:9),'(I2)') N
        WRITE(MDNAME(8:9),'(I2)') N
       ELSE IF(N.LE.999) THEN
        WRITE(FILENAME(8:10),'(I3)') N
        WRITE(MDNAME(8:10),'(I3)') N
       END IF
      END IF
* PAM2008: For MS-CASPT2 with more than one state, orbital files will
* now be numbered PT2ORB.1 ... PT2ORB.999
* depending on which CI-root of the rasscf calculation that was
* the root state of the PT2.
      LUTMP=19
      LUTMP=ISFREEUNIT(LUTMP)
* PAM 2008: Add typeindex information:
*----------------------------------------------------------------------*
*     Make typeindex information                                       *
*----------------------------------------------------------------------*
      iShift=0
      DO ISYM=1,NSYM
        IndT=0
        IndType(1+iShift)= NFRO(ISYM)
        IndT=IndT+NFRO(ISYM)
        IndType(2+iShift)= NISH(ISYM)
        IndT=IndT+NISH(ISYM)
        IndType(3+iShift)= NRAS1(ISYM)
        IndT=IndT+NRAS1(ISYM)
        IndType(4+iShift)= NRAS2(ISYM)
        IndT=IndT+NRAS2(ISYM)
        IndType(5+iShift)= NRAS3(ISYM)
        IndT=IndT+NRAS3(ISYM)
        IndType(7+iShift)= NDEL(ISYM)
        IndT=IndT+NDEL(ISYM)
        IndType(6+iShift)= NBAS(ISYM)-IndT
        iShift=iShift+7
      EndDo
      If (NSTATE.GT.1) THEN
        Write(Note,'(A41,I3,A3,f22.12)')
     &   '* CASPT2 natural orbitals for root number',N,
     &   ' E=',Energy(JSTATE)
      Else
        Note='* CASPT2 natural orbitals'
      End If

      CALL WRVEC(FILENAME,LUTMP,'COI',NSYM,NBAS,NBAS,
     &  CNAT, OCC,Dummy  ,IndType,Note)
      iUHF=0
      Call Molden_Interface(iUHF,FILENAME,MDNAME)

C Write natural orbitals to standard output.
      IF ( IPRGLB.GE.VERBOSE) THEN
       WRITE(u6,*)
       WRITE(u6,'(A)')'  The CASPT2 orbitals are computed as natural '//
     &          'orbitals of a density matrix'
       WRITE(u6,'(A)')'  defined as:'
       WRITE(u6,'(A)')'   D = (D0 + D1 + D2)/<PSI|PSI>'
       WRITE(u6,'(A)')' where D0..D2 are zeroth..2nd order'//
     &                ' contributions'
       WRITE(u6,'(A)')' and |PSI> is the total wave function.'
       WRITE(u6,'(A)')' A new RasOrb file named PT2ORB is prepared.'
       IF (PRORB) THEN
         IF ( OUTFMT.EQ.'LONG    ' ) THEN
           THRENE=2.0d0**31
           THROCC=-2.0d0**31
         ELSE IF ( OUTFMT.EQ.'DEFAULT ' ) THEN
           THRENE=Five
           THROCC=5.0d-04
         END IF
         CALL PRIMO('Output orbitals from CASPT2',
     &           .TRUE.,.FALSE.,THROCC,THRENE,NSYM,NBAS,
     &            NBAS,BNAME,DUM,OCC,CNAT,-1)
       END IF
      END IF

* compute Mulliken's orbital populations

      IF ( IPRGLB.GE.USUAL ) THEN
        WRITE(u6,*)
        WRITE(u6,*)
        WRITE(u6,'(6X,A)') 'Mulliken population Analysis:'
        WRITE(u6,'(6X,A)') '-----------------------------'

        call mma_allocate(Scr,NBAST**2,Label='Scr')
        iRc=-1
        iOpt=ibset(ibset(0,sNoOri),sNoNuc)
        iComp=1
        iSyLbl=1
        Label='Mltpl  0'
        Call RdOne(iRc,iOpt,Label,iComp,Scr,iSyLbl)
        If ( iRc.eq.0 ) then
           lSave = MSTATE(JSTATE) .eq. irlxroot
           Call Charge(nSym,nBas,bName,
     &               CNAT,OCC,Scr,2,FullMlk,lSave)
        End If
        call mma_deallocate(Scr)
      END IF

* compute one-electron properties

      IF ( IPRGLB.GE.USUAL ) THEN
        WRITE(u6,*)
        WRITE(u6,'(6X,A)') 'Expectation values of various properties:'
        WRITE(u6,'(6X,A)') '-----------------------------------------'
      END IF

      nDens=0
      Do i = 1, nSym
         nDens=nDens+nBas(i)*(nBas(i)+1)/2
      End Do
      call mma_allocate(Scr,NDENS,Label='Scr')
*
      Call DOne_Caspt2(CNAT,OCC,Scr)
      Call Put_dArray('D1ao',Scr,nDens)
*
      Note='Temporary orbital file used by prpt.'
      if (mode.eq.1) Note='var'
      LuTmp=50
      LuTmp=IsFreeUnit(LuTmp)
      Call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,
     &            CNAT,OCC,Dummy,IndType,Note)
      Call Prpt()
cnf
*
*------- ESPF analysis
      Call DecideOnESPF(Do_ESPF)
      lSave = MSTATE(JSTATE) .eq. irlxroot
      If (Do_ESPF) Call espf_analysis(lSave)
cnf
*
*---- On return from PrPt the 1-particle matrix is stored
*     in the beginning of the scratch array.
*
      call mma_deallocate(Scr)
      call mma_deallocate(CNAT)
      call mma_deallocate(OCC)

      END SUBROUTINE PRPCTL
