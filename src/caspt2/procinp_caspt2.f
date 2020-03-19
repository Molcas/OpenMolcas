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
      Subroutine ProcInp_Caspt2
CSVC: process CASPT2 input based on the data in the input table, and
C initialize global common-block variables appropriately.
      Use InputData
      Implicit None
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "ofembed.fh"
#include "para_info.fh"

      Integer iDummy

* Number of non-valence orbitals per symmetry
      Integer nCore(mxSym)
      Integer nDiff, NFI, NSD
* Geometry-determining root
      Logical Is_iRlxRoot_Set
* Cholesky stuff
      Logical REORD,DECO,timings
      Integer Algo
      COMMON /CHORAS  / REORD,DECO,ALGO
      COMMON /CHOTIME / timings

      Integer I, J, M, N
      Integer ISYM
* State selection
      Integer iGroup, IOFF
* Error condition
      Integer IERR

#include "chocaspt2.fh"

      CALL QENTER('READIN')

* Hzero and Focktype are merged together into Hzero. We keep the
* variable Focktype not to break the input keyword which is documented
* in the manual. However, eventually we will have to keep only Hzero
* and remove Focktype.
      Hzero = input%Hzero
      if (Hzero.ne.'STANDARD'.and.Hzero.ne.'CUSTOM') then
        call WarningMessage(2,
     &   'invalid 0th-order Hamiltonian: '//TRIM(Hzero))
        call Quit_OnUserError
      end if

* Choose Focktype, reset IPEA shift to 0 for non-standard fock matrices
      Focktype = input%Focktype
      if (Focktype.ne.'STANDARD') then
        if (IfChol) then
          Call WarningMessage(2,'Requested FOCKtype not possible.')
          WRITE(6,*)'Calculations using Cholesky vectors can only'
          WRITE(6,*)'be used with the standard FOCKtype!'
          Call Quit_OnUserError
        end if
* If both Hzero and Focktype are not standard, quit
        if (Hzero.ne.'STANDARD') then
          Call WarningMessage(2,'Requested combination of FOCKtype'//
     &                          ' and HZERo not possible.')
          Call Quit_OnUserError
        end if
* IPEA different from zero only for standard Focktype
        if (BSHIFT.gt.0.0d0.or.BSHIFT.lt.0.0d0) then
          BSHIFT = 0.0d0
          if (IPRGLB.ge.TERSE) then
            call WarningMessage(1,'IPEA shift reset to zero!')
          end if
        end if
      else
* user-specified IPEA shift or not?
        if (input%IPEA) then
          BSHIFT = input%BSHIFT
        else
* Set default IPEA to 0.25 Eh
          BSHIFT = 0.25d0
        end if
      end if

* Copy over to Hzero the content of Focktype, if Hzero is not CUSTOM
      if (Hzero.ne.'CUSTOM') then
        Hzero = Focktype
      end if

* print warnings if deviating from the default
      if (Hzero.ne.'STANDARD') then
        call warningmessage(1,'User-modified 0th-order Hamiltonian!')
      end if

* real/imaginary shifts
      SHIFT = Input % Shift
      SHIFTI = Input % ShiftI

* RHS algorithm selection
#ifdef _MOLCAS_MPP_
#  ifdef _GA_
C     The RHS on-demand algorithm doesn't handle serial calculations
C     because it's not adapted for use with regular Work arrays, only
C     global arrays, and needs to be switched off (using rhsall instead)
      RHSDIRECT=(Is_Real_Par().AND.Input%RHSD)
#  else
C     Without the Global Arrays library, we can't use the RHSALL2
C     and ADDRHS algorithms in parallel. Here we force the use of
C     RHS on-demand instead, depending on if the calculation is
C     really parallel or not.
      RHSDIRECT=Is_Real_Par()
#  endif
#else
      RHSDIRECT=.False.
#endif

* Cholesky: set defaults if it was not called during input
      If (.NOT.(Input % ChoI .or. Input % Chol)) Then
        Call Cho_caspt2_rdInp(.True.,iDummy)
      End If

*---  Initialize
      IDCIEX=0
      IEOF1M=0
      DO I=1,64
        IAD1M(I)=-1
      END DO
      RFpert = Input % RFPert

      NTIT=0
      OUTFMT='DEFAULT'
      G1SECIN=.FALSE.
      PRORB=.TRUE.
      PRSD=.FALSE.
      NCASES=13

      JMS = Input % JMS
      NLYROOT = Input % OnlyRoot
      NLYGROUP = 0

      DoCumulant = Input % DoCumulant

************************************************************************
*
* Root selection
*
************************************************************************
*SB: create a flat array, putting (X)MS states with increasing
* group number. Do not allow to run both MS and XMS in the same
* calculation, it can lead to catastrophic results.
* For MS, put one state per group, for XMS put all states in
* a single group.
* Example: lets say the input is:
* MULTistate = 3 3 4 5
*
* Then the arrays that are created will be:
* nstate = 3
* mstate(nstate): 3 4 5
* nGroup = 3
* nGroupState(nGroup): 1 1 1
*
* On the other hand, if the input is
* XMULtistate = 4 1 2 3 4
*
* Then the arrays that are created will be:
* nstate = 4
* mstate(nstate): 1 2 3 4
* nGroup = 1
* nGroupState(nGroup): 4
      NSTATE = 0
      MSTATE = 0
      NGROUP = 0
      NGROUPSTATE = 0
* This is the case for MS-CASPT2 and DW-CASPT2
      If (Input%MULT) Then
        If (Input%XMUL) Then
          Call WarningMessage(2,'Keyword MULTistate cannot be used '//
     &                          'together with keyword XMULtistate.')
          Call Quit_OnUserError
        End If
* Either the states were specified manually or the keyword "all"
* was used, so first we check the keyword all
        If (Input%AllMult) Then
          NSTATE = NROOTS
          MSTATE = IROOT
          NGROUP = NSTATE
          NGROUPSTATE(1:NGROUP)=1
        Else
* Save the states that need to be computed
          Do I=1,Input%nMultState
            MSTATE(I) = Input%MultGroup%State(I)
            NSTATE = NSTATE + 1
          End Do
          NGROUP = Input%nMultState
          NGROUPSTATE(1:NGROUP) = 1
        End If
      End If
      IOFF=NSTATE
* This is the case for XMS-CASPT2 and XDW-CASPT2
      if (Input%XMUL) then
        IFXMS = Input%XMUL
        if (Input%MULT) then
          call WarningMessage(2,'Keyword XMULtistate cannot be used '//
     &                          'together with keyword MULTistate.')
          call Quit_OnUserError
        end if
        IFSilPrRot = Input%SilentPrRot
        IFNOPT2=Input%IFNOPT2
* This is a XDW-CASPT2 calculation. It is actually more similar to
* a MS-CASPT2 one since we need to put one state per group and thus
* have as many groups as states. Nevertheless, it makes more sense
* from the user point of view to ask for it through XMUL and DWMS
        if (Input%DWMS) then
          if (Input%AllXMult) Then
            NSTATE = NROOTS
            MSTATE = IROOT
            NGROUP = NSTATE
            NGROUPSTATE(1:NGROUP)=1
          else
            do I=1,Input%nXMulState
              MSTATE(I) = Input%XMulGroup%State(I)
              NSTATE = NSTATE + 1
            end do
            NGROUP = Input%nXMulState
            NGROUPSTATE(1:NGROUP) = 1
          end if
* This is a XMS-CASPT2: one group with all the states
        else
          if (Input%AllXMult) Then
            NSTATE = NROOTS
            MSTATE = IROOT
            NGROUP=1
            NGROUPSTATE(1)=NSTATE
          else
            NGROUP = 1
            NGROUPSTATE(NGROUP) = Input%nXMulState
            do I=1,Input%nXMulState
              MSTATE(I) = Input%XMulGroup%State(I)
              NSTATE = NSTATE + 1
            end do
          end if
        end if
      End If
* After parsing mult or xmult, check that no two equal states where
* given in the input
      If (Input%MULT.OR.Input%XMUL) Then
        Do I=1,NSTATE
          Do J=I+1,NSTATE
            If (MSTATE(I).EQ.MSTATE(J)) Then
             Call WarningMessage(2,'The same root cannot be used '//
     &                             'twice in MULT/XMULT blocks.')
              Call Quit_OnUserError
            End If
          End Do
        End Do
      End If
* The LROOt keyword specifies a single root to be used. It should not be
* used together with either MULT or XMUL keywords.
      If (Input%LROO) Then
        If (Input%MULT.OR.Input%XMUL) Then
          Call WarningMessage(2,'Keyword LROO cannot be used together'//
     &                          'with the MULT or XMUL keywords.')
          Call Quit_OnUserError
        End If
        NSTATE=1
        MSTATE(1)=Input%SingleRoot
        NGROUP=1
        NGROUPSTATE(1)=1
      End If
* If still nothing was selected we should default to compute all the
* roots that were part of the rasscf orbital optimization.
      IF (NSTATE.EQ.0) THEN
        NSTATE=NROOTS
        MSTATE=IROOT
        NGROUP=NSTATE
        NGROUPSTATE(1:NGROUP)=1
      END IF
* Find the group number for OnlyRoot
      If (NLYROOT.ne.0) Then
        IOFF=0
        Do IGROUP=1,NGROUP
          Do I=1,NGROUPSTATE(IGROUP)
            If (IOFF+I.eq.NLYROOT) NLYGROUP=IGROUP
          End Do
          IOFF=IOFF+NGROUPSTATE(IGROUP)
        End Do
      End If
* Finally, some sanity checks.
      IF(NSTATE.LE.0.OR.NSTATE.GT.MXROOT) Then
        Call WarningMessage(2,'Number of states is <0 or too large.')
        WRITE(6,'(a,i8)')' NSTATE = ',NSTATE
        WRITE(6,*)' Check usage of keywords MULT/XMUL.'
        Call Quit_OnUserError
      END IF
* setup root to state translation
      ROOT2STATE = 0
      DO I=1,NSTATE
        ROOT2STATE(MSTATE(I))=I
      END DO
*
* Relax root selection
*
      iRlxRoot=-1
      Call Qpg_iScalar('NumGradRoot',Is_iRlxRoot_Set)
      If (Is_iRlxRoot_Set) Then
          Call Get_iScalar('NumGradRoot',iRlxRoot)
      End If
      If (Input % RlxRoot.gt.0) iRlxRoot=Input%RlxRoot
      If (iRlxRoot.eq.-1) iRlxRoot=NSTATE
      If (iRlxRoot.gt.NSTATE) Then
        If (IPRGLB.GE.TERSE) then
         Call WarningMessage(1,'Too large iRlxRoot.')
         WRITE(6,*)' Reset to NSTATE=',NSTATE
        End If
        iRlxRoot=NSTATE
      End If

************************************************************************
*
* Determine number of Frozen/Deleted orbitals
*
************************************************************************
* The number of Frozen orbitals is initially read from the reference
* wavefunction. Here, we modify that number to be the larger of what
* was in the reference and the non-valence orbitals.
      Call Get_iArray('Non valence orbitals',nCore,nSym)
      Do iSym = 1, nSym
        If (nCore(iSym).gt.nFro(iSym)) Then
          nDiff=nCore(iSym)-nFro(iSym)
          nDiff=min(nDiff,nISh(iSym))
          nFro(iSym)=nFro(iSym)+nDiff
          nISh(iSym)=nISh(iSym)-nDiff
        End If
      End Do
* If a user specified the number of frozen orbitals explicitely, then
* that number overwrites the automatically chosen numbers, but it we
* warn the user and report what the computed number was.
      IF(Input % FROZ) THEN
        IF(IPRGLB.GE.TERSE) THEN
          Call WarningMessage(1,'User changed nr of frozen orbitals.')
        END IF
        IERR=0
        DO I=1,NSYM
          NFI=NFRO(I)+NISH(I)
          IF(NFI.LT.Input%nFro(I)) THEN
            Call WarningMessage(2,'Too many frozen orbitals!')
            Call Quit_OnUserError
          ELSE
            nFro(I)=Input%nFro(I)
          END IF
          NISH(I)=NFI-nFro(I)
        END DO
      ENDIF
* Set user-specified number of deleted orbitals.
      IF(Input % DELE) THEN
        IERR=0
        DO I=1,NSYM
          NSD=NSSH(I)+NDEL(I)
          IF(NSD.LT.Input%nDel(I)) THEN
            Call WarningMessage(2,'Too many deleted orbitals!')
            Call Quit_OnUserError
          ELSE
            NDEL(I)=Input%nDel(I)
          END IF
          NSSH(I)=NSD-nDel(I)
        ENDDO
      ENDIF

* NOT TESTED
#if 0
************************************************************************
*
* Orbital-free embedding.
*
************************************************************************
      Done_OFemb=.false.
      First_OFE =.true.
      ipFMaux = -666666
      If (Input % OFEmbedding) Then
        Done_OFemb=.true.
        write(6,*)
        write(6,*)  '  --------------------------------------'
        write(6,*)  '   Orbital-Free Embedding Calculation   '
        write(6,*)  '  --------------------------------------'
        write(6,*)
      End If
#endif
************************************************************************
*
* Determine what kind of calculations are needed/requested.
*
************************************************************************
* flags for enabling/disabling sections
      IFPROP = Input % Properties
      IFDENS = Input % DENS
      IFMIX  = .NOT.Input % NoMix
      IFMSCOUP = (Input % MULT .OR. Input % XMUL)
     &           .AND.(.NOT.Input % NoMult)
      IFDW = Input % DWMS
* Set exponent for DWMS
      if (IFDW) then
        zeta = Input % zeta
      end if

      IFEFOCK = Input % EFOC
      if (Input % EFOC) then
        if (.not.(IFXMS.and.IFDW)) then
          Call WarningMessage(2,'Keyword EFOCk can only be used'//
     &                          'in (X)DW-CASPT2 calculations.')
          Call Quit_OnUserError
        else
        end if
      end if

* Choice? of preprocessing route
      ORBIN='TRANSFOR'

*---  Expectation values of H0 needed or not? Spectral resolution?
      BMATRIX='YES     '
      BTRANS ='YES     '
      BSPECT ='YES     '
*---  Overlap matrix needed? Linear dependence removal needed?
      SMATRIX='YES     '
      SDECOM ='YES     '
C Presently, we always need the ON transformation/Lindep removal.
C Consistency of these demands:
      IF(BSPECT.NE.'NO      ') BTRANS ='YES     '
      IF(BTRANS.NE.'NO      ') BMATRIX='YES     '
      IF(BTRANS.NE.'NO      ') SDECOM ='YES     '
      IF(SDECOM.NE.'NO      ') SMATRIX='YES     '

************************************************************************
*
* Thresholds
*
************************************************************************
      If (Input%THRE) Then
        If (IPRGLB.GE.TERSE) THEN
          Call WarningMessage(1,
     &          'User modified linear dependency thresholds!')
        End If
      End If
      THRSHN=Input%THRSHN
      THRSHS=Input%THRSHS
      THRCONV=Input%THRCONV
      CITHR=Input%PrWF
      THRENE=5.0d+01
      THROCC=5.0d-04
      MAXIT=Input%MaxIter
      DNMTHR=Input%DNMTHR
      CMPTHR=Input%CMPTHR
      CNTTHR=Input%CNTTHR

*---  Create the symmetry multiplication table
      MUL(1,1)=1
      M=1
      DO N=1,3
        DO I=1,M
          DO J=1,M
            MUL(I+M,J)=M+MUL(I,J)
            MUL(I,J+M)=MUL(I+M,J)
            MUL(I+M,J+M)=MUL(I,J)
          END DO
        END DO
        M=2*M
      END DO
*
*---  Exit
      CALL QEXIT('PROC_INP')
      Return
      End
