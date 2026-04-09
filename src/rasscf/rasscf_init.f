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
* Copyright (C) Per Ake Malmqvist                                      *
************************************************************************
*  RasScf_Init
*
*> @brief
*>   Initialize variables in commons, and set default values.
*>   Determine whether orbital files should be read, etc.
*> @author  P. &Aring;. Malmqvist
*>
*> @details
*> Sets values in the modules timers, rasscf_global and general_data.
************************************************************************

      Subroutine RasScf_Init()
      Use Fock_util_global, only: ALGO, Deco, DensityCheck, dmpk,
     &                            DoCholesky, DoLocK, Estimate, Nscreen,
     &                            Update
      use casvb_global, only: ifvb
      use Cholesky, only: ChFracMem, timings
      use CMS, only: iCMSOpt,CMSGiveOpt
      use UnixInfo, only: SuperName
      use gas_data, only: NGAS, NGSSH, IGSOCCX
      use timers, only: TimeAoMo, TimeCIOpt, TimeDavid, TimeDens,
     &                  TimeFock, TimeHCSCE, TimeHDiag, TimeHSel,
     &                  TimeInput, TimeOrb, TimePage, TimeRelax,
     &                  TimeSigma, TimeTotal, TimeTrans, TimeWfn
      use lucia_data, only: TDENSI, TSIGMA
      use rasscf_global, only: IROOT, CMSStartMat, CMSThreshold,
     &                         CORESHIFT, Ener, ExFac, hRoots,
     &                         iAlphaBeta, ICICH, ICICP, iCIonly,
     &                         ICIRST, ICMSIterMax, ICMSIterMin, iCMSP,
     &                         iExpand, IfCRPR, IfOrde, InOCalc,
     &                         iOrbOnly, iOrbTyp, iOrdeM, iPCMRoot,
     &                         iPhName, iPT2, iRLXRoot, iRoot, irotPsi,
     &                         iSave_Exp, iSPDen, iSupSM, itCore,
     &                         ITMAX, ITRIM, iXMSP, KSDFT,
     &                         kTight, LowMS, LRoots, LvShft, MaxIt,
     &                         MaxJT, MaxOrbOut, n_keep, NewFock,
     &                         NonEq, NQUNE, NROOTS, OutFmt1, OutFmt2,
     &                         PreThr, ProThr, PrwThr, Purify, QNSTEP,
     &                         QNUPDT, RFPert, SXSel, ThFact, Thre,
     &                         ThrEn, ThrSX, TMin, Weight, Title,
     &                         ixSym, iTri, ThrTE
      use output_ras, only: LF
      use general_data, only: SXDAMP,NSEL,LOWDIN_ON,ISPIN,STSYM,NACTEL,
     &                        NHOLE1,NELEC3,NALTER,STARTORBFILE,NASH,
     &                        NBAS,NDEL,NFRO,NISH,NRS1,NRS2,NRS3,NRS3,
     &                        NSSH
      use spinfo, only: I_ELIMINATE_GAS_MOLCAS,ISPEED
      use Molcas, only: MxOrb, MxRoot, MxSym
      use RASDim, only: MxCIIt, MxIter, MxSXIt

      Implicit None
      Integer IPRGLB_IN, IPRLOC_IN(7)
* What to do with Cholesky stuff?
      Logical, External :: Is_First_Iter
      Integer I
      Integer, External:: iPrintLevel

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* How was the program called?
*PAM 2009 Someone has put a number of possibilities here. Let it stand for now.
      IfVB=0
      If (SuperName(1:6).eq.'rasscf') Then
         ICIRST=0
*        For geometry optimizations use the old CI coefficients.
         If (.Not.Is_First_Iter()) ICIRST=1
      Else If (SuperName(1:5).eq.'casvb') Then
         IfVB=2
         ICIRST=0
      Else If (SuperName(1:6).eq.'loprop') Then
C        ICIRST=1 ! to be activated!
         ICIRST=0
      Else If (SuperName(1:11).eq.'last_energy') Then
         ICIRST=1
      Else If (SuperName(1:18).eq.'numerical_gradient') Then
         ICIRST=1
      Else
         ICIRST=0
      End If

* Initialize print levels: Module output_ras
* Global logical unit numbers for standard output
      LF=6
* Externally set default print level control. Should the program be silent?
      IPRGLB_IN=iPrintLevel(-1)
      DO I=1,7
       IPRLOC_IN(I)=IPRGLB_IN
      END DO
* Set print levels, and adjust them if needed:
      call setprlev(LF,IPRGLB_IN,IPRLOC_IN)

* Cholesky-related settings:
      Call DecideOnCholesky(DoCholesky)
      ALGO  = 1
      DensityCheck=.false.
      Deco=.true.
      timings=.false.
      DoLock=.true.
      Nscreen=10
      dmpk=1.0d-1
      Update=.true.
      Estimate=.false.
*
#ifdef _MOLCAS_MPP_
      ChFracMem=0.3d0
#else
      ChFracMem=0.0d0
#endif



      OutFmt1='DEFAULT '
      OutFmt2='DEFAULT '

* Max nr of state-specific orbital files printed:
      MAXORBOUT=100
* Default title line:
      TITLE(1)='(No title given)'
*
* iteration control
*
* maximum number of RASSCF iterations
      MAXIT=mxIter
* max number of super-CI iterations
      ITMAX=mxSxIt
* max number of iterations in Davidson diagonalization
      MAXJT=MXCIIT-2
* threshold for change in RASSCF energy
      THRE=1.D-08
*tbp, may 2013: no thre modification with Cholesky
*tbp  If (DoCholesky) then
*tbp     Call Get_dScalar('Cholesky Threshold',ThrCom)
*tbp     THRE = Max(THRE,ThrCom)
*tbp  EndIf
* threshold for max orbital rotation
*PAM2010       THRTE=1.D-04
* PAM2010: Note: This is *not* a threshold that keeps rotation down
* between iterations in order to ensure proper function of the
* optimization -- it was intended as one of the thresholds that
* determine when the calculation has converged! As such, it is
* irrelevant! The relevant threshold is the max BLB.
      THRTE=1.D-01
* threshold for max BLB matrix element
* Note: If one changes the following value, please change it in
* fock_util/cho_LK_rassi.f and fock_util/cho_LK_rassi_x.f for consistency.
      THRSX=1.D-04
* Default damping in the SXCI orbital optimization
      SXDAMP=0.0002D0
* Default thresholds used to determine convergence in CI
      THREN=1.0D-04
      THFACT=1.0D-03
* PAM 2017, Additional shift for douby occupied core states
* in order to compute core hole states. The core orbital is
* specified as one particular orbital in the input orbital set.
      CORESHIFT=0.0D0
      ITCORE=0
      IFCRPR=.false.
* PAM 2009, new default value for LVSHFT
* level shift parameter
      LVSHFT=0.5D00
* Quasi Newton update of the rotation matrix
      NQUNE=2
* only the CI calculation will be performed if iCIonly=1
      iCIonly=0
* only the orbitals from a JobIph to RasOrb if iOrbOnly=1
      iOrbOnly=0
* Default orbital type for RasOrb: Average orbitals
      iOrbTyp=1
* Root selection in the SXCI orbital optimization step.
* Values: LOWEST or HOMING.
      SXSEL='LOWEST  '
* Choose to only expand or generate information for CI-vectors if INOCALC = 1
      INOCALC = 0
* Save information on CI expansion if ISAVE_EXP = 1
      ISAVE_EXP = 0
* Expand a smaller CI vector in a larger one if IEXPAND = 1
      IEXPAND = 0
*
* wave function control bits
*
* new fock operator
      NewFock=1
* State used in response calculation
      iPCMROOT=1
* State to alaska
      iRLXROOT=0
* number of roots required in CI
      NROOTS=1
* number of roots actually used in CI-DAVIDSON
      LROOTS=1
* sequence numbers for roots in CI counted from
* lowest energy.
      Call iCopy(mxRoot,[0],0,iRoot,1)
      IROOT(1)=1
* weights used for average energy calculations
      Call dCopy_(mxRoot,[0.0D0],0,WEIGHT,1)
      WEIGHT(1)=1.0D0
* iteration energies
      Call dCopy_(mxRoot*(mxIter+2),[0.0D0],0,ENER,1)
*
      ICICH=0
* if flag is active (ICICH=1) CI roots will be selected
*             by maximum overlap with input CI function
*             ICI(NROOTS,NREF)    CSF number for each root
*             CCI(NROOTS,NREF)    corresponding CI coefficient
*             maximum number is five csf's.
*
      ISUPSM=0
* make no use of supersymmetry
      I_ELIMINATE_GAS_MOLCAS = 0
* Highly excited states are not default
      hRoots=0
* No hidden roots by default
      n_keep=0
* Number of kept vectors in Davidson chosen in ini_david by default
      IORDEM=0
* (SVC) do not force any ordering options
      IFORDE=1
* (SVC) use ordering of orbitals
* start CI Davidson with unit guess for CI vector
* use restart option if numerical gradients are computed.
      PRWTHR = 0.05D0
* threshold for printout of CI wave function

      PROTHR=-1.0D0
* occupation threshold for printout of orbitals
* ( The negative value serves to show if no user selection was made)
      PRETHR=999999.0d0
* energy threshold for printout of orbitals

      ICICP=0
* no CI coupling (not active in this version)
      NSEL=200
* Default value for explicit Hamiltonian
*
      TMIN=0.0D0
      QNSTEP='SX'
      QNUPDT=' NO'
*     Default value for tight parameter
      KTIGHT=0
*
* Default value for type of CASSCF (used for DFT)
*
      KSDFT='SCF'
      ExFac=1.0D0
** Default orthonormalization of CMOs to be with
** Gram-Schmidt
*      Lowdin_ON=.False.
* PAM Jan 12 2010, on request, Lowdin ON has been made the default.
      Lowdin_ON=.true.
*
* default for spin projection
      LOWMS=0
* default spin value (singlet)
      ISPIN=1
* default symmetry
      STSYM=1
* default number of active electrons
      NACTEL=0
* default maximum number of holes in RAS1
      NHOLE1=0
* default maximum number of electrons in RAS3
      NELEC3=0
* This run will not be the start for a CASPT2 calculation
      IPT2=0
* This key will activate pertubational reaction field
* calculations.
      RFpert=.false.
* Do compute the spin density matrix
      ISPDEN=1
* These keys will activate the calculation of the high
* frequency contribution to the reaction field
* ???
* This key controls if a non-equilibrium reaction field
* calculation is performed.
      NonEq=.False.
* This initializes nr of input orbital swaps requested:
      NAlter=0
* set default values for orbitals
*
C
      DO I=1,mxSym
        NFRO(I)=0
        NISH(I)=0
        NASH(I)=0
        NRS1(I)=0
        NRS2(I)=0
        NRS3(I)=0
        NSSH(I)=0
        NDEL(I)=0
        NBAS(I)=0
      END DO
* initialize occupation numbers for GAS
*
      NGAS=3
      NGSSH=0
      IGSOCCX=0
      DO I=1,mxOrb
        IXSYM(I)=0
      END DO
      PURIFY='NO'
*
*     Auxiliary vector ITRI(I)=I*(I-1)/2
*
      ITRI(1)=0
      DO I=2,ITRIM
       ITRI(I)=ITRI(I-1)+I-1
      END DO

* Initial guess for jobiph name to use:
      IPHNAME='JOBIPH'
* Initial guess for starting orbital file:
      StartOrbFile='INPORB'
* Initialize alpha or beta orbitals (none):
      iAlphaBeta=0
*
* Initialize speed options (turn everything that's working on)
*
      Do i = 1,size(iSpeed)
         if (i .le. 2) Then
            iSpeed(i) = 1
         else
C The rest is at the present time just to allow testing
            iSpeed(i) = 0
         end if
      End Do
*
      TimeTotal  = 0.0d0
      TimeInput  = 0.0d0
      TimeWfn    = 0.0d0
      TimeDens   = 0.0d0
      TimeSigma  = 0.0d0
      TimeHSel   = 0.0d0
      TimeHDiag  = 0.0d0
      TimeFock   = 0.0d0
      TimeAoMo   = 0.0d0
      TimeTrans  = 0.0d0
      TimeCIOpt  = 0.0d0
      TimeOrb    = 0.0d0
      TimeDavid  = 0.0d0
      TimePage   = 0.0d0
      TimeHCSCE  = 0.0d0
      TimeRelax  = 0.0d0
*
CSVC: lucia timers
      tsigma(:) = 0.0d0
      tdensi(:) = 0.0d0
*
C state rotation
      iRotPsi=0
      iXMSP=0
      iCMSP=0
      ICMSIterMax=100
      ICMSIterMin=5
      CMSThreshold=1.0d-8
      CMSStartMat='XMS'
      iCMSOpt=1
      CMSGiveOpt=.false.
      RETURN
      END Subroutine RasScf_Init
