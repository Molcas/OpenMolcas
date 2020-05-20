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
*> Sets values in common blocks in rasscf.fh, general.fh, timers.fh
************************************************************************
      Subroutine RasScf_Init_m()
      Implicit Real*8 (A-H,O-Z)
      External Get_SuperName
      Character*100 ProgName, Get_SuperName
#include "rasdim.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='RasScf_Init')
#include "rasscf.fh"
#include "casvb.fh"
#include "general.fh"
#include "gas.fh"
#include "timers.fh"
#include "lucia_ini.fh"
#include "orthonormalize.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer IPRGLB_IN, IPRLOC_IN(7)
* What to do with Cholesky stuff?
      Logical DoCholesky,timings,DensityCheck
      Logical DoLocK,Deco
      Logical Estimate,Update
      Integer ALGO,Nscreen
      Real*8  dmpk,ChFracMem
      Logical, External :: Is_First_Iter

      Common /CHLCAS / DoCholesky,ALGO
      COMMON /CHODENSITY/ DensityCheck
      COMMON /CHOTIME / timings
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen
      COMMON /CHOSCREEN/ Estimate,Update
      COMMON /CHOPAR/ ChFracMem
*----------------------------------------------------------------------*
      ProgName=Get_SuperName()
*----------------------------------------------------------------------*
* How was the program called?
*PAM 2009 Someone has put a number of possibilities here. Let it stand for now.
      IfVB=0
      If (ProgName(1:6).eq.'rasscf'.or.ProgName(1:6).eq.'mcpdft') Then
         ICIRST=0
*        For geometry optimizations use the old CI coefficients.
         If (.Not.Is_First_Iter()) ICIRST=1
      ELse If (ProgName(1:5).eq.'casvb') Then
         IfVB=2
         ICIRST=0
      ELse If (ProgName(1:6).eq.'loprop') Then
C        ICIRST=1 ! to be activated!
         ICIRST=0
      ELse If (ProgName(1:11).eq.'last_energy') Then
         ICIRST=1
      ELse If (ProgName(1:18).eq.'numerical_gradient') Then
         ICIRST=1
      Else
         ICIRST=0
      End If

* Initialize print levels: See output_ras.fh
* Global logical unit numbers for standard input and standard output
      IO=5
      LF=6
* Externally set default print level control. Should the program be silent?
      IPRGLB_IN=iPrintLevel(-1)
      DO I=1,7
       IPRLOC_IN(I)=IPRGLB_IN
      END DO
* Set print levels, and adjust them if needed:
      call setprlev_m(LF,IPRGLB_IN,IPRLOC_IN)
*
* SET UP SYMMETRY MULTIPLICATION TABLE:
      MUL(1,1)=1
      M=1
      DO  N=1,3
        DO  I=1,M
          DO  J=1,M
            MUL(I+M,J)=M+MUL(I,J)
            MUL(I,J+M)=MUL(I+M,J)
            MUL(I+M,J+M)=MUL(I,J)
          END DO
         END DO
        M=2*M
      END DO

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
#if defined (_MOLCAS_MPP_)
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
* assign ipCleanMask to dummy pointer
      ipCleanMask=ip_Dummy
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
      THRSX=1.D-04
* Default damping in the SXCI orbital optimization
      SXDAMP=0.0002D0
* Default thresholds used to determine convergence in CI
      THREN=1.0D-04
      THFACT=1.0D-03
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
*TRS
*      iRLXROOT=0
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
      KAVER=0
* not zero if density matrices are to be averaged.
*     KAVER=1 symmetries KSYM(1) and KSYM(2) averaged
*     KAVER=2 also symmetries KSYM(3) and KSYM(4) averaged
*
      ISUPSM=0
* make no use of supersymmetry
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
* Initialize KSDF coefficients (S Dong, 2018)
      CoefR = 1.0D0
      CoefX = 1.0D0
!      Write(6,*) ' Correlation energy scaling factor (init) is ',CoefR
!      Write(6,*) ' Exchange energy scaling factor (init) is ',CoefX
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
      LSYM=1
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
      ICLEAN=0
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
      Do i = 1,nSpeed
         if (i .le. 2) Then
            iSpeed(i) = 1
         else
C The rest is at the present time just to allow testing
            iSpeed(i) = 0
         end if
      End Do
*
      Ebel_3     = 0.0d0
      Eterna_3   = 0.0d0
      Rado_3     = 0.0d0
      Rolex_3    = 0.0d0
      Omega_3    = 0.0d0
      Tissot_3   = 0.0d0
      Piaget_3   = 0.0d0
      Candino_3  = 0.0d0
      Fortis_3   = 0.0d0
      Zenith_3   = 0.0d0
      Gucci_3    = 0.0d0
      Alfex_3    = 0.0d0
      WTC_3      = 0.0d0
      Longines_3 = 0.0d0
      Oris_2     = 0.0d0
      Movado_2   = 0.0d0
*
CSVC: lucia timers
      tsigma = 0.0d0
      tdensi = 0.0d0
*
      RETURN
      END
