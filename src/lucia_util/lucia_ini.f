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
      SUBROUTINE Lucia_Ini
*
      implicit real*8 (a-h,o-z)
C Input from RASSCF
#include "cstate.fh"
#include "lucia_ini.fh"
*

#include "mxpdim.fh"
#include "cecore.fh"
#include "cgas.fh"
#include "cprnt.fh"
#include "crun.fh"
#include "irat.fh"
#include "lucinp.fh"
#include "oper.fh"
#include "orbinp.fh"
#include "gasstr.fh"
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"
#include "memman.fh"
C Other definitions
      PARAMETER(MXPKW = 125)
      dimension isetkw(MXPKW)
*
C ==================================
C  Initialize marks in Lucia Memory
C ==================================
      Num_Marks = 0
C ===============================================
C  Some initialization to avoid compiler warnings
C ===============================================
c      MSCOMB_CC = 0
c      I_USE_NEWCCP = 0
C =======================
C  Some initial settings
C =======================
      INI_H0 = 1
*. Flag for compatibility with normal MOLCAS input format
*
      IEXPERT = 0
      NERROR = 0
      NWARN = 0
      EXTSPC = 0
      IECHO = 0
* No cc as default
*. Start out with normal integrals
c      I_USE_SIMTRH = 0
* Initialize flag to be used by densi_master to check whether the
* sigma-vector is on disk.
      iSigma_on_disk = 0
*
C ===================
C  Initialize isetkw
C ===================
      Do i = 1, mxpkw
         isetkw(i) = 0
      ENDDO
*
C =========================
C  Point group of orbitals
C =========================
      pntgrp = 1
*
C ==============================
C  Number of irreps of orbitals
C ==============================
         nirrep = nsym_molcas
         nsmcmp = nsym_molcas
         nsmob = nsym_molcas
         maxml = -1
         maxl = -1
         invcnt = -1
*
C ============================
C  Number of active electrons
C ============================
      nactel = nactel_Molcas
*
C =================
C  Inactive shells
C =================
      DO irrep = 1, nirrep
         ninash(irrep) = nish_molcas(irrep)
      ENDDO
*
C ==================
C  Secondary shells
C ==================

C ===========================
C  Two times spin projection
C ===========================
      ms2 = ms2_molcas

C ===================
C  Spin multiplicity
C ===================
      mults = iSpin_molcas
*
C ====================
C  Reference symmetry
C ====================
      irefsm = lsym_molcas
*
C =======
C  Roots
C =======
      nroot = nroots_molcas
      DO i = 1, nroot
         iroot(i) = iroot_molcas(i)
      ENDDO
*
C ============
C  Enviroment
C ============
      intimp = 6
      enviro(1:6) = 'RASSCF'
*
* ====================
* 27: Ms combinations
* ====================
      IF (MULTS .EQ. 1 .AND. ISPEED(1) .EQ. 1) THEN
         ISETKW(27) = 1
         PSSIGN = 1.0D0
      ELSE
         ISETKW(27) = 0
      ENDIF
*
*===========================================================
* 33 : Largest allowed number of Vectors in diagonalization
*===========================================================
      MXCIV = 2*NROOT
*
C ==========================
C  Storage mode for vectors
C ==========================
      icistr = 2
C =================================
C  Tell Lucia to use CSF-expansion
C =================================
c      noscf = 0
*
* ================================
* 39 : Define dimension of resolution matrices
* ================================
      MXINKA = 150
*
C ==========================================================
C  Generalized active space concept invoked, orbital spaces
C ==========================================================
      ngas=ngas_molcas
      ngssh=0
      do irrep = 1,nirrep
        do igas=1,ngas
          ngssh(irrep,igas)=ngssh_molcas(igas,irrep)
        end do
      end do
*
C ==================================================
C  Generalized active space occupation restrictions
C ==================================================
      ncispc = 1
      do i=1,ngas
        do j=1,2
          igsoccx(i,j,1)=igsoccx_molcas(i,j)
        end do
      end do
*
C ==========================
C  Energy convergence of CI
C ==========================
      thres_e = thre_molcas
*
C ==========================================
C  Sequency of calculations to wavefunction
C ==========================================
      ncmbspc = ncispc
      nseqci(1) = 1
      cseqci(1,1) = 'CI'
      iseqci(1,1) = itmax_molcas
*
C ======================================================
C  No calculation, save CI vector info, expand CI vector
C ======================================================
      INOCALC   = INOCALC_MOLCAS
      ISAVE_EXP = ISAVE_EXP_MOLCAS
      IEXPAND   = IEXPAND_MOLCAS
*
C ================================================
C  Length of smallest block og C an Sigma vectors
C ================================================
      lcsblk = 3000000
*
* =============================================================
* 66 : No MO-AO file
* =============================================================
*
*. No MO-AO file
      NOMOFL = 1
*
C =================================
C  Calculation of density matrices
C =================================
      idensi = 2
*
C ==================================================================
C  Isimsym: Treat all symmetryblocks with given type simultaneously
C ==================================================================

C ==============================================
C  Particle hole simplifications, default is no
C ==============================================
      IUSE_PH = 0
*
* ==============================================
* 87 : Allow the sigma routine to take advice
* ==============================================
*
      IADVICE = 1
*
C ============================================================
C  Transform CI vectors to alternative orbital representation
C ============================================================
      itraci     = 1
      itraci_cr  = 'REST'
      itraci_cn  = ' '
      If (ipt2_molcas .eq. 1) Then
         itraci_cn = 'CANO'
      Else
         itraci_cn = 'NATU'
      Endif
*
C ===========
C HEXS - Highly excited states
C DEXS - Doubly excited states
C ===========
      I_ELIMINATE_GAS = I_ELIMINATE_GAS_MOLCAS
      N_ELIMINATED_GAS = N_ELIMINATED_GAS_MOLCAS
*     WRITE(6,*) I_ELIMINATE_GAS_MOLCAS,N_ELIMINATED_GAS_MOLCAS
      DO IGAS = 1, N_ELIMINATED_GAS
        IELIMINATED_IN_GAS(IGAS) = IELIMINATED_IN_GAS_MOLCAS(IGAS)
      END DO
*
C =============
C  Printlevels
C =============
      iprcix = iprci_molcas-2
      iprorb = 0
      iprdia = 0
      iprxt = 0
      iprocc = 0
      iprden = (iprci_molcas-2)/2
      iprrsp = 0
      iprpro = 0
      iprcc = 0
*
C ==============
C  Set defaults
C ==============
*
**********************************************************************
*                                                                    *
* Part 2 : Insert defaults for missing optional keywords             *
*          and print error messages for missing mandatory keywords   *
*                                                                    *
**********************************************************************
*
      NMISS = 0
*
*
*
* 8 : Core orbitals, only of interest if EXTSPC .ne. 0
*
         CALL ISETVC(NRS0SH,0,NIRREP)
         MNHR0 = 0
*
* 9 : RAS 1 orbitals
*
         CALL ISETVC(NRSSH(1,1),0,NIRREP)
*
* 10 : RAS 2 orbitals
*
         CALL ISETVC(NRSSH(1,2),0,NIRREP)
*
* 11 : RAS 3 orbitals
*
      IMLCR3 = 0
         CALL ISETVC(NRSSH(1,3),0,NIRREP)
         IF(INTSPC.EQ.2) THEN
*. Use information from one-electron integral file to obtain
* default
            IMLCR3 = 1
         END IF
*
* 13 : Secondary space
*
         MXR4TP = 1
         DO 210 IR4TP = 1, MXR4TP
            CALL ISETVC(NRS4SH(1,IR4TP),0,NIRREP)
 210     CONTINUE
         MXER4 = 0
*
* 14 : occupation restrictions for Reference space
*
         MNRS1R = MNRS10
         MXRS3R = MXRS30
*
* 15 : Selection of active configurations
*
*. Standard is no selection
         INTSEL = 0
*
* 20 : Diagonalization routine
*
*. Standard is currently MICDV*
         IDIAG = 1
*
* 21 : Explicit Hamiltonian
*
*. Default is no explicit Hamiltonian
         MXP1 = 0
         MXP2 = 0
         MXQ  = 0
*
* 22 : Largest allowed number of CI iterations per root
*
*. Default is 20 ( not active I expect )
         MAXIT = 20
*
* 23 : Restart option
*
*. Default is no explicit Hamiltonian
         IRESTR = 0
*
* 25 : INCORE option for integrals
*
         IF(EXTSPC.EQ.0 ) THEN
            INCORE = 1
         ELSE
            INCORE = 0
         END IF
*
* 26 : DELETEd shells
*
*. If CAS + Active have been set or RAS + Ras3 have been set,
*. obtain for MOLCAS Interface from number of basis functions
         IF(INTSPC.EQ.1) THEN
            IMLCR3 = 2
         ELSE
            CALL ISETVC(NDELSH,0,NIRREP)
         END IF
*
* 27 : Ms combinations
*
      IF(ISETKW(27).EQ.0) THEN
         PSSIGN = 0.0D0
         ISETKW(27) = 2
      ELSE IF(MS2.NE.0) THEN
         WRITE(6,*) ' Spin combinations only allowed with MS2 = 0'
        WRITE(6,*) ' Your value of MS2 = ',MS2, ' differs from zero'
         WRITE(6,*) ' LUCIA will neglect your nice suggestion '
         WRITE(6,*)  ' to use spin combinations '
         PSSIGN = 0.0D0
         ISETKW(27) = 2
      END IF
*
* 28 : Ml combinations
*
         PLSIGN = 0.0D0

      IF(PSSIGN.EQ.0.0D0) THEN
         IDC = 1
      ELSE IF(PSSIGN.NE.0.0D0) THEN
         IDC = 2
      endif

C?    WRITE(6,* ) ' TEST readin IDC = ', IDC
*
* ====================================================================
* 44 : Use Minimal operatioon count method for alpha-alpha and beta-beta
* ====================================================================
* TEST OF PERFORMANCE
      MOCAA      = ISPEED(3)
      MOCAB      = ISPEED(4)


       IUSE_PA=0

*
* 33 : Number of Ci vectors in subspace
*
*
*
      IF(IDIAG.EQ.2 .AND. MXCIV.GT.2 ) THEN
        MXCIV = 2
        NWARN = NWARN + 1
        WRITE(6,*) ' Warning : You have specified TERACI '
        WRITE(6,*) '           I allow myself to set MXCIV = 2 '
        WRITE(6,*)
        WRITE(6,*) '                   Best Wishes    '
        WRITE(6,*) '                      Lucia       '
      END IF

*
* 34 : CI storage mode
*
*
* 35 : Employ CSF expansion ?
*
*. Default is no ( only possibility at the moment )
* CSF expansion must only be used when two vectors are stored in CORE
*
* 37 : Avoid any readin of integrals ( useful for obtaining
*      size of CI expansion etc.
*
        NOINT = 0
*
* 38 : Dump integrals in formatted form : Default is No
*
c        IDMPIN = 0

*
* 40 : Use CJKAKB intermediate matrices in alpha-beta loop,
*      Default is  YES !!!!!
*
        ICJKAIB = 1
*
*  41 : Initial CI in reference space, default is : No
*
         INIREF = 0
*
*  42 : Restart with CI in reference space
*
         IRESTRF = 0
*
*
* Core energy : Default is 0 / MOLCAS : Value read in !
*
         ECORE = 0.0D0
*
*. Use perturbation theory for zero order space . def is no !
*
      IF(ISETKW(47).EQ.0) THEN
        IPERT = 0
        NPERT = 0
        ISETKW(47) = 2
*. Else ensure that a CI in reference space is performed
      ELSE
        INIREF = 1
      END IF
*
*
*. 48 : Approximate Hamiltonian in reference space : NO !!
*
      IF(ISETKW(48).EQ.0) THEN
        IAPRREF = 0
        MNRS1RE = MNRS1R
        MXRS3RE = MXRS3R
        ISETKW(48) = 2
      END IF
*
*. 49 : Approximate Hamiltonian in zero order space : NO !!
*
      IF(ISETKW(49).EQ.0) THEN
        IAPRZER = 0
        MNRS1ZE = MNRS10
        MXRS3ZE = MXRS30
        ISETKW(49) = 2
      END IF
*
* 50 : GAS shells must be defined
*
*
* 52 : Combination of gasspaces : Default is just to take each  space
*      By itself
*
      IF(ISETKW(52).EQ.0) THEN
        NCMBSPC = NCISPC
        DO ICISPC = 1, NCISPC
          LCMBSPC(ICISPC) = 1
          ICMBSPC(1,ICISPC) = ICISPC
        END DO
        ISETKW(52) = 2
      END IF
*
* 53 : Convergence threshold for CI
*
*
* 54 : General sequencer : default is just straight sequence
*      of CI with default number of iterations
*
* 55 : EKT calculation : Default is no
*
      IF(ISETKW(55).EQ.0) THEN
        IEXTKOP = 0
        ISETKW(55) = 2
      END IF
*
*. 56 : Default Machine : Good old BIM machine
*
* 57 : Allow first order correction to be saved on DISC
*     (For vector free calculations )
*     Default is : NO !!
      IF(ISETKW(57).EQ.0) THEN
        IC1DSC = 0
        ISETKW(57) = 2
      END IF
*
* 58 : Restrictions on interactions of perturbation
*
*. Default is : no
      IF(ISETKW(58).EQ.0) THEN
        IH0SPC = 0
        ISETKW(58) = 2
      END IF
*
* 59 : Type of perturbation in subspaces spaces
*
* Default is specified by IPART from keyword PERTU
      IF(ISETKW(59).EQ.0) THEN
       ISETKW(59) = 2
       IF(IH0SPC.NE.0) THEN
         DO JPTSPC = 1, NPTSPC
           IH0INSPC(JPTSPC) = IPART
         END DO
       END IF
      END IF
*
* 60 : Reference Root, default is NROOT
*
*. Should be less or equal to NROOT
      IF(ISETKW(60).EQ.1) THEN
        IF(IRFROOT.GT.NROOT) THEN
          WRITE(6,*) ' Reference root (RFROOT) larger '
          WRITE(6,*) ' than total number of roots (NROOT) '
          WRITE(6,*) ' CHANGE NROOT or RFROOT '
          NMISS = NMISS + 1
        END IF
      END IF

      IF(ISETKW(60).EQ.0) THEN
       ISETKW(60) = 2
       IRFROOT = NROOT
      END IF
*
* 61 : H0EX : Orbital spaces in which exaxt Hamiltonian is used
*      No default
*.
*. Is H0EX required ( Has H0FRM = 5 been used )
      IUSED = 0
      IF(ISETKW(59).EQ.1) THEN
         IUSED = 0
         DO JPTSPC = 1, NPTSPC
           IF( IH0INSPC(JPTSPC) .EQ. 5 ) IUSED = 1
         END DO
       END IF
       IF(IUSED.EQ.0.AND.ISETKW(61).EQ.0) THEN
*. No exact spaces included and none have been defined !
         NH0EXSPC = 0
         IH0EXSPC(1) = -1
       END IF
       IF(IUSED.EQ.1.AND.ISETKW(61).EQ.0) THEN
*. Needed, but not supplied
          WRITE(6,*) ' You have specified that zero order operator'
          WRITE(6,*) ' Include exact Hamilton operator in subspace'
          WRITE(6,*) ' You should then also supply Keyword H0EX '
          NMISS = NMISS + 1
       END IF
*
*. If perturbation theory will be invoked be sure that the
*. form of perturbation theory has been specified through
* KEYWORD PERTU ( number 47 as you maybe know )
      IDOPERT = 0
      DO JCMBSPC = 1, NCMBSPC
        DO JSEQCI = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQCI,JCMBSPC).EQ.-5 ) IDOPERT = 1
        END DO
      END DO
*
      IF(IDOPERT.EQ.1 .AND. IPERT.EQ.0) THEN
        WRITE(6,*) ' Perturbation theory will be used '
        WRITE(6,*) ' Please specify form through PERTU keyword '
        NMISS = NMISS + 1
      END IF
*
*. 62 : Default Handling of degenrences of initial CI vectors
*.      Default is : No action
*
      IF(ISETKW(62).EQ.0) THEN
        INIDEG = 0
        ISETKW(62) = 2
      END IF
*
*. 63 : Use F + Lambda(H-F) as operator instead of H
*.      Default is : No i.e Lambda = 1
*
      IF(ISETKW(63).EQ.0) THEN
        XLAMBDA = 1.0D0
        ISETKW(63) = 2
      END IF
*
*. 64 : Smallest block in batch of C and sigma
*.      Default is zero
*
*
*. 66 : NO MO file : Default is no access to MO-AO file
*
*. 68 : Type of natural orbitals, default is natural orbitals
*
      IF(ISETKW(68).EQ.0) THEN
        ISETKW(68) = 2
        IFINMO = 1
      END IF
*
*. 69 : Default Threshold for individual energy correction = 0.0
*
      IF(ISETKW(69).EQ.0) THEN
        E_THRE = 0.0D0
        ISETKW(69) = 2
      END IF
*
*. 70 : Default Threshold for wave individual function corrections = 0.0
*
      IF(ISETKW(70).EQ.0) THEN
        C_THRE = 0.0D0
        ISETKW(70) = 2
      END IF
*
*. 71 : Default Threshold for total energy corrections = 0.0
*
      IF(ISETKW(71).EQ.0) THEN
        E_CONV = 0.0D0
        ISETKW(71) = 2
      END IF
*
*. 72 : Default Threshold for total wave function correction = 0.0
*
      IF(ISETKW(72).EQ.0) THEN
        C_CONV = 0.0D0
        ISETKW(72) = 2
      END IF
*
*. 73 : Perform Class selection : Default if Yes if TERACI is used
*
      IF(ISETKW(73).EQ.0) THEN
        IF(IDIAG.EQ.1) THEN
          ICLSSEL = 0
        ELSE IF (IDIAG.EQ.2) THEN
          ICLSSEL = 1
        END IF
        ISETKW(73) = 2
      END IF
*
* 74 : Calculation of density matrices : Default is
*       calculation of one-body density
*
*. If IDENSI was set to zero and properties were requested
*  overwrite input to obtain 1-el matrix
      IF(IDENSI.EQ.0.AND.ISETKW(80).EQ.1) THEN
        WRITE(6,*) ' You have specified calculation of'
        WRITE(6,*) ' one-electron properties, and this'
        WRITE(6,*) ' requires the calculation of the '
        WRITE(6,*) ' one-electron density. '
        WRITE(6,*)
        WRITE(6,*) ' You have, however, specified IDENSI=0'
        WRITE(6,*) ' corresponding  to no densities'
        WRITE(6,*)
        WRITE(6,*) ' I will allow myself to modify your'
        WRITE(6,*) ' input to allow calculation of the '
        WRITE(6,*) ' one-electron densities, so property'
        WRITE(6,*) ' calculation can proceed as planned '
        WRITE(6,*)
        WRITE(6,*)                        ' Lucia '
*. and do it
        IDENSI = 1
      END IF
*. If CC is performed, one- and two- particle densities are
*  used in current simple-minded implementation.
*. Two-electron density also needed for MCSCF
*
* 75 : Perturbation expansion of EKT, default is no
*
      IF(ISETKW(75).EQ.0) THEN
        IPTEKT = 0
        ISETKW(75) = 2
      END IF
*
* 76 : Root for zero order operator , default is NROOT
*
*. Should be less or equal to NROOT
      IF(ISETKW(76).EQ.1) THEN
        IF(IH0ROOT.GT.NROOT) THEN
          WRITE(6,*) ' Zero order operator root (H0ROOT) larger '
          WRITE(6,*) ' than total number of roots (NROOT) '
          WRITE(6,*) ' CHANGE NROOT or H0ROOT '
          NMISS = NMISS + 1
        END IF
      END IF
      IF(ISETKW(76).EQ.0) THEN
       ISETKW(76) = 2
       IH0ROOT = NROOT
      END IF
*
* 77 : NO restart from previous vectors in calculation 2
*      Deafault is NO NO, ie. restart in calc 2
*
      IF(ISETKW(77).EQ.0) THEN
        IRST2 = 1
        ISETKW(77) = 2
      END IF
*
* 78 : skip initial energy evaluations - if possible
*
      IF(ISETKW(78).EQ.0) THEN
        ISKIPEI = 1
        ISETKW(78) = 2
      END IF
*
* 79 : Symmetry of x,y,z - needed for property calculations
*
      IF(ISETKW(79).EQ.0) THEN
*. Problematic if Properties should be calculated
       IF(ISETKW(80).EQ.1.OR.ISETKW(81).EQ.1.OR.ISETKW(82).EQ.1)
     & THEN
         WRITE(6,*) ' Symmetry of X,Y,Z has not been given'
         WRITE(6,*) ' You have to specify this for property calc'
         WRITE(6,*) ' Please add KEYWORD XYZSYM '
         NMISS = NMISS + 1
         ISETKW(79) = -1
       ELSE
*. Is not needed, just supply zeroes
         DO ICOMP = 1, 3
           IXYZSYM(ICOMP) = 0
         END DO
         ISETKW(79) = 2
       END IF
      END IF
*
* 80 : Property calculation, default is no
*
      IF(ISETKW(80).EQ.0) THEN
        NPROP = 0
        ISETKW(80) = 2
      END IF
*
* 81 : Transition properties , default is no
*
      IF(ISETKW(81).EQ.0) THEN
        ITRAPRP = 0
        ISETKW(81) = 2
      END IF
*
* 82 : Response properties , default is no
*
      IF(ISETKW(82).EQ.0) THEN
        IRESPONS = 0
        ISETKW(82) = 2
        NRESP = 0
        N_AVE_OP = 0
      END IF
*. Properties should be defined if transition properties are
*. invoked
      IF(ITRAPRP.NE.0.AND.NPROP.EQ.0) THEN
        WRITE(6,*)
     &  ' You have specified transition property calculation'
        WRITE(6,*)
     &  ' (keyword TRAPRP) but no property labels have been supplied'
        WRITE(6,*)
     &  '(Keyword PROPER). Transition densities will be obtained '
      END IF
*
* 83 : Max number of iterations in linear equations
*
      IF(ISETKW(83).EQ.0) THEN
        MXITLE = 20
        ISETKW(83) = 2
      END IF
*
* 85 : Root homing, default is no
*
      IF(ISETKW(85).EQ.0) THEN
        IROOTHOMING = 0
        ISETKW(85) = 2
      END IF
*
* 90 : Perturbation expansion of Fock matrix : default is no
*
      IF(ISETKW(90).EQ.0) THEN
       IPTFOCK = 0
       ISETKW(90) = 2
      END IF
*
* 91 : Print final CI vectors : default is no
*
      IF(ISETKW(91).EQ.0) THEN
       IPRNCIV = 0
       ISETKW(91) = 2
      END IF
*
* 93 : End Calculation with integral transformation
*
      IF(ISETKW(93).EQ.0) THEN
       ITRA_FI = 0
       ISETKW(93) = 2
      END IF
*
* 94 : Initialize Calculation with integral transformation
*
      IF(ISETKW(94).EQ.0) THEN
       ITRA_IN = 0
       ISETKW(94) = 2
      END IF
*. Requires access to MO-AO file
       IF(ITRA_IN.EQ.1.AND.NOMOFL.EQ.1) THEN
         WRITE(6,*) ' Integral transformation required, '
         WRITE(6,*) ' but no mo-ao file accessible      '
         WRITE(6,*) ' REMOVE KEWORD NOMOFL '
         ISETKW(94) = -1
         NERROR = NERROR + 1
       END IF
*
* 95 : Multispace optimization in each run, default is no
*
      IF(ISETKW(95).EQ.0) THEN
        MULSPC = 0
        IFMULSPC = 0
        LPAT = 0
        ISETKW(95) = 2
      END IF
*
*
* Expert mode ( neglect mistyped keywords ) : default is no expert
*
      IF(ISETKW(97).EQ.0) THEN
        IEXPERT = 0
        ISETKW(97) = 2
      END IF
*
* Number of roots to be converged : default is total number of roots
*
      IF(ISETKW(98).EQ.0) THEN
        NCNV_RT = NROOT
        ISETKW(98) = 2
      END IF
*
* 100 : Do quantum dot calculation, default is no
*
c      IF(ISETKW(100).EQ.0) THEN
c        IDODQ = 0
c        ISETKW(100) = 2
c      END IF
*
* 101: Restrict MS2 at some intermediate level : default is no way
*
      IF(ISETKW(101).EQ.0) THEN
        I_RE_MS2_SPACE = 0
        I_RE_MS2_VALUE = 0
        ISETKW(101) = 2
      END IF
*
*
* 103 : Treat all TT blocks with given types simultaneously : Default is no
*

        ISIMSYM = 0

*
*. Thresholds only active in connection with IDIAG = 2,
*. Check and maybe issue a warning
      IF(IDIAG.EQ.2) THEN
*. Check to ensure that zero or two thresholds were  set,
        IF(ISETKW(69).NE.ISETKW(70)) THEN
          WRITE(6,*)
     &    ' Only a single threshold (E_THRE or C_THRE) '
          WRITE(6,*)
     &    ' on individual determinants given. '
          WRITE(6,*)
     &    ' One of the thresholds vanishes therefore and '
          WRITE(6,*)
     &    ' all determinants will therefore be included  '
          WRITE(6,*)
          WRITE(6,*) '                   Warns '
          WRITE(6,*)
          WRITE(6,*) '                   LUCIA  '
        END IF
      ELSE
*. Good old diagonalization, thrsholds not active
        IF(ISETKW(69).EQ.1.OR.ISETKW(70).EQ.1) THEN
          WRITE(6,*)
     &    ' Thresholds on selection of individual coefficients '
          WRITE(6,*)
     &    ' are only active in connection with keyword TERACI  '
          WRITE(6,*)
          WRITE(6,*) '                   Warns '
          WRITE(6,*)
          WRITE(6,*) '                   LUCIA  '
        END IF
      END IF
*
*
*
*
      IF(NMISS.NE.0) THEN
        WRITE(6,'(1H ,A,I9)')
     &  ' Number of missing required keyword ', NMISS
        WRITE(6,*)
     &  ' You have wounded me I give up '
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
     & '     An expert is a man who has made all the mistakes,'
        WRITE(6,*)
     &  '     which can be made, in a very narrow field        '
        WRITE(6,*)
     &  '                                                      '
        WRITE(6,*)
     &  '                                      Niels Bohr      '
        IF(IEXPERT.EQ.0) THEN
*          STOP
          CALL SYSABENDMSG('lucia_util/lucia_ini','Input error',' ')
        ELSE
          WRITE(6,*) ' Processing continues (EXPERT mode )'
        END IF
      END IF
*. Open one-electron file to obtain core energy and
*. Number of MO's and AO's
      IF(NOINT.EQ.0.AND.ENVIRO(1:4).NE.'NONE') THEN
         ecore_env = potnuc_molcas
* Manual setting of naos_env and nmos_env
         nmos_env(1) = nbas_molcas(1)
         nmos_env(2) = nbas_molcas(2)
         nmos_env(3) = nbas_molcas(4)
         nmos_env(4) = nbas_molcas(3)
         nmos_env(5) = nbas_molcas(7)
         nmos_env(6) = nbas_molcas(6)
         nmos_env(7) = nbas_molcas(8)
         nmos_env(8) = nbas_molcas(5)
         naos_env(1) = norb_molcas(1)
         naos_env(2) = norb_molcas(2)
         naos_env(3) = norb_molcas(4)
         naos_env(4) = norb_molcas(3)
         naos_env(5) = norb_molcas(7)
         naos_env(6) = norb_molcas(6)
         naos_env(7) = norb_molcas(8)
         naos_env(8) = norb_molcas(5)
         ECORE = ECORE_ENV
      ELSE
        WRITE(6,*) ' GETOBS and CHK_ORBDIM not called '
        ECORE = 0.0D0
      END IF
*. Check number of orbitals and insert occupations for ALL/REST

************************************************************
*                                                          *
* Part 3 : Print input                                     *
*                                                          *
************************************************************
*
*. Machine in use
*
*. Type of reference state
      IF (IPRCIX .GE. 100) THEN
         WRITE(6,*)
         WRITE(6,*) '*************************************'
         WRITE(6,*) '*  Symmetry and spin of CI vectors  *'
         WRITE(6,*) '*************************************'
         WRITE(6,*)
*. Point group
           WRITE(6,'(1H ,A)')
     &       '     Point group ............ D2H'
*.Spatial symmetry
           WRITE(6,'(1H ,A,I1)')
     &       '     Spatial symmetry ....... ', IREFSM
*.Spin
         WRITE(6,'(1H ,A,I2)')
     &      '     2 times spinprojection  ', MS2
*.Number of active electrons
         WRITE(6,'(1H ,A,I2)')
     &     '     Active electrons .....  ', NACTEL
         WRITE(6,*)
         WRITE(6,*) '*********************************************'
         WRITE(6,*) '*  Shell spaces and occupation constraints  *'
         WRITE(6,*) '********************************************* '
         WRITE(6,*)
*
*
      IF(XLAMBDA.NE.1.0D0) THEN
        WRITE(6,*)
        WRITE(6,'(A,F13.8)')
     &  ' Modified operator H(l) = l*F + l*(H-F) used with l =',XLAMBDA
c        IF(IUSEH0P.EQ.0) THEN
         WRITE(6,'(A)')  ' Zero-order operator without projection used '
c        ELSE
c         WRITE(6,'(A)')  ' Zero-order operator with projection used '
c        END IF
        IF(IRESTR.EQ.0) THEN
        WRITE(6,*)
     &  ' Notice : This madness starts  in second calculation'
        ELSE
         WRITE(6,*) ' You have specified a calculation with modified '
         WRITE(6,*) ' Hamiltonian (the LAMBDA option) and RESTART '
         WRITE(6,*) ' so this is what I will do '
         WRITE(6,*)
         WRITE(6,*) '   1:) Perform CI in space 1 to obtain Hamiltonian'
         WRITE(6,*) '       (no RESTART in this space )'
         WRITE(6,*) '   2:) CI calculation in space 2  with '
         WRITE(6,*) '       modified Hamiltonian and RESTART from LU21'
         WRITE(6,*) ' Space 2 should therefore correspond to the'
         WRITE(6,*) ' restarted calculation '
       END IF
      END IF
*
      WRITE(6,*)
      WRITE(6,*) '***********'
      WRITE(6,*) '*  Roots  *'
      WRITE(6,*) '*********** '
      WRITE(6,*)
      WRITE(6,'(1H ,A,I3)')
     &  '     Number of roots to be included  ', NROOT
      WRITE(6,'(1H ,A,(20I3))')
     &  '     Roots to be obtained ', (IROOT(I),I=1, NROOT )
      WRITE(6,'(1H ,A,I3)')
     &  '     Number of roots to be converged ', NCNV_RT

*
      WRITE(6,*)
      WRITE(6,*) '**************************'
      WRITE(6,*) '*  Run time definitions  *'
      WRITE(6,*) '************************** '
      WRITE(6,*)
*. Program environment
      WRITE(6,'(A,A6)')  '      Program environment... ', ENVIRO
*
*. Integral import
      IF(NOINT.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  '     No integrals will be read in       '
      ELSE IF(NOINT.EQ.0) THEN
*. Integral storage
      IF(INCORE.EQ.1) WRITE(6,'(1H ,A)')
     &  '     All integrals stored in core'
      END IF
      WRITE(6,*)
* ( END IF for NOINT
*. CSF or SD expansion
      IF(NOCSF.EQ.0) THEN
        WRITE(6,'(1H ,A)')
     &  '     CI optimization performed with CSF''s '
      ELSE
        WRITE(6,'(1H ,A)')
     &  '     CI optimization performed with SD''s '
      END IF
*. Ms,Ml combinations
      IF(ISETKW(27).EQ.1) THEN
        WRITE(6,'(1H ,A,F8.3)')
     &  '     Spin combinations used with sign ',PSSIGN
      END IF
      IF(ISETKW(28).EQ.1) THEN
        WRITE(6,'(1H ,A,F8.3)')
     &  '     ML   combinations used with sign ',PLSIGN
      END IF
*. Initial approximation to vectors
      WRITE(6,*)
      IF(IRESTR.EQ.1.AND.IRESTRF.EQ.0) THEN
         WRITE(6,'(1H ,A)')
     &  '     Restarted calculation '
      ELSE IF(IRESTRF.EQ.1) THEN
         WRITE(6,'(1H ,A)')
     &  '     Restarted calculation from REFERENCE space expansion'
      ELSE
         IF(MXP1.NE.0) THEN
           WRITE(6,'(1H ,A)')
     &  '     Initial vectors obtained from explicit Hamiltonian'
         ELSE IF(MXP1.EQ.0) THEN
           WRITE(6,'(1H ,A)')
     &  '     Initial vectors obtained from diagonal'
         END IF
      END IF
*. Handling of degenerencies of initial vectors
      IF(INIDEG.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  '     Symmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.-1) THEN
        WRITE(6,'(1H ,A)')
     &  '     Antiymmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.0) THEN
        WRITE(6,'(1H ,A)')
     &  '     No combination of degenerate initial vectors'
      END IF
*. No calculation, save CI vector information and expand CI vector
      IF(INOCALC.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  ' No calculation will be performed '
      END IF
      IF(ISAVE_EXP.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  ' Save CI vector information '
      END IF
      IF(IEXPAND.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  ' Expand shorter CI vector in longer one '
      END IF
*. Ms,Ml combinations
C     IF(ISETKW(27).EQ.1) THEN
C       WRITE(6,'(1H ,A,F8.3)')
C    &  '     Spin combinations used with sign ',PSSIGN
C     END IF
C     IF(ISETKW(28).EQ.1) THEN
C       WRITE(6,'(1H ,A,F8.3)')
C    &  '     ML   combinations used with sign ',PLSIGN
C     END IF
*. CI storage mode
      WRITE(6,*)

        WRITE(6,*)
     &  '     3 symmetry blocks will be held in core '

      IF(LCSBLK.NE.0) WRITE(6,'(A,I10)')
     &  '      Smallest allowed size of sigma- and C-batch ',LCSBLK
      WRITE(6,'(1H ,A,I4)')
     &  '     Dimension of block of resolution strings ', MXINKA
c      IF(IUSE_PH.EQ.1) THEN
c        WRITE(6,'(1H ,A)')
c     &  '     Particle-hole separation used '
c      ELSE
        WRITE(6,'(1H ,A)')
     & '      Particle-hole separation not used '
c      END IF
*
      IF(IADVICE.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  '     Advice routine call to optimize sigma generation'
      END IF
*
c      IF(IUSE_PA.EQ.1) THEN
c        WRITE(6,'(1H ,A)')
c     &  '     Strings divided into active and passive parts'
c      ELSE
        WRITE(6,'(1H ,A)')
     &  '     Strings not divided into active and passive parts'
c      END IF
c      IF(ISIMSYM.EQ.1) THEN
c        WRITE(6,'(1H ,A)')
c     &  '     ALl TTS blocks with given types treated in sigma'
c      END IF
c      IF(IUSE_HW .EQ. 1) THEN
c        WRITE(6,*) ' Hardwired routines in use '
c      END IF
*
      WRITE(6,*)
      IF(IDENSI.EQ.0) THEN
        WRITE(6,'(1H ,A)')
     &  '     No calculation of density matrices  '
      ELSE IF(IDENSI.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &  '     One-body density matrix calculated'
      ELSE IF(IDENSI.EQ.2) THEN
        WRITE(6,'(1H ,A)')
     &  '     One- and two-body density matrices  calculated'
      END IF
      WRITE(6,*)
C?    IF(MOCAA.NE.0) WRITE(6,'(1H ,A,I4)')
C?   &  '     MOC method used for alpha-alpha+beta-beta loop '
C?    IF(MOCAB.NE.0) WRITE(6,'(1H ,A,I4)')
C?   &  '     MOC method used for alpha-beta loop            '
*
*. Diagonalization information
      WRITE(6,'(1H ,A)')
     &  '     CI diagonalization : '
      WRITE(6,'(1H ,A)')
     &  '     ==================== '
*. Subspace Hamiltinian
      IF(MXP1+MXP2+MXQ .EQ.0) THEN
        WRITE(6,'(1H ,A)')
     &  '        No subspace Hamiltonian '
      ELSE
        WRITE(6,'(1H ,A,3I4)')
     &  '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
      END IF
*. Diagonalizer
      IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
        WRITE(6,'(1H ,A)')
     &    '        Diagonalizer : MINDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.2) THEN
        WRITE(6,'(1H ,A)')
     &    '        Diagonalizer : MICDV* '
      ELSE IF(IDIAG.EQ.2) THEN
      WRITE(6,'(1H ,A)')
     &    '        Diagonalizer : PICO*  '
      END IF
        WRITE(6,'(1H ,A)')
     &      '        Simple diagonal used as preconditioner  '
*. Root homing
      IF(IROOTHOMING.EQ.1) THEN
      WRITE(6,'(1H ,A)')
     &  '        Root homing will be used '
      ELSE
      WRITE(6,'(1H ,A)')
     &  '        No root homing '
      END IF
*. No restart in CI calc 2
      IF(IRST2.EQ.0) THEN
      WRITE(6,'(1H ,A)')
     &  '        No restart from previous vectors in second calc '
      END IF
      IF(ISKIPEI.EQ.1) THEN
      WRITE(6,'(1H ,A)')
     &  '        Initial energy evaluations skipped after first calc'
      WRITE(6,'(1H ,A)')
     &  '        (Only active in connection with TERACI )'
      END IF
*. Number of iterations
C     WRITE(6,'(1H ,A,I2)')
C    &  '        Allowed number of iterations    ',MAXIT
*. Number of CI vectors in subspace
      WRITE(6,'(1H ,A,I2)')
     &  '        Allowed Dimension of CI subspace',MXCIV
*
      WRITE(6,'(1H ,A,E12.5)')
     &  '        Convergence threshold for energy',THRES_E
*. Multispace (multigrid info )
c      IF(MULSPC.EQ.1) THEN
c        WRITE(6,'(1H ,A,I3)')
c     &    '        Multispace method in use from space ',
c     &             IFMULSPC
c        WRITE(6,*)
c     &    '        Pattern '
c        CALL IWRTMA(IPAT,1,LPAT,1,LPAT)
c      ELSE
        WRITE(6,'(1H ,A)')
     &    '        No multispace method in use '
c      END IF
*
      WRITE(6,*)
      IF(IDIAG.EQ.2) THEN
        WRITE(6,'(1H ,A,E12.5)')
     &   '        Individual second order energy threshold',E_THRE
        WRITE(6,'(1H ,A,E12.5)')
     &   '        Individual first order wavefunction threshold',C_THRE
        IF(ICLSSEL.EQ.1) THEN
         WRITE(6,*)
         WRITE(6,'(1H ,A)')
     &   '         Class selection will be performed : '
         WRITE(6,'(1H ,A)')
     &   '         =================================== '
         WRITE(6,'(1H ,A,E12.5)')
     &    '          Total second order energy threshold',E_CONV
         WRITE(6,'(1H ,A,E12.5)')
     &    '          Total first order wavefunction threshold',C_CONV
        ELSE
         WRITE(6,'(1H ,A)')
     &'            No class selection in iterative procedure '
        END IF
      END IF
*
      IF(IPERT.NE.0) THEN
        WRITE(6,'(1H ,A)')
     &    '     Perturbation calculation'
        WRITE(6,'(1H ,A)')
     &  '     ======================= '
        WRITE(6,*)
     &  '        Root Choosen as zero order state ', IRFROOT
        WRITE(6,*)
     &  '        Root used for zero order operator ', IH0ROOT
COLD    IF(MPORENP.EQ.1) THEN
COLD    WRITE(6,*)
COLD &  '        Moller Plesset partitioning '
COLD    ELSE IF (MPORENP.EQ.2) THEN
COLD    WRITE(6,*)
COLD &  '        Epstein-Nesbet partitioning '
COLD    ELSE IF  (MPORENP.EQ.0) THEN
COLD    WRITE(6,*)
COLD &  '        One-body Hamiltonian readin '
COLD    END IF
        IF(IE0AVEX.EQ.1) THEN
          WRITE(6,*)
     &  '        Expectation value of H0 used as zero order energy '
        ELSE IF( IE0AVEX.EQ.2) THEN
          WRITE(6,*)
     &  '        Exact energy of reference used as zero order energy'
        END IF
        WRITE(6,*)
     &  '        Correction vectors obtained through  order ', NPERT
        IF(IH0SPC.EQ.0) THEN
        WRITE(6,*)
     &  '        No restrictions on perturbation interactions '
        ELSE
        WRITE(6,*)
     &  '        Perturbation restricted to interactions in subspaces'
        END IF
*
        IF(IH0SPC.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*)
     &  '        Number of perturbation subspaces ', NPTSPC
        WRITE(6,*)
        WRITE(6,*)
     &  '        ======================== '
        WRITE(6,*)
     &  '        Perturbation subspaces : '
        WRITE(6,*)
     &  '        ======================== '
        DO JPTSPC = 1, NPTSPC
COLD      WRITE(6,'(A)')
COLD &     ' ====================================================== '
          WRITE(6,'(A)')
          WRITE(6,'(7X,A)') '         Min. occ    Max. occ '
          WRITE(6,'(7X,A)') '         ========    ======== '
          DO IGAS = 1, NGAS
            WRITE(6,'(7X,A,I2,3X,I3,9X,I3)')
     &      '   GAS',IGAS,IOCPTSPC(1,IGAS,JPTSPC)
     &                   ,IOCPTSPC(2,IGAS,JPTSPC)
          END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,'(7X,A)') ' ========================================'
        WRITE(6,'(7X,A)') ' Approximate Hamiltonian in CI subspaces '
        WRITE(6,'(7X,A)') ' ========================================'
        WRITE(6,'(7X,A)')
        WRITE(6,'(7X,A)') '    Subspace          H(apr)   '
        WRITE(6,'(7X,A)') '  ============================='
        WRITE(6,'(7X,A)')
        DO JPTSPC = 1, NPTSPC
          IF(IH0INSPC(JPTSPC).EQ.1) THEN
            WRITE(6,'(12X,I3,8X,A)')
     &      JPTSPC, ' Diagonal Fock operator '
          ELSE IF(IH0INSPC(JPTSPC).EQ.2) THEN
            WRITE(6,'(12X,I3,8X,A)')
     &      JPTSPC, ' Epstein-Nesbet operator'
          ELSE IF(IH0INSPC(JPTSPC).EQ.3) THEN
            WRITE(6,'(12X,I3,8X,A)')
     &      JPTSPC, ' Nondiagonal Fock operator '
          ELSE IF(IH0INSPC(JPTSPC).EQ.4) THEN
            WRITE(6,'(12X,I3,8X,A)')
     &      JPTSPC, ' Complete Hamiltonian  '
          ELSE IF(IH0INSPC(JPTSPC).EQ.5) THEN
            WRITE(6,'(12X,I3,8X,A)')
     &      JPTSPC, ' Mix of Fock and Exact operator '
          END IF
         END DO
         IF(ISETKW(61).GT.0) THEN
           WRITE(6,*)
           WRITE(6,'(7X,A)')
     &     ' Orbital subspaces where exact Hamiltonian is used : '
           WRITE(6,'(7X,A)')
     &      '===================================================='
           WRITE(6,*)
           WRITE(6,'(10X,10(2X,I3))') (IH0EXSPC(I),I=1, NH0EXSPC)
           WRITE(6,*)
         END IF
*
       END IF
       END IF
*
*
       IF(NPROP.EQ.0) THEN
       WRITE(6,*)
C        WRITE(6,*) '     No calculation of properties'
       ELSE
         WRITE(6,'(7X,A,I3)')
     &   ' Number of properties to be calculated', NPROP
         WRITE(6,*)
         WRITE(6,'(9X,A)')    ' Properties : '
         WRITE(6,'(9X,A)')   ' ============='
         DO IPROP = 1, NPROP
           WRITE(6,'(16X,A)') PROPER(IPROP)
         END DO
*
c         IF(IRELAX.EQ.0) THEN
           WRITE(6,'(7X,A)') ' No use of relaxed densities '
c         ELSE
c           WRITE(6,'(7X,A)')
c     &     ' Relaxed densities used for property evaluation'
C          WRITE(6,'(7X,A)') ' (implemented only for pert) '
c         END IF
       END IF
*
       IF(IEXTKOP.EQ.0.AND.IPTEKT.EQ.0) THEN
C        WRITE(6,'(5X,A)') ' No extended Koopmans'' calculations '
       ELSE IF(IEXTKOP.NE.0) THEN
         WRITE(6,'(5X,A)') ' Extended Koopmans'' calculations '
       ELSE IF(IPTEKT.NE.0) THEN
         WRITE(6,'(5X,A)') ' Perturbation expansion of EKT equations'
       END IF
*
       IF(IPTFOCK.EQ.1) THEN
         WRITE(6,*) ' Perturbation expansion of Fock matrix '
       ELSE
C        WRITE(6,*) 'No  Perturbation expansion of Fock matrix '
       END IF
*
      IF(ITRAPRP.EQ.0) THEN
C       WRITE(6,*)
C       WRITE(6,'(5X,A)')
C    &  ' No transition properties will be calculated'
      ELSE
        WRITE(6,*)
        WRITE(6,'(5X,A)')
     &  ' Transition properties will be calculated '
        WRITE(6,*)  ' Symmetry of additional states :', IEXCSYM
        WRITE(6,*)  ' Number   of additional states :', NEXCSTATE
        WRITE(6,*)
      END IF
*
*
      IF(IRESPONS.NE.0) THEN
      WRITE(6,*)
      WRITE(6,*) '**************************'
      WRITE(6,*) '*  Response Calculation  *'
      WRITE(6,*) '************************** '
      WRITE(6,*)
        WRITE(6,*)
     &  ' CI-Response will be called after each CI calculation'
        WRITE(6,*)
     &  ' Root used for response calculations (RFROOT) ',IRFROOT
        WRITE(6,*)
C       WRITE(6,*) ' Number of A-operators : ', N_AVE_OP
        WRITE(6,*) ' Labels of A-operators '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        DO IAVE = 1, N_AVE_OP
          WRITE(6,'(1H , 6X,A)') AVE_OP(IAVE)
        END DO
        WRITE(6,*)
C       WRITE(6,*) ' Number of response calculations ', NRESP
        WRITE(6,*) ' Perturbations : '
        WRITE(6,*) ' ================'
        WRITE(6,*)
        WRITE(6,*) ' Calc  Op1    Op2    Mxord1     Mxord2    Freq '
        DO IRESP = 1, NRESP
          WRITE(6,'(1H ,I2,2X,A,A,3X,I4,3X,I4,2X,F12.7)' )
     &    IRESP,RESP_OP(1,IRESP),RESP_OP(2,IRESP),MAXORD_OP(1,IRESP),
     &    MAXORD_OP(2,IRESP),RESP_W(IRESP)
        END DO
      END IF
*
C     IF(NOMOFL.EQ.0) THEN
        WRITE(6,*)
        WRITE(6,'(7X,A)') ' Final orbitals :'
        WRITE(6,'(7X,A)') ' ================'
        WRITE(6,*)
        IF(IFINMO.EQ.1) THEN
          WRITE(6,'(10X,A)') ' Natural orbitals'
        ELSE IF(IFINMO.EQ.2) THEN
          WRITE(6,'(10X,A)') ' Canonical orbitals'
        ELSE IF(IFINMO.EQ.3) THEN
          WRITE(6,'(10X,A)') ' Pseudo-natural orbitals'
          WRITE(6,'(10X,A)')
     &   ' (Density matrix diagonalized in orbital subspaces )'
        ELSE IF(IFINMO.EQ.4) THEN
          WRITE(6,'(10X,A)') ' Pseudo-canonical orbitals'
          WRITE(6,'(10X,A)')
     &   ' (FI+FA  diagonalized in orbital subspaces )'
         ELSE IF (IFINMO .EQ. 5 ) THEN
          WRITE(6,'(10X,A)')
     &   ' Pseudo-natural-canonical orbitals (sic)'
          WRITE(6,'(10X,A)')
     &   ' (Pseudo natural orbitals are first obtained'
          WRITE(6,'(10X,A)')
     &   '  by diagonalizing density matrix in orbital subpspaces.'
          WRITE(6,'(10X,A)')
     &   '  FI+FA is transformed to this basis, and the transformed'
          WRITE(6,'(10X,A)')
     &   '  matrix is block diagonalized) '
          WRITE(6,*)
          WRITE(6,'(10X,A)')
     &   ' Orbital spaces in which transformed FIFA is diagonalized'
          WRITE(6,'(10X,A)')
     &   ' ========================================================'
          DO IPSSPC = 1, NPSSPC
            WRITE(6,'(A,I2,A,10I4,6X,2I6)')
     &      '     SPACE',IPSSPC,'          ',
     &     (NPSSH(IRREP,IPSSPC),IRREP = 1, NIRREP)
          END DO
        END IF
C     END IF
*. Transformation of CI vectors
      IF(ITRACI.EQ.0) THEN
C       WRITE(6,'(5X,A)')  ' No transformation of CI vectors'
      ELSE
        WRITE(6,'(5X,A)')   ' CI vectors transformed in each run'
        WRITE(6,'(7X,A,A)')
     &        ' Complete or restricted rotations :',ITRACI_CR
        WRITE(6,'(7X,A,A)')
     &        ' Type of Final orbitals           :',ITRACI_CN
      END IF
*
* Integral Transformations
*
c      IF(ITRA_FI.EQ.1) THEN
c        WRITE(6,*) '      Integrals transformed to final MO''s '
c      END IF
      IF(ITRA_IN.EQ.1) THEN
        WRITE(6,*) '      Integrals transformed to initial  MO''s '
      END IF

*
*. Print levels
*
*
      WRITE(6,*)
      WRITE(6,'(1H ,A,E18.9)') '      Core energy : ', ECORE
*
c      IF(IDMPIN.EQ.1) THEN
c        WRITE(6,'(1H ,A)')
c        WRITE(6,*) '      Integrals written in formatted form (E22.15)'
c        WRITE(6,*) '      on file 90 '
c      END IF
*
      WRITE(6,*) ' IPART before leaving READIN = ', IPART
      END IF

C ================================================
C  Set ratio beteeen real and integer word length
C ================================================
      irat = rtoi_molcas
*
C =============================================
C  Find largest unused vector for use in Lucia
C =============================================
*
* Allocate memory for lucia
* Set up string info

      CALL lucia()
*
      return
      end

      SUBROUTINE COMBINATIONS(ICOMBI,SIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "cstate.fh"
*
      ICOMBI = 0
      SIGN   = PSSIGN
      IF (PSSIGN .NE. 0.0D0) ICOMBI = 1
*
      RETURN
      END
