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
      Use Fock_util_global, only: DoCholesky
      Use Cholesky, only: ChFracMem
      Use KSDFT_Info, Only: CoefR, CoefX
      use mcpdft_output, only:  set_print_level

      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general_mul.fh"
#include "gas.fh"
#include "timers.fh"
#include "lucia_ini.fh"
#include "WrkSpc.fh"

*----------------------------------------------------------------------*

* Set print levels, and adjust them if needed:
      call set_print_level()
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

#if defined (_MOLCAS_MPP_)
      ChFracMem=0.3d0
#else
      ChFracMem=0.0d0
#endif

* Default title line:
      TITLE(1)='(No title given)'
*
* assign ipCleanMask to dummy pointer
      ipCleanMask=ip_Dummy
*
* iteration control
*
* max number of super-CI iterations
      ITMAX=mxSxIt
* threshold for change in RASSCF energy
      THRE=1.D-08

* Choose to only expand or generate information for CI-vectors if INOCALC = 1
      INOCALC = 0
* Save information on CI expansion if ISAVE_EXP = 1
      ISAVE_EXP = 0
* Expand a smaller CI vector in a larger one if IEXPAND = 1
      IEXPAND = 0
*
* wave function control bits
*

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

      ! prethr: energy threshold for printout of orbitals
      prethr = 0.15d0

*
* Default value for type of CASSCF (used for DFT)
*
      KSDFT='SCF'
      ExFac=1.0D0
* Initialize KSDF coefficients (S Dong, 2018)
      CoefR = 1.0D0
      CoefX = 1.0D0

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

* These keys will activate the calculation of the high
* frequency contribution to the reaction field
* ???
* This key controls if a non-equilibrium reaction field
* calculation is performed.
      NonEq=.False.

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

      END
