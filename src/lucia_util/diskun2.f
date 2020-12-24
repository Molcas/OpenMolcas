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
      SUBROUTINE DISKUN2
      IMPLICIT REAL*8           (A-H,O-Z)
*
* Assign logical unit numbers for LUCIA:
*
* All file with some kind of input information  :  10 - 19
* All files containing final results            :  90 - 99
* Scratch files                                 :  40 - 59
* Internal files (retained through job)         :  60 - 69
*
#include "clunit.fh"
#include "io_util.fh"
* =========================
* Standard input and output
* =========================
*. Input file
*     LUIN = 5
*. Output file
*     LUOUT = 6
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 12
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals - MOLCAS
      LU2INT = 13
*. Input , one electron integrals - MOLCAS
      LU1INT = 14
*. Input , property one-electron integral files
      LUPRP  = 15
*. Sirius interface file
      LUSIR1 = 16
*. File containing additional states for transition densities
      LUEXC = 17
* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = IsFreeUnit(18)
      CALL DANAME_WA(LUDIA,'CIDIA')
*
*. CI vector
      LUC = IsFreeUnit(Ludia)
      CALL DANAME_WA(LUC,'LUCVECT')
*
*. Sigma vector file
      LUHC = IsFreeUnit(LuC)
      CALL DANAME_WA(LUHC,'HCFILE')
*
*. File collecting CC correction vectors, used for DIIS etc
C-jwk      LU_CCVEC = 23
*. File containing approximations to the CC solutions
C-jwk      LU_CCVECT = 24
*. File containing CC vector functions for the CC vectors on LU_CCVECT
C-jwk      LU_CCVECF = 25
*. File containing Last CC coefficients
C-jwk      LU_CCVECL = 26
*. File containing Last CC vector function
C-jwk      LU_CCVECFL = 27
* =================
* Scratch files
* =================
      LUSC1 = IsFreeUnit(LuHC)
      CALL DANAME_WA(LUSC1,'LUSC1')
*
      LUSC2 = IsFreeUnit(Lusc1)
      CALL DANAME_WA(LUSC2,'LUSC2')
*
      LUSC3 = IsFreeUnit(Lusc2)
      CALL DANAME_WA(LUSC3,'LUSC3')
*
*. Scratch space for subspace handling
      LUSC34 = IsFreeUnit(Lusc3)
      CALL DANAME_WA(LUSC34,'LUSC34')
*
      LUSC35 = IsFreeUnit(Lusc34)
      CALL DANAME_WA(LUSC35,'LUSC35')
*
      LUSC36 = IsFreeUnit(Lusc35)
      CALL DANAME_WA(LUSC36,'LUSC36')
*
      LUSC37 = IsFreeUnit(Lusc36)
      CALL DANAME_WA(LUSC37,'LUSC37')
*
      LUSC38 = IsFreeUnit(Lusc37)
      CALL DANAME_WA(LUSC38,'LUSC38')
*
      LUSC39 = IsFreeUnit(Lusc38)
      CALL DANAME_WA(LUSC39,'LUSC39')
*
      LUSC40 = IsFreeUnit(Lusc39)
      CALL DANAME_WA(LUSC40,'LUSC40')
*
*.
      LUSC51 = 51
      LUSC52 = 52
      LUSC53 = 53

* =================
* Output files
* =================
*. output file for CI-vectors
*. Not in use
      LUCIVO = 98
*. Natural orbitals in terms of input orbitals
*.
      LUMOUT = IsFreeUnit(Lusc40)
      CALL DANAME_WA(LUMOUT,'LUMOUT')
*

*. Dumping 1- and 2- electron integrals in formatted form
* ( LU90 just defined here, it is not exported )
c     LU90  = 90
*. Dumping symmmetry info, MO-AO expansion matrix and property integrals
c     LU91 = 91
*. CC amplitudes in formatted form
c      LU_CCAMP = 92
*. Result of CI=> CC conversion
c      LU_CC_FROM_CI = 93
*. Excitation operators, all symmetries
c      LU_CCEXC_OP = 94
      IDISK=0
      RETURN
      END
