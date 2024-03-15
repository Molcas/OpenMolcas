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
      SUBROUTINE FREESTR_GAS()
      use stdalloc, only: mma_deallocate
      use strbas
* allocations during strinf_gas
      use distsym
* Deallocate the memory that was set up in MEMSTR_GAS

*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "csm.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "stinf.fh"
#include "crun.fh"
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        Call mma_deallocate(OCSTR(IGRP)%I)
        CALL mma_deallocate(STREO(IGRP)%I)
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL mma_deallocate(NSTSGP(1)%I)
      CALL mma_deallocate(ISTSGP(1)%I)
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL mma_deallocate(NSTSO(ITP)%I)
        CALL mma_deallocate(ISTSO(ITP)%I)
      END DO
*
**. Lexical adressing of arrays : use array indices for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL mma_deallocate(Zmat(IGRP)%I)
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
*. IF creation is involve : Use full orbital notation
*  If only annihilation is involved, compact form will be used
        CALL mma_deallocate(STSTM(IGRP,1)%I)
        CALL mma_deallocate(STSTM(IGRP,2)%I)
      END DO
*
*. Occupation classes
*
      CALL mma_deallocate(IOCLS)
*. Annihilation/Creation map of supergroup types
      CALL mma_deallocate(SPGPAN)
      CALL mma_deallocate(SPGPCR)
*
* Allocated during strinf_gas call
      CALL mma_deallocate(ISMDFGP)
      CALL mma_deallocate(NACTSYM)
      CALL mma_deallocate(ISMSCR)
      END
