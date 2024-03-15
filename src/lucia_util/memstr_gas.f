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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE MEMSTR_GAS()
      use stdalloc, only: mma_allocate
      use strbas
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in Module STRBAS
*
* Jeppe Olsen , Winter of 1994
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
*. Start of string information
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL mma_allocate(OCSTR(IGRP)%I,LSTRIN,Label='OCSTR()')
        CALL mma_allocate(STREO(IGRP)%I,NSTRIN,Label='STREO()')
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL mma_allocate(NSTSGP(1)%I,NSMST*NGRP,Label='NSTSGP(1)')
      CALL mma_allocate(ISTSGP(1)%I,NSMST*NGRP,Label='ISTSGP(1)')
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL mma_allocate(NSTSO(ITP)%I,NSPGPFTP(ITP)*NSMST,
     &                    Label='NSTSO(ITP)')
        CALL mma_allocate(ISTSO(ITP)%I,NSPGPFTP(ITP)*NSMST,
     &                    Label='ISTSO(ITP)')
      END DO
*
**. Lexical adressing of arrays : use array indices for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL mma_allocate(Zmat(IGRP)%I,NACOB*NELFGP(IGRP),
     &                    Label='ZMat()')
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. IF creation is involve : Use full orbital notation
*  If only annihilation is involved, compact form will be used
        LENGTH=1
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN
        ELSE IF(ISTAC(IGRP,1).NE.0) THEN
*. Only annihilation map so
          LENGTH = IEL*ISTRIN
        ENDIF
        CALL mma_allocate(STSTM(IGRP,1)%I,LENGTH,LABEL='STSTM(IGRP,1)')
        CALL mma_allocate(STSTM(IGRP,2)%I,LENGTH,LABEL='STSTM(IGRP,2)')
      END DO
*
*. Occupation classes
*
      CALL mma_allocate(IOCLS,NMXOCCLS*NGAS,Label='IOCLS')
*. Annihilation/Creation map of supergroup types
      CALL mma_allocate(SPGPAN,NTSPGP*NGAS,Label='SPGPAN')
      CALL mma_allocate(SPGPCR,NTSPGP*NGAS,Label='SPGPCR')
*
      END
