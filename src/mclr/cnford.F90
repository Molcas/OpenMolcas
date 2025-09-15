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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

subroutine CNFORD(ICTSDT,ICONF,IREFSM,NORB,IPRODT,NCNFTP,NEL,IGENSG,ISGNA,ISGNB,IAGRP,IBGRP,IOOS,NORB1,NORB2,NORB3,NEL1MN,NEL3MX, &
                  NAEL,NBEL,MINOP,MAXOP,PSSIGN)
! Generate configurations in ICONF
!
! Generate determinants in configuration order and obtain
! sign array for switching between the two formats.
!
! Jeppe Olsen January 1989
!
! December 1990 : ICNFOK added
! September 1993 : Combinations added
!
! NCNFCN /= 0 indicates that additional constraints on configurations
! should be checked
! by calling CICNCH.ICNFOK(ICNF) is 1 of tests are passed, ICNFOK(ICNF)
! is zero if test fails
! IGENSG /= 0 assumes general signs of strings given in ISGNA,ISGNB

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: ICTSDT(*), ICONF(*)
integer(kind=iwp), intent(in) :: iRefSM, nOrb, IPRODT(*), NCNFTP(*), IGENSG, ISGNA(*), ISGNB(*), IAGRP, IBGRP, IOOS(*), NORB1, &
                                 NORB2, NORB3, NEL1MN, NEL3MX, NAEL, NBEL, MINOP, MAXOP
integer(kind=iwp), intent(inout) :: NEL
real(kind=wp) :: PSSIGN

! NOTE : NCNFTP IS COLUMN FOR SYMMETRY GIVEN, NOT COMPLETE MATRIX.
! Dim of IWORK : MAX(3*NORB,(MXDT+2)*NEL),
! where MXDT is the largest number of prototype determinants of
! a given type.
!
! ==============================================================
! Construct list of configurations,offset for each configuration
!  and type for each configuration
! ==============================================================

call CONFG2(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,MINOP,MAXOP,IREFSM,NEL,ICONF,NCNFTP)

! ========================================================
! Obtain determinants for each configuration and determine
! the corresponding address and phaseshift to reform into
! string form and ordering.
! ========================================================

call CNTOST(ICONF,ICTSDT,NAEL,NBEL,IPRODT,IREFSM,NORB,NEL,IGENSG,ISGNA,ISGNB,IAGRP,IBGRP,IOOS,PSSIGN)

end subroutine CNFORD
