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

use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer, intent(out) :: ICTSDT(*)
integer, intent(inout) :: ICONF(*)
integer :: iRefSM, nOrb
integer, intent(in) :: IPRODT(*)
integer, intent(in) :: NCNFTP(*)
integer :: NEL, IGENSG
integer, intent(in) :: ISGNA(*), ISGNB(*)
integer :: IAGRP, IBGRP
integer, intent(in) :: IOOS(*)
integer :: NORB1, NORB2, NORB3, NEL1MN, NEL3MX, NAEL, NBEL, MINOP, MAXOP
real*8 :: PSSIGN
! Scratch
integer, allocatable :: KL1(:), KL2(:), KL3(:)

! NOTE : NCNFTP IS COLUMN FOR SYMMETRY GIVEN, NOT COMPLETE MATRIX.
! Dim of IWORK : MAX(3*NORB,(MXDT+2)*NEL),
! where MXDT is the largest number of prototype determinants of
! a given type.
!
! ==============================================================
! Construct list of configurations,offset for each configuration
!  and type for each configuration
! ==============================================================

call mma_allocate(KL1,NORB1+NORB2+NORB3,Label='KL1')
call mma_allocate(KL2,NORB1+NORB2+NORB3,Label='KL2')
call mma_allocate(KL3,NORB1+NORB2+NORB3,Label='KL3')
call CONFG2(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,MINOP,MAXOP,IREFSM,NEL,ICONF,NCNFTP,KL1,KL2,KL3)
call mma_deallocate(KL3)
call mma_deallocate(KL2)
call mma_deallocate(KL1)

! ========================================================
! Obtain determinants for each configuration and determine
! the corresponding address and phaseshift to reform into
! string form and ordering.
! ========================================================

call CNTOST(ICONF,ICTSDT,NAEL,NBEL,IPRODT,IREFSM,NORB,NEL,IGENSG,ISGNA,ISGNB,IAGRP,IBGRP,IOOS,PSSIGN)

end subroutine CNFORD
