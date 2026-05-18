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
! Copyright (C) 2000, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2000  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine GETSGM2(ILEV,JLEV,ISYCI,CI,nCI,SGM,MSGM)
! GIVEN CI COUPLING LEVELS ILEV, JLEV, COMPUTE SGM=E(ILEV,JLEV)*CI
! ILEV,JLEV ARE IN PRINCIPLE ACTIVE ORBITAL NUMBERS, BUT POSSIBLY
! IN ANOTHER ORDER THAN THE USUAL ONE -- HERE WE USE THE ORDER
! FOLLOWED BY THE GUGA COUPLING SCHEME.
!
! THIS ROUTINE REPLACES EARLIER GETSGM, TO GET RID OF THE PACKING AND
! STORING USED EARLIER.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE!! THE EARLIER CALL GETSGM(ILEV,JLEV,IDARR,SGM) IS REPLACED BY
! GETSGM2(ILEV,JLEV,CI,SGM)!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use Symmetry_Info, only: Mul
use sguga, only: SGS, CIS, EXS
use constants, only: Zero, One
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(in) :: ILEV, JLEV, ISYCI, nCI, MSGM
real(kind=wp), intent(in) :: CI(nCI)
real(kind=wp), intent(inout) :: SGM(MSGM)
integer(kind=iwp) IS, JS, IJS, ISSG, NSGM

SGM(1:MSGM) = Zero
IS = SGS%ISM(ILEV)
JS = SGS%ISM(JLEV)
IJS = Mul(IS,JS)
ISSG = Mul(IJS,ISYCI)
NSGM = CIS%NCSF(ISSG)
if (NSGM > MSGM) then
  write(u6,*) 'GETSGM2: NSGM>MSGM'
  call Abend()
end if
if (NSGM == 0) return

call SG_Epq_Psi(SGS,CIS,EXS,ILEV,JLEV,One,ISYCI,CI,SGM)

end subroutine GETSGM2
