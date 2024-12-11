!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module StrNum
use Definitions, only: iwp
Implicit None
Private
LOGICAL EQUAL
INTEGER(KIND=iwp), PARAMETER :: MXOTYP=100
INTEGER(KIND=iwp) MAXSYM,NORB1,NORB2,NORB3,NAEL,NBEL,NAEXCI,NBEXCI,NASTR,NBSTR,   &
                  NOCTPA,NOCTPB,NL1MNA,NL1MNB,                                    &
                  IOCPTA(MXOTYP),IOCPTB(MXOTYP),                                  &
                  NSTAOA(MXOTYP),NSTAOB(MXOTYP)

Public :: EQUAL,MAXSYM,NORB1,NORB2,NORB3,NAEL,NBEL,NAEXCI,NBEXCI,NASTR,NBSTR,     &
          NOCTPA,NOCTPB,NL1MNA,NL1MNB,IOCPTA,IOCPTB,NSTAOA,NSTAOB
End Module StrNum
