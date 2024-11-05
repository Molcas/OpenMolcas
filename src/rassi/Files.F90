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
Module Files
use Cntrl, only: MXJOB
!----------------------------------------------------------------------*
!     Define files ( file names and unit numbers )                     *
!----------------------------------------------------------------------*
!
! LUIPH  - UNIT NUMBER OF JOBIPHS
! LUMCK  - UNIT NUMBER OF MCKINT FILES
! LUONE  - D:O, ONE-ELECTRON INTEGRAL FILE
! LUORD  - D:O, ORDERED TWO-ELECTRON INTEGRAL FILE
! IADR15 - TABLE OF CONTENTS, DISK ADDRESSES ON LUIPH.
! IDCMO  - Adresses to the CMO arrays on each JOBIPH
Character(LEN=8) FnOne,FnIph,FnMck,FnExt,FnOrd,FnCom,FnTDM,FnExc,FnToF,FnToM,FnEig
INTEGER LUCOM, LUEXT, LUIPH, LUMCK, LUONE, LUORD, LUEIG
INTEGER LUEXC, LUTDM, LUTOM, IDCMO(MXJOB), lIDTDM, ITOC15(30), LUDYS, LIDDYS
End Module Files
