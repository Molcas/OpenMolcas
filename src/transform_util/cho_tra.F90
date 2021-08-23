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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
!***********************************************************************
!         MODULE for TRANSFORMATION of CHOLESKY VECTORS                *
!              and GENERATION of TWO-ELECTRONS INTEGRALS               *
!***********************************************************************

module Cho_Tra

use Data_Structures, only: Alloc2DArray_Type
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MaxSym = 8, MxTCVx = 7
integer(kind=iwp) :: IAD2M(3,36*36), nAsh(MaxSym), nBas(MaxSym), nDel(MaxSym), nFro(MaxSym), nIsh(MaxSym), nOrb(MaxSym), &
                     nOsh(MaxSym), nSsh(MaxSym), nSym, NumCho(MaxSym)
logical(kind=iwp) :: DoCoul, DoExc2, DoFull, DoTCVA, IfTest, SubBlocks(3,3), TCVXist(MxTCVx,MaxSym,MaxSym)
type(Alloc2DArray_Type) :: TCVX(MxTCVx,MaxSym,MaxSym)

public :: DoCoul, DoExc2, DoFull, DoTCVA, IAD2M, IfTest, nAsh, nBas, nDel, nFro, nIsh, nOrb, nOsh, nSsh, nSym, NumCho, SubBlocks, &
          TCVX, TCVXist

end module Cho_Tra
