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

private :: v2

integer, parameter :: MaxSym = 8, MxTCVx = 7
integer IAD2M(3,36*36)
integer NSYMZ, NORBZ(MaxSym), NOSHZ(MaxSym), LUINTMZ
integer nSym, nBas(MaxSym), nFro(MaxSym), nDel(MaxSym)
integer nOrb(MaxSym), nIsh(MaxSym), nAsh(MaxSym), nOsh(MaxSym), nSsh(MaxSym)
integer NumCho(MaxSym)
logical DoTCVA, DoFull, DoCoul, DoExc2
logical SubBlocks(3,3)
logical IfTest

logical TCVXist(MxTCVx,MaxSym,MaxSym)

type v2
  real*8, allocatable :: A(:,:)
end type v2

type(v2) :: TCVX(MxTCVx,MaxSym,MaxSym)

end module Cho_Tra
