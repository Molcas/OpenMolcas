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
Module Cho_Tra

Private:: v2

Integer, Parameter:: MaxSym=8, MxTCVx=7
Integer IAD2M(3,36*36)
Integer NSYMZ,NORBZ(MaxSym),NOSHZ(MaxSym),LUINTMZ
Integer nSym, nBas(MaxSym), nFro(MaxSym), nDel(MaxSym)
Integer nOrb(MaxSym),nIsh(MaxSym),nAsh(MaxSym),nOsh(MaxSym), nSsh(MaxSym)
Integer NumCho(MaxSym)
Logical DoTCVA, DoFull, DoCoul, DoExc2
Logical SubBlocks(3,3)
Logical IfTest

Integer iMemTCVX(MxTCVx,MaxSym,MaxSym,2)
Logical  TCVXist(MxTCVx,MaxSym,MaxSym)

Type v2
  Real*8, Allocatable:: A(:,:)
End Type v2

Type (v2):: TCVX(MxTCVx,MaxSym,MaxSym)

End Module Cho_Tra
