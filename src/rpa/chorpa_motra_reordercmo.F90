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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,nFro,nOcc,      &
     &                                   nVir,nDel,CMOinp,CMOout)
      Implicit None
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
      Integer nFro(nSym)
      Integer nOcc(nSym)
      Integer nVir(nSym)
      Integer nDel(nSym)
      Real*8  CMOinp(*)
      Real*8  CMOout(*)

      Integer iSym, ip1, ip2, l

      ip1=1
      ip2=1
      Do iSym=1,nSym
         l=nBas(iSym)*nOrb(iSym)
         Call dCopy_(l,CMOinp(ip1),1,CMOout(ip2),1)
         Call fZero(CMOout(ip2+l),nBas(iSym)*nBas(iSym)-l)
         ip1=ip1+l
         ip2=ip2+nBas(iSym)*nBas(iSym)
      End Do

! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(nDel)
         Call Unused_integer_array(nFro)
         Call Unused_integer_array(nOcc)
         Call Unused_integer_array(nVir)
      End If
      End
