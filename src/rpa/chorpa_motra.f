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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine ChoRPA_MOTra(includeFrozen,includeDeleted)
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Transform Cholesky vectors to MO basis.
C
C     TODO/FIXME:
C     1. This routine computes all MO blocks (ij, ai, ab), even though
C        we may only need some of them. A more flexible interface would
C        be nice to have. Presumably not a performance issue in RPA,
C        though (remains to be verified).
C     2. For unrestricted calculations, the alpha and beta
C        transformations are done separately, which means that the AO
C        vectors are read twice. Simultaneous transformation would be
C        desirable!
C
      Implicit None
      Logical includeFrozen, includeDeleted
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam='ChoRPA_MOTra')

      Integer  RPA_iUHF
      External RPA_iUHF

      Character*6 BName(2)

      Integer nSpin
      Integer iSym
      Integer ip_lCMO, l_lCMO
      Integer ip, l
*     Integer ip_lnBas
*     Integer ip_lnOrb
      Integer ip_lnFro
      Integer ip_lnOcc
      Integer ip_Zeros
      Integer ip_lnVir
      Integer ip_lnDel
      Integer iSpin

      nSpin=RPA_iUHF()
      If (nSpin.eq.1) Then
         BName(1)='MOVECS'
         BName(2)='unused'
      Else If (nSpin.eq.2) Then
         BName(1)='MOVECa'
         BName(2)='MOVECb'
      Else
         Call RPA_Warn(3,SecNam//': illegal nSpin')
         BName(1)='unused'
         BName(2)='unused'
      End If
      l_lCMO=nBas(1)**2
      Do iSym=2,nSym
         l_lCMO=l_lCMO+nBas(iSym)**2
      End Do
      Call GetMem('locCMO','Allo','Real',ip_lCMO,l_lCMO)
      l=5*nSym
      Call GetMem('local','Allo','Inte',ip,l)
      ip_lnFro=ip
      ip_lnOcc=ip+nSym
      ip_Zeros=ip+2*nSym
      ip_lnVir=ip+3*nSym
      ip_lnDel=ip+4*nSym
      Call iZero(iWork(ip_Zeros),nSym)
      If (includeFrozen) ip_lnFro=ip_Zeros
      If (includeDeleted) ip_lnDel=ip_Zeros

      Do iSpin=1,nSpin

         ! Set orbital blocks
         If (includeFrozen) Then
            Do iSym=1,nSym
               iWork(ip_lnOcc-1+iSym)=nFro(iSym,iSpin)+nOcc(iSym,iSpin)
            End Do
         Else
            Call iCopy(nSym,nFro(1,iSpin),1,iWork(ip_lnFro),1)
            Call iCopy(nSym,nOcc(1,iSpin),1,iWork(ip_lnOcc),1)
         End If
         If (includeDeleted) Then
            Do iSym=1,nSym
               iWork(ip_lnVir-1+iSym)=nVir(iSym,iSpin)+nDel(iSym,iSpin)
            End Do
         Else
            Call iCopy(nSym,nVir(1,iSpin),1,iWork(ip_lnVir),1)
            Call iCopy(nSym,nDel(1,iSpin),1,iWork(ip_lnDel),1)
         End If
         ! Reorder CMO array
         Call ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,
     *                                nFro(1,iSpin),nOcc(1,iSpin),
     *                                nVir(1,iSpin),nDel(1,iSpin),
     *                                Work(ip_CMO(iSpin)),Work(ip_lCMO))
         ! Set base name for MO files
         ! Transform Cholesky vectors
         Call Cho_MOTra_(Work(ip_lCMO),l_lCMO,nSym,
     *                   nBas,nOrb,
     *                  iWork(ip_lnFro),iWork(ip_lnOcc),iWork(ip_Zeros),
     *                   iWork(ip_lnVir),iWork(ip_lnDel),
     *                   BName(iSpin),.false.,0,.false.)

      End Do

      Call GetMem('local','Free','Inte',ip,l)
      Call GetMem('locCMO','Free','Real',ip_lCMO,l_lCMO)

      End
************************************************************************
      Subroutine ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,nFro,nOcc,
     *                                   nVir,nDel,CMOinp,CMOout)
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

c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(nDel)
         Call Unused_integer_array(nFro)
         Call Unused_integer_array(nOcc)
         Call Unused_integer_array(nVir)
      End If
      End
