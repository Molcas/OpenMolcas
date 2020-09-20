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
* Copyright (C) 2015, Roland Lindh                                     *
************************************************************************
      SubRoutine EMFInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: to compute the electromagnetic field radiation integrals     *
*         using a complex Gauss-Hermite quadrature.                    *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : RecPrt                                                  *
*              CCrtCmp                                                 *
*              CAssmbl                                                 *
*              CVelInt                                                 *
*              CCmbnVe                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemistry - Angstrom,             *
*             University of Uppsala, Sweden. December 2015             *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Logical ABeq(3)
      Integer iStabO(0:7), iDCRT(0:7)
*
      Call EMFInt_Internal(Array)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_integer_array(iStabM)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine EMFInt_Internal(Array)
      Use Iso_C_Binding
      Real*8, Target :: Array(*)
      Complex*16, Pointer :: zAxyz(:),zBxyz(:),zQxyz(:),zVxyz(:)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 195
      iPrint = nPrint(iRout)
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1+nOrdOp) * 2
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1+nOrdOp) * 2
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1+nOrdOp)*(lb+1+nOrdOp) * 2
      If (nOrdOp.eq.1) Then
         ipVxyz = nip
         nip = nip + nZeta*6*(la+1)*(lb+1) * 2
         ipA = nip
         nip = nip + nZeta
         ipB = nip
         nip = nip + nZeta
         ipRes = nip
         nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      Else
         ipVxyz = nip
         ipA = nip
         ipB = nip
         ipRes = nip
         nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      End If
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'EMFInt: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in EMFInt'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In EMFInt: A',' ',A,1,3)
         Call RecPrt(' In EMFInt: RB',' ',RB,1,3)
         Call RecPrt(' In EMFInt: KVector',' ',CCoor,1,3)
         Call RecPrt(' In EMFInt: P',' ',P,nZeta,3)
         Write (6,*) ' In EMFInt: la,lb=',la,lb
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*     Compute the cartesian values of the basis functions angular part
*     Note that these arrays are complex.
*
      Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,
     &                  [nZeta*3*nHer*(la+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,A,zAxyz,
     &               la+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,CCoor)
      Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,
     &                  [nZeta*3*nHer*(lb+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,RB,zBxyz,
     &               lb+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,CCoor)
      Nullify(zAxyz,zBxyz)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
      Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,
     &                  [nZeta*3*nHer*(la+nOrdOp+1)])
      Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,
     &                  [nZeta*3*(la+nOrdOp+1)*(lb+nOrdOp+1)])
      Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,
     &                  [nZeta*3*nHer*(lb+nOrdOp+1)])
      Call CAssmbl(zQxyz,
     &             zAxyz,la+nOrdOp,
     &             zBxyz,lb+nOrdOp,
     &             nZeta,HerW(iHerW(nHer)),nHer)
      Nullify(zAxyz,zBxyz,zQxyz)
*
*     Compute the cartesian components for the velocity integrals.
*     The velocity components are linear combinations of overlap
*     components.
*
      If (nOrdOp.eq.1) Then
         ipAOff = ipA
         Do iBeta = 1, nBeta
            call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
            ipAOff = ipAOff + nAlpha
         End Do

         ipBOff = ipB
         Do iAlpha = 1, nAlpha
            call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
            ipBOff = ipBOff + 1
         End Do
*
         Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,
     &                     [nZeta*3*(la+1)*(lb+1)])
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,
     &                     [nZeta*3*(la+1)*(lb+1)])
         Call CVelInt(zVxyz,zQxyz,la,lb,
     &                Array(ipA),Array(ipB),nZeta)
         Nullify(zVxyz,zQxyz)
*
*     Combine the cartesian components to the full one electron
*     integral.
*
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,
     &                     [nZeta*3*(la+1)*(lb+1)])
         Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,
     &                     [nZeta*3*(la+1)*(lb+1)*2])
         Call CCmbnVe(zQxyz,nZeta,la,lb,Zeta,rKappa,
     &                Array(ipRes),nComp,zVxyz,CCoor)
         Nullify(zQxyz,zVxyz)
      Else
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,
     &                     [nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)])
         Call CCmbnMP(zQxyz,nZeta,la,lb,nOrdOp,Zeta,
     &                rKappa,Array(ipRes),nComp,CCoor)
         Nullify(zQxyz)
      End If
*
      llOper=lOper(1)
      Do iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
      End Do
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
      Do lDCRT = 0, nDCRT-1
*
*--------Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
      Return
      End SubRoutine EMFInt_internal
*
      End
