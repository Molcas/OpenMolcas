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
* Copyright (C) 1996, Per Ake Malmqvist                                *
*               1996, Roland Lindh                                     *
************************************************************************
      SubRoutine AMPInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for computing matrix elements of the          *
*         six hermitized products of two angular momentum ops          *
*                                                                      *
*     Author: Per-Ake Malmqvist, Dept. of Theoretical Chemistry,       *
*             University of Lund, SWEDEN                               *
*             November '96                                             *
*     After pattern of other SEWARD soubroutines by R. Lindh.          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local Variables
      Real*8 TC(3)
      Integer iStabO(0:7), iDCRT(0:7)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2
*
      iRout = 220
      iPrint = nPrint(iRout)
*
      nip = 1
      ipB = nip
      nip = nip + nZeta
      ipTpp = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+2)*6
      ipTp  = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*3
      ipT   = nip
      nip = nip + nZeta*nElem(la)*nElem(lb  )*6
      ipTm  = 1
      ipTmm = 1
      if(lb.gt.0) then
        ipTm  = nip
        nip = nip + nZeta*nElem(la)*nElem(lb-1)*3
        if(lb.gt.1) then
          ipTmm = nip
          nip = nip + nZeta*nElem(la)*nElem(lb-2)*6
        end if
      end if
      ipRes=nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,' AMPInt: nip-1.gt.nZeta*nArr')
         call Abend()
      End If
      ipArr = nip
      mArr = (nArr*nZeta - (nip-1))/nZeta

      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)

      ipOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
      End Do

      llOper = lOper(1)
      Do iComp = 2, nComp
         iDum=lOper(iComp)
         llOper = iOr(llOper,iDum)
      End Do

C Compute stabilizer, and then the double coset representation:
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

C Loop over the cosets of the stabilizer group:
      Do lDCRT = 0, nDCRT-1
         Call OA(iDCRT(lDCRT), Ccoor, TC)


C Generate the quadrupole integral tables:
         iComp=6
         iOrdOp = 2
         nHer = (la + (lb+2) + 2 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipTpp),nZeta,iComp,la,lb+2,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         nHer = (la +  lb    + 2 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipT  ),nZeta,iComp,la,lb  ,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         if(lb.ge.2) then
           nHer = (la + (lb-2) + 2 + 2) / 2
           Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &               Array(ipTmm),nZeta,iComp,la,lb-2,A,RB,nHer,
     &               Array(ipArr),mArr,TC,iOrdOp)
         end if
C Generate the dipole integral tables:
         iComp=3
         iOrdOp = 1
         nHer = (la + (lb+1) + 1 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipTp ),nZeta,iComp,la,lb+1,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         if(lb.ge.1) then
           nHer = (la + (lb-1) + 1 + 2) / 2
           Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &               Array(ipTm ),nZeta,iComp,la,lb-1,A,RB,nHer,
     &               Array(ipArr),mArr,TC,iOrdOp)
         end if

         if(iprint.gt.49) write(6,*)' AMPInt calling AMPr.'
         Call AMPr(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipTpp),
     &             Array(ipTp),Array(ipT),Array(ipTm),Array(ipTmm))

C Symmetry adaption:
         if(iprint.gt.49) write(6,*)' AMPInt calling SymAdO'
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp,lOper,iChO,One)
         if(iprint.gt.49) write(6,*)' Back to AMPInt.'
      End Do

      if(iprint.gt.49) write(6,*)' Leaving AMPInt.'
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nOrdOp)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(iAddPot)
      End If
      End
