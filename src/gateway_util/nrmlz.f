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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Nrmlz(Exp,nPrim,Coeff,nCntrc,iAng)
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      Real*8 Exp(nPrim), Coeff(nPrim,nCntrc)
      Real*8, Dimension(:), Allocatable :: Scrt1, Scrt2
*
      If (nPrim*nCntrc.eq.0) Return
*
      nScrt1=nPrim**2
      nScrt2=nPrim*nCntrc
      Call mma_allocate(Scrt1,nScrt1)
      Call mma_allocate(Scrt2,nScrt2)
      Call Nrmlz_Internal(Exp,nPrim,Coeff,nCntrc,Scrt1,nScrt1,
     &                                           Scrt2,nScrt2,iAng)
      Call mma_deallocate(Scrt2)
      Call mma_deallocate(Scrt1)
*
*
      Return
      End
      SubRoutine Nrmlz_Internal(Exp,nPrim,Coeff,nCntrc,
     &                          Scrt1,nScrt1,Scrt2,nScrt2,iAng)
************************************************************************
*                                                                      *
* Object: normalize the contraction coefficients with respect to the   *
*         radial overlap.                                              *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : RecPrt                                                  *
*              DGEMM_  (ESSL)                                          *
*              DnDot   (ESSL)                                          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Exp(nPrim), Coeff(nPrim,nCntrc), Scrt1(nScrt1),
     &       Scrt2(nScrt2)
#include "real.fh"
*
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) ' In Nrmlz: iAng=',iAng
      Call RecPrt(' In Nrmlz: Coefficients',' ',Coeff,nPrim,nCntrc)
      Call RecPrt(' In Nrmlz: Exponents',' ',Exp,nPrim,1)
#endif
*
*     Normalize the coefficients (only radial normalization)
*
*     Compute overlap for all primitives of this shell.
*     This formula includes the normalization constant for each
*     primitive as well as the overlap factor of between the primitives.
*     Hence, the overlap matrix elements correspond to those of the
*     normalized primitive gaussian functions.
*
      Do 200 iExp = 1, nPrim
         Do 210 jExp = 1, iExp-1
               Pro_ij  = Exp(iExp)*Exp(jExp)
               Sum_ij  = Exp(iExp)+Exp(jExp)
               iPower_i1=iAng + 1
               Power= 0.5D0*DBLE(iAng) + 0.75D0
               C       = 2.0D0/Sum_ij
*              Temp    = Sqrt(C) * C**iPower_i1 * Pro_ij**Power
               Power = DBLE(iAng) + 1.5D0
               C = 2.0D0*Sqrt(Pro_ij)/Sum_ij
               Temp = C**Power
*              Temp =  ( Two*Sqrt(Exp(iExp)*Exp(jExp)) /
*    &                   (Exp(iExp)+Exp(jExp)))**(DBLE(iAng)+Three/Two)
               Scrt1(nPrim*(iExp-1)+jExp)=Temp
               Scrt1(nPrim*(jExp-1)+iExp)=Temp
210        Continue
           Scrt1(nPrim*(iExp-1)+iExp)=One
200     Continue
*     Contract right side
      Call DGEMM_('N','N',
     &            nPrim,nCntrc,nPrim,
     &            1.0d0,Scrt1,nPrim,
     &            Coeff,nPrim,
     &            0.0d0,Scrt2,nPrim)
#ifdef _DEBUGPRINT_
      Call RecPrt(' Overlap primitives',' ',Scrt1,nPrim,nPrim)
      Call RecPrt(' Overlap PrimCon',' ',Scrt2,nPrim,nCntrc)
#endif
*
*     Compute the overlap for each contracted basis function, <i|i>
*
      Call DnDot(nCntrc,nPrim,Scrt1,1,1,Scrt2,1,nPrim,Coeff,1,nPrim)
#ifdef _DEBUGPRINT_
      Call RecPrt(' Overlap Contracted',' ',Scrt1,nCntrc,1)
#endif
*
*     Normalize coefficients, i.e. combine the normalization factor
*     of the primitive and the overlap of the unnormalized contracted
*     basis function.

      Do 300 i = 1, nCntrc
         If ( Abs(Scrt1(i)).lt.1.0D-12 ) then
            Call WarningMessage(2,
     &                  '; Error in contraction matrix, zero column'
     &                //'; ; Abend in subroutine NRMLZ')
            Call Abend()
         End If
300   Continue
*
      Rtemp=0.5D0*DBLE(iAng) + 0.75D0
      Qtemp=2.0D0**(iAng+1) * Sqrt(2.0D0) * TwoP34
      Do i = 1, nCntrc
         vRR= Scrt1(i)**(-0.5d0)
         Do j = 1, nPrim
            vR2=Exp(j)**Rtemp
            Coeff(j,i) = Coeff(j,i)*Qtemp*vRR*vR2
         End Do
      End Do
      If (nPrim.eq.1 .and. nCntrc.eq.1 .and. Exp(1).eq.Zero) Then
         Coeff(1,1)=One
      End If
#ifdef _DEBUGPRINT_
      Call Recprt(' In Nrmlz: Normalized coefficients',' ',
     &            Coeff,nPrim,nCntrc)
#endif
*
      Return
      End
