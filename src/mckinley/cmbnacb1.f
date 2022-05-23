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
* Copyright (C) 2004, Anders Bernhardsson                              *
************************************************************************
      Subroutine CmbnACB1(FA1,FB1,Final,Fact,nAlpha,nBeta,C,nC,
     &                    la,lb,iang,ifgrad,tmp,lsro,
     &                    index,mvec,idcar)
*******************************************************************************
*
*     Merges the first derivatives of ECP projection/SRO  integrals
*     for derivatives of components
*
*******************************************************************************
*
*     @parameter FA1   The first derivative of Left side , Includes no deriavtive (input)
*     @parameter FB1   The first derivative of Right side . Includes no derivative (input)
*     @parameter Final Result added up to (out)
*     @parameter Fact  Factor the reult is multiplied with bef. added up (input)
*     @parameter C     Coefficients for SRO (input)
*     @parameter nAlpha Number of exponents LS (input)
*     @parameter nBeta  Number of exponents RS (input)
*     @parameter nC     Number of exponents in SRO (input)
*     @parameter la     Angular monenta LS (input)
*     @parameter lb     Angular monenta RS (input)
*     @parameter iAng   Angular monenta SRO (input)
*     @parameter Tmp    Working Area nAlpha*nC (SRO case) (scratch)
*     @parameter lSRO   true for SRO false projection operator (input)
*     @paramaeter index array storing index for derivatives in final (out)
*     @parameter  mvec  Number of derivatives calculated (out)
*     @parameter  idcar cartesiam index for current derivative (input)
*
*******************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"

      Logical ifgrad(3,4),lsro
      Real*8 FA1(nAlpha,nC,(la+1)*(la+2)/2,(2*iang+1),*),
     &       FB1(nC,nBeta,(2*iang+1),(lb+1)*(lb+2)/2,*),
     &     Final(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,21)

      Real*8  C(*),Tmp(*)
      Integer Index(3,4)

      nelem(ixyz) = (ixyz+1)*(ixyz+2)/2

      nZeta = nAlpha*nBeta
      call dcopy_(nZeta*nElem(la)*nElem(lb)*6,
     &                    [Zero],0,Final,1)
      Call iCopy(12,[0],0,Index,1)

      mVec = 0
      Do iCent = 1, 2
         If (ifGrad(iDCar,iCent)) Then
            mVec = mVec + 1
            Index(iDcar,icent)=mvec
*
            If (iCent.eq.1) Then
               iFa = 2
               iFb = 1
            Else
               iFa = 1
               iFb = 2
            End If
*
            Do  ib = 1, nElem(lb)
               Do  ia = 1, nElem(la)

                  Do iC = 1, (2*iAng+1)
                     if (lsro) Then
                       Call DGEMM_('N','N',
     &                      nAlpha,nC,nC,
     &                      One,FA1(1,1,ia,ic,iFa),nAlpha,
     &                      C,nC,
     &                      Zero,Tmp,nAlpha)
                       Call DGEMM_('N','N',
     &                      nAlpha,nBeta,nC,
     &                      Fact,Tmp,nAlpha,
     &                      FB1(1,1,ic,ib,iFb),nC,
     &                      One,Final(1,ia,ib,mvec),nAlpha)
                     else
                       Call DGEMM_('N','N',
     &                           nAlpha,nBeta,nC,
     &                           Fact,Fa1(1,1,ia,ic,iFa),nAlpha,
     &                           Fb1(1,1,ic,ib,iFb),nC,
     &                           One,Final(1,ia,ib,mvec),
     &                           nAlpha)
                     Endif



                  End DO        ! iC
               End DO           ! iA
            End DO              ! ib

         End If
      End Do ! icent

      Return
      End
