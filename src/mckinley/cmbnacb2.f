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
      Subroutine CmbnACB2(Fa1,Fa2,Fb1,Fb2,Final,Fact,nalpha,nbeta,
     &     C,nC,la,lb,iang,jfhess,Tmp,lsro)
*******************************************************************************
*
*     Merges the second derivatives of ECP projection/SRO  integrals
*     for derivatives of components
*
*******************************************************************************
*
*     @parameter FA1   The first derivative of Left side , Includes no deriavtive
*     @parameter FA2   The second derivative of Left side
*     @parameter FB1   The first derivative of Right side . Includes no derivative
*     @parameter FB2   The second derivative of Right side
*     @parameter Final Result added up to (out)
*     @parameter Fact  Factor the reult is multiplied with bef. added up
*     @parameter C     Coefficients for SRO
*     @parameter nAlpha Number of exponents LS
*     @parameter nBeta  Number of exponents RS
*     @parameter nC     Number of exponents in SRO
*     @parameter la     Angular monenta LS
*     @parameter lb     Angular monenta RS
*     @parameter iAng   Angular monenta SRO
*     @parameter nBeta  Number of exponents RS
*     @parameter nC     Number of exponents in SRO
*     @parameter Tmp    Working Area nAlpha*nC (SRO case)
*     @parameter lSRO   true for SRO false projection operator
*
*******************************************************************************
      Implicit Real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"

      Logical jfhess(4,3,4,3),lsro

      Real*8 FA1(nAlpha,nC,(la+1)*(la+2)/2,(2*iang+1),*),
     &       FA2(nAlpha,nC,(la+1)*(la+2)/2,(2*iang+1),*),
     &       FB1(nC,nBeta,(2*iang+1),(lb+1)*(lb+2)/2,*),
     &       FB2(nC,nBeta,(2*iang+1),(lb+1)*(lb+2)/2,*),
     &       Tmp(*),c(*),
     &     Final(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,21)
*                                                                      *
************************************************************************
*                                                                      *
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*     Merge integrals with one deriavtive on each center

      Do iCar = 1, 3
         Do jCar = 1, 3
            mvec=itri(iCar+3,jcar)
            If (jfHess(2,iCar,1,jcar)) Then
               Do  ib = 1, nElem(lb)
                  Do ia = 1, nElem(la)
                     Do  iC = 1, (2*iAng+1)

                        if (lsro) then
                           call mult_sro(FA1(1,1,ia,ic,icar+1),nAlpha,
     &                                   C,nC,
     &                                   FB1(1,1,ic,ib,jcar+1),nBeta,
     &                                   Fact,Final(1,ia,ib,mVec),
     &                                   Tmp)
                        else
                           Call DGEMM_('N','N',
     &                                nAlpha,nBeta,nC,
     &                                Fact,FA1(1,1,ia,ic,icar+1),nAlpha,
     &                                FB1(1,1,ic,ib,jcar+1),nC,
     &                                One,Final(1,ia,ib,mVec),nAlpha)
                        endif
                     End do
                  End do
               End do
            End If
         End do
      End do
*                                                                      *
************************************************************************
*                                                                      *
*     Merge integrals with both derivative on center A
      Do iCar = 1, 3
         Do jCar = 1, icar
            mvec=itri(iCar,jcar)
            If (jfHess(1,iCar,1,jcar)) Then
               Do  ib = 1, nElem(lb)
                  Do ia = 1, nElem(la)
                     Do  iC = 1, (2*iAng+1)
                        if (lsro) then
                           call mult_sro(FA2(1,1,ia,ic,mvec),nAlpha,
     &                                   C,nC,
     &                                   FB1(1,1,ic,ib,1),nBeta,
     &                                   Fact,Final(1,ia,ib,mVec),
     &                                   Tmp)
                        else
                           Call DGEMM_('N','N',
     &                          nAlpha,nBeta,nC,
     &                          Fact,FA2(1,1,ia,ic,mvec),nAlpha,
     &                          FB1(1,1,ic,ib,1),nC,
     &                          One,Final(1,ia,ib,mVec),nAlpha)
                        endif
                     End do
                  End do
               End do
            End If
         End do
      End do
*                                                                      *
************************************************************************
*                                                                      *
*     Merge integrals with both derivative on center B
      Do iCar = 1, 3
         Do jCar = 1, icar
            mvec=itri(3+icar,3+jcar)
            mvecB=itri(icar,jcar)
            If (jfHess(2,iCar,2,jcar)) Then

               Do  ib = 1, nElem(lb)
                  Do ia = 1, nElem(la)
*
                     Do  iC = 1, (2*iAng+1)
                        if (lsro) then
                           call mult_sro(FA1(1,1,ia,ic,1),nAlpha,
     &                                   C,nC,
     &                                   FB2(1,1,ic,ib,mvecb),nBeta,
     &                                   Fact,Final(1,ia,ib,mVec),
     &                                   Tmp)
                        else
                           Call DGEMM_('N','N',
     &                          nAlpha,nBeta,nC,
     &                          Fact,FA1(1,1,ia,ic,1),nAlpha,
     &                          FB2(1,1,ic,ib,mvecb),nC,
     &                          One,Final(1,ia,ib,mVec),nAlpha)
                        endif
                     End do
                  End do
               End do
            End If
         End do
      End do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
