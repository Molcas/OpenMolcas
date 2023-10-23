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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
      SubRoutine OCHRR(Target,nPrim,nTrgt,la,lb,ipRs)
!***********************************************************************
! Object: this is a One Center HRR routine.                            *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************
      Implicit None
      Integer, Intent(In) :: nPrim, nTrgt, la, lb
      Integer, Intent(Out):: ipRs
      Real*8, Intent(InOut):: Target(nPrim,nTrgt)

      Integer i, ixyz, ix, iz, nElem, Ind
      Integer iout, ixb, iyb, izb, ixyzb, iybMax, ixa, iyaMax, ixab,
     &        iya, iza, izab, ixyza, iTo, iFrom, iab
!
!     Statment functions
!
      nElem(i) = (i+1)*(i+2)/2
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
      If (la.eq.0 .or. lb.eq.0) Then
         ipRs = 1
         Return
      End If
      iab = 0
      iout = nElem(la+lb)
      ipRs = iout*nPrim + 1
      Do 100 ixb = 0, lb
         iybMax = lb - ixb
         Do 110 iyb = 0, iybMax
            izb = iybMax - iyb
            ixyzb = Ind(lb,ixb,izb)
!
            Do 200 ixa = 0, la
               iyaMax = la - ixa
               ixab = ixa + ixb
               Do 210 iya = 0, iyaMax
                  iza = iyaMax - iya
                  izab = iza + izb
                  ixyza = Ind(la,ixa,iza)
                  iTo = iout + nElem(la)*(ixyzb-1) + ixyza
                  iFrom = iab + Ind(la+lb,ixab,izab)
!
                  Target(:,iTo) = Target(:,iFrom)
!
 210           Continue
 200        Continue
 110     Continue
 100  Continue
!
      Return
      End SubRoutine OCHRR
