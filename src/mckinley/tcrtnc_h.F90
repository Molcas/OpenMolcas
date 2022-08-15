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
! Copyright (C) 1990,1992,1994, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************
      SubRoutine Tcrtnc_h(Coef1,n1,m1,Coef2,n2,m2,                      &
     &                    Coef3,n3,m3,Coef4,n4,m4,                      &
     &                    ACInt,mabcd,Scrtch,nScr,ACOut,                &
     &                    IndZet,lZeta,IndEta,lEta)
!***********************************************************************
!                                                                      *
! Object: to transform the integrals from primitives to contracted     *
!         basis functions. The subroutine will do both complete and    *
!         incomplete transformations.                                  *
!                                                                      *
!         Observe that ACInt and ACOut may overlap!!!!                 *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified to back transformation, January '92.            *
!***********************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
!
!---- Cache size
!
#include "lCache.fh"
      Real*8 Coef1(n1,m1), Coef2(n2,m2), Coef3(n3,m3), Coef4(n4,m4),    &
     &       ACInt(m1*m2*m3*m4,mabcd), Scrtch(nScr),                    &
     &       ACOut(n1*n2*n3*n4,mabcd)
      Integer IndZet(lZeta), IndEta(lEta)
!
      iRout = 18
      iPrint = nPrint(iRout)
!     iPrint=99
!
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(AB|CD)',ACInt,           &
     &                                m1*m2*m3*m4*mabcd)
      If (iPrint.ge.99) Then
         Call RecPrt(' In Tcrtnc: P(ab|cd)',' ',ACInt,m1*m2,m3*m4*mabcd)
         Call RecPrt(' Coef1',' ',Coef1,n1,m1)
         Call RecPrt(' Coef2',' ',Coef2,n2,m2)
         Call RecPrt(' Coef3',' ',Coef3,n3,m3)
         Call RecPrt(' Coef4',' ',Coef4,n4,m4)
         Write (6,*) n1, n2, n3, n4
      End If
!
!---- Reduce for contraction matrix
      nCache = (3*lCache)/4 - n1*m1 - n2*m2
      lsize= m1*m2 + m2*n1
      nVec = m3*m4*mabcd
      IncVec = Min(Max(1,nCache/lsize),nVec)
      ipA3=1
      nA3=nVec*lZeta         ! This is the same for the second set!
      ipA2=ipA3 + nA3
!
      Call TncHlf_h(Coef1,m1,n1,Coef2,m2,n2,iDum,lZeta,nVec,            &
     &              IncVec,ACInt,Scrtch(ipA2),Scrtch(ipA3),IndZet)
!
      nCache = (3*lCache)/4 - n3*m3 - n4*m4
      lsize = m3*m4 + m4*n3
      nVec = mabcd*lZeta
      IncVec = Min(Max(1,nCache/lsize),nVec)
!
      lZE=lZeta*lEta
      If (mabcd.ne.1) Then
         Call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,iDum,lEta,nVec,          &
     &                 IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
         Call DGeTMO(ACOut,mabcd,mabcd,lZE,Scrtch,lZE)
         call dcopy_(mabcd*lZE,Scrtch,1,ACOut,1)
      Else
         Call TncHlf_h(Coef3,m3,n3,Coef4,m4,n4,iDum,lEta,nVec,          &
     &                 IncVec,Scrtch(ipA3),Scrtch(ipA2),ACOut,IndEta)
      End If
!
      If (iPrint.ge.59)                                                 &
     &  Call RecPrt(' In Tcrtnc: P(ab|cd) ',' ',ACOut,mabcd,lZE)
      If (iPrint.ge.19) Call WrCheck('Tcrtnc:P(ab|cd)',ACOut,           &
     &                                lZE*mabcd)
!
!     Call GetMem('Tcrtnc','CHECK','REAL',iDum,iDum)
      Return
      End
