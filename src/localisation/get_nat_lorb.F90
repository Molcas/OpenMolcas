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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************
      Subroutine Get_Nat_Lorb(Occ,FOcc,nO,nX,jOrb,Umat,iSym)

      Implicit Real*8 (a-h,o-z)
      Real*8  Occ(*), FOcc(*), Umat(*)
      Integer nO, nX, jOrb(nO), iSym
#include "WrkSpc.fh"
!
!
      If (nO.lt.1)  Return
!
      Call GetMem('eta_ik','Allo','Real',ip_eta,2*nX**2+1)
      ip_Z=ip_eta+nX**2
      ip_ZZ=ip_Z+nX
      Call FZero(Work(ip_eta),nX**2)
      Do i=1,nX
         ii=ip_eta+nX*(i-1)+i-1
         Work(ii)=Occ(i)
      End Do
      nXx=Max(1,nX)
      nOx=Max(1,nO)
      Call DGEMM_('N','N',nX,nO,nX,1.0d0,Work(ip_eta),nXx,              &
     &                                  Umat(1),nXx,                    &
     &                            0.0d0,Work(ip_Z),nXx)
      Call DGEMM_('T','N',nO,nO,nX,1.0d0,Umat(1),nXx,                   &
     &                                  Work(ip_Z),nXx,                 &
     &                            0.0d0,Work(ip_eta),nOx)

      Call Eigen_Molcas(nO,Work(ip_eta),Work(ip_Z),Work(ip_ZZ))

      call dcopy_(nO**2,Work(ip_eta),1,Umat(1),1)
      Do i=1,nO
         ii=ip_Z+i-1
         j=jOrb(i)
         FOcc(j)=Work(ii)
      End Do
      Call GetMem('eta_ik','Free','Real',ip_eta,2*nX**2+1)

      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSym)
      End
!***********************************************************************
