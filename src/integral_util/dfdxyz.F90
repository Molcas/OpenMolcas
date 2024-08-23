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
! Copyright (C) 1995, Roland Lindh                                     *
!               2004, Takashi Tsuchiya                                 *
!***********************************************************************
      Subroutine dFdxyz(mterm,mform,N,jp,ip,ixyz,ipf,jdrv)
!
      Implicit None
      Integer mTerm, mForm, jp, ip, ixyz, ipf, jdrv
      Integer N(mterm,5,mform)

      Integer nTerm, iTerm, jTerm, i
!
!
!     ipf: Phase factor in integer
!
!
!
!
      nterm=2**jdrv
      iterm=0
!
      Do jterm=1,nterm
!
!
!        downward operation
!        (derivation of angular part)
!
!
         iterm=iterm+1
         Do i=1,5
            If (i .EQ. ixyz) then
               N(iterm,ixyz,ip)=N(jterm,ixyz,jp)-1
            Else
               N(iterm,i,ip)=N(jterm,i,jp)
            End if
         End do
         N(iterm,5,ip)=N(iterm,5,ip)*N(jterm,ixyz,jp)
         N(iterm,5,ip)=N(iterm,5,ip)*ipf
!
!
!        upward operation
!        (derivation of radial  part)
!
!
         iterm=iterm+1
         Do i=1,5
            If (i .EQ. ixyz) then
               N(iterm,ixyz,ip)=N(jterm,ixyz,jp)+1
            Else
               N(iterm,i,ip)=N(jterm,i,jp)
            End if
         End do
         N(iterm,4,ip)=N(iterm,4,ip)+1
         N(iterm,5,ip)=N(iterm,5,ip)*ipf
      End do
      End Subroutine dFdxyz
