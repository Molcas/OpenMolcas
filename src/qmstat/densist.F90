!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! The RASSI-density matrix subroutine.
      Subroutine DensiSt(Dens,StVec,iS,nSt,iDim)
      Implicit Real*8 (a-h,o-z)
      Dimension Dens(*),StVec(iDim,*)

! iS        -        Which state that is occupied.
! Dens        -        The density
! StVec        -        The coefficients for how the new states are expressed
!                with the old.
      kaunt=0
      Do 101, i=1,nSt
        Do 102, j=1,i
          kaunt=kaunt+1
          Dens(kaunt)=0.0d0
102     Continue
101   Continue
      kaunt=0
      Do 112, ii=1,nSt
        Do 113, jj=1,ii
          kaunt=kaunt+1
          If(ii.eq.jj) then
            Dens(kaunt)=1.0d0*StVec(ii,iS)*StVec(jj,iS)
          Else
            Dens(kaunt)=2.0d0*StVec(ii,iS)*StVec(jj,iS)
          Endif
113     Continue
112   Continue
      Return
      End
