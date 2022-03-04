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
      subroutine getAOs2(lhigh)
      implicit real*8(a-h,o-z)
!bs   get expansions of atomic orbitals in contracted functions
#include "para.fh"
#include "param.fh"
#include "nucleus.fh"
      integer closedshells(0:LMAX),openshells(0:LMAX)
      call getocc_ao(int(charge),closedshells,openshells)
      do lrun=0,lhigh
      do irun=1,MxcontL
      do jrun=1,MxcontL
      AOcoeffs(jrun,irun,lrun)=0d0
      enddo
      enddo
      enddo
!BS   write(6,*) 'Orbitals for mean-field'
      do lrun=0,lhigh
!BS   write(6,'(A3,I3)') 'L= ',lrun
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      AOcoeffs(i,i,lrun)=1d0
      enddo
      noccorb(lrun)=closedshells(lrun)
      if (openshells(lrun).gt.0) then
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*DBLE(openshells(lrun))/DBLE(lrun+lrun+1)
      AOcoeffs(i,i,lrun)=1d0
      noccorb(lrun)=i
      endif
      if (noccorb(lrun).gt.0) then
!BS   write(6,'(A,I3)') 'number of orbitals ',noccorb(lrun)
!BS   do iorbital=1,noccorb(lrun)
!BS   write(6,'(A,8F8.4)') 'OCCUPATION: ',(occup(iorbital,lrun),
!BS  *iorbital=1,noccorb(lrun))
!BS   enddo
      endif
      enddo
      return
      end
