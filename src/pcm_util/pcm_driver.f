************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine PCM_Driver(iPrint,DMat,V,Q,nTs)
      Implicit Real*8 (A-H,O-Z)
      Real*8 V(2,nTs),Q(2,nTs),DMat(nTs,nTs)
      Data Zero,Half /0.0d0,0.5d0/
*
*     Computes PCM solvation charges given the nuclear and electronic electrostatic
*     potential on each tessera.
*     Modifies nuclear repulsion, one-electron and two electron terms.
*
      call dcopy_(2*nTs,Zero,0,Q,1)
      Do iTs = 1, nTs
        Do jTs = 1, nTs
           tmp=DMat(iTs,jTs)+DMat(jTs,iTs)
           DMat(iTs,jTs)=Half*tmp
           DMat(jTs,iTs)=Half*tmp
        End Do
      End Do

      Do iTs = 1, nTs
        Do jTs = 1, nTs
*-------- Actual nuclear charge
          Q(1,iTs) = Q(1,iTs) + DMat(iTs,jTs) * V(1,jTs)
*-------- Effective nuclear charge
c_rl      Q(1,iTs) = Q(1,iTs)
c_rl &             + Half*(DMat(iTs,jTs)+DMat(jTs,iTs)) * V(1,jTs)
*-------- Actual electronic charge
          Q(2,iTs) = Q(2,iTs) + DMat(iTs,jTs) * V(2,jTs)
        EndDo
      EndDo

!      do i=1,nTs   ! yma delete later
!        write(*,*)" == V(2,iTs) diff == ",i,V(2,i)
!      end do

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iPrint)
      End
