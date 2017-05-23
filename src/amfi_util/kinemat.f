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
      Subroutine kinemat(L,ndim,evtkin,type1,type2,Energy)
      implicit real*8 (a-h,o-z)
      parameter (fine=7.29735308D-03) !TO_BE_CHECKED
cbs   at least it's identical with Odd's valuE
      parameter (speed=1d0/fine)
      parameter (speed2=speed*speed)
      parameter (speed4=speed2*speed2)
cbs   this routine generates the kinematic A-factors=sqrt((E+mc^2)/(2E))
cbs   (type1) and   c*A/(E+mc^2) (type2)
cbs   The c in the second kinematic factor comes from Jan Almloef and
cbs   Odd Gropen in Rev in Comp.Chem. 8(1996)
      dimension evtkin(*),type1(*),type2(*),Energy(*)
c     E= sqrt(p**2 c**2 + m**2 c**4)
c     p**2= 2*m*TKIN
c     with m = 1
      do Irun=1,ndim
      if (evtkin(Irun).lt.0.0D0) call SysAbendMsg('kinemat',
     & 'strange kinetic energy ',' ')
      Energy(Irun)=(evtkin(Irun)+evtkin(Irun))*speed2+speed4
      enddo
      do Irun=1,ndim
      Energy(Irun)=sqrt(energy(irun))
      enddo
      do Irun=1,ndim
!     sqrt((E+mc^2)/(2E)):
      type1(Irun)=sqrt(0.5d0*(1d0+speed2/Energy(Irun)))
      enddo
!      c*A/(E+mc^2)
      do Irun=1,ndim
      type2(Irun)=speed*type1(Irun)/(Energy(Irun)+speed2)
      enddo
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(L)
      end
