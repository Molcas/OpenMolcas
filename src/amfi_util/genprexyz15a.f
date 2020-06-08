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
      subroutine genprexyz15a(icheckxy,icheckz,interxyz)
      implicit real*8(a-h,o-z)
      dimension icheckxy(*),icheckz(*),interxyz(16,*)
#include "para.fh"
#include "Molcas.fh"
cbs   the following M values are the ones from the cartesian
cbs   linear combinations. interxyz gives the sign sequence
cbs   for interacting spherical functions, starting with
cbs   type 1 (++++) and ending with type 16 (-++-)
      ilauf=1
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      irun=0
      if (icheckxy(ilauf)+icheckz(ilauf).gt.0) then
          if (iabs(m1+m2-m3-m4).le.1) then
          irun=irun+1
          interxyz(irun,ilauf)=1          ! + + + +
                  if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=2  ! - - - -
                  endif
          endif
          if (iabs(m1+m2-m3+m4).le.1) then
                  if (m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=3  ! + + + -
                  endif
                  if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=4  ! - - - +
                  endif
          endif
          if (iabs(m1+m2+m3-m4).le.1) then
                  if (m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=5  ! + + - +
                  endif
                  if (m1.gt.0.and.m2.gt.0.and.
     *            m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=6  ! - - + -
                  endif
          endif
          if (iabs(m1-m2-m3-m4).le.1) then
                  if (m2.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=7  ! + - + +
                  endif
                  if (m1.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=8  ! - + - -
                  endif
          endif
          if (iabs(-m1+m2-m3-m4).le.1) then
                  if (m1.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=9  ! - + + +
                  endif
                  if (m2.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=10 ! + - - -
                  endif
          endif
          if (iabs(m1+m2+m3+m4).le.1) then
                  if (m3.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=11 ! + + - -
                  endif
                  if (m1.gt.0.and.m2.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=12 ! - - + +
                  endif
          endif
          if (iabs(m1-m2-m3+m4).le.1) then
                  if (m2.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=13 ! + - + -
                  endif
                  if (m1.gt.0.and.m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=14 ! - + - +
                  endif
          endif
          if (iabs(m1-m2+m3-m4).le.1) then
                  if (m2.gt.0.and.m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=15 ! + - - +
                  endif
                  if (m1.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,ilauf)=16 ! - + + -
                  endif
          endif
      endif
      ilauf=ilauf+1
      enddo
      enddo
      enddo
      enddo
      return
      end
