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
        subroutine sumabdistt (n,idtot)
!
!        this routine distribute work for n records among
!        nprocab processors with the frequency, corresponding
!        to ideffab values
!
!        n - # of records to be distributed (I)
!        idtot - distribution vector (O)
!               (idtot(i) -  # of records to be realized by i-th node)
!
        implicit none
#include "parallel.fh"
!
        integer n
        integer idtot(1)
!
!        help parameters
!
        integer i,ntot,max,imax
!LD        integer i,j,ntot,max,imax
        real*8 sum
!
!
!1        distribute recordsc according to eff. coeficients
!
        sum=0.0d0
        do i=1,nprocab
        sum=sum+ideffab(i)
        end do
!
        do i=1,nprocab
        idtot(i)=int(((ideffab(i)*n)/sum)+0.5d0)
        end do
!
!2        do corrections, if roundoff errors caused some diferences
!
1        ntot=0
        do i=1,nprocab
        ntot=ntot+idtot(i)
        end do

        if (ntot.gt.n) then
!        ubrat treba (z najvacsieho dielu)
          max=idtot(1)
          imax=1
          do i=1,nprocab
          if (max.lt.idtot(i)) then
          max=idtot(i)
          imax=i
          end if
          end do
          idtot(imax)=idtot(imax)-1
          goto 1
        else if (ntot.lt.n) then
!        pridat treba (k najvacsiemu dielu)
          max=idtot(1)
          imax=1
          do i=1,nprocab
          if (max.lt.idtot(i)) then
          max=idtot(i)
          imax=i
          end if
          end do
          idtot(imax)=idtot(imax)+1
          goto 1
        end if
!
        return
        end
