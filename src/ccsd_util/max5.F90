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

subroutine max5(wrk,wrksize,nind,mapd,mapi,text)
! this routine finds and types:
! a) note
! b) 5 maximal elements with their indices in given vector V
! c) euclidian norm
!
! nind - number of indices in V (I)
! mapd - direct map of V (I)
! mapi - inverse map of V (I)
! text - notice (I)

#include "ccsd1.fh"
#include "wrk.fh"
integer nind
integer mapd(0:512,1:6)
integer mapi(1:8,1:8,1:8)
character*8 text
! help variables
integer nhelp1, nhelp2, i, j, a, b, it
real*8 val
integer imax(1:8,1:5)
real*8 rmax(1:5)

!0 set rmax,imax=0

do nhelp1=1,5
  rmax(nhelp1) = 0.0d0
  do nhelp2=1,8
    imax(nhelp2,nhelp1) = 0
  end do
end do

if (nind == 2) then

  !1 T1aa or T1bb amplitudes

  !1.1 find 5 max

  nhelp1 = mapd(1,1)
  do it=1,mapd(0,5)
    do i=1,dimm(mapd(0,2),mapd(it,4))
      do a=1,dimm(mapd(0,1),mapd(it,3))
        val = wrk(nhelp1)
        if (abs(val) >= abs(rmax(5))) then
          ! write this amplitude
          call max5h1(imax,rmax,mapd(it,3),0,mapd(it,4),0,a,0,i,0,val)
        end if
        nhelp1 = nhelp1+1
      end do
    end do
  end do

  !1.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,mapd,mapi,rmax,imax,text)

else if (mapd(0,6) == 0) then

  !2 T2abab amplitudes

  !2.1 find 5 max

  nhelp1 = mapd(1,1)
  do it=1,mapd(0,5)
    do j=1,dimm(mapd(0,4),mapd(it,6))
      do i=1,dimm(mapd(0,3),mapd(it,5))
        do b=1,dimm(mapd(0,2),mapd(it,4))
          do a=1,dimm(mapd(0,1),mapd(it,3))
            val = wrk(nhelp1)
            if (abs(val) >= abs(rmax(5))) then
              ! write this amplitude
              call max5h1(imax,rmax,mapd(it,3),mapd(it,4),mapd(it,5),mapd(it,6),a,b,i,j,val)
            end if
            nhelp1 = nhelp1+1
          end do
        end do
      end do
    end do
  end do

  !2.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,mapd,mapi,rmax,imax,text)

else

  !3 T2aaaa or T2bbbb amplitudes

  !3.1 find 5 max

  nhelp1 = mapd(1,1)
  do it=1,mapd(0,5)

    if (mapd(it,3) == mapd(it,4)) then
      ! case syma=symb, symi=symj

      do i=2,dimm(mapd(0,3),mapd(it,5))
        do j=1,i-1
          do a=2,dimm(mapd(0,1),mapd(it,3))
            do b=1,a-1
              val = wrk(nhelp1)
              if (abs(val) >= abs(rmax(5))) then
                ! write this amplitude
                call max5h1(imax,rmax,mapd(it,3),mapd(it,4),mapd(it,5),mapd(it,6),a,b,i,j,val)
              end if
              nhelp1 = nhelp1+1
            end do
          end do
        end do
      end do

    else
      ! case syma>symb, symi> symj
      do j=1,dimm(mapd(0,4),mapd(it,6))
        do i=1,dimm(mapd(0,3),mapd(it,5))
          do b=1,dimm(mapd(0,2),mapd(it,4))
            do a=1,dimm(mapd(0,1),mapd(it,3))
              val = wrk(nhelp1)
              if (abs(val) >= abs(rmax(5))) then
                ! write this amplitude
                call max5h1(imax,rmax,mapd(it,3),mapd(it,4),mapd(it,5),mapd(it,6),a,b,i,j,val)
              end if
              nhelp1 = nhelp1+1
            end do
          end do
        end do
      end do
    end if

  end do

  !3.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,mapd,mapi,rmax,imax,text)

end if

return

end subroutine max5
