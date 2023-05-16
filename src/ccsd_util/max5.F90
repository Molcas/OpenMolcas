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

subroutine max5(wrk,wrksize,nind,v,text)
! this routine finds and types:
! a) note
! b) 5 maximal elements with their indices in given vector V
! c) euclidian norm
!
! nind - number of indices in V (I)
! v    - map type of V (I)
! text - notice (I)

use ccsd_global, only: dimm, fullprint, Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, nind
real(kind=wp), intent(in) :: wrk(wrksize)
type(Map_Type), intent(in) :: v
character(len=8), intent(in) :: text
integer(kind=iwp) :: a, b, i, imax(8,5), it, j, nhelp1
real(kind=wp) :: rmax(5), val

!0 set rmax,imax=0

rmax(:) = Zero
imax(:,:) = 0

if (nind == 2) then

  !1 T1aa or T1bb amplitudes

  !1.1 find 5 max

  nhelp1 = v%d(1,1)
  do it=1,v%d(0,5)
    do i=1,dimm(v%d(0,2),v%d(it,4))
      do a=1,dimm(v%d(0,1),v%d(it,3))
        val = wrk(nhelp1)
        ! write this amplitude
        if (abs(val) >= abs(rmax(5))) call max5h1(imax,rmax,v%d(it,3),0,v%d(it,4),0,a,0,i,0,val)
        nhelp1 = nhelp1+1
      end do
    end do
  end do

  !1.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,v,rmax,imax,text)

else if (v%d(0,6) == 0) then

  !2 T2abab amplitudes

  !2.1 find 5 max

  nhelp1 = v%d(1,1)
  do it=1,v%d(0,5)
    do j=1,dimm(v%d(0,4),v%d(it,6))
      do i=1,dimm(v%d(0,3),v%d(it,5))
        do b=1,dimm(v%d(0,2),v%d(it,4))
          do a=1,dimm(v%d(0,1),v%d(it,3))
            val = wrk(nhelp1)
            ! write this amplitude
            if (abs(val) >= abs(rmax(5))) call max5h1(imax,rmax,v%d(it,3),v%d(it,4),v%d(it,5),v%d(it,6),a,b,i,j,val)
            nhelp1 = nhelp1+1
          end do
        end do
      end do
    end do
  end do

  !2.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,v,rmax,imax,text)

else

  !3 T2aaaa or T2bbbb amplitudes

  !3.1 find 5 max

  nhelp1 = v%d(1,1)
  do it=1,v%d(0,5)

    if (v%d(it,3) == v%d(it,4)) then
      ! case syma=symb, symi=symj

      do i=2,dimm(v%d(0,3),v%d(it,5))
        do j=1,i-1
          do a=2,dimm(v%d(0,1),v%d(it,3))
            do b=1,a-1
              val = wrk(nhelp1)
              ! write this amplitude
              if (abs(val) >= abs(rmax(5))) call max5h1(imax,rmax,v%d(it,3),v%d(it,4),v%d(it,5),v%d(it,6),a,b,i,j,val)
              nhelp1 = nhelp1+1
            end do
          end do
        end do
      end do

    else
      ! case syma>symb, symi> symj
      do j=1,dimm(v%d(0,4),v%d(it,6))
        do i=1,dimm(v%d(0,3),v%d(it,5))
          do b=1,dimm(v%d(0,2),v%d(it,4))
            do a=1,dimm(v%d(0,1),v%d(it,3))
              val = wrk(nhelp1)
              ! write this amplitude
              if (abs(val) >= abs(rmax(5))) call max5h1(imax,rmax,v%d(it,3),v%d(it,4),v%d(it,5),v%d(it,6),a,b,i,j,val)
              nhelp1 = nhelp1+1
            end do
          end do
        end do
      end do
    end if

  end do

  !3.2 write output
  if (fullprint >= 0) call max5h2(wrk,wrksize,nind,v,rmax,imax,text)

end if

return

end subroutine max5
