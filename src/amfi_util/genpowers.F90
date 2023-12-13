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

subroutine genpowers(Lhigh,powexp,coulovlp)

use AMFI_global, only: df, exponents, Lmax, MxprimL, nprimit
use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lhigh
real(kind=wp), intent(out) :: powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5)), &
                              coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)
integer(kind=iwp) :: incl1, incl2, iprim1, iprim2, irun1, irun2, L1, L2, Lrun, n1, n2
real(kind=wp) :: df1, df12, df2, fact, fact1, fact2, factor

!bs set some often used powers of exponents
do L2=0,Lhigh
  do L1=0,L2
    powexp(1:nprimit(L1),1:nprimit(L2),L1,L2,0) = One
  end do
end do
do L2=0,Lhigh
  do L1=0,L2
    do Lrun=1,(L1+L2+5)
      do irun2=1,nprimit(L2)
        do irun1=1,nprimit(L1)
          fact = sqrt(Half*(exponents(irun1,L1)+exponents(irun2,L2)))
          !bs write(u6,*) 'fact',fact,'powexp',powexp(irun1,irun2,L1,L2,Lrun-1)
          powexp(irun1,irun2,L1,L2,Lrun) = powexp(irun1,irun2,L1,L2,Lrun-1)*fact
        end do
      end do
    end do
  end do
end do
!bs generate coulovlp = overlap for normalized functions, but sometimes
!bs with shifted l-values
do l2=0,lhigh
  do incl2=-1,1
    if (l2+incl2 >= 0) then  ! do not lower l for s-functions
      n2 = l2+incl2+1
      df2 = One/sqrt(df(n2+n2-1))
      do l1=0,l2
        do incl1=-1,1
          if (l1+incl1 >= 0) then ! do not lower l for s-functions
            n1 = l1+incl1+1
            df1 = One/sqrt(df(n1+n1-1))
            df12 = df(n1+n2-1)
            do iprim2=1,nprimit(l2)
              fact2 = sqrt(powexp(iprim2,iprim2,l2,l2,n2+n2+1))
              factor = fact2*df1*df2*df12
              do iprim1=1,nprimit(l1)
                fact1 = sqrt(powexp(iprim1,iprim1,l1,l1,n1+n1+1))
                coulovlp(iprim1,iprim2,incl1,incl2,l1,l2) = fact1*factor/powexp(iprim1,iprim2,l1,l2,n1+n2+1)
                !BS write(u6,*) 'fact1',fact1,'factor ',factor,'powexp ',powexp(iprim1,iprim2,l1,l2,n1+n2+1)
              end do
            end do
          end if
        end do
      end do
    end if
  end do
end do

return

end subroutine genpowers
