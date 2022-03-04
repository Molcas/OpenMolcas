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

subroutine buildcoul(l1,l2,l3,l4,incl1,incl3,Lrun,prmints,nprim1,nprim2,nprim3,nprim4,expo1,expo2,expo3,expo4,power13,power24, &
                     quotpow1,quotpow2,coulovlp)
!bs ##################################################################
!bs  purpose: builds up the coulomb integrals
!bs  inbetween primitives and multiplies
!bs  with extra factors to correct the
!bs  normalization
!bs ##################################################################
! l1,l2,l3,l4                 : angular momenta of primitives
! incl1,incl3                 : shifts for different radial integrals
! Lrun                        : L-value for coulomb integrals
! prmints                     ! scratch for prim integral
! nprim1,nprim2,nprim3,nprim4 : number of primitives
! expo1,expo2,expo3,expo4     : arrays with the exponents

implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "dofuc.fh"
#include "real.fh"
dimension expo1(nprim1), expo2(nprim2), expo3(nprim3), expo4(nprim4), prmints(nprim1,nprim2,nprim3,nprim4), &
          power13(MxprimL,MxprimL), power24(MxprimL,MxprimL), quotpow1(nprim1,nprim2,nprim3,nprim4), &
          quotpow2(nprim1,nprim2,nprim3,nprim4), fraclist1(0:Lmax+3), fraclist2(0:Lmax+3), fact(MxprimL), frac(MxprimL), &
          cfunctx1(MxprimL), cfunctx2(MxprimL), coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)

root8ovpi = sqrt(8d0/pi)
!bs ##################################################################
!bs   prepare indices for coulint
!bs ##################################################################
n1 = l1+incl1+1
n2 = l2+1
n3 = l3+incl3+1
n4 = l4+1
n13 = n1+n3
n24 = n2+n4
index1 = N13-Lrun-1
index2 = n24+Lrun
index3 = N24-Lrun-1
index4 = n13+Lrun
do krun=0,(index1-1)/2
  fraclist1(krun) = dffrac(krun+krun+index2-1,krun+krun)*dffrac(1,index2-1)
end do
do krun=0,(index3-1)/2
  fraclist2(krun) = dffrac(krun+krun+index4-1,krun+krun)*dffrac(1,index4-1)
end do
!bs ##################################################################
!bs   common factors including double factorials
!bs ##################################################################
doff1 = dffrac(index1-1,n13-1)*dffrac(n24+Lrun-1,n24-1)
doff2 = dffrac(index3-1,n24-1)*dffrac(n13+Lrun-1,n13-1)
if (index1 == 1) then
  do irun4=1,nprim4
    do irun3=1,nprim3
      if (l2 == l4) then
        limit2 = irun4
      else
        limit2 = nprim2
      end if
      do irun2=1,limit2
        pow24inv = doff1/power24(irun4,irun2)
        if (l1 == l3) then
          limit1 = irun3
        else
          limit1 = nprim1
        end if
        do irun1=1,limit1
          prmints(irun1,irun2,irun3,irun4) = quotpow1(irun1,irun2,irun3,irun4)*sqrt(0.5d0*(expo1(irun1)+expo3(irun3)))* &
                                             power13(irun3,irun1)*pow24inv
        end do
      end do
    end do
  end do
else
  do irun4=1,nprim4
    do irun3=1,nprim3
      if (l2 == l4) then
        limit2 = irun4
      else
        limit2 = nprim2
      end if
      do irun2=1,limit2
        alpha24inv = 1d0/(expo2(irun2)+expo4(irun4))
        pow24inv = doff1/power24(irun4,irun2)
        if (l1 == l3) then
          limit1 = irun3
        else
          limit1 = nprim1
        end if
        do irun1=1,limit1
          a1324 = alpha24inv*(expo1(irun1)+expo3(irun3))
          Cfunctx1(irun1) = fraclist1(0)
          frac(irun1) = a1324/(1d0+a1324)
          fact(irun1) = frac(irun1)
        end do
        !vocl loop,repeat(Lmax+3)
        do k=1,(index1-1)/2
          do irun1=1,limit1
            Cfunctx1(irun1) = Cfunctx1(irun1)+fraclist1(k)*fact(irun1)
          end do
          do irun1=1,limit1
            fact(irun1) = fact(irun1)*frac(irun1)
          end do
        end do
        do irun1=1,limit1
          alpha13 = 0.5d0*(expo1(irun1)+expo3(irun3))
          prmints(irun1,irun2,irun3,irun4) = quotpow1(irun1,irun2,irun3,irun4)*sqrt(alpha13)*power13(irun3,irun1)*pow24inv* &
                                             Cfunctx1(irun1)
        end do
      end do
    end do
  end do
end if
if (index3 == 1) then
  do irun4=1,nprim4
    do irun3=1,nprim3
      if (l2 == l4) then
        limit2 = irun4
      else
        limit2 = nprim2
      end if
      do irun2=1,limit2
        pow24 = doff2*power24(irun4,irun2)*sqrt(0.5d0*(expo2(irun2)+expo4(irun4)))
        if (l1 == l3) then
          limit1 = irun3
        else
          limit1 = nprim1
        end if
        do irun1=1,limit1
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)+pow24*quotpow2(irun1,irun2,irun3,irun4)/ &
                                             power13(irun3,irun1)
        end do
      end do
    end do
  end do
else
  do irun4=1,nprim4
    do irun3=1,nprim3
      if (l2 == l4) then
        limit2 = irun4
      else
        limit2 = nprim2
      end if
      do irun2=1,limit2
        alpha24 = expo2(irun2)+expo4(irun4)
        pow24 = doff2*power24(irun4,irun2)*sqrt(0.5d0*alpha24)
        if (l1 == l3) then
          limit1 = irun3
        else
          limit1 = nprim1
        end if
        do irun1=1,limit1
          a2413 = alpha24/(expo1(irun1)+expo3(irun3))
          Cfunctx2(irun1) = fraclist2(0)
          frac(irun1) = a2413/(1d0+a2413)
          fact(irun1) = frac(irun1)
        end do
        !vocl loop,repeat(Lmax+3)
        do k=1,(index3-1)/2
          do irun1=1,limit1
            Cfunctx2(irun1) = Cfunctx2(irun1)+fraclist2(k)*fact(irun1)
          end do
          do irun1=1,limit1
            fact(irun1) = fact(irun1)*frac(irun1)
          end do
        end do
        do irun1=1,limit1
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)+quotpow2(irun1,irun2,irun3,irun4)*Cfunctx2(irun1)* &
                                             pow24/power13(irun3,irun1)
        end do
      end do
    end do
  end do
end if
!bs make some mirroring for identical l-values
!bs for the case that l1=l3
if (l1 == l3) then
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=1,nprim2
        do irun1=irun3+1,nprim1
          prmints(irun1,irun2,irun3,irun4) = prmints(irun3,irun2,irun1,irun4)
        end do
      end do
    end do
  end do
end if
!bs for the case that l2=l4
if (l2 == l4) then
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=irun4+1,nprim2
        do irun1=1,nprim1
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun4,irun3,irun2)
        end do
      end do
    end do
  end do
end if
!bs some factors which are the same for all cases
do irun4=1,nprim4
  do irun3=1,nprim3
    do irun2=1,nprim2
      do irun1=1,nprim1
        prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)*coulovlp(irun4,irun2,0,0,l4,l2)* &
                                           coulovlp(irun3,irun1,incl3,incl1,l3,l1)*root8ovpi
      end do
    end do
  end do
end do

!bs look for additional factors, as the
!bs coulomb integrals are calculated
!bs for normalized functions with that
!bs specific l

!bs if l was increased by one, the factor is
!bs 0.5*sqrt((2l+3)/(exponent))
!bs if l was decreased by one, the factor is
!bs 2d0*sqrt(exponent/(2l+1))

!bs  check for first function

if (incl1 == 1) then
  fact1 = 0.5d0*sqrt(dble(l1+l1+3))
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=1,nprim2
        do irun1=1,nprim1
          factor = fact1/sqrt(expo1(irun1))
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)*factor
        end do
      end do
    end do
  end do
else if (incl1 == -1) then
  fact1 = 2d0/sqrt(dble(l1+l1+1))
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=1,nprim2
        do irun1=1,nprim1
          factor = fact1*sqrt(expo1(irun1))
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)*factor
        end do
      end do
    end do
  end do
end if

!bs check for third function

if (incl3 == 1) then
  fact1 = 0.5d0*sqrt(dble(l3+l3+3))
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=1,nprim2
        do irun1=1,nprim1
          factor = fact1/sqrt(expo3(irun3))
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)*factor
        end do
      end do
    end do
  end do
else if (incl3 == -1) then
  fact1 = 2d0/sqrt(dble(l3+l3+1))
  do irun4=1,nprim4
    do irun3=1,nprim3
      do irun2=1,nprim2
        do irun1=1,nprim1
          factor = fact1*sqrt(expo3(irun3))
          prmints(irun1,irun2,irun3,irun4) = prmints(irun1,irun2,irun3,irun4)*factor
        end do
      end do
    end do
  end do
end if

return

end subroutine buildcoul
