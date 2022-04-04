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

subroutine gencoul(l1,l2,l3,l4,makemean,bonn,breit,sameorb,cont4SO,cont4OO,icont4,powexp,coulovlp)
!bs SUBROUTINE to generate all required radial
!bs integrals for the four angular momenta l1-l4

use AMFI_global, only: exponents, Lblocks, Lfirst, Llast, Lmax, Lstarter, Lvalues, MxprimL, nblock, ncontrac, nprimit
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: l1, l2, l3, l4, icont4
logical(kind=iwp), intent(in) :: makemean, bonn, breit, sameorb
real(kind=wp), intent(_OUT_) :: cont4SO(*), cont4OO(*)
real(kind=wp), intent(in) :: powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5)), coulovlp(*)
integer(kind=iwp) :: incl1, incl3, ipow1, ipow2, istart, istart2, Lanf, Lend, Lrun, nanz, nprimprod
real(kind=wp), allocatable :: Prim(:), Quot1(:), Quot2(:), QuotP1(:), QuotP2(:), Scr1(:), Scr2(:)

!bs first of all, this routine determines, for which L
!bs values the radial integrals have to be solved
!bs initialize the number of blocks for the different
!bs l-combinations
!bs no (ss|ss) contributions
if ((l1 == 0) .and. (l2 == 0) .and. (l3 == 0) .and. (l4 == 0)) return
! no integrals for <ss|ss>
if (makemean) then
  nblock = 1  ! sp sp are the first, so the first block
  Lstarter(1) = 1
else
  call SysAbendMsg('gencoul','only mean-field with this version',' ')
end if
!bs keep track of L-values for later purposes
Lvalues(1) = l1
Lvalues(2) = l2
Lvalues(3) = l3
Lvalues(4) = l4
!bs now nanz is given the new value
nanz = ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
nprimprod = nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)

call mma_allocate(Quot1,nPrimProd,Label='Quot1')
call mma_allocate(Quot2,nPrimProd,Label='Quot2')
call mma_allocate(QuotP1,nPrimProd,Label='QuotP1')
call mma_allocate(QuotP2,nPrimProd,Label='QuotP2')
call mma_allocate(Prim,nPrimProd,Label='Prim')
call mma_allocate(Scr1,nPrimProd,Label='Scr1')
call mma_allocate(Scr2,nPrimProd,Label='Scr2')

call initfrac(nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),Quot1,Quot2,exponents(:,l1),exponents(:,l2),exponents(:,l3), &
              exponents(:,l4))
!bs prepare the powers needed for cfunctx
!
! There are seven different CASES of integrals following
!   (   A  --  C)
!
! The structure is the same for all cases, therefore comments can be found only on case A

!bs ####################################################################
!bs   the (+2) cases          CASE A
!bs ####################################################################
incl1 = 1  !  Those increments define the case
incl3 = 1
!bs determine the possible L-values for the integrals by checking for triangular equation

call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)

!bs returns first and last L-values (Lanf,Lend), for which
!bs radial integrals have to be calculated
if (Lend-Lanf >= 0) then
  !bs if there are blocks
  Lblocks(1) = (Lend-Lanf)/2+1 ! L increases in steps of 2, due to parity conservation
  Lfirst(1) = Lanf
  Llast(1) = Lend
else
  Lblocks(1) = 0
end if
if (Lblocks(1) > 0) then    ! integrals have to be calculated
  !bs### check, whether integrals fit on array ################
  if (Lstarter(1)+nanz*Lblocks(1) > icont4) then
    write(u6,*) 'end at: ',Lstarter(1)+nanz*Lblocks(1)
    call SysAbendMsg('gencoul','increase icont4 in amfi.f',' ')
  end if
  !bs### check, whether integrals fit on array ################
  istart = Lstarter(1)
  ! gives the address, where to write the contracted integrals
  !bs ipow1 and ipow2 are the the numbers of powers in the prefactor
  !bs of the function Cfunct
  !bs now loop over possible L-values
  do Lrun=Lfirst(1),Llast(1),2
    ipow1 = 2+(l2+l4+Lrun)/2
    ipow2 = 2+(l1+l3+incl1+incl3+Lrun)/2
    !b those powers have to be generated...
    QuotP1(:) = Quot1(:)**(real(ipow1,kind=wp)-Half)
    !bs those powers have to be generated...
    QuotP2(:) = Quot2(:)**(real(ipow2,kind=wp)-Half)
    ! in buildcoul the radial integrals are calculated
    call buildcoul(l1,l2,l3,l4,incl1,incl3,Lrun,Prim,nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),exponents(:,l1), &
                   exponents(:,l2),exponents(:,l3),exponents(:,l4),powexp(:,:,l3,l1,lrun),powexp(:,:,l4,l2,lrun),QuotP1,QuotP2, &
                   coulovlp)
    !bs in the contcas_ routines the integrals are contracted, including exponents as prefactors...
    if (bonn .or. breit .or. sameorb) then
      call contcasASO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
    else
      call contcasASO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
      call contcasAOO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4OO)
    end if
    istart = istart+nanz! start for next block  contr integr.
  end do
end if
!bs ####################################################################
!bs   the (0) cases         CASE  B
!bs ####################################################################
incl1 = 0
incl3 = 0
call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
if (Lend-Lanf >= 0) then
  Lblocks(2) = (Lend-Lanf)/2+1
  Lfirst(2) = Lanf
  Llast(2) = Lend
  Lblocks(3) = (Lend-Lanf)/2+1
  Lfirst(3) = Lanf
  Llast(3) = Lend
else
  Lblocks(2) = 0
  Lblocks(3) = 0
end if
Lstarter(2) = Lstarter(1)+nanz*Lblocks(1)
Lstarter(3) = Lstarter(2)+nanz*Lblocks(2)
!bs primitive integrals are the same for type 2 and 3  !!!!!
if (Lblocks(2) > 0) then
  !bs### check, whether integrals fit on array ################
  if (Lstarter(2)+2*nanz*Lblocks(2) > icont4) then
    write(u6,*) 'end at: ',Lstarter(2)+2*nanz*Lblocks(2)
    call SysAbendMsg('gencoul','increase icont4 in amfi.f',' ')
  end if
  !bs### check, whether integrals fit on array ################
  istart = Lstarter(2)
  istart2 = Lstarter(3)
  do Lrun=Lfirst(2),Llast(2),2
    ipow1 = 2+(l2+l4+Lrun)/2
    ipow2 = 2+(l1+l3+incl1+incl3+Lrun)/2
    QuotP1(:) = Quot1(:)**(real(ipow1,kind=wp)-Half)
    QuotP2(:) = Quot2(:)**(real(ipow2,kind=wp)-Half)
    call buildcoul(l1,l2,l3,l4,incl1,incl3,Lrun,Prim,nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),exponents(:,l1), &
                   exponents(:,l2),exponents(:,l3),exponents(:,l4),powexp(:,:,l3,l1,lrun),powexp(:,:,l4,l2,lrun),QuotP1,QuotP2, &
                   coulovlp)
    if (bonn .or. breit .or. sameorb) then
      call contcasB1SO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,Prim,Scr1,Scr2,cont4SO)
    else
      call contcasB1SO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,Prim,Scr1,Scr2,cont4SO)
      call contcasB1OO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4OO)
      call contcasB2OO(l1,l2,l3,l4,istart2,Prim,Scr1,Scr2,cont4OO)
    end if
    istart = istart+nanz
    istart2 = istart2+nanz
  end do
end if
!bs ####################################################################
!bs   the (-2) cases      CASE C
!bs ####################################################################
if ((l1 == 0) .or. (l3 == 0)) then
  Lblocks(4) = 0
else
  incl1 = -1
  incl3 = -1
  call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
  if (Lend-Lanf >= 0) then
    Lblocks(4) = (Lend-Lanf)/2+1
    Lfirst(4) = Lanf
    Llast(4) = Lend
  else
    Lblocks(4) = 0
  end if
end if
Lstarter(4) = Lstarter(3)+nanz*Lblocks(3)
if (Lblocks(4) > 0) then
  !bs### check, whether integrals fit on array ################
  if (Lstarter(4)+nanz*Lblocks(4) > icont4) then
    write(u6,*) 'end at: ',Lstarter(4)+nanz*Lblocks(4)
    call SysAbendMsg('gencoul','increase icont4 in amfi.f',' ')
  end if
  !bs### check, whether integrals fit on array ################
  istart = Lstarter(4)
  do Lrun=Lfirst(4),Llast(4),2
    ipow1 = 2+(l2+l4+Lrun)/2
    ipow2 = 2+(l1+l3+incl1+incl3+Lrun)/2
    QuotP1(:) = Quot1(:)**(real(ipow1,kind=wp)-Half)
    QuotP2(:) = Quot2(:)**(real(ipow2,kind=wp)-Half)
    call buildcoul(l1,l2,l3,l4,incl1,incl3,Lrun,Prim,nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),exponents(:,l1), &
                   exponents(:,l2),exponents(:,l3),exponents(:,l4),powexp(:,:,l3,l1,lrun),powexp(:,:,l4,l2,lrun),QuotP1,QuotP2, &
                   coulovlp)
    if (bonn .or. breit .or. sameorb) then
      call contcasCSO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
    else
      call contcasCSO(l1,l2,l3,l4,istart,Prim,Scr1,Scr2,cont4SO)
      call contcasCOO(l1,l2,l3,l4,istart,prim,Scr1,Scr2,cont4OO)
    end if
    istart = istart+nanz
  end do
end if

call mma_deallocate(Quot2)
call mma_deallocate(Quot1)
call mma_deallocate(QuotP2)
call mma_deallocate(QuotP1)
call mma_deallocate(Prim)
call mma_deallocate(Scr2)
call mma_deallocate(Scr1)

return

end subroutine gencoul
