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

subroutine Isotop_i(n,Element,mAtoms,Temp,nTemp,Coord,double)

use Isotopes, only: ElementList, Initialize_Isotopes, Isotope, MaxAtomNum, PTab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: UtoAU
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: n, mAtoms, nTemp
character(len=2), intent(in) :: Element(mAtoms)
real(kind=wp), intent(out) :: Temp(nTemp)
real(kind=wp), intent(inout) :: Coord(3, mAtoms)
logical(kind=iwp), intent(in) :: double
integer(kind=iwp) :: AtNum, AtNum2, i, iElement, iMass, ineg, ip1, ip2, ipH, ipH2, ipVal, ipVec, IsoNum, j, k, l, m, nAt, Subs, &
                     Subs2
real(kind=wp) :: dMass, dMass1, dMass2, Mass(MxAtom), MassIso, umass(MxAtom)
logical(kind=iwp) :: Changed, EQ, Found, lmass
character(len=3) :: cmass(MxAtom)
real(kind=wp), allocatable :: DefMass(:)
integer(kind=iwp), external :: iNuclearChargeFromSymbol

!                                                                      *
!***********************************************************************
!                                                                      *
ipH = 1

call qpg_iscalar('iMass',lmass)
if (lmass) then
  call Get_iScalar('iMass',iMass)
  call Get_cArray('cmass',cmass,3*iMass)
  call Get_dArray('umass',umass,iMass)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the unsymmetrized Hessian from run file

call Get_dArray('FC-Matrix',Temp(ipH),n**2)

! Put in the initial masses and take backup copy

call Get_Mass_All(Mass,mAtoms)
call mma_allocate(DefMass,mAtoms,label='DefaultMass')
DefMass(:) = Mass(1:mAtoms)

do i=1,mAtoms
  Coord(1,i) = abs(Coord(1,i))
  Coord(2,i) = abs(Coord(2,i))
  Coord(3,i) = abs(Coord(3,i))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over common isotops

ipH2 = ipH+2*n**2
ip1 = ipH2+2*n**2
ip2 = ip1+2*n**2
ipVal = ip2+2*n**2
ipVec = ipVal+2*n**2
call Initialize_Isotopes()
!                                                                      *
!***********************************************************************
!                                                                      *
! Single substitutions

write(u6,*)
write(u6,*) ' Single substitutions:'
write(u6,*) ' -----------------------'
write(u6,*)
do i=1,mAtoms

  if (i > 1) then
    if (EQ(Coord(1,i),Coord(1,i-1))) cycle
  end if

  AtNum = iNuclearChargeFromSymbol(Element(i))
  if (lmass) then

    ! Do according to user list of isotopes.

    Subs = 0
    do nAt=1,iMass
      if (Element(i) == cmass(nAt)) then
        Mass(i) = UtoAU*umass(nAt)
        Subs = 1
      end if
    end do
  else

    ! Get list of natural isotopes

    Subs = max(ElementList(AtNum)%Natural,1)
  end if

  dMass = Mass(i)
  do k=1,Subs
    if (.not. lmass) then
      IsoNum = ElementList(AtNum)%Isotopes(k)%A
      if (IsoNum == nint(dMass/UtoAU)) cycle
      call Isotope(IsoNum,AtNum,Mass(i))
    end if

    call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
    write(u6,*) 'Masses:'
    write(u6,*) '======='
    write(u6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
    write(u6,*)
    write(u6,*)
    write(u6,*) 'Frequencies:'
    write(u6,*) '============'
    call freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
    call GFPrnt_i(Temp(ipVal),n)
  end do
  ! Put back the original mass
  Mass(i) = dMass
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Full substitutions

write(u6,*)
write(u6,*) ' Full substitutions:'
write(u6,*) ' -----------------------'
write(u6,*)
do iElement=1,MaxAtomNum

  Found = .false.
  do i=1,mAtoms
    Found = Found .or. (PTab(iElement) == Element(i))
  end do

  ! Process if element found in the molecule.

  if (Found) then

    ! Loop over all natural isotopes of this element.

    do k=1,max(ElementList(iElement)%Natural,1)
      IsoNum = ElementList(iElement)%Isotopes(k)%A
      call Isotope(IsoNum,iElement,MassIso)

      ! Substitute all instances.

      Changed = .false.
      do i=1,mAtoms
        if (PTab(iElement) == Element(i)) then
          if (IsoNum /= nint(DefMass(i)/UtoAU)) Changed = .true.
          Mass(i) = MassIso
        end if
      end do

      if (Changed) then
        call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
        write(u6,*) 'Masses:'
        write(u6,*) '======='
        write(u6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
        write(u6,*)
        write(u6,*)
        write(u6,*) 'Frequencies:'
        write(u6,*) '============'
        call Freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
        call GFPrnt_i(Temp(ipVal),n)
      end if

    end do

    ! Restore the initial mass array

    Mass(1:mAtoms) = DefMass(:)

  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Double substitutions

if (double) then
  write(u6,*)
  write(u6,*) ' Double substitutions:'
  write(u6,*) ' -----------------------'
  write(u6,*)

  do i=1,mAtoms
    AtNum = iNuclearChargeFromSymbol(Element(i))
    if (lmass) then
      Subs = 0
      do nAt=1,iMass
        if (Element(i) == cmass(nAt)) then
          Mass(i) = UtoAU*umass(nAt)
          Subs = 1
        end if
      end do
    else
      Subs = max(ElementList(AtNum)%Natural,1)
    end if

    ! All isotopes first atom
    dMass1 = Mass(i)
    do k=1,Subs
      if (.not. lmass) then
        IsoNum = ElementList(AtNum)%Isotopes(k)%A
        if (IsoNum == nint(dMass1/UtoAU)) cycle
        call Isotope(IsoNum,AtNum,Mass(i))
      end if

      ! Second atom
      do j=i+1,mAtoms
        AtNum2 = iNuclearChargeFromSymbol(Element(j))
        if (lmass) then
          Subs2 = 0
          do nAt=1,iMass
            if (Element(j) == cmass(nAt)) then
              Mass(j) = UtoAU*umass(nAt)
              Subs2 = 1
            end if
          end do
        else
          Subs2 = max(ElementList(AtNum2)%Natural,1)
        end if

        ! All isotopes for second atom
        dMass2 = Mass(j)
        do l=1,Subs2
          if (.not. lmass) then
            IsoNum = ElementList(AtNum2)%Isotopes(l)%A
            if (IsoNum == nint(dMass2/UtoAU)) cycle
            call Isotope(IsoNum,AtNum2,Mass(j))
          end if

          call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
          write(u6,*) 'Masses:'
          write(u6,*) '======='
          write(u6,'(20I4)') (nint(mass(m)/UtoAU),m=1,mAtoms)
          write(u6,*)
          write(u6,*)
          write(u6,*) 'Frequencies:'
          write(u6,*) '============'
          call freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
          call GFPrnt_i(Temp(ipVal),n)

        end do
        ! Put back the original mass
        Mass(j) = dMass2
      end do   ! End inner loop over atoms

    end do
    ! Put back the original mass
    Mass(i) = dMass1
  end do ! End outer loop over atoms
!                                                                      *
!***********************************************************************
!                                                                      *
end if
call mma_deallocate(DefMass)

return

end subroutine Isotop_i
