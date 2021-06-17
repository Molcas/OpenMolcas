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

use Isotopes

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "constants2.fh"
#include "WrkSpc.fh"
integer AtNum, IsoNum, AtNum2, Subs, Subs2
real*8 Temp(nTemp), Coord(3,mAtoms)
real*8 MassIso
real*8 Mass(MxAtom)
character*2 Element(mAtoms)
logical EQ, double, lmass, Found, Changed
character*3 cmass(MxAtom)
real*8 umass(MxAtom)

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
call GetMem('DefaultMass','ALLO','REAL',ipDefMass,mAtoms)
call dCopy_(mAtoms,Mass,1,Work(ipDefMass),1)

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

write(6,*)
write(6,*) ' Single substitutions:'
write(6,*) ' -----------------------'
write(6,*)
do i=1,mAtoms

  if (i > 1) then
    if (EQ(Coord(1,i),Coord(1,i-1))) Go To 94
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
    write(6,*) 'Masses:'
    write(6,*) '======='
    write(6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
    write(6,*)
    write(6,*)
    write(6,*) 'Frequencies:'
    write(6,*) '============'
    call freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
    call GFPrnt_i(Temp(ipVal),n)
  end do
  ! Put back the original mass
  Mass(i) = dMass
94 continue
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Full substitutions

write(6,*)
write(6,*) ' Full substitutions:'
write(6,*) ' -----------------------'
write(6,*)
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
          if (IsoNum /= nint(Work(ipDefMass+i-1)/UtoAU)) Changed = .true.
          Mass(i) = MassIso
        end if
      end do

      if (Changed) then
        call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
        write(6,*) 'Masses:'
        write(6,*) '======='
        write(6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
        write(6,*)
        write(6,*)
        write(6,*) 'Frequencies:'
        write(6,*) '============'
        call Freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
        call GFPrnt_i(Temp(ipVal),n)
      end if

    end do

    ! Restore the initial mass array

    call dCopy_(mAtoms,Work(ipDefMass),1,Mass,1)

  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Double substitutions

if (double) then
  write(6,*)
  write(6,*) ' Double substitutions:'
  write(6,*) ' -----------------------'
  write(6,*)

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
          write(6,*) 'Masses:'
          write(6,*) '======='
          write(6,'(20I4)') (nint(mass(m)/UtoAU),m=1,mAtoms)
          write(6,*)
          write(6,*)
          write(6,*) 'Frequencies:'
          write(6,*) '============'
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
call GetMem('DefaultMass','FREE','REAL',ipDefMass,mAtoms)

return

end subroutine Isotop_i
