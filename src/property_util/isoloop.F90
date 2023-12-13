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

subroutine isoloop(double)

use Isotopes, only: ElementList, Initialize_Isotopes, Isotope, MaxAtomNum, PTab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: UtoAU
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: double
integer(kind=iwp) :: AtNum, AtNum2, i, iElement, iMass, ineg, IsoNum, j, k, l, m, mAtoms, n, nAt, nTemp, Subs, Subs2
real(kind=wp) :: dMass, dMass1, dMass2, MassIso
logical(kind=iwp) :: Changed, EQ, Found, lmass
real(kind=wp), allocatable :: Coord(:,:), DefMass(:), H(:,:), H2(:,:), Mass(:), Temp(:), umass(:), Val(:), Vec(:,:)
character(len=3), allocatable :: cmass(:)
character(len=2), allocatable :: Element(:)
integer(kind=iwp), external :: iNuclearChargeFromSymbol

!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,*)
call CollapseOutput(1,'   Isotopic shifts:')
write(u6,'(3X,A)') '   ----------------'
write(u6,*)

call Get_nAtoms_All(mAtoms)
call mma_allocate(Coord,3,mAtoms,label='Coord')
call Get_Coord_All(Coord,mAtoms)
call mma_allocate(Element,mAtoms,label='Element')
call Get_Name_All(Element)

write(u6,*)
write(u6,*)
write(u6,*) '****************************************'
write(u6,*) '*                                      *'
write(u6,*) '* Isotope shifted frequencies in cm-1  *'
write(u6,*) '*                                      *'
write(u6,*) '****************************************'
write(u6,*)
n = 3*mAtoms
nTemp = 6*2*n**2
call mma_allocate(Temp,nTemp,label='ISOLOOP')
call mma_allocate(H,n,n,label='H')
call mma_allocate(H2,n,n,label='H2')
call mma_allocate(Vec,2*n,n,label='Vec')
call mma_allocate(Val,2*n,label='Val')
!                                                                      *
!***********************************************************************
!                                                                      *
call qpg_iscalar('iMass',lmass)
if (lmass) then
  call Get_iScalar('iMass',iMass)
  call mma_allocate(cmass,iMass,label='cmass')
  call mma_allocate(umass,iMass,label='umass')
  call Get_cArray('cmass',cmass,3*iMass)
  call Get_dArray('umass',umass,iMass)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the unsymmetrized Hessian from run file

call Get_dArray('FC-Matrix',H,n**2)

! Put in the initial masses and take backup copy

call mma_allocate(Mass,mAtoms,label='Mass')
call mma_allocate(DefMass,mAtoms,label='DefaultMass')
call Get_Mass_All(Mass,mAtoms)
DefMass(:) = Mass(:)

do i=1,mAtoms
  Coord(1,i) = abs(Coord(1,i))
  Coord(2,i) = abs(Coord(2,i))
  Coord(3,i) = abs(Coord(3,i))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over common isotopes

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
  Subs = 0
  if (lmass) then

    ! Do according to user list of isotopes.

    do nAt=1,iMass
      if (Element(i) == cmass(nAt)) then
        Mass(i) = UtoAU*umass(nAt)
        Subs = 1
      end if
    end do
  else if (AtNum > 0) then

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

    H2(:,:) = H(:,:)
    write(u6,*) 'Masses:'
    write(u6,*) '======='
    write(u6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
    write(u6,*)
    write(u6,*)
    write(u6,*) 'Frequencies:'
    write(u6,*) '============'
    call freq_i(n,H2,mass,Vec,Val,ineg)
    call GFPrnt_i(Val,n)
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
        H2(:,:) = H(:,:)
        write(u6,*) 'Masses:'
        write(u6,*) '======='
        write(u6,'(20I4)') (nint(mass(l)/UtoAU),l=1,mAtoms)
        write(u6,*)
        write(u6,*)
        write(u6,*) 'Frequencies:'
        write(u6,*) '============'
        call Freq_i(n,H2,mass,Vec,Val,ineg)
        call GFPrnt_i(Val,n)
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

          H2(:,:) = H(:,:)
          write(u6,*) 'Masses:'
          write(u6,*) '======='
          write(u6,'(20I4)') (nint(mass(m)/UtoAU),m=1,mAtoms)
          write(u6,*)
          write(u6,*)
          write(u6,*) 'Frequencies:'
          write(u6,*) '============'
          call freq_i(n,H2,mass,Vec,Val,ineg)
          call GFPrnt_i(Val,n)

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
if (lmass) then
  call mma_deallocate(cmass)
  call mma_deallocate(umass)
end if
call mma_deallocate(Mass)
call mma_deallocate(DefMass)
call mma_deallocate(Temp)
call mma_deallocate(H)
call mma_deallocate(H2)
call mma_deallocate(Vec)
call mma_deallocate(Val)
call mma_deallocate(Element)
call mma_deallocate(Coord)

call CollapseOutput(0,'   Isotopic shifts:')
write(u6,*)

return

end subroutine isoloop
