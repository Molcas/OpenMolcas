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

subroutine Prop(Short,qplab,cen1,cen2,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,lpole,ifallorb)
!***********************************************************************
!                                                                      *
!     purpose: preprocessing of tables for different tensor            *
!              properties                                              *
!                                                                      *
!     Short           logical option for either Short (Short           *
!                     =.true., the total electronic contribution)      *
!                     long (Short=.false., orbital contributions)      *
!                     output                                           *
!     oplab           operator label as defined in SEWARD              *
!     cen1(1:3)       coordinates of centre no.1 (operator)            *
!     cen2(1:3)       coordinates of centre no.2 (gauge)               *
!     nIrrep          the number of irreducible representations        *
!     nBas            the number of functions in each representa-      *
!     (0:nIrrep-1)    tion                                             *
!     nTot            the total number of elements supplied for        *
!                     each component of the property tensor; equal     *
!                     either to 1 (total electronic and total nuc-     *
!                     lear contributions) or to the dimension of       *
!                     the basis set.                                   *
!     Occ(1:nTot)     occupation numbers for all eigenvectors,         *
!                     a dummy for Short outputs                        *
!     ThrSV           threshold for occupation numbers; If             *
!                     Occ(i).le.ThrSV the contribution will not        *
!                     be printed                                       *
!     PrEl(1:nTot,    matrix elements for all components 1,2,...,      *
!          1:maxlab)  maxlab, nTot entries for each component          *
!                     maxlab=(lpole+1)*(lpole+2)/2                     *
!     PrNu(1:maxlab)  nuclear contributions for each component         *
!     ifallorb        logical option for whether the property of       *
!                     all orbitals are printed (and not weighted by    *
!                     occupation number)in property calculation when   *
!                     short=.false. (S.S.Dong, 2018)                   *
!                                                                      *
!     lpole is the value for l in l-pole moments                       *
!     allocate integer storage area of the size appropriate            *
!     for the actual lpole value: size as below                        *
!                                                                      *
! 2000 Dept. of Chem. Phys., Univ. of Lund, Sweden                     *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Angstrom, Debye
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Short, ifallorb
character(len=8), intent(in) :: qplab
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nTot, lpole
real(kind=wp), intent(in) :: cen1(3), cen2(3), Occ(nTot), ThrSV
real(kind=wp), intent(inout) :: PrEl(nTot,(lpole+1)*(lpole+2)/2), PrNu((lpole+1)*(lpole+2)/2)
integer(kind=iwp) :: i, icen, icen1, ilab, inp, ip_, iPL, iSt, iTol, iTol_E0, iTol_E1, ix, ixx, iy, iyy, iz, izz, j, jMax, maxlab
real(kind=wp) :: Fact, Molecular_Charge = Zero, sig, tmp, X_Coor, Y_Coor, Z_Coor
logical(kind=iwp) :: StoreInfo
integer(kind=iwp), parameter :: lmax = 16
character(len=lmax) :: lab
character(len=80) :: Line
character(len=8) :: oplab
character(len=5) :: lab5
character(len=4) :: lab4
character(len=3) :: lab3
real(kind=wp), allocatable :: PrElAug(:,:), PrNuAug(:), PrTot(:), tmat(:,:), temp(:)
character(len=lmax), allocatable :: labs(:), labsAug(:)
integer(kind=iwp), external :: Cho_X_GetTol, iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
#include "hfc_logical.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Some properties can be calculated at many, many centers,
! and writing all values through Add_Info affects performance,
! those will be deactivated and the summation is printed by the
! calling function.
! This only affects the call to Add_Info, otherwise the properties
! are calculated and printed normally
!
StoreInfo = .true.
!                                                                      *
!***********************************************************************
!                                                                      *

maxlab = (lpole+1)*(lpole+2)/2
call mma_allocate(labs,maxlab,label='labs')
call mma_allocate(tmat,maxlab,maxlab,label='tmat')
call mma_allocate(temp,maxlab,label='temp')
call mma_allocate(PrElAug,nTot,maxlab+1,label='PrElAug')
call mma_allocate(PrNuAug,maxlab+1,label='PrNuAug')
call mma_allocate(labsAug,maxlab+1,label='labsAug')

! decipher the operator label:oplab
oplab = qplab
call UpCase(oplab)
lab4 = oplab(1:4)
if (lab4 == 'MLTP') then
  if (lpole < 0) then
    call End1()
    return
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Multipole moment section ... generate labels for printing

  if (lPole > lMax) then
    write(u6,*) 'Prop: lPole.gt.lMax'
    write(u6,*) 'lPole=',lPole
    write(u6,*) 'Increase lMax and recompile!'
    call Abend()
  end if
  ilab = 0
  do ix=lpole,0,-1
    do iy=lpole-ix,0,-1
      iz = lpole-ix-iy
      ilab = ilab+1
      lab = '                '
      do izz=lmax,lmax-iz+1,-1
        lab(izz:izz) = 'Z'
      end do
      do iyy=lmax-iz,lmax-iz-iy+1,-1
        lab(iyy:iyy) = 'Y'
      end do
      do ixx=lmax-iz-iy,lmax-lpole+1,-1
        lab(ixx:ixx) = 'X'
      end do
      labs(ilab) = lab
    end do
  end do

  maxlab = ilab
  call mma_allocate(PrTot,maxlab,label='PrTot')

  ! Print cartesian moments

  !--------------------------------------------------------------------*
  if ((iPL == 2) .and. Short) then
  !--------------------------------------------------------------------*

    PrTot(:) = PrNu(1:MaxLab)-PrEl(1,1:MaxLab)

    ! New style output

    Line = ' '
    if (lPole == 0) then
      Line = 'Charge (e):'
      Fact = One
    else if (lPole == 1) then
      Line = 'Dipole Moment (Debye):'
      Fact = Debye
    else if (lPole == 2) then
      Line = 'Quadrupole Moment (Debye*Ang):'
      Fact = Debye*Angstrom
    else if (lPole == 3) then
      Line = 'Octupole Moment (Debye*Ang**2):'
      Fact = Debye*Angstrom**2
    else if (lPole == 4) then
      Line = 'Hexadecapole Moment (Debye*Ang**3):'
      Fact = Debye*Angstrom**3
    else
      Line = ''
      if (lpole <= 9) then
        write(Line,'(I1)') lPole
        iSt = 2
      else
        write(Line,'(I2)') lPole
        iSt = 3
      end if
      Line(iSt:iSt+26) = 'th-pole Moment (Debye*Ang**'
      iSt = iSt+27
      if (lpole <= 10) then
        write(Line(iSt:iSt),'(I1)') lpole-1
        iSt = iSt+1
      else
        write(Line(iSt:iSt+1),'(I2)') lpole-1
        iSt = iSt+2
      end if
      Line(iSt:iSt+2) = '):'
      Fact = Debye*Angstrom**(lPole-1)
    end if
    write(u6,'(6X,A)') trim(Line)
    if (lpole > 0) then
      write(u6,'(6X,A,3F10.4)') 'Origin of the operator (Ang)=',(cen1(i)*Angstrom,i=1,3)
    end if
    if (lPole == 0) then
      write(u6,'(6X,A,A,F10.4)') labs(1),'=',PrTot(1)*Fact
      Molecular_Charge = PrTot(1)*Fact
    else if (lPole == 1) then
      tmp = sqrt(PrTot(1)**2+PrTot(2)**2+PrTot(3)**2)
      write(u6,'(4X,4(A,A,ES12.4))') labs(1),'=',PrTot(1)*Fact,labs(2),'=',PrTot(2)*Fact,labs(3),'=',PrTot(3)*Fact, &
                                     '           Total','=',tmp*Fact
      if (abs(Molecular_Charge) > 0.9_wp) then
        write(u6,'(6X,A)') 'Center of Charge (Ang)'
        X_Coor = Angstrom*(PrTot(1)/Molecular_Charge)
        Y_Coor = Angstrom*(PrTot(2)/Molecular_Charge)
        Z_Coor = Angstrom*(PrTot(3)/Molecular_Charge)
        write(u6,'(6X,3(A,A,F14.8))') labs(1),'=',X_Coor,labs(2),'=',Y_Coor,labs(3),'=',Z_Coor
        Molecular_Charge = Zero
      end if
      call Put_DArray('Dipole moment',PrTot,3)
      !call peek_iScalar('xml opened',isopen)
      !if (isopen == 1) then
      call xml_dDump('dipole','Dipole moment','Debye',1,PrTot,3,1)
      !end if
    else if (lPole >= 2) then
      tmp = Zero
      do i=1,Maxlab
        tmp = max(tmp,abs(PrTot(i)))
      end do
      ip_ = 1
      do i=1,maxlab,4
        jMax = min(maxlab-i,3)
        write(u6,'(4X,4(A,A,ES12.4))') (labs(i+j),'=',PrTot(ip_+j)*Fact,j=0,jMax)
        ip_ = ip_+4
      end do
    end if
    if ((lpole >= 2) .and. (lpole <= 4)) then

      ! Transform cartesian moments to multipole moments

      inp = 0
      call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
      inp = 1
      call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
      PrTot(:) = PrNu(1:MaxLab)-PrEl(1,1:MaxLab)

      if (lPole >= 3) then
        write(u6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang**',lPole-1,')'
      else
        write(u6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang)'
      end if
      ip_ = 1
      do i=1,maxlab,4
        jMax = min(maxlab-i,3)
        write(u6,'(4X,4(A,A,ES12.4))') (labs(i+j),'=',PrTot(ip_+j)*Fact,j=0,jMax)
        ip_ = ip_+4
      end do

    end if

  !--------------------------------------------------------------------*
  else if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
  !--------------------------------------------------------------------*

    ! Old style output

    if (lpole == 0) lab5 = ' 0-th'
    if (lpole == 1) lab5 = ' 1-st'
    if (lpole == 2) lab5 = ' 2-nd'
    if (lpole == 3) lab5 = ' 3-rd'
    if (lpole > 3) then
      lab3 = '-th'
      write(lab5,'(i2,a3)') lpole,lab3
    end if
    write(u6,'(//6x,a5,a,3(f12.8,a))') lab5,' cartesian moments: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,76(''-''))')
    sig = -One
    call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,0,ifallorb)
    if (lpole == 1) then
      write(u6,'(6x,76(''-''))')
      write(u6,'(6x,a,3f16.8,3x,a)') 'Total             ',(PrTot(j)*Debye,j=1,3),'Debye'
      call Put_DArray('Dipole moment',PrTot,3)
    end if

    if ((lpole >= 2) .and. (lpole <= 4)) then
      write(u6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,'-pole moment: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
      write(u6,'(6x,76(''-''))')

      ! Transform cartesian moments to multipole moments

      inp = 0
      call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
      inp = 1
      call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)

      ! Print l-pole cartesian moments

      sig = -One
      call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,0,ifallorb)
    end if

  !--------------------------------------------------------------------*
  else
  !--------------------------------------------------------------------*
    PrTot(:) = PrNu(1:MaxLab)-PrEl(1,1:MaxLab)
    if (lpole == 1) call Put_DArray('Dipole moment',PrTot,3)
    if ((lpole >= 2) .and. (lpole <= 4)) then

      ! Transform cartesian moments to multipole moments

      inp = 0
      call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
      inp = 1
      call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
      PrTot(:) = PrNu(1:MaxLab)-PrEl(1,1:MaxLab)
    end if
  !--------------------------------------------------------------------*
  end if ! iPL
  !--------------------------------------------------------------------*
  ! Prop is also called in other programs where MAG_X2C could
  ! be uninitialized, if a test is required in such case
  ! please initialize MAG_X2C to false in related programs
  if (MAG_X2C) StoreInfo = .false.
else if (lab4(1:2) == 'EF') then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! electric field section

  read(oplab,'(a3,i5)') lab3,icen

  if (lab3 == 'EF0') then
    MaxLab = 1
    ! set labels
    labs(1) = '                '
    ! Print cartesian components of the electric field gradient
    ! tensor at the given centre

    if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
      write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric potential:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
      write(u6,'(6x,72(''-''))')
    else
      if (icen == 1) then
        write(u6,'(//6X,A)') 'Electric potential:'
      end if
    end if
    sig = -One
  end if
  if (lab3 == 'EF1') then
    MaxLab = 3
    ! set labels
    labs(1) = '               X'
    labs(2) = '               Y'
    labs(3) = '               Z'
    ! Print cartesian components of the electric field vector
    ! at the given centre

    if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
      write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric field:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
      write(u6,'(6x,72(''-''))')
    else
      if (icen == 1) then
        write(u6,'(//6X,A)') 'Electric field:'
        write(u6,'(5X,6A16)') (labs(i),i=1,MaxLab)
      end if
    end if
    sig = +One
  end if
  if (lab3 == 'EF2') then
    ! Seven components because we add r*r to the "normal" 6
    MaxLab = 7
    ! set labels
    labs(1) = '  (2*XX-YY-ZZ)/2'
    labs(2) = '          1.5*XY'
    labs(3) = '          1.5*XZ'
    labs(4) = '  (2*YY-ZZ-XX)/2'
    labs(5) = '          1.5*YZ'
    labs(6) = '  (2*ZZ-XX-YY)/2'
    ! Actually, r*r is already stored as the 6th component
    ! we now move it to the 7th using the "augmented" arrays
    do j=1,5
      PrElAug(:,j) = PrEl(:,j)
      PrNuAug(j) = PrNu(j)
      labsAug(j) = labs(j)
    end do
    ! Generate the actual value of the 6th component
    PrElAug(:,6) = -PrEl(:,1)-PrEl(:,4)
    PrElAug(:,7) = PrEl(:,6)
    PrNuAug(6) = PrNu(6)
    PrNuAug(7) = Zero
    labsAug(6) = labs(6)
    labsAug(7) = '     RR=XX+YY+ZZ'

    ! Print cartesian components of the electric field gradient
    ! tensor at the given centre

    if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
      write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric field gradient:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
      write(u6,'(6x,78(''-''))')
    else
      if (icen == 1) then
        write(u6,'(//6X,A)') 'Electric field gradient:'
        write(u6,'(5X,6A16)') (labsAug(i),i=1,MaxLab)
      end if
    end if
    sig = +One
  end if
  call mma_allocate(PrTot,maxlab,label='PrTot')

  ! Print the values using "augmented" arrays if needed

  if (Maxlab == 7) then
    call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrElAug,PrNuAug,maxlab,labsAug,PrTot,iPL,icen,ifallorb)
    MaxLab = 6 ! Reset so call to Add_Info is correct!

  else
    call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,icen,ifallorb)
  end if

  if (lab3 == 'EF2') then
    Tmp = PrTot(3)
    PrTot(3) = PrTot(4)
    PrTot(4) = Tmp
    if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) call Print_EigenValues(PrTot,3)
    Tmp = PrTot(3)
    PrTot(3) = PrTot(4)
    PrTot(4) = Tmp
  end if

  ! do not write the different electric field components through Add_Info
  StoreInfo = .false.
else if ((lab4 == 'DMS ') .or. (lab4(1:3) == 'PAM') .or. (lab4(1:3) == 'CNT')) then
  if (lab4 == 'DMS ') then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! diamagnetic shielding section

    read(oplab,'(a4,i2,i2)') lab4,i,icen1

    maxlab = 9
    call mma_allocate(PrTot,maxlab,label='PrTot')

    ! set labels
    labs(1) = '   (YO*YG+ZO*ZG)'
    labs(2) = '         (XO*YG)'
    labs(3) = '         (XO*ZG)'
    labs(4) = '         (YO*XG)'
    labs(5) = '   (XO*XG+ZO*ZG)'
    labs(6) = '         (YO*ZG)'
    labs(7) = '         (ZO*XG)'
    labs(8) = '         (ZO*YG)'
    labs(9) = '   (XO*XG+YO*YG)'

    ! Print cartesian components of the diamagnetic shielding
    ! tensor at the given centre (O=cen1) with the gauge origin
    ! at centre G=cen2

    write(u6,'(//6x,a,i3,1x,a,3(f12.8,a)/1x,a,3(f12.8,a))') ' Diamagnetic shielding:   centre no.',icen1,'(',cen1(1),',',cen1(2), &
                                                            ',',cen1(3),')','                       gauge origin at (',cen2(1), &
                                                            ',',cen2(2),',',cen2(3),')'
    write(u6,'(6x,78(''-''))')
    sig = One
    call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,0,ifallorb)
  end if
  if ((lab4 == 'DMS ') .or. (lab4(1:3) == 'PAM')) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! PAM section ... generate labels for printing

    ilab = 0
    do ix=lpole,0,-1
      do iy=lpole-ix,0,-1
        iz = lpole-ix-iy
        ilab = ilab+1
        lab = '                '
        do izz=lmax,lmax-iz+1,-1
          lab(izz:izz) = 'Z'
        end do
        do iyy=lmax-iz,lmax-iz-iy+1,-1
          lab(iyy:iyy) = 'Y'
        end do
        do ixx=lmax-iz-iy,lmax-lpole+1,-1
          lab(ixx:ixx) = 'X'
        end do
        labs(ilab) = lab
      end do
    end do

    maxlab = ilab
    call mma_allocate(PrTot,maxlab,label='PrTot')

    ! Print cartesian moments

    if (lpole == 0) lab5 = ' 0-th'
    if (lpole == 1) lab5 = ' 1-st'
    if (lpole == 2) lab5 = ' 2-nd'
    if (lpole == 3) lab5 = ' 3-rd'
    if (lpole > 3) then
      lab3 = '-th'
      write(lab5,'(i2,a3)') lpole,lab3
    end if
    write(u6,'(//6x,a5,a,3(f12.8,a))') lab5,' cartesian moments: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,76(''-''))')
    sig = One
    call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,0,ifallorb)

    if ((lpole >= 2) .and. (lpole <= 4)) then
      write(u6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,'-pole moment: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
      write(u6,'(6x,76(''-''))')

      inp = 0
      call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
      inp = 1
      call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)

      ! Print 0-pole PAM integrals (sig=1 in opposite multipole moments)

      sig = One
      call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,0,ifallorb)
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Contact term section

  read(oplab,'(a3,i5)') lab3,icen

  maxlab = 1
  call mma_allocate(PrTot,maxlab,label='PrTot')

  ! set labels
  labs(1) = '      Delta(R-C)'

  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
    write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Contact term:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,78(''-''))')
  else
    if (icen == 1) then
      write(u6,'(//6X,A)') 'Contact term:'
      write(u6,'(5X,6A16)') (labs(i),i=1,MaxLab)
    end if
  end if
  sig = +One
  call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,icen,ifallorb)

  ! do not write the contact term through Add_Info
  StoreInfo = .false.
else if (lab4 == 'VELO') then
  call End1()
  return
else if (lab4 == 'ANGM') then
  call End1()
  return
else

  ! invalid label supplied

  write(u6,'(//1x,a,a8,a//)') ' The label: *',oplab,'* is not a valid label ... stop'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (OpLab(1:3) == 'EF2') then
  iTol = 4
else
  iTol = 5
end if
iTol_E0 = 8
iTol_E1 = Cho_X_GetTol(iTol_E0)
iTol = int(real(iTol*iTol_E1,kind=wp)/real(iTol_E0,kind=wp))
if (StoreInfo) call Add_Info(OpLab,PrTot,maxlab,iTol)
call mma_deallocate(PrTot)
call End1()

return

contains

subroutine End1()
  call mma_deallocate(labs)
  call mma_deallocate(tmat)
  call mma_deallocate(temp)
  call mma_deallocate(PrElAug)
  call mma_deallocate(PrNuAug)
  call mma_deallocate(labsAug)
end subroutine End1

end subroutine Prop
