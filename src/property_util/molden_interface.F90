!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1999, Coen de Graaf                                    *
!               1999, Anders Bernhardsson                              *
!               1999, Roland Lindh                                     *
!***********************************************************************

subroutine Molden_Interface(iUHF,FName,filename)
!***********************************************************************
!                                                                      *
!     Object: to generate MOLDEN input file                            *
!                                                                      *
!     Authors: Coen de Graaf, Anders Bernardsson and R. Lindh, 1999    *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, MolWgh, nBas, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, lIrrep
use Sizes_of_Seward, only: S
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: iUHF
character(len=*), intent(in) :: FName, Filename
integer(kind=iwp) :: i, iAngMx_Valence, iatom, iB, iCntr, iCnttp, icontr, iD, iData, iDeg, iDummy(1), iErr, ii, iIrrep, ik, iPL, &
                     iprim, iRc = 0, iS, isegm, ishell, iWFtype, j, jData, jPL, k, kk, kk_Max, l, Lu_, mdc, MF, nAtom, nB, nData, &
                     nDeg, nOrb(8), nTest, nTot, nTot2
real(kind=wp) :: Check_CMO, Check_Energy, Check_Occupation, coeff, prim
logical(kind=iwp) :: Exists, Found, y_cart, y_sphere
character(len=40) :: VTitle
character(len=8) :: Env
integer(kind=iwp), allocatable :: Cent(:,:), Cent2(:), Cent3(:), ibas_lab(:), Phase(:,:)
real(kind=wp), allocatable :: AdCMO(:), AdCMO_ab(:), AdEor(:), AdEor_ab(:), AdOcc(:), AdOcc_ab(:), C2(:,:), C2_ab(:,:), Coor(:,:), &
                              Mull(:), r_Norm(:), V(:,:), V_ab(:,:), Znuc(:)
character(len=LenIn8+1), allocatable :: gtolabel(:)
character(len=LenIn8), allocatable :: label(:)
character(len=LenIn), allocatable :: AtomLabel(:)
character(len=8), allocatable :: MO_Label(:)
real(kind=wp), parameter :: EorbThr = 50.0_wp
character, parameter :: shelllabel(7) = ['s','p','d','f','g','h','i'], &
                        cNumber(61) = ['1','2','3','4','5','6','7','8','9','0', &
                                       'a','b','c','d','e','f','g','h','i','j', &
                                       'k','l','m','n','o','p','q','r','s','t', &
                                       'u','v','w','x','y','z','A','B','C','D', &
                                       'E','F','G','H','I','J','K','L','M','N', &
                                       'O','P','Q','R','S','T','V','W','X','Y', &
                                       'Z']
integer(kind=iwp), external :: iPrintLevel
real(kind=wp), external :: DblFac
logical(kind=iwp), external :: Reduce_Prt
!integer(kind=iwp), parameter :: MaxOrb_Molden = 400

if (iRc == 1) return

! Do nothing within numerical_gradient
if (SuperName == 'numerical_gradient') return
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the print level

iPL = iPrintLevel(-1)
jPL = iPL
if (Reduce_Prt() .and. (iPL < 3)) jPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
if (MolWgh == 1) then
  if (jPL >= 2) then
    write(u6,*) 'Molden_Interface: Unsupported normalization,Molwgh=1!'
  end if
  iRc = 1
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call GetEnvf('MOLCAS_MOLDEN',Env)
!if ((Env == ' ') .or. (Env == 'OFF')) then
if (Env == 'OFF') then
  if (jPL >= 2) then
    write(u6,*)
    write(u6,*) ' Molden files will not be produced'
    write(u6,*)
  end if
  iRC = 1
  return
end if
!VV: current version of Molden has no clear limit for MaxOrb
!if (MaxOrb > MaxOrb_Molden) then
!  if (jPL >= 2) then
!    write(u6,*)
!    write(u6,*) ' Molden_Interface: W A R N I N G !!!!'
!    write(u6,*)
!    write(u6,*) ' No Molden input file will be generated!'
!    write(u6,*)
!    write(u6,*) ' Calculation exceeds the max number of orbitals allowed for MOLDEN. To change this modify the'
!    write(u6,*) ' parameter MaxOrb_Molden in src/util/molden_interface.f and follow the instructions in Molden'
!    write(u6,*) ' on how to modify the parameter MaxOrb.'
!    write(u6,*)
!  end if
!  iRC = 1
!  return
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
Check_CMO = Zero
Check_Energy = Zero
Check_Occupation = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call f_Inquire('RUNFILE',Exists)
if (.not. Exists) then
  iRC = 1
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the characteristics of all different basis sets,
! provide each atom with a nuclear charge and establish
! a link between an atom and its basis set ---
!
! NOTICE!!!
! This call will also fill Basis_Info and Center_Info.

call mma_allocate(AtomLabel,MxAtom,label='AtomLabel')
call mma_allocate(iBas_Lab,MxAtom,label='iBas_Lab')
call mma_allocate(Coor,3,MxAtom,label='Coor')
call mma_allocate(Znuc,MxAtom,label='Znuc')
call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom)
call Qpg_iArray('nOrb',Found,nData)
if (Found) then
  call Get_iArray('nOrb',nOrb,nData)
else
  nOrb(:) = nBas(:)
end if
call mma_deallocate(iBas_Lab)
!                                                                      *
!***********************************************************************
!                                                                      *
iAngMx_Valence = 0
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    nTest = dbsc(iCnttp)%nVal-1
    iAngMx_Valence = max(iAngMx_Valence,nTest)
  end if
end do
if (iAngMx_Valence > 4) then
  if (jPL >= 2) then
    write(u6,*) 'Sorry, Molden does not know how to handle'
    write(u6,*) 'functions with angular momentum larger than g'
  end if
  call End1()
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Unnormalize contraction coefficients for the valence shells

do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    do l=0,dbsc(iCnttp)%nVal-1
      ishell = dbsc(iCnttp)%iVal+l
      if (Shells(ishell)%Transf .and. (.not. Shells(iShell)%Prjct)) then
        if (jPL >= 2) then
          write(u6,*) 'Sorry, Molden does not support contaminants'
        end if
        call End1()
        return
      end if
      call Unnrmlz(Shells(ishell)%Exp,Shells(ishell)%nExp,Shells(ishell)%pCff,Shells(ishell)%nBasis,l)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute memory requirements and allocate memory

nB = 0
do iS=0,nirrep-1
  nB = nB+nBas(is)
end do
call mma_allocate(MO_Label,nB,label='MO_Label')
call mma_allocate(gtolabel,nB,label='gtolabel')
call mma_allocate(r_Norm,nB,label='r_Norm')
call mma_allocate(Cent,8,nB,label='ICENT')
call mma_allocate(Phase,8,nB,label='IPHASE')
call mma_allocate(Cent2,nB,label='nCENT')
call mma_allocate(Cent3,nB,label='ICENTER')
call mma_allocate(C2,nB,nB,label='CMO2')
call mma_allocate(V,nB,nB,label='VECTOR')
V(:,:) = Zero
if (iUHF == 1) then
  call mma_allocate(C2_ab,nB,nB,label='CMO2')
  call mma_allocate(V_ab,nB,nB,label='VECTOR')
  V_ab(:,:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Open input file for MOLDEN

MF = 9
call molcas_open(MF,filename)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the atom information in MOLDEN format to unit MF

y_cart = .false.
y_sphere = .false.
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) cycle
  do iCntr=1,dbsc(iCnttp)%nCntr
    do l=0,dbsc(iCnttp)%nVal-1
      ! Test for the appearance of cartesian functions with l=2,3,4
      ishell = dbsc(iCnttp)%iVal+l
      if ((l >= 2) .and. (.not. y_cart)) then
        if (.not. Shells(ishell)%Transf) y_cart = .true.
      end if
      if ((l >= 2) .and. (.not. y_sphere)) then
        if (Shells(ishell)%Transf) y_sphere = .true.
      end if
      if (y_sphere .and. y_cart) then
        if (jPL >= 2) then
          write(u6,*)
          write(u6,*) 'Failed to generate input file to MOLDEN'
          write(u6,*) 'No mixing allowed of spherical and cartesian d, f, g-functions'
        end if
        call End2()
        return
      end if
    end do
  end do
end do
write(MF,'(A)') '[Molden Format]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write atomic information

write(MF,'(A)') '[N_Atoms]'
write(MF,*) natom
write(MF,'(A)') '[Atoms] (AU)'
do iatom=1,natom
  write(MF,99) AtomLabel(iatom),iatom,int(Znuc(iatom)),(coor(i,iatom),i=1,3)
end do
if (.not. y_cart) then
  if (S%iAngMx > 1) write(MF,'(A)') '[5D]'
  if (S%iAngMx > 2) write(MF,'(A)') '[7F]'
  if (S%iAngMx > 3) write(MF,'(A)') '[9G]'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out charges and dipole moments

call qpg_dArray('Mulliken Charge',Found,nData)
if (Found) then
  write(MF,'(A)') '[Charge] (Mulliken)'
  call mma_allocate(Mull,nData,label='Mull')
  call Get_dArray('Mulliken Charge',Mull,nData)

  iData = 0
  jData = 0
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%pChrg) cycle
    do iCntr=1,dbsc(iCnttp)%nCntr
      iData = iData+1
      mdc = iCntr+dbsc(iCnttp)%mdci
      nDeg = nIrrep/dc(mdc)%nStab
      do iDeg=1,nDeg
        jData = jData+1
        write(MF,*) Mull(iData)
      end do
    end do
  end do
  call mma_deallocate(Mull)
  if (iData /= nData) then
    write(u6,*) 'Molden_Interface: iData.ne.nData'
    write(u6,*) 'iData,nData=',iData,nData
    call Abend()
  end if
  if (jData /= nAtom) then
    write(u6,*) 'Molden_Interface: jData.ne.nAtom'
    write(u6,*) 'jData,nAtom=',jData,nAtom
    call Abend()
  end if
end if
!write(MF,'(A)') '[NDIPOLE]'
!write(MF,'(A)') '[DIPOLE]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write Gaussian basis set information to MOLDEN input file.

write(MF,'(A)') '[GTO] (AU)'

! Read exponents and contraction coefficients of each unique basis.
! Write the present basis set (iCnttp) to the molden.input file for
! the appropriate atoms.
! Moreover, a list is constructed which contains a label for each
! GTO (gtolabel). This list follows the MOLDEN order of GTOs.
! Later this list will be used in the transformation of sabf (the
! symmetry adapted basis functions).

iatom = 0
mdc = 0
kk = 0

do iCnttp=1,nCnttp             ! loop over unique basis sets
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) cycle

  do iCntr=1,dbsc(iCnttp)%nCntr  ! loop over sym. unique centers
    mdc = mdc+1
    nDeg = nIrrep/dc(mdc)%nStab
    do iDeg=1,nDeg             ! loop over centers
      iAtom = iAtom+1
      write(MF,'(I4)') iAtom

      do l=0,dbsc(iCnttp)%nVal-1
        ishell = dbsc(iCnttp)%iVal+l
        if (Shells(iShell)%nBasis > size(cNumber)) then
          write(u6,*) 'Interf: too many contracted functions!'
          write(u6,*) 'nBasis(iShell)=',Shells(iShell)%nBasis
          call Abend()
        end if

        ! Iterate over each contracted GTO

        do icontr=1,Shells(ishell)%nBasis

          ! Find the number of exponents with non-zero exponents

          isegm = 0
          do iprim=1,Shells(ishell)%nExp
            coeff = Shells(ishell)%pCff(iprim,icontr)
            if (coeff /= Zero) then
              isegm = isegm+1
            end if
          end do

          write(MF,'(3x,A1,I4)') shelllabel(l+1),isegm

          ! Write exponents and contraction coefficients.

          do iprim=1,Shells(ishell)%nExp
            coeff = Shells(ishell)%pCff(iprim,icontr)
            prim = Shells(ishell)%Exp(iprim)
            if (coeff /= Zero) then
              write(MF,'(ES17.9,ES17.9)') prim,coeff
            end if
          end do

          ! Construction of gtolabel
          ! Molden order: for p-functions: px(1), py(1), pz(1),
          !                                px(2), py(2), pz(2), etc.
          ! for d-, and f-functions: similar

          if (l == 0) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'01s     '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
          end if
          if (l == 1) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02px    '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02py    '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02pz    '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
          end if
          if ((l == 2) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d00   '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
          end if
          if ((l == 2) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d020000 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000200 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000002 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,2)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010100 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010001 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000101 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,1)
            Cent3(kk) = iAtom
          end if
          if ((l == 3) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f00   '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
          end if
          if ((l == 3) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f030000 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000300 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000003 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,3)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010200 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020100 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020001 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010002 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,2)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000102 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,2)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000201 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010101 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,1)
            Cent3(kk) = iAtom
          end if
          if ((l == 4) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g00   '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '//cNumber(icontr)
            r_Norm(kk) = One
            Cent3(kk) = iAtom
          end if
          if ((l == 4) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g040000 '//cNumber(icontr)
            r_Norm(kk) = CC(4,0,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000400 '//cNumber(icontr)
            r_Norm(kk) = CC(0,4,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000004 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,4)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030100 '//cNumber(icontr)
            r_Norm(kk) = CC(3,1,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030001 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010300 '//cNumber(icontr)
            r_Norm(kk) = CC(1,3,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000301 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010003 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,3)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000103 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,3)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020200 '//cNumber(icontr)
            r_Norm(kk) = CC(2,2,0)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020002 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,2)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000202 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,2)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020101 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010201 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,1)
            Cent3(kk) = iAtom
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010102 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,2)
            Cent3(kk) = iAtom
          end if
        end do
      end do
      write(MF,'(A)') ' '
    end do
  end do
end do
kk_Max = kk
if (nB > kk_max) then
  if (jPL >= 2) then
    write(u6,*) 'Molden_Interface: nB.gt.kk_max'
    write(u6,*) 'nB,kk_Max=',nB,kk_Max
  end if
  call End2()
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nTot = 0
nTot2 = 0
do iS=0,nIrrep-1
  nTot = nTot+nBas(iS)
  nTot2 = nTot2+nBas(iS)**2
end do
call mma_allocate(AdOcc,nTot,label='Occ')
call mma_allocate(AdEor,nTot,label='Eor')
call mma_allocate(AdCMO,nTot2,label='CMO')
AdOcc(:) = Zero
AdEor(:) = Zero
AdCMO(:) = Zero
if (iUHF == 1) then
  call mma_allocate(AdOcc_ab,nTot,label='Occ')
  call mma_allocate(AdEor_ab,nTot,label='Eor')
  call mma_allocate(AdCMO_ab,nTot2,label='CMO')
  AdOcc_ab(:) = Zero
  AdEor_ab(:) = Zero
  AdCMO_ab(:) = Zero
else
  call mma_allocate(AdOcc_ab,0,label='Occ')
  call mma_allocate(AdEor_ab,0,label='Eor')
  call mma_allocate(AdCMO_ab,0,label='CMO')
end if

! Read HF CMOs from file

Lu_ = 75
call RdVec_(FName,Lu_,'COE',iUHF,nIrrep,nBas,nBas,AdCMO,AdCMO_ab,AdOcc,AdOcc_ab,AdEor,AdEor_ab,iDummy,VTitle,1,iErr,iWFtype)

! Get the coeff. of sym adapted basis functions (C2)

call Dens_IF_SCF(C2,AdCMO,'F')
if (iUHF == 1) then
  call Dens_IF_SCF(C2_ab,AdCMO_ab,'F')
end if
call mma_deallocate(AdCMO)
call mma_deallocate(AdCMO_ab)
!                                                                      *
!***********************************************************************
!                                                                      *
!  Back 'transformation' of the symmetry adapted basis functions.
!  Probably somewhat clumsy, but it seems to work. If someone
!  knows a more elegant way to do it, please improve this part!
!
!  PART 1: Obtain symmetry information (soout), construct a label
!          for each sabf, which will be used in part 2 to find the
!          corresponding GTO in the MOLDEN list by comparing with
!          gtolabel
!
!  nB     --- Total number of contracted basis functions
!  cent2  --- degeneracy of a basis function
!  Cent   --- centres over which the basis function is delocalized
!  Phase  --- phase of the AO in the linear combination

call mma_allocate(label,MaxBfn+MaxBfn_Aux,label='label')
Phase(:,:) = 0
Cent(:,:) = 0
call SOout(label,Cent,Phase)
do iContr=1,nB
  Cent2(iContr) = 0
  do k=1,8
    if (Cent(k,iContr) /= 0) Cent2(iContr) = Cent2(iContr)+1
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!
! Part 2: -Take a MOLCAS symmetry functions (loop i)
!         -Find the corresponding label in the MOLDEN list (loop j)
!         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
!          by the appropriate factor (Phase),and divide by the number of
!          centres over which the sabf is delocalized (Cent3).
!         -The vectors are copied by rows!

! loop over MOLCAS symmetry functions
i = 0
ik = 0
do iIrrep=0,nIrrep-1

  do iB=1,nBas(iIrrep)
    i = i+1

    if (iB == 1) then
      ik = 1
    else
      if (label(i-1) == label(i)) then
        ik = ik+1
      else
        ik = 1
      end if
    end if
    if (ik > size(cNumber)) then
      write(u6,*) 'Molden_Interface: ik > size(cNumber)'
      write(u6,*) 'ik,size(cNumber)=',ik,size(cNumber)
      write(u6,*) 'List of labels'
      do iD=1,nB
        write(u6,'(A)') Label(i)
      end do
      call End3()
      return
    end if

    write(MO_Label(i),'(I5,A3)') iB,lirrep(iIrrep)

    do j=1,nB

      if (gtolabel(j) == label(i)//cNumber(ik)) then
        do k=1,8
          if (Cent(k,i) == Cent3(j)) then
            do ii=1,nB
              if (MolWgh == 0) then
                V(j,ii) = V(j,ii)+(C2(i,ii)*r_Norm(j))*real(Phase(k,i),kind=wp)/real(Cent2(i),kind=wp)
                if (iUHF == 1) V_ab(j,ii) = V_ab(j,ii)+(C2_ab(i,ii)*r_Norm(j))*real(Phase(k,i),kind=wp)/real(Cent2(i),kind=wp)
              else
                V(j,ii) = V(j,ii)+(C2(i,ii)*r_Norm(j))*real(Phase(k,i),kind=wp)/sqrt(real(Cent2(i),kind=wp))
                if (iUHF == 1) V_ab(j,ii) = V_ab(j,ii)+(C2_ab(i,ii)*r_Norm(j))*real(Phase(k,i),kind=wp)/sqrt(real(Cent2(i),kind=wp))
              end if
            end do
          end if
        end do
      end if
    end do
  end do
end do
call mma_deallocate(label)
!                                                                      *
!***********************************************************************
!                                                                      *
!  Dump vector in the molden.input file

write(MF,'(A)') '[MO]'
do i=1,nB
  if (AdEOr(i) <= EorbThr) then
    write(MF,'(A,A)') 'Sym= ',MO_Label(i)
    write(MF,103) AdEOr(i)
    write(MF,'(A)') 'Spin= Alpha'
    write(MF,104) AdOcc(i)
    if (AdEOr(i) < Zero) then
      Check_Energy = Check_Energy+AdEOr(i)*real(i,kind=wp)
    end if
    Check_Occupation = Check_Occupation+AdOcc(i)*real(i,kind=wp)
    do j=1,nB
      write(MF,100) j,V(j,i)
      Check_CMO = Check_CMO+V(j,i)**2
    end do
  end if

  if (iUHF == 1) then
    if (AdEOr_ab(i) <= EorbThr) then
      write(MF,'(A,A)') 'Sym= ',MO_Label(i)
      write(MF,103) AdEOr_ab(i)
      write(MF,'(A)') 'Spin= Beta'
      write(MF,104) AdOcc_ab(i)
      Check_Energy = Check_Energy+AdEOr_ab(i)*real(i,kind=wp)
      Check_Occupation = Check_Occupation+AdOcc_ab(i)*real(i,kind=wp)
      do j=1,nB
        write(MF,100) j,V_ab(j,i)
        Check_CMO = Check_CMO+V_ab(j,i)**2
      end do
    end if
  end if

end do
!                                                                      *
!***********************************************************************
!                                                                      *
!if (Env /= 'OFF') then
!  if (iUHF < 2) call Add_Info('MOLDEN_CMO',Check_CMO,1,2)
!  call Add_Info('MOLDEN_Occupation',Check_Occupation,1,2)
!  call Add_Info('MOLDEN_Energy',Check_Energy,1,2)
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%Aux) .and. (.not. dbsc(iCnttp)%Frag)) then
    do l=0,dbsc(iCnttp)%nVal-1
      ishell = dbsc(iCnttp)%iVal+l
      call Unnrmlz2(Shells(ishell)%Exp,Shells(ishell)%nExp,Shells(ishell)%pCff,Shells(ishell)%nBasis,l)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (jPL >= 2) then
  write(u6,*)
  write(u6,'(6X,A)') 'Input file to MOLDEN was generated!'
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call End3()

return

! -------------FORMATS-------------
99 format(A,2(3x,I4),5x,3F16.8,3I4)
100 format(I4,3x,F16.8)
103 format('Ene= ',F10.4)
104 format('Occup= ',F10.5)

contains

function CC(ix,iy,iz)
  real(kind=wp) :: CC
  integer(kind=iwp), intent(in) :: ix, iy, iz
  CC = sqrt(DblFac(2*ix-1)*DblFac(2*iy-1)*DblFac(2*iz-1))
end function CC

subroutine End1()
  call mma_deallocate(AtomLabel)
  call mma_deallocate(Coor)
  call mma_deallocate(Znuc)
  call ClsSew()
end subroutine End1

subroutine End2()
  call mma_deallocate(MO_Label)
  call mma_deallocate(gtolabel)
  call mma_deallocate(r_Norm)
  call mma_deallocate(Cent)
  call mma_deallocate(Phase)
  call mma_deallocate(Cent2)
  call mma_deallocate(Cent3)
  call mma_deallocate(C2)
  call mma_deallocate(V)
  if (iUHF == 1) then
    call mma_deallocate(C2_ab)
    call mma_deallocate(V_ab)
  end if
  close(MF)
  call End1()
end subroutine End2

subroutine End3()
  call mma_deallocate(AdOcc)
  call mma_deallocate(AdEor)
  call mma_deallocate(AdOcc_ab)
  call mma_deallocate(AdEor_ab)
  call End2()
end subroutine End3

end subroutine Molden_Interface
