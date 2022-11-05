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
!               2018, Jesper Norell                                    *
!***********************************************************************

subroutine Molden_DysOrb(filename,ENE,OCC,CMO,NDO,NZ)
!***********************************************************************
!                                                                      *
!     Object: to generate MOLDEN files for Dyson orbitals              *
!                                                                      *
!     Modified from molden_interface by Jesper Norell 2018             *
!     (Since requirements for e.g. number of orbitals and              *
!     treament of symmetry differs from "normal" MOs)                  *
!                                                                      *
!     Authors: Coen de Graaf, Anders Bernardsson and R. Lindh, 1999    *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, nBas, MolWgh, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, lIrrep
use Sizes_of_Seward, only: S
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: NDO, NZ
character(len=*), intent(in) :: filename
real(kind=wp), intent(in) :: ENE(NDO), OCC(NDO), CMO(NZ*NDO)
integer(kind=iwp) :: i, iAngMx_Valence, iatom, iB, iBtot, iCntr, iCnttp, icontr, iCount, iDeg, idx, iGTO, iIrrep, iLabel, iPL, &
                     iprim, iRc = 0, iS, isegm, ishell, iSymcent, jPL, k, kk, kk_Max, l, mdc, MF, nAtom, nB, nData, nDeg, nOrb(8), &
                     nTest, nTot, nTot2
real(kind=wp) :: COEF, coeff, prim
logical(kind=iwp) :: Exists, Found, y_cart, y_sphere
character(len=8) :: Env
integer(kind=iwp), allocatable :: Cent(:,:), Cent2(:), ibas_lab(:), Phase(:,:)
real(kind=wp), allocatable :: Coor(:,:), DESYM(:,:), r_Norm(:), Znuc(:)
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
!if ((Env == ' ') .or. (Env == 'OFF')) Then
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
        call End3()
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

  do iCntr=1,dbsc(iCnttp)%nCntr   ! loop over sym. unique centers
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
              write(MF,'(E17.9,E17.9)') prim,coeff
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
          end if
          if (l == 1) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02px    '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02py    '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'02pz    '//cNumber(icontr)
            r_Norm(kk) = One
          end if
          if ((l == 2) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d00   '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '//cNumber(icontr)
            r_Norm(kk) = One
          end if
          if ((l == 2) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d020000 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000200 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000002 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,2)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010100 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d010001 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'d000101 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,1)
          end if
          if ((l == 3) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f00   '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '//cNumber(icontr)
            r_Norm(kk) = One
          end if
          if ((l == 3) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f030000 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000300 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000003 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,3)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010200 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020100 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f020001 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010002 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,2)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000102 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,2)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f000201 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'f010101 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,1)
          end if
          if ((l == 4) .and. (.not. y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g00   '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '//cNumber(icontr)
            r_Norm(kk) = One
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '//cNumber(icontr)
            r_Norm(kk) = One
          end if
          if ((l == 4) .and. (y_cart)) then
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g040000 '//cNumber(icontr)
            r_Norm(kk) = CC(4,0,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000400 '//cNumber(icontr)
            r_Norm(kk) = CC(0,4,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000004 '//cNumber(icontr)
            r_Norm(kk) = CC(0,0,4)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030100 '//cNumber(icontr)
            r_Norm(kk) = CC(3,1,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g030001 '//cNumber(icontr)
            r_Norm(kk) = CC(3,0,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010300 '//cNumber(icontr)
            r_Norm(kk) = CC(1,3,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000301 '//cNumber(icontr)
            r_Norm(kk) = CC(0,3,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010003 '//cNumber(icontr)
            r_Norm(kk) = CC(1,0,3)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000103 '//cNumber(icontr)
            r_Norm(kk) = CC(0,1,3)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020200 '//cNumber(icontr)
            r_Norm(kk) = CC(2,2,0)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020002 '//cNumber(icontr)
            r_Norm(kk) = CC(2,0,2)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g000202 '//cNumber(icontr)
            r_Norm(kk) = CC(0,2,2)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g020101 '//cNumber(icontr)
            r_Norm(kk) = CC(2,1,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010201 '//cNumber(icontr)
            r_Norm(kk) = CC(1,2,1)
            kk = kk+1
            gtolabel(kk) = AtomLabel(iAtom)//'g010102 '//cNumber(icontr)
            r_Norm(kk) = CC(1,1,2)
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
!  Construct a matrix which inverses the symmetry transformation
!  of the basis fuctions (DESYM)
!
!  NOTE: Also necessary for no symmetry calculations
!  to properly keep track of i.e. px vs py vs pz

call mma_allocate(DESYM,NZ,NZ,label='DESYM')
DESYM(:,:) = Zero
iBtot = 0
iCount = 0
do iIrrep=0,nIrrep-1 ! For all the irreps of symmetrized functions

  do iB=1,nBas(iIrrep) ! For each symmetrized function in the irrep
    iBtot = iBtot+1
    iSymcent = 1

    if (iB == 1) then ! Counter within functions with the same label (e.g. H1 1s 1, H1 1s 2, ...)
      iCount = 1
    else
      if (label(iBtot-1) == label(iBtot)) then
        iCount = iCount+1
      else
        iCount = 1
      end if
    end if

    write(MO_Label(iBtot),'(A3)') lirrep(iIrrep) ! Save label for later printing

    do iGTO=1,nB ! Find matching Gaussian (GTO)

      if (gtolabel(iGTO) == label(iBtot)//cNumber(iCount)) then
        ! GTO Found
        coeff = real(Phase(iSymcent,iBtot),kind=wp)*r_Norm(iGTO)/sqrt(real(Cent2(iBtot),kind=wp))
        DESYM(iBtot,iGTO) = DESYM(iBtot,iGTO)+coeff
        iSymcent = iSymcent+1 ! Keep track of symmetry centre
      end if ! GTO finding

    end do ! iGTO=1,nB
  end do ! iB=1,nBas(iIrrep)
end do ! iIrrep=0,nIrrep-1
call mma_deallocate(label)

!                                                                      *
!***********************************************************************
!                                                                      *
!  Dump vector in the molden.input file

write(MF,'(A)') '[MO]'

! For all Dyson orbitals
do I=1,NDO

  ! Symmetry = perform a backtransformation from symmetrized to
  ! original basis functions
  ! and to properly split up e.g. px vs py vs pz compoments
  ! also for no symmetry calculations

  ! Find the correct symmetry label from the symmetrized functions
  idx = (I-1)*NB+1
  iLabel = maxloc(abs(CMO(idx:idx+nB-1)),1)
  write(MF,'(A,I0,A)') 'Sym= ',I,MO_Label(iLabel)
  write(MF,103) ENE(I)
  write(MF,'(A)') 'Spin= Alpha'
  write(MF,104) OCC(I)

  ! Now backtransform the orbitals
  do iGTO=1,NB ! For each molden GTO
    COEF = 0
    do iB=1,NB ! Gather contributions from all symmetrized basis functions
      COEF = COEF+DESYM(iB,iGTO)*CMO((I-1)*NB+iB)
    end do
    write(MF,100) iGTO,COEF
  end do

end do ! I=1,NDO (i.e. Dyson orbitals)
call mma_deallocate(DESYM)

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
  close(MF)
  call End1()
end subroutine End2

subroutine End3()
  call mma_deallocate(MO_Label)
  call mma_deallocate(gtolabel)
  call mma_deallocate(r_Norm)
  call mma_deallocate(Cent)
  call mma_deallocate(Phase)
  call mma_deallocate(Cent2)
  call End2()
end subroutine End3

end subroutine Molden_DysOrb
