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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  Qfread
!
!> @brief
!>   Read in all data that comes from external Molcas routines and prepare various quantities, such as MME
!> @author A. Ohrn
!>
!> @details
!> This subroutine handles the interaction with the rest of Molcas.
!> Here orbitals and various matrices are stored and to some extent
!> modified to our purpose. We call on this subroutine even when
!> we are running classical stuff. The reason for this is that we
!> wish to collect some numbers, but really, we could easily have
!> constructed thing differently so that this subroutine only would
!> be called when quantum-classical stuff is running. For RASSI
!> implementation, we also read in and transform the transition density
!> matrix. At the end, we also call the routines that make the
!> multicenter multipole expansion.
!>
!> @note
!> Seward is mandatory for both SCF and RASSI. For SCF also,
!> Motra, Averd; for RASSI also, RASSCF and RASSI.
!>
!> @param[out] nAtomsCC Atoms in solvent
!> @param[out] Coord    Unique coordinates of the atoms in the molecule in the QM-region
!> @param[out] nBas     Number of basis functions in QM-region
!> @param[out] nBasCC   Like nBas but for a solvent molecule
!> @param[out] nCnC_C   Like nCnC, but for solvent
!> @param[out] nntyp    Number of basis-function types.
!> @param[out] nOcc     The total number of basis functions that belong to a certain basis-function type.
!> @param[out] natyp    Number of atoms of the i:th basis-function type
!***********************************************************************

!******JoseMEP the last three variables are included to the MEP calculation
subroutine Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBasCC,nCnC_C,nOcc,natyp,nntyp)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "integral.fh"
#include "WrkSpc.fh"
#include "warnings.h"
integer(kind=iwp) :: iQ_Atoms, nAtomsCC, nBas(MxSym), nBasCC(1), nCnC_C(MxBasC), nOcc(MxBas), natyp(MxAt), nntyp
real(kind=wp) :: Coord(MxAt*3)
integer(kind=iwp) :: i, iAtom, iBas, icont, iDummy(1), iErr, iLu, ind, indold, iOe, iold, ipACC, ipC, ipC_C, ipE, ipE_C, iWarn, &
                     ix, j, jnd, k, kaunter, kk, Kmax, kold, l, lLine, m, na, nACCSizeC, nACCSizeQ, nnaa, nntypC, nSize, nSym, &
                     nSymCC
real(kind=wp) :: ChgeCC(3), CoordCC(3*3), Dummy(1)
character(len=120) :: BlLine, Line, StLine
character(len=100) :: OrbName, Title
character(len=10) :: WhatGet
integer(kind=iwp), allocatable :: iC_Icon(:,:), Icon(:,:), natypC(:), nfSh(:,:), nSh(:)
real(kind=wp), allocatable :: Chge(:), Cmo(:), Cmo_S(:), Occu(:)
integer(kind=iwp), external :: IsFreeUnit

!----------------------------------------------------------------------*
! Enter                                                                *
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
! Print the joblabel. It is obtained in get_input.f                    *
!----------------------------------------------------------------------*
if (ATitle) then
  lLine = len(Line)
  do i=1,lLine
    BlLine(i:i) = ' '
    StLine(i:i) = '*'
  end do
  write(u6,*)
  do i=1,6
    Line = BlLine
    if ((i == 1) .or. (i == 6)) Line = StLine
    if (i == 3) Line = 'Project:'
    if (i == 4) write(Line,'(A72)') JobLab
    call Center_Text(Line)
    write(u6,*) '*'//Line//'*'
  end do
  write(u6,*)
end if
write(u6,*) 'Auxiliary data being read and pre-processed.'
!----------------------------------------------------------------------*
! Collect some data from RUNFILE about the QM-region molecule.         *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
if (nSym /= 1) then !A minor restriction, no symmetry allowed, i.e. nSym=1.
  write(u6,*)
  write(u6,*) ' QmStat does not run with symmetry!'
  write(u6,*) ' The perturbation from the solvent breaks all symmetry.'
  call Quit(_RC_GENERAL_ERROR_)
end if
call Get_iScalar('Unique atoms',iQ_Atoms)
if (iQ_Atoms > MxAt) then
  write(u6,*)
  write(u6,*) 'Maximum number of atoms exceeded. Increase MxAt in maxi.fh in QmStat source directory.'
  call Quit(_RC_GENERAL_ERROR_)
end if
call mma_allocate(Chge,iQ_Atoms,label='Chge')
call Get_dArray('Nuclear charge',Chge,iQ_Atoms)
call Get_dArray('Unique Coordinates',Coord,3*iQ_Atoms)
call Get_dArray('Center of Mass',CT,3)
call Get_iArray('nBas',nBas,nSym)
if (nBas(1) > MxBas) then
  write(u6,*)
  write(u6,*) 'Maximum number of basis functions exceeded. Increase MxBas in maxi.fh in QmStat source directory.'
  call Quit(_RC_GENERAL_ERROR_)
end if

!----------------------------------------------------------------------*
! Print elementary information about molecule.                         *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,*) '     ------------------------------'
write(u6,*) '     |       QM-region data       |'
write(u6,*) '     ------------------------------'
write(u6,*)
write(u6,'(A,15X,I5)') '      Number of basis functions:',(nBas(i),i=1,nSym)
write(u6,'(A,3F10.6)') '      Centre of mass =',(CT(kk),kk=1,3)
call PrCoor()

!----------------------------------------------------------------------*
! If the Qmtype is SCF we now want the orbitals. On the other hand if  *
! we are running RASSI, another route must be taken, hence here we     *
! inquire which QM-method is used.                                     *
!----------------------------------------------------------------------*
if (QmType(1:3) == 'SCF') then
  !--------------------------------------------------------------------*
  !    SSS    CCC  FFFFFF                                              *
  !   SS    CC     FF                                                  *
  !    SS  CC      FFFF                                                *
  !     SS  CC     FF                                                  *
  !   SSS    CCC   FF                                                  *
  !--------------------------------------------------------------------*

  !--------------------------------------------------------------------*
  ! Print information about orbitals and store coefficients in new     *
  ! variable.                                                          *
  !--------------------------------------------------------------------*
  iLu = 15
  iLu = IsFreeUnit(iLu)
  write(OrbName,'(A)') 'AVEORB'
  write(WhatGet,'(A)') 'CO'
  iWarn = 1
  call mma_allocate(Cmo,MxBas**2,label='Cmo')
  call mma_allocate(Occu,MxBas,label='Occu')
  call RdVec(OrbName,iLu,WhatGet,nSym,nBas,nBas,Cmo,Occu,Dummy,iDummy,Title,iWarn,iErr)
  if (iErr /= 0) then
    write(u6,*)
    write(u6,*) 'Error when reading AVEORB'
    write(u6,*)
    call Quit(_RC_IO_ERROR_READ_)
  end if
  kaunter = 0
  call GetMem('OrbCoeffQ','Allo','Real',iV1,iOrb(1)*nBas(1))
  do j=1,iOrb(1)
    do k=1,nBas(1)
      Work(iV1+kaunter) = Cmo(k+(j-1)*nBas(1))
      kaunter = kaunter+1
    end do
  end do
  !write(u6,'(A,I4)') '      Number of Orbitals:',iOrb(1)
  call Get_dScalar('PotNuc',PotNuc)

else if (QmType(1:4) == 'RASS') then
  !--------------------------------------------------------------------*
  !   RRRR      AA      SSS   SSS   II                                 *
  !   RR  R    A  A    SS    SS     II                                 *
  !   RR R    AAAAAA    SS    SS    II                                 *
  !   RRR     AA  AA     SS    SS   II                                 *
  !   RR R    AA  AA   SSS   SSS    II                                 *
  !--------------------------------------------------------------------*
  call TdmTrans(nBas)
end if
write(u6,*)
write(u6,*)
!----------------------------------------------------------------------*
! Compute various information about system. This we use for computing  *
! integrals later. GiveMeInfo collects stuff from seward, sometime with*
! some recomputations.                                                 *
!----------------------------------------------------------------------*
call mma_allocate(Icon,MxAt,MxPrCon,label='Icon')
call mma_allocate(nSh,MxAt,label='nSh')
call mma_allocate(nfSh,MxAt,MxAngqNr,label='nfSh')
call GiveMeInfo(nntyp,natyp,BasOri,Icon,nPrimus,nBA_Q,nCBoA_Q,nBonA_Q,ipE,ipC,nSh,nfSh,nSize,iPrint,MxAt,MxPrCon,MxBas,MxAngqNr, &
                ipACC,nACCSizeQ)
iBas = 0
iAtom = 0
kold = 1
iold = 1
indold = 0
do i=1,nntyp
  nOcc(i) = 0
end do
do i=1,nntyp
  na = natyp(i)
  do j=1,na
    ind = 0
    jnd = 0
    iAtom = iAtom+1
    ChaNuc(iAtom) = Chge(iAtom)
    info_atom(iAtom) = int(Chge(iAtom))
    do k=1,nSh(i)
      nnaa = nfSh(i,k)
      do l=1,nnaa
        ibas = ibas+1
        indold = indold+1
        nOcc(i) = nOcc(i)+2*k-1
        ind = ind+1
        iQang(ibas) = k
        icont = Icon(i,ind)
        iCharOnBasQ(ibas) = int(Chge(iAtom))
        do ix=1,2*k-1  !Here we construct an array of
          if (k /= kold) then !indices which is used to put right
            if (i /= iold) then !AO-overlap in right matrix pos.
              Indold = Indold+nfSh(iold,kold)*(2*kold-2)
              iold = i
            else
              Indold = Indold+nfSh(i,kold)*(2*kold-2)
            end if
            kold = k
          end if
          iWoGehenQ(ibas,ix) = indold+nnaa*(ix-1)
        end do
        do m=1,icont
          jnd = jnd+1
          alfa(ibas,m) = Work(ipE+i-1+MxAt*(jnd-1))
          cont(ibas,m) = Work(ipC+i-1+MxAt*(jnd-1))
        end do
      end do
    end do
  end do
end do
Kmax = ibas
call dcopy_(nACCSizeQ,Work(ipACC),1,Trans,1)
! Now we do not need them, so deallocate
call mma_deallocate(Chge)
call mma_deallocate(Icon)
call GetMem('AccTransa','Free','Real',ipACC,nACCSizeQ)
call GetMem('Exponents','Free','Real',ipE,nSize*MxAt)
call GetMem('ContrCoef','Free','Real',ipC,nSize*MxAt)

!----------------------------------------------------------------------*
! Obtain and print information about solvent. This requires a renaming *
! of the runfile.                                                      *
!----------------------------------------------------------------------*
call NameRun('WRUNFIL')
call Get_iScalar('nSym',nSymCC)
if (nSymCC /= 1) then
  write(u6,*)
  write(u6,*) ' QmStat does not run with symmetry!'
  call Quit(_RC_GENERAL_ERROR_)
end if
call Get_iScalar('Unique atoms',nAtomsCC)
if (nAtomsCC /= 3) then
  write(u6,*)
  write(u6,*) 'Now now... what strange solvent molecule do you try to feed QmStat with?'
  call Quit(_RC_GENERAL_ERROR_)
end if
call Get_dArray('Nuclear charge',ChgeCC,nAtomsCC)
call Get_dArray('Unique Coordinates',CoordCC,3*nAtomsCC)
call Get_iArray('nBas',nBasCC,nSymCC)
if (nBasCC(1) > MxBasC) then
  write(u6,*)
  write(u6,*) 'Number of solvent molecule basis functions exceeded. Increase MxBasC in maxi.fh in QmStat source directory.'
  call Quit(_RC_GENERAL_ERROR_)
end if
write(u6,*)
write(u6,*) '     ------------------------------'
write(u6,*) '     |   Solvent molecule data     |'
write(u6,*) '     ------------------------------'
write(u6,*)
write(u6,'(A,15X,I5)') '      Number of basis functions:',(nBasCC(i),i=1,nSymCC)
call PrCoor()

! Collect information about the solvent orbitals.

iLu = 16
iLu = IsFreeUnit(iLu)
write(OrbName,'(A)') 'SOLORB'
write(WhatGet,'(A)') 'CE'
iWarn = 1
call GetMem('OrbitalEnergy','Allo','Real',iOe,sum(nBasCC))
call mma_allocate(Cmo_S,MxBas**2,label='Cmo_S')
call RdVec(OrbName,iLu,WhatGet,nSymCC,nBasCC,nBasCC,Cmo_S,Dummy,Work(iOe),iDummy,Title,iWarn,iErr)
do i=1,iOrb(2)
  c_orbene(i) = Work(iOe+i-1)
end do
call GetMem('OrbitalEnergy','Free','Real',iOe,sum(nBasCC))

! We should not need two solvent orbital vectors, so this should
! be removed when the orbital rotation routine is fixed.

do j=1,iOrb(2)
  do k=1,nBasCC(1)
    V3(k,j) = Cmo_S(k+(j-1)*nBasCC(1))
  end do
end do
call mma_deallocate(Cmo_S)
!write(u6,'(A,I4)') '      Number of Orbitals:',iOrb(2)
write(u6,*)
write(u6,*)
!----------------------------------------------------------------------*
! And now basis set information.                                       *
!----------------------------------------------------------------------*
call mma_allocate(natypC,MxAt,label='natypC')
call mma_allocate(iC_Icon,MxAt,MxPrCon,label='iC_Icon')
call GiveMeInfo(nntypC,natypC,SavOri,iC_Icon,mPrimus,nBA_C,nCBoA_C,nBonA_C,ipE_C,ipC_C,nSh,nfSh,nSize,iPrint,MxAt,MxPrCon,MxBas, &
                MxAngqNr,ipACC,nACCSizeC)
iBas = 0
iAtom = 0
kold = 1
iold = 1
indold = 0
do i=1,nntypC !Like the corresponding thing above for the QM-region.
  na = natypC(i)
  do j=1,na
    ind = 0
    jnd = 0
    iAtom = iAtom+1
    do k=1,nSh(i)
      nnaa = nfSh(i,k)
      do l=1,nnaa
        ibas = ibas+1
        indold = indold+1
        nCnC_C(ibas) = nnaa
        ind = ind+1
        icont = iC_Icon(i,ind)
        iqn(ibas) = k
        iCharOnBasC(ibas) = int(ChgeCC(iAtom))
        ! Here we construct an array of indices which is used to put right AO-overlap in right matrix pos.
        do ix=1,2*k-1
          if (k /= kold) then
            if (i /= iold) then
              Indold = Indold+nfSh(iold,kold)*(2*kold-2)
              iold = i
            else
              Indold = Indold+nfSh(i,kold)*(2*kold-2)
            end if
            kold = k
          end if
          iWoGehenC(ibas,ix) = indold+nnaa*(ix-1)
        end do
        do m=1,icont
          jnd = jnd+1
          beta(ibas,m) = Work(ipE_C+i-1+MxAt*(jnd-1))
          dont(ibas,m) = Work(ipC_C+i-1+MxAt*(jnd-1))
        end do
      end do
    end do
  end do
end do
Lmax = ibas
if (nACCSizeC > nACCSizeQ) then
  call dcopy_(nACCSizeC,Work(ipACC),1,Trans,1)
end if
! Now we do not need them, so deallocate.
call mma_deallocate(natypC)
call mma_deallocate(nSh)
call mma_deallocate(nfSh)
call mma_deallocate(iC_Icon)
call GetMem('AccTransa','Free','Real',ipACC,nACCSizeC)
call GetMem('Exponents','Free','Real',ipE_C,nSize*MxAt)
call GetMem('ContrCoef','Free','Real',ipC_C,nSize*MxAt)

!----------------------------------------------------------------------*
! A small test to see if max-limits are violated.                      *
!----------------------------------------------------------------------*
if ((Kmax > MxBB) .or. (Lmax > MxBB)) then
  write(u6,*)
  write(u6,*) 'ERROR! MxBB too small!'
  call Quit(_RC_INTERNAL_ERROR_)
end if
!----------------------------------------------------------------------*
! The multipoles and the Hamiltonian matrix are radically different    *
! between the QM-method alternatives, so once more an inquire.         *
!----------------------------------------------------------------------*
call NameRun('RUNFILE')
if (QmType(1:3) == 'SCF') then
  call ScfHandM(Cmo,nBas,iQ_Atoms,nOcc,natyp,nntyp,Occu)
  call mma_deallocate(Cmo)
  call mma_deallocate(Occu)
else if (QmType(1:4) == 'RASS') then
  call RassiHandM(nBas,iQ_Atoms,nOcc,natyp,nntyp)
end if

!----------------------------------------------------------------------*
! Here is the end.                                                     *
!----------------------------------------------------------------------*
return

end subroutine Qfread
