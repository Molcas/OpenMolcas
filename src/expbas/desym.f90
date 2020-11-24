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
! Copyright (C) 2009, Giovanni Li Manni                                *
!               2009, Francesco Aquilante                              *
!               2020, Oskar Weser                                      *
!***********************************************************************

module desymmetrize_mod
    use Basis_Info, only: nBas, nCnttp, dbsc, Shells, MolWgh
    use Center_Info, only: dc
    use definitions, only: wp
    use linalg_mod, only: abort_, verify_
    use Symmetry_Info, only: nIrrep, lIrrep
    use stdalloc, only: mma_allocate, mma_deallocate
    use sorting, only: swap, sort, argsort
    use sorting_funcs, only: leq_r, geq_r

    implicit none
    private
    public :: desym

contains

! symmetry-----> C1 INPORB
    Subroutine desym(iUHF)
!***********************************************************************
!                                                                      *
! Purpose: Convert INPORB file with symmetry to DeSymOrb               *
!          file with C1 symmetry.                                      *
!                                                                      *
!          G. Li Manni, F. Aquilante, University of Geneva, June 2009. *
!          Oskar Weser, MPI FKF Stuttgart,  November 2020.             *
!                                                                      *
! Example of input:                                                    *
!                                                                      *
!   &SEWARD                                                            *
!    coord                                                             *
!     file.xyz                                                         *
!    basis                                                             *
!     ano-rcc-mb                                                       *
!    oneonly                                                           *
!    nodk                                                              *
!                                                                      *
!   >>COPY $CurrDir/scf_sym.ScfOrb INPORB                              *
!                                                                      *
!   &EXPBAS                                                            *
!    NoEx                                                              *
!    Desy                                                              *
!   END                                                                *
!                                                                      *
!***********************************************************************
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "info_expbas.fh"
        integer, intent(in) :: iUHF
        !> Different kinds of orbitals
        !> f, i, 1, 2, 3, s, d
        integer, parameter :: n_orb_kinds = 7

        real(wp), parameter :: EorbThr = 50._wp
        real(wp) :: Coor(3, MxAtom), Znuc(MxAtom)
        character(len=LENIN) :: AtomLabel(MxAtom)
        character(len=512) :: FilesOrb
        character(len=LENIN8), allocatable :: label(:)
        character(len=8) :: MO_Label(maxbfn)
        integer :: ibas_lab(MxAtom), nOrb(8)
        !> This is the orbital kind for each orbital.
        integer, allocatable :: kind_per_orb(:)
        !> This is the number of orbitals for every kind.
        integer :: n_kinds(n_orb_kinds)
        integer, allocatable :: iOrdEor(:)
        character(len=LENIN8 + 1) :: gtolabel(maxbfn)
        character(len=50) :: VTitle
        character(len=128) :: SymOrbName
        logical :: Exist, y_cart, Found

        real(wp) :: check_CMO, check_energy, check_occupation
        real(wp) :: temporary
        integer :: nAtom, nData, nTest, nDeg, nTot, nTot2
        integer :: iCnttp, iAngMx_Valence
        integer :: nB, iS
        integer :: ipCent, ipCent2, ipCent3
        integer :: ipPhase, ipC2, ipV, ipC2_ab, ipV_ab
        integer :: mInd_ab
        integer :: mAdOcc, mAdEor
        integer :: mAdCMO, ipAux_ab, mAdIndt_ab, mAdOcc_ab, mAdEor_ab, mAdCMO_ab
        real(wp), allocatable :: new_orb_E(:), new_occ(:)
        real(wp), allocatable :: new_CMO(:, :)
        integer :: Lu_, iErr, notSymm
        integer :: iatom, iDeg, ishell
        integer :: iIrrep, iWfType, iWF
        integer :: iTempOrd

        integer :: iPL, jPL

        integer :: mdc, kk, i, j, ik, k, l, kk_Max, ii, iB, ipp, ic, iv
        integer :: ipc
        integer :: icontr, nBasisi, icntr

        integer :: iPrintLevel
        logical :: reduce_prt
        external :: reduce_prt, iPrintLevel

        integer, parameter :: nNumber = 61
        character(len=1), parameter :: &
            number(nNumber) = &
            ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', &
             'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', &
             'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', &
             'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', &
             'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', &
             'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', &
             'Z']
        integer, save :: iRc = 0
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        Check_CMO = 0._wp
        Check_Energy = 0._wp
        Check_Occupation = 0._wp
        y_cart = .false.

        Call f_Inquire('RUNFILE', Exist)
        call verify_(exist, 'Error finding RUNFILE')

        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !-----Read the characteristics of all different basis sets,
        !     provide each atom with a nuclear charge and establish
        !     a link between an atom and its basis set ---
        !
        !     NOTICE!!!
        !     This call will also fill info.fh and the dynamic storage in
        !     Work(ipInf)
        !
        Call Inter1(AtomLabel, iBas_Lab, Coor, Znuc, nAtom)
        Call Qpg_iArray('nOrb', Found, nData)
        If (Found) Then
            Call Get_iArray('nOrb', nOrb, nData)
        Else
            nOrb(:nIrrep) = nBas(:nIrrep)
        End If
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        iAngMx_Valence = maxval(dbsc%nVal - 1, mask=.not. (dbsc%Aux .or. dbsc%Frag))
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !     Compute memory requirements and allocate memory
        !
        nB = sum(nBas(0:nIrrep - 1))
        Call GetMem('ICENT', 'ALLO', 'INTE', ipCent, 8 * nB)
        Call GetMem('IPHASE', 'ALLO', 'INTE', ipPhase, 8 * nB)
        Call GetMem('nCENT', 'ALLO', 'INTE', ipCent2, nB)
        Call GetMem('ICENTER', 'ALLO', 'INTE', ipCent3, nB)
        Call GetMem('CMO2', 'ALLO', 'REAL', ipC2, nB**2)
        Call GetMem('VECTOR', 'ALLO', 'REAL', ipV, nB**2)
        call dcopy_(nB**2, [0._wp], 0, Work(ipV), 1)
        Call FZero(Work(ipC2), nB**2)
        If (iUHF == 1) Then
            Call GetMem('CMO2', 'ALLO', 'REAL', ipC2_ab, nB**2)
            Call GetMem('VECTOR', 'ALLO', 'REAL', ipV_ab, nB**2)
            call dcopy_(nB**2, [0._wp], 0, Work(ipV_ab), 1)
            Call FZero(Work(ipC2_ab), nB**2)
        Else
            ipC2_ab = ip_Dummy
            ipV_ab = ip_Dummy
        End If
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !     Read exponents and contraction coefficients of each unique basis.
        !     Write the present basis set (iCnttp) to the molden.input file for
        !     the appropriate atoms.
        !     Moreover, a list is constructed which contains a label for each
        !     GTO (gtolabel). This list follows the MOLDEN order of GTOs.
        !     Later this list will be used in the transformation of sabf (the
        !     symmetry adapted basis functions).
        !
        iatom = 0
        mdc = 0
        kk = 0
        !
        !            write(6,*)'nCnttp', nCnttp
        Do iCnttp = 1, nCnttp             ! loop over unique basis sets
            If (.not. (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) then
                !
                !         write(6,*)'dbsc(iCntt)%nCntr',dbsc(iCnttp)%nCntr
                Do iCntr = 1, dbsc(iCnttp)%nCntr! loop over symmetry unique centers
                    mdc = mdc + 1
                    nDeg = nIrrep / dc(mdc)%nStab
                    !            write(6,*)'nDeg', nDeg
                    Do iDeg = 1, nDeg             ! loop over centers
                        iAtom = iAtom + 1
                        !
                        If (dbsc(iCnttp)%nVal > 6) Then
                            Write (6, *) 'Desym: too high angular momentum!'
                            write (6, *) 'iCnttp and dbsc(iCnttp)%nVal= ' &
                                , iCnttp, dbsc(iCnttp)%nVal
                            Call Abend()
                        End If
                        !
                        Do l = 0, dbsc(iCnttp)%nVal - 1
                            ishell = dbsc(iCnttp)%iVal + l
                            nBasisi = Shells(iShell)%nBasis
                            !              write(6,*) 'nBasisi', Shells(iShell)%nBasis
                            If (nBasisi > nNumber) Then
                                Write (6, *) 'Desym: too many contracted functions!'
                                Write (6, *) 'nBasisi=', Shells(iShell)%nBasis
                                Call Abend()
                            End If

                            !             Iterate over each contracted GTO
                            !
                            If (l == 0) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'01s     '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If (l == 1) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02px    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02py    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'02pz    '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 2) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 3) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                            If ((l == 4) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            EndIf
                            If ((l == 5) .and. (.not. y_cart)) Then
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h04-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h03-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h02-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h01-  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h00   '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h01+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h02+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h03+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h04+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                                Do icontr = 1, nBasisi
                                    kk = kk + 1
                                    gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '// &
                                                   number(icontr)
                                    iWork(ipCent3 + kk - 1) = iAtom
                                End do
                            End If
                        End Do
                    End Do
                End Do
            end if
        End Do
        kk_Max = kk
        If (nB > kk_max) Then
            Write (6, *) 'nB > kk_max'
            Write (6, *) 'nB,kk_Max=', nB, kk_Max
            call ClsSew()
            return
        End If
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !
        nTot = 0
        nTot2 = 0
        Do iS = 0, nIrrep - 1
            nTot = nTot + nBas(iS)
            nTot2 = nTot2 + nBas(iS)**2
        End Do
        allocate (new_orb_E(0:nTot - 1))
        allocate (kind_per_orb(nTot))
        new_orb_E(:) = 0._wp
        Call GetMem('Occ', 'Allo', 'Real', mAdOcc, nTot)
        Call GetMem('Eor', 'Allo', 'Real', mAdEor, nTot)
        Call GetMem('CMO', 'Allo', 'Real', mAdCMO, nTot2)
        Call FZero(Work(mAdOcc), nTot)
        Call FZero(Work(mAdEor), nTot)
        Call FZero(Work(mAdCMO), nTot2)
        If (iUHF == 1) Then
            Call GetMem('Aux', 'ALLO', 'REAL', ipAux_ab, nTot)
            Call GetMem('INDT', 'Allo', 'Inte', mAdIndt_ab, ntot)
            Call GetMem('IndType', 'Allo', 'Inte', mInd_ab, 56)
            Call GetMem('Occ', 'Allo', 'Real', mAdOcc_ab, nTot)
            Call GetMem('Eor', 'Allo', 'Real', mAdEor_ab, nTot)
            Call GetMem('CMO', 'Allo', 'Real', mAdCMO_ab, nTot2)
            Call FZero(Work(ipAux_ab), nTot)
            Call IZero(iWork(mAdIndt_ab), nTot)
            Call IZero(iWork(mInd_ab), 56)
            Call FZero(Work(mAdOcc_ab), nTot)
            Call FZero(Work(mAdEor_ab), nTot)
            Call FZero(Work(mAdCMO_ab), nTot2)
        Else
            ipAux_ab = ip_Dummy
            mAdIndt_ab = ip_Dummy
            mInd_ab = ip_Dummy
            mAdOcc_ab = ip_Dummy
            mAdEor_ab = ip_Dummy
            mAdCMO_ab = ip_Dummy
        End If
        !
        !
        !---- Read HF CMOs from file
        !
        Lu_ = 75
        FilesOrb = EB_FileOrb
        If (len_trim(FilesOrb) == 0) FilesOrb = 'INPORB'
        if (DoExpbas) FilesOrb = 'EXPORB'
        Call RdVec_(trim(FilesOrb), Lu_, 'COEI', iUHF, nIrrep, nBas, nBas, &
                    Work(mAdCMO), Work(mAdCMO_ab), &
                    Work(mAdOcc), Work(mAdOcc_ab), &
                    Work(mAdEor), Work(mAdEor_ab), &
                    kind_per_orb, VTitle, 1, iErr, iWfType)
        if (iUHF == 1) then
            write (6, *) 'DESY keyword not implemented for DODS wf!'
            Call Abend()
        end if
        if (iErr /= 0) then
            iRc = 1
            return
        end if

        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !     Get the coeff. of sym adapted basis functions (ipC2)
        !
        Call Dens_IF_SCF(Work(ipC2), Work(mAdCMO), 'F')
        Call GetMem('CMO', 'Free', 'Real', mAdCMO, nTot2)
        If (iUHF == 1) Then
            Call Dens_IF_SCF(Work(ipC2_ab), Work(mAdCMO_ab), 'F')
            Call GetMem('CMO', 'Free', 'Real', mAdCMO_ab, nTot2)
        End If
        !        write(6,*)'write work(ipC2) in .out'
        !        write(6,*) (Work(k), k=ipC2,ipC2+nTot2-1)

        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !      Back 'transformation' of the symmetry adapted basis functions.
        !      Probably somewhat clumsy, but it seems to work.If someone
        !      knows a more elegant way to do it, please improve this part!
        !
        !      PART 1: Obtain symmetry information (soout), construct a label
        !              for each sabf, which will be used in part 2 to find the
        !              corresponding GTO in the MOLDEN list by comparing with
        !              gtolabel
        !
        !      nB       --- Total number of contracted basis functions
        !      ipcent2  --- degeneracy of a basis function
        !      ipCent   --- centres over which the basis function is
        !                   delocalized
        !      ipPhase  --- phase of the AO in the linear combination
        !
        Call mma_allocate(Label, MaxBfn + MaxBfn_Aux, label='Label')
        Call icopy(8 * nB, [0], 0, iWork(ipPhase), 1)
        Call icopy(8 * nB, [0], 0, iWork(ipCent), 1)
        Call SOout(label, iWork(ipCent), iWork(ipPhase))
        ipc = 0
        Do iContr = 1, nB
            iWork(ipCent2 + iContr - 1) = 0
            Do k = 1, 8
                If (iWork(ipCent + ipc) /= 0) &
                    iWork(ipcent2 + iContr - 1) = iWork(ipCent2 + iContr - 1) + 1
                ipc = ipc + 1
            End Do
        End Do
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !
        ! Part 2: -Take a MOLCAS symmetry functions (loop i)
        !         -Find the corresponding label in the MOLDEN list (loop j)
        !         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
        !          by the appropriate factor (ipPhase),and divide by the number of
        !          centres over which the sabf is delocalized (ipCent3).
        !         -The vectors are copied by rows!
        !
        !     loop over MOLCAS symmetry functions
        i = 0
        ik = 0
        Do iIrrep = 0, nIrrep - 1
            !
            Do iB = 1, nBas(iIrrep)
                i = i + 1
                !
                If (iB == 1) Then
                    ik = 1
                else
                    If (label(i - 1) == label(i)) Then
                        ik = ik + 1
                    else
                        ik = 1
                    End If
                End If
                !                                                                      *
                !***********************************************************************
                !                                                                      *
                Write (MO_Label(i), '(I5,A3)') iB, lirrep(iIrrep)
                !
                Do j = 1, nB
                    !
                    !
                    If (gtolabel(j) == label(i)//number(ik)) Then
                        Do k = 1, 8
                            ipc = (i - 1) * 8 + k - 1
                            ipp = ipc
                            If (iWork(ipCent + ipc) == iWork(ipcent3 + j - 1)) Then
                                Do ii = 1, nB
                                    ic = (ii - 1) * nB + (i - 1)
                                    iv = (ii - 1) * nB + (j - 1)
                                    If (MolWgh == 0) Then
                                        Work(ipV + iv) = Work(ipV + iv) &
                                                         + Work(ipC2 + ic) &
                                                         * DBLE(iWork(ipPhase + ipp)) &
                                                         / DBLE(iWork(ipcent2 + i - 1))
                                        If (iUHF == 1) Work(ipV_ab + iv) = Work(ipV_ab + iv) &
                                                                           + Work(ipC2_ab + ic) &
                                                                           * DBLE(iWork(ipPhase + ipp)) &
                                                                           / DBLE(iWork(ipcent2 + i - 1))
                                    Else
                                        Work(ipV + iv) = Work(ipV + iv) &
                                                         + Work(ipC2 + ic) &
                                                         * DBLE(iWork(ipPhase + ipp)) &
                                                         / Sqrt(DBLE(iWork(ipcent2 + i - 1)))
                                        If (iUHF == 1) Work(ipV_ab + iv) = Work(ipV_ab + iv) &
                                                                           + Work(ipC2_ab + ic) &
                                                                           * DBLE(iWork(ipPhase + ipp)) &
                                                                           / Sqrt(DBLE(iWork(ipcent2 + i - 1)))
                                    End If
                                End Do
                            End If
                        End Do
                    End If
                End Do
            End Do
        End Do
        Call mma_deallocate(Label)
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !      Dump vector in the molden.input file
        !
        ii = 0
        Do i = 0, nB - 1
        If (Work(mAdEOr + i) <= EorbThr) Then
            !          Write (MF,*) 'Sym= ',MO_Label(i+1)
            !          Write (MF,103) Work(mAdEOr+i)
            !          Write (MF,*) 'Spin= Alpha'
            !          Write (MF,104) Work(mAdOcc+i)
            If (Work(mAdEOr + i) < 0.0D0) Then
                Check_Energy = Check_Energy &
                               + Work(mAdEOr + i) * DBLE(i)
            End If
            Check_Occupation = Check_Occupation &
                               + Work(mAdOcc + i) * DBLE(i)
            Do j = 1, nB
                !            Write (MF,100) j,Work(ipV+ii+j-1)
                Check_CMO = Check_CMO &
                            + Work(ipV + ii + j - 1)**2
            End Do
        End If
        !
        If (iUHF == 1) Then
            If (Work(mAdEOr_ab + i) <= EorbThr) Then
                !             Write (MF,*) 'Sym= ',MO_Label(i+1)
                !             Write (MF,103) Work(mAdEOr_ab+i)
                !             Write (MF,*) 'Spin= Beta'
                !             Write (MF,104) Work(mAdOcc_ab+i)
                Check_Energy = Check_Energy &
                               + Work(mAdEOr_ab + i) * DBLE(i)
                Check_Occupation = Check_Occupation &
                                   + Work(mAdOcc_ab + i) * DBLE(i)
                Do j = 1, nB
                    !                Write (MF,100) j,Work(ipV_ab+ii+j-1)
                    Check_CMO = Check_CMO &
                                + Work(ipV_ab + ii + j - 1)**2
                End Do
            End If
        End If
        !
        ii = ii + nB
        End Do
        !**************************** START SORTING ****************************

        allocate (iOrdEor(0:nTot - 1))
        allocate (new_CMO(0:nTot - 1, 0:nTot - 1))
        allocate (new_occ(0:nTot - 1))

        n_kinds(:) = 0

        do i = lbound(kind_per_orb, 1), ubound(kind_per_orb, 1)
            n_kinds(kind_per_orb(i)) = n_kinds(kind_per_orb(i)) + 1
        end do

        iOrdEor(:) = argsort(Work(mAdEor:mAdEor + nTot - 1), leq_r) - 1

        new_orb_E(:) = Work(mAdEor + iOrdEor)

        do i = 0, nTot - 1
            do k = 0, nTot - 1
                new_CMO(k, i) = work(ipV + nTot * iOrdEor(i) + k)
            end do
        end do

        new_occ(:) = Work(mAdOcc + iOrdEor(:))

        iOrdEor(:) = argsort(new_occ, geq_r) - 1
        new_occ(:) = new_occ(iOrdEor)
        new_CMO(:, :) = new_CMO(:, iOrdEor)

        !> write

        notSymm = 1
        iWF = 9
        SymOrbName = 'DESORB'
        VTitle = 'Basis set desymmetrized orbital file DESORB'
        Call WrVec_(SymOrbName, iWF, 'COEI', iUHF, notSymm, [nTot], [nTot], &
                    new_CMO, Work(ipV_ab), &
                    new_occ, Work(mAdOcc_ab), &
                    new_orb_E, Work(ipAux_ab), &
                    n_kinds, VTitle, iWFtype)
        call Add_Info('desym CMO', new_CMO, 999, 8)
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        iPL = iPrintLevel(-1)
        jPL = iPL
        If (Reduce_Prt() .and. iPL < 3) jPL = 0
        If (jPL >= 3) Then
            Write (6, *)
            Write (6, *) ' INPORB_C1 file was generated!'
            Write (6, *)
        End If
        !                                                                      *
        !***********************************************************************
        !                                                                      *
        deallocate (new_occ)
        deallocate (new_CMO)
        Call GetMem('Eor', 'Free', 'Real', mAdEor, nTot)
        Call GetMem('Occ', 'Free', 'Real', mAdOcc, nTot)
        deallocate (new_orb_E)
        deallocate (kind_per_orb)
        If (iUHF == 1) Then
            Call GetMem('Eor', 'Free', 'Real', mAdEor_ab, nTot)
            Call GetMem('Occ', 'Free', 'Real', mAdOcc_ab, nTot)
            Call GetMem('Aux', 'FREE', 'REAL', ipAux_ab, nTot)
            Call GetMem('INDT', 'Free', 'Inte', mAdIndt_ab, nTot)
            Call GetMem('IndType', 'Free', 'Inte', mInd_ab, 56)
        End If
        Call GetMem('ICENT', 'FREE', 'INTE', ipCent, 8 * nB)
        Call GetMem('IPHASE', 'FREE', 'INTE', ipPhase, 8 * nB)
        Call GetMem('nCENT', 'FREE', 'INTE', ipCent2, nB)
        Call GetMem('ICENTER', 'FREE', 'INTE', ipCent3, nB)
        If (iUHF == 1) Then
            Call GetMem('CMO2', 'FREE', 'REAL', ipC2_ab, nB**2)
            Call GetMem('VECTOR', 'FREE', 'REAL', ipV_ab, nB**2)
        End If
        Call GetMem('CMO2', 'FREE', 'REAL', ipC2, nB**2)
        Call GetMem('VECTOR', 'FREE', 'REAL', ipV, nB**2)

        call ClsSew()

    end subroutine
end module
