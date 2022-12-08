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
! Copyright (C) 2017,2020, Roland Lindh                                *
!***********************************************************************

subroutine Fix_FockOp(LuRd)
!***********************************************************************
!                                                                      *
!    Objective: To compute the fock operator for basis sets which do   *
!               not carry this information explicitly in the basis set *
!               file or are inline basis sets.                         *
!                                                                      *
!    The Fock operator is defined as F=\sum_{k,m} |k>e_{k,m}<m|.       *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
! Calling    : GetBS                                                   *
!                                                                      *
!     Author:  Roland Lindh (thanks to Ben Swerts)                     *
!                                                                      *
!     Modified for muons by R. Lindh January 2017                      *
!***********************************************************************

use Her_RW, only: nPrp
use Basis_Info, only: dbsc, nCnttp, Shells
use Sizes_of_Seward, only: S
use Gateway_Info, only: UnNorm, Do_FckInt, FNMC
use Isotopes, only: PTab
use Index_Functions, only: nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Six, Eight, Ten, Twelve
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuRd
#include "itmax.fh"
#include "Molcas.fh"
integer(kind=iwp) :: BasisTypes(4), i, iAng, iAngMax_Proj, iAtom, iB, iBF, iC, iCmp_a, iCmp_r, iCnttp, iComp, iFerm, iFrom, ijB, &
                     ijC, ijTri, Indx, iShll, iShll_a, iShll_Proj_r, iShll_r, iTo, jB, jBF, jShll, kEval, Last, List(0:iTabMx), &
                     List_Add(0:iTabMx), List_AE(0:iTabMx), lSTDINP, mCnttp, MemNA, MmKnEP, MmMltp, naa, nBF, nCntrc_a, &
                     nCntrc_Proj, nCntrc_r, nCntrc_t, nHer, nOrdOp, nPrim_a, nPrim_r, nRemove, nSAA, nSAR, nSBB, nSCC, nScr1, &
                     nScr2, nScr3, nSRR
real(kind=wp) :: A(4), C_ik, C_jk, Charge_Actual, Charge_Effective, Check, D, e, e12i, qTest, Test_Charge, Tmp, xFactor, xMass
logical(kind=iwp) :: Do_Cycle, lPP, Try_Again
character(len=256) :: Basis_lib, Fname
character(len=180) :: Ref(2)
character(len=80) :: Bsl_, BSLbl
real(kind=wp), allocatable :: C(:,:), E_R(:), EVal(:), EVec(:,:), FockOp_t(:,:), FPrim(:,:), Hm1(:,:), KnE(:), NAE(:), Ovr(:,:), &
                              Ovrlp(:), S12i(:,:), S_AA(:), S_AR(:), SAA(:), SAR(:), Scr1(:), Scr2(:), Scr3(:), Temp(:,:), &
                              Tmp1(:), Tmp2(:), Tmp3(:)
character(len=180), allocatable :: STDINP(:) ! CGGn
character(len=*), parameter :: DefNm = 'basis_library'
real(kind=wp), external :: DDot_
external :: KnEPrm, MltPrm, NAPrm

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
nPrint(114) = 99
nPrint(116) = 99
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
lPP = .false.
do i=1,nCnttp
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate a dummy center. This is fine since we will only do
! 1-center overlap integrals here.

A(:) = Zero
call mma_allocate(STDINP,mxAtom*2,label='STDINP')

nOrdOp = 2
iComp = 1

nPrp = max(4,S%nMltpl)

List(:) = 0
List_AE(:) = 0
BasisTypes(:) = 0
lSTDINP = 0

! Loop over all valence shell with a non-funtional FockOp

mCnttp = nCnttp   ! to be restored at the end
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
do iCnttp=1,mCnttp
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
  iFerm = 1
  if (dbsc(iCnttp)%fMass /= One) iFerm = 2

  if (dbsc(iCnttp)%FOp .and. (dbsc(iCnttp)%Charge == Zero)) then
    do iAng=0,dbsc(iCnttp)%nVal-1
      iShll_a = dbsc(iCnttp)%iVal+iAng
      Shells(iShll_a)%FockOp(:,:) = Zero
    end do
  end if

  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nFragType > 0) .or. dbsc(iCnttp)%FOp) then
    cycle
  end if

  ! Special treatment for muonic basis sets

  if (iFerm == 2) then

    iShll = S%Mx_Shll-1
    jShll = iShll

    ! The Fock operator will simply be the one-particle
    ! Hamiltonian (kinetic + nuclear-attraction operator)

    xFactor = One/dbsc(iCnttp)%fMass
    if (FNMC) then
      iAtom = dbsc(iCnttp)%AtmNr
      ! Get the atom mass in au (me=1)
      xMass = dbsc(iCnttp)%CntMass
      ! Substract the electron mass to get the nuclear mass.
      xMass = xMass-real(iAtom,kind=wp)
      xfactor = xfactor+One/xMass
    end if

    do iAng=0,dbsc(iCnttp)%nVal-1

      iShll_a = dbsc(iCnttp)%iVal+iAng
      nPrim_a = Shells(iShll_a)%nExp
      if (nPrim_a == 0) cycle
      nCntrc_a = Shells(iShll_a)%nBasis_C
      iCmp_a = (iAng+1)*(iAng+2)/2
      if (Shells(iShll_a)%Prjct) iCmp_a = 2*iAng+1
      naa = nTri_Elem1(iAng)*nTri_Elem1(iAng)
      nScr1 = max(nPrim_a,nPrim_a)*max(nCntrc_a,nCntrc_a)*naa
      nScr2 = max(nCntrc_a,nCntrc_a)**2*naa
      call mma_allocate(Scr1,nScr1,Label='Scr1')
      call mma_allocate(Scr2,nScr2,Label='Scr2')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute the kinetic integrals

      nOrdOp = 2
      nSAA = nCntrc_a**2*naa

      call KnEMmP(nHer,MmKnEP,iAng,iAng,nOrdOp)
      nScr3 = nPrim_a**2*MmKnEP
      call mma_allocate(Scr3,nScr3,Label='Scr1')

      call mma_Allocate(KnE,NSAA,Label='KnE')
      call One_Int(KnEPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,Scr1,nScr1,Scr2,nScr2,naa,KnE,nSAA,iShll_a,nPrim_a,Shells(iShll_a)%Exp, &
                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,iShll_a,nPrim_a,Shells(iShll_a)%Exp,nCntrc_a, &
                   Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
      call mma_deallocate(Scr3)
#     ifdef _DEBUGPRINT_
      call DScal_(nCntrc_a**2*iCmp_a**2,xFactor,KnE,1)
      call RecPrt('Kinetric Energy Integrals',' ',KnE,nCntrc_a**2,iCmp_a**2)
      call DScal_(nCntrc_a**2*iCmp_a**2,One/xFactor,KnE,1)
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute the nuclear-attraction integrals

      nOrdOp = 0
      A(4) = real(iCnttp,kind=wp) ! Dirty tweak
      nSBB = nCntrc_a**2*naa
      call mma_Allocate(NAE,nSBB,Label='NAE')

      call NAMem(nHer,MemNA,iAng,iAng,nOrdOp)
      nScr3 = nPrim_a**2*MemNA
      call mma_allocate(Scr3,nScr3,Label='Scr3')

      call One_Int(NAPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,Scr1,nScr1,Scr2,nScr2,naa,NAE,nSBB,iShll_a,nPrim_a,Shells(iShll_a)%Exp, &
                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,iShll_a,nPrim_a,Shells(iShll_a)%Exp,nCntrc_a, &
                   Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
      call mma_deallocate(Scr3)
#     ifdef _DEBUGPRINT_
      call RecPrt('Nuclear-attraction Integrals',' ',NAE,nCntrc_a**2,iCmp_a**2)
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Add together the kinetic and nuclear-attraction

      call DaXpY_(nCntrc_a**2*iCmp_a**2,xFactor,KnE,1,NAE,1)
      call mma_deallocate(KnE)

      ! Change to proper order (nCntrc_a * iCmp_a)

      call mma_allocate(Hm1,nCntrc_a**2,iCmp_a**2,Label='Hm1')
      call Reorder_GW(NAE,Hm1,nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute the overlap integrals

      nOrdOp = 0
      nSCC = nCntrc_a**2*naa
      call mma_allocate(Ovrlp,nSCC,Label='Ovrlp')

      call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
      nScr3 = nPrim_a**2*MmMltp
      call mma_allocate(Scr3,nScr3,Label='Scr3')

      call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,Scr1,nScr1,Scr2,nScr2,naa,Ovrlp,nSCC,iShll_a,nPrim_a,Shells(iShll_a)%Exp, &
                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,iShll_a,nPrim_a,Shells(iShll_a)%Exp,nCntrc_a, &
                   Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
      call mma_deallocate(Scr3)
#     ifdef _DEBUGPRINT_
      call RecPrt('Overlap Integrals',' ',Ovrlp,nCntrc_a**2,iCmp_a**2)
#     endif

      ! Change to proper order (nCntrc_a * iCmp_a)

      nBF = nCntrc_a*iCmp_a
      call mma_allocate(Ovr,nBF,nBF,Label='Ovr')
      call Reorder_GW(Ovrlp,Ovr,nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
      call mma_deallocate(Ovrlp)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Now we need to convert it to the Fock operator!
      !
      ! Solve F C = e S C
      !
      ! S^(-1/2) F S^(-1/2) S^(1/2) C = e S^(1/2) C
      ! Set F' = S^(-1/2) F S^(-1/2)
      !     C' = S^(1/2) C
      !
      ! Solve F' C' = e C' , generate C = S^-(1/2) C'

      call mma_Allocate(S12i,nBF,nBF,Label='S')
      S12i(:,:) = Zero

      ! 1) Compute the eigenvectors and eigenvalues of the overlap matrix

      call mma_allocate(EVal,nBF*(nBF+1)/2,Label='EVal')
      call mma_allocate(EVec,nBF,nBF,Label='EVec')
      call unitmat(EVec,nBF)
      do iBF=1,nBF
        do jBF=1,iBF
          ijTri = (iBF-1)*iBF/2+jBF
          EVal(ijTri) = Ovr(iBF,jBF)
        end do
      end do
      call mma_deallocate(Ovr)
      call NIDiag_new(EVal,EVec,nBF,nBF)

      ! 2) Construct S^(1/2) and S^(-1/2)

      do kEval=1,nBF
        e = EVal(kEval*(kEval+1)/2)
        e12i = One/sqrt(e)
        do iBF=1,nBF
          C_ik = EVec(iBF,kEVal)
          do jBF=1,nBF
            C_jk = EVec(jBF,kEVal)
            S12i(iBF,jBF) = S12i(iBF,jBF)+C_ik*e12i*C_jk
          end do
        end do
      end do

      ! 3) Form F' =  S^(-1/2) F S^(-1/2)

      call mma_allocate(FPrim,nBF,nBF,Label='FPrim')
      FPrim(:,:) = Zero
      call mma_allocate(Temp,nBF,nBF,Label='Temp')
      call DGEMM_('N','N',nBF,nBF,nBF,One,S12i,nBF,Hm1,nBF,Zero,Temp,nBF)
      call DGEMM_('N','N',nBF,nBF,nBF,One,Temp,nBF,S12i,nBF,Zero,FPrim,nBF)

      ! 4) Compute C' and the eigenvalues

      call unitmat(EVec,nBF)
      do iBF=1,nBF
        do jBF=1,iBF
          ijTri = (iBF-1)*iBF/2+jBF
          EVal(ijTri) = FPrim(iBF,jBF)
        end do
      end do
      call mma_deallocate(Temp)
      call mma_deallocate(FPrim)
      call NIDiag_new(EVal,EVec,nBF,nBF)

      ! 5) Form C = S^(-1/2) C'

      call mma_allocate(C,nBF,nBF,Label='C')
      C(:,:) = Zero
      call DGEMM_('N','N',nBF,nBF,nBF,One,S12i,nBF,EVec,nBF,Zero,C,nBF)
#     ifdef _DEBUGPRINT_
      call RecPrt('Cs for F',' ',C,nBF,nBF)
#     endif

      ! 6) Form the matrix representation of the Fock operator

      call mma_deallocate(Hm1)
      call mma_allocate(Hm1,nBF,nBF,Label='Hm1')
      Hm1(:,:) = Zero
      do kEval=1,nBF
        e = EVal(kEval*(kEval+1)/2)
        do iBF=1,nBF
          C_ik = C(iBF,kEVal)
          do jBF=1,nBF
            C_jk = C(jBF,kEVal)
            Hm1(iBF,jBF) = Hm1(iBF,jBF)+C_ik*e*C_jk
          end do
        end do
      end do
      call mma_deallocate(C)

      call Reorder_GW(Hm1,NAE,nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
      call mma_deallocate(Hm1)

      ! Make result isotropic and distribute

      do iB=1,nCntrc_a
        do jB=1,nCntrc_a
          ijB = (jB-1)*nCntrc_a+iB
          Tmp = Zero
          do iC=1,iCmp_a
            ijC = (iC-1)*iCmp_a+iC
            iFrom = (ijC-1)*nCntrc_a**2+ijB
            Tmp = Tmp+NAE(iFrom)
          end do
          Shells(iShll_a)%FockOp(iB,jB) = Tmp/real(iCmp_a,kind=wp)
        end do
      end do
      call mma_deallocate(NAE)
#     ifdef _DEBUGPRINT_
      call RecPrt('Actual Fock operator',' ',Shells(iShll_a)%FockOp,nCntrc_a,nCntrc_a)
#     endif
      call mma_deallocate(EVal)
      call mma_deallocate(EVec)
      call mma_deallocate(S12i)
      call mma_deallocate(Scr1)
      call mma_deallocate(Scr2)
    end do

    dbsc(iCnttp)%FOp = .true.
    cycle
  end if

  ! create a new basis set index (temporary)

  nCnttp = mCnttp+1
  if (nCnttp > Mxdbsc) then
    call WarningMessage(2,'Fix_FockOp: Increase Mxdbsc')
    call Abend()
  end if

  ! create the temporary basis set label for this element to
  ! read the corresponding ANO-RCC basis set.

  BSLbl = ' '
  BSLbl = PTab(dbsc(iCnttp)%AtmNr)

  if (BSLbl(1:1) == ' ') then
    BSLbl = BSLbl(2:2)//'.ANO-RCC.....'
  else
    BSLbl = BSLbl(1:2)//'.ANO-RCC.....'
  end if

  Last = len_trim(BSLbl)
  Indx = index(BSLbl,'/')

  Bsl_ = ' '
  if (Indx == 0) then
    call WhichMolcas(Basis_lib)
    if (Basis_lib(1:1) /= ' ') then
      ib = index(Basis_lib,' ')-1
      if (ib < 1) call SysAbendMsg('Read_ANO_RCC','Too long PATH to MOLCAS',' ')
      Fname = Basis_lib(1:ib)//'/basis_library'
    else
      Fname = DefNm
    end if
    Indx = Last+1
    Bsl_ = BSLbl
  else
    Fname = BSLbl(Indx+2:Last)
    if (Fname == ' ') then
      call WarningMessage(2,' No basis set library specified for BSLbl='//BSLbl//',Fname='//Fname)
      call Quit_OnUserError()
    end if
    Fname = adjustl(Fname)
    Bsl_ = BSLbl(1:Indx-1)
  end if

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*)
  write(u6,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLbl(1:Indx-1)
  write(u6,'(1X,A,A)') 'Basis set is read from library:',Fname
# endif

  ! Let's get the reference basis set (ANO-RCC).

  iShll = S%Mx_Shll-1
  jShll = iShll
  call GetBS(Fname,Bsl_,iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,lSTDINP,.false.,.true.,' ')

  if (.not. dbsc(nCnttp)%FOp) then
    write(u6,*) 'Fix_FockOp: reference basis doesn''t contain a proper Fock operator'
    cycle
  end if
  Shells(jShll+1)%Transf = .false.
  Shells(jShll+1)%Prjct = .false.
  Shells(jShll+2)%Transf = .false.
  Shells(jShll+2)%Prjct = .false.
  dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal+dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nPP
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Start processing shells of iCnttp and mCnttp. Loop only over
  ! the shells of iCnttp (mCnttp might be larger!)

  Try_Again = .true.
  List_Add(:) = 0

  Do_Cycle = .true.
  do while (Do_Cycle)
    Test_Charge = Zero
    do iAng=0,dbsc(iCnttp)%nVal-1

      ! Pointers to the actuall shell

      iShll_a = dbsc(iCnttp)%iVal+iAng
      nPrim_a = Shells(iShll_a)%nExp
      if (nPrim_a == 0) cycle
      nCntrc_a = Shells(iShll_a)%nBasis_C
      iCmp_a = (iAng+1)*(iAng+2)/2
      if (Shells(iShll_a)%Prjct) iCmp_a = 2*iAng+1

      ! Pointers to the reference shell

      iShll_r = dbsc(nCnttp)%iVal+iAng
      nPrim_r = Shells(iShll_r)%nExp
      if (nPrim_r == 0) then
        write(u6,*) 'GuessOrb option turned off!'
        dbsc(iCnttp)%FOp = .false.
        exit
      end if
      nCntrc_r = Shells(iShll_r)%nBasis_C
      iCmp_r = (iAng+1)*(iAng+2)/2
      if (Shells(iShll_r)%Prjct) iCmp_r = 2*iAng+1

      !                                                                *
      !*****************************************************************
      !                                                                *
      if (dbsc(iCnttp)%ECP) then
#       ifdef _DEBUGPRINT_
        if (lPP) then
          write(u6,*) 'Reference is ECP (Pseudo Potential)'
        else
          write(u6,*) 'Reference is ECP (Huzinaga type)'
        end if
        call RecPrt('Reference Exponents',' ',Shells(iShll_r)%Exp,1,nPrim_r)
        call RecPrt('Reference Coefficients',' ',Shells(iShll_r)%Cff_c(1,1,1),nPrim_r,nCntrc_r)
        call RecPrt('Reference Fock operator',' ',Shells(iShll_r)%FockOp,nCntrc_r,nCntrc_r)
#       endif
        call OrbType(dbsc(nCnttp)%AtmNr,List_AE,31)
        call ECP_Shells(dbsc(iCnttp)%AtmNr,List)
        if (lPP .or. (dbsc(iCnttp)%nM1 == 0)) then

          ! Pseud potential case

          nRemove = List_AE(iAng)-List(iAng)

        else

          ! Huzinaga type, remove according to the number of projected shells.

          iAngMax_Proj = dbsc(iCnttp)%nPrj
          if (iAng <= iAngMax_Proj) then
            iShll_Proj_r = dbsc(iCnttp)%iPrj+iAng
            nCntrc_Proj = Shells(iShll_Proj_r)%nBasis
            nRemove = nCntrc_Proj
          else
            nRemove = 0
          end if

          ! If too many try the default

          if (nRemove > nCntrc_r) then
            nRemove = List_AE(iAng)-List(iAng)
          end if

        end if ! lPP
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nRemove=',nRemove
        write(u6,*) 'List_Add(iAng)=',List_Add(iAng)
#       endif
        nRemove = nRemove-List_Add(iAng)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nRemove=',nRemove
#       endif
        Test_Charge = Test_Charge+real(2*(2*iAng+1)*nRemove,kind=wp)

        ! Update pointers in case of ECP

        ! Update the number of contracted functions of ref.
        nCntrc_t = nCntrc_r-nRemove
        ! Pick up relevant parts of the FockOp matrix of ref.
        call mma_allocate(FockOp_t,nCntrc_t,nCntrc_t)
        do i=1,nCntrc_t
          FockOp_t(:,i) = Shells(iShll_r)%FockOp(nRemove+1:nRemove+nCntrc_t,nRemove+i)
        end do
        nCntrc_r = nCntrc_t
      else
        nRemove = 0
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *

#     ifdef _DEBUGPRINT_
      call RecPrt('Actual Exponents',' ',Shells(iShll_a)%Exp,1,nPrim_a)
      call RecPrt('Actual Coefficients',' ',Shells(iShll_a)%Cff_c(1,1,1),nPrim_a,nCntrc_a)
      call RecPrt('Reference Exponents',' ',Shells(iShll_r)%Exp,1,nPrim_r)
      call RecPrt('Reference Coefficients',' ',Shells(iShll_r)%Cff_c(1,nRemove+1,1),nPrim_r,nCntrc_r)
      if (allocated(FockOp_t)) then
        call RecPrt('Reference Fock operator',' ',FockOp_t,nCntrc_r,nCntrc_r)
      else
        call RecPrt('Reference Fock operator',' ',Shells(iShll_r)%FockOp,nCntrc_r,nCntrc_r)
      end if
#     endif
      if (allocated(FockOp_t)) then
        Check = DDot_(nCntrc_r**2,FockOp_t,1,FockOp_t,1)
      else
        Check = DDot_(nCntrc_r**2,Shells(iShll_r)%FockOp,1,Shells(iShll_r)%FockOp,1)
      end if
      if ((Check == Zero) .or. (dbsc(iCnttp)%Charge == Zero)) then
        if (allocated(FockOp_t)) call mma_deallocate(FockOp_t)
        cycle
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      naa = nTri_Elem1(iAng)*nTri_Elem1(iAng)
      nScr1 = max(nPrim_a,nPrim_r)*max(nCntrc_a,nCntrc_r)*naa
      nScr2 = max(nCntrc_a,nCntrc_r)**2*naa
      call mma_allocate(Scr1,nScr1,Label='Scr1')
      call mma_allocate(Scr2,nScr2,Label='Scr2')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute S_AA

      nOrdOp = 0
      nSAA = nCntrc_a**2*naa
      call mma_allocate(SAA,nSAA,Label='SAA')

      call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
      nScr3 = nPrim_a**2*MmMltp
      call mma_allocate(Scr3,nScr3,Label='Scr3')

      call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,Scr1,nScr1,Scr2,nScr2,naa,SAA,nSAA,iShll_a,nPrim_a,Shells(iShll_a)%Exp, &
                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,iShll_a,nPrim_a,Shells(iShll_a)%Exp,nCntrc_a, &
                   Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
      call mma_deallocate(Scr3)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute S_AR

      nOrdOp = 0
      nSAR = nCntrc_a*nCntrc_r*naa
      call mma_allocate(SAR,nSAR,Label='SAR')

      call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
      nScr3 = nPrim_a*nPrim_r*MmMltp
      call mma_allocate(Scr3,nScr3,Label='Scr3')

      call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,Scr1,nScr1,SCr2,nScr2,naa,SAR,nSAR,iShll_a,nPrim_a,Shells(iShll_a)%Exp, &
                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,iShll_r,nPrim_r,Shells(iShll_r)%Exp,nCntrc_r, &
                   Shells(iShll_r)%Cff_c(1,1+nRemove,1),iCmp_a)
      call mma_deallocate(Scr3)

      nSRR = nCntrc_r**2*naa
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Reorder and compute the inverse of SAA

      call mma_allocate(S_AA,nSAA,Label='S_AA')
      call Reorder_GW(SAA,S_AA,nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
#     ifdef _DEBUGPRINT_
      call RecPrt('Reordered SAA',' ',S_AA,nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#     endif
      call MInv(S_AA,SAA,D,nCntrc_a*iCmp_a)
      call mma_deallocate(S_AA)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Det=',D
      call RecPrt('Inverse of SAA',' ',SAA,nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#     endif

      ! Reorder SAR
      call mma_allocate(S_AR,nSAR,Label='S_AR')
      call Reorder_GW(SAR,S_AR,nCntrc_a,nCntrc_r,iCmp_a,iCmp_r)
      call mma_deallocate(SAR)
#     ifdef _DEBUGPRINT_
      call RecPrt('Reordered SAR',' ',S_AR,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#     endif

      ! Expand and reorder the reference fock operator

      call mma_allocate(E_R,nSRR,Label='E_R')
      call mma_allocate(Tmp1,nSRR,Label='Tmp1')
      Tmp1(:) = Zero
      if (allocated(FockOp_t)) then
        do iB=1,nCntrc_r
          do jB=1,nCntrc_r
            ijB = (jB-1)*nCntrc_r+iB
            do iC=1,iCmp_r
              ijC = (iC-1)*iCmp_r+iC
              iTo = (ijC-1)*nCntrc_r**2+ijB
              Tmp1(iTo) = FockOp_t(iB,jB)
            end do
          end do
        end do
      else
        do iB=1,nCntrc_r
          do jB=1,nCntrc_r
            ijB = (jB-1)*nCntrc_r+iB
            Tmp = Shells(iShll_r)%FockOp(iB,jB)
            do iC=1,iCmp_r
              ijC = (iC-1)*iCmp_r+iC
              iTo = (ijC-1)*nCntrc_r**2+ijB
              Tmp1(iTo) = Tmp
            end do
          end do
        end do
      end if
#     ifdef _DEBUGPRINT_
      call RecPrt('Expanded ER',' ',Tmp1,nCntrc_r*nCntrc_r,iCmp_r*iCmp_r)
#     endif
      call Reorder_GW(Tmp1,E_R,nCntrc_r,nCntrc_r,iCmp_r,iCmp_r)
#     ifdef _DEBUGPRINT_
      call RecPrt('Reordered ER',' ',E_R,nCntrc_r*iCmp_r,nCntrc_r*iCmp_r)
#     endif
      call mma_deallocate(Tmp1)

      ! Form (SAA)-1 SAR

      call mma_allocate(Tmp1,nSAR,Label='Tmp1')
      call DGEMM_('N','N',nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_a*iCmp_a,One,SAA,nCntrc_a*iCmp_a,S_AR,nCntrc_a*iCmp_a,Zero,Tmp1, &
                  nCntrc_a*iCmp_a)
#     ifdef _DEBUGPRINT_
      call RecPrt('(SAA)^-1 SAR',' ',Tmp1,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#     endif
      call mma_deallocate(S_AR)

      ! Form (SAA)-1 SAR ER

      call mma_allocate(Tmp2,nSAR,Label='Tmp2')
      call DGEMM_('N','N',nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_r*iCmp_r,One,Tmp1,nCntrc_a*iCmp_a,E_R,nCntrc_r*iCmp_r,Zero,Tmp2, &
                  nCntrc_a*iCmp_a)
#     ifdef _DEBUGPRINT_
      call RecPrt('(SAA)^-1 SAR ER',' ',Tmp2,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#     endif
      call mma_deallocate(E_R)

      ! Form (SAA)-1 SAR ER (SAR)^T (SAA)-1

      call DGEMM_('N','T',nCntrc_a*iCmp_a,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,One,Tmp2,nCntrc_a*iCmp_a,Tmp1,nCntrc_a*iCmp_a,Zero,SAA, &
                  nCntrc_a*iCmp_a)
#     ifdef _DEBUGPRINT_
      call RecPrt('EA',' ',SAA,nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#     endif
      call mma_deallocate(Tmp2)
      call mma_deallocate(Tmp1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Now we just need to reorder and put it into place!

      call mma_allocate(Tmp3,nSAA,Label='Tmp3')
      call Reorder_GW(SAA,Tmp3,nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
      call mma_deallocate(SAA)
#     ifdef _DEBUGPRINT_
      call RecPrt('Reordered EA',' ',Tmp3,nCntrc_a*nCntrc_a,iCmp_a*iCmp_a)
#     endif

      do iB=1,nCntrc_a
        do jB=1,nCntrc_a
          ijB = iB+(jB-1)*nCntrc_a
          Tmp = Zero
          do iC=1,iCmp_a
            ijC = iC+(iC-1)*iCmp_a
            iFrom = ijB+(ijC-1)*nCntrc_a**2
            Tmp = Tmp+Tmp3(iFrom)
          end do
          Shells(iShll_a)%FockOp(iB,jB) = Tmp/real(iCmp_a,kind=wp)
        end do
      end do
      if (allocated(FockOp_t)) call mma_deallocate(FockOp_t)
#     ifdef _DEBUGPRINT_
      call RecPrt('Actual Fock operator',' ',Shells(iShll_a)%FockOp,nCntrc_a,nCntrc_a)
#     endif
      call mma_deallocate(Tmp3)
      call mma_deallocate(Scr1)
      call mma_deallocate(Scr2)
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do  ! iAng
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Deallocate the memory for the reference Fock operator

    do iShll_r=jShll+1,iShll
      if (allocated(Shells(iShll_r)%Exp)) call mma_deallocate(Shells(iShll_r)%Exp)
      Shells(iShll_r)%nExp = 0
      if (allocated(Shells(iShll_r)%FockOp)) call mma_deallocate(Shells(iShll_r)%FockOp)
      Shells(iShll_r)%nFockOp = 0
      if (allocated(Shells(iShll_r)%pCff)) call mma_deallocate(Shells(iShll_r)%pCff)
      if (allocated(Shells(iShll_r)%Cff_c)) call mma_deallocate(Shells(iShll_r)%Cff_c)
      if (allocated(Shells(iShll_r)%Cff_p)) call mma_deallocate(Shells(iShll_r)%Cff_p)
      Shells(iShll_r)%nExp = 0
      Shells(iShll_r)%nBasis = 0
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *

    Charge_Actual = real(dbsc(iCnttp)%AtmNr,kind=wp)
    Charge_Effective = dbsc(iCnttp)%Charge
    qTest = Test_Charge-(Charge_Actual-Charge_Effective)
    !write(u6,*) 'qtest, Test_Charge = ',qtest,Test_Charge
    !write(u6,*) 'Charge_Actual,Charge_Effective = ',Charge_Actual,Charge_Effective
    if ((qTest == Zero) .or. (dbsc(iCnttp)%Charge == Zero)) then
      dbsc(iCnttp)%FOp = .true.
    else if (Try_Again) then
      if (qTest == Two) then
        ! s
        List_Add(0) = 1
      else if (qTest == Six) then
        ! p
        List_Add(1) = 1
      else if (qTest == Ten) then
        ! d
        List_Add(2) = 1
      else if (qTest == Eight) then
        ! s,p
        List_Add(0) = 1
        List_Add(1) = 1
      else if (qTest == Twelve) then
        ! s,d
        List_Add(0) = 1
        List_Add(2) = 1
      else if (qTest == 16.0_wp) then
        ! p,d
        List_Add(1) = 1
        List_Add(2) = 1
      else if (qTest == 18.0_wp) then
        ! s,p,d
        List_Add(0) = 1
        List_Add(1) = 1
        List_Add(2) = 1
      else if (qTest == 26.0_wp) then
        ! 2s,2p,d
        List_Add(0) = 2
        List_Add(1) = 2
        List_Add(2) = 1
      end if
      Try_Again = .false.
      cycle
    else
      write(u6,*) 'GuessOrb option turned off!'
      dbsc(iCnttp)%FOp = .false.
    end if
    Do_Cycle = .false.
  end do
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
end do ! iCnttp
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Restore the correct nCnttp value

nCnttp = mCnttp

#ifdef _INSANE_DEBUGPRINT_
nPrint(113) = 5
nPrint(114) = 5
nPrint(116) = 5
nPrint(122) = 5
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Check if we can activate the computation of FckInt!

Do_FckInt = .true.
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nFragType > 0) .or. dbsc(iCnttp)%FOp) cycle

  Do_FckInt = Do_FckInt .and. dbsc(iCnttp)%FOp ! To be activated!

end do
call mma_deallocate(STDINP)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Fix_FockOp
