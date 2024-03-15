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
! Copyright (C) 2012, Roland Lindh                                     *
!***********************************************************************

subroutine Mk_aCD_acCD_Shells(iCnttp,W2L)
!***********************************************************************
!                                                                      *
!    Objective: To generate aCD auxiliary basis sets on-the-fly.       *
!                                                                      *
! Called from: Mk_RICD_Shells                                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry - Angstrom              *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem, nTri_Elem1
use RI_procedures, only: Drv2El_Atomic_NoSym, Fix_Exponents
use SOAO_Info, only: iAOtSO, nSOInf, SOAO_Info_Free, SOAO_Info_Init
use Basis_Info, only: dbsc, Extend_Shells, Max_Shells, nCnttp, Shells
use Sizes_of_Seward, only: S
use RICD_Info, only: Do_acCD_Basis, Skip_High_AC, Thrshld_CD
use Integral_interfaces, only: Int_PostProcess, Integral_RICD
use define_af, only: iTabMx
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCnttp
logical(kind=iwp), intent(in) :: W2L
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAng, iAngMax, iAngMin, iAO, iBS, iCho_c, iCho_p, iCmp, iCntrc, iDum, iExp_k, iExp_l, ijS, ijS_req, ijSO, &
                     ijT, ik, ikl, il, Indx, iOff, iOff_Ak, iOff_Qk, ip_Exp, iRC, iSeed, iShell, iShll, iShll_, iSO, iSph, &
                     istatus, iTheta, iTheta_full, iVal, iZ, j, jAng, jAngMax, jAngMin, jCho_p, jCnttp, jkl, jp_Exp, jp_Exp_Max, &
                     jShll, jShll_, jTheta, jTheta_full, kAng, kC, Keep_All, Keep_Shell, kShll, lAng, lC, LinDep, lScr, lShll, &
                     Lu_A, Lu_B, Lu_lib, mData, mdc, mPrim, mSOInf, n, nBS, nCmp, nCmpi, nCmpj, nCnt, nCntrc, nCntrc_Max, &
                     nCnttp_Start, nExpi, nExpk, nExpl, nk, nl, nn, nPhi, nPhi_All, npi, npj, npk, npl, nPrim, nPrim_Max, nSO, &
                     nSO_p, nTest, nTheta, nTheta_All, nTheta_Full, nTInt_c, nTInt_p, nTri, NumCho_c, NumCho_p
real(kind=wp) :: Coeff_, Coeff_k, Coeff_kk, Coeff_kl, Coeff_l, Coeff_lk, Coeff_ll, Dummy(1), Exp_i, Exp_j, Fact, Thr_aCD, ThrAO, &
                 Thrs, Thrshld_CD_p
logical(kind=iwp) :: Diagonal, Found, Hit, In_Core, Keep_Basis
character(len=80) :: atom, author, Aux, basis, BSLbl, btype, CGTO, Label
integer(kind=iwp), allocatable :: Con(:), ConR(:,:), iD_c(:), iD_p(:), iList2_c(:,:), iList2_p(:,:), Indkl(:), Indkl_p(:), &
                                  LTP(:,:), Prm(:)
real(kind=wp), allocatable :: A(:), ADiag(:), C(:), Q(:), QTmp(:), Scr(:), Temp(:), TInt_c(:), TInt_p(:), Tmp(:), TP(:), tVp(:), &
                              tVt(:), tVtF(:), Vec(:), Wg(:), Z(:)
#ifdef _DEBUGPRINT_
real(kind=wp) :: Det
real(kind=wp), allocatable :: H(:), tVtInv(:), U(:)
#endif
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
ThrAO = Zero
mData = 4
nCnttp_Start = nCnttp
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Loop now over all unique valence basis sets and generate the
! corresponding aCD auxiliary basis sets. Note that there are two
! different types of aCD auxiliary basis sets, aCD and acCD.

nSO_p = 0
nTheta_All = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up the threshold for the CD procedure. Note that basis
! sets might have individual accuracy!

mdc = dbsc(iCnttp)%mdci
Thr_aCD = dbsc(iCnttp)%aCD_Thr*Thrshld_CD
!
nTest = dbsc(iCnttp)%nVal-1
!                                                                      *
!***********************************************************************
!                                                                      *
if (Skip_High_AC) then

  ! Pick up the angular index of the highest valence shell

  if (dbsc(iCnttp)%AtmNr <= 2) then
    iVal = 0
  else if (dbsc(iCnttp)%AtmNr <= 10) then
    iVal = 1
  else if (dbsc(iCnttp)%AtmNr <= 18) then
    iVal = 1
  else if (dbsc(iCnttp)%AtmNr <= 36) then
    iVal = 2
  else if (dbsc(iCnttp)%AtmNr <= 54) then
    iVal = 2
  else if (dbsc(iCnttp)%AtmNr <= 86) then
    iVal = 3
  else
    iVal = 3
  end if
  Keep_All = 2*nTest
  ! Find the number of polarization shells
  iZ = max(0,nTest-iVal)
  ! Reduce the product basis from excessive shells
  Keep_Shell = Keep_All-iZ
else
  Keep_Shell = iTabmx
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Define some parameters to facilitate the atomic calculation

iShell = dbsc(iCnttp)%nVal
S%nShlls = iShell
!                                                                      *
!***********************************************************************
!                                                                      *
! Use the name of the old valence basis

Label = dbsc(iCnttp)%Bsl_old

Hit = .true.
call Decode(Label,atom,1,Hit)
Hit = .true.
call Decode(Label,btype,2,Hit)
Hit = .true.
call Decode(Label,author,3,Hit)
Hit = .true.
call Decode(Label,basis,4,Hit)
Hit = .true.
call Decode(Label,CGTO,5,Hit)
Hit = .false.
call Decode(Label,Aux,6,Hit)
if (.not. Hit) Aux = ' '

n = index(Atom,' ')-1
Label = ' '
Label(1:n+1) = atom(1:n)//'.'
nn = n+1

n = index(btype,' ')-1
if (Do_acCD_Basis) then
  Label(nn+1:nn+n+23) = btype(1:n)//'....acCD-aux-basis.'
else
  Label(nn+1:nn+n+22) = btype(1:n)//'....aCD-aux-basis.'
end if

Indx = index(Label,' ')
BSLbl = ' '
BSLbl(1:Indx-1) = Label(1:Indx-1)

! Make a temporary setup of the SOAO_Info arrays for the
! atomic auxiliary basis set.
! Note that the auxiliary basis set will carry an angular value
! which at most is twice that of valence basis set.

mSOInf = 0

do iAng=0,2*nTest
  nCmp = nTri_Elem1(iAng)
  mSOInf = mSOInf+nCmp
end do
call SOAO_Info_Init(mSOInf,1)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! C O N T R A C T E D    S E C T I O N
!
! Run in contracted mode to generate the auxiliary basis for the
! aCD primitive product basis.

call Flip_Flop(.false.)
!                                                                      *
!***********************************************************************
!                                                                      *
! Define AOtSO

iAO = 0
iSO = 0
nSO = 0
do iAng=0,nTest
  iShll_ = dbsc(iCnttp)%iVal+iAng
  nCmp = nTri_Elem1(iAng)
  if (Shells(iShll_)%Prjct) nCmp = 2*iAng+1
  iSO = 0
  if (Shells(iShll_)%nBasis_C*Shells(iShll_)%nExp == 0) cycle
  do iCmp=1,nCmp
    iAO = iAO+1
    if (iAO > nSOInf) then
      write(u6,*) 'mk_acd_accd_shells: iAO>nSOInf (1)'
      write(u6,*) 'iAO=',iAO
      write(u6,*) 'nSOInf=',nSOInf
      call Abend()
    end if
    iAOtSO(iAO,0) = iSO+1
    iSO = iSO+Shells(iShll_)%nBasis
  end do
  nSO = nSO+iSO
end do

! Generate list

nPhi_All = nTri_Elem(nSO)
call mma_allocate(iList2_c,mData*2,nPhi_All,label='iList2_c')
call Mk_List2(iList2_c,nPhi_All,mData,nSO,iCnttp,nTest,0)
!                                                                      *
!***********************************************************************
!                                                                      *
! If the full product basis is used no need for decomposition!

if (Thr_aCD == Zero) then
  nTInt_c = nPhi_All
  call mma_allocate(iD_c,nTInt_c,label='iD_c')
  do i=1,nTInt_c
    iD_c(i) = i
  end do
  NumCho_c = nTInt_c
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate atomic two-electron integrals to decompose.

  ijS_req = 0
  Int_PostProcess => Integral_RICD
  call Drv2El_Atomic_NoSym(ThrAO,iCnttp,iCnttp,TInt_c,nTInt_c,In_Core,ADiag,Lu_A,ijS_req,Keep_Shell)
  Int_PostProcess => null()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Let us now decompose and retrieve the most important
  ! contracted products, indicies stored in iD_c

  call mma_allocate(iD_c,nTInt_c,label='iD_c')

  ! Temporary code for weights to be used in a MS-aCD/acCD
  ! scheme. Currently set to unit giving the convential
  ! all purpose aCD/acCD auxiliary basis sets.

  call mma_allocate(Wg,nTInt_c,label='Wg')
  Wg(:) = One

  if (In_Core) then
#   ifdef _DEBUGPRINT_
    call RecPrt('TInt_c',' ',TInt_c,nTInt_c,nTInt_c)
#   endif
    call mma_allocate(Vec,nTInt_c**2,label='Vec')

    call CD_InCore_p_w(TInt_c,nTInt_c,Wg,Vec,nTInt_c,iD_c,NumCho_c,Thr_aCD,iRC)

    if (iRC /= 0) then
      call WarningMessage(2,'Error in Mk_RICD_Shells')
      write(u6,*) 'Mk_aCD_Shells: CD_InCore_p(c) failed!'
      call Abend()
    end if
#   ifdef _DEBUGPRINT_
    call RecPrt('Vec',' ',Vec,nTInt_c,NumCho_c)
#   endif
    call mma_deallocate(TInt_c)
    call mma_deallocate(Vec)

  else    ! out-of-core part

    call mma_maxDBLE(lScr)
    lScr = min(lScr-2*nTInt_c,nTInt_c**2+3*nTInt_c)
    call mma_Allocate(Scr,lScr,label='Scr')

    iSeed = Lu_A+1
    Lu_B = IsFreeUnit(iSeed)
    call DaName_MF_WA(Lu_B,'AVEC1')

    call Get_Pivot_idx_w(ADiag,Wg,nTInt_c,NumCho_c,Lu_A,Lu_B,iD_c,Scr,lScr,Thr_aCD)

    call mma_deallocate(Scr)
    call mma_deallocate(ADiag)
    call DaEras(Lu_B)
    call DaEras(Lu_A)

  end if

  call mma_deallocate(Wg)

end if

if (NumCho_c < 1) then
  call WarningMessage(2,'Error in Mk_RICD_Shells')
  write(u6,*) 'Mk_aCD_Shells: NumCho_c < 1 is illegal!'
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Thr_aCD:',Thr_aCD
write(u6,*) 'NumCho_c:',NumCho_c
call iVcPrt('iD_c',' ',iD_c,NumCho_c)
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Define AOtSO for primitive integral calculations.

if (Do_acCD_Basis) then
  iAO = 0
  iSO = 0
  nSO_p = 0
  do iAng=0,nTest
    iShll_ = dbsc(iCnttp)%iVal+iAng
    nCmp = nTri_Elem1(iAng)
    if (Shells(iShll_)%Prjct) nCmp = 2*iAng+1
    iSO = 0
    do iCmp=1,nCmp
      iAO = iAO+1
      if (iAO > nSOInf) then
        write(u6,*) 'mk_acd_accd_shells: iAO>nSOInf (2)'
        call Abend()
      end if
      iAOtSO(iAO,0) = iSO+1
      iSO = iSO+Shells(iShll_)%nExp
    end do
    nSO_p = nSO_p+iSO
  end do
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
! Loop through angular products. Note that all the products
! of an atom require multiple basis sets since Seward is not
! structured to handle more than one shell of a specific
! angular at the time, i.e. a basis set contains only, for
! example, one d-shell. For an atomic basis spd we will have
! the p*p and d*s resulting in two independent shells with
! the same total angular momentum, d.

iShll = S%Mx_Shll-1

! Start now looping over the products and analys the result
! of the CD. Note the very peculiar loop structure over
! iBS, iAng, and jAng. This to reduce the number of
! created basis sets.

nBS = (nTest+2)/2
do iBS=0,nBS-1
  iAngMin = iBS
  iAngMax = nTest-iBS

  nCnttp = nCnttp+1
  Keep_Basis = .false.

  if (nCnttp > Mxdbsc) then
    call WarningMessage(2,'Error in Mk_RICD_Shells')
    write(u6,*) 'Mk_RICD_Shells: Increase Mxdbsc'
    call Abend()
  end if

  ! Some generic setting of information

  dbsc(nCnttp)%Bsl = Label
  dbsc(nCnttp)%Bsl_old = dbsc(nCnttp)%Bsl
  dbsc(nCnttp)%pChrg = dbsc(iCnttp)%pChrg
  dbsc(nCnttp)%Fixed = dbsc(iCnttp)%Fixed
  dbsc(nCnttp)%Parent_iCnttp = iCnttp
  dbsc(nCnttp)%iVal = iShll+1
  dbsc(nCnttp)%Aux = .true.
  dbsc(nCnttp)%aCD_Thr = dbsc(iCnttp)%aCD_Thr
  dbsc(nCnttp)%fMass = dbsc(iCnttp)%fMass
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Loop over shell pairs

  jShll = iShll
  do iAng=0,iAngMax
    jAngMax = min(iAng,iAngMin)
    iShll_ = dbsc(iCnttp)%iVal+iAng
    if (iAng == iAngMax) jAngMax = iAngMax
    if (iAng < iAngMin) jAngMax = 0
    jAngMin = iAngMin
    if (iAng <= iAngMin) jAngMin = 0
    do jAng=jAngMin,jAngMax
      jShll_ = dbsc(iCnttp)%iVal+jAng

      iShll = iShll+1
#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'iAng,jAng=',iAng,jAng
      write(u6,*) 'iAngMax=',iAngMax
#     endif
      if (iShll > size(Shells)) call Extend_Shells()
      Diagonal = iAng == jAng

      ! Examine if any contracted products of these two shells
      ! survived the CD procedure, or that it is an empty shell.

      Found = .false.
      kShll = -1
      lShll = -1
      do iCho_c=1,NumCho_c
        ijSO = iD_c(iCho_c)
        kAng = iList2_c(1,ijSO)
        lAng = iList2_c(2,ijSO)
        if ((iAng == kAng) .and. (jAng == lAng)) then
          kShll = iList2_c(7,ijSO)
          lShll = iList2_c(8,ijSO)
          Found = .true.
          exit
        end if
      end do

      ! Fake Found=.FALSE. for shells which should explicitly be empty.

      Found = Found .and. (jAng >= iAngMin) .and. (iAng >= iAngMin) .and. (iAng+jAng <= Keep_Shell)
      Keep_Basis = Found .or. Keep_Basis
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Found,kShll,lShll=',Found,kShll,lShll
#     endif
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! P R I M I T I V E   S E C T I O N
      !
      ! Run in uncontracted mode to produce a SLIM
      ! primitive product  basis.

      if (Do_acCD_Basis .and. Found) then

        call Flip_Flop(.true.)

        ! Generate list

        npi = Shells(iShll_)%nExp
        nCmpi = nTri_Elem1(iAng)
        if (Shells(iShll_)%Prjct) nCmpi = 2*iAng+1
        npj = Shells(jShll_)%nExp
        nCmpj = nTri_Elem1(jAng)
        if (Shells(jShll_)%Prjct) nCmpj = 2*jAng+1
        if (iAng == jAng) then
          nTheta_All = nTri_Elem(npi*nCmpi)
        else
          nTheta_All = npi*nCmpi*npj*nCmpj
        end if

        call mma_allocate(iList2_p,2*mData,nTheta_all,label='iList2_p')

        ijS_Req = iTri(iAng+1,jAng+1)

        call Mk_List2(iList2_p,nTheta_All,mData,nSO_p,iCnttp,nTest,ijS_Req)
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Generate atomic two-electron integrals

        Int_PostProcess => Integral_RICD
        call Drv2El_Atomic_NoSym(ThrAO,iCnttp,iCnttp,TInt_p,nTInt_p,In_Core,ADiag,Lu_A,ijS_Req,Keep_Shell)
        Int_PostProcess => null()

        if (.not. In_Core) then
          call WarningMessage(2,'Error in Mk_RICD_Shells')
          write(u6,*) 'Out-of-core acCD not implemented!'
          call Abend()
        end if
#       ifdef _DEBUGPRINT_
        call RecPrt('TInt_p','(5G20.11)',TInt_p,nTInt_p,nTInt_p)
#       endif
        call Flip_Flop(.false.)

      end if
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! Now mimic the procedure of GetBS!
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Working on the CONTRACTED functions.
      !
      ! This section is identical for acCD and aCD auxiliary basis sets!
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (Found) then

        lAng = iAng+jAng

        ! Now figure out how many and which!

        nk = Shells(kShll)%nBasis_C
        nl = Shells(lShll)%nBasis_C
        if (Diagonal) then
          nCntrc_Max = nTri_Elem(nk)
        else
          nCntrc_Max = nk*nl
        end if
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nCntrc_Max=',nCntrc_Max
#       endif
        call mma_allocate(Con,nCntrc_Max,label='Con')
        call mma_allocate(ConR,2,nCntrc_Max,label='ConR')
        Con(:) = 0
        ConR(:,:) = 0
        nCntrc = 0
        do iCho_c=1,NumCho_c
          ijSO = iD_c(iCho_c)
          kAng = iList2_c(1,ijSO)
          lAng = iList2_c(2,ijSO)
          if ((kAng == iAng) .and. (lAng == jAng)) then

            ! Pick up the radial index!

            ik = iList2_c(5,ijSO)
            il = iList2_c(6,ijSO)

            if (Diagonal) then
              ikl = iTri(ik,il)
            else
              ikl = (il-1)*nk+ik
            end if

            ! Note that this migh be done several time since several
            ! angular pairs might have the same radial function!

            if (Con(ikl) == 0) then
              nCntrc = nCntrc+1
              Con(ikl) = 1
              ConR(1,nCntrc) = ik
#             ifdef _DEBUGPRINT_
              write(u6,*) 'iCho_c,  ijSO=',iCho_c+1,ijSO
#             endif
              ConR(2,nCntrc) = il
            end if
          end if
        end do    !  iCho_c
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nCntrc=',nCntrc
        call iVcPrt('Con',' ',Con,nCntrc_Max)
        call iVcPrt('ConR',' ',ConR,2*nCntrc)
        write(u6,*)
        write(u6,*) 'ConR'
        write(u6,'(30I3)') (ConR(1,i),i=1,nCntrc)
        write(u6,'(30I3)') (ConR(2,i),i=1,nCntrc)
#       endif

      else

        ! Let us put in an empty shell!

        nk = 0
        nl = 0
        nCntrc = 0
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Work on the PRIMITIVE products!
      !
      ! Here the work is trivial in case of the aCD basis
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (Found) then
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Produce the SLIM primitive products
        !                                                              *
        !***************************************************************
        !                                                              *
        if (Do_acCD_Basis) then

          ! Now figure out how many and which!

          npk = Shells(kShll)%nExp
          npl = Shells(lShll)%nExp
          if (Diagonal) then
            nPrim_Max = nTri_Elem(npk)
          else
            nPrim_Max = npk*npl
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) 'nPrim_Max:',nPrim_Max
#         endif
          call mma_allocate(Prm,nPrim_Max,label='Prm')
          Prm(:) = 0

          ! Pick up the diagonal elements from TInt_p
          ! corresponding to this shell pair. We sum over
          ! the angular parts identical to those of the contracted.

          ! First make a list from the contracted which angular products to include.

          call mma_allocate(TP,nPrim_Max**2,label='TP')
          call mma_allocate(LTP,2,nPrim_Max,label='LTP')
          call Mk_TInt_P(TInt_p,nTheta_All,TP,nPrim_Max,iList2_p,nTheta_All,2*mData,iAng,jAng,npk,LTP)

#         ifdef _DEBUGPRINT_
          call RecPrt('TIntP','(5G20.10)',TP,nPrim_Max,nPrim_Max)
          call iVcPrt('List_TP',' ',LTP,2*nPrim_Max)
#         endif
          ! Let us now decompose and retrieve the most
          ! important primitive products, indicies stored in iD_p

          call mma_allocate(iD_p,nPrim_Max,label='iD_p')
          call mma_allocate(Vec,nPrim_Max**2,label='Vec')

          Thrshld_CD_p = Thr_aCD*0.2_wp
          do
            call CD_InCore_p(TP,nPrim_Max,Vec,nPrim_Max,iD_p,NumCho_p,Thrshld_CD_p,iRC)
            if (NumCho_p < 1) then
              call WarningMessage(2,'Error in Mk_RICD_Shells')
              write(u6,*) 'Mk_aCD_Shells: NumCho_p < 1 is illegal!'
              write(u6,*) 'iAng,jAng=',iAng,jAng
              write(u6,*) 'nPrim_Max=',nPrim_Max
              write(u6,*) 'NumCho_p=',NumCho_p
              write(u6,*) 'iRC=',iRC
              call Abend()
            end if

#           ifdef _DEBUGPRINT_
            write(u6,*) 'Thrshld_CD_p:',Thrshld_CD_p
            write(u6,*) 'NumCho_p    :',NumCho_p
            call iVcPrt('iD_p',' ',iD_p,NumCho_p)
            call RecPrt('Vec',' ',Vec,nPrim_Max,NumCho_p)
#           endif
            if (NumCho_p >= nCntrc) exit
            write(u6,*) 'W a r n i n g!'
            write(u6,*) 'Fewer primitive functions than contracted functions!'
            write(u6,*) 'NumCho_p=',NumCho_p
            write(u6,*) '  nCntrc=',nCntrc
            Thrshld_CD_p = Thrshld_CD_p*Half
            if (Thrshld_CD_p <= 1.0e-14_wp) then
              call WarningMessage(2,'Error in Mk_RICD_Shells')
              write(u6,*) 'Thrshld_CD_p is too low!'
              write(u6,*) 'iAng, jAng:',iAng,jAng
              call Abend()
            end if
            call Mk_TInt_P(TInt_p,nTheta_All,TP,nPrim_Max,iList2_p,nTheta_All,2*mData,iAng,jAng,npk,LTP)
          end do
          call mma_deallocate(TP)
          call mma_deallocate(Vec)

          do iCho_p=1,NumCho_p
            ijSO = iD_p(iCho_p)
            Prm(ijSO) = 1
          end do
          nPrim = NumCho_p
          call mma_allocate(Shells(iShll)%Exp,nPrim,Label='ExpacCD')
          Shells(iShll)%nExp = nPrim

#         ifdef _DEBUGPRINT_
          write(u6,*) 'nPrim=',nPrim
          call iVcPrt('Prm',' ',Prm,nPrim_Max)
#         endif
          call mma_allocate(Indkl_p,nPrim_Max,label='Indkl_p')
          call Mk_Indkl(Prm,Indkl_p,nPrim_Max)

          ! Observe that the exponents are ordered according
          ! to their importance as given by the CD.

          do iCho_p=1,NumCho_p
            iTheta = iD_p(iCho_p)
            ik = LTP(1,iTheta)
            il = LTP(2,iTheta)
            Exp_i = Shells(kShll)%Exp(ik)
            Exp_j = Shells(lShll)%Exp(il)
            Shells(iShll)%Exp(iCho_p) = Exp_i+Exp_j
          end do
#         ifdef _DEBUGPRINT_
          call RecPrt('SLIM Exponents',' ',Shells(iShll)%Exp,1,nPrim)
#         endif
          !                                                            *
          !*************************************************************
          !                                                            *
        else ! Do_aCD_Basis
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Put in the aCD set of exponents, i.e. all unique sums.

          nExpk = Shells(kShll)%nExp
          nExpl = Shells(lShll)%nExp
          if (Diagonal) then
            nPrim = nTri_Elem(nExpk)
          else
            nPrim = nExpk*nExpl
          end if
          call mma_allocate(Shells(iShll)%Exp,nPrim,Label='ExpaCD')
          Shells(iShll)%nExp = nPrim

          iOff = 0
          do ip_Exp=1,nExpk
            jp_Exp_Max = nExpl
            if (Diagonal) jp_Exp_Max = ip_Exp
            do jp_Exp=1,jp_Exp_Max
              iOff = iOff+1
              Shells(iShll)%Exp(iOff) = Shells(kShll)%Exp(ip_Exp)+Shells(lShll)%Exp(jp_Exp)
            end do
          end do

          if (iOff /= nPrim) then
            call WarningMessage(2,'Error in Mk_RICD_Shells')
            write(u6,*) 'Mk_aCD_Shell: iOff /= iEnd'
            call Abend()
          end if

#         ifdef _DEBUGPRINT_
          if (Diagonal) then
            call TriPrt('aCD Exponents',' ',Shells(iShll)%Exp,nExpk)
          else
            call RecPrt('aCD Exponents',' ',Shells(iShll)%Exp,nExpk,nExpl)
          end if
#         endif
        end if

      else
        !                                                              *
        !***************************************************************
        !                                                              *
        ! An empty shell

        nPrim = 0
        Shells(iShll)%nExp = nPrim
        !                                                              *
        !***************************************************************
        !                                                              *
      end if ! Found
      !
      !*****************************************************************
      !
      lAng = iAng+jAng
      S%iAngMx = max(S%iAngMx,lAng)
      S%MaxPrm(lAng) = max(S%MaxPrm(lAng),nPrim)

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'iShll=',iShll
      write(u6,*) 'nPrim,nCntrc=',nPrim,nCntrc
      write(u6,*) 'lAng=',lAng
      write(u6,*) 'S%MaxPrm(lAng)=',S%MaxPrm(lAng)
#     endif

      Shells(iShll)%nBasis_c = nCntrc
      !                                                                *
      !*****************************************************************
      !                                                                *
      call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
      call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
      Shells(iShll)%nBasis = nCntrc
      call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! C O N T R A C T I O N    C O E F F I C I E N T S

      if (Found) then
        !                                                              *
        !***************************************************************
        !                                                              *
        ! For SLIM basis sets
        !                                                              *
        !***************************************************************
        !                                                              *
        if (Do_acCD_Basis) then

          ! Alright this is a bit more elaborate than for
          ! the aCD basis set. Surprise!

          ! Some care has to be taken here. There might be
          ! different angular products, for example, px*px
          ! and px*py, which carry the same radial part but
          ! have different angular part! To overcome this
          ! possible source of redundancy we use the sum of
          ! such terms in the fitting procedure!

          nTheta = nPrim
          nExpk = Shells(kShll)%nExp
          nExpl = Shells(lShll)%nExp
          if (iAng == jAng) then
            nTheta_Full = nTri_Elem(nExpk)
          else
            nTheta_Full = nExpk*nExpl
          end if
          nPhi = nCntrc

          ! Generate the (theta'|V|theta') matrix in the
          ! SLIM primitive product basis.

          call mma_allocate(tVt,nTheta**2,label='tVt')
          call Mk_tVt(TInt_p,nTInt_p,tVt,nTheta,iList2_p,2*mData,Prm,nPrim_Max,iAng,jAng,nExpk,Indkl_p,nPrim_Max)

#         ifdef _DEBUGPRINT_
          write(u6,*)
          write(u6,*) 'tVt(Diag)'
          write(u6,*) (tVt(i),i=1,nTheta**2,nTheta+1)
          call RecPrt('tVt',' ',tVt,nTheta,nTheta)
          call iVcPrt('iD_p',' ',iD_p,NumCho_p)
#         endif

          ! Generate (theta'|V|theta')^{-1}

          ! Let's do a Cholesky decomposition with pivoting
          ! according to the previous CD.

          nTri = nTri_Elem(nTheta)
          call mma_allocate(Q,nTri,label='Q')
          call mma_allocate(A,nTri,label='A')
          call mma_allocate(Z,nTheta,label='Z')
          do iCho_p=1,NumCho_p
            iTheta_full = iD_p(iCho_p)
            iTheta = Indkl_p(iTheta_full)
            do jCho_p=1,iCho_p
              jTheta_full = iD_p(jCho_p)
              jTheta = Indkl_p(jTheta_full)
              ijT = iTri(iCho_p,jCho_p)
              ijS = (jTheta-1)*nTheta+iTheta
              A(ijT) = tVt(ijS)
            end do
          end do
#         ifdef _DEBUGPRINT_
          call TriPrt('A',' ',A,nTheta)

          call mma_allocate(H,nTri,label='H')
          call mma_allocate(U,nTri,label='U')
          H(:) = A
          call unitmat(U,nTheta)
          call Jacob(H,U,nTheta,nTheta)
          call TriPrt('H','(10G20.10)',H,nTheta)
          call RecPrt('U',' ',U,nTheta,nTheta)
          call mma_deallocate(H)
          call mma_deallocate(U)
#         endif

#         ifdef _DEBUGPRINT_
          call mma_allocate(tVtInv,nTheta**2,label='tVtInv')
          Det = Zero
          call MInv(tVt,tVtInv,Det,nTheta)
          write(u6,*) 'Det=',Det
#         endif
          call mma_deallocate(tVt)

          do iTheta=1,nTheta
            iOff_Ak = nTri_Elem(iTheta-1)+1
            iOff_Qk = nTri_Elem(iTheta-1)+1
            Thrs = Thrshld_CD_p
            !Thrs = 1.0e-12_wp
            call Inv_Cho_Factor(A(iOff_Ak),iTheta,A,Q,iTheta,iDum,iDum,Dummy,0,Z,Dummy,Thrs,Q(iOff_Qk),LinDep)
            if (LinDep == 1) then
              call WarningMessage(2,'Error in Mk_RICD_Shells')
              write(u6,*) 'Mk_aCD_Shells: linear dependence found in tVt!'
              write(u6,*) 'Found for vector:',iTheta
              call Abend()
            end if
          end do
          call mma_deallocate(Z)
          call mma_deallocate(A)
#         ifdef _DEBUGPRINT_
          call TriPrt('Q','(9G10.3)',Q,nTheta)
#         endif

          ! Generate the (theta'|V|theta) matrix in the mixed product basis.

          call mma_allocate(tVp,nTheta*nPhi,label='tVp')
          call mma_allocate(tVtF,nTheta*nTheta_Full,label='tVtF')
          call Mk_tVtF(TInt_p,nTInt_p,tVtF,nTheta,iList2_p,2*mData,Prm,nPrim_Max,iAng,jAng,nExpk,Indkl_p,nPrim_Max,nTheta_Full)
#         ifdef _DEBUGPRINT_
          call RecPrt('tVtF',' ',tVtF,nTheta,nTheta_Full)
#         endif

          ! Pick up the contraction coefficients of the aCD
          ! basis set. Be careful what this means in the
          ! case that the shells are identical!

          call mma_allocate(Indkl,nCntrc_Max,label='Indkl')
          call Mk_Indkl(Con,Indkl,nCntrc_Max)
          call mma_allocate(C,nTheta_Full*nPhi,label='C')
          call Mk_Coeffs(Shells(kShll)%Cff_c(1,1,1),nExpk,Shells(kShll)%nBasis_C,Shells(lShll)%Cff_c(1,1,1),nExpl, &
                         Shells(lShll)%nBasis_C,C,nTheta_Full,nPhi,iD_c,NumCho_c,iList2_c,2*mData,nPhi_All,Indkl,nCntrc_Max, &
                         Shells(kShll)%nBasis_C,iAng,jAng,Shells(kShll)%Cff_p(1,1,1),Shells(lShll)%Cff_p(1,1,1))
          call mma_deallocate(Indkl)
#         ifdef _DEBUGPRINT_
          call RecPrt('C',' ',C,nTheta_Full,nPhi)
#         endif

          ! Generate the (theta'|V|phi') matrix.

          call DGEMM_('N','N',nTheta,nPhi,nTheta_Full,One,tVtF,nTheta,C,nTheta_Full,Zero,tVp,nTheta)
          call mma_deallocate(tVtF)
#         ifdef _DEBUGPRINT_
          call RecPrt('tVp',' ',tVp,nTheta,nPhi)
#         endif
          call mma_deallocate(C)

          ! Generate the contraction coefficients of the
          ! SLIM contracted product basis in terms of the
          ! SLIM primitive product basis as
          ! Sum(nu') (mu'|V|nu')^-1  (nu'|V|i')=C_{mu',i'}

          ! To simplify life I will put the Q matrix into square storage.

          call mma_Allocate(Temp,nTheta**2,label='Temp')
          Temp(:) = Zero
          do iTheta=1,nTheta
            do jTheta=1,iTheta
              ijT = iTri(iTheta,jTheta)
              ijS = (iTheta-1)*nTheta+jTheta
              Temp(ijS) = Q(ijT)
            end do
          end do
          call mma_deallocate(Q)
#         ifdef _DEBUGPRINT_
          call RecPrt('Q',' ',Temp,nTheta,nTheta)
#         endif

          ! Resort the external index back to original
          ! order. The column index is external.

          call mma_allocate(QTmp,nTheta**2,label='QTmp')
          do iCho_p=1,NumCho_p
            iTheta_Full = iD_p(iCho_p)
            iTheta = Indkl_p(iTheta_Full)
            call dcopy_(nTheta,Temp(iCho_p),nTheta,QTmp(iTheta),nTheta)
          end do
          call mma_deallocate(Temp)
#         ifdef _DEBUGPRINT_
          call RecPrt('Q',' ',QTmp,nTheta,nTheta)
          call RecPrt('tVp',' ',tVp,nTheta,nPhi)
#         endif
          ! Q(T) tVp
          call mma_allocate(Scr,nTheta*nPhi,label='Scr')
          Scr(:) = Zero
          call DGEMM_('T','N',nTheta,nPhi,nTheta,One,QTmp,nTheta,tVp,nTheta,Zero,Scr,nTheta)
          ! QQ(T) tVp
          call DGEMM_('N','N',nTheta,nPhi,nTheta,One,QTmp,nTheta,Scr,nTheta,Zero,Shells(iShll)%Cff_c(1,1,1),nTheta)
#         ifdef _DEBUGPRINT_
          call RecPrt('SLIM coeffcients',' ',Shells(iShll)%Cff_c(1,1,1),nTheta,nPhi)
          Scr(:) = Zero
          call DGEMM_('N','N',nTheta,nPhi,nTheta,One,tVtInv,nTheta,tVp,nTheta,Zero,Scr,nTheta)
          call RecPrt('SLIM coeffcients2',' ',Scr,nTheta,nPhi)
          call mma_deallocate(tVtInv)
#         endif
          call mma_deallocate(tVp)

          ! Now reorder the coefficients to the CD order of the exponents.

          call mma_allocate(Tmp,nTheta*nPhi,label='Tmp')
          call dcopy_(nTheta*nPhi,Shells(iShll)%Cff_c(:,:,1),1,Tmp,1)
          do iCho_p=1,NumCho_p
            iTheta_full = iD_p(iCho_p)
            iTheta = Indkl_p(iTheta_full)
            call dcopy_(nPhi,Tmp(iTheta),nTheta,Shells(iShll)%Cff_c(iCho_p,1,1),nTheta)
          end do
          call mma_deallocate(Tmp)

          ! Modify from coefficients for normalized
          ! Gaussians to unnormalized Gaussians.

          do iCho_p=1,NumCho_p
            iTheta = iD_p(iCho_p)
            ik = LTP(1,iTheta)
            il = LTP(2,iTheta)
            Fact = Shells(kShll)%Cff_p(ik,ik,1)*Shells(lShll)%Cff_p(il,il,1)
            Shells(iShll)%Cff_c(iCho_p,:,1) = Fact*Shells(iShll)%Cff_c(iCho_p,:,1)
          end do
          call mma_deallocate(LTP)

          call mma_deallocate(iD_p)
          call mma_deallocate(Indkl_p)
          call mma_deallocate(Scr)
          call mma_deallocate(QTmp)
#         ifdef _DEBUGPRINT_
          call RecPrt('SLIM coeffcients',' ',Shells(iShll)%Cff_c(:,:,1),nTheta,nPhi)
#         endif
          !                                                            *
          !*************************************************************
          !                                                            *
        else ! Do_aCD_Basis
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Put in the selected set of coeffients. Note
          ! again that the order should be that according
          ! to the CD in order to prepivot, since the CD
          ! itself is implemented without pivoting.

          do iCntrc=1,nCntrc
            kC = ConR(1,iCntrc)
            lC = ConR(2,iCntrc)
#           ifdef _DEBUGPRINT_
            write(u6,*) 'kC,lC=',kC,lC
#           endif
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Form the unnormalized coefficients!

            jkl = 0
            if (Diagonal) then
              do iExp_k=1,Shells(kShll)%nExp

                Coeff_kk = Shells(kShll)%Cff_c(iExp_k,kC,1)
                Coeff_lk = Shells(lShll)%Cff_c(iExp_k,lC,1)

                do iExp_l=1,iExp_k

                  Coeff_ll = Shells(lShll)%Cff_c(iExp_l,lC,1)
                  Coeff_kl = Shells(kShll)%Cff_c(iExp_l,kC,1)
                  Coeff_ = Coeff_ll*Coeff_kk+Coeff_kl*Coeff_lk
                  if (iExp_k == iExp_l) then
                    Coeff_ = Coeff_*Half
                  end if
                  jkl = jkl+1

                  Shells(iShll)%Cff_c(jkl,iCntrc,1) = Coeff_

                end do
              end do
            else
              do iExp_k=1,Shells(kShll)%nExp

                Coeff_k = Shells(kShll)%Cff_c(iExp_k,kC,1)

                do iExp_l=1,Shells(lShll)%nExp

                  Coeff_l = Shells(lShll)%Cff_c(iExp_l,lC,1)

                  Coeff_kl = Coeff_l*Coeff_k

                  jkl = jkl+1
                  Shells(iShll)%Cff_c(jkl,iCntrc,1) = Coeff_kl
                end do
              end do
            end if

          end do ! iCntrc
#         ifdef _DEBUGPRINT_
          call RecPrt('aCD Coefficients','(6G20.12)',Shells(iShll)%Cff_c(1,1,1),nPrim,nCntrc)
#         endif
          !                                                            *
          !*************************************************************
          !                                                            *
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        Shells(iShll)%Cff_c(:,:,2) = Shells(iShll)%Cff_c(:,:,1)

        call mma_deallocate(Con)
        call mma_deallocate(ConR)
        if (Do_acCD_Basis) call mma_deallocate(Prm)

        ! Put in unit matrix of uncontracted set

        Shells(iShll)%Cff_p(:,:,1) = Zero
        do i=1,nPrim
          Shells(iShll)%Cff_p(i,i,1) = One
        end do

        Shells(iShll)%Cff_p(:,:,2) = Shells(iShll)%Cff_p(:,:,1)
        call Nrmlz(Shells(iShll)%Exp,nPrim,Shells(iShll)%Cff_p(1,1,1),nPrim,lAng)
#       ifdef _DEBUGPRINT_
        call RecPrt('uncon1',' ',Shells(iShll)%Cff_p(:,:,1),nPrim,nPrim)
        call RecPrt('uncon2',' ',Shells(iShll)%Cff_p(:,:,2),nPrim,nPrim)
#       endif

        ! OK let's do the correction now!

#       ifdef _DEBUGPRINT_
        call RecPrt('Coefficients 10',' ',Shells(iShll)%Cff_c(:,:,1),nPrim,nCntrc)
        iOff = nPrim*nCntrc
        call RecPrt('Coefficients 20',' ',Shells(iShll)%Cff_c(:,:,2),nPrim,nCntrc)
#       endif
        iOff = nPrim*nCntrc
        call Fix_Coeff(nPrim,nCntrc,Shells(iShll)%Cff_c(:,:,2),Shells(iShll)%Cff_p(:,:,1),'F')
#       ifdef _DEBUGPRINT_
        call RecPrt('Coefficients 1',' ',Shells(iShll)%Cff_c(:,:,1),nPrim,nCntrc)
        iOff = nPrim*nCntrc
        call RecPrt('Coefficients 2','(6G20.13)',Shells(iShll)%Cff_c(:,:,2),nPrim,nCntrc)
#       endif

        ! Now remove any primitives with all zero coefficents!

        call Fix_Exponents(nPrim,mPrim,nCntrc,Shells(iShll)%Exp,Shells(iShll)%Cff_c,Shells(iShll)%Cff_p)
        nPrim = mPrim
        Shells(iShll)%nExp = nPrim
#       ifdef _DEBUGPRINT_
        call RecPrt('Coefficients 1',' ',Shells(iShll)%Cff_c(:,:,1),nPrim,nCntrc)
        iOff = nPrim*nCntrc
        call RecPrt('Coefficients 2',' ',Shells(iShll)%Cff_c(:,:,2),nPrim,nCntrc)
#       endif
      end if

      Shells(iShll)%nBasis = Shells(iShll)%nBasis_c
      if ((jAng == 0) .and. Found) then
        Shells(iShll)%Transf = Shells(kShll)%Transf
        Shells(iShll)%Prjct = Shells(kShll)%Prjct
      else
        Shells(iShll)%Transf = .true.
        Shells(iShll)%Prjct = .false.
      end if
      Shells(iShll)%Aux = .true.

      if (Do_acCD_Basis .and. Found) then
        call mma_deallocate(iList2_p)
        call mma_deallocate(TInt_p)
      end if
      Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)

    end do ! jAng
  end do ! iAng

  dbsc(nCnttp)%nVal = iShll-jShll
  dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  if (Keep_Basis) then
    if (Show .and. (nPrint(2) >= 6)) then
      write(u6,*)
      write(u6,*)
      write(u6,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLbl(1:Indx-1)
      write(u6,'(1X,A)') 'On-the-fly basis set generation'
    end if

    ! Transfer the coordinate information

    nCnt = dbsc(iCnttp)%nCntr
    dbsc(nCnttp)%nCntr = nCnt
    dbsc(nCnttp)%mdci = mdc
    ! Create a pointer to the actual coordinates
    dbsc(nCnttp)%Coor => dbsc(iCnttp)%Coor(1:3,1:nCnt)

    ! Compute the number of elements stored in the dynamic memory so far.
    S%Mx_Shll = iShll+1
    Max_Shells = S%Mx_Shll
    S%Mx_mdc = mdc

  else

    ! If all the shells are empty, skip the whole basis set!

    nCnttp = nCnttp-1
  end if

end do   ! iBS

! Done for this valence basis set.

S%Mx_Shll = iShll+1
Max_Shells = S%Mx_Shll
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate

call mma_deallocate(iD_c)
call mma_deallocate(iList2_c)
!                                                                      *
!***********************************************************************
!                                                                      *
! Let us now Gram-Schmidt orthonormalize the auxiliary basis for
! better numerics and balance.

do jCnttp=nCnttp_start+1,nCnttp
  call Renorm2(jCnttp)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Optionally add auxiliary basis set to the end of the
! temporary auxiliary basis set library.

if (W2L) then
  Lu_lib = 17
  Lu_lib = IsFreeUnit(Lu_lib)
  call molcas_open(Lu_lib,'RICDLIB')
  rewind(Lu_lib)
  istatus = 0
  do while (istatus == 0)
    read(Lu_lib,*,iostat=istatus)
  end do
  backspace(Lu_lib)

  do jCnttp=nCnttp_start+1,nCnttp
    if (jCnttp == nCnttp_start+1) then
      write(Lu_lib,'(A)') '/'//Label
    else
      write(Lu_lib,'(A)') Label
    end if
    if (jCnttp == nCnttp_start+1) then
      write(Lu_lib,'(F6.2,2I10)') dbsc(jCnttp)%Charge,dbsc(jCnttp)%nVal-1,nCnttp-nCnttp_start
    else
      write(Lu_lib,'(F6.2, I10)') dbsc(jCnttp)%Charge,dbsc(jCnttp)%nVal-1
    end if
    write(Lu_lib,*) ' Dummy reference line.'
    write(Lu_lib,*) ' Dummy reference line.'
    do iAng=0,dbsc(jCnttp)%nVal-1
      iShll_ = dbsc(jCnttp)%iVal+iAng
      nExpi = Shells(iShll_)%nExp
      iSph = 0
      if (Shells(iShll_)%Prjct) iSph = 1
      if (Shells(iShll_)%Transf) iSph = iSph+2
      write(Lu_lib,'(3I10)') nExpi,Shells(iShll_)%nBasis,iSph

      ! Skip if the shell is empty.

      if (nExpi == 0) cycle

      ! Write out the exponents

      write(Lu_lib,'(5(1X,ES20.13))') (Shells(iShll_)%Exp(i),i=1,nExpi)

      ! Write out the contraction coefficients

      do i=1,nExpi
        write(Lu_lib,'(5(1X,ES20.13))') (Shells(iShll_)%Cff_c(i,j,1),j=1,Shells(iShll_)%nBasis)
      end do

    end do
  end do
  close(Lu_lib)
end if

call SOAO_Info_Free()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_aCD_acCD_Shells
