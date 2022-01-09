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

subroutine SOCtl_Seward(Mamn,nMamn)

use Basis_Info
use Center_Info
use Symmetry_Info, only: iChTbl, iOper, iChBas, lIrrep, lBsFnc, iSkip, nIrrep
use SOAO_Info, only: SOAO_Info_Init, nSOInf, iSOInf, iAOtSO, iOffSO
use real_spherical, only: iSphCr, LblCBs, LblSBs
use Temporary_Parameters, only: Primitive_Pass
use Sizes_of_Seward, only: S

implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "rinfo.fh"
#include "real.fh"
#include "print.fh"
#include "stdalloc.fh"
logical lFAIEMP
character ChOper(0:7)*3, ChTemp*8, Mamn(nMamn)*(LENIN8)
character LP_Names(MxAtom)*(LENIN4)
character*60 Fmt
logical type(0:7), lSkip, kECP, TstFnc, output, Get_BasisType
logical IsBasisAE
logical IsBasisANO
logical IsBasisUNK
integer Occ, Vir
parameter(Occ=1,Vir=0)
integer List(0:iTabMx), nFCore(0:7), nCore_Sh(0:iTabMx), List_AE(0:iTabMx)
integer jOffSO(0:7)
integer, dimension(:), allocatable :: Index, Index2, IndC, iCI, jCI, iOT, LPA, LPMM
real*8, dimension(:), allocatable :: LPQ
real*8, dimension(:,:), allocatable :: SM, LPC
character*(LENIN8) Clean_BName, ChTmp
external Clean_BName
!SVC: the basis ids are tuples (c,n,l,m) with c the center index,
!     n the shell index, l the angmom value, and m the angmom component.
!     the angmom components of p are mapped (x,y,z) -> (1,-1,0)
!     examples: 3d1+ on atom 1: (1,3,2,1); 2py on atom 5: (5,2,1,-1)
!IFG: for Cartesian shells, l -> -l, m -> T(ly+lz)-(lx+ly), where T(n)=n*(n+1)/2
integer :: llab, mlab
integer, allocatable :: basis_ids(:,:), desym_basis_ids(:,:)
integer, allocatable :: fermion_type(:)
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! LVAL end MVAL dimensioned for L = iTabMx
dimension LVAL((iTabMx+1)*(iTabMx+1))
dimension MVAL((iTabMx+1)*(iTabMx+1))

!                                                                      *
!***********************************************************************
!                                                                      *
IsBasisAE = .false.
IsBasisANO = .false.
IsBasisUNK = .false.
iRout = 2
iPrint = nPrint(iRout)
!vv LP_NAMES was used later without initialization.
do i=1,MxAtom
  LP_NAMES(i)(1:LENIN) = 'crap'
  LP_NAMES(i)(LENIN1:LENIN4) = 'crap'
end do
lFAIEMP = .false.
do i=1,nCnttp
  lFAIEMP = lFAIEMP .or. dbsc(i)%Frag
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute iBas, iBas_Aux, and iBas_Frag used for double checking
! in SOCtl.
! Compute cdMax, EtMax, and S%nShlls.

call Misc_Seward(iBas,iBas_Aux,iBas_Frag)
call SOAO_Info_Init(iBas+iBas_Frag+iBas_Aux,nIrrep)
!                                                                      *
!***********************************************************************
!                                                                      *
! initialize LVAL and MVAL
! (note: this is wrong for Cartesian shells)

k = 0
do i=0,iTabMx
  do j=-i,i
    k = k+1
    lval(k) = i
    mval(k) = j
  end do
end do
!write(6,*) ' lval',k,(iTabMx+1)**2
! correct mval order for p-functions
mval(2) = 1
mval(3) = -1
mval(4) = 0
call ICopy(MxAO,[-99],0,iCent,1)
call ICopy(MxAO,[-99],0,lnAng,1)
!write(6,'(20i4)') (lval(i),i=1,k)
!write(6,*) ' lval',k
!write(6,'(20i4)') (mval(i),i=1,k)

call ICopy(1+iTabMx,[0],0,List,1)
call ICopy(1+iTabMx,[0],0,List_AE,1)

isymunit = isfreeunit(58)
call molcas_open(isymunit,'SYMINFO')
rewind isymunit
write(isymunit,'(A)') 'Symmetry information from seward'
write(isymunit,'(A)') '#of funct, unique centre, L, M , # of sym.ad.functions , Phases'
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of symmetry adapted or petite list basis functions
!                                                                      *
!***********************************************************************
!                                                                      *
iSO = 0
iSO_Aux = iBas+iBas_Frag
iSO_Frag = iBas
iSO_Tot = 0
S%n2Tot = 0
S%nDim = 0
iAO = 0
lSkip = .false.

call ICopy(8,[0],0,nFCore,1)

call mma_Allocate(iCI,iBas,label='iCI')        ! Stuff for LoProp
call mma_Allocate(jCI,iBas,label='jCI')        ! Stuff for LocalDKH/X2C/BSS
call mma_Allocate(iOT,iBas,label='iOT')        ! Stuff for LoProp
call mma_Allocate(LPC,3,S%mCentr,label='LPC')  ! Stuff (not just) for LoProp
call mma_Allocate(LPQ,S%mCentr,label='LPQ')    ! Stuff (not just) for LoProp
call mma_Allocate(LPMM,S%mCentr,label='LPMM')  ! Stuff (not just) for LoProp
call mma_Allocate(LPA,S%mCentr,label='LPA')
call mma_allocate(basis_ids,4,maxbfn+maxbfn_aux)
call mma_allocate(desym_basis_ids,4,maxbfn+maxbfn_aux)
call mma_allocate(fermion_type,maxbfn+maxbfn_aux)

IsBasisAE = Get_BasisType('AE_')
IsBasisANO = Get_BasisType('ANO')
IsBasisUNK = Get_BasisType('UNK')
if (Show .and. (iPrint >= 6)) then
  write(6,*)
  call CollapseOutput(1,'   SO/AO info:')
  write(6,'(3X,A)') '   -----------'
end if
if (nIrrep == 1) Go To 199
!                                                                      *
!***********************************************************************
!                                                                      *
! Symmetry case.

if (Show .and. (iPrint >= 6)) then
  write(6,*)
  write(6,'(19x,a)') ' **************************************************'
  write(6,'(19x,a)') ' ******** Symmetry adapted Basis Functions ********'
  write(6,'(19x,a)') ' **************************************************'
  write(6,*)
end if

call mma_allocate(Index,5*iBas,label='Index')
call mma_allocate(Index2,5*iBas,label='Index2')
call ICopy(5*iBas,[0],0,Index,1)
call ICopy(5*iBas,[0],0,Index2,1)
iCounter = 0
jCounter = 0
call mma_Allocate(SM,iBas,iBas,label='SM')
call FZero(SM,iBas**2)
call mma_Allocate(IndC,2*S%mCentr)
iAtoms = 0

! Loop over irreducible representations and symmetry operations,
! respectively, for SO and Petite list, respectively.

do iIrrep=0,7
  jOffSO(iIrrep) = 0
end do
kIrrep = 0
do iIrrep=0,nIrrep-1
  iOffSO(iIrrep) = iSO_Tot
  jOffSO(iIrrep) = iSO
  iAO = 0
  jSO = 0
  nBas(iIrrep) = 0
  nBas_Aux(iIrrep) = 0
  nBas_Frag(iIrrep) = 0
  type(iIrrep) = .true.

  ! Loop over distinct shell types

  mc = 1
  iShell = 0
  if (iSkip(iIrrep) /= 0) then
    write(6,*)
    write(6,*) ' All basis functions of Irrep',iIrrep+1,' are removed!'
    write(6,*)
    lSkip = .true.
    Go To 2011
  end if
  iCnttp = 0
  do jCnttp=1,nCnttp

    ! Make sure that we process the dummy shell last

    if ((jCnttp == iCnttp_Dummy) .and. (jCnttp /= nCnttp)) then
      iCnttp = iCnttp+2
    else if ((jCnttp == nCnttp) .and. (iCnttp == jCnttp)) then
      iCnttp = iCnttp_Dummy
    else
      iCnttp = iCnttp+1
    end if

    output = show .and. (iPrint >= 6)
    if (dbsc(iCnttp)%Aux) output = output .and. (iPrint >= 10) .and. (iCnttp /= iCnttp_Dummy)
    if (dbsc(iCnttp)%Frag) output = output .and. (iPrint >= 10)
    kECP = dbsc(iCnttp)%ECP
    lMax = dbsc(iCnttp)%nVal-1

    call OrbType(dbsc(iCnttp)%AtmNr,List_AE,31)
    if (kECP) then

      ! ECP case

      call ECP_Shells(dbsc(iCnttp)%AtmNr,list)

      ! No core to freeze!

      call ICopy(lMax+1,[0],0,nCore_Sh,1)
    else

      ! Non-ECP case

      ! Pick up the number of occupied orbitals in each shell type.

      call ICopy(1+iTabMx,List_AE,1,List,1)

      ! Pick up which orbitals should be frozen as default.

      if (dbsc(iCnttp)%Charge /= Zero) then
        call Freeze_Default(dbsc(iCnttp)%AtmNr,nCore_Sh,lMax)
      else

        ! If there charge is zero we presume that these are
        ! ghost orbitals or something else. In any case
        ! we do not freeze any orbitals!

        call Freeze_Default(0,nCore_Sh,lMax)
      end if
    end if

    ! Loop over distinct centers

    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = iCnt+dbsc(iCnttp)%mdci

      ! Loop over shells associated with this center
      ! Start with s type shells

      kComp = 0
      kculf = 0
      iSh = dbsc(iCnttp)%iVal-1
      if (dbsc(iCnttp)%nVal < 1) then
        do iCo=0,nIrrep/dc(mdc)%nStab-1
          iyy = Index_Center(mdc,iCo,IndC,iAtoms,S%mCentr)
          iR = NrOpr(dc(mdc)%iCoSet(iCo,0))

          LPC(1:3,iyy) = dbsc(iCnttp)%Coor(1:3,iCnt)
          if (iand(iOper(iR),1) /= 0) LPC(1,iyy) = -LPC(1,iyy)
          if (iand(iOper(iR),2) /= 0) LPC(2,iyy) = -LPC(2,iyy)
          if (iand(iOper(iR),4) /= 0) LPC(3,iyy) = -LPC(3,iyy)
          LPQ(iyy) = dbsc(iCnttp)%Charge
          LPA(iyy) = dbsc(iCnttp)%AtmNr
          LPMM(iyy) = dbsc(iCnttp)%IsMM
          LP_Names(iyy) = dc(mdc)%LblCnt(1:LENIN)//':'//ChOper(iOper(iR))
        end do
      end if
      do iAng=0,dbsc(iCnttp)%nVal-1
        nCore = nCore_Sh(iAng)
        iSh = iSh+1
        iShell = iShell+1
        nExpi = Shells(iSh)%nExp
        nBasisi = Shells(iSh)%nBasis
        if (Shells(iSh)%Prjct) then
          jComp = 2*iAng+1
        else
          jComp = (iAng+1)*(iAng+2)/2
        end if
        if (nExpi == 0) Go To 2033
        if (nBasisi == 0) Go To 2033

        do iComp=1,jComp
          iAO = iAO+1
          if (iAO > MxAO) then
            call ErrTra()
            write(6,*) ' Increase MxAO'
            call Abend()
          end if
          lComp = kComp+iComp
          lculf = kculf+icomp
          ! Get character of basis function
          iChBs = iChBas(lComp)
          if (Shells(iSh)%Transf) iChBs = iChBas(iSphCr(lComp))

          ! Skip if function not a basis of irreps.

          if (.not. TstFnc(dc(mdc)%iCoSet,iIrrep,iChBs,dc(mdc)%nStab)) Go To 204
          if (.not. (Shells(iSh)%Frag .or. dbsc(iCnttp)%Aux)) nFCore(iIrrep) = nFCore(iIrrep)+nCore
          if (output .and. type(iIrrep)) then
            write(6,*)
            write(6,'(10X,A,A)') ' Irreducible representation : ',lIrrep(iIrrep)
            write(6,'(10X,2A)') ' Basis function(s) of irrep: ',lBsFnc(iIrrep)
            write(6,*)
            write(6,'(A)') ' Basis Label        Type   Center Phase'
            type(iIrrep) = .false.
          end if

          if (S%MaxBas(iAng) > 0) iAOtSO(iAO,iIrrep) = jSO+1
          S%m2Max = max(S%m2Max,nExpi**2)
          do iCntrc=1,nBasisi
            iSO_Tot = iSO_Tot+1
            if (Shells(iSh)%Aux) then
              iSO_Aux = iSO_Aux+1
              iSO_ = iSO_Aux
              nBas_Aux(iIrrep) = nBas_Aux(iIrrep)+1
            else if (Shells(iSh)%Frag) then
              iSO_Frag = iSO_Frag+1
              iSO_ = iSO_Frag
              nBas_Frag(iIrrep) = nBas_Frag(iIrrep)+1
            else
              iSO = iSO+1
              iSO_ = iSO
              nBas(iIrrep) = nBas(iIrrep)+1
            end if
            if (iSO_ > nMamn) then
              write(6,*) ' iSO_ > nMamn'
              write(6,*) 'nMamn=',nMamn
              call Abend()
            end if
            jSO = jSO+1

            ChTemp = LblCBs(lComp)
            if (Shells(iSh)%Transf) ChTemp = LblSbs(lComp)

            call Name_to_lm(ChTemp,llab,mlab)

            ! Introduce a somewhat better labelling. Thnx LG!

            if (IsBasisAE) then
              if (IsBasisANO) then
                write(ChTemp(1:2),'(I2.2)') iAng+iCntrc
              else
                if (nExpi == nBasisi) then
                  write(ChTemp(1:1),'(A1)') '*'
                  if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
                else if (iCntrc <= list(iAng)) then
                  write(ChTemp(1:2),'(I2.2)') iAng+iCntrc
                else
                  write(ChTemp(1:1),'(A1)') '*'
                  if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
                end if
              end if
            else if (.not. IsBasisUNK) then
              if (nExpi == nBasisi) then
                write(ChTemp(1:1),'(A1)') '*'
                if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
              else if (iCntrc <= list(iAng)) then
                write(ChTemp(1:2),'(I2.2)') iAng+iCntrc+(List_AE(iAng)-List(iAng))
              else
                write(ChTemp(1:1),'(A1)') '*'
                if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
              end if

            end if
            ChTmp = Clean_BName(ChTemp,0)

            if (output) write(6,'(I5,3X,A8,4X,A8,8(I3,4X,I2,4X))') iSO_,dc(mdc)%LblCnt,ChTmp, &
              (mc+iCo,iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
               iCo=0,nIrrep/dc(mdc)%nStab-1)

            if (iSO_ > nSOInf) then
              write(6,*) 'iSO_ > nSOInf'
              call Abend()
            end if
            iSOInf(1,iSO_) = iCnttp
            iSOInf(2,iSO_) = iCnt
            iSOInf(3,iSO_) = iAng

            if (Shells(iSh)%Aux .or. Shells(iSh)%Frag) Go To 205

            if (.not. Primitive_Pass) then
              write(isymunit,'(13(I4,4X))') iSO_,mdc,LVAL(lculf),MVAL(lculf),nIrrep/dc(mdc)%nStab, &
                (iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))),iCo=0,nIrrep/dc(mdc)%nStab-1)
            end if
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Stuff (not just) for LoProp

            do iCo=0,nIrrep/dc(mdc)%nStab-1
              ixxx = Index_NoSym(iCntrc,iComp,iAng,mdc,iCo,Index,iCounter,iBas)
              jxxx = Index_NoSym(iCntrc,iComp,iAng,mdc,iirrep,Index2,jCounter,iBas)
              fact = dble(iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))))

              FacN = One/dble(nIrrep/dc(mdc)%nStab)
              if (MolWgh == 1) then
                FacN = One
              else if (MolWgh == 2) then
                FacN = sqrt(FacN)
              end if
              SM(ixxx,iSO) = Fact*FacN
              iyy = Index_Center(mdc,iCo,IndC,iAtoms,S%mCentr)

              iCI(ixxx) = iyy
              jCI(jxxx) = icnt

              if (iCntrc <= list(iAng)) then
                iOT(ixxx) = Occ
              else
                iOT(ixxx) = Vir
              end if

              iR = NrOpr(dc(mdc)%iCoSet(iCo,0))
              LPC(1:3,iyy) = dbsc(iCnttp)%Coor(1:3,iCnt)
              if (iand(iOper(iR),1) /= 0) LPC(1,iyy) = -LPC(1,iyy)
              if (iand(iOper(iR),2) /= 0) LPC(2,iyy) = -LPC(2,iyy)
              if (iand(iOper(iR),4) /= 0) LPC(3,iyy) = -LPC(3,iyy)

              LPQ(iyy) = dbsc(iCnttp)%Charge
              LPMM(iyy) = dbsc(iCnttp)%IsMM
              LPA(iyy) = dbsc(iCnttp)%AtmNr

              LP_Names(iyy) = dc(mdc)%LblCnt(1:LENIN)//':'//ChOper(iOper(iR))
              desym_basis_ids(1,ixxx) = iyy
              desym_basis_ids(2,ixxx) = iCntrc
              desym_basis_ids(3,ixxx) = llab
              desym_basis_ids(4,ixxx) = mlab
            end do
            !                                                          *
            !***********************************************************
            !                                                          *
            Mamn(iSO) = dc(mdc)%LblCnt(1:LENIN)//ChTemp(1:8)
            basis_ids(1,iSO) = mdc
            basis_ids(2,iSO) = iCntrc
            basis_ids(3,iSO) = llab
            basis_ids(4,iSO) = mlab
            fermion_type(iSO) = 0
            if (dbsc(iCnttp)%fMass /= 1.0d0) fermion_type(iSO) = 1
            if (.not. Primitive_Pass) then
              kIrrep = kIrrep+1
              icent(kIrrep) = mdc
              lnang(kIrrep) = lval(lculf)
              lmag(kIrrep) = mval(lculf)
              lant(kIrrep) = nIrrep/dc(mdc)%nStab
            end if
205         continue
          end do

204       continue
        end do
2033    kComp = kComp+(iAng+1)*(iAng+2)/2
        kculf = kculf+2*iAng+1
      end do ! iAng
      mc = mc+nIrrep/dc(mdc)%nStab
    end do ! iCnt

  end do ! jCnttp
2011 continue
  !ulf
  nrSym = nIrrep
  nrBas(iIrrep+1) = nBas(iIrrep)
  !write(6,*) ' nBas(iIrrep)', iIrrep, nBas(iIrrep)
  S%nDim = S%nDim+nBas(iIrrep)
  S%n2Tot = S%n2Tot+nBas(iIrrep)**2
end do ! iIrrep
!if (lSkip) S%nDim = iBas
if ((iBas /= iSO) .and. (iBas_Aux /= iSO_Aux-iSO) .and. (.not. lSkip)) then
  write(6,*) 'iBas=',iBas
  write(6,*) 'iBas_Aux=',iBas_Aux
  write(6,*) 'iSO=',iSO
  write(6,*) 'iSO_Aux=',iSO_Aux-iSO
  write(6,*) 'iSO_Tot=',iSO_Tot
  call ErrTra()
  call Abend()
end if
! redefine iOffSO array in case of Fragment AIEMP
if (lFAIEMP) then
  do iIrrep=0,nIrrep-1
    iOffSO(iIrrep) = jOffSO(iIrrep)
  end do
end if
#ifdef _DEBUGPRINT_
call RecPrt('Symmetrization Matrix','(20F5.2)',SM,iBas,iBas)
#endif
call Put_dArray('SM',SM,iBas**2)

!SVC: basis IDs of both symmetric and non-symmetric case
if (.not. Primitive_Pass) then
  call Put_iArray('Fermion IDs',fermion_type,iSO)
  call Put_iArray('Basis IDs',basis_ids,4*iSO)
  call Put_iArray('Desym Basis IDs',desym_basis_ids,4*iBas)
end if

call mma_deallocate(IndC)
call mma_deallocate(SM)
call mma_deallocate(Index)
call mma_deallocate(Index2)
Go To 198
!                                                                      *
!***********************************************************************
!                                                                      *
! No symmetry case.

199 continue
if (Show .and. (iPrint >= 6)) then
  write(6,*)
  write(6,'(19x,a)') ' **************************************************'
  write(6,'(19x,a)') ' ********** Petite list Basis Functions ***********'
  write(6,'(19x,a)') ' **************************************************'
  write(6,*)
end if

kIrrep = 0
do iIrrep=0,nIrrep-1
  iOffSO(iIrrep) = iSO_Tot
  iAO = 0
  jSO = 0
  nBas(iIrrep) = 0
  nBas_Aux(iIrrep) = 0
  nBas_Frag(iIrrep) = 0
  type(iIrrep) = .true.

  ! Loop over distinct shell types

  mc = 1
  iShell = 0
  iCnttp = 0
  do jCnttp=1,nCnttp

    ! Make sure that we process the dummy shell last

    if ((jCnttp == iCnttp_Dummy) .and. (jCnttp /= nCnttp)) then
      iCnttp = iCnttp+2
    else if ((jCnttp == nCnttp) .and. (iCnttp == jCnttp)) then
      iCnttp = iCnttp_Dummy
    else
      iCnttp = iCnttp+1
    end if

    output = show .and. (iPrint >= 6)
    if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) output = output .and. (iPrint >= 10) .and. (iCnttp /= iCnttp_Dummy)
    kECP = dbsc(iCnttp)%ECP
    lMax = dbsc(iCnttp)%nVal-1
    call OrbType(dbsc(iCnttp)%AtmNr,List_AE,31)
    if (kECP) then
      call ECP_Shells(dbsc(iCnttp)%AtmNr,list)
      call ICopy(lmax+1,[0],0,nCore_Sh,1)
    else
      call ICopy(1+iTabMx,List_AE,1,List,1)
      if (dbsc(iCnttp)%Charge /= Zero) then
        call Freeze_Default(dbsc(iCnttp)%AtmNr,nCore_Sh,lMax)
      else
        call Freeze_Default(0,nCore_Sh,lMax)
      end if
    end if

    ! Loop over distinct centers

    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = iCnt+dbsc(iCnttp)%mdci

      ! Loop over shells associated with this center
      ! Start with s type shells

      kComp = 0
      kculf = 0
      iSh = dbsc(iCnttp)%iVal-1
      if (dbsc(iCnttp)%nVal < 1) then
        LPC(1:3,mdc) = dbsc(iCnttp)%Coor(1:3,iCnt)
        LPQ(mdc) = dbsc(iCnttp)%Charge
        LPMM(mdc) = dbsc(iCnttp)%IsMM
        LPA(mdc) = dbsc(iCnttp)%AtmNr
        LP_Names(mdc) = dc(mdc)%LblCnt(1:LENIN)//'    '
      end if
      do iAng=0,dbsc(iCnttp)%nVal-1
        nCore = nCore_Sh(iAng)
        iSh = iSh+1
        iShell = iShell+1
        nExpi = Shells(iSh)%nExp
        nBasisi = Shells(iSh)%nBasis
        if (Shells(iSh)%Prjct) then
          jComp = 2*iAng+1
        else
          jComp = (iAng+1)*(iAng+2)/2
        end if
        if (nExpi == 0) Go To 3033
        if (nBasisi == 0) Go To 3033
        do iComp=1,jComp
          iAO = iAO+1
          if (iAO > MxAO) then
            call ErrTra()
            write(6,*) ' Increase MxAO'
            call Abend()
          end if
          lComp = kComp+iComp
          lculf = kculf+iComp

          ! Skip if symmetry operator is not in the coset of this center.

          do imc=0,(nIrrep/dc(mdc)%nStab)-1
            if (dc(mdc)%iCoSet(imc,0) == iOper(iIrrep)) Go To 307
          end do
          Go To 304
307       continue
          if (output .and. type(iIrrep)) then
            write(6,*)
            write(6,'(10X,2A)') ' Basis functions generated by ',ChOper(iIrrep)
            write(6,*)
            write(6,'(A)') ' Basis Label        Type   Center'
            type(iIrrep) = .false.
          end if

          if (S%MaxBas(iAng) > 0) iAOtSO(iAO,iIrrep) = jSO+1
          S%m2Max = max(S%m2Max,nExpi**2)
          if (.not. (Shells(iSh)%Frag .or. dbsc(iCnttp)%Aux)) nFCore(0) = nFCore(0)+nCore

          ! Loop over contracted basis functions

          do iCntrc=1,nBasisi
            iSO_Tot = iSO_Tot+1
            if (Shells(iSh)%Aux) then
              iSO_Aux = iSO_Aux+1
              iSO_ = iSO_Aux
              nBas_Aux(iIrrep) = nBas_Aux(iIrrep)+1
            else if (Shells(iSh)%Frag) then
              iSO_Frag = iSO_Frag+1
              iSO_ = iSO_Frag
              nBas_Frag(iIrrep) = nBas_Frag(iIrrep)+1
            else
              iSO = iSO+1
              iSO_ = iSO
              nBas(iIrrep) = nBas(iIrrep)+1
            end if
            if (iSO_ > nMamn) then
              write(6,*) ' iSO_ > nMamn'
              write(6,*) 'nMamn=',nMamn
              call Abend()
            end if
            jSO = jSO+1

            ChTemp = LblCBs(lComp)
            if (Shells(iSh)%Transf) ChTemp = LblSbs(lComp)

            call Name_to_lm(ChTemp,llab,mlab)

            ! Introduce a somewhat better labelling. Thnx LG!

            if (IsBasisAE) then
              if (IsBasisANO) then
                write(ChTemp(1:2),'(I2.2)') iAng+iCntrc
              else
                if (nExpi == nBasisi) then
                  write(ChTemp(1:1),'(A1)') '*'
                  if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
                else if (iCntrc <= list(iAng)) then
                  write(ChTemp(1:2),'(I2.2)') iAng+iCntrc
                else
                  write(ChTemp(1:1),'(A1)') '*'
                  if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
                end if
              end if
            else if (.not. IsBasisUNK) then
              if (nExpi == nBasisi) then
                write(ChTemp(1:1),'(A1)') '*'
                if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
              else if (iCntrc <= list(iAng)) then
                write(ChTemp(1:2),'(I2.2)') iAng+iCntrc+(List_AE(iAng)-List(iAng))
              else
                write(ChTemp(1:1),'(A1)') '*'
                if (llab >= 0) write(ChTemp(2:2),'(A1)') '0'
              end if
            end if
            ChTmp = Clean_BName(ChTemp,0)

            if (output) write(6,'(I5,2X,A8,5X,A8,I3)') iSO_,dc(mdc)%LblCnt,ChTmp,mc+imc

            if (iSO_ > nSOInf) then
              write(6,*) 'iSO_ > nSOInf'
              call Abend()
            end if
            iSOInf(1,iSO_) = iCnttp
            iSOInf(2,iSO_) = iCnt
            iSOInf(3,iSO_) = iAng

            if (Shells(iSh)%Aux .or. Shells(iSh)%Frag) Go To 305
            write(isymunit,'(13(I4,4X))') iSO,mdc,LVAL(lculf),MVAL(lculf),nIrrep/dc(mdc)%nStab, &
              (iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))),iCo=0,nIrrep/dc(mdc)%nStab-1)
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Stuff (not just) for LoProp

            iCI(iSO) = mdc
            jCI(iSO) = mdc
            if (iCntrc <= list(iAng)) then
              iOT(iSO) = Occ
            else
              iOT(iSO) = Vir
            end if
            LPC(1:3,mdc) = dbsc(iCnttp)%Coor(1:3,iCnt)
            LPQ(mdc) = dbsc(iCnttp)%Charge
            LPMM(mdc) = dbsc(iCnttp)%IsMM
            LPA(mdc) = dbsc(iCnttp)%AtmNr
            LP_Names(mdc) = dc(mdc)%LblCnt(1:LENIN)//'    '
            !                                                          *
            !***********************************************************
            !                                                          *
            Mamn(iSO) = dc(mdc)%LblCnt(1:LENIN)//ChTemp(1:8)
            basis_ids(1,iSO) = mdc
            basis_ids(2,iSO) = iCntrc
            basis_ids(3,iSO) = llab
            basis_ids(4,iSO) = mlab
            fermion_type(iSO) = 0
            if (dbsc(iCnttp)%fMass /= 1.0d0) fermion_type(iSO) = 1
            if (.not. Primitive_Pass) then
              kIrrep = kIrrep+1
              icent(kIrrep) = mdc
              lnang(kIrrep) = lval(lculf)
              lmag(kIrrep) = mval(lculf)
              lant(kIrrep) = nIrrep/dc(mdc)%nStab
            end if
305         continue
          end do

304       continue
        end do
3033    kComp = kComp+(iAng+1)*(iAng+2)/2
        kculf = kculf+2*iAng+1
      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do

  end do
  nrSym = nIrrep
  nrBas(iIrrep+1) = nBas(iIrrep)
  S%nDim = S%nDim+nBas(iIrrep)
  S%n2Tot = S%n2Tot+nBas(iIrrep)**2
end do

!SVC: basis IDs of non-symmetric case
if (.not. Primitive_Pass) then
  call Put_iArray('Fermion IDs',fermion_type,iBas)
  call Put_iArray('Basis IDs',basis_ids,4*iBas)
end if
!
do jAO=1,iAO
  do jIrrep=0,nIrrep-1
    iAOtSO(jAO,jIrrep) = iAOtSO(jAO,jIrrep)+iOffSO(jIrrep)
  end do
end do

! Fix index list such that redundant operators will have the
! same AO index.

mc = 1
iShell = 0
iAO = 0
do iCnttp=1,nCnttp
  kECP = dbsc(iCnttp)%ECP

  ! Loop over distinct centers

  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = iCnt+dbsc(iCnttp)%mdci
    iChxyz = dc(mdc)%iChCnt

    ! Loop over shells associated with this center
    ! Start with s type shells

    kComp = 0
    iSh = dbsc(iCnttp)%iVal-1
    do iAng=0,dbsc(iCnttp)%nVal-1
      iSh = iSh+1
      iShell = iShell+1
      nExpi = Shells(iSh)%nExp
      nBasisi = Shells(iSh)%nBasis
      if (nExpi == 0) Go To 4033
      if (nBasisi == 0) Go To 4033
      jComp = (iAng+1)*(iAng+2)/2
      if (Shells(iSh)%Prjct) jComp = 2*iAng+1
      do iComp=1,jComp
        iAO = iAO+1
        if (iAO > MxAO) then
          call ErrTra()
          write(6,*) ' Increase MxAO'
          call Abend()
        end if
        lComp = kComp+iComp
        if (S%MaxBas(iAng) > 0) then

          do imc=0,(nIrrep/dc(mdc)%nStab)-1
            itest1 = iand(dc(mdc)%iCoSet(imc,0),iChxyz)
            Nr = NrOpr(dc(mdc)%iCoSet(imc,0))
            do jIrrep=0,nIrrep-1
              itest2 = iand(iOper(jIrrep),iChxyz)
              if (itest1 == itest2) iAOtSO(iAO,jIrrep) = iAOtSO(iAO,Nr)
            end do
          end do
        end if
      end do
4033  kComp = kComp+(iAng+1)*(iAng+2)/2
    end do
    mc = mc+nIrrep/dc(mdc)%nStab
  end do

end do

198 continue
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  if (iPrint >= 6) then ! std print has been moved to seward.f

    ! Print out basis set information

    Fmt = '(6X,A,T30,8I4)'
    write(6,*)
    write(6,'(6X,A)') 'Basis set specifications :'
    write(6,'(6X,A,T30,8(1X,A))') 'Symmetry species',(lIrrep(i),i=0,nIrrep-1)
    write(6,Fmt) 'Basis functions',(nBas(i),i=0,nIrrep-1)

  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write info (not just) for LoProp

if (.not. Primitive_Pass) then
  call Put_cArray('LP_L',LP_Names(1),(LENIN4)*S%mCentr)
  call Put_iArray('LP_A',LPA,S%mCentr)
  call Put_dArray('LP_Q',LPQ,S%mCentr)
  call Put_dArray('LP_Coor',LPC,3*S%mCentr)
  call Put_iScalar('LP_nCenter',S%mCentr)
  call Put_iArray('IsMM Atoms',LPMM,S%mCentr)
  call Put_iArray('Center Index',iCI,iBas)
  call Put_iArray('Orbital Type',iOT,iBas)
  call Put_iArray('Non valence orbitals',nFCore,nIrrep)
else
  call Put_iArray('Ctr Index Prim',jCI,iBas)

end if
call mma_deallocate(fermion_type)
call mma_deallocate(desym_basis_ids)
call mma_deallocate(basis_ids)
call mma_deallocate(LPA)
call mma_deallocate(LPMM)
call mma_deallocate(LPQ)
call mma_deallocate(LPC)
call mma_deallocate(iCI)
call mma_deallocate(jCI)
call mma_deallocate(iOT)

if (Show .and. (iPrint >= 6)) then
  call CollapseOutput(0,'   SO/AO info:')
  write(6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*) ' *** iAOtSO ***'
do jAO=1,iAO
  write(6,*) (iAOtSO(jAO,jIrrep),jIrrep=0,nIrrep-1)
end do
#endif

write(isymunit,'(A3)') 'END'
close(isymunit)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine SOCtl_Seward
