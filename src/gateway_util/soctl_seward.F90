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

use Basis_Info, only: dbsc, iCnttp_Dummy, nBas, nBas_Aux, nBas_Frag, nCnttp, MolWgh, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iChBas, iChTbl, iOper, iSkip, lBsFnc, lIrrep, nIrrep
use SOAO_Info, only: iAOtSO, iSOInf, iOffSO, nSOInf, SOAO_Info_Init
use real_spherical, only: iSphCr, LblCBs, LblSBs
use Gateway_global, only: Primitive_Pass
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nMamn
character(len=LenIn8), intent(out) :: Mamn(nMamn)
#include "itmax.fh"
#include "rinfo.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAng, iAO, iAtoms, iBas, iBas_Aux, iBas_Frag, iChBs, iChxyz, iCnt, iCntrc, iCnttp, iCo, iComp, iCounter, &
                     iIrrep, imc, iPrint, iR, iRout, iSh, iShell, iSO, iSO_, iSO_Aux, iSO_Frag, iSO_Tot, isymunit, itest1, itest2, &
                     ixxx, iyy, j, jAO, jCnttp, jComp, jCounter, jIrrep, jOffSO(0:7), jSO, jxxx, k, kComp, kculf, kIrrep, lComp, &
                     lculf, llab, lMax, mc, mdc, mlab, nBasisi, nCore, nExpi, nFCore(0:7), Nr
real(kind=wp) :: FacN, fact
logical(kind=iwp) :: IsBasisAE, IsBasisANO, IsBasisUNK, kECP, lFAIEMP, lSkip, output, TstFnc, bType(0:7)
character(len=LenIn8) :: ChTmp
character(len=8) :: ChTemp
character(len=60) :: Frmt
!SVC: the basis ids are tuples (c,n,l,m) with c the center index,
!     n the shell index, l the angmom value, and m the angmom component.
!     the angmom components of p are mapped (x,y,z) -> (1,-1,0)
!     examples: 3d1+ on atom 1: (1,3,2,1); 2py on atom 5: (5,2,1,-1)
!IFG: for Cartesian shells, l -> -l, m -> T(ly+lz)-(lx+ly), where T(n)=n*(n+1)/2
integer(kind=iwp), allocatable :: basis_ids(:,:), desym_basis_ids(:,:), fermion_type(:), iCI(:), IndC(:), Index1(:), Index2(:), &
                                  iOT(:), jCI(:), List(:), List_AE(:), LPA(:), LPMM(:), LVAL(:), MVAL(:), nCore_Sh(:)
real(kind=wp), allocatable :: LPC(:,:), LPQ(:), SM(:,:)
character(len=LenIn4), allocatable :: LP_Names(:)
integer(kind=iwp), parameter :: Occ = 1, Vir = 0
character(len=*), parameter :: ChOper(0:7) = ['E  ', &
                                              'x  ', &
                                              'y  ', &
                                              'xy ', &
                                              'z  ', &
                                              'xz ', &
                                              'yz ', &
                                              'xyz']
integer(kind=iwp), external :: Index_Center, Index_Nosym, iPrmt, isfreeunit, NrOpr
logical(kind=iwp), external :: Get_BasisType
character(len=LenIn8), external :: Clean_BName

!                                                                      *
!***********************************************************************
!                                                                      *
IsBasisAE = .false.
IsBasisANO = .false.
IsBasisUNK = .false.
iRout = 2
iPrint = nPrint(iRout)
!vv LP_NAMES was used later without initialization.
call mma_allocate(LP_Names,MxAtom,label='LP_Names')
LP_NAMES(:)(1:LenIn) = 'crap'
LP_NAMES(:)(LenIn1:) = 'crap'
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

call mma_allocate(LVAL,(iTabMx+1)**2,label='LVAL')
call mma_allocate(MVAL,(iTabMx+1)**2,label='MVAL')
k = 0
do i=0,iTabMx
  do j=-i,i
    k = k+1
    lval(k) = i
    mval(k) = j
  end do
end do
!write(u6,*) ' lval',k,(iTabMx+1)**2
! correct mval order for p-functions
mval(2) = 1
mval(3) = -1
mval(4) = 0
iCent(:) = -99
lnAng(:) = -99
!write(u6,'(20i4)') (lval(i),i=1,k)
!write(u6,*) ' lval',k
!write(u6,'(20i4)') (mval(i),i=1,k)

call mma_allocate(nCore_Sh,[0,iTabMx],label='nCore_Sh')
call mma_allocate(List,[0,iTabMx],label='List')
call mma_allocate(List_AE,[0,iTabMx],label='List_AE')
List(:) = 0
List_AE(:) = 0

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

nFCore(:) = 0

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
  write(u6,*)
  call CollapseOutput(1,'   SO/AO info:')
  write(u6,'(3X,A)') '   -----------'
end if
if (nIrrep /= 1) then

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Symmetry case.

  if (Show .and. (iPrint >= 6)) then
    write(u6,*)
    write(u6,'(19x,a)') ' **************************************************'
    write(u6,'(19x,a)') ' ******** Symmetry adapted Basis Functions ********'
    write(u6,'(19x,a)') ' **************************************************'
    write(u6,*)
  end if

  call mma_allocate(Index1,5*iBas,label='Index1')
  call mma_allocate(Index2,5*iBas,label='Index2')
  Index1(:) = 0
  Index2(:) = 0
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
    bType(iIrrep) = .true.

    ! Loop over distinct shell types

    mc = 1
    iShell = 0
    if (iSkip(iIrrep) /= 0) then
      write(u6,*)
      write(u6,*) ' All basis functions of Irrep',iIrrep+1,' are removed!'
      write(u6,*)
      lSkip = .true.
    else
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

          nCore_Sh(0:lMax) = 0
        else

          ! Non-ECP case

          ! Pick up the number of occupied orbitals in each shell type.

          List(:) = List_AE

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
              if (btest(iOper(iR),0)) LPC(1,iyy) = -LPC(1,iyy)
              if (btest(iOper(iR),1)) LPC(2,iyy) = -LPC(2,iyy)
              if (btest(iOper(iR),2)) LPC(3,iyy) = -LPC(3,iyy)
              LPQ(iyy) = dbsc(iCnttp)%Charge
              LPA(iyy) = dbsc(iCnttp)%AtmNr
              LPMM(iyy) = dbsc(iCnttp)%IsMM
              LP_Names(iyy) = dc(mdc)%LblCnt(1:LenIn)//':'//ChOper(iOper(iR))
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
            if ((nExpi /= 0) .and. (nBasisi /= 0)) then

              do iComp=1,jComp
                iAO = iAO+1
                if (iAO > MxAO) then
                  write(u6,*) ' Increase MxAO'
                  call Abend()
                end if
                lComp = kComp+iComp
                lculf = kculf+icomp
                ! Get character of basis function
                iChBs = iChBas(lComp)
                if (Shells(iSh)%Transf) iChBs = iChBas(iSphCr(lComp))

                ! Skip if function not a basis of irreps.

                if (.not. TstFnc(dc(mdc)%iCoSet,iIrrep,iChBs,dc(mdc)%nStab)) cycle
                if (.not. (Shells(iSh)%Frag .or. dbsc(iCnttp)%Aux)) nFCore(iIrrep) = nFCore(iIrrep)+nCore
                if (output .and. bType(iIrrep)) then
                  write(u6,*)
                  write(u6,'(10X,A,A)') ' Irreducible representation : ',lIrrep(iIrrep)
                  write(u6,'(10X,2A)') ' Basis function(s) of irrep: ',lBsFnc(iIrrep)
                  write(u6,*)
                  write(u6,'(A)') ' Basis Label        Type   Center Phase'
                  bType(iIrrep) = .false.
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
                    write(u6,*) ' iSO_ > nMamn'
                    write(u6,*) 'nMamn=',nMamn
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

                  if (output) write(u6,'(I5,3X,A8,4X,A8,8(I3,4X,I2,4X))') iSO_,dc(mdc)%LblCnt,ChTmp, &
                    (mc+iCo,iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
                     iCo=0,nIrrep/dc(mdc)%nStab-1)

                  if (iSO_ > nSOInf) then
                    write(u6,*) 'iSO_ > nSOInf'
                    call Abend()
                  end if
                  iSOInf(1,iSO_) = iCnttp
                  iSOInf(2,iSO_) = iCnt
                  iSOInf(3,iSO_) = iAng

                  if (Shells(iSh)%Aux .or. Shells(iSh)%Frag) cycle

                  if (.not. Primitive_Pass) then
                    write(isymunit,'(13(I4,4X))') iSO_,mdc,LVAL(lculf),MVAL(lculf),nIrrep/dc(mdc)%nStab, &
                      (iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
                       iCo=0,nIrrep/dc(mdc)%nStab-1)
                  end if
                  !                                                    *
                  !*****************************************************
                  !                                                    *
                  ! Stuff (not just) for LoProp

                  do iCo=0,nIrrep/dc(mdc)%nStab-1
                    ixxx = Index_NoSym(iCntrc,iComp,iAng,mdc,iCo,Index1,iCounter,iBas)
                    jxxx = Index_NoSym(iCntrc,iComp,iAng,mdc,iirrep,Index2,jCounter,iBas)
                    fact = real(iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))),kind=wp)

                    FacN = One/real(nIrrep/dc(mdc)%nStab,kind=wp)
                    if (MolWgh == 1) then
                      FacN = One
                    else if (MolWgh == 2) then
                      FacN = sqrt(FacN)
                    end if
                    SM(ixxx,iSO) = Fact*FacN
                    iyy = Index_Center(mdc,iCo,IndC,iAtoms,S%mCentr)

                    iCI(ixxx) = iyy
                    jCI(jxxx) = icnt

                    if (iCntrc <= List(iAng)) then
                      iOT(ixxx) = Occ
                    else
                      iOT(ixxx) = Vir
                    end if

                    iR = NrOpr(dc(mdc)%iCoSet(iCo,0))
                    LPC(1:3,iyy) = dbsc(iCnttp)%Coor(1:3,iCnt)
                    if (btest(iOper(iR),0)) LPC(1,iyy) = -LPC(1,iyy)
                    if (btest(iOper(iR),1)) LPC(2,iyy) = -LPC(2,iyy)
                    if (btest(iOper(iR),2)) LPC(3,iyy) = -LPC(3,iyy)

                    LPQ(iyy) = dbsc(iCnttp)%Charge
                    LPMM(iyy) = dbsc(iCnttp)%IsMM
                    LPA(iyy) = dbsc(iCnttp)%AtmNr

                    LP_Names(iyy) = dc(mdc)%LblCnt(1:LenIn)//':'//ChOper(iOper(iR))
                    desym_basis_ids(1,ixxx) = iyy
                    desym_basis_ids(2,ixxx) = iCntrc
                    desym_basis_ids(3,ixxx) = llab
                    desym_basis_ids(4,ixxx) = mlab
                  end do
                  !                                                    *
                  !*****************************************************
                  !                                                    *
                  Mamn(iSO) = dc(mdc)%LblCnt(1:LenIn)//ChTemp(1:8)
                  basis_ids(1,iSO) = mdc
                  basis_ids(2,iSO) = iCntrc
                  basis_ids(3,iSO) = llab
                  basis_ids(4,iSO) = mlab
                  fermion_type(iSO) = 0
                  if (dbsc(iCnttp)%fMass /= One) fermion_type(iSO) = 1
                  if (.not. Primitive_Pass) then
                    kIrrep = kIrrep+1
                    icent(kIrrep) = mdc
                    lnang(kIrrep) = lval(lculf)
                    lmag(kIrrep) = mval(lculf)
                    lant(kIrrep) = nIrrep/dc(mdc)%nStab
                  end if
                end do

              end do
            end if
            kComp = kComp+(iAng+1)*(iAng+2)/2
            kculf = kculf+2*iAng+1
          end do ! iAng
          mc = mc+nIrrep/dc(mdc)%nStab
        end do ! iCnt

      end do ! jCnttp
    end if
    !ulf
    nrSym = nIrrep
    nrBas(iIrrep+1) = nBas(iIrrep)
    !write(u6,*) ' nBas(iIrrep)', iIrrep, nBas(iIrrep)
    S%nDim = S%nDim+nBas(iIrrep)
    S%n2Tot = S%n2Tot+nBas(iIrrep)**2
  end do ! iIrrep
  !if (lSkip) S%nDim = iBas
  if ((iBas /= iSO) .and. (iBas_Aux /= iSO_Aux-iSO) .and. (.not. lSkip)) then
    write(u6,*) 'iBas=',iBas
    write(u6,*) 'iBas_Aux=',iBas_Aux
    write(u6,*) 'iSO=',iSO
    write(u6,*) 'iSO_Aux=',iSO_Aux-iSO
    write(u6,*) 'iSO_Tot=',iSO_Tot
    call Abend()
  end if
  ! redefine iOffSO array in case of Fragment AIEMP
  if (lFAIEMP) then
    do iIrrep=0,nIrrep-1
      iOffSO(iIrrep) = jOffSO(iIrrep)
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('Symmetrization Matrix','(20F5.2)',SM,iBas,iBas)
# endif
  call Put_dArray('SM',SM,iBas**2)

  !SVC: basis IDs of both symmetric and non-symmetric case
  if (.not. Primitive_Pass) then
    call Put_iArray('Fermion IDs',fermion_type,iSO)
    call Put_iArray('Basis IDs',basis_ids,4*iSO)
    call Put_iArray('Desym Basis IDs',desym_basis_ids,4*iBas)
  end if

  call mma_deallocate(IndC)
  call mma_deallocate(SM)
  call mma_deallocate(Index1)
  call mma_deallocate(Index2)

else

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! No symmetry case.

  if (Show .and. (iPrint >= 6)) then
    write(u6,*)
    write(u6,'(19x,a)') ' **************************************************'
    write(u6,'(19x,a)') ' ********** Petite list Basis Functions ***********'
    write(u6,'(19x,a)') ' **************************************************'
    write(u6,*)
  end if

  kIrrep = 0
  do iIrrep=0,nIrrep-1
    iOffSO(iIrrep) = iSO_Tot
    iAO = 0
    jSO = 0
    nBas(iIrrep) = 0
    nBas_Aux(iIrrep) = 0
    nBas_Frag(iIrrep) = 0
    bType(iIrrep) = .true.

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
        nCore_Sh(0:lmax) = 0
      else
        List(:) = List_AE
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
          LP_Names(mdc) = dc(mdc)%LblCnt(1:LenIn)//'    '
        end if
        do iAng=0,dbsc(iCnttp)%nVal-1
          nCore = nCore_Sh(iAng)
          iSh = iSh+1
          iShell = iShell+1
          nExpi = Shells(iSh)%nExp
          nBasisi = Shells(iSh)%nBasis
          if ((nExpi /= 0) .and. (nBasisi /= 0)) then
            if (Shells(iSh)%Prjct) then
              jComp = 2*iAng+1
            else
              jComp = (iAng+1)*(iAng+2)/2
            end if
            do iComp=1,jComp
              iAO = iAO+1
              if (iAO > MxAO) then
                write(u6,*) ' Increase MxAO'
                call Abend()
              end if
              lComp = kComp+iComp
              lculf = kculf+iComp

              ! Skip if symmetry operator is not in the coset of this center.

              do imc=0,(nIrrep/dc(mdc)%nStab)-1
                if (dc(mdc)%iCoSet(imc,0) /= iOper(iIrrep)) cycle
                if (output .and. bType(iIrrep)) then
                  write(u6,*)
                  write(u6,'(10X,2A)') ' Basis functions generated by ',ChOper(iIrrep)
                  write(u6,*)
                  write(u6,'(A)') ' Basis Label        Type   Center'
                  bType(iIrrep) = .false.
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
                    write(u6,*) ' iSO_ > nMamn'
                    write(u6,*) 'nMamn=',nMamn
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

                  if (output) write(u6,'(I5,2X,A8,5X,A8,I3)') iSO_,dc(mdc)%LblCnt,ChTmp,mc+imc

                  if (iSO_ > nSOInf) then
                    write(u6,*) 'iSO_ > nSOInf'
                    call Abend()
                  end if
                  iSOInf(1,iSO_) = iCnttp
                  iSOInf(2,iSO_) = iCnt
                  iSOInf(3,iSO_) = iAng

                  if (Shells(iSh)%Aux .or. Shells(iSh)%Frag) cycle
                  write(isymunit,'(13(I4,4X))') iSO,mdc,LVAL(lculf),MVAL(lculf),nIrrep/dc(mdc)%nStab, &
                    (iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
                     iCo=0,nIrrep/dc(mdc)%nStab-1)
                  !                                                    *
                  !*****************************************************
                  !                                                    *
                  ! Stuff (not just) for LoProp

                  iCI(iSO) = mdc
                  jCI(iSO) = mdc
                  if (iCntrc <= List(iAng)) then
                    iOT(iSO) = Occ
                  else
                    iOT(iSO) = Vir
                  end if
                  LPC(1:3,mdc) = dbsc(iCnttp)%Coor(1:3,iCnt)
                  LPQ(mdc) = dbsc(iCnttp)%Charge
                  LPMM(mdc) = dbsc(iCnttp)%IsMM
                  LPA(mdc) = dbsc(iCnttp)%AtmNr
                  LP_Names(mdc) = dc(mdc)%LblCnt(1:LenIn)//'    '
                  !                                                    *
                  !*****************************************************
                  !                                                    *
                  Mamn(iSO) = dc(mdc)%LblCnt(1:LenIn)//ChTemp(1:8)
                  basis_ids(1,iSO) = mdc
                  basis_ids(2,iSO) = iCntrc
                  basis_ids(3,iSO) = llab
                  basis_ids(4,iSO) = mlab
                  fermion_type(iSO) = 0
                  if (dbsc(iCnttp)%fMass /= One) fermion_type(iSO) = 1
                  if (.not. Primitive_Pass) then
                    kIrrep = kIrrep+1
                    icent(kIrrep) = mdc
                    lnang(kIrrep) = lval(lculf)
                    lmag(kIrrep) = mval(lculf)
                    lant(kIrrep) = nIrrep/dc(mdc)%nStab
                  end if
                end do
              end do
            end do
          end if
          kComp = kComp+(iAng+1)*(iAng+2)/2
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

      iSh = dbsc(iCnttp)%iVal-1
      do iAng=0,dbsc(iCnttp)%nVal-1
        iSh = iSh+1
        iShell = iShell+1
        nExpi = Shells(iSh)%nExp
        nBasisi = Shells(iSh)%nBasis
        if ((nExpi /= 0) .and. (nBasisi /= 0)) then
          if (Shells(iSh)%Prjct) then
            jComp = 2*iAng+1
          else
            jComp = (iAng+1)*(iAng+2)/2
          end if
          do iComp=1,jComp
            iAO = iAO+1
            if (iAO > MxAO) then
              write(u6,*) ' Increase MxAO'
              call Abend()
            end if
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
        end if
      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do

  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  if (iPrint >= 6) then ! std print has been moved to seward.f

    ! Print out basis set information

    Frmt = '(6X,A,T30,8I4)'
    write(u6,*)
    write(u6,'(6X,A)') 'Basis set specifications :'
    write(u6,'(6X,A,T30,8(1X,A))') 'Symmetry species',(lIrrep(i),i=0,nIrrep-1)
    write(u6,Frmt) 'Basis functions',(nBas(i),i=0,nIrrep-1)

  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write info (not just) for LoProp

if (.not. Primitive_Pass) then
  call Put_cArray('LP_L',LP_Names(1),(LenIn4)*S%mCentr)
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
call mma_deallocate(LP_Names)
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
call mma_deallocate(List)
call mma_deallocate(List_AE)
call mma_deallocate(LVAL)
call mma_deallocate(MVAL)
call mma_deallocate(nCore_Sh)

if (Show .and. (iPrint >= 6)) then
  call CollapseOutput(0,'   SO/AO info:')
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) ' *** iAOtSO ***'
do jAO=1,iAO
  write(u6,*) (iAOtSO(jAO,jIrrep),jIrrep=0,nIrrep-1)
end do
#endif

write(isymunit,'(A3)') 'END'
close(isymunit)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine SOCtl_Seward
