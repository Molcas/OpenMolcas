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
! Copyright (C) 2011, Daoling Peng                                     *
!               2021, Rulin Feng                                       *
!***********************************************************************

subroutine DKRelint_DP()
! modified by D. Peng, ETH Zurich, October 2011
!
! Interface/Driver routine for scalar relativistic
!       arbitrary-order DKH method,
!       exact decoupling X2C method &
!       exact decoupling BSS method.

use Basis_Info, only: dbsc, nBas, ncnttp
use DKH_Info, only: cLightAU, iRelae, LDKroll, radiLD
use Symmetry_Info, only: nIrrep
use Gateway_Info, only: lMXTC
use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: dkhorder, dkhparam, i, i_Dim, iAngr, iBas, icnt, iCnttp, iComp, idbg, idum(1), iExp, iibas, iMEF, iOpt, &
                     iPrint, iProps, iRC, iRout, iSize, iSizec, iSizep, iSizes, iTemp, jCent, jExp, k, kAng, kC, kCof, kCofi, &
                     kCofj, kExp, kExpi, kExpj, ks, kz, L, lOper, lOper_save, Lu_One, n, n_Int, nAtoms, nBas_cont(8), &
                     nBas_prim(8), nbl, nComp, nSym, numb_props, Op, relmethod, xorder
real(kind=wp) :: rCofi, rCofj, rEpsilon, rExpi, rExpj, rI, rNorm, rSum, VELIT
logical(kind=iwp) :: DoFullLT
!character(len=3) :: paramtype
character(len=8) :: Label, pXpLbl
integer(kind=iwp), allocatable :: indx(:), Loc(:), Map(:), Ind(:,:)
real(kind=wp), allocatable :: iK(:), SS(:), V(:), pVp(:), K_Save(:), K_Done(:), U_L(:), U_S(:), X(:), pXp(:), Prop(:), Pmag(:), &
                              Y(:), P(:), G(:), Ev2(:,:), Eig(:,:), Sinv(:,:), Ew(:), E(:), Aa(:), Rr(:), Tt(:), Re1r(:,:), &
                              Auxi(:,:), Tmp(:), Even1(:,:), Pvpt(:), Bu(:), H(:), H_nr(:), H_temp(:)
logical(kind=iwp), parameter :: Debug = .false.
integer(kind=iwp), external :: nProp_Int
#include "Molcas.fh"
#include "rinfo.fh"
#include "print.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 77
iPrint = nPrint(iRout)

if (Debug) then
  idbg = 6
else
  idbg = -1
end if

!                                                                      *
!***********************************************************************
!                                                                      *
! Save basis set info from contracted run

if (iprint >= 10) write(u6,*) ' In DKRelInt',ncnttp
kCof = 0
kAng = 0
kExp = 0
kC = 0

! Normalize coefficients

do iCnttp=1,nCnttp

  ! The non-valence type shells come at the end.
  ! When this block is encountered stop the procedure.

  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nFragType > 0)) exit

  do icnt=1,dbsc(iCnttp)%nCntr
    kC = kC+1
    do iAngr=0,nAngr(kC)
      rI = real(iAngr,kind=wp)+One+Half
      kAng = kAng+1
      do iBas=1,nBasisr(kAng)
        rSum = Zero
        kExpi = kExp
        kCofi = kCof
        do iExp=1,nPrimr(kAng)
          kExpi = kExpi+1
          kCofi = kCofi+1
          rExpi = rExp(kExpi)
          !write(u6,'(a11,f20.8)') ' Exponents',rExpi
          rCofi = rCof(kCofi)
          kExpj = kExp
          kCofj = kCof
          do jExp=1,nPrimr(kAng)
            kExpj = kExpj+1
            kCofj = kCofj+1
            rExpj = rExp(kExpj)
            rCofj = rCof(kCofj)
            rSum = rSum+rCofi*rCofj*(Two*sqrt(rExpi*rExpj)/(rExpi+rExpj))**rI
          end do
        end do
        rNorm = One/sqrt(rSum)
        if (iprint >= 10) write(u6,*) ' rNorm',kAng,rNorm
        do iExp=1,nPrimr(kAng)
          rCof(kCof+iExp) = rCof(kCof+iExp)*rNorm
          if (iprint >= 10) then
            write(u6,'(a24,f20.6)') ' normalized coefficients',rCof(kCof+iExp)
          end if
        end do
        kCof = kCof+nPrimr(kAng)
      end do
      kExp = kExp+nPrimr(kAng)
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 10) then
  i = 0
  do L=1,nrSym
    write(u6,*) ' Irreducible representation',L
    do ibas=1,nrBas(L)
      i = i+1
      write(u6,'(20i4)') i,icent(i),lnang(i),lmag(i)
    end do
  end do
end if

nBas_Cont(:) = nBas
nSym = nIrrep
!                                                                      *
!***********************************************************************
!                                                                      *
! Close ONEINT

iOpt = 0
call ClsOne(iRC,iOpt)
!                                                                      *
!***********************************************************************
!                                                                      *
! HFC magnetic integrals

if (lMXTC) then
  call Get_nAtoms_All(nAtoms)
  call copy_mag_ints(nAtoms)
end if

! Open ONEREL

iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
if (iRC /= 0) call Error()

call OneBas('PRIM')
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iArray('nBas_Prim',nBas,nSym)
nBas_prim(:) = nBas
if (iPrint >= 10) then
  write(u6,'(a,8i5)') ' Symmetries          ',nSym
  write(u6,'(a,8i5)') ' Primitive basis fcns',(nBas(i),i=0,nSym-1)
end if

! Allocate memory for relativistic part

VELIT = cLightAU
iSizep = 0
iSizes = 0
iSizec = 0
iibas = 0
do L=0,nSym-1
  iSizep = iSizep+nBas(L)*(nBas(L)+1)/2
  iSizes = iSizes+nBas(L)*nBas(L)
  iSizec = iSizec+nrBas(L+1)*(nrBas(L+1)+1)/2
  iibas = iibas+nbas(L)
end do
if (iPrint >= 10) write(u6,*) ' iSizep',iSizep

call mma_allocate(iK,iSizep+4,label='Kin')
call mma_allocate(SS,iSizep+4,label='SS')
call mma_allocate(V,iSizep+4,label='V')
call mma_allocate(pVp,iSizep+4,label='pVp')

Label = 'Mltpl  0'
iComp = 1
iOpt = 0
iRC = -1
call RdOne(iRC,iOpt,Label,iComp,SS,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
nComp = 1
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,[1],SS)
Label = 'Attract '
iRC = -1
call RdOne(iRC,iOpt,Label,iComp,V,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,[1],V)
Label = 'Kinetic '
iRC = -1
call RdOne(iRC,iOpt,Label,iComp,iK,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,[1],iK)
Label = 'pVp     '
iRC = -1
call RdOne(iRC,iOpt,Label,iComp,pVp,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,[1],pVp)

iOpt = 0
call ClsOne(iRC,iOpt)
!                                                                      *
!***********************************************************************
!                                                                      *
if (IRELAE >= 100) then
  if (IRELAE >= 1000) then

    ! IRELAE codes the method (DKH=1), the order (01 to 99), and the
    ! parametrization (1 to 5), i.e. 3rd order DKH with sqrt
    ! parametrization: 1-03-3, i.e. IRELAE=1033

    relmethod = 1

    ! relmethod = 1  arbitrary order DKH
    !           = 2  exact decoupling X2C
    !           = 3  exact decoupling BSS

    xOrder = iRELAE/10000
    iTemp = iRELAE-1000-xOrder*10000
    dkhorder = iTemp/10
    iTemp = iTemp-dkhorder*10
    dkhparam = iTemp
    if (xorder > dkhorder) then
      write(u6,'(a,i3,a,i3)') ' xorder was reduced from ',xorder,' to ',dkhorder
      xorder = dkhorder
    end if
    !if (iTemp == 1) then
    !  paramtype = 'OPT'
    !else if (iTemp == 2) then
    !  paramtype = 'EXP'
    !else if (iTemp == 3) then
    !  paramtype = 'SQR'
    !else if (iTemp == 4) then
    !  paramtype = 'MCW'
    !else if (iTemp == 5) then
    !  paramtype = 'CAY'
    !else
    if ((iTemp < 1) .or. (iTemp > 5)) then
      write(u6,*) 'dkrelint: Illegal parametrization!'
      call Abend()
    end if
    !write(u6,'(A)') ' Computing the arbitrary order DKH Hamiltonian'
  else if (IRELAE == 101) then
    relmethod = 2
    xorder = 1
    !write(u6,'(A)') ' Computing the exact decoupling X2C Hamiltonian'
  else if (IRELAE == 102) then
    relmethod = 3
    xorder = 1
    !write(u6,'(A)') ' Computing the exact decoupling BSS Hamiltonian'
  else
    write(u6,*) 'dkrelint: Unknown method !'
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_allocate(K_Save,iSizep+4,label='K_Save')
  call mma_allocate(K_Done,iSizep+4,label='K_Done')
  call dcopy_(iSizep+4,iK,1,K_Save,1)

  call mma_allocate(U_L,iSizes+4,label='U_L')
  call mma_allocate(U_S,iSizes+4,label='U_S')

  ! Read block information if do local transformation

  if (LDKroll) then
    call mma_allocate(indx,iibas,label='Index')
    call xdr_indx(iibas,indx)
    DoFullLT = .true.
    if (radiLD == Zero) DoFullLT = .false.
    if (DoFullLT) then
      if ((relmethod == 1) .and. (xorder == 0)) then
        xorder = dkhorder
      end if
      write(u6,'(A)') '   DLU Local Transformation'
    else
      write(u6,'(A)') '   DLH Local Transformation'
    end if
  end if

  ! Do the Hamiltonian separately

  k = 1
  ks = 1
  kz = 1

  do L=0,nSym-1
    n = nBas(L)
    iSize = n*(n+1)/2
    if (iSize == 0) cycle
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (LDKroll) then
      call mma_allocate(Loc,n,label='InfoLoc')
      call mma_allocate(Map,n,label='InfoMap')
      call xdr_info_local(n,indx(kz),nbl,Loc,Map)
      !DP write(u6,'(a,i1,i5,a,99i4)') '   Sym: ',L+1,n,'  = Local ',(Loc(i),i=1,nbl)
      call XDR_Local_Ham(n,isize,n*n,relmethod,dkhparam,dkhorder,xorder,SS(k),iK(k),V(k),pVp(k),U_L(ks),U_S(ks),nbl,Loc,Map, &
                         DoFullLT,cLightAU)
      call mma_deallocate(Loc)
      call mma_deallocate(Map)
    else
      call XDR_Ham(n,isize,n*n,relmethod,dkhparam,dkhorder,xorder,SS(k),iK(k),V(k),pVp(k),U_L(ks),U_S(ks),cLightAU)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ks = ks+n*n
    kz = kz+n
    k = k+isize
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call dcopy_(iSizep+4,iK,1,K_Done,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (xOrder > 0) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Pick up the number of property integrals to process.

#   ifndef MOLPRO
    numb_props = nProp_Int(.false.,idum,0)
    call mma_allocate(Ind,4,numb_props,label='Ind')
    numb_props = nProp_Int(.true.,Ind,numb_props)
    do iProps=1,numb_props
      Op = Ind(1,iProps)
      iMEF = Ind(2,iProps)
      iComp = Ind(3,iProps)
      jCent = Ind(4,iProps)

      call dcopy_(iSizep+4,K_Save,1,iK,1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Read the property integrals
      !
      ! Open ONEREL

      iOpt = 0
      iRC = -1
      Lu_One = 2
      call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
      if (iRC /= 0) call Error()

      call OneBas('PRIM')

      if (Op == 1) then
        write(Label,'(a,i2)') 'MLTPL ',iMEF
      else if (Op == 2) then
        write(Label,'(a,i1,i5)') 'EF',iMEF,jCent
      else if (Op == 3) then
        write(Label,'(a,i5)') 'Cnt',jCent
      else if (Op == 4) then
        write(Label,'(A,I3)') 'MAGXP',jCent
      else
        write(u6,*) 'DKRelInt: illegal property!'
        call Abend()
      end if
      !write(u6,*)
      !write(u6,*) 'Label=',Label
      !write(u6,*) 'iComp=',iComp
      !write(u6,*)

      iOpt = ibset(0,sOpSiz)
      iRC = -1
      lOper = -1
      call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
      if (iRC == 0) n_Int = idum(1)
      !write(u6,*) 'lOper=',lOper
      call mma_allocate(X,n_Int+4,label='X')
      iRC = -1
      iOpt = 0
      call RdOne(iRC,iOpt,Label,iComp,X,lOper)
      if (iRC /= 0) then
        write(u6,*) 'DKRelInt: Error reading from ONEREL'
        write(u6,'(A,A)') 'Label=',Label
        write(u6,'(A,A)') 'iRC=',iRC
        call Abend()
      end if
      call CmpInt(X,n_Int,nBas_Prim,nSym,lOper)
      if (n_Int == 0) then
        iOpt = 0
        call ClsOne(iRC,iOpt)
      else

        if (Op == 1) then
          write(pXpLbl,'(A,I2)') 'pMp   ',iMEF
        else if (Op == 2) then
          write(pXpLbl,'(A,I1,I5)') 'PP',iMEF,jCent
        else if (Op == 3) then
          write(pXpLbl,'(A,I2)') 'pCp   ',jCent
        else if (Op == 4) then
          write(pXpLbl,'(A,I3)') 'MAGPX',jCent
        end if
        iOpt = ibset(0,sOpSiz)
        iRC = -1
        call iRdOne(iRC,iOpt,pXpLbl,iComp,idum,lOper)
        if (iRC == 0) n_Int = idum(1)
        call mma_allocate(pXp,n_Int+4,label='pXp')
        iOpt = 0
        iRC = -1
        call RdOne(iRC,iOpt,pXpLbl,iComp,pXp,lOper)
        if (iRC /= 0) then
          write(u6,*) 'DKRelInt: Error reading from ONEREL'
          write(u6,'(A,A)') 'pXpLbl=',pXpLbl
          write(u6,'(A,A)') 'iRC=',iRC
          call Abend()
        end if
        call CmpInt(pXp,n_Int,nBas_Prim,nSym,lOper)

        iOpt = 0
        call ClsOne(iRC,iOpt)

        !call mma_maxDBLE(Mem_Available)
        !write(u6,*) 'Mem_Available=',Mem_Available
        k = 1
        ks = 1
        do L=0,nSym-1
          n = nBas(L)
          iSize = n*(n+1)/2
          if (iSize == 0) cycle

          ! Skip if the propetry operator does not have a total
          ! symmetric component!

          if (btest(lOper,0)) then
            !                                                          *
            !***********************************************************
            !                                                          *
            call XDR_Prop(n,isize,n*n,relmethod,dkhparam,xorder,SS(k),iK(k),V(k),pVp(k),X(k),pXp(k),U_L(ks),U_S(ks),cLightAU, &
                          Label,iComp,iSizec)
            ks = ks+n*n
          end if
          k = k+isize
        end do
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Put the picture change corrected integral back to the
        ! ONEINT file. Primitives in X.
        !
        ! First contract the result, store in Prop

        call mma_allocate(Prop,iSizec+4,label='Prop')
        call repmat(idbg,X,Prop,.true.)
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Read the contracted property integrals from OneInt

        iOpt = 0
        iRC = -1
        Lu_One = 2
        call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
        if (iRC /= 0) call Error()

        if (Op == 4) then
          call mma_allocate(Pmag,iSizec+4,label='Pmag')
          call repmat(idbg,X,Pmag,.false.)
          lOper_save = lOper
          lOper = 255
          call WrOne(iRC,iOpt,Label,iComp,Pmag,lOper)
          call WrOne(iRC,iOpt,pXpLbl,iComp,Pmag,lOper)
          iOpt = 0
          call ClsOne(iRC,iOpt)
          call mma_deallocate(X)
          call mma_deallocate(pXp)
          call mma_deallocate(Prop)
          call mma_deallocate(Pmag)
          lOper = lOper_save
          cycle
        end if

        iOpt = ibset(0,sOpSiz)
        iRC = -1
        lOper = -1
        call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
        if (iRC == 0) n_Int = idum(1)
        call mma_allocate(Y,n_Int+4,label='Y')
        iRC = -1
        iOpt = 0
        call RdOne(iRC,iOpt,Label,iComp,Y,lOper)
        !write(u6,*) 'Y1=',DDot_(n_Int,Y,1,One,0)
        if (iRC /= 0) then
          write(u6,*) 'DKRelInt: Error reading from ONEINT'
          write(u6,'(A,A)') 'Label=',Label
          call Abend()
        end if

        ! Put the picture change corrected blocks in. Note that this
        ! is just the diagonal symmetry blocks.

        call Cp_Prop_Int(Y,n_Int,Prop,iSizec,nrBas,nIrrep,lOper)

        ! Now write it back to disc

        iOpt = 0
        call WrOne(iRC,iOpt,Label,iComp,Y,lOper)

        iOpt = 0
        call ClsOne(iRC,iOpt)

        call mma_deallocate(Prop)
        call mma_deallocate(Y)

        call mma_deallocate(pXp)
      end if
      call mma_deallocate(X)
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do ! iProps

    call mma_deallocate(Ind)
#   endif
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call dcopy_(iSizep+4,K_Done,1,iK,1)
  call mma_deallocate(K_Save)
  call mma_deallocate(K_Done)
  call mma_deallocate(U_L)
  call mma_deallocate(U_S)
  if (LDKroll) then
    call mma_deallocate(indx)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  !      The Old code, no option for property integrals.
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  !      Loop over the symmetry blocks
  !
  rEpsilon = 1.0e-10_wp
  k = 1
  do L=0,nSym-1
    n = nBas(L)
    iSize = n*(n+1)/2
    if (iSize == 0) cycle
    ! Allocate

    call mma_allocate(P,iSize+4,label='P')
    call mma_allocate(G,iSize+4,label='G')
    call mma_allocate(Ev2,n,n,label='Ev2')
    call mma_allocate(Eig,n,n,label='Eig')
    call mma_allocate(Sinv,n,n,label='Sinv')
    call mma_allocate(Ew,n,label='Ew')
    call mma_allocate(E,n,label='E')
    call mma_allocate(Aa,n,label='Aa')
    call mma_allocate(Rr,n,label='Rr')
    call mma_allocate(Tt,n,label='Tt')
    call mma_allocate(Re1r,n,n,label='Re1r')
    call mma_allocate(Auxi,n,n,label='Auxi')
    call mma_allocate(Tmp,iSize,label='Tmp')
    if (IRELAE == 0) then
      call mma_allocate(Even1,1,1,label='Even1')
      call mma_allocate(Pvpt,1,label='Pvpt')
      call mma_allocate(Bu,1,label='Bu')
      i_Dim = 1
    else
      call mma_allocate(Even1,n,n,label='Even1')
      call mma_allocate(Pvpt,iSize,label='Pvpt')
      call mma_allocate(Bu,iSize,label='Bu')
      i_Dim = n
    end if

    ! call to package relsew

    call SCFCLI(idbg,rEpsilon,SS(k),iK(k),V(k),pVp(k),n,iSize,VELIT,Bu,P,G,Ev2,Eig,Sinv,Ew,E,Aa,Rr,Tt,Pvpt,Even1,Re1r,Auxi,Tmp, &
                i_Dim)

    call mma_deallocate(P)
    call mma_deallocate(G)
    call mma_deallocate(Ev2)
    call mma_deallocate(Eig)
    call mma_deallocate(Sinv)
    call mma_deallocate(Ew)
    call mma_deallocate(E)
    call mma_deallocate(Aa)
    call mma_deallocate(Rr)
    call mma_deallocate(Tt)
    call mma_deallocate(Re1r)
    call mma_deallocate(Auxi)
    call mma_deallocate(Tmp)
    call mma_deallocate(Even1)
    call mma_deallocate(Pvpt)
    call mma_deallocate(Bu)
    k = k+isize
  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate arrays in contracted basis

call mma_allocate(H,iSizec+4,label='H')
call mma_allocate(H_nr,iSizec+4,label='H_o')
call mma_allocate(H_temp,iSizec+4,label='H_temp')
H(:) = Zero
H_nr(:) = Zero
H_temp(:) = Zero

#ifdef MOLPRO
! store relativistic H
call repmat(idbg,ik,H_temp,.true.)
call fperm(H_temp,H)
call writem(H,iSizec+2,1,1200,0,'H0')
call writem(H,iSizec+2,1,1210,0,'H01')
! store V=H-T
call lesw(SS,iSizec,1,1400,0)
call daxpy_(iSizec,-One,SS,1,H,1)
call writem(H,iSizec+2,1,1410,0,'POT')
! reset contracted basis size
nBas(:) = nBas_Cont
call mma_deallocate(iK)
call mma_deallocate(SS)
call mma_deallocate(V)
#else

! Note: in combination with ECPs V is only based on the effective
!       charges of the atoms. In the primitive basis, however, we
!       have temporarily introduced the actual atomic charges. We
!       have to fix this now. Hence the somewhat strange way in
!       which the DKH corrected Hamiltonian is computed.
!
! Compute stripped non-relativistic H H_temp
!                                                                      *
!***********************************************************************
!                                                                      *
! Open ONEREL

iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
if (iRC /= 0) call Error()

call OneBas('PRIM')

Label = 'Kinetic '
iComp = 1
call RdOne(iRC,iOpt,Label,iComp,SS,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'Attract '
call RdOne(iRC,iOpt,Label,iComp,V,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
! Add Kinetic and Attraction term
call DaXpY_(iSizep+4,One,SS,1,V,1)
if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Prim,1)
  call PrMtrx('Attract+Kinetic (prim)',[lOper],nComp,[1],V)
  call PrMtrx('Kinetic (prim)',[lOper],nComp,[1],SS)
  call iSwap(8,nBas,1,nBas_Prim,1)
end if
H_temp(iSizec+1:) = Zero
! Contract and store in H_temp
call repmat(idbg,V,H_temp,.true.)

call mma_deallocate(SS)
call mma_deallocate(V)

! Close ONEREL and re-open ONEINT

iOpt = 0
iRC = -1
call ClsOne(iRC,iOpt)
if (iRC /= 0) call Error()
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
if (iRC /= 0) call Error()

if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Cont,1)
  call PrMtrx('H_temp (cont)',[lOper],nComp,[1],H_temp)
  call iSwap(8,nBas,1,nBas_Cont,1)
end if

! Transform DKH Hamiltonian to contracted basis (H)

if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Prim,1)
  call PrMtrx('iK (prim)',[lOper],nComp,[1],iK)
  call iSwap(8,nBas,1,nBas_Prim,1)
end if
call repmat(idbg,iK,H,.true.)
if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Cont,1)
  call PrMtrx('H (cont)',[lOper],nComp,[1],H)
  call iSwap(8,nBas,1,nBas_Cont,1)
end if

call mma_deallocate(iK)

iOpt = 0
iRC = -1
Label = 'OneHam 0'
call RdOne(iRC,iOpt,Label,iComp,H_nr,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iOpt = 0
iRC = -1

! final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
! where (s) is stripped and with full charge

call DaXpY_(iSizec+4,-One,H_temp,1,H,1)
call DaXpY_(iSizec+4,One,H_nr,1,H,1)

call Get_iArray('nBas',nBas,nSym)
if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
end if
Label = 'OneHam 0'
lOper = 1
nComp = 1
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,[1],H)

! Replace 1-el Hamiltonian on ONEINT

iRC = -1
call WrOne(iRC,iOpt,Label,1,H,lOper)
Label = 'OneHam  '
call WrOne(iRC,iOpt,Label,1,H,lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

if (IRELAE == 23) then   ! IORA

  ! Replace overlap on ONEINT

  call repmat(idbg,pVp,H,.true.)
  iRC = -1
  Label = 'Mltpl  0'
  call WrOne(iRC,iOpt,Label,1,H,lOper)
  if (iRC /= 0) then
    write(u6,*) 'DKInt: Error reading from ONEINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
end if
#endif

call mma_deallocate(pVp)
call mma_deallocate(H)
call mma_deallocate(H_nr)
call mma_deallocate(H_temp)

return

contains

subroutine Error()
  write(u6,*) ' *** Error in subroutine DKRelInt ***'
  write(u6,*) '     Abend in subroutine OpnOne'
  call Abend()
end subroutine Error

end subroutine DKRelint_DP
