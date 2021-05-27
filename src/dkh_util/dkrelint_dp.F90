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
use DKH_Info, only: LDKroll, radiLD
use Symmetry_Info, only: nIrrep
use Logical_Info, only: lMXTC
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: dkhorder, dkhparam, i, i_Dim, iAa, iAngr, iAuxi, iBas, iBu, icnt, iCnttp, iComp, idbg, idum(1), iE, iEig, &
                     iEv2, iEven1, iEw, iExp, iG, iH, iH_nr, iH_temp, iibas, iK, iK_Done, iK_Save, iLoc, iMap, iMEF, indx, iOpt, &
                     iP, ip_MEF, ip_Pmag, ip_Prop, ipaddr(3), ipiComp, ipInd, ipjCent, ipOp, iPrint, iProps, ipVp, iPvpt, ipXp, &
                     iRC, iRe1r, iRout, iRr, iSinv, iSize, iSizec, iSizep, iSizes, iSS, iTemp, iTt, iTwrk4, iU_L, iU_S, iV, iX, &
                     iY, jCent, jExp, k, kAng, kC, kCof, kCofi, kCofj, kExp, kExpi, kExpj, ks, kz, L, Length, Length2, lOper, &
                     lOper_save, Lu_One, Mem_Available, n, n_Int, nAtoms, nBas_cont(8), nBas_prim(8), nbl, nComp, nrSym, nSym, &
                     numb_props, relmethod, xorder
real(kind=wp) :: rCofi, rCofj, rEpsilon, rExpi, rExpj, rI, rNorm, rSum, VELIT
logical(kind=iwp) :: DoFullLT
!character(len=3) :: paramtype
character(len=8) :: Label, pXpLbl
logical(kind=iwp), parameter :: Debug = .false.
integer(kind=iwp), external :: nProp_Int
#include "Molcas.fh"
#include "rinfo.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "RelLight.fh"
#include "relae.fh"

if (IRFLAG1 == 1) then
  call DKRelint()
  return
end if
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

  ! The none valence type shells comes at the end. When this block
  ! is encountered stop the procedure.

  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%nFragType > 0) Go To 999

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
999 continue
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

call iCopy(8,nBas,1,nBas_Cont,1)
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
if (iRC /= 0) Go To 9999

call OneBas('PRIM')
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iArray('nBas_Prim',nBas,nSym)
call iCopy(8,nBas,1,nBas_prim,1)
if (iPrint >= 10) then
  write(u6,'(a,8i5)') ' Symmetries          ',nSym
  write(u6,'(a,8i5)') ' Primitive basis fcns',(nBas(i),i=0,nSym-1)
end if

! Allocate memory for relativistic part

VELIT = CLightAU
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

call GetMem('Kin     ','ALLO','REAL',iK,iSizep+4)
call GetMem('SS      ','ALLO','REAL',iSS,iSizep+4)
call GetMem('V       ','ALLO','REAL',iV,iSizep+4)
call GetMem('pVp     ','ALLO','REAL',ipVp,iSizep+4)

if (iprint >= 20) write(u6,*) '  indices',iss,ik,iv,ipvp
Label = 'Mltpl  0'
iComp = 1
iOpt = 0
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
nComp = 1
ipaddr(1) = iSS
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
Label = 'Attract '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
ipaddr(1) = iV
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
Label = 'Kinetic '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(iK),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
ipaddr(1) = iK
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
Label = 'pVp     '
iRC = -1
call RdOne(iRC,iOpt,Label,1,Work(ipVp),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
ipaddr(1) = ipVp
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)

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
    !   paramtype='OPT'
    !else if (iTemp == 2) then
    !   paramtype='EXP'
    !else if (iTemp == 3) then
    !   paramtype='SQR'
    !else if (iTemp == 4) then
    !   paramtype='MCW'
    !else if (iTemp == 5) then
    !   paramtype='CAY'
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
  call Allocate_Work(iK_Save,iSizep+4)
  call Allocate_Work(iK_Done,iSizep+4)
  call dcopy_(iSizep+4,Work(iK),1,Work(iK_Save),1)

  call Allocate_Work(iU_L,iSizes+4)
  call Allocate_Work(iU_S,iSizes+4)

  ! Read block information if do local transformation

  if (LDKroll) then
    call GetMem('Index  ','ALLO','INTE',indx,iibas+4)
    call xdr_indx(iibas,iWork(indx))
    DoFullLT = .true.
    if (radiLD == Zero) DoFullLT = .false.
    if (DoFullLT) then
      if (relmethod == 1 .and. xorder == 0) then
        xorder = dkhorder
      end if
      write(u6,'(A)') '   DLU Local Transformation'
    else
      write(u6,'(A)') '   DLH Local Transformation'
    end if
  end if

  ! Do the Hamiltonian separately

  k = 0
  ks = 0
  kz = 0

  do L=0,nSym-1
    n = nBas(L)
    iSize = n*(n+1)/2
    if (iSize == 0) Go To 911
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (LDKroll) then
      call GetMem('InfoLoc','ALLO','INTE',iLoc,n+4)
      call GetMem('MapLoc ','ALLO','INTE',iMap,n+4)
      call xdr_info_local(n,iWork(indx+kz),nbl,iWork(iLoc),iWork(iMap))
      !DP write(u6,'(a,i1,i5,a,99i4)') '   Sym: ',L+1,n,'  = Local ',(iWork(iLoc+i),i=0,nbl-1)
      call XDR_Local_Ham(n,isize,n*n,relmethod,dkhparam,dkhorder,xorder,Work(iSS+k),Work(iK+k),Work(iV+k),Work(ipVp+k), &
                         Work(iU_L+ks),Work(iU_S+ks),iWork(indx+kz),nbl,iWork(iLoc),iWork(iMap),DoFullLT,clightau)
      call GetMem('InfoLoc','FREE','INTE',iLoc,n+4)
      call GetMem('MapLoc ','FREE','INTE',iMap,n+4)
    else
      call XDR_Ham(n,isize,n*n,relmethod,dkhparam,dkhorder,xorder,Work(iSS+k),Work(iK+k),Work(iV+k),Work(ipVp+k),Work(iU_L+ks), &
                   Work(iU_S+ks),clightau)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ks = ks+n*n
    kz = kz+n
911 continue
    k = k+isize
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call dcopy_(iSizep+4,Work(iK),1,Work(iK_Done),1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (xOrder <= 0) Go To 912
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pick up the number of property integrals to process.

# ifndef MOLPRO
  numb_props = nProp_Int(.false.,iWork(ip_iDummy),0)
  call Allocate_iWork(ipInd,4*numb_props)
  numb_props = nProp_Int(.true.,iWork(ipInd),numb_props)
  do iProps=1,numb_props
    ipOp = ipInd+(iProps-1)*4
    ip_MEF = ipInd+(iProps-1)*4+1
    ipiComp = ipInd+(iProps-1)*4+2
    ipjCent = ipInd+(iProps-1)*4+3

    call dcopy_(iSizep+4,Work(iK_Save),1,Work(iK),1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the property integrals
    !
    ! Open ONEREL

    iOpt = 0
    iRC = -1
    Lu_One = 2
    call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
    if (iRC /= 0) Go To 9999

    call OneBas('PRIM')

    iMEF = iWork(ip_MEF)
    if (iWork(ipOp) == 1) then
      write(Label,'(a,i2)') 'MLTPL ',iMEF
    else if (iWork(ipOp) == 2) then
      jCent = iWork(ipjCent)
      write(Label,'(a,i1,i5)') 'EF',iMEF,jCent
    else if (iWork(ipOp) == 3) then
      jCent = iWork(ipjCent)
      write(Label,'(a,i5)') 'Cnt',jCent
    else if (iWork(ipOp) == 4) then
      jCent = iWork(ipjCent)
      write(Label,'(A,I3)') 'MAGXP',jCent
    else
      write(u6,*) 'DKRelInt: illegal property!'
      call Abend()
    end if
    iComp = iWork(ipiComp)
    !write(u6,*)
    !write(u6,*) 'Label=',Label
    !write(u6,*) 'iComp=',iComp
    !write(u6,*)

    iOpt = 1
    iRC = -1
    lOper = -1
    call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
    if (iRC == 0) n_Int = idum(1)
    !write(u6,*) 'lOper=',lOper
    call GetMem('X       ','ALLO','REAL',iX,n_Int+4)
    iRC = -1
    iOpt = 0
    call RdOne(iRC,iOpt,Label,iComp,Work(iX),lOper)
    if (iRC /= 0) then
      write(u6,*) 'DKRelInt: Error reading from ONEREL'
      write(u6,'(A,A)') 'Label=',Label
      write(u6,'(A,A)') 'iRC=',iRC
      call Abend()
    end if
    call CmpInt(Work(iX),n_Int,nBas_Prim,nSym,lOper)
    if (n_Int == 0) then
      iOpt = 0
      call ClsOne(iRC,iOpt)
      Go To 666
    end if

    if (iWork(ipOp) == 1) then
      write(pXpLbl,'(A,I2)') 'pMp   ',iMEF
    else if (iWork(ipOp) == 2) then
      write(pXpLbl,'(A,I1,I5)') 'PP',iMEF,jCent
    else if (iWork(ipOp) == 3) then
      write(pXpLbl,'(A,I2)') 'pCp   ',jCent
    else if (iWork(ipOp) == 4) then
      write(pXpLbl,'(A,I3)') 'MAGPX',jCent
    end if
    iOpt = 1
    iRC = -1
    call iRdOne(iRC,iOpt,pXpLbl,iComp,idum,lOper)
    if (iRC == 0) n_Int = idum(1)
    call GetMem('pXp     ','ALLO','REAL',ipXp,n_Int+4)
    iOpt = 0
    iRC = -1
    call RdOne(iRC,iOpt,pXpLbl,iComp,Work(ipXp),lOper)
    if (iRC /= 0) then
      write(u6,*) 'DKRelInt: Error reading from ONEREL'
      write(u6,'(A,A)') 'pXpLbl=',pXpLbl
      write(u6,'(A,A)') 'iRC=',iRC
      call Abend()
    end if
    call CmpInt(Work(ipXp),n_Int,nBas_Prim,nSym,lOper)

    iOpt = 0
    call ClsOne(iRC,iOpt)

    call GetMem('Core','Max','Real',iDum(1),Mem_Available)
    !write(u6,*) 'Mem_Available=',Mem_Available
    k = 0
    ks = 0
    kz = 0
    do L=0,nSym-1
      n = nBas(L)
      iSize = n*(n+1)/2
      if (iSize == 0) Go To 91

      ! Skip if the propetry operator does not have a total
      ! symmetric component!

      if (iand(1,lOper) == 0) Go To 91
      !                                                                *
      !*****************************************************************
      !                                                                *
      call XDR_Prop(n,isize,n*n,relmethod,dkhparam,dkhorder,xorder,Work(iSS+k),Work(iK+k),Work(iV+k),Work(ipVp+k),Work(iX+k), &
                    Work(ipXp+k),Work(iU_L+ks),Work(iU_S+ks),clightau,Label,iComp,iSizec)
      ks = ks+n*n
      kz = kz+n
91    continue
      k = k+isize
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Put the picture change corrected integral back to the
    ! ONEINT file. Primitives in Work(iX).
    !
    ! First contract the result, store in Work(ip_Prop)

    call Allocate_Work(ip_Prop,iSizec+4)
    call repmat(idbg,Work(iX),Work(ip_Prop),.true.)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the contracted property integrals from OneInt

    iOpt = 0
    iRC = -1
    Lu_One = 2
    call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
    if (iRC /= 0) Go To 9999

    if (iWork(ipOp) == 4) then
      call Allocate_Work(ip_Pmag,iSizec+4)
      call repmat(idbg,Work(iX),Work(ip_Pmag),.false.)
      lOper_save = lOper
      lOper = 255
      call WrOne(iRC,iOpt,Label,iComp,Work(ip_Pmag),lOper)
      call WrOne(iRC,iOpt,pXpLbl,iComp,Work(ip_Pmag),lOper)
      iOpt = 0
      call ClsOne(iRC,iOpt)
      call Free_Work(ip_Prop)
      call Free_Work(ip_Pmag)
      call GetMem('pXp     ','FREE','REAL',ipXp,iSizep+4)
      call GetMem('X       ','FREE','REAL',iX,iSizep+4)
      lOper = lOper_save
      cycle
    end if

    iOpt = 1
    iRC = -1
    lOper = -1
    call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
    if (iRC == 0) n_Int = idum(1)
    call GetMem('Y       ','ALLO','REAL',iY,n_Int+4)
    iRC = -1
    iOpt = 0
    call RdOne(iRC,iOpt,Label,iComp,Work(iY),lOper)
    !write(u6,*) 'Y1=',DDot_(n_Int,Work(iY),1,One,0)
    if (iRC /= 0) then
      write(u6,*) 'DKRelInt: Error reading from ONEINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if

    ! Put the picture change corrected blocks in. Note that this
    ! is just the diagonal symmetry blocks.

    call Cp_Prop_Int(Work(iY),n_Int,Work(ip_Prop),iSizec,nrBas,nIrrep,lOper)

    ! Now write it back to disc

    iOpt = 0
    call WrOne(iRC,iOpt,Label,iComp,Work(iY),lOper)

    iOpt = 0
    call ClsOne(iRC,iOpt)

    call Free_Work(iY)
    call Free_Work(ip_Prop)

    call GetMem('pXp     ','FREE','REAL',ipXp,iSizep+4)
666 continue
    call GetMem('X       ','FREE','REAL',iX,iSizep+4)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do ! iProps

  call Free_iWork(ipInd)
# endif
912 continue
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call dcopy_(iSizep+4,Work(iK_Done),1,Work(iK),1)
  call Free_Work(iK_Done)
  call Free_Work(iK_Save)
  call Free_Work(iU_L)
  call Free_Work(iU_S)
  if (LDKroll) then
    call GetMem('Index  ','FREE','INTE',indx,iibas+4)
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
  k = 0
  do L=0,nSym-1
    n = nBas(L)
    iSize = n*(n+1)/2
    if (iSize == 0) goto 9
    ! Allocate

    call GetMem('P       ','ALLO','REAL',iP,isize+4)
    call GetMem('G       ','ALLO','REAL',iG,isize+4)
    call GetMem('Ev2     ','ALLO','REAL',iEv2,n*n+4)
    call GetMem('Eig     ','ALLO','REAL',iEig,n*n+4)
    call GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
    call GetMem('Ew      ','ALLO','REAL',iEw,n+4)
    call GetMem('E       ','ALLO','REAL',iE,n+4)
    call GetMem('Aa      ','ALLO','REAL',iAa,n+4)
    call GetMem('Rr      ','ALLO','REAL',iRr,n+4)
    call GetMem('Tt      ','ALLO','REAL',iTt,n+4)
    call GetMem('Re1r    ','ALLO','REAL',iRe1r,n*n+4)
    call GetMem('Auxi    ','ALLO','REAL',iAuxi,n*n+4)
    call GetMem('Twrk4   ','ALLO','REAL',iTwrk4,n*200+4)
    Length = N*N+4
    Length2 = iSize+4
    i_Dim = N
    if (IRELAE == 0) then
      Length = 1
      Length2 = 1
      i_Dim = 1
    end if
    call GetMem('Even1   ','ALLO','REAL',iEven1,Length)
    call GetMem('Pvpt    ','ALLO','REAL',iPvpt,Length2)
    call GetMem('Bu      ','ALLO','REAL',iBu,Length2)

    ! call to package relsew

    call SCFCLI(idbg,rEpsilon,Work(iSS+k),Work(iK+k),Work(iV+k),Work(ipVp+k),n,iSize,VELIT,Work(iBu),Work(iP),Work(iG),Work(iEv2), &
                Work(iEig),Work(iSinv),Work(iEw),Work(iE),Work(iAa),Work(iRr),Work(iTt),Work(iPvpt),Work(iEven1),Work(iRe1r), &
                Work(iAuxi),Work(iTwrk4),i_Dim)

    call GetMem('Bu      ','FREE','REAL',iBu,Length2)
    call GetMem('P       ','FREE','REAL',iP,isize+4)
    call GetMem('G       ','FREE','REAL',iG,isize+4)
    call GetMem('Ev2     ','FREE','REAL',iEv2,n*n+4)
    call GetMem('Eig     ','FREE','REAL',iEig,n*n+4)
    call GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
    call GetMem('Ew      ','FREE','REAL',iEw,n+4)
    call GetMem('E       ','FREE','REAL',iE,n+4)
    call GetMem('Aa      ','FREE','REAL',iAa,n+4)
    call GetMem('Rr      ','FREE','REAL',iRr,n+4)
    call GetMem('Tt      ','FREE','REAL',iTt,n+4)
    call GetMem('Pvpt    ','FREE','REAL',iPvpt,Length2)
    call GetMem('Even1   ','FREE','REAL',iEven1,Length)
    call GetMem('Re1r    ','FREE','REAL',iRe1r,n*n+4)
    call GetMem('Auxi    ','FREE','REAL',iAuxi,n*n+4)
    call GetMem('Twrk4   ','FREE','REAL',iTwrk4,n*200+4)
9   continue
    k = k+isize
  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate arrays in contracted basis

call GetMem('H       ','ALLO','REAL',iH,iSizec+4)
call FZero(Work(iH),iSizec+4)
call GetMem('H_o     ','ALLO','REAL',iH_nr,iSizec+4)
call FZero(Work(iH_nr),iSizec+4)
call GetMem('H_temp  ','ALLO','REAL',iH_temp,iSizec+4)
call FZero(Work(iH_temp),iSizec+4)

#ifdef MOLPRO
! store relativistic H
call repmat(idbg,Work(ik),Work(iH_temp),.true.)
call fperm(Work(iH_temp),Work(ih))
call writem(Work(iH),isizec+2,1,1200,0,'H0')
call writem(Work(iH),isizec+2,1,1210,0,'H01')
! store V=H-T
call lesw(Work(iss),iSizec,1,1400,0)
call daxpy_(iSizec,-one,Work(iss),1,Work(iH),1)
call writem(Work(iH),isizec+2,1,1410,0,'POT')
! reset contracted basis size
call iCopy(8,nBas_Cont,1,nBas,1)
call GetMem('V       ','FREE','REAL',iV,iSizep+4)
call GetMem('SS      ','FREE','REAL',iSS,iSizep+4)
call GetMem('Kin     ','FREE','REAL',iK,iSizep+4)
#else

! Note: in combination with ECPs V is only based on the effective
!       charges of the atoms. In the primitive basis, however, we
!       have temporarily introduced the actual atomic charges. We
!       have to fix this now. Hence the somewhat strange way in
!       which the DKH corrected Hamiltonian is computed.
!
! Compute stripped non-relativistic H (iH_Temp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Open ONEREL

iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
if (iRC /= 0) Go To 9999

call OneBas('PRIM')

Label = 'Kinetic '
call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
Label = 'Attract '
call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
! Add Kinetic and Attraction term
call DaXpY_(iSizep+4,One,Work(iSS),1,Work(iV),1)
if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Prim,1)
  call PrMtrx('Attract+Kinetic (prim)',[lOper],nComp,[iV],Work)
  call PrMtrx('Kinetic (prim)',[lOper],nComp,[iSS],Work)
  call iSwap(8,nBas,1,nBas_Prim,1)
end if
call dcopy_(4,[Zero],0,Work(iH_Temp+iSizec),1)
! Contract and store in iH_temp
call repmat(idbg,Work(iV),Work(iH_temp),.true.)

call GetMem('V       ','FREE','REAL',iV,iSizep+4)
call GetMem('SS      ','FREE','REAL',iSS,iSizep+4)

! Close ONEREL and re-open ONEINT

iOpt = 0
iRC = -1
call ClsOne(iRC,iOpt)
if (iRC /= 0) Go To 9999
iOpt = 0
iRC = -1
Lu_One = 2
call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
if (iRC /= 0) Go To 9999

if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Cont,1)
  call PrMtrx('iH_temp (cont)',[lOper],nComp,[iH_temp],Work)
  call iSwap(8,nBas,1,nBas_Cont,1)
end if

! Transform DKH Hamiltonian to contracted basis (iH)

if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Prim,1)
  call PrMtrx('iK (prim)',[lOper],nComp,[iK],Work)
  call iSwap(8,nBas,1,nBas_Prim,1)
end if
call repmat(idbg,Work(iK),Work(iH),.true.)
if (iPrint >= 20) then
  call iSwap(8,nBas,1,nBas_Cont,1)
  call PrMtrx('iH (cont)',[lOper],nComp,[iH],Work)
  call iSwap(8,nBas,1,nBas_Cont,1)
end if

call GetMem('Kin     ','FREE','REAL',iK,iSizep+4)

iOpt = 0
iRC = -1
Label = 'OneHam 0'
call RdOne(iRC,iOpt,Label,1,Work(iH_nr),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iOpt = 0
iRC = -1

! final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
! where (s) is stripped and with full charge

call DaXpY_(iSizec+4,-One,Work(iH_temp),1,Work(iH),1)
call DaXpY_(iSizec+4,One,Work(iH_nr),1,Work(iH),1)

call Get_iArray('nBas',nBas,nSym)
if (iPrint >= 10) then
  write(u6,'(a11,10i5)') ' Symmetries',nSym
  write(u6,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
end if
Label = 'OneHam 0'
lOper = 1
nComp = 1
ipaddr(1) = iH
if (iPrint >= 20) call PrMtrx(Label,[lOper],nComp,ipaddr,Work)

! Replace 1-el Hamiltonian on ONEINT

iRC = -1
call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
Label = 'OneHam  '
call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
if (iRC /= 0) then
  write(u6,*) 'DKRelInt: Error reading from ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

if (IRELAE == 23) then   ! IORA

  ! Replace overlap on ONEINT

  call repmat(idbg,Work(ipVp),Work(iH),.true.)
  iRC = -1
  Label = 'Mltpl  0'
  call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
  if (iRC /= 0) then
    write(u6,*) 'DKInt: Error reading from ONEINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
end if
#endif

call GetMem('OneHam  ','FREE','REAL',iH,iSizec+4)
call GetMem('H_o     ','FREE','REAL',iH_nr,iSizec+4)
call GetMem('H_temp  ','FREE','REAL',iH_temp,iSizec+4)
call GetMem('pVp     ','FREE','REAL',ipVp,iSizep+4)

return

9999 continue
write(u6,*) ' *** Error in subroutine DKRelInt ***'
write(u6,*) '     Abend in subroutine OpnOne'
call Abend()

end subroutine DKRelint_DP
