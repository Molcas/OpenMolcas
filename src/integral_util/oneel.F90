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
! Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine OneEl(Kernel,KrnlMm,Label,ip,lOper,nComp,CoorO,nOrdOp,rNuc,rHrmt,iChO,PtChrg,nGrid,iAddPot)

use Basis_Info, only: nBas
use PrpPnt, only: nDen, nOcc, nVec, Occ, Vec
use Gateway_global, only: IfAllOrb, PrPrt, Short
use Sizes_of_Seward, only: S
use Gateway_Info, only: Thrs
use Symmetry_Info, only: nIrrep
use Integral_interfaces, only: int_kernel, int_mem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), nOrdOp, iChO(nComp), nGrid, iAddPot
integer(kind=iwp), intent(out) :: ip(nComp)
real(kind=wp), intent(in) :: CoorO(3,nComp), rNuc(nComp), rHrmt, PtChrg(nGrid)
integer(kind=iwp) :: iAdr, iComp, iComp_, iDIsk, iEF, ii, iInd1, iInd2, iIrrep, iOpt, iPAMCount, ipC2, ipEl, ipNuc, ipOut, iRC, &
                     iSmLbl, iStabO(0:7), jComp, LenInt, LenTot, llOper, lPole, LuTmp, mDim, n_Int, nIC, nStabO
real(kind=wp) :: rSum
character(len=8) :: L_Temp
character(len=4) :: LBL
real(kind=wp), allocatable :: Array(:), El(:), Nuc(:), Out_(:)
integer(kind=iwp), external :: n2Tri

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) ' In OneEl: Label',Label
write(u6,*) ' In OneEl: nComp'
write(u6,'(1X,8I5)') nComp
write(u6,*) ' In OneEl: lOper'
write(u6,'(1X,8I5)') lOper
write(u6,*) ' In OneEl: n2Tri'
do iComp=1,nComp
  ip(iComp) = n2Tri(lOper(iComp))
end do
write(u6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
call RecPrt(' CoorO',' ',CoorO,3,nComp)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of blocks from each component of the operator
! and the irreps it will span.

nIC = 0
llOper = 0
do iComp=1,nComp
  llOper = ior(llOper,lOper(iComp))
  do iIrrep=0,nIrrep-1
    if (btest(lOper(iComp),iIrrep)) nIC = nIC+1
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' nIC =',nIC
#endif

if (nIC == 0) return

call SOS(iStabO,nStabO,llOper)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

ip(:) = -1
LenTot = 0
do iComp=1,nComp
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(Array,LenTot,label='Array')
Array(:) = Zero
ip(1) = 1
iadr = ip(1)
do iComp=1,nComp
  LenInt = n2Tri(lOper(iComp))
  ip(icomp) = iadr
  iadr = iadr+LenInt+4
  ! Copy center of operator to work area.
  Array(ip(iComp)+LenInt:ip(iComp)+LenInt+2) = CoorO(:,iComp)
  ! Copy nuclear contribution to work area.
  Array(ip(iComp)+LenInt+3) = rNuc(iComp)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute all SO integrals for all components of the operator.

call OneEl_Inner(Kernel,KrnlMm,Label,ip,lOper,nComp,CoorO,nOrdOp,rHrmt,iChO,iStabO,nStabO,nIC,PtChrg,nGrid,iAddPot,Array,LenTot)
!                                                                      *
!***********************************************************************
!                                                                      *
!                    P O S T P R O C E S S I N G                       *
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call PrMtrx(Label,lOper,nComp,ip,Array)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Make a square sum on all the integrals for verification
!call VrfMtrx(Label,lOper,nComp,ip,Array)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute properties or write integrals to disc.

do iComp=1,nComp
  iSmLbl = lOper(iComp)
  if (Prprt) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute properties directly from integrals

    ! Allocate some memory

    if (iComp == 1) then
      if (short) then
        mDim = 1
      else
        mDim = S%nDim
      end if
      call mma_allocate(Out_,mDim*nComp,label='Out_')
      Out_(:) = Zero
      ipOut = 1
      call mma_allocate(Nuc,nComp,label='Nuc')
      Nuc(:) = Zero
      ipNuc = 1
    end if
    n_Int = n2Tri(iSmLbl)
    if (n_Int /= 0) call CmpInt(Array(ip(iComp)),n_Int,nBas,nIrrep,iSmLbl)
    Nuc(ipNuc+(iComp-1)) = Array(ip(iComp)+n_Int+3)
    if (n_Int /= 0) call XProp(Short,ifallorb,nIrrep,nBas,nVec,Vec,nOcc,Occ,nDen,Array(ip(iComp)),Out_(ipOut+(iComp-1)*mDim))

    if (Label(1:3) == 'PAM') then
      !open(unit=28,file='R_vect',access='append')
      endfile(28)
      if (Short) then
        write(28,'(a8,2x,f20.14)') Label,Out_(ipOut+(iComp-1)*mDim)
      else
        rSum = Zero
        do ii=1,nOcc
          rSum = rSum+Out_(ipOut+ii-1+(iComp-1)*mDim)
        end do
        write(28,'(a8,2x,f20.14)') Label,-rSum
      end if
      !close(28)
    end if

    ! Once all components have been computed print them.

    if (iComp == nComp) then
      LBL = Label(1:4)
      call UpCase(LBL)
      lpole = 0
      if (LBL == 'MLTP') then
        read(Label,'(5X,I3)') lpole
      else if (LBL == 'PAM ') then
        read(Label,'(5X,I1)') lpole
      else if (LBL == 'L_MP') then
        read(Label,'(5X,I1)') lpole
      else if (LBL(1:2) == 'EF') then
        read(Label,'(2X,I1)') lpole
      else if (LBL == 'DMS ') then
        lpole = 3
      else if (LBL == 'VELO') then
        lpole = 1
      end if
      if (nComp == 1) then
        ipC2 = 1 ! dummy
      else
        ipC2 = 2 ! Used only for diamagnetic shielding.
      end if
      call Prop(Short,Label,CoorO(1,1),CoorO(1,ipC2),nIrrep,nBas,mDim,Occ,Thrs,Out_,Nuc,lpole,ifallorb)

      ! For a properties calculation, save the values of EF or CNT operators,
      ! they will be used to write the sum through Add_Info in Drv1El

      if (PrPrt .and. ((LBL(1:2) == 'EF') .or. (LBL(1:3) == 'CNT'))) then
        call mma_allocate(El,nComp,label='El')
        El(:) = Zero
        ipEl = 1
        ! Compute the sum of all orbital components
        do jComp=0,nComp-1
          iInd1 = ipEl+jComp
          iInd2 = ipOut+jComp*mDim
          El(iInd1) = El(iInd1)+sum(Out_(iInd2:iInd2+mDim-1))
        end do
        ! Write electronic and nuclear components in temp file
        LuTmp = 10
        call DaName(LuTmp,'TMPPRP')
        read(Label(4:8),*) iEF
        iDisk = (iEF-1)*2
        call dDaFile(LuTmp,1,El,nComp,iDisk)
        call dDaFile(LuTmp,1,Nuc,nComp,iDisk)
        call DaClos(LuTmp)
        call mma_deallocate(El)
      end if

      call mma_deallocate(Nuc)
      call mma_deallocate(Out_)
    end if
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Write integrals to disc

    iOpt = 0
    iRC = -1
    if (Label(1:3) == 'PAM') then
      write(L_Temp,'(A5,I3.3)') 'PAM  ',iPAMcount
      iPAMcount = iPAMcount+1
      iComp_ = 1
    else if (Label(1:5) == 'EMFR0') then
      iComp_ = 1
      if (iComp == 1) then
        L_Temp = 'EMFR0  R'
      else
        L_Temp = 'EMFR0  I'
      end if
    else if (Label(1:5) == 'EMFR ') then
      iComp_ = mod(iComp+2,3)+1
      if ((iComp+2)/3 == 1) then
        L_Temp = 'EMFR  RS'
      else if ((iComp+2)/3 == 2) then
        L_Temp = 'EMFR  RA'
      else if ((iComp+2)/3 == 3) then
        L_Temp = 'EMFR  IS'
      else if ((iComp+2)/3 == 4) then
        L_Temp = 'EMFR  IA'
      end if
    else if (Label(1:5) == 'TMOM0') then
      if (iComp == 1) then
        L_Temp = 'TMOM0  R'
        iComp_ = 1
      else
        L_Temp = 'TMOM0  I'
        iComp_ = 1
      end if
    else if (Label(1:5) == 'TMOM2') then
      if (iComp == 1) then
        L_Temp = 'TMOM2  R'
        iComp_ = 1
      else
        L_Temp = 'TMOM2  I'
        iComp_ = 1
      end if
    else if (Label(1:5) == 'TMOM ') then
      iComp_ = mod(iComp+2,3)+1
      if ((iComp+2)/3 == 1) then
        L_Temp = 'TMOM  RS'
      else if ((iComp+2)/3 == 2) then
        L_Temp = 'TMOM  RA'
      else if ((iComp+2)/3 == 3) then
        L_Temp = 'TMOM  IS'
      else if ((iComp+2)/3 == 4) then
        L_Temp = 'TMOM  IA'
      end if
    else
      L_Temp = Label
      iComp_ = iComp
    end if
    call WrOne(iRC,iOpt,L_Temp,iComp_,Array(ip(iComp)),iSmLbl)

    if (iRC /= 0) then
      call WarningMessage(2,' *** Error in subroutine ONEEL ***,     Abend in subroutine WrOne')
      call Abend()
    end if
  end if
end do  ! iComp

if (Label == 'Attract ') call Add_info('SEWARD_ATTRACT',Array(ip(1)),1,5)
if (Label == 'Kinetic ') call Add_info('SEWARD_KINETIC',Array(ip(1)),1,5)
if (Label == 'Mltpl  1') call Add_info('SEWARD_MLTPL1X',Array(ip(1)),1,5)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory for integral

call mma_deallocate(Array)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine OneEl
