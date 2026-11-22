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

! Warm start of the SA Z-vector equations (SA-MCSCF and MC-PDFT) from the solution
! of the previous call (previous geometry step) with the same target.
! Active inside a geometry loop (MOLCAS_ITER >= 1) or a molecular dynamics
! run (Dynamix data on the runfile), except on the displaced
! geometries of a numerical differentiation in Slapaf; MOLCAS_MCLR_WARM=0
! switches it off.
!
! The solution of the previous call is stored together with its MOs and
! reference CI vectors. It is mapped onto the current MOs with
! T = C_new^T S C_old (kappa' = T kappa T^T) and onto the current
! reference states with O_ij = <c_new_i|c_old_j>.

module ZVec_Guess

use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp), parameter :: nKey = 9

public :: ZG_Active, ZG_Load, ZG_Save

contains

function ZG_Active()
  use input_mclr, only: ColdStart
  logical(kind=iwp) :: ZG_Active
  integer(kind=iwp) :: iErr, Iter
  character(len=80) :: Val
  logical(kind=iwp) :: MD
  ! Only inside a geometry loop (MOLCAS_ITER >= 1) or a molecular dynamics run
  ! (Dynamix has stored 'MD_Time'); the keyword COLD or MOLCAS_MCLR_WARM=0 switches it off
  ZG_Active = .false.
  if (ColdStart) then
    write(u6,*) 'ZVec_Guess: keyword COLD, cold start'
    return
  end if
  Val = ' '
  call get_environment_variable('MOLCAS_MCLR_WARM',Val)
  if (trim(adjustl(Val)) == '0') return
  Val = ' '
  call get_environment_variable('MOLCAS_ITER',Val)
  read(Val,*,iostat=iErr) Iter
  if (iErr /= 0) Iter = 0
  call Qpg_dScalar('MD_Time',MD)
  if ((Iter < 1) .and. (.not. MD)) return
  ZG_Active = .not. NumDiff_Step()
end function ZG_Active

! Slapaf sets 'Slapaf NumDiff' = 1 on the displaced geometries of a numerical
! differentiation: there the errors of a warm-started solution do not cancel
! between displacements, so start cold.
function NumDiff_Step()
  logical(kind=iwp) :: NumDiff_Step
  integer(kind=iwp) :: iFlag
  logical(kind=iwp) :: Found
  NumDiff_Step = .false.
  call Qpg_iScalar('Slapaf NumDiff',Found)
  if (Found) then
    call Get_iScalar('Slapaf NumDiff',iFlag)
    NumDiff_Step = (iFlag == 1)
  end if
  if (NumDiff_Step) write(u6,*) 'ZVec_Guess: numerical differentiation step, cold start'
end function NumDiff_Step

! kind of equations: 0 SA-MCSCF, 1 SA-MCSCF with a PT2 Lagrangian, 2 MC-PDFT
function iEqs(Method)
  use input_mclr, only: PT2
  integer(kind=iwp) :: iEqs
  integer(kind=iwp), intent(in) :: Method
  iEqs = merge(2,merge(1,0,PT2),Method == 1)
end function iEqs

! one file per target: ZV<equations><NAC><root or NAC state 1><NAC state 2>
function FName(Method)
  use MCLR_Data, only: IRLXROOT, ISNAC, NACSTATES
  character(len=8) :: FName
  integer(kind=iwp), intent(in) :: Method
  if (isNAC) then
    write(FName,'(A2,2I1,2I2.2)') 'ZV',iEqs(Method),1,mod(NACstates(1),100),mod(NACstates(2),100)
  else
    write(FName,'(A2,2I1,2I2.2)') 'ZV',iEqs(Method),0,mod(irlxroot,100),0
  end if
end function FName

subroutine Make_Key(Key,Method)
  use MCLR_Data, only: IRLXROOT, ISNAC, NACSTATES, nConf1, nDensC, nCMO
  use input_mclr, only: nRoots
  integer(kind=iwp), intent(out) :: Key(nKey)
  integer(kind=iwp), intent(in) :: Method
  Key(1) = merge(1,0,isNAC)
  Key(2) = NACstates(1)
  Key(3) = NACstates(2)
  Key(4) = irlxroot
  Key(5) = iEqs(Method)
  Key(6) = nDensC
  Key(7) = nConf1
  Key(8) = nRoots
  Key(9) = nCMO
end subroutine Make_Key

! Method: 0 (default) SA-MCSCF (wfctl_sa), 1 MC-PDFT (wfctl_pdft)
subroutine ZG_Save(Kappa,CIT,Xvec,CI,Method)
  use MCLR_Data, only: CMO, nConf1, nDensC, nCMO
  use input_mclr, only: nRoots
  ! inout because dDaFile's buffer is intent(inout) (read and write); the arrays are only written to disk
  real(kind=wp), intent(inout) :: Kappa(nDensC), CIT(nConf1*nRoots), Xvec(nRoots**2), CI(nConf1*nRoots)
  integer(kind=iwp), intent(in), optional :: Method
  integer(kind=iwp) :: iDisk, iMeth, Key(nKey), Lu
  iMeth = 0
  if (present(Method)) iMeth = Method
  Lu = 21
  call Make_Key(Key,iMeth)
  call DaName(Lu,FName(iMeth))
  iDisk = 0
  call iDaFile(Lu,1,Key,nKey,iDisk)
  call dDaFile(Lu,1,CMO,nCMO,iDisk)
  call dDaFile(Lu,1,CI,nConf1*nRoots,iDisk)
  call dDaFile(Lu,1,Kappa,nDensC,iDisk)
  call dDaFile(Lu,1,CIT,nConf1*nRoots,iDisk)
  call dDaFile(Lu,1,Xvec,nRoots**2,iDisk)
  call DaClos(Lu)
end subroutine ZG_Save

subroutine ZG_Load(Kappa,CIT,Xvec,CI,Found,Method)
  use Index_Functions, only: nTri_Elem
  use MCLR_Data, only: CMO, ipMat, nConf1, nDens, nDensC, nCMO
  use input_mclr, only: nBas, nOrb, nRoots
  use Orbital_Spaces, only: nSym
  use OneDat, only: sNoNuc, sNoOri
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  real(kind=wp), intent(out) :: Kappa(nDensC), CIT(nConf1*nRoots), Xvec(nRoots**2)
  real(kind=wp), intent(in) :: CI(nConf1*nRoots)
  logical(kind=iwp), intent(out) :: Found
  integer(kind=iwp), intent(in), optional :: Method
  integer(kind=iwp) :: iMeth, i, iComp, iDisk, iOff, iOffC, iOpt, irc, iSym, iSymLbl, j, k, Key(nKey), KeyOld(nKey), Lu, nB, nO
  real(kind=wp) :: dev, q
  logical(kind=iwp) :: Exists
  character(len=8) :: Label
  real(kind=wp), allocatable :: CIOld(:), CITOld(:), CMOOld(:), Kf(:), Kn(:), KapOld(:), O(:,:), S(:), Sq(:), T(:), Tmp(:), &
                                Tmp2(:), XOld(:,:), X(:,:)
  real(kind=wp), external :: DDot_

  iMeth = 0
  if (present(Method)) iMeth = Method
  Found = .false.
  Kappa(:) = Zero
  CIT(:) = Zero
  Xvec(:) = Zero
  inquire(file=trim(FName(iMeth)),exist=Exists)
  if (.not. Exists) then
    call WarningMessage(1,'ZVec_Guess: no stored Z-vector solution to restart from; cold start')
    return
  end if

  Lu = 21
  call DaName(Lu,FName(iMeth))
  iDisk = 0
  call iDaFile(Lu,2,KeyOld,nKey,iDisk)
  call Make_Key(Key,iMeth)
  if (any(Key /= KeyOld)) then
    call DaClos(Lu)
    write(u6,*) 'ZVec_Guess: stored solution belongs to another target, cold start'
    return
  end if
  call mma_allocate(CMOOld,nCMO,Label='CMOOld')
  call mma_allocate(CIOld,nConf1*nRoots,Label='CIOld')
  call mma_allocate(KapOld,nDensC,Label='KapOld')
  call mma_allocate(CITOld,nConf1*nRoots,Label='CITOld')
  call mma_allocate(XOld,nRoots,nRoots,Label='XOld')
  call dDaFile(Lu,2,CMOOld,nCMO,iDisk)
  call dDaFile(Lu,2,CIOld,nConf1*nRoots,iDisk)
  call dDaFile(Lu,2,KapOld,nDensC,iDisk)
  call dDaFile(Lu,2,CITOld,nConf1*nRoots,iDisk)
  call dDaFile(Lu,2,XOld,nRoots**2,iDisk)
  call DaClos(Lu)

  ! AO overlap at the current geometry
  nB = 0
  do iSym=1,nSym
    nB = nB+nTri_Elem(nBas(iSym))
  end do
  call mma_allocate(S,nB,Label='S')
  iSymLbl = 1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  Label = 'Mltpl  0'
  iComp = 1
  call RdOne(irc,iOpt,Label,iComp,S,iSymLbl)

  ! Orbital part: kappa' = T kappa T^T, T = C_new^T S C_old
  call mma_allocate(Kf,nDens,Label='Kf')
  call mma_allocate(Kn,nDens,Label='Kn')
  call UnCompress(KapOld,Kf,1)
  Kn(:) = Zero
  dev = Zero
  iOff = 0
  iOffC = 1
  do iSym=1,nSym
    nB = nBas(iSym)
    nO = nOrb(iSym)
    if (nO > 0) then
      call mma_allocate(Sq,nB*nB,Label='Sq')
      call mma_allocate(T,nO*nO,Label='T')
      call mma_allocate(Tmp,nB*nO,Label='Tmp')
      call mma_allocate(Tmp2,nO*nO,Label='Tmp2')
      k = iOff
      do i=1,nB
        do j=1,i
          k = k+1
          Sq(j+(i-1)*nB) = S(k)
          Sq(i+(j-1)*nB) = S(k)
        end do
      end do
      call dGemm_('N','N',nB,nO,nB,One,Sq,nB,CMOOld(iOffC),nB,Zero,Tmp,nB)
      call dGemm_('T','N',nO,nO,nB,One,CMO(iOffC),nB,Tmp,nB,Zero,T,nO)
      ! deviation of T from orthogonality, as a diagnostic
      call dGemm_('T','N',nO,nO,nO,One,T,nO,T,nO,Zero,Tmp2,nO)
      do i=1,nO
        Tmp2(i+(i-1)*nO) = Tmp2(i+(i-1)*nO)-One
      end do
      dev = max(dev,maxval(abs(Tmp2)))
      call dGemm_('N','N',nO,nO,nO,One,T,nO,Kf(ipMat(iSym,iSym)),nO,Zero,Tmp2,nO)
      call dGemm_('N','T',nO,nO,nO,One,Tmp2,nO,T,nO,Zero,Kn(ipMat(iSym,iSym)),nO)
      call mma_deallocate(Tmp2)
      call mma_deallocate(Tmp)
      call mma_deallocate(T)
      call mma_deallocate(Sq)
    end if
    iOff = iOff+nTri_Elem(nB)
    iOffC = iOffC+nB*nO
  end do
  call Compress(Kn,Kappa,1)
  call mma_deallocate(Kn)
  call mma_deallocate(Kf)
  call mma_deallocate(S)

  ! CI part: rotate within the reference space, then project it out
  call mma_allocate(O,nRoots,nRoots,Label='O')
  do i=1,nRoots
    do j=1,nRoots
      O(i,j) = DDot_(nConf1,CI(1+(i-1)*nConf1),1,CIOld(1+(j-1)*nConf1),1)
    end do
  end do
  q = sum(O**2)/real(nRoots,kind=wp)
  if (q > 0.9_wp) then
    do i=1,nRoots
      do j=1,nRoots
        CIT(1+(i-1)*nConf1:i*nConf1) = CIT(1+(i-1)*nConf1:i*nConf1)+O(i,j)*CITOld(1+(j-1)*nConf1:j*nConf1)
      end do
    end do
    do i=1,nRoots
      do k=1,nRoots
        q = DDot_(nConf1,CI(1+(k-1)*nConf1),1,CIT(1+(i-1)*nConf1),1)
        CIT(1+(i-1)*nConf1:i*nConf1) = CIT(1+(i-1)*nConf1:i*nConf1)-q*CI(1+(k-1)*nConf1:k*nConf1)
      end do
    end do
    call mma_allocate(X,nRoots,nRoots,Label='X')
    call dGemm_('N','T',nRoots,nRoots,nRoots,One,XOld,nRoots,O,nRoots,Zero,X,nRoots)
    call dGemm_('N','N',nRoots,nRoots,nRoots,One,O,nRoots,X,nRoots,Zero,XOld,nRoots)
    Xvec(:) = reshape(XOld,[nRoots**2])
    call mma_deallocate(X)
    write(u6,'(A,ES10.3,A,F8.5)') ' ZVec_Guess: warm start; max|T^T T-1| =',dev,', reference overlap =',sum(O**2)/nRoots
  else
    write(u6,'(A,ES10.3,A,F8.5)') ' ZVec_Guess: warm start, orbital part only; max|T^T T-1| =',dev,', reference overlap =',q
  end if
  call mma_deallocate(O)

  call mma_deallocate(XOld)
  call mma_deallocate(CITOld)
  call mma_deallocate(KapOld)
  call mma_deallocate(CIOld)
  call mma_deallocate(CMOOld)
  Found = .true.

end subroutine ZG_Load

end module ZVec_Guess
