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

subroutine Compute_Xhole_Int(nBasLop,nSym,nB1,SqMom,Func,nSize)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBasLop(nSym), nB1
real(kind=wp), intent(out) :: SqMom(nB1,nB1), Func
integer(kind=iwp), intent(out) :: nSize
integer(kind=iwp) :: i, iDum(1), iO1, iO2, iOpt, irc, iSmLbl, kaunt1, kaunt2, Lu_One, nB, nCMO, nDens, nDiff, nOrb, nTri
logical(kind=iwp) :: DSCF !, Do_Gamma, Do_Grad, On_Top, Do_Tau, Do_MO, Do_TwoEl
character(len=4) :: DFTFOCK, KSDFT
real(kind=wp), allocatable :: CMO(:), D1ao(:), MatEl(:), Mult1(:), MultSq(:,:), OrbDip(:,:), TEMP(:,:)
integer(kind=iwp), external :: ip_of_Work, IsFreeUnit
#include "nq_info.fh"

! Check symmetry

if (nSym /= 1) then
  write(u6,*)
  write(u6,*) ' You should not run LoProp with symmetry!'
  call Abend()
end if

! Set a lot of numbers and labels. See below for more help.

nB = nBasLop(1)     ! number of basis functions
nTri = nB*(nB+1)/2  ! obvious!
Func = Zero         ! initialize
!IFG: uncomment when DrvNQ call fixed (see below)
!Dens = Zero         ! initialize
!nFckDim = 1         ! for open-shell the fock-matrix has two versions, but this implementation is not for UHF (but for CASSCF)
!nD = 1              ! much like nFckDim
!Do_Gamma = .false.  ! optional stuff in DFT-theory; for our purpose, just scheiss.
!Do_Grad = .false.   ! see above.
!On_Top = .false.    ! see above.
!Do_Tau = .false.    ! see above.
!Do_MO = .true.      ! in our integration kernel, we need MOs.
!Do_TwoEl = .false.  ! scheiss in our functional.

write(DFTFOCK,'(A)') 'XHOL'          ! Tell the routines that we wish to run a xhole calculation.
write(KSDFT,'(A)') 'LDA'             ! Just to get the routines to take the fastest path, no other meaning.
call mma_allocate(D1ao,nTri)
nDens = nTri
call Get_D1ao(D1ao,nDens)            ! The density matrix.
Functional_type = LDA_type           ! Number from nq_info.fh.
call Put_dScalar('EThr',1.0e-9_wp)   ! A threshold for energy accuracy in DFT/SCF. Dummy.
call Seward_Init()                   ! Initialize a lot of shit.
nDiff = 0                            ! like above
DSCF = .false.                       ! like above
call GetInf(DSCF,nDiff)              ! like above
call SetUp_iSD()                     ! like above
call Get_iScalar('nSym',mIrrep)      ! Stupid number in nq_info needed to fool do_mo stuff,
                                     ! since they normally are only activated by CASDFT
call Get_iArray('nBas',mBas,mIrrep)  ! Like above

! Then we need orbital density dipoles.

nCMO = nB**2
call mma_allocate(CMO,nCMO,Label='CMO')
nOrb = int(sqrt(real(nCMO,kind=wp)))
call Get_CMO(CMO,nCMO)
nB = mBas(0)
call mma_allocate(MultSq,nB,nB,label='MultSq')
call mma_allocate(TEMP,nOrb,nB,label='TEMP')
call mma_allocate(Mult1,nOrb**2+4,label='MultiKulti')
call mma_allocate(OrbDip,nOrb*(nOrb+1)/2,3,label='OrbDips')
ip_OrbDip(1) = ip_of_Work(OrbDip(1,1))
ip_OrbDip(2) = ip_of_Work(OrbDip(1,2))
ip_OrbDip(3) = ip_of_Work(OrbDip(1,3))
irc = -1
Lu_One = 49
Lu_One = IsFreeUnit(Lu_One)
call OpnOne(irc,0,'ONEINT',Lu_One)
if (irc /= 0) then
  write(u6,*)
  write(u6,*) 'ERROR! Could not open one-electron integral file.'
  call Abend()
end if
do i=1,3
  irc = -1
  iOpt = 1
  iSmLbl = 1
  call iRdOne(irc,iOpt,'Mltpl  1',i,iDum,iSmLbl)
  nSize = iDum(1)
  irc = -1
  iOpt = 0
  iSmLbl = 0
  call RdOne(irc,iOpt,'Mltpl  1',i,Mult1,iSmLbl)
  call Square(Mult1,MultSq,1,nB,nB)
  call DGEMM_('T','N',nOrb,nB,nB,One,CMO,nB,MultSq,nB,Zero,TEMP,nOrb)
  call DGEMM_('N','N',nOrb,nOrb,nB,One,TEMP,nOrb,CMO,nB,Zero,Mult1,nOrb)
  kaunt1 = 0
  kaunt2 = 0
  do iO1=1,nOrb
    do iO2=1,nOrb
      kaunt2 = kaunt2+1
      if (iO1 >= iO2) then
        kaunt1 = kaunt1+1
        OrbDip(kaunt1,i) = Mult1(kaunt2)
      end if
    end do
  end do
end do
call ClsOne(irc,Lu_One)
call mma_deallocate(MultSq)
call mma_deallocate(TEMP)
call mma_deallocate(Mult1)
call mma_deallocate(CMO)

! Anders' little helper on DrvNQ:
!   Argument (1): The name on the integration kernel. Is a routine
!                 in src/nq_util/ directory.
!            (2): Through some "fooling" of the integration routines
!                 this argument will upon return contain the matrix
!                 elements.
!            (3): nFckDim is described above. Just one.
!            (4): Here the functional "energy" comes, in other words
!                 the molecular expectation value
!            (5): Just shit.
!            (6): The one-electron density matrix.
!            (7): Dimension on all one-electron matrices.
!            (8): Like (3).
!            (9-17): Total scheiss from our perspective.
!            (18): Label to signal that we will need MOs.
!            (19): We use no two-electron stuff, hence false.
!            (20): Usually, this tells what type of Fock-matrix
!                  we use, but now we "rob" this variable for our
!                  purpose and let the value XHOL signify that we
!                  are computing the xhole dipole.
call mma_allocate(MatEl,nTri,label='X-Dipole elements')
!IFG: The call to DrvNQ below is completely messed up, please fix
call WarningMessage(2,'There is surely a bug here!')
!call DrvNQ(Do_XHoleDip,MatEl,nFckDim,Func,Dens,D1ao,nTri,nD,Do_Gamma,Do_Grad,Dummy,iDummy,Dummy,Dummy,iDummy,On_Top,Do_Tau,Do_MO, &
!           Do_TwoEl,DFTFOCK)
!FFF = ddot_(nDens,D1ao,1,MatEl,1)
!write(u6,*) 'YYY:',nDens,FFF,Func

! Put the second-moments in square form.

call Square(MatEl,SqMom,1,nB,nB)

! Deallocate
call mma_deallocate(D1ao)
call Free_iSD()
call mma_deallocate(MatEl)
call mma_deallocate(OrbDip)
call Free_HerRW()

return

end subroutine Compute_Xhole_Int
