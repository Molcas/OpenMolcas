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

subroutine MP2Dens_drv(E2BJAI,REFC)
!***********************************************************************
!                                                                      *
! called from: Mp2_driver                                              *
!                                                                      *
!***********************************************************************

use MBPT2_Global, only: CMO, Density, DiaA, EMP2, iPoVec, MP2Lagr, VECL2, WDensity
use Data_Structures, only: Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: E2BJAI, REFC
integer(kind=iwp) :: i, iA, iAdr, iI, iSym, iSymIA, iSymJB, Iter, iVecOff(8), j, l_TriDens, lVec, nIter, nOccAll(8), nOrbAll(8)
real(kind=wp) :: Eps, res, TotLagr
logical(kind=iwp) :: Done
real(kind=wp), allocatable :: AP(:), AOTriDens(:), Mult(:), MultN(:), P(:), PN(:), R(:), RN(:), WAOTriDens(:), Z(:), ZN(:)
#include "corbinf.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Start

Eps = 1.0e-8_wp
Done = .false.
nIter = 100
!                                                                      *
!***********************************************************************
!                                                                      *
call MP2gDens_setup()
call rhs_mp2()
call Mp2Diag()
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
do iSym=1,nSym
  write(u6,*) 'Symmetry nr',iSym
  call RecPrt('InvDia','',DiaA%SB(iSym)%A1,nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Now we are going to solve Ax = b using the pcg-method
! Reference here to the PCG-paper with the right symbols.
! We need an initial guess for x, we will use 0
! this will give r_1 = b - A*x_1 = b so r_1 is our
! rhs, that is the MP2-lagrangian.
!
! Since we are using a preconditioner we need also a
! vector z_1 = inv(A~)*r_1 where A~ is our preconditioner.
! We use the diagonal of A as preconditioner and we have
! what we need to calculate z
! We also need a vector p that is chosen to r_0 in the first
! iteration
!
! All vectors we use are separated into one block for each symmetry.
! To know where a certain symmetry begins iPoVec(iSym) is added to
! the vectors memory adress. nVec is the total length of the vector.

iAdr = 1
iPoVec(1) = 0
lVec = size(Mp2Lagr%A0)
do iSym=1,nSym
  iPoVec(iAdr+1) = iPoVec(iAdr)+(nFro(iSym)+nOcc(iSym))*(nExt(iSym)+nDel(iSym))
  iAdr = iAdr+1
end do

! Allocate the vectors needed for the PCG

call mma_allocate(Z,lVec,label='z_vector')
call mma_allocate(ZN,lVec,label='z_next')
call mma_allocate(R,lVec,label='r_vector')
call mma_allocate(RN,lVec,label='r_next')
call mma_allocate(P,lVec,label='p_vector')
call mma_allocate(PN,lVec,label='p_next')
call mma_allocate(AP,lVec,label='Ap_vector')
call mma_allocate(Mult,lVec,label='LagrMult')
call mma_allocate(MultN,lVec,label='LagrMult_next')

! Initialize vectors to zero.

ZN(:) = Zero
RN(:) = Zero
PN(:) = Zero
Mult(:) = Zero
MultN(:) = Zero

iVecOff(1) = 0
do iSym=2,nSym
  iVecOff(iSym) = iVecOff(iSym-1)+(nFro(iSym-1)+nOcc(iSym-1))*(nExt(iSym-1)+nDel(iSym-1))
end do

! Calculate initial values for some of the vectors

# ifdef _DEBUGPRINT_
do iSym=1,nSym
  call RecPrt('(ia|ia)',' ',DiaA%SB(iSym)%A1,nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
  call RecPrt('MP2Lagr',' ',Mp2Lagr%SB(iSym)%A1,nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
end do
# endif
Z(:) = Mp2Lagr%A0(:)*DiaA%A0(:)
P(:) = Z(:)
R(:) = Mp2Lagr%A0(:)

! Check if the mp2-lagrangian is zero, in that case the cphf-solution
! is trivial and P_ia = 0 for all i and a, in either case
! MP2Lagr should be deallocated.

TotLagr = Zero
do i=1,size(Mp2Lagr%A0)
  TotLagr = TotLagr+Mp2Lagr%A0(i)
end do
call Deallocate_DT(Mp2Lagr)
if (abs(TotLagr) < 1.0e-12_wp) then
  Done = .true.
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now we have all the initial stuff and should make the PCG-loop.
  ! The maximum number of iterations are the one defined globally
  ! in MCLR and not MP2-specific.
  do Iter=1,nIter

#   ifdef _DEBUGPRINT_
    write(u6,*) 'P ITER:',Iter
    do i=1,lVec
      write(u6,*) P(i)
    end do
#   endif

    ! The reason for not doing the whole CG to a black box routine
    ! is that the quantity A*p is method dependent since A is too
    ! big to store on disk so we calculate this outside for each
    ! iteration.

    AP(:) = Zero
    do iSymIA=1,nSym
      do iSymJB=1,iSymIA
        if ((nOrb(iSymIA)+nDel(iSymIA))*(nOrb(iSymJB)+nDel(iSymJB)) /= 0) then
          call MP2Ap(iSymIA,iSymJB,AP,P)
        end if
      end do
    end do
#   ifdef _DEBUGPRINT_
    !write(u6,*) 'MP2Ap'
    !do i=1,lVec
    !  write(u6,*) Ap(i)
    !end do
#   endif
    ! Makes a call to a routine that makes one CG-update and checks convergence.
    call Conj_Grad(Done,lVec,DiaA%A0,Mult,MultN,R,RN,P,PN,Z,ZN,AP,Eps,res)
    if (Done) exit
  end do
end if

if (.not. Done) then
  write(u6,*) '***************WARNING************************'
  write(u6,*) ''
  write(u6,*) 'Too many iterations, this is what you get after 50'
  write(u6,*) 'The residual is',res,'and not',Eps
  write(u6,*) '**********************************************'
end if

! The PCG is done and we now have the Lagrange multipliers we sought.
! This is a good time to release all the memory related to PCG.

do iSym=1,nSym
  do iI=1,nFro(iSym)+nOcc(iSym)
    do iA=1,nExt(iSym)+nDel(iSym)
      Density%SB(iSym)%A2(iI,iA+nFro(iSym)+nOcc(iSym)) = Mult(iVecOff(iSym)+iI+(nFro(iSym)+nOcc(iSym))*(iA-1))
    end do
  end do
end do

! Now we have the upper half of the one particle density matrix and only need
! to copy down off diagonal terms

do iSym=1,nSym
  do i=1,nOrb(iSym)+nDel(iSym)
    do j=1,i-1
      Density%SB(iSym)%A2(i,j) = Density%SB(iSym)%A2(j,i)
    end do
  end do
end do

#ifdef _DEBUGPRINT_
do iSym=1,nSym
  write(u6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Density%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
  call RecPrt('MP2WDensity','',WDensity%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Z)
call mma_deallocate(ZN)
call mma_deallocate(R)
call mma_deallocate(RN)
call mma_deallocate(P)
call mma_deallocate(PN)
call mma_deallocate(AP)
call mma_deallocate(Mult)
call mma_deallocate(MultN)

! Set up Density matrices for MO and AO and one for triangular
! stored AO.

l_TriDens = 0
do iSym=1,nSym
  l_TriDens = l_TriDens+(nOrb(iSym)+nDel(iSym))*(nOrb(iSym)+nDel(iSym)+1)/2
end do

call mma_allocate(AOTriDens,l_TriDens,label='AOTriDens')
call mma_allocate(WAOTriDens,l_TriDens,label='WAOTriDens')
do iSym=1,8
  nOrbAll(iSym) = nOrb(iSym)+nDel(iSym)
  nOccAll(iSym) = nFro(iSym)+nOcc(iSym)
end do

call Finish_WDensity()

do iSym=1,nSym
  do i=1,nOccAll(iSym)
    Density%SB(iSym)%A2(i,i) = Density%SB(iSym)%A2(i,i)+Two
  end do
end do

! use the old interface for now ... (RL)

call Build_Mp2Dens_Old(AOTriDens,Density,CMO,nSym,nOrbAll,.true.)
call Build_Mp2Dens_Old(WAOTriDens,WDensity,CMO,nSym,nOrbAll,.false.)

#ifdef _DEBUGPRINT_
write(u6,*) 'Normal Dens'
do i=1,l_TriDens
  write(u6,*) AOTriDens(i)
end do
write(u6,*) 'WDens'
do i=1,l_TriDens
  write(u6,*) WAOTriDens(i)
end do
#endif

call Put_dArray('D1aoVar',AOTriDens,l_TriDens)
!call Put_dArray('D1ao',AOTriDens,l_TriDens)
call Put_dArray('FockOcc',WAOTriDens,l_TriDens)

call mma_deallocate(AOTriDens)
call mma_deallocate(WAOTriDens)
! We now have the density matrix for both the MO-basis and the AO-basis in the
! compact form it is supposed to have on the runfile so what is left is to write
! it to the runfile.
!
! Overwrite nonvariational density to fool LoProp. (should not be done
! this way)
#ifdef _DEBUGPRINT_
write(u6,*) 'EMP2 is ',EMP2
write(u6,*) ' '
do iSym=1,nSym
  write(u6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Density%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(u6,*) 'WDensity matrix for Symm:',iSym
  call RecPrt('MP2WDensity','',WDensity%SB(iSym)%A1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
#endif
call Deallocate_DT(Density)
call Deallocate_DT(WDensity)
call Deallocate_DT(DiaA)

E2BJAI = EMP2
REFC = VECL2

return

end subroutine MP2Dens_drv
