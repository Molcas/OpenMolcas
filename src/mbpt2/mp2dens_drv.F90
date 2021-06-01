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

implicit real*8(a-h,o-z)
logical Done
integer iVecOff(8), nOccAll(8), nOrbAll(8)
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
!                                                                      *
!***********************************************************************
!                                                                      *
! Statement functions
iMult(i,j,k) = ip_mult+iVecOff(k)+j-1+(nFro(k)+nOcc(k))*(i-1)
iDensVirOcc(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Start

Eps = 1.0D-8
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
  write(6,*) 'Symmetry nr',iSym
  call RecPrt('InvDia','',Work(ip_DiaA(iSym)),nFro(iSym)+nOcc(iSym),nExt(iSym)+nDel(iSym))
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
lVec = l_Mp2Lagr
do iSym=1,nSym
  iPoVec(iAdr+1) = iPoVec(iAdr)+(nFro(iSym)+nOcc(iSym))*(nExt(iSym)+nDel(iSym))
  iAdr = iAdr+1
end do

! Allocate the vectors needed for the PCG

call GetMem('z_vector','Allo','Real',ip_Z,lVec)
call GetMem('z-next','Allo','Real',ip_ZN,lVec)
call GetMem('r_vector','Allo','Real',ip_R,lVec)
call GetMem('r_next','Allo','Real',ip_RN,lVec)
call GetMem('p_vector','Allo','Real',ip_P,lVec)
call GetMem('p_next','Allo','Real',ip_PN,lVec)
call GetMem('Ap_vector','Allo','Real',ip_AP,lVec)
call GetMem('LagrMult','Allo','Real',ip_Mult,lVec)
call GetMem('LagrMult_next','Allo','Real',ip_MultN,lVec)

! Initiate all vectors to zero.

call FZero(Work(ip_Z),lVec)
call FZero(Work(ip_ZN),lVec)
call FZero(Work(ip_P),lVec)
call FZero(Work(ip_PN),lVec)
call FZero(Work(ip_R),lVec)
call FZero(Work(ip_RN),lVec)
call FZero(Work(ip_Mult),lVec)
call FZero(Work(ip_MultN),lVec)

iVecOff(1) = 0
do iSym=2,nSym
  iVecOff(iSym) = iVecOff(iSym-1)+(nFro(iSym-1)+nOcc(iSym-1))*(nExt(iSym-1)+nDel(iSym-1))
end do

! Calculate initial values for some of the vectors

do iSym=1,nSym
  nI = nFro(iSym)+nOcc(iSym)
  nA = nExt(iSym)+nDel(iSym)
# ifdef _DEBUGPRINT_
  call RecPrt('(ia|ia)',' ',Work(ip_DiaA(iSym)),nI,nA)
  call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),nI,nA)
# endif
  do i=1,nI*nA
    Work(ip_Z+iVecOff(iSym)+i-1) = Work(ip_Mp2Lagr(iSym)+i-1)*Work(ip_DiaA(iSym)+i-1)
    Work(ip_P+iVecOff(iSym)+i-1) = Work(ip_Z+iVecOff(iSym)+i-1)
    Work(ip_R+iVecOff(iSym)+i-1) = Work(ip_Mp2Lagr(iSym)+i-1)
  end do
end do

! Check if the mp2-lagrangian is zero, in that case the cphf-solution
! is trivial and P_ia = 0 for all i and a, in either case
! MP2Lagr should be deallocated.

TotLagr = 0.0d0
do iSym=1,nSym
  do i=1,(nFro(iSym)+nOcc(iSym))*(nExt(iSym)+nDel(iSym))
    TotLagr = TotLagr+Work(ip_Mp2Lagr(iSym)+i-1)
  end do
end do
if (abs(TotLagr) < 1.0d-12) then
  call GetMem('MP2Lagr','Free','Real',ip_First_Mp2Lagr,l_Mp2Lagr)
  Go To 100
else
  call GetMem('MP2Lagr','Free','Real',ip_First_Mp2Lagr,l_Mp2Lagr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Now we have all the initial stuff and should make the PCG-loop.
! The maximum number of iterations are the one defined globally
! in MCLR and not MP2-specific.
do Iter=1,nIter

# ifdef _DEBUGPRINT_
  write(6,*) 'P ITER:',Iter
  do i=0,lVec-1
    write(6,*) Work(ip_P+i)
  end do
# endif

  ! The reason for not doing the whole CG to a black box routine
  ! is that the quantity A*p is method dependent since A is too
  ! big to store on disk so we calculate this outside for each
  ! iteration.

  call FZero(Work(ip_AP),lVec)
  do iSymIA=1,nSym
    do iSymJB=1,iSymIA
      if ((nOrb(iSymIA)+nDel(iSymIA))*(nOrb(iSymJB)+nDel(iSymJB)) /= 0) then
        call MP2Ap(iSymIA,iSymJB,ip_AP,ip_P)
      end if
    end do
  end do
# ifdef _DEBUGPRINT_
  !write(6,*) 'MP2Ap'
  !do i=0,lVec-1
  !  write(6,*) Work(ip_Ap+i)
  !end do
# endif
  ! Makes a call to a routine that makes one CG-update and checks convergence.
  call Conj_Grad(Done,lVec,Work(ip_DiaA(1)),Work(ip_Mult),Work(ip_MultN),Work(ip_R),Work(ip_RN),Work(ip_P),Work(ip_PN),Work(ip_Z), &
                 Work(ip_ZN),Work(ip_AP),Eps,res)
  if (Done) goto 100
end do

write(6,*) '***************WARNING************************'
write(6,*) ''
write(6,*) 'Too many iterations, this is what you get after 50'
write(6,*) 'The residual is',res,'and not',Eps
write(6,*) '**********************************************'

! The PCG is done and we now have the Lagrange multipliers we sought.
! This is a good time to release all the memory related to PCG.

100 continue
do iSym=1,nSym
  do iI=1,nFro(iSym)+nOcc(iSym)
    do iA=1,nExt(iSym)+nDel(iSym)
      Work(iDensVirOcc(iA,iI,iSym)) = Work(iMult(iA,iI,iSym))
    end do
  end do
end do

! Now we have the upper half of the one particle density matrix and only need
! to copy down off diagonal terms

do iSym=1,nSym
  do i=1,nOrb(iSym)+nDel(iSym)
    do j=1,i-1
      Work(ip_Density(iSym)+i-1+(j-1)*(nOrb(iSym)+nDel(iSym))) = Work(ip_Density(iSym)+j-1+(i-1)*(nOrb(iSym)+nDel(iSym)))
    end do
  end do
end do

#ifdef _DEBUGPRINT_
do iSym=1,nSym
  write(6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Work(ip_Density(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
  call RecPrt('MP2WDensity','',Work(ip_WDensity(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call GetMem('z_vector','Free','Real',ip_Z,lVec)
call GetMem('z-next','Free','Real',ip_ZN,lVec)
call GetMem('r_vector','Free','Real',ip_R,lVec)
call GetMem('r_next','Free','Real',ip_RN,lVec)
call GetMem('p_vector','Free','Real',ip_P,lVec)
call GetMem('p_next','Free','Real',ip_PN,lVec)
call GetMem('Ap_vector','Free','Real',ip_AP,lVec)
call GetMem('LagrMult_next','Free','Real',ip_MultN,lVec)
call GetMem('LagrMult','Free','Real',ip_Mult,lVec)

! Set up Density matrices for MO and AO and one for triangular
! stored AO.

l_TriDens = 0
do iSym=1,nSym
  l_TriDens = l_TriDens+(nOrb(iSym)+nDel(iSym))*(nOrb(iSym)+nDel(iSym)+1)/2
end do

call GetMem('AOTriDens','Allo','Real',ip_AOTriDens,l_TriDens)
call GetMem('AOWTriDens','Allo','Real',ip_WAOTriDens,l_TriDens)
do iSym=1,8
  nOrbAll(iSym) = nOrb(iSym)+nDel(iSym)
  nOccAll(iSym) = nFro(iSym)+nOcc(iSym)
end do

call Finish_WDensity()

do iSym=1,nSym
  do i=1,nOccAll(iSym)
    Work(ip_density(iSym)+i-1+(nOrbAll(iSym))*(i-1)) = Work(ip_density(iSym)+i-1+(nOrbAll(iSym))*(i-1))+2.0d0
  end do
end do

! use the old interface for now ... (RL)

call Build_Mp2Dens_Old(ip_AOTriDens,ip_Density,Work(ipCMO),nSym,nOrbAll,nOccAll,.true.)
call Build_Mp2Dens_Old(ip_WAOTriDens,ip_WDensity,Work(ipCMO),nSym,nOrbAll,nOccAll,.false.)

#ifdef _DEBUGPRINT_
write(6,*) 'Normal Dens'
do i=0,l_TriDens-1
  write(6,*) Work(ip_AOTriDens+i)
end do
write(6,*) 'WDens'
do i=0,l_TriDens-1
  write(6,*) Work(ip_WAOTriDens+i)
end do
#endif

call Put_D1ao_Var(Work(ip_AOTriDens),l_TriDens)
!call Put_D1ao(Work(ip_AOTriDens),l_TriDens)
call Put_Fock_Occ(Work(ip_WAOTriDens),l_TriDens)

call GetMem('AOTriDens','Free','Real',ip_AOTriDens,l_TriDens)
call GetMem('AOWTriDens','Free','Real',ip_WAOTriDens,l_TriDens)
! We now have the density matrix for both the MO-basis and the AO-basis in the
! compact form it is supposed to have on the runfile so what is left is to write
! it to the runfile.
!
! Overwrite nonvariational density to fool LoProp. (should not be done
! this way)
#ifdef _DEBUGPRINT_
write(6,*) 'EMP2 is ',EMP2
write(6,*) ' '
do iSym=1,nSym
  write(6,*) 'Density matrix for Symm:',iSym
  call RecPrt('MP2Density','',Work(ip_Density(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
do iSym=1,nSym
  write(6,*) 'WDensity matrix for Symm:',iSym
  call RecPrt('MP2WDensity','',Work(ip_WDensity(iSym)),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
end do
#endif
call GetMem('MP2Density','Free','Real',ip_First_Density,l_Density)
call GetMem('MP2WDensity','Free','Real',ip_First_WDensity,l_Density)
call GetMem('MP2DiaA','Free','Real',ip_First_DiaA,l_DiaA)

E2BJAI = EMP2
REFC = VECL2

return

end subroutine MP2Dens_drv
