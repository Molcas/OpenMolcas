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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  RassiHandM
!
!> @brief
!>   Make a multicenter multipole expansion of the various densities in the RASSI-state Hamiltonian to be.
!>   We also construct the gas-phase RASSI-state Hamiltonian
!> @author A. Ohrn
!>
!> @details
!> First construct the unperturbed RASSI-state Hamiltonian. The \c RasEne
!> are given in input (could be changed later). Then we obtain the
!> MME of the densities of each unique pair of AO-basis functions,
!> which we with the transition density matrices transform to their
!> RASSI-state counterparts---a process that requires some knowledge
!> about that matrix. For example, it should be noted that the matrix
!> is triangularily stored *with* corrections made for the difference
!> between diagonal and non-diagonal elements; therefore we do not
!> need to treat them differently, like we have to in the subroutine
!> ::scfhandm. If requested, we compute total charges and dipoles of
!> every state and print.
!>
!> @note
!> Requires Qfread and of course a RASSI-computation and also MPPROP.
!>
!> @param[in] nBas     Number of contracted AO-basis functions
!> @param[in] iQ_Atoms
!> @param[in] nOcc     Number of basis functions on the \f$ i \f$ th atom-type.
!> @param[in] natyp    Number of atoms of the \f$ i \f$ th atom-type
!> @param[in] nntyp    Number of atom-types
!***********************************************************************

subroutine RassiHandM(nBas,iQ_Atoms,nOcc,natyp,nntyp)

use qmstat_global, only: AddExt, ChaNuc, ChargedQM, CT, HmatSOld, HmatState, iPrint, lSlater, MoAveRed, MxMltp, MxSymQ, nMlt, &
                         nRedMO, nState, outxyzRAS, qTot, RasCha, RasDip, RasQua
use Index_Functions, only: nTri3_Elem, nTri_Elem
use Data_Structures, only: Alloc1DArray_Type, Allocate_DT, Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas(MxSymQ), iQ_Atoms, nntyp, nOcc(nntyp), natyp(nntyp)
integer(kind=iwp) :: i, iCi, j, k, kaunter, kk, nTyp
real(kind=wp) :: D1, D2, D3, dipx, dipx0, dipy, dipy0, dipz, dipz0, dQxx, dQxy, dQxz, dQyy, dQyz, dQzz, dTox, dToy, dToz, Q, qEl, &
                 quaDxx, quaDxy, quaDxz, quaDyy, quaDyx, quaDyz, quaDzx, quaDzy, quaDzz, quaQxx, quaQxy, quaQxz, quaQyy, quaQyz, &
                 quaQzz, quaxx, quaxy, quaxz, quayy, quayz, quazz, Tra, Trace1, Trace2, Trace3
integer(kind=iwp), allocatable :: Dum(:,:), iCent(:)
type(Alloc1DArray_Type), allocatable :: MME(:)

! A modest entrance.

iCi = nTri_Elem(iQ_Atoms)
call mma_allocate(RasCha,nTri_Elem(nState),iCi,label='RasCha')
call mma_allocate(RasDip,nTri_Elem(nState),3,iCi,label='RasDip')
call mma_allocate(RasQua,nTri_Elem(nState),6,iCi,label='RasQua')

! Zeros.

RasCha(:,:) = Zero
RasDip(:,:,:) = Zero
RasQua(:,:,:) = Zero

! Construct H_0 with external perturbation if requested. Construct a copy also.

write(u6,*)
if (.not. AddExt) then
  write(u6,*) '     Constructs H_0.'
else
  write(u6,*) '     Constructs H_0 with external perturbation.'
end if
call Chk_OneHam(nBas)
call RasH0(nBas(1))
HMatSOld(:) = HmatState
write(u6,*) '     ...Done!'
write(u6,*)

! Here the MME is performed. First the expansion is performed in the
! AO-basis. Then depending on whether we use a reduced MO-basis or
! not, the proper approach is chosen.

! Say what we do.

write(u6,*)
write(u6,*) '     Multicenter multipole expanding the charge density expressed in RASSI eigenstates.'

! First obtain MME in AO-basis.

call mma_allocate(iCent,nTri_Elem(nBas(1)),label='iCent')
call mma_allocate(outxyzRAS,3,nTri_Elem(iQ_Atoms),label='outxyzRAS')
call mma_allocate(Dum,nBas(1),nBas(1),label='Dummy')
call Allocate_DT(MME,[1,nTri3_Elem(MxMltp)],label='MME')
call MultiNew(iQ_Atoms,nBas(1),nOcc,natyp,nntyp,MME,iCent,Dum,nMlt,outxyzRAS,lSlater)
call mma_deallocate(Dum)

! Set nTyp, which is number of unique multipole components.

nTyp = 0
do i=1,nMlt
  nTyp = nTyp+nTri_Elem(i)
end do

! Transform to State-basis. The logical flag MoAveRed decides which path to go in this subroutine.

call StateMME(MoAveRed,nBas(1),nRedMO,nState,nTyp,MME,iCent,RasCha,RasDip,RasQua)

! Deallocate the MME in AO-basis.

call mma_deallocate(iCent)

call Deallocate_DT(MME)

! Buckinghamification of the quadrupoles.

RasQua(:,:,:) = RasQua*OneHalf
kaunter = 0
do i=1,nState
  do j=1,i
    kaunter = kaunter+1
    do k=1,iCi
      Tra = (RasQua(kaunter,1,k)+RasQua(kaunter,3,k)+RasQua(kaunter,6,k))/Three
      RasQua(kaunter,1,k) = RasQua(kaunter,1,k)-Tra
      RasQua(kaunter,3,k) = RasQua(kaunter,3,k)-Tra
      RasQua(kaunter,6,k) = RasQua(kaunter,6,k)-Tra
    end do
  end do
end do

! And do some printing if asked for.

if (iPrint >= 10) then
  write(u6,*)
  write(u6,*) '    Distributed multipoles for each state'
  do i=1,nState
    k = nTri_Elem(i)
    write(u6,*) '     State ',i
    do j=1,iCi
      write(u6,*) '        Center ',j
      Q = -RasCha(k,j)
      if (j <= iQ_Atoms) Q = Q+ChaNuc(j)
      write(u6,*) '          ',Q
      D1 = -RasDip(k,1,j)
      D2 = -RasDip(k,2,j)
      D3 = -RasDip(k,3,j)
      write(u6,*) '          ',D1,D2,D3
    end do
  end do
end if
if (iPrint >= 5) then
  write(u6,*)
  write(u6,*) '    Summed multipoles for each state (not state-overlaps)'
  write(u6,*) '              Charge  Dipole(x)   Dipole(y)   Dipole(z)   '// &
              'Quadrup(xx) Quadrup(xy) Quadrup(xz) Quadrup(yy) Quadrup(yz) Quadrup(zz)'
end if
do i=1,nState !Total charge and dipole.
  k = nTri_Elem(i)
  qEl = Zero
  dipx = Zero
  dipy = Zero
  dipz = Zero
  dipx0 = Zero
  dipy0 = Zero
  dipz0 = Zero
  quaxx = Zero
  quaxy = Zero
  quaxz = Zero
  quayy = Zero
  quayz = Zero
  quazz = Zero
  quaDxx = Zero
  quaDxy = Zero
  quaDxz = Zero
  quaDyx = Zero
  quaDyy = Zero
  quaDyz = Zero
  quaDzx = Zero
  quaDzy = Zero
  quaDzz = Zero
  quaQxx = Zero
  quaQxy = Zero
  quaQxz = Zero
  quaQyy = Zero
  quaQyz = Zero
  quaQzz = Zero
  qtot = Zero
  dTox = Zero
  dToy = Zero
  dToz = Zero
  dQxx = Zero
  dQxy = Zero
  dQxz = Zero
  dQyy = Zero
  dQyz = Zero
  dQzz = Zero
  Trace1 = Zero
  Trace2 = Zero
  Trace3 = Zero
  do j=1,iCi
    qEl = qEl+RasCha(k,j)
    dipx = dipx+RasDip(k,1,j)
    dipy = dipy+RasDip(k,2,j)
    dipz = dipz+RasDip(k,3,j)
    dipx0 = dipx0+RasCha(k,j)*outxyzRAS(1,j)
    dipy0 = dipy0+RasCha(k,j)*outxyzRAS(2,j)
    dipz0 = dipz0+RasCha(k,j)*outxyzRAS(3,j)
    quaxx = quaxx+RasQua(k,1,j)
    quaxy = quaxy+RasQua(k,2,j)
    quaxz = quaxz+RasQua(k,4,j)
    quayy = quayy+RasQua(k,3,j)
    quayz = quayz+RasQua(k,5,j)
    quazz = quazz+RasQua(k,6,j)
    quaDxx = quaDxx+RasDip(k,1,j)*(outxyzRAS(1,j)-CT(1))
    quaDxy = quaDxy+RasDip(k,1,j)*(outxyzRAS(2,j)-CT(2))
    quaDxz = quaDxz+RasDip(k,1,j)*(outxyzRAS(3,j)-CT(3))
    quaDyx = quaDyx+RasDip(k,2,j)*(outxyzRAS(1,j)-CT(1))
    quaDyy = quaDyy+RasDip(k,2,j)*(outxyzRAS(2,j)-CT(2))
    quaDyz = quaDyz+RasDip(k,2,j)*(outxyzRAS(3,j)-CT(3))
    quaDzx = quaDzx+RasDip(k,3,j)*(outxyzRAS(1,j)-CT(1))
    quaDzy = quaDzy+RasDip(k,3,j)*(outxyzRAS(2,j)-CT(2))
    quaDzz = quaDzz+RasDip(k,3,j)*(outxyzRAS(3,j)-CT(3))
    quaQxx = quaQxx+RasCha(k,j)*(outxyzRAS(1,j)-CT(1))*(outxyzRAS(1,j)-CT(1))
    quaQxy = quaQxy+RasCha(k,j)*(outxyzRAS(1,j)-CT(1))*(outxyzRAS(2,j)-CT(2))
    quaQxz = quaQxz+RasCha(k,j)*(outxyzRAS(1,j)-CT(1))*(outxyzRAS(3,j)-CT(3))
    quaQyy = quaQyy+RasCha(k,j)*(outxyzRAS(2,j)-CT(2))*(outxyzRAS(2,j)-CT(2))
    quaQyz = quaQyz+RasCha(k,j)*(outxyzRAS(2,j)-CT(2))*(outxyzRAS(3,j)-CT(3))
    quaQzz = quaQzz+RasCha(k,j)*(outxyzRAS(3,j)-CT(3))*(outxyzRAS(3,j)-CT(3))
  end do
  do kk=1,iQ_Atoms
    qtot = qtot+ChaNuc(kk)
    dTox = dTox+ChaNuc(kk)*outxyzRAS(1,kk)
    dToy = dToy+ChaNuc(kk)*outxyzRAS(2,kk)
    dToz = dToz+ChaNuc(kk)*outxyzRAS(3,kk)
    dQxx = dQxx+ChaNuc(kk)*(outxyzRAS(1,kk)-CT(1))*(outxyzRAS(1,kk)-CT(1))
    dQxy = dQxy+ChaNuc(kk)*(outxyzRAS(1,kk)-CT(1))*(outxyzRAS(2,kk)-CT(2))
    dQxz = dQxz+ChaNuc(kk)*(outxyzRAS(1,kk)-CT(1))*(outxyzRAS(3,kk)-CT(3))
    dQyy = dQyy+ChaNuc(kk)*(outxyzRAS(2,kk)-CT(2))*(outxyzRAS(2,kk)-CT(2))
    dQyz = dQyz+ChaNuc(kk)*(outxyzRAS(2,kk)-CT(2))*(outxyzRAS(3,kk)-CT(3))
    dQzz = dQzz+ChaNuc(kk)*(outxyzRAS(3,kk)-CT(3))*(outxyzRAS(3,kk)-CT(3))
  end do
  ! Observe! qTot is not just a check. It is used later as a sign to see if QM-region is charged.
  qtot = qtot-qEl
  dTox = dTox-dipx-dipx0
  dToy = dToy-dipy-dipy0
  dToz = dToz-dipz-dipz0
  Trace1 = dQxx+dQyy+dQzz
  Trace2 = -quaDxx-quaDyy-quaDzz
  Trace3 = -quaQxx-quaQyy-quaQzz
  Trace1 = Trace1/Three
  Trace2 = Two*Trace2/Three
  Trace3 = Trace3/Three
  dQxx = OneHalf*(dQxx-Two*quaDxx-quaQxx-Trace1-Trace2-Trace3)-quaxx
  dQyy = OneHalf*(dQyy-Two*quaDyy-quaQyy-Trace1-Trace2-Trace3)-quayy
  dQzz = OneHalf*(dQzz-Two*quaDzz-quaQzz-Trace1-Trace2-Trace3)-quazz
  dQxy = OneHalf*(dQxy-quaDxy-quaDyx-quaQxy)-quaxy
  dQxz = OneHalf*(dQxz-quaDxz-quaDzx-quaQxz)-quaxz
  dQyz = OneHalf*(dQyz-quaDyz-quaDzy-quaQyz)-quayz
  if (iPrint >= 5) then
    write(u6,9001) '      State ',i,'  ',qtot,dTox,dToy,dToz,dQxx,dQxy,dQxz,dQyy,dQyz,dQzz
  end if
  if (abs(qtot) > 1.0e-4_wp) ChargedQM = .true.
end do
write(u6,*) '     ...Done!'

!----------------------------------------------------------------------*
! The end has come.                                                    *
!----------------------------------------------------------------------*
return

!Jose This format has problems to print anions
!9001 format(A,i2,A,F5.3,9(F12.8))
9001 format(A,i2,A,F5.1,9(F12.8))

end subroutine RassiHandM
