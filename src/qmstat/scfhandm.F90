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
!  ScfHandM
!
!> @brief
!>   (i) construct the multicenter multipole expansion for the MO-Hamiltonian,
!>   (ii) read in the one-electron part of the Hamiltonian and
!>   (iii) construct the super-matrix. The last two are MO-transformed
!> @author A. Ohrn
!>
!> @details
!> First call on the MME-routine. It returns things in AO-basis
!> which we need to put into MO-form. Observe that the quadrupoles
!> in Qmstat are ordered differently compared to Molcas. Qmstat:
!> \f$ xx \f$ \f$ xy \f$ \f$ yy \f$ \f$ xz \f$ \f$ yz \f$ \f$ zz \f$; Molcas:
!> \f$ xx \f$ \f$ xy \f$ \f$ xz \f$ \f$ yy \f$ \f$ yz \f$ \f$ zz \f$. Then we
!> read in parts of the unperturbed Hamiltonian and construct the
!> super-matrix.
!>
!> @param[in] Cmo   Orbital coeff.
!> @param[in] nBas  Number of contracted basis functions
!> @param[in] nOcc  Number of contracted basis functions of a certain atom-type
!> @param[in] natyp Number of atoms of a certain atom-type (for water, hydrogen is 2)
!> @param[in] nntyp Number of atom-types in molecule
!> @param[in] Occu  Orbital occupation numbers
!***********************************************************************

subroutine ScfHandM(Cmo,nBas,iQ_Atoms,nOcc,natyp,nntyp,Occu)

use qmstat_global, only: Cha, ChaNuc, ChargedQM, DipMy, iOrb, iPrint, lSlater, Mp2DensCorr, nMlt, outxyz, qTot, Quad, SlExpQ
use Index_Functions, only: iTri, nTri3_Elem, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
real(kind=wp) :: Cmo(MxBas**2), Occu(MxBas)
integer(kind=iwp) :: nBas(MxSym), iQ_Atoms, nOcc(MxBas), natyp(MxAt), nntyp
integer(kind=iwp) :: i, i1, i2, iB1, iB2, iCi, iDum, ii, iMME(nTri3_Elem(MxMltp)), iMtot, indMME, iO1, iO2, ipO, iTyp, j, k, &
                     kaunta, kaunter, kk, nTyp
real(kind=wp) :: cProd, dipx, dipx0, dipy, dipy0, dipz, dipz0, dTox, dToy, dToz, qEl, Tra
character(len=20) :: MMElab
character(len=2) :: ChCo
integer(kind=iwp), allocatable :: iCent(:,:)

!----------------------------------------------------------------------*
! Zeros.                                                               *
!----------------------------------------------------------------------*
iCi = nTri_Elem(iQ_Atoms)
do i=1,nTri_Elem(iOrb(1))
  do j=1,iCi
    Cha(i,j) = Zero
    DipMy(i,1,j) = Zero
    DipMy(i,2,j) = Zero
    DipMy(i,3,j) = Zero
    Quad(i,1,j) = Zero
    Quad(i,2,j) = Zero
    Quad(i,3,j) = Zero
    Quad(i,4,j) = Zero
    Quad(i,5,j) = Zero
    Quad(i,6,j) = Zero
  end do
end do

! MulticenterMultipoleExpansion

! First get the expansion in AO-format.

call mma_allocate(iCent,nBas(1),nBas(1),label='iCent')
call GetMem('Dummy','Allo','Inte',iDum,nBas(1)**2)
call MultiNew(iQ_Atoms,nBas(1),nOcc,natyp,nntyp,iMME,iWork(iDum),iCent,nMlt,outxyz,SlExpQ,lSlater)
call GetMem('Dummy','Free','Inte',iDum,nBas(1)**2)

! If MP2 density correction is requested, go here. This option is
! not working nicely, alas, hence do not use!

if (Mp2DensCorr) then
  call Mbpt2Corr(nBas(1),Cmo)
end if

! Then transform this to MO-format since that is the basis we use here!

nTyp = 0
do i=1,nMlt
  nTyp = nTyp+nTri_Elem(i)
end do
call GetMem('OnTheWay','Allo','Real',ipO,nTyp)
kaunter = 0
do iO1=1,iOrb(1)
  do iO2=1,iO1
    kaunter = kaunter+1
    do iB1=1,nBas(1)
      do iB2=1,nBas(1)
        kaunta = iCent(iB2,iB1)
        indMME = iTri(iB1,iB2)
        cProd = Cmo(iB1+(iO1-1)*nBas(1))*Cmo(iB2+(iO2-1)*nBas(1))
        do iTyp=1,nTyp
          Work(ipO+iTyp-1) = cProd*Work(iMME(iTyp)+indMME-1)
        end do
        Cha(kaunter,kaunta) = Cha(kaunter,kaunta)+Work(ipO)
        DipMy(kaunter,1,kaunta) = DipMy(kaunter,1,kaunta)+Work(ipO+1)
        DipMy(kaunter,2,kaunta) = DipMy(kaunter,2,kaunta)+Work(ipO+2)
        DipMy(kaunter,3,kaunta) = DipMy(kaunter,3,kaunta)+Work(ipO+3)
        Quad(kaunter,1,kaunta) = Quad(kaunter,1,kaunta)+Work(ipO+4)
        Quad(kaunter,2,kaunta) = Quad(kaunter,2,kaunta)+Work(ipO+5)
        Quad(kaunter,3,kaunta) = Quad(kaunter,3,kaunta)+Work(ipO+7) !Why seven? See @details
        Quad(kaunter,4,kaunta) = Quad(kaunter,4,kaunta)+Work(ipO+6)
        Quad(kaunter,5,kaunta) = Quad(kaunter,5,kaunta)+Work(ipO+8)
        Quad(kaunter,6,kaunta) = Quad(kaunter,6,kaunta)+Work(ipO+9)
      end do
    end do
  end do
end do
call GetMem('OnTheWay','Free','Real',ipO,nTyp)
call mma_deallocate(iCent)

! Deallocate the AO-multipoles.

do i=1,nTyp
  write(ChCo,'(I2.2)') i
  write(MMElab,*) 'MME'//ChCo
  call GetMem(MMElab,'Free','Real',iMME(i),nTri_Elem(nBas(1)))
end do

! Put quadrupoles in Buckingham form.

kaunter = 0
do i1=1,iOrb(1)
  do i2=1,i1
    kaunter = kaunter+1
    do k=1,iCi
      do j=1,6
        Quad(kaunter,j,k) = Quad(kaunter,j,k)*OneHalf
      end do
      Tra = Quad(kaunter,1,k)+Quad(kaunter,3,k)+Quad(kaunter,6,k)
      Tra = Tra/Three
      Quad(kaunter,1,k) = Quad(kaunter,1,k)-Tra
      Quad(kaunter,3,k) = Quad(kaunter,3,k)-Tra
      Quad(kaunter,6,k) = Quad(kaunter,6,k)-Tra
    end do
  end do
end do
!----------------------------------------------------------------------*
! To conclude, we sum up. This serves two purposes: (1) to make a check*
! if things have proceeded nicely, (2) to deduce if the QM-molecule is *
! charged.                                                             *
!----------------------------------------------------------------------*
qEl = Zero
dipx = Zero
dipy = Zero
dipz = Zero
dipx0 = Zero
dipy0 = Zero
dipz0 = Zero
qtot = Zero
dTox = Zero
dToy = Zero
dToz = Zero
call GetMem('TotMME','Allo','Real',iMtot,10*nTri_Elem(iQ_Atoms))
call dCopy_(10*nTri_Elem(iQ_Atoms),[Zero],0,Work(iMtot),1)
do ii=1,iOrb(1)
  i = ntri_Elem(ii)
  do j=1,nTri_Elem(iQ_Atoms)
    qEl = qEl+Cha(i,j)*Occu(ii)
    Work(iMtot+10*(j-1)) = Work(iMtot+10*(j-1))+Cha(i,j)*Occu(ii)
    dipx = dipx+DipMy(i,1,j)*Occu(ii)
    dipy = dipy+DipMy(i,2,j)*Occu(ii)
    dipz = dipz+DipMy(i,3,j)*Occu(ii)
    Work(iMtot+10*(j-1)+1) = Work(iMtot+10*(j-1)+1)+DipMy(i,1,j)*Occu(ii)
    Work(iMtot+10*(j-1)+2) = Work(iMtot+10*(j-1)+2)+DipMy(i,2,j)*Occu(ii)
    Work(iMtot+10*(j-1)+3) = Work(iMtot+10*(j-1)+3)+DipMy(i,3,j)*Occu(ii)
    dipx0 = dipx0+Cha(i,j)*outxyz(j,1)*Occu(ii)
    dipy0 = dipy0+Cha(i,j)*outxyz(j,2)*Occu(ii)
    dipz0 = dipz0+Cha(i,j)*outxyz(j,3)*Occu(ii)
    Work(iMtot+10*(j-1)+4) = Work(iMtot+10*(j-1)+4)+Quad(i,1,j)*Occu(ii)
    Work(iMtot+10*(j-1)+5) = Work(iMtot+10*(j-1)+5)+Quad(i,2,j)*Occu(ii)
    Work(iMtot+10*(j-1)+6) = Work(iMtot+10*(j-1)+6)+Quad(i,3,j)*Occu(ii)
    Work(iMtot+10*(j-1)+7) = Work(iMtot+10*(j-1)+7)+Quad(i,4,j)*Occu(ii)
    Work(iMtot+10*(j-1)+8) = Work(iMtot+10*(j-1)+8)+Quad(i,5,j)*Occu(ii)
    Work(iMtot+10*(j-1)+9) = Work(iMtot+10*(j-1)+9)+Quad(i,6,j)*Occu(ii)
  end do
end do
if (iPrint >= 10) then
  write(u6,*)
  write(u6,*) '    Distributed multipole in each centre'
  write(u6,*) '    (Compare with output from MpProp.)'
  do j=1,nTri_Elem(iQ_Atoms)
    if (j <= iQ_Atoms) then
      Work(iMtot+10*(j-1)) = Work(iMtot+10*(j-1))-ChaNuc(j)
    end if
    write(u6,*) '      Center: ',j
    write(u6,*) '      Charge: ',-Work(iMtot+10*(j-1))
    write(u6,*) '      Dipole: ',(-Work(iMtot+10*(j-1)+kk),kk=1,3)
    write(u6,*)
  end do
end if
call GetMem('TotMME','Free','Real',iMtot,10*nTri_Elem(iQ_Atoms))
do i=1,iQ_Atoms
  qtot = qtot+ChaNuc(i)
  dTox = dTox+ChaNuc(i)*outxyz(i,1)
  dToy = dToy+ChaNuc(i)*outxyz(i,2)
  dToz = dToz+ChaNuc(i)*outxyz(i,3)
end do
qtot = qtot-qEl
dTox = dTox-dipx-dipx0
dToy = dToy-dipy-dipy0
dToz = dToz-dipz-dipz0
if (iPrint >= 5) then
  write(u6,*)
  write(u6,*) '    Summed multipoles for unperturbed w.f.'
  write(u6,*) '      Charge: ',qtot
  write(u6,*) '      Dipole: ',dTox,',',dToy,',',dToz
end if
if (abs(qtot) > 1.0e-4_wp) ChargedQM = .true.

! Make a check of the one-electron matrix: is it pure?

call Chk_OneHam(nBas)

! So read integrals and construct super-matrix.

call ScfH0(nBas)

! The End...

return

end subroutine ScfHandM
