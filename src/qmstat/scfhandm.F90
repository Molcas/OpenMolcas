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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "WrkSpc.fh"
#include "tratoc.fh"
dimension Cmo(MxBas**2), Occu(MxBas), nOcc(MxBas), natyp(MxAt)
dimension nBas(MxSym), iCent(MxBas*MxBas)
dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6)
character MMElab*20, ChCo*2

!----------------------------------------------------------------------*
! Zeros.                                                               *
!----------------------------------------------------------------------*
iCi = iQ_Atoms*(iQ_Atoms+1)/2
do i=1,iOrb(1)*(iOrb(1)+1)/2
  do j=1,iCi
    Cha(i,j) = 0
    DipMy(i,1,j) = 0
    DipMy(i,2,j) = 0
    DipMy(i,3,j) = 0
    Quad(i,1,j) = 0
    Quad(i,2,j) = 0
    Quad(i,3,j) = 0
    Quad(i,4,j) = 0
    Quad(i,5,j) = 0
    Quad(i,6,j) = 0
  end do
end do

! MulticenterMultipoleExpansion

! First get the expansion in AO-format.

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
  nTyp = nTyp+i*(i+1)/2
end do
call GetMem('OnTheWay','Allo','Real',ipO,nTyp)
kaunter = 0
do iO1=1,iOrb(1)
  do iO2=1,iO1
    kaunter = kaunter+1
    kaunta = 0
    do iB1=1,nBas(1)
      do iB2=1,nBas(1)
        kaunta = kaunta+1
        iX1 = max(iB1,iB2)
        iX2 = min(iB1,iB2)
        indMME = iX2+iX1*(iX1-1)/2
        cProd = Cmo(iB1+(iO1-1)*nBas(1))*Cmo(iB2+(iO2-1)*nBas(1))
        do iTyp=1,nTyp
          Work(ipO+iTyp-1) = cProd*Work(iMME(iTyp)+indMME-1)
        end do
        Cha(kaunter,iCent(kaunta)) = Cha(kaunter,iCent(kaunta))+Work(ipO)
        DipMy(kaunter,1,iCent(kaunta)) = DipMy(kaunter,1,iCent(kaunta))+Work(ipO+1)
        DipMy(kaunter,2,iCent(kaunta)) = DipMy(kaunter,2,iCent(kaunta))+Work(ipO+2)
        DipMy(kaunter,3,iCent(kaunta)) = DipMy(kaunter,3,iCent(kaunta))+Work(ipO+3)
        Quad(kaunter,1,iCent(kaunta)) = Quad(kaunter,1,iCent(kaunta))+Work(ipO+4)
        Quad(kaunter,2,iCent(kaunta)) = Quad(kaunter,2,iCent(kaunta))+Work(ipO+5)
        Quad(kaunter,3,iCent(kaunta)) = Quad(kaunter,3,iCent(kaunta))+Work(ipO+7) !Why seven? See @details
        Quad(kaunter,4,iCent(kaunta)) = Quad(kaunter,4,iCent(kaunta))+Work(ipO+6)
        Quad(kaunter,5,iCent(kaunta)) = Quad(kaunter,5,iCent(kaunta))+Work(ipO+8)
        Quad(kaunter,6,iCent(kaunta)) = Quad(kaunter,6,iCent(kaunta))+Work(ipO+9)
      end do
    end do
  end do
end do
call GetMem('OnTheWay','Free','Real',ipO,nTyp)

! Deallocate the AO-multipoles.

do i=1,nTyp
  write(ChCo,'(I2.2)') i
  write(MMElab,*) 'MME'//ChCo
  call GetMem(MMElab,'Free','Real',iMME(i),nBas(1)*(nBas(1)+1)/2)
end do

! Put quadrupoles in Buckingham form.

kaunter = 0
do i1=1,iOrb(1)
  do i2=1,i1
    kaunter = kaunter+1
    do k=1,ici
      do j=1,6
        Quad(kaunter,j,k) = Quad(kaunter,j,k)*1.5
      end do
      Tra = Quad(kaunter,1,k)+Quad(kaunter,3,k)+Quad(kaunter,6,k)
      Tra = Tra/3
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
qEl = 0
dipx = 0
dipy = 0
dipz = 0
dipx0 = 0
dipy0 = 0
dipz0 = 0
qtot = 0
dTox = 0
dToy = 0
dToz = 0
call GetMem('TotMME','Allo','Real',iMtot,10*iQ_Atoms*(iQ_Atoms+1)/2)
call dCopy_(10*iQ_Atoms*(iQ_Atoms+1)/2,[0.0d0],0,Work(iMtot),1)
do ii=1,iOrb(1)
  i = ii*(ii+1)/2
  do j=1,iQ_Atoms*(iQ_Atoms+1)/2
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
  write(6,*)
  write(6,*) '    Distributed multipole in each centre'
  write(6,*) '    (Compare with output from MpProp.)'
  do j=1,iQ_Atoms*(iQ_Atoms+1)/2
    if (j <= iQ_Atoms) then
      Work(iMtot+10*(j-1)) = Work(iMtot+10*(j-1))-Chanuc(j)
    end if
    write(6,*) '      Center: ',j
    write(6,*) '      Charge: ',-Work(iMtot+10*(j-1))
    write(6,*) '      Dipole: ',(-Work(iMtot+10*(j-1)+kk),kk=1,3)
    write(6,*)
  end do
end if
call GetMem('TotMME','Free','Real',iMtot,4*iQ_Atoms*(iQ_Atoms+1)/2)
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
  write(6,*)
  write(6,*) '    Summed multipoles for unperturbed w.f.'
  write(6,*) '      Charge: ',qtot
  write(6,*) '      Dipole: ',dTox,',',dToy,',',dToz
end if
if (abs(qtot) > 0.0001) ChargedQM = .true.

! Make a check of the one-electron matrix: is it pure?

call Chk_OneHam(nBas)

! So read integrals and construct super-matrix.

call ScfH0(nBas)

! The End...

return

end subroutine ScfHandM
