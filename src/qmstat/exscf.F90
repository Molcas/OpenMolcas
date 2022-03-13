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

subroutine ExScf(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms,nAtomsCC,Ax,Ay,Az,itri,Smat,SmatPure,InCutOff,ipAOSum)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
dimension Smat(MxOT), CorTemp(3), nCnC_C(MxBasC)
dimension SmatPure(MxOT)
dimension Inside(MxAt,3)
logical Inside, NearBy, InCutOff

!----------------------------------------------------------------------*
! Deduce how much the QM-molecule is translated from its position as   *
! definied in Seward.                                                  *
!----------------------------------------------------------------------*
Ax = Cordst(1,1)-outxyz(1,1)
Ay = Cordst(1,2)-outxyz(1,2)
Az = Cordst(1,3)-outxyz(1,3)
!----------------------------------------------------------------------*
! Rotate solvent orbitals and make AO integration.                     *
!----------------------------------------------------------------------*
Cut_ExSq1 = Cut_Ex1**2
Cut_ExSq2 = Cut_Ex2**2
nOrbSize = iOrb(1)*iOrb(2)
nV2size = iOrb(2)*nBaseC
nAObaseSize = nBaseQ*nBaseC
nStorlek = iOrb(1)*nBaseC
call GetMem('RotOrb','Allo','Real',iV2,nV2size)
call GetMem('Sint','Allo','Real',ipAOint,nAObaseSize)
call GetMem('Sintpar','Allo','Real',ipAOintpar,nAObaseSize)
call GetMem('OvlMO','Allo','Real',iOvlMO,nOrbSize)
call GetMem('Intermed','Allo','Real',iInte,nStorlek)
call GetMem('OvlMOpure','Allo','Real',iOPure,nOrbSize)
call GetMem('OvlMOene','Allo','Real',iOvlMOE,nOrbSize)
call GetMem('AUX','Allo','Real',ipAUX,iOrb(1)**2)
call GetMem('AUXp','Allo','Real',ipAUXp,iOrb(1)**2)
call GetMem('AUXtri','Allo','Real',ipAUXtri,iTri)
!Jose****************************************************************
if (lExtr(8)) then
  nAOqMOcl = nBaseQ*iOrb(2)
  iAOAOTri = nBaseQ*(nBaseQ+1)/2
  call GetMem('qAOclMOOvl','Allo','Real',iAOMOOvl,nAOqMOcl)
  call GetMem('qAOclMOOvlE','Allo','Real',iAOMOOvlE,nAOqMOcl)
  call GetMem('AuxAOp','Allo','Real',ipAOAUX,nBaseQ**2)
  call GetMem('AuxAOpTri','Allo','Real',ipAOAUXtri,iAOAOTri)
end if
!call GetMem('SumOvlAOQ','Allo','Real',ipAOSum,iAOAOTri)
!********************************************************************
InCutOff = .false.
do i=1,iTri
  Smat(i) = 0
  SmatPure(i) = 0
end do
nInsideCut = 0
do N=iCStart-1,nCent*(nPart-1),nCent

  ! Initialize.

  dist_sw = 1d20
  r3 = 1d20
  do i=1,MxAt
    Inside(i,1) = .false.
    Inside(i,2) = .false.
    Inside(i,3) = .false.
  end do
  NearBy = .false.
  ! Loop over atoms.
  do inwm=1,iQ_Atoms
    do k=1,3
      CorTemp(k) = (Cordst(N+1,k)-Cordst(inwm,k))**2
    end do
    r2 = CorTemp(1)+CorTemp(2)+CorTemp(3)
    dist_sw = min(dist_sw,r2)
    DH1 = 0.0d0 !Distances for the inner cut-off. Also include the hydrogens.
    DH2 = 0.0d0
    do k=1,3
      DH1 = DH1+(Cordst(N+2,k)-Cordst(inwm,k))**2
      DH2 = DH2+(Cordst(N+3,k)-Cordst(inwm,k))**2
    end do
    r3temp1 = min(DH1,DH2)
    r3temp2 = min(r3temp1,r2)
    r3 = min(r3,r3temp2)
    ! Check if this atom-atom pair inside
    if (r2 < Cut_ExSq1) then
      Inside(inwm,1) = .true.
      NearBy = .true.
    end if
    if (DH1 < Cut_ExSq1) then
      Inside(inwm,2) = .true.
      NearBy = .true.
    end if
    if (DH2 < Cut_ExSq1) then
      Inside(inwm,3) = .true.
      NearBy = .true.
    end if
  end do

  ! Now make the cut-off test.

  if (.not. NearBy) Go To 501
  if (r3 < Cut_ExSq2) then !Inner cut-off.
    InCutOff = .true.
  end if
  nInsideCut = nInsideCut+1

  ! Make the AO-AO overlap integration.

  call AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax,Inside)

  ! Transform to MO-MO overlap.

  call Dgemm_('T','N',iOrb(1),nBaseC,nBaseQ,ONE,Work(iV1),nBaseQ,Work(ipAOint),nBaseQ,ZERO,Work(iInte),iOrb(1))
  call Dgemm_('N','N',iOrb(1),iOrb(2),nBaseC,ONE,Work(iInte),iOrb(1),Work(iV2),nBaseC,ZERO,Work(iOvlMO),iOrb(1))
  call dcopy_(nOrbSize,[ZERO],iZERO,Work(iOvlMOE),iONE)
  do i=1,iOrb(2)
    ind = iOrb(1)*(i-1)
    call DaxPy_(iOrb(1),c_orbene(i),Work(iOvlMO+ind),iONE,Work(iOvlMOE+ind),iONE)
  end do

  !**Jose
  if (lExtr(8)) then
    call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,ONE,Work(ipAOint),nBaseQ,Work(iV2),nBaseC,ZERO,Work(iAOMOOvl),nBaseQ)
    call dcopy_(nAOqMOcl,[ZERO],iZERO,Work(iAOMOOvlE),iONE)
    do i=1,iOrb(2)
      ind = nBaseQ*(i-1)
      call DaxPy_(nBaseQ,c_orbene(i),Work(iAOMOOvl+ind),iONE,Work(iAOMOOvlE+ind),iONE)
    end do
  end if
  !*******************

  ! If you are interested, print some bla bla bla.

  if (iPrint >= 29) then
    write(6,*)
    write(6,*) 'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE ',N/nCent
    write(6,*) 'QM-MO  SOLV-MO  OVERLAP'
    call dcopy_(iOrb(1)*iOrb(2),Work(iOvlMO),iONE,Work(iOPure),iONE)
    do i=0,iOrb(1)-1
      do j=0,iOrb(2)-1
        write(6,8888) i+1,j+1,Work(iOPure+i+j*iOrb(1))
      end do
    end do
8888 format(I3,'    ',I3,'       ',F12.10)
  end if

  ! Construct the perturbation and accumulate pure overlap for
  ! subsequent construction of higher order term.

  call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),exrep2,Work(iOvlMO),iOrb(1),Work(iOvlMOE),iOrb(1),ZERO,Work(ipAUX),iOrb(1))
  call SqToTri_Q(Work(ipAUX),Work(ipAUXtri),iOrb(1))
  call DaxPy_(iTri,ONE,Work(ipAUXtri),iONE,Smat,iONE)
  call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),ONE,Work(iOvlMO),iOrb(1),Work(iOvlMO),iOrb(1),ZERO,Work(ipAUXp),iOrb(1))
  call SqToTri_Q(Work(ipAUXp),Work(ipAUXtri),iOrb(1))
  call DaxPy_(iTri,ONE,Work(ipAUXtri),iONE,SmatPure,iONE)

  !Jose*********************************
  if (lExtr(8)) then
    call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2,Work(iAOMOOvl),nBaseQ,Work(iAOMOOvlE),nBaseQ,ZERO,Work(ipAOAUX),nBaseQ)
    call SqToTri_Q(Work(ipAOAUX),Work(ipAOAUXtri),nBaseQ)
    call DaxPy_(iAOAOTri,ONE,Work(ipAOAUXtri),iONE,Work(ipAOSum),iONE)
  end if
  !*************************************

  ! This solvent molecule ends now!

501 continue
end do

call GetMem('RotOrb','Free','Real',iV2,nV2size)
call GetMem('Sint','Free','Real',ipAOint,nAObaseSize)
call GetMem('Sintpar','Free','Real',ipAOintpar,nAObaseSize)
call GetMem('OvlMO','Free','Real',iOvlMO,nOrbSize)
call GetMem('Intermed','Free','Real',iInte,nStorlek)
call GetMem('OvlMOpure','Free','Real',iOPure,nOrbSize)
call GetMem('OvlMOene','Free','Real',iOvlMOE,nOrbSize)
call GetMem('AUX','Free','Real',ipAUX,iOrb(1)**2)
call GetMem('AUXp','Free','Real',ipAUXp,iOrb(1)**2)
call GetMem('AUXtri','Free','Real',ipAUXtri,iTri)
!Jose****************************************************************
if (lExtr(8)) then
  call GetMem('qAOclMOOvl','Free','Real',iAOMOOvl,nAOqMOcl)
  call GetMem('qAOclMOOvlE','Free','Real',iAOMOOvlE,nAOqMOcl)
  call GetMem('AuxAOp','Free','Real',ipAOAUX,nBaseQ**2)
  call GetMem('AuxAOpTri','Free','Real',ipAOAUXtri,iAOAOTri)
end if
!********************************************************************

return

end subroutine ExScf
