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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************
Subroutine validate_wfn()
  use definitions,only:iwp,u6
  use rasscf_global,only:lRoots,NAC,NIN,nRoots,iRoot
  use general_data,only:nbas,ndel,nssh,nrs3,nrs2,nrs1,nash,nish,nfro,nbas,norb,stsym,ntot,nsym,nelec3, &
                         ispin,nrs3t,nrs2t,nactel,nhole1,nrs1t
  implicit none

#include "rasdim.fh"

  integer(kind=iwp) :: ierr,i,ia0,ib0,ic0,ierr1
  integer(kind=iwp) :: ierr2,iSym

  IERR = 0
  If(NTOT > mxOrb) Then
    write(u6,*)
    write(u6,*) '************ ERROR ******************'
    Call WarningMessage(2,'Too many orbitals.')
    write(u6,'(1X,A,I8)') 'Too many orbitals NTOT=',NTOT
    write(u6,'(1X,A,I8)') 'Limit is MXORB=',MXORB
    write(u6,*) '*************************************'
    call Quit_OnUserError()
  EndIf
  If(NAC > mxAct) then
    write(u6,*)
    write(u6,*) '*************** ERROR ***************'
    Call WarningMessage(2,'Too many active orbitals.')
    write(u6,'(1X,A,I8)') 'Too many active orbitals NAC=',NAC
    write(u6,'(1X,A,I8)') 'Limit is MXACT=',MXACT
    write(u6,*) '*************************************'
    call Quit_OnUserError()
  Endif
  If(NIN > mxIna) then
    write(u6,*)
    write(u6,*) '*************** ERROR ***************'
    Call WarningMessage(2,'Too many inactive orbitals.')
    write(u6,'(1X,A,I8)') 'Too many inactive orbitals NIN=',NIN
    write(u6,'(1X,A,I8)') 'Limit is MXINA=',MXINA
    write(u6,*) '*************************************'
    call Quit_OnUserError()
  Endif
  If(NACTEL > 2*NAC) then
    write(u6,*)
    write(u6,*) '********************* ERROR **********************'
    Call WarningMessage(2,'Too many active electrons.')
    write(u6,'(1X,A,I6)') 'Too many active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of active orbitals=',2*NAC
    write(u6,*) '**************************************************'
    call Quit_OnUserError()
  Endif
  If(NHOLE1 > 2*NRS1T) then
    write(u6,*)
    write(u6,*) '******************** WARNING *********************'
    Call WarningMessage(2,'Too many holes in Ras1.')
    write(u6,'(1X,A,I6)') 'You allow too many holes in Ras1 NHOLE1=',NHOLE1
    write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras1 orbitals=',2*NRS1T
    write(u6,*) '**************************************************'
    call Quit_OnUserError()
  Endif
  If(NELEC3 > 2*NRS3T) then
    write(u6,*)
    write(u6,*) '******************** WARNING *********************'
    Call WarningMessage(1,'Too many electrons in Ras3.')
    write(u6,'(1X,A,I6)') 'You allow too many electrons in Ras3 NELEC3=',NELEC3
    write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras3 orbitals=',2*NRS3T
    NELEC3 = 2*NRS3T
    write(u6,'(1X,A,I6)') 'NELEC3 has been reset to ',NELEC3
    write(u6,*) '**************************************************'
  Endif
  If(NACTEL > 2*(NRS1T+NRS2T)+NELEC3) then
    write(u6,*)
    write(u6,*) '********************* ERROR **********************'
    Call WarningMessage(2,'Too many active electrons.')
    write(u6,'(1X,A,I8)') 'Too many active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(u6,*) '**************************************************'
    call Quit_OnUserError()
  EndIf
  If(NACTEL < 2*NRS1T-NHOLE1) then
    write(u6,*)
    write(u6,*) '********************* ERROR **********************'
    Call WarningMessage(2,'Too few active electrons.')
    write(u6,'(1X,A,I8)') 'Too few active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(u6,*) '**************************************************'
    call Quit_OnUserError()
  EndIf

  If(NSYM /= 1 .and. NSYM /= 2 .and. NSYM /= 4 .and. NSYM /= 8) Then
    write(u6,*)
    Call WarningMessage(2,'Wrong nr of symmetries.')
    write(u6,*) '*************** ERROR ***************'
    write(u6,'(1X,A,I8)') 'Nr of symmetries NSYM=',NSYM
    write(u6,*) ' Only possible values are 1,2,4 or 8.'
    write(u6,*) '*************************************'
    call Quit_OnUserError()
  EndIf

  IERR1 = 0
  Do iSym = 1,nSym
    If(nBas(iSym) < 0) IERR1 = 1
    If(nFro(iSym) < 0) IERR1 = 1
    If(nDel(iSym) < 0) IERR1 = 1
    If(nOrb(iSym) < 0) IERR1 = 1
    If(nISh(iSym) < 0) IERR1 = 1
    If(nASh(iSym) < 0) IERR1 = 1
    If(nSSh(iSym) < 0) IERR1 = 1
    If(nRS1(iSym) < 0) IERR1 = 1
    If(nRS2(iSym) < 0) IERR1 = 1
    If(nRS3(iSym) < 0) IERR1 = 1
  EndDo
  IERR2 = 0
  Do iSym = 1,nSym
    If(nBas(iSym) > mxBas) IERR2 = 1
    If(nFro(iSym) > mxBas) IERR2 = 1
    If(nDel(iSym) > mxBas) IERR2 = 1
    If(nOrb(iSym) > mxBas) IERR2 = 1
    If(nISh(iSym) > mxBas) IERR2 = 1
    If(nASh(iSym) > mxBas) IERR2 = 1
    If(nSSh(iSym) > mxBas) IERR2 = 1
    If(nRS1(iSym) > mxBas) IERR2 = 1
    If(nRS2(iSym) > mxBas) IERR2 = 1
    If(nRS3(iSym) > mxBas) IERR2 = 1
  EndDo
  If(IERR1+IERR2 > 0) Then
    write(u6,*)
    write(u6,*) '****************** ERROR *******************'
    Call WarningMessage(2,'Erroneous nr of orbitals.')
    write(u6,*) 'Inappropriate nr of orbitals. One or more of'
    write(u6,*) 'these orbital counts is wrong or too large.'
    If(IERR1 > 0) write(u6,*) ' Negative values.'
    If(IERR2 > 0) write(u6,*) ' Extremely large values.'
    write(u6,'(1X,A,8I4)') '   All orbitals:',(NBAS(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '         Frozen:',(NFRO(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '       Inactive:',(NISH(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '         Active:',(NASH(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '          RAS-1:',(NRS1(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '          RAS-2:',(NRS2(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '          RAS-3:',(NRS3(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '      Secondary:',(NSSH(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') '        Deleted:',(NDEL(ISYM),ISYM=1,NSYM)
    write(u6,'(1X,A,8I4)') 'Basis functions:',(NBAS(ISYM),ISYM=1,NSYM)
    If(IERR1 > 0) Then
      write(u6,*) ' Have you used a too small basis set?'
    EndIf
    write(u6,*) '********************************************'
    Call Quit_OnUserError()
  EndIf

  IERR = 0
  IB0 = ISPIN-1
  IA0 = (NACTEL-IB0)/2
  IC0 = (NAC-IA0-IB0)
  If((2*IA0+IB0) /= NACTEL) IERR = 1
  If(IA0 < 0) IERR = 1
  If(IB0 < 0) IERR = 1
  If(IC0 < 0) IERR = 1
  If(IERR == 1) Then
    write(u6,*)
    write(u6,*) '******************** ERROR *********************'
    Call WarningMessage(2,'No such wave function.')
    write(u6,*) 'The following combined specifications are wrong.'
    write(u6,'(1X,A,8I4)') 'Nr of active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,8I4)') 'Nr of active orbitals     NAC=',NAC
    write(u6,'(1X,A,8I4)') 'Spin degeneracy         ISPIN=',ISPIN
    write(u6,*) 'There can be no such wave function.'
    write(u6,*) '************************************************'
    Call Quit_OnUserError()
  EndIf

  IERR = 0
  If(NROOTS > mxRoot) IERR = 1
  If(LROOTS > mxRoot) IERR = 1
  If(IERR == 1) Then
    write(u6,*)
    write(u6,*) '***************** ERROR *****************'
    Call WarningMessage(2,'Max roots exceeded.')
    write(u6,'(1X,A,I6)') 'Input Error: Max roots exceeded.',mxRoot
    write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
    write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
    write(u6,*) '************************************************'
    Call Quit_OnUserError()
  EndIf
  Do i = 1,NROOTS
    If(IROOT(i) < 0 .or. IROOT(i) > LROOTS) IERR = 1
  EndDo
  If(IERR == 1) Then
    write(u6,*)
    write(u6,*) '***************** ERROR *****************'
    Call WarningMessage(2,'Root specifications are wrong.')
    write(u6,*) 'CHKINP Error: Root specifications are wrong.'
    write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
    write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
    write(u6,*) '************************************************'
    call Quit_OnUserError()
  EndIf

  If(STSYM > NSYM) then
    write(u6,*)
    write(u6,*) '***************** ERROR *****************'
    Call WarningMessage(2,'Wrong symmetry.')
    write(u6,'(1X,A,I8)') 'CHKINP Error: Wrong symmetry.'
    write(u6,'(1X,A,I8)') 'State symmetry   STSYM=',STSYM
    write(u6,'(1X,A,I8)') 'Point group order NSYM=',NSYM
    write(u6,*) '************************************************'
    call Quit_OnUserError()
  EndIf

  IERR = 0
  If(NHOLE1 < 0 .and. NRS1T /= 0) IERR = 1
  If(NELEC3 < 0 .and. NRS3T /= 0) IERR = 1
  If(NACTEL < 0 .and. NRS2T /= 0) IERR = 1
  If(IERR == 1) Then
    write(u6,*)
    write(u6,*) '***************** ERROR *****************'
    Call WarningMessage(2,'Error in RAS specification.')
    write(u6,*) 'Error in RAS specification.'
    write(u6,'(1X,A,I8)') 'Max holes in Ras1,     NHOLE1=',NHOLE1
    write(u6,'(1X,A,I8)') 'Max electrons in Ras3, NELEC3=',NELEC3
    write(u6,'(1X,A,I8)') 'Nr of active electrons NACTEL=',NACTEL
    write(u6,*) '*****************************************'
    call Quit_OnUserError()
  EndIf
End
