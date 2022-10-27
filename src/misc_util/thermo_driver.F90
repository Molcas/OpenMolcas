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

subroutine Thermo_Driver(UserT,UserP,nUserPT,nsRot,EVal,nFreq,lSlapaf)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
integer nUserPT, nsRot, nFreq
real*8 UserT(64), UserP, Eval(*)
logical lSlapaf ! If .True. then Thermo_Driver called by SLAPAF
logical lTest

lTest = .false.

if (lSlapaf) then
  call Get_iScalar('NSYM',nSym)
  if (nSym /= 1) then
    write(6,'(A)') 'WARNING: No thermochemistry analysis conducted for numerical frequencies unless no symmetry is used!'
    return
  end if
end if

write(6,*)
call CollapseOutput(1,'Thermochemistry')
write(6,*)
write(6,'(1X,A)') '*********************'
write(6,'(1X,A)') '*                   *'
write(6,'(1X,A)') '*  THERMOCHEMISTRY  *'
write(6,'(1X,A)') '*                   *'
write(6,'(1X,A)') '*********************'
write(6,*)

if (lTest) then
  write(6,*) '----------------------------------------------------'
  write(6,*) '[Thermo_Driver] Input Data:'
  write(6,*) '    UserP=',UserP,'  nsRot=',nsRot,'nUserPT=',nUserPT
  write(6,*) '    UserT(1-5)==',(UserT(i),i=1,5)
  write(6,'(A,I3,A,256F8.2)') '  nFreq=',nFreq,'  Freq(i)==',(EVal(i),i=1,nFreq)
  write(6,*) '----------------------------------------------------'
  call XFlush(6)
end if

call Rotation(TotalM,TRotA,TRotB,TRotC,nsRot,nAtom,lSlapaf)
call Get_iScalar('Multiplicity',iMult)

if (lTest) then
  write(6,*) ' Calling ThermoChem,  iMult=',iMult
  write(6,*) ' UserP=',UserP,'  nsRot=',nsRot,'  nAtom=',nAtom
  write(6,*) ' TotalM,TRotA,TRotB,TRotC==',TotalM,TRotA,TRotB,TRotC
  write(6,'(A,I3,A,256F8.2)') ' nFreq=',nFreq,'  Freq(i)==',(EVal(i),i=1,nFreq)
  call XFlush(6)
end if

call ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,nUserPT,nsRot,iMult,EVal,nFreq,lSlapaf)
call CollapseOutput(0,'Thermochemistry')

return

end subroutine Thermo_Driver
