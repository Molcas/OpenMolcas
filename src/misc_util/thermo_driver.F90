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

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: UserT(64), UserP
real(kind=wp), intent(in) :: Eval(*)
integer(kind=iwp), intent(inout) :: nUserPT, nsRot
integer(kind=iwp), intent(in) :: nFreq
logical(kind=iwp), intent(in) :: lSlapaf ! If .true. then Thermo_Driver called by SLAPAF
integer(kind=iwp) :: iMult, nAtom, nSym
real(kind=wp) :: TotalM, TRotA, TRotB, TRotC
logical(kind=iwp) :: lTest

lTest = .false.

if (lSlapaf) then
  call Get_iScalar('NSYM',nSym)
  if (nSym /= 1) then
    write(u6,'(A)') 'WARNING: No thermochemistry analysis conducted for numerical frequencies unless no symmetry is used!'
    return
  end if
end if

write(u6,*)
call CollapseOutput(1,'Thermochemistry')
write(u6,*)
write(u6,'(1X,A)') '*********************'
write(u6,'(1X,A)') '*                   *'
write(u6,'(1X,A)') '*  THERMOCHEMISTRY  *'
write(u6,'(1X,A)') '*                   *'
write(u6,'(1X,A)') '*********************'
write(u6,*)

if (lTest) then
  write(u6,*) '----------------------------------------------------'
  write(u6,*) '[Thermo_Driver] Input Data:'
  write(u6,*) '    UserP=',UserP,'  nsRot=',nsRot,'nUserPT=',nUserPT
  write(u6,*) '    UserT(1-5)==',UserT(1:5)
  write(u6,'(A,I3,A,256F8.2)') '  nFreq=',nFreq,'  Freq(i)==',EVal(1:nFreq)
  write(u6,*) '----------------------------------------------------'
  call XFlush(u6)
end if

call Rotation(TotalM,TRotA,TRotB,TRotC,nsRot,nAtom,lSlapaf)
call Get_iScalar('Multiplicity',iMult)

if (lTest) then
  write(u6,*) ' Calling ThermoChem_,  iMult=',iMult
  write(u6,*) ' UserP=',UserP,'  nsRot=',nsRot,'  nAtom=',nAtom
  write(u6,*) ' TotalM,TRotA,TRotB,TRotC==',TotalM,TRotA,TRotB,TRotC
  write(u6,'(A,I3,A,256F8.2)') ' nFreq=',nFreq,'  Freq(i)==',EVal(1:nFreq)
  call XFlush(u6)
end if

call ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,nUserPT,nsRot,iMult,EVal,nFreq,lSlapaf)
call CollapseOutput(0,'Thermochemistry')

return

end subroutine Thermo_Driver
