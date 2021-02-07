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

subroutine RdInp_Dynamix(LuSpool,Task,nTasks,mTasks)

#ifdef _HDF5_
use mh5, only: mh5_put_dset
#endif

implicit real*8(a-h,o-z)
#include "MD.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
integer Task(nTasks), maxHop
integer VelVer, VV_First, VV_Second, Gromacs, VV_Dump
parameter(VelVer=1,VV_First=2,VV_Second=3,Gromacs=4,VV_Dump=5)
#ifdef _HDF5_
real*8, allocatable :: Mass(:)
#endif
logical lHop
character Title*72
character*180 Key, Line
character*180 Get_Ln
external Get_Ln

mTasks = 0

! Start of input

rewind(LuSpool)
call RdNLst(LuSpool,'Dynamix')

999 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)
!
if (Line(1:4) == 'TITL') goto 1100
if (Line(1:4) == 'PRIN') goto 1101
if (Line(1:4) == 'VV_F') goto 1102
if (Line(1:4) == 'VV_S') goto 1103
if (Line(1:4) == 'VV_D') goto 1104
if (Line(1:4) == 'THER') goto 1105
if (Line(1:4) == 'VELO') goto 1106
if (Line(1:2) == 'DT') goto 1107
if (Line(1:4) == 'GROM') goto 1108
if (Line(1:4) == 'TIME') goto 1109
if (Line(1:4) == 'VELV') goto 1110
if (Line(1:3) == 'HOP') goto 1111
if (Line(1:4) == 'REST') goto 1112
if (Line(1:4) == 'TEMP') goto 1113
if (Line(1:4) == 'ISOT') goto 1114
if (Line(1:4) == 'H5RE') goto 1115
if (Line(1:3) == 'OUT') goto 1116
if (Line(1:2) == 'IN') goto 1117
if (Line(1:3) == 'END') goto 9000

!>>>>>>>>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1100 continue
Line = Get_Ln(LuSpool)
call Get_S(1,Title,72)
goto 999
!>>>>>>>>>>>>>>>>>>>> PRIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1101 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,iPrint)
!>>>>>>>>>>>>>>>>>>>> VV_First <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1102 continue
write(6,*) ' VV_First 1'
mTasks = mTasks+1
Task(mTasks) = VV_First
write(6,*) ' VV_First 2'
goto 999
!>>>>>>>>>>>>>>>>>>>> VV_Second <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1103 continue
mTasks = mTasks+1
Task(mTasks) = VV_Second
goto 999
!>>>>>>>>>>>>>>>>>>>> VV_Dump <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1104 continue
mTasks = mTasks+1
Task(mTasks) = VV_Dump
goto 999
!>>>>>>>>>>>>>>>>>>>> THERmostat <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1105 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading THERMO.'
#endif
Line = Get_Ln(LuSpool)
call Get_I1(1,THERMO)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading THERMO.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> VELOcities <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1106 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading VELO.'
#endif
Line = Get_Ln(LuSpool)
call Get_I1(1,VELO)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading VELO.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> DT   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1107 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading DT.'
#endif
Line = Get_Ln(LuSpool)
call Get_F1(1,DT)
call Put_dScalar('Timestep',DT)
#ifdef _HDF5_
call mh5_put_dset(dyn_dt,DT)
#endif
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading DT.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> GROM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1108 continue
mTasks = mTasks+1
Task(mTasks) = Gromacs
goto 999
!>>>>>>>>>>>>>>>>>>>> TIME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1109 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,TIME)
goto 999
!>>>>>>>>>>>>>>>>>>>> VelVer <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1110 continue
!     This is the keyword for Velocity Verlet algorithm
mTasks = mTasks+1
Task(mTasks) = VelVer
goto 999
!>>>>>>>>>>>>>>>>>>>> Hop    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1111 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,maxHop)
lHop = .false.
call qpg_iScalar('MaxHops',lHop)
if (.not. lHop) then
  call Put_iScalar('MaxHops',maxHop)
end if
#ifdef _DEBUGPRINT_
write(6,*) ' lHop = ',lHop,'maxHop = ',maxHop
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> Restart <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1112 continue
Line = Get_Ln(LuSpool)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading RESTART.'
#endif
call Get_F1(1,RESTART)
goto 999
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading RESTART.'
#endif

!>>>>>>>>>>>>>>>>>>>> TEMPERATURE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1113 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading Temperature.'
#endif
Line = Get_Ln(LuSpool)
call Get_F1(1,TEMP)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading Temperature.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> Isotope <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1114 continue
write(6,*) 'ISOTope keyword is obsolete in DYNAMIX,'
write(6,*) 'use it in GATEWAY to specify isotopes/masses'
call Abend()
goto 999
!      Write (6,*) 'Unknown keyword:', Key
!      CALL Abend()
!>>>>>>>>>>>>>>>>>>>> Restart from HDF5 file <<<<<<<<<<<<<<<<<<<<<<<<<
1115 continue
#ifdef _HDF5_
lH5Restart = .true.
Line = Get_Ln(LuSpool)
call Get_S(1,FILE_H5RES,1)
#else
write(6,*) 'The user asks to restart the dynamics calculation '
write(6,*) 'from a HDF5 file, but this is not supported in this'
write(6,*) 'installation.'
call Quit_OnUserError()
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> project OUT some coordinates <<<<<<<<<<<<<<<<<<<<<<<
1116 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading OUT.'
#endif
Line = Get_Ln(LuSpool)
call Get_I1(1,POUT)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading OUT.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> keep IN only some coordinates <<<<<<<<<<<<<<<<<<<<<<<
1117 continue
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix starts reading IN.'
#endif
Line = Get_Ln(LuSpool)
call Get_I1(1,PIN)
#ifdef _DEBUGPRINT_
write(6,*) ' Dynamix ends reading IN.'
#endif
goto 999
!>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
9000 continue
write(6,*)

#ifdef _HDF5_
call Get_nAtoms_All(natom)
call mma_allocate(Mass,natom)
call Get_Mass_All(Mass,natom)
call mh5_put_dset(dyn_mass,Mass)
call mma_deallocate(Mass)
#endif

return

end subroutine RdInp_Dynamix
