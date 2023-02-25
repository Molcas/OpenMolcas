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
use Dynamix_Globals, only: dyn_dt, dyn_mass, File_H5Res, lH5Restart
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Dynamix_Globals, only: DT, iPrint, PIN, POUT, RESTART, TEMP, THERMO, VELO, VelVer, VV_First, VV_Second, Gromacs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool, nTasks
integer(kind=iwp), intent(inout) :: Task(nTasks)
integer(kind=iwp), intent(out) :: mTasks
integer(kind=iwp) :: maxhop
real(kind=wp) :: TIME
logical(kind=iwp) :: lHop
character(len=72) :: Title
character(len=180) :: Key, Line
character(len=180), external :: Get_Ln
#ifdef _HDF5_
integer(kind=iwp) :: natom
real(kind=wp), allocatable :: Mass(:)
#endif

mTasks = 0

! Start of input

rewind(LuSpool)
call RdNLst(LuSpool,'Dynamix')

do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)

  if (Line(1:4) == 'TITL') then
    !>>>>>>>>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line = Get_Ln(LuSpool)
    call Get_S(1,Title,72)
  else if (Line(1:4) == 'PRIN') then
    !>>>>>>>>>>>>>>>>>>>> PRIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line = Get_Ln(LuSpool)
    call Get_I1(1,iPrint)
  else if (Line(1:4) == 'VV_F') then
    !>>>>>>>>>>>>>>>>>>>> VV_First <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    write(u6,*) ' VV_First 1'
    mTasks = mTasks+1
    Task(mTasks) = VV_First
    write(u6,*) ' VV_First 2'
  else if (Line(1:4) == 'VV_S') then
    !>>>>>>>>>>>>>>>>>>>> VV_Second <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    mTasks = mTasks+1
    Task(mTasks) = VV_Second
  else if (Line(1:4) == 'THER') then
    !>>>>>>>>>>>>>>>>>>>> THERmostat <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading THERMO.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_I1(1,THERMO)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading THERMO.'
#   endif
  else if (Line(1:4) == 'VELO') then
    !>>>>>>>>>>>>>>>>>>>> VELOcities <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading VELO.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_I1(1,VELO)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading VELO.'
#   endif
  else if (Line(1:2) == 'DT') then
    !>>>>>>>>>>>>>>>>>>>> DT   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading DT.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_F1(1,DT)
    call Put_dScalar('Timestep',DT)
#   ifdef _HDF5_
    call mh5_put_dset(dyn_dt,DT)
#   endif
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading DT.'
#   endif
  else if (Line(1:4) == 'GROM') then
    !>>>>>>>>>>>>>>>>>>>> GROM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    mTasks = mTasks+1
    Task(mTasks) = Gromacs
  else if (Line(1:4) == 'TIME') then
    !>>>>>>>>>>>>>>>>>>>> TIME <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line = Get_Ln(LuSpool)
    call Get_F1(1,TIME)
  else if (Line(1:4) == 'VELV') then
    !>>>>>>>>>>>>>>>>>>>> VelVer <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! This is the keyword for Velocity Verlet algorithm
    mTasks = mTasks+1
    Task(mTasks) = VelVer
  else if (Line(1:3) == 'HOP') then
    !>>>>>>>>>>>>>>>>>>>> Hop    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line = Get_Ln(LuSpool)
    call Get_I1(1,maxHop)
    lHop = .false.
    call qpg_iScalar('MaxHops',lHop)
    if (.not. lHop) then
      call Put_iScalar('MaxHops',maxHop)
    end if
#   ifdef _DEBUGPRINT_
    write(u6,*) ' lHop = ',lHop,'maxHop = ',maxHop
#   endif
  else if (Line(1:4) == 'REST') then
    !>>>>>>>>>>>>>>>>>>>> Restart <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Line = Get_Ln(LuSpool)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading RESTART.'
#   endif
    call Get_F1(1,RESTART)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading RESTART.'
#   endif
  else if (Line(1:4) == 'TEMP') then
    !>>>>>>>>>>>>>>>>>>>> TEMPERATURE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading Temperature.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_F1(1,TEMP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading Temperature.'
#   endif
  else if (Line(1:4) == 'ISOT') then
    !>>>>>>>>>>>>>>>>>>>> Isotope <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    write(u6,*) 'ISOTope keyword is obsolete in DYNAMIX,'
    write(u6,*) 'use it in GATEWAY to specify isotopes/masses'
    call Abend()
  else if (Line(1:4) == 'H5RE') then
    !>>>>>>>>>>>>>>>>>>>> Restart from HDF5 file <<<<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _HDF5_
    lH5Restart = .true.
    Line = Get_Ln(LuSpool)
    call Get_S(1,File_H5Res,1)
#   else
    write(u6,*) 'The user asks to restart the dynamics calculation '
    write(u6,*) 'from a HDF5 file, but this is not supported in this'
    write(u6,*) 'installation.'
    call Quit_OnUserError()
#   endif
  else if (Line(1:3) == 'OUT') then
    !>>>>>>>>>>>>>>>>>>>> project OUT some coordinates <<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading OUT.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_I1(1,POUT)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading OUT.'
#   endif
  else if (Line(1:2) == 'IN') then
    !>>>>>>>>>>>>>>>>>>>> keep IN only some coordinates <<<<<<<<<<<<<<<<<<<<<<<
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix starts reading IN.'
#   endif
    Line = Get_Ln(LuSpool)
    call Get_I1(1,PIN)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Dynamix ends reading IN.'
#   endif
  else if (Line(1:3) == 'END') then
    !>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    write(u6,*)
    exit
  else
    write(u6,*) 'Unknown keyword:',Key
    call Abend()
  end if
end do

#ifdef _HDF5_
call Get_nAtoms_All(natom)
call mma_allocate(Mass,natom)
call Get_Mass_All(Mass,natom)
call mh5_put_dset(dyn_mass,Mass)
call mma_deallocate(Mass)
#endif

return

end subroutine RdInp_Dynamix
