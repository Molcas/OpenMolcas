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

subroutine rdinp_surfacehop()

use Tully_variables, only: tullyL, DECO, decoherence, NSUBSTEPS, Ethreshold, RandThreshold, fixedrandL, FixedRand, InitSeed, &
                           iseedL, tullySubVerb
#ifdef _HDF5_
use Surfacehop_globals, only: lH5Restart, File_H5Res
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
character(len=180) :: Key, Line
logical(kind=iwp) :: Found
integer(kind=iwp) :: LuSpool, NSTATE, i, j, ndata, maxHop
real(kind=wp), allocatable :: temp(:), AmatrixVR(:), AmatrixVI(:)
character(len=180), external :: Get_Ln
integer(kind=iwp), external :: IsFreeUnit

LuSpool = IsFreeUnit(21)
call SpoolInp(LuSpool)

call RdNLst(LuSpool,'surfacehop')

read_input: do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)

  select case(Line(1:4))
    case('TULL')
      tullyL=.true.
    case('DECO')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,DECO)
      decoherence = .true.
    case('SUBS')
      Line = Get_Ln(LuSpool)
      call Get_I1(1,NSUBSTEPS)
    case('ETHR')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,Ethreshold)
    case('RTHR')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,RandThreshold)
    case('PSUB')
      tullySubVerb = .true.
    case('DMTX')
      Line = Get_Ln(LuSpool)
      call Get_I1(1,NSTATE)
      call mma_allocate(AmatrixVR,NSTATE*NSTATE,label='AmatrixVR')
      call mma_allocate(AmatrixVI,NSTATE*NSTATE,label='AmatrixVI')
      call mma_allocate(temp,NSTATE,label='temp')
      do i=1,NSTATE
        Line = Get_Ln(LuSpool)
        call Get_F(1,temp,NSTATE)
        do j=1,NSTATE
          AmatrixVR((i-1)*NSTATE+j) = temp(j)
        end do
      end do
      do i=1,NSTATE
         Line = Get_Ln(LuSpool)
         call Get_F(1,temp,NSTATE)
         do j=1,NSTATE
            AmatrixVI((i-1)*NSTATE+j) = temp(j)
         end do
      end do
      call mma_deallocate(temp)
      call Qpg_zArray('AmatrixV',Found,ndata)
      if (.not.Found) then
        call Put_dArray('RAmatrixV',AmatrixVR,NSTATE*NSTATE)
        call Put_dArray('IAmatrixV',AmatrixVI,NSTATE*NSTATE)
      else if (ndata /= NSTATE*NSTATE) then
        call warningmessage(2,'The A matrix dimension on the input is different from the one in runfile')
        call Abend()
      end if
      call mma_deallocate(AmatrixVR)
      call mma_deallocate(AmatrixVI)
    case('FRAN')
      Line = Get_Ln(LuSpool)
      call Get_F1(1,FixedRand)
      fixedrandL = .true.
    case('ISEE') ! initial seed number
      Line = Get_Ln(LuSpool)
      call Get_I1(1,InitSeed)
      iseedL = .true.
    case('MAXH')
      Line = Get_Ln(LuSpool)
      call Get_I1(1,maxHop)
      call Put_iScalar('MaxHopsTully',maxHop)
!     write(u6,*) 'MaxHops set to ', maxHop
    case('H5RE')
#ifdef _HDF5_
      lH5Restart = .true.
      Line = Get_Ln(LuSpool)
      call Get_S(1,File_H5Res,1)
#else
      write(u6,*) 'The user asks to restart the dynamics calculation from a HDF5 file,'
      write(u6,*) 'but this is not supported in this installation.'
      call Quit_OnUserError()
#endif
    case('END ')
      exit
    case default
      write (u6,*) 'Unknown keyword: ', trim(Key)
      call Abend()
  end select
end do read_input

close(LuSpool)

return

end subroutine rdinp_surfacehop
