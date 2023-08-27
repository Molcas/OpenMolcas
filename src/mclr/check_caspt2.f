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

      subroutine check_caspt2(mode)
      use MckDat, only: sNew
      use Arrays, only: CMO, G2t, G1t

      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "sa.fh"
#include "warnings.h"

      Character*72 Line
      Character*8 Method
      character(len=128) :: FileName
      character(len=16) :: StdIn, mstate1
      Logical Found,Exists,NeedGrdt
C
C     Call RdInp_MCLR()  ! Read in input
C     write (*,*) "isnac = ", isnac

      iRlxRoot    = 0
      iRlxRootPT2 = 0
      call Get_iScalar('SA ready',iGo)
      !! Requested root for gradient
      Call Get_iScalar('Relax CASSCF root',iRlxRoot)
      if (iGo.eq.1) then
        !! CASPT2 density has been computed
        !! Check the root of the density
        Call Get_iScalar('Relax original root',iRlxRootPT2)
      end if
C     write (*,*) "iGo = ", iGo
C     write (*,*) "irlxroot = ", irlxroot
C     write (*,*) "irlxrootPT2 = ", irlxrootPT2
C
      !! Check this is NAC or not
      isNAC = .false.
      call Get_cArray('MCLR Root',mstate1,16)
      if (index(mstate1,'@') /= 0) then
        !! Requested root for NAC by ALASKA
        read(mstate1,'(1X,I7,1X,I7)') NACStates(1),NACStates(2)
        if (NACStates(1) /= 0) isNAC = .true.
        if (NACStates(1) == 0) iRlxRoot = NACStates(2)
      else if (iGo.eq.1 .and. iRlxRoot.ne.iRlxRootPT2) then
        !! this means CASPT2 density has been computed for the states
        !! specified by NAC in &CASPT2
        !! in this case, perform MCLR anyway(?)
        NACStates(1) = iRlxRoot
        NACStates(2) = iRlxRootPT2
        isNAC = .true.
        override = .true.
      end if

      if (mode.eq.1) return

C     write (*,*) "isnac = ", isnac
C     if (isnac) then
C       write (*,*) "requested NAC:", nacstates(1),nacstates(2)
C       write (*,*) "computed  NAC:", irlxroot,irlxrootpt2
C     else
C       write (*,*) "requested GRD:", irlxroot
C       write (*,*) "computed  GRD:", irlxrootpt2
C     endif

      !! If CASPT2 density has been computed, and the root of
      !! the density is the desired one in ALASKA, go for MCLR
      if (isNAC) then
        iRoot1req = MAX(NACStates(1),NACStates(2))
        iRoot2req = MIN(NACStates(1),NACStates(2))
        iRoot1com = MAX(iRlxRoot,iRlxRootPT2)
        iRoot2com = MIN(iRlxRoot,iRlxRootPT2)
        if (iGo.ne.0 .and. iRoot1req.eq.iRoot1com
     *               .and. iRoot2req.eq.iRoot2com) return
      else
        if (iGo.ne.0 .and. iRlxRoot.eq.iRlxRootPT2) return
      end if
      !! Otherwise, compute CASPT2 density for the state specified
      !! by ALASKA (iRlxRoot)

      !! This most likely means that density for the state
      !! has not been computed
      LuInput = 11
      LuInput = IsFreeUnit(LuInput)
      call StdIn_Name(StdIn)
      call Molcas_open(LuInput,StdIn)

      write(LuInput,'(A)') '>ECHO OFF'
      write(LuInput,'(A)') '>export MCLR_OLD_TRAP=$MOLCAS_TRAP'
      write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

      FileName = 'CASPTINP'
      call f_inquire(Filename,Exists)

      if (Exists) then
        LuSpool2 = 77
        LuSpool2 = IsFreeUnit(LuSpool2)
        call Molcas_Open(LuSpool2,Filename)

        NeedGrdt = isStructure() /= 1
        do
          read(LuSpool2,'(A)',iostat=istatus) Line
C         write (*,'(a)') line
          if (istatus > 0) call Abend()
          if (istatus < 0) exit
          write(LuInput,'(A)') Line
          Call UpCase(Line)
          if (Line(1:4).eq.'GRDT') NeedGrdt = .false.
        end do

        if (NeedGrdt) then
          backspace LuInput
          write (LuInput,'(A)') 'GRDT'
        end if

        close(LuSpool2)

      else

        write(6,'(A)') "this cannot happen, ig"
        call abend()

      end if

      FileName = 'MCLRINP'
      call f_inquire(Filename,Exists)

      if (Exists) then
        LuSpool2 = 77
        LuSpool2 = IsFreeUnit(LuSpool2)
        call Molcas_Open(LuSpool2,Filename)

        do
          read(LuSpool2,'(A)',iostat=istatus) Line
C         write (*,'(a)') line
          if (istatus > 0) call Abend()
          if (istatus < 0) exit
          write(LuInput,'(A)') Line
        end do

        close(LuSpool2)

      else

        write(LuInput,'(A)') ' &Mclr &End'
        write(LuInput,'(A)') 'End of Input'

      end if

!     write(LuInput,'(A)') '>RM -FORCE $Project.MckInt'
      write(LuInput,'(A)') '>export MOLCAS_TRAP=$MCLR_OLD_TRAP'
      write(LuInput,'(A)') '>ECHO ON'

C     rewind luinput
C     write (*,*)
C     write (*,*) "show luinput"
C     write (*,*)
C     do
C       read(LuInput,'(A)',iostat=istatus) Line
C       write (*,'(a)') line
C       if (istatus > 0) call Abend()
C       if (istatus < 0) exit
C     end do

      close(LuInput)

      call Finish(_RC_INVOKED_OTHER_MODULE_)

      return

      end subroutine check_caspt2
