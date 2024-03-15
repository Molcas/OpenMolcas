************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2023, Yoshio Nishimoto                                 *
************************************************************************

      subroutine check_caspt2(mode)
C
C     Check the roots to be considered in CASPT2 gradient
C
C     With mode = 0, this subroutine should first try to decide the root
C     for which CASPT2 density and MCLR are performed. If we do not have
C     CASPT2 density, i.e. iGo /= 3, call CASPT2 using the root
C     specified by ALASKA (or 'Relax CASSCF root'). Othewise (iGo = 3),
C     we are going to perform MCLR next. If ALASKA has not specified
C     roots, perform MCLR for the roots specified by CASPT2 ('Relax
C     original root' for gradients, or that and 'Relax CASSCF root' for
C     NAC). If ALASKA has specified, the roots are obtained from 'MCLR
C     Root', and CASPT2 is then called to compute the density etc.
C
C     With mode = 1, this subroutine obtains the character in 'MCLR
C     Root' and determine the roots for NAC calculation.
C
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
      character(len=128) :: FileName
      character(len=16) :: StdIn, mstate1
      Logical Exists,NeedGrdt
C
      iRlxRoot    = 0
      iRlxRootPT2 = 0
      call Get_iScalar('SA ready',iGo)
      !! Requested root for gradient
      Call Get_iScalar('Relax CASSCF root',iRlxRoot)
      if (iGo == 3) then
        !! iGo=3 means that CASPT2 density has been computed
        !! Check the root of the density
        Call Get_iScalar('Relax original root',iRlxRootPT2)
      end if
C
      !! Check this is NAC or not
      isNAC = .false.
      call Get_cArray('MCLR Root',mstate1,16)
      if (index(mstate1,'@') /= 0) then
        !! Requested root for NAC by ALASKA
        read(mstate1,'(1X,I7,1X,I7)') NACStates(1),NACStates(2)
        if (NACStates(1) /= 0) isNAC = .true.
        if (NACStates(1) == 0) iRlxRoot = NACStates(2)
      else if ((iGo == 3) .and. (iRlxRoot /= iRlxRootPT2)) then
        !! This means CASPT2 density has been computed for the states
        !! specified by the NAC option in &CASPT2 (either specified by
        !! the original input or the call below).
        !! In this case, perform MCLR anyway(?)
        !! The states can be different from those ALASKA requests.
        !! If different, ALASKA will call MCLR then CASPT2 again with
        !! the correct states.
        NACStates(1) = iRlxRoot
        NACStates(2) = iRlxRootPT2
        isNAC = .true.
        override = .true.
      end if

      !! With mode = 1, just set NACStates
      if (mode == 1) return

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
        if ((iGo /= 0) .and. (iRoot1req == iRoot1com)
     *                 .and. (iRoot2req == iRoot2com))   return
      else
        if ((iGo /= 0) .and. (iRlxRoot  == iRlxRootPT2)) return
      end if

      !! Otherwise, compute CASPT2 density for the state specified
      !! by ALASKA (iRlxRoot)

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

        NeedGrdt = (isStructure() /= 1)
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
        write(6,'(A)') "CASPT2 gradient without &CASPT2?"
        write(6,'(A)') "this cannot happen, ig"
        call abend()
      end if

      FileName = 'MCLRINP'
      call f_inquire(Filename,Exists)

      !! NAC states are obtained from "MCLR Roots"
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
