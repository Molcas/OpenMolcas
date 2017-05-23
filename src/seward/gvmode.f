************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GvMode(isGvMode)
c
c Make a decision about GV mode, based on MOLCAS_GV variable
c
c if MOLCAS_GV can be set to Yes/GuessOrb or SCF
c
c  Return:
c  -1 or 0 : normal run
c   1 - use GuessOrb
c   2 - use SCF
c
      character Value*256
      isGvMode=-1
      Value=' '
      Call getenvf('MOLCAS_GV',Value)
      if(Value.eq.' ') then
         isGvMode=0
         return
      endif
      if(Value(1:1).eq.'y'.or.Value(1:1).eq.'Y'.or.
     &   Value(1:1).eq.'g'.or.Value(1:1).eq.'G') then
      isGvMode=1
      else
        if(Value(1:1).eq.'s'.or.Value(1:1).eq.'S') then
          isGvMode=2
        else
          Write(6,*) 'Unknown value of MOLCAS_GV'
          isGvMode=-1
        endif
      endif
      return
      end
c
      Subroutine DoGvMode(isGvMode)
      character*16 StdIn
      if(isGvMode.eq.1) then
*                                                                      *
************************************************************************
*                                                                      *
      Write (6,*)
      Write (6,*) ' Seward requests the grid_it to be computed!'
      Write (6,*)
*
      LuInput=11
      LuInput=IsFreeUnit(LuInput)
      Call StdIn_Name(StdIn)
      Call Molcas_Open(LuInput,StdIn)
      Write (LuInput,'(A)') ' &Grid_It &End'
c
c vv: temporary fix - use ASCII instead of PACK
c      Write (LuInput,'(A)') 'Pack'

      Write (LuInput,'(A)') 'ASCII'
c VV
      Write (LuInput,'(A)') 'Spar'
      Write (LuInput,'(A)') 'End of Input'
      Write (LuInput,'(A)') '>> exit 0'
      Close(LuInput)
*                                                                      *
************************************************************************
*                                                                      *
      endif
      return
      end
