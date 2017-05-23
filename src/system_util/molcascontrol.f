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
* Copyright (C) 2000-2016, Valera Veryazov                             *
************************************************************************
******************************************************************************
*                                                                            *
* Author:   Valera Veryazov 2000-2016                                        *
*           Theoretical Chemistry                                            *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
******************************************************************************

      Subroutine MolcasControl(Label,Value)
************************************************************
*
*   <DOC>
*     <Name>MolcasControl</Name>
*     <Syntax>Call MolcasControl(Label,Value)</Syntax>
*     <Arguments>
*       \Argument{Label}{Query string}{Character}{in}
*       \Argument{Value}{Returned value}{Character}{out}
*     </Arguments>
*     <Purpose>Query a string from the control file</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V.Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>file molcas.control</Side_Effects>
*     <Description>
*
*      Only lines started from ! are in use.
*
*       if user modified molcas.control file
*       (by placing $!$ instead of $\#$)
*       return a value of the field (as a string)
*       and mark the label as a comment
*       else - return blank value
*
*
*       Usage:
*            Call MolcasControl('SHUTDOWN',Value)
*            if(Value.eq.'YES') Call abend
*
*     </Description>
*    </DOC>
*
************************************************************
      parameter(nLines=20)
      character filename*16
      character*80 Line(nLines)
      character*(*) Label
      character*(*) Value
      Logical Exist
      Logical Modify
      Integer StrnLn

      filename='molcas.control'
      Value=' '
      Modify=.false.
      Value=' '
      call f_inquire(filename,Exist)
      if(.Not.Exist) return
      Lu=1
      Lu=isfreeunit(Lu)
      open(Lu,File=filename)
      iLine=1
1     read(Lu,'(a)',end=100,err=100) Line(iLine)
      if(Line(iLine)(1:1).eq.'!') Modify=.true.
      iLine=iLine+1
      if(iLine.lt.nLines) goto 1

100   continue
      close(Lu)
c

      if(.not.Modify) return


      open(Lu,File=filename)
      do ic=1,iLine-1
      if(Line(ic)(1:1).eq.'!') then
        i=index(Line(ic)(2:),'=')
        if(i.gt.0) then
          if(Line(ic)(2:i).eq.Label) then
             Line(ic)(1:1)='#'
             Modify=.true.
             Value=Line(ic)(i+2:)
          endif
        endif
      endif
      i=StrnLn(Line(ic))
      write(Lu,'(a)') Line(ic)(1:i)
      enddo
      close(Lu)

      Return
      End
c
      Subroutine MolcasControlInit(label)
************************************************************
*
*   <DOC>
*     <Name>MolcasControlInit</Name>
*     <Syntax>Call MolcasControlInit(label)</Syntax>
*     <Arguments>
*       \Argument{label}{Coma separated list of fields}{Character}{in}
*     </Arguments>
*     <Purpose>Initiate molcas.control file</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V.Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Create a dummy file molcas.control
*        Not more than 20 entries are permitted
*        Initial value (as well as = sign) can bi omitted
*        Usage:
*            Call MolcasControlInit('SHUTDOWN=YES,ITER')
*     </Description>
*    </DOC>
*
************************************************************
      parameter(nLines=20)
      character*(*) Label
      character*512 tmp
      character*32 My
      character filename*16
      Integer StrnLn


      filename='molcas.control'
      iC=0
      tmp=Label(1:len(Label))
      Lu=1
      Lu=isfreeunit(Lu)
      open(Lu,File=filename)
      write(Lu,'(a)')'# Molcas control file: change # to ! to activate.'
      islast=0

10    i=index(tmp,',')
      My=' '
      if(i.gt.0) then
           My(1:i-1)=tmp(1:i-1)
           tmp=tmp(i+1:)
      else
           My=trim(tmp)
           islast=1
      endif
      iC=iC+1
      if(ic.gt.nLines) call abend()
      i=StrnLn(My)
      if(index(My,'=').eq.0) then
       i=i+1
       My(i:i)='='
      endif
      write(Lu,'(a,a)') '#',My(1:i)
      if(islast.eq.0) goto 10
      close(Lu)
      Return
      end
