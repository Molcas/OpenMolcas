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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
      Subroutine Add_Info(Label,Array,nArray,iToll)
************************************************************************
*                                                                      *
*     written by:                                                      *
*     V.Veryazov                                                       *
*                                                                      *
*    parameters:                                                       *
*       Label - text label                                             *
*       Array(nArray) - values to check                                *
*       iToll - tolerance:                                             *
*                if positive:10**(-iToll)                              *
*                else                                                  *
*                           : -iToll                                   *
*                                                                      *
************************************************************************
*
      Character*(*) Label
      Character*120 Line
      Real*8 Array(nArray)
      Character*32 File_Name
c     Logical Exist
      Character*24 Junk
      Character*8 Toll
c     Integer irecl
c     Logical is_error
      Character*256 STRING
      Character*256 STMP
      Character*256  STRING2
      Character*256 Collect
      Logical Found
*--- Some variables for Geo environment //Jonas 2011
      Integer iGeoInfo(2)
      Character*15 Energy_File
      Character*13 GeoDataF
      Integer iGeoData,iuGeoData
      Integer nIntCoord
#include "para_info.fh"
*------------------------------------------------
c If this is a fake parallel run (e.g. inside the parallel loop of CASPT2_gradient,
c then do not add info - just return immidiately
#ifdef _MOLCAS_MPP_
      If ((.not.King()).and.(.not.Is_Real_Par())) Return
#endif
c Number - is a number of exported variables from an array.
      Number=20
      File_Name='molcas_info'
      if(iToll.eq.0) iToll=8
*
*     Call qEnter('Add_Info')
      Lu_Info=99
*
*---------------------------------------------------------------------*
*     Check the file status                                           *
*---------------------------------------------------------------------*
      call open_molcas_info
c      call f_Inquire (File_name, Exist)
c*----------------------------------------------------------------------*
c*     Open existing file and position the record pointer to the end    *
c*----------------------------------------------------------------------*
c      If ( Exist ) then
c      call molcas_open_Ext2(Lu_info,file_name,'sequential','formatted',
c     &         ios,.false.,irecl,'unknown',is_error)
c         If ( ios.ne.0 ) Then
c            Write (6,*) 'Add_Info: can not create info file'
c            Write (6,*) 'Check file permissions!'
c            Call Abend
c         End If
c         nLines = 0
c         Do while (.True.)
c           Read (Lu_Info,'(A)',End=100,iostat=ios) Line
c           nLines = nLines+1
c         End Do
c100      Continue
c         Rewind(Lu_Info)
c         Do iLine = 1,nLines
c           Read (Lu_Info,'(A)') Line
c         End Do
c#ifdef NAGFOR
cC FIXME: ugly hack to make NAG compiler happy
c         close(Lu_Info)
c         open(Lu_Info,file=file_name,position='append')
c#endif
c*----------------------------------------------------------------------*
c*     Open new file                                                    *
c*----------------------------------------------------------------------*
c      Else
c      call molcas_open_Ext2(Lu_info,file_name,'sequential','formatted',
c     &         ios,.false.,irecl,'unknown',is_error)
c         If ( ios.ne.0 ) Then
c            Write (6,*) 'Add_Info: can not create a new file'
c            Write (6,*) 'Check file permissions!'
c            Call Abend
c         End If
c         Do 20 i=1,Len(Line)
c            Line(i:i)='#'
c20       Continue
c         Write (Lu_Info,'(A)') Line
c         Write (Lu_Info,'(A)') '# MOLCAS-Info_File Vers.No. 1.1'
c         Write (Lu_Info,'(A)') Line
c      End If
*----------------------------------------------------------------------*
*     Append new information                                           *
*----------------------------------------------------------------------*
      write(Toll,'(i8)')iToll
      nlabel=Len(label)
      Line=label
      n0=nlabel
      Do i = 1, nlabel
         If (label(i:i).eq.' ') Line(i:i)='_'
      End Do
      call upcase(Line)
*------------------ Some code for Geo-Environment //Jonas 2011
*----------------------------------------------------------------------*
*     Check if we are in the GEO-Environment and append if energy      *
*----------------------------------------------------------------------*
c666   Continue
      Call qpg_iArray('GeoInfo',Found,length)
      If(Found) Then
         Call get_iArray('GeoInfo',iGeoInfo,2)
         If((nArray.eq.1) .and. (iGeoInfo(1) .eq. 1) .and.
     &      Label(1:2) .eq. 'E_') Then
            Write(Energy_File,'(A,I4.4)') 'disp.energy',iGeoInfo(2)
            LuDispEn = 1
            LuDispEn = isfreeunit(LuDispEn)
            Call Molcas_Open(LuDispEn,Energy_File)
            Write(LuDispEn,'(F16.8)') Array(nArray)
            Close(LuDispEn)
            GeoDataF='GEODATA'
            iGeoData = 0
            iuGeoData = 10
            iuGeoData = isFreeUnit(iuGeoData)
            Call DaName_WA(iuGeoData,GeoDataF)
            Call iDaFile(iuGeoData,2,nIntCoord,1,iGeoData)
            iGeoData=iGeoInfo(2)*(nIntCoord+1)+1
            Call dDaFile(iuGeoData,1,Array(nArray),1,iGeoData)
            Call DaClos(iuGeoData)
         End If
      End If
      If (MyRank.ne.0) GoTo 888
*-----------------------------------------------------------------------


C If this label should not be checked, then just return immediately:
      STRING=' '
      CALL GETENVF('MOLCAS_NOCHECK',STRING)
      call upcase(STRING)
      STMP=STRING
      isfound=0
999   icoma=index(STMP,',')
      if(icoma.ne.0) then
        STRING=' '
        STRING=STMP(1:icoma-1)
        STMP=STMP(icoma+1:)
      else
        STRING=STMP
        STMP=' '
      endif
      k=0
      l=len(STRING)
      do i=1,l
        if(STRING(i:i).eq.' ') then
          if(k.gt.0) then
            if(STRING2(1:k).eq.Line(1:k)) then
               isfound=1
               goto 898
            end if
          end if
          k=0
        else
          k=k+1
          STRING2(k:k)=STRING(i:i)
        end if
      end do
898   continue
      if(STMP.ne.' ') goto 999
      if(isfound.eq.1) goto 888

      Do iArray = 1, nArray
      nlabel=n0
       if(nArray.ne.1) then
         write(Junk,'(a,i3,a)') '[',iArray-1,']'
         do j=1,5
           if(Junk(j:j).ne.' ') then
           nlabel=nlabel+1
           Line(nlabel:nlabel)=Junk(j:j)
           endif
         enddo
       endif
         nlabel=nlabel+1
         Line(nlabel:nlabel)='='
         nlabel=nlabel+1
         Line(nlabel:nlabel)='"'
         ia=int(Array(iArray)+0.3d0)
         if(abs(Array(iArray)-ia).lt.0.0000001.and.ia.ne.0) then
          write(Junk,'(I24)')ia
          else
            if(abs(Array(iArray)).gt.1D-14)then
              write(Junk,'(F24.12)') Array(iArray)
            else
           Junk='0.0'
           endif
         endif
         do j=1,24
           if(Junk(j:j).ne.' ') then
           nlabel=nlabel+1
           Line(nlabel:nlabel)=Junk(j:j)
           endif
         enddo
         nlabel=nlabel+1
         Line(nlabel:nlabel)='"'
        if(iArray.lt.Number) then
*----------------------------------------------------------------------*
*     Export only head of Array                                        *
*----------------------------------------------------------------------*
c        write(Lu_Info,*)
        collect=Line(1:nlabel)
        call add_molcas_info(collect,nlabel)
c        write (Lu_Info,'(a)') Line(1:nlabel)
        if (iArray.eq.nArray) then
          collect='export '//Line(1:n0)
          call add_molcas_info(collect,n0+7)
c          write (Lu_Info,'(a,a)') 'export ',Line(1:n0)
        endif
        endif
      ik=0
      do i=1,8
      if(Toll(i:i).ne.' ') then
        ik=ik+1
        Junk(ik:ik)=Toll(i:i)
      endif
      enddo
      collect='#> '// Line(1:nlabel)// '/'//Junk(1:ik)
      call add_molcas_info(collect,nlabel+ik+4)
c      write (Lu_Info,'(a,a,a,a)') '#> ', Line(1:nlabel),'/',Junk(1:ik)
      enddo
*----------------------------------------------------------------------*
*     Close file                                                       *
*----------------------------------------------------------------------*
888   continue
      call close_molcas_info
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
*     Call qExit('Add_Info')
      Return
      End
