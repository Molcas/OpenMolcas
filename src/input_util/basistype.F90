!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001-2005, Valera Veryazov                             *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************
      Subroutine BasisType(Filename,inline,BasisTypes)
!
!  translate basis file names
!
      character *(*) Filename
      Integer StrnLn
      External StrnLn
      character *256 DirName, OrigName
      character *256 Line, LineComp
      Character *4 tmp
      integer BasisTypes(4),irecl
      logical is_error,found
      character*16 CONT, ALLE, RELA, NUCL
      character*3 vCONT, vALLE, vRELA, vNUCL
#include "basistype.fh"
      CONT='#Contraction '
      ALLE='#AllElectron '
      RELA='#Hamiltonian '
      NUCL='#Nucleus '
      vCONT='UNK'
      vALLE='UNK'
      vRELA='UNK'
      vNUCL='UNK'
      isKnown=0
      If (inline.eq.1) Then
         BasisTypes(1)=-1
         BasisTypes(2)=-1
         BasisTypes(3)=-1
         BasisTypes(4)=-1
         call SysWarnMsg('BasisType','inline basis is used',            &
     & 'assuming all defaults for the basis types')
         Return
      End If
      BasisTypes(3)=0
      iunit=20
      iunit=isfreeunit(iunit)
      ileft=StrnLn(FileName)
      i=StrnLn(FileName)
      found=.false.

      do while((.not.found).and.(i.gt.1))
        if(FileName(i:i).eq.'/') then
           ileft=i
           found=.true.
        endif
        i=i-1
      enddo

      if(.not.found) then
        ileft=StrnLn(FileName)
        i=StrnLn(FileName)
        do while((.not.found).and.(i.gt.1))
           if(FileName(i:i).eq.'_') then
           ileft=i
           found=.true.
           endif
        i=i-1
        enddo
      endif
      DirName=Filename(1:ileft)
      i=index(Filename,' ')
      if(i.eq.0) i=len(Filename)
      OrigName=Filename(ileft+1:i)
      LenOrig=i-ileft
!
! first, check tbl.
      call molcas_open_ext2(iunit,DirName(1:ileft)//'basistype.tbl',    &
     & 'sequential','formatted',iostat,.false.,irecl,'old',is_error)
!      open(iunit,file=DirName(1:ileft)//'basistype.tbl',
!     &        Form='FORMATTED',IOSTAT=IOStat)
         if(Iostat.ne.0) then
         close(iunit)
          call molcas_open_ext2(iunit,'BASLIB_basistype.tbl',           &
     & 'sequential','formatted',iostat,.false.,irecl,'old',is_error)
!         open(iunit,file='BASLIB_basistype.tbl',
!     &        Form='FORMATTED',IOSTAT=IOStat)
        endif
        If (IOStat.ne.0) Then
          write(6,*) 'basistype.tbl is not found'
         close(iunit)
       goto 40
      endif
20    read(iunit,'(a)',end=30, err=30) Line
      if(Line(1:1).eq.'#') goto 20
      i=index(Line,OrigName(1:LenOrig-1))
      if(i.ne.1.or.Line(LenOrig:LenOrig).ne.' ') goto 20
       i=LenOrig+1
       iL=1
       iflag=0
       LineComp=' '
       Do jj=i,256
       if(Line(jj:jj).eq.' ') then
         if(iflag.eq.0) then
           iflag=1
           LineComp(iL:iL)=':'
           iL=iL+1
         endif
       else
         iflag=0
         LineComp(iL:iL)=Line(jj:jj)
         iL=iL+1
       endif
       enddo
       isKnown=1
       vCONT=LineComp(2:4)
       vALLE=LineComp(6:8)
       vRELA=LineComp(10:12)
       vNUCL=LineComp(14:16)
       if(vCONT.eq.' ') vCONT='UNK'
       if(vALLE.eq.' ') vALLE='UNK'
       if(vRELA.eq.' ') vRELA='UNK'
       if(vNUCL.eq.' ') vNUCL='UNK'
30     close(iunit)


! now let's check the header of the basis file
!
      if(OrigName(1:LenOrig-1).eq.' ') goto 40
      Call f_Inquire(DirName(1:ileft)//OrigName(1:LenOrig-1),Found)
      if(.not.Found) goto 40
!      print *, 'open >',DirName(1:ileft)//OrigName(1:LenOrig-1),'<'
      call molcas_open(iunit,DirName(1:ileft)//OrigName(1:LenOrig-1))
777   read(iunit,'(a)', end=35, err=35) Line
      if(Line(1:1).eq.'/') then
!  we already reached body of the file.
       LineComp=':'//vCONT//':'//vALLE//':'//vRELA//':'//vNUCL//':'
             isKnown=1
             goto 35
      endif
!     finish reading, let's create the return string

!       print *,'line=',Line
      if(index(Line,CONT(1:index(CONT,' '))).eq.1) then
         i=index(Line,' ')
         do j=i,len(Line)
           if(Line(j:j).ne.' ') goto 701
         enddo
701         vCONT=Line(j:j+3)
      endif
      if(index(Line,ALLE(1:index(ALLE,' '))).eq.1) then
         i=index(Line,' ')
         do j=i,len(Line)
           if(Line(j:j).ne.' ') goto 702
         enddo
702        vALLE=Line(j:j+3)
      endif
      if(index(Line,RELA(1:index(RELA,' '))).eq.1) then
         i=index(Line,' ')
         do j=i,len(Line)
           if(Line(j:j).ne.' ') goto 703
         enddo
703        vRELA=Line(j:j+3)
      endif

      if(index(Line,NUCL(1:index(NUCL,' '))).eq.1) then
         i=index(Line,' ')
         do j=i,len(Line)
           if(Line(j:j).ne.' ') goto 704
         enddo

704       vNUCL=Line(j:j+3)
      endif

      goto 777

!
! well. now use tbl
!
35     close(iunit)
40     if(isKnown.eq.0) LineComp=':UNK:UNK:UNK:UNK:'
!       print *,'DEBUG=',LineComp
          tmp=LineComp(2:5)
          i1=index(BasTypeCon,tmp)
          if(i1.eq.0.or.tmp.eq.'UNK:') Then
             BasisTypes(1)=-1
          Else
             i1=i1/4+1
             BasisTypes(1)=i1
          End If
!
          tmp=LineComp(6:9)
          i1=index(BasTypeAll,tmp)
          if(i1.eq.0.or.tmp.eq.'UNK:') Then
             BasisTypes(2)=-1
          Else
             i1=i1/4+1
             BasisTypes(2)=i1
! hack to map YES to AE_, and NO_ to NAE
             if(BasisTypes(2).eq.3) BasisTypes(2)=1
             if(BasisTypes(2).eq.4) BasisTypes(2)=2
          EndIf
!
          tmp=LineComp(10:13)
          i1=index(BasTypeRel,tmp)
          if(i1.eq.0.or.tmp.eq.'UNK:') Then
             BasisTypes(3)=-1
          Else
             i1=i1/4+1
             BasisTypes(3)=i1
          EndIf
!
          tmp=LineComp(14:17)
          i1=index(BasTypeNuc,tmp)
          if(i1.eq.0.or.tmp.eq.'UNK:') Then
             BasisTypes(4)=-1
          Else
             i1=i1/4+1
             BasisTypes(4)=i1
          EndIf

       return
      end
