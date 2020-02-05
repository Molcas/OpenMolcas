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
      Subroutine Put_ln(In_Line)
************************************************************************
* This function replaces function getln                                *
*                                                                      *
* It reads, broadcasts, and parses an input line                       *
*                                                                      *
* Blank lines or lines containing star (*) in column 1 are skipped     *
* Lines staring with exclaimation (!) in column 1 are skipped          *
*                                                                      *
* After this routine has been called, data can be retrieved using      *
* the subroutines                                                      *
*                                                                      *
*   Get_F(icol,array,n)  (for floating point values)                   *
*   Get_I(icol,iarry,n)  (for integer values)                          *
*   Get_S(icol,strgs,n)  (for character strings)                       *
*                                                                      *
* where icol is the first non-blank work to be taken, and n is the     *
* number of data.                                                      *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      Character*(*) In_line
      Character*180 line
* mxn should be len(line)/2+1
      parameter (mxn=91)
      common/cgetlc/ line
      common/cgetln/ ncol,istrt(mxn),iend(mxn)
      Line=In_Line
      l=len(line)
      Do i = 1, l
         if(ichar(line(i:i)).eq.9) line(i:i)=' '
         If (line(i:i).eq.';') line(i:l)=' '
      End Do
      ncol=0
      j=1
 10   icom=0
      do i=j,l
        if(line(i:i).eq.',') then
          icom=icom+1
          if(icom.ne.1) goto 20
        else
          if(line(i:i).ne.' ') goto 20
        end if
      end do
      Return
 20   do j=i,l
        if(line(j:j).eq.' '.or.line(j:j).eq.',') goto 30
      end do
      j=l+1
 30   ncol=ncol+1
      istrt(ncol)=i
      iend(ncol)=j-1
      goto 10
      End
      character*180 function get_ln(lunit)
      logical Quit_On_Error
      common /getlnQOE/ Quit_On_Error
      character*180 get_ln_quit
      character*4 keya(45)
      character*4 tmp
      integer i,ik
      common /ln_stack/keya
      common /ln_stack_count/ia,ik
      data ifirst/0/
      save ifirst
      if(ifirst.eq.0) then
c      print *,'DEBUG', 'first init'
        ifirst=1
        ia=0
	ik=0
      endif
      if(ia.eq.1.and.ik.gt.0) then
      
        tmp=keya(1)
        do i=1,ik-1
         keya(i)=keya(i+1)
        enddo
        keya(ik)=' '
        ik=ik-1
        if(ik.eq.0) ia=0
        get_ln=tmp

      else
      get_ln=get_ln_quit(lunit,1)
c      print *,'DEBUG', 'get line=',get_ln
      if(Quit_On_Error) Then
        Call WarningMessage(2,'Error in Get_Ln')
        Call Quit_OnUserError()
      End If
      endif
      Return
      End
c VV
      subroutine activate_ln_stack()
      common /ln_stack/keya
      common /ln_stack_count/ia,ik
      character*4 keya(45)
      integer ia,ik
      ia=1
      return
      end
c syntax: call set_ln_stack("xxx;yyy;zzz")
c no spaces are allowed
      subroutine set_ln_stack(keys)
      character *180 keys, keyst, keystmp
      character*4 keya(45)
      character*4 tmp
      common /ln_stack/keya
      common /ln_stack_count/ia,ik
      integer i,ia,ik
      ik=0
      ia=0
      if(keys.eq.' ') return
      keyst=' '
      keyst=keys
100   i=index(keyst,';')
      if(keyst.eq.' ') return
      if(i.eq.0) then
      ik=ik+1
      tmp=keyst(1:4)
      keya(ik)=tmp
      else
      ik=ik+1
      tmp=keyst(1:4)
      if(i.le.4) tmp(i:)=' '
      keya(ik)=tmp
      keystmp=' '
      keystmp(1:)=keyst(i+1:)
      keyst=keystmp
      goto 100
      endif
      return
      end
      
      
      character*180 function get_ln_EOF(lunit)
      logical Quit_On_Error
      common /getlnQOE/ Quit_On_Error
      character*180 get_ln_quit
      get_ln_EOF=get_ln_quit(lunit,0)
      if(Quit_On_Error) get_ln_EOF='EOF'
      End

      character*180 function get_ln_quit(lunit,icritical)
************************************************************************
* This function replaces function getln                                *
*                                                                      *
* It reads, broadcasts, and parses an input line                       *
*                                                                      *
* Blank lines or lines containing star (*) in column 1 are skipped     *
* Lines staring with exclaimation (!) in column 1 are skipped          *
*                                                                      *
* After this routine has been called, data can be retrieved using      *
* the subroutines                                                      *
*                                                                      *
*   Get_F(icol,array,n)  (for floating point values)                   *
*   Get_I(icol,iarry,n)  (for integer values)                          *
*   Get_S(icol,strgs,n)  (for character strings)                       *
*                                                                      *
* where icol is the first non-blank work to be taken, and n is the     *
* number of data.                                                      *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      Character*180 line
      Character*256 filename
* mxn should be len(line)/2+1
      parameter (mxn=91)
      common/cgetlc/ line
      common/cgetln/ ncol,istrt(mxn),iend(mxn)
      common /igetline/ igetline, myunit
      logical Quit_On_Error
      common /getlnQOE/ Quit_On_Error
      Quit_On_Error=.false.
      myunit=lunit
1     read(lunit,'(A)',err=100,end=200) line
        igetline=igetline+1
      if(line.eq.' '.or.line(1:1).eq.'*'.or.line(1:1).eq.'!') goto 1
      l=len(line)
      Do i = 1, l
         If (ichar(line(i:i)).eq.9) line(i:i)=' '
         If (line(i:i).eq.';') line(i:l)=' '
      End Do
      ncol=0
      j=1
 10   icom=0
      do i=j,l
        if(line(i:i).eq.',') then
          icom=icom+1
          if(icom.ne.1) goto 20
        else
          if(line(i:i).ne.' ') goto 20
        end if
      end do
      get_ln_quit=line
      return
 20   do j=i,l
        if(line(j:j).eq.' ' .or.  line(j:j).eq.',') goto 30
      end do
      j=l+1
 30   ncol=ncol+1
      istrt(ncol)=i
      iend(ncol)=j-1
      goto 10
*
100   filename=' '
      inquire(unit=lunit,name=filename)
      if(filename.ne.' ') then
        write(6,'(a,a)') 'Error reading file=',filename
      else
        write(6,'(a,i8)') 'Error reading unit=',lunit
      endif
        write(6,'(a)') 'Line: ',line(1:80)
      Quit_On_Error=.true.
200   if(icritical.eq.0) then
       Quit_On_Error=.true.
      return
      endif
c
      filename=' '
      inquire(unit=lunit,name=filename)
      if(filename.ne.' ') then
        write(6,'(a,a)') 'EOF reached for file=',filename
      else
        write(6,'(a,i8)') 'EOF reached for unit=',lunit
      endif

      Quit_On_Error=.true.
      end
*
      subroutine Get_F(icol,val,n)
      implicit real*8 (a-h,o-z)
      Character*180 line
      Character*80 string
* mxn should be len(line)/2+1
      parameter (mxn=91)
      common/cgetlc/ line
      common/cgetln/ ncol,istrt(mxn),iend(mxn)
      common /igetline/ igetline, myunit
      dimension val(n)
      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          i1=istrt(ic)
          i2=iend(ic)
          if(i1.le.i2) then
            string=' '
            string(80+i1-i2:80)=line(i1:i2)
            read(string,'(F80.0)',err=600,end=600) val(i)
          else
            val(i)=0.0D0
          end if
          ic=ic+1
        else
          write(6,110) icol+n-1,line
110       format(/' ERROR IN GET_F: TRYING TO READ',i4,' VALUES'/1x,a)
          Call FindErrorLine
          Call WarningMessage(2,'Error in Get_F')
          Call Quit_OnUserError()
        end if
      end do
      return
*
600       Call FindErrorLine
          Call WarningMessage(2,'Error in Get_F')
          Call Quit_OnUserError()
      end

      subroutine Get_I(icol,ival,n)
      implicit real*8 (a-h,o-z)
      Character*180 line
      Character*80 string
* mxn should be len(line)/2+1
      parameter (mxn=91)
      common/cgetlc/ line
      common/cgetln/ ncol,istrt(mxn),iend(mxn)
      common /igetline/ igetline, myunit
      dimension ival(n)
      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          i1=istrt(ic)
          i2=iend(ic)
          if(i1.le.i2) then
            string=' '
            string(80+i1-i2:80)=line(i1:i2)
            read(string,'(i80)',err=600,end=600) ival(i)
          else
            ival(i)=0
          end if
          ic=ic+1
        else
          write(6,210) icol+n-1,line
210       format(/' ERROR IN GET_I: TRYING TO READ',i4,' VALUES'/1x,a)
          Call FindErrorLine
          Call WarningMessage(2,'Error in Get_I')
          Call Quit_OnUserError()
        end if
      end do
      return
600       Call FindErrorLine
          Call WarningMessage(2,'Error in Get_I')
          Call Quit_OnUserError()
      end
*
      Subroutine Get_S(icol,str,n)
      character*(*) str(n)
      Character*180 line
      parameter (mxn=91)
      common/cgetlc/ line
      common/cgetln/ ncol,istrt(mxn),iend(mxn)
      common /igetline/ igetline, myunit
      ic=icol
      do i=1,n
        if(ic.le.ncol) then
          if(iend(ic).ge.istrt(ic)) then
            str(i)=line(istrt(ic):iend(ic))
          else
            str(i)=' '
          end if
          ic=ic+1
        else
          write(6,210) icol+n-1,line
210       format(/' ERROR IN GET_S: TRYING TO READ',i4,' STRINGS'/1x,a)
          Call FindErrorLine
          Call WarningMessage(2,'Error in Get_S')
          Call Quit_OnUserError()
        end if
      end do
      end
      subroutine Read_v(lunit,work,istrt,iend,inc,ierr)
      implicit real*8 (a-h,o-z)
      dimension work(iend)
      ierr=0
      read(lunit,*,err=100) (work(istrt+i),i=0,iend-istrt,inc)
      goto 110
100   ierr=1
110   continue
      return
      end
      subroutine Read_iv(lunit,iwork,istrt,iend,inc,ierr)
      implicit real*8 (a-h,o-z)
      dimension iwork(iend)
      ierr=0
      read(lunit,*,err=100) (iwork(istrt+i),i=0,iend-istrt,inc)
      goto 110
100   ierr=1
110   continue
      return
      end
       subroutine FindErrorLine
       character *180 line
       common /igetline/igetline,myunit
       lunit=myunit
       isave=igetline
       rewind (lunit)
 2     read(lunit,'(a)', end=300) Line
       Call LeftAd(Line)
       Call UpCase(Line)
       if(Line(1:1).eq.'&') then
         line=line(2:)
          goto 3
       endif
       goto 2
 3      igetline=0
        write(6,'(a,a,a)') ' >>>>> Input file for module ',
     *   line(1:index(line,' ')),' <<<<<'
 1     read(lunit,'(A)',err=100,end=200) line
        igetline=igetline+1
        if(igetline.eq.isave) then
         write (6,*) '******   Error  *******'
         write (6,'(a)') line
         write (6,'(a)')
         Call WarningMessage(2,'Error in FindErrorLine')
         call Quit_OnUserError()
        endif
        if(isave-igetline.le.50) then
         write (6,'(a)') line
        endif
        goto 1
c        write(6,'(a)') ' >>>>> Input error <<<<<'
c        rewind(lunit)
c        igetline=0
c        goto 1
100    continue
200    continue
300    continue
       Call WarningMessage(1,'FindErrorLine:'//
     &  ' Error in input was not located;'//
     &  '  Please, check it manually!')
       return
       end

       subroutine ResetErrorLine
       common /igetline/igetline,myunit
       igetline=0
       return
       end
