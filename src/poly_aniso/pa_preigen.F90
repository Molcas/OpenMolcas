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

subroutine pa_preigen(exch,lmax,ibas,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,JITO_exchange,WLIN,WDIP,WKEX,WITO,W,Z,iPrint)
! this function prints the energies and eigenvectors of the interaction Hamiltonians:

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: exch, lmax, ibas(exch,lmax), iPrint
logical(kind=iwp), intent(in) :: Dipol, AnisoLines1, AnisoLines3, AnisoLines9, KE, JITO_exchange
real(kind=wp), intent(in) :: WLIN(exch), WDIP(exch), WKEX(exch), WITO(exch), W(exch)
complex(kind=wp), intent(in) :: Z(exch,exch)
integer(kind=iwp) :: i, ipar, iss, j, jEnd, m, nf
character(len=80) :: fmtline

write(u6,*)
write(u6,'(a)') repeat('%',100)
write(u6,'(10x,a)') 'EigenValues of the Magnetic Interaction'
write(u6,'(a)') repeat('%',100)

#ifdef _DEBUGPRINT_
write(u6,*) 'pa_preigen:   AnisoLines1 =',AnisoLines1
write(u6,*) 'pa_preigen:   AnisoLines3 =',AnisoLines3
write(u6,*) 'pa_preigen:   AnisoLines9 =',AnisoLines9
write(u6,*) 'pa_preigen:         Dipol =',Dipol
write(u6,*) 'pa_preigen:            KE =',KE
write(u6,*) 'pa_preigen: JITO_exchange =',JITO_exchange
write(u6,*) WITO
#endif

!----------------------------------------------------------------------!
if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |',' Dipole-Dipole Inter. |','   Kinetic Exchange   |', &
                   '         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |', &
                   '  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'

  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WDIP(i),'|',WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |',' Dipole-Dipole Inter. |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WDIP(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |','   Kinetic Exchange   |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|',' Dipole-Dipole Inter. |','   Kinetic Exchange   |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WDIP(i),'|',WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|',' Dipole-Dipole Inter. |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WDIP(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6a)') 'Coupled|','   Kinetic Exchange   |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WKEX(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'

! JITO_exchange active below:
!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','   Kinetic Exchange   |','     ITO Exchange     |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WKEX(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','     ITO Exchange     |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|',' Dipole-Dipole Inter. |','     ITO Exchange     |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WDIP(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|',' Dipole-Dipole Inter. |','   Kinetic Exchange   |','     ITO Exchange     |', &
                   '         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |', &
                   '  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WDIP(i),'|',WKEX(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |','     ITO Exchange     |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |',' Dipole-Dipole Inter. |','     ITO Exchange     |', &
                   '         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |', &
                   '  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WDIP(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |','   Kinetic Exchange   |','     ITO Exchange     |', &
                   '         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |', &
                   '  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',WLIN(i),'|',WKEX(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '---------------------------------------------|'

!----------------------------------------------------------------------!
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE .and. JITO_exchange) then
  write(u6,'(6a)') 'Coupled|','     Lines Model      |',' Dipole-Dipole Inter. |','   Kinetic Exchange   |', &
                   '     ITO Exchange     |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |','  absolute val., cm-1 |', &
                   '  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,6(f20.12,2x,A))') i,'|',WLIN(i),'|',WDIP(i),'|',WKEX(i),'|',WITO(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','----------------------|','----------------------|', &
                   '----------------------|','---------------------------------------------|'

!----------------------------------------------------------------------!
else

  write(u6,'(6a)') 'Coupled|','  Magnetic Exchange   |','         Total Magnetic Interaction          |'
  write(u6,'(6a)') 'State  |','  absolute val., cm-1 |','  absolute val., cm-1 |  relative val., cm-1 |'
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
  do i=1,exch
    write(u6,'(i5,2x,A,5(f20.12,2x,A))') i,'|',W(i),'|',W(i),'|',W(i)-W(1),'|'
  end do
  write(u6,'(6a)') '-------|','----------------------|','---------------------------------------------|'
end if

!cc  eigenvectors
if (iPrint >= 3) then

  write(u6,*)
  write(u6,'(a)') repeat('%',100)
  write(u6,'(10x,a)') 'EigenVectors of the Total Magnetic Interaction'
  write(u6,'(a)') repeat('%',100)

  if (exch >= 4) then
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('---------------------------|',i=1,4)
    write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,39x,A,38x,A)'
    write(u6,fmtline) ('   ',m=1,lmax),'eigenvectors of the exchange matrix','|'
  else if (exch == 3) then
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('---------------------------|',i=1,3)
    write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,25x,A,24x,A)'
    write(u6,fmtline) ('   ',m=1,lmax),'eigenvectors of the exchange matrix','|'
  else if (exch == 2) then
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('---------------------------|',i=1,2)
    write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,11x,A,10x,A)'
    write(u6,fmtline) ('   ',m=1,lmax),'eigenvectors of the exchange matrix','|'
  else if (exch == 1) then
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|','---------------------------|'
    write(fmtline,'(A,i2,A)') '(1x,',lmax,'A,5x,A,3x,A)'
    write(u6,fmtline) ('   ',m=1,lmax),'exchange eigenvector','|'
  end if

  do J=1,exch,4
    jEnd = min(exch,J+3)
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('---------------------------|',i=j,jEnd)

    if (lmax == 1) then

      write(u6,'(A,5A)') 'Exch|',('      exchange state       |',i=j,jEnd)

      write(u6,'(A,6A)') ' bas|',('                           |',i=j,jEnd)
      write(u6,'(A,6(a,i5,a))') '    |',('         ',i,'             |',i=j,jEnd)
      write(fmtline,'(A,i2,A)') '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
      write(u6,fmtline) ' ',(i,i=1,lmax),'|',('  E =',w(i)-w(1),' |',i=j,jEnd)

    else if (lmax == 2) then

      write(u6,'(A,6A)') 'Exch.  |',('      exchange state       |',i=j,jEnd)

      write(u6,'(A,6A)') 'basis  |',('                           |',i=j,jEnd)
      write(u6,'(A,6(a,i5,a))') 'on site|',('         ',i,'             |',i=j,jEnd)
      write(fmtline,'(A,i2,A)') '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
      write(u6,fmtline) ' ',(i,i=1,lmax),'|',('  E =',w(i)-w(1),' |',i=j,jEnd)

    else if (lmax == 3) then

      write(fmtline,'(A,i2,A,i2,A)') '(A,A,4A)'
      write(u6,fmtline) ' Exchange |',('      exchange state       |',i=j,jEnd)

      write(u6,'(A,6A)') ' basis on |',('                           |',i=j,jEnd)
      write(u6,'(A,6(a,i5,a))') '    sites |',('         ',i,'             |',i=j,jEnd)
      write(fmtline,'(A,i2,A)') '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
      write(u6,fmtline) ' ',(i,i=1,lmax),'|',('  E =',w(i)-w(1),' |',i=j,jEnd)

    else

      ipar = mod(lmax,2)
      nf = (3*lmax-8+ipar)/2
      write(fmtline,'(A,i2,A,i2,A)') '(',nf,'x,A,',nf,'x,A,4A)'
      write(u6,fmtline) ' Exchange','|',('      exchange state       |',i=j,jEnd)
      write(u6,fmtline) ' basis on','|',('                           |',i=j,jEnd)
      write(fmtline,'(A,i2,A,i2,A)') '(',nf,'x,A,',nf,'x,A,4(a,i5,a))'
      write(u6,fmtline) '    sites','|',('         ',i,'             |',i=j,jEnd)
      write(fmtline,'(A,i2,A)') '(A,',lmax,'i3,a,6(a,f19.12,2x,a))'
      write(u6,fmtline) ' ',(i,i=1,lmax),'|',('  E =',w(i)-w(1),' |',i=j,jEnd)

    end if

    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('----- Real ------- Imag ---|',i=j,jEnd)

    write(fmtline,'(A,i2,A)') '(A,',lmax,'i3,A,5(2F13.9,1x,a))'
    do iss=1,exch
      write(u6,fmtline) '<',(ibas(iss,m)+1,m=1,lmax),'|',(Z(iss,i),'|',i=j,jEnd)
    end do  !iss
    write(u6,'(A,16A)') '-',('---',m=1,lmax),'|',('---------------------------|',i=j,jEnd)
    write(u6,*)
  end do ! j

else
  write(u6,'(10x,a)') 'Printing of EigenVectors of the Total Magnetic Interaction is Suppressed.'
  write(u6,'(10x,a)') 'Set PRLV >=3 to activate printing of this information.'
end if

return

end subroutine pa_preigen
