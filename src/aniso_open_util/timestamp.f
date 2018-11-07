      Subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  ModIfied:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      Implicit None

      Character * ( 8 ) ampm
      Integer d
      Character * ( 8 ) date
      Integer h
      Integer m
      Integer mm
      Character * ( 9 ) month(12)
      Integer n
      Integer s
      Character * ( 10 ) time
      Integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      Call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      If ( h .lt. 12 ) Then
        ampm = 'AM'
      Else If ( h .eq. 12 ) Then
        If ( n .eq. 0 .and. s .eq. 0 ) Then
          ampm = 'Noon'
        Else
          ampm = 'PM'
        End If
      Else
        h = h - 12
        If ( h .lt. 12 ) Then
          ampm = 'PM'
        Else If ( h .eq. 12 ) Then
          If ( n .eq. 0 .and. s .eq. 0 ) Then
            ampm = 'Midnight'
          Else
            ampm = 'AM'
          End If
        End If
      End If

      Write ( 6,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      Return
      End

