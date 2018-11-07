      Subroutine pa_sort (a, n)
      Implicit None
      Integer n, a(n), temp, k, l
      Do k = 1, n - 1
        Do l = k+1, n
         If (a(k).gt.a(l)) Then
         temp = a(k)
         a(k) = a(l)
         a(l) = temp
         End If
        End Do
      End Do
      Return
      End

