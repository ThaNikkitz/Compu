program GSeidel

	implicit none

	integer, parameter :: n = 10
	integer :: stat, i, j, iter, centinela, imax, k
!	double precision a(np, np), y(np, np), b(n), d, a1(np, np), num, c(np,np)
	double precision :: a(n, n), b(n), x(n), dummy, suma, es, ea, lambda, old
!	Definición Matriz A

    OPEN(UNIT=11,FILE='mat_diag.txt',IOSTAT=stat,ACTION='READ')
    IF (stat /= 0) THEN
        WRITE(*,*) 'Opening Unit 10 Failed with iostat ', stat, '.'
    END IF
    
    DO i = 1,n
        READ(11,*) (a(i,j), j = 1,n)
    END DO
    READ(11,*) 
    DO i = 1,n
        READ(11,*) b(i)
    END DO
    
    CLOSE(UNIT=11,IOSTAT=stat)
    IF (stat /= 0) THEN
        WRITE(*,*) 'Closing Unit 10 Failed with iostat ', stat, '.'
    END IF

01  format(1x,1p6e15.7)

!	Presentación Matriz

!    write(*,*) 'La matriz A='
!    do i=1,n
!    write(*,01) (a(i,j), j = 1,n)
!    end do
!    write(*,*)
!    write(*,*) 'El vector b='
!    do i=1,n
!    write(*,01) b(i)
!    end do

!write(*,01) x
!write(*,*)

	call GSeid(a, b, n, x, imax, es, lambda)

write(*,*) x

end program

subroutine Gseid(a, b, n, x, imax, es, lambda)

	implicit none

	integer :: i, j, n, iter, centinela, imax
	double precision :: a(n, n), b(n), x(n), dummy, suma, es, ea, lambda, old

	es = 0.0001
	lambda = 0.001
	imax = 1000000

	Do i = 1, n
		dummy = a(i,i)
		Do j = 1, n
			a(i,j) = a(i,j)/dummy
		End Do
		b(i) = b(i)/dummy
	End Do
	Do i = 1, n
		suma = b(i)
		Do j = 1, n
			If (i /= j) Then
				suma = suma - a(i,j)*x(j)
			End If
		End Do
		x(i) = suma
	End Do
	iter = 1
	centinela = 0
	Do While(centinela /= 1 .and. iter < imax)
		old = x(i)
		centinela = 1
		Do i = 1, n
			suma = b(i)
			Do j = 1,n
				If (i /= j) Then
					suma = suma - a(i,j)*x(j)
				End If
			End Do
			x(i) = lambda*suma + (1.0d0-lambda)*old
			If ((centinela == 1) .and. (x(i) /= 0.0d0)) Then 
				ea = Abs((x(i) - old)/x(i))*100
				If (ea > es) Then 
					centinela = 0
				End If
			End If
		End Do
		iter = iter + 1
	End Do
End Subroutine
