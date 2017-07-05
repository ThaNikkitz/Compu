program Dif_Finitas
implicit none

integer :: i, j, n, c, k, l, m, imax
double precision :: delta, es, lambda
double precision, allocatable :: qx(:), qy(:), T(:), A(:,:), b(:) 

write(*,*) 'Ingrese el tamaño del paso'
read(*,*) delta

n = 50/delta

allocate(qx(1:n), qy(1:n), T(1:n**2), A(1:n**2,1:n**2), b(1:n**2))

Open(Unit = 10, file = 'input.txt')
read(10,*), (b(k), k = 1, n**2)
 Close(Unit = 10)

 c = 1

Do j = 1,n
	Do i = 1,n
		If(i == 1 .and. j == 1) then		!nodo superior izquierdo
			A(c,c) = -4
			A(c,c+1) = 2
			A(c,c+n) = 1
		Elseif(i == n .and. j == 1) then	!nodo superior derecho
			A(c,c) = -4
			A(c,c-1) = 1
			A(c,c+n) = 1
		Elseif(i == n .and. j == n) then	!nodo inferior derecho
			A(c,c) = -4
			A(c,c-1) = 1
			A(c,c-n) = 2
		Elseif(i == 1 .and. j == n) then	!nodo inferior izquierdo
			A(c,c) = -4
			A(c,c-n) = 2
			A(c,c+1) = 2
		Elseif(j == 1) then			!nodos borde superior
			A(c,c) = -4
			A(c,c+n) = 1
			A(c,c+1) = 1
			A(c,c-1) = 1		
		Elseif(i == 1) then			!nodos borde izquierdo
			A(c,c) = -4
			A(c,c+n) = 1
			A(c,c+1) = 2
			A(c,c-n) = 1
		Elseif(i == n) then			!nodos borde derecho
			A(c,c) = -4
			A(c,c-1) = 1
			A(c,c+n) = 1
			A(c,c-n) = 1
		Elseif(j == n) then			!nodos borde inferior
			A(c,c) = -4
			A(c,c+1) = 1
			A(c,c-1) = 1
			A(c,c-n) = 2
		Else					!nodos centrales
			A(c,c) = -4
			A(c,c+1) = 1
			A(c,c-1) = 1
			A(c,c+n) = 1
			A(c,c-n) = 1
		End If
		c = c + 1
	End Do
End Do

Open(Unit = 11, file = 'ouput.txt')

Do l = 1,n**2
	write(11,"(16F5.2)") A(l,:)
End Do

 Close(Unit = 11)

01	format(1x, 1p3e15.7)

Call Gseid(A, b, n, T, imax, es, lambda)

!Open(Unit = 11, file = 'ouput.txt')
!Do k = 1,n
!	Do l = 1,n
!		write(11,01) T(k*l)
!	End Do
!End Do
 !Close(Unit = 11)

write(*,01) T(57:64)

deallocate(qx, qy, T, A, b)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Gseid(a, b, n, x, imax, es, lambda)

	implicit none

	integer :: ii, jj, n, iter, centinela, imax
	double precision :: a(n**2, n**2), b(n**2), x(n**2), dummy, suma, es, ea, lambda, old

	es = 0.0001
	lambda = 0.01
	imax = 1000000

	Do ii = 1, n**2
		dummy = a(ii,ii)
		Do jj = 1, n**2
			a(ii,jj) = a(ii,jj)/dummy
		End Do
		b(ii) = b(ii)/dummy
	End Do

	Do ii = 1, n**2
		suma = b(ii)
		Do jj = 1, n**2
			If (ii /= jj) Then
				suma = suma - a(ii,jj)*x(jj)
			End If
		End Do
		x(ii) = suma
	End Do

	iter = 1
	centinela = 0
	Do While(centinela /= 1 .and. iter < imax)
		centinela = 1
		Do ii = 1, n**2
			old = x(ii)
			suma = b(ii)
			Do jj = 1,n**2
				If (ii /= jj) Then
					suma = suma - a(ii,jj)*x(jj)
				End If
			End Do
			x(ii) = lambda*suma + (1.0d0-lambda)*old
			If ((centinela == 1) .and. (x(ii) /= 0.0d0)) Then 
				ea = Abs((x(ii) - old)/x(ii))*100
				If (ea > es) Then 
					centinela = 0
				End If
			End If
		End Do
		iter = iter + 1
	End Do
	write(*,*) iter
End Subroutine


