program Dif_Finitas
implicit none

integer :: i, j, n, c, k, l, m, imax
double precision :: delta, es, lambda, kappa, posx, posy
double precision, allocatable :: qx(:,:), qy(:,:), T(:), A(:,:), b(:), Temperaturas(:,:), Calores(:,:)

write(*,*) 'Ingrese el tamaño del paso'
read(*,*) delta

n = 50/delta
posx = 0.0d0
posy = 50.0d0

kappa = 0.49!*4.184

allocate(qx(1:n,1:n), qy(1:n,1:n), T(1:n**2), A(1:n**2,1:n**2), b(1:n**2), Temperaturas(1:n+1,1:n+1), Calores(1:n,1:n))

If (n == 4) then
Open(Unit = 10, file = 'input1.txt')
read(10,*), (b(k), k = 1, n**2)
 Close(Unit = 10)

Elseif (n == 8) then
Open(Unit = 11, file = 'input2.txt')
read(11,*), (b(k), k = 1, n**2)
 Close(Unit = 11)
End If
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

If (n == 4) Then
Open(Unit = 13, file = 'Mat1.txt')
Do m = 1, n**2
	write(13,'(99(1x,f4.0))') (A(m,j), j = 1,n**2)
End Do

Elseif(n == 8) Then
Open(Unit = 13, file = 'Mat2.txt')
Do m = 1, n**2
	write(13,'(99(1x,f4.0))') (A(m,j), j = 1,n**2)
End Do
 Close(Unit = 13)
End If

Call Gseid(A, b, n, T, imax, es, lambda)

Do i = 2,n+1
	Temperaturas(i,:) = T((i-2)*n+1:(i-1)*n)
End Do
Temperaturas(1,:) = 100
Temperaturas(:,n+1) = 50
Temperaturas(1,n+1) = 150

qx(:,2:n) = -kappa*(Temperaturas(2:n+1,3:n+1) - Temperaturas(2:n+1,1:n))/(2*delta)
qy(1:n,:) = -kappa*(Temperaturas(1:n,1:n) - Temperaturas(3:n+1,1:n))/(2*delta)
qx(:,1) = 0
qy(n,:) = 0

If (n == 4) then
Open(Unit = 12, file = 'ouput1.txt')
Do l = 1,n+1
	write(12,"(15F10.5)") Temperaturas(l,:)
End Do
 Close(Unit = 12)

Elseif(n == 8) then
Open(Unit = 12, file = 'ouput2.txt')
Do l = 1,n+1
	write(12,"(15F10.5)") Temperaturas(l,:)
End Do
 Close(Unit = 12)
End If

01	format(1x, 1p16e15.7)

!Call Gseid(A, b, n, T, imax, es, lambda)

!Do i = 2,n+1
	!Temperaturas(i,:) = T((i-2)*n+1:(i-1)*n)
!End Do
!Temperaturas(1,:) = 100
!Temperaturas(:,n+1) = 50
!Temperaturas(1,n+1) = 150


!qx(:,2:n) = -kappa*(Temperaturas(2:n+1,3:n+1) - Temperaturas(2:n+1,1:n))/(2*delta)
!qy(1:n,:) = -kappa*(Temperaturas(1:n,1:n) - Temperaturas(3:n+1,1:n))/(2*delta)
!qx(:,1) = 0
!qy(n,:) = 0

!Do i = 1, n+1
		!write(*,01) Temperaturas(i,:)
!End Do

!write(*,*) Temperaturas(5,1)

Open(Unit = 14, file = 'Calor.txt')
write(14,*) 'Componentes en x del flujo de calor'
Do k = 1,n
	write(14,01) qx(k,:)
End Do
write(14,*) 'Componentes en y del flujo de calor'
Do l = 1,n
	write(14,01) qy(l,:)
End Do
write(14,*) 'Magnitud del flujo de calor'
Do m = 1,n
	write(14,01) Sqrt(qx(m,:)**2 + qy(m,:)**2)
End Do
 Close(Unit = 14)

deallocate(qx, qy, T, A, b)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Gseid(a, b, n, x, imax, es, lambda)

	implicit none

	integer :: ii, jj, n, iter, centinela, imax
	double precision :: a(n**2, n**2), b(n**2), x(n**2), dummy, suma, es, ea, lambda, old

	es = 1.0d-10
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
!	write(*,*) iter
End Subroutine


