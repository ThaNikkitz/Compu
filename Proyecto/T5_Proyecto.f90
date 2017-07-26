program Dif_Finitas
implicit none

integer :: i, j, n, c, l, m, imax, k
double precision :: delta, es, lambda, kappa, ro, Cp, alpha, q
double precision, allocatable :: qx(:,:), qy(:,:), T(:), A(:,:), b(:), Temperaturas(:,:), Calores(:,:)

write(*,*) 'Ingrese el número de nodos en "x" y en "y" como nx,ny. Debe ser par'
read(*,*) n, m

delta = 2.0d0/n
kappa = 100.0d0
ro = 8862.0d0
 Cp = 421.0d0
alpha = kappa/(ro*Cp)
q = 1.0d5

allocate(qx(1:n,1:m), qy(1:n,1:m), T(1:n*m), A(1:n*m,1:n*m), b(1:n*m), Temperaturas(1:m,1:n+2), Calores(1:n,1:m))

Do i = 1,m
	Do j = 1,n
		If (i == 1 .and. j <= n/2) Then
			b((j-1)*m+i) = -2.0d0*delta*q/kappa
		Elseif (i == m .and. j <= n/2) Then
			b((j-1)*m+i) = -2.0d0*delta*q/kappa
		Else
			b((j-1)*m+i) = 0.0d0
		End If
	End Do
End Do

 c = 1
! Se resuelve para un vector (T_1,1 ; T_1,2 ; T_1,3 ... T_n,m-1 ; T_n,m)
Do j = 1,n
	Do i = 1,m
		If(i == 1 .and. j == 1) then		!nodo superior izquierdo
			A(c,c) = -4
			A(c,c+1) = 2
			A(c,c+m) = 1
		Elseif(i == 1 .and. j == n) then	!nodo superior derecho
			A(c,c) = -4
			A(c,c+1) = 2
			A(c,c-m) = 1
		Elseif(i == m .and. j == n) then	!nodo inferior derecho
			A(c,c) = -4
			A(c,c-1) = 2
			A(c,c-m) = 1
		Elseif(i == m .and. j == 1) then	!nodo inferior izquierdo
			A(c,c) = -4
			A(c,c+m) = 1
			A(c,c-1) = 2
		Elseif(i == 1) then			!nodos borde superior
			A(c,c) = -4
			A(c,c+m) = 1
			A(c,c-m) = 1
			A(c,c+1) = 2		
		Elseif(j == 1) then			!nodos borde izquierdo
			A(c,c) = -4
			A(c,c+1) = 1
			A(c,c+m) = 1
			A(c,c-1) = 1
		Elseif(j == n) then			!nodos borde derecho
			A(c,c) = -4
			A(c,c-m) = 1
			A(c,c+1) = 1
			A(c,c-1) = 1
		Elseif(i == m) then			!nodos borde inferior
			A(c,c) = -4
			A(c,c+m) = 1
			A(c,c-1) = 2
			A(c,c-m) = 1
		Else					!nodos centrales
			A(c,c) = -4
			A(c,c+1) = 1
			A(c,c-1) = 1
			A(c,c+m) = 1
			A(c,c-m) = 1
		End If
		c = c + 1
	End Do
End Do

Open(Unit = 13, file = 'A.txt')
Do i = 1, n*m
	write(13,'(99(1x,f4.0))') (A(i,j), j = 1,n*m)
End Do
write(13,*)
write(13,*) b
 Close(Unit = 13)

Call Gseid(A, b, n, m, T, imax, es, lambda)

Do j = 2,n+1
	Temperaturas(:,j) = T((j-2)*m+1:(j-1)*m)
End Do

!write(*,*) Temperaturas(1,n+2)

!Temperaturas(:,1) = 0.0d0
!Temperaturas(:,n+2) = 0.0d0

!write(*,*) 'Todo bien por aca'

!qx(:,2:n) = -kappa*(Temperaturas(2:n+1,3:n+1) - Temperaturas(2:n+1,1:n))/(2*delta)
!qy(1:n,:) = -kappa*(Temperaturas(1:n,1:n) - Temperaturas(3:n+1,1:n))/(2*delta)
!qy(:,1) = 0
!qy(n,:) = 0

01	format(1x, 1p16e15.7)

Open(Unit = 12, file = 'output.txt')
Do l = 1,m
	write(12,01) Temperaturas(l,:)
End Do
 Close(Unit = 12)

!Acá borré mucho, ojo :P

deallocate(qx, qy, T, A, b, Calores, Temperaturas)

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Gseid(a, b, n, m, x, imax, es, lambda)

	implicit none

	integer :: ii, jj, n, m, iter, centinela, imax
	double precision :: a(n*m, n*m), b(n*m), x(n*m), dummy, suma, es, ea, lambda, old

	es = 1.0d-10
	lambda = 0.01
	imax = 1000000

	Do ii = 1,n*m
		dummy = a(ii,ii)
		Do jj = 1, n*m
			a(ii,jj) = a(ii,jj)/dummy
		End Do
		b(ii) = b(ii)/dummy
	End Do

	Do ii = 1, n*m
		suma = b(ii)
		Do jj = 1, n*m
			If (ii /= jj) Then
				suma = suma - a(ii,jj)*x(jj)
			End If
		End Do
		x(ii) = suma
	End Do

write(*,*) 'checkpoint!'

	iter = 1
	centinela = 0
	Do While(centinela /= 1 .and. iter < imax)
		centinela = 1
		Do ii = 1, n*m
			old = x(ii)
			suma = b(ii)
			Do jj = 1,n*m
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
		If (Mod(iter,1500) == 0) Then
			write(*,*) 'Iteracion número: ', iter
		End If
		If (Mod(iter,300000) == 1) Then
			write(*,*) 'Anda a tomarte un tecito, le falta harto'
			write(*,*) 'Para que te hagas una idea, el error va en: ', ea
		End If
		If (Mod(iter,400000) == 1) Then
			write(*,*) 'Y tu pc ya empieza a oler a quemado...'
		End If
	End Do
	write(*,*) iter
End Subroutine


