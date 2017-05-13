!Cálculo de un sistema lineal de ecuaciones mediante el Método de Sustitución de Gauss con pivote parcial.

program GAUSS
implicit none

integer, parameter :: n = 10
integer i, j, k, er, stat
double precision a(n,n), e(n,n), b(n), x(n), x1(n), c(n), r(n), tol, suma
character(150) message

tol = 1.0d-10

!Definición de la matriz A

OPEN(UNIT = 10, FILE = 'matrix3.txt', IOSTAT = stat, ACTION = 'READ')
DO i = 1,n
	read(10,*) (a(i,j), j = 1,n)
		do j=1,n
		e(i,j) = a(i,j)
		end do
end do

read(10,*)
do i = 1,n
	read(10,*) b(i)
end do

	CLOSE(UNIT=0, IOSTAT = stat)

!Presentación de la Matriz

write(*,*) 'La matriz A='
do i = 1,n
	write(*,01) (a(i,j), j = 1,n)
end do

write(*,*) 'El vector b='
do i = 1,n
	write(*,01) b(i)
end do

01	format(1x, 1p3e15.7)
	
call Gauss_method(a, b, n, x, tol, er)

OPEN(unit=2, file = 'vector_sol.dat')

write(*,*) 'Vector solución:'
do i = 1, n
	write(2,*) x(i)
	x1(i) = x(i)
	do j = 1, n
		suma = 0.d0
		suma = suma + e(i,j)*x1(j)
	end do
!	c(i) = suma
!	write(*, 01) c(i)
end do

	CLOSE(unit=2)

!Subrutina Gauss---------------------------------

contains

SUBROUTINE Gauss_method(a, b, n, x, tol, er)
integer n, i, j, k, er
double precision a(n,n), b(n), x(n), tol, s(n)
er = 0.d0
do i= 1, n
	s(i) = abs(a(i,1))
	do j= 2, n
		if (abs(a(i,j))>s(i)) then
			s(i) = abs(a(i,j))
		end if
	end do
end do
call Eliminacion(a, b, s, n, tol, er)
if (er /= -1) then
	call Sustitucion(a, b, n, x)
end if
end SUBROUTINE

!Fin Gauss---------------------------------------

!Subrutina Eliminacion---------------------------

subroutine Eliminacion(a, b, s, n, tol, er)
integer er, k, n, j, i
double precision a(n,n), b(n), s(n), tol, f
do k = 1, n - 1
	call Pivote(a, b, s, n, k)
	if (abs(a(k,k)/s(k)) < tol) then
		er = -1
	end if
	do i = k + 1, n
		f = a(i,k)/a(k,k)
		do j = k + 1, n
			a(i,j) = a(i,j) - f*a(k,j)
		end do
	b(i) = b(i) - f*b(k)
	end do
end do
if (abs(a(k,k)/s(k)) < tol) then
		er = -1
end if
end subroutine

!Fin Eliminacion---------------------------------

!Subrutina Pivote--------------------------------

subroutine Pivote(a, b, s, n, k)
integer k, n, ii, jj, p, er
double precision big, dummy, a(n,n), b(n), s(n)
p = k
big = abs(a(k,k)/s(k))
do ii = k + 1, n
	dummy = abs(a(ii,k)/s(ii))
	if (dummy > big) then
		big = dummy
		p = ii
	end if
end do
if (p /= k) then
	do jj = k, n
		dummy = a(p,jj)
		a(p,jj) = a(k,jj)
		a(k,jj) = dummy
	end do
	dummy = b(p)
	b(p) = b(k)
	b(k) = dummy
	dummy = s(p)
	s(p) = s(k)
	s(k) = dummy
end if
end subroutine

!Fin Pivote--------------------------------------

!Subrutina Sustitucion---------------------------

subroutine Sustitucion(a, b, n, x)
integer n, k, i, j
double precision a(n,n), b(n), x(n), su
x(n) = b(n)/a(n,n)
do i = n - 1, 1, -1
	su = 0
	do j = i + 1, n
		su = su + a(i,j)*x(j)
	end do
	x(i) = (b(i) - su)/a(i,i)
end do
end subroutine

!Fin Sustitucion---------------------------------

end program
