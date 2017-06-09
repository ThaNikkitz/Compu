program EDO

implicit none

integer :: i, j, k, n, m
double precision :: y, xi, xf, dx, xout, x, xend, h
double precision, allocatable :: xp(:), yp(:)


write(*,*) 'Ingrese, en el siguiente orden, los parámetros correspondientes a: "Valor inicial de la variable dependiente", "Valor inicial de la variable independiente", "Valor final de la variable independiente", "Tamaño del paso", "Intervalo de salida"'
read(*,*) y, xi, xf, dx, xout
write(*,*) 'Gracias'

x = xi
n = (xf - xi)/dx + 1
m = 0

allocate(xp(0:n), yp(0:n))

xp(0) = x
yp(0) = y

Do
	xend = x + xout
	If (xend >  xf) Then
		xend = xf
	End If
	h = dx	
	Call Integrator(x, y, h, xend)
	m = m + 1
	xp(m) = x
	yp(m) = y
	If(x => xf) Then Exit
	End If

deallocate(xp,yp)

end program
