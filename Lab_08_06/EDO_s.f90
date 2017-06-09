program EDO

implicit none

integer :: i, j, k, n, m
double precision :: y, xi, xf, dx, xout, x, xend, h
double precision, allocatable :: xp(:), yp(:)


write(*,*) "Ingrese, en el siguiente orden, los parámetros correspondientes a: 'Valor inicial de la variable dependiente', 'Valor inicial de la variable independiente', 'Valor final de la variable independiente', 'Tamaño del paso', 'Intervalo de salida'"
read(*,*) y, xi, xf, dx, xout
write(*,*) 'Gracias'

x = xi
n = (xf - xi)/dx + 1
m = 1

allocate(xp(1:n), yp(1:n))

xp(1) = x
yp(1) = y

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
	If(x >= xf) Then
		Exit
	End If
End Do

deallocate(xp,yp)

end program

!-------------------------------------------------------------------------------------------------------------------

Subroutine Integrator(x, y, h, xend)

Do
	If((xend - x) < h) Then
		h = xend - x
	End If
	Call Euler(x, y, h, ynew)
	yy = ynew
	If(x >= xend) Then
		Exit
	End If
End Do
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine Euler(x, y, h, ynew)
	Call Derivs(x, y, dydx)
	ynew = y + dydx*h
	x = x + h
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine Derivs(x, y, dydx)
dydx = -2*x^3 + 12*x^2 - 20*x + 8.5
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine Heun(x, y, h, ynew)
	Call Derivs(x, y, h, dy1dx)
	ye = y + dy1dx*h
	Call Derivs(x + h, ye, dy2dx)
	Slope = (dy1dx + dy2dx)/2.0d0
	ynew = y + Slope*h
	x = x + h
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine Midpoint(x, y, h, ynew)
	Call Derivs(x, y, dydx)
	ym = y + dydx*h/2
	Call Derivs(x + h/2, ym, dymdx)
	ynew = y + dymdx*h
	x = x + h
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine HeunIter(x, y, h, ynew)
	es = 0.01
	maxit = 20
	Call Derivs(x, y, dy1dx)
	ye = y + dy1dx*h
	iter = 0
	Do
		yeold = ye
		Call Derivs(x + h, ye, dy2dx)
		slope = (dy1dx + dy2dx)/2
		ye = y + slope*h
		ea = abs(ye-yeold)*100/ye
		If(ea <= es .or. iter > maxit) Then
			Exit
		End If
	End Do
ynew = ye
x = x + h
End Subroutine

!-------------------------------------------------------------------------------------------------------------------

Subroutine RK4(x, y, h, ynew)
	Call Derivs(x, y, k1)
	ym = y + k1*h/2
	Call Derivs(x + h/2, ym, k2)
	ym = y + k2*h/2
	Call Derivs(x + h/2, ym, k3)
	ye = y + k3*h
	Call Derivs(x + h, ye, k4)
	slope = (k1 + 2*(k2 + k3) + k4)/6
	ynew = y + slope*h
	x = x + h
End Subroutine
