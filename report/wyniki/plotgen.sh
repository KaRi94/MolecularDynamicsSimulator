#!/bin/bash

# {/Symbol \260} stopien

#$XLABEL = "Temperatura [{/Symbol \260}C]"
#$YLABEL = "Napiecie [V]"
#$FILE = "fioletowy.dat"
#set key left top

gnuplot <<- EOF
	set xtics font "Verdana,20"
	set xlabel font "Verdana,25"
	set xlabel "czas [s]"
	set ytics font "Verdana,20"
	set ylabel font "Verdana,25"
	set ylabel "Energia [J]"
	set pointsize 2

	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "2particles-collision.eps"
	plot 'particles-collision.dat' using 1:5 w l title 'Energia kinetyczna', 'particles-collision.dat' using 1:6 w l title 'Energia potencjalna', 'particles-collision.dat' using 1:7 w l title 'Energia calkowita'



	set xlabel "czas [s]"
	set ylabel "x [m]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "pos-particles-collision.eps"
	plot 'particles-collision.dat' using 1:4 title 'Polozenie czastki 1'



	set xlabel "x [m]"
	set ylabel "y [m]"
	set autoscale xy

	set output "one-particle-xy.eps"
	plot 'one-particle.dat' using 2:3 title 'Polozenie czastki xy'



	set xlabel "y [m]"
	set ylabel "z [m]"
	set autoscale xy

	set output "one-particle-yz.eps"
	plot 'one-particle.dat' using 3:4 title 'Polozenie czastki yz'



	set xlabel "czas [s]"
	set ylabel "Energia [J]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "oscilating-particles.eps"
	plot 'oscilating-particles.dat' using 1:5 w l title 'Energia kinetyczna', 'oscilating-particles.dat' using 1:6 w l title 'Energia potencjalna', 'oscilating-particles.dat' using 1:7 w l title 'Energia calkowita'


	
	set xlabel "czas [s]"
	set ylabel "Energia [J]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "border-oscilating-particles.eps"
	plot 'border-oscilating-particles.dat' using 1:5 w l title 'Energia kinetyczna', 'border-oscilating-particles.dat' using 1:6 w l title 'Energia potencjalna', 'border-oscilating-particles.dat' using 1:7 w l title 'Energia calkowita'

	set xlabel "czas [s]"
	set ylabel "Energia [J]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "border-oscilating-particles2.eps"
	plot 'border-oscilating-particles2.dat' using 1:5 w l title 'Energia kinetyczna', 'border-oscilating-particles2.dat' using 1:6 w l title 'Energia potencjalna', 'border-oscilating-particles2.dat' using 1:7 w l title 'Energia calkowita'

	set xlabel "czas [s]"
	set ylabel "x [m]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "border-oscilating-particles2-pos.eps"
	plot 'border-oscilating-particles2.dat' using 1:2 title 'czastka 1', 'border-oscilating-particles2.dat' using 1:3 title 'czastka 2'

	set xlabel "czas [s]"
	set ylabel "x [m]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "verlet-radius-out.eps"
	plot 'verlet-radius-out.dat' using 1:2 title 'czastka 1', 'verlet-radius-out.dat' using 1:3 title 'czastka 2'

	
	set xlabel "czas [s]"
	set ylabel "x [m]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "verlet-radius-in.eps"
	plot 'verlet-radius-in.dat' using 1:2 title 'czastka 1', 'verlet-radius-in.dat' using 1:3 title 'czastka 2'

	set xlabel "czas [s]"
	set ylabel "Energia [J]"
	set autoscale xy

	set term postscript eps size 8,6 enhanced color font 'Helvetica,20' linewidth 1
	set output "simulation.eps"
	plot 'simulation-param.dat' using 1:2 w l title 'Energia kinetyczna', 'simulation-param.dat' using 1:3 w l title 'Energia potencjalna', 'simulation-param.dat' using 1:4 w l title 'Energia calkowita'

	set xlabel "czas [s]"
	set ylabel "Temperatura [K]"
	set autoscale xy

	set output "simulation-temp.eps"
	plot 'simulation-param.dat' using 1:5 w l title 'Temperatura'

	set xlabel "czas [s]"
	set ylabel "Cisnienie [Pa]"
	set xrange [1e-11:1e-10]
	set yrange [750000:1e6]
	set format x "%11.2e"

	set output "simulation-cis.eps"
	plot 'simulation-param.dat' using 1:6 w l title 'Temperatura'

EOF
