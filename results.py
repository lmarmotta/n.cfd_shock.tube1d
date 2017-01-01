#!/usr/bin/python

# -----------------------------------------------------------------------------
#  Post processing script CC299 Prj 1 and 2.
#  Author: Leonardo Motta
#  Date: 10/10/2016
# -----------------------------------------------------------------------------

import os
import fileinput
import sys
import shutil
import glob
import time

os.chdir(os.getcwd())

# Builds my config file.
def write_config(p1,rho1,scheme,dissip):

    f = open('input.in','w')

    f.write("&PAR_physical\n")
    f.write("    p1 = " + str(p1) + "d0\n")
    f.write("    p4 = 1.0d0\n")
    f.write("    rho1 = " + str(rho1) + "d0\n")
    f.write("    rho4 = 1.0d0\n")
    f.write("    fgamma = 1.4d0\n")
    f.write("    R_const = 287.0d0\n")
    f.write("    F_Cp = 1004.5\n")
    f.write("    F_Cv = 717.5\n")
    f.write("/\n")
    f.write("")
    f.write("&PAR_geometry\n")
    f.write("")
    f.write("    total_mesh_points = 1001\n")
    f.write("    start_mesh_point = -5.0d0\n")
    f.write("    final_mesh_point = 5.0d0\n")
    f.write("    print_step = 100\n")
    f.write("/\n")
    f.write("\n")
    f.write("&PAR_numeric\n")
    f.write("    scheme = " + str(scheme) + "\n")
    f.write("    iterations = 10000\n")
    f.write("    time_step = 0.0001d0\n")
    f.write("/\n")
    f.write("")
    f.write("&PAR_dissip\n")
    f.write("    dissp_scheme = " + str(dissip) + "\n")
    f.write("    dissip_omega = 0.5d0\n")
    f.write("/\n")
    f.write("\n")

    f.close()


# Gets the discretization scheme name.
def get_sheme(number):
    if number == 1:
        return "Simple Centered Scheme"
    elif number == 2:
        return "Lax-Wendroff Scheme"
    elif number == 3:
        return "MacCormack Scheme"
    elif number == 4:
        return "Beam-Warming Scheme"
    elif number == 5:
        return "Steger-Warming 1st Order Scheme"
    elif number == 6:
        return "Steger-Warming 2nd Order Scheme"
    elif number == 7:
        return "Van Leer 1st Order Scheme"
    elif number == 8:
        return "Van Leer 2nd Order Scheme"
    elif number == 9:
        return "Roe Scheme"
    elif number == 10:
        return "AUSM+ Scheme"
    elif number == 11:
        return "Implicit Steger-Warming 1st Order"


# Gets the dissipation scheme name.
def get_dissip(number):
    if number == 1:
        return "Pullian Non-Linear dissipation"
    elif number == 2:
        return "Second difference dissipation"
    elif number == 3:
        return "Fourth difference dissipation"
    elif number == 0:
        return "No dissipation"


# Generate gnuplot files for epslatex outputs.
def gnuplot_eps(variable,pressure_ratio,scheme,dissip_scheme):

    # Building name of the output file.
    outFileName = variable.lower() + '_ratio' + str(int(pressure_ratio)) + '_' + scheme.replace(" ","") + '_dissp' + dissip_scheme.replace(" ","") + '.tex'
    gnuFileName = variable.lower() + '_ratio' + str(int(pressure_ratio)) + '_' + scheme.replace(" ","") + '_dissp' + dissip_scheme.replace(" ","") + '.gnu'

    # Build title of the plot
    title = 'Analytical vs Numerical | ' + scheme

    f = open(gnuFileName,'w')

    f.write("set term cairolatex monochrome size 15.0cm, 8cm\n")
    f.write('set output "' + outFileName + '"\n')
    f.write("set grid\n")
    f.write('set xtics font "Times-Roman, 10"\n')
    f.write('set ytics font "Times-Roman, 10"\n')
    f.write('set xlabel "x position" center\n')
    f.write('set ylabel "' + variable + '" center\n')
    f.write('set title "Analytical vs Numerical | ' + variable + ' | ' + scheme + ' | ' +  dissip_scheme + '" \n')
    f.write('set pointsize 0.5\n')
    f.write('set key font ",10"\n')
    f.write('set border linewidth 3\n')
    f.write('set offset graph 0.10,0.10,0.10,0.10\n')
    

    fortran_out_analytical = 'a' + variable.lower() + '.out'
    fortran_out_numerical  = variable.lower() + 'Output.out'

    f.write('plot "' + fortran_out_analytical +'" u 1:2 with linespoints lt -1 lw 1 pt 4 title "Analytical",\\\n')
    f.write( '"' + fortran_out_numerical + '" u 1:2 with lines lw 5 title "Numerical"\n')
    f.close()


# Generate Gnuplots for simple pdf outputs.
def gnuplot_pdf(variable,pressure_ratio,scheme,dissip_scheme):

    # Building name of the output file.
    outFileName = variable.lower() + '_ratio' + str(int(pressure_ratio)) + '_' + scheme.replace(" ","") + '_dissp' + dissip_scheme.replace(" ","") + '.tex'
    gnuFileName = variable.lower() + '_ratio' + str(int(pressure_ratio)) + '_' + scheme.replace(" ","") + '_dissp' + dissip_scheme.replace(" ","") + '.gnu'

    # Build title of the plot
    title = 'Analytical vs Numerical | ' + scheme

    f = open(gnuFileName,'w')

    f.write("set term cairolatex monochrome size 15.0cm, 8cm\n")
    f.write('set output "' + outFileName + '"\n')
    f.write("set grid\n")
    f.write('set xtics font "Times-Roman, 10"\n')
    f.write('set ytics font "Times-Roman, 10"\n')
    f.write('set xlabel "x position" center\n')
    f.write('set ylabel "' + variable + '" center\n')
    f.write('set title "Analytical vs Numerical | ' + variable + ' | ' + scheme + ' | ' +  dissip_scheme + '" \n')
    f.write('set pointsize 0.5\n')
    f.write('set key font ",10"\n')
    f.write('set border linewidth 3\n')
    f.write('set offset graph 0.10,0.10,0.10,0.10\n')
    

    fortran_out_analytical = 'a' + variable.lower() + '.out'
    fortran_out_numerical  = variable.lower() + 'Output.out'

    f.write('plot "' + fortran_out_analytical +'" u 1:2 with linespoints lt -1 lw 1 pt 4 title "Analytical",\\\n')
    f.write( '"' + fortran_out_numerical + '" u 1:2 with lines lw 5 title "Numerical"\n')
    f.close()

# Generate latex code.
def generate_latex(text_image_file,caption,scheme,flag):

    if scheme == 1:
        wrk = "Centrado simples"
    elif scheme == 2:
        wrk = "Lax-Wendroff"
    elif scheme == 3:
        wrk = "MacCormack"
    elif scheme == 4:
        wrk = "Beam-Warminga Impl\\'icito"
    elif scheme == 5:
        wrk = "Steger-Warming 1st order scheme"
    elif scheme == 6:
        wrk = "Steger-Warming 2st order scheme"
    elif scheme == 7:
        wrk = "Van Leer 1st order scheme"
    elif scheme == 8:
        wrk = "Van Leer 2st order scheme"
    elif scheme == 9:
        wrk = "Roe 1st order scheme"
    elif scheme == 10:
        wrk = "AUSM+"
    elif scheme == 11:
        wrk = "Implicit Steger-Warming 1st Order"
    elif scheme == 12:
        wrk = "Implicit Van-Leer 1st Order"
    elif scheme == 13:
        wrk = "Harten TVD 1st Order - Roe avg"
    elif scheme == 14:
        wrk = "Harten TVD 2nd Order - Roe avg"
    elif scheme == 15:
        wrk = "Harten TVD 1st Order - Artm avg"
    elif scheme == 16:
        wrk = "Harten TVD 2nd Order - Artm avg"
        

    if flag == 1:

        latex.write("\\subsection{" + wrk + "}\n")

        latex.write("\\begin{figure}[H]\n")
        latex.write("    \centering\n")
        latex.write("    \input{" + text_image_file + "}\n")
        latex.write("    \caption{"+ caption + "}\n")
        latex.write("    \label{fig:digraph}\n")
        latex.write("\\end{figure}\n")
        latex.write("\n\n")

    elif flag == 0:

        latex.write("\\begin{figure}[H]\n")
        latex.write("    \centering\n")
        latex.write("    \input{" + text_image_file + "}\n")
        latex.write("    \caption{"+ caption + "}\n")
        latex.write("    \label{fig:digraph}\n")
        latex.write("\\end{figure}\n")
        latex.write("\n\n")



# ----------------------------------------------------------------------- #
#                            MAIN LOOP.                                   #
# ----------------------------------------------------------------------- #

# ===================== INPUTS OF THE SCRIPT ============================ #
# Remember here to look at the main loop order since the sequence I use was the
# sequence of results I wanted in my report. Look at the latex generator functi
# on limitations to think about the organization.

# Pressure ratios desired.
pressure_ratios = [5.0]

# Schemes desired
schemes = [9]

# Artificial dissipation desired.
dissips = [0]


# Automated zone.... No user inputs from now on...
# Open Latex export file.
latex = open("bizu.txt",'w')

i = 0

# Loop through pressure ratios.
for i in xrange(len(pressure_ratios)):

    print "---------------------------------------------"
    print " + Configuring File for pressure ratio: " + str(pressure_ratios[i])
    print "--------------------------ii-----------------\n"

    # Loopt through artificial dissipations.
    for kkk in xrange(len(dissips)):

        # Loop through numerical schemes.
        for jj in xrange(len(schemes)):

            # In case of an explosion lets have some files of gnuplot.
            os.system("touch pressureOutput.out")
            os.system("touch densityOutput.out")

            # It's usefull to warn the user abouth the current state.
            print "------------------------------------------------------------"
            print " + Configuring file for scheme: " + get_sheme(schemes[jj]) + "\n"
            print " + Configuring file for dissip: " + get_dissip(dissips[kkk])
            print "------------------------------------------------------------\n"

            # Time to read the state.
            time.sleep(2)

            # Select the inputs for config file (fortran code).
            p_1 = pressure_ratios[i]
            r_1 = pressure_ratios[i]
            sch = schemes[jj]
            dis = dissips[kkk]

            # Write the config file.
            write_config(p_1,r_1,sch,dis) 


            ### Calling program ### 
            os.system("./sod")


            # Configuring gnuplot scripts for pressure and density.
            gnuplot_pdf('Density',pressure_ratios[i],get_sheme(schemes[jj]),get_dissip(dissips[kkk]))
            gnuplot_pdf('Pressure',pressure_ratios[i],get_sheme(schemes[jj]),get_dissip(dissips[kkk]))

            # Generate latex file for density
            outFileNamed = 'density_ratio' + str(int(pressure_ratios[i])) + '_' + get_sheme(schemes[jj]) + '_dissp' + get_dissip(dissips[kkk]) + '.tex'
            captiond = "Gr\\'afico de densidade | " + get_sheme(schemes[jj]) + " at $p_{4}/p_{1} = " + str(pressure_ratios[i]) + "$" + " usando " + get_dissip(dissips[kkk])

            generate_latex(outFileNamed.replace(" ",""),captiond,schemes[jj],1)

            # Generate latex file for pressure
            outFileNamep = 'pressure_ratio' + str(int(pressure_ratios[i])) + '_' + get_sheme(schemes[jj]) + '_dissp' + get_dissip(dissips[kkk]) + '.tex'
            captionp = "Gr\\'fico de press\~oes " + get_sheme(schemes[jj]) + " at $p_{4}/p_{1} = " + str(pressure_ratios[i]) + "$" + " usando " + get_dissip(dissips[kkk])

            generate_latex(outFileNamep.replace(" ",""),captionp,schemes[jj],0)

            # plotting
            to_be_ploted = 'density_ratio' + str(int(pressure_ratios[i])) + '_' + get_sheme(schemes[jj]).replace(" ","") + '_dissp' + get_dissip(dissips[kkk]).replace(" ","") + '.gnu'
            to_be_plotep = 'pressure_ratio' + str(int(pressure_ratios[i])) + '_' + get_sheme(schemes[jj]).replace(" ","") + '_dissp' + get_dissip(dissips[kkk]).replace(" ","") + '.gnu'

            os.system("gnuplot " + to_be_ploted)
            os.system("gnuplot " + to_be_plotep)

            # Cleaning useless stuff.
            os.system("rm *.gnu")
            os.system("rm *Output.out")


latex.close()
