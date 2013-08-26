
        set term post eps color enhanced 12
        set out "likelihood_best_section_for_mass6.eps"
        set xlabel "{/=28 mass:6}"
        set title "{/=28 likelihood: best section for mass:6}" 
        a=173.07
        b=0.0
        c=1696.35678673
        f(x)=1/(2*b**2)*(x-a)**2+c
        fit f(x) 'likelihood_best_section_for_mass6.dat' u 1:2:3 via a,b,c
        plot 'likelihood_best_section_for_mass6.dat' u 1:2:3 with errorbar title "" , 1/(2*b**2)*(x-a)**2+c title "fit" 
