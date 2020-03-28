import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.stats
import math

# Error function for initial first-pass fitting of thermal unfolding curves
def single_thermal_curve_errfxn(xxx_todo_changeme, fluorescence, temperatures, Cp ):

    (Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope) = xxx_todo_changeme
    if (len(fluorescence) != len(temperatures)):
        print("Mismatch in fluorescence and T array lengths")
        quit()

    sq_err = 0
    for index in range(len(temperatures)):
        T = temperatures[index] + 273.15
        R = 1.987/1000
        dG = dH*(1-T/(Tm+273.15)) - Cp*(Tm+273.15-T + T*np.log(T/(Tm+273.15)))
        Ku = np.exp(-dG/(R*T))
        Y = (Ku/(1+Ku))*(unfolded_slope*T + unfolded_intercept) + (1/(1+Ku))*(folded_slope*T + folded_intercept)
        err = Y - fluorescence[index]
        sq_err += err*err

    return sq_err

# Carry out initial first-pass fitting of thermal unfolding curves
def fit_single_thermal_curve(temperatures, fluorescence, Cp):
    
    init_dH = 150
    window_size = int(len(temperatures)/10)
    ( init_bottom_slope, init_bottom_intercept, junk1, junk2, junk3 ) = scipy.stats.linregress(temperatures[0:window_size]+273.15, fluorescence[0:window_size])
    tse=len(temperatures)-1
    ( init_top_slope, init_top_intercept, junk1, junk2, junk3 ) = scipy.stats.linregress(temperatures[tse-window_size:tse]+273.15, fluorescence[tse-window_size:tse])
    fluo_midpoint = ( np.amax(fluorescence) - np.amin(fluorescence) ) / 2
    diff = abs(fluo_midpoint*10)
    init_Tm = 0
    for index in range(len(temperatures)):
        curr_diff = abs( fluorescence[index] - fluo_midpoint)
        if (curr_diff < diff) :
            init_Tm = temperatures[index]
            diff = curr_diff
    res = scipy.optimize.minimize( single_thermal_curve_errfxn, (init_Tm, init_dH, init_top_intercept, init_bottom_intercept, init_top_slope, init_bottom_slope), args=(fluorescence, temperatures, Cp), method='Nelder-Mead' )
    return res.x

def fit_all_individual_thermal_curves( concs, temperatures, all_dat, Cp, output_file):
    num_datasets = len(concs)
    single_fit_params = np.empty([6,num_datasets])
    ft = open(output_file,'w')
    ft.write('Conc\tT_m\n')
    for index in range(num_datasets):
        fit_params = fit_single_thermal_curve(temperatures,all_dat[:,index],Cp)
#        single_fit_params[:,index] = fit_params
        single_fit_params[0,index] = fit_params[0]
        single_fit_params[1,index] = fit_params[1]
        single_fit_params[2,index] = fit_params[2]
        single_fit_params[3,index] = fit_params[3]
        single_fit_params[4,index] = fit_params[4]
        single_fit_params[5,index] = fit_params[5]
        min_Tm = fit_params[0]
        min_dH = fit_params[1]
        ft.write(('%.3f uM\t%.2f oC.\n' % (concs[index], min_Tm)))
    ft.close()
    return single_fit_params

# Build the (isothermal) fraction folded curve at a specified temperature
def calculate_fraction_unfolded(binding_temperature, concs, fitting_params, Cp):
    
    num_datasets = len(concs)
    fraction_unfolded = np.empty((concs).shape)
    for c_index in range(num_datasets):
        Tm = fitting_params[0,c_index]   # Tm
        dH = fitting_params[1,c_index]   # dH

        T = binding_temperature + 273.15
        R = 1.987/1000
        dG = dH*(1-T/(Tm+273.15)) - Cp*(Tm+273.15-T + T*np.log(T/(Tm+273.15)))
        Ku = np.exp(-dG/(R*T))
        fraction_unfolded[c_index] = (Ku/(1+Ku))
        
    return fraction_unfolded


# Use Ku and Kd (and protein+ligand concs) to get the fraction unfolded
def calculate_fitted_isothermal_point(ligand_conc, Ku, Kd, protein_conc):
    Kd = abs(Kd)
    Ku = abs(Ku)
    b = protein_conc + Kd*(1+Ku) - ligand_conc
    c = -1.0 * ligand_conc * Kd * (1+Ku)
    L_free = (-b+math.sqrt(b*b-4*c))/2
    # assume that L_free = L_tot (if we don't want to use the quadratic equation above)
    #L_free = concs[index]
    # then use L_free to get fraction unfolded
    fit_fraction_unfolded = 1 / (1 + (1/Ku)*(1+L_free/Kd))
    return fit_fraction_unfolded


            
# Fitting function for binding curves
def binding_curve_errfxn(xxx_todo_changeme2, fraction_unfolded, concs, protein_conc ):
    
    (Ku, Kd) = xxx_todo_changeme2
    if (len(fraction_unfolded) != len(concs)):
        print("Mismatch in fraction_unfolded and concs array lengths")
        quit()

    sq_err = 0
    for index in range(len(concs)):
        fit_fraction_unfolded = calculate_fitted_isothermal_point(concs[index], Ku, Kd, protein_conc)
        err = fit_fraction_unfolded - fraction_unfolded[index]
        sq_err += err*err
    return sq_err


# Fit the (isothermal) binding curves
def fit_isothermal_Kd( fraction_unfolded, concs, unique_concs, protein_conc ):

    if (len(fraction_unfolded) != len(concs)):
        print("Mismatch in fraction_unfolded and concs array lengths")
        quit()

    max_Fu = 0
    for index in range(len(fraction_unfolded)):
        if ( fraction_unfolded[index] > max_Fu ):
            max_Fu = fraction_unfolded[index]
    init_Ku = max_Fu/(1.0-max_Fu)
    init_midpoint = max_Fu / 2
    init_EC50 = 0
    diff = 1
    for index in range(len(fraction_unfolded)):
        curr_diff = abs( fraction_unfolded[index] - init_midpoint)
        if (curr_diff < diff) :
            init_EC50 = concs[index]
            diff = curr_diff
    init_Kd = max_Fu * ( init_EC50 - protein_conc / 2 )

    res = scipy.optimize.minimize( binding_curve_errfxn, (init_Ku, init_Kd), args=(fraction_unfolded, concs, protein_conc), method='Nelder-Mead' )

    Ku = abs(res.x[0])
    Kd = abs(res.x[1])

    fit_values = np.empty(unique_concs.shape)
    for c_index in range(len(unique_concs)):
        fit_values[c_index] = calculate_fitted_isothermal_point(unique_concs[c_index], Ku, Kd, protein_conc)

    return Ku, Kd, fit_values


    
# Make plots of the fits compared to the experimental thermal unfolding curves
def plot_fits_to_thermal_unfolding_data(concs, temperatures, fitting_params, all_dat, Cp, save_location ):
    
    num_datasets = len(concs)
    # figure out good dimensions for the subplot step (in case we're using it)
    subplot_horiz = math.ceil(math.sqrt(num_datasets))
    subplot_vert = math.ceil(num_datasets / subplot_horiz)

    plt.figure(5, figsize=(18,11))

    for c_index in range(num_datasets):
        Tm = fitting_params[0,c_index]   # Tm
        dH = fitting_params[1,c_index]   # dH
        ti = fitting_params[2,c_index]   # unfolded intercept
        bi = fitting_params[3,c_index]   # folded intercept
        ts = fitting_params[4,c_index]   # unfolded slope
        bs = fitting_params[5,c_index]   # folded slope

        ax = plt.subplot(subplot_vert, subplot_horiz, c_index+1)

        fluorescence = all_dat[:,c_index]
        fit_fluo = np.empty((fluorescence).shape)
        for t_index in range(len(temperatures)):
            T = temperatures[t_index] + 273.15
            R = 1.987/1000
            dG = dH*(1-T/(Tm+273.15)) - Cp*(Tm+273.15-T + T*np.log(T/(Tm+273.15)))
            Ku = np.exp(-dG/(R*T))
            fit_fluo[t_index] = (Ku/(1+Ku))*(ts*T + ti) + (1/(1+Ku))*(bs*T + bi)
    
        plt.plot(temperatures, fluorescence, 'bs', temperatures, fit_fluo, 'r')

        if ( (c_index % subplot_horiz) == 0 ):
            # print the y-axis label only on the left-most column
            plt.ylabel('Fluorescence')
        # take the numbers off all columns (not just the middle ones)
        ax.set_yticklabels([])
        if ( ( c_index + subplot_horiz >= num_datasets ) ):
            # print the x-axis label only on the bottom row
            plt.xlabel('Temperature (oC)')
        else:
            # take the numbers off all other rows
            ax.set_xticklabels([])
        plt.title('Ligand conc. '+str(concs[c_index])+' uM')
            
    # make the plot
    plt.savefig(save_location+'thermal_unfolding_curves.png', dpi=100)
    #plt.show(block=False)
    return 0


# Make plots of the fits compared to the experimental thermal unfolding curves
def plot_Tm_shifts(concs, Tms, save_location):
    # Plot the Tm shifts
    plt.figure(4, figsize=(8,6))
    plt.semilogx(concs, Tms, 'ro')
    plt.xlabel('Ligand concentration (uM)')
    plt.ylabel('Tm (oC)')
    plt.title('Tm-shifts')
    plt.savefig(save_location+'Tm_shifts.png', dpi=100)
    #plt.show(block=False)

    
# Make plots of the fits compared to the isothermal fraction folded values
def plot_isothermal_binding_fits(isothermal_Tlist, concs, all_isothermal_data, unique_concs, all_isothermal_fits, save_location ):
    
    num_plots = len(isothermal_Tlist)
    # figure out good dimensions for the subplot step (in case we're using it)
    subplot_horiz = math.ceil(math.sqrt(num_plots))
    subplot_vert = math.ceil(num_plots / subplot_horiz)

    plt.figure(3, figsize=(14,12))

    for t_index in range(len(isothermal_Tlist)):
        temperature_string = isothermal_Tlist[t_index]                
        ax = plt.subplot(subplot_vert, subplot_horiz, t_index+1)
        plt.semilogx(concs, all_isothermal_data[:,t_index], 'bs', unique_concs, all_isothermal_fits[:,t_index], 'r')

        if ( (t_index % subplot_horiz) == 0 ):
            # print the y-axis label only on the left-most column
            plt.ylabel('Fraction unfolded')
        if ( ( t_index + subplot_horiz >= num_plots ) ):
            # print the x-axis label only on the bottom row
            plt.xlabel('Ligand concentration (uM)')

        plt.title('Isothermal data at '+temperature_string+' oC')
            
    # make the plot (with many panels)
    plt.savefig(save_location+'isothermal_curves.single.png', dpi=100)
    #plt.show(block=False)

    # make another plot with all the curves on a single plot
    if ( num_plots > 1 ):
        plt.figure(2, figsize=(8,6))
        for t_index in range(len(isothermal_Tlist)):
            temperature_string = isothermal_Tlist[t_index]                
            plt.semilogx(concs, all_isothermal_data[:,t_index], 'bs', unique_concs, all_isothermal_fits[:,t_index], 'r')
            plt.xlabel('Ligand concentration (uM)')
            plt.ylabel('Fraction unfolded')
            plt.title('Isothermal data at all requested temperatures')
            plt.savefig(save_location+'isothermal_curves.all.png', dpi=100)
            #plt.show(block=False)
    
    return 0

def process_single_protein(protein):
    print('\n-------------'+protein+'-------------')
    protein_conc = np.float(input('Enter protein concentration for '+protein+' :'))
    Cp = np.float(input('Enter Cp for '+protein+':'))
    fluo_datafile = '../Output/' + protein + '/' + protein + '_fluo.txt'
    conc_datafile = '../Output/' + protein + '/' + protein + '_conc.txt'
    parameter_datafile = '../Output/' + protein + '/' + protein + '_parameters.txt'
    Tm_datafile = '../Output/' + protein + '/' + protein + '_Tm_values.txt'
    KuKd_datafile = '../Output/' + protein + '/' + protein + '_KuKd_values.txt'
    isothermal_output_tag= '../Output/' + protein + '/' + protein + '_isothermal'
    save_location = '../Output/' + protein + '/'


    print("Reading input data")
    concs = np.loadtxt(conc_datafile)
    temp_fluo_data = np.loadtxt(fluo_datafile, skiprows=1, delimiter="\t")

    # Temperatures are in the first column, make these a separate array
    temperatures = temp_fluo_data[:,0]
    fluo_data = np.delete(temp_fluo_data, 0, 1)

    if (len(concs) != len(fluo_data[0, :])):
        print("Mismatch in concs and all_dat array lengths")
        quit()

    print("Proceeding with", len(temperatures), "temperatures")
    print("Proceeding with", len(concs), "ligand concentrations")

    print("Fitting thermal curves individually")
    fitting_params = fit_all_individual_thermal_curves(concs, temperatures, fluo_data, Cp, Tm_datafile)
    print("Writing parameters ")
    np.savetxt(parameter_datafile, fitting_params.transpose(), fmt='%15.5f',
               header="        Tm (oC)     dH (kcal/mol)       intercept_U       intercept_F           slope_U           slope_F",
               delimiter="   ", comments="")
    print("Plotting fits to thermal unfolding data")
    plot_fits_to_thermal_unfolding_data(concs, temperatures, fitting_params, fluo_data, Cp, save_location)
    print("Plotting Tm shifts")
    plot_Tm_shifts(concs, fitting_params[0, :], save_location)

    isothermal_temp_string = input('Enter comma separated temperatures for isothermals :')
    isothermal_temp_list = isothermal_temp_string.split(',')

    unique_concs = np.unique(concs)
    all_isothermal_data = np.empty([len(concs), len(isothermal_temp_list)])
    all_isothermal_fits = np.empty([len(unique_concs), len(isothermal_temp_list)])
    ft = open(KuKd_datafile, 'w')
    ft.write('Temp\tKu\tkD\n')

    for t_index in range(len(isothermal_temp_list)):
        binding_temperature = float(isothermal_temp_list[t_index])
        fraction_unfolded = calculate_fraction_unfolded(binding_temperature, concs, fitting_params, Cp)
        all_isothermal_data[:, t_index] = fraction_unfolded.transpose()
        (Ku, Kd, fit_values) = fit_isothermal_Kd(fraction_unfolded, concs, unique_concs, protein_conc)
        all_isothermal_fits[:, t_index] = fit_values.transpose()

        ft.write('%5.2foC\t%4.4f\t%4.4fuM\n' % (binding_temperature, Ku, Kd))

        fname_isothermal_data = isothermal_output_tag + "_T=" + isothermal_temp_list[t_index] + ".data.txt"
        fname_isothermal_fit = isothermal_output_tag + "_T=" + isothermal_temp_list[t_index] + ".fit.txt"

        fd = open(fname_isothermal_data, 'w')
        ff = open(fname_isothermal_fit, 'w')

        fd.write("    conc (uM)\tfraction_unfolded\n")
        ff.write("    conc (uM)\tfraction_unfolded\n")
        for c_index in range(len(concs)):
            fd.write(
                '{:13.5f}'.format(concs[c_index]) + "\t" + '{:17.5f}'.format(fraction_unfolded[c_index]) + "\n")
        for c_index in range(len(unique_concs)):
            ff.write(
                '{:13.5f}'.format(unique_concs[c_index]) + "\t" + '{:17.5f}'.format(fit_values[c_index]) + "\n")
        fd.close()
        ff.close()
    ft.close()

    print("Plotting isothermals")
    plot_isothermal_binding_fits(isothermal_temp_list, concs, all_isothermal_data, unique_concs, all_isothermal_fits, save_location)
    plt.close('all')
    print('----------------------------------\n')

