def BE_adaptive(th=0.0017, Pd_0=0.0226867, Pd_1=0.0146142, G=5, F=3, R=1e4, K=4, t_pi=140, IF=150e6, tau0=25,
                T2=3000, window_size=350, majority=False, adaptive=True, testTime=False, check_iteration=False,
                cut_theta_loop=True, likelihood=0, msb=1, delta=0.1, amplitude=0.5, precision=0.1, theta_num=5):
    
    tau_vec = [int(arg) for arg in tau0 * (2 ** np.linspace(K, 0, K + 1))]  # tau_vec goes from short to long

    meas_len = 2000
    f = IF
    pi_time = t_pi
    Th = th * R
    alpha = (Pd_0 + Pd_1) / 2
    v_alpha = (-Pd_0 + Pd_1) / 2  # This is V*alpha without the exponent of tau
    fB_min = float(-(1 / (2 * tau0)))  # After multiplication by units of tau
    fB_max = float((1 / (2 * tau0)))  # After multiplication by units of tau
    dfB = float(delta * 1e-3)  # After multiplication by units of tau
    fB_vec = np.arange(fB_min, fB_max + dfB / 2, dfB)
    one_over_R = 1 / R
    norm_msb = msb
    Pf_py = [1 / len(fB_vec)] * int(len(fB_vec))
    Pf_estimator_py = np.dot(Pf_py, fB_vec)

    with program() as BE_adaptive:

        ind_pf = declare(int)
        ind_tau = declare(int, value=-1)
        tau = declare(int)
        n = declare(int)  # index for collecting counts
        i = declare(int)  # index for not exceeding the number of phases per tau - G+(K-k)F
        j = declare(int)  # index for iterating over all R and for saving Pf
        u = declare(bool)  # for majority 0/1 compared to Th
        U = declare(fixed)  # after adding 0.5 to u
        theta = declare(fixed)
        last_theta = declare(fixed, size=theta_num)  # vector holding all 5 last thetas to do the cut_theta option
        theta_i = declare(int)  # index to run over thetas in last_theta vec
        photoluminescence_1 = declare(int)
        rel_photoluminescence_1 = declare(int)
        result1 = declare(int, size=4)
        result_len1 = declare(int)
        result2 = declare(int, size=4)
        result_len2 = declare(int)
        one_over_2m = declare(fixed, value=[1.0 / (2 * (G + (K - k) * F)) for k in range(K + 1)])  # phase addition to non-adaptive option
        one_over_4m = declare(fixed, value=[1.0 / (4 * (G + (K - k) * F)) for k in range(K + 1)])
        
        k = declare(int, value=K)

        B = declare(fixed, value=[arg for arg in np.exp(-(np.array(tau_vec) / (T2 / 4)) ** 2)])  # the exponent part of V
        C = declare(fixed)
        fB = declare(fixed)
        Pf = declare(fixed, value=[1 / len(fB_vec)] * int(len(fB_vec)))  # w2pi
        one_over_norm = declare(fixed, value=1.0)
        a = declare(fixed)  # amplitude of MW
        assign(a, amplitude)
        cut_flag = declare(int, value=1)  # for the cu_theta option

        if majority:
            msb = declare(int)  # for majority normalization
            right_shift = declare(int)
            left_shift = declare(int)

        else:
            Pd1 = declare(fixed, size=len(fB_vec))
            ##

            sigma = declare(int)

            one_over_sigma = declare(fixed)

            ##
            int_Pf = declare(int, value=[int(np.log2(1 / len(fB_vec)))] * int(len(fB_vec)))  # int part of log2 of Pf
            argument = declare(int, value=0)  # the argument of the exp of Pf (practically after log)
            squared = declare(int, value=0)  # argument squared
            max_value = declare(int, value=-2 ** 31)  # for normalization
            max_idx = declare(int, value=0)  # for normalization
            peak_counter = declare(int, value=0)  # for normalization after exp
            one_over_peak_counter = declare(fixed)  # for normalization after exp
            one = declare(int, value=1)

        if adaptive:
            fB_vec_qua = declare(fixed, value=[fB_vec[x] for x in range(len(fB_vec))])
            Pf_estimator = declare(fixed, value=Pf_estimator_py)
            norm = declare(fixed, value=1.0)
            r_over_R = declare(fixed)
            one_minus_two_rR = declare(fixed)
            # Initial average values:
            a_two_pi = declare(fixed, value=(alpha ** 2 + (1 - 2 * alpha) * alpha) * 2 * np.pi)
            b = declare(fixed, value=(1 - 2 * alpha) * v_alpha)
        if not majority:  # Declare variables for 20.12 fixed calculation
            zero = declare(fixed, value=0.0)
            two = declare(fixed, value=2.0)
            int_part = declare(int)
            pow2_int_part = declare(fixed)

        update_frequency('qubit', f)
        ###############
        # Tau loop:
        ###############
        
        with for_(tau, tau0* (2 ** K), tau >= tau0 , tau >> 1):

            if not adaptive:
                assign(theta, (-1 * one_over_2m[k]))
            assign(ind_tau, ind_tau + 1)

            assign(i, 0)
            ###############
            # zero last_theta_vec:
            ###############
            if (cut_theta_loop and adaptive):
                with for_(theta_i, 0, theta_i < theta_num, theta_i + 1):
                    assign(last_theta[theta_i], 0)
            
            ###############
            # Theta loop:
            ###############
            with while_(i < (G + (K - k) * F)): # loop of phases
                assign(photoluminescence_1, 0)
                assign(rel_photoluminescence_1, 0)

                ###############
                # Theta Update:
                ###############
                if adaptive:
                    assign(theta,
                           Cast.mul_fixed_by_int(Pf_estimator, tau) * norm + 0.25 - 0.5 * Math.div(b, a_two_pi) * B[ind_tau])
                else:
                    assign(theta, theta + one_over_2m[k])
                save(theta, 'theta')

                ###############
                # calc Pd1:
                ###############

                if not majority:  # Initialize Pd1 vector
                    assign(ind_pf, 0)
                    with for_(fB, fB_min, fB < fB_max + dfB / 2, fB + dfB):
                        assign(C, Math.cos2pi(Cast.mul_fixed_by_int(fB, tau) - theta + likelihood)) 
                        assign(Pd1[ind_pf], (alpha + v_alpha * B[ind_tau] * C))

                        if check_iteration:
                            save(Pd1[ind_pf],'Pd1')
                        assign(ind_pf, ind_pf + 1)

                ############
                # Test Time:
                ############
                if testTime:
                    frame_rotation_2pi(theta, 'qubit')
                    play('MW_pi', 'qubit', duration=((pi_time // 2) // 4))  # pi/2 pulse
                    reset_frame('qubit')

                ###############
                # R loop:
                ###############
                align('qubit', 'readout1', 'readout2', 'laser')
                with for_(j, 0, j < R, j + 1):  # loop over R repetitions

                    #########
                    # Ramsey:
                    #########
                    play('MW_pi' * amp(a), 'qubit', duration=((pi_time // 2) // 4))  # pi/2 pulse
                    wait(tau, 'qubit')  # wait time tau
                    frame_rotation_2pi(theta, 'qubit')
                    play('MW_pi' * amp(a), 'qubit', duration=((pi_time // 2) // 4))  # pi/2 pulse rotated by theta
                    reset_frame('qubit')

                    ##############
                    # Measurement:
                    ##############
                    align('qubit', 'readout1', 'readout2', 'laser')
                    wait(100 // 4, 'laser')
                    play('laser_pulse', 'laser', duration=(meas_len - 400) // 4)
                    measure('readout', 'readout1', None, time_tagging.analog(result1, meas_len, targetLen=result_len1))
                    measure('readout', 'readout2', None, time_tagging.analog(result2, meas_len, targetLen=result_len2))

                    ###########
                    # collect PL:
                    ###########
                    with for_(n, 0, n < result_len1, n + 1):
                        assign(photoluminescence_1, Util.cond((result1[n] > 220) & (result1[n] < (220 + window_size)),
                                                              photoluminescence_1 + 1, photoluminescence_1))

                    with for_(n, 0, n < result_len2, n + 1):
                        assign(photoluminescence_1, Util.cond((result2[n] > 220) & (result2[n] < (220 + window_size)),
                                                              photoluminescence_1 + 1, photoluminescence_1))

                    with for_(n, 0, n < result_len1, n + 1):
                        assign(rel_photoluminescence_1,
                               Util.cond((result1[n] > (1700 - window_size)) & (result1[n] < 1700),
                                         rel_photoluminescence_1 + 1, rel_photoluminescence_1))

                    with for_(n, 0, n < result_len2, n + 1):
                        assign(rel_photoluminescence_1,
                               Util.cond((result2[n] > (1700 - window_size)) & (result2[n] < 1700),
                                         rel_photoluminescence_1 + 1, rel_photoluminescence_1))

                ######################
                # Bayesian Estimation:
                ######################
                align('qubit', 'readout1', 'readout2', 'laser')
                ############
                # Test Time:
                ############
                if testTime:
                    frame_rotation_2pi(theta, 'qubit')
                    play('MW_pi', 'qubit', duration=((pi_time // 2) // 4))  # pi/2 pulse
                    reset_frame('qubit')

                ############
                # Majority BE:
                ############
                if majority:
                    assign(u, (photoluminescence_1) <= Th)
                    assign(U, Cast.to_fixed(u) - 0.5)
                    assign(ind_pf, 0)
                    with for_(fB, fB_min, fB < fB_max + dfB / 2, fB + dfB):
                        assign(C, Math.cos2pi(Cast.mul_fixed_by_int(fB, tau) - theta + likelihood))
                        assign(Pf[ind_pf], (0.5 + U * B[ind_tau] * C) * Pf[ind_pf])
                        assign(ind_pf, ind_pf + 1)
                    if check_iteration:
                        with for_(j, 0, j < Pf.length(), j + 1):
                            save(Pf[j], "Pf_vec_before_norm")

                    ############
                    # normalization of majority:
                    ############
                    assign(one_over_norm, Math.sum(Pf))
                    assign(msb, Math.msb(one_over_norm))
                    assign(right_shift, Util.cond(msb > norm_msb, msb - norm_msb, 0))
                    assign(left_shift, Util.cond(msb < (norm_msb - 1), norm_msb - msb - 1, 0))

                    with for_(ind_pf, 0, ind_pf < Pf.length(), ind_pf + 1):
                        assign(Pf[ind_pf], Pf[ind_pf] >> right_shift)
                        assign(Pf[ind_pf], Pf[ind_pf] << left_shift)
                

                ############
                # batching BE:
                ############
                else:

                    assign(ind_pf, 0)
                    assign(max_value, -2 ** 31)  # Lowest value possible
                    assign(max_idx, 0)  # Lowest value possible

                    with for_(fB, fB_min, fB < fB_max + dfB / 2, fB + dfB):
                        # These shifts, left and right, change the number to be represented by a fixed 20.12 instead of
                        # 4.28, this should give good precision for what follows.
                        assign(argument,((photoluminescence_1 << 12) - Cast.mul_int_by_fixed((int(R) << 12), Pd1[ind_pf])) >> 12)
                        assign(squared, argument * argument)
                        ##

                        assign(sigma,Cast.mul_int_by_fixed(int(R),Pd1[ind_pf]))
                        save(sigma,'sigma')

                        assign(sigma,Cast.mul_int_by_fixed(sigma,(1-Pd1[ind_pf])))
                        save(sigma,'sigma')

                        assign(sigma,2*sigma)

                        assign(one_over_sigma,Math.div(one,sigma))

                        ##
                        ##

                        assign(squared, Cast.mul_int_by_fixed(squared,one_over_sigma))

                        ##
                        assign(int_Pf[ind_pf], -Math.abs(squared) + int_Pf[ind_pf])
                        assign(max_idx, Util.cond((int_Pf[ind_pf] > max_value), ind_pf, max_idx))
                        assign(max_value, int_Pf[max_idx])
                        assign(ind_pf, ind_pf + 1)
                        if check_iteration:
                            save(argument, 'argument')
                            save(squared, 'squared')
                            save(sigma,'sigma')
                            save(one_over_sigma,'one_over_sigma')
                    if check_iteration:
                        with for_(j, 0, j < int_Pf.length(), j + 1):
                            save(int_Pf[j], "int_log_vec")

                    ############
                    # Normalization of log Pf:
                    ############
                    assign(peak_counter, 0)
                    with for_(ind_pf, 0, ind_pf < int_Pf.length(), ind_pf + 1):
                        assign(int_Pf[ind_pf], int_Pf[ind_pf] - max_value)
                        assign(peak_counter, Util.cond(Math.abs(int_Pf[ind_pf]) == 0, peak_counter + 1, peak_counter))
                        if check_iteration:
                            save(max_value,'int_max')
                            save(max_idx,'ind_max')
                        

                    if check_iteration:
                        with for_(j, 0, j < int_Pf.length(), j + 1):
                            save(int_Pf[j], "int_log_vec_norm")

                    ############
                    # get Pf from log and Normalization of exp Pf:
                    ############
                    if check_iteration:
                        save(peak_counter,'peak_counter')
                    assign(one_over_peak_counter, Math.div(one, peak_counter))
                    with for_(ind_pf, 0, ind_pf < int_Pf.length(), ind_pf + 1):
                        # This part converts log2(Pf), which is represented by 20.12, to Pf, which is fixed 4.28
                        
                        assign(int_part, int_Pf[ind_pf])
                        assign(int_part, Util.cond(Math.abs(int_part) > 28, 28, Math.abs(int_part)))
                        assign(pow2_int_part, two >> (int_part + 1))
                        assign(Pf[ind_pf], pow2_int_part * one_over_peak_counter)
                        if check_iteration:
                            save(int_part,'int_part')
                            save(pow2_int_part,'pow2_int_part')

                if check_iteration:
                    with for_(j, 0, j < Pf.length(), j + 1):
                        save(Pf[j], "Pf_vec")

                ###############
                # Theta Update:
                ###############
                if adaptive:
                    assign(r_over_R, Cast.mul_fixed_by_int(one_over_R, photoluminescence_1))
                    assign(one_minus_two_rR, 1 - 2 * r_over_R)
                    assign(a_two_pi, ((r_over_R * r_over_R + (
                                one_minus_two_rR * alpha)) * 2 * np.pi))  # times 2 and not 4 because of overflows
                    assign(b, one_minus_two_rR * v_alpha)
                    assign(one_over_norm, Math.sum(Pf))
                    assign(norm, Math.inv(one_over_norm))
                    assign(Pf_estimator, Math.dot(fB_vec_qua, Pf))
                    if check_iteration:
                        save(a_two_pi,'a')
                        save(b,'b')
                        save(one_over_norm,'one_over_norm')
                        save(norm,'norm')
                        save(Pf_estimator,'Pf_estimator')
                    ############
                    # check if theta converged and update i for the theta_loop:
                    ############
                    if cut_theta_loop:
                        assign(theta_i, 0)
                        assign(cut_flag, 1)
                        with for_(theta_i, 0, theta_i < theta_num, theta_i+1):
                            assign(cut_flag, Util.cond(Math.abs(theta - last_theta[theta_i]) > precision, 0, 1))
                        with if_(cut_flag == 1):
                            assign(i, (G + (K - k) * F))
                        with else_():
                            assign(i, i + 1)
                            with for_(theta_i, 0, theta_i < theta_num - 1, theta_i + 1):
                                assign(last_theta[theta_i], last_theta[theta_i + 1])
                            assign(last_theta[theta_num - 1], theta)
                    else:
                        assign(i, i + 1)

                else:
                    assign(i, i + 1)

                save(photoluminescence_1, 'PL_1')
                save(rel_photoluminescence_1, 'PL_rel_1')
                save(tau, 'tau')

            
            assign(k, k - 1)


        with for_(j, 0, j < Pf.length(), j + 1):
            save(Pf[j], "Pf")

    return BE_adaptive
