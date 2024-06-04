from pydream.core import run_dream
from pydream.parameters import SampledParam
from pydream.convergence import Gelman_Rubin
from pysb.simulator import ScipyOdeSimulator
from scipy.stats import norm, uniform
import numpy as np
import matplotlib.pyplot as plt
import re
import seaborn as sns
import glob


class ParameterCalibration(object):

    def __init__(self, model, exp_data_file, sim_protocols, priors=None, no_sample=None, default_prior=('norm', 2.0),
                 param_expts_map=None):

        self.model = model
        # if only one simulation protocol is given, put it inside a list of length 1
        self.sim_protocols = [sim_protocols] if len(np.array(sim_protocols).shape) == 0 else sim_protocols

        # read in experimental data
        self.raw_data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")

        # determine how many experiments there are and what the time points and observables are for each
        self.experiments = np.unique([d['expt_id'] for d in self.raw_data])
        self.n_experiments = len(self.experiments)
        self.observables = [np.unique([d['observable'] for d in self.raw_data if d['expt_id'] == expt])
                            for expt in self.experiments]
        self.tdata = [np.unique([d['time'] for d in self.raw_data if d['expt_id'] == expt])
                      for expt in self.experiments]

        # for a given experiment, there may be missing data points, so create masks that indicate for each observable in
        # each experiment which time points we have data for
        self.tspan_masks = []
        for expt, obs, tspan in zip(self.experiments, self.observables, self.tdata):
            self.tspan_masks.append({})
            for name in obs:
                self.tspan_masks[-1][name] = [False] * len(tspan)
                tdata = [d['time'] for d in self.raw_data if d['observable'] == name and d['expt_id'] == expt]
                for i in range(len(tspan)):
                    if tspan[i] in tdata:
                        self.tspan_masks[-1][name][i] = True

        # create the list of parameters to be sampled and save their indices
        priors = {} if priors is None else priors
        no_sample = [] if no_sample is None else no_sample
        param_expts_map = {} if param_expts_map is None else param_expts_map
        self.sampled_params_list = []
        self.parameter_idxs = []
        self.samples_idxs = [[] for n in range(self.n_experiments)]  # this is for when params vary between expts
        for i, p in enumerate(self.model.parameters):
            if p.name not in no_sample:
                self.parameter_idxs.append(i)
                p_value = self.model.parameters[p.name].value
                prior = default_prior if p.name not in list(priors.keys()) else priors[p.name]
                print("Will sample parameter {} with {} prior around log10({}) = {}".format(p.name, prior,
                                                                                            p_value, np.log10(p_value)))
                # loop over groups of experiments (default is one group with all experiments)
                for group in param_expts_map.get(p.name, [[i for i in range(self.n_experiments)]]):
                    if prior[0] == 'uniform':
                        self.sampled_params_list.append(SampledParam(uniform, loc=np.log10(p_value) - 0.5 * prior[1],
                                                                     scale=prior[1]))
                    elif prior[0] == 'norm':
                        self.sampled_params_list.append(SampledParam(norm, loc=np.log10(p_value), scale=prior[1]))
                    else:
                        raise ValueError("Prior shape {} not recognized".format(prior[0]))
                    # save index of sampled_params_list for each experiment in the group
                    for n in range(self.n_experiments):
                        if n in group:
                            self.samples_idxs[n].append(len(self.sampled_params_list) - 1)

        # create normal distributions around each data point to use in the likelihood function
        self.like_data = []
        for expt, obs, tspan in zip(self.experiments, self.observables, self.tdata):
            self.like_data.append({})
            for name in obs:
                data = [d for d in self.raw_data if d['observable'] == name and d['expt_id'] == expt]
                avg = []
                se = []
                for t in tspan:  # add data in same order as the time points in tspan
                    if t in [d['time'] for d in data]:
                        avg.append([d['average'] for d in data if d['time'] == t][0])
                        se.append([d['stderr'] for d in data if d['time'] == t][0])
                self.like_data[-1][name] = norm(loc=avg, scale=se)

        # store the full set of parameter values
        self.param_values = np.array([p.value for p in self.model.parameters])

        # create the array for storing the likelihood values
        # creating the array here and reusing it, rather than creating a new one every time the likelihood function is
        # accessed, should save some time
        self.logp_data = [0.] * self.n_experiments

        # create the simulator
        self.solver = ScipyOdeSimulator(self.model)

    # likelihood function
    def likelihood(self, position):
        y = np.copy(position)
        for n in range(self.n_experiments):
            self.logp_data[n] = 0.  # reinitialize this to zero
            self.param_values[self.parameter_idxs] = 10 ** y[self.samples_idxs[n]]
            # run simulation
            output = self.sim_protocols[n](self.solver, self.tdata[n], self.param_values)
            # calculate log-likelihood
            for obs in self.like_data[n].keys():
                self.logp_data[n] += np.sum(self.like_data[n][obs].logpdf(output[obs][self.tspan_masks[n][obs]]))
            if np.isnan(self.logp_data[n]):
                self.logp_data[n] = -np.inf
        return sum(self.logp_data)

    # def likelihood(position):
    #     y = np.copy(position)
    #     logp_data = [0] * n_experiments
    #     for n in range(n_experiments):
    #         param_values[rates_mask] = 10 ** y
    #         # equilibration
    #         equil = solver.run(tspan=np.linspace(-500, 0, 2), param_values=param_values)
    #         # add tumor cells
    #         initials = equil.species[-1]
    #         idx_tumor = [str(sp) for sp in model.species].index('Tumor()')  # get index of Tumor species
    #         initials[idx_tumor] = 1  # fM
    #         sim = solver.run(tspan=tspan[n], param_values=param_values, initials=initials).all
    #         # calculate log-likelihood
    #         for sp in like_data[n].keys():
    #             logp_data[n] += np.sum(like_data[n][sp].logpdf(sim[sp][tspan_mask[n][sp]]))
    #         if np.isnan(logp_data[n]):
    #             logp_data[n] = -np.inf
    #     return sum(logp_data)

    def run(self, niterations=50000, nchains=3, multitry=False, gamma_levels=4,  adapt_gamma=True,
            history_thin=1, verbose=True, plot_results=True):

        sampled_params, log_ps = run_dream(parameters=self.sampled_params_list,
                                           likelihood=self.likelihood,
                                           niterations=niterations,
                                           nchains=nchains,
                                           multitry=multitry,
                                           gamma_levels=gamma_levels,
                                           adapt_gamma=adapt_gamma,
                                           history_thin=history_thin,
                                           model_name='dreamzs_%dchain' % nchains,
                                           verbose=verbose)
        total_iterations = niterations
        burnin = int(total_iterations / 2)
        # Save sampling output (sampled parameter values and their corresponding logps).
        for chain in range(len(sampled_params)):
            np.save('dreamzs_%dchain_sampled_params_chain_%d_%d' %
                    (nchains, chain, total_iterations), sampled_params[chain])
            np.save('dreamzs_%dchain_logps_chain_%d_%d' % (nchains, chain, total_iterations), log_ps[chain])
        old_samples = sampled_params

        # Check convergence and continue sampling if not converged
        GR = Gelman_Rubin(sampled_params)
        print('At iteration: ', total_iterations, ' GR = ', GR)
        np.savetxt('dreamzs_%dchain_GelmanRubin_iteration_%d.txt' % (nchains, total_iterations), GR)
        if np.any(GR > 1.2):
            starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
            converged = False
            while not converged:
                total_iterations += niterations
                burnin += niterations
                sampled_params, log_ps = run_dream(parameters=self.sampled_params_list,
                                                   likelihood=self.likelihood,
                                                   niterations=niterations,
                                                   nchains=nchains,
                                                   start=starts,
                                                   multitry=False,
                                                   gamma_levels=4,
                                                   adapt_gamma=True,
                                                   history_thin=1,
                                                   model_name='dreamzs_%dchain' % nchains,
                                                   verbose=True,
                                                   restart=True)
                for chain in range(len(sampled_params)):
                    np.save('dreamzs_%dchain_sampled_params_chain_%d_%d' %
                            (nchains, chain, total_iterations), sampled_params[chain])
                    np.save('dreamzs_%dchain_logps_chain_%d_%d' % (nchains, chain, total_iterations), log_ps[chain])
                old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
                GR = Gelman_Rubin(old_samples)
                print('At iteration: ', total_iterations, ' GR = ', GR)
                np.savetxt('dreamzs_%dchain_GelmanRubin_iteration_%d.txt' % (nchains, total_iterations), GR)
                if np.all(GR < 1.2):
                    converged = True

        if plot_results:

            print('Plotting log-likelihoods')
            logps_files = glob.glob('dreamzs*logps*')
            self.plot_log_likelihood(logps_files, cutoff=2)

            print('Plotting parameter distributions')
            samples_files = glob.glob('dreamzs*params*')
            param_samples = (
                self.plot_param_dist(samples_files, labels=[self.model.parameters[i].name for i in self.parameter_idxs],
                                     cutoff=2))

            print('Plotting time courses')
            # get the time units for the x-label
            time_units = np.unique([d['time_units'] for d in self.raw_data])
            if len(time_units) > 1:
                print("WARNING: Multiple time units included in the data file. Using default 'time' for xlabel.")
            xlabel = 'time (%s)' % time_units[0] if len(time_units) == 1 else 'time'
            # get the amount units for the y-labels, which can be different for different observables
            ylabels = []
            for n in range(self.n_experiments):
                ylabels.append([])
                for obs in self.observables[n]:
                    amount_units = np.unique([d['amount_units'] for d in self.raw_data if d['observable'] == obs
                                              and d['expt_id'] == self.experiments[n]])
                    if len(amount_units) > 1:
                        print("WARNING: Multiple amount units included in the data file for observable %s (expt id %s)."
                              " Using default 'amount' for ylabel.")
                    ylabel = 'amount (%s)' % amount_units[0] if len(amount_units) == 1 else 'amount'
                    ylabels[-1].append(ylabel)
            # increase the number of time points for the simulations by a factor of 10
            tspans = [np.linspace(self.tdata[n][0], self.tdata[n][-1],
                                  int((self.tdata[n][-1] - self.tdata[n][0]) * 10 + 1))
                      for n in range(self.n_experiments)]
            self.plot_timecourses(self.model, self.raw_data, self.sim_protocols, tspans, param_samples,
                                  self.parameter_idxs, samples_idxs=self.samples_idxs, xlabel=xlabel, ylabels=ylabels)

    @staticmethod
    def plot_log_likelihood(logps_files, cutoff=None, show_plot=False, save_plot=True):

        # get chains and iterations from file names
        chains = []
        iterations = []
        for file in logps_files:
            m = re.search(r'chain_(\d+)_(\d+).npy$', file)
            chains.append(int(m.group(1)))
            iterations.append(int(m.group(2)))
        chains = np.unique(chains)
        iterations = np.unique(iterations)

        # read the likelihoods from the files
        log_ps = []
        for chain in chains:
            log_ps.append(np.concatenate(
                tuple(np.load('dreamzs_%dchain_logps_chain_%d_%d.npy' % (len(chains), chain, it)).flatten()
                      for it in iterations)))

        # plot the likelihoods
        plt.figure(constrained_layout=True)
        # calculate mean and variance for last half of steps from last iteration
        burnin = int(max(iterations)) - int(min(iterations) / 2)
        log_ps_max = -np.inf
        log_ps_mean = 0
        log_ps_var = 0
        for chain in chains:
            plt.plot(range(len(log_ps[chain])), log_ps[chain], label='chain %d' % chain)
            log_ps_max = np.max(log_ps[chain]) if log_ps_max < np.max(log_ps[chain]) else log_ps_max
            log_ps_mean += np.mean(log_ps[chain][burnin:]) / len(chains)
            log_ps_var += np.var(log_ps[chain][burnin:]) / len(chains)  # this is the mean of the variances, but that's fine
        # plt.axvline()
        top = np.ceil(log_ps_mean + 5 * np.sqrt(log_ps_var))
        bottom = np.floor(log_ps_mean - 20 * np.sqrt(log_ps_var))
        print('max: %g, mean: %g, sdev: %g, top: %g, bottom: %g' %
              (log_ps_max, log_ps_mean, np.sqrt(log_ps_var), top, bottom))
        if cutoff is not None:
            plt.axhline(log_ps_mean - cutoff * np.sqrt(log_ps_var), color='k', ls='--', lw=2)
        plt.ylim(bottom=bottom, top=top)
        plt.xlabel('iteration')
        plt.ylabel('log-likelihood')
        plt.legend(loc=0)

        if show_plot:
            plt.show()
        if save_plot is not None:
            filename = 'fig_PyDREAM_log_ps' if save_plot is True else save_plot
            plt.savefig(filename)

    @staticmethod
    def plot_param_dist(sample_files, labels, groups=None, cutoff=None, show_plot=False, save_plot=True, **kwargs):

        # get chains and iterations from file names
        chains = []
        iterations = []
        for file in sample_files:
            m = re.search(r'chain_(\d+)_(\d+).npy$', file)
            chains.append(int(m.group(1)))
            iterations.append(int(m.group(2)))
        chains = np.unique(chains)
        iterations = np.unique(iterations)

        # read in parameter values and likelihoods from files
        burnin = int(min(iterations) / 2)
        samples = np.concatenate(
            tuple(np.load('dreamzs_%dchain_sampled_params_chain_%d_%d.npy' %
                          (len(chains), chain, max(iterations)))[burnin:] for chain in chains))

        # if a likelihood cutoff is defined, load the likelihoods and remove samples that fall below the cutoff
        if cutoff is not None:
            log_ps = np.concatenate(
                tuple(np.load('dreamzs_%dchain_logps_chain_%d_%d.npy' % (len(chains), chain, max(iterations)))[burnin:]
                      for chain in chains))
            log_ps_mean = np.mean(log_ps)
            log_ps_sdev = np.std(log_ps)
            keep_idxs = [i for i in range(len(samples)) if log_ps[i] > log_ps_mean - cutoff * log_ps_sdev]
            samples = samples[keep_idxs]

        # plot histograms
        if groups is None:
            groups = [[i for i in range(len(samples[0]))]]
            labels = [labels]
        for n, label, group in zip(range(len(labels)), labels, groups):
            ndims = len(group)  # number of dimensions (i.e., parameters)
            # set plot parameters
            fscale = np.ceil(ndims / 16)
            figsize = kwargs.get('figsize', fscale * np.array([6.4, 4.8]))
            labelsize = kwargs.get('labelsize', 10 * max(1, (2/5 * fscale)))
            fontsize = kwargs.get('fontsize', 10 * max(1, (3/5 * fscale)))
            ncols = kwargs.get('ncols', int(np.ceil(np.sqrt(ndims))))
            nrows = int(np.ceil(ndims/ncols))
            # create figure
            colors = sns.color_palette(n_colors=ndims)
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=False, constrained_layout=True,
                                    figsize=figsize)
            row = 0
            col = 0
            for dim in range(ndims):
                print(group[dim], end=' ')
                sns.distplot(samples[:, group[dim]], color=colors[dim], norm_hist=True, ax=axs[row][col])
                axs[row][col].set_yticklabels([])
                axs[row][col].set_ylabel(None)
                axs[row][col].set_title(label[dim], fontsize=labelsize)
                axs[row][col].tick_params(axis='x', labelsize=labelsize)
                col += 1
                if col % ncols == 0:
                    col = 0
                    row += 1
            print()
            fig.supxlabel(r'log$_{10}$ value', fontsize=fontsize)
            fig.supylabel('Density', fontsize=fontsize)
            # delete extra plots
            if col > 0:
                while col < ncols:
                    fig.delaxes(axs[row][col])
                    col += 1
            # save plots
            if save_plot is not None:
                filename = 'fig_PyDREAM_histograms' if save_plot is True else save_plot
                suffix = '' if len(groups) == 1 else '_group_%d' % n
                plt.savefig(filename + suffix)

        if show_plot:
            plt.show()

        return samples


    @staticmethod
    def plot_timecourses(model, exp_data, sim_protocols, tspans, param_samples, parameter_idxs, samples_idxs=None,
                         show_plot=False, save_plot=True, **kwargs):

        # if there's only one experiment, put the following objects inside a list of size 1, since the following code
        # is set up to loop over the number of experiments
        if len(np.array(sim_protocols).shape) == 0:
            sim_protocols = [sim_protocols]
        if len(np.array(tspans).shape) == 1:
            tspans = [tspans]
        if samples_idxs is not None and len(np.array(samples_idxs).shape) == 1:
            samples_idxs = [samples_idxs]

        # experimental data
        raw_data = exp_data if isinstance(exp_data, np.ndarray) \
            else np.genfromtxt(exp_data, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
        experiments = np.unique([d['expt_id'] for d in raw_data])
        n_experiments = len(experiments)
        observables = [np.unique([d['observable'] for d in raw_data if d['expt_id'] == expt]) for expt in experiments]

        # process kwargs
        fill_between = kwargs.get('fill_between', (5, 95))
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
        colors = kwargs.get('colors', [[cycle[i % 10] for i in range(len(observables[n]))]
                                       for n in range(n_experiments)])
        locs = kwargs.get('locs', [[0] * len(observables[n]) for n in range(n_experiments)])
        xlabel = kwargs.get('xlabel', 'time')
        ylabels = kwargs.get('ylabels', [['amount'] * len(observables[n]) for n in range(n_experiments)])
        # make sure arrays are 2D
        if len(np.array(colors).shape) == 1:
            colors = [colors]
        if len(np.array(locs).shape) == 1:
            locs = [locs]
        if len(np.array(ylabels).shape) == 1:
            ylabels = [ylabels]

        # store original parameter values and create simulator
        param_values = np.array([p.value for p in model.parameters])
        solver = ScipyOdeSimulator(model, verbose=False)

        # if sample indices are not provided, all parameters are used in all experiments
        if samples_idxs is None:
            samples_idxs = [[i for i in range(len(parameter_idxs))] for expt in experiments]

        # loop over experiments and run simulations
        # only use unique parameter samples
        samples_unique, counts = np.unique(param_samples, return_counts=True, axis=0)
        print('Running %d simulations' % len(samples_unique))
        for n in range(n_experiments):
            print('Experiment %s' % experiments[n])
            outputs = []
            # run simulations
            for i, sample in enumerate(samples_unique):
                print(i, end=' ')
                param_values[parameter_idxs] = 10 ** sample[samples_idxs[n]]
                outputs.append(sim_protocols[n](solver, tspans[n], param_values))
            print()
            # use 'counts' to generate full set of simulation outputs for correct weighting for plots
            outputs = np.repeat(outputs, counts, axis=0)
            # plot results
            for i, obs_name in enumerate(observables[n]):
                plt.figure(num='%s_%s' % (obs_name, experiments[n]), constrained_layout=True)
                # plot simulated data as a percent envelope
                yvals = np.array([output[obs_name] for output in outputs])
                yvals_min = np.percentile(yvals, fill_between[0], axis=0)
                yvals_max = np.percentile(yvals, fill_between[1], axis=0)
                plt.fill_between(tspans[n], yvals_min, yvals_max, alpha=0.25, color=colors[n][i], label=obs_name)
                # plot experimental data
                time = [d['time'] for d in raw_data if d['observable'] == obs_name and d['expt_id'] == experiments[n]]
                avg = [d['average'] for d in raw_data if d['observable'] == obs_name and d['expt_id'] == experiments[n]]
                stderr = [d['stderr'] for d in raw_data if
                          d['observable'] == obs_name and d['expt_id'] == experiments[n]]
                label = 'experiment' if n_experiments == 1 else 'experiment %d' % n
                plt.errorbar(time, avg, yerr=stderr, capsize=6, fmt='o', ms=8, mfc=colors[n][i], mec=colors[n][i],
                             ecolor='k', label=label)
                plt.xlabel(xlabel)
                plt.ylabel(ylabels[n][i])
                plt.legend(loc=locs[n][i])

                if save_plot is not None:
                    filename = 'fig_PyDREAM_tc_%s_expt_%s' % (obs_name, experiments[n]) if save_plot is True \
                        else save_plot
                    plt.savefig(filename)

        if show_plot:
            plt.show()
