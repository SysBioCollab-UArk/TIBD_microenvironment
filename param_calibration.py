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
from itertools import combinations
import os


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
                      for expt in self.experiments]  # TODO: change this to xvar, rather than time

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
                for group in param_expts_map.get(p.name, [[expt for expt in self.experiments]]):
                    if prior[0] == 'uniform':
                        self.sampled_params_list.append(SampledParam(uniform, loc=np.log10(p_value) - 0.5 * prior[1],
                                                                     scale=prior[1]))
                    elif prior[0] == 'norm':
                        self.sampled_params_list.append(SampledParam(norm, loc=np.log10(p_value), scale=prior[1]))
                    else:
                        raise ValueError("Prior shape {} not recognized".format(prior[0]))
                    # save index of sampled_params_list for each experiment in the group
                    for n in range(self.n_experiments):
                        if self.experiments[n] in group:
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
        # NOTE: creating the array here and reusing it, rather than creating a new one every time the likelihood
        # function is accessed, should save some time
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

    def run(self, niterations=50000, nchains=3, multitry=False, gamma_levels=4,  adapt_gamma=True, history_thin=1,
            verbose=True, obs_labels=None, plot_results=True, plot_ll_args=None, plot_pd_args=None, plot_tc_args=None):

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
                                                   multitry=multitry,
                                                   gamma_levels=gamma_levels,
                                                   adapt_gamma=adapt_gamma,
                                                   history_thin=history_thin,
                                                   model_name='dreamzs_%dchain' % nchains,
                                                   verbose=verbose,
                                                   start=starts,
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
            logps_files = glob.glob('dreamzs*logps*')
            samples_files = glob.glob('dreamzs*params*')
            self.create_figures(logps_files, samples_files, obs_labels=obs_labels, plot_ll_args=plot_ll_args,
                                plot_pd_args=plot_pd_args, plot_tc_args=plot_tc_args)

    def create_figures(self, logps_files, samples_files, obs_labels=None, show_plots=False, plot_ll_args=None,
                       plot_pd_args=None, plot_tc_args=None):

        # process plotting function arguments
        _plot_ll_args = {'cutoff': 2}
        _plot_pd_args = {'labels': None, 'groups': None, 'cutoff': 2, 'save_plot': None}
        _plot_tc_args = {'tspans': None, 'xlabel': None, 'ylabels': None, 'leg_labels': None, 'separate_plots': True}
        if plot_ll_args is not None:
            _plot_ll_args.update(plot_ll_args)
        if plot_pd_args is not None:
            _plot_pd_args.update(plot_pd_args)
        if plot_tc_args is not None:
            _plot_tc_args.update(plot_tc_args)

        print('Plotting log-likelihoods')
        self.plot_log_likelihood(logps_files, **_plot_ll_args)

        print('Plotting parameter distributions')
        # get groups of parameters common for different sets of experiments
        param_groups = []
        filenames = []  # save info about which experiments each group is associated with in the filename
        expt_idxs = [n for n in range(self.n_experiments)]
        for n in expt_idxs:
            for combo in combinations(expt_idxs, n + 1):  # all possible combinations of expt groups
                param_list_1 = [set(self.samples_idxs[i]) for i in combo]
                param_list_2 = [set(self.samples_idxs[i]) for i in expt_idxs if i not in combo]
                if len(param_list_2) == 0:
                    param_list_2 = [set()]  # empty set
                # get parameters unique to this group of expts (might be none)
                param_diff = set.intersection(*param_list_1) - set.intersection(*param_list_2)
                if len(param_diff) > 0:
                    param_groups.append(sorted(list(param_diff)))
                    filename = 'fig_PyDREAM_histograms_EXPTS'
                    for i in combo:
                        filename += '_%s' % str(self.experiments[i])  # cast to a string in case an integer
                    filenames.append(filename)
        # get the parameter names to label the plots in the histogram
        param_labels = []
        for group in param_groups:
            param_labels.append([])
            for i in group:
                for s_idx_list in self.samples_idxs:
                    if i in s_idx_list:
                        param_labels[-1].append(
                            self.model.parameters[self.parameter_idxs[s_idx_list.index(i)]].name)
        # make the plots
        if _plot_pd_args['labels'] is None:
            _plot_pd_args['labels'] = param_labels
        if _plot_pd_args['groups'] is None:
            _plot_pd_args['groups'] = param_groups
        if _plot_pd_args['save_plot'] is None:
            _plot_pd_args['save_plot'] = filenames
        param_samples = self.plot_param_dist(samples_files, **_plot_pd_args)

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
        # create observable legend labels if obs_labels exists
        leg_labels = None if obs_labels is None else [[obs_labels.get(obs_name, obs_name) for obs_name
                                                       in self.observables[n]] for n in range(self.n_experiments)]
        # increase the number of time points for the simulations by a factor of 10
        tspans = _plot_tc_args.pop('tspans')  # pop tspans out of the dictionary since it's not passed as a kwarg below
        if tspans is None:
            tspans = [np.linspace(
                self.tdata[n][0], self.tdata[n][-1],
                int((self.tdata[n][-1] - self.tdata[n][0]) * 10 + 1))
                for n in range(self.n_experiments)]
        # make the plots
        # _plot_tc_args = {'tspans': None, 'xlabel': None, 'ylabels': None}
        if _plot_tc_args['xlabel'] is None:
            _plot_tc_args['xlabel'] = xlabel
        if _plot_tc_args['ylabels'] is None:
            _plot_tc_args['ylabels'] = ylabels
        if _plot_tc_args['leg_labels'] is None:
            _plot_tc_args['leg_labels'] = leg_labels
        self.plot_timecourses(self.model, tspans, self.sim_protocols, param_samples, self.parameter_idxs,
                              exp_data=self.raw_data, samples_idxs=self.samples_idxs, **_plot_tc_args)

        if show_plots:
            plt.show()

    @staticmethod
    def plot_log_likelihood(logps_files, cutoff=None, show_plot=False, save_plot=True):

        # get the path from the first file
        path, file = os.path.split(logps_files[0])

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
            log_ps.append(np.concatenate(tuple(np.load(
                os.path.join(path, 'dreamzs_%dchain_logps_chain_%d_%d.npy' % (len(chains), chain, it))).flatten()
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
            log_ps_var += np.var(log_ps[chain][burnin:]) / len(chains)  # mean of the variances, but that's fine
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

        # get the path from the first file
        path, file = os.path.split(sample_files[0])

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
        samples = np.concatenate(tuple(np.load(
            os.path.join(path, 'dreamzs_%dchain_sampled_params_chain_%d_%d.npy' %
                         (len(chains), chain, max(iterations))))[burnin:] for chain in chains))

        # if a likelihood cutoff is defined, load the likelihoods and remove samples that fall below the cutoff
        if cutoff is not None:
            log_ps = np.concatenate(tuple(np.load(
                os.path.join(path, 'dreamzs_%dchain_logps_chain_%d_%d.npy' %
                             (len(chains), chain, max(iterations))))[burnin:] for chain in chains))
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
                if save_plot is True:
                    filename = 'fig_PyDREAM_histograms'
                    suffix = '' if len(groups) == 1 else '_group_%d' % n
                    filename += suffix
                else:
                    if isinstance(save_plot, str):  # a string prefix has been passed
                        filename = save_plot
                        suffix = '' if len(groups) == 1 else '_group_%d' % n
                        filename += suffix
                    else:  # the last possibility is an array of filenames
                        filename = save_plot[n]
                plt.savefig(filename)

        if show_plot:
            plt.show()

        return samples


    @staticmethod
    def plot_timecourses(model, tspans, sim_protocols, param_samples, parameter_idxs, observables=None, exp_data=None,
                         samples_idxs=None, show_plot=False, save_plot=True, separate_plots=True, **kwargs):

        # if there's only one simulation to run, put the following objects inside a list of size 1, since the code
        # is set up to loop over the number of simulations
        if len(np.array(tspans).shape) == 1:
            tspans = [tspans]
        if len(np.array(sim_protocols).shape) == 0:
            sim_protocols = [sim_protocols]
        if observables is not None and len(np.array(observables).shape) == 1:
            observables = [observables]
        if samples_idxs is not None and len(np.array(samples_idxs).shape) == 1:
            samples_idxs = [samples_idxs]

        # experimental data
        n_experiments = 0
        if exp_data is not None:
            raw_data = exp_data if isinstance(exp_data, np.ndarray) \
                else np.genfromtxt(exp_data, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
            experiments = np.unique([d['expt_id'] for d in raw_data])
            n_experiments = len(experiments)
            if observables is None:
                observables = [np.unique([d['observable'] for d in raw_data if d['expt_id'] == expt])
                               for expt in experiments]
        elif observables is None:
            raise Exception('No experimental data or observables provided')

        # if sample indices are not provided, all parameters are used in all simulations
        if samples_idxs is None:
            samples_idxs = [[i for i in range(len(parameter_idxs))] for sim in range(len(tspans))]

        # error check
        if len(np.unique([len(tspans), len(sim_protocols), len(observables), len(samples_idxs)])) != 1:
            raise Exception(
                "The following arrays must all be equal length: 'tspans', 'sim_protocols', 'observables', "
                "'sample_idxs")

        # number of simulations
        n_sims = len(tspans)

        # process kwargs
        fill_between = kwargs.get('fill_between', (5, 95))
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
        colors = kwargs.get('colors',
                            [[cycle[i % 10] for i in range(len(obs))] for obs in observables] if separate_plots
                            else
                            [[cycle[(i + len(observables[n]) * n) % 10] for i in range(len(observables[n]))]
                             for n in range(n_sims)])  # different colors if sims plotted on same plot
        locs = kwargs.get('locs', [[0] * len(observables[n]) for n in range(n_sims)])
        xlabel = kwargs.get('xlabel', 'time')
        ylabels = kwargs.get('ylabels', [['amount'] * len(observables[n]) for n in range(n_sims)])
        leg_labels = kwargs.get('leg_labels', [[obs_name for obs_name in observables[n]] for n in range(n_sims)])
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

        # run simulations using only unique parameter samples
        samples_unique, counts = np.unique(param_samples, return_counts=True, axis=0)
        print('Running %d simulations' % len(samples_unique))

        # loop over simulations (experiments + perturbations)
        for n in range(n_sims):
            if n < n_experiments:
                print("Experiment '%s' (%d of %d)" % (str(experiments[n]), n+1, n_experiments))
            else:
                print("Perturbation %d of %d" % (n - n_experiments + 1, n_sims - n_experiments))
            outputs = []
            # run simulations
            for i, sample in enumerate(samples_unique):
                print(i, end=' ')
                param_values[parameter_idxs] = 10 ** sample[samples_idxs[n]]
                # TODO: let sim_protocol return a legend label for the simulated experiment
                outputs.append(sim_protocols[n](solver, tspans[n], param_values))
                if (i + 1) % 20 == 0:
                    print()
            print()
            # use 'counts' to generate full set of simulation outputs for correct weighting for plots
            outputs = np.repeat(outputs, counts, axis=0)
            # plot results
            for i, obs_name in enumerate(observables[n]):
                figname = '%s_sim_%s' % (obs_name, n) if separate_plots else obs_name
                plt.figure(num=figname, constrained_layout=True)
                # plot simulated data as a percent envelope
                yvals = np.array([output[obs_name] for output in outputs])
                yvals_min = np.percentile(yvals, fill_between[0], axis=0)
                yvals_max = np.percentile(yvals, fill_between[1], axis=0)
                plt.fill_between(tspans[n], yvals_min, yvals_max, alpha=0.25, color=colors[n][i],
                                 label=leg_labels[n][i])
                # plot experimental data
                if n < n_experiments:
                    time = [d['time'] for d in raw_data if d['observable'] == obs_name
                            and d['expt_id'] == experiments[n]]
                    avg = [d['average'] for d in raw_data if d['observable'] == obs_name
                           and d['expt_id'] == experiments[n]]
                    stderr = [d['stderr'] for d in raw_data if d['observable'] == obs_name
                              and d['expt_id'] == experiments[n]]
                    label = 'experiment' if n_experiments == 1 else 'experiment %s' % str(experiments[n])
                    plt.errorbar(time, avg, yerr=stderr, capsize=6, fmt='o', ms=8, mfc=colors[n][i], mec=colors[n][i],
                                 ecolor=colors[n][i], label=label)
                plt.xlabel(xlabel)
                plt.ylabel(ylabels[n][i])
                plt.legend(loc=locs[n][i])

                if save_plot is not None and separate_plots:
                    filename = 'fig_PyDREAM_tc_%s' % figname if save_plot is True else save_plot
                    plt.savefig(filename)

        if save_plot is not None and not separate_plots:
            # flatten the observables matrix and then pull out the unique names so we can loop over them
            obs_names = np.unique(np.array(observables).flatten())
            for obs_name in obs_names:
                # need to loop over the experiments in order to get the correct legend location
                for n in range(n_sims):
                    # check if the observable is in the list for this experiment. if it is, reorder the legend, save the
                    # file, and break out of the loop
                    if obs_name in observables[n]:
                        filename = 'fig_PyDREAM_tc_%s' % obs_name if save_plot is True else save_plot
                        # fix legend order
                        plt.figure(num=obs_name)
                        handles, labels = plt.gca().get_legend_handles_labels()
                        idx_order = []
                        for i in range(n_experiments):
                            idx_order += [j for j in range(i, len(handles), n_sims)]
                        idx_order += [j for j in range(len(handles)) if j not in idx_order]
                        plt.legend([handles[j] for j in idx_order], [labels[j] for j in idx_order],
                                   loc=locs[n][list(observables[n]).index(obs_name)])
                        plt.savefig(filename)
                        break

        if show_plot:
            plt.show()
