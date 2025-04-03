from pydream.core import run_dream
from pydream.parameters import SampledParam
from pydream.convergence import Gelman_Rubin
from scipy.stats import norm, uniform
import numpy as np
import matplotlib.pyplot as plt
import re
import seaborn as sns
import glob
from itertools import combinations
import os
import csv
import multiprocessing
import math


# Function for timeout option
def run_simulation(sim_protocol, tspan, p_values):
    return sim_protocol.run(tspan, p_values)


class SimulationProtocol(object):
    def __init__(self, solver, t_equil=None):
        self.solver = solver
        self.t_equil = t_equil

    # default simulation protocol (can be overwritten)
    def run(self, tspan, param_values):
        if self.t_equil is not None:
            out = self.solver.run(tspan=np.linspace(-self.t_equil + tspan[0], tspan[0], 2),
                                  param_values=param_values)
            initials = out.species[-1]
        else:
            initials = None
        output = self.solver.run(tspan=tspan, param_values=param_values, initials=initials).all
        return output


class ParameterCalibration(object):

    def __init__(self, model, exp_data_file, sim_protocols, priors=None, no_sample=None, default_prior=('norm', 2.0),
                 param_expts_map=None):

        self.model = model
        # if only one simulation protocol is given, put it inside a list of length 1
        self.sim_protocols = [sim_protocols] if len(np.array(sim_protocols).shape) == 0 else sim_protocols

        # read in experimental data
        self.raw_data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")

        # determine how many experiments there are and what the time points and observables are for each
        self.experiments = list(dict.fromkeys([d['expt_id'] for d in self.raw_data])) # preserves order of appearance
        self.n_experiments = len(self.experiments)
        self.observables = [np.unique([d['observable'] for d in self.raw_data if d['expt_id'] == expt])
                            for expt in self.experiments]
        # all time points are extracted for all observables in each experiment
        self.tdata = [list(np.unique([d['time'] for d in self.raw_data if d['expt_id'] == expt]))
                      for expt in self.experiments]
        # If the 'alt_expt_id' column doesn't exist, create a new array and add it with values equal to 'expt_id'
        if 'alt_expt_id' not in self.raw_data.dtype.names:
            expt_id_dtype = self.raw_data.dtype['expt_id']
            new_dtype = np.dtype(self.raw_data.dtype.descr + [('alt_expt_id', expt_id_dtype)])
            new_data = np.zeros(self.raw_data.shape, dtype=new_dtype)
            for col in self.raw_data.dtype.names:
                new_data[col] = self.raw_data[col]
            new_data['alt_expt_id'] = self.raw_data['expt_id']
            self.raw_data = new_data
        # get alternative experiment IDs for all experiments, retaining order of appearance
        # (allows multiple experiments to be merged under one ID)
        alt_experiments = dict(
            zip(
                self.experiments,
                [list(dict.fromkeys([d['alt_expt_id'] for d in self.raw_data if d['expt_id'] == expt]))
                 for expt in self.experiments]
            )
        )

        # for a given experiment, there may be missing data points, so create masks that indicate for each observable in
        # each experiment which time points we have data for
        self.tspan_masks = []
        for expt, obs, tspan in zip(self.experiments, self.observables, self.tdata):
            self.tspan_masks.append({})
            for name in obs:
                # create 2D masks for each observable for each experiment based on the alt_expt IDs
                self.tspan_masks[-1][name] = [[False for t in tspan] for x in alt_experiments[expt]]
                tdata = [[d['time'] for d in self.raw_data if d['observable'] == name and d['expt_id'] == expt and
                         d['alt_expt_id'] == alt_expt] for alt_expt in alt_experiments[expt]]
                for i in range(len(tspan)):
                    for j in range(len(tdata)):
                        if tspan[i] in tdata[j]:
                            self.tspan_masks[-1][name][j][i] = True

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
        for n, (expt, obs, tspan) in enumerate(zip(self.experiments, self.observables, self.tdata)):
            self.like_data.append({})
            for name in obs:
                data = [d for d in self.raw_data if d['observable'] == name and d['expt_id'] == expt]
                tdata = [d['time'] for d in data]
                n_data_pts = [len([m for m in mask if m is True]) for mask in self.tspan_masks[n][name]]
                avg = []
                se = []
                start = 0
                for n_pts in n_data_pts:
                    for t in tspan:  # add data in same order as the time points in tspan
                        if t in tdata[start: start + n_pts]:
                            avg.append([d['average'] for d in data[start: start + n_pts] if d['time'] == t][0])
                            se.append([d['stderr'] for d in data[start: start + n_pts] if d['time'] == t][0])
                    start += n_pts
                self.like_data[-1][name] = norm(loc=avg, scale=se)
                # now flatten the tspan_mask, so nothing will change for the simple case of one alt_expt_id per expt_id
                self.tspan_masks[n][name] = np.array(self.tspan_masks[n][name]).flatten()

        # store the full set of parameter values
        self.param_values = np.array([p.value for p in self.model.parameters])

        # create the array for storing the likelihood values
        # NOTE: creating the array here and reusing it, rather than creating a new one every time the likelihood
        # function is accessed, should save some time
        self.logp_data = [0.] * self.n_experiments

    # likelihood function
    def likelihood(self, position):
        y = np.copy(position)
        for n in range(self.n_experiments):
            self.logp_data[n] = 0.  # reinitialize this to zero
            self.param_values[self.parameter_idxs] = 10 ** y[self.samples_idxs[n]]
            # run simulation
            output = self.sim_protocols[n].run(self.tdata[n], self.param_values)
            # calculate log-likelihood
            for obs in self.like_data[n].keys():
                if np.any(np.isnan(output[obs])):  # return -inf if simulation returns NaNs
                    return -np.inf
                self.logp_data[n] += np.sum(self.like_data[n][obs].logpdf(output[obs][self.tspan_masks[n][obs]]))
            if np.isnan(self.logp_data[n]):  # return -inf if logp is NaN
                return -np.inf
        return sum(self.logp_data)

    def run(self, nchains=3, niterations=50000, start=None, restart=False, verbose=True, adapt_gamma=True,
            history_thin=1, gamma_levels=4, multitry=False, obs_labels=None, plot_results=True, plot_ll_args=None,
            plot_pd_args=None, plot_tc_args=None):

        sampled_params, log_ps = run_dream(parameters=self.sampled_params_list,
                                           likelihood=self.likelihood,
                                           nchains=nchains,
                                           niterations=niterations,
                                           start=start,
                                           restart=restart,
                                           verbose=verbose,
                                           adapt_gamma=adapt_gamma,
                                           history_thin=history_thin,
                                           gamma_levels=gamma_levels,
                                           multitry=multitry,
                                           model_name='dreamzs_%dchain' % nchains)
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
                                                   nchains=nchains,
                                                   niterations=niterations,
                                                   start=starts,
                                                   restart=True,
                                                   verbose=verbose,
                                                   adapt_gamma=adapt_gamma,
                                                   history_thin=history_thin,
                                                   gamma_levels=gamma_levels,
                                                   multitry=multitry,
                                                   model_name='dreamzs_%dchain' % nchains)
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
                       plot_pd_args=None, plot_tc_args=None, which_plots='all', max_time_points=1000):

        # which plots should we create? (NOTE: user can pass an integer too)
        which_plots = \
            3 if which_plots == 'all' or which_plots == 'll_pd_tc' else \
            2 if which_plots == 'll_pd' else \
            1 if which_plots == 'll' else \
            0 if which_plots == 'none' else \
            which_plots

        # error check
        if not isinstance(which_plots, int):
            raise Exception("'which_plots' must be an integer or one of the following strings: 'all', " +
                            "'ll_pd_tc', 'll_pd', 'll', 'none'")

        # process plotting function arguments
        _plot_ll_args = {'cutoff': 2}
        _plot_pd_args = {'labels': None, 'groups': None, 'save_plot': None, 'sharex': 'all', 'sharey': 'none'}
        _plot_tc_args = {'tspans': None, 'xlabel': None, 'ylabels': None, 'leg_labels': None, 'separate_plots': True}
        if plot_ll_args is not None:
            _plot_ll_args.update(plot_ll_args)
        if plot_pd_args is not None:
            _plot_pd_args.update(plot_pd_args)
        if plot_tc_args is not None:
            _plot_tc_args.update(plot_tc_args)

        # List of objects returned from each function
        # NOTE: only `plot_param_dist()` returns anything right now. The other functions return `None`
        return_objects = [None] * 3

        if which_plots > 0:
            print('Plotting log-likelihoods')
            return_objects[0] = self.plot_log_likelihood(logps_files, **_plot_ll_args)

        if which_plots > 1:
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
                    param_diff = sorted(list(set.intersection(*param_list_1) - set.union(*param_list_2)))
                    if len(param_diff) > 0:
                        param_groups.append(param_diff)
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
                            break
            # make the plots
            if _plot_pd_args['labels'] is None:
                _plot_pd_args['labels'] = param_labels
            if _plot_pd_args['groups'] is None:
                _plot_pd_args['groups'] = param_groups
            _plot_pd_args['cutoff'] = _plot_ll_args['cutoff']
            if _plot_pd_args['save_plot'] is None:
                _plot_pd_args['save_plot'] = filenames
            return_objects[1] = param_samples = self.plot_param_dist(samples_files, **_plot_pd_args)

        if which_plots > 2:
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
            # create observable legend labels
            leg_labels = [[obs_name for obs_name in self.observables[n]] for n in range(self.n_experiments)] \
                if obs_labels is None else [[obs_labels.get(obs_name, obs_name) for obs_name in self.observables[n]]
                                            for n in range(self.n_experiments)]
            # get output time points for simulations
            tspans = _plot_tc_args.pop('tspans')  # pop tspans out of the dict since it's not passed as a kwarg below
            # if tspans not defined by user, increase the number of time points for the simulations by a factor of 10
            if tspans is None:
                tspans = [np.linspace(
                    self.tdata[n][0], self.tdata[n][-1],
                    min(int((self.tdata[n][-1] - self.tdata[n][0]) * 10 + 1), max_time_points))
                    for n in range(self.n_experiments)]
                # add in the experimental time points and filter them out if they overlap with the points above
                for n in range(self.n_experiments):
                    tspans[n] = np.unique(np.append(tspans[n], self.tdata[n][1:-1]))
            # make the plots
            if _plot_tc_args['xlabel'] is None:
                _plot_tc_args['xlabel'] = xlabel
            if _plot_tc_args['ylabels'] is None:
                _plot_tc_args['ylabels'] = ylabels
            if _plot_tc_args['leg_labels'] is None:
                _plot_tc_args['leg_labels'] = leg_labels
            # pop separate_plots out of the dict since it's not passed as a kwarg below
            separate_plots = _plot_tc_args.pop('separate_plots')
            return_objects[2] = self.plot_timecourses(
                self.model, tspans, self.sim_protocols, param_samples, self.parameter_idxs, exp_data=self.raw_data,
                samples_idxs=self.samples_idxs, separate_plots=separate_plots, **_plot_tc_args)

        if show_plots and which_plots > 0:
            plt.show()

        return tuple(return_objects)


    @staticmethod
    def plot_log_likelihood(logps_files, cutoff=None, show_plot=False, save_plot=True):

        # get the path from the first file
        path, file = os.path.split(logps_files[0])
        m = re.search(r'(.+)_chain', file)
        prefix = m.group(1)

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
                os.path.join(path, '%s_chain_%d_%d.npy' % (prefix, chain, it))).flatten() for it in iterations)))
        log_ps = np.array(log_ps)

        # plot the likelihoods
        plt.figure(constrained_layout=True)
        # calculate mean and variance for last half of steps
        burnin = int(max(iterations) / 2)
        log_ps_max = -np.inf
        log_ps_mean = 0
        log_ps_var = 0
        for i, chain in enumerate(chains):
            plt.plot(range(len(log_ps[i])), log_ps[i], label='chain %d' % chain)
            log_ps_max = np.max(log_ps[i]) if log_ps_max < np.max(log_ps[i]) else log_ps_max
            log_ps_mean += np.mean(log_ps[i][burnin:]) / len(chains)
            log_ps_var += np.var(log_ps[i][burnin:]) / len(chains)  # mean of the variances, but that's fine
        # plot the mean over all chains
        steps = list(np.arange(0, len(log_ps[0]), len(log_ps[0]) // 100))
        plt.plot(steps, np.mean(log_ps, axis=0)[steps], 'k', lw=2, label='average')
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
        if save_plot is not False:
            filename = 'fig_PyDREAM_log_ps' if save_plot is True else save_plot
            plt.savefig(filename)

        return None  # this is a placeholder, in case we want to return something in the future

    @staticmethod
    def plot_param_dist(sample_files, labels=None, groups=None, cutoff=None, show_plot=False, save_plot=True, **kwargs):


        # Helper function for getting the optimal number of columns for the histogram figure
        def get_ncols(ndims):
            if not isinstance(ndims, int) or ndims <= 0:
                raise ValueError("'ndims' must be a positive integer")
            r1 = round(math.sqrt(ndims))  # round() returns an int
            r2 = math.ceil(math.sqrt(ndims))  # math.ceil() also returns an int
            while r1 * r2 >= ndims:
                r1 -= 1
                r2 += 1
            return min(r1 + 1, r2 - 1)  # the smaller of the two integers is the # of columns


        # get chains and iterations from file names
        chains = []
        iterations = []
        for file in sample_files:
            m = re.search(r'chain_(\d+)_(\d+).npy$', file)
            chains.append(int(m.group(1)))
            iterations.append(int(m.group(2)))
        chains = np.unique(chains)
        iterations = np.unique(iterations)
        burnin = int(max(iterations) / 2)  # discard first 50% of samples

        # read in parameter values from files
        samples_chain = []
        for chain in chains:
            # get files and sort them numerically by number of iterations (using key=lambda function)
            files = sorted([file for file in sample_files if re.search(r'chain_%d' % chain, file)],
                           key=lambda f: int(re.search(r'(\d+).npy$', f).group(1)))
            samples_chain.append(np.concatenate(tuple(np.load(file) for file in files)))

        # if a likelihood cutoff is defined, load the likelihoods and remove samples that fall below the cutoff
        if cutoff is not None:
            log_ps_chain = []
            for chain in chains:
                # get files and sort them numerically by number of iterations (using key=lambda function)
                files = sorted([file.replace('sampled_params', 'logps')
                                for file in sample_files if re.search(r'chain_%d' % chain, file)],
                               key=lambda f: int(re.search(r'(\d+).npy$', f).group(1)))
                log_ps_chain.append(np.concatenate(tuple(np.load(file) for file in files)).flatten())
            log_ps = np.concatenate(tuple(log_ps_chain[chain][burnin:] for chain in range(len(chains))))
            log_ps_mean = np.mean(log_ps)
            log_ps_sdev = np.std(log_ps)
            # get indices for samples with log-likelihood > cutoff
            for chain in range(len(chains)):
                keep_idxs = [i for i in range(burnin, len(samples_chain[chain]))
                             if log_ps_chain[chain][i] > log_ps_mean - cutoff * log_ps_sdev]
                samples_chain[chain] = samples_chain[chain][keep_idxs]

        # combine all samples together
        samples = np.concatenate(tuple(samples_chain[chain] for chain in range(len(chains))))

        # plot histograms
        print("Number of samples: %d (of %d total)" % (len(samples), burnin * len(chains)))
        if groups is None:
            groups = [[i for i in range(len(samples[0]))]]
        if labels is None:
            labels = [['p_%d_%d' % (i, j) for j in range(len(groups[i]))] for i in range(len(groups))]
        fig_axs_list = []
        for n, label, group in zip(range(len(labels)), labels, groups):
            ndims = len(group)  # number of dimensions (i.e., parameters)
            # set plot parameters
            ncols = kwargs.get('ncols', get_ncols(ndims))
            nrows = math.ceil(ndims/ncols)
            figsize = kwargs.get('figsize', (0.65 * ncols * 6.4, 0.5 * nrows * 4.8))
            labelsize = kwargs.get('labelsize', 10 * max(1, (2/5 * np.ceil(nrows / 2))))
            fontsize = kwargs.get('fontsize', 10 * max(1, (3/5 * np.ceil(nrows / 2))))
            sharex = kwargs.get('sharex', 'all')
            sharey = kwargs.get('sharey', 'none')
            bw_adjust = kwargs.get('bw_adjust', 1)
            # create figure
            colors = sns.color_palette(n_colors=ndims)
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey, constrained_layout=True,
                                    figsize=figsize)
            if nrows == 1 and ncols == 1:
                axs = np.array([axs])  # make it an np.ndarray
            # save the fig and axes for later
            fig_axs_list.append((fig, axs.flatten()))
            ##### TODO
            # fig2, axs2 = plt.subplots(nrows=nrows, ncols=ncols, sharex='none', sharey=sharey, constrained_layout=True,
            #                           figsize=figsize)
            #####
            row = 0
            col = 0
            for dim in range(ndims):
                print(group[dim], end=' ')
                sns.kdeplot(samples[:, group[dim]], color=colors[dim], fill=True, common_norm=False, ax=axs[row][col],
                            bw_adjust=bw_adjust)
                axs[row][col].set_yticklabels([])
                axs[row][col].set_ylabel(None)
                axs[row][col].set_title(label[dim], fontsize=labelsize)
                axs[row][col].tick_params(axis='x', labelsize=labelsize)
                ##### TODO
                # for chain in range(len(chains)):
                #     sns.kdeplot(samples_chain[chain][:, group[dim]], fill=None, common_norm=False, ax=axs2[row][col])
                #     #, label='chain %d' % chains[chain])
                # # axs2[row][col].legend(loc=0)
                # axs2[row][col].set_yticklabels([])
                # axs2[row][col].set_ylabel(None)
                # axs2[row][col].set_title(label[dim], fontsize=labelsize)
                # axs2[row][col].tick_params(axis='x', labelsize=labelsize)
                #####
                col += 1
                if col % ncols == 0:
                    col = 0
                    row += 1
            print()
            fig.supxlabel(r'log$_{10}$ value', fontsize=fontsize)
            fig.supylabel('Density', fontsize=fontsize)
            ##### TODO
            # fig2.supxlabel(r'log$_{10}$ value', fontsize=fontsize)
            # fig2.supylabel('Density', fontsize=fontsize)
            #####
            # delete extra plots
            if col > 0:
                while col < ncols:
                    fig.delaxes(axs[row][col])
                    ##### TODO
                    # fig2.delaxes(axs[row][col])
                    #####
                    col += 1

        # Determine common x-limits for all figures, if requested (only applicable if there are multiple param groups)
        if len(fig_axs_list) > 1 and kwargs.get('sharex', 'all') == 'all':
            # Extract x-limits across all subplots
            x_min = np.inf
            x_max = -np.inf
            for fig, axs in fig_axs_list:
                for ax in axs:
                    x_min = min(x_min, ax.get_xlim()[0])
                    x_max = max(x_max, ax.get_xlim()[1])
            # Apply the common x-limits to all figures
            for fig, axs in fig_axs_list:
                for ax in axs:
                    ax.set_xlim(x_min, x_max)

        # save plots
        if save_plot is not False:
            for n, (fig, axs) in enumerate(fig_axs_list):
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
                fig.savefig(filename)

        if show_plot:
            plt.show()

        return samples  # these samples have been trimmed based on log-likelihood values


    @staticmethod
    def plot_timecourses(model, tspans, sim_protocols, param_samples, parameter_idxs, observables=None, exp_data=None,
                         samples_idxs=None, show_plot=False, save_plot=True, separate_plots=True, timeout_seconds=None,
                         **kwargs):

        # if there's only one simulation to run, put the following objects inside a list of size 1, since the code
        # is set up to loop over the number of simulations
        # NOTE: checking first elements of 'tspans', 'observables', and 'samples_idxs' because they could be jagged, in
        # which case calling 'np.array' would raise an error
        if len(np.array(tspans[0]).shape) == 0:
            tspans = [tspans]
        if len(np.array(sim_protocols).shape) == 0:
            sim_protocols = [sim_protocols]
        if observables is not None and len(np.array(observables[0]).shape) == 0:
            observables = [observables]
        if samples_idxs is not None and len(np.array(samples_idxs[0]).shape) == 0:
            samples_idxs = [samples_idxs]

        # experimental data
        n_experiments = 0
        if exp_data is not None:
            raw_data = exp_data if isinstance(exp_data, np.ndarray) \
                else np.genfromtxt(exp_data, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
            experiments = list(dict.fromkeys([d['expt_id'] for d in raw_data]))  # preserves order of appearance
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
                "'samples_idxs")

        # number of simulations
        n_sims = len(tspans)

        # flatten the observables matrix and pull out the unique names so we can loop over them later
        obs_names = np.unique(np.concatenate([obs for obs in observables]))

        # process kwargs
        fill_between = kwargs.get('fill_between', (5, 95))
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # standard colors
        colors = kwargs.get('colors',
                            # each observable has a different color
                            np.array([[cycle[i % 10] for j in range(n_sims)] for i in range(len(obs_names))]).T
                            if separate_plots else
                            # each experiment/simulation has a different color
                            np.array([[cycle[j % 10] for j in range(n_sims)] for name in obs_names]).T)

        locs = kwargs.get('locs', [[0] * len(observables[n]) for n in range(n_sims)])
        xlabel = kwargs.get('xlabel', 'time')
        ylabels = kwargs.get('ylabels', [['amount'] * len(observables[n]) for n in range(n_sims)])
        leg_labels = kwargs.get('leg_labels', [[obs_name for obs_name in observables[n]] for n in range(n_sims)])
        save_sim_data = kwargs.get('save_sim_data', False)
        # make sure arrays are 2D
        # NOTE: checking first elements of the 'colors', 'locs', and 'ylabels' arrays because they could be jagged, in
        # which case calling 'np.array' would raise an error
        if len(np.array(colors[0]).shape) == 0:
            colors = [colors]
        if len(np.array(locs[0]).shape) == 0:
            locs = [locs]
        if len(np.array(ylabels[0]).shape) == 0:
            ylabels = [ylabels]

        # store original parameter values
        param_values = np.array([p.value for p in model.parameters])

        # run simulations using only unique parameter samples
        samples_unique, counts = np.unique(param_samples, return_counts=True, axis=0)
        # FOR DEBUGGING ### TODO
        # samples_unique = samples_unique[:10]
        # counts = counts[:10]
        # #################
        # save simulation data, if requested
        if save_sim_data:
            csvfile = open("SIM_DATA.csv", 'w')
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(['observable', 'time', 'yval_min', 'yval_max', 'sim_id'])
        # loop over simulations (experiments + perturbations)
        print('Running %d simulations (of %d samples)' % (len(samples_unique), len(param_samples)))
        for n in range(n_sims):
            if n < n_experiments:
                print("Experiment '%s' (%d of %d)..." % (str(experiments[n]), n+1, n_experiments))
            else:
                print("Perturbation %d of %d..." % (n - n_experiments + 1, n_sims - n_experiments))
            outputs = []
            counts_n = counts.copy()  # don't change 'counts' array, make a copy
            # run simulations
            for i, sample in enumerate(samples_unique):
                print(i, end=' ')
                param_values[parameter_idxs] = 10 ** sample[samples_idxs[n]]
                if timeout_seconds is None:
                    outputs.append(sim_protocols[n].run(tspans[n], param_values))
                else:  # use timeout option (if PyDREAM ran as expected, this shouldn't be necessary)
                    output = None
                    with multiprocessing.Pool(processes=1) as pool:
                        future_output = pool.apply_async(run_simulation, (sim_protocols[n], tspans[n], param_values,))
                        try:
                            output = future_output.get(timeout=timeout_seconds)  # Wait for completion
                        except multiprocessing.TimeoutError:
                            print(f"Skipping parameter sample {i} due to timeout.")
                            pool.terminate()  # Kill the process if it exceeds the timeout
                            counts_n = np.delete(counts_n, i, axis=0)  # remove counts for this parameter set
                    if output is not None:
                        outputs.append(output)
                if (i + 1) % 20 == 0:
                    print()
            print('DONE')
            # remove any simulations that produced NaNs or ended prematurely
            idx_remove = [i for i in range(len(outputs)) if np.any(np.isnan(outputs[i][observables[n][0]]))
                          or len(outputs[i]) != len(tspans[n])]
            if len(idx_remove) > 0:
                print('Removing simulations', idx_remove, 'that produced NaNs')
            outputs = np.delete(outputs, idx_remove, axis=0)
            counts_n = np.delete(counts_n, idx_remove, axis=0)
            # use 'counts' to generate full set of simulation outputs for correct weighting of plots
            outputs = np.repeat(outputs, counts_n, axis=0)
            # plot results
            for i, obs_name in enumerate(observables[n]):
                print(obs_name)
                figname = '%s_exp_%s' % (obs_name, str(experiments[n])) if separate_plots and n < n_experiments \
                    else '%s_sim_%s' % (obs_name, n - n_experiments + 1) if separate_plots and n >= n_experiments \
                    else obs_name
                plt.figure(num=figname, constrained_layout=True)
                # plot simulated data as a percent envelope
                print('   Simulated data...', end='')
                yvals = np.array([output[obs_name] for output in outputs])
                # reshape yvals_min and yvals_max to handle the case where multiple experiments are merged into one
                # using the ParallelExperiments simulation protocol class
                yvals_min = np.percentile(yvals, fill_between[0], axis=0).reshape(-1, len(tspans[n]))
                yvals_max = np.percentile(yvals, fill_between[1], axis=0).reshape(-1, len(tspans[n]))
                # loop over each sub-experiment (if there's only 1 expt, will only loop once)
                hatch_patterns = [None, '/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
                for j, (ymin, ymax) in enumerate(zip(yvals_min, yvals_max)):
                    hatch = hatch_patterns[j % len(hatch_patterns)]
                    plt.fill_between(tspans[n], ymin, ymax, alpha=0.25, color=colors[n][i], hatch=hatch,
                                     hatch_linewidth=3, label=leg_labels[n][i])
                    # save simulation data, if requested
                    if save_sim_data:
                        sim_id = experiments[n] if n < n_experiments else n
                        for line in zip(tspans[n], ymin, ymax):
                            csvwriter.writerow([obs_name] + list(line) + [sim_id])
                print('DONE')
                # plot experimental data
                if n < n_experiments:
                    label = 'experiment' if n_experiments == 1 else 'experiment %s' % str(experiments[n])
                    suffix = ['']
                    if len(yvals_min) == 1:
                        times = [[d['time'] for d in raw_data if d['observable'] == obs_name
                                and d['expt_id'] == experiments[n]]]
                        avgs = [[d['average'] for d in raw_data if d['observable'] == obs_name
                               and d['expt_id'] == experiments[n]]]
                        stderrs = [[d['stderr'] for d in raw_data if d['observable'] == obs_name
                                  and d['expt_id'] == experiments[n]]]
                    else: # need to use if/else because 'alt_expt_id' column may not exist if len(yvals_min) == 1
                        # get alt_expt_ids, retaining order of appearance
                        alt_expt_ids = list(dict.fromkeys(
                            [d['alt_expt_id'] for d in raw_data if d['expt_id'] == experiments[n]]))
                        suffix = [' (alt %d)' % (id + 1) for id in range(len(alt_expt_ids))]
                        times = [[d['time'] for d in raw_data if d['observable'] == obs_name
                                and d['expt_id'] == experiments[n] and d['alt_expt_id'] == alt_expt_id]
                                for alt_expt_id in alt_expt_ids]
                        avgs = [[d['average'] for d in raw_data if d['observable'] == obs_name
                               and d['expt_id'] == experiments[n] and d['alt_expt_id'] == alt_expt_id]
                               for alt_expt_id in alt_expt_ids]
                        stderrs = [[d['stderr'] for d in raw_data if d['observable'] == obs_name
                                  and d['expt_id'] == experiments[n] and d['alt_expt_id'] == alt_expt_id]
                                  for alt_expt_id in alt_expt_ids]
                    markers = ['o', '^', 's', 'v', '<', '>', 'D', 'p', 'H', '*']
                    for j, (time, avg, stderr) in enumerate(zip(times, avgs, stderrs)):
                        marker = markers[j % len(markers)]
                        plt.errorbar(time, avg, yerr=stderr, capsize=6, fmt=marker, ms=8, mfc=colors[n][i],
                                     mec=colors[n][i], ecolor=colors[n][i], label=label+suffix[j])

                    # To avoid confusion, only print the following message if there is actual expt data
                    # NOTE: to get the correct legend labels, we still need to make the errorbar plot above, even if
                    # there is no expt data
                    if len(times[0]) > 0:
                        print('   Experimental data...DONE')
                plt.xlabel(xlabel)
                plt.ylabel(ylabels[n][i])
                plt.legend(loc=locs[n][i])

                if save_plot is not None and separate_plots:
                    filename = 'fig_PyDREAM_tc_%s' % figname if save_plot is True else save_plot
                    plt.savefig(filename)

        if save_sim_data:
            csvfile.close()

        if save_plot is not None and not separate_plots:
            # 'obs_names' is an array of unique observables (defined above)
            for obs_name in obs_names:  # each plot is an observable, so this is actually a loop over all plots
                # need to loop over the experiments in order to get the correct legend location
                for n in range(n_sims):
                    # check if the observable is in the list for this experiment. if it is, reorder the legend, save the
                    # file, and break out of the loop
                    if obs_name in observables[n]:
                        filename = 'fig_PyDREAM_tc_%s' % obs_name if save_plot is True else save_plot
                        # fix legend order
                        plt.figure(num=obs_name)
                        handles, labels = plt.gca().get_legend_handles_labels()
                        N = len([label for label in labels if 'experiment' in label])  # No. of expts in this plot
                        M = len(labels) - N  # total number of sims in this plot

                        new_handles = [tuple([handles[j] for j in range(i, len(handles), M)]) for i in range(N)] + \
                                      [handles[i + N] for i in range(M-N)]
                        new_labels = [': '.join([labels[j] for j in range(i, len(labels), M)]) for i in range(N)] + \
                                     [labels[i + N] for i in range(M-N)]

                        # if there are no expt data points, remove the second handle from the 'new_handles' tuples
                        for i in range(N):
                            if len(new_handles[i][1].lines[0].get_data()[0]) == 0:
                                new_handles[i] = new_handles[i][0]

                        plt.legend(new_handles, new_labels, loc=locs[n][list(observables[n]).index(obs_name)])
                        plt.savefig(filename)
                        break

        if show_plot:
            plt.show()

        return None  # this is a placeholder, in case we want to return something in the future
