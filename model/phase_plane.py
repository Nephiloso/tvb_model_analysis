# -*- coding: utf-8 -*-
#
#
#  Modified from TheVirtualBrain-Scientific Package.
#
# (c) 2012-2022, Baycrest Centre for Geriatric Care ("Baycrest") and others
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#   CITATION:
# When using The Virtual Brain for scientific publications, please cite it as follows:
#
#   Paula Sanz Leon, Stuart A. Knock, M. Marmaduke Woodman, Lia Domide,
#   Jochen Mersmann, Anthony R. McIntosh, Viktor Jirsa (2013)
#       The Virtual Brain: a simulator of primate brain network dynamics.
#   Frontiers in Neuroinformatics (7:10. doi: 10.3389/fninf.2013.00010)
#
#

"""
An phase-plane plot generated from a Model object of TVB.

Optionally an Integrator object from TVB can be specified, this will be used to
generate sample trajectories -- not the phase-plane. This is mainly interesting
for visualising the effect of noise on a trajectory.

Demo::

    from phase_plane import PhasePlaneInteractive
    ppi_fig = PhasePlaneInteractive()
    ppi_fig.show()


Example specifying a Model and stochastic sample trajectories::

    import tvb.simulator.models
    from tvb.simulator.plot.phase_plane_interactive import PhasePlaneInteractive
    MODEL = tvb.simulator.models.JansenRit()
    import tvb.simulator.integrators
    INTEGRATOR = tvb.simulator.integrators.HeunStochastic(dt=2**-5)
    ppi_fig = PhasePlaneInteractive(model=MODEL, integrator=INTEGRATOR)
    ppi_fig.show()


.. moduleauthor:: Stuart A. Knock <Stuart@tvb.invalid>

"""

import numpy
import pylab
import colorsys
# import matplotlib.widgets as widgets

from tvb.simulator.common import get_logger
import tvb.simulator.models as models_module
import tvb.simulator.integrators as integrators_module
from tvb.basic.neotraits.api import HasTraits, Attr, NArray, List

LOG = get_logger(__name__)

# Define a colour theme... see: matplotlib.colors.cnames.keys()
BACKGROUNDCOLOUR = "lightgray"
EDGECOLOUR = "darkslateblue"
AXCOLOUR = "steelblue"

# Set the resolution of the phase-plane and sample trajectories.
NUMBEROFGRIDPOINTS = 42


def get_color(num_colours):
    for hue in range(num_colours):
        hue = 1.0 * hue / num_colours
        col = [int(x) for x in colorsys.hsv_to_rgb(hue, 1.0, 230)]
        yield "#{0:02x}{1:02x}{2:02x}".format(*col)


class PhasePlaneInteractive(HasTraits):
    """
    The GUI for the interactive phase-plane viewer provides sliders for setting:
        - The value of all parameters of the Model.
        - The extent of the axes.
        - A fixed value for the state-variables which aren't currently selected.
        - The noise strength, if a stocahstic integrator is specified.

    and radio buttons for selecting:
        - Which state-variables to show on each axis.
        - Which mode to show, if the Model has them.

    Clicking on the phase-plane will generate a sample trajectory, originating
    from where you clicked.

    """

    model = Attr(
        field_type=models_module.Model,
        label="Model",
        default=models_module.Generic2dOscillator(),
        doc="""An instance of the local dynamic model to be investigated with
        PhasePlaneInteractive.""")

    integrator = Attr(
        field_type=integrators_module.Integrator,
        label="Integrator",
        default=integrators_module.RungeKutta4thOrderDeterministic(),
        doc="""The integration scheme used to for generating sample
        trajectories on the phase-plane. NOTE: This is not used for generating
        the phase-plane itself, ie the vector field and nulclines.""")

    exclude_sliders = List(
        of=str)

    def __init__(self, traj_steps = 4096, **kwargs):
        """
        Initialise based on provided keywords or their traited defaults. Also,
        initialise the place-holder attributes that aren't filled until the
        show() method is called.
        """
        super(PhasePlaneInteractive, self).__init__(**kwargs)
        LOG.debug(str(kwargs))

        # figure
        self.ipp_fig = None
        self.traj_steps = traj_steps

        # p hase-plane
        self.pp_ax = None
        self.X = None
        self.Y = None
        self.U = None
        self.V = None
        self.UVmag = None
        self.nullcline_x = None
        self.nullcline_y = None
        self.pp_quivers = None

        # Current state
        self.svx = None
        self.svy = None
        self.default_sv = None
        self.no_coupling = None
        self.mode = None

        self.initial_states = [[0.2, 0.4], [0.4, 0.85], [0.85, 0.2]]
        self.traj = None

    def show(self):
        """ Generate the interactive phase-plane figure. """
        model_name = self.model.__class__.__name__
        msg = "Generating an interactive phase-plane plot for %s"
        LOG.info(msg % model_name)

        # Make sure the model is fully configured...
        self.model.configure()

        # Setup the inital(current) state
        self.svx = self.model.state_variables[0]  # x-axis: 1st state variable
        self.svy = self.model.state_variables[1]  # y-axis: 2nd state variable
        self.mode = 0
        self.set_state_vector()

        # Make the figure:
        self.create_figure()
        if isinstance(self.integrator, integrators_module.IntegratorStochastic):
            if self.integrator.noise.ntau > 0.0:
                self.integrator.noise.configure_coloured(self.integrator.dt,
                                                         (1, self.model.nvar, 1,
                                                          self.model.number_of_modes))
            else:
                self.integrator.noise.configure_white(self.integrator.dt,
                                                      (1, self.model.nvar, 1,
                                                       self.model.number_of_modes))

        # Calculate the phase plane
        self.set_mesh_grid()
        self.calc_phase_plane()

        # Plot phase plane
        self.plot_phase_plane()

        for i in range(len(self.initial_states)):
            if i==0:
                self.traj = self.plot_trajectory(self.initial_states[i][0], self.initial_states[i][1])
            else:
                self.plot_trajectory(self.initial_states[i][0], self.initial_states[i][1], False)

        pylab.show()

    ##------------------------------------------------------------------------##
    ##----------------- Functions for building the figure --------------------##
    ##------------------------------------------------------------------------##

    def create_figure(self):
        """ Create the figure and phase-plane axes. """
        # Figure and main phase-plane axes
        model_name = self.model.__class__.__name__
        integrator_name = self.integrator.__class__.__name__
        figsize = 10, 5
        try:
            figure_window_title = "Phase-plane: " + model_name
            figure_window_title += "   --   %s" % integrator_name
            self.ipp_fig = pylab.figure(num=figure_window_title,
                                        figsize=figsize,
                                        facecolor=BACKGROUNDCOLOUR,
                                        edgecolor=EDGECOLOUR)
        except ValueError:
            LOG.info("My life would be easier if you'd update your PyLab...")
            self.ipp_fig = pylab.figure(num=42, figsize=figsize,
                                        facecolor=BACKGROUNDCOLOUR,
                                        edgecolor=EDGECOLOUR)

        self.pp_ax = self.ipp_fig.add_axes([0.265, 0.2, 0.5, 0.75])

        self.pp_splt = self.ipp_fig.add_subplot(212)
        self.ipp_fig.subplots_adjust(left=0.265, bottom=0.02, right=0.765,
                                     top=0.3, wspace=0.1, hspace=None)
        self.pp_splt.set_prop_cycle(color=get_color(self.model.nvar))
        self.pp_splt.plot(numpy.arange(self.traj_steps + 1) * self.integrator.dt,
                          numpy.zeros((self.traj_steps + 1, self.model.nvar)))
        if hasattr(self.pp_splt, 'autoscale'):
            self.pp_splt.autoscale(enable=True, axis='y', tight=True)
        self.pp_splt.legend(self.model.state_variables)
        
        # self.pp_sprm = self.ipp_fig.add_subplot(222)


    def set_mesh_grid(self):
        """
        Generate the phase-plane gridding based on currently selected 
        state-variables and their range values.
        """
        xlo = self.model.state_variable_range[self.svx][0]
        xhi = self.model.state_variable_range[self.svx][1]
        ylo = self.model.state_variable_range[self.svy][0]
        yhi = self.model.state_variable_range[self.svy][1]

        self.X = numpy.mgrid[xlo:xhi:(NUMBEROFGRIDPOINTS * 1j)]
        self.Y = numpy.mgrid[ylo:yhi:(NUMBEROFGRIDPOINTS * 1j)]

    def set_state_vector(self):
        """
        Set up a vector containing the default state-variable values and create
        a filler(all zeros) for the coupling arg of the Model's dfun method.
        This method is called once at initialisation (show()).
        """
        # import pdb; pdb.set_trace()
        sv_mean = numpy.array([self.model.state_variable_range[key].mean() for key in self.model.state_variables])
        sv_mean = sv_mean.reshape((self.model.nvar, 1, 1))
        self.default_sv = sv_mean.repeat(self.model.number_of_modes, axis=2)
        self.no_coupling = numpy.zeros((self.model.nvar, 1,
                                        self.model.number_of_modes))

    def calc_phase_plane(self):
        """ Calculate the vector field. """
        svx_ind = self.model.state_variables.index(self.svx)
        svy_ind = self.model.state_variables.index(self.svy)

        # Calculate the vector field discretely sampled at a grid of points
        grid_point = self.default_sv.copy()
        self.U = numpy.zeros((NUMBEROFGRIDPOINTS, NUMBEROFGRIDPOINTS,
                              self.model.number_of_modes))
        self.V = numpy.zeros((NUMBEROFGRIDPOINTS, NUMBEROFGRIDPOINTS,
                              self.model.number_of_modes))
        for ii in range(NUMBEROFGRIDPOINTS):
            grid_point[svy_ind] = self.Y[ii]
            for jj in range(NUMBEROFGRIDPOINTS):
                # import pdb; pdb.set_trace()
                grid_point[svx_ind] = self.X[jj]

                d = self.model.dfun(grid_point, self.no_coupling)

                for kk in range(self.model.number_of_modes):
                    self.U[ii, jj, kk] = d[svx_ind, 0, kk]
                    self.V[ii, jj, kk] = d[svy_ind, 0, kk]

        # Colours for the vector field quivers
        # self.UVmag = numpy.sqrt(self.U**2 + self.V**2)

        # import pdb; pdb.set_trace()
        if numpy.isnan(self.U).any() or numpy.isnan(self.V).any():
            LOG.error("NaN")

    def plot_phase_plane(self):
        """ Plot the vector field and its nullclines. """
        # Set title and axis labels
        model_name = self.model.__class__.__name__
        self.pp_ax.set(title=model_name + " mode " + str(self.mode))
        self.pp_ax.set(xlabel="State Variable " + self.svx)
        self.pp_ax.set(ylabel="State Variable " + self.svy)

        # import pdb; pdb.set_trace()
        # Plot a discrete representation of the vector field
        if numpy.all(self.U[:, :, self.mode] + self.V[:, :, self.mode] == 0):
            self.pp_ax.set(title=model_name + " mode " + str(self.mode) + ": NO MOTION IN THIS PLANE")
            X, Y = numpy.meshgrid(self.X, self.Y)
            self.pp_quivers = self.pp_ax.scatter(X, Y, s=8, marker=".", c="k")
        else:
            self.pp_quivers = self.pp_ax.quiver(self.X, self.Y,
                                                self.U[:, :, self.mode],
                                                self.V[:, :, self.mode],
                                                # self.UVmag[:, :, self.mode],
                                                width=0.001, headwidth=8)

        # Plot the nullclines
        self.nullcline_x = self.pp_ax.contour(self.X, self.Y,
                                              self.U[:, :, self.mode],
                                              [0], colors="r")
        self.nullcline_y = self.pp_ax.contour(self.X, self.Y,
                                              self.V[:, :, self.mode],
                                              [0], colors="g")
        pylab.draw()
# !
    def plot_trajectory(self, x, y, plot_splt=True):
        """
        Plot a sample trajectory, starting at the position x,y in the
        phase-plane. This method is called as a result of a mouse click on the 
        phase-plane.
        """
        svx_ind = self.model.state_variables.index(self.svx)
        svy_ind = self.model.state_variables.index(self.svy)

        # Calculate an example trajectory
        state = self.default_sv.copy()
        self.integrator.clamped_state_variable_indices = numpy.setdiff1d(
            numpy.r_[:len(self.model.state_variables)], numpy.r_[svx_ind, svy_ind])
        # import pdb; pdb.set_trace()
        self.integrator.clamped_state_variable_values = self.default_sv[self.integrator.clamped_state_variable_indices]
        state[svx_ind] = x
        state[svy_ind] = y
        scheme = self.integrator.scheme
        traj = numpy.zeros((self.traj_steps + 1, self.model.nvar, 1,
                            self.model.number_of_modes))
        traj[0, :] = state
        for step in range(self.traj_steps):
            # import pdb; pdb.set_trace()
            state = scheme(state, self.model.dfun, self.no_coupling, 0.0, 0.0)
            traj[step + 1, :] = state

        self.pp_ax.scatter(x, y, s=42, c='g', marker='o', edgecolor=None)
        self.pp_ax.plot(traj[:, svx_ind, 0, self.mode],
                        traj[:, svy_ind, 0, self.mode])

        # Plot the selected state variable trajectories as a function of time
        if plot_splt:
            self.pp_splt.plot(numpy.arange(self.traj_steps + 1) * self.integrator.dt,
                              traj[:, :, 0, self.mode])

        pylab.draw()

        return traj