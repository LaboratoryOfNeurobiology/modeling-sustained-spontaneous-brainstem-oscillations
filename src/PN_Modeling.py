#!/home/zupanc-lab/anaconda3/bin/python
from __future__ import division
from neuron import h, nrn, gui
import matplotlib
from matplotlib import pyplot
import numpy as np
import math
import os
import ast

matplotlib.use('TkAgg')
h.load_file("nrngui.hoc")


###############################################################################

# CLASSES
class Synapse(object):

    def __init__(self, pre_syn, post_syn, conduct_range):
        self.pre = pre_syn
        self.post = post_syn
        self.pre_syn_ident = self.pre.get_ident()
        self.post_syn_ident = self.post.get_ident()
        self.create_halfgaps(conduct_range)

    def create_halfgaps(self, conduct_range):
        self.gaps = [0, 0]
        (pre_soma_dx, pre_axon_dx) = self.pre.give_dxs()
        (post_soma_dx, post_axon_dx) = self.post.give_dxs()
        # Create HalfGap Objects
        self.gaps[0] = h.HalfGap(self.pre.axon(1 - (pre_axon_dx / 2)))
        self.gaps[1] = h.HalfGap(self.post.soma(0 + (post_soma_dx / 2)))
        # Identify cathode vs anode
        self.gaps[0].isanode = 1
        self.gaps[1].isanode = -1
        # Set voltage pointers
        h.setpointer(self.post.soma(0 + (post_soma_dx) / 2)._ref_v,
                     'vgap', self.gaps[0])
        h.setpointer(self.pre.axon(1 - (pre_axon_dx) / 2)._ref_v,
                     'vgap', self.gaps[1])
        # Set conductance min and max
        self.gaps[0].gmin = conduct_range[0]
        self.gaps[1].gmin = conduct_range[0]
        self.gaps[0].gmax = conduct_range[1]
        self.gaps[1].gmax = conduct_range[1]

    def give_gaps(self):
        return self.gaps


class Graph(object):
    """ GRAPH DICT REPRESENTATION
        KEY:
        RANGE(0, NUM_PACEMAKERS)
            represents PACEMAKER cells
        RANGE(NUM_PACEMAKERS, NUM_PACEMAKERS + NUM_RELAYS)
            represents RELAY cells

        Each key value representing a cell has a list [] of cells that it
        connects to in the network.
    """

    def __init__(self, graph_dict=None):
        """ initializes a graph object
            If no dictionary or None is given,
            an empty dictionary will be used
        """
        if graph_dict is None:
            graph_dict = {}
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def edges_of_vertex(self, vertex):
        return self.__graph_dict[vertex]

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary.
            Otherwise nothing has to be done.
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list;
            between two vertices can be multiple edges!
        """
        (vertex1, vertex2) = edge
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def __generate_edges(self):
        """ A static method generating the edges of the
            graph "graph". Edges are represented as sets
            with one (a loop back to the vertex) or two
            vertices
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res


class PacemakerCell(object):
    """Two-section cell: A soma with active channels and
    a axon with active channels."""

    def __init__(self, ps_pars, pa_pars, identifier):
        self.ident = identifier
        self.x = self.y = self.z = 0
        self.connections = []  # List of cell identifiers this cell projects to
        self.synapses = []  # List of Synapse Objects corr. to identifiers
        self.create_sections()
        self.build_subsets()
        self.define_geometry(ps_pars, pa_pars)
        self.define_biophysics(ps_pars, pa_pars)
        self.build_topology()
        self.set_recording_vectors(ps_pars[6], pa_pars[6])

    def create_sections(self):
        """Create the sections of the cell."""
        # NOTE: cell=self is required to tell NEURON of this object.
        self.soma = h.Section(name='soma', cell=self)
        self.axon = h.Section(name='axon', cell=self)

    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        self.axon.connect(self.soma(1))

    def define_geometry(self, ps_pars, pa_pars):
        """Set the 3D geometry of the cell."""
        self.soma.L = self.soma.diam = 30  # microns
        self.axon.L = 45  # microns
        self.axon.diam = 8  # microns
        self.soma.nseg = ps_pars[6]
        self.axon.nseg = pa_pars[6]
        self.soma_dx = (1 / ps_pars[6])
        self.axon_dx = (1 / pa_pars[6])
        self.shape_3D()

    def define_biophysics(self, ps_pars, pa_pars):
        """Assign the membrane properties across the cell."""
        for sec in self.all:  # 'all' defined in build_subsets
            sec.Ra = 100  # Axial resistance in Ohm * cm
            sec.cm = 1  # Membrane capacitance in micro Farads / cm^2

        # Insert active Hodgkin-Huxley current in the soma
        self.soma.insert('hh')
        self.soma.insert('na_ion')
        self.soma.insert('k_ion')
        for seg in self.soma:
            seg.hh.gnabar = ps_pars[3]  # Sodium conductance in S/cm2
            seg.hh.gkbar = ps_pars[4]  # Potassium conductance in S/cm2
            seg.ena = ps_pars[1]
            seg.ek = ps_pars[0]
            seg.hh.gl = ps_pars[5]  # Leak conductance in S/cm2
            seg.hh.el = ps_pars[2]  # Reversal potential in mV

        # Insert active current in the axon
        self.axon.insert('hh')
        self.axon.insert('na_ion')
        self.axon.insert('k_ion')
        for seg in self.axon:
            seg.hh.gnabar = pa_pars[3]  # Sodium conductance in S/cm2
            seg.hh.gkbar = pa_pars[4]  # Potassium conductance in S/cm2
            seg.ena = pa_pars[1]
            seg.ek = pa_pars[0]
            seg.hh.gl = pa_pars[5]  # Leak conductance in S/cm2
            seg.hh.el = pa_pars[2]  # Reversal potential in mV

    def build_subsets(self):
        """Build subset lists. For now we define 'all'."""
        self.all = h.SectionList()
        self.all.append(sec=self.soma)
        self.all.append(sec=self.axon)

    def shape_3D(self):
        """
        Set the default shape of the cell in 3D coordinates.
        Set soma(0) to the origin (0,0,0) and axon extending along
        the X-axis.
        """
        len1 = self.soma.L
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(0, 0, 0, self.soma.diam, sec=self.soma)
        h.pt3dadd(len1, 0, 0, self.soma.diam, sec=self.soma)
        len2 = self.axon.L
        h.pt3dclear(sec=self.axon)
        h.pt3dadd(len1, 0, 0, self.axon.diam, sec=self.axon)
        h.pt3dadd(len1 + len2, 0, 0, self.axon.diam, sec=self.axon)

    def set_position(self, x, y, z):
        """
        Set the base location in 3D and move all other
        parts of the cell relative to that location.
        """
        for sec in self.all:
            # note: iterating like this changes the context for all NEURON
            # functions that depend on a section, so no need to specify sec=
            for i in range(sec.n3d()):
                h.pt3dchange(i,
                             x - self.x + sec.x3d(i),
                             y - self.y + sec.y3d(i),
                             z - self.z + sec.z3d(i),
                             sec.diam3d(i), sec=sec)
        self.x, self.y, self.z = x, y, z

    def set_recording_vectors(self, J, M):
        self.soma_v = h.Vector()
        self.axon_v = h.Vector()
        self.time_vec = h.Vector()
        self.soma_spike_times = h.Vector()
        self.axon_spike_times = h.Vector()

        self.soma_v.record(self.soma(self.soma_dx * (J-1) + (self.soma_dx / 2))._ref_v)
        self.soma_v.record(self.soma(self.axon_dx * (M-1) + (self.axon_dx / 2))._ref_v)
        self.time_vec.record(h._ref_t)
        self.soma_spike_detector = h.NetCon(self.soma(1)._ref_v, None,
                                            sec=self.soma)
        self.soma_spike_detector.threshold = 10
        self.axon_spike_detector = h.NetCon(self.axon(1)._ref_v, None,
                                            sec=self.axon)
        self.axon_spike_detector.threshold = 10
        self.soma_spike_detector.record(self.soma_spike_times)
        self.axon_spike_detector.record(self.axon_spike_times)

    def rotateZ(self, theta):
        """Rotate the cell about the Z axis."""
        rel_origin = (self.soma.x3d(0), self.soma.y3d(0))
        for sec in self.all:
            for i in range(sec.n3d()):
                point = (sec.x3d(i), sec.y3d(i))
                dx = (rel_origin[0] + math.cos(theta)
                      * (point[0] - rel_origin[0]) - math.sin(theta)
                      * (point[1] - rel_origin[1]))
                dy = (rel_origin[1] + math.sin(theta)
                      * (point[0] - rel_origin[0]) + math.cos(theta)
                      * (point[1] - rel_origin[1]))
                new_point = (dx, dy)
                h.pt3dchange(i, dx, dy, 0, sec.diam3d(i), sec=sec)

    def give_dxs(self):
        return (self.soma_dx, self.axon_dx)

    def get_ident(self):
        return self.ident

    def give_synapses(self):
        return self.synapses

    def give_spikes(self):
        return self.soma_spike_times, self.axon_spike_times, self.time_vec

    def give_vec_halves(self):
        return self.soma_v, self.axon_v

    def give_len(self, str):
        if (str == "soma"):
            return self.soma.L
        elif (str == "axon"):
            return self.axon.L
        else:
            return 0

    def add_synapse(self, target_cell, conduct_range):
        # target_cell is a cell object
        if target_cell.get_ident() not in self.connections:
            self.connections.append(target_cell.get_ident())
            syn = Synapse(self, target_cell, conduct_range)
            self.synapses.append(syn)


class RelayCell(object):
    """Two-section cell: A soma with active channels and
    an axon with active channels."""

    def __init__(self, rs_pars, ra_pars, identifier):
        self.ident = identifier
        self.x = self.y = self.z = 0
        self.connections = []  # List of cell identifiers this cell projects to
        self.synapses = []  # List of Synapse Objects corr. to identifiers
        self.create_sections()
        self.build_subsets()
        self.define_geometry(rs_pars, ra_pars)
        self.define_biophysics(rs_pars, ra_pars)
        self.build_topology()
        self.set_recording_vectors(rs_pars[6], ra_pars[6])

    def create_sections(self):
        """Create the sections of the cell."""
        # NOTE: cell=self is required to tell NEURON of this object.
        self.soma = h.Section(name='soma', cell=self)
        self.axon = h.Section(name='axon', cell=self)

    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        self.axon.connect(self.soma(1))

    def define_geometry(self, rs_pars, ra_pars):
        """Set the 3D geometry of the cell."""
        self.soma.L = self.soma.diam = 60  # microns
        self.axon.L = 40  # microns
        self.axon.diam = 7  # microns
        self.soma.nseg = rs_pars[6]
        self.axon.nseg = ra_pars[6]
        self.soma_dx = (1 / rs_pars[6])
        self.axon_dx = (1 / ra_pars[6])
        self.shape_3D()

    # raparams = [ra_EK, ra_ENA, ra_EL, ra_gNa, ra_gK, ra_gL, N]
    def define_biophysics(self, rs_pars, ra_pars):
        """Assign the membrane properties across the cell."""
        for sec in self.all:
            if (sec == self.axon):
                sec.Ra = 500  # Axial resistance in Ohm * cm
            else:
                sec.Ra = 100  # Axial resistance in Ohm * cm
            sec.cm = 1  # Membrane capacitance in micro Farads / cm^2

        # Insert active channels in the soma
        self.soma.insert('hh')
        self.soma.insert('na_ion')
        self.soma.insert('k_ion')
        for seg in self.soma:
            seg.hh.gnabar = rs_pars[3]
            seg.hh.gkbar = rs_pars[4]
            seg.ena = rs_pars[1]
            seg.ek = rs_pars[0]
            seg.hh.gl = rs_pars[5]  # Leak conductance in S/cm2
            seg.hh.el = rs_pars[2]  # Reversal potential in mV

        # Insert active Hodgkin-Huxley current in the axon
        self.axon.insert('hh')
        self.axon.insert('na_ion')
        self.axon.insert('k_ion')
        for seg in self.axon:
            seg.hh.gnabar = ra_pars[3]  # Sodium conductance in S/cm2
            seg.hh.gkbar = ra_pars[4]  # Potassium conductance in S/cm2
            seg.ena = ra_pars[1]
            seg.ek = ra_pars[0]
            seg.hh.gl = ra_pars[5]  # Leak conductance in S/cm2
            seg.hh.el = ra_pars[2]  # Reversal potential in mV

    def build_subsets(self):
        """Build subset lists. For now we define 'all'."""
        self.all = h.SectionList()
        self.all.append(sec=self.soma)
        self.all.append(sec=self.axon)

    def shape_3D(self):
        """
        Set the default shape of the cell in 3D coordinates.
        Set soma(0) to the origin (0,0,0) and axon extending along
        the X-axis.
        """
        len1 = self.soma.L
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(0, 0, 0, self.soma.diam, sec=self.soma)
        h.pt3dadd(len1, 0, 0, self.soma.diam, sec=self.soma)
        len2 = self.axon.L
        h.pt3dclear(sec=self.axon)
        h.pt3dadd(len1, 0, 0, self.axon.diam, sec=self.axon)
        h.pt3dadd(len1 + len2, 0, 0, self.axon.diam, sec=self.axon)

    def set_position(self, x, y, z):
        """
        Set the base location in 3D and move all other
        parts of the cell relative to that location.
        """
        for sec in self.all:
            # note: iterating like this changes the context for all NEURON
            # functions that depend on a section, so no need to specify sec=
            for i in range(sec.n3d()):
                h.pt3dchange(i,
                             x - self.x + sec.x3d(i),
                             y - self.y + sec.y3d(i),
                             z - self.z + sec.z3d(i),
                             sec.diam3d(i), sec=sec)
        self.x, self.y, self.z = x, y, z

    def set_recording_vectors(self, K, N):
        self.soma_v = h.Vector()
        self.axon_v = h.Vector()
        self.time_vec = h.Vector()
        self.soma_spike_times = h.Vector()
        self.axon_spike_times = h.Vector()

        self.soma_v.record(self.soma(self.soma_dx * (K-1) + (self.soma_dx / 2))._ref_v)
        self.soma_v.record(self.soma(self.axon_dx * (N-1) + (self.axon_dx / 2))._ref_v)
        self.time_vec.record(h._ref_t)
        self.soma_spike_detector = h.NetCon(self.soma(1)._ref_v, None,
                                            sec=self.soma)
        self.soma_spike_detector.threshold = 10
        self.axon_spike_detector = h.NetCon(self.axon(1)._ref_v, None,
                                            sec=self.axon)
        self.axon_spike_detector.threshold = 10
        self.soma_spike_detector.record(self.soma_spike_times)
        self.axon_spike_detector.record(self.axon_spike_times)

    def rotateZ(self, theta):
        """Rotate the cell about the Z axis."""
        rel_origin = (self.soma.x3d(0), self.soma.y3d(0))
        for sec in self.all:
            for i in range(sec.n3d()):
                point = (sec.x3d(i), sec.y3d(i))
                dx = (rel_origin[0] + math.cos(theta)
                      * (point[0] - rel_origin[0]) - math.sin(theta)
                      * (point[1] - rel_origin[1]))
                dy = (rel_origin[1] + math.sin(theta)
                      * (point[0] - rel_origin[0]) + math.cos(theta)
                      * (point[1] - rel_origin[1]))
                new_point = (dx, dy)
                h.pt3dchange(i, dx, dy, 0, sec.diam3d(i), sec=sec)

    def give_dxs(self):
        return self.soma_dx, self.axon_dx

    def get_ident(self):
        return self.ident

    def give_synapses(self):
        return self.synapses

    def give_spikes(self):
        return self.soma_spike_times, self.axon_spike_times, self.time_vec

    def give_vec_halves(self):
        return self.soma_v, self.axon_v

    def give_len(self, str):
        if str == "soma":
            return self.soma.L
        elif str == "axon":
            return self.axon.L
        else:
            return 0

    def add_synapse(self, target_cell, conduct_range):
        # target_cell is a cell object
        if target_cell.get_ident() not in self.connections:
            self.connections.append(target_cell.get_ident())
            syn = Synapse(self, target_cell, conduct_range)
            self.synapses.append(syn)


# UTILITY FUNCTIONS
def simulate(tstop):
    """Initialize and run a simulation.

    :param tstop: Duration of the simulation. Currently set to 9.12 to make it
    easier to visually determine a.p. frequency. In this sim. time frame, we
    need 8 full action potentials for a target frequency of 880Hz.
    """
    h.tstop = tstop
    h.run()


def raster(cells, show_fig, dir, X, p_s_f, p_a_f, r_s_f, r_a_f):
    cells = cells
    show_fig = show_fig
    dir = dir
    x = X
    p_s_f = p_s_f
    p_a_f = p_a_f
    r_s_f = r_s_f
    r_a_f = r_a_f

    pyplot.figure()
    ax = pyplot.axes()
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('cell identifier (#)')
    lbl_ps = f"P_soma F={p_s_f} Hz"
    lbl_pa = f"P_axon F={p_a_f} Hz"
    lbl_rs = f"R_soma F={r_s_f} Hz"
    lbl_ra = f"R_axon F={r_a_f} Hz"
    label_set1 = False
    label_set2 = False

    for i, cell in enumerate(cells):
        soma_v, axon_v, t_v = cell.give_spikes()
        if i < 87:
            soma = 'black'
            axon = 'red'
            if label_set1:
                pyplot.vlines(soma_v, i + 0.3, i + 1.3, colors=soma, alpha=0.6)
                pyplot.vlines(axon_v, i + 0.3, i + 1.3, colors=axon, alpha=0.6)
            else:
                pyplot.vlines(soma_v, i + 0.3, i + 1.3, colors=soma, alpha=0.6, label=lbl_ps)
                pyplot.vlines(axon_v, i + 0.3, i + 1.3, colors=axon, alpha=0.6, label=lbl_pa)
                label_set1 = True
        else:
            soma = 'blue'
            axon = 'green'
            if label_set2:
                pyplot.vlines(soma_v, i + 0.3, i + 1.3, colors=soma, alpha=0.6)
                pyplot.vlines(axon_v, i + 0.3, i + 1.3, colors=axon, alpha=0.6)
            else:
                pyplot.vlines(soma_v, i + 0.3, i + 1.3, colors=soma, alpha=0.6, label=lbl_rs)
                pyplot.vlines(axon_v, i + 0.3, i + 1.3, colors=axon, alpha=0.6, label=lbl_ra)
                label_set2 = True
    title = f"{x}_Raster"
    ax.set_title(title)
    ax.legend(loc='upper right', bbox_to_anchor=(0.9, 0.7))
    pyplot.savefig(f"{dir}/{title}.pdf", format='pdf')
    pyplot.close()


def frequency_vs_EK(eks, freqs, conf):
    pts = zip(freqs, eks)
    scatter = pyplot.figure()
    ax = scatter.gca()
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlabel('EK (mV)')
    for pt in pts:
        ax.scatter(pt[1], pt[0], c="blue")
    ax.set_title(f"{conf}_Network_EK_vs_F")
    pyplot.savefig(f"{conf}_Network_EK_vs_F.eps", format='eps')
    pyplot.show()
    pyplot.close()


def cellular_potentials(cells, dir, X, p_f, r_f):
    """ List of Cell -> Figure of Cellular Potentials vs Time vs Cell #.

    Receives a LoCell and graphs their membrane potential vs time vs Cell #.
    """
    cells = cells
    dir = dir
    x = X
    p_f = p_f
    r_f = r_f
    mem_potentials_fig = pyplot.figure()
    ax = pyplot.axes(projection='3d')
    ax.set_xlabel('Simulation Time (ms)')
    ax.set_ylabel('Cell Identifier (#)')
    ax.set_zlabel('Membrane Potential (mV)')
    lbl_p = f"P_cell F={p_f} Hz"
    lbl_r = f"R_cell F={r_f} Hz"

    for i in range(0, len(cells)):
        cell = cells[i]
        ident = cell.get_ident()
        soma_vecs, axon_vecs, time_vec = cell.give_vectors()
        s_half, a_half = cell.give_vec_halves()
        if i < 87:
            color = 'red'
            if i == 0:
                ax.plot3D(time_vec, i * np.ones(len(time_vec)), s_half, c=color, alpha=.4, label=lbl_p)
            else:
                ax.plot3D(time_vec, i * np.ones(len(time_vec)), s_half, c=color, alpha=.4)
        else:
            color = 'green'
            if i == 87:
                ax.plot3D(time_vec, i * np.ones(len(time_vec)), s_half, c=color, alpha=.4, label=lbl_r)
            else:
                ax.plot3D(time_vec, i * np.ones(len(time_vec)), s_half, c=color, alpha=.4)

    title = f"{x}_Mem_Pots"
    ax.set_title(title)
    ax.legend(loc='upper right', bbox_to_anchor=(0.9, 0.7))
    pyplot.savefig(f"{dir}/{title}.pdf", format='pdf')
    pyplot.close()


def post_sim_analysis(path):
    # Loop through all dirs in current dir.
    ct = 0
    for subdir in next(os.walk(path))[1]:
        if subdir.endswith("]"):
            print(f"Found directory {subdir}")
            ct += 1
            with open(f"{path}/{subdir}/{subdir}_spiketimes.txt") as f:
                labels = []
                t_vec = []
                somas = []
                axons = []
                soma_volts = []
                axon_volts = []
                # Parse and store data
                for line in f:
                    line = ast.literal_eval(line.rstrip('\n'))
                    case = line[0]
                    if (case == 'time'):
                        t_vec = line[1]
                    elif (case == 'soma'):
                        somas.append([line[1], line[2]])
                    elif (case == 'axon'):
                        axons.append([line[1], line[2]])
                    elif (case == 'freqs'):
                        labels.append(f"P_soma F={line[1]} Hz")
                        labels.append(f"P_axon F={line[2]} Hz")
                        labels.append(f"R_soma F={line[3]} Hz")
                        labels.append(f"R_axon F={line[4]} Hz")
                    elif (case == 'volt'):
                        soma_volts.append([line[1], line[2]])
                        axon_volts.append([line[1], line[3]])

                if len(labels) == 0:
                    continue

                # Plot raster data
                raster = pyplot.figure()
                ax = pyplot.axes()
                ax.set_xlabel('time (ms)')
                ax.set_ylabel('cell identifier (#)')
                for soma, axon in zip(somas, axons):
                    idents = [soma[0], axon[0]]
                    vecs = [soma[1], axon[1]]

                    for idx, (idt, vec) in enumerate(zip(idents, vecs)):
                        if idx == 0:  # soma
                            pace_c = "black"
                            relay_c = "blue"
                            lbls = [labels[0], labels[2]]
                        else:  # axon
                            pace_c = "red"
                            relay_c = "green"
                            lbls = [labels[1], labels[3]]
                        if idt == 0:
                            pyplot.vlines(vec, idt + 0.3, idt + 1.3, colors=pace_c, alpha=0.6, label=lbls[0])
                        elif idt == 87:
                            pyplot.vlines(vec, idt + 0.3, idt + 1.3, colors=relay_c, alpha=0.6, label=lbls[1])
                        else:
                            if idt < 87:
                                pyplot.vlines(vec, idt + 0.3, idt + 1.3, colors=pace_c, alpha=0.6)
                            else:
                                pyplot.vlines(vec, idt + 0.3, idt + 1.3, colors=relay_c, alpha=0.6)
                string = subdir.split("=", 1)[0]
                title = f"{string}_Raster"
                ax.set_title(title)
                ax.legend(loc='upper right', bbox_to_anchor=(0.9, 0.7))
                pyplot.savefig(f"{subdir}/{title}.pdf", format='pdf')
                pyplot.close()
                print(f"Raster processed for {subdir}.")

                # Plot membrane potential data
                mem_potentials = pyplot.figure()
                ax = pyplot.axes(projection='3d')
                ax.set_xlabel('Simulation Time (ms)')
                ax.set_ylabel('Cell Identifier (#)')
                ax.set_zlabel('Membrane Potential (mV)')
                for soma_v, axon_v in zip(soma_volts, axon_volts):
                    idents = [soma_v[0], axon_v[0]]
                    vecs = [soma_v[1], axon_v[1]]

                    for idx, (idt, vec) in enumerate(zip(idents, vecs)):
                        if idx == 0:  # soma
                            pace_c = "black"
                            relay_c = "blue"
                            lbls = [labels[0], labels[2]]
                        else:  # axon
                            pace_c = "red"
                            relay_c = "green"
                            lbls = [labels[1], labels[3]]
                        if idt == 0:
                            ax.plot3D(t_vec, idt * np.ones(len(t_vec)), vec, c=pace_c, alpha=.4, label=lbls[0])
                        elif idt == 87:
                            ax.plot3D(t_vec, idt * np.ones(len(t_vec)), vec, c=relay_c, alpha=.4, label=lbls[1])
                        else:
                            if idt < 87:
                                ax.plot3D(t_vec, idt * np.ones(len(t_vec)), vec, c=pace_c, alpha=.4)
                            else:
                                ax.plot3D(t_vec, idt * np.ones(len(t_vec)), vec, c=relay_c, alpha=.4)
                string = subdir.split("=", 1)[0]
                title = f"{string}_Mem_Pots"
                ax.set_title(title)
                ax.legend(loc='upper right', bbox_to_anchor=(0.9, 0.7))
                pyplot.savefig(f"{subdir}/{title}.pdf", format='pdf')
                pyplot.close()
                print(f"Membrane potentials plot processed for {subdir}.")
    print(f"Found {ct} matching directories and processed their raw data.")


def attach_current_clamp(cell, delay, dur, amp=1, loc=0, freq=1):
    """Attach a current Clamp to a cell.

    :param cell: Cell object to attach the current clamp.
    :param delay: Onset of the injected current.
    :param dur: Duration of the stimulus.
    :param amp: Magnitude of the current.
    :param loc: Location on the soma where the stimulus is placed.
    """
    stim = h.IClamp(cell.soma(loc))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp
    return stim
