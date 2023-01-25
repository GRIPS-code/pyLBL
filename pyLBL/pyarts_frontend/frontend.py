from logging import getLogger
from os.path import join
from pathlib import Path
from subprocess import run
from sys import stderr, stdout

from numpy import asarray

import pyarts


info = getLogger(__name__).info


def absorption_line(line):
    """ Translates a single pyLBL Transition object to ARTS AbsorptionSingleLine
    
    Input
    -----
        line : Transition
            A pyLBL transition object
    
    Output
    ------
        q_key : QuantumIdentifier
            The pyarts ID of the absorption species
        ls : AbsorptionSingleLine
            A single ARTS absorption line
    """
    iso = line.local_iso_id
    if iso == 11:
        iso = 'A'
    elif iso == 12:
        iso = 'B'
    iso = str(iso)
    
    r = pyarts.arts.hitran.ratio(line.molecule_id, iso)
    qkey = pyarts.arts.hitran.quantumidentity(line.molecule_id, iso)

    slf = pyarts.arts.LineShapeSingleSpeciesModel(
        G0=pyarts.arts.LineShapeModelParameters("T1",
                                                pyarts.arts.convert.kaycm_per_atm2hz_per_pa(
                                                    line.gamma_self),
                                                line.n_air),
        D0=pyarts.arts.LineShapeModelParameters("T0",
                                                pyarts.arts.convert.kaycm_per_atm2hz_per_pa(
                                                    line.delta_air)
                                                ))

    air = pyarts.arts.LineShapeSingleSpeciesModel(
        G0=pyarts.arts.LineShapeModelParameters("T1",
                                                pyarts.arts.convert.kaycm_per_atm2hz_per_pa(
                                                    line.gamma_air),
                                                line.n_air),
        D0=pyarts.arts.LineShapeModelParameters("T0",
                                                pyarts.arts.convert.kaycm_per_atm2hz_per_pa(
                                                    line.delta_air)
                                                ))

    sl = pyarts.arts.AbsorptionSingleLine(
        F0=pyarts.arts.convert.kaycm2freq(line.nu),
        I0=pyarts.arts.convert.kaycm_per_cmsquared2hz_per_msquared(line.sw / r),
        E0=pyarts.arts.convert.kaycm2joule(line.elower),
        lineshape=pyarts.arts.LineShapeModel([slf, air]))
    
    return qkey, sl


def absorption_lines(lines):
    """ Translates a list of pyLBL Transition object to ARTS ArrayOfAbsorptionLines
    
    Input
    -----
        lines : List[Transition]
            A list of the database entries for all the lines
    
    Output
    ------
        lines : ArrayOfAbsorptionLines
            As pyarts abs_lines
        
    """
    data = {}

    for line in lines:
        qkey, sl = absorption_line(line)
        key = str(qkey)
        if key in data:
            data[key].append(sl)
        else:
            data[key] = [sl]

    aal = pyarts.arts.ArrayOfAbsorptionLines()
    for x in data:
        aal.append(
            pyarts.arts.AbsorptionLines(
                selfbroadening=True,
                bathbroadening=True,
                cutoff="None",
                mirroring="None",
                population="LTE",
                normalization="SFS",
                lineshapetype="SplitVP",
                quantumidentity=x,
                broadeningspecies=[x.split('-')[0], "Bath"],
                T0=296,
                lines=data[x]
            ))
    return aal

class PyArtsGas(object):
    def __init__(self, lines_database, formula):
        self.ws = pyarts.workspace.Workspace()
        self.ws.abs_speciesSet(species=[formula])
        self.ws.abs_lines_per_species = [absorption_lines(lines_database.gas(formula)[2])]
        
        self.ws.jacobianOff()
        self.ws.Touch(self.ws.rtp_nlte)
        self.ws.Touch(self.ws.rtp_mag)
        self.ws.Touch(self.ws.rtp_los)
        self.ws.propmat_clearsky_agendaAuto()
        self.ws.lbl_checkedCalc()
        self.ws.stokes_dim = 1

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               remove_pedestal=False, cut_off=25):
        """Calculates absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber [cm-1].

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        # Configure spectral grid.
        self.ws.f_grid = pyarts.arts.convert.kaycm2freq(grid)

        # Configure the atmosphere.
        self.ws.rtp_pressure = pressure
        self.ws.rtp_temperature = temperature
        self.ws.rtp_vmr = [volume_mixing_ratio]

        # Calculate the absorption coefficient.
        self.ws.AgendaExecute(a=self.ws.propmat_clearsky_agenda)
        x = pyarts.arts.physics.number_density(pressure, temperature) * volume_mixing_ratio
        return self.ws.propmat_clearsky.value.data.value.flatten() / x
