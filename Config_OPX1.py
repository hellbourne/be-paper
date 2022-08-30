from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import numpy as np
# from qm import SimulationConfig, ConfigHelper
from qm import LoopbackInterface
import matplotlib.pyplot as plt
import time


IF_freq = 150e6
rad_freq = 200e6
LO_freq = 1.5e9-IF_freq

config = {

    'version': 1,

    'controllers': {

        'finkler1': {
            'type': 'opx1',
            'analog_outputs': {
                1: {'offset': 0.0},  # I
                2: {'offset': 0.0},  # Q
                3: {'offset': 0.01407},  # I
                4: {'offset': 0.01286},  # Q
            },
            'digital_outputs': {
                1: {'offset': 0.0},
                2: {'offset': 0.0},
                7: {'offset': 0.0},
                8: {'offset': 0.0},
            },  # laser
            'analog_inputs': {
                1: {'offset': -0.032},  # readout1
                2: {'offset': -0.037},  # readout2

            }
        }
    },

    'elements': {

        'qubit': {
            'mixInputs': {
                'I': ('finkler1', 1),
                'Q': ('finkler1', 2),
                'lo_frequency': LO_freq,
                'mixer': 'mixer_qubit'
            },
            'intermediate_frequency': IF_freq,
            'digitalInputs': {
                "ttl_in": {
                    "port": ("finkler1", 7),
                    "delay": 136,
                    "buffer": 15,
                },
                "ttl_check": {
                    "port": ("finkler1", 8),
                    "delay": 136,
                    "buffer": 15,
                },
            },
            'operations': {
                'MW_pi': 'MW_pi_pulse',
            }
        },
        'qubit2': {
            'mixInputs': {
                'I': ('finkler1', 3),
                'Q': ('finkler1', 4),
                'lo_frequency': LO_freq,
                'mixer': 'mixer_qubit2'
            },
            'intermediate_frequency': IF_freq,
            'operations': {
                'MW_pi': 'MW_pi_pulse2',
            }
        },

        'radical': {
            'mixInputs': {
                'I': ('finkler1', 1),
                'Q': ('finkler1', 2),
                'lo_frequency': LO_freq,
                'mixer': 'mixer_qubit'
            },
            'intermediate_frequency': rad_freq,
            'digitalInputs': {
                "ttl_in": {
                    "port": ("finkler1", 7),
                    "delay": 136,
                    "buffer": 15,
                },
            },
            'operations': {
                'MW_pi': 'MW_pi_pulse',
            }
        },


        "laser": {
            'digitalInputs': {
                "laser_in": {
                    "port": ("finkler1", 1),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'laser_pulse': 'laser_trig_pulse',
            }
        },

        'readout1': {

            # open: fake
            "singleInput": {
                "port": ("finkler1", 1)
            },


            'operations': {
                'readout': 'readout_pulse',
            },
            "outputs": {
                'out1': ("finkler1", 1)
            },
            'time_of_flight': 28,
            'smearing': 0,
            
            'outputPulseParameters': {
                'signalThreshold': -300,
                'signalPolarity': 'Descending',
                'derivativeThreshold': -50,
                'derivativePolarity': 'Descending'
            }
        },

         'readout2': {

             # open: fake
             "singleInput": {
                 "port": ("finkler1", 1)
             },


             'operations': {
                 'readout': 'readout_pulse',
             },
             "outputs": {
                 'out1': ("finkler1", 2)
             },
             'time_of_flight': 28,
             'smearing': 0,
             
             'outputPulseParameters': {
                'signalThreshold': -300,
                'signalPolarity': 'Descending',
                'derivativeThreshold': -50,
                'derivativePolarity': 'Descending'}
         },

        "NI": {
            'digitalInputs': {
                "NI_in": {
                    "port": ("finkler1", 2),
                    "delay": 0,
                    "buffer": 0,
                },
            },
            'operations': {
                'follow': 'follow_pulse',
            }
        },
    },

    "pulses": {
        'MW_pi_pulse': {
            'operation': "control",
            'length': 24*4,
            'waveforms': {
                "I": "const_wf",
                "Q": "zero_wf"
            },
            'digital_marker': 'ON'
        },

        'MW_pi_pulse2': {
            'operation': "control",
            'length': 24*4,
            'waveforms': {
                "I": "const_wf",
                "Q": "zero_wf"
            },
        },

        'Radical_pi_pulse': {
            'operation': 'control',
            'length': 24*4,
            'waveforms': {
                'single': 'const_radical_wf',
            },
        },


        'laser_trig_pulse': {
            'operation': "control",
            'length': 100,
            'digital_marker': 'ON'
        },

        'follow_pulse': {
            'operation': "control",
            'length': 80,
            'digital_marker': 'ON'
        },

        'readout_pulse': {
            'operation': "measurement",
            "length": 2000,
            "waveforms": {
                "single": "zero_wf"  # fake!
            },
            'digital_marker': 'ON',
        },
    },

    'waveforms': {
        "const_wf": {
            "type": "constant",  # arbitrary give list og amps
            "sample": 0.4  # need to give samples
        },

        "zero_wf": {
            "type": "constant",
            "sample": 0.0
        },

        'const_radical_wf': {
            "type": "constant",  # arbitrary give list og amps
            "sample": 0.49  # need to give samples
        },
    },

    "digital_waveforms": {
        "ON": {
            "samples": [(1, 0)]  # (on/off, ns)
        }
    },

    "mixers": {
        'mixer_qubit2': [
            
            {'intermediate_frequency': rad_freq, 'lo_frequency': LO_freq, 'correction': [1.18,  0.19,  0.24,  0.96]},  # SSB (Both (-) left (+) right)
            {'intermediate_frequency': IF_freq, 'lo_frequency': LO_freq, 'correction': [1.18,  0.19,  0.24,  0.96]}  # SSB (Both (-) left (+) right)
        ],
        'mixer_qubit': [
            {'intermediate_frequency': rad_freq, 'lo_frequency': LO_freq, 'correction': [1.18,  0.19,  0.24,  0.96]},  # SSB (Both (-) left (+) right)
            {'intermediate_frequency': IF_freq, 'lo_frequency': LO_freq, 'correction': [1.0,  0.0,  0.0,  1.0]}  # SSB (Both (-) left (+) right)

        ]
    },
}
